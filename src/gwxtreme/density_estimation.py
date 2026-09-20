# Copyright (C) 2026 Joseph David Quinn-Vitabile
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""Flexible and general wrapper classes for various probability density estimator implementations.

These implementations allow enforcement of boundary conditions on the data space by transforming data to unbounded space.

Implemented classes:
    ``BoundedKDE``
        Wrapper class for a Scipy Gaussian kernel density estimator

    ``EnsembleNormalizingFlow``
        Wrapper class for a collection of several PyTorch/Zuko-based MAF models
"""

import os
import pathlib
import threading
from collections.abc import Sequence

import matplotlib.pyplot as plt
import numpy as np
import onnxruntime
import scipy.stats
import torch
from zuko.flows import MAF

# Process-wide cache of ONNX inference sessions, keyed by
# (pid, flow files, intra_op_num_threads). See ``_get_onnx_sessions``.
_onnx_session_cache: dict[tuple, list[onnxruntime.InferenceSession]] = {}
_onnx_session_cache_lock = threading.Lock()


def _get_onnx_sessions(flow_files: tuple[str, ...], intra_op_num_threads: int = 1) -> list[onnxruntime.InferenceSession]:
    """Get the ONNX inference sessions for the given flow files, creating them if this
    process has not already done so.

    ``onnxruntime.InferenceSession`` objects cannot be pickled, so they are held here
    rather than on the ``EnsembleNormalizingFlow`` instances that use them.

    The cache key includes the process ID because inference sessions are not fork safe:
    a child process that inherits this cache through ``fork`` must build its own
    sessions rather than reuse the parent's.

    Parameters
    ----------
    flow_files
        Paths to the .onnx files of the ensemble
    intra_op_num_threads
        Number of threads each session may use to parallelize the execution of an
        individual operator

    Returns
    -------
        List of inference sessions, one per file in ``flow_files``
    """

    key = (os.getpid(), flow_files, intra_op_num_threads)

    with _onnx_session_cache_lock:
        sessions = _onnx_session_cache.get(key)

        if sessions is None:
            options = onnxruntime.SessionOptions()
            options.intra_op_num_threads = intra_op_num_threads
            options.inter_op_num_threads = 1
            providers = ["CPUExecutionProvider"]

            sessions = [onnxruntime.InferenceSession(flow_file, options, providers=providers) for flow_file in flow_files]
            _onnx_session_cache[key] = sessions

        return sessions


def to_latent_space(x: np.ndarray, bounds: Sequence[tuple[float, float]]) -> np.ndarray:
    """Transform points to the (unbounded) latent space used by the underlying density estimator.

    This transformation operates on each parameter of the data (iterating over the second dimension
    of x) to compute the latent representation z of each parameter as follows:

        If x is unbounded,              z = x
        If x has a lower bound x_min,   z = log(x - x_min)
        If x has an upper bound x_max,  z = log(x_max - x)
        If x has both bounds,           z = logit( (x - x_min) / (x_max - x_min) )

    Parameters
    ----------
    x
        Array of points in the data space with shape (N_points, N_dim)

    Returns
    -------
        Array of points in the latent space with shape (N_points, N_dim)
    """
    assert x.ndim == 2, "x should be an (N, D) shaped array"

    transformed = []
    with np.errstate(divide="ignore"):
        for dim in range(x.shape[-1]):
            inf_bounds = np.logical_not(np.isfinite(bounds[dim]))

            # already unbounded, no transformation
            if np.all(inf_bounds):
                z = x[:, dim]

            # upper bound is inf, bounded below only
            elif inf_bounds[1]:
                z = np.log(x[:, dim] - bounds[dim][0])

            # lower bound is -inf, bounded above only
            elif inf_bounds[0]:
                z = np.log(bounds[dim][1] - x[:, dim])

            # no inf bounds, bounded above and below
            else:
                a = (x[:, dim] - bounds[dim][0]) / (bounds[dim][1] - bounds[dim][0])
                z = np.log(a / (1 - a))

            transformed.append(z)

    return np.stack(transformed, axis=-1)


def to_data_space(z: np.ndarray, bounds: Sequence[tuple[float, float]]) -> np.ndarray:
    """Transform points to the (bounded) data space.

    This transformation inverts the transformation done by _to_latent_space.

    Parameters
    ----------
    z
        Array of points in the latent space with shape (N_points, N_dim)

    Returns
    -------
        Array of points in the data space with shape (N_points, N_dim)
    """
    assert z.ndim == 2, "z should be an (N, D) shaped array"

    transformed = []
    for dim in range(z.shape[-1]):
        inf_bounds = np.logical_not(np.isfinite(bounds[dim]))

        # already unbounded, no transformation
        if np.all(inf_bounds):
            x = z[:, dim]

        # upper bound is inf, bounded below only
        elif inf_bounds[1]:
            x = np.exp(z[:, dim]) + bounds[dim][0]

        # lower bound is -inf, bounded above only
        elif inf_bounds[0]:
            x = bounds[dim][1] - np.exp(z[:, dim])

        # no inf bounds, bounded above and below
        else:
            a = z[:, dim] * (bounds[dim][1] - bounds[dim][0]) + bounds[dim][0]
            x = 1 / (1 + np.exp(-a))

        transformed.append(x)

    return np.stack(transformed, axis=-1)


def get_log_abs_det_jacobian(x: np.ndarray, bounds: Sequence[tuple[float, float]]) -> np.ndarray:
    """Compute the log of the absolute value of the determinant of the jacobian (ladj) of the
    data-to-latent-space transformation.

    The returned ladj values are the sum of those for individual dimensions of x (since the
    Jacobian matrix for the data transformation is diagonal), which are each computed as follows:

        If x is unbounded,              ladj = 0
        If x has a lower bound x_min,   ladj = -log(x - x_min)
        If x has an upper bound x_max,  ladj = -log(x_max - x)
        If x has both bounds,           ladj = log(x_max - x_min) - log(x - x_min) - log(x_max - x)

    Parameters
    ----------
    x
        Array of points in the data space with shape (N_points, N_dim)

    Returns
    -------
        Array of ladj values with shape (N_points,) corresponding to the given points x
    """

    assert x.ndim == 2, "x should be an (N, D) shaped array"

    ladj = np.zeros_like(x[:, 0])

    with np.errstate(divide="ignore"):
        for dim in range(x.shape[-1]):
            inf_bounds = np.logical_not(np.isfinite(bounds[dim]))

            # already unbounded, no transformation jacobian
            if np.all(inf_bounds):
                continue

            # upper bound is inf, bounded below only
            elif inf_bounds[1]:
                ladj += -np.log(x[:, dim] - bounds[dim][0])

            # lower bound is -inf, bounded above only
            elif inf_bounds[0]:
                ladj += -np.log(bounds[dim][1] - x[:, dim])

            # no inf bounds, bounded above and below
            else:
                ladj += np.log(bounds[dim][1] - bounds[dim][0]) - np.log(x[:, dim] - bounds[dim][0]) - np.log(bounds[dim][1] - x[:, dim])

    return ladj


class BoundedKDE:
    """Kernel Density Estimator, applicable to bounded data using transformation."""

    def __init__(self, posterior_samples: np.ndarray, bounds: Sequence[tuple[float, float]], downsample: bool = True):
        """
        Parameters
        ----------
        posterior_samples
            Array of samples (with shape (N_samples, N_dim)) from the distribution on which density estimation is being performed
        bounds
            List of tuples of the form (lower bound, upper bound) for each of the N_dim parameters appearing in posterior_samples.
            For infinite bounds (unbounded), use -np.inf or np.inf.
        downsample
            If true, and if the number of posterior samples N is > 25,000, downsample the data by keeping every floor(N/25,000)-th
            sample.
        """

        if downsample and len(posterior_samples) > 25000:
            keep_every = int(np.floor(len(posterior_samples) / 25000))
            self.posterior_samples = posterior_samples[::keep_every]
        else:
            self.posterior_samples = posterior_samples

        self.bounds = bounds

        Z = to_latent_space(self.posterior_samples, self.bounds)
        finite_valued_points = np.all(np.isfinite(Z), axis=1)

        self.base_kde = scipy.stats.gaussian_kde(Z[finite_valued_points].T)

    def log_pdf(self, x: np.ndarray, resample: bool = False) -> np.ndarray:
        """Compute log probability densities of individual points.

        Parameters
        ----------
        x
            Array of points in the data space with shape (N_points, N_dim)
        resample
            Whether to resample the KDE and use the resampled instance for this computation.
            Resampling the KDE means drawing a sample from it (with size equal to the original
            data size) and instantiating a new KDE from this sample.
            Default: False

        Returns
        -------
            Array of log probability densities with shape (N_points,).
        """

        if resample:
            Z = self.base_kde.resample(len(self.posterior_samples))
            kde = scipy.stats.gaussian_kde(Z)
        else:
            kde = self.base_kde

        z = to_latent_space(x, self.bounds)
        ladj = get_log_abs_det_jacobian(x, self.bounds)

        finite_valued_points = np.all(np.isfinite(z), axis=1)

        lp = np.empty(x.shape[0], dtype=np.float64)
        lp[finite_valued_points] = kde.logpdf(z[finite_valued_points].T).T + ladj[finite_valued_points]
        lp[~finite_valued_points] = -np.inf
        return lp

    def pdf(self, x: np.ndarray, resample: bool = False) -> np.ndarray:
        """Compute probability densities of individual points.

        Parameters
        ----------
        x
            Array of points in the data space with shape (N_points, N_dim)
        resample
            Whether to resample the KDE and use the resampled instance for this computation.
            Resampling the KDE means drawing a sample from it (with size equal to the original
            data size) and instantiating a new KDE from this sample.
            Default: False

        Returns
        -------
            Array of probability densities with shape (N_points,)
        """

        return np.exp(self.log_pdf(x, resample=resample))


class EnsembleNormalizingFlow:
    """Density estimator constructed from an ensemble of trained normalizing flow models."""

    def __init__(
        self,
        flows_dir: str,
        bounds: Sequence[tuple[float, float]],
    ):
        """
        Parameters
        ----------
        flows_dir
            Path to a directory containing ONNX files for the flow ensemble
        bounds
            List of tuples of the form (lower bound, upper bound) for each of the N_dim parameters appearing in posterior_samples.
            For infinite bounds (unbounded), use -np.inf or np.inf.

        Notes
        -----
        The ONNX inference sessions are created lazily on the first density evaluation and cached per process
        (see ``_get_onnx_sessions``), so instances of this class can be pickled and sent to worker processes.
        """

        self.bounds = bounds

        # Sorted so that the ensemble order is the same in every process
        self.flow_files = []
        for file in sorted(pathlib.Path(flows_dir).iterdir()):
            if file.suffix == ".onnx":
                self.flow_files.append(str(file))
        self.flow_files = tuple(self.flow_files)

        # Fail now rather than on first use
        missing = [flow_file for flow_file in self.flow_files if not os.path.isfile(flow_file)]
        if missing:
            raise FileNotFoundError(f"Missing flow model files in {flows_dir}: {missing}")

    @property
    def sessions(self) -> list[onnxruntime.InferenceSession]:
        """ONNX inference sessions for the flows in the ensemble, created in this process on first access."""
        return _get_onnx_sessions(self.flow_files)

    def log_pdf(self, x: np.ndarray) -> list[np.ndarray]:
        """Compute log probability densities of individual points.

        Parameters
        ----------
        x
            Array of points in the data space with shape (N_points, N_dim)

        Returns
        -------
            List of arrays, one for each flow in the ensemble. Each array contains the log densities for that model with shape (N_points,).
        """

        z = to_latent_space(x, self.bounds)
        ladj = get_log_abs_det_jacobian(x, self.bounds)

        finite_valued_points = np.all(np.isfinite(z), axis=1)

        z = np.array(z, dtype=np.float32)

        lps = []
        for flow in self.sessions:
            lp: np.ndarray = flow.run([flow.get_outputs()[0].name], {flow.get_inputs()[0].name: z})[0] + ladj
            lp[~finite_valued_points] = -np.inf
            lps.append(lp)

        return lps

    def pdf(self, x: np.ndarray) -> list[np.ndarray]:
        """Compute probability densities of individual points.

        Parameters
        ----------
        x
            Array of points in the data space with shape (N_points, N_dim)

        Returns
        -------
            List of arrays, one for each flow in the ensemble. Each array contains the densities for that model with shape (N_points,).
        """

        p = [np.exp(lp) for lp in self.log_pdf(x)]
        return p
