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

    ``BayesianNormalizingFlow``
        Wrapper class for a PyTorch/Zuko-based Bayesian MAF model

    ``EnsembleNormalizingFlow``
        Wrapper class for a collection of several PyTorch/Zuko-based MAF models
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats


def to_latent_space(x: np.ndarray, bounds: list[tuple[float, float]]) -> np.ndarray:
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

    # To prevent transforming values at the bounds to -inf
    eps = 1e-4  # (forces the distance of x from any of its bounds to be >= 1e-4)

    transformed = []
    for dim in range(x.shape[-1]):
        inf_bounds = np.logical_not(np.isfinite(bounds[dim]))

        # already unbounded, no transformation
        if np.all(inf_bounds):
            z = x[:, dim]

        # upper bound is inf, bounded below only
        elif inf_bounds[1]:
            z = np.log(x[:, dim] - bounds[dim][0] + eps)

        # lower bound is -inf, bounded above only
        elif inf_bounds[0]:
            z = np.log(bounds[dim][1] - x[:, dim] + eps)

        # no inf bounds, bounded above and below
        else:
            a = (x[:, dim] - bounds[dim][0]) / (bounds[dim][1] - bounds[dim][0])
            z = np.log(a / (1 - a) + eps)

        transformed.append(z)

    return np.stack(transformed, axis=-1)


def to_data_space(z: np.ndarray, bounds: list[tuple[float, float]]) -> np.ndarray:
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


def get_log_abs_det_jacobian(x: np.ndarray, bounds: list[tuple[float, float]]) -> np.ndarray:
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

    # To prevent Jacobian determinant of +inf
    eps = 1e-4  # (forces the distance of x from one of its bounds to be >= 1e-4)

    ladj = np.zeros_like(x[:, 0])

    for dim in range(x.shape[-1]):
        inf_bounds = np.logical_not(np.isfinite(bounds[dim]))

        # already unbounded, no transformation jacobian
        if np.all(inf_bounds):
            continue

        # upper bound is inf, bounded below only
        elif inf_bounds[1]:
            ladj += -np.log(x[:, dim] - bounds[dim][0] + eps)

        # lower bound is -inf, bounded above only
        elif inf_bounds[0]:
            ladj += -np.log(bounds[dim][1] - x[:, dim] + eps)

        # no inf bounds, bounded above and below
        else:
            ladj += np.log(bounds[dim][1] - bounds[dim][0]) - np.log(x[:, dim] - bounds[dim][0] + eps) - np.log(bounds[dim][1] - x[:, dim] + eps)

    return ladj


class BoundedKDE:
    """Kernel Density Estimator, applicable to bounded data using transformation."""

    def __init__(self, posterior_samples: np.ndarray, bounds: list[tuple[float, float]]):
        """
        Parameters
        ----------
        posterior_samples
            Array of samples (with shape (N_samples, N_dim)) from the distribution on which density estimation is being performed
        bounds
            List of tuples of the form (lower bound, upper bound) for each of the N_dim parameters appearing in posterior_samples.
            For infinite bounds (unbounded), use -np.inf or np.inf.
        """

        self.posterior_samples = posterior_samples
        self.bounds = bounds

        Z = to_latent_space(self.posterior_samples, self.bounds)

        # Scipy KDE fails if there are inf or -inf values
        Z = np.nan_to_num(Z, posinf=1e4, neginf=-1e4)

        self.base_kde = scipy.stats.gaussian_kde(Z.T)

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

        # Scipy KDE fails if there are inf or -inf values
        z = np.nan_to_num(z, posinf=1e4, neginf=-1e4)

        lp = kde.logpdf(z.T).T + ladj
        lp = np.nan_to_num(lp, nan=-np.inf)
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
