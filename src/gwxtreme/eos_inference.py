# Copyright (C) 2026 Joseph David Quinn-Vitabile
# Copyright (C) 2022 Shaon Ghosh, Michael Camilo, Xiaoshu Liu
# Copyright (C) 2021 Anarya Ray
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

"""Neutron star equation of state inference algorithms.

This module implements the following classes:
- ``ModelSelector``
- ``JointModelSelector``
- ``ParameterizedEoSSampler``

``ModelSelector`` is used for single-event inference, and can be
directly used to compute Bayes factors between EOS models.

``JointModelSelector`` utilizes several ``ModelSelector`` s to
perform joint inference across several events.

``ParameterizedEoSSampler`` is used to infer the posterior
distribution of the parameters of a parameterized EOS model,
such as the 4-parameter spectral decomposition model, by
running MCMC stochastic sampling. The likelihood function used
for this sampling is the joint evidence for the proposed EOS,
which comes from ``JointModelSelector.compute_joint_eos_evidence_ratio``.
The parameterized inference can therefore be done across multiple
events.
"""

import json
import logging
import os
import pathlib
from typing import Literal

import emcee
import h5py
import numpy as np
import ray

from gwxtreme.density_estimation import BoundedKDE
from gwxtreme.utils import get_q_range, read_prior_or_posterior_file

try:
    from gwxtreme.flow_density_estimation import BayesianNormalizingFlow, EnsembleNormalizingFlow
except ImportError:
    BayesianNormalizingFlow = None  # type: ignore
    EnsembleNormalizingFlow = None  # type: ignore

from gwxtreme.eos_interpolator import EOSInterpolator, convert_masses
from gwxtreme.eos_prior import is_valid_eos
from gwxtreme.utils import (
    get_gw_event_pe_posterior_samples,
    get_mean_mchirp_for_cbc_event,
    get_nicer_pulsar_pe_posterior_samples,
)

logger = logging.getLogger(__name__)


@ray.remote
def _distributed_eos_evidence_integration(
    density_estimator: BayesianNormalizingFlow | BoundedKDE, points: np.ndarray, event_type: Literal["gw-2d", "gw-3d", "gw-4d", "psr"]
):
    if event_type == "gw-4d":
        m1_arr = points[0, :, 0]
        m2_arr = points[:, 0, 1]

        n_grid = points.shape[0]
        points = points.reshape(n_grid**2, 4)

        density = density_estimator.pdf(points, resample=True)
        density = density.reshape((n_grid, n_grid))

        # Integrate over m1 and m2
        integral_over_m2 = np.trapezoid(density, m2_arr, axis=0)
        evidence = np.trapezoid(integral_over_m2, m1_arr)

    else:
        density = density_estimator.pdf(points, resample=True)

        # integrate over: ...
        integrate_dim = {
            "gw-2d": 1,  # q
            "gw-3d": 1,  # q
            "psr": 0,  # m
        }.get(event_type)
        evidence = np.trapezoid(y=density, x=points[:, integrate_dim])

    return evidence


class ModelSelector:
    """Approximate Bayesian model selection of the neutron star equation of state.

    This class is used for inference on single observed events with available parameter
    estimation results. Events may be gravitational wave detections of binary neutron star
    or neutron star-black hole mergers, or mass-radius measurements of pulsars.
    """

    def __init__(
        self,
        posterior_file: str,
        event_type: Literal["gw-2d", "gw-3d", "gw-4d", "psr"],
        density_est_method: Literal["kde", "bayes-flow", "ensemble-flow"] = "kde",
        flow_file: str | None = None,
        prior_file: str | None = None,
        q_min: float | None = None,
        q_max: float | None = None,
    ):
        """
        Parameters
        ----------
        posterior_file
            Path to a .json, .h5/.hdf5, or .dat file containing posterior samples of the necessary
            EOS-dependent parameters of the event.
            The needed parameter samples depends on ``event_type``:

            - "gw-2d" requires samples for mass ratio (q) and dominant tidal deformability (LambdaTilde)
            - "gw-3d" requires samples for mass ratio (q) and the two tidal deformabilities (Lambda1, Lambda2)
            - "gw-4d" requires samples for both masses and both tidal deformabilities (m1, m2, Lambda1, Lambda2)
            - "psr" requires samples for mass (in solar masses) and compactness

        event_type
            Type of the event/observation from which EOS inference is being conducted and the associated variant
            of the approximation scheme to use. Must be one of:

            - "gw-2d" : GW detection, 2 dimensional inference approximation scheme
            - "gw-3d" : GW detection, 3 dimensional inference approximation scheme
            - "gw-4d" : GW detection, 4 dimensional inference approximation scheme
            - "psr" : Pulsar mass-radius measurement

        density_est_method
            Choice of the type of density estimator to use during the inference. The density estimator is used to
            convert the data from the posterior file into a probability density function that can be integrated along
            EOS lines. Must be one of:

            - "kde" : (default) Gaussian kernel density estimator from Scipy, wrapped with gwxtreme.density_estimation.BoundedKDE
            - "bayes-flow" : Bayesian normalizing flow PyTorch/Zuko model, pre-trained on event data and wrapped with gwxtreme.density_estimation.BayesianNormalizingFlow. If chosen, ``flow_file`` must also be passed.
            - "ensemble-flow" : Set of normalizing flow PyTorch/Zuko models, trained on event data identically and only differing due to random weight initializations. This approach is designed to support a reproducible alternative to the Bayesian flow approach, with uncertainty estimation coming from the variance in density estimates from the ensemble of models. If chosen, ``flow-file`` must be passed.

        flow_file
            If ``density_est_method`` is "bayes-flow", provide a path to a .pt file containing the weights and configuration for a PyTorch/Zuko-based Bayesian Normalizing Flow model.
            If ``density_est_method`` is "ensemble-flow", provide a path to a directory containing the ensemble of PyTorch/Zuko-based Normalizing Flow models in the form of .pt files (1 per model).
            In either case, the .pt model files should correspond to models trained on the (transformed) posterior sample data for the given event, with the model parameterization matching that required based on ``event_type``.

        prior_file
            Path to a .json, .h5/.hdf5, or .dat file containing prior samples of the necessary
            EOS-dependent parameters of the event. This data is used to determine the bounds of integration for the EOS
            evidence integral. The required parameter(s) depends on the event type:

            - "gw-2d" : q_min and q_max are taken from the given prior samples
            - "gw-3d" : q_min and q_max are taken from the given prior samples
            - "gw-4d" : m1_min, m1_max, m2_min, and m2_max are taken from the given prior samples
            - "psr" : m_min and m_max are taken from the given prior samples

            If None, the bounds on the respective integration parameter(s) will be taken from the min and max of the posterior samples.

        q_min, q_max
            Prior bounds on q that will be used for EOS evidence integration for "gw-2d" and "gw-3d" events. Note that these
            arguments will be ignored if ``prior_file`` is supplied, in which case the bounds on q will be taken as the min and
            max values of q from the prior samples. If neither these parameters nor ``prior_file`` is given, the bounds on q will be
            taken from the min and max values of q in the posterior samples.
        """

        assert event_type in ["gw-2d", "gw-3d", "gw-4d", "psr"]
        assert density_est_method in ["kde", "bayes-flow", "ensemble-flow"]

        logger.info(
            "Creating ModelSelector with\nevent_type={}\ndensity_est_method={}\nflow_file={}\nposterior_file={}\nprior_file={}".format(
                event_type, density_est_method, flow_file, posterior_file, prior_file
            )
        )

        self.event_type = event_type
        self.density_est_method = density_est_method

        if self.event_type == "gw-2d":
            self.mean_mchirp = get_mean_mchirp_for_cbc_event(posterior_file, cbc_dim=2)

            logger.info("Mean chirp mass from CBC posterior = {}".format(self.mean_mchirp))

            if prior_file is not None:
                self.q_min, self.q_max = get_q_range(prior_file, cbc_dim=2)
            elif q_min is not None and q_max is not None:
                self.q_min, self.q_max = q_min, q_max
            else:
                self.q_min, self.q_max = get_q_range(posterior_file, cbc_dim=2)

            logger.info("q integration range = ({}, {})".format(self.q_min, self.q_max))

            posterior_samples = get_gw_event_pe_posterior_samples(posterior_file, cbc_dim=2)
            parameter_bounds = [(0.0, np.inf), (0.0, 1.0)]

        elif self.event_type == "gw-3d":
            self.mean_mchirp = get_mean_mchirp_for_cbc_event(posterior_file, cbc_dim=3)

            if prior_file is not None:
                self.q_min, self.q_max = get_q_range(prior_file, cbc_dim=3)
            elif q_min is not None and q_max is not None:
                self.q_min, self.q_max = q_min, q_max
            else:
                self.q_min, self.q_max = get_q_range(posterior_file, cbc_dim=3)

            logger.info("q integration range = ({}, {})".format(self.q_min, self.q_max))

            posterior_samples = get_gw_event_pe_posterior_samples(posterior_file, cbc_dim=3)
            parameter_bounds = [(0.0, np.inf), (0.0, 1.0), (0.0, np.inf)]

        elif self.event_type == "gw-4d":
            posterior_samples = get_gw_event_pe_posterior_samples(posterior_file, cbc_dim=4)
            parameter_bounds = [
                (0.0, np.inf),
                (0.0, np.inf),
                (0.0, np.inf),
                (0.0, np.inf),
            ]

            if prior_file is not None:
                prior_samples = get_gw_event_pe_posterior_samples(prior_file, cbc_dim=4)
                self.m1_min = np.min(prior_samples[:, 0])
                self.m1_max = np.max(prior_samples[:, 0])
                self.m2_min = np.min(prior_samples[:, 1])
                self.m2_max = np.max(prior_samples[:, 1])
            else:
                self.m1_min = np.min(posterior_samples[:, 0])
                self.m1_max = np.max(posterior_samples[:, 0])
                self.m2_min = np.min(posterior_samples[:, 1])
                self.m2_max = np.max(posterior_samples[:, 1])

        elif self.event_type == "psr":
            posterior_samples = get_nicer_pulsar_pe_posterior_samples(posterior_file)
            parameter_bounds = [(0.0, np.inf), (0.0, np.inf)]

            if prior_file is not None:
                prior_samples = get_nicer_pulsar_pe_posterior_samples(prior_file)
                self.m_min = np.min(prior_samples[:, 0])
                self.m_max = np.max(prior_samples[:, 0])
            else:
                self.m_min = np.min(posterior_samples[:, 0])
                self.m_max = np.max(posterior_samples[:, 0])

        logger.info("Posterior samples shape = {}".format(posterior_samples.shape))
        logger.info("Posterior density estimation bounds = {}".format(parameter_bounds))

        # Check if flow features are available in this version of the package. If not, code will terminate here with error.
        if density_est_method in ["bayes-flow", "ensemble-flow"] and (BayesianNormalizingFlow is None or EnsembleNormalizingFlow is None):
            raise NotImplementedError(
                "'bayes-flow' and 'ensemble-flow' density estimation methods are not implemented in this version of GWXtreme. \
                    Please use the available KDE implementation by specifying density_est_method='kde' (the default value)."
            )

        if density_est_method == "bayes-flow":
            assert flow_file is not None, (
                "To use 'bayes-flow' density estimator, must provide a path to the flow model file via the ``flow_file`` argument"
            )
            self.density_estimator = BayesianNormalizingFlow(
                posterior_samples=posterior_samples,
                bounds=parameter_bounds,
                flow_file=flow_file,
            )
        elif density_est_method == "kde":
            self.density_estimator = BoundedKDE(posterior_samples=posterior_samples, bounds=parameter_bounds)
        elif density_est_method == "ensemble-flow":
            assert flow_file is not None, (
                "To use 'ensemble-flow' density estimator, must provide a path to the directory of flow files via the ``flow_file`` argument"
            )
            self.density_estimator = EnsembleNormalizingFlow(posterior_samples=posterior_samples, bounds=parameter_bounds, flows_dir=flow_file)

    def compute_eos_evidence_ratio(
        self,
        target_eos_name: str | None = None,
        reference_eos_name: str | None = None,
        target_eos_mass_lambda_file: str | None = None,
        reference_eos_mass_lambda_file: str | None = None,
        target_eos_mass_radius_k_file: str | None = None,
        reference_eos_mass_radius_k_file: str | None = None,
        n_grid: int = 200,
        n_resamplings: int = 0,
        n_jobs: int = 1,
        save_file: str | None = None,
    ) -> np.ndarray:
        """Compute evidence ratios (Bayes factors) between EOS models.

        Bayes factors are computed as Target EOS Evidence / Reference EOS Evidence.

        Supply the two EOS models either as LALSuite model names, mass-lambda files, or mass-radius-tidal Love number files.
        The target and reference models may be supplied in different forms, but one ``target`` and one ``reference`` argument should be given.

        Parameters
        ----------
        target_eos_name, reference_eos_name
            Name of an EOS as implemented in LALSuite

        target_eos_mass_lambda_file, reference_eos_mass_lambda_file
            .txt file containing mass and tidal deformability (λ) values

            The data should be in the following format (without titles)::

                (mass column)       (lambda column)
                min_mass            ...
                ...                 ...
                ...                 ...
                max_mass            ...

            The values of masses should be in units of solar masses. The
            tidal deformability λ should be supplied in SI units.

            Note: Supplying an EOS in this form is not sufficient for inferences
            with 'psr'-type events, due to the inability to compute radii
            from mass and tidal deformability alone.

        target_eos_mass_radius_k_file, reference_eos_mass_radius_k_file
            .txt file containing mass, radii, and tidal Love number values

            The data should be in the following format (without titles)::

                (mass column)       (radius column)     (kappa column)
                min_mass            ...                 ...
                ...                 ...                 ...
                ...                 ...                 ...
                max_mass            ...                 ...

            The values of masses should be in units of solar masses. The radius should
            be supplied in meters.

        n_grid
            Number of points to use when computing each evidence integral, by default 200.
            Note: For 'gw-4d' events, the evidence integral is 2-dimensional over (m1, m2), so the square
            of the given number of points will be used. It is recommended to thus use a smaller
            number of points if using the 4D method, to prevent intractable computational time.

        n_resamplings
            Number of Bayes factor re-computations to perform by resampling the density estimator
            and re-integrating the probability density along the EOS line. These re-computed Bayes
            factor values are returned in an array along with the original Bayes factor. Default: 0
            NOTE: This parameter has no effect when the "ensemble-flow" density estimation method is used; the
            number of re-computed Bayes factors will be equal to the number of models provided in the ensemble.

        n_jobs
            Determines whether to use Ray for multi-core parallel processing of evidence re-computations,
            (the number of which are specified via the ``n_resamplings`` argument).
            By default (n_jobs = 1), the computation will be serial. Changing this to use parallel
            processing is only recommended for n_resamplings > 100.
            Options:

                - n_jobs = 1 (default): Serial execution; Ray not used.
                - n_jobs > 1 : Ray will be allocated the given number of CPU cores on the machine
                - n_jobs = -1 : Ray will be allocated 95% of the available CPU cores on the machine

            NOTE: This parameter has no effect when the "ensemble-flow" density estimation method is used.

        save_file
            Optional path to a json file that will be created and used to store Bayes factor
            results, with this schema:
            {
            "target_eos": ``target_eos_name`` (or file path if given),
            "reference_eos": ``reference_eos_name`` (or file path if given),
            "bayes_factors": [<original Bayes factor>, <n resampled Bayes factors>...]
            }

        Returns
        -------
            Array of Bayes factors (target EOS evidences / reference EOS evidences) with size ``n_resamplings`` + 1, structured like [<original Bayes factor>, <n resampled Bayes factors>...].
        """

        target_eos_evidences = self.compute_eos_evidence(
            eos_name=target_eos_name,
            eos_mass_lambda_file=target_eos_mass_lambda_file,
            eos_mass_radius_k_file=target_eos_mass_radius_k_file,
            n_grid=n_grid,
            n_resamplings=n_resamplings,
            n_jobs=n_jobs,
        )

        reference_eos_evidences = self.compute_eos_evidence(
            eos_name=reference_eos_name,
            eos_mass_lambda_file=reference_eos_mass_lambda_file,
            eos_mass_radius_k_file=reference_eos_mass_radius_k_file,
            n_grid=n_grid,
            n_resamplings=n_resamplings,
            n_jobs=n_jobs,
        )

        bayes_factors = target_eos_evidences / reference_eos_evidences

        if save_file is not None:
            logger.info("Saving Bayes factors to {}".format(save_file))
            results = {
                "target_eos": next(
                    (eos for eos in [target_eos_name, target_eos_mass_lambda_file, target_eos_mass_radius_k_file] if eos is not None), None
                ),
                "reference_eos": next(
                    (eos for eos in [reference_eos_name, reference_eos_mass_lambda_file, reference_eos_mass_radius_k_file] if eos is not None), None
                ),
                "bayes_factors": bayes_factors.tolist(),
            }
            with open(save_file, "w") as f:
                json.dump(results, f, indent=4)

        return bayes_factors

    def compute_eos_evidence(
        self,
        eos_name: str | None = None,
        eos_mass_lambda_file: str | None = None,
        eos_mass_radius_k_file: str | None = None,
        n_grid: int = 200,
        n_resamplings: int = 0,
        n_jobs: int = 1,
    ) -> np.ndarray:
        """Compute evidence for a single EOS.

        Supply EOS model either as a LALSuite model name, mass-lambda file, or mass-radius-tidal Love number file.

        Note that this evidence value is meaningless on its own because it is
        not normalized; this function should only ever be used when taking the
        evidence ratio (Bayes factor) between two different EOS models. ``ModelSelector.compute_eos_evidence_ratio``
        can be used directly for this purpose. This function is available for
        convenience and performance savings when computing Bayes factors between
        many EOS models and a single reference model (e.g. to prevent unnecessary
        repeated computations of the evidence for the reference model).

        Parameters
        ----------
        eos_name
            Name of an EOS as implemented in LALSuite

        eos_mass_lambda_file
            .txt file containing mass and tidal deformability (λ) values

            The data should be in the following format (without titles)::

                (mass column)       (lambda column)
                min_mass            ...
                ...                 ...
                ...                 ...
                max_mass            ...

            The values of masses should be in units of solar masses. The
            tidal deformability λ should be supplied in SI units.

            Note: Supplying an EOS in this form is not sufficient for inferences
            with 'psr'-type events, due to the inability to compute radii
            from mass and tidal deformability alone.

        eos_mass_radius_k_file
            .txt file containing mass, radii, and tidal Love number values

            The data should be in the following format (without titles)::

                (mass column)       (radius column)     (kappa column)
                min_mass            ...                 ...
                ...                 ...                 ...
                ...                 ...                 ...
                max_mass            ...                 ...

            The values of masses should be in units of solar masses. The radius should
            be supplied in meters.

        n_grid
            Number of points to use when computing each evidence integral, by default 200.
            Note: For 'gw-4d' events, the evidence integral is 2-dimensional over (m1, m2), so the square
            of the given number of points will be used. It is recommended to thus use a smaller
            number of points if using the 4D method, to prevent intractable computational time.

        n_resamplings
            Number of evidence re-computations to perform by resampling the density estimator
            and re-integrating the probability density along the EOS line. These re-computed evidence
            values are returned in an array along with the original value. Default: 0
            NOTE: This parameter has no effect when the "ensemble-flow" density estimation method is used; the
            number of re-computed evidences will be equal to the number of models provided in the ensemble.

        n_jobs
            Determines whether to use Ray for multi-core parallel processing of evidence re-computations,
            (the number of which are specified via the ``n_resamplings`` argument).
            By default (n_jobs = 1), the computation will be serial. Changing this to use parallel
            processing is only recommended for n_resamplings > 100.
            Options:

                - n_jobs = 1 (default): Serial execution; Ray not used.
                - n_jobs > 1 : Ray will be allocated the given number of CPU cores on the machine
                - n_jobs = -1 : Ray will be allocated 95% of the available CPU cores on the machine

            NOTE: This parameter has no effect when the "ensemble-flow" density estimation method is used.

        Returns
        -------
            Array of evidences with size ``n_resamplings`` + 1, structured like [<original evidence>, <n resampled evidences>...].
        """

        logger.info(
            "Computing evidence for EOS:\neos_name={}\neos_mass_lambda_file={}\neos_mass_radius_k_file={}".format(
                eos_name, eos_mass_lambda_file, eos_mass_radius_k_file
            )
        )

        interpolator = EOSInterpolator(eos_name=eos_name, mass_lambda_file=eos_mass_lambda_file, mass_radius_k_file=eos_mass_radius_k_file)
        eos_path = self.get_path_for_eos(interpolator, n_grid)

        logger.info("EOS path has shape {}".format(eos_path.shape))

        eos_evidences = self._path_evidence(eos_path, n_resamplings, n_jobs)

        return eos_evidences

    def compute_parameterized_eos_evidence(
        self,
        parameters: np.ndarray | tuple | list,
        parameterization: Literal["spectral", "polytrope"],
        parameter_bounds: list[tuple[float, float]],
        n_grid: int = 200,
    ) -> float:
        """Compute the evidence for a parameterized EOS model.

        Parameters
        ----------
        parameters
            Array of values for parameters characterizing the EOS
        parameterization
            Must be one of "spectral" (4-parameter spectral
            decomposition model) or "polytrope" (4-parameter
            piecewise-polytrope model)
        parameter_bounds
            Bounds for a uniform prior distribution of the EOS parameters, structured like
            [(g1_min, g1_max), (g2_min, g2_max), (g3_min, g3_max), (g4_min, g4_max)]

        Returns
        -------
            Evidence of the EOS constructed from the given parameters and parameterization
        """

        interpolator = EOSInterpolator(
            eos_parameters=np.array(parameters),
            parameterization=parameterization,
            parameter_bounds=parameter_bounds,
        )

        path = self.get_path_for_eos(interpolator, n_grid=n_grid)
        evidence = self._path_evidence(path)

        # Need to average the array of evidences corresponding to the
        # ensemble of models, if applicable.
        if isinstance(self.density_estimator, EnsembleNormalizingFlow):
            return np.mean(evidence).item()

        return evidence.item()

    def get_path_for_eos(
        self,
        eos_interpolator: EOSInterpolator,
        n_grid: int = 200,
    ) -> np.ndarray:
        """Use the given EOS interpolator to compute points of the
        line (or surface) that the EOS corresponds to in the event
        parameter space.

        The returned points are computed depending on the model
        selector ``event_type`` as follows:

        - "gw-2d" : A range of values for the NS component masses
                    are computed from a range of the mass ratio q
                    (from 0 to 1) and the event mean chirp mass.
                    These masses are then used to compute the
                    dominant tidal deformability (LambdaTilde).
                    A line of points (q, LambdaTilde) is returned
                    (with shape (n_grid, 2)).
        - "gw-3d" : The component masses are computed as described
                    above, and are then used to compute the tidal
                    deformabilities (Lambda1 and Lambda2). A line
                    of points (Lambda1, q, Lambda2) is returned (
                    with shape (n_grid, 3)).
        - "cbd-4d" : Component masses and tidal deformabilities are
                    computed as described above. A surface of points
                    (m1, m2, Lambda1, Lambda2) is returned (with
                    shape (n_grid, n_grid, 4)).
        - "psr" : The EOS-predicted NS radius is computed for a
                    uniform range of masses from 0.8 to 3.0
                    solar masses. A line of points (mass, radius) is
                    returned (with shape (n_grid, 2)).

        Parameters
        ----------
        eos_interpolator
            EOSInterpolator object encapsulating the EOS model which
            will be used to interpolate tidal deformabilities (or
            radii) from masses
        n_grid
           Number of points to return for the EOS line, by default 200.
           (For "gw-4d" event type, the number of returned points will be n_grid^2.)

        Returns
        -------
            Array of points with content and shape described above
        """

        if self.event_type == "gw-2d":
            q = np.linspace(self.q_min, self.q_max, n_grid)
            m1, m2 = convert_masses(q, self.mean_mchirp)
            m1, m2, q = eos_interpolator.apply_bns_mass_constraint(m1, m2, q)  # type: ignore

            lambdat = eos_interpolator.get_lambda_tilde(m1, m2)
            points = np.stack((lambdat, q), axis=-1)

        elif self.event_type == "gw-3d":
            q = np.linspace(self.q_min, self.q_max, n_grid)
            m1, m2 = convert_masses(q, self.mean_mchirp)
            m1, m2, q = eos_interpolator.apply_bns_mass_constraint(m1, m2, q)  # type: ignore

            lambda1 = eos_interpolator.get_lambda(m1)
            lambda2 = eos_interpolator.get_lambda(m2)
            points = np.stack((lambda1, q, lambda2), axis=-1)

        elif self.event_type == "gw-4d":
            m1 = np.linspace(self.m1_min, self.m1_max, n_grid)
            m2 = np.linspace(self.m2_min, self.m2_max, n_grid)
            m1_grid, m2_grid = np.meshgrid(m1, m2)

            lambda1 = eos_interpolator.get_lambda(m1_grid.reshape((n_grid**2))).reshape((n_grid, n_grid))
            lambda2 = eos_interpolator.get_lambda(m2_grid.reshape((n_grid**2))).reshape((n_grid, n_grid))

            points = np.stack((m1_grid, m2_grid, lambda1, lambda2), axis=-1)

        elif self.event_type == "psr":
            mass = np.linspace(self.m_min, self.m_max, n_grid)

            radius = eos_interpolator.get_radius(mass)
            points = np.stack((mass, radius), axis=-1)

        return points

    def _path_evidence(self, points: np.ndarray, n_resamplings: int = 0, n_jobs: int = 1) -> np.ndarray:
        """Compute the evidence(s) of the EOS line characterized by
        ``points`` over the event posterior density.

        Parameters
        ----------
        points
            Array of points as returned by ``ModelSelector.get_path_for_eos``

        n_resamplings
            Number of evidence recomputations to perform by resampling
            the posterior density estimator, by default 0.
            NOTE: This parameter has no effect when the "ensemble-flow" density estimation method is used; the
            number of re-computed evidences will be equal to the number of models provided in the ensemble.

        n_jobs
            Determines whether to use Ray for multi-core parallel processing of evidence re-computations,
            (the number of which are specified via the ``n_resamplings`` argument).
            By default (n_jobs = 1), the computation will be serial. Changing this to use parallel
            processing is only recommended for n_resamplings > 100.
            Options:

                - n_jobs = 1 (default): Serial execution; Ray not used.
                - n_jobs > 1 : Ray will be allocated the given number of CPU cores on the machine, up to a maximum of 95% of the total available cores
                - n_jobs = -1 : Ray will be allocated 95% of the available CPU cores on the machine
                - Any other given value will fall back to the default option

            NOTE: This parameter has no effect when the "ensemble-flow" density estimation method is used.

        Returns
        -------
            Array of evidences with size n_resamplings + 1
        """

        original_evidence = self._integrate_eos_path(points)
        if n_resamplings > 0 and not isinstance(self.density_estimator, EnsembleNormalizingFlow):
            logger.info("Re-computing evidence over {} re-samplings of the density estimator".format(n_resamplings))

            evidences = np.empty(n_resamplings + 1)
            evidences[0] = original_evidence

            # Serial execution (do not use Ray)
            #   n_jobs < -1 is incorrect input
            #   os.cpu_count() == None means cpu count is indeterminable
            cpu_count = os.cpu_count()
            if n_jobs == 1 or n_jobs == 0 or n_jobs < -1 or cpu_count is None:
                logger.info("Using serial execution")

                for i in range(1, n_resamplings + 1):
                    evidences[i] = self._integrate_eos_path(points, resample=True)

            # Parallel execution (use Ray)
            else:
                # Determine number of cores to use
                max_cores = np.floor(0.95 * cpu_count)

                if n_jobs == -1:
                    n_cores = max_cores
                else:
                    n_cores = min(n_jobs, max_cores)

                logger.info("Using parallel execution with {} cores".format(n_cores))

                # Spin up Ray with the specific cpu ceiling if not already running
                if not ray.is_initialized():
                    logger.info("Initializing Ray")
                    ray.init(num_cpus=n_cores, log_to_driver=False, logging_level=logging.WARNING)
                else:
                    logger.info("Available Ray cluster already exists; connecting to it")

                # Put these (somewhat) large objects in shared memory
                density_estimator_ref = ray.put(self.density_estimator)
                points_ref = ray.put(points)

                # Launch parallel tasks
                logger.info("Dispatching {} EOS evidence integration tasks to Ray cluster".format(n_resamplings))
                results = []
                for _ in range(n_resamplings):
                    results.append(
                        _distributed_eos_evidence_integration.options(enable_task_events=False).remote(
                            density_estimator_ref, points_ref, self.event_type
                        )
                    )

                evidences[1:] = ray.get(results)

                logger.info("Ray tasks returned")

            return evidences

        else:
            return original_evidence

    def _integrate_eos_path(self, points: np.ndarray, resample: bool = False) -> np.ndarray:
        """Compute the EOS evidence by integrating the EOS line
        characterized by ``points`` over the event posterior density.

        The integral is performed numerically using ``np.trapezoid``.

        Parameters
        ----------
        points
            Array of points as returned by ``ModelSelector.get_path_for_eos``
        resample
            Whether to resample the density estimator before
            integrating, by default False

        Returns
        -------
            Single-element array containing the evidence value
        """

        if not isinstance(self.density_estimator, EnsembleNormalizingFlow):
            if self.event_type == "gw-4d":
                m1_arr = points[0, :, 0]
                m2_arr = points[:, 0, 1]

                n_grid = points.shape[0]
                points = points.reshape(n_grid**2, 4)

                density = self.density_estimator.pdf(points, resample)
                density = density.reshape((n_grid, n_grid))

                # Integrate over m1 and m2
                integral_over_m2 = np.trapezoid(density, m2_arr, axis=0)
                evidence = np.trapezoid(integral_over_m2, m1_arr)

            else:
                density = self.density_estimator.pdf(points, resample)

                # integrate over: ...
                integrate_dim = {
                    "gw-2d": 1,  # q
                    "gw-3d": 1,  # q
                    "psr": 0,  # m
                }.get(self.event_type)
                evidence = np.trapezoid(y=density, x=points[:, integrate_dim])

            return evidence

        else:
            if self.event_type == "gw-4d":
                m1_arr = points[0, :, 0]
                m2_arr = points[:, 0, 1]

                n_grid = points.shape[0]
                points = points.reshape(n_grid**2, 4)

                ensemble_densities = self.density_estimator.pdf(points)

                ensemble_evidences = []
                for density in ensemble_densities:
                    density = density.reshape((n_grid, n_grid))

                    # Integrate over m1 and m2
                    integral_over_m2 = np.trapezoid(density, m2_arr, axis=0)
                    ensemble_evidences.append(np.trapezoid(integral_over_m2, m1_arr))

            else:
                ensemble_densities = self.density_estimator.pdf(points)

                # integrate over: ...
                integrate_dim = {
                    "gw-2d": 1,  # q
                    "gw-3d": 1,  # q
                    "psr": 0,  # m
                }.get(self.event_type)

                ensemble_evidences = []
                for density in ensemble_densities:
                    ensemble_evidences.append(np.trapezoid(y=density, x=points[:, integrate_dim]))

            return np.array(ensemble_evidences)

    def plot_eos_paths(
        self, eos_interpolators: list[EOSInterpolator], labels: list[str], colors: list[str], n_resamplings: int = 100, grid_size: int = 100
    ):
        paths = [self.get_path_for_eos(interp) for interp in eos_interpolators]

        assert len(eos_interpolators) == len(labels) == len(colors), "Must provide a line label and color for each EOS"

        plabels = {
            "gw-2d": [r"$\tilde{\Lambda}$", "$q$"],
            "gw-3d": [r"$\Lambda_1$", "$q$", r"$\Lambda_2$"],
            "gw-4d": [r"$m_1$", r"$m_2$", r"$\Lambda_1$", r"$\Lambda_2$"],
            "psr": [r"$m (M_{\odot})$", "$R (km)$"],
        }[self.event_type]

        pbounds = {
            "gw-2d": [(0.01, 3000), (0.001, 0.999)],
            "gw-3d": [(0.01, 3000), (0.001, 0.999), (0.01, 3000)],
            "gw-4d": [(1.0, 3.0), (1.0, 3.0), (0.0, 1500), (0.0, 1500)],
            "psr": [(0.8, 3.0), (8.0, 30.0)],
        }[self.event_type]

        fig, axes = corner_plot(
            self.density_estimator,
            parameter_labels=plabels,
            parameter_plot_bounds=pbounds,
            n_resamplings=n_resamplings,
            grid_size=grid_size,
        )

        ndim = len(plabels)
        for row in range(ndim):
            for col in range(ndim):
                if col < row:  # 2D marginal plots
                    for i, path in enumerate(paths):
                        x = path[:, col]
                        y = path[:, row]
                        axes[row, col].plot(x, y, label=labels[i], color=colors[i])
                    axes[row, col].legend()

        return fig, axes


class JointModelSelector:
    """Joint model selection of the neutron star equation-of-state.

    This class performs joint inference on several events by composing
    several ``ModelSelector`` instances, one for each event. Joint evidences
    are the product of the individual event evidences.
    """

    def __init__(
        self,
        posterior_files: list[str],
        event_types: list[str],
        density_est_method: Literal["kde", "bayes-flow", "ensemble-flow"] = "kde",
        flow_files: list[str] | None = None,
        prior_files: list[str] | None = None,
        q_mins: list[float] | None = None,
        q_maxes: list[float] | None = None,
    ):
        """
        Parameters
        ----------
        posterior_files
            Paths to .json, .h5/.hdf5, or .dat files containing posterior samples of the necessary
            EOS-dependent parameters of each event.
            The needed parameter samples for each file depends on the corresponding ``event_type``:

            - "gw-2d" requires samples for mass ratio (q) and dominant tidal deformability (LambdaTilde)
            - "gw-3d" requires samples for mass ratio (q) and the two tidal deformabilities (Lambda1, Lambda2)
            - "gw-4d" requires samples for both masses and both tidal deformabilities (m1, m2, Lambda1, Lambda2)
            - "psr" requires samples for mass (in solar masses) and compactness

        event_types
            Types of the events/observations from which EOS inference is being conducted and the associated variants
            of the approximation scheme to use. Must be a list with the same length as ``posterior_files``,
            where each element is one of:

            - "gw-2d" : GW detection, 2 dimensional inference approximation scheme
            - "gw-3d" : GW detection, 3 dimensional inference approximation scheme
            - "gw-4d" : GW detection, 4 dimensional inference approximation scheme
            - "psr" : Pulsar mass-radius measurement

        density_est_method
            Choice of the type of density estimator to use during the inference. The density estimator is used to
            convert the data from each posterior file into a probability density function that can be integrated along
            EOS lines. Must be one of:

            - "kde" : (default) Gaussian kernel density estimator from Scipy, wrapped with gwxtreme.density_estimation.BoundedKDE
            - "bayes-flow" : Bayesian normalizing flow PyTorch/Zuko model, pre-trained on event data and wrapped with gwxtreme.density_estimation.BayesianNormalizingFlow. If chosen, ``flow_files`` must also be passed.
            - "ensemble-flow" : Set of normalizing flow PyTorch/Zuko models, trained on event data identically and only differing due to random weight initializations. This approach is designed to support a reproducible alternative to the Bayesian flow approach, with uncertainty estimation coming from the variance in density estimates from the ensemble of models. If chosen, ``flow-files`` must be passed.

        flow_files
            If ``density_est_method`` is "bayes-flow", provide a list of paths (1 per event) to .pt files containing the weights and configurations for PyTorch/Zuko-based Bayesian Normalizing Flow models.
            If ``density_est_method`` is "ensemble-flow", provide a list of paths (1 per event) to directories containing the ensembles of PyTorch/Zuko-based Normalizing Flow models in the form of .pt files (1 per model).
            In either case, the .pt model files should correspond to models trained on the (transformed) posterior sample data for the given events, with the model parameterization matching that required based on ``event_types``.

        prior_files
            Paths to .json, .h5/.hdf5, or .dat files containing prior samples of the necessary
            EOS-dependent parameters for each event. This data is used to determine the bounds of integration for the EOS
            evidence integral for each event. The required parameter(s) depends on the event type:

            - "gw-2d" : q_min and q_max are taken from the given prior samples
            - "gw-3d" : q_min and q_max are taken from the given prior samples
            - "gw-4d" : m1_min, m1_max, m2_min, and m2_max are taken from the given prior samples
            - "psr" : m_min and m_max are taken from the given prior samples

            If None, the bounds on the respective integration parameter(s) will be taken from the min and max of the posterior samples.

        q_mins, q_maxes
            Prior bounds on q for each event that will be used for EOS evidence integration (only applicable to 'gw-2d' and 'gw-3d' events). Note that these
            arguments will be ignored if ``prior_files`` are supplied, in which case the bounds on q will be taken as the min and
            max values of q from the prior samples for each event. If neither these parameters nor ``prior_files`` is given, the bounds on q will be
            taken from the min and max values of q in each event's posterior samples.
        """

        logger.info(
            "Creating JointModelSelector with\nevent_types={}\ndensity_est_method={}\nflow_files={}\nposterior_files={}\nprior_files={}".format(
                event_types, density_est_method, flow_files, posterior_files, prior_files
            )
        )

        if density_est_method != "kde":
            assert flow_files is not None, (
                "If using 'bayes-flow' or 'ensemble-flow' density_est_method, must pass a set of model files or ensemble directory paths; see class __init__ docs"
            )
            assert len(posterior_files) == len(flow_files), (
                "Number of posterior_files should match the number of given flow_files when 'bayes-flow' or 'ensemble-flow' are chosen for density_est_method"
            )
        else:
            flow_files = [None] * len(posterior_files)  # type: ignore

        if prior_files is not None:
            assert len(posterior_files) == len(prior_files), "Number of posterior_files should match the number of prior_files"
        else:
            prior_files = [None] * len(posterior_files)  # type: ignore

        if q_mins is not None and q_maxes is not None:
            assert len(q_mins) == len(q_maxes) == len(posterior_files)
        else:
            q_mins = [None] * len(posterior_files)  # type: ignore
            q_maxes = [None] * len(posterior_files)  # type: ignore

        # For the parameterized EOS evidence calculation (used by the sampler class), define the number
        # of grid points to use: 100 points for all event types except the 4D method, which uses a 2-dimensional
        # surface integral, so 200 * 200 points will be used (40,000 KDE evaluations)
        self.parameterized_evidence_integral_n_grids = [200 if et == "gw-4d" else 100 for et in event_types]

        self.model_selectors = []
        for post_file, event_type, flow_file, prior_file, q_min, q_max in zip(posterior_files, event_types, flow_files, prior_files, q_mins, q_maxes):
            self.model_selectors.append(
                ModelSelector(
                    posterior_file=post_file,
                    event_type=event_type,
                    density_est_method=density_est_method,
                    flow_file=flow_file,
                    prior_file=prior_file,
                    q_min=q_min,
                    q_max=q_max,
                )
            )

    def compute_joint_eos_evidence_ratio(
        self,
        target_eos_name: str | None = None,
        reference_eos_name: str | None = None,
        target_eos_mass_lambda_file: str | None = None,
        reference_eos_mass_lambda_file: str | None = None,
        target_eos_mass_radius_k_file: str | None = None,
        reference_eos_mass_radius_k_file: str | None = None,
        n_grid: int = 200,
        n_resamplings: int = 0,
        n_jobs: int = 1,
        save_file: str | None = None,
    ) -> np.ndarray:
        """Compute joint evidence ratios (Bayes factors) between EOS models.

        Bayes factors are computed as Target EOS Evidence / Reference EOS Evidence.
        Joint bayes factors are obtained by multiplying the Bayes factors computed
        from individual events.

        Supply the two EOS models either as LALSuite model names, mass-lambda files, or mass-radius-tidal Love number files.
        The target and reference models may be supplied in different forms, but one ``target`` and one ``reference`` argument should be given.

        Parameters
        ----------
        target_eos_name, reference_eos_name
            Name of an EOS as implemented in LALSuite

        target_eos_mass_lambda_file, reference_eos_mass_lambda_file
            .txt file containing mass and tidal deformability (λ) values

            The data should be in the following format (without titles)::

                (mass column)       (lambda column)
                min_mass            ...
                ...                 ...
                ...                 ...
                max_mass            ...

            The values of masses should be in units of solar masses. The
            tidal deformability λ should be supplied in SI units.

            Note: Supplying an EOS in this form is not sufficient for inferences
            with 'psr'-type events, due to the inability to compute radii
            from mass and tidal deformability alone.

        target_eos_mass_radius_k_file, reference_eos_mass_radius_k_file
            .txt file containing mass, radii, and tidal Love number values

            The data should be in the following format (without titles)::

                (mass column)       (radius column)     (kappa column)
                min_mass            ...                 ...
                ...                 ...                 ...
                ...                 ...                 ...
                max_mass            ...                 ...

            The values of masses should be in units of solar masses. The radius should
            be supplied in meters.

        n_grid
            Number of points to use when computing each evidence integral, by default 200.
            Note: For 'gw-4d' events, the evidence integral is 2-dimensional over (m1, m2), so the square
            of the given number of points will be used. It is recommended to thus use a smaller
            number of points if using the 4D method, to prevent intractable computational time.

        n_resamplings
            Number of Bayes factor re-computations to perform by resampling the density estimator
            and re-integrating the probability density along the EOS line. These re-computed Bayes
            factor values are returned in an array along with the original Bayes factor. Default: 0.
            NOTE: This parameter has no effect when the "ensemble-flow" density estimation method is used; the
            number of re-computed Bayes factors will be equal to the number of models provided in the ensemble.

        n_jobs
            Determines whether to use Ray for multi-core parallel processing of evidence re-computations,
            (the number of which are specified via the ``n_resamplings`` argument).
            By default (n_jobs = 1), the computation will be serial. Changing this to use parallel
            processing is only recommended for n_resamplings > 100.
            Options:

                - n_jobs = 1 (default): Serial execution; Ray not used.
                - n_jobs > 1 : Ray will be allocated the given number of CPU cores on the machine, up to a maximum of 95% of the total available cores
                - n_jobs = -1 : Ray will be allocated 95% of the available CPU cores on the machine
                - Any other given value will fall back to the default option

            NOTE: This parameter has no effect when the "ensemble-flow" density estimation method is used.

        save_file
            Optional path to a json file that will be created and used to store Bayes factor
            results, with this schema:
            {
            "target_eos": ``target_eos_name`` (or file path if provided),
            "reference_eos": ``reference_eos_name`` (or file path if provided),
            "bayes_factors": [<original joint Bayes factor>, <resampled joint Bayes factors>...],
            "per_event_bayes_factors": [ [Bayes factors for event 1], [Bayes factors for event 2], ...]
            }

        Returns
        -------
            Array of joint Bayes factors (target EOS evidences / reference EOS evidences) with size ``n_resamplings`` + 1, structured like [<original Bayes factor>, <n resampled Bayes factors>...].
        """

        target_eos_joint_evidences, target_eos_per_event_evidences = self.compute_joint_eos_evidence(
            eos_name=target_eos_name,
            eos_mass_lambda_file=target_eos_mass_lambda_file,
            eos_mass_radius_k_file=target_eos_mass_radius_k_file,
            n_grid=n_grid,
            n_resamplings=n_resamplings,
            n_jobs=n_jobs,
        )

        reference_eos_joint_evidences, reference_eos_per_event_evidences = self.compute_joint_eos_evidence(
            eos_name=reference_eos_name,
            eos_mass_lambda_file=reference_eos_mass_lambda_file,
            eos_mass_radius_k_file=reference_eos_mass_radius_k_file,
            n_grid=n_grid,
            n_resamplings=n_resamplings,
            n_jobs=n_jobs,
        )

        per_event_bayes_factors = [
            (target_eos_event_ev / reference_eos_event_ev).tolist()
            for target_eos_event_ev, reference_eos_event_ev in zip(target_eos_per_event_evidences, reference_eos_per_event_evidences)
        ]
        joint_bayes_factors = np.prod(np.array(per_event_bayes_factors), axis=0)

        if save_file is not None:
            logger.info("Saving joint Bayes factors to {}".format(save_file))
            results = {
                "target_eos": next(
                    (eos for eos in [target_eos_name, target_eos_mass_lambda_file, target_eos_mass_radius_k_file] if eos is not None), None
                ),
                "reference_eos": next(
                    (eos for eos in [reference_eos_name, reference_eos_mass_lambda_file, reference_eos_mass_radius_k_file] if eos is not None), None
                ),
                "bayes_factors": joint_bayes_factors.tolist(),
                "per_event_bayes_factors": per_event_bayes_factors,
            }

            with open(save_file, "w") as f:
                json.dump(results, f, indent=4)

        return joint_bayes_factors

    def compute_joint_eos_evidence(
        self,
        eos_name: str | None = None,
        eos_mass_lambda_file: str | None = None,
        eos_mass_radius_k_file: str | None = None,
        n_grid: int = 200,
        n_resamplings: int = 0,
        n_jobs: int = 1,
    ) -> tuple[np.ndarray, list[np.ndarray]]:
        """Compute joint evidence for a single EOS.

        Supply EOS model either as a LALSuite model name, mass-lambda file, or mass-radius-tidal Love number file.

        Note that this evidence value is meaningless on its own because it is
        not normalized; this function should only ever be used when taking the
        evidence ratio (Bayes factor) between two different EOS models. ``JointModelSelector.compute_joint_eos_evidence_ratio``
        can be used directly for this purpose. This function is available for
        convenience and performance savings when computing Bayes factors between
        many EOS models and a single reference model (e.g. to prevent unnecessary
        repeated computations of the evidence for the reference model).

        Parameters
        ----------
        eos_name
            Name of an EOS as implemented in LALSuite

        eos_mass_lambda_file
            .txt file containing mass and tidal deformability (λ) values

            The data should be in the following format (without titles)::

                (mass column)       (lambda column)
                min_mass            ...
                ...                 ...
                ...                 ...
                max_mass            ...

            The values of masses should be in units of solar masses. The
            tidal deformability λ should be supplied in SI units.

            Note: Supplying an EOS in this form is not sufficient for inferences
            with 'psr'-type events, due to the inability to compute radii
            from mass and tidal deformability alone.

        eos_mass_radius_k_file
            .txt file containing mass, radii, and tidal Love number values

            The data should be in the following format (without titles)::

                (mass column)       (radius column)     (kappa column)
                min_mass            ...                 ...
                ...                 ...                 ...
                ...                 ...                 ...
                max_mass            ...                 ...

            The values of masses should be in units of solar masses. The radius should
            be supplied in meters.

        n_grid
            Number of points to use when computing each evidence integral, by default 200.
            Note: For 'gw-4d' events, the evidence integral is 2-dimensional over (m1, m2), so the square
            of the given number of points will be used. It is recommended to thus use a smaller
            number of points if using the 4D method, to prevent intractable computational time.

        n_resamplings
            Number of evidence re-computations to perform by resampling the density estimator
            and re-integrating the probability density along the EOS line. These re-computed evidence
            values are returned in an array along with the original value. Default: 0.
            NOTE: This parameter has no effect when the "ensemble-flow" density estimation method is used; the
            number of re-computed Bayes factors will be equal to the number of models provided in the ensemble.

        n_jobs
            Determines whether to use Ray for multi-core parallel processing of evidence re-computations,
            (the number of which are specified via the ``n_resamplings`` argument).
            By default (n_jobs = 1), the computation will be serial. Changing this to use parallel
            processing is only recommended for n_resamplings > 100.
            Options:

                - n_jobs = 1 (default): Serial execution; Ray not used.
                - n_jobs > 1 : Ray will be allocated the given number of CPU cores on the machine, up to a maximum of 95% of the total available cores
                - n_jobs = -1 : Ray will be allocated 95% of the available CPU cores on the machine
                - Any other given value will fall back to the default option

            NOTE: This parameter has no effect when the "ensemble-flow" density estimation method is used.

        Returns
        -------
            Tuple containing an array of values for the joint evidence with size ``n_resamplings`` + 1, structured like [<original evidence>, <n resampled evidences>...] and a list of evidence arrays (with the same size/structure) for each individual event
        """

        joint_evidences = np.ones(n_resamplings + 1, dtype=np.float32)
        per_event_evidences = []

        for model_selector in self.model_selectors:
            evidences = model_selector.compute_eos_evidence(
                eos_name=eos_name,
                eos_mass_lambda_file=eos_mass_lambda_file,
                eos_mass_radius_k_file=eos_mass_radius_k_file,
                n_grid=n_grid,
                n_resamplings=n_resamplings,
                n_jobs=n_jobs,
            )
            per_event_evidences.append(evidences)
            joint_evidences *= evidences

        return joint_evidences, per_event_evidences

    def compute_parameterized_eos_joint_evidence(
        self, parameters, parameterization: Literal["spectral", "polytrope"], parameter_bounds
    ) -> tuple[float, np.ndarray]:
        """Compute the joint evidence for a parameterized EOS model.

        The joint evidence is the product of the evidences from individual events.

        Parameters
        ----------
        parameters
            Array of values for parameters characterizing the EOS
        parameterization
            Must be one of "spectral" (4-parameter spectral
            decomposition model) or "polytrope" (4-parameter
            piecewise-polytrope model)
        parameter_bounds
            Bounds for a uniform prior distribution of the EOS parameters, structured like
            [(g1_min, g1_max), (g2_min, g2_max), (g3_min, g3_max), (g4_min, g4_max)]

        Returns
        -------
            Tuple where the first element is the joint evidence and the second element is an array containing the individual event evidences
        """

        per_event_evidences = []

        for i, model_selector in enumerate(self.model_selectors):
            per_event_evidences.append(
                model_selector.compute_parameterized_eos_evidence(
                    parameters, parameterization, parameter_bounds, self.parameterized_evidence_integral_n_grids[i]
                )
            )

        joint_evidence = np.prod(per_event_evidences)
        return float(joint_evidence), np.array(per_event_evidences)


class ParameterizedEoSSampler:
    """Performs MCMC parameter estimation for a parameterized EOS model using ``emcee``.

    The sampling likelihood is the joint evidence of the EOS from several provided events,
    as computed by ``JointModelSelector.compute_parameterized_eos_joint_evidence``.
    """

    def __init__(
        self,
        posterior_files: list[str],
        event_types: list[str],
        eos_prior_bounds: list[tuple],
        largest_observed_ns_mass: float = 1.97,
        density_est_method: Literal["kde", "bayes-flow", "ensemble-flow"] = "kde",
        flow_files: list[str] | None = None,
        parameterization: Literal["spectral", "polytrope"] = "spectral",
        prior_files: list[str] | None = None,
        q_mins: list[float] | None = None,
        q_maxes: list[float] | None = None,
    ):
        """
        Parameters
        ----------
        posterior_files
            Paths to .json, .h5/.hdf5, or .dat files containing posterior samples of the necessary
            EOS-dependent parameters of each event.
            The needed parameter samples for each file depends on the corresponding ``event_type``:

            - "gw-2d" requires samples for mass ratio (q) and dominant tidal deformability (LambdaTilde)
            - "gw-3d" requires samples for mass ratio (q) and the two tidal deformabilities (Lambda1, Lambda2)
            - "gw-4d" requires samples for both masses and both tidal deformabilities (m1, m2, Lambda1, Lambda2)
            - "psr" requires samples for mass (in solar masses) and compactness

        event_types
            Types of the events/observations from which EOS inference is being conducted and the associated variants
            of the approximation scheme to use. Must be a list with the same length as ``posterior_files``,
            where each element is one of:

            - "gw-2d" : GW detection, 2 dimensional inference approximation scheme
            - "gw-3d" : GW detection, 3 dimensional inference approximation scheme
            - "gw-4d" : GW detection, 4 dimensional inference approximation scheme
            - "psr" : Pulsar mass-radius measurement

        eos_prior_bounds
            Bounds for a uniform prior distribution of the EOS parameters, structured like
            [(g1_min, g1_max), (g2_min, g2_max), (g3_min, g3_max), (g4_min, g4_max)]

        largest_observed_ns_mass
            Mass of the heaviest observed NS, in solar masses. During sampling, EOSs will be required to
            support a maximum mass at least this large to have non-zero likelihood.

        density_est_method
            Choice of the type of density estimator to use during the inference. The density estimator is used to
            convert the data from each posterior file into a probability density function that can be integrated along
            EOS lines. Must be one of:

            - "kde" : (default) Gaussian kernel density estimator from Scipy, wrapped with gwxtreme.density_estimation.BoundedKDE
            - "bayes-flow" : Bayesian normalizing flow PyTorch/Zuko model, pre-trained on event data and wrapped with gwxtreme.density_estimation.BayesianNormalizingFlow. If chosen, ``flow_files`` must also be passed.
            - "ensemble-flow" : Set of normalizing flow PyTorch/Zuko models, trained on event data identically and only differing due to random weight initializations. This approach is designed to support a reproducible alternative to the Bayesian flow approach, with uncertainty estimation coming from the variance in density estimates from the ensemble of models. If chosen, ``flow-files`` must be passed.

        flow_files
            If ``density_est_method`` is "bayes-flow", provide a list of paths (1 per event) to .pt files containing the weights and configurations for PyTorch/Zuko-based Bayesian Normalizing Flow models.
            If ``density_est_method`` is "ensemble-flow", provide a list of paths (1 per event) to directories containing the ensembles of PyTorch/Zuko-based Normalizing Flow models in the form of .pt files (1 per model).
            In either case, the .pt model files should correspond to models trained on the (transformed) posterior sample data for the given events, with the model parameterization matching that required based on ``event_types``.

        parameterization
            Must be one of "spectral" (4-parameter spectral decomposition model) or "polytrope" (4-parameter piecewise-polytrope model)

        prior_files
            Paths to .json, .h5/.hdf5, or .dat files containing prior samples of the necessary
            EOS-dependent parameters for each event. This data is used to determine the bounds of integration for the EOS
            evidence integral for each event. The required parameter(s) depends on the event type:

            - "gw-2d" : q_min and q_max are taken from the given prior samples
            - "gw-3d" : q_min and q_max are taken from the given prior samples
            - "gw-4d" : m1_min, m1_max, m2_min, and m2_max are taken from the given prior samples
            - "psr" : m_min and m_max are taken from the given prior samples

            If None, the bounds on the respective integration parameter(s) will be taken from the min and max of the posterior samples.

        q_mins, q_maxes
            Prior bounds on q for each event that will be used for EOS evidence integration (only applicable to 'gw-2d' and 'gw-3d' events). Note that these
            arguments will be ignored if ``prior_files`` are supplied, in which case the bounds on q will be taken as the min and
            max values of q from the prior samples for each event. If neither these parameters nor ``prior_files`` is given, the bounds on q will be
            taken from the min and max values of q in each event's posterior samples.
        """
        logger.info(
            "Creating ParameterizedEoSSampler with\nevent_types={}\ndensity_est_method={}\nflow_files={}\nposterior_files={} \
                \nprior_files={}\nparameterization={}\neos_prior_bounds={}".format(
                event_types, density_est_method, flow_files, posterior_files, prior_files, parameterization, eos_prior_bounds
            )
        )

        self.eos_prior_bounds = eos_prior_bounds
        self.parameterization = parameterization
        self.largest_observed_ns_mass = largest_observed_ns_mass

        self.joint_selector = JointModelSelector(
            posterior_files=posterior_files,
            event_types=event_types,
            density_est_method=density_est_method,
            flow_files=flow_files,
            prior_files=prior_files,
            q_mins=q_mins,
            q_maxes=q_maxes,
        )

    def log_post(self, parameters: np.ndarray) -> float:
        """Compute the log posterior density (MCMC sampling likelihood) for
        the given ``parameters``.

        If the given parameter vector describes a valid EOS (i.e. adheres to
        uniform prior bounds given when instantiating this class, adheres to
        causality, thermodynamic stability, and observational consistency with
        the most massive observed neutron star) then the log likelihood given by
        this function is simply the log of the multi-event joint evidence for the EOS
        as computed by ``JointModelSelector.compute_parameterized_eos_joint_evidence``.

        If the EOS is not valid, returns -inf.

        Parameters
        ----------
        parameters
            Array of values for parameters characterizing the EOS

        Returns
        -------
            Joint evidence of the EOS
        """

        # Given parameter vector should lie within specified uniform prior
        if not all([parameters[i] >= self.eos_prior_bounds[i][0] and parameters[i] <= self.eos_prior_bounds[i][1] for i in range(4)]):
            return -np.inf

        # Checking for physical and observational consistency
        if not is_valid_eos(parameters, self.parameterization, largest_ns_mass=self.largest_observed_ns_mass):
            return -np.inf

        joint_evidence, _ = self.joint_selector.compute_parameterized_eos_joint_evidence(parameters, self.parameterization, self.eos_prior_bounds)
        log_evidence = np.log(joint_evidence)
        return log_evidence

    def _initialize_walkers(self, nwalkers: int) -> list[np.ndarray]:
        """Initializes the walkers for MCMC by constructing
        an initial state in the EOS parameter space for each walker.

        This initial state must satisfy the prior and physical
        constraints on the EOS.

        Parameters
        ----------
        nwalkers
            Number of MCMC walkers to use; this value is passed
            directly to the ``emcee`` sampler
        """

        logger.info("Initializing {} MCMC walkers".format(nwalkers))

        n_valid_walkers = 0
        state0 = []

        param_lower_bounds = [bound[0] for bound in self.eos_prior_bounds]
        param_upper_bounds = [bound[1] for bound in self.eos_prior_bounds]

        while n_valid_walkers < nwalkers:
            params = np.random.uniform(param_lower_bounds, param_upper_bounds)

            if is_valid_eos(params, self.parameterization, largest_ns_mass=self.largest_observed_ns_mass):
                logger.info("Attempted EOS with parameter vector {} is valid; appending to initial walker state".format(params))

                state0.append(params)
                n_valid_walkers += 1

        return state0

    def run_sampler(self, nsteps: int, nwalkers: int, save_file: str, reset: bool = True) -> None:
        """Runs MCMC to sample EOS parameters from their joint posterior across
        all provided events.

        The ``emcee.EnsembleSampler`` is used; the below parameters are identical to
        how they are defined by the ``emcee`` API.

        Note: During the run, the progress of the sampler will be saved in the given file.
        This can be used to resume previous runs (see the ``reset`` argument). You should
        NOT try to open this file to inspect the chain while the sampler is running; this
        will crash the sampler.

        Parameters
        ----------
        nsteps
            Number of steps to run. Note: this value is a maximum number of steps; if the
            sampler converges before reaching ``nsteps``, it will terminate early. Convergence
            is judged by looking at the estimated autocorrelation time of the chain as computed
            by ``emcee``. If the chain is longer than 100 * the autocorr time and if the estimate
            of the autocorr time changes by less than 1% between checks, then the chain is
            considered to have converged.
        nwalkers
            Number of MCMC walkers to use.
        save_file
            Path to an hdf5 file that will be used to store the sample chain and the associated log
            probability densities. This file can be read after the run with the ``load_samples`` function.
        reset
            Whether to overwrite the samples currently stored in the given ``save_file``, such
            as from a previous run with the same data/method. If False, the sampler will continue
            from where the previous run left off, using the last sample in the stored chain as the
            initial state of the walkers. If True (default), anything in ``save_file`` will be
            overwritten with a fresh run. See https://emcee.readthedocs.io/en/stable/tutorials/monitor/
        """

        logger.info("Running MCMC for {} EOS with {} walkers for {} steps".format(self.parameterization, nwalkers, nsteps))

        # Won't be able to create the file if parent dir doesn't exist beforehand
        if not pathlib.Path(save_file).parent.exists():
            raise FileNotFoundError("Parent directory of given save_file doesn't exist: {}".format(save_file))

        logger.info("Saving samples and logp chains to {}".format(save_file))
        backend = emcee.backends.HDFBackend(save_file)

        backend_is_empty = False
        try:
            backend.get_chain()
        except AttributeError:
            backend_is_empty = True

        if reset or backend_is_empty:
            backend.reset(nwalkers=nwalkers, ndim=4)
            initial_state = self._initialize_walkers(nwalkers)
        else:
            initial_state = None

        sampler = emcee.EnsembleSampler(nwalkers=nwalkers, ndim=4, log_prob_fn=self.log_post, backend=backend)

        # Track average autocorrelation time to determine convergence
        old_tau = np.inf

        for sample in sampler.sample(initial_state if initial_state else sampler.get_last_sample(), iterations=nsteps, progress=True):
            # Check for convergence every 100 steps
            if sampler.iteration % 100:
                continue

            # Compute autocorrelation time
            tau = sampler.get_autocorr_time(tol=0)

            logger.info(f"Estimated autocorrelation time at iteration {sampler.iteration}: {tau}")

            # Check for convergence (chain is longer than 100 *
            # the autocorr time and if the estimated autocorr
            # time has changed by less than 1%)
            converged = np.all(tau * 100 < sampler.iteration) and np.all(np.abs(old_tau - tau) / tau < 0.01)
            if converged:
                logger.info(f"Sampler converged at iteration {sampler.iteration} with estimated autocorrelation time {tau}")
                break

            old_tau = tau


def load_samples(
    samples_file: str,
    burn_in_frac: float = 0.5,
    thin: int | None = None,
) -> np.ndarray:
    """Read the given HDF5 samples chain file and return the samples
    as a flat chain (all walkers combined).

    Optionall discard the first ``burn_in_frac`` percentage of ``samples``
    and return every ``thin``'th sample from the remaining array.

    See https://emcee.readthedocs.io/en/stable/tutorials/autocorr/ for some documentation
    on choosing thinning and burn-in options.

    Parameters
    ----------
    samples_file
        Path to hdf5 file containing the samples chain
    burn_in_frac
        Percentage of samples to discard from the beginning of the chain.
        Default: 0.5
    thin
        Take only every ``thin`` samples from the chain. By default,
        1/2 of the estimated maximum integrated autocorrelation time will be used,
        if the chain is long enough for ``emcee`` to compute this value. If
        not, ``thin`` will be taken as 1/50 th of the chain length. For no
        thinning, provide ``thin = 1``.

    Returns
    -------
        Samples chain array with shape (nsteps * n_walkers, 4)
    """

    reader = emcee.backends.HDFBackend(samples_file)
    samples = reader.get_chain(flat=False)

    burn_in = int(samples.shape[0] * burn_in_frac)

    if not thin:
        try:
            thin = int(max(np.array(emcee.autocorr.integrated_time(samples))) / 2.0)
        except emcee.autocorr.AutocorrError as e:
            logger.exception(e)
            thin = int(samples.shape[0] / 50)

    thin = max(thin, 1)

    return reader.get_chain(flat=True, discard=burn_in, thin=thin)
