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

"""Neutron star equation of state inference algorithms."""

import atexit
import json
import logging
import os
import pathlib
import shutil
from collections.abc import Sequence
from typing import Literal

import emcee
import h5py
import lal
import numpy as np
import ray

from gwxtreme.density_estimation import BoundedKDE, EnsembleNormalizingFlow
from gwxtreme.eos_interpolator import EOSInterpolator, convert_masses
from gwxtreme.eos_prior import is_valid_eos

logger = logging.getLogger(__name__)


_owns_ray_session = False
_ray_session_dir: str | None = None


def _shutdown_owned_ray():
    if _owns_ray_session and _ray_session_dir is not None:
        ray.shutdown()
        shutil.rmtree(_ray_session_dir, ignore_errors=True)


@ray.remote
def _distributed_eos_evidence_integration(density_estimator: BoundedKDE, points: np.ndarray, event_type: str):
    if event_type == "gw-4d":
        q_arr = points[0, :, 0]
        mchirp_arr = points[:, 0, 1]

        n_grid = points.shape[0]
        points = points.reshape(n_grid**2, 4)

        density = density_estimator.pdf(points, resample=True)
        density = density.reshape((n_grid, n_grid))

        integral_over_mchirp = np.trapezoid(density, mchirp_arr, axis=0)
        evidence = np.trapezoid(integral_over_mchirp, q_arr)

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

    This class is used for inference on neutron star merger PE datasets or NICER pulsar
    mass-radius samples.
    """

    def __init__(
        self,
        posterior_files: Sequence[str],
        event_types: Sequence[str],
        density_est_method: Literal["kde", "flow"] = "kde",
        flow_directories: Sequence[str | pathlib.Path] | None = None,
        integration_bounds: Sequence[tuple[float, float]] | Sequence[Sequence[tuple[float, float]]] | None = None,
    ):
        """
        Parameters
        ----------
        posterior_files
            Paths to .json, .h5/.hdf5, or .dat files containing posterior samples of the necessary
            EOS-dependent parameters of each event.
            The needed parameter samples for each file depends on the corresponding ``event_type``:

            - "gw-2d" requires samples for mass ratio (q), chirp mass, and dominant tidal deformability (LambdaTilde)
            - "gw-3d" and "gw-4d" require samples for mass ratio (q), chirp mass, and the two tidal deformabilities (Lambda1, Lambda2)
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
            - "flow" : Set of normalizing flow PyTorch/Zuko models, trained on event data identically and only differing due to random weight initializations. This approach is designed to support a reproducible alternative to the Bayesian flow approach, with uncertainty estimation coming from the variance in density estimates from the ensemble of models. If chosen, ``flow-files`` must be passed.

        flow_directories
            If ``density_est_method`` is "flow", provide a list of paths (1 per event) to directories containing the ensembles of Normalizing Flow models in the form of .onnx files (1 per model).

        integration_bounds
            Bounds for the EOS evidence integral.

            These should correspond to:

            - (q_min, q_max) for "gw-2d" and "gw-3d" events
            - [(q_min, q_max), (mchirp_min, mchirp_max)] for "gw-4d" events
            - (mass_min, mass_max) for "psr" events (with mass in solar masses)

            If None, the bounds on the respective integration parameter(s) will be taken from the min and max of the posterior samples.
        """

        data_space_bounds = {
            "gw-2d": [(0.0, np.inf), (0.0, 1.0)],
            "gw-3d": [(0.0, np.inf), (0.0, 1.0), (0.0, np.inf)],
            "gw-4d": [
                (0.0, 1.0),
                (0.0, np.inf),
                (0.0, np.inf),
                (0.0, np.inf),
            ],
            "psr": [(0.0, np.inf), (0.0, np.inf)],
        }

        self.density_estimators = []
        self.integration_bounds = []
        self.mean_mchirps = []
        self.event_types = event_types

        for i in range(len(posterior_files)):
            posterior_samples = get_posterior_samples(posterior_files[i], event_types[i])

            if event_types[i] == "gw-2d" or event_types[i] == "gw-3d":
                if event_types[i] == "gw-2d":
                    samples_array = np.stack((posterior_samples["lambdat"], posterior_samples["q"]), axis=-1, dtype=np.float32)
                else:
                    samples_array = np.stack(
                        (posterior_samples["lambda1"], posterior_samples["q"], posterior_samples["lambda2"]), axis=-1, dtype=np.float32
                    )

                self.mean_mchirps.append(np.mean(posterior_samples["mchirp"]))

                if integration_bounds is not None:
                    self.integration_bounds.append(integration_bounds)
                else:
                    self.integration_bounds.append((np.min(posterior_samples["q"]), np.max(posterior_samples["q"])))

            elif event_types[i] == "gw-4d":
                samples_array = np.stack(
                    (posterior_samples["q"], posterior_samples["mchirp"], posterior_samples["lambda1"], posterior_samples["lambda2"]),
                    axis=-1,
                    dtype=np.float32,
                )

                self.mean_mchirps.append(None)

                if integration_bounds is not None:
                    if not hasattr(integration_bounds[0], "__len__"):
                        raise ValueError(
                            f"For 'gw-4d' events, must provide integration bounds for both q and mchirp, such as [(0.3, 1.0), (1.1, 1.4)] - got {integration_bounds}"
                        )
                    self.integration_bounds.append(integration_bounds)
                else:
                    self.integration_bounds.append(
                        (
                            (np.min(posterior_samples["q"]), np.max(posterior_samples["q"])),
                            (np.min(posterior_samples["mchirp"]), np.max(posterior_samples["mchirp"])),
                        )
                    )

            elif event_types[i] == "psr":
                samples_array = np.stack((posterior_samples["mass"], posterior_samples["radius"]), axis=-1, dtype=np.float32)

                self.mean_mchirps.append(None)

                if integration_bounds is not None:
                    self.integration_bounds.append(integration_bounds)
                else:
                    self.integration_bounds.append((np.min(posterior_samples["mass"]), np.max(posterior_samples["mass"])))

            if density_est_method == "kde":
                self.density_estimators.append(BoundedKDE(posterior_samples=samples_array, bounds=data_space_bounds[event_types[i]]))
            elif density_est_method == "flow":
                self.density_estimators.append(EnsembleNormalizingFlow(flows_dir=str(flow_directories[i]), bounds=data_space_bounds[event_types[i]]))

        # For the parameterized EOS evidence calculation (used by the sampler class), define the number
        # of grid points to use: 100 points for all event types except the 4D method, which uses a 2-dimensional
        # surface integral, so 200 * 200 points will be used (40,000 KDE evaluations)
        self.parameterized_evidence_integral_n_grids = [200 if et == "gw-4d" else 100 for et in event_types]

    def evidence_ratio(
        self,
        target_eos_name: str | None = None,
        reference_eos_name: str | None = None,
        target_eos_mass_lambda_file: str | None = None,
        reference_eos_mass_lambda_file: str | None = None,
        target_eos_mass_radius_k_file: str | None = None,
        reference_eos_mass_radius_k_file: str | None = None,
        n_grid: int = 100,
        n_resamplings: int = 0,
        n_jobs: int = 1,
        ray_address: str | None = None,
        save_file: str | None = None,
    ) -> np.ndarray:
        """Compute evidence ratios (Bayes factors) between EOS models.

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
            .txt file containing mass and dimensionless tidal deformability values

            The data should be in the following format (without titles)::

                (mass column)       (Lambda column)
                min_mass            ...
                ...                 ...
                ...                 ...
                max_mass            ...

            The values of masses should be in units of solar masses. The
            tidal deformability values should be dimensionless.

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
            Note: For 'gw-4d' events, the evidence integral is 2-dimensional over (q, mchirp), so the square
            of the given number of points will be used. It is recommended to thus use a smaller
            number of points if using the 4D method, to prevent intractable computational time.

        n_resamplings
            Number of Bayes factor re-computations to perform by resampling the density estimator
            and re-integrating the probability density along the EOS line. These re-computed Bayes
            factor values are returned in an array along with the original Bayes factor. Default: 0.
            NOTE: This parameter has no effect when the "flow" density estimation method is used; the
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

            NOTE: This parameter has no effect when the "flow" density estimation method is used.

        ray_address
            Address of Ray cluster to use for parallel processing when ``n_jobs`` is set. By default,
            any existing (already initialized) Ray cluster will be used, or one will be initialized
            if none exists.

        save_file
            Optional path to a json file that will be created and used to store Bayes factor
            results, with this schema:
            {
            "target_eos": ``target_eos_name`` (or file path if provided),
            "reference_eos": ``reference_eos_name`` (or file path if provided),
            "bayes_factors": [<original Bayes factor>, <resampled Bayes factors>...],
            "per_event_bayes_factors": [ [Bayes factors for event 1], [Bayes factors for event 2], ...]
            }

        Returns
        -------
            Array of Bayes factors (target EOS evidences / reference EOS evidences) with size ``n_resamplings`` + 1, structured like [<original Bayes factor>, <n resampled Bayes factors>...].
        """

        _, target_eos_per_event_evidences = self._joint_evidence(
            eos_name=target_eos_name,
            eos_mass_lambda_file=target_eos_mass_lambda_file,
            eos_mass_radius_k_file=target_eos_mass_radius_k_file,
            n_grid=n_grid,
            n_resamplings=n_resamplings,
            n_jobs=n_jobs,
            ray_address=ray_address,
        )

        _, reference_eos_per_event_evidences = self._joint_evidence(
            eos_name=reference_eos_name,
            eos_mass_lambda_file=reference_eos_mass_lambda_file,
            eos_mass_radius_k_file=reference_eos_mass_radius_k_file,
            n_grid=n_grid,
            n_resamplings=n_resamplings,
            n_jobs=n_jobs,
            ray_address=ray_address,
        )

        per_event_bayes_factors = [
            (target_eos_event_ev / reference_eos_event_ev).tolist()
            for target_eos_event_ev, reference_eos_event_ev in zip(target_eos_per_event_evidences, reference_eos_per_event_evidences)
        ]
        joint_bayes_factors = np.prod(np.array(per_event_bayes_factors), axis=0)

        if save_file is not None:
            logger.info(f"Saving Bayes factors to {save_file}")
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

    def _joint_evidence(
        self,
        eos_name: str | None = None,
        eos_mass_lambda_file: str | None = None,
        eos_mass_radius_k_file: str | None = None,
        n_grid: int = 200,
        n_resamplings: int = 0,
        n_jobs: int = 1,
        ray_address: str | None = None,
    ) -> tuple[np.ndarray, list[np.ndarray]]:
        """Compute joint evidence for a single EOS.

        Supply EOS model either as a LALSuite model name, mass-lambda file, or mass-radius-tidal Love number file.

        Note that this evidence value is meaningless on its own because it is
        not normalized; this function should only ever be used when taking the
        evidence ratio (Bayes factor) between two different EOS models. ``ModelSelector.evidence_ratio``
        can be used directly for this purpose. This function is available for
        convenience and performance savings when computing Bayes factors between
        many EOS models and a single reference model (e.g. to prevent unnecessary
        repeated computations of the evidence for the reference model).

        Parameters
        ----------
        eos_name
            Name of an EOS as implemented in LALSuite

        eos_mass_lambda_file
            .txt file containing mass and dimensionless tidal deformability values

            The data should be in the following format (without titles)::

                (mass column)       (Lambda column)
                min_mass            ...
                ...                 ...
                ...                 ...
                max_mass            ...

            The values of masses should be in units of solar masses. The
            tidal deformability values should be dimensionless.

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
            Note: For 'gw-4d' events, the evidence integral is 2-dimensional over (q, mchirp), so the square
            of the given number of points will be used. It is recommended to thus use a smaller
            number of points if using the 4D method, to prevent intractable computational time.

        n_resamplings
            Number of evidence re-computations to perform by resampling the density estimator
            and re-integrating the probability density along the EOS line. These re-computed evidence
            values are returned in an array along with the original value. Default: 0.
            NOTE: This parameter has no effect when the "flow" density estimation method is used; the
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

            NOTE: This parameter has no effect when the "flow" density estimation method is used.

        ray_address
            Address of Ray cluster to use for parallel processing when ``n_jobs`` is set. By default,
            any existing (already initialized) Ray cluster will be used, or one will be initialized
            if none exists.

        Returns
        -------
            Tuple containing an array of values for the joint evidence with size ``n_resamplings`` + 1, structured like [<original evidence>, <n resampled evidences>...] and a list of evidence arrays (with the same size/structure) for each individual event
        """

        interpolator = EOSInterpolator(
            eos_name=eos_name,
            mass_lambda_file=eos_mass_lambda_file,
            mass_radius_k_file=eos_mass_radius_k_file,
        )
        per_event_evidences = []
        for i in range(len(self.density_estimators)):
            evidences = self._single_event_evidence(
                interpolator,
                self.density_estimators[i],
                self.event_types[i],
                self.integration_bounds[i],
                self.mean_mchirps[i],
                n_grid=n_grid,
                n_resamplings=n_resamplings,
                n_jobs=n_jobs,
                ray_address=ray_address,
            )
            per_event_evidences.append(evidences)

        joint_evidences = np.prod(per_event_evidences, axis=0)
        return joint_evidences, per_event_evidences

    def _single_event_evidence(
        self,
        interpolator: EOSInterpolator,
        density_estimator: BoundedKDE | EnsembleNormalizingFlow,
        event_type: str,
        integration_bounds,
        mean_mchirp: float | None,
        n_grid: int = 100,
        n_resamplings: int = 0,
        n_jobs: int = 1,
        ray_address: str | None = None,
    ):
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
        n_grid
            Number of points to use when computing each evidence integral, by default 200.
            Note: For 'gw-4d' events, the evidence integral is 2-dimensional over (q, mchirp), so the square
            of the given number of points will be used. It is recommended to thus use a smaller
            number of points if using the 4D method, to prevent intractable computational time.

        n_resamplings
            Number of evidence re-computations to perform by resampling the density estimator
            and re-integrating the probability density along the EOS line. These re-computed evidence
            values are returned in an array along with the original value. Default: 0
            NOTE: This parameter has no effect when the "flow" density estimation method is used; the
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

            NOTE: This parameter has no effect when the "flow" density estimation method is used.

        ray_address
            Address of Ray cluster to use for parallel processing when ``n_jobs`` is set. By default,
            any existing (already initialized) Ray cluster will be used, or one will be initialized
            if none exists.

        Returns
        -------
            Array of evidences with size ``n_resamplings`` + 1, structured like [<original evidence>, <n resampled evidences>...].
        """

        points = self._create_eos_contour(interpolator, event_type, integration_bounds, mean_mchirp, n_grid)

        original_evidence = self._integrate_eos_contour(points, density_estimator, event_type)
        if n_resamplings > 0 and isinstance(density_estimator, BoundedKDE):
            logger.info(f"Re-computing evidence over {n_resamplings} re-samplings of the density estimator")

            evidences = np.empty(n_resamplings + 1)
            evidences[0] = original_evidence

            # Serial execution (do not use Ray)
            #   n_jobs < -1 is incorrect input
            #   os.cpu_count() == None means cpu count is indeterminable
            cpu_count = os.cpu_count()
            if n_jobs == 1 or n_jobs == 0 or n_jobs < -1 or cpu_count is None:
                logger.info("Using serial execution")

                for i in range(1, n_resamplings + 1):
                    evidences[i] = self._integrate_eos_contour(points, density_estimator, event_type, resample=True)

            # Parallel execution (use Ray)
            else:
                # Determine number of cores to use
                max_cores = np.floor(0.95 * cpu_count)

                if n_jobs == -1:
                    n_cores = max_cores
                else:
                    n_cores = min(n_jobs, max_cores)

                logger.info(f"Using parallel execution with {n_cores} cores")

                # Spin up Ray with the specific cpu ceiling if not already running
                if not ray.is_initialized():
                    logger.info("Initializing Ray")

                    if ray_address is not None:
                        # Connecting to a user-specified existing cluster: resource
                        # arguments can't be combined with `address`, and this process
                        # must never manage the lifecycle of a cluster it didn't create.
                        ray.init(address=ray_address, log_to_driver=False, logging_level=logging.WARNING)
                    else:
                        ctx = ray.init(num_cpus=n_cores, log_to_driver=False, logging_level=logging.WARNING)

                        # Register shutdown procedure to kill Ray processes and remove the
                        # session dir it creates in /tmp/ray. Only done when this call is
                        # the one that started the cluster.
                        global _owns_ray_session, _ray_session_dir
                        _owns_ray_session = True
                        _ray_session_dir = ctx.address_info["session_dir"]  # type: ignore[attr-defined]
                        atexit.register(_shutdown_owned_ray)

                else:
                    logger.info("Available Ray cluster already exists; connecting to it")

                # Put these (somewhat) large objects in shared memory
                density_estimator_ref = ray.put(density_estimator)
                points_ref = ray.put(points)

                # Launch parallel tasks
                logger.info(f"Dispatching {n_resamplings} EOS evidence integration tasks to Ray cluster")
                results = []
                for _ in range(n_resamplings):
                    results.append(
                        _distributed_eos_evidence_integration.options(enable_task_events=False).remote(density_estimator_ref, points_ref, event_type)
                    )

                evidences[1:] = ray.get(results)

                logger.info("Ray tasks returned")

            return evidences

        else:
            return original_evidence

    def _create_eos_contour(
        self,
        eos_interpolator: EOSInterpolator,
        event_type: str,
        integration_bounds,
        mean_mchirp: float | None,
        n_grid: int = 100,
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
                    (q, mchirp, Lambda1, Lambda2) is returned (with
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

        if event_type == "gw-2d":
            q = np.linspace(integration_bounds[0], integration_bounds[1], n_grid)
            m1, m2 = convert_masses(q, mean_mchirp)
            m1, m2, q = eos_interpolator.apply_bns_mass_constraint(m1, m2, q)  # type: ignore

            lambdat = eos_interpolator.get_lambda_tilde(m1, m2)
            points = np.stack((lambdat, q), axis=-1)

        elif event_type == "gw-3d":
            q = np.linspace(integration_bounds[0], integration_bounds[1], n_grid)
            m1, m2 = convert_masses(q, mean_mchirp)
            m1, m2, q = eos_interpolator.apply_bns_mass_constraint(m1, m2, q)  # type: ignore

            lambda1 = eos_interpolator.get_lambda(m1)
            lambda2 = eos_interpolator.get_lambda(m2)
            points = np.stack((lambda1, q, lambda2), axis=-1)

        elif event_type == "gw-4d":
            q = np.linspace(integration_bounds[0][0], integration_bounds[0][1], n_grid)
            mchirp = np.linspace(integration_bounds[1][0], integration_bounds[1][1], n_grid)

            q_grid, mchirp_grid = np.meshgrid(q, mchirp)

            # make same size 2D grid in m1, m2 in order to compute
            # Lambda1 and Lambda2 grids
            m1, m2 = convert_masses(q, mchirp)
            m1_grid, m2_grid = np.meshgrid(m1, m2)

            lambda1 = eos_interpolator.get_lambda(m1_grid.reshape(n_grid**2)).reshape((n_grid, n_grid))
            lambda2 = eos_interpolator.get_lambda(m2_grid.reshape(n_grid**2)).reshape((n_grid, n_grid))

            points = np.stack((q_grid, mchirp_grid, lambda1, lambda2), axis=-1)

        elif event_type == "psr":
            mass = np.linspace(integration_bounds[0], integration_bounds[1], n_grid)
            mass = eos_interpolator.apply_ns_mass_constraint(mass)

            radius = eos_interpolator.get_radius(mass)
            points = np.stack((mass, radius), axis=-1)

        return points

    def _integrate_eos_contour(
        self, points: np.ndarray, density_estimator: BoundedKDE | EnsembleNormalizingFlow, event_type: str, resample: bool = False
    ) -> np.ndarray:
        """Compute the EOS evidence by integrating the EOS line
        characterized by ``points`` over the event posterior density.

        The integral is performed numerically using ``np.trapezoid``.

        Parameters
        ----------
        points
            Array of points as returned by ``ModelSelector._create_eos_contour``
        resample
            Whether to resample the density estimator before
            integrating, by default False

        Returns
        -------
            Single-element array containing the evidence value
        """

        if isinstance(density_estimator, BoundedKDE):
            if event_type == "gw-4d":
                q_arr = points[0, :, 0]
                mchirp_arr = points[:, 0, 1]

                n_grid = points.shape[0]
                points = points.reshape(n_grid**2, 4)

                density = density_estimator.pdf(points, resample)
                density = density.reshape((n_grid, n_grid))

                integral_over_mchirp = np.trapezoid(density, mchirp_arr, axis=0)
                evidence = np.trapezoid(integral_over_mchirp, q_arr)

            else:
                density = density_estimator.pdf(points, resample)

                # integrate over: ...
                integrate_dim = {
                    "gw-2d": 1,  # q
                    "gw-3d": 1,  # q
                    "psr": 0,  # m
                }.get(event_type)
                evidence = np.trapezoid(y=density, x=points[:, integrate_dim])

            return evidence

        else:
            if event_type == "gw-4d":
                q_arr = points[0, :, 0]
                mchirp_arr = points[:, 0, 1]

                n_grid = points.shape[0]
                points = points.reshape(n_grid**2, 4)

                ensemble_densities = density_estimator.pdf(points)

                ensemble_evidences = []
                for density in ensemble_densities:
                    density = density.reshape((n_grid, n_grid))

                    integral_over_mchirp = np.trapezoid(density, mchirp_arr, axis=0)
                    ensemble_evidences.append(np.trapezoid(integral_over_mchirp, q_arr))

            else:
                ensemble_densities = density_estimator.pdf(points)

                # integrate over: ...
                integrate_dim = {
                    "gw-2d": 1,  # q
                    "gw-3d": 1,  # q
                    "psr": 0,  # m
                }.get(event_type)

                ensemble_evidences = []
                for density in ensemble_densities:
                    ensemble_evidences.append(np.trapezoid(y=density, x=points[:, integrate_dim]))

            return np.array(ensemble_evidences)

    def parameterized_eos_evidence(self, parameters, parameterization: Literal["spectral", "polytrope"]) -> float:
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

        Returns
        -------
            Evidence of the EOS constructed from the given parameters and parameterization
        """

        interpolator = EOSInterpolator(
            eos_parameters=np.array(parameters),
            parameterization=parameterization,
        )

        per_event_evidences = []
        for i in range(len(self.density_estimators)):
            points = self._create_eos_contour(
                interpolator,
                self.event_types[i],
                self.integration_bounds[i],
                self.mean_mchirps[i],
                n_grid=self.parameterized_evidence_integral_n_grids[i],
            )

            evidence = self._integrate_eos_contour(points, self.density_estimators[i], self.event_types[i])

            if isinstance(self.density_estimators[i], EnsembleNormalizingFlow):
                evidence = np.mean(evidence).item()

            per_event_evidences.append(evidence)

        joint_evidence = np.prod(per_event_evidences)
        return float(joint_evidence)


class ParameterizedEoSSampler:
    """Performs MCMC parameter estimation for a parameterized EOS model using ``emcee``.

    The sampling likelihood is the joint evidence of the EOS from several provided events,
    as computed by ``ModelSelector.parameterized_eos_evidence``.
    """

    def __init__(
        self,
        posterior_files: Sequence[str],
        event_types: Sequence[str],
        eos_prior_bounds: Sequence[tuple],
        largest_observed_ns_mass: float = 1.97,
        density_est_method: Literal["kde", "flow"] = "kde",
        flow_directories: Sequence[str] | None = None,
        parameterization: Literal["spectral", "polytrope"] = "spectral",
        integration_bounds: Sequence[tuple[float, float]] | Sequence[Sequence[tuple[float, float]]] | None = None,
    ):
        """
        Parameters
        ----------
        posterior_files
            Paths to .json, .h5/.hdf5, or .dat files containing posterior samples of the necessary
            EOS-dependent parameters of each event.
            The needed parameter samples for each file depends on the corresponding ``event_type``:

            - "gw-2d" requires samples for mass ratio (q), chirp mass, and dominant tidal deformability (LambdaTilde)
            - "gw-3d" and "gw-4d" require samples for mass ratio (q), chirp mass, and the two tidal deformabilities (Lambda1, Lambda2)
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
            - "flow" : Set of normalizing flow PyTorch/Zuko models, trained on event data identically and only differing due to random weight initializations. This approach is designed to support a reproducible alternative to the Bayesian flow approach, with uncertainty estimation coming from the variance in density estimates from the ensemble of models. If chosen, ``flow-files`` must be passed.

        flow_directories
            If ``density_est_method`` is "flow", provide a list of paths (1 per event) to directories containing the ensembles of PyTorch/Zuko-based Normalizing Flow models in the form of .onnx files (1 per model).

        parameterization
            Must be one of "spectral" (4-parameter spectral decomposition model) or "polytrope" (4-parameter piecewise-polytrope model)

        integration_bounds
            Bounds for the EOS evidence integral.

            These should correspond to:

            - (q_min, q_max) for "gw-2d" and "gw-3d" events
            - [(q_min, q_max), (mchirp_min, mchirp_max)] for "gw-4d" events
            - (mass_min, mass_max) for "psr" events (with mass in solar masses)

            If None, the bounds on the respective integration parameter(s) will be taken from the min and max of the posterior samples.
        """

        self.eos_prior_bounds = eos_prior_bounds
        self.parameterization = parameterization
        self.largest_observed_ns_mass = largest_observed_ns_mass

        self.joint_selector = ModelSelector(
            posterior_files=posterior_files,
            event_types=event_types,
            density_est_method=density_est_method,
            flow_directories=flow_directories,
            integration_bounds=integration_bounds,
        )

    def _log_post(self, parameters: np.ndarray) -> float:
        """Compute the log posterior density (MCMC sampling likelihood) for
        the given ``parameters``.

        If the given parameter vector describes a valid EOS (i.e. adheres to
        uniform prior bounds given when instantiating this class, adheres to
        causality, thermodynamic stability, and observational consistency with
        the most massive observed neutron star) then the log likelihood given by
        this function is simply the log of the multi-event joint evidence for the EOS
        as computed by ``ModelSelector.parameterized_eos_evidence``.

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

        try:
            joint_evidence = self.joint_selector.parameterized_eos_evidence(parameters, self.parameterization)
        except RuntimeError as e:
            logger.error(f"RuntimeError in _log_post({parameters}): {e}")
            return -np.inf

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

        logger.info(f"Initializing {nwalkers} MCMC walkers. EOS prior bounds: {self.eos_prior_bounds}")

        n_valid_walkers = 0
        state0 = []

        param_lower_bounds = [bound[0] for bound in self.eos_prior_bounds]
        param_upper_bounds = [bound[1] for bound in self.eos_prior_bounds]

        while n_valid_walkers < nwalkers:
            params = np.random.uniform(param_lower_bounds, param_upper_bounds)

            if is_valid_eos(params, self.parameterization, largest_ns_mass=self.largest_observed_ns_mass):
                logger.info(f"Walker {n_valid_walkers + 1} initial state: {params}")

                state0.append(params)
                n_valid_walkers += 1

        return state0

    def run_sampler(self, nsteps: int, nwalkers: int, save_file: str, reset: bool = True, pool=None) -> None:
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
            Number of MCMC walkers to use. Recommended: 300.
        save_file
            Path to an hdf5 file that will be used to store the sample chain and the associated log
            probability densities. This file can be read after the run with the ``load_samples`` function.
        reset
            Whether to overwrite the samples currently stored in the given ``save_file``, such
            as from a previous run with the same data/method. If False, the sampler will continue
            from where the previous run left off, using the last sample in the stored chain as the
            initial state of the walkers. If True (default), anything in ``save_file`` will be
            overwritten with a fresh run. See https://emcee.readthedocs.io/en/stable/tutorials/monitor/
        pool
            A multiprocessing pool object that will be used by emcee to parallelize the sampling
            steps over the walkers.

            Example:

            ``with multiprocessing.Pool(processes=96) as pool:``
            ``    sampler.run_sampler(100, nwalkers=500, save_file=save_file, pool=pool)``
        """

        logger.info(f"Running MCMC for {self.parameterization} EOS with {nwalkers} walkers for {nsteps} steps")

        # Won't be able to create the file if parent dir doesn't exist beforehand
        if not pathlib.Path(save_file).parent.exists():
            raise FileNotFoundError(f"Parent directory of given save_file doesn't exist: {save_file}")

        logger.info(f"Saving samples and logp chains to {save_file}")
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
            logger.info(f"Continuing existing MCMC chain from {save_file}. Initial state: {backend.get_chain()[-1]}")

        sampler = emcee.EnsembleSampler(
            nwalkers=nwalkers,
            ndim=4,
            log_prob_fn=self._log_post,
            backend=backend,
            moves=[
                (emcee.moves.StretchMove(), 0.2),
                (emcee.moves.DEMove(), 0.6),
                (emcee.moves.DESnookerMove(), 0.2),
            ],
            pool=pool,
        )

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
            converged = np.all(tau * 100 < (sampler.iteration * nwalkers)) and np.all(np.abs(old_tau - tau) / tau < 0.01)
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


def get_posterior_samples(posterior_file: str, event_type: str) -> dict:
    """Convenience utility to read a posterior file and return
    a stacked array containing the necessary parameter samples
    for GWXtreme inference based on the variant of the approximation
    specified by ``cbc_dim``.

    Parameters
    ----------
    posterior_file
        Required contents and/or keys depend on given ``event_type``.
        If ``event_type`` is "gw-2d", the file must include:

            q, mc_source, lambdat

            Possible alternative names:
            mass_ratio, chirp_mass_source, lambda_tilde

            Must be in one of these formats: .h5/.hdf5, .json, .dat

        If ``event_type`` is "gw-3d" or "gw-4d", the file must include:

            q, mc_source, lambda_1, lambda_2

            Possible alternative names:
            mass_ratio, chirp_mass_source, lambda1, lambda2

            Must be in one of these formats: .h5/.hdf5, .json, .dat

        If ``event_type`` is "psr", the file must be a text file
        containing two columns: mass in solar masses
        and compactness.

    event_type
        "gw-2d", "gw-3d", "gw-4d", or "psr"

    Returns
    -------
        Contents of returned dict depend on ``event_type``:
        "gw-2d": {"q": ..., "lambdat": ...}
        "gw-3d" or "gw-4d": {"q": ..., "mchirp": ..., "lambda1": ..., "lambda2": ...}
        "psr": {"mass": ..., "radius": ...}
    """

    if event_type == "psr":
        samples = np.loadtxt(posterior_file, dtype=np.float32)

        # Convert compactness to radius in km
        mass_in_sm, compactness = samples[:, 0], samples[:, 1]
        mass_in_kg = lal.MSUN_SI * mass_in_sm
        radius_in_km = (lal.G_SI * mass_in_kg / (lal.C_SI**2 * compactness)) / 1000

        return {"mass": mass_in_sm, "radius": radius_in_km}

    else:
        ext = pathlib.Path(posterior_file).suffix

        if ext == ".h5" or ext == ".hdf5":
            with h5py.File(posterior_file) as f:
                samples = np.array(f["posterior_samples"])

        elif ext == ".json":
            with open(posterior_file) as f:
                samples = json.load(f)["posterior"]["content"]

        elif ext == ".dat":
            samples = np.genfromtxt(posterior_file, names=True)

        try:
            q = np.array(samples["q"])
        except KeyError:
            q = np.array(samples["mass_ratio"])

        try:
            mchirp = np.array(samples["mc_source"])
        except KeyError:
            mchirp = np.array(samples["chirp_mass_source"])

        if event_type == "gw-2d":
            try:
                lambdat = np.array(samples["lambdat"])
            except KeyError:
                lambdat = np.array(samples["lambda_tilde"])

            return {"q": q, "lambdat": lambdat, "mchirp": mchirp}

        elif event_type in ["gw-3d", "gw-4d"]:
            try:
                lambda1 = np.array(samples["lambda_1"])
                lambda2 = np.array(samples["lambda_2"])
            except KeyError:
                lambda1 = np.array(samples["lambda1"])
                lambda2 = np.array(samples["lambda2"])

            return {"q": q, "lambda1": lambda1, "lambda2": lambda2, "mchirp": mchirp}

        else:
            raise NotImplementedError()
