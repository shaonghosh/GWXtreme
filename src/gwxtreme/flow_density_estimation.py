import pathlib

import numpy as np
import torch
import zuko
import zuko.bayesian

from gwxtreme.density_estimation import get_log_abs_det_jacobian, to_data_space, to_latent_space

# By default, PyTorch parallelizes tensor operations and sets of operations across several threads.
# For our machine with 128 logical CPUs / hardware threads, PyTorch tries by default to use 64 and
# 128 threads for intra-op and inter-op parallelization, respectively.
# Changing these settings to a smaller number has resulted in a ~ 16 X speedup, as allowing the
# default use of too many threads results in CPU oversubscription.
torch.set_num_threads(max(2, int(torch.get_num_threads() / 8)))
torch.set_num_interop_threads(max(4, int(torch.get_num_interop_threads() / 8)))


class EnsembleNormalizingFlow:
    """Density estimator constructed from an ensemble of pre-trained normalizing flow models."""

    def __init__(
        self,
        posterior_samples: np.ndarray,
        bounds: list[tuple[float, float]],
        flows_dir: str,
    ):
        """
        Parameters
        ----------
        posterior_samples
            Array of samples (with shape (N_samples, N_dim)) from the distribution on which density estimation is being performed
        bounds
            List of tuples of the form (lower bound, upper bound) for each of the N_dim parameters appearing in posterior_samples.
            For infinite bounds (unbounded), use -np.inf or np.inf.
        flows_dir
            Path to a directory containing .pt files containg the weights (state dicts) and configurations for a set of trained PyTorch/Zuko-based MAF models.
        """

        self.posterior_samples = posterior_samples
        self.bounds = bounds
        self.flows = load_flow_ensemble(flows_dir)

    def log_pdf(self, x: np.ndarray) -> list[np.ndarray]:
        """Compute log probability densities of individual points.

        Parameters
        ----------
        x
            Array of points in the data space with shape (N_points, N_dim)

        Returns
        -------
            List of arrays, 1 for each flow in the ensemble. Each array contains the log densities for that model with shape (N_points,).
        """

        z = to_latent_space(x, self.bounds)
        ladj = get_log_abs_det_jacobian(x, self.bounds)
        lps = []
        for flow in self.flows:
            with torch.no_grad():
                lp = flow().log_prob(torch.tensor(z, dtype=torch.float32)).detach().numpy() + ladj
            lp = np.nan_to_num(lp, nan=-np.inf)
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
            List of arrays, 1 for each flow in the ensemble. Each array contains the densities for that model with shape (N_points,).
        """

        p = [np.exp(lp) for lp in self.log_pdf(x)]
        return p


class BayesianNormalizingFlow:
    """Density estimator constructed from a pre-trained Bayesian Normalizing Flow model."""

    def __init__(
        self,
        posterior_samples: np.ndarray,
        bounds: list[tuple[float, float]],
        flow_file: str,
    ):
        """
        Parameters
        ----------
        posterior_samples
            Array of samples (with shape (N_samples, N_dim)) from the distribution on which density estimation is being performed
        bounds
            List of tuples of the form (lower bound, upper bound) for each of the N_dim parameters appearing in posterior_samples.
            For infinite bounds (unbounded), use -np.inf or np.inf.
        flow_file
            Path to a .pt file containg the weights and configuration for a PyTorch/Zuko-based Bayesian Normalizing Flow model.
        """

        self.posterior_samples = posterior_samples
        self.bounds = bounds
        self.bayesian_flow = load_bayesian_flow(flow_file)
        self.bayesian_flow.eval()
        self.flow = self.bayesian_flow.sample_model()

    def log_pdf(self, x: np.ndarray, resample: bool = False) -> np.ndarray:
        """Compute log probability densities of individual points.

        Parameters
        ----------
        x
            Array of points in the data space with shape (N_points, N_dim)
        resample
            Whether to sample a new instance of the Bayesian flow and use the
            resampled instance for this computation. Sampling a new instance
            means drawing a new set of model weights from the posterior distribution
            of weights that was inferred during training.
            Default: False

        Returns
        -------
            Array of log probability densities with shape (N_points,).
        """

        if resample:
            flow = self.bayesian_flow.sample_model()
        else:
            flow = self.flow

        z = to_latent_space(x, self.bounds)
        ladj = get_log_abs_det_jacobian(x, self.bounds)
        with torch.no_grad():
            lp = flow().log_prob(torch.tensor(z, dtype=torch.float32)).detach().numpy() + ladj
        lp = np.nan_to_num(lp, nan=-np.inf)
        return lp

    def pdf(self, x: np.ndarray, resample: bool = False) -> np.ndarray:
        """Compute probability densities of individual points.

        Parameters
        ----------
        x
            Array of points in the data space with shape (N_points, N_dim)
        resample
            Whether to sample a new instance of the Bayesian flow and use the
            resampled instance for this computation. Sampling a new instance
            means drawing a new set of model weights from the posterior distribution
            of weights that was inferred during training.
            Default: False

        Returns
        -------
            Array of probability densities with shape (N_points,)
        """

        return np.exp(self.log_pdf(x, resample=resample))


def load_bayesian_flow(path: str) -> zuko.bayesian.BayesianModel:
    """Load a stored ``zuko.bayesian.BayesianModel``
    wrapping a ``zuko.flows.MAF`` masked autoregressive
    flow model, to be used for density estimation.

    Parameters
    ----------
    path
        .pt file containing the trained model weights along with
        the configuration arguments that were used to instantiate
        the stored BayesianModel and MAF objects before training.
        Must contain the following arguments in this schema:

        {
        "state_dict": <PyTorch state dict of the model>,
        "flow_config": {
        "features": <dimensionality of the density estimator>,
        "transforms": <number of layers of the MAF>,
        "randperm": <whether the parameter order is permuted with each layer>
        }
        "init_logvar": <argument passed to BayesianModel upon creation>,
        "include_params": <argument passed to BayesianModel upon creation>,
        "exclude_params": <argument passed to BayesianModel upon creation>,
        }

    Returns
    -------
        The Zuko Bayesian MAF model
    """

    ckpt = torch.load(path, map_location="cpu", weights_only=False)

    flow = zuko.flows.MAF(**ckpt["flow_config"])
    bayes_flow = zuko.bayesian.BayesianModel(
        flow,
        init_logvar=ckpt["init_logvar"],
        include_params=ckpt["include_params"],
        exclude_params=ckpt["exclude_params"],
    )
    bayes_flow.load_state_dict(ckpt["state_dict"])
    bayes_flow.eval()
    return bayes_flow


def load_flow_ensemble(flows_dir: str) -> list[zuko.flows.MAF]:
    """Load an ensemble of trained ``zuko.flows.MAF`` masked
    autoregressive flow models from the given directory on disk.

    Parameters
    ----------
    flows_dir
        Path to directory containing the ensemble of models,
        each in a .pt file containing the model configuration
        arguments and the state dict (trained weights).
        Must contain the following arguments in this schema:

        {
        "state_dict": <PyTorch state dict of the model>,
        "flow_config": {
        "features": <dimensionality of the density estimator>,
        "transforms": <number of layers of the MAF>,
        "randperm": <whether the parameter order is permuted with each layer>
        }
        }

    Returns
    -------
        List of the Zuko MAF models
    """

    flows = []
    for flow_file in pathlib.Path(flows_dir).iterdir():
        ckpt = torch.load(flow_file, map_location="cpu", weights_only=False)
        flow = zuko.flows.MAF(**ckpt["flow_config"])
        flow.load_state_dict(ckpt["state_dict"])
        flow.eval()
        flows.append(flow)
    return flows
