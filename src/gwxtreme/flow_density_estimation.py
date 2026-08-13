import pathlib

import numpy as np
import onnxruntime

from gwxtreme.density_estimation import get_log_abs_det_jacobian, to_latent_space


class EnsembleNormalizingFlow:
    """Density estimator constructed from an ensemble of trained normalizing flow models."""

    def __init__(
        self,
        flows_dir: str,
        bounds: list[tuple[float, float]],
    ):
        """
        Parameters
        ----------
        flows_dir
            Path to a directory containing .onnx files, one for each flow included in the ensemble
        bounds
            List of tuples of the form (lower bound, upper bound) for each of the N_dim parameters appearing in posterior_samples.
            For infinite bounds (unbounded), use -np.inf or np.inf.
        """

        self.bounds = bounds

        options = onnxruntime.SessionOptions()
        options.intra_op_num_threads = 1
        options.inter_op_num_threads = 1
        providers = ["CPUExecutionProvider"]

        self.sessions = []
        for flow_file in pathlib.Path(flows_dir).iterdir():
            if flow_file.suffix in [".onnx"]:
                self.sessions.append(onnxruntime.InferenceSession(flow_file, options, providers=providers))

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
        for flow in self.sessions:
            lp = flow.run([flow.get_outputs()[0].name], {flow.get_inputs()[0].name: z})[0]
            lp += ladj
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
