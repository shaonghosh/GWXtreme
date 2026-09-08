import logging
import pathlib
import time
from collections.abc import Sequence
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import onnxruntime
import torch
import zuko
import optuna

from gwxtreme.density_estimation import get_log_abs_det_jacobian, to_latent_space
from gwxtreme.flow_density_estimation import EnsembleNormalizingFlow
from gwxtreme.flow_training import density_corner_plot, evaluate_over_grid, flow_training_checkpoint_to_onnx, train_flow
from gwxtreme.utils import get_gw_event_pe_posterior_samples, get_nicer_pulsar_pe_posterior_samples
from flow_training import sample_optimizer_kwargs, sample_flow_kwargs

# By default, PyTorch parallelizes tensor operations and sets of operations across several threads.
# This can be problematic. E.g. For a machine with 128 logical CPUs / hardware threads, PyTorch tries
# by default to use 64 and 128 threads for intra-op and inter-op parallelization, respectively.
# Changing these settings to a smaller number has resulted in a ~ 16 X speedup, as allowing the
# default use of too many threads results in CPU oversubscription.
torch.set_num_threads(max(2, int(torch.get_num_threads() / 8)))
torch.set_num_interop_threads(max(4, int(torch.get_num_interop_threads() / 8)))


class TransformBoundedFlow:
    def __init__(
        self,
        flow_file: str,
        bounds: Sequence[tuple[float, float]],
    ):
        self.bounds = bounds

        options = onnxruntime.SessionOptions()
        options.intra_op_num_threads = 1
        options.inter_op_num_threads = 1
        providers = ["CPUExecutionProvider"]

        self.session = onnxruntime.InferenceSession(flow_file, options, providers=providers)

    def log_pdf(self, x: np.ndarray) -> np.ndarray:
        z = to_latent_space(x, self.bounds)
        ladj = get_log_abs_det_jacobian(x, self.bounds)

        finite_valued_points = np.all(np.isfinite(z), axis=1)

        lp = self.session.run([self.session.get_outputs()[0].name], {self.session.get_inputs()[0].name: z})[0] + ladj
        lp[~finite_valued_points] = -np.inf
        return lp

    def pdf(self, x: np.ndarray) -> np.ndarray:
        return np.exp(self.log_pdf(x))


def single_flow_testing(outdir,run_optuna=False,n_trials=5,optimizer_param_space=None):
    #########################################
    data_dim = 2
    posterior_file = "/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/cbc_pe_samples/GW170817/posterior_samples/GW170817_posterior_samples_TaylorF2_narrow_spin_uniform_lambda_tildes_prior.dat"
    flow_class = zuko.flows.MAF
    flow_kwargs = {"features": data_dim, "transforms": 5}
    optimizer_class = torch.optim.Adam
    optimizer_kwargs = {"lr": 5e-4, "betas": (0.9, 0.999)}
    #########################################

    # Bounds on parameters for 2D, 3D, 4D events
    bounds = {
        2: [(0.0, np.inf), (0.0, 1.0)],  # LambdaTilde, q
        3: [(0.0, np.inf), (0.0, 1.0), (0.0, np.inf)],  # Lambda1, q, Lambda2
        4: [(0.0, 1.0), (0.0, np.inf), (0.0, np.inf), (0.0, np.inf)],  # q, m_chirp, Lambda1, Lambda2
    }

    # Ranges of values of parameters to use for plotting/testing
    test_bounds = {
        2: [(0.0, 2500.0), (0.0, 1.0)],
        3: [(0.0, 2500.0), (0.0, 1.0), (0.0, 2500.0)],
        4: [(0.0, 1.0), (0.8, 2.0), (0.0, 2500.0), (0.0, 2500.0)],
    }

    parameter_labels = {
        2: [r"$\tilde{\Lambda}$", "$q$"],
        3: [r"$\Lambda_1$", "$q$", r"$\Lambda_2$"],
        4: ["$q$", r"\mathcal{M}", r"$\Lambda_1$", r"$\Lambda_2$"],
    }
    default_optuna_params = {
        "lr" : {"low": 1e-5, "high": 1e-4},
        "beta1": {"low": 0.80, "high": 0.999},
        "beta2": {"low": 0.90, "high": 0.9999},
        "transforms": {"low": 3, "high": 8}
    }


    data = get_gw_event_pe_posterior_samples(posterior_file, cbc_dim=data_dim)
    latent_data = to_latent_space(data, bounds[data_dim])
    latent_data = torch.from_numpy(latent_data)

    # can add other labeling here to distinguish runs
    run_label = f"{flow_class.__name__}_{data_dim}D"  # MAF_2D
    savedir = pathlib.Path(outdir) / run_label

    if run_optuna:
        param_space = optimizer_param_space or default_optuna_params
        def objective(trial):
            trial_kwargs = sample_optimizer_kwargs(trial, optimizer_kwargs, param_space)
            trial_flow_kwargs = sample_flow_kwargs(trial, flow_kwargs, param_space)
            trial_save_dir = savedir / f"trial_{trial.number}"
            filename_prefix = f"trial_{trial.number}_"
            train_losses, val_losses = train_flow(
                zuko_flow_class=flow_class,
                flow_kwargs=trial_flow_kwargs,
                optimizer_class=optimizer_class,
                optimizer_kwargs=trial_kwargs,
                data=latent_data,
                savedir=str(trial_save_dir),
                n_epochs=20,
                batch_size=200,
                max_norm=25.0,
                early_stopping=10,
                filename_prefix = filename_prefix
                
            )

            trial_ckpt_path = trial_save_dir / f"{filename_prefix}best_checkpoint.pth"
            trial_checkpoint = torch.load(trial_ckpt_path)
            
            #  # The best_checkpoint.pth file is created by the train_flow function in the given savedir.
            # checkpoint_file = savedir / "best_checkpoint.pth"
            onnx_file = trial_ckpt_path.with_suffix(".onnx")
            flow_training_checkpoint_to_onnx(str(trial_ckpt_path), str(onnx_file))
        
            flow = TransformBoundedFlow(
                flow_file=str(onnx_file),
                bounds=bounds[data_dim],
            )
        
            density, grids = evaluate_over_grid(
                f=flow.pdf,
                grid_bounds=test_bounds[data_dim],
                grid_size=100,  # this can be increased for 2D and maybe 3D, but not 4D events
            )
        
            fig, axes = density_corner_plot(density=density, grids=grids, samples=data, parameter_labels=parameter_labels[data_dim])
            plt.close()
            fig.savefig(str(trial_save_dir / f"{filename_prefix}data_density.png"))

        
            return min(val_losses)  



        savedir.mkdir(parents=True, exist_ok=True)
        study = optuna.create_study(direction = "minimize", study_name = run_label,
                                   storage=f"sqlite:///{savedir / 'optuna_study.db'}", load_if_exists=True)
        study.optimize(objective, n_trials=n_trials)
        best_ckpt_path = savedir / f"trial_{study.best_trial.number}" / f"trial_{study.best_trial.number}_best_checkpoint.pth"
        best_checkpoint = torch.load(best_ckpt_path)
        print(f"best trial #{study.best_trial.number}, "
              f"val_loss={study.best_value:.6f}, params={study.best_params}")


        
    else:
        train_flow(
            zuko_flow_class=flow_class,
            flow_kwargs=flow_kwargs,
            optimizer_class=optimizer_class,
            optimizer_kwargs=optimizer_kwargs,
            data=latent_data,
            savedir=str(savedir),
            n_epochs=20,
            batch_size=200,
            max_norm=25.0,  # maximum norm of the gradient for each training parameter, it will be clipped to this value if above
            early_stopping=10,  # epochs with no decrease in validation loss to continue before quitting
        )

        #The best_checkpoint.pth file is created by the train_flow function in the given savedir.
        best_ckpt_path = savedir / "best_checkpoint.pth"  # Path throughout
        best_checkpoint = torch.load(best_ckpt_path)
    
        onnx_file = best_ckpt_path.with_suffix(".onnx")   # works now — best_ckpt_path is a real Path
        flow_training_checkpoint_to_onnx(str(best_ckpt_path), str(onnx_file))

        ############### Load and evaluate ###############
        flow = TransformBoundedFlow(
            flow_file=str(onnx_file),
            bounds=bounds[data_dim],
        )
    
        # This uses a helper function just to evaluate a given function (here, flow.pdf())
        # over a grid defined by given grid_bounds. The density grid can then be plotted easily.
        density, grids = evaluate_over_grid(
            f=flow.pdf,
            grid_bounds=test_bounds[data_dim],
            grid_size=100,  # this can be increased for 2D and maybe 3D, but not 4D events
        )
    
        fig, axes = density_corner_plot(density=density, grids=grids, samples=data, parameter_labels=parameter_labels[data_dim])
        fig.savefig(str(savedir / "data_density.png"))
