import logging
import pathlib
import time
from collections import Counter
from typing import Literal

import numpy as np
import torch
import tqdm
import zuko
import zuko.bayesian

from gwxtreme.density_estimation import to_latent_space
from gwxtreme.utils import get_gw_event_pe_posterior_samples, get_nicer_pulsar_pe_posterior_samples

logger = logging.getLogger(__name__)


def train_masked_autoregressive_flow(
    data: torch.Tensor,
    n_epochs: int,
    batch_size: int,
    model_label: str,
    save_dir: str,
    max_norm=100.0,
    n_transforms=3,
    num_epochs_no_improvement_before_stop: int = 10,
) -> list:
    assert len(data.shape) == 2  # (N, n_dim)
    assert len(data) >= batch_size, f"need at least {batch_size} samples, got {len(data)}"
    assert pathlib.Path(save_dir).exists(), "given save_dir directory doesn't exist"

    flow = zuko.flows.MAF(features=data.shape[-1], transforms=n_transforms, randperm=False)

    optimizer = torch.optim.Adam(params=flow.parameters(), lr=1e-3, betas=(0.9, 0.999))
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=0.5,
        patience=7,
        min_lr=1e-8,
    )

    train_dataset, test_dataset = torch.utils.data.random_split(data, [0.80, 0.20])

    train_loader = torch.utils.data.DataLoader(dataset=train_dataset, drop_last=True, batch_size=batch_size, shuffle=True)

    model_save_path = pathlib.Path(save_dir).joinpath(model_label)

    logger.info(
        f"data.shape: {data.shape}\nflow: {flow}\noptimizer: {optimizer}\nN_epochs: {n_epochs}\nbatch_size: {batch_size}\nmax_norm: {max_norm}"
    )

    epoch_mean_losses = []
    minimum_epoch_mean_loss = torch.inf
    best_epoch = -1

    # Set model to train mode
    flow.train()

    start = time.perf_counter()
    for epoch in tqdm.tqdm(range(n_epochs)):
        # Collect losses per train batch for this epoch
        losses = []

        for d in train_loader:
            batch = d[0]

            loss = -flow().log_prob(batch).mean()
            loss.backward()

            preclip_norm = torch.nn.utils.clip_grad_norm_(flow.parameters(), max_norm=max_norm)
            if preclip_norm > max_norm:
                logger.warning(f"[{epoch:6d} / {n_epochs}]\tPre-clip norm exceeded clip threshold: {preclip_norm:.2f} > {max_norm}")

            optimizer.step()
            optimizer.zero_grad(set_to_none=True)
            losses.append(loss.detach())

        # Show top 10 largest gradients
        by_param = Counter({name: p.grad.norm().item() for name, p in flow.named_parameters() if p.grad is not None})
        for name, norm in sorted(by_param.items(), key=lambda x: -x[1])[:10]:
            logger.info(f"[{epoch:6d} / {n_epochs}]\tPer-parameter gradient norms\t{name}: {norm:.2f}")

        epoch_mean_loss = torch.stack(losses).mean().item()
        epoch_mean_losses.append(epoch_mean_loss)

        logger.info(f"[{epoch:6d} / {n_epochs}]\tavg. loss = {epoch_mean_loss:3.4f}")

        prev_lr = optimizer.param_groups[0]["lr"]
        scheduler.step(epoch_mean_loss)
        current_lr = optimizer.param_groups[0]["lr"]
        if current_lr != prev_lr:
            logger.info(f"[{epoch:6d} / {n_epochs}]\tLR reduced: {prev_lr:.2e} -> {current_lr:.2e}")

        if epoch_mean_loss < minimum_epoch_mean_loss:
            minimum_epoch_mean_loss = epoch_mean_loss
            best_epoch = epoch
            torch.save(
                {
                    "state_dict": flow.state_dict(),
                    "flow_config": {"features": data.shape[-1], "transforms": n_transforms, "randperm": False},
                },
                model_save_path.with_suffix(".pt"),
            )
        else:
            if epoch - best_epoch >= num_epochs_no_improvement_before_stop:
                logger.info(f"Stopping early - no improvement in loss after {num_epochs_no_improvement_before_stop} epochs.")
                break

    end = time.perf_counter()

    flow.eval()
    avg_log_likelihood_on_trainset = flow().log_prob(train_dataset.dataset[train_dataset.indices]).mean().item()
    avg_log_likelihood_on_testset = flow().log_prob(test_dataset.dataset[test_dataset.indices]).mean().item()

    logger.info(f"Training finished: time = {(end - start) / 60:.2f} minutes")
    logger.info(f"Avg LL (train set) = {avg_log_likelihood_on_trainset:.3f}\nAvg LL (test set) = {avg_log_likelihood_on_testset:.3f}")

    return epoch_mean_losses


def train_bayesian_masked_autoregressive_flow(
    data: torch.Tensor,
    n_epochs: int,
    batch_size: int,
    model_label: str,
    save_dir: str,
    init_logvar: float = -6.0,
    max_norm=100.0,
    n_transforms=3,
    num_epochs_no_improvement_before_stop: int = 10,
) -> list:
    assert len(data.shape) == 2  # (N, n_dim)
    assert len(data) >= batch_size, f"need at least {batch_size} samples, got {len(data)}"
    assert pathlib.Path(save_dir).exists(), "given save_dir directory doesn't exist"

    flow = zuko.flows.MAF(features=data.shape[-1], transforms=n_transforms, randperm=False)

    include_params = ["transform.transforms.*.hyper.4.weight"]
    exclude_params = ["**.bias"]
    bayes_flow = zuko.bayesian.BayesianModel(flow, init_logvar=init_logvar, include_params=include_params, exclude_params=exclude_params)

    optimizer = torch.optim.Adam(params=bayes_flow.parameters(), lr=1e-3, betas=(0.9, 0.999))
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=0.5,
        patience=7,
        min_lr=1e-8,
    )

    train_dataset, test_dataset = torch.utils.data.random_split(data, [0.80, 0.20])

    train_loader = torch.utils.data.DataLoader(dataset=train_dataset, drop_last=True, batch_size=batch_size, shuffle=True)

    model_save_path = pathlib.Path(save_dir).joinpath(model_label)

    logger.info(
        f"data.shape: {data.shape}\nflow: {flow}\noptimizer: {optimizer}\nN_epochs: {n_epochs}\nbatch_size: {batch_size}\ninit_logvar: {init_logvar}\nmax_norm: {max_norm}"
    )

    # Bayesian model parameters
    logger.info(f"{len(bayes_flow.means)} variational parameter groups:")
    for key in bayes_flow.means:
        logger.info(f"  {key.replace('-', '.')}  shape={tuple(bayes_flow.means[key].shape)}")

    epoch_mean_losses = []
    minimum_epoch_mean_loss = torch.inf
    best_epoch = -1

    # Set model to train mode
    bayes_flow.train()

    start = time.perf_counter()
    for epoch in tqdm.tqdm(range(n_epochs)):
        # Collect losses per train batch for this epoch
        losses = []

        for d in train_loader:
            batch = d[0]
            kl_per_sample = bayes_flow.kl_divergence() / len(data)

            with bayes_flow.reparameterize() as flow_rep:
                nll = -flow_rep().log_prob(batch).mean()
                loss = nll + kl_per_sample
                loss.backward()

            preclip_norm = torch.nn.utils.clip_grad_norm_(bayes_flow.parameters(), max_norm=max_norm)
            if preclip_norm > max_norm:
                logger.warning(f"[{epoch:6d} / {n_epochs}]\tPre-clip norm exceeded clip threshold: {preclip_norm:.2f} > {max_norm}")

            optimizer.step()
            optimizer.zero_grad(set_to_none=True)
            losses.append(loss.detach())

        # Show top 10 largest gradients
        by_param = Counter({name: p.grad.norm().item() for name, p in bayes_flow.named_parameters() if p.grad is not None})
        for name, norm in sorted(by_param.items(), key=lambda x: -x[1])[:10]:
            logger.info(f"[{epoch:6d} / {n_epochs}]\tPer-parameter gradient norms\t{name}: {norm:.2f}")

        epoch_mean_loss = torch.stack(losses).mean().item()
        epoch_mean_losses.append(epoch_mean_loss)

        logger.info(f"[{epoch:6d} / {n_epochs}]\tavg. loss = {epoch_mean_loss:3.4f}")

        prev_lr = optimizer.param_groups[0]["lr"]
        scheduler.step(epoch_mean_loss)
        current_lr = optimizer.param_groups[0]["lr"]
        if current_lr != prev_lr:
            logger.info(f"[{epoch:6d} / {n_epochs}]\tLR reduced: {prev_lr:.2e} -> {current_lr:.2e}")

        if epoch_mean_loss < minimum_epoch_mean_loss:
            minimum_epoch_mean_loss = epoch_mean_loss
            best_epoch = epoch
            torch.save(
                {
                    "state_dict": bayes_flow.state_dict(),
                    "flow_config": {"features": data.shape[-1], "transforms": n_transforms, "randperm": False},
                    "init_logvar": init_logvar,
                    "include_params": include_params,
                    "exclude_params": exclude_params,
                },
                model_save_path.with_suffix(".pt"),
            )
        else:
            if epoch - best_epoch >= num_epochs_no_improvement_before_stop:
                logger.info(f"Stopping early - no improvement in loss after {num_epochs_no_improvement_before_stop} epochs.")
                break

    end = time.perf_counter()

    bayes_flow.eval()
    flow_instance = bayes_flow.sample_model()
    avg_log_likelihood_on_trainset = flow_instance().log_prob(train_dataset.dataset[train_dataset.indices]).mean().item()
    avg_log_likelihood_on_testset = flow_instance().log_prob(test_dataset.dataset[test_dataset.indices]).mean().item()

    logger.info(f"Training finished: time = {(end - start) / 60:.2f} minutes")
    logger.info(f"Avg LL (train set) = {avg_log_likelihood_on_trainset:.3f}\nAvg LL (test set) = {avg_log_likelihood_on_testset:.3f}")

    return epoch_mean_losses


def train_flow_for_event(
    posterior_file: str, event_type: Literal["gw-2d", "gw-3d", "gw-4d", "psr"], model_label: str, save_dir: str, bayesian: bool = True
):
    if event_type == "gw-2d":
        posterior_samples = get_gw_event_pe_posterior_samples(posterior_file, cbc_dim=2)
        parameter_bounds = [(0.0, np.inf), (0.0, 1.0)]

    elif event_type == "gw-3d":
        posterior_samples = get_gw_event_pe_posterior_samples(posterior_file, cbc_dim=3)
        parameter_bounds = [(0.0, np.inf), (0.0, 1.0), (0.0, np.inf)]

    elif event_type == "gw-4d":
        posterior_samples = get_gw_event_pe_posterior_samples(posterior_file, cbc_dim=4)
        parameter_bounds = [
            (0.0, 1.0),
            (0.0, np.inf),
            (0.0, np.inf),
            (0.0, np.inf),
        ]

    elif event_type == "psr":
        posterior_samples = get_nicer_pulsar_pe_posterior_samples(posterior_file)
        parameter_bounds = [(0.0, np.inf), (0.0, np.inf)]

    latent_space_samples = to_latent_space(posterior_samples, parameter_bounds)

    if bayesian:
        train_bayesian_masked_autoregressive_flow(
            data=torch.tensor(latent_space_samples, dtype=torch.float32),
            n_epochs=500,
            batch_size=100,
            model_label=model_label,
            save_dir=save_dir,
            init_logvar=-6.0,
            max_norm=100.0,
            n_transforms=15,
            num_epochs_no_improvement_before_stop=75,
        )
    else:
        train_masked_autoregressive_flow(
            data=torch.tensor(latent_space_samples, dtype=torch.float32),
            n_epochs=500,
            batch_size=100,
            model_label=model_label,
            save_dir=save_dir,
            max_norm=100.0,
            n_transforms=15,
            num_epochs_no_improvement_before_stop=75,
        )
