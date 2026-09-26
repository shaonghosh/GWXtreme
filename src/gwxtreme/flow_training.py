import logging
import pathlib
import time
from collections import Counter
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import torch
import tqdm
import zuko

logger = logging.getLogger(__name__)


class FlowSavingWrapper(torch.nn.Module):
    """Wraps the trained flow object so that it can be saved
    to ONNX appropriately.

    Torch models saved to ONNX will have inputs passed to their
    forward(x) methods, but for the Zuko flows, the density evaluation
    call is done through flow().log_prob(x). This class converts the latter
    into a simple forward() call.
    """

    def __init__(self, flow):
        super().__init__()
        self.flow = flow

    def forward(self, z):
        return self.flow().log_prob(z)


def flow_training_checkpoint_to_onnx(checkpoint_file: str, savefile: str):
    ckpt = torch.load(checkpoint_file, map_location="cpu", weights_only=False)

    flow_class = getattr(zuko.flows, ckpt["flow_class"])
    flow_kwargs = ckpt["flow_kwargs"]

    flow: torch.nn.Module = flow_class(**flow_kwargs)
    flow.load_state_dict(ckpt["flow_state_dict"])
    flow.eval()
    model = FlowSavingWrapper(flow).eval()

    save_as_onnx(model, d=ckpt["flow_kwargs"]["features"], savefile=savefile)


def save_as_onnx(model: torch.nn.Module, d: int, savefile: str):
    model.eval()
    torch.onnx.export(
        model=model,
        args=(torch.zeros(1, d),),
        f=pathlib.Path(savefile).with_suffix(".onnx"),
        input_names=["z"],
        output_names=["log_density"],
        dynamic_shapes={"z": {0: "batch_size"}},
        verify=True,
    )


def train_flow(
    zuko_flow_class,
    flow_kwargs: dict,
    optimizer_class,
    optimizer_kwargs: dict,
    data: torch.Tensor,
    savedir: str,
    filename_prefix: str = "",
    n_epochs: int = 200,
    batch_size: int = 256,
    max_norm: float = 20.0,
    early_stopping: int = 30,
    density_plot_grid_size: int = 100,
) -> tuple[list[float], list[float]]:
    assert len(data.shape) == 2  # (N, n_dim)
    assert len(data) >= batch_size, f"need at least {batch_size} samples, got {len(data)}"
    assert data.shape[-1] == flow_kwargs["features"]

    pathlib.Path(savedir).mkdir(parents=True, exist_ok=True)

    flow: torch.nn.Module = zuko_flow_class(**flow_kwargs)
    optimizer: torch.optim.Optimizer = optimizer_class(params=flow.parameters(), **optimizer_kwargs)

    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer,
        T_max=n_epochs,
        eta_min=1e-8,
    )

    dataset = torch.utils.data.TensorDataset(data)
    train_dataset, val_dataset, test_dataset = torch.utils.data.random_split(dataset, [0.80, 0.10, 0.10])

    train_loader = torch.utils.data.DataLoader(dataset=train_dataset, drop_last=True, batch_size=batch_size, shuffle=True)

    model_save_path = pathlib.Path(savedir).joinpath(f"{filename_prefix}best_checkpoint.pth")
    trainset_path = model_save_path.with_name(f"{filename_prefix}train_data.pth")
    valset_path = model_save_path.with_name(f"{filename_prefix}val_data.pth")
    testset_path = model_save_path.with_name(f"{filename_prefix}test_data.pth")
    log_path = model_save_path.with_name(f"{filename_prefix}log.txt")

    torch.save(dataset[train_dataset.indices][0], trainset_path)
    torch.save(dataset[val_dataset.indices][0], valset_path)
    torch.save(dataset[test_dataset.indices][0], testset_path)

    handler = logging.FileHandler(log_path, "w", "utf-8")
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    logger.info(
        f"flow:\n{flow}\n\noptimizer:\n{optimizer}\n\nscheduler:\n{scheduler}\n\nn_epochs: {n_epochs}\nbatch_size: {batch_size}\nmax_norm: {max_norm}"
    )
    logger.info(f"data.shape: {data.shape}\ntrain size: {len(train_dataset)}\nvalidation size: {len(val_dataset)}\ntest size: {len(test_dataset)}")

    mean_epoch_train_losses = []
    mean_epoch_val_losses = []
    min_epoch_val_loss = torch.inf
    best_epoch = -1

    flow.train()

    start = time.perf_counter()
    for epoch in tqdm.tqdm(range(n_epochs)):
        epoch_train_loss = 0.0

        for batch in train_loader:
            (z,) = batch

            # Negative log likelihood
            loss = -flow().log_prob(z).mean()

            optimizer.zero_grad(set_to_none=True)
            loss.backward()

            preclip_norm = torch.nn.utils.clip_grad_norm_(flow.parameters(), max_norm=max_norm)
            if preclip_norm > max_norm:
                logger.warning(f"[{epoch:6d} / {n_epochs}]\tPre-clip norm exceeded clip threshold: {preclip_norm:.2f} > {max_norm}")

            optimizer.step()

            epoch_train_loss += loss.detach().item()

        # Show top 5 largest gradients
        by_param = Counter({name: p.grad.norm().item() for name, p in flow.named_parameters() if p.grad is not None})
        for name, norm in sorted(by_param.items(), key=lambda x: -x[1])[:5]:
            logger.info(f"[{epoch:6d} / {n_epochs}]\tTop 5 per-parameter gradient norms\t{name}: {norm:.2f}")

        mean_epoch_train_loss = epoch_train_loss / len(train_loader)
        mean_epoch_train_losses.append(mean_epoch_train_loss)

        logger.info(f"[{epoch:6d} / {n_epochs}]\ttrain loss = {mean_epoch_train_loss:3.4f}")

        prev_lr = optimizer.param_groups[0]["lr"]
        scheduler.step()
        current_lr = optimizer.param_groups[0]["lr"]
        if current_lr != prev_lr:
            logger.info(f"[{epoch:6d} / {n_epochs}]\tLR reduced: {prev_lr:.2e} -> {current_lr:.2e}")

        flow.eval()
        with torch.no_grad():
            mean_epoch_val_loss = -flow().log_prob(dataset[val_dataset.indices][0]).mean().item()
        flow.train()

        mean_epoch_val_losses.append(mean_epoch_val_loss)
        logger.info(f"[{epoch:6d} / {n_epochs}]\tval loss = {mean_epoch_val_loss:3.4f}")

        if mean_epoch_val_loss < min_epoch_val_loss:
            min_epoch_val_loss = mean_epoch_val_loss
            best_epoch = epoch

            logger.info(f"[{epoch:6d} / {n_epochs}]\tNew minimum validation loss reached; saving checkpoint")
            torch.save(
                {
                    "epoch": epoch,
                    "validation_loss": mean_epoch_val_loss,
                    "flow_class": zuko_flow_class.__name__,
                    "flow_kwargs": flow_kwargs,
                    "flow_state_dict": flow.state_dict(),
                    "optimizer_class": optimizer_class.__name__,
                    "optimizer_kwargs": optimizer_kwargs,
                    "optimizer_state_dict": optimizer.state_dict(),
                },
                model_save_path,
            )
        else:
            if epoch - best_epoch == early_stopping:
                logger.info(f"[{epoch:6d} / {n_epochs}]\tStopping early - no improvement in validation loss after {early_stopping} epochs.")
                break

    end = time.perf_counter()

    flow.eval()
    with torch.no_grad():
        avg_log_likelihood_on_trainset = flow().log_prob(dataset[train_dataset.indices][0]).mean().item()
        avg_log_likelihood_on_valset = flow().log_prob(dataset[val_dataset.indices][0]).mean().item()
        avg_log_likelihood_on_testset = flow().log_prob(dataset[test_dataset.indices][0]).mean().item()

    logger.info(f"Training finished. Took {(end - start) / 60:.2f} minutes")
    logger.info(f"Mean log likelihood (train set) = {avg_log_likelihood_on_trainset:3.4f}")
    logger.info(f"Mean log likelihood (val set)   = {avg_log_likelihood_on_valset:3.4f}")
    logger.info(f"Mean log likelihood (test set)  = {avg_log_likelihood_on_testset:3.4f}")

    logger.removeHandler(handler)

    plt.clf()
    plt.plot(mean_epoch_train_losses, "r--", label="train loss")
    plt.plot(mean_epoch_val_losses, "g--", label="validation loss")
    plt.legend()
    plt.xlabel("epoch")
    plt.ylabel("negative log likelihood")
    plt.savefig(model_save_path.with_name(f"{filename_prefix}losses.png"))

    def wrap_call(x: np.ndarray) -> np.ndarray:
        with torch.no_grad():
            density = torch.exp(flow().log_prob(torch.from_numpy(x))).numpy()
        return density

    density, grids = evaluate_over_grid(
        f=wrap_call,
        grid_bounds=[(torch.min(data[:, i]).item(), torch.max(data[:, i]).item()) for i in range(data.shape[-1])],
        grid_size=density_plot_grid_size,
        ##
    )

    fig, _ = density_corner_plot(
        density=density,
        grids=grids,
        samples=data.numpy(),
        ##
    )
    fig.savefig(model_save_path.with_name(f"{filename_prefix}density.png"))

    return mean_epoch_train_losses, mean_epoch_val_losses


def evaluate_over_grid(f, grid_bounds: list[tuple[float, float]], grid_size: int = 50) -> tuple[np.ndarray, list[np.ndarray]]:
    ndim = len(grid_bounds)

    if grid_size**ndim > 1e7:
        raise UserWarning(
            f"density_corner_plot:\nTrying to compute density over too many points (grid_size ** ndim = {grid_size**ndim:,} - limit = 10,000,000). Decrease grid_size for higher dimensional distributions."
        )

    grids = [np.linspace(lo, hi, grid_size, dtype=np.float32) for lo, hi in grid_bounds]

    mesh = np.meshgrid(*grids, indexing="ij")
    x = np.stack([m.ravel() for m in mesh], axis=-1)

    density = f(x).reshape((grid_size,) * ndim)
    return density, grids


def density_corner_plot(
    density: np.ndarray,
    grids: list[np.ndarray],
    samples: np.ndarray,
    parameter_labels: list[str] | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    ndim = samples.shape[-1]

    parameter_plot_bounds = [(np.min(g), np.max(g)) for g in grids]

    if parameter_labels is None:
        parameter_labels = [f"$z_{i}$" for i in range(1, ndim + 1)]

    def marginalize(arr, keep_dims):
        # arr shape: (N, ..., N); integrate out all dims not in keep_dims
        drop = sorted([d for d in range(ndim) if d not in keep_dims], reverse=True)
        for d in drop:
            arr = np.trapezoid(arr, grids[d], axis=d)
        return arr  # (N,) or (N, N)

    max_points_to_show = int(min(samples.shape[0], 3000))

    fig, axes = plt.subplots(ndim, ndim, figsize=(5 * ndim, 5 * ndim))

    for row in range(ndim):
        for col in range(ndim):
            ax = axes[row, col]
            if col > row:
                ax.set_visible(False)
                continue

            if col == row:  # 1D marginal
                marginal_1d = marginalize(density, keep_dims=[col])
                ax.hist(
                    samples[:, col],
                    bins=30,
                    density=True,
                    color="#593791",
                    alpha=0.4,
                    histtype="stepfilled",
                )

                ax.plot(grids[col], marginal_1d, color="black")

                ax.set_ylim(bottom=0.0)
                ax.set_xlim(parameter_plot_bounds[col])

            else:  # 2D marginal
                marginal_2d = marginalize(density, keep_dims=[col, row])

                qcs = ax.contourf(grids[col], grids[row], marginal_2d.T, levels=25, cmap="inferno")
                fig.colorbar(qcs, ax=ax)

                ax.scatter(
                    samples[:max_points_to_show, col],
                    samples[:max_points_to_show, row],
                    s=1.0,
                    color="gray",
                    alpha=0.5,
                    linewidths=0,
                )

                ax.set_xlim(parameter_plot_bounds[col])
                ax.set_ylim(parameter_plot_bounds[row])

            if row == ndim - 1:
                ax.set_xlabel(parameter_labels[col])
            if col == 0 and row != 0:
                ax.set_ylabel(parameter_labels[row])

    fig.tight_layout()
    return fig, axes
