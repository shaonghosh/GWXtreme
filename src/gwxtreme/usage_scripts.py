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

import json
import pathlib
from typing import Literal

import lalsimulation
import matplotlib.pyplot as plt
import numpy as np

from gwxtreme.eos_inference import (
    JointModelSelector,
    ModelSelector,
)
from gwxtreme.eos_interpolator import get_eos_lambdas_from_masses, get_log_pressure, get_parameterized_eos


def compute_bayes_factors_for_several_named_eos(
    model_selector: ModelSelector | JointModelSelector,
    named_eos_list: list[str],
    reference_eos_name: str,
    n_grid: int,
    n_resamplings: int,
    n_jobs: int,
    save_file: str | None = None,
) -> dict:
    results = {}

    if isinstance(model_selector, ModelSelector):
        reference_eos_evidences = model_selector.compute_eos_evidence(
            eos_name=reference_eos_name, n_grid=n_grid, n_resamplings=n_resamplings, n_jobs=n_jobs
        )

        for eos in named_eos_list:
            evidences = model_selector.compute_eos_evidence(eos_name=eos, n_grid=n_grid, n_resamplings=n_resamplings, n_jobs=n_jobs)

            bfs = evidences / reference_eos_evidences

            results[eos] = bfs.tolist()

    else:
        reference_eos_joint_evidences, ref_per_event_evidences = model_selector.compute_joint_eos_evidence(
            eos_name=reference_eos_name, n_grid=n_grid, n_resamplings=n_resamplings, n_jobs=n_jobs
        )

        per_event_results = [{} for _ in range(len(ref_per_event_evidences))]

        for eos in named_eos_list:
            joint_evidences, per_event_evidences = model_selector.compute_joint_eos_evidence(
                eos_name=eos, n_grid=n_grid, n_resamplings=n_resamplings, n_jobs=n_jobs
            )

            per_event_bfs = [
                event_evidence / ref_event_evidence for event_evidence, ref_event_evidence in zip(per_event_evidences, ref_per_event_evidences)
            ]

            for i, event_bfs in enumerate(per_event_bfs):
                per_event_results[i][eos] = event_bfs

            results[eos] = np.prod(np.array(per_event_bfs), axis=0).tolist()

    if save_file:
        with open(save_file, "w") as f:
            json.dump(results, f, indent=4)

        # Additionally save per-event Bayes factors if using multi-event selector
        if isinstance(model_selector, JointModelSelector):
            save_path = pathlib.Path(save_file)
            for i, event_bfs in enumerate(per_event_results):
                results = {eos: event_bfs[eos].tolist() for eos in named_eos_list}

                with open(save_path.with_stem(f"{save_path.stem}__event_{i}"), "w") as f:
                    json.dump(results, f, indent=4)

    return results


def compute_eos_constraints_from_parameter_samples(
    samples: np.ndarray,
    parameterization: Literal["spectral", "polytrope"],
    save_file: str | None = None,
):
    # Turn into confidence interval data
    log_pressure = []
    density = np.logspace(17.1, 18.25, 1000)

    for params in samples:
        eos = get_parameterized_eos(params, parameterization)
        log_pressure.append(get_log_pressure(eos, density))

    log_pressure = np.array(log_pressure)

    # Compute median and 90% credible interval of log pressure values for each density
    logp_CIup = np.quantile(log_pressure, 0.95, axis=0)
    logp_CIlow = np.quantile(log_pressure, 0.05, axis=0)
    logp_med = np.quantile(log_pressure, 0.50, axis=0)
    logp_min = np.min(log_pressure, axis=0)
    logp_max = np.max(log_pressure, axis=0)

    out = np.array([density, logp_CIlow, logp_med, logp_CIup, logp_min, logp_max]).T

    if save_file is not None:
        # Save credible interval data
        np.savetxt(save_file, out)

    return out


def compute_eos_lambdas_from_parameter_samples(
    samples: np.ndarray,
    parameterization: Literal["spectral", "polytrope"],
    reference_ns_mass: float = 1.4,
    save_file: str | None = None,
) -> np.ndarray:
    lambdas = []

    for params in samples:
        eos = get_parameterized_eos(params, parameterization)
        fam = lalsimulation.CreateSimNeutronStarFamily(eos)

        _, lambda_ = get_eos_lambdas_from_masses(eos_fam=fam, masses=np.array([reference_ns_mass]))

        lambdas.append(lambda_)

    lambdas = np.array(lambdas)

    if save_file is not None:
        np.save(pathlib.Path(save_file).with_suffix(".npy"), lambdas)

    return lambdas


def plot_bayes_factors_bar_chart(
    bayes_factor_bar_heights: dict[str, dict[str, np.ndarray | list]],
    bayes_factor_error_bars: dict[str, dict[str, np.ndarray | list]],
    save_file: str | None = None,
    yscale: Literal["linear", "log"] = "linear",
) -> tuple[plt.Figure, plt.Axes]:
    num_methods = len(bayes_factor_bar_heights.keys())
    if num_methods > 4:
        raise UserWarning("More than 4 BF sets / methods not supported for plotting.")

    colors = ["#f94b42", "#c8c0ff", "#ffa551", "#0d741b"]
    spacing_options = {
        1: [0.0],
        2: [-0.10, 0.10],
        3: [-0.20, 0.0, 0.20],
        4: [-0.30, -0.10, 0.10, 0.30],
    }
    spacing = spacing_options[num_methods]

    fig, axes = plt.subplots(figsize=(16, 10))

    x_axis = np.arange(len(list(bayes_factor_bar_heights.values())[0].keys()))

    for i, ((run_label, bfs), bf_errors) in enumerate(zip(bayes_factor_bar_heights.items(), bayes_factor_error_bars.values())):
        axes.bar(
            x=x_axis + spacing[i],
            height=list(bfs.values()),
            width=0.20,
            label=run_label,
            color=colors[i],
        )
        axes.errorbar(
            x=x_axis + spacing[i],
            y=list(bfs.values()),
            yerr=list(bf_errors.values()),
            ls="none",
            ecolor="black",
        )

        # needed for x-ticks
        named_eos_list = list(bfs.keys())

    if yscale == "log":
        axes.set_yscale("log")
        ylim_bottom, ylim_top = axes.get_ylim()
        axes.set_ylim(bottom=max(ylim_bottom, 1e-20), top=min(ylim_top, 1e6))
    else:
        axes.set_ylim(bottom=0.0)

    axes.set_xticks(x_axis, named_eos_list, rotation=90, ha="right", fontsize=14)
    axes.axhline(1.0, color="k", linestyle="--", alpha=0.2)
    axes.set_ylabel("Bayes Factor w.r.t. SLY", fontsize=18)
    axes.legend(loc="upper right", fontsize=16)

    if save_file:
        fig.savefig(save_file, bbox_inches="tight")

    return fig, axes


def plot_parameterized_eos_constraints(
    constraints_files: list[str],
    labels: list[str],
    save_file: str,
    named_eos_list: list[str] | None = None,
    prior_constraints_file: str | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    colors = ["#64ACDC", "#c06161", "#6b9e64", "#82ca84"]
    hatches = ["", "|", "\\", "/"]
    named_eos_styles = ["k--", "g--", "r--", "b--"]

    fig, axes = plt.subplots(figsize=(8, 8))
    plt.rc("xtick", direction="out", color="black")
    plt.rc("ytick", direction="out", color="black")
    plt.rc("lines", linewidth=2)

    for file, label, color, hatch in zip(constraints_files, labels, colors, hatches):
        constraints = np.loadtxt(file)
        density = constraints[:, 0]
        logp_CIlow = constraints[:, 1]
        logp_CIup = constraints[:, 3]

        log10_density = np.log10(density)

        axes.fill_between(
            log10_density,
            logp_CIlow,
            logp_CIup,
            color=color,
            alpha=0.45,
            label=label,
            zorder=1.0,
            hatch=hatch,
        )

    if named_eos_list is not None:
        for i, eos in enumerate(named_eos_list):
            named_eos_logp = get_log_pressure(lalsimulation.SimNeutronStarEOSByName(eos), density)

            axes.plot(
                log10_density,
                named_eos_logp,
                named_eos_styles[i],
                linewidth=2.0,
                label=eos,
                alpha=0.45,
            )

    if prior_constraints_file is not None:
        constraints = np.loadtxt(prior_constraints_file)
        density = constraints[:, 0]
        logp_CIlow = constraints[:, 1]
        logp_CIup = constraints[:, 3]
        logp_min = constraints[:, 4]
        logp_max = constraints[:, 5]

        log10_density = np.log10(density)

        axes.plot(log10_density, logp_CIlow, color="red", linestyle="dotted", lw=1.3)
        axes.plot(log10_density, logp_CIup, color="red", linestyle="dotted", label="prior 90%", lw=1.3)

        axes.plot(log10_density, logp_min, color="black", linestyle="dotted", lw=1.3)
        axes.plot(log10_density, logp_max, color="black", linestyle="dotted", label="prior extrema", lw=1.3)

    axes.set_xlim(left=np.min(log10_density), right=np.max(log10_density))
    axes.set_xlabel(r"$\log_{10}(\frac{\rho}{\mathrm{g \, cm^{-3}}})$", fontsize=17)
    axes.set_ylabel(r"$\log_{10}(\frac{p}{\mathrm{dyn \, cm^{-2}}})$", fontsize=17)
    axes.legend(fontsize=12)
    fig.savefig(save_file, bbox_inches="tight", dpi=300)

    return fig, axes


def plot_lambda_dist(
    lambdas_files: list[str],
    labels: list[str],
    save_file: str,
    reference_ns_mass: float | None = 1.4,
    true_eos: str | None = None,
):
    colors = ["#64ACDC", "#c06161", "#6b9e64", "#82ca84"]

    fig, axes = plt.subplots(figsize=(8, 8))
    plt.rc("xtick", direction="out", color="black")
    plt.rc("ytick", direction="out", color="black")
    plt.rc("lines", linewidth=2)

    for file, color, method_label in zip(lambdas_files, colors, labels):
        Lambdas = np.load(file)
        axes.hist(
            Lambdas,
            label=method_label,
            alpha=1.0,
            fill=False,
            bins=25,
            density=True,
            color=color,
            histtype="step",
        )

    if true_eos:
        eos = lalsimulation.SimNeutronStarEOSByName(true_eos)
        fam = lalsimulation.CreateSimNeutronStarFamily(eos)

        _, true_lambda = get_eos_lambdas_from_masses(fam, np.array([reference_ns_mass]))

        axes.axvline(x=true_lambda[0].item(), label=true_eos, color="black")

    axes.set_xlabel(r"$\Lambda(1.4$ $\mathrm{M}_{\odot})$", fontsize=17)
    axes.set_yticks([])
    axes.legend(fontsize=12)
    fig.savefig(save_file, bbox_inches="tight", dpi=300)

    return fig, axes
