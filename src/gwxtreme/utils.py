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

import h5py
import lal
import numpy as np


def read_prior_or_posterior_file(file: str, cbc_dim: int) -> dict:
    """Return samples for the intrinsic BNS/NSBH parameters
    read from the given samples file.

    Parameters
    ----------
    file
        Must be in one of these formats: .h5/.hdf5, .json, .dat

        Required contents and keys depend on given ``cbc_dim``.
        If ``cbc_dim`` is 2, the file must include:

            q, mc_source, lambdat

            Possible alternative names:
            mass_ratio, chirp_mass_source, lambda_tilde

        If ``cbc_dim`` is 3 or 4, the file must include:

            m1_source, m2_source, q, mc_source, lambda_1, lambda_2

            Possible alternative names:
            mass_1_source, mass_2_source, mass_ratio, chirp_mass_source, lambda1, lambda2

    cbc_dim
        Must be 2, 3, or 4. This corresponds to the event type used in the
        main inference codes of gw-2d, gw-3d, or gw-4d.

    Returns
    -------
        Dictionary with keys ["m1_source", "m2_source", "q", "mc_source", "lambdat", "lambda1", "lambda2"]
    """

    file_ = pathlib.Path(file)
    ext = file_.suffix

    m1, m2, q, mc, lambda1, lambda2, lambdat = None, None, None, None, None, None, None

    if ext == ".h5" or ext == ".hdf5":
        with h5py.File(file_) as f:
            data = np.array(f["posterior_samples"])

    elif ext == ".json":
        with open(file_) as f:
            data = json.load(f)["posterior"]["content"]

    elif ext == ".dat":
        data = np.genfromtxt(file_, names=True)

    if cbc_dim == 2:
        try:
            q = np.array(data["q"])
        except KeyError:
            q = np.array(data["mass_ratio"])

        try:
            mc = np.array(data["mc_source"])
        except KeyError:
            mc = np.array(data["chirp_mass_source"])

        try:
            lambdat = np.array(data["lambdat"])
        except KeyError:
            lambdat = np.array(data["lambda_tilde"])

    else:
        try:
            m1 = np.array(data["m1_source"])
            m2 = np.array(data["m2_source"])
        except KeyError:
            m1 = np.array(data["mass_1_source"])
            m2 = np.array(data["mass_2_source"])

        try:
            q = np.array(data["q"])
        except KeyError:
            q = np.array(data["mass_ratio"])

        try:
            mc = np.array(data["mc_source"])
        except KeyError:
            mc = np.array(data["chirp_mass_source"])

        try:
            lambda1 = np.array(data["lambda_1"])
            lambda2 = np.array(data["lambda_2"])
        except KeyError:
            lambda1 = np.array(data["lambda1"])
            lambda2 = np.array(data["lambda2"])

    return {
        "m1_source": m1,
        "m2_source": m2,
        "q": q,
        "mc_source": mc,
        "lambdat": lambdat,
        "lambda1": lambda1,
        "lambda2": lambda2,
    }


def get_gw_event_pe_posterior_samples(posterior_file: str, cbc_dim: int) -> np.ndarray:
    """Convenience utility to read a posterior file and return
    a stacked array containing the necessary parameter samples
    for GWXtreme inference based on the variant of the approximation
    specified by ``cbc_dim``.

    Parameters
    ----------
    posterior_file
        Must be in one of these formats: .h5/.hdf5, .json, .dat

        Required contents and keys depend on given ``cbc_dim``.
        If ``cbc_dim`` is 2, the file must include:

            q, mc_source, lambdat

            Possible alternative names:
            mass_ratio, chirp_mass_source, lambda_tilde

        If ``cbc_dim`` is 3 or 4, the file must include:

            m1_source, m2_source, q, mc_source, lambda_1, lambda_2

            Possible alternative names:
            mass_1_source, mass_2_source, mass_ratio, chirp_mass_source, lambda1, lambda2

    cbc_dim
        Must be 2, 3, or 4. This corresponds to the event type used in the
        main inference codes of gw-2d, gw-3d, or gw-4d.

    Returns
    -------
        Array of shape (n_samples, cbc_dim). The order of parameters is (lambdat, q) for 2D, (lambda1, q, lambda2) for 3D, and (m1, m2, lambda1, lambda2) for 4D.
    """

    samples = read_prior_or_posterior_file(posterior_file, cbc_dim)

    if cbc_dim == 2:
        return np.stack(
            (
                samples["lambdat"],
                samples["q"],
            ),
            axis=-1,
        )
    elif cbc_dim == 3:
        return np.stack(
            (
                samples["lambda1"],
                samples["q"],
                samples["lambda2"],
            ),
            axis=-1,
        )
    elif cbc_dim == 4:
        return np.stack(
            (
                samples["m1_source"],
                samples["m2_source"],
                samples["lambda1"],
                samples["lambda2"],
            ),
            axis=-1,
        )
    else:
        raise NotImplementedError()


def get_nicer_pulsar_pe_posterior_samples(posterior_file: str) -> np.ndarray:
    """Read mass, compactness samples and convert to mass, radius (in km)

    Parameters
    ----------
    posterior_file
        .txt file containing 2 columns: mass in solar masses
        and compactness

    Returns
    -------
        Array of shape (n_samples, 2). The first column contains the mass values as given in solar masses. The second column contains radii in km.
    """

    if posterior_file[-4:] == ".txt":
        samples = np.loadtxt(posterior_file)

        # Convert compactness to radius in km
        mass_in_sm, compactness = samples[:, 0], samples[:, 1]
        mass_in_kg = lal.MSUN_SI * mass_in_sm
        radius_in_km = (lal.G_SI * mass_in_kg / (lal.C_SI**2 * compactness)) / 1000
        return np.stack((mass_in_sm, radius_in_km), axis=-1)
    else:
        raise NotImplementedError()


def get_mean_mchirp_for_cbc_event(posterior_file: str, cbc_dim: int) -> float:
    """Return the mean chirp mass from the given samples file"""

    samples = read_prior_or_posterior_file(posterior_file, cbc_dim=cbc_dim)
    return np.mean(samples["mc_source"]).item()


def get_q_range(samples_file: str, cbc_dim: int) -> tuple[float, float]:
    """Return the min and max values of q from the given samples file"""

    samples = read_prior_or_posterior_file(samples_file, cbc_dim=cbc_dim)
    q = samples["q"]
    return np.min(q).item(), np.max(q).item()
