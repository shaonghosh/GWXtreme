# Copyright (C) 2026 Joseph David Quinn-Vitabile
# Copyright (C) 2017 Daniel Wysocki
# Assimilated from Wysocki D., O'Shaughnessy R., Bayesian Parametric Population Models,
# 2017–, bayesian-parametric-population-models.readthedocs.io
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

from typing import Literal

import lal
import lalsimulation
import numpy as np

from gwxtreme.eos_interpolator import get_parameterized_eos


def spectral_eos_adiabatic_index(spectral_params: list | tuple | np.ndarray, log10_pressure: np.ndarray) -> np.ndarray:
    """Compute the adiabatic index at a particular value of log pressure,
    using the spectral EOS, given a particular choice of spectral
    parameters.

    Parameters
    ----------
    spectral_params
        Vector of 4 parameters (gammas) characterizing the spectral EOS
    log10_pressure
        Values at which the adiabatic index will be computed

    Returns
    -------
        Array of adiabatic index values
    """

    log_gamma = (
        spectral_params[0]
        + spectral_params[1] * log10_pressure
        + spectral_params[2] * log10_pressure * log10_pressure
        + spectral_params[3] * log10_pressure * log10_pressure * log10_pressure
    )
    return np.exp(log_gamma)


def is_valid_adiabatic_index(spectral_params: list | tuple | np.ndarray) -> bool:
    """Confirm that the adiabatic index is within a tolerable range
    (0.6, 4.5) for the given spectral EOS.

    Parameters
    ----------
    spectral_params
        Vector of 4 parameters (gammas) characterizing the spectral EOS

    Returns
    -------
        Bool indicating if the adiabatic index is valid
    """

    log10_pressure_grid = np.linspace(0.0, 12.3081, 500)
    adiabatic_index = spectral_eos_adiabatic_index(spectral_params, log10_pressure_grid)
    return bool(np.all(adiabatic_index > 0.6)) and bool(np.all(adiabatic_index < 4.5))


def has_enough_points(eos_fam) -> bool:
    """Confirm that the TOV solver returns increasing
    masses for increasing pressure for at least 8 points
    in the pressure grid. If this is not true, ``lalsimulation``
    may return interpolation errors.

    Parameters
    ----------
    eos_fam
        ``lalsimulation`` EOS family object

    Returns
    -------
        Bool indicating if the EOS is monotonically increasing for 8 points
    """

    min_points = 8

    logpmin = 75.5
    logpmax = np.log(lalsimulation.SimNeutronStarEOSMaxPressure(eos_fam))

    dlogp = (logpmax - logpmin) / 100

    m_prev = 0.0
    for i in range(min_points):
        p = np.exp(logpmin + i * dlogp)
        r, m, k = lalsimulation.SimNeutronStarTOVODEIntegrate(p, eos_fam)

        if m <= m_prev:
            return False
        m_prev = m

    return True


def eos_max_sound_speed(eos, eos_fam, max_mass_kg: float | None = None) -> float:
    """Compute the maximum sound speed (at the center of the NS)
    predicted by the given EOS.

    Parameters
    ----------
    eos
        ``lalsimulation`` EOS object
    eos_fam
        ``lalsimulation`` EOS family object
    max_mass_kg
        Optionally provide the max mass in kg for the given
        EOS to reduce computation if this is known already

    Returns
    -------
        Max sound speed
    """

    # Maximum allowed mass
    # allow passing max_mass_kg to avoid a bit of redundant computation when calling
    # this from is_valid_eos()
    if not max_mass_kg:
        max_mass_kg = lalsimulation.SimNeutronStarMaximumMass(eos_fam)

    # Central pressure
    p_max = lalsimulation.SimNeutronStarCentralPressure(max_mass_kg, eos_fam)
    # Pseudo-enthalpy at the core.
    h_max = lalsimulation.SimNeutronStarEOSPseudoEnthalpyOfPressure(p_max, eos)
    # Maximum sound speed should occur at the max pseudo-enthalpy for a typical
    # TOV neutron star.
    max_sound_speed = lalsimulation.SimNeutronStarEOSSpeedOfSoundGeometerized(h_max, eos)

    return max_sound_speed


def is_valid_eos(eos_params: np.ndarray | tuple | list, parameterization: Literal["spectral", "polytrope"], largest_ns_mass: float = 1.97) -> bool:
    """Check if the given parameterized EOS is physical and stable, i.e.
    it must be thermally stable (monotonically increasing in density(pressure)),
    causal (maximum speed of sound < 1.10 c), and result in a maximum NS mass
    that is larger than the mass of the most massive observed NS.


    Parameters
    ----------
    eos_params
        Vector of 4 parameters characterizing the EOS
    parameterization
        "spectral" for 4-parameter spectral model or "polytrope"
        for 4-parameter piecewise polytrope model
    largest_ns_mass
        Mass of the heaviest known NS in solar masses, by default 1.97

    Returns
    -------
        Bool indicating if the given EOS satisfies all the mentioned conditions
    """

    if parameterization == "spectral":
        if not is_valid_adiabatic_index(eos_params):
            return False

    eos = get_parameterized_eos(eos_params, parameterization)

    # to avoid interpolation errors from LAL
    if not has_enough_points(eos):
        return False

    eos_fam = lalsimulation.CreateSimNeutronStarFamily(eos)
    eos_m_max = lalsimulation.SimNeutronStarMaximumMass(eos_fam)

    if eos_m_max < largest_ns_mass * lal.MSUN_SI:
        return False

    # final condition: causality
    return eos_max_sound_speed(eos, eos_fam, max_mass_kg=eos_m_max) < 1.10
