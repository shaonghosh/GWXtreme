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

"""This module contains functions for interfacing with ``lalsimulation`` to construct
and compute with neutron star equation of state models.

The ``EOSInterpolator`` class serves as a convenient utility for creating named or
parameterized EOS models and using them to interpolate tidal deformability or radius
values from given masses. Mass constraints can also be applied from this class based
on the minimum NS mass supported by the created EOS.
"""

from typing import Literal

import lal
import lalsimulation
import numpy as np
import scipy.interpolate


def get_lal_named_eos(eos_name: str):
    """Retrieve the ``lalsimulation`` neutron star EOS corresponding to
    ``eos_name``.

    Parameters
    ----------
    eos_name
        Name of the equation of state, e.g. "APR4_EPP", "H4", etc.

    Returns
    -------
        ``lalsimulation.SimNeutronStarEOSByName(eos_name)``
    """

    return lalsimulation.SimNeutronStarEOSByName(eos_name)


def get_parameterized_eos(
    params: np.ndarray | tuple | list, parameterization: Literal["spectral", "polytrope"]
) -> tuple[scipy.interpolate.interp1d, float, float]:
    """Retrieve the ``lalsimulation`` parameterized EOS model with the given
    parameterization and parameter values.

    Parameters
    ----------
    params
        Parameter vector characterizing the EOS
    parameterization
        One of "spectral" (4-parameter spectral decomposition) or "polytrope"
        (4-parameter piecewise polytrope)

    Returns
    -------
        The parameterized EOS model
    """

    if parameterization == "polytrope":
        # params are log_p1_SI, g1, g2, g3
        eos = lalsimulation.SimNeutronStarEOS4ParameterPiecewisePolytrope(params[0], params[1], params[2], params[3])
    else:
        # params are g0, g1, g2, g3
        eos = lalsimulation.SimNeutronStarEOS4ParameterSpectralDecomposition(params[0], params[1], params[2], params[3])

    return eos


def read_eos_mass_lambda_file(mass_lambda_file: str) -> tuple[np.ndarray, np.ndarray]:
    """Read neutron star mass, tidal deformability values from the given .txt file, representing the mass-radius relation prescribed by some equation of state.

    The data should be in the following format (without titles)::

        (mass column)       (lambda column)
        min_mass            ...
        ...                 ...
        ...                 ...
        max_mass            ...

    The values of masses should be in units of solar masses. The
    tidal deformability λ should be supplied in SI units.

    Parameters
    ----------
    mass_lambda_file
        .txt file containing mass and tidal deformability (λ) values

    Returns
    -------
        masses (as given), dimensionless tidal deformability
    """

    masses, lambdas = np.loadtxt(mass_lambda_file, unpack=True)
    lambdas = lal.G_SI * lambdas * (1 / (lal.MRSUN_SI * masses) ** 5)

    return masses, lambdas


def read_eos_mass_radius_k_file(mass_radius_k_file: str):
    """Read neutron star mass, radius, and tidal Love number values from the given .txt file, representing the mass-radius-kappa relation prescribed by some equation of state.

    The data should be in the following format (without titles)::

        (mass column)       (radius column)     (kappa column)
        min_mass            ...                 ...
        ...                 ...                 ...
        ...                 ...                 ...
        max_mass            ...                 ...

    The values of masses should be in units of solar masses. The radius should
    be supplied in meters.

    Parameters
    ----------
    mass_radius_k_file
        .txt file containing mass, radii, and tidal Love number values

    Returns
    -------
        masses (as given), dimensionless tidal deformabilities, radii (in km)
    """

    masses, radii, kappas = np.loadtxt(mass_radius_k_file, unpack=True)
    compactness = masses * lal.MRSUN_SI / radii
    lambdas = (2 / 3) * kappas / (compactness**5)

    return masses, lambdas, radii * 1e-3  # convert radii to km for the interpolator


def get_eos_lambdas_from_masses(eos_fam, masses: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute the tidal deformabilities corresponding to the given
    ``masses`` for the ``lalsimulation`` EOS family given by ``eos_fam``.

    Note that some of the provided mass values may induce runtime
    errors in ``lalsimulation`` functions, in which case these masses
    will be skipped. The output mass and tidal deformability arrays
    may therefore be smaller than the size of the input ``masses``.

    Parameters
    ----------
    eos_fam
        ``lalsimulation`` EOS family to use
    masses
        Array of neutron star masses in solar masses

    Returns
    -------
        Tuple of 2 arrays: the NS masses and the associated tidal deformabilities
    """

    lambdas = []
    valid_masses = []
    for m in masses:
        try:
            r = lalsimulation.SimNeutronStarRadius(m * lal.MSUN_SI, eos_fam)
            k = lalsimulation.SimNeutronStarLoveNumberK2(m * lal.MSUN_SI, eos_fam)
            c = m * lal.MRSUN_SI / r
            lambdas.append((2 / 3) * k / (c**5))
            valid_masses.append(m)
        except RuntimeError:
            continue

    return np.array(valid_masses), np.array(lambdas)


def get_eos_radii_from_masses(eos_fam, masses: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute the radii corresponding to the given
    ``masses`` for the ``lalsimulation`` EOS family given by ``eos_fam``.

    Note that some of the provided mass values may induce runtime
    errors in ``lalsimulation`` functions, in which case these masses
    will be skipped. The output mass and radius arrays may therefore
    be smaller than the size of the input ``masses``.

    Parameters
    ----------
    eos_fam
        ``lalsimulation`` EOS family to use
    masses
        Array of neutron star masses in solar masses

    Returns
    -------
        Tuple of 2 arrays: the NS masses and the associated radii
    """

    radii = []
    valid_masses = []
    for m in masses:
        try:
            r = lalsimulation.SimNeutronStarRadius(m * lal.MSUN_SI, eos_fam)
            radii.append(r)
            valid_masses.append(m)
        except RuntimeError:
            continue

    return np.array(valid_masses), np.array(radii) / 1000  # convert to km


def convert_masses(q: float | np.ndarray, mc: float | np.ndarray) -> tuple[float | np.ndarray, float | np.ndarray]:
    """Convert compact binary mass ratio and chirp mass parameters into
    individual component masses.

    Parameters
    ----------
    q
        Mass ratio
    mc
        Chirp mass

    Returns
    -------
        Tuple of 2 arrays (m1, m2) containing the component object masses
    """

    m1 = mc * (1 + q) ** (1 / 5) * (q) ** (-3 / 5)
    m2 = mc * (1 + q) ** (1 / 5) * (q) ** (2 / 5)
    return m1, m2


def convert_lambdas(m1: np.ndarray, m2: np.ndarray, Lambda1: np.ndarray, Lambda2: np.ndarray) -> np.ndarray:
    """Convert compact binary masses and tidal deformabilities
    to the dominant tidal deformability parameter (LambdaTilde).

    Parameters
    ----------
    m1
        Mass 1
    m2
        Mass 2
    Lambda1
        Tidal deformability 1
    Lambda2
        Tidal deformability 2

    Returns
    -------
        Dominant binary tidal deformability parameter lambda tilde
    """

    LambdaTilde = (2**5) / (26 * (m1 + m2) ** 5)
    LambdaTilde *= (m1**5 + 12 * m2 * m1**4) * Lambda1 + (m2**5 + 12 * m1 * m2**4) * Lambda2
    return LambdaTilde


def get_log_pressure(eos, density: np.ndarray) -> np.ndarray:
    """Compute log10 pressure for a range of densities for any EOS.

    Parameters
    ----------
    eos
        ``lalsimulation`` EOS model, named or parameterized
    density
        Array of densities in g cm^-3

    Returns
    -------
        Array of log10 pressure values corresponding to the given densities
    """

    min_pressure = 5.0e31
    max_pressure = min(6.0e36, lalsimulation.SimNeutronStarEOSMaxPressure(eos))
    log10_pressure_grid = np.linspace(np.log10(min_pressure), np.log10(max_pressure), 256)

    pressure_grid = np.power(10.0, log10_pressure_grid)
    density_grid = np.empty_like(pressure_grid)

    for i, pressure in np.ndenumerate(pressure_grid):
        h = lalsimulation.SimNeutronStarEOSPseudoEnthalpyOfPressure(pressure, eos)
        density_grid[i] = lalsimulation.SimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(h, eos)

    log10_density_grid = np.log10(density_grid)

    log10_p_out = np.interp(
        np.log10(density),
        log10_density_grid,
        log10_pressure_grid,
        left=-np.inf,
        right=-np.inf,
    )
    return log10_p_out


class EOSInterpolator:
    """Combination of neutron star equation of state interpolators for mass-radius and mass-tidal deformability."""

    def __init__(
        self,
        eos_name: str | None = None,
        eos_parameters: np.ndarray | None = None,
        parameterization: Literal["spectral", "polytrope"] | None = None,
        parameter_bounds: list[tuple] | None = None,
        mass_lambda_file: str | None = None,
        mass_radius_k_file: str | None = None,
    ):
        """Construct the combined interpolator from either a LALSuite named EOS, a mass-lambda file,
        a mass-radius-tidal Love number file, or a parameterized model from the 4-parameter spectral
        decomposition or 4-parameter piecewise polytrope parameterizations.

        Provide either:
            ``eos_name``

            **OR**

            ``mass_lambda_file``

            **OR**

            ``mass_radius_k_file``

            **OR**

            ``eos_parameters``, ``parameterization``, and ``parameter_bounds``

        Parameters
        ----------
        eos_name
            Name of an EOS as implemented in ``lalsimulation``

        eos_parameters
            Parameter vector characterizing the EOS

        parameterization
            Must be one of "spectral" or "polytrope"; required if passing ``eos_parameters``

        parameter_bounds
            Bounds for a uniform prior distribution of EOS parameters, structured like
            [(g1_min, g1_max), (g2_min, g2_max), (g3_min, g3_max), (g4_min, g4_max)];
            required if passing ``eos_parameters``

        mass_lambda_file
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

        mass_radius_k_file
            .txt file containing mass, radii, and tidal Love number values

            The data should be in the following format (without titles)::

                (mass column)       (radius column)     (kappa column)
                min_mass            ...                 ...
                ...                 ...                 ...
                ...                 ...                 ...
                max_mass            ...                 ...

            The values of masses should be in units of solar masses. The radius should
            be supplied in meters.


        Raises
        ------
        UserWarning
            If necessary arguments are missing
        """

        if mass_lambda_file is not None:
            masses, lambdas = read_eos_mass_lambda_file(mass_lambda_file)

            self.mass_lambda_interpolant = scipy.interpolate.interp1d(masses, lambdas)

            self.min_mass = np.min(masses)
            self.max_mass = np.max(masses)

        elif mass_radius_k_file is not None:
            masses, lambdas, radii = read_eos_mass_radius_k_file(mass_radius_k_file)

            self.mass_lambda_interpolant = scipy.interpolate.interp1d(masses, lambdas)
            self.mass_radius_interpolant = scipy.interpolate.interp1d(masses, radii)
            self.min_mass = np.min(masses)
            self.max_mass = np.max(masses)

        else:
            if eos_name is not None:
                assert eos_name in list(lalsimulation.SimNeutronStarEOSNames), (
                    f'Given EOS name "{eos_name}" is not available in lalsimulation\nAllowed EOS are :\n' + str(lalsimulation.SimNeutronStarEOSNames)
                )
                self.eos = get_lal_named_eos(eos_name)

            elif eos_parameters is not None and parameterization is not None and parameter_bounds is not None:
                self.eos_parameters = eos_parameters
                self.parameterization = parameterization
                self.parameter_bounds = parameter_bounds
                self.eos = get_parameterized_eos(eos_parameters, parameterization)

            else:
                raise UserWarning("Must pass either eos_name or eos_parameters, parameterization, and eos_parameter_bounds.")

            self.eos_fam = lalsimulation.CreateSimNeutronStarFamily(self.eos)

            self.min_mass = lalsimulation.SimNeutronStarFamMinimumMass(self.eos_fam) / lal.MSUN_SI
            self.max_mass = lalsimulation.SimNeutronStarMaximumMass(self.eos_fam) / lal.MSUN_SI

            self.min_mass = int(self.min_mass * 1000 + 1) / 1000
            self.max_mass = int(self.max_mass * 1000) / 1000

            n_points = 200
            masses = np.linspace(self.min_mass, self.max_mass, n_points)
            masses = masses[masses <= self.max_mass]

            grav_masses, radii = get_eos_radii_from_masses(self.eos_fam, masses)
            self.mass_radius_interpolant = scipy.interpolate.interp1d(grav_masses, radii)

            grav_masses, lambdas = get_eos_lambdas_from_masses(self.eos_fam, masses)
            self.mass_lambda_interpolant = scipy.interpolate.interp1d(grav_masses, lambdas)

    def get_lambda_tilde(self, m1: np.ndarray, m2: np.ndarray) -> np.ndarray:
        """Compute the compact binary dominant tidal deformability parameter
        from component object masses and the EOS mass-tidal deformability interpolator.

        Parameters
        ----------
        m1
            Mass 1
        m2
            Mass 2

        Returns
        -------
            Array of lambda tilde values corresponding to given masses
        """

        Lambda1 = self.get_lambda(m1)
        Lambda2 = self.get_lambda(m2)

        # compute binary tidal deformability
        LambdaT = convert_lambdas(m1, m2, Lambda1, Lambda2)

        return LambdaT

    def get_lambda(self, mass: np.ndarray) -> np.ndarray:
        """Compute neutron star tidal deformability from mass using the EOS
        mass-tidal deformability interpolator.

        Parameters
        ----------
        mass
            Array of neutron star masses in solar masses

        Returns
        -------
            Array of neutron star tidal deformabilities
        """

        kerr_cases = mass >= self.max_mass
        Lambda = np.zeros_like(mass)

        Lambda[~kerr_cases] = self.mass_lambda_interpolant(mass[~kerr_cases])
        return Lambda

    def get_radius(self, mass: np.ndarray) -> np.ndarray:
        """Compute neutron star radius from mass using the EOS mass-radius interpolator.

        Parameters
        ----------
        mass
            Array of neutron star masses in solar masses

        Returns
        -------
            Array of neutron star radii in kilometers
        """

        radius = self.mass_radius_interpolant(mass)
        return radius

    def apply_bns_mass_constraint(
        self,
        m1: np.ndarray,
        m2: np.ndarray,
        q: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the sub-arrays of the given m1, m2, and q values
        for which neither m1 nor m2 are less than this EOS's minimum
        mass.

        Parameters
        ----------
        m1
            Array of mass values in solar masses
        m2
            Array of mass values in solar masses
        q
            Array of mass ratio values

        Returns
        -------
            Tuple (m1, m2, q) containing the respective sub-arrays of given values after applying the EOS minimum mass constraint
        """

        min_mass_violation_1 = m1 < self.min_mass
        min_mass_violation_2 = m2 < self.min_mass

        mass_violation = min_mass_violation_1 + min_mass_violation_2

        m1 = m1[~mass_violation]
        m2 = m2[~mass_violation]
        q = q[~mass_violation]
        return m1, m2, q

    def apply_ns_mass_constraint(self, mass: np.ndarray) -> np.ndarray:
        """Return the sub-array of the given mass values that are
        less than this EOS's minimum mass.

        Parameters
        ----------
        mass
            Array of mass values in solar masses

        Returns
        -------
            Sub-array of valid masses after applying the EOS minimum mass constraint
        """

        mass = mass[mass > self.min_mass]
        mass = mass[mass < self.max_mass]
        return mass
