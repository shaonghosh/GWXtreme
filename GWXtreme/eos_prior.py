# Copyright (C) 2017 Daniel Wysocki
# Assimilated from Wysocki D., O'Shaughnessy R., Bayesian Parametric Population Models, 
# 2017–, bayesian-parametric-population-models.readthedocs.io
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


import typing
import numpy
import scipy.stats
import h5py


largest_ns_mass = 1.97

def eos_p_of_rho(rho,eos):
    import lalsimulation
    FAM=lalsimulation.CreateSimNeutronStarFamily(eos)
    p_max_i = min(6e35, lalsimulation.SimNeutronStarEOSMaxPressure(eos))
    log10_p_grid = numpy.linspace(*numpy.log10([5e31, p_max_i]), 128)
    p_grid = numpy.power(10.0, log10_p_grid)
    rho_grid = numpy.empty_like(p_grid)
    for j, p in numpy.ndenumerate(p_grid):
                h = lalsimulation.SimNeutronStarEOSPseudoEnthalpyOfPressure(p, eos)
                rho_grid[j] = (
                    lalsimulation.SimNeutronStarEOSRestMassDensityOfPseudoEnthalpy(
                        h, eos,
                    )
                )


    log10_rho_grid = numpy.log10(rho_grid)
    log10_p_out = numpy.interp(
                numpy.log10(rho), log10_rho_grid, log10_p_grid,
                left=numpy.NINF, right=numpy.NINF,
            )
    return log10_p_out


def spectral_eos(eos_parameters: tuple) -> typing.Any:
    import lalsimulation

    
    gamma1, gamma2, gamma3, gamma4 = eos_parameters
    eos = lalsimulation.SimNeutronStarEOS4ParameterSpectralDecomposition(
        gamma1, gamma2, gamma3, gamma4,
    )
    

    return eos

def polytrope_eos(eos_parameters: tuple) -> typing.Any:
    import lalsimulation

    logP1, gamma1, gamma2, gamma3 = eos_parameters

    return lalsimulation.SimNeutronStarEOS4ParameterPiecewisePolytrope(
        logP1,
        gamma1, gamma2, gamma3,
    )


def spectral_eos_adiabatic_index(x,spectral_parameters):
    x_sq = x * x
    x_cu = x_sq * x

    gamma1, gamma2, gamma3, gamma4 = spectral_parameters

    log_gamma = (
        gamma1 +
        gamma2 * x +
        gamma3 * x_sq +
        gamma4 * x_cu
    )

    return numpy.exp(log_gamma)


_x_min = 0.0
_x_max = 12.3081
_x_grid = numpy.linspace(_x_min, _x_max, 500)



def is_valid_adiabatic_index(spectral_parameters: tuple):
    # Confirm that the adiabatic index is within a tolerable range.
    Gamma = spectral_eos_adiabatic_index(_x_grid, spectral_parameters)
    return (0.6 < Gamma).all() and (Gamma < 4.5).all()




def has_enough_points(eos: typing.Any) -> bool:
    import lalsimulation

    min_points = 8

    logpmin = 75.5
    logpmax = numpy.log(lalsimulation.SimNeutronStarEOSMaxPressure(eos))

    dlogp = (logpmax - logpmin) / 100

    m_prev = 0.0

    for i in range(min_points):
        p = numpy.exp(logpmin + i*dlogp)

        r, m, k = lalsimulation.SimNeutronStarTOVODEIntegrate(p, eos)

        if m <= m_prev:
            return False

        m_prev = m

    return True

def eos_max_sound_speed(eos: typing.Any, eos_fam: typing.Any) -> float:
    import lalsimulation

    # Maximum allowed mass
    m_max_kg = lalsimulation.SimNeutronStarMaximumMass(eos_fam)
    # Central pressure
    p_max = lalsimulation.SimNeutronStarCentralPressure(m_max_kg, eos_fam)
    # Pseudo-enthalpy at the core.
    h_max = lalsimulation.SimNeutronStarEOSPseudoEnthalpyOfPressure(p_max, eos)
    # Maximum sound speed should occur at the max pseudo-enthalpy for a typical
    # TOV neutron star.
    return lalsimulation.SimNeutronStarEOSSpeedOfSoundGeometerized(h_max, eos)


def is_causal_eos(eos: typing.Any, eos_fam: typing.Any) -> bool:
    c_max = eos_max_sound_speed(eos, eos_fam)

    # Confirm that the sound speed is less than speed of light, with an added
    # buffer (10%) to account for imperfect modeling.
    c_buffer = 0.1
    return c_max < (1.0 + c_buffer)


def eos_mass_range(eos_fam: typing.Any) -> typing.Tuple[float,float]:
    import lal
    import lalsimulation

    m_min = lalsimulation.SimNeutronStarFamMinimumMass(eos_fam)
    m_max = lalsimulation.SimNeutronStarMaximumMass(eos_fam)

    return m_min / lal.MSUN_SI, m_max / lal.MSUN_SI



#########################
#  Combined EoS Priors  #
#########################




def is_valid_eos(
        parameters, prior_settings,
        eos_coordinates="spectral", largest_ns_mass=largest_ns_mass,
        require_mass_ranges=None,
    ):
    import lalsimulation
    if(eos_coordinates=="spectral"):
        
        gamma1, gamma2, gamma3, gamma4 = parameters['gamma1'],parameters['gamma2'],parameters['gamma3'],parameters['gamma4']

        params_shape = gamma1.shape

        valid = (
        (gamma1 >= prior_settings["gamma1"]["params"]["min"]) &
        (gamma1 <= prior_settings["gamma1"]["params"]["max"]) &
        (gamma2 >= prior_settings["gamma2"]["params"]["min"]) &
        (gamma2 <= prior_settings["gamma2"]["params"]["max"]) &
        (gamma3 >= prior_settings["gamma3"]["params"]["min"]) &
        (gamma3 <= prior_settings["gamma3"]["params"]["max"]) &
        (gamma4 >= prior_settings["gamma4"]["params"]["min"]) &
        (gamma4 <= prior_settings["gamma4"]["params"]["max"])
    )

        for i in numpy.ndindex(*params_shape):
            # Skip already invalidated samples
            if not valid[i]:
                continue

            parameters_at_i = {
                param_name : values[i]
                for param_name, values in parameters.items()
        }
            spectral_parameters_at_i = (gamma1[i], gamma2[i], gamma3[i], gamma4[i])

            

            try:
                if not is_valid_adiabatic_index(spectral_parameters_at_i):
                    G = spectral_eos_adiabatic_index(
                        _x_grid, spectral_parameters_at_i,
                    )
                    
                    valid[i] = False
                    continue

                eos_at_i = spectral_eos(spectral_parameters_at_i)

                if not has_enough_points(eos_at_i):
                    
                    valid[i] = False
                    continue

                
                eos_fam_at_i = lalsimulation.CreateSimNeutronStarFamily(eos_at_i)
                

                if not is_causal_eos(eos_at_i, eos_fam_at_i):
                    
                    valid[i] = False
                    continue

                eos_m_min, eos_m_max = eos_mass_range(eos_fam_at_i)

                if eos_m_max < largest_ns_mass:
                    
                    valid[i] = False
                    continue

                # Check that a series of mass ranges are included within the EOS's
                # allowed mass range.
                if require_mass_ranges is not None:
                    for m_mins, m_maxs in require_mass_ranges:
                        if m_mins[i] >= eos_m_max:
                            
                            valid[i] = False
                            continue
                        if m_maxs[i] <= eos_m_min:
                            
                            valid[i] = False
                            continue

            except RuntimeError as e:
                if str(e) != "Generic failure":
                    raise
                else:
                    
                    valid[i] = False
    elif(eos_coordinates=="polytropic"):
        
        logP, gamma1, gamma2, gamma3 =  parameters["logP"],\
            parameters['gamma1'],\
            parameters['gamma2'],\
            parameters['gamma3']

        params_shape = logP.shape

        valid = (
        (logP >= prior_settings["logP"]["params"]["min"]) &
        (logP <= prior_settings["logP"]["params"]["max"]) &
        (gamma1 >= prior_settings["gamma1"]["params"]["min"]) &
        (gamma1 <= prior_settings["gamma1"]["params"]["max"]) &
        (gamma2 >= prior_settings["gamma2"]["params"]["min"]) &
        (gamma2 <= prior_settings["gamma2"]["params"]["max"]) &
        (gamma3 >= prior_settings["gamma3"]["params"]["min"]) &
        (gamma3 <= prior_settings["gamma3"]["params"]["max"])
    )

        for i in numpy.ndindex(*params_shape):
            # Skip already invalidated samples
            if not valid[i]:
                continue

            parameters_at_i = {
                param_name : values[i]
                for param_name, values in parameters.items()
        }
            polytropic_parameters_at_i = (logP[i], gamma1[i], gamma2[i], gamma3[i])

            

            try:
                

                eos_at_i = polytrope_eos(polytropic_parameters_at_i)

                if not has_enough_points(eos_at_i):
                    
                    valid[i] = False
                    continue

                
                eos_fam_at_i = lalsimulation.CreateSimNeutronStarFamily(eos_at_i)
                

                if not is_causal_eos(eos_at_i, eos_fam_at_i):
                    
                    valid[i] = False
                    continue

                eos_m_min, eos_m_max = eos_mass_range(eos_fam_at_i)

                if eos_m_max < largest_ns_mass:
                    
                    valid[i] = False
                    continue

                # Check that a series of mass ranges are included within the EOS's
                # allowed mass range.
                if require_mass_ranges is not None:
                    for m_mins, m_maxs in require_mass_ranges:
                        if m_mins[i] >= eos_m_max:
                            
                            valid[i] = False
                            continue
                        if m_maxs[i] <= eos_m_min:
                            
                            valid[i] = False
                            continue

            except RuntimeError as e:
                if str(e) != "Generic failure":
                    raise
                else:
                    
                    valid[i] = False
    else:
         raise ValueError("Does not have a known EOS coordinate system.")
    return valid


