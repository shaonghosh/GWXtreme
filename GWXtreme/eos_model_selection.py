# Copyright (C) 2018 Shaon Ghosh, Xiaoshu Liu
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


from __future__ import division, print_function

import os
import sys

import numpy as np
from scipy.interpolate import interp1d

import lal
import lalsimulation as lalsim

from .bounded_2d_kde import Bounded_2d_kde


def getMasses(q, mc):
    '''
    Given chirp-mass and mass ratio, compute the individual masses
    '''
    m1 = mc * (1 + q)**(1/5) * (q)**(-3/5)
    m2 = mc * (1 + q)**(1/5) * (q)**(2/5)
    return (m1, m2)


def getLambdaT(m1, m2, Lambda1, Lambda2):
    '''
    This function converts the Lambda1, Lambda2, mass1, mass2
    to Lambda tilde.
    '''
    LambdaTilde = (2**5)/(26*(m1 + m2)**5)
    LambdaTilde *= (m1**5 + 12*m2*m1**4)*Lambda1 +\
                   (m2**5 + 12*m1*m2**4)*Lambda2
    return LambdaTilde


class Model_selection:
    def __init__(self, posteriorFile, priorFile=None):
        '''
        Initiates the Bayes factor calculator with the posterior
        samples from the uniform LambdaT, dLambdaT parameter
        estimation runs.

        posteriorFile :: The full path to the posterior_samples.dat
                         file

        priorFile     :: The full path to the priors file (optional).
                         If the prior file is supplied, the mass
                         boundaries for the KDE computation is
                         obtained from the prior file. If this is not
                         supplied, the posterior samples will be used
                         to determine the bounds.
        '''
        self.data = np.recfromtxt(posteriorFile, names=True)

        if priorFile:
            self.prior = np.recfromtxt(priorFile, names=True)
            self.minMass = np.min(self.prior['m2_source'])
            self.maxMass = np.max(self.prior['m1_source'])
            self.q_max = np.max(self.prior['q'])
            self.q_min = np.min(self.prior['q'])
        else:
            self.prior = None
            self.minMass = np.min(self.data['m2_source'])  # min posterior mass
            self.maxMass = np.max(self.data['m1_source'])  # max posterior mass
            self.q_max = np.max(self.data['q'])
            self.q_min = np.min(self.data['q'])

        # store useful parameters
        self.mc_mean = np.mean(self.data['mc_source'])

        # whiten data
        self.var_LambdaT = np.std(self.data['lambdat'])
        self.var_q = np.std(self.data['q'])

        self.q_max /= self.var_q
        self.q_min /= self.var_q
        self.yhigh = 1.0/self.var_q  # For reflection boundary condition

        self.margPostData = np.vstack((self.data['lambdat']/self.var_LambdaT,
                                       self.data['q']/self.var_q)).T
        self.bw = len(self.margPostData)**(-1/6.)  # Scott's bandwidth factor

        # Compute the KDE for the marginalized posterior distribution #
        self.kde = Bounded_2d_kde(self.margPostData,
                                  xlow=0.0,
                                  xhigh=None,
                                  ylow=0.0,
                                  yhigh=self.yhigh,
                                  bw=self.bw)

    def getEoSInterp(self, eosname, m_min=1.0, N=100):
        '''
        This method accepts one of the NS native equations of state
        and uses that to return a list [s, mass, Λ, max_mass] where
        s is the interpolation function for the mass and the tidal
        deformability.

        eosname     :: Equation of state native to LALsuite

        m_min       :: The minimum mass of the NS from which value
                       the interpolant will be constructed
                       (default = 1.0).

        N           :: Number of points that will be used for the
                       construction of the interpolant.
        '''

        allowedEoS = ['ALF1', 'ALF2', 'ALF3', 'ALF4', 'AP1', 'AP2', 'AP3',
                      'AP4', 'APR4_EPP', 'BBB2', 'BGN1H1', 'BPAL12', 'BSK19',
                      'BSK20', 'BSK21', 'ENG', 'FPS', 'GNH3', 'GS1', 'GS2',
                      'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'MPA1', 'MS1B',
                      'MS1B_PP', 'MS1_PP', 'MS1', 'MS2', 'PAL6', 'PCL2', 'PS',
                      'QMC700', 'SLY4', 'SLY', 'SQM1', 'SQM2', 'SQM3', 'WFF1',
		      'WFF2', 'WFF3', 'APR', 'BHF_BBB2', 'KDE0V', 'KDE0V1',
                      'RS', 'SK255', 'SK272', 'SKA', 'SKB', 'SKI2', 'SKI3',
                      'SKI4', 'SKI5', 'SKI6', 'SKMP', 'SKOP', 'SLY2',
                      'SLY230A', 'SLY9', 'HQC18']

        try:
            assert eosname in allowedEoS
        except AssertionError:
            print('EoS family is not available in lalsimulation')
            print('Allowed EoS are :\n' + str(allowedEoS))
            sys.exit(1)
        eos = lalsim.SimNeutronStarEOSByName(eosname)
        fam = lalsim.CreateSimNeutronStarFamily(eos)
        max_mass = lalsim.SimNeutronStarMaximumMass(fam)/lal.MSUN_SI

        # This is necessary so that interpolant is computed over the full range
        # Keeping number upto 3 decimal places
        # Not rounding up, since that will lead to RuntimeError
        max_mass = int(max_mass*1000)/1000
        masses = np.linspace(m_min, max_mass, N)
        masses = masses[masses <= max_mass]
        Lambdas = []
        gravMass = []
        for m in masses:
            try:
                rr = lalsim.SimNeutronStarRadius(m*lal.MSUN_SI, fam)
                kk = lalsim.SimNeutronStarLoveNumberK2(m*lal.MSUN_SI, fam)
                cc = m*lal.MRSUN_SI/rr
                Lambdas = np.append(Lambdas, (2/3)*kk/(cc**5))
                gravMass = np.append(gravMass, m)
            except RuntimeError:
                break
        Lambdas = np.array(Lambdas)
        gravMass = np.array(gravMass)
        s = interp1d(gravMass, Lambdas)
        return [s, gravMass, Lambdas, max_mass]

    def getEoSInterpFromMLambdaFile(self, tidalFile):
        '''
        This method accepts the data from a file that have the
        tidal deformability information in the following format:

        #mass		λ

        ...	    	...

        ...		    ...

        max_mass	...

        The values of masses should be in units of solar masses. The
        tidal deformability λ should be supplied in SI unit.

        The method computes the dimensionless tidal deformabiliy Λ and
        returns a list [s, mass, Λ, max_mass] where s is the interpolation
        function for the mass and the tidal deformability.
        '''
        masses, lambdas = np.loadtxt(tidalFile, unpack=True)
        self.minMass = np.min(masses)
        Lambdas = lal.G_SI*lambdas*(1/(lal.MRSUN_SI*masses)**5)
        s = interp1d(masses, Lambdas)
        max_mass = np.max(masses)
        return [s, masses, Lambdas, max_mass]


    def getEoSInterpFromMRFile(self, MRFile):
        '''
        This method accepts the data from a file that have the
        mass-radius-love deformability information in the following format:

        #mass		Radius       love_num

        ...	    	...          ...

        ...		    ...          ...

        max_mass	...          ...

        The values of masses should be in units of solar masses. The
        tidal deformability radius should be supplied in meters.

        The method computes the dimensionless tidal deformabiliy Λ and
        returns a list [s, mass, Λ, max_mass] where s is the interpolation
        function for the mass and the tidal deformability.
        '''
        masses, radius, kappa = np.loadtxt(MRFile, unpack=True)
        self.minMass = np.min(masses)
        compactness = masses*lal.MRSUN_SI/radius
        Lambdas = (2/3)*kappa/(compactness**5)
        s = interp1d(masses, Lambdas)
        max_mass = np.max(masses)
        return [s, masses, Lambdas, max_mass]


    def computeEvidenceRatio(self, EoS1, EoS2,
                             gridN=1000, trials=0):
        '''
        This method computes the ratio of evidences for two
        tabulated EoS. It first checks if a file exists with
        the name associated with the strings EoS1 and EoS2.
        If it does, then use method getEoSInterpFromFile,
        else use method getEoSInterp. This computation is
        conducted for multiple trials to get an estimation of
        the uncertainty.

        EoS1    :: The name of the first tabulated equation of
                   state or the name of the file from which the
                   EoS data is to be read.
        EoS2    :: The name of the second tabulated equation of
                   state or the name of the file from which the
                   EoS data is to be read.
        gridN   :: Number of grid points over which the
                   line-integral is computed. (Default = 1000)
        trials  :: Number of trials for estimating the
                   uncertainty in the Bayes-factor.
        '''

        # generate interpolators for both EOS

        if os.path.exists(EoS1):
            print('Trying m-R-k file to compute EoS interpolant')
            try:
                [s1, _, _,
                 max_mass_eos1] = self.getEoSInterpFromMRFile(EoS1)
            except ValueError:
                print('Trying m-λ file to compute EoS interpolant')
                [s1, _, _,
                 max_mass_eos1] = self.getEoSInterpFromMLambdaFile(EoS1)
        else:
            [s1, _, _,
             max_mass_eos1] = self.getEoSInterp(EoS1,
                                                m_min=self.minMass)

        if os.path.exists(EoS2):
            print('Trying m-R-k file to compute EoS interpolant')
            try:
                [s2, _, _,
                 max_mass_eos2] = self.getEoSInterpFromMRFile(EoS2)
            except ValueError:
                print('Trying m-λ file to compute EoS interpolant')
                [s2, _, _,
                 max_mass_eos2] = self.getEoSInterpFromMLambdaFile(EoS2)
        else:
            [s2, _, _,
             max_mass_eos2] = self.getEoSInterp(EoS2,
                                                m_min=self.minMass)

        # compute support
        [lambdat_eos1,
         q_eos1, support2D1] = self.integrator(self.q_min, self.q_max,
                                               self.mc_mean, s1,
                                               max_mass_eos1, self.kde,
                                               gridN=gridN,
                                               var_LambdaT=self.var_LambdaT,
                                               var_q=self.var_q)

        [lambdat_eos2,
         q_eos2, support2D2] = self.integrator(self.q_min, self.q_max,
                                               self.mc_mean, s2,
                                               max_mass_eos2, self.kde,
                                               gridN=gridN,
                                               var_LambdaT=self.var_LambdaT,
                                               var_q=self.var_q)

        # iterate to determine uncertainty via re-drawing from
        # smoothed distribution
        # NOTE: this is known to introduce a non-zero bias into
        # the mean and variance estimate!
        support2D1_list = []
        support2D2_list = []
        if trials == 0:
            return (support2D1/support2D2)
        for ii in range(trials):

            # generate new (synthetic) data
            new_margPostData = np.array([])
            counter = 0
            while len(new_margPostData) < len(self.margPostData):
                prune_adjust_factor = 1.1 + counter/10.
                N_resample = int(len(self.margPostData)*prune_adjust_factor)
                new_margPostData = self.kde.resample(size=N_resample).T
                unphysical = (new_margPostData[:, 0] < 0) +\
                             (new_margPostData[:, 1] > self.yhigh) +\
                             (new_margPostData[:, 1] < 0)
                new_margPostData = new_margPostData[~unphysical]
                counter += 1
            indices = np.arange(len(new_margPostData))
            chosen = np.random.choice(indices, len(self.margPostData))
            new_margPostData = new_margPostData[chosen]

            # generate a new kde
            new_kde = Bounded_2d_kde(new_margPostData, xlow=0.0,
                                     xhigh=None, ylow=0.0,
                                     yhigh=self.yhigh,
                                     bw=self.bw)

            # integrate to get support
            [this_lambdat_eos1, this_q_eos1,
             this_support2D1] = self.integrator(self.q_min, self.q_max,
                                                self.mc_mean, s1,
                                                max_mass_eos1, new_kde,
                                                gridN=gridN,
                                                var_LambdaT=self.var_LambdaT,
                                                var_q=self.var_q)
            [this_lambdat_eos2, this_q_eos2,
             this_support2D2] = self.integrator(self.q_min, self.q_max,
                                                self.mc_mean, s2,
                                                max_mass_eos2, new_kde,
                                                gridN=gridN,
                                                var_LambdaT=self.var_LambdaT,
                                                var_q=self.var_q)
            # store the result
            support2D1_list.append(this_support2D1)
            support2D2_list.append(this_support2D2)

        sup_array = np.array(support2D1_list)/np.array(support2D2_list)
        return [support2D1/support2D2, sup_array]

    # The integrator function #
    def integrator(self, q_min, q_max, mc, eosfunc,
                   max_mass_eos, postfunc,
                   gridN=1000, var_LambdaT=1.0, var_q=1.0):
        '''
        This function numerically integrates the KDE along the
        EoS curve.

        q_min  	:: Minimum value of mass-ratio for the EoS curve

        q_max  	:: Maximum value of mass-ratio for the EoS curve

        mc	:: Chirp mass (fixed to mean value of posterior)

        eosfunc	 :: interpolation function of Λ = eosfunc(m)

        max_mass_eos :: Maximum mass allowed by the EoS.

        postfunc :: K(Λ, q) KDE of the posterior distr

        gridN  :: Number of steps for the integration (default=1K)

        var_LambdaT :: Standard deviation of the LambdT

        var_q  :: Standard deviation of the mass-ratio

        If for the choice of mass-ratio and mc, the masses of one
        or both the object goes above the maximum mass of NS
        allowed by the EoS, then the object(s) is(are) treated as
        BH (Λ=0). If the masses are below the minimum mass, the
        points are excludeds from the integral.

        '''
        # scale these appropriately to evaluate prior bounds
        q_min *= var_q
        q_max *= var_q

        # get values for line integral
        q = np.linspace(q_min, q_max, gridN)
        m1, m2 = getMasses(q, mc)

        # apply constraints
        min_mass_violation_1 = m1 < self.minMass
        min_mass_violation_2 = m2 < self.minMass
        min_mass_violation = min_mass_violation_1 + min_mass_violation_2
        m1 = m1[~min_mass_violation]
        m2 = m2[~min_mass_violation]
        q = q[~min_mass_violation]

        # determine when bodies are black holes
        kerr_cases_1 = m1 >= max_mass_eos
        kerr_cases_2 = m2 >= max_mass_eos

        Lambda1 = np.zeros_like(m1)
        Lambda2 = np.zeros_like(m2)

        # interpolate from known curves to obtain tidal
        # deformabilities as a function of mass for the
        # rest of the points
        Lambda1[~kerr_cases_1] = eosfunc(m1[~kerr_cases_1])
        Lambda2[~kerr_cases_2] = eosfunc(m2[~kerr_cases_2])

        # compute chirp tidal deformability
        LambdaT = getLambdaT(m1, m2, Lambda1, Lambda2)

        # scale things back so they make sense with the KDE
        LambdaT_scaled, q_scaled = LambdaT/var_LambdaT, q/var_q

        # perform integration via trapazoidal approximation
        dq = np.diff(q)
        f = postfunc.evaluate(np.vstack((LambdaT_scaled, q_scaled)).T)
        f_centers = 0.5*(f[1:] + f[:-1])
        int_element = f_centers * dq

        return [LambdaT_scaled, q_scaled, np.sum(int_element)]
