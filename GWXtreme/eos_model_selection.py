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
import json

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


def get_LambdaT_for_eos(m1, m2, max_mass_eos, eosfunc):
    '''
    This function accepts the masses and an equation of state interpolant
    with its maximum allowed mass, and return the values of LambdaT.
    '''
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

    return LambdaT


def join_json_files(list_of_jsons, nametag="model"):
    '''
    This helper function joins the JSON files from the output of the Bayes
    factor computation methods in Model_selection and Stacking into one single
    JSON file that can be directly used for further analysis.

    list_of_jsons :: A list of all the JSON files that needs to be combined.
    nametag :: Use this option to create a name-tag for the resulting JSON
    files
    '''
    alldata = []
    reference_models = []
    for json_file in list_of_jsons:
        with open(json_file, 'r') as f:
            thisdata = json.load(f)
            alldata.append(thisdata)
            reference_models.append(thisdata['ref_eos'])

    # Keep unique names of the equation of state models
    reference_models = np.unique(np.array(reference_models)).tolist()

    for model in reference_models:
        combined_dict = {}
        for data in alldata:
            if model == data['ref_eos']:
                value = {'joint_bf': data['joint_bf'],
                         'joint_bf_array': data['joint_bf_array'],
                         'all_bf': data['all_bf'],
                         'all_bf_err': data['all_bf_err']}
                combined_dict[data['target_eos']] = value

        filename = 'bayes_factors_against_' + model + '_' + nametag + '.json'
        with open(filename, 'w') as f:
            json.dump(combined_dict, f, indent=2, sort_keys=True)


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

    def getEoSInterp(self, eosname=None, m_min=1.0, N=100):
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

        if eosname is None:
            print('Allowed equation of state models are:')
            print(lalsim.SimNeutronStarEOSNames)
            print('Pass the model name as a string')
            return None
        try:
            assert eosname in list(lalsim.SimNeutronStarEOSNames)
        except AssertionError:
            print('EoS family is not available in lalsimulation')
            print('Allowed EoS are :\n' + str(lalsim.SimNeutronStarEOSNames))
            print('Make sure that if you are passing a custom file, it exists')
            print('in the path that you have provided...')
            sys.exit(0)

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

        #mass    	λ

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

    def computeEvidenceRatio(self, EoS1, EoS2, gridN=1000,
                             save=None, trials=0):
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
             max_mass_eos1] = self.getEoSInterp(eosname=EoS1,
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
             max_mass_eos2] = self.getEoSInterp(eosname=EoS2,
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
        if save:
            bf_dict = {}
            bf_dict['ref_eos'] = EoS2
            bf_dict['target_eos'] = EoS1
            bf_dict['bf'] = support2D1/support2D2
            bf_dict['bf_array'] = sup_array.tolist()
            # Making sure that the file extension is json
            if (save.split('.')[-1] != 'json') and (save.split('.')[-1] != 'JSON'):
                save += '.json'
            with open(save, 'w') as f:
                json.dump(bf_dict, f, indent=2, sort_keys=True)
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

        m1, m2, q = self.apply_mass_constraint(m1, m2, q)
        LambdaT = get_LambdaT_for_eos(m1, m2, max_mass_eos, eosfunc)

        # scale things back so they make sense with the KDE
        LambdaT_scaled, q_scaled = LambdaT/var_LambdaT, q/var_q

        # perform integration via trapazoidal approximation
        dq = np.diff(q)
        f = postfunc.evaluate(np.vstack((LambdaT_scaled, q_scaled)).T)
        f_centers = 0.5*(f[1:] + f[:-1])
        int_element = f_centers * dq

        return [LambdaT_scaled, q_scaled, np.sum(int_element)]

    def apply_mass_constraint(self, m1, m2, q):
        '''
        Apply constraints on masses based on the prior or posterior sample
        spread.
        '''
        min_mass_violation_1 = m1 < self.minMass
        min_mass_violation_2 = m2 < self.minMass
        min_mass_violation = min_mass_violation_1 + min_mass_violation_2
        m1 = m1[~min_mass_violation]
        m2 = m2[~min_mass_violation]
        q = q[~min_mass_violation]
        return (m1, m2, q)

    def plot_func(self, eos_list, gridN=1000, filename='posterior_support.png',
                  full_mc_dist=False):
        '''
        This method takes as input a list of equation of state models
        and creates a plot where these equation of state models are
        overlayed on the 2D posterior samples and their corresponding
        KDE.

        eos_list :: A list of equation state models. The members of
                    list could be either one of the named equation of
                    state in LALSimulation, or text files with columns
                    giving the mass and tidal deformability information,
                    or text files with columns giving mass, radius and
                    tidal Love number. The method also accepts a string
                    if the user wishes to plot a single equation of state.
        '''
        import pylab as pl

        pl.clf()
        pl.rcParams.update({'font.size': 18})
        pl.figure(figsize=(15, 10))

        lambdat_grid = np.linspace(0, np.max(self.data['lambdat']), 100)
        q_grid = np.linspace(np.min(self.data['q']), 1.0, 100)
        L_GRID, Q_GRID = np.meshgrid(lambdat_grid, q_grid)
        grid2D = np.array([L_GRID, Q_GRID]).T
        a, b, c = np.shape(grid2D)
        grid2D_reshaped = grid2D.reshape(a*b, c)
        sample_data = np.vstack((self.data['lambdat'], self.data['q'])).T
        kde = Bounded_2d_kde(sample_data, xlow=0.0, xhigh=None, ylow=0.0,
                             yhigh=1.0, bw=self.bw)

        support2Dgrid = kde.evaluate(grid2D_reshaped)

        support2D_matrix = support2Dgrid.reshape(len(lambdat_grid),
                                                 len(q_grid))
        pl.pcolormesh(L_GRID, Q_GRID, support2D_matrix.T)
        pl.colorbar()
        pl.scatter(self.data['lambdat'], self.data['q'], marker='.', c='k',
                   s=1, alpha=0.1)

        q_min = self.q_min*self.var_q
        q_max = self.q_max*self.var_q
        mc = np.mean(self.data['mc_source'])
        if full_mc_dist:
            mc_low = np.min(self.data['mc_source'])
            mc_hi = np.max(self.data['mc_source'])

        q = np.linspace(q_min, q_max, gridN)
        m1, m2 = getMasses(q, mc)
        if full_mc_dist:
            m1_low, m2_low = getMasses(q, mc_low)
            m1_low, m2_low, q_low = self.apply_mass_constraint(m1_low, m2_low,
                                                               q)
            m1_hi, m2_hi = getMasses(q, mc_hi)
            m1_hi, m2_hi, q_hi = self.apply_mass_constraint(m1_hi, m2_hi, q)
            q_fill = np.intersect1d(q_low, q_hi)
            m1_hi = m1_hi[np.in1d(q_hi, q_fill)]
            m2_hi = m2_hi[np.in1d(q_hi, q_fill)]
            m1_low = m1_low[np.in1d(q_low, q_fill)]
            m2_low = m2_low[np.in1d(q_low, q_fill)]
        m1, m2, q = self.apply_mass_constraint(m1, m2, q)

        assert (type(eos_list) == str or type(eos_list) == list)
        if type(eos_list) == str:
            eos_list = [eos_list]
        for eos in eos_list:
            if os.path.exists(eos):
                print('Trying m-R-k file to compute EoS interpolant')
                try:
                    [s, _, _,
                     max_mass_eos] = self.getEoSInterpFromMRFile(eos)
                except ValueError:
                    print('Trying m-λ file to compute EoS interpolant')
                    [s, _, _,
                     max_mass_eos] = self.getEoSInterpFromMLambdaFile(eos)
            else:
                [s, _, _,
                 max_mass_eos] = self.getEoSInterp(eosname=eos,
                                                   m_min=self.minMass)

            LambdaT = get_LambdaT_for_eos(m1, m2, max_mass_eos, s)
            if full_mc_dist:
                LambdaT_low = get_LambdaT_for_eos(m1_low, m2_low,
                                                  max_mass_eos, s)
                LambdaT_hi = get_LambdaT_for_eos(m1_hi, m2_hi, max_mass_eos, s)

            if full_mc_dist:
                p = pl.plot(LambdaT, q, linewidth=1, label=eos)
            else:
                p = pl.plot(LambdaT, q, linewidth=3, label=eos)
            color = p[0].get_color()
            if full_mc_dist:
                pl.fill_betweenx(q_fill, LambdaT_low, LambdaT_hi,
                                 facecolor=color, alpha=0.5)
            pl.xlabel('$\\tilde{\\Lambda}$')
            pl.ylabel('$q$')
            pl.xlim([np.min(self.data['lambdat']),
                     np.max(self.data['lambdat'])])
            pl.ylim([np.min(self.data['q']),
                     np.max(self.data['q'])])
            pl.legend()

        pl.title('EoS = {}'.format(eos_list))
        pl.savefig(filename, bbox_inches='tight')


class Stacking():
    def __init__(self, event_list, event_priors=None, labels=None):
        '''
        This class takes as input a list of posterior-samples files for
        various events. Optionally, prior samples files can also be
        supplied and allows us to compute the various quantities related
        to each of the posterior samples.
        '''
        if type(event_list) != list:  # event_list must be a list
            print('All arguments for Stacking must be a list of file-names')
            sys.exit(0)

        if event_priors:
            if type(event_priors) != list:
                print('All arguments for Stacking must be lists of file-names')
                sys.exit(0)

        if labels is None:
            labels = [None]*len(event_list)

        # Loop over the list and make sure all the paths exists.
        # Keep only those events whose file exits.

        sanitized_event_list = []
        self.labels = []
        for event, label in zip(event_list, labels):
            if os.path.exists(event):
                sanitized_event_list.append(event)
                self.labels.append(label)
            else:
                print('Could not file {}. Skipping event'.format(event))

        self.event_list = sanitized_event_list

        if event_priors:
            sanitized_event_priors = []
            for event_prior in event_priors:
                if os.path.exists(event_prior):
                    sanitized_event_priors.append(event_prior)
                else:
                    print('Could not file {}. Skipping'.format(event_prior))

            self.event_priors = sanitized_event_priors

        else:
            self.event_priors = [None]*len(self.event_list)

        # Right now the method demands a unique prior file for each event
        # This may be changed later.
        if len(self.event_priors) != len(self.event_list):
            print('Number of prior and posterior files should be same')
            sys.exit(0)

    def stack_events(self, EoS1, EoS2, trials=0, gridN=1000, save=None):
        '''
        Loop through each event and compute the joint Bayes-factor.
        Each individual event's Bayes-factor can be accessed from the Stacking
        object, using Stacking.all_bayes_factors. Uncertainty for each case can
        be accessed by using Stacking.all_bayes_factors_errors.

        EoS1 :: The name of the first equation of state model. This can be
                either one of the named equation of state models from LALSuite,
                or a file containing the information of the equation of state,
                (m, λ) or (m, r, κ).

        EoS2 :: The name of the second equation of state model. This can be
                either one of the named equation of state models from LALSuite,
                or a file containing the information of the equation of state,
                (m, λ) or (m, r, κ).

        trials :: Number of trials to be used to computed the uncertainty in
                  the Bayes-factor.

        gridN :: Number of grid points over which the line-integral is
                 computed (Default = 1000).

        save :: Use this option to save the results into json files. If nothing
                is provided, then output will not be saved. If a name is
                provided then the output will be saved to a file with that
                name.
        '''
        joint_bf = 1.0
        self.all_bayes_factors = []  # To be populated by B.Fs from all events
        if trials > 0:
            joint_bf_array = np.ones(trials)
            self.all_bayes_factors_errors = []
        for prior_file, event_file in zip(self.event_priors, self.event_list):
            modsel = Model_selection(posteriorFile=event_file,
                                     priorFile=prior_file)
            bayes_factor = modsel.computeEvidenceRatio(EoS1, EoS2,
                                                       gridN=gridN,
                                                       trials=trials)
            if type(bayes_factor) == np.float64:
                joint_bf *= bayes_factor
                self.all_bayes_factors.append(bayes_factor)
            elif type(bayes_factor) == list:
                joint_bf *= bayes_factor[0]
                self.all_bayes_factors.append(bayes_factor[0])
                this_event_error = 2*np.std(bayes_factor[-1])
                self.all_bayes_factors_errors.append(this_event_error)
                joint_bf_array *= bayes_factor[-1]

        if save is None:
            if trials > 0:
                joint_bf = [joint_bf, joint_bf_array]
            return joint_bf

        stack_dict = {}
        stack_dict['ref_eos'] = EoS2
        stack_dict['target_eos'] = EoS1

        stack_dict['joint_bf'] = joint_bf
        if trials > 0:
            stack_dict['joint_bf_array'] = joint_bf_array.tolist()
        else:
            stack_dict['joint_bf_array'] = None

        stack_dict['all_bf'] = self.all_bayes_factors
        if trials > 0:
            stack_dict['all_bf_err'] = self.all_bayes_factors_errors
        else:
            stack_dict['all_bf_err'] = None

        # Making sure that the file extension is json
        if (save.split('.')[-1] != 'json') and (save.split('.')[-1] != 'JSON'):
            save += '.json'

        with open(save, 'w') as f:
            json.dump(stack_dict, f, indent=2, sort_keys=True)

        if trials > 0:
            joint_bf = [joint_bf, joint_bf_array]

        return joint_bf

    def plot_stacked_bf(self, eos_list=None, ref_eos='SLY', trials=0,
                        gridN=1000, filename='stacked_bf.pdf'):
        '''
        This method makes bar plots for bayes-factor between various EoS
        and a reference EoS. It does this for multiple events. The bar plots
        are generated for each event. The results are stacked and the then
        a combined bayes-factor bar-plot is generated. Alternatively, this
        method can also be used to make plots directly from data files.

        eos_list :: List of strings, either named equation of state models from
                   LALSuite, or names of files contining the equation of state
                   information, (m, λ) or (m, r, κ). If no list is provided, a
                   default list will be used.

        ref_eos :: Reference equation of state The equation of state against
                   which the Bayes-factor is to be computed. If no reference
                   model is used, the SLY model from LALSuite will be used.

        trials :: Number of trials to be used to computed the uncertainty in
                  the Bayes-factor.

        gridN :: Number of grid points over which the line-integral
                 computed (Default = 1000).

        filename :: Name of the file where the plot will be saved.
        '''

        import pylab as pl
        pl.clf()

        if eos_list is None:
            eos_list = ['APR4_EPP', 'BHF_BBB2', 'H4', 'HQC18',
                        'KDE0V', 'KDE0V1', 'MPA1', 'MS1B_PP',
                        'MS1_PP', 'RS', 'SK255', 'SK272',
                        'SKI2', 'SKI3', 'SKI4', 'SKI5', 'SKI6',
                        'SKMP', 'SKOP', 'SLY9', 'WFF1']

        N = len(eos_list)
        ind = np.arange(N)
        width = 0.10

        bf_combined = []
        d_bf_combined = []
        bf_all_events = []
        d_bf_all_events = []
        for eos in eos_list:
            print('Stacking events for model: {}'.format(eos))
            this_eos_bf = self.stack_events(eos, ref_eos, trials=trials)
            if trials > 0:
                bf_combined.append(this_eos_bf[0])
                d_bf_combined.append(2*np.std(this_eos_bf[-1]))
            else:
                bf_combined.append(this_eos_bf)

            bf_all_events.append(self.all_bayes_factors)
            if trials > 0:
                d_bf_all_events.append(self.all_bayes_factors_errors)

        bf_all_events = np.array(bf_all_events).T
        d_bf_all_events = np.array(d_bf_all_events).T

        if trials > 0:
            pl.bar(ind, bf_combined, width, yerr=d_bf_combined,
                   error_kw=dict(lw=3, capsize=6, capthick=2),
                   label="Joint Bayes' factor")
        else:
            pl.bar(ind, bf_combined, width, yerr=None,
                   label="Joint Bayes' factor")

        shift = 1
        if trials > 0:
            for bf, dbf, ll in zip(bf_all_events, d_bf_all_events, self.labels):  # noqa E501
                pl.bar(ind + shift*width, bf, width, yerr=dbf, capsize=6,
                       label=ll, alpha=0.4)
                pl.xticks(ind + 2.5*width, eos_list, rotation=50)
                shift += 1

        else:
            for bf, ll in zip(bf_all_events, self.labels):
                pl.bar(ind + shift*width, bf, width, label=ll, alpha=0.4)
                pl.xticks(ind + 2.5*width, eos_list, rotation=50)
                shift += 1

        pl.legend(loc='best')
        ax = pl.gca()
        pl.ylim([0, 1.9])
        pl.ylabel('Bayes-Factor w.r.t {}'.format(ref_eos))
        pl.savefig(filename, bbox_inches='tight')
