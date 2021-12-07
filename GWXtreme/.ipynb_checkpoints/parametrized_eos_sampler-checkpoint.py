# Copyright (C) 2021 Anarya Ray
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


import numpy as np
from .eos_model_selection import Stacking
import h5py
import corner
from multiprocessing import cpu_count, Pool
import time
import emcee as mc
from pop_models.applications.bns.iid_gaussian_mass_with_eos.prior import is_valid_eos
from pop_models.astro_models import eos






class mcmc_sampler():
    def __init__(self, posterior_files, prior_bounds, outfile, nwalkers=100, Nsamples=10000, ndim=4, spectral=True,npool=1):
        '''
        Initiates Parametric EoS mcmc Sampler Class
        that also stacks over multiple events,from the
        single event, uniform in LambdaT, dLambdaT parameter
        estimation runs. Parametrization chosen: 4 parameter 
        Spectal decomposition of Adiabatic Index in terms of
        Pressure or 4 parameter piec-wise polytrope.
            
        prior_bounds :: dictionary containining prior bounds
                            of spectral parameters
                              
        posterior_files :: Array containing full paths to 
                                single event PE posterior files
                                    
        nwalkers  ::  Number of walkers to use for mcmc
            
        Nsamples  ::   Numper of samples to draw using mcmc
                        to infer the joint posterior of spectral
                        parameters
                           
            
        out_file  ::   .h5 File name that will store samples
            
        
        npool     ::    Number of cpu's to use for parallelization
        '''
        
        self.posteriorfiles=posterior_files
        self.priorbounds=prior_bounds
        self.outfile=outfile+'.h5'
        self.nwalkers=nwalkers
        self.nsamples=Nsamples
        self.ndim=ndim
        self.spectral=spectral
        self.eosmodel=Stacking(posterior_files,spectral=spectral)
        self.npool=npool
    def log_post(self,p):
        '''
        This mathod accepts an array of spectral parameters
        and returns their log posterior given gw data from all
        the events provided while initializing the class.
        
        p  :: array of spectral parameters.
        
        gridN :: Number of grid points to use to perform integral
                over q
        '''
        
        params={"gamma{}".format(i):np.array([p[i-1]]) for i in range(1,len(p)+1)}
        tr=is_valid_eos(params,self.priorbounds)
        if not tr:
            return -np.inf
        return np.nan_to_num(np.log(self.eosmodel.joint_evidence(p,gridN=100)))
    
    def initialize_walkers(self):
        
        '''
        This method initializes the walkers for mcmc 
        (to be run on the spectral parameters posterior
        given GW data from all the events)inside the prior 
        region.
        '''
        
        n=0
        p0=[]
        while True:
            g=np.array([np.random.uniform(self.spectral_bounds["gamma{}".format(i)]["params"]["min"],self.spectral_bounds["gamma{}".format(i)]["params"]["max"]) for i in range(1,5)])
            params={"gamma{}".format(i):np.array([g[i-1]]) for i in range(1,len(g)+1)}
    
            if(is_valid_eos(params,self.priorbounds)):
                p0.append(g)
                n+=1
            if(n>=self.nwalkers):
                break
            
        self.p0=p0
    
    
    def run_sampler(self):
        '''
        runs mcmc sampler to draw samples of 
        the spectral parameters from their
        posterior given GW data from all the 
        events.
        '''
        
        if (self.npool>1):
            
            with Pool(min(cpu_count(),self.npool)) as pool:
                
                sampler=mc.EnsembleSampler(self.nwalkers,self.ndim,self.log_post,pool=pool)
                sampler.run_mcmc(self.p0,self.nsamples,progress=True)
        
                self.samples=sampler.get_chain()
                self.logp=sampler.get_log_prob()
        else:
              sampler=mc.EnsembleSampler(self.nwalkers,self.ndim,self.log_post)
              sampler.run_mcmc(self.p0,self.nsamples,progress=True)
        
              self.samples=sampler.get_chain()
              self.logp=sampler.get_log_prob() 
    
    def save_data(self):
        '''
        This method saves the posterior samples 
        of the spectral parameters and there corresponding
        log probablities
        '''
        f=h5py.File(self.outfile,'w')
        f.create_dataset('chains',data=np.array(self.samples))
        f.create_dataset('logp',data=np.array(self.logp))
        f.close()
        

        
    def plot(self,cornerplot={'plot':False,'true vals':None},p_vs_rho=False):
        '''
        This method plots the posterior of the spectral
        parameters in a corner plot and also the pressure
        density credible intervals inferred from these
        posteriors in a separate plot. It reurns an array
        containing matplotlib.pyplot.figure and axes objects
        corresponding to each plot
        
        cornerplot :: dictionary saying whether to plot corner or
                      not and whether or not to mark true values
                      (if any) on the plot.
                      deafult is
                      {'plot':False,'true vals':None}
        
        p_vs_rho   :: Mentions whether or not to plot 
                      Pressure Density Credible intervals
                      (default is False)
        '''
        fig={'corner':None,'p_vs_rho':None}
        Ns=self.samples.shape
        burn_in=int(Ns[0]/2.)
        samples=[]
        thinning=int(Ns[0]/50.)
        
        try:
            thinning=int(max(mc.autocorr.integrated_time(self.samples))/2.)
        except mc.autocorr.AutocorrError as e:
            print(e)
        for i in range(burn_in,Ns[0],thinning):
            for j in range(Ns[1]):
                samples.append(self.samples[i,j,:])

        samples=np.array(samples)
        
        if cornerplot['plot']:
                                                                     
            fig_corner=corner.corner(samples,labels=[r'$\gamma_0$',r'$\gamma_1$',r'$\gamma_2$',r'$\gamma_3$'],smooth=1.2,quantiles=[0.1,0.5,0.9],color='b',show_titles=True,title_kwargs={"fontsize": 12,"color":'b'},use_math_text=True,label='GWXtreme',hist_kwargs={"density":True})
            if(cornerplot['true vals'] is not None):
                ndim=4
                axes=np.array(fig_corner.axes).reshape((ndim,ndim))

                Tr=cornerplot['true vals']
        
                for i in range(ndim):
                    ax=axes[i,i]
                    ax.axvline(x=Tr[i],color='orange')
    
    
                for yi in range(ndim):
                    for xi in range(yi):
                        ax=axes[yi,xi]
                        ax.axvline(x=Tr[xi],color='orange')
                        ax.axhline(y=Tr[yi],color='orange')
                        ax.plot(Tr[xi],Tr[yi],color='orange')
                fig['corner']=fig_corner
                        
        if(p_vs_rho):
            logp=[]
            rho=np.logspace(17.25,18.25,1000)
            

            for s in samples:
                params=(s[0],s[1],s[2],s[3])
                lp=np.zeros(len(rho))
                p=eos.spectral_eos_p_of_rho(rho,params)
                arg1=np.where(p>0.)
                lp[arg1]=np.log10(p[arg1])
                logp.append(lp)
    
            logp=np.array(logp)
            logp_CIup=np.array([np.quantile(logp[:,i],0.95) for i in range(len(rho))])
            logp_CIlow=np.array([np.quantile(logp[:,i],0.05) for i in range(len(rho))])
            logp_med=np.array([np.quantile(logp[:,i],0.5) for i in range(len(rho))])
            fig_eos,ax_eos
            ax_eos.errorbar(np.log10(rho),logp_med,color='cyan',yerr=[logp_med-logp_CIlow,logp_CIup-logp_med],elinewidth=2.0,capsize=1.5,ecolor='cyan',fmt='')
            ax_eos.set_xlabel(r'$\log10{\frac{\rho}{g cm^-3}}$',fontsize=20)
            ax_eos.set_ylabel(r'$log10(\frac{p}{dyne cm^{-2}})$',fontsize=20)
            fig['p_vs_rho']=(fig_eos,ax_eos)
        return fig