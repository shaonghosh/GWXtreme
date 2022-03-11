from GWXtreme.parametrized_eos_sampler import mcmc_sampler

#Array Containing list of paths to the .dat files  containing the posterior samples for the events:
fnames=[]



#Name of/ Path to file in which EoS parameter posterior samples will be saved:
outname='Ap4_O3_injections'


#Initialize Sampler Object:

"""For SPectral"""

sampler=mcmc_sampler(fnames, {'gamma1':{'params':{"min":0.2,"max":2.00}},'gamma2':{'params':{"min":-1.6,"max":1.7}},'gamma3':{'params':{"min":-0.6,"max":0.6}},'gamma4':{'params':{"min":-0.02,"max":0.02}}}, outname, spectral=True, nwalkers=100, Nsamples=10000, ndim=4, spectral=True,npool=100)


"""OR"""

"""For Piece wise polytrope"""

sampler=mcmc_sampler(fnames, {'logP':{'params':{"min":33.6,"max":34.5}},'gamma1':{'params':{"min":2.0,"max":4.5}},'gamma2':{'params':{"min":1.1,"max":4.5}},'gamma3':{'params':{"min":1.1,"max":4.5}}}, '/home/anarya.ray/gwxtreme-project/repos/'+out, spectral=False, nwalkers=100, Nsamples=10000, ndim=4, spectral=True,npool=100)

#Run, Save , Plot

sampler.initialize_walkers()
sampler.run_sampler()
sampler.save_data()

fig=sampler.plot(cornerplot={'plot':True,'true vals':None},p_vs_rho={'plot':True,'true_eos':'AP4'})
fig['corner'].savefig('corner4_O3.png')
fig['p_vs_rho'][0].savefig('eos4_O3.png')




