Parameterized EOS Inference
===========================

GWXtreme uses the ``emcee`` implementation of MCMC sampling to sample the posterior distribution of EOS parameters given gravitational wave event or pulsar mass-radius data.

The below script demonstrates usage of the ``ParameterizedEoSSampler`` class to conduct this inference for the 4-parameter spectral decomposition EOS model.
Data from GW170817 and GW190425 are supplied. The 3-dimensional approximation method is chosen (gw-3d).

For these inferences, the bounds for a uniform prior distribution of the EOS parameters must be supplied. GWXtreme additionally checks during the sampling that every proposed EOS is physically valid.
For more details, consult `Ray et al. 2023 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.107.043035>`__.

.. code:: python

   import logging
   from gwxtreme.eos_inference import ParameterizedEoSSampler, load_samples

   # It is worth collecting logs from gwxtreme, especially when using the sampler, because
   # they will help you track convergence of the MCMC chain
   logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
   
   GW170817_pe_samples_file = "./GW170817_posterior_samples_with_uniform_lambda1_lambda2_prior.dat"
   GW190425_pe_samples_file = "./GW190425_posterior_samples_with_uniform_lambda1_lambda2_prior.dat"
   
   sampler = ParameterizedEoSSampler(
      posterior_files=[GW170817_pe_samples_file, GW190425_pe_samples_file],
      event_types=['gw-3d', 'gw-3d'],
      eos_prior_bounds=[(0.2, 2.0), (-1.6, 1.7), (-0.6, 0.6), (-0.02, 0.02)],
      parameterization='spectral'
   )

   sampler.run_sampler(nsteps=5_000, nwalkers=30, save_file='./samples_chain.hdf5')
   samples = load_samples('./samples_chain.hdf5', burn_in_frac=0.5, thin=5)

To continue sampling from a previous run, simply repeat the above but pass ``reset=False`` when running the sampler. The given HDF5 file will then be used to determine
where sampling should resume (i.e. the initial state for the new run will be the last sample of the existing, stored chain). New samples will be appended to the chain
in the same file.

.. code:: python
   
   sampler.run_sampler(nsteps=1_000, nwalkers=30, save_file='./samples_chain.hdf5', reset=False)

We provide convenience functions to make plots of the 90% credible level of EOS curves in pressure-density space from the posterior samples of EOS parameters.

.. code:: python

   from gwxtreme.usage_scripts import compute_eos_constraints_from_parameter_samples, plot_parameterized_eos_constraints

   samples_file = './samples_chain.hdf5'
   constraints_file = './eos_constraints.txt'

   samples = load_samples(samples_file, burn_in_frac=0.5, thin=5)

   compute_eos_constraints_from_parameter_samples(samples, 'spectral', save_file=constraints_file)

   plot_parameterized_eos_constraints(
      constraints_files=[constraints_file],  # files/labels from several runs may be provided
      labels=['GWXtreme spectral EOS constraints'],
      save_file='./constraints_plot.png',
      named_eos_list=['APR4_EPP']   # optionally provide LALSuite EOS names to include in the plot
   )