EOS Model Selection
===================

Single-Event Inference
----------------------

The below example demonstrates how to compute Bayes factors between ``lalsimulation`` tabulated equation of state models - in this case, APR4_EPP and SLY.
Parameter estimation data from the binary neutron star merger GW170817 is used for this example inference.

The necessary parameter samples depend on the variant of the inference approximation scheme being used.
Here, the 2-dimensional posterior approximation method is chosen, which expects mass ratio (q) and binary tidal deformability (Lambda tilde) samples.

The method choice is set via the ``event_type`` argument when creating a ``ModelSelector`` - "gw-2d" in this example.

.. code:: python

   from gwxtreme.eos_inference import ModelSelector
   
   GW170817_pe_samples_file = "/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/cbc_pe_samples/GW170817/posterior_samples/GW170817_posterior_samples_broad_spin_prior.dat"
   
   ms = ModelSelector(posterior_file=GW170817_pe_samples_file, event_type='gw-2d')
   
   bayes_factor = ms.compute_eos_evidence_ratio('APR4_EPP', 'SLY')

Computing the evidence of a given EOS model from the event data involves integrating a probability density function over the event parameter space.
By default, a simple (but slow) kernel density estimator is used to construct the necessary density function from the given parameter estimation samples.

Faster density estimation (1 - 2 orders of magnitude) can be achieved using a normalizing flow, which is a deep neural network architecture suited for rapid density estimation.
GWXtreme supports the use of Bayesian normalizing flows, which have the added benefit of learning a distribution of model weights during training, allowing for uncertainty estimation by drawing many instances of the model from this distribution and averaging the results.

The below code snippet shows how a Bayesian normalizing flow can be used during EOS inference with GWXtreme. 
A PyTorch tensor file containing the model weights and configuration for a pre-trained flow must be provided. These model files, pre-trained for various GW events and pulsar mass-radius datasets, can be found ...

.. code:: python

   ms = ModelSelector(
      posterior_file=GW170817_pe_samples_file, 
      event_type='gw-2d',
      density_est_method='flow',
      flow_file='/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/density_estimators/GW170817/GW170817_gw-2d_TaylorF2.pt'
   )

Multi-Event Joint Inference
---------------------------

Multi-event EOS inferences can be conducted using ``JointModelSelector``. Provide several posterior samples files and a corresponding event type for each, as shown below.
If the normalizing flow density estimation method is used, several flow model files must be provided, corresponding to the several events.

.. code:: python

   from gwxtreme.eos_inference import JointModelSelector
   
   GW170817_pe_samples_file = "/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/cbc_pe_samples/GW170817/posterior_samples/GW170817_posterior_samples_broad_spin_prior.dat"
   GW190425_pe_samples_file = "/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/cbc_pe_samples/GW190425/posterior_samples/posterior_samples_GW190425_broad_spin_prior.dat"
   
   jms = JointModelSelector(
      posterior_files=[GW170817_pe_samples_file, GW190425_pe_samples_file],
      event_types=['gw-2d', 'gw-2d'],
      density_est_method='flow',
      flow_files=[
         '/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/density_estimators/GW170817/GW170817_gw-2d_TaylorF2.pt',
         '/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/density_estimators/GW190425/GW190425_gw-2d_TaylorF2.pt'
      ]
   )
   
   joint_bayes_factor = jms.compute_joint_eos_evidence_ratio('APR4_EPP', 'SLY')

Mult-Messenger Joint Inference
------------------------------

Multi-messenger inference can easily be performed by supplying the model selector with both gravitational wave event data and pulsar mass-radius measurement data
in the form (mass in solar masses, compactness). If using a flow density estimator, it should be trained on (mass in solar masses, radius in km) data.
These pulsar events are passed with event type "psr".

.. code:: python

   from gwxtreme.eos_inference import JointModelSelector
   
   GW170817_pe_samples_file = "/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/cbc_pe_samples/GW170817/posterior_samples/GW170817_posterior_samples_broad_spin_prior.dat"
   PSRJ0030_mass_compactness_samples_file = "/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/nicer_pulsar_samples/N50k_J0030_3spot_CM.txt"
   
   jms = JointModelSelector(
      posterior_files=[GW170817_pe_samples_file, PSRJ0030_mass_compactness_samples_file],
      event_types=['gw-2d', 'psr'],
      density_est_method='flow',
      flow_files=[
         '/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/density_estimators/GW170817/GW170817_gw-2d_TaylorF2.pt',
         '/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/density_estimators/PSRJ0030/PSRJ0030.pt'
      ]
   )
   
   joint_bayes_factor = jms.compute_joint_eos_evidence_ratio('APR4_EPP', 'SLY')

Model Selection with a Set of Equations
---------------------------------------

The script below will compute Bayes factors for a set of tabulated LALSuite models with respect to the SLY model.
To allow for estimating uncertainties in the results, the Bayes factor calculation for each EOS is repeated 1000 times (specified with ``n_resamplings = 1000``).
This is accomplished by re-sampling the density estimator and re-integrating the resulting density. For a KDE, re-sampling entails sampling points from the KDE and then constructing a new KDE from this sample.
For a Bayesian flow, re-sampling simply involves drawing a new set of model weights from the posterior learned during training.

.. code:: python

   from gwxtreme.eos_inference import ModelSelector
   from gwxtreme.usage_scripts import compute_bayes_factors_for_several_named_eos

   GW170817_pe_samples_file = "/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/cbc_pe_samples/GW170817/posterior_samples/GW170817_posterior_samples_broad_spin_prior.dat"

   ms = ModelSelector(
      posterior_file=GW170817_pe_samples_file, 
      event_type='gw-2d',
      density_est_method='flow',
      flow_file='/home/joseph/LocalProjects/gwxtreme_work/using_gwxtreme/density_estimators/GW170817/GW170817_gw-2d_TaylorF2.pt'
   )

   eos_list = [
      "APR4_EPP",
      "HQC18",
      "SKOP",
      "MPA1",
      "SKI4",
      "SKI6",
      "SKMP",
      "SK272",
      "SK255",
      "RS",
      "SKI3",
      "SKI2",
      "SKI5",
      "H4",
      "MS1B_PP",
      "MS1_PP",
   ]

   compute_bayes_factors_for_several_named_eos(
      model_selector=ms,
      named_eos_list=eos_list,
      reference_eos_name='SLY',
      n_grid=1000,
      n_resamplings=1000,
      save_file='./results.json',
   )

The below code can be used to plot a bar chart of the computed Bayes factors. The determination of the height of these bars, which capture the "best" value for the Bayes factor for each equation, depends on the type of density estimator used.

For a Bayesian flow (this example), the best value for the Bayes factor is simply the mean of all computed values, including the original and the re-sampled ones (and the uncertainty can be taken as the standard deviation of all values).
For a KDE, the best value would be the original, or first computed Bayes factor, since the re-samplings of the KDE introduce bias due to over-smoothing. The uncertainty could then come from the standard deviation of just the re-sampled values.

.. code:: python

   import json
   from gwxtreme.usage_scripts import plot_bayes_factors_bar_chart

   with open('./results.json') as f:
      bayes_factors = json.load(f)

   # An outer label ("GW170817") is given to match the format which
   # the below plotting function expects. For comparisons between 
   # several runs/events, results for each event can be added to this 
   # dict in a similar fashion.
   bar_chart_heights = {"GW170817": {}}
   bar_chart_error_bars = {"GW170817": {}}

   for eos, bfs in bayes_factors.items():
      bf_best = np.mean(bfs)
      bf_error = np.std(bfs)

      bar_chart_heights["GW170817"][eos] = bf_best
      bar_chart_error_bars["GW170817"][eos] = bf_error

   plot_bayes_factors_bar_chart(
      bayes_factor_bar_heights=bar_chart_heights,
      bayes_factor_error_bars=bar_chart_error_bars,
      save_file='./GW170817_results.png',
   )