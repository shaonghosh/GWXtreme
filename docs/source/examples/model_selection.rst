EOS Model Selection
===================

Single-Event Inference
----------------------

The below example demonstrates how to compute Bayes factors between ``lalsimulation`` tabulated equation of state models - in this case, APR4_EPP and SLY.
Parameter estimation data from the binary neutron star merger GW170817 is used for this example inference.

The necessary parameter samples depend on the variant of the inference approximation scheme being used. (See :doc:`inference methods </inference_methods>`.)
Here, the 3-dimensional posterior approximation method is chosen, which expects mass ratio (q) and dimensionless tidal deformability (Lambda1 and Lambda2) samples.

The method choice is set via the ``event_types`` argument when creating a ``ModelSelector`` - "gw-3d" in this example.

.. code:: python

   from gwxtreme.eos_inference import ModelSelector
   
   GW170817_pe_samples_file = "./GW170817_posterior_samples_with_uniform_lambda1_lambda2_prior.dat"
   
   ms = ModelSelector(posterior_files=[GW170817_pe_samples_file], event_types=['gw-3d'])
   
   bayes_factor = ms.evidence_ratio('APR4_EPP', 'SLY')

Either or both of the EOS models can also be supplied using mass-tidal deformability files or mass-radius-tidal Love number files.
A different format can be used for the *target* and *reference* models, if desired (Bayes factor = *target* model evidence / *reference* model evidence).

.. code:: python

   bayes_factor = ms.evidence_ratio(
      target_eos_mass_lambda_file='./model_A_mass_lambda_file.txt',
      reference_eos_mass_radius_k_file='./model_B_mass_radius_k_file.txt'
   )

To estimate uncertainties in computed Bayes factors, you can specify a number of repeated computations using the ``n_resamplings`` argument.
For these recomputations, the kernel density estimator (KDE) that is constructed over the event PE data and integrated to compute evidences will be resampled, producing
a slightly different value for the EOS evidences and resulting Bayes factor.

We recommended using the first, or original value for the Bayes factor as the 'best' value, rather than the mean of the recomputed values. This is because the KDE is subject to over-smoothing bias when resampling.
Likewise, because of this over-smoothing, we recommend taking the uncertainty in the estimate to be twice the standard deviation of the recomputed values.

.. code:: python

   import numpy as np

   # This may take a few minutes to run
   bayes_factors = ms.evidence_ratio(
      'H4', 'SLY', 
      n_resamplings=1000, 
      n_jobs=12,
      save_file='./H4_SLY_bayes_factors.json'
   )

   BF_best = bayes_factors[0]
   BF_error = 2 * np.std(bayes_factors[1:])

The above example also demonstrates how to run these recomputations in parallel by specifying a number of virtual CPU cores to use with the ``n_jobs`` parameter.
We use `Ray <https://docs.ray.io/en/latest/ray-overview/getting-started.html>`_ for parallel processing.

Multi-Event Joint Inference
---------------------------

Multi-event EOS inferences can be conducted by providing several posterior files and a corresponding event type for each, as shown below.

.. code:: python

   from gwxtreme.eos_inference import ModelSelector
   
   GW170817_pe_samples_file = "./GW170817_posterior_samples_with_uniform_lambda1_lambda2_prior.dat"
   GW190425_pe_samples_file = "./GW190425_posterior_samples_with_uniform_lambda_tilde_delta_lambda_tilde_prior.dat"
   
   ms = ModelSelector(
      posterior_files=[GW170817_pe_samples_file, GW190425_pe_samples_file],
      event_types=['gw-3d', 'gw-2d']
   )
   
   joint_bayes_factor = ms.evidence_ratio('APR4_EPP', 'SLY', save_file='./APR4-EPP_SLY_bayes_factors.json')

Per-event Bayes factors will also be available to view in the given JSON file.

Mult-Messenger Joint Inference
------------------------------

Multi-messenger inference can easily be performed by supplying the model selector with both gravitational wave event data and pulsar mass-radius measurement data
in the form (mass in solar masses, compactness). These pulsar events are specified using an event type of "psr".

.. code:: python
   
   GW170817_pe_samples_file = "./GW170817_posterior_samples_with_uniform_lambda1_lambda2_prior.dat"
   PSRJ0030_mass_compactness_samples_file = "./J0030_mass_compactness_samples.txt"
   
   ms = ModelSelector(
      posterior_files=[GW170817_pe_samples_file, PSRJ0030_mass_compactness_samples_file],
      event_types=['gw-3d', 'psr'],
   )
   
   joint_bayes_factor = ms.evidence_ratio('APR4_EPP', 'SLY', save_file='./APR4-EPP_SLY_bayes_factors.json')