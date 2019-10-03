.. GWXtreme documentation master file, created by
   sphinx-quickstart on Tue Jul 23 14:30:07 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation of GWXtreme
=========================

Introduction:
-------------
The GWXtreme package hosts the tools to compute the Bayes' factor for two equations of 
state of neutron star from gravitational-wave (GW) data. As an input this requires 
posterior samples of masses and tidal-deformability obtained from parameter estimation 
of GW data for.

In this current implementation, the code requires a careful choice of priors for the 
parameter estimation (PE) runs to be conducted. Once the posterior samples are available 
the GWXtreme allows the user to use one of the native neutron star equations of state, 
or a supply their favorite equation of state via a text file to compute the Bayes' factor 
for the various models.

The Prior:
----------
Currently, this package is equipped to handle posterior samples from only a specific type 
of prior in the PE runs. It demands that the PE is conducted using a uniform prior in the
modified tidal-deformability parameters :math:`(\tilde{\Lambda}, \delta\tilde{\Lambda})`,
where these are defined as:

.. math::
   :nowrap:

   \begin{eqnarray}
      \tilde{\Lambda} & = & \frac{8}{13}\left[(1 + 7\eta - 31\eta^2)(\Lambda_1 + \Lambda_2) + \sqrt{1 - 4\eta}\,(1 + 9\eta - 11\eta^2)(\Lambda_1 - \Lambda_2)\right] \\
      \delta\tilde{\Lambda} & = & \frac{1}{2}\left[\sqrt{1 - 4\eta}\left(1 - \frac{13272}{1319}\eta + \frac{8944}{1319}\eta^2\right)(\Lambda_1 + \Lambda_2) + \left(1 - \frac{15910}{1319}\eta + \frac{32850}{1319}\eta^2 + \frac{3380}{1319}\eta^3\right)(\Lambda_1 - \Lambda_2)\right]
   \end{eqnarray}

The choice of this prior has an implication. It leads to some of the points having unphysical
tidal-deformabilities :math:`(\Lambda_1, \Lambda_2)`. This requires careful studies to 
quantify its effect on the computed Bayes' factor. More details on this is available below.


The Method:
-----------
Bayesian model selection is a technique of comparing two competing models by taking ratios of 
marginalized likelihoods, also known as the Bayes' factor. In this case the competing models 
are the various tabulated neutron star equations of state. The likelihood function marginalized 
over the part of the parameter space consistent with any given model, :math:`{\rm M_i}` is

.. math::
   :nowrap:

   \begin{eqnarray}
      \mathcal{L}_{marg}({\rm M_i}) = \int p\left(x|\mathcal{M}, q, \tilde{\Lambda}(\mathcal{M}, q | {\rm M_i}), \delta\tilde{\Lambda}(\mathcal{M}, q| {\rm M_i})\right) \,p(\mathcal{M}, q)\,d\mathcal{M}\,dq\,,
   \end{eqnarray}

where, :math:`x` is the gravitational-wave data, :math:`\mathcal{M}` is the chirp-mass defined as 
:math:`\frac{(m_1 m_2)^{3/5}}{(m_1 + m_2)^{1/5}}`, :math:`q` is the mass ratio :math:`(m_2/m_1)`,
and :math:`(m_1, m_2)` being the mass of the heavier and the lighter object respectively.

In the above expression we are assuming that we have marginalized over all the other parameters
on which the signal morphology of the gravitational-wave depend. The equation of state model 
being deterministic, fixes the values of :math:`\tilde{\Lambda}` and 
:math:`\delta\tilde{\Lambda}` for a choice of :math:`(\mathcal{M}, q)` thus, restricting the 
possible values of the former, which we write as follows

.. math::
   :nowrap:

   \begin{eqnarray}
	\mathcal{L}_{marg}({\rm M_i})  = \int p\left(x|\mathcal{M}, q, \tilde{\Lambda}, \delta\tilde{\Lambda}\right)\delta\left(\tilde{\Lambda} - \tilde{\Lambda}_{\rm M_i}(\mathcal{M}, q)\right) \delta\left(\delta\tilde{\Lambda} - \delta\tilde{\Lambda}_{\rm eos}(\mathcal{M}, q)\right) p(\mathcal{M}, q) d\mathcal{M} dq d\tilde{\Lambda}d\delta\tilde{\Lambda}\,.
   \end{eqnarray}
  
Of course the first term in the integral is immediately identified to be the likelihood function 
in :math:`\left(\mathcal{M}, q, \tilde{\Lambda}, \delta\tilde{\Lambda}\right)`.

.. math::
   :nowrap:

   \begin{equation}
	\mathcal{L}_{marg}({\rm M_i})  = \int \mathcal{L}\left(\mathcal{M}, q, \tilde{\Lambda}, \delta\tilde{\Lambda}\right)\,\delta\left(\tilde{\Lambda} - \tilde{\Lambda}_{\rm M_i}(\mathcal{M}, q)\right)\, \delta\left(\delta\tilde{\Lambda} - \delta\tilde{\Lambda}_{\rm M_i}(\mathcal{M}, q)\right) \,p(\mathcal{M}, q)\,d\mathcal{M}\,dq\,d\tilde{\Lambda}\,d\delta\tilde{\Lambda}\,.
   \end{equation}

Thus, all we need now is to compute this likelihood function. This is where the results from the
parameter estimation is used. We conduct a Bayesian MCMC analysis of the gravitational wave data 
using uniform priors on :math:`\tilde{\Lambda}` and :math:`\delta\tilde{\Lambda}`.The likelihood 
function can be extracted from this analysis by simply invoking Bayes' theorem. 

.. math::
   :nowrap:

   \begin{equation}
	\mathcal{L}\left(\mathcal{M}, q, \tilde{\Lambda}, \delta\tilde{\Lambda}\right) = \frac{p\left(\mathcal{M}, q, \tilde{\Lambda}, \delta\tilde{\Lambda}\, | \,x\right)}{p(\mathcal{M}, q)\, p(\tilde{\Lambda}, \delta\tilde{\Lambda})}\,.
   \end{equation}

The choice of the uniform prior on :math:`\tilde{\Lambda}` and :math:`\delta\tilde{\Lambda}`
simplifies the aboves expression of likelihood such that the marginalized likelihood becomes

.. math::
   :nowrap:

   \begin{equation}
	\mathcal{L}_{marg}({\rm M_i})  = \int p\left(\mathcal{M}, q, \tilde{\Lambda}, \delta\tilde{\Lambda}\, | \,x\right)\,\delta\left(\tilde{\Lambda} - \tilde{\Lambda}_{\rm eos}(\mathcal{M}, q\right)\, \delta\left(\delta\tilde{\Lambda} - \delta\tilde{\Lambda}_{\rm eos}(\mathcal{M}, q\right) \,d\mathcal{M}\,dq\,d\tilde{\Lambda}\,d\delta\tilde{\Lambda}\,.
   \end{equation}

But this choice has implications. Independently varying the values of 
:math:`\tilde{\Lambda}` and :math:`\delta\tilde{\Lambda}` may lead to negative values in 
:math:`\Lambda_1` and :math:`\Lambda_2`, the tidal-deformability of the individual 
objects. Of course this is unphysical and one has to take careful steps to make sure that
this does not affect our results. We are saved by the fact that this posterior samples
are truly an intermediate step and final results is only obtained after integrating along
the equation of state curves, which presumably will pass through points that are never
going to be unphysical. Thus even though we are computing the posterior samples for all
the points which includes the unphysical points in the parameter space, we are however
only computing the integral at points that constitute physical waveform. 


The Approximation:
------------------

We now make two approximations motivated by physical observations. We note that the 
likelihood function is extremely sharply peaked in :math:`\mathcal{M}`, and hence can be 
approximated with a delta-function. We also noted that the likelihood is actually 
independent of :math:`\delta\tilde{\Lambda}`. Under the choice of the aforementioned 
uniform prior, this can be written as 

.. math::
   :nowrap:

   \begin{equation}
	p\left(\mathcal{M}, q, \tilde{\Lambda}, \delta\tilde{\Lambda}|x\right) \propto p\left(q, \tilde{\Lambda}|x\right) \delta(\mathcal{M} - \mathcal{M}_0)\,,
   \end{equation}

Performing the integrals in the evidence, one has

.. math::
   :nowrap:

   \begin{equation}
	\mathcal{L}_{marg}({\rm M_i}) \propto \int p\left(q, \tilde{\Lambda}_{\rm M_i}(q, \mathcal{M}_0)\, | \,x\right)\,dq\,,
   \end{equation}

where, :math:`\mathcal{M}_0` is the value of chirp-mass where the likelihood is 
peaked. This line integral is computed by approximating the likelihood by a kernel 
density approximation. Finally, the Bayes factor between two competing models is 
simply :math:`\mathcal{L}_{marg}({\rm M_i})/\mathcal{L}_{marg}({\rm M_j})`.


Unphysicality in approximation scheme:
--------------------------------------
The approximation scheme allows us to reduce the dimensionality of the problem from
four to two. However, in doing so, we can no longer claim that the equation of state
curves don't pass through the unphysical points. This is because we have essentially
suppressed the entire :math:`\delta\tilde{\Lambda}` coordinate. Within the validity
of this approximation, however we are rescued by the fact that the waveform generation
code simply uses the :math:`\tilde{\Lambda}` to compute the waveforms. This means, 
for all practical purpose we can fix :math:`\delta\tilde{\Lambda}` to some value that
keeps :math:`(\Lambda_1, \Lambda_2)` in the physical-space and let :math:`\tilde{\Lambda}`
vary uniformly, we will still get the same results because likelihood is practically
independent of :math:`\delta\tilde{\Lambda}`.

Nevertheless, one can modify the above method to restrict the computation in the 
physical-space. To do this let us define a step-function 
:math:`V(q, \tilde{\Lambda}, \delta\tilde{\Lambda})` which is zero when the combination
of :math:`(q, \tilde{\Lambda}, \delta\tilde{\Lambda})` gives negative values of 
:math:`(\Lambda_1` or :math:`\Lambda_2)`, and one elsewhere. Using this step-function
we can write

.. math::
   :nowrap:

   \begin{equation}
	p(q, \tilde{\Lambda} | x) = \frac{1}{N}\int p(q, \tilde{\Lambda}, \delta\tilde{\Lambda} | x)\, V(q, \tilde{\Lambda}, \delta\tilde{\Lambda})\,d\delta\tilde{\Lambda}\,,
   \end{equation}

where, 

.. math::
   :nowrap:

   \begin{equation}
        N(q, \tilde{\Lambda}) = \int V(q, \tilde{\Lambda}, \delta\tilde{\Lambda})\,d\delta\tilde{\Lambda}\,.
   \end{equation}

All we have to do now is to evaluate the weighted KDE to compute the above integral, 
where the weight for each 2D posterior sample is the quantity 
:math:`N(q, \tilde{\Lambda})`.




.. toctree::
   :maxdepth: 2
   :caption: Contents:

   code

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
