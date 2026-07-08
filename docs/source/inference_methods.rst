Inference Methods
=================

GWXtreme supports the use of several variants of its core approximation scheme for gravitational wave BNS/NSBH events that correspond
to using different parameterizations of the PE data for these events. In any case, the parameterization used includes only a subset
of the PE parameters, corresponding to those that are most sensitive to the NS equation of state. All variants are approximations because
they ignore NS spins, focusing only on masses and (first-order) tidal parameters.

These approximation variants are:

* 2D method: uses :math:`(q, \tilde{\Lambda})` PE samples, as well as the mean chirp mass from the PE posterior
* 3D method: uses :math:`(q, \Lambda_1, \Lambda_2)` PE samples, as well as the mean chirp mass from the PE posterior
* 4D method: uses :math:`(m_1, m_2, \Lambda_1, \Lambda_2)` PE samples

**The recommended variant is the 3D method**, because standard uniform priors over :math:`(\Lambda_1, \Lambda_2)` can be used in obtaining the PE data,
which means any GW waveform model, including inspiral-merger-ringdown models, can also be used for the PE.

The 2D method is similar in performance to the 3D, however it requires that the PE data was obtained using non-standard uniform priors over
:math:`(\tilde{\Lambda}, \delta \tilde{\Lambda})`. This prior leads to unphysical, negative values for the individual tidal deformabilities, which are
not supported by most waveform models. TaylorF2 supports this method, which was used in `Ghosh et al. 2021 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.083003>`__.

The 4D method is also similar in performance to the 3D approach, and it is compatible with the same priors and waveforms. However, it is ~ 100 times slower
to compute with and is therefore not practical for multi-event or parameterized EOS inference. (See :doc:`the examples </examples/index>`.)