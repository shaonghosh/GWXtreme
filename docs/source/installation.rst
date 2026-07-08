Installation and Data
=====================

Installation
------------

GWXtreme can be installed with pip:

.. code::
    
    pip install gwxtreme

For an editable install directly from source:

.. code::
    
    git clone https://github.com/shaonghosh/GWXtreme.git
    cd GWXtreme
    pip install -e .

Compact Binary Coalescence (LIGO GW event) Parameter Estimation Data
--------------------------------------------------------------------

GWXtreme requires EOS-agnostic parameter estimation data to infer the EOS from gravitational wave events (i.e. BNS or NSBH detections).
This data can be obtained from the regular PE run results conducted for GW events, or from your own PE runs, if desired.
To run GWXtreme on the PE data utilized in `Ghosh et al. 2021 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.083003>`__, in which the TaylorF2 waveform was used for sampling, and which is compatible with the 2D variant
of the GWXtreme approximation as described in `Ghosh et al. 2021 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.104.083003>`__ and `Ray et al. 2023 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.107.043035>`__, you may download those original PE results from https://zenodo.org/records/4679013.

In general, PE data supplied to GWXtreme must include samples for the subset of EOS-sensitive parameters which the inference procedure relies on, which depends
on the variant of the approximation scheme used. It is required that uniform priors over the tidal deformability parameters were used in producing the PE samples (see the linked papers above for more details).
To summarize the requirements:

    * The ``gw-3d`` method (recommended) requires PE samples for :math:`(q, \Lambda_1, \Lambda_2)`, with sampling conducted using uniform priors over :math:`(\Lambda_1, \Lambda_2)`
    * The ``gw-2d`` method requires PE samples for :math:`(q, \tilde{\Lambda})`, with sampling conducted using uniform priors over :math:`(\tilde{\Lambda}, \delta \tilde{\Lambda})`
    * The ``gw-4d`` method requires PE samples for :math:`(m_1, m_2, \Lambda_1, \Lambda_2)`, with sampling conducted using uniform priors over :math:`(\Lambda_1, \Lambda_2)`

Pulsar Mass-Radius (NICER observation) Data
-----------------------------------------------------

To conduct EOS inference from simultaneous mass-radius measurements of pulsars, GWXtreme requires data files containing mass-compactness samples.