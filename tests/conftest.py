import pytest
from astropy.io.votable.tree import Param
from gwxtreme.eos_inference import ModelSelector, ParameterizedEoSSampler


@pytest.fixture
def model_selector_instance():
    ms = ModelSelector(
        posterior_file="./tests/GW170817_posterior_samples_IMRPhenomNRT_uniform_lambdas_prior.json",
        event_type="gw-3d",
        density_est_method="kde",
    )
    return ms


@pytest.fixture
def sampler_instance():
    sampler = ParameterizedEoSSampler(
        posterior_files=["./tests/GW170817_posterior_samples_IMRPhenomNRT_uniform_lambdas_prior.json"],
        event_types=["gw-3d"],
        eos_prior_bounds=[(0.2, 2.0), (-1.6, 1.7), (-0.6, 0.6), (-0.02, 0.02)],
        largest_observed_ns_mass=1.97,
        density_est_method="kde",
        parameterization="spectral",
    )
    return sampler
