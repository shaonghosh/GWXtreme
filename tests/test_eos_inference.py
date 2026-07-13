import json

import numpy as np
import pytest
from gwxtreme.eos_inference import ModelSelector, ParameterizedEoSSampler, load_samples


def test_model_selector(model_selector_instance: ModelSelector, tmp_path):
    savepath = tmp_path / "model_selector_test_output.json"
    result = model_selector_instance.compute_eos_evidence_ratio("APR4_EPP", "SLY", n_resamplings=10, save_file=str(savepath))
    with open(savepath) as f:
        data = json.load(f)

    assert result is not None
    assert data is not None


def test_sampler(sampler_instance: ParameterizedEoSSampler, tmp_path):
    savepath = tmp_path / "sampler_test_output.hdf5"
    sampler_instance.run_sampler(nsteps=10, nwalkers=10, save_file=str(savepath))
    samples = load_samples(savepath, burn_in_frac=0, thin=1)

    assert samples is not None


@pytest.mark.parametrize(
    "eos_name,eos_mass_lambda_file,eos_mass_radius_k_file",
    [
        ("APR4_EPP", "./tests/APR4_EPP_mass_lambda.txt", None),
        ("APR4_EPP", None, "./tests/APR4_EPP_mass_radius_k.txt"),
    ],
)
def test_model_selection_from_files(model_selector_instance: ModelSelector, eos_name, eos_mass_lambda_file, eos_mass_radius_k_file):
    result = model_selector_instance.compute_eos_evidence_ratio(
        target_eos_name=eos_name,
        reference_eos_mass_lambda_file=eos_mass_lambda_file,
        reference_eos_mass_radius_k_file=eos_mass_radius_k_file,
    )
    assert result.item() == pytest.approx(1.0, abs=1e-5)
