from pathlib import Path

import pytest

from GWXtreme import eos_model_selection as ems

sample_filename = Path(__file__).parent / 'data/test_samples.dat'
sample_filename = str(sample_filename.absolute())


@pytest.mark.parametrize(
    'ns_value,eos_1,eos_2,evidence_ratio',
    [(3, 'MS1', 'SLY', 0.0),
     (5, 'SLY', 'SLY', 1.0),
     (9, 'SLY', 'SLY', 1.0),
     (None, 'SLY', 'SLY', 1.0)]
)
def test_model_selection(ns_value, eos_1, eos_2, evidence_ratio):
    modsel_1 = ems.Model_selection(
        posteriorFile=sample_filename,
        priorFile=None,
        Ns=ns_value
    )
    assert modsel_1.computeEvidenceRatio(
        eos_1, eos_2) == pytest.approx(evidence_ratio, abs=1e-3)
