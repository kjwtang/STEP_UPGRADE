from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from test_reference_stream import dataset, options
from run_reference_stream import run
from summarize_reference_catalog import summarize


def test_area_mass_fractions_conserve_labeled_catalog(tmp_path):
    pytest.importorskip('matplotlib')
    source = dataset(tmp_path / 'data')
    run(options(source, tmp_path / 'run'))
    report = summarize(tmp_path / 'run', tmp_path / 'screening')
    assert sum(r['object_fraction'] for r in report['area_groups']) == 1
    assert sum(r['labeled_rain_fraction'] for r in report['area_groups']) == 1
    assert report['area_groups'][1]['object_fraction'] == 1
    assert report['area_groups'][1]['labeled_rain_mm_cells'] == 250
    assert report['branch_span_hours_max'] == 6
    assert report['one_snapshot_branch_fraction'] == 0
