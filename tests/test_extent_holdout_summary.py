import csv
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from summarize_extent_holdouts import CHECKS, read_run, run


def test_completion_guard_and_seed_area_bins(tmp_path):
    with pytest.raises(ValueError, match='Incomplete'):
        read_run(tmp_path)
    (tmp_path / 'SUCCESS').write_text('completed')
    summary = {key: True for key in CHECKS}
    summary['extents'] = {'0.1': {'assigned_rain': 10}}
    (tmp_path / 'summary.json').write_text(json.dumps(summary))
    (tmp_path / 'configuration.json').write_text(json.dumps({'preset': {'id': 'unchanged'}}))
    rows = [dict(threshold_mm_h=.1, core_cells=1, envelope_rain=2,
                 envelope_to_core_cell_ratio=100, centroid_drift_cells=10),
            dict(threshold_mm_h=.1, core_cells=1000, envelope_rain=8,
                 envelope_to_core_cell_ratio=2, centroid_drift_cells=1)]
    with (tmp_path / 'extent_objects.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = read_run(tmp_path)
    tiny, _, _, large = result['seed_area_bins']
    assert tiny['assigned_rain_fraction'] == .2
    assert large['assigned_rain_fraction'] == .8
    assert tiny['area_ratio_max'] == 100
    summary['extents']['0.1']['assigned_rain'] = 11
    (tmp_path / 'summary.json').write_text(json.dumps(summary))
    with pytest.raises(AssertionError, match='sum'):
        read_run(tmp_path)
    summary['exact_seed_hashes_equal'] = False
    (tmp_path / 'summary.json').write_text(json.dumps(summary))
    with pytest.raises(ValueError, match='Integrity'):
        read_run(tmp_path)


def test_no_inputs_or_output_overwrite(tmp_path):
    with pytest.raises(ValueError, match='At least one'):
        run([], tmp_path / 'empty')
    with pytest.raises(FileExistsError):
        run([], tmp_path)
