import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from compare_frozen_extents import digest_frame
from compare_rain_thresholds import sha_file
from plot_frozen_extent_case import run
from run_npy_validation import disk
from step.identification import identify_frame


def test_full_grid_case_replay_and_seed_hash_guard(tmp_path):
    rain = np.zeros((8, 20), np.float32)
    rain[2, 2], rain[2, 3:9] = 2, .2
    seeds = identify_frame(rain, disk(0), min_size=1, threshold=1)
    source = tmp_path / 'rain.npy'
    np.save(source, rain[None])
    metadata = tmp_path / 'metadata.json'
    metadata.write_text(json.dumps({'timestamps': ['1996-07-15T01:00:00']}))
    extents = tmp_path / 'extents'
    extents.mkdir()
    (extents / 'SUCCESS').write_text('complete')
    (extents / 'configuration.json').write_text(json.dumps(dict(
        prepared_input=str(source), metadata=str(metadata), input_sha256=sha_file(source),
        metadata_sha256=sha_file(metadata), preset={'identification': dict(
            bridge_radius_cells=0, min_size_cells=1)})))
    (extents / 'seed_frame_hashes.json').write_text(json.dumps([digest_frame(seeds)]))
    case = dict(frame=0, core_label=1, threshold_mm_h=.1, core_cells=1,
                envelope_cells=7, core_node_id=1, branch_id=1)
    (extents / 'extreme_objects.json').write_text(json.dumps({'area_ratio': [case]}))
    result = run(extents, tmp_path / 'case')
    assert result['core_max_rain_mm_h'] == 2
    assert result['seed_labels_in_touched_weak_components'] == [1]
    assert (tmp_path / 'case' / 'case.png').exists()
    case_high = dict(case, threshold_mm_h=.5, envelope_cells=1)
    (extents / 'extreme_objects.json').write_text(json.dumps(
        {'area_ratio': [case], 'centroid_drift': [case_high]}))
    high_result = run(extents, tmp_path / 'high_case', 'centroid_drift')
    assert high_result['threshold_mm_h'] == .5
    assert high_result['extent_0p5_cells'] == 1
    assert high_result['extent_0p1_cells'] == 7
    assert high_result['extent_mean_rain_threshold_mm_h'] == .1
    with pytest.raises(FileExistsError):
        run(extents, tmp_path / 'case')
    (extents / 'seed_frame_hashes.json').write_text(json.dumps(['corrupt']))
    with pytest.raises(AssertionError, match='raster'):
        run(extents, tmp_path / 'bad')
    assert not (tmp_path / 'bad').exists()
