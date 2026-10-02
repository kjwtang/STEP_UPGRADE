import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from compare_frozen_extents import extent_frame, run, support_change
from compare_rain_thresholds import DEFAULT_PRESET, run as paired_run

pytest.importorskip('skimage', reason='Install envelope extra')


def prepared_reference(tmp_path, rain):
    source = tmp_path / 'rain.npy'
    np.save(source, rain)
    stamps = ['1996-06-01T%02d:00:00' % (i+1) for i in range(len(rain))]
    metadata = tmp_path / 'metadata.json'
    metadata.write_text(json.dumps(dict(timestamps=stamps,
        grid=dict(shape=list(rain.shape[1:]), DX=4000, DY=4000))))
    preset = json.loads(DEFAULT_PRESET.read_text())
    preset['identification']['bridge_radius_cells'] = 0
    preset_path = tmp_path / 'preset.json'
    preset_path.write_text(json.dumps(preset))
    paired_run(source, metadata, tmp_path / 'paired', preset_path,
               workers=1, chunk_frames=2, plots=False)
    return source, metadata, tmp_path / 'paired' / 'A_1mm'


def test_frozen_extent_dry_barriers_weak_bridge_and_chunk_invariance(tmp_path):
    rain = np.zeros((4, 7, 12), dtype=np.float32)
    rain[:, 2, 1:9] = .2
    rain[:, 2, 1], rain[:, 2, 8] = 2, 3
    rain[:, 5, 10] = .2  # Coreless rain must not be assigned.
    rain[:, 0, 0] = np.nan
    source, meta, reference = prepared_reference(tmp_path, rain)
    original = DEFAULT_PRESET.read_bytes()
    first = run(source, meta, reference, tmp_path / 'one', chunk_frames=1, plots=False)
    second = run(source, meta, reference, tmp_path / 'three', chunk_frames=3, plots=False)
    assert first['extents'] == second['extents']
    assert first['lower_threshold_support_changes'] == second['lower_threshold_support_changes']
    assert first['previous_0p5_extent_comparison']['geometry_equal']
    assert first['reference_structure']['nodes'] == 8
    assert first['extents']['0.1']['multi_seed_weak_components'] == 4
    assert first['extents']['0.1']['unassigned_cells'] == 4
    assert first['extents']['0.1']['assigned_cells'] == 32
    assert first['extents']['0.5']['assigned_cells'] == 8
    assert first['lower_threshold_support_changes']['added_cells'] == 24
    assert first['graph_source_files_unchanged'] and first['exact_seed_hashes_equal']
    assert (tmp_path / 'one' / 'SUCCESS').exists()
    assert DEFAULT_PRESET.read_bytes() == original
    with pytest.raises(FileExistsError):
        run(source, meta, reference, tmp_path / 'one', plots=False)
    changed = json.loads(meta.read_text())
    changed['timestamps'][0] = '2046-06-01T01:00:00'
    meta.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match='checksum'):
        run(source, meta, reference, tmp_path / 'mismatch', plots=False)
    assert not (tmp_path / 'mismatch').exists()


def test_dry_reference_and_corrupted_seed_hash_fail_closed(tmp_path):
    rain = np.zeros((2, 4, 6), dtype=np.float32)
    rain[:, 0, 0] = np.nan
    source, meta, reference = prepared_reference(tmp_path, rain)
    result = run(source, meta, reference, tmp_path / 'dry', plots=False)
    assert result['extents']['0.1']['assigned_rain_fraction'] is None
    assert result['extents']['0.1']['largest_seed_count_in_weak_component'] == 0
    assert 'NaN' not in (tmp_path / 'dry' / 'summary.json').read_text()
    hashes = reference / 'chunk_000000' / 'frame_hashes.json'
    record = json.loads(hashes.read_text())
    record['identified_frame_sha256'][0] = 'corrupted'
    hashes.write_text(json.dumps(record))
    with pytest.raises(AssertionError, match='seed raster'):
        run(source, meta, reference, tmp_path / 'bad', plots=False)
    assert not (tmp_path / 'bad' / 'SUCCESS').exists()


def test_support_change_distinguishes_reassignment_from_new_support():
    rain = np.array([[2., .6, .2, .8, 3.]])
    seeds = np.array([[1, 0, 0, 0, 2]])
    high = np.array([[1, 1, 0, 0, 2]])
    low = np.array([[1, 2, 2, 2, 2]])
    stats = support_change(rain, seeds, high, low, 0)
    assert stats['reassigned_existing_cells'] == 1
    assert stats['reassigned_existing_rain'] == pytest.approx(.6)
    assert stats['added_below_0p5_cells'] == 1
    assert stats['added_existing_ge_0p5_cells'] == 1
    assert stats['added_rain'] == pytest.approx(1.)
    low[0, 0] = 0
    with pytest.raises(AssertionError, match='lost assigned'):
        support_change(rain, seeds, high, low, 0)


def test_invalid_rate_input_never_marks_success(tmp_path):
    rain = np.zeros((2, 4, 6), dtype=np.float32)
    source, meta, reference = prepared_reference(tmp_path, rain)
    rain[0, 1, 1] = -1
    np.save(source, rain)
    # Re-signing a malformed input does not bypass explicit negative-rate checks.
    from compare_rain_thresholds import sha_file
    config_path = reference.parent / 'configuration.json'
    config = json.loads(config_path.read_text())
    config['input_sha256'] = sha_file(source)
    config_path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match='Negative hourly'):
        run(source, meta, reference, tmp_path / 'bad', plots=False)
    assert not (tmp_path / 'bad' / 'SUCCESS').exists()
