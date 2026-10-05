import json
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from compare_frozen_extents import run as extent_run
from compare_rain_thresholds import reference_equivalence
from compare_tracking_geometry import run
from summarize_tracking_geometry import (fixed_seed_metrics, event_seed_sizes,
                                         edge_change_seed_sizes, run as summarize)
from plot_tracking_geometry_case import select_case, run as plot_case
from test_frozen_extents import prepared_reference

pytest.importorskip('skimage', reason='Install envelope extra')


def fixture(tmp_path):
    rain = np.zeros((4, 7, 12), dtype=np.float32)
    rain[:, 2, 1:9] = .2
    rain[:, 2, 1], rain[:, 2, 8] = 2, 3
    rain[:, 5, 10] = .2  # No seed: no new object may enter.
    rain[:, 0, 0] = np.nan
    source, meta, reference = prepared_reference(tmp_path, rain)
    extents = tmp_path / 'extents'
    extent_run(source, meta, reference, extents, plots=False)
    return extents


def test_fixed_entry_exact_control_and_chunk_invariance(tmp_path):
    extents = fixture(tmp_path)
    first = run(extents, tmp_path / 'one', chunk_frames=1, workers=1)
    second = run(extents, tmp_path / 'three', chunk_frames=3, workers=2)
    assert all(first['baseline_control'].values())
    assert first['source_files_unchanged']
    for name in ('seed_control', 'extent_0p5', 'extent_0p1'):
        assert all(reference_equivalence(tmp_path / 'one' / name,
                                         tmp_path / 'three' / name).values())
    assert first['comparisons']['extent_0p1']['fixed_node_identities_equal']
    assert first['comparisons']['extent_0p1']['after']['nodes'] == 8
    core = json.loads((tmp_path / 'one/seed_control/chunk_000000/frame_hashes.json').read_text())
    expanded = json.loads((tmp_path / 'one/extent_0p1/chunk_000000/frame_hashes.json').read_text())
    assert core['identified_frame_sha256'] == expanded['identified_frame_sha256']
    assert core['tracking_geometry_frame_sha256'] != expanded['tracking_geometry_frame_sha256']
    assert (tmp_path / 'one/SUCCESS').exists()
    report = summarize([tmp_path / 'one'], tmp_path / 'report')
    assert report[0]['variants']['extent_0p1']['fixed_seed_metrics']['seed_rain_denominator'] == 15
    assert (tmp_path / 'report/SUCCESS').exists()
    with pytest.raises(ValueError, match='nonempty distinct'):
        summarize([tmp_path / 'one', tmp_path / 'one'], tmp_path / 'duplicate')
    with pytest.raises(FileExistsError):
        run(extents, tmp_path / 'one')


def test_corrupt_saved_seed_hash_and_deadline_fail_closed(tmp_path):
    extents = fixture(tmp_path)
    with pytest.raises(TimeoutError, match='Deadline reserve'):
        run(extents, tmp_path / 'late', finish_before=datetime.now(timezone.utc))
    assert not (tmp_path / 'late/SUCCESS').exists()
    hashes = extents / 'seed_frame_hashes.json'
    record = json.loads(hashes.read_text())
    record[0] = 'corrupted'
    hashes.write_text(json.dumps(record))
    with pytest.raises(AssertionError, match='seed hash'):
        run(extents, tmp_path / 'corrupt')
    assert not (tmp_path / 'corrupt/SUCCESS').exists()


def test_changed_input_and_invalid_execution_parameters_rejected(tmp_path):
    extents = fixture(tmp_path)
    for kwargs in ({'workers': 0}, {'chunk_frames': 0}):
        with pytest.raises(ValueError, match='Positive'):
            run(extents, tmp_path / 'invalid', **kwargs)
    assert not (tmp_path / 'invalid').exists()
    config = json.loads((extents / 'configuration.json').read_text())
    metadata = Path(config['metadata'])
    record = json.loads(metadata.read_text())
    record['timestamps'][0] = '2046-06-01T01:00:00'
    metadata.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='Prepared input changed'):
        run(extents, tmp_path / 'changed')
    assert not (tmp_path / 'changed').exists()


def test_association_weights_use_original_seed_rain_not_expanded_mass():
    seeds = [dict(node_id='1', time='0', local_label='1', precipitation_sum='1'),
             dict(node_id='2', time='1', local_label='1', precipitation_sum='2'),
             dict(node_id='3', time='1', local_label='2', precipitation_sum='8')]
    expanded = [dict(r, branch_id=str(i+1), precipitation_sum='1000') for i, r in enumerate(seeds)]
    stats = fixed_seed_metrics(seeds, expanded, [dict(child_node_id='2')], 2)
    assert stats['incoming_seed_rain_fraction'] == .2
    assert stats['incoming_seed_node_fraction'] == .5
    assert stats['seed_rain_denominator'] == 10
    expanded[-1]['local_label'] = '3'
    with pytest.raises(ValueError, match='Seed identity'):
        fixed_seed_metrics(seeds, expanded, [], 2)


def test_event_size_bins_and_explicit_case_selection():
    seeds = [dict(node_id='1', time='0', local_label='1', area_cells='1000', precipitation_sum='10'),
             dict(node_id='2', time='1', local_label='1', area_cells='1', precipitation_sum='1'),
             dict(node_id='3', time='1', local_label='2', area_cells='2000', precipitation_sum='20')]
    events = [dict(parent_node_ids='1', child_node_ids='2;3')]
    assert event_seed_sizes(seeds, events)['events_by_smallest_original_seed_cells']['1'] == 1
    changes = [dict(pair=[0, 1, 1, 1]), dict(pair=[0, 1, 1, 2])]
    assert edge_change_seed_sizes(seeds, changes) == {'1': 1, '2-15': 0, '16-999': 0, '>=1000': 1}
    assert select_case(seeds, changes)['pair'] == [0, 1, 1, 2]
    with pytest.raises(ValueError, match='No changes'):
        select_case(seeds, [])


def test_incomplete_plot_source_fails_without_output(tmp_path):
    with pytest.raises(ValueError, match='Successful experiment'):
        plot_case(tmp_path / 'missing', tmp_path / 'plot')
    assert not (tmp_path / 'plot').exists()
