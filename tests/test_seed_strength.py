import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_seed_strength import (THRESHOLDS, branch_duration_context, longest_run,
                                 policies, run, strong_patches, verify_manifest_coverage)
from compare_frozen_extents import run as extent_run
from test_frozen_extents import prepared_reference


def test_raw_connected_patch_not_scattered_or_dilated_area():
    rain = np.zeros((3, 10), np.float32)
    seeds = np.zeros(rain.shape, np.int32)
    rain[0, :3], seeds[0, :3] = [1, 2, 5], 2
    rain[0, 8], seeds[0, 8] = 5, 2  # Same STEP object, disconnected raw patch.
    rain[2, 6:9], seeds[2, 6:9] = 2, 7
    original = seeds.copy()
    count, largest, patches = strong_patches(rain, seeds, 2)
    assert (count[2], largest[2], patches[2]) == (3, 2, 2)
    assert (count[7], largest[7], patches[7]) == (3, 3, 1)
    count, largest, patches = strong_patches(rain, seeds, 5)
    assert (count[2], largest[2], patches[2]) == (2, 1, 2)
    assert largest[7] == 0
    assert np.array_equal(seeds, original)
    for invalid in (.5, np.nan):
        with pytest.raises(ValueError):
            strong_patches(rain, seeds, invalid)
    with pytest.raises(AssertionError, match='crosses'):
        strong_patches(np.ones((1, 2)), np.array([[2, 7]]), 1)


def test_consecutive_hours_break_on_gaps_and_weak_snapshots():
    assert longest_run([(3, True), (2, True), (0, True)]) == 2
    assert longest_run([(0, True), (1, False), (2, True)]) == 1
    assert longest_run([(0, True), (1, True), (2, True)]) == 3
    assert longest_run([]) == 0


def test_legacy_manifest_never_ignores_changed_or_new_graph_files(tmp_path):
    reference = tmp_path / 'paired' / 'A_1mm'
    graph = str((reference / 'chunk_000000' / 'edges.csv').resolve())
    old = {graph: 'same'}
    optional = str((reference.parent / 'c_envelope_objects.csv').resolve())
    assert verify_manifest_coverage(old, {**old, optional: 'new-signature'}, reference) == [optional]
    for current in ({graph: 'changed'}, {}, {**old, '/unexpected/chunk/objects.csv': 'new'}):
        with pytest.raises(ValueError, match='source changed'):
            verify_manifest_coverage(old, current, reference)


def test_branch_phase_is_retrospective_not_every_snapshot_or_family():
    records = []
    for threshold in THRESHOLDS:
        for node, branch, frame, area in [(1, 10, 0, 0), (2, 10, 1, 4), (3, 10, 2, 0),
                                         (4, 20, 0, 4), (5, 20, 2, 4), (6, 20, 3, 4)]:
            records.append(dict(node_id=node, branch_id=branch, frame=frame,
                threshold_mm_h=threshold, largest_patch_cells=area))
    result = policies(records, {t: {n: 1. for n in range(1, 7)} for t in (.5, .1)},
                      10., {'0.5': 6., '0.1': 6.})
    select = lambda p, d: next(r for r in result if r['threshold_mm_h'] == 2 and
        r['min_patch_cells'] == 4 and r['selection'] == p and r['min_consecutive_hours'] == d)
    assert select('snapshot', None)['retained_snapshots'] == 4
    assert select('branch_phase', 1)['retained_snapshots'] == 6
    assert select('branch_phase', 3)['retained_snapshots'] == 0


def test_duration_context_flags_overlap_without_causal_attribution():
    objects = [dict(node_id=n, branch_id=b, time=t) for n, b, t in (
        (1, 10, 0), (2, 10, 1), (3, 20, 2), (4, 20, 4),
        (5, 30, 1), (6, 30, 2), (7, 30, 3))]
    edges = [dict(parent_node_id=2, child_node_id=3, event='split'),
             dict(parent_node_id=3, child_node_id=4, event='gap_continue')]
    result = branch_duration_context(objects, edges, {n: 1 for n in range(1, 8)}, 5)
    groups = result['subsets']
    assert groups['all_short_consecutive_branches']['branches'] == 2
    assert groups['touching_window_boundary']['assigned_0p1_rain_fraction'] == pytest.approx(4/7)
    assert groups['touching_event_endpoint']['assigned_0p1_rain_fraction'] == pytest.approx(4/7)
    assert groups['containing_gap']['assigned_0p1_rain_fraction'] == pytest.approx(2/7)
    assert groups['without_these_flags']['branches'] == 0


def test_integrity_control_and_changed_source_fail_closed(tmp_path):
    rain = np.zeros((4, 6, 8), np.float32)
    rain[:, 2, 1:5] = 3
    rain[:, 2, 5:7] = .2
    source, metadata, reference = prepared_reference(tmp_path, rain)
    extent_run(source, metadata, reference, tmp_path / 'extents', plots=False)
    result = run(tmp_path / 'extents', tmp_path / 'strength')
    assert result['exact_seed_hashes_equal'] and result['graph_and_extent_sources_unchanged']
    assert result['baseline_preserves_all_assigned_rain']
    selected = next(r for r in result['policies'] if r['threshold_mm_h'] == 2 and
        r['min_patch_cells'] == 4 and r['selection'] == 'branch_phase' and r['min_consecutive_hours'] == 3)
    assert selected['retained_snapshots'] == 4
    assert selected['retained_fraction_of_assigned_rain']['0.1'] == 1
    assert (tmp_path / 'strength' / 'SUCCESS').is_file()
    with pytest.raises(TimeoutError, match='reserve'):
        run(tmp_path / 'extents', tmp_path / 'cutoff',
            datetime.now(timezone.utc)+timedelta(minutes=5))
    assert (tmp_path / 'cutoff' / 'PLAN.json').is_file()
    assert not (tmp_path / 'cutoff' / 'SUCCESS').exists()
    with pytest.raises(FileExistsError):
        run(tmp_path / 'extents', tmp_path / 'strength')
    config = tmp_path / 'extents' / 'seed_frame_hashes.json'
    hashes = json.loads(config.read_text())
    hashes[0] = 'changed'
    config.write_text(json.dumps(hashes))
    with pytest.raises(AssertionError, match='seed hash'):
        run(tmp_path / 'extents', tmp_path / 'bad')
    assert not (tmp_path / 'bad' / 'SUCCESS').exists()
    summary_path = tmp_path / 'extents' / 'summary.json'
    summary = json.loads(summary_path.read_text())
    summary['extents']['0.5']['assigned_rain'] = 0
    summary_path.write_text(json.dumps(summary))
    with pytest.raises(ValueError, match='positive assigned'):
        run(tmp_path / 'extents', tmp_path / 'empty')
    assert not (tmp_path / 'empty').exists()
