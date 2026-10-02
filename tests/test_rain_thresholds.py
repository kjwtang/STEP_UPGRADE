import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from compare_rain_thresholds import (DEFAULT_PRESET, accounting, core_memberships,
                                    projected_edges, reference_equivalence, run)
from step.envelope import identify_envelope
from audit_reference_stream import inspect

pytest.importorskip('skimage', reason='Install envelope extra for partition sensitivity')


def test_memberships_do_not_assume_equal_ids_and_keep_merged_seeds():
    core = np.array([[1, 1, 0, 2, 2], [0, 0, 3, 0, 0]])
    lower = np.array([[4, 4, 4, 4, 4], [0, 0, 7, 0, 0]])
    assert core_memberships(core, lower, 5) == [
        dict(frame=5, core_label=1, lower_label=4, shared_core_pixels=2),
        dict(frame=5, core_label=2, lower_label=4, shared_core_pixels=2),
        dict(frame=5, core_label=3, lower_label=7, shared_core_pixels=1)]
    lower[0, 0] = 0
    with pytest.raises(AssertionError, match='lost'):
        core_memberships(core, lower, 5)


def test_partition_mass_ownership_and_coreless_weak_rain():
    rain = np.zeros((7, 12), dtype=np.float32)
    rain[2, 1:9] = .5
    rain[2, 1] = 2
    rain[2, 8] = 3
    rain[5, 10] = .7
    rain[0, 0] = np.nan
    core, partition = identify_envelope(rain, 1, .5, mode='partition')
    _, system = identify_envelope(rain, 1, .5, mode='system')
    lower = (np.isfinite(rain) & (rain >= .5)).astype(np.int32)
    row = accounting(rain, core, lower, partition, system)
    assert row['missing_cells'] == 1
    assert row['unassigned_lower_cells'] == 1
    assert row['unassigned_lower_rain'] == pytest.approx(.7)
    assert row['c_rain'] == pytest.approx(8)
    assert partition[2, 1] != partition[2, 8]
    assert system[2, 1] == system[2, 8]


def test_projected_edges_exclude_multicore_and_weak_only_endpoints():
    a = [dict(node_id=str(n), time=str(t), local_label=str(l))
         for n, t, l in ((1, 0, 1), (2, 0, 2), (3, 1, 1), (4, 1, 2))]
    b = [dict(node_id=str(n), time=str(t), local_label=str(l))
         for n, t, l in ((11, 0, 8), (12, 0, 9), (13, 1, 8), (14, 1, 9))]
    ae = [dict(parent_node_id='1', child_node_id='3', event='continue'),
          dict(parent_node_id='2', child_node_id='4', event='continue')]
    be = [dict(parent_node_id='11', child_node_id='13', event='continue'),
          dict(parent_node_id='12', child_node_id='14', event='continue')]
    mapping = [dict(frame=t, core_label=l, lower_label=8, shared_core_pixels=1)
               for t in (0, 1) for l in (1, 2)]
    summary, changes = projected_edges(a, ae, b, be, mapping)
    assert summary['bijectively_mapped_core_snapshots'] == 0
    assert summary['comparable_a_edges'] == summary['comparable_b_edges'] == 0
    assert summary['excluded_a_edges'] == summary['excluded_b_edges'] == 2
    assert not changes['added'] and not changes['removed']
    mapping = [dict(frame=t, core_label=l, lower_label=l+7, shared_core_pixels=1)
               for t in (0, 1) for l in (1, 2)]
    summary, _ = projected_edges(a, ae, b, be, mapping)
    assert summary['retained'] == 2 and summary['bijectively_mapped_core_snapshots'] == 4


def test_synthetic_paired_run_chunk_equality_and_no_default_mutation(tmp_path):
    rain = np.zeros((4, 8, 32), dtype=np.float32)
    rain[:, 3, 2:23] = .5
    rain[:, 3, 2] = 2
    rain[:, 3, 22] = 3
    rain[:, 6, 29] = .7
    np.save(tmp_path / 'rain.npy', rain)
    stamps = ['1996-06-01T0%d:00:00' % h for h in (1, 2, 3, 4)]
    meta = dict(timestamps=stamps, grid=dict(shape=[8, 32], DX=4000, DY=4000))
    (tmp_path / 'metadata.json').write_text(json.dumps(meta))
    original = DEFAULT_PRESET.read_bytes()
    preset = json.loads(original)
    preset['identification']['bridge_radius_cells'] = 0
    (tmp_path / 'preset.json').write_text(json.dumps(preset))
    args = (tmp_path / 'rain.npy', tmp_path / 'metadata.json')
    first = run(*args, tmp_path / 'one', tmp_path / 'preset.json', chunk_frames=1, plots=False)
    second = run(*args, tmp_path / 'three', tmp_path / 'preset.json', chunk_frames=3, plots=False)
    assert first['accounting'] == second['accounting']
    assert first['accounting']['b_objects_with_multiple_a_seeds'] == 4
    assert first['accounting']['b_objects_without_a_seed'] == 4
    assert first['accounting']['c_rain'] < first['accounting']['b_rain']
    for variant in ('A_1mm', 'B_0p5mm'):
        _, o, e, v, r = inspect(tmp_path / 'one' / variant)
        _, oo, ee, vv, rr = inspect(tmp_path / 'three' / variant)
        assert (o, e, v, r) == (oo, ee, vv, rr)
    assert DEFAULT_PRESET.read_bytes() == original
    assert (tmp_path / 'one' / 'SUCCESS').exists()
    ref = tmp_path / 'three' / 'A_1mm'
    assert all(reference_equivalence(tmp_path / 'one' / 'A_1mm', ref).values())
    stamp_path = ref / 'chunk_000000' / 'metadata.json'
    saved = json.loads(stamp_path.read_text())
    saved['timestamps'][0] = '2046-06-01T01:00:00'
    stamp_path.write_text(json.dumps(saved))
    with pytest.raises(AssertionError, match='saved_timestamps_equal'):
        reference_equivalence(tmp_path / 'one' / 'A_1mm', ref)
    with pytest.raises(FileExistsError):
        run(*args, tmp_path / 'one', tmp_path / 'preset.json', plots=False)
    meta['timestamps'][1] = meta['timestamps'][0]
    (tmp_path / 'metadata.json').write_text(json.dumps(meta))
    with pytest.raises(ValueError, match='hourly timestamps'):
        run(*args, tmp_path / 'bad', tmp_path / 'preset.json', plots=False)
    assert not (tmp_path / 'bad').exists()


def test_dry_window_has_no_nan_summary_and_missing_cells_remain_barriers(tmp_path):
    rain = np.zeros((2, 4, 6), dtype=np.float32)
    rain[:, 0, 0] = np.nan
    np.save(tmp_path / 'rain.npy', rain)
    meta = dict(timestamps=['1996-06-01T01:00:00', '1996-06-01T02:00:00'],
                grid=dict(shape=[4, 6], DX=4000, DY=4000))
    (tmp_path / 'metadata.json').write_text(json.dumps(meta))
    result = run(tmp_path / 'rain.npy', tmp_path / 'metadata.json', tmp_path / 'dry',
                 workers=1, chunk_frames=1, plots=False)
    assert result['a_structure']['nodes'] == result['b_structure']['nodes'] == 0
    assert result['accounting']['a_rain_fraction'] is None
    assert result['c_centroid_drift_cells_quantiles']['1.0'] is None
    assert 'NaN' not in (tmp_path / 'dry' / 'summary.json').read_text()
    rain[0, 1, 1] = -1
    np.save(tmp_path / 'negative.npy', rain)
    with pytest.raises(ValueError, match='Negative hourly rate'):
        run(tmp_path / 'negative.npy', tmp_path / 'metadata.json', tmp_path / 'negative',
            workers=1, plots=False)
    assert not (tmp_path / 'negative' / 'SUCCESS').exists()
