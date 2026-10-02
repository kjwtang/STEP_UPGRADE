import csv
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pytest
xr = pytest.importorskip('xarray')
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import run_reference_stream as stream


def dataset(root, count=7):
    root.mkdir()
    cumulative = np.zeros((50, 60), dtype=np.float32)
    for i in range(count):
        rain = np.zeros_like(cumulative)
        if i not in (0, 3):
            rain[20:25, 10 + i:15 + i] = 2
        cumulative += rain
        attrs = {key: 1. for key in ('MAP_PROJ', 'CEN_LAT', 'CEN_LON',
            'TRUELAT1', 'TRUELAT2', 'STAND_LON', 'MOAD_CEN_LAT', 'POLE_LAT', 'POLE_LON')}
        attrs.update(DX=4000., DY=4000.)
        ds = xr.Dataset({'RAINNC': (('Time', 'south_north', 'west_east'),
            cumulative[None].copy(), {'units': 'mm'})}, attrs=attrs)
        ds.to_netcdf(root / f'cstm_d01_1996-01-01_0{i}:00:00.rain')
    return root


def options(source, output, **overrides):
    settings = dict(input=source, output_dir=output, preset=stream.DEFAULT_PRESET,
        start=0, hours=6, sequence_id='test_reference_v1', workers=1,
        chunk_frames=2, save_labels=True, resume=False, stop_after_chunks=None,
        no_progress=True)
    settings.update(overrides)
    return SimpleNamespace(**settings)


def rows(root, name):
    result = []
    for part in sorted(root.glob('chunk_*')):
        with (part / name).open() as handle:
            result.extend(csv.DictReader(handle))
    return result


def frame_hashes(root, key):
    return [value for part in sorted(root.glob('chunk_*'))
            for value in json.loads((part / 'COMMITTED.json').read_text())[key]]


def test_resume_different_workers_and_chunks_preserves_ids_and_graph(tmp_path):
    source = dataset(tmp_path / 'data')
    whole, resumed = tmp_path / 'whole', tmp_path / 'resumed'
    expected = stream.run(options(source, whole, chunk_frames=6))
    paused = stream.run(options(source, resumed, stop_after_chunks=1))
    assert paused == dict(status='paused', committed_frames=2)
    assert not (resumed / 'SUCCESS').exists()
    actual = stream.run(options(source, resumed, resume=True, workers=2, chunk_frames=3,
                               id_executor='persistent'))
    assert actual['status'] == 'completed'
    for key in ('nodes', 'edges', 'gap_edges', 'events', 'historical_branches'):
        assert expected[key] == actual[key]
    for name in ('edges.csv', 'events.csv'):
        assert rows(whole, name) == rows(resumed, name)
    # Per-chunk as-of family roots can differ after future merges; retain all
    # immutable fields and compare final canonical roots separately if needed.
    for root in (whole, resumed):
        assert len(rows(root, 'objects.csv')) == expected['nodes']
    assert [{k: v for k, v in r.items() if k != 'family_id'} for r in rows(whole, 'objects.csv')] == [
        {k: v for k, v in r.items() if k != 'family_id'} for r in rows(resumed, 'objects.csv')]
    for key in ('identified_frame_sha256', 'tracked_frame_sha256'):
        assert frame_hashes(whole, key) == frame_hashes(resumed, key)
    assert stream.run(options(source, resumed, resume=True))['frames'] == 6
    with pytest.raises(FileExistsError):
        stream.run(options(source, resumed))


@pytest.mark.parametrize('change', ['sequence', 'preset', 'input', 'output', 'code'])
def test_resume_rejects_contract_changes(tmp_path, monkeypatch, change):
    source = dataset(tmp_path / 'data')
    out = tmp_path / 'run'
    stream.run(options(source, out, stop_after_chunks=1))
    a = options(source, out, resume=True)
    if change == 'sequence':
        a.sequence_id = 'different_sequence'
    elif change == 'preset':
        value = json.loads(a.preset.read_text())
        value['identification']['threshold_mm_h'] = .1
        a.preset = tmp_path / 'changed.json'
        a.preset.write_text(json.dumps(value))
    elif change == 'input':
        path = sorted(source.glob('*.rain'))[-1]
        path.touch()
    elif change == 'output':
        a.save_labels = False
    else:
        real_digest = stream.digest
        monkeypatch.setattr(stream, 'digest', lambda path: 'changed' if path.name == 'tracking.py' else real_digest(path))
    with pytest.raises(ValueError, match='Resume contract differs'):
        stream.run(a)
    assert not (out / 'SUCCESS').exists()


def test_corrupt_published_output_rejected(tmp_path):
    source = dataset(tmp_path / 'data')
    out = tmp_path / 'run'
    stream.run(options(source, out, stop_after_chunks=1))
    (out / 'chunk_000000' / 'edges.csv').write_text('corrupted')
    with pytest.raises(ValueError, match='Published output changed'):
        stream.run(options(source, out, resume=True))


def test_interrupted_publication_is_replayed_not_skipped(tmp_path, monkeypatch):
    source = dataset(tmp_path / 'data')
    out, control = tmp_path / 'run', tmp_path / 'control'
    expected = stream.run(options(source, control))
    stream.run(options(source, out, stop_after_chunks=1))
    real_replace = stream.os.replace
    def interruption(src, dest):
        if Path(dest).name == 'chunk_000002':
            raise RuntimeError('simulated kill before chunk publication')
        return real_replace(src, dest)
    with monkeypatch.context() as patch:
        patch.setattr(stream.os, 'replace', interruption)
        with pytest.raises(RuntimeError, match='simulated kill'):
            stream.run(options(source, out, resume=True))
    assert not (out / 'chunk_000002').exists()
    pending = list(out.glob('.*.pending'))
    assert pending and not (out / 'SUCCESS').exists()
    actual = stream.run(options(source, out, resume=True))
    assert actual['edges'] == expected['edges']
    assert frame_hashes(out, 'tracked_frame_sha256') == frame_hashes(control, 'tracked_frame_sha256')
    assert all(path.exists() for path in pending)  # no destructive cleanup


def test_malformed_chunks_do_not_resume(tmp_path):
    source = dataset(tmp_path / 'data')
    out = tmp_path / 'run'
    stream.run(options(source, out, stop_after_chunks=1))
    (out / 'chunk_000000').rename(out / 'chunk_000002')
    with pytest.raises(ValueError, match='not contiguous'):
        stream.run(options(source, out, resume=True))


def test_grid_preset_and_hourly_plan_checks(tmp_path):
    source = dataset(tmp_path / 'data')
    path = sorted(source.glob('*.rain'))[-1]
    path.rename(source / 'cstm_d01_1996-01-01_08:00:00.rain')
    with pytest.raises(ValueError, match='consecutive hourly'):
        stream.run(options(source, tmp_path / 'bad_time'))
    assert not (tmp_path / 'bad_time').exists()


def test_reuses_manifest_without_rescanning_each_chunk(tmp_path, monkeypatch):
    source = dataset(tmp_path / 'data')
    calls = []
    real_manifest = stream.manifest
    def once(*args, **kwargs):
        calls.append(1)
        return real_manifest(*args, **kwargs)
    monkeypatch.setattr(stream, 'manifest', once)
    stream.run(options(source, tmp_path / 'run'))
    assert calls == [1]


def test_single_writer_lock_releases_after_failure(tmp_path):
    out = tmp_path / 'run'
    out.mkdir()
    with stream.execution_lock(out):
        with pytest.raises(RuntimeError, match='Another process'):
            with stream.execution_lock(out):
                pass
    with pytest.raises(RuntimeError, match='injected failure'):
        with stream.execution_lock(out):
            raise RuntimeError('injected failure')
    with stream.execution_lock(out):
        pass


def test_resume_rejects_dependency_change(tmp_path, monkeypatch):
    source = dataset(tmp_path / 'data')
    out = tmp_path / 'run'
    stream.run(options(source, out, stop_after_chunks=1))
    monkeypatch.setattr(stream, 'dependency_versions', lambda: {'numpy': 'other'})
    with pytest.raises(ValueError, match='Resume contract differs'):
        stream.run(options(source, out, resume=True))


def test_corrupt_success_metadata_does_not_claim_completion(tmp_path):
    source = dataset(tmp_path / 'data')
    out = tmp_path / 'run'
    stream.run(options(source, out))
    (out / 'SUCCESS').write_text('{}')
    with pytest.raises(ValueError, match='SUCCESS/summary'):
        stream.run(options(source, out, resume=True))
