import csv
import json
from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_reference_stream import audit, inspect
from test_reference_stream import dataset, options
from run_reference_stream import run


def test_cross_chunk_audit_exact_equivalence(tmp_path):
    source = dataset(tmp_path / 'data')
    first, second = tmp_path / 'first', tmp_path / 'second'
    run(options(source, first, chunk_frames=2))
    run(options(source, second, chunk_frames=3))
    summary = audit(second, tmp_path / 'audit', reference=first)
    assert all(summary['reference_equivalence'].values())
    assert summary['structural_checks_pass']


def test_audit_detects_duplicate_nodes_and_missing_edge(tmp_path):
    source = dataset(tmp_path / 'data')
    out = tmp_path / 'run'
    run(options(source, out))
    path = out / 'chunk_000000' / 'objects.csv'
    with path.open() as handle:
        rows = list(csv.reader(handle))
    with path.open('w') as handle:
        csv.writer(handle).writerows(rows + rows[1:2])
    with pytest.raises(ValueError, match='Duplicate node'):
        inspect(out)


def test_identical_graphs_cannot_hide_different_source_dates(tmp_path):
    source = dataset(tmp_path / 'data')
    first, second = tmp_path / 'first', tmp_path / 'second'
    run(options(source, first))
    run(options(source, second))
    meta = second / 'chunk_000000' / 'metadata.json'
    value = json.loads(meta.read_text())
    value['timestamps'][0] = '2046-01-01T01:00:00'
    meta.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='timestamps differ'):
        audit(second, tmp_path / 'audit', reference=first)


def test_missing_saved_timestamps_rejected(tmp_path):
    source = dataset(tmp_path / 'data')
    out = tmp_path / 'run'
    run(options(source, out))
    meta = out / 'chunk_000000' / 'metadata.json'
    value = json.loads(meta.read_text())
    value['timestamps'] = value['timestamps'][:-1]
    meta.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='complete stream extent'):
        audit(out, tmp_path / 'audit')
