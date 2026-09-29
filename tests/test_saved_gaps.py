import csv
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_saved_gaps import audit, footprint_candidates


def write_rows(path, data):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def test_gap_across_chunks_with_intermediate_object(tmp_path):
    source = tmp_path / 'source'
    for time in range(3):
        part = source / f'chunk_{time:06d}'
        part.mkdir(parents=True)
        labels = np.zeros((1, 5, 5), dtype=np.int32)
        labels[0, 2, 2] = 1
        np.save(part / 'identified_labels.npy', labels)
        write_rows(part / 'objects.csv', [dict(node_id=time+1, time=time, local_label=1, area_cells=1)])
        # Header-only files for the first two chunks.
        edge = dict(parent_node_id=1, child_node_id=3, time=2, gap=1, event='gap_continue')
        with (part / 'edges.csv').open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(edge))
            writer.writeheader()
            if time == 2:
                writer.writerow(edge)
    result = audit(source, tmp_path / 'audit', plot=True)
    assert result[0]['intermediate'][0]['candidates'][0]['node_id'] == 2
    assert len(list((tmp_path / 'audit').glob('*.png'))) == 1
    with pytest.raises(FileExistsError):
        audit(source, tmp_path / 'audit', plot=False)


def test_no_overlap_is_not_a_candidate():
    labels = np.array([[0, 1], [0, 0]])
    mask = np.array([[True, False], [False, False]])
    assert footprint_candidates(labels, mask, {}) == []
