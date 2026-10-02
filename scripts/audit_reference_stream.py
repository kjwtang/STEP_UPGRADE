#!/usr/bin/env python3
"""Audit saved stream structure and exact execution equivalence across chunks.

Accepts diagnostic validate_real_data streams and run_reference_stream outputs.
Does not establish independent meteorological truth or raw-input correctness.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from step.tracking import load_tracking_state, _find


def table(root, name):
    rows = []
    for part in sorted(root.glob('chunk_*')):
        with (part / name).open(newline='') as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def frame_hashes(root, key):
    hashes = []
    filename = {'identified_frame_sha256': 'identified_labels.npy',
                'tracked_frame_sha256': 'tracked_labels.npy'}[key]
    for part in sorted(root.glob('chunk_*')):
        record = part / 'COMMITTED.json'
        if record.exists():
            hashes.extend(json.loads(record.read_text())[key])
        elif (part / filename).exists():
            array = np.load(part / filename, mmap_mode='r')
            hashes.extend(hashlib.sha256(memoryview(frame)).hexdigest() for frame in array)
        else:
            return None
    return hashes


def inspect(root):
    if not (root / 'SUCCESS').exists():
        raise ValueError(f'Incomplete stream: {root}')
    parts = sorted(root.glob('chunk_*'))
    if not parts:
        raise ValueError('No committed stream chunks')
    state = load_tracking_state(parts[-1] / 'state.json')
    objects = table(root, 'objects.csv')
    edges = table(root, 'edges.csv')
    events = table(root, 'events.csv')
    by_node = {int(row['node_id']): row for row in objects}
    if len(by_node) != len(objects):
        raise ValueError('Duplicate node ID')
    if len({(int(r['time']), int(r['local_label'])) for r in objects}) != len(objects):
        raise ValueError('Duplicate frame/local-label identity')
    if len({(int(r['time']), int(r['branch_id'])) for r in objects}) != len(objects):
        raise ValueError('Branch reused within one frame')
    if any(not 0 <= int(r['time']) <= state.last_time for r in objects):
        raise ValueError('Object time lies outside checkpoint extent')
    edge_keys = set()
    same_branch_edges = set()
    for row in edges:
        p, c = int(row['parent_node_id']), int(row['child_node_id'])
        if p not in by_node or c not in by_node or (p, c) in edge_keys:
            raise ValueError('Missing edge endpoint or duplicate edge')
        edge_keys.add((p, c))
        parent, child = by_node[p], by_node[c]
        delta = int(child['time']) - int(parent['time'])
        if delta != int(row['gap']) + 1 or delta <= 0 or int(row['time']) != int(child['time']):
            raise ValueError('Invalid forward edge time or gap')
        same_branch = parent['branch_id'] == child['branch_id']
        if same_branch != (row['event'] in ('continue', 'gap_continue')):
            raise ValueError('Branch identity does not match edge type')
        if row['event'] == 'gap_continue' and int(row['gap']) < 1:
            raise ValueError('Gap continuation has no gap')
        if same_branch:
            same_branch_edges.add((p, c))
    groups = defaultdict(list)
    for row in objects:
        groups[int(row['branch_id'])].append((int(row['time']), int(row['node_id'])))
    for branch, values in groups.items():
        values.sort()
        if any((a[1], b[1]) not in same_branch_edges for a, b in zip(values, values[1:])):
            raise ValueError(f'Disconnected/reused branch: {branch}')
        if branch not in state.family_parent:
            raise ValueError('Observed branch absent from checkpoint history')
    # IDs must be allocated monotonically on first appearance, not compacted
    # after old branches disappear. Multiple births within a frame are sorted.
    births = sorted((min(v)[0], branch) for branch, v in groups.items())
    if [branch for _, branch in births] != sorted(groups):
        raise ValueError('Branch IDs are not monotonically allocated')
    if state.next_branch_id <= max(groups, default=0) or state.next_node_id <= max(by_node, default=0):
        raise ValueError('Checkpoint next-ID counters would recycle an observed ID')
    event_ids = [int(r['event_id']) for r in events]
    if len(event_ids) != len(set(event_ids)) or state.next_event_id <= max(event_ids, default=0):
        raise ValueError('Event ID duplicated or next counter would recycle it')
    for row in events:
        parents = [int(v) for v in row['parent_node_ids'].split(';') if v]
        children = [int(v) for v in row['child_node_ids'].split(';') if v]
        if any(n not in by_node for n in parents + children):
            raise ValueError('Event has missing endpoint')
        # Complex components may contain only some Cartesian endpoint pairs.
        if any(not any((p, c) in edge_keys for c in children) for p in parents) or any(
                not any((p, c) in edge_keys for p in parents) for c in children):
            raise ValueError('Event endpoints are not represented in graph edges')
    roots = {b: _find(state.family_parent, b) for b in sorted(state.family_parent)}
    canonical_objects = [{**row, 'family_id': str(roots[int(row['branch_id'])])} for row in objects]
    summary = dict(frames=state.last_time + 1, nodes=len(objects), edges=len(edges),
        branches=len(groups), families=len(set(roots.values())),
        events=dict(Counter(r['event'] for r in events)),
        edge_types=dict(Counter(r['event'] for r in edges)),
        structural_checks_pass=True, final_checkpoint_bytes=(parts[-1] / 'state.json').stat().st_size,
        total_checkpoint_bytes=sum((p / 'state.json').stat().st_size for p in parts),
        sequence_id=state.sequence_id,
        note='Graph/catalog structure only; not meteorological accuracy or independent ingestion validation.')
    return summary, canonical_objects, edges, events, roots


def audit(source, output, reference=None, prefix_reference=None):
    if output.exists():
        raise FileExistsError(output)
    summary, objects, edges, events, roots = inspect(source)
    if reference is not None:
        other, o, e, v, r = inspect(reference)
        checks = dict(canonical_objects_equal=objects == o, full_edges_and_scores_equal=edges == e,
            events_equal=events == v, final_family_roots_equal=roots == r)
        for key in ('identified_frame_sha256', 'tracked_frame_sha256'):
            left, right = frame_hashes(source, key), frame_hashes(reference, key)
            checks[key + '_equal'] = None if left is None or right is None else left == right
        if any(value is False for value in checks.values()):
            raise AssertionError(f'Execution equivalence failed: {checks}')
        summary['reference_equivalence'] = checks
        summary['reference_structure'] = other
    if prefix_reference is not None:
        other, *_ = inspect(prefix_reference)
        hashes = {}
        for key in ('identified_frame_sha256', 'tracked_frame_sha256'):
            left, right = frame_hashes(source, key), frame_hashes(prefix_reference, key)
            if left is None or right is None:
                raise ValueError('Prefix comparison requires hashes or saved rasters')
            hashes[key + '_equal'] = left[:len(right)] == right
            if not hashes[key + '_equal']:
                raise AssertionError('Prefix raster hashes differ')
        summary['prefix_raster_equivalence'] = dict(frames=other['frames'], **hashes)
    output.mkdir(parents=True)
    (output / 'summary.json').write_text(json.dumps(summary, indent=2))
    (output / 'REPORT.md').write_text('# Stream integrity and execution equivalence\n\n'
        'Structural audit passed. Exact comparisons are listed below; null means unavailable, not passed. '
        'Family roots are normalized to the final checkpoint before catalog comparison. '
        'Input correctness and meteorological truth are not independently certified.\n\n'
        '```json\n' + json.dumps(summary, indent=2) + '\n```\n')
    (output / 'SUCCESS').write_text('Requested structure and exact-equivalence audits passed\n')
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--reference', type=Path)
    p.add_argument('--prefix-reference', type=Path,
                   help='Earlier shorter stream starting at the identical snapshot')
    a = p.parse_args()
    print(json.dumps(audit(a.source, a.output_dir, a.reference, a.prefix_reference), indent=2))
