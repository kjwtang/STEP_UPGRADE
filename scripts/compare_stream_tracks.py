#!/usr/bin/env python3
"""Compare saved stream graph edges by immutable frame/local-label identity."""
import argparse
from collections import Counter
import json
from pathlib import Path
import numpy as np
from replay_stream_pairs import read_rows, DEFAULT_PAIRS


def read_catalog(root):
    nodes, edges = {}, {}
    parts = sorted(root.glob('chunk_*'))
    if not parts:
        raise ValueError(f'No saved chunks in {root}')
    for part in parts:
        for row in read_rows(part/'objects.csv'):
            nid = int(row['node_id'])
            if nid in nodes:
                raise ValueError('Duplicate node ID')
            nodes[nid] = row
    identities = {nid:(int(r['time']),int(r['local_label'])) for nid,r in nodes.items()}
    if len(set(identities.values())) != len(identities):
        raise ValueError('Duplicate frame/local-label identity')
    for part in sorted(root.glob('chunk_*')):
        for row in read_rows(part/'edges.csv'):
            key = identities[int(row['parent_node_id'])]+identities[int(row['child_node_id'])]
            if key in edges:
                raise ValueError('Duplicate edge')
            edges[key] = row['event']
    measurements = {identities[nid]:{k:v for k,v in r.items()
                    if k not in ('node_id','branch_id','family_id')} for nid,r in nodes.items()}
    return identities, measurements, edges


def compare(before, after, output, reference=None):
    if output.exists():
        raise FileExistsError(output)
    ids, old_objects, old = read_catalog(before)
    _, new_objects, new = read_catalog(after)
    if old_objects != new_objects:
        raise ValueError('Object measurements differ: not a tracking-only comparison')
    old_parts, new_parts = sorted(before.glob('chunk_*')), sorted(after.glob('chunk_*'))
    if [p.name for p in old_parts] != [p.name for p in new_parts]:
        raise ValueError('This comparison requires the same chunk boundaries')
    for a,b in zip(old_parts,new_parts):
        if not np.array_equal(np.load(a/'identified_labels.npy',mmap_mode='r'),
                              np.load(b/'identified_labels.npy',mmap_mode='r')):
            raise ValueError('Identification masks differ')
    added = [dict(pair=k,event=new[k]) for k in sorted(new.keys()-old.keys())]
    removed = [dict(pair=k,event=old[k]) for k in sorted(old.keys()-new.keys())]
    changed = [dict(pair=k,before=old[k],after=new[k]) for k in sorted(old.keys() & new.keys())
               if old[k]!=new[k]]
    focus=[]
    for p,c in DEFAULT_PAIRS:
        if p in ids and c in ids:
            key=ids[p]+ids[c]
            focus.append(dict(source_node_pair=[p,c],object_pair=key,
                              before=old.get(key),after=new.get(key)))
    summary=dict(identification_and_measurements_equal=True,
        before=str(before.resolve()),after=str(after.resolve()),
        before_edges=len(old),after_edges=len(new),added=len(added),removed=len(removed),
        changed_event=len(changed),before_edge_types=dict(Counter(old.values())),
        after_edge_types=dict(Counter(new.values())),focus_pairs=focus,
        note='More edges or fewer gaps is not proof of better tracking. Compare event maps and false links. Focus node IDs refer to the specified before run.')
    if reference is not None:
        _, ref_objects, ref_edges = read_catalog(reference)
        if ref_objects != old_objects:
            raise ValueError('Reference object measurements differ')
        ref_parts=sorted(reference.glob('chunk_*'))
        if [p.name for p in ref_parts] != [p.name for p in old_parts]:
            raise ValueError('Reference chunk boundaries differ')
        for a,b in zip(ref_parts,old_parts):
            if not np.array_equal(np.load(a/'identified_labels.npy',mmap_mode='r'),
                                  np.load(b/'identified_labels.npy',mmap_mode='r')):
                raise ValueError('Reference identification differs')
        ref_continue={k for k,v in ref_edges.items() if v=='continue'}
        previously_missing=ref_continue-old.keys()
        summary['reference_recovery']=dict(
            reference=str(reference.resolve()),previously_missing_continue=len(previously_missing),
            recovered_pairs=[dict(pair=k,event=new[k]) for k in sorted(previously_missing & new.keys())],
            still_missing_pairs=sorted(previously_missing-new.keys()),
            newly_missing_pairs=sorted((ref_continue & old.keys())-new.keys()),
            note='Recovered pair may now be an event edge; recovery is not proof of correctness.')
    output.mkdir(parents=True)
    (output/'summary.json').write_text(json.dumps(summary,indent=2))
    (output/'edge_changes.json').write_text(json.dumps(dict(added=added,removed=removed,changed_event=changed),indent=2))
    print(json.dumps(summary,indent=2))
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('before',type=Path)
    p.add_argument('after',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--reference',type=Path,help='Optional original baseline to quantify recovery of missing continue edges')
    a=p.parse_args()
    compare(a.before,a.after,a.output_dir,a.reference)
