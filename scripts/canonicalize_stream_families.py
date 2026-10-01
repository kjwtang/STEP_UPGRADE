#!/usr/bin/env python3
"""Export final-known family roots without rewriting branch IDs or source files."""
import argparse
import csv
import json
from pathlib import Path

from step.tracking import load_tracking_state, _find


def export(source, output):
    parts = sorted(source.glob('chunk_*'))
    if not parts or not (source/'SUCCESS').exists():
        raise ValueError('Require a successfully completed saved stream')
    if output.exists():
        raise FileExistsError(output)
    state = load_tracking_state(parts[-1]/'state.json')
    roots = {branch: _find(state.family_parent, branch)
             for branch in sorted(state.family_parent)}
    output.mkdir(parents=True)
    count = 0
    fields = None
    with (output/'objects.csv').open('w', newline='') as handle:
        writer = None
        for part in parts:
            with (part/'objects.csv').open(newline='') as rows:
                reader = csv.DictReader(rows)
                if 'canonical_family_id' in (reader.fieldnames or []):
                    raise ValueError('Input is already canonicalized')
                if fields is None:
                    fields = reader.fieldnames
                    if not fields or not {'branch_id','family_id','node_id'} <= set(fields):
                        raise ValueError('Missing required object fields')
                    writer = csv.DictWriter(handle,fieldnames=fields+['canonical_family_id'])
                    writer.writeheader()
                elif fields != reader.fieldnames:
                    raise ValueError('Chunk catalog headers differ')
                for row in reader:
                    branch = int(row['branch_id'])
                    if branch not in roots:
                        raise ValueError('Branch absent from final checkpoint')
                    row['canonical_family_id'] = roots[branch]
                    writer.writerow(row)
                    count += 1
    with (output/'family_map.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['branch_id','canonical_family_id'])
        writer.writerows(roots.items())
    summary = dict(objects=count,branches=len(roots),families=len(set(roots.values())),
                   sequence_id=state.sequence_id,through_time=state.last_time,
                   note='Final-known lineage roots, not immutable storm identities. Original family_id and branch_id are retained.')
    (output/'summary.json').write_text(json.dumps(summary,indent=2))
    (output/'SUCCESS').write_text('Canonical family export completed\n')
    return summary


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(export(args.source,args.output_dir),indent=2))
