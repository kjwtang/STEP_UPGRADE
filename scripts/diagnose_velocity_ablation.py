#!/usr/bin/env python3
"""Paired, exact tracking replay for velocity reset on/off ablation.

Select reference continuations newly lost after disabling resets, up to three
recovered controls, and adjacent edges newly classified as split/merge/complex.
Selection uses frame/local-label identities, never shared branch or node IDs.
"""
import argparse
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from compare_stream_tracks import compare, read_catalog
from replay_stream_pairs import replay
from step.tracking import load_tracking_state


def validate_settings(before, after):
    settings=[]
    for root in (before,after):
        parts=sorted(root.glob('chunk_*'))
        if not parts:
            raise ValueError(f'No chunks in {root}')
        state=load_tracking_state(parts[-1]/'state.json')
        c=dict(state.tracking_config)
        reset=c.pop('velocity_reset_morphology_v1',
                    float(bool(c.get('adjacent_overlap_first_v1'))))
        settings.append((c,reset,state.sequence_id))
    if settings[0][0]!=settings[1][0] or settings[0][2]!=settings[1][2]:
        raise ValueError('Runs must differ only in effective velocity-reset policy')
    if settings[0][1]!=1 or settings[1][1]!=0:
        raise ValueError('Expected reset-enabled before and reset-disabled after')


def select_pairs(old,new,ref):
    original={k for k,v in ref.items() if v=='continue'}
    lost=sorted((original & old.keys())-new.keys())
    recovered=sorted((original-old.keys()) & new.keys())[:3]
    events=sorted(k for k,v in new.items() if v in ('split','merge','complex')
                  and old.get(k)!=v and k[2]==k[0]+1)
    selected=sorted(set(lost+recovered+events))
    return selected,dict(newly_lost_reference=lost,recovered_controls=recovered,
                         new_event_pairs=events)


def diagnose(before,after,reference,output):
    if output.exists():
        raise FileExistsError(output)
    validate_settings(before,after)
    print('Checking identical inputs and object measurements across all three runs',flush=True)
    with TemporaryDirectory() as tmp:
        # Reuse the strict three-way identification/measurement check. Do not
        # leave a successful top-level report if either replay later fails.
        with redirect_stdout(io.StringIO()):
            compare(before,after,Path(tmp)/'check',reference=reference)
    old_ids,_,old=read_catalog(before)
    new_ids,_,new=read_catalog(after)
    _,_,ref=read_catalog(reference)
    keys,groups=select_pairs(old,new,ref)
    if not keys:
        raise ValueError('No selected adjacent pairs for this comparison')
    output.mkdir(parents=True)
    reports=[]
    for name,root,ids in (('reset_on',before,old_ids),('reset_off',after,new_ids)):
        lookup={identity:nid for nid,identity in ids.items()}
        pairs=[(lookup[k[:2]],lookup[k[2:]]) for k in keys]
        print(f'\n{name}: exact replay of {len(pairs)} selected pairs',flush=True)
        reports.append(replay(root,output/name,pairs))
    maps=[{tuple(r['object_pair']):r for r in report['pairs']} for report in reports]
    paired=[dict(object_pair=k,reset_on=maps[0][k],reset_off=maps[1][k]) for k in keys]
    summary=dict(both_replays_equal=True,selected_pairs=len(keys),groups=groups,
                 source_before=str(before.resolve()),source_after=str(after.resolve()),
                 pairs=paired,note='Diagnostics, not automatic causal or meteorological truth. Motion histories can diverge after earlier decisions. Raw overlap is unchanged input evidence; advected overlap depends on estimated velocity.')
    (output/'paired_diagnosis.json').write_text(json.dumps(summary,indent=2))
    lines=['# Velocity reset paired diagnosis','',
           'Both full-stream raster and catalog replays passed.',
           'Pair notation: (parent frame, local label) -> (child frame, local label).','']
    for item in paired:
        k=item['object_pair']; a=item['reset_on']; b=item['reset_off']
        lines.append(f'## {k[:2]} -> {k[2:]}')
        lines.append(f"Areas {a['parent_area']} -> {a['child_area']}; raw intersection {a['raw_intersection_cells']} cells; "
                     f"raw coverage parent/child={a['raw_parent_coverage']:.3f}/{a['raw_child_coverage']:.3f}.")
        for name,r in (('ON',a),('OFF',b)):
            edges=','.join(e['event'] for e in r['accepted_edges']) or 'none'
            lines.append(f"- {name}: {r['decision']} ({edges}); score={r['score_if_evaluated']:.4f}; "
                f"velocity(y,x)=({r['velocity_y_cells_per_frame']:.2f},{r['velocity_x_cells_per_frame']:.2f}); "
                f"prediction error={r['prediction_error_cells']:.2f}; adv IoU={r['advected_iou']:.4f}; "
                f"strong raw={r['strong_raw_overlap']}; event eligible={r['event_candidate']}.")
        lines.append('')
    lines.extend([summary['note'],''])
    (output/'REPORT.md').write_text('\n'.join(lines))
    (output/'SUCCESS').write_text('Both exact replays completed.\n')
    print(f'\nCompleted: {output}/REPORT.md',flush=True)
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('before',type=Path)
    p.add_argument('after',type=Path)
    p.add_argument('--reference',required=True,type=Path)
    p.add_argument('--output-dir',required=True,type=Path)
    a=p.parse_args()
    diagnose(a.before,a.after,a.reference,a.output_dir)
