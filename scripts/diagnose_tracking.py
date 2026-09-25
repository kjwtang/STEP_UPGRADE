#!/usr/bin/env python3
"""Read saved comparison data and explain tracker decisions without changing them."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from step.tracking import Tracker, Candidate, _objects, _pixel_overlap


def evidence(tracker, parent, child, time):
    dt=time-parent.node.time
    velocity=np.asarray(tracker._velocity(parent))
    previous=np.asarray(parent.node.object.centroid)
    predicted=previous+velocity*dt
    actual=float(np.linalg.norm(np.asarray(child.centroid)-previous))
    distance=float(np.linalg.norm(np.asarray(child.centroid)-predicted))
    radius=tracker.max_displacement*dt
    raw,rpc,rcc=_pixel_overlap(parent.node.object,child)
    shift=tuple(int(round(v*dt)) for v in velocity)
    adv,apc,acc=_pixel_overlap(parent.node.object,child,shift)
    intensity=min(parent.node.object.mean_intensity,child.mean_intensity)/max(parent.node.object.mean_intensity,child.mean_intensity,1e-12)
    score=.6*max(raw,adv)+.25*np.exp(-distance/max(radius,1.))+.15*intensity
    return dict(actual_displacement_cells=actual,prediction_error_cells=distance,search_radius_cells=radius,
        velocity_y_cells_per_frame=float(velocity[0]),velocity_x_cells_per_frame=float(velocity[1]),
        raw_iou=raw,advected_iou=adv,parent_coverage=max(rpc,apc),child_coverage=max(rcc,acc),
        intensity_ratio=intensity,score_if_evaluated=float(score),
        passes_score=bool(score>=tracker.tau),
        passes_event_overlap=bool(max(rpc,apc,rcc,acc)>=tracker.event_overlap),
        parent_area=parent.node.object.area,child_area=child.area)


def old_id(frame,mask):
    values=np.unique(frame[mask])
    if len(values)!=1 or values[0]<=0:
        raise ValueError("Original tracking does not have exactly one positive ID per identified object")
    return int(values[0])


def assignment_example():
    scores=[(0,0,.9),(0,1,.34),(1,0,.8),(1,1,0.)]
    candidates=[Candidate(i,j,s,0,0,0,0,0) for i,j,s in scores]
    got=Tracker._assignment(candidates,[0,1],[0,1],.35)
    return dict(tau=.35,scores=scores,actual_selected=[(x.parent_index,x.child_index,x.score) for x in got],
        feasible_higher_score_selection=[(0,0,.9)],
        explanation="The .34 subthreshold edge participates in assignment before filtering, displacing .90 in favor of .80. This reproduces a code issue, not proof of its role in the RCC case.")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("comparison",type=Path,help="Completed original_vs_upgrade_600_mem24 directory")
    p.add_argument("--output-dir",required=True,type=Path)
    p.add_argument("--frames",nargs="+",type=int,default=[1,5],help="Child frame indices to audit; default 00->01 and 04->05")
    a=p.parse_args()
    root=a.comparison.resolve()
    config=json.loads((root/"configuration.json").read_text())
    metadata=json.loads((root/"input_metadata.json").read_text())
    statuses=json.loads((root/"stage_results.json").read_text())
    for stage in ["new_id","new_track","new_id_original_track"]:
        if statuses[stage]["status"]!="completed":
            raise ValueError(f"Need completed {stage}; use the mem24 run")
    rain=np.load(root/"rain_mm_h.npy",mmap_mode="r")
    labels=np.load(root/"new_id/labels.npy",mmap_mode="r")
    saved=np.load(root/"new_track/labels.npy",mmap_mode="r")
    original=np.load(root/"new_id_original_track/labels.npy",mmap_mode="r")
    if not rain.shape==labels.shape==saved.shape==original.shape:
        raise ValueError("Saved shapes differ")
    if hashlib.sha256(rain.tobytes()).hexdigest()!=config["input_sha256"]:
        raise ValueError("Saved rain differs from comparison input hash")
    if any(t<1 or t>=len(rain) for t in a.frames):
        raise ValueError("Frames must be child-frame indices between 1 and N-1")
    spec=json.loads((root/"new_track/spec.json").read_text())
    tracker=Tracker(tau=spec["new_tau"],max_displacement=spec["new_displacement"],
                    max_gap=1,gap_tau=.45,gap_ambiguity=.05,event_overlap=.10,sequence_id="original_comparison")
    rows=[]
    frame_checks=[]
    for t in range(len(rain)):
        data=np.where(rain[t]>=config["threshold"],rain[t],0).astype(np.float32)
        parents=list(tracker.state.active)
        children=_objects(labels[t],data) if t in a.frames else []
        candidates=tracker._candidates(parents,children,t) if children else []
        candidate_map={(c.parent_index,c.child_index):c for c in candidates}
        pending=[]
        if t in a.frames:
            print("Auditing transition",t-1,"->",t,flush=True)
            for pi,parent in enumerate(parents):
                before=old_id(original[t-1],labels[t-1]==parent.node.object.label)
                for ci,child in enumerate(children):
                    after=old_id(original[t],labels[t]==child.label)
                    e=evidence(tracker,parent,child,t)
                    admitted=(pi,ci) in candidate_map
                    if admitted and not np.isclose(e["score_if_evaluated"],candidate_map[(pi,ci)].score,rtol=0,atol=1e-12):
                        raise ValueError("Diagnostic score does not reproduce actual candidate score")
                    pending.append(dict(frame=t,timestamp=metadata["timestamps"][t],
                        parent_node_id=parent.node.node_id,parent_local_label=parent.node.object.label,
                        child_local_label=child.label,parent_branch=parent.node.branch_id,
                        old_parent_id=before,old_child_id=after,original_same_id=before==after,
                        admitted_to_candidates=admitted,**e))
        raster,graph=tracker.update(t,labels[t],data)
        if not np.array_equal(raster,saved[t]):
            raise ValueError(f"Replay differs from saved labels at frame {t}. Do not interpret diagnostics from a different algorithm/configuration.")
        frame_checks.append(t)
        nodes={n.object.label:n for n in graph.objects}
        edges={(e.parent_id,e.child_id):e for e in graph.edges}
        parents_used={e.parent_id for e in graph.edges}
        children_used={e.child_id for e in graph.edges}
        for row in pending:
            child=nodes[row["child_local_label"]]
            edge=edges.get((row["parent_node_id"],child.node_id))
            row.update(child_node_id=child.node_id,child_branch=child.branch_id,
                       accepted_edge=edge.event if edge else "",same_branch=row["parent_branch"]==child.branch_id)
            if edge:
                reason="accepted_"+edge.event
            elif not row["admitted_to_candidates"]:
                reason="outside_prediction_gate"
            elif not row["passes_score"]:
                reason="below_score_threshold"
            elif row["parent_node_id"] in parents_used or child.node_id in children_used:
                reason="endpoint_used_by_other_edge"
            else:
                reason="eligible_but_not_assigned"
            row["decision"]=reason
            row["other_edges_from_parent"]=";".join(f"{e.child_id}:{e.event}" for e in graph.edges if e.parent_id==row["parent_node_id"] and e.child_id!=child.node_id)
            row["other_edges_to_child"]=";".join(f"{e.parent_id}:{e.event}" for e in graph.edges if e.child_id==child.node_id and e.parent_id!=row["parent_node_id"])
            rows.append(row)
    # Write only after all frames reproduce the saved result exactly.
    out=a.output_dir.resolve(); out.mkdir(parents=True,exist_ok=False)
    with (out/"pair_decisions.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    focus=[r for r in rows if (r["frame"],r["parent_branch"],r["child_branch"]) in {(1,12,26),(5,36,76)}]
    disagreements=[r for r in rows if r["original_same_id"] and not r["accepted_edge"]]
    summary=dict(replay_equal=True,frames_replayed=len(frame_checks),frames_audited=a.frames,
        comparison=str(root),tracking_parameters=tracker.state.tracking_config,
        original_same_id_without_new_edge_counts=dict(Counter(r["decision"] for r in disagreements)),
        focus_pairs=focus,assignment_counterexample=assignment_example(),
        interpretation="Original same-ID relations are not ground truth and can be many-to-many. Scores outside gate are counterfactual only. Rejection reasons are pair-level; endpoint competition requires inspecting other edges. No algorithm or tracking parameters were changed.")
    (out/"diagnosis.json").write_text(json.dumps(summary,indent=2))
    lines=["# Tracking diagnosis", "",f"Exact raster replay passed for {len(frame_checks)} frames.",
           "", "## Focus transitions", ""]
    for r in focus:
        lines.append(f"- Frame {r['frame']-1}->{r['frame']}, branch {r['parent_branch']}->{r['child_branch']}: "
            f"{r['decision']}; actual displacement={r['actual_displacement_cells']:.3f} cells, "
            f"prediction error={r['prediction_error_cells']:.3f}, gate={r['search_radius_cells']:.1f}; "
            f"raw IoU={r['raw_iou']:.4f}, advected IoU={r['advected_iou']:.4f}, "
            f"score={r['score_if_evaluated']:.4f}, tau={tracker.tau:.3f}.")
    if not focus:
        lines.append("Requested screenshot branch pairs absent; inspect pair_decisions.csv rather than assuming IDs.")
    lines += ["",summary["interpretation"],"", "A synthetic assignment counterexample is recorded separately. It establishes a post-filtering issue in the current code, not its contribution to these real-data failures."]
    (out/"REPORT.md").write_text("\n".join(lines))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=="__main__":
    main()
