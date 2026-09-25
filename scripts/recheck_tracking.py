#!/usr/bin/env python3
"""Re-run only upgraded tracking; preserve and reuse all old comparison outputs."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from compare_original_step import run_stage, make_report, dump
from validate_real_data import provenance
from step.tracking import track_with_graph, save_tracking_state, load_tracking_state


def edge_map(directory):
    with (directory/"objects.csv").open() as f:
        nodes={int(r["node_id"]):(int(r["time"]),int(r["local_label"])) for r in csv.DictReader(f)}
    with (directory/"edges.csv").open() as f:
        return {(*nodes[int(r["parent_node_id"])],*nodes[int(r["child_node_id"])]):r["event"] for r in csv.DictReader(f)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("comparison",type=Path)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--chunk-frames",type=int,default=6)
    p.add_argument("--plot-frames",type=int,default=6)
    a=p.parse_args()
    if a.chunk_frames<1 or a.plot_frames<1:
        p.error("Chunk and plot counts must be positive")
    source=a.comparison.resolve()
    config=json.loads((source/"configuration.json").read_text())
    metadata=json.loads((source/"input_metadata.json").read_text())
    old_status=json.loads((source/"stage_results.json").read_text())
    for name in ("original_id","new_id","original_track","new_id_original_track","new_track"):
        if old_status[name]["status"]!="completed":
            raise ValueError("Use the completed mem24 comparison directory")
    rain=np.load(source/"rain_mm_h.npy",mmap_mode="r")
    if hashlib.sha256(rain.tobytes()).hexdigest()!=config["input_sha256"]:
        raise ValueError("Saved precipitation hash changed")
    out=a.output_dir.resolve(); out.mkdir(parents=True,exist_ok=False)
    # Reused arrays/stages are read-only inputs. Do not overwrite or move sources.
    for name in ("rain_mm_h.npy","input_metadata.json","original_id","new_id","original_track","new_id_original_track"):
        (out/name).symlink_to(source/name,target_is_directory=(source/name).is_dir())
    config.update(new_version=provenance(a),rerun_from=str(source),rain=str(out/"rain_mm_h.npy"))
    dump(out/"configuration.json",config)
    results={k:dict(v,reused=True) for k,v in old_status.items() if k!="new_track"}
    results["new_track"]=run_stage("new_track",config,out,out/"new_id/labels.npy")
    dump(out/"stage_results.json",results)
    if results["new_track"]["status"]!="completed":
        raise SystemExit("New tracking failed; see new_track/run.log")
    print("Verifying chunk/checkpoint equivalence",flush=True)
    identified=np.load(out/"new_id/labels.npy",mmap_mode="r")
    whole=np.load(out/"new_track/labels.npy",mmap_mode="r")
    state=None
    combined=None
    kwargs=dict(tau=config["new_tau"],km=config["new_displacement"],max_gap=1,
                gap_tau=.45,gap_ambiguity=.05,event_overlap=.10,sequence_id="original_comparison")
    for start in range(0,len(rain),a.chunk_frames):
        end=min(len(rain),start+a.chunk_frames)
        precip=np.where(rain[start:end]>=config["threshold"],rain[start:end],0).astype(np.float32)
        raster,graph,state=track_with_graph(identified[start:end],precip,state=state,return_state=True,**kwargs)
        if not np.array_equal(raster,whole[start:end]):
            raise ValueError(f"Chunk equivalence failed at {start}:{end}")
        if combined is None:
            combined=type(graph)()
        combined.objects.extend(graph.objects); combined.edges.extend(graph.edges); combined.events.extend(graph.events)
        combined.family_map=graph.family_map
        save_tracking_state(state,out/"chunk_state.json")
        state=load_tracking_state(out/"chunk_state.json")
    # Write chunk graph in a distinct directory to compare exact CSV contents.
    from run_npy_validation import write_graph
    chunks=out/"chunk_verification"; chunks.mkdir()
    write_graph(combined,chunks)
    for filename in ("objects.csv","edges.csv","events.csv","family_map.csv"):
        if (chunks/filename).read_bytes()!=(out/"new_track"/filename).read_bytes():
            raise ValueError(f"Chunk graph differs: {filename}")
    old_edges=edge_map(source/"new_track")
    new_edges=edge_map(out/"new_track")
    old_raster=np.load(source/"new_track/labels.npy",mmap_mode="r")
    focus=[]
    for t,pb,cb in ((1,12,26),(5,36,76)):
        if t>=len(rain):
            continue
        parents=np.unique(identified[t-1][old_raster[t-1]==pb])
        children=np.unique(identified[t][old_raster[t]==cb])
        for parent in parents[parents>0]:
            for child in children[children>0]:
                key=(t-1,int(parent),t,int(child))
                focus.append(dict(object_pair=key,old_branch_pair=[pb,cb],old_edge=old_edges.get(key),new_edge=new_edges.get(key)))
    def records(keys,mapping):
        return [dict(object_pair=k,event=mapping[k]) for k in sorted(keys)]
    summary=dict(source=str(source),chunk_rasters_and_graph_equal=True,
        old_edges=len(old_edges),new_edges=len(new_edges),
        added=records(new_edges.keys()-old_edges.keys(),new_edges),
        removed=records(old_edges.keys()-new_edges.keys(),old_edges),
        changed_event=[dict(object_pair=k,old=old_edges[k],new=new_edges[k]) for k in sorted(old_edges.keys()&new_edges.keys()) if old_edges[k]!=new_edges[k]],
        new_edge_types=dict(Counter(new_edges.values())),focus_pairs=focus,
        note="Pairs use immutable frame/local-label identities, not changed branch IDs. Added links are not automatically correct. Original stages and identification were reused, not rerun.")
    dump(out/"edge_changes.json",summary)
    make_report(out,results,metadata,a.plot_frames)
    with (out/"REPORT.md").open("a") as f:
        f.write("\n\n## Tracking revision recheck\n\nOriginal stages and identification were reused from "+str(source)+
                ". Only new tracking was rerun. Raster and graph chunk/checkpoint equivalence passed. "
                "See edge_changes.json for added, removed and changed-event edges by frame/local-label identity.\n")
    (out/"SUCCESS").touch()
    print(json.dumps({k:v for k,v in summary.items() if k not in ("added","removed","changed_event")},indent=2),flush=True)
    print("Results:",out,flush=True)


if __name__=="__main__":
    main()
