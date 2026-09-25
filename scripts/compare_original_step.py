#!/usr/bin/env python3
"""Pinned RDCEP STEP vs STEP_UPGRADE on one identical bounded rain cube."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

import numpy as np

PIN = "ef75c083addcddbcabf6209dc95b1dff89055b0c"
ROOT = Path(__file__).resolve().parents[1]


def dump(path, value):
    Path(path).write_text(json.dumps(value,indent=2,default=str))


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def worker(spec_path):
    s = json.loads(Path(spec_path).read_text())
    # Hard per-process virtual-memory ceiling on RCC/Linux; RSS watchdog also applies.
    if sys.platform.startswith("linux"):
        limit = int(s["memory_gb"]*1024**3)
        resource.setrlimit(resource.RLIMIT_AS,(limit,limit))
    from run_npy_validation import disk, write_graph
    source = np.load(s["rain"],mmap_mode="r")
    precip = np.where(source>=s["threshold"],source,0).astype(np.float32)
    stage = s["stage"]
    start = time.perf_counter()
    if stage=="new_id":
        from step.identification import identify
        result = identify(precip,disk(s["radius"]),workers=s["workers"],min_size=1)
    elif stage=="original_id":
        module = load_module(Path(s["original"])/"step/identification.py","original_step_identification")
        result = module.identify(precip,disk(s["radius"]))
    elif stage=="new_track":
        from step.tracking import track_with_graph
        labels = np.load(s["labels"],mmap_mode="r")
        result, graph = track_with_graph(labels,precip,tau=s["new_tau"],km=s["new_displacement"],
            max_gap=1,gap_tau=.45,gap_ambiguity=.05,event_overlap=.10,sequence_id="original_comparison")
    else:
        module = load_module(Path(s["original"])/"step/tracking.py","original_step_tracking")
        labels = np.load(s["labels"],mmap_mode="r")
        result = module.track(labels,precip,tau=s["original_tau"],phi=s["original_phi"],km=s["original_km"],test=False)
    core = time.perf_counter()-start
    if stage=="new_track":
        write_graph(graph,Path(s["stage_dir"]))
    if result.shape!=precip.shape or not np.isfinite(result).all() or np.any(result<0):
        raise ValueError("Algorithm returned invalid label array")
    if not np.equal(result,np.floor(result)).all():
        raise ValueError("Algorithm returned noninteger labels")
    np.save(s["result"],result)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=="darwin" else 1024)
    dump(Path(s["stage_dir"])/"worker_metrics.json",dict(core_seconds=core,parent_peak_rss_mib=rss))


def run_stage(name, config, out, labels=None):
    import psutil
    directory = out/name
    directory.mkdir()
    spec = dict(config,stage=name,stage_dir=str(directory),result=str(directory/"labels.npy"))
    if labels is not None:
        spec["labels"] = str(labels)
    dump(directory/"spec.json",spec)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
    start = time.perf_counter()
    peak = 0
    errors = set()
    status = "running"
    print(f"START {name}: timeout={config['timeout_seconds']}s, memory={config['memory_gb']}GiB",flush=True)
    with (directory/"run.log").open("w") as log:
        proc = subprocess.Popen([sys.executable,str(Path(__file__).resolve()),"--worker",str(directory/"spec.json")],
            stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        try:
            while proc.poll() is None:
                rss = 0
                try:
                    parent = psutil.Process(proc.pid)
                    for p in [parent]+parent.children(recursive=True):
                        try:
                            rss += p.memory_info().rss
                        except psutil.NoSuchProcess:
                            pass
                except (psutil.Error,OSError) as e:
                    errors.add(type(e).__name__)
                peak = max(peak,rss)
                elapsed = time.perf_counter()-start
                if elapsed>config["timeout_seconds"]:
                    status="timeout"
                elif rss>config["memory_gb"]*1024**3:
                    status="memory_limit"
                if status!="running":
                    try:
                        os.killpg(proc.pid,signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    break
                time.sleep(.2)
            code = proc.wait()
        finally:
            if proc.poll() is None:
                try:
                    os.killpg(proc.pid,signal.SIGKILL)
                except ProcessLookupError:
                    pass
                proc.wait()
    if status=="running":
        status="completed" if code==0 and (directory/"labels.npy").exists() else "failed"
    metrics = dict(stage=name,status=status,returncode=code,wall_seconds=time.perf_counter()-start,
                   sampled_tree_peak_rss_mib=peak/1024**2 if peak else None,memory_sampling_errors=sorted(errors))
    if status=="completed":
        metrics.update(json.loads((directory/"worker_metrics.json").read_text()))
    dump(directory/"status.json",metrics)
    print(json.dumps(metrics),flush=True)
    return metrics


def partition_metrics(a,b):
    both=(a>0)&(b>0)
    union=(a>0)|(b>0)
    av,bv=a[both],b[both]
    n=len(av)
    ari=None
    if n>=2:
        _,ac=np.unique(av,return_counts=True)
        _,bc=np.unique(bv,return_counts=True)
        _,joint=np.unique(np.stack([av,bv],axis=1),axis=0,return_counts=True)
        pairs=lambda counts: float(np.sum(counts.astype(float)*(counts-1)/2))
        total=n*(n-1)/2
        expected=pairs(ac)*pairs(bc)/total
        denom=(pairs(ac)+pairs(bc))/2-expected
        ari=(pairs(joint)-expected)/denom if denom else 1.
    return dict(original_objects=int(np.count_nonzero(np.unique(a))),
                new_objects=int(np.count_nonzero(np.unique(b))),
                wet_mask_iou=float(both.sum()/union.sum()) if union.any() else 1.,
                common_wet_pixels=n,ari_common_wet=ari)


def same_id_links(identified,tracked):
    """Observed adjacent-frame same-ID associations, not inferred truth."""
    previous={}
    links=set()
    reused_within_frame=0
    for t in range(len(identified)):
        current={}
        for label in np.unique(identified[t]):
            if not label:
                continue
            values=np.unique(tracked[t][identified[t]==label])
            if len(values)!=1 or values[0]==0:
                raise ValueError("Tracked ID does not map uniquely to identified object")
            current.setdefault(int(values[0]),[]).append(int(label))
        reused_within_frame+=sum(len(v)>1 for v in current.values())
        for track,labels in current.items():
            for before in previous.get(track,[]):
                for after in labels:
                    links.add((t,int(before),int(after)))
        previous=current
    return links,reused_within_frame


def make_report(out, results, meta, plot_frames):
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    rain=np.load(out/"rain_mm_h.npy",mmap_mode="r")
    maps={name:np.load(out/name/"labels.npy",mmap_mode="r") for name,r in results.items() if r["status"]=="completed"}
    metrics=[]
    if "original_id" in maps and "new_id" in maps:
        metrics=[dict(frame=t,timestamp=meta["timestamps"][t],**partition_metrics(maps["original_id"][t],maps["new_id"][t])) for t in range(len(rain))]
        with (out/"identification_comparison.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(metrics[0])); w.writeheader(); w.writerows(metrics)
    tracking_stats={}
    for stage,id_stage in [("original_track","original_id"),("new_id_original_track","new_id"),("new_track","new_id")]:
        if stage in maps:
            links,duplicates=same_id_links(maps[id_stage],maps[stage])
            tracking_stats[stage]=dict(adjacent_same_id_object_links=len(links),
                shared_id_frame_instances=duplicates,unique_raster_ids=int(np.count_nonzero(np.unique(maps[stage]))))
    if "new_id_original_track" in maps and "new_track" in maps:
        old,_=same_id_links(maps["new_id"],maps["new_id_original_track"])
        new,_=same_id_links(maps["new_id"],maps["new_track"])
        dump(out/"tracking_link_disagreements.json",dict(
            definition="Same identified objects; compare adjacent same-ID associations. New split/merge edges end branches and are not counted here; see new_track/edges.csv.",
            common=len(old&new),original_only=sorted(old-new),new_only=sorted(new-old)))
    dump(out/"tracking_summary.json",tracking_stats)
    cmap=ListedColormap(np.random.default_rng(12).uniform(.15,.95,(256,3)))
    cmap.set_bad("white")
    positive=rain[rain>0]
    vmax=max(2.,float(np.percentile(positive,99))) if positive.size else 2.
    for t in range(min(plot_frames,len(rain))):
        for kind,stages in [("identification",["original_id","new_id"]),
                            ("tracking",["original_track","new_id_original_track","new_track"])]:
            fig,axes=plt.subplots(1,len(stages)+1,figsize=(5*(len(stages)+1),5),layout="constrained")
            im=axes[0].imshow(rain[t],origin="lower",cmap="Blues",vmin=0,vmax=vmax)
            axes[0].set_title("Rain (mm/h)"); fig.colorbar(im,ax=axes[0])
            for ax,name in zip(axes[1:],stages):
                ax.set_title(name.replace("_"," "),fontsize=11)
                if name not in maps:
                    ax.text(.5,.5,results[name]["status"],ha="center",transform=ax.transAxes)
                else:
                    frame=maps[name][t]
                    ax.imshow(np.ma.masked_where(frame==0,frame%256),origin="lower",cmap=cmap,vmin=0,vmax=255)
                    for label in np.unique(frame):
                        if label:
                            yy,xx=np.nonzero(frame==label)
                            ax.text(xx.mean(),yy.mean(),str(int(label)),fontsize=6)
                ax.set_xlabel("x (same cropped cells)")
            fig.suptitle(f"{meta['timestamps'][t]} | {kind}; IDs/colors independent across algorithms")
            fig.savefig(out/f"{kind}_{t:03d}.png",dpi=120); plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4),layout="constrained")
    if metrics:
        axes[0].plot([r["original_objects"] for r in metrics],label="Original identification")
        axes[0].plot([r["new_objects"] for r in metrics],label="New identification")
        axes[0].legend()
    axes[0].set(xlabel="Frame",ylabel="Identified objects")
    names=list(results)
    axes[1].bar(names,[results[n].get("core_seconds",np.nan) for n in names])
    axes[1].tick_params(axis="x",rotation=35)
    axes[1].set_ylabel("Core seconds (completed stages only)")
    fig.savefig(out/"comparison_summary.png",dpi=140); plt.close(fig)
    report=["# Original STEP comparison", "", "## Stage outcomes", ""]
    for name,r in results.items():
        report.append(f"- {name}: {r['status']}; wall={r.get('wall_seconds',0):.2f}s; core={r.get('core_seconds','N/A')}s")
    report += ["", "The original is a baseline, not ground truth. Shared inputs use identical cumulative differencing, "
        "crop, threshold and disk footprint. Original tutorial tracking scores are not calibrated for this hourly dataset; "
        "tau and distance semantics differ between implementations. Same numerical tau would not ensure equivalence. "
        "Branch IDs/colors are independent between algorithms. New branches end at split/merge; do not interpret every "
        "ID change as a failure. Within-frame shared original IDs can represent multiple objects; check label maps. "
        "Same-ID link disagreements isolate tracking on new identification, but do not establish physical correctness. "
        "Timeout/failed stages have no accuracy or speedup claim. Read run.log for exceptions. "
        "Memory sampling sums RSS and can double-count shared pages. Linux workers additionally have an RLIMIT_AS ceiling. "
        "Statistics use the entire selected sample; preview PNGs show the first requested frames only."]
    (out/"REPORT.md").write_text("\n".join(report))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",type=Path,help="CSTM hourly directory")
    p.add_argument("--original-dir",required=True,type=Path)
    p.add_argument("--output-dir",required=True,type=Path)
    p.add_argument("--start",type=int,default=0)
    p.add_argument("--hours",type=int,default=24)
    p.add_argument("--crop",type=int,default=600)
    p.add_argument("--threshold",type=float,default=1.)
    p.add_argument("--radius",type=int,default=9)
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--timeout-seconds",type=float,default=900,help="Per algorithm stage, including startup")
    p.add_argument("--memory-gb",type=float,default=8,help="Per-worker virtual memory and sampled process-tree RSS limits in GiB")
    p.add_argument("--original-tau",type=float,default=.7,help="Original tutorial starting point, not calibrated")
    p.add_argument("--original-phi",type=float,default=.003,help="Original tutorial starting point, not calibrated")
    p.add_argument("--original-km",type=float,default=30.,help="120km / 4km grid =30 cells; original code has directional fallback")
    p.add_argument("--new-tau",type=float,default=.35)
    p.add_argument("--new-displacement",type=float,default=20.)
    p.add_argument("--plot-frames",type=int,default=6)
    a=p.parse_args()
    if a.crop<1 or a.crop>600:
        p.error("Comparison is limited to the current <=600-cell domain; no full-domain original tracking")
    if a.hours<2 or a.workers<1 or a.radius<0 or a.start<0 or a.plot_frames<1:
        p.error("Invalid frame/worker/radius/start/plot count")
    for key in ["threshold","timeout_seconds","memory_gb","original_tau","original_phi","original_km","new_tau","new_displacement"]:
        value=getattr(a,key)
        if not np.isfinite(value) or value<=0:
            p.error(f"{key} must be finite and positive")
    original=a.original_dir.resolve()
    commit=subprocess.check_output(["git","-C",str(original),"rev-parse","HEAD"],text=True).strip()
    if commit!=PIN:
        p.error(f"Original commit must be {PIN}; got {commit}")
    for file in ["step/identification.py","step/tracking.py"]:
        committed=subprocess.check_output(["git","-C",str(original),"show",f"{PIN}:{file}"])
        if (original/file).read_bytes()!=committed:
            p.error(f"Original source has local modifications: {file}")
    out=a.output_dir.resolve(); out.mkdir(parents=True,exist_ok=False)
    from types import SimpleNamespace
    from wrf_rain_input import read_wrf
    from validate_real_data import provenance
    start=time.perf_counter()
    rain,meta=read_wrf(SimpleNamespace(input=a.input,inspect=False,variable="RAINNC",rain_kind="cumulative",
        start=a.start,hours=a.hours,dt_hours=1.,units=None,crop=a.crop,negative_tolerance_mm=0.))
    if not np.isfinite(rain).all():
        raise ValueError("Original STEP lacks missing-data barriers: comparison refuses missing pixels rather than filling zero")
    np.save(out/"rain_mm_h.npy",rain)
    dump(out/"input_metadata.json",meta)
    config=dict(vars(a),original=str(original),original_commit=commit,rain=str(out/"rain_mm_h.npy"),
        source_sha256={name:hashlib.sha256((original/"step"/name).read_bytes()).hexdigest() for name in ["identification.py","tracking.py"]},
        input_sha256=hashlib.sha256(rain.tobytes()).hexdigest(),read_seconds=time.perf_counter()-start,
        new_version=provenance(a),original_parameter_source="https://github.com/relttira/STEP/wiki/Tutorial; tau/phi example only, km adapted to 4km grid")
    dump(out/"configuration.json",config)
    del rain
    results={}
    # Produce upgraded output even if original identification/tracking cannot finish.
    for name,dependency in [("new_id",None),("new_track","new_id"),("original_id",None),
                            ("original_track","original_id"),("new_id_original_track","new_id")]:
        if dependency and results[dependency]["status"]!="completed":
            results[name]=dict(stage=name,status="dependency_failed")
        else:
            results[name]=run_stage(name,config,out,None if dependency is None else out/dependency/"labels.npy")
        dump(out/"stage_results.json",results)
    make_report(out,results,meta,a.plot_frames)
    complete=all(r["status"]=="completed" for r in results.values())
    (out/("SUCCESS" if complete else "PARTIAL")).touch()
    print("Comparison report:",out/"REPORT.md",flush=True)
    if not complete:
        raise SystemExit(2)


if __name__=="__main__":
    if len(sys.argv)==3 and sys.argv[1]=="--worker":
        worker(sys.argv[2])
    else:
        main()
