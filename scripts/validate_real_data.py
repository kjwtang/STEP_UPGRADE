#!/usr/bin/env python3
"""Small-sample scientific diagnostics; does not certify scientific accuracy."""
import argparse
from collections import Counter
from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import time
import resource
import sys
import os
import platform

import numpy as np

from step.identification import identify
from step.tracking import track_with_graph, save_tracking_state, load_tracking_state
from run_npy_validation import disk, write_graph


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path, help="One NetCDF or (time,y,x) NPY file")
    p.add_argument("--inspect", action="store_true", help="Print metadata without loading rain or running STEP")
    p.add_argument("--variable", help="Explicit NetCDF precipitation variable")
    p.add_argument("--dims", nargs=3, metavar=("TIME", "Y", "X"), help="NetCDF dimension order; defaults to variable order")
    p.add_argument("--time-coordinate", help="Decoded one-dimensional time coordinate, if different from time dimension")
    p.add_argument("--rain-kind", choices=["interval", "rate"], help="interval=mm per interval; rate=mm/h; cumulative rain unsupported")
    p.add_argument("--units", choices=["mm", "mm/h"], help="Explicit confirmed units, overrides metadata")
    p.add_argument("--dt-hours", type=float, help="Required for NPY or undecodable time; checked against decoded time")
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--hours", type=int, default=72, help="Number of FRAMES (72 hours only for hourly data)")
    p.add_argument("--crop", type=int, default=300, help="Central square width; 0 uses full spatial domain")
    p.add_argument("--output-dir", type=Path)
    p.add_argument("--sequence-id", default="validation_case01")
    p.add_argument("--threshold", type=float, default=1., help="mm/h after normalization")
    p.add_argument("--bridge-radius", type=int, default=9)
    p.add_argument("--min-size", type=int, default=1)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--tau", type=float, default=.35)
    p.add_argument("--max-displacement", type=float, default=20., help="Grid cells per frame")
    p.add_argument("--max-gap", type=int, default=1)
    p.add_argument("--gap-tau", type=float, default=.45)
    p.add_argument("--gap-ambiguity", type=float, default=.05)
    p.add_argument("--event-overlap", type=float, default=.10)
    p.add_argument("--chunk-frames", type=int, default=24)
    p.add_argument("--grid-km", type=float, help="Square grid spacing; enables km/h and km2 statistics")
    p.add_argument("--gif", action="store_true", help="Also write all-frame GIF (slower)")
    p.add_argument("--stream", action="store_true", help="High-resolution benchmark: process one chunk at a time; no whole-run equivalence check")
    p.add_argument("--save-labels", action="store_true", help="Streaming mode: save per-chunk label arrays (extra disk I/O)")
    return p.parse_args()


def load_input(a):
    """Slice before loading, preserve missing pixels, reject ambiguous inputs."""
    dt, metadata = a.dt_hours, {}
    if a.input.suffix == ".npy":
        source = np.load(a.input, mmap_mode="r")
        if a.inspect:
            print({"shape": source.shape, "dtype": str(source.dtype)})
            return None
        if source.ndim != 3:
            raise ValueError("NPY must have shape (time,y,x)")
        shape = source.shape
        stamps = None
        units = a.units
        ds = None
    else:
        import xarray as xr
        ds = xr.open_dataset(a.input)
        if a.inspect:
            print(ds)
            for name, var in ds.data_vars.items():
                print(name, var.dims, var.shape, dict(var.attrs))
            ds.close()
            return None
        if not a.variable or a.variable not in ds:
            ds.close()
            raise ValueError("Choose --variable from --inspect output")
        if a.variable.upper() in {"RAINNC", "RAINC", "RAINSH"}:
            ds.close()
            raise ValueError("Cumulative WRF rain needs restart-aware preprocessing first")
        var = ds[a.variable]
        dims = tuple(a.dims or var.dims)
        if len(dims) != 3 or set(dims) != set(var.dims):
            ds.close()
            raise ValueError("Select a 3-D rain variable and specify --dims TIME Y X")
        source = var.transpose(*dims)
        shape = source.shape
        coordinate = ds.get(a.time_coordinate or dims[0])
        stamps = None
        if coordinate is not None and coordinate.dims == (dims[0],):
            vals = coordinate.values
            if np.issubdtype(vals.dtype, np.datetime64) or (len(vals) and hasattr(vals[0], "calendar")):
                stamps = vals
        units = a.units or var.attrs.get("units")
        metadata = {"variable": a.variable, "dims": dims, "source_units": var.attrs.get("units"),
                    "units_override": a.units, "source_attributes": dict(var.attrs)}
    try:
        if not a.rain_kind:
            raise ValueError("Confirm --rain-kind interval or rate; cumulative input is not supported")
        normalized = str(units).lower().replace(" ", "")
        allowed = {"interval": {"mm"}, "rate": {"mm/h", "mmhr-1", "mmh-1", "mmhour-1"}}
        if normalized not in allowed[a.rain_kind]:
            raise ValueError("Units missing or incompatible: confirm --units mm or mm/h and --rain-kind")
        end = a.start + a.hours
        if a.start < 0 or end > shape[0]:
            raise ValueError("Requested frames extend outside input; adjust --start/--hours")
        selected = None if stamps is None else stamps[a.start:end]
        if selected is not None:
            cadence_times = stamps[max(0,a.start-1):min(len(stamps),end+1)]
            if np.issubdtype(cadence_times.dtype, np.datetime64):
                if np.isnat(cadence_times).any():
                    raise ValueError("Missing timestamps")
                intervals = np.diff(cadence_times) / np.timedelta64(1, "h")
            else:
                intervals = np.array([(b-a).total_seconds()/3600 for a,b in zip(cadence_times[:-1], cadence_times[1:])])
            if not len(intervals):
                raise ValueError("Need at least two timestamps to verify cadence")
            if np.any(intervals <= 0) or not np.allclose(intervals, intervals[0]):
                raise ValueError("Duplicate, decreasing or irregular timestamps: do not bridge missing observations")
            measured = float(intervals[0])
            if dt is not None and not np.isclose(dt, measured):
                raise ValueError("--dt-hours disagrees with timestamps")
            dt = measured
        if dt is None or not np.isfinite(dt) or dt <= 0:
            raise ValueError("No decoded time: supply confirmed --dt-hours; cadence cannot be independently checked")
        width = a.crop or max(shape[1:])
        y = max(0, (shape[1]-width)//2)
        x = max(0, (shape[2]-width)//2)
        slices = (slice(a.start, end), slice(y, y+width), slice(x, x+width))
        data = np.array(source[slices], dtype=np.float32, copy=True)
        if np.any(data[np.isfinite(data)] < 0):
            raise ValueError("Negative precipitation: decode fill values or fix preprocessing first")
        if np.any(~np.isfinite(data).any(axis=(1,2))):
            raise ValueError("Entire missing observation frame: explicit missing-frame policy required")
        if a.rain_kind == "interval":
            data /= dt
        metadata.update(source_shape=list(shape), crop_origin_yx=[y,x], dt_hours=dt,
                        time_verified=selected is not None, normalized_units="mm/h",
                        timestamps=None if selected is None else [str(t) for t in selected])
        return data, metadata
    finally:
        if ds is not None:
            ds.close()


def signature(graph):
    nodes = []
    for n in graph.objects:
        item = asdict(n)
        item["family_id"] = graph.family_map[n.branch_id]
        nodes.append(item)
    return {"nodes": nodes, "edges": [asdict(e) for e in graph.edges],
            "events": [asdict(e) for e in graph.events], "family_map": graph.family_map}


def plots(data, labels, tracked, graph, dt, out, gif):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.animation import FuncAnimation, PillowWriter
    rng = np.random.default_rng(12)
    colors = rng.uniform(.15, .95, (256, 3))
    cmap = ListedColormap(colors)
    cmap.set_bad("white")
    positive = data[np.isfinite(data) & (data > 0)]
    vmax = max(2., float(np.percentile(positive, 99))) if positive.size else 2.
    nodes = {n.node_id: n for n in graph.objects}
    stamps = json.loads((out / "metadata.json").read_text())["timestamps"]
    def draw(t, axes):
        for ax in axes:
            ax.clear()
        im = axes[0].imshow(data[t], origin="lower", vmin=0, vmax=vmax, cmap="Blues")
        if np.any(labels[t]):
            from skimage.segmentation import find_boundaries
            border = find_boundaries(labels[t], mode="inner")
            axes[0].imshow(np.ma.masked_where(~border, border), origin="lower", cmap="autumn", vmin=0, vmax=1)
        axes[0].set_title("Rain (mm/h) + object boundaries")
        coded = np.ma.masked_where(tracked[t] == 0, tracked[t] % 256)
        axes[1].imshow(coded, origin="lower", cmap=cmap, vmin=0, vmax=255)
        for n in graph.objects:
            if n.time == t:
                axes[1].text(n.object.centroid[1], n.object.centroid[0], str(n.branch_id), fontsize=6)
        axes[1].set_title("Branch IDs (colors repeat every 256)")
        for ax in axes:
            ax.set_xlabel("x (cropped grid cells)")
            ax.set_ylabel("y (cropped grid cells)")
        axes[0].figure.suptitle(stamps[t] if stamps else f"Frame {t}; elapsed {t*dt:g} h")
        return im
    selected = set(np.linspace(0,len(data)-1,min(6,len(data)),dtype=int).tolist())
    # Include first example of each lineage event plus a gap's two endpoints.
    examples = {}
    for e in graph.events:
        examples.setdefault(e.event, e.time)
    for t in examples.values():
        selected.update([max(0,t-1),t])
    for e in graph.edges:
        if e.gap:
            selected.update([nodes[e.parent_id].time, nodes[e.child_id].time])
            break
    for t in sorted(selected):
        fig, axes = plt.subplots(1,2,figsize=(12,5),layout="constrained")
        im = draw(t, axes)
        fig.colorbar(im, ax=axes[0], label="mm/h")
        fig.savefig(out / f"frame_{t:03d}.png", dpi=140)
        plt.close(fig)
    if gif:
        fig, axes = plt.subplots(1,2,figsize=(10,4))
        ani = FuncAnimation(fig, lambda t: draw(t,axes), frames=len(data))
        ani.save(out / "tracking.gif", writer=PillowWriter(fps=3), dpi=90)
        plt.close(fig)


def peak_rss_mb():
    # Linux ru_maxrss is KiB; macOS uses bytes. Parent only, not worker sum.
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**2 if sys.platform == "darwin" else 1024)


def provenance(a):
    settings = vars(a).copy()
    settings.update(python=sys.version,host=platform.node(),platform=platform.platform(),
                    slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                    slurm_cpus_per_task=os.environ.get("SLURM_CPUS_PER_TASK"),
                    numpy_version=np.__version__)
    try:
        settings["git_commit"] = subprocess.check_output(["git","rev-parse","HEAD"],
            cwd=Path(__file__).resolve().parents[1],text=True).strip()
        settings["git_dirty"] = bool(subprocess.check_output(["git","status","--porcelain"],
            cwd=Path(__file__).resolve().parents[1],text=True).strip())
    except (OSError,subprocess.CalledProcessError):
        settings["git_commit"] = "unknown"
    return settings


def stream_benchmark(a):
    """Bound resident rain/label storage by chunk size, preserve tracking state."""
    import copy
    import csv
    import threading
    import psutil
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not a.output_dir:
        raise ValueError("--output-dir required")
    out = a.output_dir
    out.mkdir(parents=True, exist_ok=False)
    (out / "settings.json").write_text(json.dumps(provenance(a),indent=2,default=str))
    peak_tree = [0.]
    memory_errors = set()
    stop_monitor = threading.Event()
    process = psutil.Process()
    def monitor():
        while not stop_monitor.is_set():
            total = 0
            try:
                processes = [process] + process.children(recursive=True)
            except (psutil.Error, OSError) as error:
                memory_errors.add(type(error).__name__)
                processes = [process]
            for proc in processes:
                try:
                    total += proc.memory_info().rss
                except (psutil.Error, OSError) as error:
                    memory_errors.add(type(error).__name__)
            peak_tree[0] = max(peak_tree[0],total/1024**2)
            stop_monitor.wait(.2)
    thread = threading.Thread(target=monitor,daemon=True)
    thread.start()
    rows, state = [], None
    frame_stats = []
    reference_dt = None
    time_verified = True
    events, n_nodes, n_edges, n_gaps = Counter(), 0, 0, 0
    began = time.perf_counter()
    cpu_started = time.process_time()
    kwargs = dict(tau=a.tau, km=a.max_displacement, max_gap=a.max_gap, gap_tau=a.gap_tau,
                  gap_ambiguity=a.gap_ambiguity,event_overlap=a.event_overlap,sequence_id=a.sequence_id)
    try:
        for offset in range(0,a.hours,a.chunk_frames):
            b = copy.copy(a)
            b.start, b.hours = a.start+offset, min(a.chunk_frames,a.hours-offset)
            t = time.perf_counter()
            data, meta = load_input(b)
            read_s = time.perf_counter()-t
            if reference_dt is not None and not np.isclose(reference_dt,meta["dt_hours"]):
                raise ValueError("Cadence changes between chunks")
            reference_dt = meta["dt_hours"]
            time_verified &= meta["time_verified"]
            # Independently read the two frames around each file slice boundary.
            if offset and meta["time_verified"]:
                boundary = copy.copy(a)
                boundary.start, boundary.hours, boundary.crop = b.start-1, 2, 1
                check_data, check_meta = load_input(boundary)
                if not np.isclose(reference_dt,check_meta["dt_hours"]):
                    raise ValueError("Missing or irregular timestamp at chunk boundary")
                del check_data
                read_s = time.perf_counter()-t
            part_dir = out / f"chunk_{offset:06d}"
            part_dir.mkdir()
            (part_dir / "metadata.json").write_text(json.dumps(meta,indent=2,default=str))
            print(f"Chunk {offset}:{offset+b.hours}, shape={data.shape}: identification",flush=True)
            t = time.perf_counter()
            labels = identify(data,disk(a.bridge_radius),workers=a.workers,threshold=a.threshold,min_size=a.min_size)
            identify_s = time.perf_counter()-t
            print("Tracking",flush=True)
            t = time.perf_counter()
            tracked, graph, state = track_with_graph(labels,data,state=state,return_state=True,**kwargs)
            track_s = time.perf_counter()-t
            t = time.perf_counter()
            write_graph(graph,part_dir)
            save_tracking_state(state,part_dir / "state.json")
            state = load_tracking_state(part_dir / "state.json")
            if a.save_labels:
                np.save(part_dir / "identified_labels.npy",labels)
                np.save(part_dir / "tracked_labels.npy",tracked)
            write_s = time.perf_counter()-t
            t = time.perf_counter()
            if offset == 0:
                plots(data,labels,tracked,graph,meta["dt_hours"],part_dir,a.gif)
            plot_s = time.perf_counter()-t
            counts = Counter(n.time for n in graph.objects)
            for i in range(len(data)):
                frame_stats.append(dict(frame=a.start+offset+i,objects=counts[offset+i],
                    mean_rain_mm_h=float(np.nanmean(data[i])),max_rain_mm_h=float(np.nanmax(data[i])),
                    wet_fraction=float(np.mean(labels[i]>0)),missing_fraction=float(np.mean(~np.isfinite(data[i])))))
            n_nodes += len(graph.objects)
            n_edges += len(graph.edges)
            n_gaps += sum(e.gap>0 for e in graph.edges)
            events.update(e.event for e in graph.events)
            row = dict(start_frame=b.start,frames=b.hours,ny=data.shape[1],nx=data.shape[2],
                       read_seconds=read_s,identification_seconds=identify_s,tracking_seconds=track_s,
                       output_seconds=write_s,plot_seconds=plot_s,nodes=len(graph.objects),
                       max_objects_per_frame=max(counts.values(),default=0),
                       parent_peak_rss_mb=peak_rss_mb(),sampled_tree_peak_rss_mb=peak_tree[0] or None)
            rows.append(row)
            with (out / "performance.csv").open("w",newline="") as handle:
                writer = csv.DictWriter(handle,fieldnames=list(row))
                writer.writeheader()
                writer.writerows(rows)
            print(json.dumps(row),flush=True)
            del data, labels, tracked, graph
        total = time.perf_counter()-began
        core = sum(r["identification_seconds"]+r["tracking_seconds"] for r in rows)
        timed_keys = ["read_seconds","identification_seconds","tracking_seconds","output_seconds","plot_seconds"]
        totals = {key:sum(r[key] for r in rows) for key in timed_keys}
        summary = dict(mode="stream_benchmark",frames=a.hours,spatial_shape=[rows[0]["ny"],rows[0]["nx"]],
                       total_wall_seconds=total,parent_cpu_seconds=time.process_time()-cpu_started,
                       stage_seconds=totals,core_seconds_per_frame=core/a.hours,
                       core_million_cells_per_second=a.hours*rows[0]["ny"]*rows[0]["nx"]/max(core,1e-12)/1e6,
                       parent_peak_rss_mb=peak_rss_mb(),sampled_tree_peak_rss_mb=peak_tree[0] or None,
                       memory_sampling_errors=sorted(memory_errors),
                       nodes=n_nodes,edges=n_edges,gap_edges=n_gaps,events=dict(events),time_verified=time_verified,
                       chunk_equivalence_checked=False,
                       memory_note="Tree RSS sampled every 0.2s; shared fork pages are counted multiple times. Parent HWM excludes workers. Confirm job memory with Slurm accounting.")
        (out / "summary.json").write_text(json.dumps(summary,indent=2))
        with (out / "frame_statistics.csv").open("w",newline="") as handle:
            writer = csv.DictWriter(handle,fieldnames=list(frame_stats[0]))
            writer.writeheader()
            writer.writerows(frame_stats)
        fig,axes = plt.subplots(2,1,figsize=(10,6),layout="constrained")
        axes[0].plot([r["frame"] for r in frame_stats],[r["objects"] for r in frame_stats])
        axes[0].set_ylabel("Objects per frame")
        axes[1].plot([r["frame"] for r in frame_stats],[r["mean_rain_mm_h"] for r in frame_stats])
        axes[1].set(xlabel="Input frame",ylabel="Domain mean rain (mm/h)")
        fig.savefig(out / "statistics.png",dpi=150)
        plt.close(fig)
        fig,axes = plt.subplots(1,2,figsize=(12,4),layout="constrained")
        axes[0].bar([k.replace("_seconds","") for k in timed_keys],list(totals.values()))
        axes[0].tick_params(axis="x",rotation=25)
        axes[0].set_ylabel("Wall seconds")
        axes[1].plot([r["start_frame"] for r in rows],[r["parent_peak_rss_mb"] for r in rows],label="Parent high water")
        axes[1].plot([r["start_frame"] for r in rows],[r["sampled_tree_peak_rss_mb"] or np.nan for r in rows],label="Sampled tree RSS sum (if available)")
        axes[1].set(xlabel="Input frame",ylabel="Cumulative peak RSS (MiB)")
        axes[1].legend()
        fig.savefig(out / "performance.png",dpi=150)
        plt.close(fig)
        (out / "REPORT.md").write_text("# High-resolution streaming benchmark\n\n"
            "Completed chunked processing with checkpoint reloads. No whole-run equivalence comparison was performed in this mode. "
            "Run small-sample validation separately. See performance.csv, performance.png and summary.json. "
            "Preview maps are from the first chunk. Each chunk contains graph tables and state. "
            "Use finalize_catalog.py with all chunk directories to canonicalize families. "
            "Timing depends on rain-object density as well as grid size; do not extrapolate linearly to a full season. "
            "RSS includes shared pages more than once; inspect Slurm accounting for actual job usage.\n")
        (out / "SUCCESS").touch()
        print(json.dumps(summary,indent=2),flush=True)
    finally:
        stop_monitor.set()
        thread.join(timeout=2)


def main():
    a = arguments()
    if a.hours < 2 or a.crop < 0 or a.chunk_frames < 1 or a.workers < 1 or a.bridge_radius < 0:
        raise ValueError("Need >=2 frames and positive crop/chunk/workers; radius >=0")
    if a.grid_km is not None and (not np.isfinite(a.grid_km) or a.grid_km <= 0):
        raise ValueError("--grid-km must be positive")
    if a.stream and not a.inspect:
        stream_benchmark(a)
        return
    read_started = time.perf_counter()
    loaded = load_input(a)
    read_seconds = time.perf_counter()-read_started
    if a.inspect:
        return
    if a.output_dir is None:
        raise ValueError("--output-dir required")
    out = a.output_dir
    out.mkdir(parents=True, exist_ok=False)
    data, meta = loaded
    dt = meta["dt_hours"]
    (out / "metadata.json").write_text(json.dumps(meta, indent=2, default=str))
    settings = provenance(a)
    (out / "settings.json").write_text(json.dumps(settings, indent=2, default=str))
    started = time.perf_counter()
    print("Identifying", data.shape, flush=True)
    valid = np.isfinite(data)
    labels = identify(data, disk(a.bridge_radius), workers=a.workers, threshold=a.threshold,
                      min_size=a.min_size, valid_mask=valid)
    identification_seconds = time.perf_counter()-started
    kwargs = dict(tau=a.tau, km=a.max_displacement, max_gap=a.max_gap, gap_tau=a.gap_tau,
                  gap_ambiguity=a.gap_ambiguity, event_overlap=a.event_overlap, sequence_id=a.sequence_id)
    print("Tracking whole sample", flush=True)
    track_started = time.perf_counter()
    tracked, graph, state = track_with_graph(labels, data, return_state=True, **kwargs)
    tracking_seconds = time.perf_counter()-track_started
    elapsed = time.perf_counter()-started
    write_graph(graph, out)
    np.save(out / "identified_labels.npy", labels)
    np.save(out / "tracked_labels.npy", tracked)
    save_tracking_state(state, out / "state.json")
    print("Checking independent chunk identification and JSON checkpoint resume", flush=True)
    comparison_started = time.perf_counter()
    combined = type(graph)()
    resumed = None
    raster_equal = True
    id_equal = True
    chunk_size = min(a.chunk_frames, max(1,len(data)//2))
    for start in range(0,len(data),chunk_size):
        stop = min(len(data),start+chunk_size)
        local = identify(data[start:stop], disk(a.bridge_radius), workers=a.workers,
                         threshold=a.threshold, min_size=a.min_size, valid_mask=valid[start:stop])
        id_equal &= np.array_equal(local, labels[start:stop])
        raster, part, resumed = track_with_graph(local,data[start:stop],state=resumed,return_state=True,**kwargs)
        raster_equal &= np.array_equal(raster, tracked[start:stop])
        combined.objects.extend(part.objects)
        combined.edges.extend(part.edges)
        combined.events.extend(part.events)
        combined.family_map = part.family_map
        save_tracking_state(resumed, out / "chunk_state.json")
        resumed = load_tracking_state(out / "chunk_state.json")
    lookup = {n.node_id:n for n in graph.objects}
    comparison_seconds = time.perf_counter()-comparison_started
    checks = {"identification_chunk_equal": bool(id_equal), "branch_raster_chunk_equal": bool(raster_equal),
              "graph_chunk_equal": signature(graph)==signature(combined),
              "unique_node_ids": len(lookup)==len(graph.objects),
              "unique_branch_per_frame": len({(n.time,n.branch_id) for n in graph.objects})==len(graph.objects),
              "edge_references_and_forward_time": all(e.parent_id in lookup and e.child_id in lookup and
                  lookup[e.parent_id].time < lookup[e.child_id].time for e in graph.edges)}
    branches = {}
    for n in graph.objects:
        branches.setdefault(n.branch_id, []).append(n)
    spans = [(ns[-1].time-ns[0].time+1)*dt for ns in branches.values()]
    speed = [float(np.linalg.norm(np.subtract(lookup[e.child_id].object.centroid,lookup[e.parent_id].object.centroid))) /
             ((lookup[e.child_id].time-lookup[e.parent_id].time)*dt)
             for e in graph.edges if lookup[e.child_id].branch_id==lookup[e.parent_id].branch_id]
    families = Counter(graph.family_map[n.branch_id] for n in graph.objects)
    summary = {"checks": checks, "shape": list(data.shape), "runtime_identification_tracking_seconds": elapsed,
               "read_seconds": read_seconds,"identification_seconds": identification_seconds,
               "tracking_seconds": tracking_seconds,"chunk_comparison_seconds": comparison_seconds,
               "parent_peak_rss_mb": peak_rss_mb(),
               "nodes": len(graph.objects), "branches": len(branches), "families": len(families),
               "events": dict(Counter(e.event for e in graph.events)), "gap_edges": sum(e.gap>0 for e in graph.edges),
               "edges": len(graph.edges), "largest_family_nodes": max(families.values(),default=0),
               "missing_pixel_fraction": float(1-valid.mean()), "wet_pixel_fraction": float(np.mean(labels>0)),
               "branch_span_hours_median": float(np.median(spans)) if spans else None,
               "branch_span_hours_max": max(spans,default=None),
               "speed_cells_per_hour_median": float(np.median(speed)) if speed else None,
               "time_or_crop_boundary_branches": sum(any(n.time in (0,len(data)-1) or
                   n.object.bbox[0][0]==0 or n.object.bbox[1][0]==0 or
                   n.object.bbox[0][1]==data.shape[1] or n.object.bbox[1][1]==data.shape[2] for n in ns) for ns in branches.values())}
    if a.grid_km:
        summary["speed_km_per_hour_median"] = float(np.median(speed))*a.grid_km if speed else None
    (out / "summary.json").write_text(json.dumps(summary,indent=2))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2,2,figsize=(11,8),layout="constrained")
    counts = Counter(n.time for n in graph.objects)
    axes[0,0].plot(np.arange(len(data))*dt,[counts[t] for t in range(len(data))])
    axes[0,0].set(xlabel="Elapsed hours",ylabel="Objects per frame")
    axes[0,1].hist(spans,bins=20)
    axes[0,1].set(xlabel="Observed branch span (hours; gaps included)",ylabel="Branches")
    axes[1,0].hist([n.object.area*(a.grid_km**2 if a.grid_km else 1) for n in graph.objects],bins=25)
    axes[1,0].set(xlabel="Object area (km2)" if a.grid_km else "Object area (cells)",ylabel="Object snapshots")
    axes[1,1].hist(np.array(speed)*(a.grid_km or 1),bins=25)
    axes[1,1].set(xlabel="Within-branch speed (km/h)" if a.grid_km else "Within-branch speed (cells/h)",ylabel="Links")
    fig.savefig(out / "statistics.png",dpi=150)
    plt.close(fig)
    print("Rendering sample and event frames", flush=True)
    plots(data,labels,tracked,graph,dt,out,a.gif)
    (out / "REPORT.md").write_text("# STEP real-data validation\n\n"+
        ("PASS" if all(checks.values()) else "FAIL")+": computational checks only.\n\n"+
        "See summary.json, statistics.png, frame_*.png and the node/edge/event tables. "
        "All rainfall was normalized to mm/h before identification and tracking. "
        "Object precipitation_sum is a sum over cells, not a water volume. "
        "Branch spans include gaps and are censored by sample/crop boundaries and split/merge events; "
        "they are not complete storm lifetimes. Review visual matching and event plausibility before production. "
        "This run does not compare to legacy STEP or calibrate thresholds. "
        "Array coordinates are cropped grid cells, not longitude/latitude.\n")
    print(json.dumps(summary,indent=2), flush=True)
    if not all(checks.values()):
        raise SystemExit("FAILED computational validation; inspect summary.json")
    (out / "SUCCESS").touch()
    print("Completed:",out,flush=True)


if __name__ == "__main__":
    main()
