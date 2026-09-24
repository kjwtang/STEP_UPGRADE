#!/usr/bin/env python3
"""Run a reproducible STEP_UPGRADE validation on a 3-D NumPy precipitation cube.

The input must have shape (time, y, x).  This intentionally handles .npy
first; adapting the loading block to a site's NetCDF/Zarr convention is a
small, isolated follow-up once the data layout is known.
"""

import argparse
import csv
import json
import os
from pathlib import Path
import shutil

import numpy as np

from step.identification import identify
from step.tracking import (
    load_tracking_state,
    save_tracking_state,
    track_with_graph,
)


def disk(radius):
    yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return (xx * xx + yy * yy) <= radius * radius


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Path to a .npy array with shape (time, y, x).")
    parser.add_argument("--output-dir", required=True, help="Directory for validation outputs.")
    parser.add_argument("--start", type=int, default=0, help="First input time index (inclusive).")
    parser.add_argument("--stop", type=int, default=None, help="Last input time index (exclusive).")
    parser.add_argument("--threshold", type=float, default=0.0, help="Keep precipitation >= this value.")
    parser.add_argument("--bridge-radius", type=int, default=1, help="Binary dilation radius in grid cells.")
    parser.add_argument("--min-size", type=int, default=1, help="Minimum rain-object area, in cells.")
    parser.add_argument("--workers", type=int, default=0, help="0 uses SLURM_CPUS_PER_TASK; otherwise explicit count.")
    parser.add_argument("--tau", type=float, default=0.35, help="Minimum sparse tracking score.")
    parser.add_argument("--max-displacement", type=float, default=20.0, help="Maximum hourly motion, in grid cells.")
    parser.add_argument("--max-gap", type=int, default=1, help="Maximum missing object frames to reconnect.")
    parser.add_argument("--gap-tau", type=float, default=None, help="Minimum gap reconnection score.")
    parser.add_argument("--gap-ambiguity", type=float, default=0.05, help="Reject gap matches when the two best scores differ by less than this value.")
    parser.add_argument("--event-overlap", type=float, default=0.10, help="Minimum parent/child coverage for split/merge evidence.")
    parser.add_argument("--state-input", default=None, help="Optional JSON checkpoint from the preceding chunk/month.")
    parser.add_argument("--state-output", default=None, help="Optional JSON checkpoint path for the following chunk/month.")
    parser.add_argument("--absolute-start", type=int, default=None, help="Absolute integer time for the first processed frame.")
    parser.add_argument("--sequence-id", default="default", help="Climate/year/member namespace checked when resuming.")
    return parser.parse_args()


def write_graph(graph, output_dir):
    with (output_dir / "objects.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "node_id", "branch_id", "family_id", "local_label", "area_cells", "centroid_y", "centroid_x", "mean_intensity", "max_intensity", "precipitation_sum", "y_start", "y_stop", "x_start", "x_stop"])
        for node in graph.objects:
            obj = node.object
            writer.writerow([node.time, node.node_id, node.branch_id, graph.family_map[node.branch_id], obj.label, obj.area, obj.centroid[0], obj.centroid[1], obj.mean_intensity, obj.max_intensity, obj.precipitation_sum,
                             obj.bbox[0][0], obj.bbox[0][1], obj.bbox[1][0], obj.bbox[1][1]])
    with (output_dir / "edges.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "parent_node_id", "child_node_id", "score", "event", "gap", "raw_iou", "advected_iou", "distance", "parent_coverage", "child_coverage"])
        for edge in graph.edges:
            writer.writerow([edge.time, edge.parent_id, edge.child_id, edge.score, edge.event, edge.gap, edge.raw_iou, edge.advected_iou, edge.distance, edge.parent_coverage, edge.child_coverage])
    with (output_dir / "events.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["event_id", "time", "event", "parent_node_ids", "child_node_ids"])
        for event in graph.events:
            writer.writerow([event.event_id, event.time, event.event,
                             ";".join(map(str, event.parent_ids)),
                             ";".join(map(str, event.child_ids))])
    with (output_dir / "family_map.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["branch_id", "family_id"])
        writer.writerows(sorted(graph.family_map.items()))


def main():
    args = parse_args()
    workers = args.workers or int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    if workers < 1:
        raise ValueError("workers must be positive")
    if args.bridge_radius < 0:
        raise ValueError("bridge-radius cannot be negative")

    source = np.load(args.input, mmap_mode="r")
    if source.ndim != 3:
        raise ValueError("input must have shape (time, y, x)")
    stop = args.stop if args.stop is not None else source.shape[0]
    if not 0 <= args.start < stop <= source.shape[0]:
        raise ValueError("invalid --start/--stop time window")

    # Copy only the requested time window to working memory.  float32 halves
    # the memory of common float64 source files without changing object labels.
    precip = np.array(source[args.start:stop], dtype=np.float32, copy=True)
    valid = np.isfinite(precip)
    precip[~valid] = 0
    structure = disk(args.bridge_radius)

    labels = identify(
        precip, structure, workers=workers, min_size=args.min_size,
        threshold=args.threshold, valid_mask=valid,
    )
    state = load_tracking_state(args.state_input) if args.state_input else None
    tracked, graph, state = track_with_graph(
        labels, precip, tau=args.tau, km=args.max_displacement,
        max_gap=args.max_gap, gap_tau=args.gap_tau,
        gap_ambiguity=args.gap_ambiguity,
        event_overlap=args.event_overlap, state=state,
        start_time=args.absolute_start, return_state=True,
        sequence_id=args.sequence_id,
    )

    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"output directory already exists: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_dir.with_name(f".{output_dir.name}.{os.getpid()}.tmp")
    if temporary.exists():
        raise FileExistsError(f"temporary output already exists: {temporary}")
    temporary.mkdir()
    settings = vars(args).copy()
    settings.update({"resolved_workers": workers, "input_shape": list(source.shape), "processed_shape": list(precip.shape)})
    try:
        np.save(temporary / "identified_labels.npy", labels)
        np.save(temporary / "tracked_labels.npy", tracked)
        write_graph(graph, temporary)
        with (temporary / "run_settings.json").open("w") as handle:
            json.dump(settings, handle, indent=2, sort_keys=True)
        if args.state_output:
            requested_state = Path(args.state_output).resolve()
            try:
                relative_state = requested_state.relative_to(output_dir)
            except ValueError as error:
                raise ValueError(
                    "state-output must be inside output-dir so data and state publish together"
                ) from error
            save_tracking_state(state, temporary / relative_state)
        os.replace(temporary, output_dir)
    except Exception:
        if temporary.exists():
            shutil.rmtree(temporary)
        raise

    print("Processed shape:", precip.shape)
    print("Rain-object snapshots:", len(graph.objects))
    print("Track edges:", len(graph.edges))
    print("Maximum branch ID in chunk:", int(tracked.max()))
    print("Outputs:", output_dir)


if __name__ == "__main__":
    main()
