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

import numpy as np

from step.identification import identify
from step.tracking import track_with_graph


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
    parser.add_argument("--bridge-radius", type=int, default=1, help="Maximum gap to bridge, in grid cells.")
    parser.add_argument("--min-size", type=int, default=1, help="Minimum rain-object area, in cells.")
    parser.add_argument("--workers", type=int, default=0, help="0 uses SLURM_CPUS_PER_TASK; otherwise explicit count.")
    parser.add_argument("--tau", type=float, default=0.35, help="Minimum sparse tracking score.")
    parser.add_argument("--max-displacement", type=float, default=20.0, help="Maximum hourly motion, in grid cells.")
    return parser.parse_args()


def write_graph(graph, output_dir):
    with (output_dir / "objects.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "storm_id", "local_label", "area_cells", "centroid_y", "centroid_x", "mean_intensity", "y_start", "y_stop", "x_start", "x_stop"])
        for time, storm_id, obj in graph.objects:
            writer.writerow([time, storm_id, obj.label, obj.area, obj.centroid[0], obj.centroid[1], obj.mean_intensity,
                             obj.bbox[0].start, obj.bbox[0].stop, obj.bbox[1].start, obj.bbox[1].stop])
    with (output_dir / "edges.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", "parent_id", "child_id", "score", "event"])
        for edge in graph.edges:
            writer.writerow([edge.time, edge.parent_id, edge.child_id, edge.score, edge.event])


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
    precip = np.asarray(source[args.start:stop], dtype=np.float32)
    precip[~np.isfinite(precip)] = 0
    precip[precip < args.threshold] = 0
    structure = disk(args.bridge_radius)

    labels = identify(precip, structure, workers=workers, min_size=args.min_size)
    tracked, graph = track_with_graph(labels, precip, tau=args.tau, km=args.max_displacement)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    np.save(output_dir / "identified_labels.npy", labels)
    np.save(output_dir / "tracked_labels.npy", tracked)
    write_graph(graph, output_dir)
    settings = vars(args).copy()
    settings.update({"resolved_workers": workers, "input_shape": list(source.shape), "processed_shape": list(precip.shape)})
    with (output_dir / "run_settings.json").open("w") as handle:
        json.dump(settings, handle, indent=2, sort_keys=True)

    print("Processed shape:", precip.shape)
    print("Rain-object snapshots:", len(graph.objects))
    print("Track edges:", len(graph.edges))
    print("Persistent storm IDs:", int(tracked.max()))
    print("Outputs:", output_dir)


if __name__ == "__main__":
    main()
