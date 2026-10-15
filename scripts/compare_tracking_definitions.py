#!/usr/bin/env python3
"""Compare storm-event definitions on one hourly rain-rate sequence.

Input: .npz (array ``rain_mm_h``) or .npy, shape (time, y, x), mm/h, uniform
grid. Runs, on identical input:

  A   paper S1 rules on 3-hourly mean rates (first multiple of 3 hours)
  B   paper S1 rules on hourly frames
  C1  hourly, selection = most_similar
  C2  hourly, most_similar, displacement threshold = speed x dt (default 40 km)
  D   STEP_UPGRADE reference preset (1 mm/h, 9-cell disk, lineage tracker)
  D2  STEP_UPGRADE lineage tracker on the paper (B) segments

Writes results.json (factors, N, sum_l, timings, link diagnostics) and the
per-segment tables. Reports rain-volume shares removed by each threshold.
Event counts in short windows are censored at both ends; compare runs
relative to each other on the same window.
"""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from step import chang2016 as c
from step import track_with_graph
from step.identification_fast import identify_disk


def load(path):
    if path.suffix == ".npz":
        with np.load(path) as data:
            return np.asarray(data["rain_mm_h"], dtype=np.float32)
    return np.load(path).astype(np.float32)


def link_diagnostics(rows):
    linked = [r for r in rows if r["rule"] != "new"]
    if not linked:
        return {}
    disp = np.array([r["displacement_km"] for r in linked])
    area = {(r["time"], r["segment_label"]): r["area_cells"] for r in rows}
    far = [r for r in linked if r["displacement_km"] >= 300]
    return dict(
        links=len(linked),
        rules={k: sum(r["rule"] == k for r in rows) for k in ("new", "distance", "angle")},
        displacement_km_p50_p90_p99_max=[*np.percentile(disp, [50, 90, 99]).round(1).tolist(),
                                         float(disp.max())],
        links_ge_300km=len(far),
        median_child_area_cells_ge_300km=float(np.median([r["area_cells"] for r in far])) if far else None,
        median_parent_area_cells_ge_300km=float(np.median(
            [area[(r["time"] - 1, r["parent_label"])] for r in far])) if far else None,
    )


def paper_run(name, labels, rain, dt, grid_km, cell_area, tau, out, **kw):
    start = time.perf_counter()
    events, rows, tracker = c.track(labels, rain, grid_km=grid_km, tau=tau, **kw)
    seconds = time.perf_counter() - start
    f = c.factorize(rows, cell_area, dt, len(labels))
    with (out / f"segments_{name}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return dict(definition="paper events", dt_hours=dt, frames=len(labels),
                track_s_per_frame=seconds / len(labels),
                N=f["number_of_events"], sum_l=len({(r["event_id"], r["time"]) for r in rows}),
                intensity_mm_h=f["intensity_mm_per_h"], size_1e4_km2=f["size_factor_km2"] / 1e4,
                duration_h=f["duration_factor_h"], options=kw,
                pairs=dict(exact=tracker.pairs_evaluated,
                           certified=tracker.pairs_certified + tracker.pairs_block_certified,
                           bound_rejected=tracker.pairs_block_rejected,
                           pruned=tracker.pairs_pruned),
                links=link_diagnostics(rows)), rows


def lineage_run(labels, rain, preset, cell_area, dt):
    p = dict(preset["tracking"])
    p["km"] = p.pop("max_displacement")
    start = time.perf_counter()
    _, graph = track_with_graph(labels.astype(np.int32), rain, **p)
    seconds = time.perf_counter() - start
    wet = sum(n.object.area for n in graph.objects)
    rate = sum(n.object.precipitation_sum for n in graph.objects)
    result = dict(track_s_per_frame=seconds / len(labels), objects=len(graph.objects),
                  topology={})
    for e in graph.events:
        result["topology"][e.event] = result["topology"].get(e.event, 0) + 1
    for key, attr in (("branch", "branch_id"), ("family", "family_id")):
        ids = {getattr(n, attr) for n in graph.objects}
        sum_l = len({(getattr(n, attr), n.time) for n in graph.objects})
        result[key] = dict(N=len(ids), sum_l=sum_l, intensity_mm_h=rate / wet,
                           size_1e4_km2=cell_area * wet / sum_l / 1e4,
                           duration_h=dt * sum_l / len(ids))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("rain", type=Path)
    p.add_argument("--grid-km", type=float, default=4.0)
    p.add_argument("--season", choices=("summer", "winter"), default="summer")
    p.add_argument("--speed-kmh", type=float, default=40.0, help="C2 displacement per hour")
    p.add_argument("--runs", default="A,B,C1,C2,D,D2")
    p.add_argument("--fft-workers", type=int, default=1)
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=False)
    runs = set(a.runs.split(","))
    root = Path(__file__).resolve().parents[1]
    preset = json.loads((root / "configs/rainfall_lineage_4km_1h_reference_v1.json").read_text())
    paper = c.PAPER_PRESETS[a.season]
    cell_area = a.grid_km ** 2
    rain = load(a.rain)
    total = float(rain.sum(dtype=np.float64))
    results = dict(input=str(a.rain), shape=list(rain.shape), season=a.season, grid_km=a.grid_km,
                   rain_share_removed=dict(
                       paper_0p033_mm_h=float(rain[rain <= 0.1 / 3].sum(dtype=np.float64) / total),
                       preset_1_mm_h=float(rain[rain < 1].sum(dtype=np.float64) / total)),
                   timings={}, runs={})
    kw = dict(fft_workers=a.fft_workers)

    t = time.perf_counter()
    seg = c.identify(rain, radius_km=paper["radius_km"], grid_km=a.grid_km)
    results["timings"]["paper_identification_s_per_frame"] = (time.perf_counter() - t) / len(rain)
    np.save(a.output_dir / "segments_paper_1h.npy", seg)

    if "A" in runs and len(rain) >= 3:
        r3 = c.aggregate_rate(rain[:len(rain) // 3 * 3], 3)
        t = time.perf_counter()
        seg3 = c.identify(r3, radius_km=paper["radius_km"], grid_km=a.grid_km)
        results["timings"]["paper_identification_3h_s_per_frame"] = (time.perf_counter() - t) / len(r3)
        results["runs"]["A"], _ = paper_run("A", seg3, r3, 3, a.grid_km, cell_area, paper["tau"],
                                            a.output_dir, **kw)
    specs = dict(B=dict(), C1=dict(selection="most_similar"),
                 C2=dict(selection="most_similar", max_displacement_km=a.speed_kmh))
    for name, opts in specs.items():
        if name in runs:
            results["runs"][name], rows = paper_run(name, seg, rain, 1, a.grid_km, cell_area,
                                                    paper["tau"], a.output_dir, **opts, **kw)
            if name == "B":
                rows_b = rows
            if name == "C1" and "B" in runs:
                pb = {(r["time"], r["segment_label"]): r["parent_label"] for r in rows_b}
                pc = {(r["time"], r["segment_label"]): r["parent_label"] for r in rows}
                both = [k for k in pb if pb[k] and pc[k]]
                results["B_vs_C1_parent_disagreement"] = dict(
                    linked_in_both=len(both), different_parent=sum(pb[k] != pc[k] for k in both))
    if "D" in runs:
        ident = preset["identification"]
        t = time.perf_counter()
        seg_u = identify_disk(rain, ident["bridge_radius_cells"], threshold=ident["threshold_mm_h"],
                              min_size=ident["min_size_cells"])
        results["timings"]["upgrade_identification_fast_s_per_frame"] = (time.perf_counter() - t) / len(rain)
        results["runs"]["D"] = lineage_run(seg_u, rain, preset, cell_area, 1)
    if "D2" in runs:
        results["runs"]["D2"] = lineage_run(seg, rain, preset, cell_area, 1)
    (a.output_dir / "results.json").write_text(json.dumps(results, indent=2, default=float))
    print(json.dumps(results, indent=2, default=float))


if __name__ == "__main__":
    main()
