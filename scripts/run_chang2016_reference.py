#!/usr/bin/env python3
"""Run the Chang et al. (2016) S1 reference algorithm on a rain-rate array.

Input: a NumPy file (time, y, x) of rain rates in mm/h on a uniform grid.
Hourly input is averaged to 3-hourly rates by the paper presets. Output:
segment and event label arrays, a segment table, and the paper's
four-factor decomposition. Existing STEP_UPGRADE outputs are not modified.
"""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from step import chang2016 as c


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("rain", type=Path, help=".npy (time, y, x) rain rates in mm/h")
    p.add_argument("--preset", type=Path, required=True)
    p.add_argument("--grid-km", type=float, required=True,
                   help="uniform grid spacing in km (paper: 12 model, 4 Stage IV)")
    p.add_argument("--input-dt-hours", type=float, default=1.0)
    p.add_argument("--cell-area-km2", type=float,
                   help="default: grid_km squared (the paper's convention)")
    p.add_argument("--fft-workers", type=int, default=1)
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args(argv)
    a.output_dir.mkdir(parents=True, exist_ok=False)
    preset = json.loads(a.preset.read_text())
    ident, trk = preset["identification"], preset["tracking"]
    rain = np.load(a.rain, mmap_mode="r")
    factor = int(preset["input_contract"].get("aggregate_hourly_input", 1))
    dt = a.input_dt_hours
    if factor > 1:
        if a.input_dt_hours != 1:
            raise SystemExit("aggregation presets expect hourly input")
        rain = c.aggregate_rate(rain, factor)
        dt = float(factor)
    expected = float(preset["input_contract"]["dt_hours"])
    if dt != expected:
        raise SystemExit(f"preset expects dt={expected} h after aggregation, got {dt} h")
    start = time.perf_counter()
    labels = c.identify(np.asarray(rain), radius_km=ident["radius_km"], grid_km=a.grid_km,
                        threshold=ident["threshold_mm_h"], strict=ident["strict"],
                        erosion_border=ident["erosion_border"])
    t_id = time.perf_counter() - start
    events, rows, tracker = c.track(
        labels, np.asarray(rain), grid_km=a.grid_km, tau=trk["tau"],
        phi_per_km=trk["phi_per_km"], max_displacement_km=trk["max_displacement_km"],
        max_angle_deg=trk["max_angle_deg"], no_history=trk["no_history"],
        fft_workers=a.fft_workers)
    t_track = time.perf_counter() - start - t_id
    np.save(a.output_dir / "segment_labels.npy", labels)
    np.save(a.output_dir / "event_labels.npy", events)
    with (a.output_dir / "segments.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["time"])
        writer.writeheader()
        writer.writerows(rows)
    area = a.cell_area_km2 if a.cell_area_km2 else a.grid_km ** 2
    factors = c.factorize(rows, cell_area_km2=area, dt_hours=dt, n_time_steps=len(labels))
    summary = dict(preset=preset, grid_km=a.grid_km, cell_area_km2=area, dt_hours=dt,
                   frames=int(len(labels)), factors=factors,
                   units=dict(intensity="mm/h", size="km2", duration="h",
                              total="mm km2 (1 mm km2 = 1000 m3)"),
                   seconds=dict(identification=t_id, tracking=t_track),
                   pairs_evaluated=tracker.pairs_evaluated, pairs_pruned=tracker.pairs_pruned)
    (a.output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["factors"], indent=2))


if __name__ == "__main__":
    main()
