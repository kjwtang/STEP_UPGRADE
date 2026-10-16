"""Hourly-sampling variant of the Chang et al. (2016) segment tracker.

The paper's motion rule (centroid displacement below 120 km, or an angle below
120 degrees to the parent's previous displacement) was designed for 3-hourly
data. At hourly sampling it accepts links of several hundred km between a
small cell and a large system, and the centroid of a segment jumps whenever
segments are regrouped. This variant keeps the paper's event semantics (one
parent per segment, several segments may share a parent, so splits share an
event) and replaces three rules:

* Candidate parents must satisfy ``similarity > tau`` (the paper's kernel
  similarity, unchanged) and lie within ``max_distance_km`` of the child
  *mask to mask* (shortest distance between any two pixels), not centroid to
  centroid. Touching or overlapping segments have distance zero.
* Among the qualifying parents the most similar one is chosen (ties: smaller
  mask distance, then larger area), not the one with the most grid cells.
* A child without an adjacent-frame parent may link to a segment whose event
  ended one frame earlier (``allow_gap``), within ``gap_distance_factor`` times
  the distance. A gap link keeps the event id; the missing frame is not counted
  as a time step of the event.

No existing module is changed; rows are compatible with
:func:`step.chang2016.factorize`.
"""
from dataclasses import dataclass, field
import math
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

from .chang2016 import (PAPER_PHI_PER_KM, Chang2016Tracker, Segment,
                        segments_from_labels)


def mask_distance_cells(a: Segment, b: Segment, limit: float) -> float:
    """Shortest pixel-to-pixel distance in cells; ``inf`` if above ``limit``."""
    if a.area > b.area:
        a, b = b, a
    pad = int(math.ceil(limit))
    keep = ((b.ys >= a.bbox[0] - pad) & (b.ys <= a.bbox[1] + pad)
            & (b.xs >= a.bbox[2] - pad) & (b.xs <= a.bbox[3] + pad))
    if not keep.any():
        return math.inf
    tree = cKDTree(np.column_stack((a.ys, a.xs)))
    dist, _ = tree.query(np.column_stack((b.ys[keep], b.xs[keep])),
                         distance_upper_bound=limit * (1.0 + 1e-9) + 1e-9)
    return float(dist.min())


def _bbox_gap(a: Segment, b: Segment) -> float:
    dy = max(0, a.bbox[0] - b.bbox[1], b.bbox[0] - a.bbox[1])
    dx = max(0, a.bbox[2] - b.bbox[3], b.bbox[2] - a.bbox[3])
    return math.hypot(dy, dx)


@dataclass
class HourlyLink:
    time: int
    label: int
    event_id: int
    parent_label: int          # 0 for a new event
    similarity: float
    displacement_km: float     # centroid to centroid, kept for diagnostics
    mask_distance_km: float
    gap: int                   # 0 adjacent frame, 1 skipped one frame
    rule: str                  # 'overlap', 'distance', 'gap' or 'new'


@dataclass
class HourlyChang2016Tracker:
    grid_km: float
    tau: float = 0.05
    phi_per_km: float = PAPER_PHI_PER_KM
    max_distance_km: float = 40.0
    allow_gap: bool = True
    gap_distance_factor: float = 2.0
    fft_workers: int = 1
    previous: List[Segment] = field(default_factory=list)
    ended: List[Segment] = field(default_factory=list)
    next_event_id: int = 1
    last_time: Optional[int] = None

    def __post_init__(self):
        if self.grid_km <= 0 or self.phi_per_km <= 0 or self.max_distance_km < 0:
            raise ValueError("grid_km and phi_per_km must be positive, distance non-negative")
        if self.gap_distance_factor < 1:
            raise ValueError("gap_distance_factor must be at least one")
        self._sim = Chang2016Tracker(grid_km=self.grid_km, tau=self.tau,
                                     phi_per_km=self.phi_per_km,
                                     fft_workers=self.fft_workers)

    def _candidates(self, current: List[Segment], parents: List[Segment],
                    only: List[int], limit_cells: float) -> Dict[int, List[Tuple[int, float, float]]]:
        """Per child: (parent index, mask distance in cells, similarity > tau)."""
        pairs: List[Tuple[int, int]] = []
        dists: Dict[Tuple[int, int], float] = {}
        for ci in only:
            for pi, parent in enumerate(parents):
                if _bbox_gap(current[ci], parent) > limit_cells:
                    continue
                d = mask_distance_cells(current[ci], parent, limit_cells)
                if math.isfinite(d):
                    pairs.append((ci, pi))
                    dists[(ci, pi)] = d
        self._sim.previous = parents
        sims = self._sim._exact(current, pairs)
        out: Dict[int, List[Tuple[int, float, float]]] = {}
        for (ci, pi), s in sims.items():
            if s > self.tau:
                out.setdefault(ci, []).append((pi, dists[(ci, pi)], s))
        return out

    @staticmethod
    def _best(options):
        return max(options, key=lambda o: (o[2], -o[1], o[3].area, -o[0]))

    def update(self, time: int, labels: np.ndarray, rain: np.ndarray):
        time = int(time)
        if self.last_time is not None and time != self.last_time + 1:
            raise ValueError("time steps must be consecutive")
        current = segments_from_labels(labels, rain)
        limit = self.max_distance_km / self.grid_km
        adjacent = self._candidates(current, self.previous, list(range(len(current))), limit)
        chosen: Dict[int, tuple] = {}
        for ci, options in adjacent.items():
            pi, d, s = self._best([(pi, d, s, self.previous[pi]) for pi, d, s in options])[:3]
            chosen[ci] = (self.previous[pi], d, s, 0)
        if self.allow_gap and self.ended:
            rest = [ci for ci in range(len(current)) if ci not in chosen]
            gapped = self._candidates(current, self.ended, rest, limit * self.gap_distance_factor)
            for ci, options in gapped.items():
                pi, d, s = self._best([(pi, d, s, self.ended[pi]) for pi, d, s in options])[:3]
                chosen[ci] = (self.ended[pi], d, s, 1)
        links: List[HourlyLink] = []
        for ci, seg in enumerate(current):
            if ci not in chosen:
                seg.event_id = self.next_event_id
                self.next_event_id += 1
                links.append(HourlyLink(time, seg.label, seg.event_id, 0, 0.0, 0.0, 0.0, 0, "new"))
                continue
            parent, d, s, gap = chosen[ci]
            seg.event_id = parent.event_id
            seg.previous_centroid = parent.centroid
            disp = float(np.hypot(*(np.subtract(seg.centroid, parent.centroid) * self.grid_km)))
            rule = "gap" if gap else ("overlap" if d == 0 else "distance")
            links.append(HourlyLink(time, seg.label, seg.event_id, parent.label, s, disp,
                                    d * self.grid_km, gap, rule))
        raster = np.zeros(np.shape(labels), dtype=np.int64)
        for seg in current:
            raster[seg.ys, seg.xs] = seg.event_id
        alive = {s.event_id for s in current}
        self.ended = [p for p in self.previous if p.event_id not in alive]
        self.previous = current
        self.last_time = time
        return raster, links, current


def track_hourly(labels: np.ndarray, rain: np.ndarray, grid_km: float, tau: float = 0.05, **kwargs):
    """Track a (time, y, x) label array; rows match :func:`step.chang2016.track`."""
    labels, rain = np.asarray(labels), np.asarray(rain)
    if labels.shape != rain.shape or labels.ndim != 3:
        raise ValueError("labels and rain need identical (time, y, x) shapes")
    tracker = HourlyChang2016Tracker(grid_km=grid_km, tau=tau, **kwargs)
    events = np.zeros(labels.shape, dtype=np.int64)
    rows: List[dict] = []
    for t in range(labels.shape[0]):
        raster, links, segments = tracker.update(t, labels[t], rain[t])
        events[t] = raster
        by_label = {s.label: s for s in segments}
        for link in links:
            s = by_label[link.label]
            rows.append(dict(
                time=t, segment_label=link.label, event_id=link.event_id,
                parent_label=link.parent_label, rule=link.rule, gap=link.gap,
                similarity=link.similarity, displacement_km=link.displacement_km,
                mask_distance_km=link.mask_distance_km, angle_deg=None,
                area_cells=s.area, mean_intensity=s.rain_sum / s.area,
                precipitation_sum=s.rain_sum,
                centroid_y=s.centroid[0], centroid_x=s.centroid[1]))
    return events, rows, tracker
