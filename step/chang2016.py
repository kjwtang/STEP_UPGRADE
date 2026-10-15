"""Reference implementation of the Chang et al. (2016) S1 algorithm.

This module is additive: it does not change ``step.identification`` or
``step.tracking``. It follows Supplemental Material S1 of
Chang, Stein, Wang, Kotamarthi and Moyer (2016), J. Climate 29, 8355-8376,
doi:10.1175/JCLI-D-15-0844.1, with physical units made explicit.

Paper settings (S1.1, S1.2):

* structuring disks: radius 60 km (summer) and 120 km (winter);
* similarity sim = sum_s sum_u w_i(s) w_j(u) exp(-phi d(s, u)), d in km,
  phi = 1/120 km^-1, w = rain / sum(rain) within each segment;
* tau = 0.05 (summer) and 0.01 (winter); a link needs sim > tau;
* displacement between intensity-weighted centroids < 120 km, or angle
  between the new and the previous displacement < 120 degrees;
* among qualifying previous segments, choose the one with most grid cells;
* several current segments may share a predecessor (splits share an event);
  a current segment links to exactly one predecessor (merges continue one).

The paper used 3-hourly data and a cutoff of > 0.1 mm per 3 hours. Rain
intensity is mm/hour and size is km^2; the two are never combined into a
mm/hour/km^2 quantity. Choices the supplement does not specify are explicit
keyword arguments and are documented in docs/CHANG2016_REFERENCE.md.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import ndimage
from scipy import fft as sp_fft
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

PAPER_THRESHOLD_MM_PER_3H = 0.1
PAPER_PRESETS = {
    "summer": dict(radius_km=60.0, tau=0.05),
    "winter": dict(radius_km=120.0, tau=0.01),
}
PAPER_PHI_PER_KM = 1.0 / 120.0
PAPER_MAX_DISPLACEMENT_KM = 120.0
PAPER_MAX_ANGLE_DEG = 120.0

_EIGHT = np.ones((3, 3), dtype=bool)


# ---------------------------------------------------------------------------
# Morphology on binary masks (exact disk operations via Euclidean distance).
# ---------------------------------------------------------------------------

def disk(radius_cells: float) -> np.ndarray:
    """Disk structuring element {(y, x): y^2 + x^2 <= r^2}."""
    r = float(radius_cells)
    n = int(math.floor(r))
    yy, xx = np.ogrid[-n:n + 1, -n:n + 1]
    return (yy * yy + xx * xx) <= r * r + 1e-9


def _squared_distance(mask_false_is_target: np.ndarray) -> np.ndarray:
    """Exact squared Euclidean distance (cells) to the nearest False cell."""
    d = ndimage.distance_transform_edt(mask_false_is_target)
    return np.rint(d * d)


def dilate_disk(mask: np.ndarray, radius_cells: float) -> np.ndarray:
    """Binary dilation by ``disk(radius_cells)``; cost independent of radius."""
    mask = np.asarray(mask, dtype=bool)
    if radius_cells <= 0 or not mask.any():
        return mask.copy()
    return _squared_distance(~mask) <= radius_cells ** 2 + 1e-9


def erode_disk(mask: np.ndarray, radius_cells: float,
               border: str = "ignore") -> np.ndarray:
    """Binary erosion by ``disk(radius_cells)``.

    ``border='ignore'`` treats cells outside the domain as rain, which is the
    behaviour of skimage.morphology.erosion used by the original Python STEP.
    ``border='background'`` treats them as dry (scipy binary_erosion default).
    """
    mask = np.asarray(mask, dtype=bool)
    if radius_cells <= 0:
        return mask.copy()
    if not mask.any():
        return np.zeros_like(mask)
    if border not in ("ignore", "background"):
        raise ValueError("border must be 'ignore' or 'background'")
    pad = int(math.ceil(radius_cells)) + 1
    padded = np.pad(mask, pad, constant_values=(border == "ignore"))
    d2 = _squared_distance(padded)[pad:-pad, pad:-pad]
    return d2 > radius_cells ** 2 + 1e-9


def _stable_relabel(labels: np.ndarray) -> np.ndarray:
    flat = labels.ravel()
    positions = np.flatnonzero(flat > 0)
    if not len(positions):
        return np.zeros(labels.shape, dtype=np.int32)
    present, first = np.unique(flat[positions], return_index=True)
    ordered = present[np.argsort(first)]
    lookup = np.zeros(int(labels.max()) + 1, dtype=np.int32)
    lookup[ordered] = np.arange(1, len(ordered) + 1, dtype=np.int32)
    return lookup[labels]


# ---------------------------------------------------------------------------
# S1.1 identification.
# ---------------------------------------------------------------------------

def identify_frame(rain: np.ndarray, radius_km: float, grid_km: float,
                   threshold: float = PAPER_THRESHOLD_MM_PER_3H / 3.0,
                   strict: bool = True, valid_mask: Optional[np.ndarray] = None,
                   erosion_border: str = "ignore") -> np.ndarray:
    """Identify rainstorm segments C_t,1..m at one time step (S1.1).

    ``rain`` is a rate (normally mm/h). ``threshold`` uses the same unit; the
    default is the paper's 0.1 mm/3 h expressed as a mean rate. ``strict``
    selects ``rain > threshold`` (paper wording) instead of ``>=``.
    Returns int32 labels; 0 is background. Labels are local to the frame.
    """
    rain = np.asarray(rain, dtype=float)
    if rain.ndim != 2:
        raise ValueError("rain must be two-dimensional")
    if grid_km <= 0 or radius_km < 0:
        raise ValueError("grid_km must be positive and radius_km non-negative")
    valid = np.isfinite(rain)
    if valid_mask is not None:
        valid &= np.asarray(valid_mask, dtype=bool)
    wet = valid & ((rain > threshold) if strict else (rain >= threshold)) & (rain > 0)
    if not wet.any():
        return np.zeros(rain.shape, dtype=np.int32)
    r = radius_km / grid_km

    # Step 1: contiguous regions B_t,1..r_t (8-connectivity).
    regions, n_regions = ndimage.label(wet, structure=_EIGHT)

    # Large regions keep >= 1 cell after erosion. Erosion of the union equals
    # the union of per-region erosions because a disk is connected.
    eroded = erode_disk(wet, r, erosion_border)
    is_large = np.zeros(n_regions + 1, dtype=bool)
    is_large[np.unique(regions[eroded])] = True
    is_large[0] = False

    group = np.zeros(n_regions + 1, dtype=np.int64)  # region -> segment id
    next_group = 1

    # Step 2: opening O = dilate(erode(B_large)), then almost-connected
    # labelling of O (dilate and connect). Large regions are grouped through
    # their opened cells; a region whose opening falls in several clusters
    # joins them, because the region itself is contiguous.
    if is_large.any():
        opened = dilate_disk(eroded, r) & wet
        almost = dilate_disk(opened, r) & valid
        clusters, n_clusters = ndimage.label(almost, structure=_EIGHT)
        reg = regions[opened]
        clu = clusters[opened]
        large_ids = np.flatnonzero(is_large)
        index = np.zeros(n_regions + 1, dtype=np.int64)
        index[large_ids] = np.arange(len(large_ids))
        pairs = np.unique(np.stack([reg, clu], axis=1), axis=0)
        n_nodes = len(large_ids) + n_clusters
        graph = coo_matrix((np.ones(len(pairs)),
                            (index[pairs[:, 0]], len(large_ids) + pairs[:, 1] - 1)),
                           shape=(n_nodes, n_nodes))
        _, component = connected_components(graph, directed=False)
        roots = component[:len(large_ids)]
        _, compact = np.unique(roots, return_inverse=True)
        group[large_ids] = compact + 1
        next_group = int(compact.max()) + 2

    # Step 3: each small region is dilated individually; it joins the large
    # segment sharing the most cells with its dilation (ties: lowest id).
    large_raster = group[regions] * is_large[regions]
    small_ids = np.flatnonzero(~is_large)
    small_ids = small_ids[small_ids > 0]
    unassigned = []
    if len(small_ids):
        slices = ndimage.find_objects(regions)
        pad = int(math.floor(r))
        H, W = wet.shape
        for rid in small_ids:
            sl = slices[rid - 1]
            y0, y1 = max(0, sl[0].start - pad), min(H, sl[0].stop + pad)
            x0, x1 = max(0, sl[1].start - pad), min(W, sl[1].stop + pad)
            local = regions[y0:y1, x0:x1] == rid
            hit = large_raster[y0:y1, x0:x1][dilate_disk(local, r)]
            hit = hit[hit > 0]
            if len(hit):
                counts = np.bincount(hit)
                group[rid] = int(np.argmax(counts))  # argmax = lowest id on ties
            else:
                unassigned.append(rid)

    # Step 4: almost-connected labelling of the remaining small regions.
    if unassigned:
        unassigned = np.asarray(unassigned)
        remaining = np.isin(regions, unassigned)
        almost = dilate_disk(remaining, r) & valid
        clusters, _ = ndimage.label(almost, structure=_EIGHT)
        slices = ndimage.find_objects(regions)
        for rid in unassigned:
            sl = slices[rid - 1]
            cid = clusters[sl][regions[sl] == rid][0]
            group[rid] = next_group + int(cid) - 1

    labels = group[regions]
    labels[~wet] = 0
    return _stable_relabel(labels)


def identify(rain: np.ndarray, radius_km: float, grid_km: float, **kwargs) -> np.ndarray:
    """Apply :func:`identify_frame` to every frame of a (time, y, x) array."""
    rain = np.asarray(rain)
    if rain.ndim != 3:
        raise ValueError("rain must have shape (time, y, x)")
    valid = kwargs.pop("valid_mask", None)
    return np.stack([
        identify_frame(rain[t], radius_km, grid_km,
                       valid_mask=None if valid is None else valid[t], **kwargs)
        for t in range(rain.shape[0])
    ])


# ---------------------------------------------------------------------------
# S1.2 similarity.
# ---------------------------------------------------------------------------

@dataclass
class Segment:
    label: int
    ys: np.ndarray
    xs: np.ndarray
    weights: np.ndarray
    centroid: Tuple[float, float]   # intensity-weighted, grid cells
    area: int
    rain_sum: float
    bbox: Tuple[int, int, int, int]  # ymin, ymax, xmin, xmax (inclusive)
    event_id: int = 0
    previous_centroid: Optional[Tuple[float, float]] = None
    spread: float = 0.0  # rain-weighted mean distance to centroid, cells
    _blocks: Optional[tuple] = field(default=None, repr=False, compare=False)


def segments_from_labels(labels: np.ndarray, rain: np.ndarray) -> List[Segment]:
    labels, rain = np.asarray(labels), np.asarray(rain, dtype=float)
    out = []
    for label, sl in enumerate(ndimage.find_objects(labels), start=1):
        if sl is None:
            continue
        local = labels[sl] == label
        yy, xx = np.nonzero(local)
        ys, xs = yy + sl[0].start, xx + sl[1].start
        values = rain[ys, xs]
        total = float(values.sum())
        if total > 0:
            weights = values / total
        else:  # cannot occur for thresholded positive rain; kept for safety
            weights = np.full(len(values), 1.0 / len(values))
        centroid = (float(np.dot(ys, weights)), float(np.dot(xs, weights)))
        spread = float(np.dot(weights, np.hypot(ys - centroid[0], xs - centroid[1])))
        out.append(Segment(label, ys, xs, weights, centroid, len(ys), total,
                           (int(ys.min()), int(ys.max()), int(xs.min()), int(xs.max())),
                           spread=spread))
    return out


def similarity_bruteforce(a: Segment, b: Segment, phi_per_cell: float) -> float:
    """Direct O(n_a n_b) evaluation of the S1.2 similarity (for testing)."""
    dy = a.ys[:, None] - b.ys[None, :]
    dx = a.xs[:, None] - b.xs[None, :]
    k = np.exp(-phi_per_cell * np.sqrt(dy * dy + dx * dx))
    return float(a.weights @ k @ b.weights)


def _direct_similarity(a: Segment, b: Segment, phi_per_cell: float,
                       block: int = 2048) -> float:
    total = 0.0
    for start in range(0, a.area, block):
        ys, xs = a.ys[start:start + block], a.xs[start:start + block]
        dy = ys[:, None] - b.ys[None, :]
        dx = xs[:, None] - b.xs[None, :]
        total += float(a.weights[start:start + block]
                       @ np.exp(-phi_per_cell * np.sqrt(dy * dy + dx * dx)) @ b.weights)
    return total


_KERNEL_CACHE: Dict[Tuple[float, int, int], Tuple[np.ndarray, Tuple[int, int]]] = {}
_KERNEL_CACHE_LIMIT = 64


def _bucket(n: int) -> int:
    """Round a domain length up to a small set of sizes (<= ~19% larger)."""
    size = 16
    while size < n:
        size = int(math.ceil(size * 1.19))
    return size


def _kernel_spectrum(phi_per_cell: float, Hb: int, Wb: int, workers: int):
    key = (phi_per_cell, Hb, Wb)
    if key not in _KERNEL_CACHE:
        if len(_KERNEL_CACHE) >= _KERNEL_CACHE_LIMIT:
            _KERNEL_CACHE.pop(next(iter(_KERNEL_CACHE)))
        shape = (sp_fft.next_fast_len(3 * Hb - 2, True),
                 sp_fft.next_fast_len(3 * Wb - 2, True))
        oy, ox = np.ogrid[-(Hb - 1):Hb, -(Wb - 1):Wb]
        kernel = np.exp(-phi_per_cell * np.sqrt(oy * oy + ox * ox))
        _KERNEL_CACHE[key] = (sp_fft.rfft2(kernel, shape, workers=workers), shape)
    return _KERNEL_CACHE[key]


def _smoothed_field(b: Segment, box: Tuple[int, int, int, int],
                    phi_per_cell: float, workers: int = 1) -> np.ndarray:
    """(K * w_b) on ``box`` with the full, untruncated kernel (exact).

    A linear (not circular) convolution: the FFT length is at least
    3*bucket-2, so no wrap-around occurs. Kernel spectra are cached by
    bucketed domain shape.
    """
    y0, y1, x0, x1 = box
    H, W = y1 - y0 + 1, x1 - x0 + 1
    Hb, Wb = _bucket(H), _bucket(W)
    spectrum, shape = _kernel_spectrum(phi_per_cell, Hb, Wb, workers)
    field_b = np.zeros((H, W))
    field_b[b.ys - y0, b.xs - x0] = b.weights
    full = sp_fft.irfft2(sp_fft.rfft2(field_b, shape, workers=workers) * spectrum,
                         shape, workers=workers)
    return full[Hb - 1:Hb - 1 + H, Wb - 1:Wb - 1 + W]


def _union_box(segments) -> Tuple[int, int, int, int]:
    return (min(s.bbox[0] for s in segments), max(s.bbox[1] for s in segments),
            min(s.bbox[2] for s in segments), max(s.bbox[3] for s in segments))


def similarity(a: Segment, b: Segment, phi_per_cell: float,
               direct_limit: int = 250_000) -> float:
    """Exact S1.2 similarity without forming the n_a x n_b matrix in memory.

    Small pairs use blocked direct sums. Large pairs use one linear FFT
    convolution of b's weight field with the full (untruncated) kernel over
    the union bounding box, which is exact up to floating-point rounding.
    """
    if a.area * b.area <= direct_limit:
        return _direct_similarity(a, b, phi_per_cell)
    box = _union_box((a, b))
    smoothed = _smoothed_field(b, box, phi_per_cell)
    return float(np.dot(a.weights, smoothed[a.ys - box[0], a.xs - box[2]]))


def similarities_to_parent(parent: Segment, children: Sequence[Segment],
                           phi_per_cell: float, direct_limit: int = 250_000,
                           workers: int = 1) -> List[float]:
    """Similarity of several children to one parent, sharing one FFT field."""
    if sum(ch.area for ch in children) * parent.area <= direct_limit:
        return [_direct_similarity(ch, parent, phi_per_cell) for ch in children]
    box = _union_box((parent, *children))
    smoothed = _smoothed_field(parent, box, phi_per_cell, workers)
    return [float(np.dot(ch.weights, smoothed[ch.ys - box[0], ch.xs - box[2]]))
            for ch in children]


def _block_summary(seg: Segment, k: int):
    """Per k x k block: weight, weighted centroid, spread, occupied extent."""
    cache = seg._blocks if seg._blocks is not None else {}
    if k in cache:
        return cache[k]
    by, bx = seg.ys // k, seg.xs // k
    key = by * (int(bx.max()) + 1) + bx
    _, inv = np.unique(key, return_inverse=True)
    n = int(inv.max()) + 1
    w = np.bincount(inv, seg.weights, n)
    cy = np.bincount(inv, seg.weights * seg.ys, n) / w
    cx = np.bincount(inv, seg.weights * seg.xs, n) / w
    r = np.bincount(inv, seg.weights * np.hypot(seg.ys - cy[inv], seg.xs - cx[inv]), n) / w
    ext = np.empty((n, 4))
    ext[:, 0] = np.full(n, np.inf); ext[:, 2] = np.full(n, np.inf)
    ext[:, 1] = -np.inf; ext[:, 3] = -np.inf
    np.minimum.at(ext[:, 0], inv, seg.ys); np.maximum.at(ext[:, 1], inv, seg.ys)
    np.minimum.at(ext[:, 2], inv, seg.xs); np.maximum.at(ext[:, 3], inv, seg.xs)
    summary = (w, cy, cx, r, ext)
    cache[k] = summary
    seg._blocks = cache
    return summary


def similarity_bounds(a: Segment, b: Segment, phi_per_cell: float, k: int = 8,
                      chunk: int = 4096) -> Tuple[float, float]:
    """Certified [lower, upper] bounds on the S1.2 similarity from k x k blocks.

    Lower: Jensen within each block pair plus the triangle inequality.
    Upper: the minimum distance between occupied block extents.
    """
    wa, ya, xa, ra, ea = _block_summary(a, k)
    wb, yb, xb, rb, eb = _block_summary(b, k)
    lower = upper = 0.0
    for lo in range(0, len(wa), chunk):
        sl = slice(lo, lo + chunk)
        dc = np.hypot(ya[sl, None] - yb[None, :], xa[sl, None] - xb[None, :])
        ww = wa[sl, None] * wb[None, :]
        lower += float(np.sum(ww * np.exp(-phi_per_cell * (dc + ra[sl, None] + rb[None, :]))))
        dy = np.maximum(0, np.maximum(ea[sl, None, 0] - eb[None, :, 1], eb[None, :, 0] - ea[sl, None, 1]))
        dx = np.maximum(0, np.maximum(ea[sl, None, 2] - eb[None, :, 3], eb[None, :, 2] - ea[sl, None, 3]))
        upper += float(np.sum(ww * np.exp(-phi_per_cell * np.hypot(dy, dx))))
    return lower, upper


def _bbox_gap_cells(a: Segment, b: Segment) -> float:
    dy = max(0, a.bbox[0] - b.bbox[1], b.bbox[0] - a.bbox[1])
    dx = max(0, a.bbox[2] - b.bbox[3], b.bbox[2] - a.bbox[3])
    return math.hypot(dy, dx)


def _angle_deg(u: np.ndarray, v: np.ndarray) -> Optional[float]:
    nu, nv = float(np.linalg.norm(u)), float(np.linalg.norm(v))
    if nu == 0 or nv == 0:
        return None
    c = float(np.dot(u, v) / (nu * nv))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


# ---------------------------------------------------------------------------
# S1.2 tracking.
# ---------------------------------------------------------------------------

@dataclass
class Link:
    time: int
    label: int
    event_id: int
    parent_label: int          # 0 for a new event
    similarity: float
    displacement_km: float
    angle_deg: Optional[float]
    rule: str                  # 'distance', 'angle' or 'new'


@dataclass
class Chang2016Tracker:
    """Sequential S1.2 tracker; call :meth:`update` once per time step.

    ``no_history``: when the displacement is >= max_displacement_km and the
    previous segment has no predecessor, the angle cannot be evaluated. The
    supplement does not specify this case; 'reject' (default) refuses the
    link, 'accept' reproduces the original Python STEP behaviour at t=2.

    ``selection``: S1.2 chooses the qualifying previous segment with the most
    grid cells ('largest', default). The main text describes the "most
    morphologically similar" event; 'most_similar' is the sensitivity option.
    """

    grid_km: float
    tau: float
    phi_per_km: float = PAPER_PHI_PER_KM
    max_displacement_km: float = PAPER_MAX_DISPLACEMENT_KM
    max_angle_deg: float = PAPER_MAX_ANGLE_DEG
    no_history: str = "reject"
    selection: str = "largest"
    previous: List[Segment] = field(default_factory=list)
    next_event_id: int = 1
    last_time: Optional[int] = None
    pairs_evaluated: int = 0
    pairs_pruned: int = 0
    direct_limit: int = 250_000
    fft_workers: int = 1
    exact_similarity: bool = False
    block_levels: Tuple[int, ...] = (32, 8)
    pairs_certified: int = 0
    pairs_block_certified: int = 0
    pairs_block_rejected: int = 0

    def __post_init__(self):
        if self.grid_km <= 0 or self.phi_per_km <= 0:
            raise ValueError("grid_km and phi_per_km must be positive")
        if self.no_history not in ("reject", "accept"):
            raise ValueError("no_history must be 'reject' or 'accept'")
        if self.selection not in ("largest", "most_similar"):
            raise ValueError("selection must be 'largest' or 'most_similar'")

    @property
    def phi_per_cell(self) -> float:
        return self.phi_per_km * self.grid_km

    def _prune_gap_cells(self) -> float:
        # sim <= exp(-phi * d_min) <= exp(-phi * bbox_gap); sim > tau is then
        # impossible when bbox_gap >= ln(1/tau) / phi. Exact, not heuristic.
        if self.tau <= 0:
            return math.inf
        if self.tau >= 1:
            return 0.0
        return math.log(1.0 / self.tau) / self.phi_per_cell

    def _motion(self, seg: Segment, parent: Segment):
        """S1.2 motion rule; returns (rule, displacement_km, angle) or None."""
        e = np.subtract(seg.centroid, parent.centroid) * self.grid_km
        disp = float(np.hypot(e[0], e[1]))
        if disp < self.max_displacement_km:
            return "distance", disp, None
        if parent.previous_centroid is None:
            return ("angle", disp, None) if self.no_history == "accept" else None
        e_star = np.subtract(parent.centroid, parent.previous_centroid) * self.grid_km
        angle = _angle_deg(e, e_star)
        if angle is None or not angle < self.max_angle_deg:
            return None
        return "angle", disp, angle

    def _exact(self, current: List[Segment], pairs) -> Dict[Tuple[int, int], float]:
        """Exact similarities; small pairs direct, others on a shared hub field."""
        sims: Dict[Tuple[int, int], float] = {}
        hubs: Dict[Tuple[str, int], List[Tuple[int, int]]] = {}
        for ci, pi in pairs:
            seg, parent = current[ci], self.previous[pi]
            self.pairs_evaluated += 1
            if seg.area * parent.area <= self.direct_limit:
                sims[(ci, pi)] = _direct_similarity(seg, parent, self.phi_per_cell)
            elif parent.area >= seg.area:
                hubs.setdefault(("p", pi), []).append((ci, pi))
            else:
                hubs.setdefault(("c", ci), []).append((ci, pi))
        for (side, index), group in hubs.items():
            hub = self.previous[index] if side == "p" else current[index]
            partners = [current[ci] if side == "p" else self.previous[pi] for ci, pi in group]
            sims.update(zip(group, similarities_to_parent(
                hub, partners, self.phi_per_cell, direct_limit=0, workers=self.fft_workers)))
        return sims

    def update(self, time: int, labels: np.ndarray, rain: np.ndarray):
        time = int(time)
        if self.last_time is not None and time != self.last_time + 1:
            raise ValueError("time steps must be consecutive")
        current = segments_from_labels(labels, rain)
        cutoff = self._prune_gap_cells()
        phi = self.phi_per_cell
        # Certified decisions (exact, not approximate). For any pair:
        #   sim <= exp(-phi * bbox_gap)                       (upper bound)
        #   sim >= exp(-phi * (|dc| + spread_i + spread_j))   (lower bound;
        #          Jensen, exp(-phi d) is convex, plus the triangle inequality)
        # so sim > tau is decided without the kernel sum unless
        # bbox_gap < cutoff <= |dc| + spread_i + spread_j, cutoff = ln(1/tau)/phi.
        nc, npv = len(current), len(self.previous)
        qualified: Dict[int, List[Tuple[int, str, float, Optional[float], Optional[float]]]] = {}
        uncertain: List[Tuple[int, int]] = []
        motion: Dict[Tuple[int, int], Tuple[str, float, Optional[float]]] = {}
        if nc and npv:
            cb = np.array([s_.bbox for s_ in current], dtype=float)
            pb = np.array([s_.bbox for s_ in self.previous], dtype=float)
            cc = np.array([s_.centroid for s_ in current])
            pc = np.array([s_.centroid for s_ in self.previous])
            cs = np.array([s_.spread for s_ in current])
            ps = np.array([s_.spread for s_ in self.previous])
            for lo in range(0, nc, 512):
                hi = min(nc, lo + 512)
                dy = np.maximum(0, np.maximum(cb[lo:hi, None, 0] - pb[None, :, 1],
                                              pb[None, :, 0] - cb[lo:hi, None, 1]))
                dx = np.maximum(0, np.maximum(cb[lo:hi, None, 2] - pb[None, :, 3],
                                              pb[None, :, 2] - cb[lo:hi, None, 3]))
                gap = np.hypot(dy, dx)
                dc = np.hypot(cc[lo:hi, None, 0] - pc[None, :, 0],
                              cc[lo:hi, None, 1] - pc[None, :, 1])
                bound = dc + cs[lo:hi, None] + ps[None, :]
                keep = gap < cutoff
                self.pairs_pruned += int((~keep).sum())
                for i, j in zip(*np.nonzero(keep)):
                    ci, pi = lo + int(i), int(j)
                    rule = self._motion(current[ci], self.previous[pi])
                    if rule is None:
                        continue
                    motion[(ci, pi)] = rule
                    if bound[i, j] < cutoff:
                        self.pairs_certified += 1
                        qualified.setdefault(ci, []).append((pi, *rule, None))
                    else:
                        uncertain.append((ci, pi))
        for k in self.block_levels:  # coarse to fine certified bounds
            still = []
            for ci, pi in uncertain:
                lower, upper = similarity_bounds(current[ci], self.previous[pi], phi, k)
                if upper <= self.tau:
                    self.pairs_block_rejected += 1
                elif lower > self.tau:
                    self.pairs_block_certified += 1
                    qualified.setdefault(ci, []).append((pi, *motion[(ci, pi)], None))
                else:
                    still.append((ci, pi))
            uncertain = still
        sims = self._exact(current, uncertain)
        for (ci, pi), value in sims.items():
            if value > self.tau:
                qualified.setdefault(ci, []).append((pi, *motion[(ci, pi)], value))
        if self.selection == "most_similar" or self.exact_similarity:
            need = [(ci, pi) for ci, items in qualified.items()
                    for pi, _, _, _, v in items if v is None]
            exact = self._exact(current, need)
            qualified = {ci: [(pi, r, d, a, exact.get((ci, pi), v)) for pi, r, d, a, v in items]
                         for ci, items in qualified.items()}
        links: List[Link] = []
        for ci, seg in enumerate(current):
            best = None
            for pi, rule, disp, angle, sim in qualified.get(ci, ()):
                parent = self.previous[pi]
                simv = -1.0 if sim is None else sim
                key = ((parent.area, simv, -pi) if self.selection == "largest"
                       else (simv, parent.area, -pi))
                if best is None or key > best[0]:
                    best = (key, parent, rule, sim, disp, angle)
            if best is None:
                seg.event_id = self.next_event_id
                self.next_event_id += 1
                links.append(Link(time, seg.label, seg.event_id, 0, 0.0, 0.0, None, "new"))
            else:
                _, parent, rule, sim, disp, angle = best
                seg.event_id = parent.event_id
                seg.previous_centroid = parent.centroid
                links.append(Link(time, seg.label, seg.event_id, parent.label,
                                  float("nan") if sim is None else sim, disp, angle, rule))
        raster = np.zeros(np.shape(labels), dtype=np.int64)
        for seg in current:
            raster[seg.ys, seg.xs] = seg.event_id
        self.previous = current
        self.last_time = time
        return raster, links, current


def track(labels: np.ndarray, rain: np.ndarray, grid_km: float, tau: float, **kwargs):
    """Track a (time, y, x) label array; returns event raster and link rows."""
    labels, rain = np.asarray(labels), np.asarray(rain)
    if labels.shape != rain.shape or labels.ndim != 3:
        raise ValueError("labels and rain need identical (time, y, x) shapes")
    tracker = Chang2016Tracker(grid_km=grid_km, tau=tau, **kwargs)
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
                parent_label=link.parent_label, rule=link.rule,
                similarity=link.similarity, displacement_km=link.displacement_km,
                angle_deg=link.angle_deg, area_cells=s.area,
                mean_intensity=s.rain_sum / s.area, precipitation_sum=s.rain_sum,
                centroid_y=s.centroid[0], centroid_x=s.centroid[1]))
    return events, rows, tracker


# ---------------------------------------------------------------------------
# Paper Section 4b factorization.
# ---------------------------------------------------------------------------

def factorize(rows: Sequence[dict], cell_area_km2: float, dt_hours: float,
              n_time_steps: Optional[int] = None) -> Dict[str, float]:
    """Total = Intensity x Size x Duration x Number (paper Section 4b).

    ``rows`` are segment rows from :func:`track`; rain sums must be rates in
    mm/h. An event's size at a time step is the sum of its segments' cells.

    * intensity_mm_per_h = sum_{i,t} a(S_i,t) s(S_i,t) / sum s(S_i,t)
    * size_factor_km2    = cell_area_km2 * sum s(S_i,t) / sum_i l(S_i)
    * duration_factor_h  = dt_hours * sum_i l(S_i) / N
    * number_of_events   = N
    The product equals sum(rain rate * cell area * dt) in mm km^2.
    """
    cells: Dict[Tuple[int, int], int] = {}
    sums: Dict[Tuple[int, int], float] = {}
    for r in rows:
        key = (int(r["event_id"]), int(r["time"]))
        cells[key] = cells.get(key, 0) + int(r["area_cells"])
        sums[key] = sums.get(key, 0.0) + float(r["precipitation_sum"])
    if not cells:
        return dict(intensity_mm_per_h=float("nan"), size_factor_km2=float("nan"),
                    duration_factor_h=float("nan"), number_of_events=0,
                    total_mm_km2=0.0)
    lifetimes: Dict[int, int] = {}
    for event, _ in cells:
        lifetimes[event] = lifetimes.get(event, 0) + 1
    total_cells = sum(cells.values())
    total_rate = sum(sums.values())
    sum_l = sum(lifetimes.values())
    n = len(lifetimes)
    result = dict(
        intensity_mm_per_h=total_rate / total_cells,
        size_factor_km2=cell_area_km2 * total_cells / sum_l,
        duration_factor_h=dt_hours * sum_l / n,
        number_of_events=n,
        total_mm_km2=total_rate * cell_area_km2 * dt_hours,
    )
    if n_time_steps:
        result["events_per_hour"] = n / (n_time_steps * dt_hours)
    return result


def aggregate_rate(rain_mm_h: np.ndarray, factor: int) -> np.ndarray:
    """Mean rate (mm/h) over consecutive blocks of ``factor`` frames.

    Use factor=3 on hourly rates to reproduce the paper's 3-hourly sampling;
    an incomplete trailing block is rejected rather than padded.
    """
    rain = np.asarray(rain_mm_h, dtype=float)
    if rain.ndim != 3 or factor < 1:
        raise ValueError("rain must be (time, y, x) and factor >= 1")
    if rain.shape[0] % factor:
        raise ValueError("number of frames must be a multiple of factor")
    return rain.reshape(rain.shape[0] // factor, factor, *rain.shape[1:]).mean(axis=1)
