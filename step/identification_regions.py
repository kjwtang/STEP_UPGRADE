"""Opt-in large/small region identification for STEP_UPGRADE (Chang et al. 2016 S1 style).

The default STEP_UPGRADE identification dilates the whole rain mask and labels
the result, so any chain of cells closer than the bridge radius becomes one
object. This module separates regions by size first, as in S1:

1. contiguous rain regions (8-connectivity) at ``threshold``;
2. a region is *large* if it keeps at least one cell after erosion by a disk
   of ``radius_cells`` (the core radius);
3. large regions are grouped (``large_join``):
   'almost' (S1) dilates their openings by the radius and connects;
   'region' keeps every large region separate;
4. each small region is dilated individually and joins the large segment it
   touches most; small regions never connect two large segments;
5. small regions that touch no large segment (``small_join``):
   'almost' (S1) are almost-connected among themselves;
   'region' stay separate contiguous regions.

Steps 3-5 dilate by ``bridge_radius_cells`` (default: ``radius_cells``, as in
S1). A smaller core radius lets more regions count as large, so that fewer
small regions are left for step 5.

``edge_cells`` treats a band of that width along the domain border as invalid
(no rain is labelled there and no bridge crosses it); default 0.

Parameters are in grid cells and use ``rain >= threshold`` like
``step.identification.identify``. The defaults reproduce
``step.chang2016.identify_frame`` with ``strict=False`` exactly. No existing
module is changed.
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .chang2016 import _EIGHT, _stable_relabel, dilate_disk, erode_disk


def edge_mask(shape, edge_cells: int) -> np.ndarray:
    """True inside the domain, False within ``edge_cells`` of its border.

    WRF output near the lateral boundaries is affected by the relaxation zone;
    the mask can also be passed as ``valid_mask`` to ``identify_disk``.
    """
    if int(edge_cells) != edge_cells or edge_cells < 0:
        raise ValueError("edge_cells must be a non-negative integer")
    mask = np.ones(shape, dtype=bool)
    e = int(edge_cells)
    if e:
        mask[:e, :] = mask[-e:, :] = False
        mask[:, :e] = mask[:, -e:] = False
    return mask


def identify_regions_frame(rain: np.ndarray, radius_cells: float, threshold: float = 1.0,
                           large_join: str = "almost", small_join: str = "almost",
                           min_size: int = 1, valid_mask: Optional[np.ndarray] = None,
                           erosion_border: str = "ignore",
                           bridge_radius_cells: Optional[float] = None,
                           edge_cells: int = 0) -> np.ndarray:
    if large_join not in ("almost", "region") or small_join not in ("almost", "region"):
        raise ValueError("large_join and small_join must be 'almost' or 'region'")
    if radius_cells < 0 or threshold < 0 or min_size < 1:
        raise ValueError("radius_cells, threshold must be non-negative; min_size >= 1")
    rain = np.asarray(rain, dtype=float)
    if rain.ndim != 2:
        raise ValueError("rain must be two-dimensional")
    valid = np.isfinite(rain) & edge_mask(rain.shape, edge_cells)
    if valid_mask is not None:
        valid &= np.asarray(valid_mask, dtype=bool)
    wet = valid & (rain >= threshold) & (rain > 0)
    if not wet.any():
        return np.zeros(rain.shape, dtype=np.int32)
    r = float(radius_cells)
    rb = r if bridge_radius_cells is None else float(bridge_radius_cells)
    if rb < 0:
        raise ValueError("bridge_radius_cells must be non-negative")
    regions, n_regions = ndimage.label(wet, structure=_EIGHT)
    eroded = erode_disk(wet, r, erosion_border)
    is_large = np.zeros(n_regions + 1, dtype=bool)
    is_large[np.unique(regions[eroded])] = True
    is_large[0] = False

    group = np.zeros(n_regions + 1, dtype=np.int64)
    next_group = 1
    large_ids = np.flatnonzero(is_large)
    if len(large_ids):
        if large_join == "region":
            group[large_ids] = np.arange(1, len(large_ids) + 1)
            next_group = len(large_ids) + 1
        else:
            opened = dilate_disk(eroded, r) & wet
            almost = dilate_disk(opened, rb) & valid
            clusters, n_clusters = ndimage.label(almost, structure=_EIGHT)
            reg, clu = regions[opened], clusters[opened]
            index = np.zeros(n_regions + 1, dtype=np.int64)
            index[large_ids] = np.arange(len(large_ids))
            pairs = np.unique(np.stack([reg, clu], axis=1), axis=0)
            n_nodes = len(large_ids) + n_clusters
            graph = coo_matrix((np.ones(len(pairs)),
                                (index[pairs[:, 0]], len(large_ids) + pairs[:, 1] - 1)),
                               shape=(n_nodes, n_nodes))
            _, component = connected_components(graph, directed=False)
            _, compact = np.unique(component[:len(large_ids)], return_inverse=True)
            group[large_ids] = compact + 1
            next_group = int(compact.max()) + 2

    large_raster = group[regions] * is_large[regions]
    small_ids = np.flatnonzero(~is_large)
    small_ids = small_ids[small_ids > 0]
    unassigned = []
    if len(small_ids):
        unassigned = _attach_small(regions, large_raster, small_ids, rb, wet.shape, group)

    if unassigned:
        unassigned = np.asarray(unassigned)
        if small_join == "region":
            group[unassigned] = next_group + np.arange(len(unassigned))
        else:
            remaining = np.isin(regions, unassigned)
            clusters, _ = ndimage.label(dilate_disk(remaining, rb) & valid, structure=_EIGHT)
            slices = ndimage.find_objects(regions)
            for rid in unassigned:
                sl = slices[rid - 1]
                group[rid] = next_group + int(clusters[sl][regions[sl] == rid][0]) - 1

    labels = group[regions]
    labels[~wet] = 0
    labels = _stable_relabel(labels)
    if min_size > 1 and labels.max() > 0:
        counts = np.bincount(labels.ravel())
        keep = counts >= min_size
        keep[0] = False
        labels[~keep[labels]] = 0
        labels = _stable_relabel(labels)
    return labels


def identify_regions(data: np.ndarray, radius_cells: float, threshold: float = 1.0, **kwargs) -> np.ndarray:
    """Apply :func:`identify_regions_frame` to every frame of a (time, y, x) array."""
    data = np.asarray(data)
    if data.ndim != 3:
        raise ValueError("data must have shape (time, y, x)")
    valid = kwargs.pop("valid_mask", None)
    return np.stack([identify_regions_frame(data[t], radius_cells, threshold,
                                            valid_mask=None if valid is None else valid[t], **kwargs)
                     for t in range(data.shape[0])])


def _attach_small(regions, large_raster, small_ids, rb, shape, group):
    """Dilate each small region and join the large segment it covers most."""
    unassigned = []
    slices = ndimage.find_objects(regions)
    pad = int(math.floor(rb))
    H, W = shape
    for rid in small_ids:
        sl = slices[rid - 1]
        y0, y1 = max(0, sl[0].start - pad), min(H, sl[0].stop + pad)
        x0, x1 = max(0, sl[1].start - pad), min(W, sl[1].stop + pad)
        hit = large_raster[y0:y1, x0:x1][dilate_disk(regions[y0:y1, x0:x1] == rid, rb)]
        hit = hit[hit > 0]
        if len(hit):
            group[rid] = int(np.argmax(np.bincount(hit)))
        else:
            unassigned.append(rid)
    return unassigned

