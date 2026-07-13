"""Sparse, object-based storm tracking.

This module intentionally never forms a pixel-by-pixel distance matrix.  It
only compares objects whose predicted centres can plausibly meet, then scores
their overlap, distance, and intensity continuity.
"""

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage
from skimage.segmentation import relabel_sequential


@dataclass
class StormObject:
    label: int
    area: int
    centroid: tuple
    mean_intensity: float
    bbox: tuple


@dataclass
class TrackEdge:
    time: int
    parent_id: int
    child_id: int
    score: float
    event: str


@dataclass
class TrackGraph:
    """Compact event graph; raster labels alone cannot represent lineage."""
    edges: list = field(default_factory=list)
    objects: list = field(default_factory=list)


def _objects(labels, precip):
    labels = np.asarray(labels)
    boxes = ndimage.find_objects(labels)
    result = []
    for label, box in enumerate(boxes, start=1):
        if box is None:
            continue
        local = labels[box] == label
        area = int(local.sum())
        if not area:
            continue
        values = precip[box][local]
        weight = float(values.sum())
        yy, xx = np.nonzero(local)
        if weight > 0:
            cy = float(box[0].start + np.dot(yy, values) / weight)
            cx = float(box[1].start + np.dot(xx, values) / weight)
        else:
            cy = float(box[0].start + yy.mean())
            cx = float(box[1].start + xx.mean())
        result.append(StormObject(label, area, (cy, cx), float(values.mean()), box))
    return result


def _iou(previous_map, previous, current_map, current):
    y0 = max(previous.bbox[0].start, current.bbox[0].start)
    y1 = min(previous.bbox[0].stop, current.bbox[0].stop)
    x0 = max(previous.bbox[1].start, current.bbox[1].start)
    x1 = min(previous.bbox[1].stop, current.bbox[1].stop)
    if y0 >= y1 or x0 >= x1:
        return 0.0
    p = previous_map[y0:y1, x0:x1] == previous.label
    c = current_map[y0:y1, x0:x1] == current.label
    intersection = int(np.logical_and(p, c).sum())
    return intersection / float(previous.area + current.area - intersection)


def _score(previous_map, previous, current_map, current, max_displacement):
    dy = current.centroid[0] - previous.centroid[0]
    dx = current.centroid[1] - previous.centroid[1]
    distance = float(np.hypot(dy, dx))
    # A modest size-dependent margin handles expanding systems without opening
    # comparisons to all objects in the domain.
    reach = max_displacement + 0.5 * (np.sqrt(previous.area) + np.sqrt(current.area))
    if distance > reach:
        return None
    overlap = _iou(previous_map, previous, current_map, current)
    distance_score = np.exp(-distance / max(max_displacement, 1.0))
    intensity_score = min(previous.mean_intensity, current.mean_intensity) / max(previous.mean_intensity, current.mean_intensity, 1e-12)
    return 0.60 * overlap + 0.25 * distance_score + 0.15 * intensity_score


def track_with_graph(labeled_maps, precip_data, tau=0.35, phi=None, km=20.0, workers=1):
    """Track labelled rain objects and return ``(persistent_labels, graph)``.

    ``phi`` is retained only for STEP call compatibility. ``km`` has the old
    meaning: maximum one-timestep displacement expressed in grid cells.  The
    feature extraction is linear in active pixels; matching is sparse.
    """
    labels = np.asarray(labeled_maps)
    precip = np.asarray(precip_data)
    if labels.shape != precip.shape or labels.ndim != 3:
        raise ValueError("labeled_maps and precip_data must have identical (time, y, x) shapes")
    if not 0.0 <= tau <= 1.0:
        raise ValueError("tau must be in [0, 1]")
    if km <= 0:
        raise ValueError("km must be positive")

    result = np.zeros(labels.shape, dtype=np.int32)
    graph = TrackGraph()
    next_id = 1
    previous_objects = []
    previous_ids = {}

    for time in range(labels.shape[0]):
        current_map = labels[time].astype(np.int32, copy=False)
        current_objects = _objects(current_map, precip[time])
        current_ids = {}
        if time == 0:
            for obj in current_objects:
                current_ids[obj.label] = next_id
                graph.objects.append((time, next_id, obj))
                next_id += 1
        else:
            candidates = []
            for prev in previous_objects:
                for curr in current_objects:
                    score = _score(labels[time - 1], prev, current_map, curr, km)
                    if score is not None and score >= tau:
                        candidates.append((score, prev, curr))
            candidates.sort(key=lambda item: item[0], reverse=True)
            by_child = {}
            for score, prev, curr in candidates:
                by_child.setdefault(curr.label, []).append((score, prev))
            for curr in current_objects:
                parents = by_child.get(curr.label, [])
                if parents:
                    score, primary = parents[0]
                    current_ids[curr.label] = previous_ids[primary.label]
                    event = "continue" if len(parents) == 1 else "merge"
                    for parent_score, parent in parents:
                        graph.edges.append(TrackEdge(time, previous_ids[parent.label], current_ids[curr.label], parent_score, event))
                else:
                    current_ids[curr.label] = next_id
                    next_id += 1
                graph.objects.append((time, current_ids[curr.label], curr))

        for obj in current_objects:
            result[time][current_map == obj.label] = current_ids[obj.label]
        previous_objects = current_objects
        previous_ids = current_ids

    return relabel_sequential(result)[0].astype(np.int32, copy=False), graph


def track(labeled_maps, precip_data, tau, phi, km, test=False, workers=1):
    """Backward-compatible STEP entry point returning only persistent labels."""
    if test:
        print("STEP_UPGRADE uses sparse object matching; phi is ignored.")
    result, _ = track_with_graph(labeled_maps, precip_data, tau=tau, phi=phi, km=km, workers=workers)
    return result
