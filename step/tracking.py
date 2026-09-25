"""Deterministic, stateful precipitation-object lineage tracking.

The tracker distinguishes observed object nodes, uninterrupted branches, and
split/merge families. It can resume across chunk and month boundaries without
changing identifiers or graph semantics.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree


@dataclass
class StormObject:
    label: int
    area: int
    centroid: Tuple[float, float]
    mean_intensity: float
    bbox: Tuple[Tuple[int, int], Tuple[int, int]]
    max_intensity: float
    precipitation_sum: float
    pixels: List[Tuple[int, int]] = field(repr=False)


@dataclass
class TrackNode:
    time: int
    node_id: int
    branch_id: int
    family_id: int
    object: StormObject


@dataclass
class TrackEdge:
    time: int
    parent_id: int
    child_id: int
    score: float
    event: str
    gap: int = 0
    raw_iou: float = 0.0
    advected_iou: float = 0.0
    distance: float = 0.0
    parent_coverage: float = 0.0
    child_coverage: float = 0.0


@dataclass
class TrackEvent:
    event_id: int
    time: int
    event: str
    parent_ids: List[int]
    child_ids: List[int]


@dataclass
class Candidate:
    parent_index: int
    child_index: int
    score: float
    raw_iou: float
    advected_iou: float
    distance: float
    parent_coverage: float
    child_coverage: float


@dataclass
class Terminal:
    node: TrackNode
    previous_centroid: Optional[Tuple[float, float]] = None
    previous_time: Optional[int] = None


@dataclass
class TrackGraph:
    """Lineage graph produced for one call or streaming chunk."""

    objects: List[TrackNode] = field(default_factory=list)
    edges: List[TrackEdge] = field(default_factory=list)
    events: List[TrackEvent] = field(default_factory=list)
    family_map: Dict[int, int] = field(default_factory=dict)

    @property
    def nodes(self):
        return self.objects


@dataclass
class TrackingState:
    """Minimal state required to resume a time-ordered sequence."""

    last_time: Optional[int] = None
    next_node_id: int = 1
    next_branch_id: int = 1
    next_event_id: int = 1
    active: List[Terminal] = field(default_factory=list)
    dormant: List[Terminal] = field(default_factory=list)
    family_parent: Dict[int, int] = field(default_factory=dict)
    sequence_id: Optional[str] = None
    grid_shape: Optional[Tuple[int, int]] = None
    tracking_config: Dict[str, float] = field(default_factory=dict)


def _storm_to_dict(obj: StormObject) -> dict:
    value = asdict(obj)
    value["centroid"] = list(obj.centroid)
    value["bbox"] = [list(obj.bbox[0]), list(obj.bbox[1])]
    value["pixels"] = [list(pixel) for pixel in obj.pixels]
    return value


def _storm_from_dict(value: dict) -> StormObject:
    return StormObject(
        label=int(value["label"]), area=int(value["area"]),
        centroid=tuple(value["centroid"]),
        mean_intensity=float(value["mean_intensity"]),
        bbox=(tuple(value["bbox"][0]), tuple(value["bbox"][1])),
        max_intensity=float(value["max_intensity"]),
        precipitation_sum=float(value["precipitation_sum"]),
        pixels=[tuple(pixel) for pixel in value["pixels"]],
    )


def _node_to_dict(node: TrackNode) -> dict:
    return {
        "time": node.time, "node_id": node.node_id,
        "branch_id": node.branch_id, "family_id": node.family_id,
        "object": _storm_to_dict(node.object),
    }


def _node_from_dict(value: dict) -> TrackNode:
    return TrackNode(
        time=int(value["time"]), node_id=int(value["node_id"]),
        branch_id=int(value["branch_id"]), family_id=int(value["family_id"]),
        object=_storm_from_dict(value["object"]),
    )


def _terminal_to_dict(terminal: Terminal) -> dict:
    return {
        "node": _node_to_dict(terminal.node),
        "previous_centroid": (
            list(terminal.previous_centroid)
            if terminal.previous_centroid is not None else None
        ),
        "previous_time": terminal.previous_time,
    }


def _terminal_from_dict(value: dict) -> Terminal:
    centroid = value.get("previous_centroid")
    return Terminal(
        node=_node_from_dict(value["node"]),
        previous_centroid=tuple(centroid) if centroid is not None else None,
        previous_time=value.get("previous_time"),
    )


def state_to_dict(state: TrackingState) -> dict:
    return {
        "schema_version": 1, "last_time": state.last_time,
        "next_node_id": state.next_node_id,
        "next_branch_id": state.next_branch_id,
        "next_event_id": state.next_event_id,
        "active": [_terminal_to_dict(item) for item in state.active],
        "dormant": [_terminal_to_dict(item) for item in state.dormant],
        "family_parent": {str(k): v for k, v in state.family_parent.items()},
        "sequence_id": state.sequence_id,
        "grid_shape": list(state.grid_shape) if state.grid_shape is not None else None,
        "tracking_config": state.tracking_config,
    }


def state_from_dict(value: dict) -> TrackingState:
    if value.get("schema_version") != 1:
        raise ValueError("unsupported tracking-state schema")
    shape = value.get("grid_shape")
    return TrackingState(
        last_time=value.get("last_time"),
        next_node_id=int(value["next_node_id"]),
        next_branch_id=int(value["next_branch_id"]),
        next_event_id=int(value["next_event_id"]),
        active=[_terminal_from_dict(item) for item in value.get("active", [])],
        dormant=[_terminal_from_dict(item) for item in value.get("dormant", [])],
        family_parent={int(k): int(v) for k, v in value.get("family_parent", {}).items()},
        sequence_id=value.get("sequence_id"),
        grid_shape=tuple(shape) if shape is not None else None,
        tracking_config={
            str(k): float(v) for k, v in value.get("tracking_config", {}).items()
        },
    )


def save_tracking_state(state: TrackingState, path) -> None:
    """Atomically save a portable JSON checkpoint."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w") as handle:
        json.dump(state_to_dict(state), handle, sort_keys=True)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def load_tracking_state(path) -> TrackingState:
    with Path(path).open() as handle:
        return state_from_dict(json.load(handle))


def _objects(labels, precip) -> List[StormObject]:
    labels, precip = np.asarray(labels), np.asarray(precip)
    result = []
    for label, slices in enumerate(ndimage.find_objects(labels), start=1):
        if slices is None:
            continue
        local = labels[slices] == label
        area = int(local.sum())
        if not area:
            continue
        values = np.asarray(precip[slices][local], dtype=float)
        yy, xx = np.nonzero(local)
        gy, gx = yy + slices[0].start, xx + slices[1].start
        weight = float(values.sum())
        if weight > 0:
            centroid = (float(np.dot(gy, values) / weight),
                        float(np.dot(gx, values) / weight))
        else:
            centroid = (float(gy.mean()), float(gx.mean()))
        result.append(StormObject(
            label, area, centroid, float(values.mean()),
            ((slices[0].start, slices[0].stop),
             (slices[1].start, slices[1].stop)),
            float(values.max()), weight,
            list(zip(gy.tolist(), gx.tolist())),
        ))
    return result


def _object_without_pixels(obj: StormObject) -> StormObject:
    """Return the table representation; masks remain only in live state."""
    return StormObject(
        obj.label, obj.area, obj.centroid, obj.mean_intensity, obj.bbox,
        obj.max_intensity, obj.precipitation_sum, [],
    )


def _find(parent: Dict[int, int], value: int) -> int:
    parent.setdefault(value, value)
    root = value
    while parent[root] != root:
        root = parent[root]
    while parent[value] != value:
        following = parent[value]
        parent[value] = root
        value = following
    return root


def _union(parent: Dict[int, int], values: Iterable[int]) -> int:
    roots = sorted({_find(parent, value) for value in values})
    if not roots:
        raise ValueError("cannot union an empty family")
    for other in roots[1:]:
        parent[other] = roots[0]
    return roots[0]


def _pixel_overlap(previous, current, shift=(0, 0)):
    p = {(y + shift[0], x + shift[1]) for y, x in previous.pixels}
    c = set(current.pixels)
    intersection = len(p & c)
    union = len(p) + len(c) - intersection
    return (
        intersection / float(union) if union else 1.0,
        intersection / float(previous.area),
        intersection / float(current.area),
    )


class Tracker:
    """Stateful tracker for one climate/year/member sequence."""

    def __init__(self, tau=0.35, max_displacement=20.0, max_gap=1,
                 gap_tau=None, gap_ambiguity=0.05, event_overlap=0.10,
                 state=None, sequence_id=None):
        if not 0 <= tau <= 1 or not 0 <= event_overlap <= 1:
            raise ValueError("tau and event_overlap must be in [0, 1]")
        if max_displacement <= 0 or max_gap < 0:
            raise ValueError("invalid displacement or gap")
        self.tau = float(tau)
        self.max_displacement = float(max_displacement)
        self.max_gap = int(max_gap)
        self.gap_tau = float(gap_tau if gap_tau is not None else max(tau, 0.45))
        if gap_ambiguity < 0:
            raise ValueError("gap_ambiguity must be non-negative")
        self.gap_ambiguity = float(gap_ambiguity)
        self.event_overlap = float(event_overlap)
        self.state = state if state is not None else TrackingState()
        requested_sequence = sequence_id or self.state.sequence_id or "default"
        if self.state.sequence_id not in (None, requested_sequence):
            raise ValueError("checkpoint belongs to a different sequence_id")
        self.state.sequence_id = requested_sequence
        config = {
            "algorithm_revision": 3.0,
            "tau": self.tau,
            "max_displacement": self.max_displacement,
            "max_gap": float(self.max_gap),
            "gap_tau": self.gap_tau,
            "gap_ambiguity": self.gap_ambiguity,
            "event_overlap": self.event_overlap,
        }
        if self.state.tracking_config and self.state.tracking_config != config:
            raise ValueError("tracking parameters differ from the checkpoint")
        self.state.tracking_config = config

    def _new_branch(self):
        branch = self.state.next_branch_id
        self.state.next_branch_id += 1
        self.state.family_parent[branch] = branch
        return branch

    @staticmethod
    def _velocity(terminal):
        if terminal.previous_centroid is None or terminal.previous_time is None:
            return 0.0, 0.0
        dt = terminal.node.time - terminal.previous_time
        if dt <= 0:
            return 0.0, 0.0
        return tuple((terminal.node.object.centroid[i] - terminal.previous_centroid[i]) / dt
                     for i in range(2))

    def _candidates(self, parents, children, time):
        if not parents or not children:
            return []
        tree = cKDTree(np.asarray([child.centroid for child in children]))
        # Bounding-box centres/radii conservatively index potentially overlapping
        # masks, independently of intensity-weighted centroid movement.
        boxes = np.asarray([child.bbox for child in children], dtype=float)
        centres = boxes.mean(axis=2)
        radii = np.linalg.norm((boxes[:, :, 1] - boxes[:, :, 0]) / 2, axis=1)
        overlap_tree = cKDTree(centres)
        output = []
        for pi, terminal in enumerate(parents):
            dt = time - terminal.node.time
            vy, vx = self._velocity(terminal)
            predicted = (terminal.node.object.centroid[0] + vy * dt,
                         terminal.node.object.centroid[1] + vx * dt)
            radius = self.max_displacement * dt
            proximity = set(tree.query_ball_point(predicted, radius))
            proximity.update(tree.query_ball_point(terminal.node.object.centroid, radius))
            shift = (int(round(vy * dt)), int(round(vx * dt)))
            overlap_indices = set()
            parent_box = np.asarray(terminal.node.object.bbox, dtype=float)
            parent_radius = np.linalg.norm((parent_box[:, 1] - parent_box[:, 0]) / 2)
            for translation in ((0, 0), shift):
                moved = parent_box + np.asarray(translation)[:, None]
                possible = overlap_tree.query_ball_point(moved.mean(axis=1), parent_radius + radii.max())
                for ci in possible:
                    if np.all(np.minimum(moved[:, 1], boxes[ci, :, 1]) >
                              np.maximum(moved[:, 0], boxes[ci, :, 0])):
                        overlap_indices.add(ci)
            for ci in sorted(proximity | overlap_indices):
                child = children[ci]
                distance = float(np.hypot(child.centroid[0] - predicted[0],
                                          child.centroid[1] - predicted[1]))
                raw_iou, raw_pc, raw_cc = _pixel_overlap(terminal.node.object, child)
                adv_iou, adv_pc, adv_cc = _pixel_overlap(terminal.node.object, child, shift)
                # Box intersection alone never admits a distant candidate.
                coverage = max(raw_pc, raw_cc, adv_pc, adv_cc)
                if ci not in proximity and not (coverage > 0 and coverage >= self.event_overlap):
                    continue
                overlap = max(raw_iou, adv_iou)
                intensity = min(terminal.node.object.mean_intensity,
                                child.mean_intensity) / max(
                                    terminal.node.object.mean_intensity,
                                    child.mean_intensity, 1e-12)
                score = float(0.60 * overlap +
                              0.25 * np.exp(-distance / max(radius, 1.0)) +
                              0.15 * intensity)
                output.append(Candidate(pi, ci, score, raw_iou, adv_iou,
                                        distance, max(raw_pc, adv_pc),
                                        max(raw_cc, adv_cc)))
        return output

    @staticmethod
    def _components(candidates):
        p2c, c2p = {}, {}
        for item in candidates:
            p2c.setdefault(item.parent_index, set()).add(item.child_index)
            c2p.setdefault(item.child_index, set()).add(item.parent_index)
        seen, output = set(), []
        for start in sorted(p2c):
            if start in seen:
                continue
            parents, children, pending = set(), set(), [start]
            while pending:
                parent = pending.pop()
                if parent in parents:
                    continue
                parents.add(parent)
                for child in p2c.get(parent, ()):
                    children.add(child)
                    pending.extend(c2p.get(child, ()))
            seen.update(parents)
            output.append((sorted(parents), sorted(children)))
        return output

    @staticmethod
    def _assignment(candidates, parent_indices, child_indices, minimum):
        parents, children = set(parent_indices), set(child_indices)
        eligible = [c for c in candidates if c.parent_index in parents and
                    c.child_index in children and c.score >= minimum]
        # Independent connected components avoid one full-domain dense matrix.
        by_parent = {}
        for c in eligible:
            by_parent.setdefault(c.parent_index, []).append(c)
        selected = []
        for ps, cs in Tracker._components(eligible):
            ppos, cpos = {p:i for i,p in enumerate(ps)}, {c:i for i,c in enumerate(cs)}
            matrix = np.full((len(ps), len(cs) + len(ps)), -1e6)
            matrix[:, len(cs):] = 0.0  # One available unmatched slot per parent.
            lookup = {}
            for pi in ps:
                for item in by_parent[pi]:
                    key = (ppos[pi], cpos[item.child_index])
                    if item.score > matrix[key]:
                        matrix[key], lookup[key] = item.score, item
            rows, columns = linear_sum_assignment(-matrix)
            selected.extend(lookup[(r,c)] for r,c in zip(rows,columns) if (r,c) in lookup)
        return sorted(selected, key=lambda item: (item.child_index, item.parent_index))

    @staticmethod
    def _edge(time, parent, child_id, item, event, gap=0):
        return TrackEdge(time, parent.node.node_id, child_id, item.score, event,
                         gap, item.raw_iou, item.advected_iou, item.distance,
                         item.parent_coverage, item.child_coverage)

    def update(self, time, labeled_map, precip):
        time = int(time)
        if self.state.last_time is not None and time <= self.state.last_time:
            raise ValueError("tracking times must be strictly increasing")
        if self.state.last_time is not None and time != self.state.last_time + 1:
            raise ValueError(
                "tracking times must be contiguous; represent a below-threshold "
                "hour with an empty valid frame"
            )
        labeled_map, precip = np.asarray(labeled_map), np.asarray(precip)
        if labeled_map.shape != precip.shape or labeled_map.ndim != 2:
            raise ValueError("labels and precipitation must share a 2-D shape")
        if self.state.grid_shape not in (None, tuple(labeled_map.shape)):
            raise ValueError("grid shape differs from the checkpoint")
        self.state.grid_shape = tuple(labeled_map.shape)

        objects, active = _objects(labeled_map, precip), list(self.state.active)
        graph = TrackGraph()
        raster = np.zeros(labeled_map.shape, dtype=np.int64)
        node_ids = list(range(self.state.next_node_id,
                              self.state.next_node_id + len(objects)))
        self.state.next_node_id += len(objects)
        candidates = self._candidates(active, objects, time)
        event_candidates = [c for c in candidates if c.score >= self.tau and
                            max(c.parent_coverage, c.child_coverage) >= self.event_overlap]
        used_p, used_c, child_branch, child_parent = set(), set(), {}, {}

        for parents, children in self._components(event_candidates):
            if len(parents) == len(children) == 1:
                continue
            event = "split" if len(parents) == 1 else (
                "merge" if len(children) == 1 else "complex")
            branches = [self._new_branch() for _ in children]
            parent_branches = [active[i].node.branch_id for i in parents]
            family = _union(self.state.family_parent, parent_branches + branches)
            for ci, branch in zip(children, branches):
                child_branch[ci] = branch
                self.state.family_parent[branch] = family
            graph.events.append(TrackEvent(
                self.state.next_event_id, time, event,
                [active[i].node.node_id for i in parents],
                [node_ids[i] for i in children]))
            self.state.next_event_id += 1
            for item in event_candidates:
                if item.parent_index in parents and item.child_index in children:
                    graph.edges.append(self._edge(
                        time, active[item.parent_index], node_ids[item.child_index],
                        item, event))
            used_p.update(parents)
            used_c.update(children)

        remaining_p = [i for i in range(len(active)) if i not in used_p]
        remaining_c = [i for i in range(len(objects)) if i not in used_c]
        for item in self._assignment(candidates, remaining_p, remaining_c, self.tau):
            parent = active[item.parent_index]
            child_branch[item.child_index] = parent.node.branch_id
            child_parent[item.child_index] = parent
            used_p.add(item.parent_index)
            used_c.add(item.child_index)
            graph.edges.append(self._edge(time, parent, node_ids[item.child_index],
                                           item, "continue"))

        dormant = [d for d in self.state.dormant
                   if time - d.node.time <= self.max_gap + 1]
        dormant.extend(active[i] for i in range(len(active)) if i not in used_p)
        gap_children = [i for i in range(len(objects)) if i not in used_c]
        gap_candidates = [
            item for item in self._candidates(dormant, objects, time)
            if 1 < time - dormant[item.parent_index].node.time <= self.max_gap + 1
        ]
        by_child = {}
        for item in gap_candidates:
            by_child.setdefault(item.child_index, []).append(item)
        ambiguous_children = {
            child for child, items in by_child.items()
            if len(items) > 1
            and sorted((item.score for item in items), reverse=True)[0]
            - sorted((item.score for item in items), reverse=True)[1]
            < self.gap_ambiguity
        }
        gap_candidates = [
            item for item in gap_candidates
            if item.child_index not in ambiguous_children
        ]
        used_dormant = set()
        for item in self._assignment(gap_candidates, range(len(dormant)),
                                     gap_children, self.gap_tau):
            elapsed = time - dormant[item.parent_index].node.time
            if elapsed <= 1 or elapsed > self.max_gap + 1:
                continue
            parent = dormant[item.parent_index]
            child_branch[item.child_index] = parent.node.branch_id
            child_parent[item.child_index] = parent
            used_c.add(item.child_index)
            used_dormant.add(item.parent_index)
            graph.edges.append(self._edge(time, parent, node_ids[item.child_index],
                                           item, "gap_continue", elapsed - 1))

        new_active = []
        for ci, obj in enumerate(objects):
            branch = child_branch.get(ci)
            if branch is None:
                branch = self._new_branch()
            family = _find(self.state.family_parent, branch)
            node = TrackNode(time, node_ids[ci], branch, family, obj)
            parent = child_parent.get(ci)
            terminal = Terminal(node) if parent is None else Terminal(
                node, parent.node.object.centroid, parent.node.time)
            new_active.append(terminal)
            graph.objects.append(TrackNode(
                node.time, node.node_id, node.branch_id, node.family_id,
                _object_without_pixels(node.object),
            ))
            raster[labeled_map == obj.label] = branch

        self.state.active = new_active
        self.state.dormant = [d for i, d in enumerate(dormant)
                              if i not in used_dormant and
                              time - d.node.time <= self.max_gap]
        self.state.last_time = time
        graph.family_map = {b: _find(self.state.family_parent, b)
                            for b in sorted(self.state.family_parent)}
        for node in graph.objects:
            node.family_id = graph.family_map[node.branch_id]
        return raster, graph


def track_with_graph(labeled_maps, precip_data, tau=0.35, phi=None, km=20.0,
                     workers=1, max_gap=1, gap_tau=None, gap_ambiguity=0.05,
                     event_overlap=0.10, state=None, start_time=None,
                     return_state=False, sequence_id=None):
    """Track a chunk and optionally return resumable state.

    ``km`` retains the legacy name and is grid cells per frame. ``phi`` and
    ``workers`` remain accepted for compatibility. Pass returned state into
    the next chunk/month, or persist it with :func:`save_tracking_state`.
    """
    labels, precip = np.asarray(labeled_maps), np.asarray(precip_data)
    if labels.shape != precip.shape or labels.ndim != 3:
        raise ValueError("labels and precipitation need identical (time, y, x) shapes")
    if start_time is None:
        start_time = 0 if state is None or state.last_time is None else state.last_time + 1
    tracker = Tracker(
        tau, km, max_gap, gap_tau, gap_ambiguity, event_overlap, state,
        sequence_id
    )
    result, combined = np.zeros(labels.shape, dtype=np.int64), TrackGraph()
    for offset in range(labels.shape[0]):
        raster, graph = tracker.update(start_time + offset, labels[offset], precip[offset])
        result[offset] = raster
        combined.objects.extend(graph.objects)
        combined.edges.extend(graph.edges)
        combined.events.extend(graph.events)
        combined.family_map = graph.family_map
    if return_state:
        return result, combined, tracker.state
    return result, combined


def track(labeled_maps, precip_data, tau, phi, km, test=False, workers=1, **kwargs):
    """Backward-compatible entry point returning branch labels only."""
    if test and phi is not None:
        print("STEP_UPGRADE accepts phi for compatibility; lineage scoring ignores it.")
    result, _ = track_with_graph(labeled_maps, precip_data, tau=tau, phi=phi,
                                 km=km, workers=workers, **kwargs)
    return result
