"""Opt-in coverage scoring and shared parents for the lineage tracker.

Two independent options on top of :class:`step.tracking.Tracker`; with both at
their defaults the results are identical to ``Tracker``.

``overlap_metric='coverage'``
    The overlap term of the adjacent score uses the overlap coefficient,
    intersection / min(parent area, child area), instead of IoU. IoU penalizes
    pairs of very different size: a 50-cell child lying entirely inside a
    3000-cell parent has IoU ~0.02 and cannot reach ``tau``. Split/merge and
    raw-overlap rules are unchanged (they reconstruct intersections from IoU).
    Caution: a small parent lying inside a large child also has coefficient 1,
    so the large child can inherit the small parent's branch.

``overlap_metric='child_coverage'``
    As above but with intersection / child area only, so only the small-child
    side of a size mismatch is favoured.

``share_parent=True``
    After the normal update, a child that received no edge but has at least
    ``share_min_coverage`` of its cells (and at least ``share_min_pixels``
    cells) inside one object of the previous frame joins that object's branch
    with a ``share`` edge, as in Chang et al. (2016) S1, where several segments
    may have the same predecessor. A branch can then hold several nodes at one
    time. The joined child carries no velocity (its centroid is not a motion
    estimate of the parent).

``rescue_overlap=True``
    Before sharing, unlinked children and parents that received no edge are
    paired when each is the other's largest raw overlap (at least
    ``rescue_min_pixels`` cells): first with the previous frame
    (``overlap_continue``), then with objects last seen two frames earlier
    that were not continued (``gap_overlap``, gap 1). This keeps one-to-one
    branch semantics. It addresses small, fast-changing cells whose IoU-based
    score stays below ``tau`` although they overlap, and the rarely used gap
    rule (``gap_tau`` and ``gap_min_pixels`` of the preset).

No existing module is changed.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .tracking import (Candidate, TrackEdge, TrackGraph, Tracker, _find,
                       _set_overlap, coherent_score)

_EXTRA_KEYS = ("overlap_metric_coverage_v1", "share_parent_v1",
               "share_min_coverage", "share_min_pixels", "rescue_overlap_v1",
               "rescue_min_pixels")


class CoverageTracker(Tracker):
    def __init__(self, *args, overlap_metric="iou", share_parent=False,
                 share_min_coverage=0.5, share_min_pixels=1, rescue_overlap=False,
                 rescue_min_pixels=1, **kwargs):
        if overlap_metric not in ("iou", "coverage", "child_coverage"):
            raise ValueError("overlap_metric must be 'iou', 'coverage' or 'child_coverage'")
        if not 0 < share_min_coverage <= 1 or int(share_min_pixels) != share_min_pixels \
                or share_min_pixels < 1:
            raise ValueError("share_min_coverage in (0, 1], share_min_pixels a positive integer")
        if int(rescue_min_pixels) != rescue_min_pixels or rescue_min_pixels < 1:
            raise ValueError("rescue_min_pixels must be a positive integer")
        extra = {}
        if overlap_metric != "iou":
            extra["overlap_metric_coverage_v1"] = 1.0 if overlap_metric == "coverage" else 2.0
        if share_parent:
            extra.update(share_parent_v1=1.0, share_min_coverage=float(share_min_coverage),
                         share_min_pixels=float(share_min_pixels))
        if rescue_overlap:
            extra.update(rescue_overlap_v1=1.0, rescue_min_pixels=float(rescue_min_pixels))
        state = kwargs.get("state")
        saved = {}
        if state is not None and state.tracking_config:
            saved = {k: state.tracking_config.pop(k) for k in _EXTRA_KEYS
                     if k in state.tracking_config}
        super().__init__(*args, **kwargs)
        if self.state.last_time is not None and saved != extra:
            raise ValueError("coverage/share options differ from the checkpoint")
        self.state.tracking_config.update(extra)
        self.overlap_metric = overlap_metric
        self.share_parent = bool(share_parent)
        self.share_min_coverage = float(share_min_coverage)
        self.share_min_pixels = int(share_min_pixels)
        self.rescue_overlap = bool(rescue_overlap)
        self.rescue_min_pixels = int(rescue_min_pixels)

    # Copy of Tracker._candidates; only the overlap term differs under 'coverage'.
    def _candidates(self, parents, children, time):
        if self.overlap_metric == "iou":
            return super()._candidates(parents, children, time)
        if not parents or not children:
            return []
        tree = cKDTree(np.asarray([child.centroid for child in children]))
        boxes = np.asarray([child.bbox for child in children], dtype=float)
        centres = boxes.mean(axis=2)
        radii = np.linalg.norm((boxes[:, :, 1] - boxes[:, :, 0]) / 2, axis=1)
        overlap_tree = cKDTree(centres)
        output = []
        child_sets = {}
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
            indices = sorted(proximity | overlap_indices)
            if not indices:
                continue
            parent_pixels = set(terminal.node.object.pixels)
            moved_pixels = ({(y + shift[0], x + shift[1]) for y, x in parent_pixels}
                            if shift != (0, 0) else parent_pixels)
            for ci in indices:
                child = children[ci]
                distance = float(np.hypot(child.centroid[0] - predicted[0],
                                          child.centroid[1] - predicted[1]))
                if ci not in child_sets:
                    child_sets[ci] = set(child.pixels)
                pixels = child_sets[ci]
                raw_iou, raw_pc, raw_cc = _set_overlap(parent_pixels, pixels,
                                                      terminal.node.object.area, child.area)
                adv_iou, adv_pc, adv_cc = (_set_overlap(moved_pixels, pixels,
                    terminal.node.object.area, child.area) if shift != (0, 0)
                    else (raw_iou, raw_pc, raw_cc))
                coverage = max(raw_pc, raw_cc, adv_pc, adv_cc)
                if ci not in proximity and not (coverage > 0 and coverage >= self.event_overlap):
                    continue
                if self.overlap_metric == "coverage":
                    raw_ov, adv_ov = max(raw_pc, raw_cc), max(adv_pc, adv_cc)
                else:  # child_coverage
                    raw_ov, adv_ov = raw_cc, adv_cc
                intensity = min(terminal.node.object.mean_intensity,
                                child.mean_intensity) / max(
                                    terminal.node.object.mean_intensity,
                                    child.mean_intensity, 1e-12)
                score = float(0.60 * max(raw_ov, adv_ov) +
                              0.25 * np.exp(-distance / max(radius, 1.0)) +
                              0.15 * intensity)
                if self.score_policy == 'coherent_adjacent' and dt == 1:
                    raw_distance = float(np.linalg.norm(np.asarray(child.centroid) -
                                                       terminal.node.object.centroid))
                    score = coherent_score(raw_ov, adv_ov, raw_distance, distance, radius, intensity)
                output.append(Candidate(pi, ci, score, raw_iou, adv_iou,
                                        distance, max(raw_pc, adv_pc),
                                        max(raw_cc, adv_cc)))
        return output

    @staticmethod
    def _overlaps(pool, labeled_map):
        """{child label: {pool index (1-based): intersection cells}}."""
        index = np.zeros(labeled_map.shape, dtype=np.int64)
        for i, terminal in enumerate(pool, start=1):
            if terminal.node.object.pixels:
                ys, xs = zip(*terminal.node.object.pixels)
                index[np.asarray(ys), np.asarray(xs)] = i
        both = (labeled_map > 0) & (index > 0)
        out = {}
        if both.any():
            key = labeled_map[both].astype(np.int64) * (len(pool) + 1) + index[both]
            keys, counts = np.unique(key, return_counts=True)
            for k, n in zip(keys.tolist(), counts.tolist()):
                label, pi = divmod(k, len(pool) + 1)
                out.setdefault(label, {})[pi] = n
        return out

    def _join(self, time, node, parent, event, gap, inter, raster, labeled_map, graph, terminals):
        old, new = node.branch_id, parent.branch_id
        family = _find(self.state.family_parent, new)
        self.state.family_parent[old] = family
        node.branch_id, node.family_id = new, family
        live = terminals[node.node_id].node
        live.branch_id, live.family_id = new, family
        raster[labeled_map == node.object.label] = new
        union = parent.object.area + node.object.area - inter
        coverage = inter / node.object.area
        graph.edges.append(TrackEdge(
            int(time), parent.node_id, node.node_id, coverage, event, gap,
            inter / union, inter / union,
            float(np.hypot(node.object.centroid[0] - parent.object.centroid[0],
                           node.object.centroid[1] - parent.object.centroid[1])),
            inter / parent.object.area, coverage))

    def _rescue(self, time, pool, gap, raster, labeled_map, graph, terminals):
        """Mutual-best raw overlap between unlinked children and free parents."""
        if not pool:
            return
        linked = {e.child_id for e in graph.edges}
        children = {n.object.label: n for n in graph.objects if n.node_id not in linked}
        overlaps = {lab: v for lab, v in self._overlaps(pool, labeled_map).items() if lab in children}
        parent_best = {}
        for lab, per in overlaps.items():
            for pi, n in per.items():
                if pi not in parent_best or n > parent_best[pi][1] or \
                        (n == parent_best[pi][1] and lab < parent_best[pi][0]):
                    parent_best[pi] = (lab, n)
        joined = set()
        for lab, per in sorted(overlaps.items()):
            pi, n = max(per.items(), key=lambda item: (item[1], -item[0]))
            if n < self.rescue_min_pixels or parent_best[pi][0] != lab:
                continue
            parent = pool[pi - 1].node
            self._join(time, children[lab], parent, "gap_overlap" if gap else "overlap_continue",
                       gap, n, raster, labeled_map, graph, terminals)
            joined.add(parent.node_id)
        if joined:
            self.state.dormant = [d for d in self.state.dormant if d.node.node_id not in joined]

    def update(self, time, labeled_map, precip):
        time = int(time)
        previous = [t for t in self.state.active if t.node.time == time - 1]
        gap_pool = ([d for d in self.state.dormant if d.node.time == time - 2]
                    if self.rescue_overlap else [])
        raster, graph = super().update(time, labeled_map, precip)
        if not (self.share_parent or self.rescue_overlap):
            return raster, graph
        labeled_map = np.asarray(labeled_map)
        terminals = {t.node.node_id: t for t in self.state.active}
        if self.rescue_overlap:
            used = {e.parent_id for e in graph.edges}
            self._rescue(time, [t for t in previous if t.node.node_id not in used], 0,
                         raster, labeled_map, graph, terminals)
            used = {e.parent_id for e in graph.edges}
            self._rescue(time, [d for d in gap_pool if d.node.node_id not in used], 1,
                         raster, labeled_map, graph, terminals)
        if self.share_parent and previous:
            best = {lab: max(per.items(), key=lambda item: (item[1], -item[0]))
                    for lab, per in self._overlaps(previous, labeled_map).items()}
            linked = {e.child_id for e in graph.edges}
            for node in graph.objects:
                if node.node_id in linked or node.object.label not in best:
                    continue
                pi, inter = best[node.object.label]
                if inter / node.object.area < self.share_min_coverage or inter < self.share_min_pixels:
                    continue
                self._join(time, node, previous[pi - 1].node, "share", 0, inter,
                           raster, labeled_map, graph, terminals)
        branches = (self.state.family_parent if self.family_map_scope == 'global'
                    else {n.branch_id for n in graph.objects})
        graph.family_map = {b: _find(self.state.family_parent, b) for b in sorted(branches)}
        for node in graph.objects:
            node.family_id = graph.family_map[node.branch_id]
        return raster, graph


def track_with_graph_coverage(labeled_maps, precip_data, tau=0.35, km=20.0, max_gap=1,
                              gap_tau=None, gap_ambiguity=0.05, event_overlap=0.10,
                              state=None, start_time=None, return_state=False,
                              sequence_id=None, overlap_metric="iou", share_parent=False,
                              share_min_coverage=0.5, share_min_pixels=1,
                              rescue_overlap=False, rescue_min_pixels=1, **options):
    """Like :func:`step.tracking.track_with_graph` with the two extra options.

    ``options`` are the remaining ``Tracker`` keyword arguments (event_policy,
    adjacent_policy, score_policy, ...). ``km`` is grid cells per frame.
    """
    labels, precip = np.asarray(labeled_maps), np.asarray(precip_data)
    if labels.shape != precip.shape or labels.ndim != 3:
        raise ValueError("labels and precipitation need identical (time, y, x) shapes")
    if start_time is None:
        start_time = 0 if state is None or state.last_time is None else state.last_time + 1
    tracker = CoverageTracker(tau=tau, max_displacement=km, max_gap=max_gap, gap_tau=gap_tau,
                              gap_ambiguity=gap_ambiguity, event_overlap=event_overlap,
                              state=state, sequence_id=sequence_id,
                              overlap_metric=overlap_metric, share_parent=share_parent,
                              share_min_coverage=share_min_coverage,
                              share_min_pixels=share_min_pixels,
                              rescue_overlap=rescue_overlap,
                              rescue_min_pixels=rescue_min_pixels, **options)
    result, combined = np.zeros(labels.shape, dtype=np.int64), TrackGraph()
    for offset in range(labels.shape[0]):
        raster, graph = tracker.update(start_time + offset, labels[offset], precip[offset])
        result[offset] = raster
        combined.objects.extend(graph.objects)
        combined.edges.extend(graph.edges)
        combined.events.extend(graph.events)
        if tracker.family_map_scope == 'global':
            combined.family_map = graph.family_map
        else:
            combined.family_map.update(graph.family_map)
    if tracker.family_map_scope == 'observed':
        combined.family_map = {b: _find(tracker.state.family_parent, b)
                               for b in sorted(combined.family_map)}
    if return_state:
        return result, combined, tracker.state
    return result, combined
