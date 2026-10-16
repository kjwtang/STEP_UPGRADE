"""Paper-style event ids from a STEP_UPGRADE lineage graph.

``branch_id`` ends at every split and merge, and ``family_id`` unites every
branch that ever split or merged, so neither is the event of Chang et al.
(2016) S1, where each segment has exactly one predecessor and several segments
may share a predecessor. This post-processing assigns that event:

* a node without incoming edges starts a new event;
* a node with one incoming edge (continue, gap, share, rescue, or one parent
  of a split) takes the parent's event, so split children share it;
* a node with several incoming edges (merge, complex) takes the event of the
  parent with the most grid cells, as S1 chooses the largest predecessor
  (ties: lowest parent node id).

The graph is not modified. No existing module is changed.
"""
from __future__ import annotations

from typing import Dict


def paper_event_ids(graph) -> Dict[int, int]:
    """Return {node_id: event_id}; event ids start at 1 in time order."""
    nodes = sorted(graph.objects, key=lambda n: (n.time, n.node_id))
    by_id = {n.node_id: n for n in nodes}
    parents: Dict[int, list] = {}
    for e in graph.edges:
        if e.parent_id in by_id:
            parents.setdefault(e.child_id, []).append(e.parent_id)
    event: Dict[int, int] = {}
    next_event = 1
    for n in nodes:
        candidates = [p for p in parents.get(n.node_id, ()) if p in event]
        if not candidates:
            event[n.node_id] = next_event
            next_event += 1
            continue
        best = max(candidates, key=lambda p: (by_id[p].object.area, -p))
        event[n.node_id] = event[best]
    return event


def event_factors(graph, events: Dict[int, int], cell_area_km2: float, dt_hours: float):
    """Size, duration and number factors (paper Section 4b) for ``events``."""
    pairs = {(events[n.node_id], n.time) for n in graph.objects}
    wet = sum(n.object.area for n in graph.objects)
    n_events = len(set(events.values()))
    return dict(N=n_events, sum_l=len(pairs),
                size_km2=cell_area_km2 * wet / len(pairs),
                duration_h=dt_hours * len(pairs) / n_events)
