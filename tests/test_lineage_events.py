import json
from pathlib import Path

import numpy as np
from scipy import ndimage

from step import track_with_graph
from step.lineage_events import event_factors, paper_event_ids
from step.tracking import StormObject, TrackEdge, TrackGraph, TrackNode

PRESET = json.loads((Path(__file__).resolve().parents[1] /
                     "configs/rainfall_lineage_4km_1h_reference_v1.json").read_text())


def _node(time, node_id, area):
    obj = StormObject(node_id, area, (0.0, 0.0), 1.0, ((0, 1), (0, 1)), 1.0, float(area), [])
    return TrackNode(time, node_id, node_id, node_id, obj)


def _graph(nodes, edges):
    g = TrackGraph()
    g.objects = [_node(*n) for n in nodes]
    g.edges = [TrackEdge(t, p, c, 1.0, ev) for t, p, c, ev in edges]
    return g


def test_split_shares_and_merge_continues_largest():
    # t0: 1 (big), 2 (small) -> t1: 3 merges both -> t2: 4, 5 split from 3
    g = _graph([(0, 1, 100), (0, 2, 10), (1, 3, 120), (2, 4, 60), (2, 5, 50)],
               [(1, 1, 3, "merge"), (1, 2, 3, "merge"), (2, 3, 4, "split"), (2, 3, 5, "split")])
    ev = paper_event_ids(g)
    assert ev[3] == ev[1] != ev[2]
    assert ev[4] == ev[5] == ev[1]
    f = event_factors(g, ev, 16.0, 1.0)
    assert f["N"] == 2 and f["sum_l"] == 4 and f["duration_h"] == 2.0


def test_merge_tie_goes_to_lowest_node_id():
    g = _graph([(0, 1, 50), (0, 2, 50), (1, 3, 90)],
               [(1, 2, 3, "merge"), (1, 1, 3, "merge")])
    ev = paper_event_ids(g)
    assert ev[3] == ev[1]


def test_events_on_real_tracker_graph():
    shape = (80, 120)
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    blob = lambda cy, cx, r: np.where((yy - cy) ** 2 + (xx - cx) ** 2 <= r * r, 2.0, 0.0)
    rain = np.stack([blob(40, 60, 12), blob(31, 60, 8) + blob(49, 60, 8)])
    labels = np.stack([ndimage.label(f > 0)[0] for f in rain]).astype(np.int32)
    p = dict(PRESET["tracking"]); p["km"] = p.pop("max_displacement")
    _, g = track_with_graph(labels, rain, **p)
    assert {e.event for e in g.edges} == {"split"}
    assert len({n.branch_id for n in g.objects}) == 3
    assert len(set(paper_event_ids(g).values())) == 1
