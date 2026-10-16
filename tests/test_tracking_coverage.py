import json
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

from step import track_with_graph
from step.tracking_coverage import track_with_graph_coverage

PRESET = json.loads((Path(__file__).resolve().parents[1] /
                     "configs/rainfall_lineage_4km_1h_reference_v1.json").read_text())


def _params():
    p = dict(PRESET["tracking"])
    p["km"] = p.pop("max_displacement")
    return p


def _blob(shape, cy, cx, r, value=2.0):
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    return np.where((yy - cy) ** 2 + (xx - cx) ** 2 <= r * r, value, 0.0)


def _labels(rain):
    return np.stack([ndimage.label(f > 0, structure=np.ones((3, 3)))[0] for f in rain]).astype(np.int32)


def _edges(graph):
    return sorted((e.time, e.parent_id, e.child_id, e.event, round(e.score, 12)) for e in graph.edges)


def test_defaults_reproduce_tracker_exactly():
    rng = np.random.default_rng(0)
    shape = (90, 120)
    rain = []
    for t in range(8):
        f = sum(_blob(shape, 20 + 8 * k + rng.integers(-2, 3), 15 + 12 * k + 3 * t, 4 + k % 4,
                      1.0 + rng.random()) for k in range(7))
        f = f * (ndimage.gaussian_filter(rng.random(shape), 2) > 0.45)
        rain.append(f)
    rain = np.stack(rain)
    labels = _labels(rain)
    r0, g0 = track_with_graph(labels, rain, **_params())
    r1, g1 = track_with_graph_coverage(labels, rain, **_params())
    assert np.array_equal(r0, r1)
    assert _edges(g0) == _edges(g1)
    assert [(n.node_id, n.branch_id, n.family_id) for n in g0.objects] == \
           [(n.node_id, n.branch_id, n.family_id) for n in g1.objects]


def _inside_case():
    shape = (80, 80)
    t0 = _blob(shape, 40, 40, 20)
    t1 = _blob(shape, 40, 36, 11) + _blob(shape, 40, 55, 3)   # small cell inside old footprint
    return np.stack([t0, t1])


def test_small_child_inside_parent_is_new_by_default_and_shared_with_option():
    rain = _inside_case()
    labels = _labels(rain)
    small = labels[1, 40, 55]
    big = labels[1, 40, 36]
    raster, graph = track_with_graph(labels, rain, **_params())
    assert raster[1, 40, 55] != raster[0, 40, 40]
    assert raster[1, 40, 36] == raster[0, 40, 40]
    raster, graph = track_with_graph_coverage(labels, rain, share_parent=True, **_params())
    assert raster[1, 40, 55] == raster[1, 40, 36] == raster[0, 40, 40]
    share = [e for e in graph.edges if e.event == "share"]
    assert len(share) == 1 and share[0].child_coverage == 1.0
    assert small != big


def test_coverage_metric_links_a_lone_small_remnant():
    shape = (80, 80)
    rain = np.stack([_blob(shape, 40, 40, 20), _blob(shape, 40, 52, 3)])
    labels = _labels(rain)
    raster, _ = track_with_graph(labels, rain, **_params())
    assert raster[1, 40, 52] != raster[0, 40, 40]
    raster, _ = track_with_graph_coverage(labels, rain, overlap_metric="coverage", **_params())
    assert raster[1, 40, 52] == raster[0, 40, 40]


def test_checkpoint_options_must_match():
    rain = _inside_case()
    labels = _labels(rain)
    _, _, state = track_with_graph_coverage(labels[:1], rain[:1], share_parent=True,
                                            return_state=True, **_params())
    with pytest.raises(ValueError):
        track_with_graph_coverage(labels[1:], rain[1:], state=state, **_params())


def test_resume_matches_single_call():
    rain = _inside_case()
    rain = np.concatenate([rain, rain[::-1]])
    labels = _labels(rain)
    kw = dict(share_parent=True, overlap_metric="coverage", **_params())
    r_all, g_all = track_with_graph_coverage(labels, rain, **kw)
    r_a, g_a, state = track_with_graph_coverage(labels[:2], rain[:2], return_state=True, **kw)
    r_b, g_b = track_with_graph_coverage(labels[2:], rain[2:], state=state, **kw)
    assert np.array_equal(r_all, np.concatenate([r_a, r_b]))
    assert _edges(g_all) == sorted(_edges(g_a) + _edges(g_b))


def test_rescue_links_small_low_iou_overlap_with_free_parent():
    shape = (40, 60)
    t0 = np.zeros(shape); t0[20, 20:24] = 2.0; t0[21, 20] = 2.0      # 5 cells
    t1 = np.zeros(shape); t1[20, 23:31] = 1.0; t1[19:23, 28:31] = 1.0  # shifted, 1-cell overlap
    rain = np.stack([t0, t1])
    labels = _labels(rain)
    raster, _ = track_with_graph_coverage(labels, rain, **_params())
    assert raster[1, 20, 23] != raster[0, 20, 20]
    raster, graph = track_with_graph_coverage(labels, rain, rescue_overlap=True, **_params())
    assert raster[1, 20, 23] == raster[0, 20, 20]
    assert [e.event for e in graph.edges] == ["overlap_continue"]


def test_rescue_bridges_one_missing_frame_by_overlap():
    shape = (40, 60)
    a = np.zeros(shape); a[18:23, 18:23] = 1.2
    b = np.zeros(shape); b[19:24, 20:25] = 1.2
    rain = np.stack([a, np.zeros(shape), b])
    labels = _labels(rain)
    p = _params()
    p["gap_min_pixels"] = 100     # disable the built-in gap rule for this object size
    raster, _ = track_with_graph_coverage(labels, rain, **p)
    assert raster[2, 20, 22] != raster[0, 20, 20]
    raster, graph = track_with_graph_coverage(labels, rain, rescue_overlap=True, **p)
    assert raster[2, 20, 22] == raster[0, 20, 20]
    assert [(e.event, e.gap) for e in graph.edges] == [("gap_overlap", 1)]


def test_rescue_keeps_one_to_one_and_resumes():
    rain = _inside_case()
    rain = np.concatenate([rain, rain[::-1]])
    labels = _labels(rain)
    kw = dict(share_parent=True, rescue_overlap=True, overlap_metric="child_coverage", **_params())
    r_all, g_all = track_with_graph_coverage(labels, rain, **kw)
    r_a, g_a, state = track_with_graph_coverage(labels[:2], rain[:2], return_state=True, **kw)
    r_b, g_b = track_with_graph_coverage(labels[2:], rain[2:], state=state, **kw)
    assert np.array_equal(r_all, np.concatenate([r_a, r_b]))
    rescued = [e for e in g_all.edges if e.event in ("overlap_continue", "gap_overlap")]
    assert len({e.parent_id for e in rescued}) == len(rescued)
