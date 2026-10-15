import math

import numpy as np
import pytest
from scipy import ndimage

from step import chang2016 as c


def _blob(shape, cy, cx, r, value=2.0):
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    return np.where((yy - cy) ** 2 + (xx - cx) ** 2 <= r * r, value, 0.0)


@pytest.mark.parametrize("radius", [1, 2.5, 4, 7])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_disk_morphology_matches_scipy(radius, seed):
    rng = np.random.default_rng(seed)
    mask = ndimage.gaussian_filter(rng.standard_normal((80, 90)), 3) > 0.1
    s = c.disk(radius)
    assert np.array_equal(c.dilate_disk(mask, radius), ndimage.binary_dilation(mask, s))
    assert np.array_equal(c.erode_disk(mask, radius, "background"),
                          ndimage.binary_erosion(mask, s, border_value=0))
    assert np.array_equal(c.erode_disk(mask, radius, "ignore"),
                          ndimage.binary_erosion(mask, s, border_value=1))


def test_paper_threshold_is_strict_and_in_mm_per_hour():
    rain = np.zeros((20, 20))
    rain[5, 5] = 0.1 / 3.0          # exactly 0.1 mm per 3 h: not > cutoff
    rain[10, 10] = 0.04
    labels = c.identify_frame(rain, radius_km=0, grid_km=4)
    assert labels[5, 5] == 0 and labels[10, 10] == 1


def test_small_regions_join_nearest_large_segment_without_chaining():
    shape = (60, 120)
    # Two large blobs 40 cells apart, joined by a string of small specks.
    rain = _blob(shape, 30, 20, 10) + _blob(shape, 30, 100, 10)
    for x in range(33, 90, 4):
        rain[30, x] = 1.0
    labels = c.identify_frame(rain, radius_km=12, grid_km=4)  # r = 3 cells
    left, right = labels[30, 20], labels[30, 100]
    assert left and right and left != right
    # specks near the left blob join it; mid specks form their own segment(s)
    assert labels[30, 33] == left and labels[30, 89] == right
    # plain almost-connected labelling would chain everything together
    plain = ndimage.label(c.dilate_disk(rain > 0, 3), structure=np.ones((3, 3)))[0]
    assert plain[30, 20] == plain[30, 100]


def test_nearby_large_regions_group_into_one_segment():
    shape = (60, 80)
    rain = _blob(shape, 30, 20, 8) + _blob(shape, 30, 42, 8)  # gap of 5 cells
    labels = c.identify_frame(rain, radius_km=12, grid_km=4)
    assert labels[30, 20] == labels[30, 42] != 0
    labels_far = c.identify_frame(_blob(shape, 30, 15, 6) + _blob(shape, 30, 65, 6),
                                  radius_km=12, grid_km=4)
    assert labels_far[30, 15] != labels_far[30, 65]


def test_labels_cover_exactly_the_wet_cells():
    rng = np.random.default_rng(3)
    rain = np.clip(ndimage.gaussian_filter(rng.standard_normal((120, 150)), 4) * 20, 0, None)
    labels = c.identify_frame(rain, radius_km=20, grid_km=4)
    assert np.array_equal(labels > 0, rain > 0.1 / 3)
    assert set(np.unique(labels)) == set(range(int(labels.max()) + 1))


@pytest.mark.parametrize("limit", [10 ** 9, 0])
def test_similarity_matches_bruteforce_for_direct_and_fft_paths(limit):
    rng = np.random.default_rng(4)
    shape = (70, 90)
    a_rain = _blob(shape, 30, 30, 9) * (1 + rng.random(shape))
    b_rain = _blob(shape, 36, 41, 7) * (1 + rng.random(shape))
    a = c.segments_from_labels((a_rain > 0).astype(int), a_rain)[0]
    b = c.segments_from_labels((b_rain > 0).astype(int), b_rain)[0]
    phi = 4 / 120
    expected = c.similarity_bruteforce(a, b, phi)
    assert c.similarity(a, b, phi, direct_limit=limit) == pytest.approx(expected, rel=1e-10)
    assert c.similarity(a, a, phi) <= 1.0


def test_bbox_pruning_is_exact():
    tracker = c.Chang2016Tracker(grid_km=4, tau=0.05)
    cutoff = tracker._prune_gap_cells()
    shape = (10, 400)
    a_rain, b_rain = np.zeros(shape), np.zeros(shape)
    a_rain[5, 5] = 1
    b_rain[5, 5 + int(math.ceil(cutoff))] = 1
    a = c.segments_from_labels((a_rain > 0).astype(int), a_rain)[0]
    b = c.segments_from_labels((b_rain > 0).astype(int), b_rain)[0]
    assert c.similarity_bruteforce(a, b, tracker.phi_per_cell) <= tracker.tau


def _moving_frames(T=5, shape=(80, 160), step=6, r=6):
    return np.stack([_blob(shape, 40, 20 + step * t, r, 3.0) for t in range(T)])


def test_moving_storm_keeps_one_event():
    rain = _moving_frames()
    labels = c.identify(rain, radius_km=12, grid_km=4)
    events, rows, _ = c.track(labels, rain, grid_km=4, tau=0.05)
    assert len({r["event_id"] for r in rows}) == 1
    assert [r["rule"] for r in rows] == ["new"] + ["distance"] * 4


def test_split_shares_event_and_merge_continues_largest():
    shape = (80, 120)
    t0 = _blob(shape, 40, 60, 10)
    t1 = _blob(shape, 25, 60, 6) + _blob(shape, 55, 60, 6)       # split
    labels = c.identify(np.stack([t0, t1]), radius_km=4, grid_km=4)
    events, rows, _ = c.track(labels, np.stack([t0, t1]), grid_km=4, tau=0.05)
    assert len({r["event_id"] for r in rows}) == 1
    big, small = _blob(shape, 40, 30, 10), _blob(shape, 40, 90, 5)
    merged = _blob(shape, 40, 60, 12)
    stack = np.stack([big + small, merged])
    labels = c.identify(stack, radius_km=4, grid_km=4)
    events, rows, _ = c.track(labels, stack, grid_km=4, tau=0.01)
    assert events[1, 40, 60] == events[0, 40, 30] != events[0, 40, 90]


def test_angle_rule_and_no_history_policy():
    shape = (40, 400)
    # 20 cells (80 km), then 40 cells (160 km) per step at 4 km spacing
    rain = np.stack([_blob(shape, 20, x, 8, 3.0) for x in (20, 40, 80)])
    labels = c.identify(rain, radius_km=4, grid_km=4)
    _, rows, _ = c.track(labels, rain, grid_km=4, tau=1e-6)
    assert [r["rule"] for r in rows] == ["new", "distance", "angle"]
    # a fast first step has no predecessor history: rejected by default
    fast = np.stack([_blob(shape, 20, x, 8, 3.0) for x in (20, 60, 100)])
    labels = c.identify(fast, radius_km=4, grid_km=4)
    _, rows, _ = c.track(labels, fast, grid_km=4, tau=1e-6)
    assert [r["rule"] for r in rows] == ["new", "new", "new"]
    _, rows, _ = c.track(labels, fast, grid_km=4, tau=1e-6, no_history="accept")
    assert [r["rule"] for r in rows] == ["new", "angle", "angle"]
    # reversing direction beyond 120 km is refused
    back = np.stack([_blob(shape, 20, x, 8, 3.0) for x in (200, 220, 180)])
    labels = c.identify(back, radius_km=4, grid_km=4)
    _, rows, _ = c.track(labels, back, grid_km=4, tau=1e-6)
    assert [r["rule"] for r in rows] == ["new", "distance", "new"]


def test_factorization_identity_and_units():
    rain = _moving_frames(T=4)
    labels = c.identify(rain, radius_km=12, grid_km=4)
    _, rows, _ = c.track(labels, rain, grid_km=4, tau=0.05)
    f = c.factorize(rows, cell_area_km2=16, dt_hours=1, n_time_steps=4)
    product = (f["intensity_mm_per_h"] * f["size_factor_km2"] *
               f["duration_factor_h"] * f["number_of_events"])
    wet = rain > 0.1 / 3
    assert product == pytest.approx(f["total_mm_km2"])
    assert f["total_mm_km2"] == pytest.approx(rain[wet].sum() * 16 * 1)
    assert f["intensity_mm_per_h"] == pytest.approx(rain[wet].mean())
    assert f["duration_factor_h"] == 4 and f["number_of_events"] == 1


def test_aggregate_rate_mean_and_rejects_partial_blocks():
    rain = np.arange(6 * 4, dtype=float).reshape(6, 2, 2)
    agg = c.aggregate_rate(rain, 3)
    assert agg.shape == (2, 2, 2)
    assert np.allclose(agg[0], rain[:3].mean(axis=0))
    with pytest.raises(ValueError):
        c.aggregate_rate(rain[:5], 3)


def test_shared_parent_field_matches_bruteforce():
    rng = np.random.default_rng(5)
    shape = (90, 120)
    parent_rain = _blob(shape, 45, 40, 15) * (1 + rng.random(shape))
    kids = [_blob(shape, 30, 70, 6) * (1 + rng.random(shape)),
            _blob(shape, 60, 75, 9) * (1 + rng.random(shape))]
    parent = c.segments_from_labels((parent_rain > 0).astype(int), parent_rain)[0]
    children = [c.segments_from_labels((k > 0).astype(int), k)[0] for k in kids]
    phi = 4 / 120
    got = c.similarities_to_parent(parent, children, phi, direct_limit=0)
    want = [c.similarity_bruteforce(ch, parent, phi) for ch in children]
    assert np.allclose(got, want, rtol=1e-10, atol=0)


def test_reference_runner_end_to_end(tmp_path):
    import json
    import subprocess
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    hourly = np.repeat(_moving_frames(T=4, step=2), 3, axis=0)  # 12 hourly frames
    np.save(tmp_path / "rain.npy", hourly)
    out = tmp_path / "out"
    subprocess.run([sys.executable, str(root / "scripts/run_chang2016_reference.py"),
                    str(tmp_path / "rain.npy"), "--preset",
                    str(root / "configs/chang2016_s1_summer_3h.json"),
                    "--grid-km", "4", "--output-dir", str(out)], check=True,
                   capture_output=True, text=True)
    summary = json.loads((out / "summary.json").read_text())
    assert summary["dt_hours"] == 3 and summary["frames"] == 4
    assert summary["factors"]["number_of_events"] == 1
    assert summary["factors"]["duration_factor_h"] == 12
    assert np.load(out / "event_labels.npy").shape == (4, 80, 160)


def test_selection_largest_versus_most_similar():
    shape = (60, 140)
    t0 = _blob(shape, 30, 40, 14) + _blob(shape, 30, 95, 6)
    t1 = _blob(shape, 30, 88, 6)                      # sits on the small one
    stack = np.stack([t0, t1])
    labels = c.identify(stack, radius_km=2, grid_km=2)
    ev, _, _ = c.track(labels, stack, grid_km=2, tau=1e-4)
    assert ev[1, 30, 88] == ev[0, 30, 40]             # largest qualifying
    ev, _, _ = c.track(labels, stack, grid_km=2, tau=1e-4, selection="most_similar")
    assert ev[1, 30, 88] == ev[0, 30, 95]


def _random_storms(seed, T=4, shape=(160, 200)):
    rng = np.random.default_rng(seed)
    base = ndimage.gaussian_filter(rng.standard_normal((shape[0] + 40, shape[1] + 40)), 5)
    return np.stack([np.clip((base[3 * t:3 * t + shape[0], 2 * t:2 * t + shape[1]] - 0.05)
                             * 300, 0, None) for t in range(T)])


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("tau", [0.05, 0.3, 0.7])
def test_certified_decisions_equal_exact_decisions(seed, tau):
    rain = _random_storms(seed)
    labels = c.identify(rain, radius_km=20, grid_km=4)
    fast, rows_fast, tf = c.track(labels, rain, grid_km=4, tau=tau)
    slow, rows_slow, ts = c.track(labels, rain, grid_km=4, tau=tau, exact_similarity=True)
    assert np.array_equal(fast, slow)
    assert [r["rule"] for r in rows_fast] == [r["rule"] for r in rows_slow]


def test_jensen_lower_bound_holds():
    rain = _random_storms(7, T=2)
    labels = c.identify(rain, radius_km=12, grid_km=4)
    a = c.segments_from_labels(labels[0], rain[0])
    b = c.segments_from_labels(labels[1], rain[1])
    phi = 4 / 120
    for s in a[:8]:
        for t in b[:8]:
            exact = c.similarity_bruteforce(s, t, phi)
            dc = math.hypot(s.centroid[0] - t.centroid[0], s.centroid[1] - t.centroid[1])
            assert exact >= math.exp(-phi * (dc + s.spread + t.spread)) - 1e-12
            assert exact <= math.exp(-phi * c._bbox_gap_cells(s, t)) + 1e-12


@pytest.mark.parametrize("k", [1, 4, 8])
def test_block_bounds_bracket_exact(k):
    rain = _random_storms(11, T=2)
    labels = c.identify(rain, radius_km=12, grid_km=4)
    a = c.segments_from_labels(labels[0], rain[0])
    b = c.segments_from_labels(labels[1], rain[1])
    phi = 4 / 120
    for s in a[:6]:
        for t in b[:6]:
            exact = c.similarity_bruteforce(s, t, phi)
            lo, hi = c.similarity_bounds(s, t, phi, k)
            assert lo - 1e-12 <= exact <= hi + 1e-12


@pytest.mark.parametrize("levels", [(), (4,), (32, 8)])
def test_block_certification_keeps_decisions(levels):
    rain = _random_storms(3)
    labels = c.identify(rain, radius_km=20, grid_km=4)
    ref, _, _ = c.track(labels, rain, grid_km=4, tau=0.3, exact_similarity=True)
    got, _, _ = c.track(labels, rain, grid_km=4, tau=0.3, block_levels=levels)
    assert np.array_equal(ref, got)


def test_compare_tracking_definitions_smoke(tmp_path):
    import json
    import subprocess
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    rain = _random_storms(2, T=6, shape=(120, 150)).astype(np.float32)
    np.savez_compressed(tmp_path / "rain.npz", rain_mm_h=rain)
    out = tmp_path / "cmp"
    subprocess.run([sys.executable, str(root / "scripts/compare_tracking_definitions.py"),
                    str(tmp_path / "rain.npz"), "--output-dir", str(out)],
                   check=True, capture_output=True, text=True)
    r = json.loads((out / "results.json").read_text())
    assert set(r["runs"]) == {"A", "B", "C1", "C2", "D", "D2"}
    assert r["runs"]["B"]["intensity_mm_h"] == pytest.approx(r["runs"]["C1"]["intensity_mm_h"])
