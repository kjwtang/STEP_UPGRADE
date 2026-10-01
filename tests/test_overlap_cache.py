import numpy as np

from step.tracking import Tracker, _objects, _pixel_overlap, _set_overlap


def test_cached_overlap_matches_reference_for_translation_and_disjoint_masks():
    rng = np.random.default_rng(42)
    for shift in ((0, 0), (3, -7), (-4, 2), (100, 100)):
        a = (rng.random((20, 30)) > .7).astype(int)
        b = (rng.random((20, 30)) > .8).astype(int)
        parent = _objects(a, a.astype(float))[0]
        child = _objects(b, b.astype(float))[0]
        moved = {(y + shift[0], x + shift[1]) for y, x in parent.pixels}
        assert _set_overlap(moved, set(child.pixels), parent.area, child.area) == \
            _pixel_overlap(parent, child, shift)


def test_candidate_cache_does_not_survive_mask_change():
    a = np.zeros((20, 30), dtype=int)
    a[3:8, 5:10] = 1
    tracker = Tracker(max_displacement=100)
    tracker.update(0, a, a.astype(float))
    children = _objects(a, a.astype(float))
    assert tracker._candidates(tracker.state.active, children, 1)[0].raw_iou == 1
    children[0].pixels = [(y, x + 10) for y, x in children[0].pixels]
    assert tracker._candidates(tracker.state.active, children, 1)[0].raw_iou == 0
