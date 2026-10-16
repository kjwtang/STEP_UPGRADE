import numpy as np
import pytest
from scipy import ndimage

from step import chang2016 as c
from step.identification import identify
from step.identification_regions import identify_regions, identify_regions_frame


def _blob(shape, cy, cx, r, value=2.0):
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    return np.where((yy - cy) ** 2 + (xx - cx) ** 2 <= r * r, value, 0.0)


def _chain_field():
    """Two large blobs 60 cells apart joined by specks every 4 cells."""
    shape = (60, 160)
    rain = _blob(shape, 30, 20, 10) + _blob(shape, 30, 140, 10)
    for x in range(33, 128, 4):
        rain[30, x] = 2.0
    return rain


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_defaults_match_chang2016_reference(seed):
    rng = np.random.default_rng(seed)
    rain = ndimage.gaussian_filter(rng.random((90, 110)), 2) * 4
    rain[rain < 1.6] = 0
    expected = c.identify_frame(rain, radius_km=4.0, grid_km=1.0, threshold=1.0, strict=False)
    assert np.array_equal(identify_regions_frame(rain, 4, 1.0), expected)


def test_preset_identification_chains_but_regions_do_not():
    rain = _chain_field()
    structure = c.disk(3)
    preset = identify(rain[None], morph_structure=structure, threshold=1.0)[0]
    assert preset[30, 20] == preset[30, 140] != 0          # one chained object
    labels = identify_regions_frame(rain, 3, 1.0)
    assert labels[30, 20] != labels[30, 140]                 # large blobs stay apart
    assert labels[30, 33] == labels[30, 20]                  # nearby speck joins its blob


def test_small_join_region_keeps_isolated_specks_apart():
    shape = (40, 120)
    rain = np.zeros(shape)
    for x in range(10, 110, 4):
        rain[20, x] = 2.0                                   # only small regions
    almost = identify_regions_frame(rain, 3, 1.0)
    assert len(np.unique(almost[almost > 0])) == 1
    separate = identify_regions_frame(rain, 3, 1.0, small_join="region")
    assert len(np.unique(separate[separate > 0])) == 25


def test_large_join_region_separates_nearby_large_regions():
    shape = (60, 100)
    rain = _blob(shape, 30, 30, 10) + _blob(shape, 30, 64, 10)   # gap of 13 cells < 2r
    joined = identify_regions_frame(rain, 8, 1.0)
    assert joined[30, 30] == joined[30, 64]
    apart = identify_regions_frame(rain, 8, 1.0, large_join="region")
    assert apart[30, 30] != apart[30, 64]


def test_threshold_is_inclusive_and_output_only_on_rain():
    rain = np.zeros((20, 20))
    rain[5, 5] = 1.0
    rain[10, 10] = 0.99
    labels = identify_regions_frame(rain, 2, 1.0)
    assert labels[5, 5] > 0 and labels[10, 10] == 0
    assert np.all(labels[rain < 1.0] == 0)


def test_stack_and_min_size():
    rain = np.stack([_chain_field(), np.zeros((60, 160))])
    labels = identify_regions(rain, 3, 1.0, min_size=2)
    assert labels.shape == rain.shape and labels[1].max() == 0
    assert labels[0][30, 128] == 0 or labels[0][30, 128] == labels[0][30, 140]


def test_small_core_radius_with_wide_bridge():
    rain = _chain_field()
    labels = identify_regions_frame(rain, 1, 1.0, bridge_radius_cells=3)
    assert labels[30, 20] != labels[30, 140]
    assert labels[30, 33] == labels[30, 20] and labels[30, 125] != labels[30, 140]  # 5-cell gap > 3
    narrow = identify_regions_frame(rain, 1, 1.0, bridge_radius_cells=1)
    assert narrow[30, 33] != narrow[30, 20]                 # 3-cell gap needs a bridge of 3



def test_edge_cells_removes_boundary_rain_only():
    from step.identification_regions import edge_mask
    rain = np.zeros((50, 60))
    rain[1:4, 5:55] = 2.0            # strip along the south edge
    rain[20:30, 20:30] = 2.0         # interior blob
    labels = identify_regions_frame(rain, 2, 1.0, edge_cells=5)
    assert labels[2, 30] == 0 and labels[25, 25] > 0
    assert np.array_equal(labels, identify_regions_frame(rain, 2, 1.0, valid_mask=edge_mask(rain.shape, 5)))
    assert identify_regions_frame(rain, 2, 1.0)[2, 30] > 0
