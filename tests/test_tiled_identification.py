import numpy as np
import pytest
from step.identification import identify


@pytest.mark.parametrize('structure',[np.ones((7,7)),np.ones((4,6)),np.eye(9)])
def test_halo_dilation_and_global_connectivity_match_frame_parallel(structure):
    rng=np.random.default_rng(15)
    rain=np.where(rng.random((3,63,81))>.93,2.,0.)
    rain[:,30:34,5:75]=2
    rain[:,10:55,40]=np.nan
    rain[:,52:55,40]=0 # distant opening: connectivity must be global
    expected=identify(rain,structure,threshold=1,min_size=4)
    tiled=identify(rain,structure,workers=4,threshold=1,min_size=4,parallel_mode='tiles')
    framed=identify(rain,structure,workers=2,threshold=1,min_size=4)
    np.testing.assert_array_equal(expected,tiled)
    np.testing.assert_array_equal(expected,framed)
    np.testing.assert_array_equal(expected[:1],identify(rain[:1],structure,workers=8,
                                 threshold=1,min_size=4,parallel_mode='tiles'))
