import numpy as np
from step.identification import _stable_relabel


def test_vectorized_relabel_matches_first_rain_pixel_reference():
    rng=np.random.default_rng(8)
    for count in (0,1,4,100):
        values=np.arange(1,count+1)*3
        lab=rng.choice(np.concatenate([[0],values]),size=(40,50)).astype(int)
        present=[v for v in np.unique(lab) if v>0]
        ordered=sorted(present,key=lambda v:int(np.flatnonzero(lab==v)[0]))
        expected=np.zeros(lab.shape,dtype=np.int32)
        for new,old in enumerate(ordered,1):
            expected[lab==old]=new
        np.testing.assert_array_equal(_stable_relabel(lab),expected)
