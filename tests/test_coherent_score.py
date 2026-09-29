import numpy as np
import pytest
from step.tracking import coherent_score, Tracker, track_with_graph, save_tracking_state, load_tracking_state


def test_bad_prediction_cannot_penalize_stationary_hypothesis():
    expected=.6*.31+.25*np.exp(-3.6/20)+.15*.9
    assert coherent_score(.31,0,3.6,75,20,.9)==pytest.approx(expected)


def test_no_mixing_raw_overlap_with_predicted_distance():
    raw,adv,actual,predicted=.1148409894,.0241131463,65.16,10.90
    mixed=.6*raw+.25*np.exp(-predicted/20)+.15*.9
    expected=max(.6*raw+.25*np.exp(-actual/20)+.15*.9,
                 .6*adv+.25*np.exp(-predicted/20)+.15*.9)
    assert coherent_score(raw,adv,actual,predicted,20,.9)==pytest.approx(expected)
    assert expected<mixed


def test_valid_translation_still_uses_motion_hypothesis():
    expected=.6*.8+.25*np.exp(-1/20)+.15
    assert coherent_score(.1,.8,30,1,20,1)==pytest.approx(expected)


def test_checkpoint_and_gap_scoring_unchanged(tmp_path):
    lab=np.zeros((4,12,30),dtype=int)
    lab[0,3:7,4:10]=1
    lab[1,3:7,6:12]=1
    lab[3,3:7,8:14]=1
    rain=(lab>0).astype(float)*2
    opts=dict(score_policy='coherent_adjacent',adjacent_policy='overlap_first',velocity_reset='off')
    whole,g=track_with_graph(lab,rain,**opts)
    a,ga,state=track_with_graph(lab[:2],rain[:2],return_state=True,**opts)
    save_tracking_state(state,tmp_path/'state.json')
    restored=load_tracking_state(tmp_path/'state.json')
    b,gb=track_with_graph(lab[2:],rain[2:],state=restored,**opts)
    np.testing.assert_array_equal(whole,np.concatenate([a,b]))
    assert g.edges==ga.edges+gb.edges and g.objects==ga.objects+gb.objects
    assert any(e.gap for e in g.edges)
    # Same state and same pair: multi-frame candidate scoring remains legacy.
    from step.tracking import _objects
    children=_objects(lab[3],rain[3])
    coh=Tracker(state=load_tracking_state(tmp_path/'state.json'),**opts)
    old=Tracker()
    assert coh._candidates(coh.state.active,children,3)==old._candidates(coh.state.active,children,3)
    with pytest.raises(ValueError,match='parameters differ'):
        Tracker(state=load_tracking_state(tmp_path/'state.json'),adjacent_policy='overlap_first',velocity_reset='off')
