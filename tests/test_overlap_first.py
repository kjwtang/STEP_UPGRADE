import numpy as np
import pytest

from step.tracking import Tracker, track_with_graph, save_tracking_state, load_tracking_state


def morphing():
    lab = np.zeros((3, 12, 120), dtype=int)
    lab[:, 3:7, 5:15] = 1
    # Same identified object gains a distant lobe for one frame. Its weighted
    # centroid jumps, but the persistent core never goes away.
    lab[1, 3:7, 80:105] = 1
    return lab, (lab > 0).astype(float)*2


def test_morphology_rescue_keeps_id_and_resets_velocity():
    lab, rain = morphing()
    _, baseline = track_with_graph(lab, rain, km=5)
    assert not any(e.event == 'continue' for e in baseline.edges)
    audit = []
    raster, graph = track_with_graph(lab, rain, km=5, adjacent_policy='overlap_first',
                                     overlap_callback=audit.extend)
    assert [e.event for e in graph.edges] == ['continue','continue']
    assert len(set(raster[:,4,6])) == 1
    assert len([a for a in audit if a['action']=='reset_velocity']) == 2
    assert any(a.get('below_tau') for a in audit)


def test_split_merge_not_forced_into_one_to_one():
    lab = np.zeros((3,12,50),dtype=int)
    lab[0,3:7,5:35] = lab[2,3:7,5:35] = 1
    lab[1,3:7,5:17] = 1
    lab[1,3:7,23:35] = 2
    rain=(lab>0).astype(float)*2
    _, graph = track_with_graph(lab,rain,tau=.99,adjacent_policy='overlap_first')
    assert [e.event for e in graph.events] == ['split','merge']
    assert len(graph.edges)==4
    assert len({n.branch_id for n in graph.objects})==4
    assert len(set(graph.family_map.values()))==1


@pytest.mark.parametrize('case',['disjoint','tiny','sliver'])
def test_no_rescue_for_unrelated_or_insignificant_contacts(case):
    lab=np.zeros((2,20,60),dtype=int)
    lab[0,3:7,5:25]=1
    if case=='disjoint':
        lab[1,8:12,5:25]=1
    elif case=='tiny':
        lab[1,3:5,5:7]=1
    else:
        lab[1,3:7,24:44]=1
    _,graph=track_with_graph(lab,(lab>0).astype(float)*2,tau=.99,adjacent_policy='overlap_first')
    assert not graph.edges


def test_empty_frame_still_gets_gap():
    lab,rain=morphing()
    lab[1]=0
    rain[1]=0
    _,graph=track_with_graph(lab,rain,adjacent_policy='overlap_first',
                             gap_conflict_policy='endpoint_overlap')
    assert [e.event for e in graph.edges]==['gap_continue']


@pytest.mark.parametrize('reset',['auto','off','morphology'])
def test_chunk_resume_rasters_graph_and_diagnostics(tmp_path,reset):
    lab,rain=morphing()
    opts=dict(km=5,adjacent_policy='overlap_first',gap_conflict_policy='endpoint_overlap',velocity_reset=reset)
    audits=[]
    whole,gw=track_with_graph(lab,rain,overlap_callback=audits.extend,**opts)
    rasters=[]; edges=[]; nodes=[]; partial=[]; state=None
    for t in range(3):
        rr,gg,state=track_with_graph(lab[t:t+1],rain[t:t+1],state=state,
            return_state=True,overlap_callback=partial.extend,**opts)
        rasters.append(rr); edges.extend(gg.edges); nodes.extend(gg.objects)
        save_tracking_state(state,tmp_path/'state.json')
        state=load_tracking_state(tmp_path/'state.json')
    np.testing.assert_array_equal(whole,np.concatenate(rasters))
    assert gw.edges==edges and gw.objects==nodes and partial==audits
    assert gw.family_map==gg.family_map
    with pytest.raises(ValueError,match='parameters differ'):
        Tracker(state=state,max_displacement=5,gap_conflict_policy='endpoint_overlap')


def test_default_is_exactly_explicit_score():
    lab,rain=morphing()
    a,ga=track_with_graph(lab,rain)
    b,gb=track_with_graph(lab,rain,adjacent_policy='score')
    np.testing.assert_array_equal(a,b)
    assert ga==gb


def test_velocity_reset_ablation_keeps_overlap_rescue_and_legacy_config():
    lab,rain=morphing()
    outputs={}
    for reset in ('auto','morphology','off'):
        audit=[]
        raster,graph,state=track_with_graph(lab,rain,km=5,adjacent_policy='overlap_first',
            velocity_reset=reset,return_state=True,overlap_callback=audit.extend)
        outputs[reset]=(raster,graph,state,audit)
        assert [e.event for e in graph.edges]==['continue','continue']
        assert len(set(raster[:,4,6]))==1
    a,b=outputs['auto'],outputs['morphology']
    np.testing.assert_array_equal(a[0],b[0])
    assert a[1]==b[1] and a[2].tracking_config==b[2].tracking_config and a[3]==b[3]
    assert not any(r['action']=='reset_velocity' for r in outputs['off'][3])
    assert Tracker._velocity(outputs['off'][2].active[0])!=(0.,0.)
    assert Tracker._velocity(a[2].active[0])==(0.,0.)
    with pytest.raises(ValueError,match='parameters differ'):
        Tracker(state=a[2],max_displacement=5,adjacent_policy='overlap_first',velocity_reset='off')
    with pytest.raises(ValueError,match='parameters differ'):
        Tracker(state=outputs['off'][2],max_displacement=5,adjacent_policy='overlap_first')
