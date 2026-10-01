import numpy as np
import pytest
from step.tracking import track_with_graph, save_tracking_state, load_tracking_state


OPTIONS=dict(adjacent_policy='overlap_balanced',velocity_reset='off',
             score_policy='coherent_adjacent',gap_conflict_policy='endpoint_overlap',gap_min_pixels=16)


def test_balanced_overlap_keeps_moderate_symmetric_overlap_despite_low_centroid_score():
    lab=np.zeros((2,20,100),dtype=int)
    lab[0,5:10,5:45]=1;lab[1,5:10,27:67]=1
    rain=(lab>0).astype(float)*2
    _,old=track_with_graph(lab,rain,tau=.99,adjacent_policy='overlap_first')
    _,new=track_with_graph(lab,rain,tau=.99,**OPTIONS)
    assert not old.edges
    assert len(new.edges)==1 and new.edges[0].event=='continue'


def test_balanced_event_does_not_attach_tiny_contained_fragment():
    lab=np.zeros((2,30,50),dtype=int)
    lab[0,5:15,5:35]=1
    lab[1,5:15,5:33]=1;lab[1,5,34]=2
    _,graph=track_with_graph(lab,(lab>0).astype(float)*2,tau=.99,**OPTIONS)
    assert not graph.events
    assert len(graph.edges)==1 and graph.edges[0].event=='continue'


def test_balanced_split_merge_gap_and_checkpoint(tmp_path):
    lab=np.zeros((5,30,50),dtype=int)
    lab[0,5:15,5:35]=1
    lab[1,5:15,5:19]=1;lab[1,5:15,21:35]=2
    lab[2,5:15,5:35]=1;lab[4,5:15,5:35]=1
    rain=(lab>0).astype(float)*2
    whole,graph=track_with_graph(lab,rain,**OPTIONS)
    assert [e.event for e in graph.events]==['split','merge']
    assert sum(e.gap>0 for e in graph.edges)==1
    first,g1,state=track_with_graph(lab[:3],rain[:3],return_state=True,**OPTIONS)
    save_tracking_state(state,tmp_path/'state.json')
    restored=load_tracking_state(tmp_path/'state.json')
    last,g2=track_with_graph(lab[3:],rain[3:],state=restored,**OPTIONS)
    np.testing.assert_array_equal(whole,np.concatenate([first,last]))
    assert graph.edges==g1.edges+g2.edges and graph.events==g1.events+g2.events
    with pytest.raises(ValueError,match='parameters differ'):
        track_with_graph(lab[3:],rain[3:],state=restored,adjacent_policy='overlap_first')


def test_gap_size_floor_keeps_observed_tiny_objects_but_not_gap_links():
    lab=np.zeros((3,10,10),dtype=int);lab[0,5,5]=lab[2,5,5]=1
    _,graph=track_with_graph(lab,(lab>0).astype(float)*2,**OPTIONS)
    assert len(graph.objects)==2 and not graph.edges


def test_four_way_split_is_not_misreported_as_a_continuation():
    lab=np.zeros((2,40,60),dtype=int);lab[0,10:20,5:49]=1
    for j in range(4):
        lab[1,10:20,5+j*11:15+j*11]=j+1
    _,graph=track_with_graph(lab,(lab>0).astype(float)*2,**OPTIONS)
    assert [e.event for e in graph.events]==['split']
    assert len(graph.edges)==4 and all(e.event=='split' for e in graph.edges)


def test_observed_family_scope_matches_global_roots_across_chunk_merge(tmp_path):
    from step.tracking import _find
    lab=np.zeros((4,30,60),dtype=int)
    lab[:2,5:15,5:19]=1;lab[:2,5:15,21:35]=2
    lab[2:,5:15,5:35]=1
    rain=(lab>0).astype(float)*2
    expected,g,state=track_with_graph(lab,rain,return_state=True,**OPTIONS)
    actual,h,sparse=track_with_graph(lab,rain,return_state=True,family_map_scope='observed',**OPTIONS)
    np.testing.assert_array_equal(expected,actual)
    assert h.edges==g.edges and h.events==g.events and h.objects==g.objects
    assert all(h.family_map[n.branch_id]==g.family_map[n.branch_id] for n in g.objects)
    a,ga,s=track_with_graph(lab[:2],rain[:2],return_state=True,family_map_scope='observed',**OPTIONS)
    save_tracking_state(s,tmp_path/'state.json')
    b,gb,last=track_with_graph(lab[2:],rain[2:],state=load_tracking_state(tmp_path/'state.json'),
                            return_state=True,family_map_scope='observed',**OPTIONS)
    np.testing.assert_array_equal(actual,np.concatenate([a,b]))
    assert h.edges==ga.edges+gb.edges
    for branch in state.family_parent:
        assert _find(last.family_parent,branch)==_find(state.family_parent,branch)
