import numpy as np
from step.tracking import track_with_graph


OPTIONS=dict(adjacent_policy='overlap_ranked',velocity_reset='off',score_policy='coherent_adjacent')


def test_ranked_isolated_weak_overlap_continues_without_lowering_score_tau():
    lab=np.zeros((2,20,100),dtype=int)
    lab[0,5:10,5:45]=1;lab[1,5:10,37:77]=1
    rain=(lab>0).astype(float)*2
    _,balanced=track_with_graph(lab,rain,tau=.99,adjacent_policy='overlap_balanced')
    raster,ranked=track_with_graph(lab,rain,tau=.99,**OPTIONS)
    assert not balanced.edges
    assert [e.event for e in ranked.edges]==['continue']
    assert raster[0,5,5]==raster[1,5,37]


def test_ranked_ambiguous_equal_overlaps_do_not_invent_a_continuation():
    lab=np.zeros((2,20,100),dtype=int);lab[0,5:10,25:65]=1
    lab[1,5:10,0:33]=1;lab[1,5:10,57:90]=2
    _,graph=track_with_graph(lab,(lab>0).astype(float)*2,tau=.99,**OPTIONS)
    assert not graph.edges and not graph.events


def test_ranked_keeps_significant_four_way_event_priority():
    lab=np.zeros((2,30,60),dtype=int);lab[0,10:20,5:49]=1
    for j in range(4):lab[1,10:20,5+j*11:15+j*11]=j+1
    _,graph=track_with_graph(lab,(lab>0).astype(float)*2,**OPTIONS)
    assert [e.event for e in graph.events]==['split']
    assert len(graph.edges)==4 and all(e.event=='split' for e in graph.edges)
