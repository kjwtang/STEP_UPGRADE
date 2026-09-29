import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from replay_stream_pairs import replay
from compare_stream_tracks import compare
from run_npy_validation import write_graph
from step.tracking import track_with_graph, save_tracking_state


def test_reference_reports_partial_recovery(tmp_path):
    from dataclasses import replace
    lab=np.ones((3,4,4),dtype=int)
    _,graph=track_with_graph(lab,lab.astype(float))
    assert len(graph.edges)==2
    for name,edges in [('reference',graph.edges),('before',[]),('after',graph.edges[:1])]:
        part=tmp_path/name/'chunk_000000'
        part.mkdir(parents=True)
        write_graph(replace(graph,edges=edges),part)
        np.save(part/'identified_labels.npy',lab)
    result=compare(tmp_path/'before',tmp_path/'after',tmp_path/'report',tmp_path/'reference')
    recovery=result['reference_recovery']
    assert recovery['previously_missing_continue']==2
    assert len(recovery['recovered_pairs'])==1
    assert len(recovery['still_missing_pairs'])==1
    assert not recovery['newly_missing_pairs']


@pytest.mark.parametrize('guard',['off','endpoint_overlap'])
@pytest.mark.parametrize('adjacent',['score','overlap_first'])
@pytest.mark.parametrize('reset',['auto','off'])
def test_replay_saved_statistics_exact_and_reject_changed_mask(tmp_path,guard,adjacent,reset):
    lab=np.zeros((3,20,20),dtype=int)
    lab[:,5:10,5:10]=1
    lab[1,2:14,2:14]=1
    rain=(lab>0).astype(float)*2
    # Nonuniform intensity exercises saved weighted centroid reconstruction.
    rain[:,6,6]=7
    source=tmp_path/'source'
    state=None
    for i in range(3):
        part=source/f'chunk_{i:06d}'
        part.mkdir(parents=True)
        tracked,graph,state=track_with_graph(lab[i:i+1],rain[i:i+1],state=state,
            tau=.99,return_state=True,gap_conflict_policy=guard,adjacent_policy=adjacent,velocity_reset=reset)
        write_graph(graph,part)
        save_tracking_state(state,part/'state.json')
        np.save(part/'identified_labels.npy',lab[i:i+1])
        np.save(part/'tracked_labels.npy',tracked)
    report=replay(source,tmp_path/'out',[(1,2),(2,3)])
    assert report['replay_equal']
    assert report['frames_replayed']==3
    assert all(r['decision']==('accepted' if adjacent=='overlap_first' else 'below_score_threshold')
               for r in report['pairs'])
    assert bool(report['gap_conflicts'])==(guard!='off' and adjacent=='score')
    comparison=compare(source,source,tmp_path/'comparison',reference=source)
    assert comparison['added']==comparison['removed']==comparison['changed_event']==0
    assert comparison['reference_recovery']['previously_missing_continue']==0
    np.save(source/'chunk_000002'/'tracked_labels.npy',np.zeros_like(lab[:1]))
    with pytest.raises(ValueError,match='Raster replay mismatch'):
        replay(source,tmp_path/'bad',[(1,2)])
    assert not (tmp_path/'bad').exists()
