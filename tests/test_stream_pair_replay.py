import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from replay_stream_pairs import replay
from run_npy_validation import write_graph
from step.tracking import track_with_graph, save_tracking_state


@pytest.mark.parametrize('guard',['off','endpoint_overlap'])
def test_replay_saved_statistics_exact_and_reject_changed_mask(tmp_path,guard):
    lab=np.zeros((3,20,20),dtype=int)
    lab[:,5:10,5:10]=1
    lab[1,2:15,2:15]=1
    rain=(lab>0).astype(float)*2
    # Nonuniform intensity exercises saved weighted centroid reconstruction.
    rain[:,6,6]=7
    source=tmp_path/'source'
    state=None
    for i in range(3):
        part=source/f'chunk_{i:06d}'
        part.mkdir(parents=True)
        tracked,graph,state=track_with_graph(lab[i:i+1],rain[i:i+1],state=state,
            tau=.99,return_state=True,gap_conflict_policy=guard)
        write_graph(graph,part)
        save_tracking_state(state,part/'state.json')
        np.save(part/'identified_labels.npy',lab[i:i+1])
        np.save(part/'tracked_labels.npy',tracked)
    report=replay(source,tmp_path/'out',[(1,2),(2,3)])
    assert report['replay_equal']
    assert report['frames_replayed']==3
    assert all(r['decision']=='below_score_threshold' for r in report['pairs'])
    assert bool(report['gap_conflicts'])==(guard!='off')
    np.save(source/'chunk_000002'/'tracked_labels.npy',np.zeros_like(lab[:1]))
    with pytest.raises(ValueError,match='Raster replay mismatch'):
        replay(source,tmp_path/'bad',[(1,2)])
    assert not (tmp_path/'bad').exists()
