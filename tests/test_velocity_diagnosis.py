import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from diagnose_velocity_ablation import diagnose, select_pairs
from run_npy_validation import write_graph
from step.tracking import track_with_graph, save_tracking_state


def test_selects_lost_controls_and_new_event_types():
    a,b,c=(0,1,1,1),(1,1,2,1),(2,1,3,1)
    old={a:'continue',c:'continue'}
    new={b:'continue',c:'split'}
    ref={a:'continue',b:'continue',c:'continue'}
    keys,groups=select_pairs(old,new,ref)
    assert keys==[a,b,c]
    assert groups['newly_lost_reference']==[a]
    assert groups['recovered_controls']==[b]
    assert groups['new_event_pairs']==[c]


def test_paired_replay_reports_motion_difference_and_checks_policy(tmp_path):
    lab=np.zeros((3,10,80),dtype=int)
    for t,x in enumerate((5,20,40)):
        lab[t,3:7,x:x+20]=1
    rain=(lab>0).astype(float)*2
    for name,reset in [('on','morphology'),('off','off')]:
        root=tmp_path/name
        state=None
        for t in range(3):
            part=root/f'chunk_{t:06d}'
            part.mkdir(parents=True)
            raster,graph,state=track_with_graph(lab[t:t+1],rain[t:t+1],state=state,
                return_state=True,tau=.2,km=5,adjacent_policy='overlap_first',velocity_reset=reset)
            write_graph(graph,part)
            save_tracking_state(state,part/'state.json')
            np.save(part/'identified_labels.npy',lab[t:t+1])
            np.save(part/'tracked_labels.npy',raster)
    out=tmp_path/'audit'
    result=diagnose(tmp_path/'on',tmp_path/'off',tmp_path/'off',out)
    assert result['both_replays_equal']
    assert result['selected_pairs']==1
    pair=result['pairs'][0]
    assert pair['reset_on']['accepted_edges']==[]
    assert pair['reset_off']['accepted_edges'][0]['event']=='continue'
    assert pair['reset_on']['velocity_x_cells_per_frame']==0
    assert pair['reset_off']['velocity_x_cells_per_frame']==15
    assert (out/'SUCCESS').exists()
    with pytest.raises(ValueError,match='reset-enabled'):
        diagnose(tmp_path/'off',tmp_path/'on',tmp_path/'off',tmp_path/'bad')
    assert not (tmp_path/'bad').exists()
