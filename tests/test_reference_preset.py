import json
from pathlib import Path
import numpy as np
from step.identification import identify
from step.tracking import Tracker,track_with_graph


def test_frozen_reference_preset_is_executable_and_keeps_gap_identity():
    root=Path(__file__).resolve().parents[1]
    preset=json.loads((root/'configs/rainfall_lineage_4km_1h_reference_v1.json').read_text())
    assert preset['input_contract']['grid_km']==4
    assert preset['input_contract']['dt_hours']==1
    radius=preset['identification']['bridge_radius_cells']
    y,x=np.ogrid[-radius:radius+1,-radius:radius+1]
    rain=np.zeros((4,60,60),dtype=np.float32)
    rain[0,20:30,20:30]=2;rain[1,20:30,21:31]=2;rain[3,20:30,21:31]=2
    labels=identify(rain,x*x+y*y<=radius*radius,
        threshold=preset['identification']['threshold_mm_h'],
        min_size=preset['identification']['min_size_cells'])
    Tracker(**preset['tracking'])
    kwargs=dict(preset['tracking'])
    kwargs['km']=kwargs.pop('max_displacement')
    tracked,graph,state=track_with_graph(labels,rain,return_state=True,**kwargs)
    assert tracked[0,25,25]==tracked[1,25,25]==tracked[3,25,25]>0
    assert [e.event for e in graph.edges]==['continue','gap_continue']
    assert state.tracking_config['adjacent_overlap_ranked_v2']==1
