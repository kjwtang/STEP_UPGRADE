import csv
from pathlib import Path
import numpy as np
import pytest
from step.tracking import track_with_graph,save_tracking_state


def test_final_family_export_retains_branch_and_asof_family_ids(tmp_path):
    scripts=Path(__file__).resolve().parents[1]/'scripts'
    import sys
    sys.path.insert(0,str(scripts))
    from run_npy_validation import write_graph
    from canonicalize_stream_families import export
    labels=np.zeros((3,20,30),dtype=int)
    labels[0,3:8,3:10]=1;labels[0,3:8,12:19]=2
    labels[1:,3:8,3:19]=1
    rain=(labels>0).astype(float)*2
    source=tmp_path/'stream';source.mkdir()
    first,g,s=track_with_graph(labels[:1],rain[:1],return_state=True,
                              adjacent_policy='overlap_balanced',family_map_scope='observed')
    p=source/'chunk_000000';p.mkdir();write_graph(g,p);save_tracking_state(s,p/'state.json')
    before=(p/'objects.csv').read_bytes()
    _,h,s=track_with_graph(labels[1:],rain[1:],state=s,return_state=True,
                          adjacent_policy='overlap_balanced',family_map_scope='observed')
    p=source/'chunk_000001';p.mkdir();write_graph(h,p);save_tracking_state(s,p/'state.json')
    (source/'SUCCESS').write_text('done')
    summary=export(source,tmp_path/'canonical')
    assert summary['objects']==4 and summary['families']==1
    rows=list(csv.DictReader((tmp_path/'canonical/objects.csv').open()))
    assert [int(r['branch_id']) for r in rows[:2]]==[1,2]
    assert [int(r['family_id']) for r in rows[:2]]==[1,2]
    assert {int(r['canonical_family_id']) for r in rows}=={1}
    assert (source/'chunk_000000/objects.csv').read_bytes()==before
    with pytest.raises(FileExistsError):export(source,tmp_path/'canonical')
