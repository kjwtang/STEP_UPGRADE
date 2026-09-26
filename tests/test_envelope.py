import numpy as np
import pytest
from step.envelope import identify_envelope

pytest.importorskip('skimage', reason='Install the envelope extra for experimental tests')


@pytest.mark.parametrize('mode',['system','partition'])
def test_unseeded_weak_rain_and_missing_barrier(mode):
    rain=np.full((7,12),.2)
    rain[:,6]=np.nan
    rain[3,2]=2
    core,env=identify_envelope(rain,mode=mode)
    assert np.all(env[:,:6]>0)
    assert not env[:,6:].any()
    assert np.all(env[core>0]>0)


def test_weak_bridge_has_two_explicit_interpretations():
    rain=np.zeros((9,15))
    rain[4,2:13]=.1
    rain[4,2:4]=2
    rain[4,11:13]=2
    core,partition=identify_envelope(rain)
    _,system=identify_envelope(rain,mode='system')
    assert core.max()==partition.max()==2
    assert system.max()==1
    assert np.array_equal(partition[core>0],core[core>0])
    assert np.array_equal(partition>0,system>0)


def test_seed_size_threshold_and_determinism():
    rain=np.full((6,6),.1)
    rain[2,2]=1
    assert not identify_envelope(rain,min_core_cells=2)[1].any()
    assert np.array_equal(identify_envelope(rain)[1],identify_envelope(rain)[1])
    with pytest.raises(ValueError):
        identify_envelope(rain,envelope_threshold=2)


def test_connectivity_and_dry_barrier():
    rain=np.zeros((5,5)); rain[1,1]=2; rain[2,2]=.1
    assert identify_envelope(rain,connectivity=4)[1][2,2]==0
    assert identify_envelope(rain,connectivity=8)[1][2,2]>0
    assert not identify_envelope(rain)[1][rain==0].any()


def test_saved_cube_experiment_cli(tmp_path):
    import hashlib
    import json
    from pathlib import Path
    import subprocess
    import sys
    from step.identification import identify
    rain=np.zeros((2,8,12),dtype=np.float32)
    rain[0,3,2:10]=.2; rain[0,3,2]=2; rain[0,3,9]=2
    source=tmp_path/'input'; (source/'new_id').mkdir(parents=True)
    np.save(source/'rain_mm_h.npy',rain)
    np.save(source/'new_id/labels.npy',identify(rain,np.ones((1,1)),threshold=1))
    (source/'configuration.json').write_text(json.dumps(dict(radius=0,threshold=1,
        input_sha256=hashlib.sha256(rain.tobytes()).hexdigest())))
    out=tmp_path/'out'
    subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'scripts/compare_envelopes.py'),
        str(source),'--output-dir',str(out)],check=True,capture_output=True,text=True)
    assert (out/'SUCCESS').exists()
    assert json.loads((out/'summary.json').read_text())['baseline_seeds_verified']
    assert np.load(out/'system_0.1/labels.npy')[0].max()==1
    assert np.load(out/'partition_0.1/labels.npy')[0].max()==2
    assert not np.load(out/'partition_0.1/labels.npy')[1].any()
