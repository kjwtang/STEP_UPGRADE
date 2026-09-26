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
    assert (out/'partition_0.1/objects.csv').exists()
    assert np.load(out/'core_labels.npy')[0].max()==2
    bounded=tmp_path/'bounded'
    subprocess.run([sys.executable,str(Path(__file__).resolve().parents[1]/'scripts/compare_envelopes.py'),
        str(source),'--output-dir',str(bounded),'--max-distance-km','1',
        '--spacing-km','1','1','--cell-area-km2','2','--min-core-area-km2','2'],
        check=True,capture_output=True,text=True)
    assert (bounded/'SUCCESS').exists()
    assert not np.load(bounded/'system_0.1/labels.npy')[0,3,5]
    import csv
    rows=list(csv.DictReader((bounded/'system_0.1/objects.csv').open()))
    assert float(rows[0]['core_area_km2'])==2


@pytest.mark.parametrize('mode',['partition','system'])
def test_distance_cap_is_physical_and_does_not_cross_missing(mode):
    rain=np.full((9,15),.2); rain[4,3]=2; rain[:,7]=np.nan
    cores,env=identify_envelope(rain,mode=mode,max_distance_km=4,spacing_km=(2,1))
    assert env[4,6]>0 and env[6,3]>0
    assert not env[7,3] and not env[:,8:].any()
    assert np.all(env[cores>0]>0)
    _,zero=identify_envelope(rain,mode=mode,max_distance_km=0,spacing_km=(2,1))
    assert np.array_equal(zero>0,cores>0)


def test_explicit_physical_area_not_just_pixel_count():
    rain=np.zeros((5,8)); rain[1,1:3]=2; rain[3,5:7]=2
    area=np.ones(rain.shape); area[3,5:7]=3
    cores,_=identify_envelope(rain,min_core_area_km2=5,cell_area_km2=area)
    assert not cores[1].any() and cores[3].max()==1
    with pytest.raises(ValueError):
        identify_envelope(rain,min_core_area_km2=5)
    with pytest.raises(ValueError):
        identify_envelope(rain,max_distance_km=4)
    with pytest.raises(ValueError):
        identify_envelope(rain,cell_area_km2=-1)


def test_statistics_mass_membership_and_censoring():
    from step.envelope_statistics import envelope_statistics
    rain=np.array([[2.,.2,0],[.1,.1,np.nan],[0,0,0]])
    core=np.array([[1,0,0],[0,0,0],[0,0,0]])
    env=np.array([[7,7,0],[7,7,0],[0,0,0]])
    rows,members=envelope_statistics(rain,core,env,cell_area_km2=4)
    r=rows[0]
    assert r['envelope_label']==7
    assert r['envelope_area_km2']==16 and r['core_area_km2']==4
    assert r['rain_volume_rate_m3_h']==pytest.approx(9600)
    assert r['noncore_rain_rate_fraction']==pytest.approx(1/6)
    assert r['touches_missing'] and r['touches_domain_boundary']
    assert members==[dict(envelope_label=7,core_label=1,core_cells=1)]
    noarea,_=envelope_statistics(rain,core,env)
    assert noarea[0]['rain_volume_rate_m3_h'] is None


def test_statistics_many_to_many_disconnected_seed():
    from step.envelope_statistics import envelope_statistics
    rain=np.array([[2.,0,2.]])
    cores=np.array([[1,0,1]])
    env=np.array([[1,0,2]])
    rows,members=envelope_statistics(rain,cores,env)
    assert len(rows)==len(members)==2
    assert {m['core_label'] for m in members}=={1}
