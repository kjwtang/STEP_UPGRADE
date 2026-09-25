"""Data contract and end-to-end diagnostics, including streaming final short chunk."""
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

xr = pytest.importorskip("xarray")
pytest.importorskip("matplotlib")
pytest.importorskip("psutil")
SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from validate_real_data import arguments, load_input


def make_source(tmp_path, irregular=False):
    rain = np.zeros((5,16,16),dtype=np.float32)
    for t in range(5):
        rain[t,4:7,3+t:6+t] = 4
    rain[:,0,0] = np.nan
    times = np.datetime64("2005-07-01") + np.arange(5).astype("timedelta64[h]")
    if irregular:
        times[3:] += np.timedelta64(1,"h")
    ds = xr.Dataset({"PRECIP": (("time","y","x"),rain,{"units":"mm"})},coords={"time":times})
    path = tmp_path / "rain.nc"
    ds.to_netcdf(path)
    return path


def options(monkeypatch,path,*extra):
    monkeypatch.setattr(sys,"argv",["validate",str(path),"--variable","PRECIP","--rain-kind","interval","--hours","5",*extra])
    return arguments()


def test_loader_units_missing_and_full_domain(tmp_path,monkeypatch):
    path = make_source(tmp_path)
    a = options(monkeypatch,path,"--crop","0")
    data,meta = load_input(a)
    assert data.shape==(5,16,16)
    assert np.isnan(data[:,0,0]).all()
    assert meta["dt_hours"]==1 and meta["time_verified"]
    a.dt_hours = 2
    with pytest.raises(ValueError,match="disagrees"):
        load_input(a)


def test_loader_rejects_irregular_and_cumulative(tmp_path,monkeypatch):
    path = make_source(tmp_path,irregular=True)
    a = options(monkeypatch,path)
    with pytest.raises(ValueError,match="irregular"):
        load_input(a)
    with xr.open_dataset(path) as ds:
        renamed = ds.rename({"PRECIP":"RAINNC"}).load()
    raw = tmp_path / "raw.nc"
    renamed.to_netcdf(raw)
    a.input, a.variable = raw,"RAINNC"
    with pytest.raises(ValueError,match="Cumulative"):
        load_input(a)


@pytest.mark.parametrize("stream",[False,True])
def test_diagnostic_cli(tmp_path,stream):
    source = make_source(tmp_path)
    out = tmp_path / "results"
    args = [sys.executable,str(SCRIPTS / "validate_real_data.py"),str(source),
            "--variable","PRECIP","--rain-kind","interval","--hours","5","--crop","0",
            "--bridge-radius","0","--chunk-frames","2","--output-dir",str(out)]
    if stream:
        args += ["--stream","--save-labels"]
    else:
        args += ["--gif"]
    subprocess.run(args,check=True,capture_output=True,text=True,timeout=120)
    summary = json.loads((out / "summary.json").read_text())
    assert (out / "SUCCESS").exists()
    assert (out / "statistics.png").stat().st_size>1000
    if stream:
        assert summary["nodes"]==5
        assert summary["stage_seconds"]["tracking_seconds"]>0
        assert len(list(out.glob("chunk_*")))==3
        subprocess.run([sys.executable,str(SCRIPTS / "finalize_catalog.py"),
                        *map(str,sorted(out.glob("chunk_*"))),"--output-dir",str(tmp_path / "catalog")],
                       check=True,capture_output=True,text=True)
    else:
        assert all(summary["checks"].values())
        assert (out / "tracking.gif").exists()


def wrf_source(tmp_path, issue=None):
    directory = tmp_path / "wrf"
    directory.mkdir()
    yy,xx = np.indices((12,12),dtype=np.float32)
    for i in range(6):
        stamp = f"2005-07-01_{i:02d}:00:00"
        internal = stamp if issue!="time" or i!=2 else "2012-08-01_02:00:00"
        rain = np.zeros((12,12),dtype=np.float32)
        rain[3:6,4:7] = 2*i
        if issue=="reset" and i==3:
            rain[:] = 0
        if issue=="missing" and i==3:
            continue
        coords_x = xx+1 if issue=="grid" and i==3 else xx
        ds = xr.Dataset({"RAINNC":(("south_north","west_east"),rain,{"units":"mm"}),
                         "Times":(("DateStrLen",),np.array(list(internal),dtype="S1")),
                         "XLAT":(("south_north","west_east"),yy),
                         "XLONG":(("south_north","west_east"),coords_x)})
        ds.to_netcdf(directory / f"cstm_d01_{stamp}.nc")
    return directory


def test_wrf_difference_and_chunk_boundary(tmp_path,monkeypatch):
    source = wrf_source(tmp_path)
    a = options(monkeypatch,source,"--wrf-cumulative","--rain-kind","cumulative","--variable","RAINNC")
    whole,meta = load_input(a)
    assert whole.shape==(5,12,12)
    assert np.all(whole[:,3:6,4:7]==2)
    assert len(meta["input_files"])==6
    a.hours = 2
    first,_ = load_input(a)
    a.start,a.hours = 2,3
    second,_ = load_input(a)
    np.testing.assert_array_equal(whole,np.concatenate([first,second]))
    out = tmp_path / "stream"
    subprocess.run([sys.executable,str(SCRIPTS / "validate_real_data.py"),str(source),
        "--wrf-cumulative","--rain-kind","cumulative","--hours","5","--chunk-frames","2",
        "--bridge-radius","0","--stream","--output-dir",str(out)],
        check=True,capture_output=True,text=True,timeout=120)
    assert (out / "SUCCESS").exists()


@pytest.mark.parametrize("issue,match",[("time","Filename/Times"),("reset","decreased"),
                                      ("missing","missing/irregular"),("grid","Coordinates change")])
def test_wrf_rejects_bad_sequences(tmp_path,monkeypatch,issue,match):
    source = wrf_source(tmp_path,issue)
    a = options(monkeypatch,source,"--wrf-cumulative","--rain-kind","cumulative","--variable","RAINNC")
    if issue=="missing":
        a.hours=4
    with pytest.raises(ValueError,match=match):
        load_input(a)


def test_partition_metrics_label_permutation_and_split():
    from compare_original_step import partition_metrics
    a=np.array([[0,1,1],[2,2,0]])
    b=np.array([[0,9,9],[4,4,0]])
    assert partition_metrics(a,b)["ari_common_wet"]==1
    b[0,2]=8
    assert partition_metrics(a,b)["ari_common_wet"]<1


def test_original_comparison_optional(tmp_path):
    import os
    reference=os.environ.get("STEP_ORIGINAL_TEST_DIR")
    if not reference:
        pytest.skip("Set STEP_ORIGINAL_TEST_DIR to the pinned original checkout for integration test")
    source=wrf_source(tmp_path)
    out=tmp_path/"comparison"
    subprocess.run([sys.executable,str(SCRIPTS/"compare_original_step.py"),str(source),
        "--original-dir",reference,"--output-dir",str(out),"--hours","5","--crop","12",
        "--radius","0","--workers","1","--plot-frames","2","--timeout-seconds","30"],
        check=True,capture_output=True,text=True,timeout=120)
    assert (out/"SUCCESS").exists()
    assert (out/"tracking_001.png").stat().st_size>1000
    assert all(v["status"]=="completed" for v in json.loads((out/"stage_results.json").read_text()).values())


def test_comparison_timeout_preserves_status(tmp_path):
    from compare_original_step import run_stage
    result=run_stage("new_id",dict(timeout_seconds=1e-9,memory_gb=8),tmp_path)
    assert result["status"]=="timeout"
    assert (tmp_path/"new_id/status.json").exists()
