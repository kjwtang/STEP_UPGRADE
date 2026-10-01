import sys
from pathlib import Path
import numpy as np
import pytest
xr=pytest.importorskip('xarray')
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare_rain_only import prepare


def test_filename_only_rain_and_reset_rejection(tmp_path):
    attrs={k:1. for k in ('DX','DY','MAP_PROJ','CEN_LAT','CEN_LON','TRUELAT1',
                          'TRUELAT2','STAND_LON','MOAD_CEN_LAT','POLE_LAT','POLE_LON')}
    for i,value in enumerate((10,12,15)):
        ds=xr.Dataset({'RAINNC':(('Time','south_north','west_east'),
            np.full((1,3,4),value,dtype=np.float32),{'units':'mm'})},attrs=attrs)
        ds.to_netcdf(tmp_path/f'cstm_d01_1996-01-02_0{i}:00:00.rain')
    meta=prepare(tmp_path,'1996-01-02T00:00:00',2,tmp_path/'prepared')
    assert not meta['internal_time_verified'] and not meta['coordinates_verified']
    assert meta['filename_time_verified'] and meta['cadence_hours']==1
    result=np.load(tmp_path/'prepared/rain_mm_h.npy')
    assert np.all(result[0]==2) and np.all(result[1]==3)
    with pytest.raises(FileExistsError):
        prepare(tmp_path,'1996-01-02T00:00:00',2,tmp_path/'prepared')
    ds['RAINNC'].values[:]=1
    ds.to_netcdf(tmp_path/'cstm_d01_1996-01-02_02:00:00.rain')
    with pytest.raises(ValueError,match='decreased'):
        prepare(tmp_path,'1996-01-02T00:00:00',2,tmp_path/'bad')
    assert not (tmp_path/'bad').exists()
