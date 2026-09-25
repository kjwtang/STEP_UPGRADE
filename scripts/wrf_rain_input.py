"""Read hourly 2-D CSTM RAINNC snapshots without loading 3-D fields.

Restart/reset intervals are rejected, never guessed or replaced with zero.
"""
from pathlib import Path
import re

import numpy as np
import xarray as xr


STAMP = re.compile(r"(\d{4}-\d{2}-\d{2})_(\d{2}:\d{2}:\d{2})\.nc$")


def manifest(directory):
    entries = []
    for path in Path(directory).rglob("cstm_d01_*.nc"):
        match = STAMP.search(path.name)
        if not match:
            raise ValueError(f"Unrecognized timestamp filename: {path}")
        entries.append((np.datetime64(match[1]+"T"+match[2],"s"),path))
    entries.sort(key=lambda item:(item[0],str(item[1])))
    if len(entries)<2:
        raise ValueError("Need at least two cstm_d01_*.nc snapshots for one rain interval")
    stamps = np.array([t for t,_ in entries])
    if np.any(np.diff(stamps)<=np.timedelta64(0,"s")):
        raise ValueError("Duplicate timestamps in directory; select one climate/member/month sequence")
    return entries


def internal_time(ds,path):
    if "Times" not in ds:
        raise ValueError(f"Missing Times: {path}")
    values = np.asarray(ds["Times"].values).ravel()
    text = "".join(v.decode("ascii") if isinstance(v,(bytes,np.bytes_)) else str(v) for v in values).strip("\x00 ")
    if len(text)!=19:
        raise ValueError(f"Expected one WRF Times timestamp: {path}: {text!r}")
    return np.datetime64(text.replace("_","T"),"s")


def read_wrf(a):
    entries = manifest(a.input)
    if a.inspect:
        print(f"Snapshots: {len(entries)}; possible intervals: {len(entries)-1}")
        print("First:",entries[0],"Last:",entries[-1])
        with xr.open_dataset(entries[0][1]) as ds:
            print("First internal Times:",internal_time(ds,entries[0][1]))
            print("RAINNC:",ds["RAINNC"].dims,ds["RAINNC"].shape,dict(ds["RAINNC"].attrs))
            print("Source attributes:",dict(ds.attrs))
        return None
    if a.variable not in (None,"RAINNC") or a.rain_kind != "cumulative":
        raise ValueError("WRF directory requires --wrf-cumulative --rain-kind cumulative (RAINNC only)")
    end = a.start+a.hours+1
    if a.start<0 or end>len(entries):
        raise ValueError(f"Need {a.hours+1} snapshots starting at index {a.start}; found {len(entries)}")
    selected = entries[a.start:end]
    times = np.array([t for t,_ in selected])
    intervals = np.diff(times)/np.timedelta64(1,"h")
    if not np.all(intervals==1):
        raise ValueError("CSTM hourly sequence has missing/irregular snapshots; cannot difference across missing hours")
    if a.dt_hours is not None and not np.isclose(a.dt_hours,1):
        raise ValueError("CSTM input is hourly; --dt-hours must be 1")
    previous = None
    reference_grid = None
    reference_shape = None
    source_refs = []
    warnings = []
    rain = None
    tiny_negatives = 0
    for index,(stamp,path) in enumerate(selected):
        with xr.open_dataset(path) as ds:
            actual = internal_time(ds,path)
            if actual!=stamp:
                raise ValueError(f"Filename/Times mismatch: {path.name} has internal Times={actual}; resolve provenance before proceeding")
            if "RAINNC" not in ds or ds["RAINNC"].dims != ("south_north","west_east"):
                raise ValueError(f"Expected 2-D RAINNC(south_north,west_east): {path}")
            units = a.units or ds["RAINNC"].attrs.get("units","")
            if units.strip().lower()!="mm":
                raise ValueError(f"RAINNC must be mm: {path}")
            shape = ds["RAINNC"].shape
            if reference_shape is not None and shape!=reference_shape:
                raise ValueError("Grid shape changes across files")
            reference_shape = shape
            width = a.crop or max(shape)
            y,x = max(0,(shape[0]-width)//2),max(0,(shape[1]-width)//2)
            selection = {"south_north":slice(y,y+width),"west_east":slice(x,x+width)}
            grid = []
            for coordinate in ("XLAT","XLONG"):
                if coordinate not in ds or ds[coordinate].dims != ("south_north","west_east"):
                    raise ValueError(f"Missing or invalid {coordinate}: {path}")
                grid.append(ds[coordinate].isel(selection).values)
            if reference_grid is None:
                reference_grid = grid
            elif not all(np.array_equal(v,r,equal_nan=True) for v,r in zip(grid,reference_grid)):
                raise ValueError(f"Coordinates change across snapshots: {path}")
            current = np.array(ds["RAINNC"].isel(selection).values,dtype=np.float32,copy=True)
            if np.any(current[np.isfinite(current)]<0):
                raise ValueError(f"Negative cumulative rain: {path}")
            source_ref = str(ds.attrs.get("wrfout_ref",""))
            source_refs.append(source_ref)
            year = str(stamp)[:4]
            if source_ref and year not in source_ref:
                warnings.append(f"{path.name}: wrfout_ref does not contain timestamp year {year}: {source_ref}")
        if previous is not None:
            valid = np.isfinite(current)&np.isfinite(previous)
            if not valid.any():
                raise ValueError(f"Missing whole rain interval ending {stamp}")
            increment = np.full(current.shape,np.nan,dtype=np.float32)
            np.subtract(current,previous,out=increment,where=valid)
            bad = valid & (increment < -a.negative_tolerance_mm)
            if bad.any():
                raise ValueError(f"Cumulative RAINNC decreased at {stamp}: {int(bad.sum())} pixels, min={float(np.nanmin(increment))} mm. Possible restart/bucket reset; supply correct preprocessing, do not clip.")
            tiny = valid & (increment<0)
            tiny_negatives += int(tiny.sum())
            increment[tiny] = 0
            if rain is None:
                rain = np.empty((a.hours,*current.shape),dtype=np.float32)
            rain[index-1] = increment
        previous = current
    metadata = dict(source_type="CSTM hourly cumulative RAINNC",source_shape=[len(entries),*reference_shape],
        crop_origin_yx=[y,x],dt_hours=1.,time_verified=True,normalized_units="mm/h",source_units="mm",
        timestamps=[str(t) for t in times[1:]],interval_start_times=[str(t) for t in times[:-1]],
        input_files=[str(p) for _,p in selected],wrfout_refs=source_refs,
        provenance_warnings=warnings,negative_tolerance_mm=a.negative_tolerance_mm,
        tiny_negative_pixels_clipped=tiny_negatives,
        precipitation_scope="Grid-scale RAINNC only; does not include RAINC if parameterized convective rain exists")
    if warnings:
        print("PROVENANCE WARNING:",warnings[0],flush=True)
    return rain,metadata
