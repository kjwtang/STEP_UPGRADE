"""Bounded hourly precipitation reads across single- or multi-time NetCDF files.

Classification uses declared metadata and precipitation names, never numerical
monotonicity. Index-only files have an explicitly assumed, not verified, timeline.
"""
from pathlib import Path
import re

import numpy as np
import xarray as xr


GRID_KEYS = ('DX', 'DY', 'MAP_PROJ', 'CEN_LAT', 'CEN_LON', 'TRUELAT1',
             'TRUELAT2', 'STAND_LON', 'MOAD_CEN_LAT', 'POLE_LAT', 'POLE_LON')
AMOUNT_UNITS = {'mm': 1., 'millimeter': 1., 'millimeters': 1.,
                'kg m-2': 1., 'kg m^-2': 1., 'kg m**-2': 1., 'kg/m2': 1., 'm': 1000.}
RATE_UNITS = {'mm/h': 1., 'mm/hr': 1., 'mm h-1': 1., 'mm hr-1': 1.,
              'mm h^-1': 1., 'mm h**-1': 1., 'mm hour-1': 1.,
              'mm/s': 3600., 'kg m-2 s-1': 3600., 'kg/m2/s': 3600.,
              'kg m^-2 s^-1': 3600., 'kg m**-2 s**-1': 3600.,
              'm/s': 3600000.}


def input_paths(inputs):
    """Preserve explicit file order; sort directory contents lexicographically."""
    inputs = [inputs] if isinstance(inputs, (str, Path)) else inputs
    paths = []
    for item in inputs:
        path = Path(item)
        if path.is_dir():
            paths.extend(sorted(p for p in path.rglob('*')
                                if p.is_file() and p.suffix.lower() in ('.nc', '.rain')))
        elif path.is_file():
            paths.append(path)
        else:
            raise ValueError(f'Input does not exist: {path}')
    if not paths or len({p.resolve() for p in paths}) != len(paths):
        raise ValueError('Supply at least one NetCDF file, without repeated paths')
    return paths


def classify(variable, override='auto', units_override=None):
    units = str(units_override or variable.attrs.get('units', '')).strip().lower()
    if units not in AMOUNT_UNITS and units not in RATE_UNITS:
        raise ValueError(f'Unsupported or missing precipitation units: {units!r}; use --units')
    declared = str(variable.attrs.get('rain_kind', '')).lower()
    description = ' '.join(str(variable.attrs.get(k, '')) for k in
                           ('description', 'long_name', 'standard_name')).lower()
    if override != 'auto':
        kind, evidence = override, 'explicit --rain-kind'
    elif declared in ('cumulative', 'interval', 'rate'):
        kind, evidence = declared, 'rain_kind attribute'
    elif units in RATE_UNITS:
        kind, evidence = 'rate', 'rate units'
    elif 'cumulat' in description or 'accumulated total' in description:
        kind, evidence = 'cumulative', 'cumulative metadata'
    elif str(variable.name).upper() in ('RAINNC', 'RAINC', 'RAINSH'):
        kind, evidence = 'cumulative', 'WRF cumulative variable name'
    elif re.search(r'time\s*:\s*sum', str(variable.attrs.get('cell_methods', '')).lower()):
        kind, evidence = 'interval', 'cell_methods time: sum'
    else:
        raise ValueError('Ambiguous precipitation amount; specify --rain-kind cumulative or interval')
    if (kind == 'rate') != (units in RATE_UNITS):
        raise ValueError(f'{kind} is inconsistent with units {units!r}')
    return kind, units, (RATE_UNITS if kind == 'rate' else AMOUNT_UNITS)[units], evidence


def axes(variable, dims=None):
    if variable.ndim != 3:
        raise ValueError('Multi-time precipitation must have exactly three dimensions')
    if dims:
        result = tuple(dims)
    else:
        aliases = ({'time', 't'}, {'y', 'south_north', 'lat', 'latitude'},
                   {'x', 'west_east', 'lon', 'longitude'})
        result = tuple(next((d for d in variable.dims if d.lower() in names), '')
                       for names in aliases)
    if len(set(result)) != 3 or set(result) != set(variable.dims):
        raise ValueError(f'Cannot identify (time, y, x) axes in {variable.dims}; use --dims TIME Y X')
    return result


def choose_variable(ds, name=None):
    if name:
        if name not in ds.data_vars:
            raise ValueError(f'Missing precipitation variable: {name}')
        return ds[name]
    wrf = [n for n in ds.data_vars if n.upper() == 'RAINNC']
    if len(wrf) == 1:
        return ds[wrf[0]]
    candidates = [n for n, v in ds.data_vars.items() if v.ndim == 3 and
                  (str(v.attrs.get('units', '')).strip().lower() in RATE_UNITS or
                   any(s in n.lower() for s in ('rain', 'precip')) or
                   'precipitation' in str(v.attrs.get('standard_name', '')).lower())]
    if len(candidates) != 1:
        raise ValueError(f'Cannot select one precipitation field from {candidates}; use --variable')
    return ds[candidates[0]]


class NetCDFRainSource:
    """Plan time/shape metadata once; load only requested rain slices per chunk."""

    def __init__(self, inputs, options):
        self.paths = input_paths(inputs)
        self.files = {}
        self.records = []
        self.initial_stats = {}
        self.kind = None
        self.calendar = None
        self.assumptions = []
        origin = getattr(options, 'time_origin', None)
        for path in self.paths:
            before = (path.stat().st_size, path.stat().st_mtime_ns)
            with xr.open_dataset(path) as ds:
                v = choose_variable(ds, getattr(options, 'variable', None))
                dims = axes(v, getattr(options, 'dims', None))
                kind, units, factor, evidence = classify(v, getattr(options, 'rain_kind', 'auto'),
                                                       getattr(options, 'units', None))
                attributes = {k: float(ds.attrs[k]) for k in GRID_KEYS if k in ds.attrs}
                spacing = getattr(options, 'grid_km', None)
                grid_spacing_source = 'source_attributes'
                if spacing is not None:
                    if not np.isfinite(spacing) or spacing <= 0:
                        raise ValueError('--grid-km must be finite and positive')
                    if not all(k in attributes for k in ('DX', 'DY')):
                        grid_spacing_source = 'explicit_grid_km'
                    for k in ('DX', 'DY'):
                        if k in attributes and attributes[k] != spacing * 1000:
                            raise ValueError('--grid-km disagrees with source grid attributes')
                        attributes[k] = spacing * 1000
                if not all(k in attributes for k in ('DX', 'DY')):
                    raise ValueError('Missing DX/DY; provide actual grid metadata or explicit --grid-km')
                if not all(np.isfinite(x) for x in attributes.values()):
                    raise ValueError('Nonfinite grid attributes')
                shape = [v.sizes[d] for d in dims[1:]]
                if min(shape) < 1:
                    raise ValueError('Empty spatial dimension')
                spatial_coordinates = {d: np.asarray(ds[d].values) for d in dims[1:]
                                       if d in ds.coords and ds[d].dims == (d,)}
                if self.kind is not None:
                    if kind != self.kind or attributes != self.attributes or shape != self.shape:
                        raise ValueError('Precipitation kind or grid changes across files')
                    if set(spatial_coordinates) != set(self.spatial_coordinates) or any(
                            not np.array_equal(value, self.spatial_coordinates[d], equal_nan=True)
                            for d, value in spatial_coordinates.items()):
                        raise ValueError('Spatial coordinates change across files')
                self.kind, self.attributes, self.shape = kind, attributes, shape
                self.spatial_coordinates = spatial_coordinates
                time_name = getattr(options, 'time_coordinate', None) or dims[0]
                if not getattr(options, 'time_coordinate', None) and time_name not in ds and 'Times' in ds:
                    time_name = 'Times'
                count = v.sizes[dims[0]]
                if count < 1:
                    raise ValueError(f'Empty time dimension: {path}')
                if time_name in ds:
                    coordinate = ds[time_name]
                    values = np.asarray(coordinate.values)
                    if coordinate.ndim == 2 and coordinate.dims[0] == dims[0] and values.dtype.kind in 'SU':
                        values = np.array([''.join(x.decode('ascii') if isinstance(x, bytes) else str(x)
                                                  for x in row).strip('\x00 ') for row in values])
                    elif coordinate.dims != (dims[0],):
                        raise ValueError('Time coordinate must be one-dimensional on the time axis, or WRF Times characters')
                    if values.dtype.kind in 'SU' or values.dtype.kind == 'O' and all(
                            isinstance(x, (str, bytes)) for x in values):
                        text = [x.decode('ascii') if isinstance(x, bytes) else str(x) for x in values]
                        try:
                            values = np.array([np.datetime64(x.replace('_', 'T')) for x in text])
                        except ValueError as error:
                            raise ValueError('Invalid string timestamps') from error
                elif getattr(options, 'time_coordinate', None):
                    raise ValueError(f'Missing time coordinate: {time_name}')
                else:
                    coordinate, values = None, np.arange(count)
                if np.issubdtype(values.dtype, np.datetime64):
                    times = values.astype('datetime64[s]')
                    calendar = True
                    if np.isnat(times).any() or not np.array_equal(times, values):
                        raise ValueError('Missing or subsecond calendar timestamps')
                elif np.issubdtype(values.dtype, np.number) or np.issubdtype(values.dtype, np.timedelta64):
                    durations = np.issubdtype(values.dtype, np.timedelta64)
                    numeric = (values / np.timedelta64(1, 'h')).astype(float) if durations else values.astype(float)
                    time_units = '' if durations else str(coordinate.attrs.get('units', '') if coordinate is not None else '').lower()
                    scale = {'': 1., 'hour': 1., 'hours': 1., 'h': 1.,
                             'seconds': 1/3600, 's': 1/3600, 'days': 24., 'd': 24.}.get(time_units)
                    if scale is None:
                        raise ValueError(f'Unsupported index time units: {time_units!r}')
                    if not np.isfinite(numeric).all() or not np.allclose(np.diff(numeric) * scale, 1., rtol=0, atol=1e-9):
                        raise ValueError('Index time values must be consecutive hourly observations')
                    calendar = False
                    # Local indices often restart at zero in each monthly file.
                    # Without calendar metadata, file order is the only available
                    # boundary evidence; record this assumption explicitly.
                    times = np.arange(len(self.records), len(self.records) + count)
                    self.assumptions.append(f'{path}: hourly index cadence; file order assumes continuous boundaries')
                    if self.records and numeric[0] != 0:
                        previous = self.files[self.records[-1][1]]['last_numeric_hour']
                        if not np.isclose(numeric[0] * scale - previous, 1., rtol=0, atol=1e-9):
                            raise ValueError('Numeric time gap between files; use continuous indices or zero-based local indices')
                else:
                    raise ValueError('Use CF-decoded Gregorian time or numeric hourly indices')
                if self.calendar is not None and calendar != self.calendar:
                    raise ValueError('Cannot mix calendar timestamps and index-only time files')
                self.calendar = calendar
                self.files[path] = dict(variable=str(v.name), dims=dims, kind=kind, units=units,
                    factor=factor, evidence=evidence, last_numeric_hour=None if calendar else numeric[-1] * scale)
                self.files[path].update(grid_spacing_source=grid_spacing_source, time_coordinate=time_name)
                self.records.extend((t, path, i) for i, t in enumerate(times))
            if (path.stat().st_size, path.stat().st_mtime_ns) != before:
                raise ValueError('Input changed during planning; possible incomplete transfer')
            self.initial_stats[path] = before
        if self.calendar:
            self.records.sort(key=lambda r: r[0])
        if origin:
            if self.calendar:
                raise ValueError('--time-origin applies only to index-only input')
            base = np.datetime64(origin, 's')
            if np.isnat(base):
                raise ValueError('Invalid --time-origin')
            self.records = [(base + np.timedelta64(int(t), 'h'), p, i) for t, p, i in self.records]
            self.assumptions.append('Calendar origin supplied explicitly by --time-origin')
        self.dated = self.calendar or bool(origin)
        self.entries = [(t, p) for t, p, _ in self.records]
        self.extra = int(self.kind == 'cumulative')
        self.available = len(self.records) - self.extra
        self.time_source = 'decoded_coordinate' if self.calendar else (
            'index_with_explicit_origin' if origin else 'index_assumed_hourly')

    def validate_plan(self, start, hours):
        selected = self.records[start:start + hours + self.extra]
        if start < 0 or hours < 1 or len(selected) != hours + self.extra:
            raise ValueError(f'Need {hours + self.extra} input samples; available intervals: {self.available}')
        delta = np.timedelta64(1, 'h') if self.dated else 1
        if any(b[0] - a[0] != delta for a, b in zip(selected, selected[1:])):
            raise ValueError('Planned samples must be consecutive hourly observations; duplicates/gaps are not allowed')

    def inspect(self):
        return dict(files=[str(p) for p in self.paths], samples=len(self.records),
            available_hours=self.available, rain_kind=self.kind, time_source=self.time_source,
            time_verified=bool(self.calendar), default_dt_hours=1., spatial_shape=self.shape,
            first_sample=str(self.records[0][0]), last_sample=str(self.records[-1][0]),
            ordering='calendar_time' if self.calendar else 'file_order',
            grid_attributes=self.attributes, assumptions=list(dict.fromkeys(self.assumptions)),
            file_fields={str(p): {k: list(v) if k == 'dims' else v
                                 for k, v in info.items() if k != 'last_numeric_hour'}
                         for p, info in self.files.items()})

    def read(self, a, entries=None):
        self.validate_plan(a.start, a.hours)
        selected = self.records[a.start:a.start + a.hours + self.extra]
        paths = list(dict.fromkeys(p for _, p, _ in selected))
        before = {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in paths}
        if any(before[p] != self.initial_stats[p] for p in paths):
            raise ValueError('Input changed since planning; possible incomplete transfer')
        blocks = {}
        for path in paths:
            info = self.files[path]
            indices = [i for _, p, i in selected if p == path]
            low, high = min(indices), max(indices) + 1
            with xr.open_dataset(path) as ds:
                v = ds[info['variable']].isel({info['dims'][0]: slice(low, high)}).transpose(*info['dims'])
                blocks[path] = (low, np.array(v.values, dtype=np.float32, copy=True) * info['factor'])
        samples = np.stack([blocks[p][1][i - blocks[p][0]] for _, p, i in selected])
        del blocks
        if np.isinf(samples).any() or np.any(samples[np.isfinite(samples)] < 0):
            raise ValueError('Negative or infinite precipitation values')
        if self.kind == 'cumulative':
            rain = np.diff(samples, axis=0)
            if np.any(rain[np.isfinite(rain)] < 0):
                raise ValueError('Cumulative precipitation decreased; possible restart/reset, including at a file boundary')
        else:
            rain = samples
        if any(not np.isfinite(frame).any() for frame in rain):
            raise ValueError('Missing whole rain interval')
        if any((p.stat().st_size, p.stat().st_mtime_ns) != before[p] for p in paths):
            raise ValueError('Input changed during read; possible incomplete transfer')
        ends = selected[self.extra:]
        metadata = dict(source_type='NetCDF hourly precipitation sequence',
            source_shape=[len(self.records), *self.shape], crop_origin_yx=[0, 0],
            dt_hours=1., normalized_units='mm/h', time_verified=bool(self.calendar),
            time_source=self.time_source, grid_source='attributes', grid_attributes=self.attributes,
            timestamps=[str(t) for t, _, _ in ends] if self.dated else None,
            interval_start_times=[str(t - np.timedelta64(1, 'h')) for t, _, _ in ends] if self.dated else None,
            sample_indices=list(range(a.start + self.extra, a.start + self.extra + a.hours)),
            input_files=[str(p) for p in paths], rain_kind=self.kind,
            coordinates_verified=False, assumptions=list(dict.fromkeys(self.assumptions)),
            input_detection=self.inspect()['file_fields'],
            precipitation_scope='Selected variable only; other precipitation components are not added')
        return np.ascontiguousarray(rain, dtype=np.float32), metadata
