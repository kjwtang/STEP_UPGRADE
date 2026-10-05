"""Input-only regressions: file boundaries must not change scientific results."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from netcdf_rain_input import NetCDFRainSource
import run_reference_stream as stream
from test_reference_stream import dataset, options, rows, frame_hashes


def source_options(**overrides):
    values = dict(rain_kind='auto', variable=None, dims=None, units=None,
                  time_coordinate=None, time_origin=None, grid_km=None)
    values.update(overrides)
    return SimpleNamespace(**values)


def write(path, values, times=None, name='RAINNC', units='mm', attrs=None,
          dims=('time', 'y', 'x'), variable_attrs=None):
    grid = dict(DX=4000., DY=4000.) if attrs is None else attrs
    coords = {dims[0]: np.arange(len(values)) if times is None else times,
              dims[1]: np.arange(values.shape[1]), dims[2]: np.arange(values.shape[2])}
    xr.Dataset({name: (dims, values, dict(units=units, **(variable_attrs or {})))},
               coords=coords, attrs=grid).to_netcdf(path)
    return path


def read(source, start=0, hours=None):
    return source.read(SimpleNamespace(start=start, hours=source.available-start if hours is None else hours))


def test_single_numeric_month_defaults_hourly_and_keeps_first_predecessor(tmp_path):
    values = np.arange(7, dtype=np.float32)[:, None, None] * np.ones((1, 4, 5))
    path = write(tmp_path / 'month.nc', values)
    source = NetCDFRainSource(path, source_options())
    rain, meta = read(source)
    assert source.available == 6
    np.testing.assert_array_equal(rain, np.ones((6, 4, 5)))
    assert meta['timestamps'] is None and meta['time_verified'] is False
    assert meta['dt_hours'] == 1 and meta['sample_indices'] == list(range(1, 7))
    assert meta['time_source'] == 'index_assumed_hourly'


def test_multiple_numeric_months_cross_boundary_without_losing_first_hour(tmp_path):
    a = write(tmp_path / 'jan.nc', np.arange(3)[:, None, None].astype('f4'))
    b = write(tmp_path / 'feb.nc', np.arange(3, 6)[:, None, None].astype('f4'))
    source = NetCDFRainSource([a, b], source_options())
    assert source.available == 5
    np.testing.assert_array_equal(read(source)[0], np.ones((5, 1, 1)))
    np.testing.assert_array_equal(read(source, start=1, hours=3)[0], np.ones((3, 1, 1)))
    assert len(source.inspect()['assumptions']) == 2


@pytest.mark.parametrize('kind,units,name,attr,factor', [
    ('rate', 'mm/h', 'rain', {}, 1),
    ('rate', 'kg m-2 s-1', 'pr', {'standard_name': 'precipitation_flux'}, 3600),
    ('interval', 'mm', 'precip', {'cell_methods': 'time: sum'}, 1),
    ('interval', 'm', 'rain', {'rain_kind': 'interval'}, 1000),
])
def test_declared_rate_and_interval_are_not_differenced(tmp_path, kind, units, name, attr, factor):
    path = write(tmp_path / 'rain.nc', np.ones((3, 2, 2), dtype='f4'), name=name,
                 units=units, variable_attrs=attr)
    source = NetCDFRainSource(path, source_options())
    assert source.kind == kind and source.available == 3
    np.testing.assert_array_equal(read(source)[0], np.full((3, 2, 2), factor))


def test_ambiguous_amount_requires_explicit_kind(tmp_path):
    path = write(tmp_path / 'rain.nc', np.ones((3, 2, 2)), name='precip')
    with pytest.raises(ValueError, match='Ambiguous'):
        NetCDFRainSource(path, source_options())
    source = NetCDFRainSource(path, source_options(rain_kind='interval'))
    assert source.available == 3


def test_rainnc_can_be_explicitly_declared_prepared_interval(tmp_path):
    path = write(tmp_path / 'rain.nc', np.ones((3, 2, 2)))
    source = NetCDFRainSource(path, source_options(rain_kind='interval'))
    np.testing.assert_array_equal(read(source)[0], np.ones((3, 2, 2)))


def test_calendar_files_are_sorted_by_internal_time_not_names(tmp_path):
    times = np.datetime64('1996-01-31T22:00:00') + np.arange(6) * np.timedelta64(1, 'h')
    a = write(tmp_path / 'z_jan.nc', np.arange(3)[:, None, None].astype('f4'), times[:3])
    b = write(tmp_path / 'a_feb.nc', np.arange(3, 6)[:, None, None].astype('f4'), times[3:])
    source = NetCDFRainSource(tmp_path, source_options())
    rain, meta = read(source)
    np.testing.assert_array_equal(rain, np.ones((5, 1, 1)))
    assert meta['time_verified'] and meta['timestamps'][2] == '1996-02-01T01:00:00'


@pytest.mark.parametrize('hours', [[0, 1, 3], [0, 1, 1]])
def test_irregular_or_duplicate_calendar_times_rejected(tmp_path, hours):
    times = np.datetime64('1996-01-01') + np.array(hours) * np.timedelta64(1, 'h')
    path = write(tmp_path / 'rain.nc', np.arange(3)[:, None, None].astype('f4'), times)
    source = NetCDFRainSource(path, source_options())
    with pytest.raises(ValueError, match='consecutive hourly'):
        read(source)


@pytest.mark.parametrize('times', [[0, 1, 3], [0, 0, 1], [0, 2, 4]])
def test_irregular_numeric_time_is_not_silently_hourly(tmp_path, times):
    path = write(tmp_path / 'rain.nc', np.zeros((3, 2, 2)), times)
    with pytest.raises(ValueError, match='consecutive hourly'):
        NetCDFRainSource(path, source_options())


def test_file_boundary_reset_rejected_not_clipped(tmp_path):
    a = write(tmp_path / 'jan.nc', np.array([0, 10], dtype='f4')[:, None, None])
    b = write(tmp_path / 'feb.nc', np.array([0, 2], dtype='f4')[:, None, None])
    with pytest.raises(ValueError, match='decreased'):
        read(NetCDFRainSource([a, b], source_options()))


def test_bounded_read_does_not_inspect_unrequested_rain_values(tmp_path):
    path = write(tmp_path / 'rain.nc', np.array([0, 1, 2, -10], dtype='f4')[:, None, None])
    source = NetCDFRainSource(path, source_options())
    np.testing.assert_array_equal(read(source, hours=2)[0], np.ones((2, 1, 1)))
    with pytest.raises(ValueError, match='Negative'):
        read(source)


def test_axes_can_be_reordered(tmp_path):
    path = write(tmp_path / 'rain.nc', np.arange(4)[:, None, None] * np.ones((1, 3, 2)))
    with xr.open_dataset(path) as ds:
        changed = ds.load().transpose('x', 'time', 'y')
    other = tmp_path / 'transposed.nc'
    changed.to_netcdf(other)
    source = NetCDFRainSource(other, source_options())
    assert read(source)[0].shape == (3, 3, 2)


def test_grid_override_is_explicit_and_cannot_contradict_attributes(tmp_path):
    path = write(tmp_path / 'rain.nc', np.zeros((3, 2, 2)), attrs={})
    with pytest.raises(ValueError, match='Missing DX/DY'):
        NetCDFRainSource(path, source_options())
    assert NetCDFRainSource(path, source_options(grid_km=4)).attributes['DX'] == 4000
    path = write(tmp_path / 'known.nc', np.zeros((3, 2, 2)))
    with pytest.raises(ValueError, match='disagrees'):
        NetCDFRainSource(path, source_options(grid_km=1))


def test_explicit_origin_records_assumption_not_verified_calendar(tmp_path):
    path = write(tmp_path / 'rain.nc', np.arange(3)[:, None, None].astype('f4'))
    source = NetCDFRainSource(path, source_options(time_origin='1996-01-31T23:00:00'))
    _, meta = read(source)
    assert meta['timestamps'][0] == '1996-02-01T00:00:00'
    assert not meta['time_verified']


def test_snapshot_vs_single_file_vs_multifile_resume_exact_graph_and_ids(tmp_path):
    directory = dataset(tmp_path / 'snapshots')
    arrays = []
    for path in sorted(directory.glob('*.rain')):
        with xr.open_dataset(path) as ds:
            arrays.append(ds.RAINNC.values[0])
            attrs = dict(ds.attrs)
    values = np.stack(arrays)
    times = np.datetime64('1996-01-01') + np.arange(7) * np.timedelta64(1, 'h')
    month = write(tmp_path / 'month.nc', values, times, attrs=attrs)
    first = write(tmp_path / 'part1.nc', values[:3], times[:3], attrs=attrs)
    second = write(tmp_path / 'part2.nc', values[3:], times[3:], attrs=attrs)
    roots = [tmp_path / s for s in ('snap_result', 'single_result', 'multi_result')]
    baseline = stream.run(options(directory, roots[0], chunk_frames=6))
    single = stream.run(options(month, roots[1], hours=None, chunk_frames=4))
    paused = stream.run(options([second, first], roots[2], stop_after_chunks=1))
    assert paused['committed_frames'] == 2
    multiple = stream.run(options([second, first], roots[2], resume=True, chunk_frames=3))
    assert baseline['frames'] == single['frames'] == multiple['frames'] == 6
    for name in ('edges.csv', 'events.csv', 'objects.csv'):
        assert rows(roots[0], name) == rows(roots[1], name) == rows(roots[2], name)
    for key in ('identified_frame_sha256', 'tracked_frame_sha256'):
        assert frame_hashes(roots[0], key) == frame_hashes(roots[1], key) == frame_hashes(roots[2], key)


@pytest.mark.parametrize('kind', ['cumulative', 'interval'])
def test_numeric_sequence_resume_and_input_change_guard(tmp_path, kind):
    path = write(tmp_path / 'rain.nc', np.arange(7)[:, None, None] * np.ones((1, 4, 5)))
    out = tmp_path / 'run'
    a = options(path, out, hours=None, rain_kind=kind, stop_after_chunks=1)
    assert stream.run(a)['status'] == 'paused'
    assert stream.run(options(path, out, hours=None, rain_kind=kind, resume=True))['frames'] == (6 if kind == 'cumulative' else 7)
    path.touch()
    with pytest.raises(ValueError, match='Resume contract differs'):
        stream.run(options(path, out, hours=None, rain_kind=kind, resume=True))


def test_inspect_needs_no_output_or_sequence(tmp_path):
    path = write(tmp_path / 'rain.nc', np.zeros((741, 1, 1)))
    result = stream.run(options(path, None, hours=None, sequence_id=None, inspect=True))
    assert result['samples'] == 741 and result['available_hours'] == 740
    assert not list(tmp_path.glob('chunk_*'))


def test_numeric_time_units_seconds_and_decoded_durations(tmp_path):
    path = write(tmp_path / 'rain.nc', np.arange(3)[:, None, None].astype('f4'),
                 np.arange(3) * np.timedelta64(1, 'h'))
    source = NetCDFRainSource(path, source_options())
    np.testing.assert_array_equal(read(source)[0], np.ones((2, 1, 1)))


def test_numeric_cross_file_gap_rejected(tmp_path):
    a = write(tmp_path / 'jan.nc', np.zeros((3, 1, 1)), times=[0, 1, 2])
    b = write(tmp_path / 'feb.nc', np.ones((3, 1, 1)), times=[4, 5, 6])
    with pytest.raises(ValueError, match='Numeric time gap'):
        NetCDFRainSource([a, b], source_options())


def test_duplicate_boundary_calendar_time_rejected(tmp_path):
    times = np.datetime64('1996-01-31T23:00:00') + np.arange(5) * np.timedelta64(1, 'h')
    a = write(tmp_path / 'jan.nc', np.arange(3)[:, None, None].astype('f4'), times[:3])
    b = write(tmp_path / 'feb.nc', np.arange(2, 5)[:, None, None].astype('f4'), times[2:])
    with pytest.raises(ValueError, match='duplicates/gaps'):
        read(NetCDFRainSource([a, b], source_options()))


def test_grid_or_coordinates_changed_across_files_rejected(tmp_path):
    a = write(tmp_path / 'jan.nc', np.zeros((3, 2, 2)))
    b = write(tmp_path / 'feb.nc', np.zeros((3, 2, 2)), attrs=dict(DX=4000., DY=5000.))
    with pytest.raises(ValueError, match='grid changes'):
        NetCDFRainSource([a, b], source_options())
    b = write(tmp_path / 'feb2.nc', np.zeros((3, 2, 2)))
    with xr.open_dataset(b) as ds:
        changed = ds.load().assign_coords(x=[2, 3])
    other = tmp_path / 'shifted.nc'
    changed.to_netcdf(other)
    with pytest.raises(ValueError, match='Spatial coordinates'):
        NetCDFRainSource([a, other], source_options())


def test_units_kind_conflict_and_incomplete_transfer_rejected(tmp_path):
    path = write(tmp_path / 'rain.nc', np.arange(3)[:, None, None].astype('f4'), units='mm/h')
    with pytest.raises(ValueError, match='inconsistent'):
        NetCDFRainSource(path, source_options(rain_kind='cumulative'))
    source = NetCDFRainSource(path, source_options())
    path.touch()
    with pytest.raises(ValueError, match='changed since planning'):
        read(source)


def test_nan_barriers_and_whole_missing_interval(tmp_path):
    values = np.arange(3)[:, None, None] * np.ones((1, 2, 2))
    values[1, 0, 0] = np.nan
    path = write(tmp_path / 'rain.nc', values)
    rain, _ = read(NetCDFRainSource(path, source_options()))
    assert np.isnan(rain[:, 0, 0]).all() and np.isfinite(rain[:, 1, 1]).all()
    values[1] = np.nan
    path = write(tmp_path / 'missing.nc', values)
    with pytest.raises(ValueError, match='Missing whole rain interval'):
        read(NetCDFRainSource(path, source_options()))


def test_wrf_multitime_character_times_are_not_ignored(tmp_path):
    path = tmp_path / 'wrf_month.nc'
    times = ['1996-01-31_23:00:00', '1996-02-01_00:00:00', '1996-02-01_01:00:00']
    xr.Dataset({'RAINNC': (('Time', 'south_north', 'west_east'),
        np.arange(3, dtype='f4')[:, None, None], {'units': 'mm'}),
        'Times': (('Time', 'DateStrLen'), np.array([list(t) for t in times], dtype='S1'))},
        attrs=dict(DX=4000., DY=4000.)).to_netcdf(path)
    source = NetCDFRainSource(path, source_options())
    rain, meta = read(source)
    assert meta['time_verified'] and meta['timestamps'][0] == '1996-02-01T00:00:00'
    assert source.inspect()['file_fields'][str(path)]['time_coordinate'] == 'Times'
    np.testing.assert_array_equal(rain, np.ones((2, 1, 1)))
