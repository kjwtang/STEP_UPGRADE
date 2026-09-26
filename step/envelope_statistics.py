"""Per-frame core/envelope associations, not convective-storm classification."""
import numpy as np
from scipy import ndimage
from .envelope import area_weights

OBJECT_FIELDS = ['envelope_label', 'core_count', 'envelope_cells', 'core_cells',
                 'envelope_to_core_cell_ratio', 'peak_mm_h', 'core_peak_mm_h',
                 'mean_mm_h', 'rain_rate_sum', 'core_rain_rate_sum',
                 'noncore_rain_rate_fraction', 'touches_domain_boundary',
                 'touches_missing', 'envelope_area_km2', 'core_area_km2',
                 'rain_volume_rate_m3_h']


def envelope_statistics(rain, cores, envelopes, cell_area_km2=None):
    """Return object rows and many-to-many core membership rows.

    Core labels may be spatially disconnected under baseline morphology, so a
    core label can occur in several system envelopes. Preserve that relation.
    Rain volume rate is m3/hour, not accumulated volume: mm/hour * km2 * 1000.
    """
    rain, cores, envelopes = map(np.asarray, (rain, cores, envelopes))
    if rain.ndim != 2 or not rain.shape == cores.shape == envelopes.shape:
        raise ValueError('Rain and labels must have identical 2-D shapes')
    for labels in (cores, envelopes):
        if not np.issubdtype(labels.dtype, np.integer) or np.any(labels < 0):
            raise ValueError('Labels must be nonnegative integers')
    if np.any((envelopes > 0) & (~np.isfinite(rain) | (rain < 0))):
        raise ValueError('Assigned rain must be finite and nonnegative')
    if np.any((cores > 0) & (envelopes == 0)):
        raise ValueError('Every retained core pixel must belong to an envelope')
    area = area_weights(cell_area_km2, rain.shape)
    ids, inverse = np.unique(envelopes, return_inverse=True)
    inverse = inverse.reshape(rain.shape)
    n = len(ids)
    def totals(weights):
        return np.bincount(inverse.ravel(), weights=np.asarray(weights).ravel(), minlength=n)
    cells = totals(np.ones(rain.shape))
    core_cells = totals(cores > 0)
    rates = np.where(np.isfinite(rain), np.maximum(rain, 0), 0)
    sums, core_sums = totals(rates), totals(rates * (cores > 0))
    peaks = ndimage.maximum(rates, inverse, np.arange(n))
    core_peaks = ndimage.maximum(rates * (cores > 0), inverse, np.arange(n))
    boundary = np.zeros(rain.shape, dtype=bool)
    boundary[[0,-1],:] = True; boundary[:,[0,-1]] = True
    missing_neighbour = ndimage.binary_dilation(~np.isfinite(rain), structure=np.ones((3,3)))
    touches_boundary, touches_missing = totals(boundary), totals(missing_neighbour)
    env_area = totals(area) if area is not None else None
    core_area = totals(area * (cores > 0)) if area is not None else None
    volume = totals(rates * area) * 1000 if area is not None else None
    pairs, pair_counts = np.unique(np.column_stack((envelopes[cores > 0], cores[cores > 0])),
                                  axis=0, return_counts=True)
    memberships = [dict(envelope_label=int(e), core_label=int(c), core_cells=int(count))
                   for (e,c),count in zip(pairs,pair_counts)]
    seed_counts = dict(zip(*np.unique(pairs[:,0], return_counts=True))) if len(pairs) else {}
    rows = []
    for i, label in enumerate(ids):
        if label == 0:
            continue
        rows.append(dict(envelope_label=int(label), core_count=int(seed_counts.get(label,0)),
            envelope_cells=int(cells[i]), core_cells=int(core_cells[i]),
            envelope_to_core_cell_ratio=float(cells[i]/core_cells[i]) if core_cells[i] else None,
            peak_mm_h=float(peaks[i]), core_peak_mm_h=float(core_peaks[i]),
            mean_mm_h=float(sums[i]/cells[i]), rain_rate_sum=float(sums[i]),
            core_rain_rate_sum=float(core_sums[i]),
            noncore_rain_rate_fraction=float((sums[i]-core_sums[i])/sums[i]) if sums[i] else 0.,
            touches_domain_boundary=bool(touches_boundary[i]),
            touches_missing=bool(touches_missing[i]),
            envelope_area_km2=float(env_area[i]) if area is not None else None,
            core_area_km2=float(core_area[i]) if area is not None else None,
            rain_volume_rate_m3_h=float(volume[i]) if area is not None else None))
    return rows, memberships
