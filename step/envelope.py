"""Experimental seeded rain envelopes; production identification is unchanged."""
import numpy as np
from scipy import ndimage

from .identification import identify_frame


def area_weights(cell_area_km2, shape):
    """Validate explicit cell areas; never infer physical area from array indices."""
    if cell_area_km2 is None:
        return None
    area = np.asarray(cell_area_km2, dtype=float)
    if area.ndim and area.shape != shape:
        raise ValueError('cell_area_km2 must be a scalar or match the frame')
    if not np.all(np.isfinite(area) & (area > 0)):
        raise ValueError('cell areas must be finite and positive')
    return np.broadcast_to(area, shape)


def identify_envelope(data, core_threshold=1.0, envelope_threshold=.1,
                      mode='partition', min_core_cells=1, connectivity=8,
                      morph_structure=None, min_core_area_km2=None,
                      cell_area_km2=None, max_distance_km=None, spacing_km=None):
    """Return (core labels, envelope labels).

    partition: rain-topography marker watershed; cores retain independent IDs.
    system: whole connected weak-rain components containing a retained core.
    Neither mode crosses dry/missing cells. Thresholds are rates in mm/hour.
    Core is an operational seed, not a diagnosis of deep convection.
    Optional distance cap is Euclidean distance to ANY retained core pixel,
    before segmentation. It is not an own-label or wet-path distance bound.
    spacing_km is (dy, dx) on a uniform Cartesian grid; not lon/lat degrees.
    """
    data = np.asarray(data)
    if data.ndim != 2:
        raise ValueError('Expected a 2-D rain-rate frame')
    if not (np.isfinite(core_threshold) and np.isfinite(envelope_threshold)
            and 0 < envelope_threshold <= core_threshold):
        raise ValueError('Require 0 < envelope_threshold <= core_threshold')
    if mode not in ('partition', 'system') or connectivity not in (4, 8):
        raise ValueError('Invalid mode or connectivity')
    if int(min_core_cells) != min_core_cells or min_core_cells < 1:
        raise ValueError('min_core_cells must be a positive integer')
    area = area_weights(cell_area_km2, data.shape)
    if min_core_area_km2 is not None:
        if area is None or not np.isfinite(min_core_area_km2) or min_core_area_km2 <= 0:
            raise ValueError('Positive min_core_area_km2 requires explicit cell areas')
    if max_distance_km is not None:
        spacing = np.asarray(spacing_km, dtype=float)
        if (not np.isfinite(max_distance_km) or max_distance_km < 0 or
                spacing.shape != (2,) or not np.all(np.isfinite(spacing) & (spacing > 0))):
            raise ValueError('Nonnegative distance cap requires positive (dy, dx) spacing_km')
    valid = np.isfinite(data)
    footprint = ndimage.generate_binary_structure(2, 1 if connectivity == 4 else 2)
    if morph_structure is None:
        cores, _ = ndimage.label(valid & (data >= core_threshold), footprint)
        counts = np.bincount(cores.ravel())
        keep = counts >= min_core_cells
        keep[0] = False
        cores[~keep[cores]] = 0
        ids = np.unique(cores)
        lookup = np.zeros(int(cores.max())+1, dtype=np.int32)
        retained = ids[ids > 0]
        lookup[retained] = np.arange(1,len(retained)+1)
        cores = lookup[cores]
    else:
        # Optional compatibility seeds use the existing STEP joining rule.
        cores = identify_frame(data, morph_structure=morph_structure,
                               min_size=min_core_cells, threshold=core_threshold)
    if min_core_area_km2 is not None:
        sizes = np.bincount(cores.ravel(), weights=area.ravel())
        keep = sizes >= min_core_area_km2
        keep[0] = False
        cores[~keep[cores]] = 0
        ids = np.unique(cores[cores > 0])
        lookup = np.zeros(int(cores.max()) + 1, dtype=np.int32)
        lookup[ids] = np.arange(1, len(ids) + 1)
        cores = lookup[cores]
    wet = valid & (data >= envelope_threshold)
    if not cores.any():
        return cores, np.zeros(data.shape,dtype=np.int32)
    if max_distance_km is not None:
        wet &= ndimage.distance_transform_edt(cores == 0, sampling=spacing) <= max_distance_km
    if mode == 'system':
        components, _ = ndimage.label(wet, footprint)
        seeded = np.unique(components[cores > 0])
        lookup = np.zeros(int(components.max())+1,dtype=np.int32)
        lookup[seeded] = np.arange(1,len(seeded)+1)
        return cores, lookup[components]
    try:
        from skimage.segmentation import watershed
    except ImportError as exc:
        raise ImportError('Partition experiments require pip install -e ".[envelope]"') from exc
    surface = -np.where(valid, data, 0).astype(np.float64)
    envelopes = watershed(surface, markers=cores, mask=wet, connectivity=footprint)
    return cores, envelopes.astype(np.int32)
