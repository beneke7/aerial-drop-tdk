"""Band interpolation of a contour map: the central reading of Chapter 2 (§2.2).

A published coverage map gives only its contour levels C_1 < C_2 < ... A stepped field (each cell at
the highest level whose contour encloses it) is the lower bound. Inside the band between contour k
(level C_k) and the next inner contour (level C_{k+1}) this module interpolates by distance:

    C = C_k + (C_{k+1} - C_k) d_out / (d_out + d_in),

with d_out the distance to the band's outer contour and d_in the distance to the next inner contour.
The value stays inside the band the map gives, so thresholding the result at any contour level gives
back that contour exactly. Inside a contour that has no higher contour in it (the top contour, or a
local top) the level itself is kept, because the peak is not observed; below the lowest contour the
field stays zero (unobserved).

Distances are between cell centres on the raster, less half a cell, so they measure to the contour
that runs between two cells. `d_in` is measured to the inner contour of the same connected region.
"""

import numpy as np
from scipy import ndimage


def band_interp(C, cell):
    """Band-interpolated field from a stepped field C[ix, iy] (any unit, zero = below the lowest
    contour) on square cells of side `cell` (m). The levels are the positive values present in C."""
    C = np.pad(np.asarray(C, float), 1)  # a ring of zeros, so every region has an outside
    out = C.copy()
    levels = np.unique(C[C > 0])
    for lo, hi in zip(levels[:-1], levels[1:]):
        sup, inner = C >= lo, C >= hi
        if not (sup & ~inner).any():
            continue
        d_out = ndimage.distance_transform_edt(sup) * cell - cell / 2
        lab, _ = ndimage.label(sup)
        for i, sl in enumerate(ndimage.find_objects(lab), start=1):
            own = lab[sl] == i
            core = own & inner[sl]
            band = own & ~inner[sl]
            if not core.any() or not band.any():
                continue  # a local top: its peak is unknown, keep the level
            d_in = ndimage.distance_transform_edt(~core) * cell - cell / 2
            f = d_out[sl][band] / (d_out[sl][band] + d_in[band])
            o = out[sl]
            o[band] = lo + (hi - lo) * f
    return out[1:-1, 1:-1]
