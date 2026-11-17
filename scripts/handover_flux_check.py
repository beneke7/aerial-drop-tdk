#!/usr/bin/env python3
"""Exact liquid flux through the handover surface S against the band crossing estimator.

    .venv/bin/python scripts/handover_flux_check.py <bands dir> [--xh 5 --yh -10]

<bands dir> holds postProcessing/{liquidBand,SFlux,nozzleFlux,*OutletFlux,waterVolume} of a
Bambi v2c run (CASE_RECORD.md section 13.1). The face zone SZone is B's whole boundary, so the
liquid leaving through S per step is sum_SZone(alphaPhi_) - sum_nozzle(alphaPhi_). For every band
write interval the script compares that volume with the volume the crossing test hands over
(alpha >= 0.05, and all crossing liquid).
"""

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import pyvista as pv


def series(base: Path, name: str, col: int):
    rows = np.vstack([np.loadtxt(f, comments="#", ndmin=2) for f in sorted(glob.glob(str(base / name / "*" / "*.dat")))])
    rows = rows[np.argsort(rows[:, 0], kind="stable")]
    _, i = np.unique(rows[:, 0], return_index=True)  # restarts repeat a few steps
    return rows[i, 0], rows[i, col]


def cumulative(t, q):
    return np.r_[0.0, np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bands", type=Path)
    ap.add_argument("--xh", type=float, default=5.0)
    ap.add_argument("--yh", type=float, default=-10.0)
    ap.add_argument("--amin", type=float, default=0.05)
    a = ap.parse_args()
    pp = a.bands / "postProcessing"
    tz, qz = series(pp, "SFlux", 2)
    tn, qn = series(pp, "nozzleFlux", 2)
    qs = qz - np.interp(tz, tn, qn)  # out through S, m3/s
    S = cumulative(tz, qs)
    inj = cumulative(tn, -qn)
    files = json.loads(next((pp / "liquidBand").glob("*.series")).read_text())["files"]
    rows, prev = [], 0.0
    for f in files:
        t = f["time"]
        g = pv.read(pp / "liquidBand" / f["name"]).combine().compute_cell_sizes(length=False, area=False)
        al = np.asarray(g.cell_data["alpha.water"])
        v = al * np.asarray(g.cell_data["Volume"])
        x, u = np.asarray(g.cell_centers().points), np.asarray(g.cell_data["U"])
        dt = t - prev

        def inside(p):
            return (p[:, 0] < a.xh) & (p[:, 1] > a.yh)

        cr = inside(x - u * dt) & ~inside(x)
        exact = float(np.interp(t, tz, S) - np.interp(prev, tz, S))
        rows.append({"t": t, "dt": dt, "exact_m3": exact, "est_all_m3": float(v[cr].sum()),
                     "est_amin_m3": float(v[cr & (al >= a.amin)].sum()), "band_liquid_m3": float(v.sum())})
        prev = t
    ex, ea, em = (np.array([r[k] for r in rows]) for k in ("exact_m3", "est_all_m3", "est_amin_m3"))
    out = {"injected_m3": float(inj[-1]), "exact_through_S_m3": float(S[-1]),
           "estimator_all_alpha_m3": float(ea.sum()), "estimator_alpha_ge_amin_m3": float(em.sum()),
           "ratio_all_to_exact": float(ea.sum() / ex.sum()), "ratio_amin_to_exact": float(em.sum() / ex.sum()),
           "per_write": rows}
    (a.bands / "flux_check.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "per_write"}, indent=1))


if __name__ == "__main__":
    main()
