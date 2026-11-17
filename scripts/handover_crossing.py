#!/usr/bin/env python3
"""Continuous VOF -> drops handover at a fixed surface, from band files written during a run.

Bambi v2 (cases/bambi_bucket/CASE_RECORD.md section 13). Aircraft frame: nozzle at the origin, y up,
air arriving along +x at U. Region B = {x < xh, y > yh}. Each band file (vtkWrite, time t_w) holds
the cells just outside B. A liquid cell (alpha >= amin) is handed over at t_w if, by its own
velocity, it was inside B one write interval earlier: inB(x - u dt) and not inB(x). Connected
crossing cells form one parcel with their exact liquid volume, centroid and momentum.

Ground frame: air at rest, the aircraft moving along -x at U, the nozzle at the origin at t = 0.
A parcel starts at (x - U t_w, y, z) with velocity u - (U, 0, 0) and flies to the ground H below the
nozzle with handover_drops.fly (drag, gravity, secondary breakup; the air is steady, so each
parcel's own start time does not matter). Output: the ground deposit (kg per cell) on 0.5 m and
0.25 m cells, the mass ledger, the drop spectrum and per-write counts.

    handover_crossing.py <band.vtu.series> --U 21.607 --H 42.062 [--xh 5 --yh -10 --wind 0
                         --fragments hf --tag NAME --out DIR]

Deterministic; input hashes and arguments go into the JSON.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import handover_drops as hd
import numpy as np
import pyvista as pv


def crossing_parcels(path: Path, dt: float, xh: float, yh: float, amin: float, per_cell: bool = False):
    """Volume, centroid and velocity of the connected liquid (alpha >= amin) that left B in the last
    dt, and the volume of crossing liquid below amin (not handed over)."""
    g = pv.read(path)
    if isinstance(g, pv.MultiBlock):  # vtkWrite writes .vtm with an internal.vtu block
        g = g.combine()
    if g.n_cells == 0:
        return np.zeros(0), np.zeros((0, 3)), np.zeros((0, 3)), 0.0
    g = g.clean(tolerance=1e-6).compute_cell_sizes(length=False, area=False)
    x, u = np.asarray(g.cell_centers().points), np.asarray(g.cell_data["U"], float)
    a = np.asarray(g.cell_data["alpha.water"], float)

    def inside(p):
        return (p[:, 0] < xh) & (p[:, 1] > yh)

    crossed = inside(x - u * dt) & ~inside(x)
    vcell = np.asarray(g.cell_data["Volume"], float) * a
    diffuse = float(vcell[crossed & (a < amin)].sum())  # crossing liquid not handed over (ledger)
    sel = np.where((a >= amin) & crossed)[0]
    if not len(sel):
        return np.zeros(0), np.zeros((0, 3)), np.zeros((0, 3)), diffuse
    sub = g.extract_cells(sel)
    if per_cell:  # each crossing cell its own parcel: the fine limit (VOF resolves nothing smaller)
        lab = np.arange(sub.n_cells)
    else:  # connected crossing cells form one sphere: the coarse limit
        lab = np.asarray(sub.connectivity(extraction_mode="all").cell_data["RegionId"])
    v = np.asarray(sub.cell_data["Volume"], float) * np.asarray(sub.cell_data["alpha.water"], float)
    xs, us = np.asarray(sub.cell_centers().points), np.asarray(sub.cell_data["U"], float)
    n = lab.max() + 1
    vol = np.bincount(lab, v, n)
    cen = np.stack([np.bincount(lab, v * xs[:, i], n) for i in range(3)], 1) / vol[:, None]
    vel = np.stack([np.bincount(lab, v * us[:, i], n) for i in range(3)], 1) / vol[:, None]
    return vol, cen, vel, diffuse


def ground_map(land: np.ndarray, mass: np.ndarray, cell: float) -> dict:
    """Mass (kg) per square cell of side `cell` on a support covering every landing point."""
    lo = np.floor(land.min(0) / cell) * cell - cell
    hi = np.ceil(land.max(0) / cell) * cell + cell
    ex, ez = (np.arange(lo[i], hi[i] + cell / 2, cell) for i in (0, 2))
    m, _, _ = np.histogram2d(land[:, 0], land[:, 2], bins=[ex, ez], weights=mass)
    return {"x_edges": ex, "z_edges": ez, "mass_kg": m}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("series", type=Path, help="vtkWrite .vtm.series (or .vtu.series) of the band writer")
    ap.add_argument("--U", type=float, required=True, help="flight (crossflow) speed, m/s")
    ap.add_argument("--H", type=float, required=True, help="nozzle height above ground, m")
    ap.add_argument("--xh", type=float, default=5.0)
    ap.add_argument("--yh", type=float, default=-10.0)
    ap.add_argument("--amin", type=float, default=0.05)
    ap.add_argument("--wind", type=float, default=0.0, help="wind speed, m/s (crosswind along z unless --wind-dir)")
    ap.add_argument("--wind-dir", type=float, default=90.0,
                    help="direction the wind blows TO, degrees from +x (the air-frame downstream direction): "
                         "0 = headwind (the aircraft flies -x), 90 = crosswind +z, 180 = tailwind")
    ap.add_argument("--fragments", choices=("stable", "hf"), default="hf")
    ap.add_argument("--flux-check", type=Path, default=None,
                    help="flux_check.json (scripts/handover_flux_check.py): scale each write's parcel volumes so the "
                         "handed-over volume equals the exact liquid flux through S in that interval")
    ap.add_argument("--per-cell", action="store_true",
                    help="hand over every crossing cell as its own parcel (fine limit) instead of connected regions")
    ap.add_argument("--fill", type=float, default=1.0,
                    help="fill fraction f: keep the liquid crossing S after the first (1-f) V of the exact S flux "
                         "(prismatic Torricelli: a partial fill is the tail of the full drain; needs --flux-check)")
    ap.add_argument("--zref", type=float, default=None,
                    help="log wind profile through --wind at this height (m); default: uniform wind")
    ap.add_argument("--z0", type=float, default=0.03, help="roughness length (m) of the log profile")
    ap.add_argument("--turb", action="store_true", help="surface-layer turbulence (needs --zref)")
    ap.add_argument("--obukhov", type=float, default=None,
                    help="Obukhov length L (m; <0 unstable, >0 stable); default neutral (needs --zref)")
    ap.add_argument("--meander-deg", type=float, default=0.0, help="sigma of the slow wind-direction walk, deg (ASSUMED)")
    ap.add_argument("--meander-tau", type=float, default=60.0, help="time scale of the direction walk, s (ASSUMED)")
    ap.add_argument("--hold-track", action="store_true",
                    help="the aircraft holds its GROUND track: the release points stay on the line z = 0 instead of "
                         "drifting with the wind (default: fixed heading in the air, z += w_z t). Only the cross-track "
                         "release position changes (z of every parcel -= w_z t_w, + z toward the wind); velocities, "
                         "timing and the crab angle (body axis vs track) are not modified")
    ap.add_argument("--hf-classes", type=int, default=None,
                    help="number of HF child size classes (default 5 fixed root-normal classes)")
    ap.add_argument("--hf-jitter", action="store_true",
                    help="random quantile offset per split within each size class (stratified continuous sampling)")
    ap.add_argument("--hf-smd-scale", type=float, default=1.0,
                    help="multiply the Hsiang-Faeth child SMD by this factor (sensitivity study; default 1 = unchanged)")
    ap.add_argument("--scale-count", action="store_true",
                    help="apply the --flux-check ledger scaling to the number of drops per parcel, not to their diameter "
                         "(the missed water is taken as more drops of the same size; 17 Nov handover diagnostic)")
    ap.add_argument("--range-check", action="store_true",
                    help="add a breakup_range_check block (mass outside the HF range, as handover_drops.py) to the JSON")
    ap.add_argument("--cache", type=Path, default=None,
                    help="npz of the crossing parcels before scaling: read if present, else written")
    ap.add_argument("--snapshots", action="store_true",
                    help="also write snapshots_<tag>.npz: parcel positions every 0.1 s of flight (visualization)")
    ap.add_argument("--air-field", type=Path, default=None,
                    help="gridded time-mean VOF air field (scripts/air_field_grid.py): parcels fly through it "
                         "(rotor wake, air dragged by the liquid), added to the wind; zero outside the VOF domain")
    ap.add_argument("--rh", type=float, default=None, help="relative humidity (fraction); switches Denner evaporation on")
    ap.add_argument("--tair", type=float, default=None, help="air temperature, K (with --rh)")
    ap.add_argument("--tdrop", type=float, default=None, help="initial drop temperature, K (default: --tair)")
    ap.add_argument("--drag", choices=hd.DRAG_LAWS, default="hf",
                    help="drag law (handover_drops.drag_coefficient): hf = shock-tube deformation (default, all maps before "
                         "17 Nov 2026); tab = TAB equilibrium distortion, free-fall terminal velocities within 6 %% of Gunn-Kinzer")
    ap.add_argument("--turb-seed", type=int, default=1, help="seed of the turbulence stream (default 1 = existing maps)")
    ap.add_argument("--tag", default="crossing")
    ap.add_argument("--out", type=Path, default=Path("results/runs/handover"))
    a = ap.parse_args()
    ev = None
    if a.rh is not None:
        if a.tair is None:
            ap.error("--rh needs --tair")
        ev = hd.load_evaporation(a.rh, a.tair, a.tdrop if a.tdrop is not None else a.tair)

    if a.hf_smd_scale != 1.0:
        _smd = hd.hf_smd
        hd.hf_smd = lambda d, u: a.hf_smd_scale * _smd(d, u)
    files = json.loads(a.series.read_text())["files"]
    times = np.array([f["time"] for f in files], float)
    gaps = np.diff(times)
    if a.cache and a.cache.exists():
        z = np.load(a.cache)
        vol, cen, vel, tw = (z[k] for k in ("vol", "cen", "vel", "tw"))
        per_write, hashes = json.loads(str(z["per_write"])), json.loads(str(z["hashes"]))
    else:
        vol, cen, vel, tw, per_write, hashes = [], [], [], [], [], {}
        for k, f in enumerate(files):
            path = a.series.parent / f["name"]
            dt = gaps[k - 1] if k else (np.median(gaps) if len(gaps) else times[0])
            v, c, u, diffuse = crossing_parcels(path, float(dt), a.xh, a.yh, a.amin, a.per_cell)
            vol.append(v), cen.append(c), vel.append(u), tw.append(np.full(len(v), times[k]))
            per_write.append({"t": float(times[k]), "parcels": int(len(v)), "volume_m3": float(v.sum()),
                              "diffuse_below_amin_m3": diffuse})
            hashes[f["name"]] = hashlib.sha256(path.read_bytes()).hexdigest()
        vol, cen, vel, tw = np.concatenate(vol), np.vstack(cen), np.vstack(vel), np.concatenate(tw)
        if a.cache:
            np.savez(a.cache, vol=vol, cen=cen, vel=vel, tw=tw, per_write=json.dumps(per_write), hashes=json.dumps(hashes))
    vol_unscaled = vol.copy()
    if a.flux_check:  # per-interval ledger correction (declared assumption: the missed liquid is
        # distributed like the counted liquid of the same interval)
        exact = {round(r["t"], 6): r["exact_m3"] for r in json.loads(a.flux_check.read_text())["per_write"]}
        ex = np.array([exact.get(round(t, 6), 0.0) for t in times])
        drop = (1 - a.fill) * ex.sum()  # partial fill: the first (1 - f) V through S is not released
        keep = np.clip((np.cumsum(ex) - drop) / np.where(ex > 0, ex, 1), 0, 1)
        for k in range(len(times)):
            sel = tw == times[k]
            got = vol[sel].sum()
            scale = ex[k] * keep[k] / got if got > 0 else 0.0
            vol[sel] *= scale
            per_write[k]["ledger_scale"] = float(scale)
    else:
        assert a.fill == 1.0, "--fill needs --flux-check"
    ok = vol > 0
    vol, cen, vel, tw, vol_unscaled = vol[ok], cen[ok], vel[ok], tw[ok], vol_unscaled[ok]
    # air frame (the VOF) -> ground frame: the aircraft moves at w - (U, 0, 0), the air at w (at the bucket)
    th = np.radians(a.wind_dir)
    atm = hd.atmosphere(a.wind, a.wind_dir, a.zref, a.z0, a.turb, L=a.obukhov, meander_deg=a.meander_deg,
                         meander_tau=a.meander_tau, meander_seed=int(a.wind_dir)) if a.zref else None
    assert a.zref or not a.turb, "--turb needs --zref"
    w_bucket = a.wind * (np.log(a.H / a.z0) / np.log(a.zref / a.z0) if atm and a.obukhov is None else 1.0)
    if atm and a.obukhov is not None:
        w_bucket = float(hd.wind_speed(atm, a.H))
    w = w_bucket * np.array([np.cos(th), 0.0, np.sin(th)])
    cen[:, 0] += (w[0] - a.U) * tw
    vg_z = 0.0 if a.hold_track else w[2]  # ground cross-track speed of the aircraft
    cen[:, 2] += vg_z * tw
    vel += w - np.array([a.U, 0.0, 0.0])
    m = hd.RHO_L * vol
    info: dict = {}
    rec = {"dt": 0.1} if a.snapshots else None
    field = None
    if atm is not None:
        atm["t0"] = tw  # global release time of each input parcel (meander)
    if a.air_field:  # the aircraft frame moves at w - (U, 0, 0) over the ground, as for cen above
        z = np.load(a.air_field)
        field = {"F": {"origin": z["origin"], "h": float(z["h"]), "du": z["du"]}, "t0": tw,
                 "v_frame": np.array([w[0] - a.U, 0.0, vg_z])}
        hashes["air_field"] = hashlib.sha256(a.air_field.read_bytes()).hexdigest()
    land, tl, d, n, frac, parent = hd.fly((6 * (vol_unscaled if a.scale_count else vol) / np.pi) ** (1 / 3), cen.copy(), vel.copy(), a.H, True,
                                          w, fragments=a.fragments, evap=ev, info=info, atm=atm,
                                          record=rec, field=field,
                                          hf_sampling=hd.hf_class_sampler(a.hf_classes or 5, a.hf_jitter)
                                          if (a.hf_classes or a.hf_jitter) else None, drag=a.drag, turb_seed=a.turb_seed)
    mg = m[parent] * frac * info["w"]
    ok = np.isfinite(tl)
    a.out.mkdir(parents=True, exist_ok=True)
    maps = {f"{c:g}": ground_map(land[ok], mg[ok], c) for c in (0.5, 0.25)}
    if rec:  # flight-time snapshots; a parcel's global time is tw[parent] + flight time
        sn = rec["snaps"]
        np.savez_compressed(a.out / f"snapshots_{a.tag}.npz", t=np.array([s[0] for s in sn]), n=np.array([len(s[1]) for s in sn]),
                            parent=np.concatenate([s[1] for s in sn]), pos=np.concatenate([s[2] for s in sn]),
                            d=np.concatenate([s[3] for s in sn]), share=np.concatenate([s[4] for s in sn]),
                            tw=tw, m=m, land=land, tland=tl, mg=mg, land_parent=parent)
    np.savez_compressed(a.out / f"ground_map_{a.tag}.npz",
                        **{f"{k}_{c}": v for c, g in maps.items() for k, v in g.items()})
    res = {
        "args": {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(a).items()},
        "frame": "ground: x along track (aircraft moves -x at U), z cross-track, nozzle at the origin at t = 0",
        "handover_surface": {"xh_m": a.xh, "yh_m": a.yh, "region_B": "x < xh and y > yh (aircraft frame)"},
        "parcels": int(len(vol)), "handed_over_m3": float(vol.sum()),
        "crossing_below_amin_m3": float(sum(w["diffuse_below_amin_m3"] for w in per_write)), "per_write": per_write,
        "mass_ledger": hd.book(m, parent, frac, tl, info),
        "map_integral_check_kg": {c: float(g["mass_kg"].sum()) for c, g in maps.items()},
        "ground_spectrum": hd.quantiles(d[ok], mg[ok]) if ok.any() else {},
        "ground_mass_fraction_d_below_0p3mm": float(mg[ok][d[ok] < 3e-4].sum() / mg[ok].sum()) if ok.any() else 0.0,
        "input_sha256": hashes,
    }
    if a.range_check:
        evs = info["events"]
        me = evs["mass_frac"] * m[evs["parent"].astype(int)]
        fi, out = evs["first"].astype(bool), evs["out"].astype(bool)
        res["breakup_range_check"] = {
            "range": f"We < {hd.HF_WE_MAX:g}, Oh < {hd.HF_OH_MAX}",
            "breakup_events": int(len(me)),
            "mass_outside_range_fraction_all_events": float(me[out].sum() / me.sum()) if len(me) else None,
            "first_breakups": int(fi.sum()),
            "mass_outside_range_fraction_first_breakups": float(me[fi & out].sum() / me[fi].sum()) if fi.any() else None,
            "first_breakup_We_mass_median": hd.wmedian(evs["we"][fi], me[fi]) if fi.any() else None,
            "ground_mass_fraction_with_outside_range_ancestor": float(mg[ok & info["oor"]].sum() / mg[ok].sum()),
            "handed_over_mass_kg": float(m.sum()),
        }
    if ev:
        res["evaporation"] = {"status": "ASSUMED weather, not measured", "model": "handover_drops.load_evaporation (Denner 2026)",
                              "rh": a.rh, "T_air_K": a.tair, "T_drop_K": ev["tdrop"], "wet_bulb_K": ev["twb"]}
    (a.out / f"crossing_{a.tag}.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("parcels", "handed_over_m3", "mass_ledger", "map_integral_check_kg",
                                          "ground_spectrum")}, indent=1))


if __name__ == "__main__":
    main()
