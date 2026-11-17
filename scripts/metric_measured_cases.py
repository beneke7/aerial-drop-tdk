#!/usr/bin/env python3
"""Line yield on 20 digitized water-pattern panels (cup, infrared and canopy).

    .venv/bin/python scripts/metric_measured_cases.py   -> results/06_hold_metric/measured/

Each map is a set of superlevel polygons per coverage level (lower and upper envelopes; the
Canadian maps also a central set). The dose field is the step-down field: the highest level whose
polygon contains the cell, so it never exceeds the published contours. η is normalized by the
reported released volume (water, 1000 kg/m³) where reported, otherwise by the deposited water of
the field itself. The Black Hawk report gives no released volume, so those maps use an assumed
80 % load (user decision, 6 Nov 2026). Canopy maps are throughfall under a pine stand, not
open-field deposits; they are kept but flagged.
"""

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402

from aerial_drop.hold import line_yield  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CX16 = ROOT / "data" / "measured_maps"  # digitised measured maps (not included)
OUT = ROOT / "results/06_hold_metric/measured"
GPC = 1000 * 0.003785411784 / (100 * 0.3048**2)  # kg/m2, exact unit conversions
POINTS = {"grass": (1.0 * GPC, 3.75), "shrub": (2.4, 9.0)}  # (C_req kg/m2, W m), as in chapter 2
C_GRID = np.linspace(0.05, 4.1, 82)
W_GRID = np.linspace(0, 60, 121)
CELL = 0.5
BAMBI_LOAD_M3 = 0.8 * 660 * 0.003785411784  # [ASSUMED] 80 % of the 660 US gal capacity; the report states no load


def water_cases():
    """The digitised water maps, with the assumed Black Hawk load filled in and labelled."""
    cases = [c for c in json.load(open(CX16 / "data/cases.json"))["cases"] if c.get("liquid") == "water"]
    for c in cases:
        if c["case_id"].startswith("blackhawk_bambi660"):
            c["released_volume_m3"], c["basis"] = BAMBI_LOAD_M3, "assumed 80 % load"
    return cases


def basis_of(c):
    return c.get("basis") or ("released" if c.get("released_volume_m3") else "deposited (step-down)")


def variant_of(p):
    tag = (p.get("interpretation") or p.get("variant") or "").lower()
    return "lower" if "lower" in tag else "upper" if "upper" in tag else "central" if "central" in tag else None


def level_of(p):
    return p.get("coverage_gpc", p.get("level_gpc"))


def fields(geojson):
    """Step-down dose fields (kg/m2) per variant on one common grid; x along track."""
    sets = {}
    for f in geojson["features"]:
        p, g = f["properties"], f["geometry"]
        v, lv = variant_of(p), level_of(p)
        if v is None or lv is None or not g or not g.get("coordinates"):
            continue
        polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
        dose = p["depth_mm"] if p.get("depth_mm") is not None else lv * GPC
        sets.setdefault(v, []).append((dose, polys))
    pts = np.vstack([np.asarray(r)[:, :2] for rings in sets.values() for _, ps in rings for poly in ps for r in poly])
    lo, hi = pts.min(0) - 2, pts.max(0) + 2
    xs, ys = np.arange(lo[0], hi[0], CELL) + CELL / 2, np.arange(lo[1], hi[1], CELL) + CELL / 2
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    q = np.c_[X.ravel(), Y.ravel()]
    out = {}
    for v, rings in sets.items():
        C = np.zeros(len(q))
        for dose, polys in sorted(rings, key=lambda t: t[0]):
            inside = np.zeros(len(q), bool)
            for poly in polys:
                m = MPath(np.asarray(poly[0])[:, :2]).contains_points(q)
                for hole in poly[1:]:
                    m &= ~MPath(np.asarray(hole)[:, :2]).contains_points(q)
                inside |= m
            C[inside] = dose
        out[v] = C.reshape(X.shape)
    return out, xs, ys


# Real scenarios: USFS coverage level by fuel (Black Hawk report MTDC 0157-2804, Table 1, p. 5;
# retardant advice transferred to water here as an assumed scenario) x intensity -> Byram flame
# length L = 0.0775 I^0.46 (m, kW/m) -> required width 2.5 L (Restas 2023).
FUELS = {"grass (1 GPC)": 1, "timber litter, sagebrush (2)": 2, "green brush, light slash (3)": 3,
         "heavy litter (4)": 4, "cured brush, chaparral (6)": 6}
INTENSITIES = (100, 500, 1000, 2000, 4000)  # kW/m


def width_req(i_kw_m):
    return 2.5 * 0.0775 * i_kw_m**0.46


def tagged_spacing(C, c_req, w_req):
    """Largest along-track spacing (m) at which a train of identical drops gives a continuous held
    line: the map is folded modulo the spacing (periodic superposition of shifted copies), and every
    column of the folded field must hold. None if no spacing works (wider than one pass)."""
    nx = C.shape[0]
    for k in range(nx, 0, -1):  # spacing in cells, largest first
        f = np.pad(C, ((0, (-nx) % k), (0, 0))).reshape(-1, k, C.shape[1]).sum(axis=0)
        if line_yield(f, CELL, CELL, c_req, w_req, 1.0)["line_m"] >= k * CELL - 1e-9:
            return k * CELL
    return None


def scenario_table(cases):
    rows = []
    for c in cases:
        if "canopy" in c["case_id"] or c["case_id"].startswith("bell_s205"):
            continue  # throughfall, or unresolved contour labels
        fs, _, _ = fields(json.load(open(CX16 / c["geojson_path"])))
        v = "central" if "central" in fs else "lower"
        C = fs[v]
        litres = 1000 * c["released_volume_m3"] if c.get("released_volume_m3") else float(C.sum() * CELL * CELL)
        Cup = step_up(C)  # optimistic: dose read at the next contour above
        # water basis consistent with the field read: at least the water that field contains
        litres_up = max(1000 * (c.get("released_volume_m3") or 0), float(Cup.sum() * CELL * CELL))
        for fuel, gpc in FUELS.items():
            for i in INTENSITIES:
                s = tagged_spacing(C, gpc * GPC, width_req(i))
                su = tagged_spacing(Cup, gpc * GPC, width_req(i))
                rows.append({"case_id": c["case_id"], "variant": v, "fuel": fuel, "coverage_gpc": gpc,
                             "intensity_kw_m": i, "width_req_m": round(width_req(i), 2),
                             "spacing_m": s, "drops_per_100m": round(100 / s, 2) if s else None,
                             "litres_per_100m": round(100 / s * litres) if s else None,
                             "litres_per_100m_optimistic": round(100 / su * litres_up) if su else None,
                             "litres_basis": basis_of(c)})
    return rows


# Hotspotting (McFayden et al. 2023, Fig. 8-3, Section 8.3): cover a burning patch entirely to the
# depth that cools the combustion zone. Depth by intensity from their Table 8-3 header (C-2 fuel,
# p. 29): kW/m -> mm (= kg/m2).
HOTSPOT_DEPTH = {300: 0.5, 1300: 0.75, 3200: 1.5, 7700: 2.5}


def step_up(C):
    """The other bound of the step field: each cell at the next published level above its own."""
    lv = np.unique(C)
    nxt = dict(zip(lv, np.r_[lv[1:], lv[-1]]))
    return np.vectorize(nxt.get)(C)


def tile_spacing(C, c_req, step=2):
    """Largest tile area (m2) at which a lattice of identical drops (spacing sx along track, sy
    across) covers every point to c_req: the map is folded modulo (sx, sy). None if impossible."""
    nx, ny = C.shape
    best = None
    for kx in range(nx, 0, -step):
        for ky in range(ny, 0, -step):
            if best and kx * ky * CELL * CELL <= best[0]:
                break
            f = np.pad(C, ((0, (-nx) % kx), (0, (-ny) % ky))).reshape(-1, kx, (ny + (-ny) % ky) // ky, ky).sum(axis=(0, 2))
            if f.min() >= c_req:
                best = (kx * ky * CELL * CELL, kx * CELL, ky * CELL)
                break
    return best


def hotspot_table(cases):
    rows = []
    for c in cases:
        if "canopy" in c["case_id"] or c["case_id"].startswith("bell_s205"):
            continue
        fs, _, _ = fields(json.load(open(CX16 / c["geojson_path"])))
        C = fs["central" if "central" in fs else "lower"]
        for i, depth in HOTSPOT_DEPTH.items():
            row = {"case_id": c["case_id"], "intensity_kw_m": i, "depth_mm": depth,
                   "litres_basis": "max(released, water in the field read)"}
            for tag, field in (("certain", C), ("possible", step_up(C))):
                t = tile_spacing(field, depth)  # additive: overlapping fringes add up
                basis = max(1000 * (c.get("released_volume_m3") or 0), float(field.sum() * CELL * CELL))
                row[f"drops_per_ha_{tag}"] = round(1e4 / t[0], 1) if t else None
                row[f"litres_per_ha_{tag}"] = round(1e4 / t[0] * basis) if t else None
            # McFayden et al. (p. 23): water below the effective depth does not add; drops are
            # stacked exactly (k per spot) and their effective areas tile perfectly (a lower bound).
            areas = {k: float((C >= depth / k).sum()) * CELL * CELL for k in (1, 2, 3, 4)}
            per = {k: k * 1e4 / a for k, a in areas.items() if a > 0}
            k = min(per, key=per.get) if per else None
            row["drops_per_ha_nonadditive"] = round(per[k], 1) if k else None
            row["stack_nonadditive"] = k
            rows.append(row)
    return rows


THREE = {"blackhawk_bambi660_fig5": "Bucket: Black Hawk + Bambi 660, 42 kt*",
         "canada_bell212_fast_open_mwir_figA2_6_p72": "Belly tank: Bell 212, 40 kt (1400 L)",
         "usfs_s2t_2006_fig2_water": "Fixed wing: S-2T, 118 kt (4542 L)"}


def three_panel(sc):
    """Litres per 100 m of continuous held line, fuel x intensity, for one system per release type."""
    import matplotlib.colors as mcolors

    fuels, ints = list(FUELS), list(INTENSITIES)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
    norm = mcolors.LogNorm(vmin=500, vmax=20000)
    for ax, (cid, title) in zip(axes, THREE.items()):
        grid = np.full((len(fuels), len(ints)), np.nan)
        for r in sc:
            if r["case_id"] == cid and r["litres_per_100m"]:
                grid[fuels.index(r["fuel"]), ints.index(r["intensity_kw_m"])] = r["litres_per_100m"]
        im = ax.imshow(grid, cmap="Blues", norm=norm, aspect="auto")
        for (a_, b_), v in np.ndenumerate(grid):
            ax.text(b_, a_, "—" if np.isnan(v) else f"{v / 1000:.1f}k", ha="center", va="center", fontsize=7.5,
                    color="white" if v > 6000 else "black")
        ax.set_xticks(range(len(ints)), [f"{i}\n({width_req(i):.1f} m)" for i in ints], fontsize=7)
        ax.set_title(title, fontsize=8, loc="left")
        ax.set_xlabel("fireline intensity, kW/m (required width)", fontsize=8)
    axes[0].set_yticks(range(len(fuels)), fuels, fontsize=7.5)
    fig.colorbar(im, ax=axes, shrink=0.9, label="litres per 100 m of held line")
    fig.savefig(OUT / "three_systems_scenarios.png", dpi=160, bbox_inches="tight")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cases = water_cases()
    rows, planes = [], {}
    for c in cases:
        fs, xs, ys = fields(json.load(open(CX16 / c["geojson_path"])))
        row = {k: c.get(k) for k in ("case_id", "aircraft", "release_type", "speed_m_s", "height_m", "released_volume_m3")}
        row["canopy"] = "canopy" in c["case_id"]
        for v, C in fs.items():
            dep = float(C.sum() * CELL * CELL)  # kg (water)
            m = 1000 * c["released_volume_m3"] if c.get("released_volume_m3") else dep
            row[f"norm_{v}"] = basis_of(c)
            row[f"deposited_kg_{v}"] = round(dep, 1)
            eta = np.array([[line_yield(C, CELL, CELL, cq, w, m)["eta"] for w in W_GRID] for cq in C_GRID])
            i, j = np.unravel_index(eta.argmax(), eta.shape)
            row[f"eta_max_{v}"] = round(float(eta.max()), 3)
            row[f"argmax_C_kg_m2_{v}"] = round(float(C_GRID[i]), 2)
            row[f"argmax_W_m_{v}"] = float(W_GRID[j])
            for name, (cq, w) in POINTS.items():
                r = line_yield(C, CELL, CELL, cq, w, m)
                row[f"eta_{name}_{v}"] = round(r["eta"], 3)
                row[f"line_m_{name}_{v}"] = round(r["line_m"], 1)
            planes.setdefault(c["case_id"], {})[v] = eta
        rows.append(row)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in rows[0], k))
    with open(OUT / "measured_cases.csv", "w", newline="") as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        w.writerows(rows)
    (OUT / "measured_cases.json").write_text(json.dumps(rows, indent=1))

    n = len(rows)
    cols = 5
    fig, axes = plt.subplots(int(np.ceil(n / cols)), cols, figsize=(13, 2.6 * np.ceil(n / cols)), sharex=True, sharey=True)
    for ax, r in zip(axes.ravel(), rows):
        v = "central" if "central" in planes[r["case_id"]] else "lower"
        im = ax.pcolormesh(W_GRID, C_GRID, np.clip(planes[r["case_id"]][v], 0, 1), cmap="Blues", vmin=0, vmax=0.3, shading="auto")
        for (cq, wq), mk in zip(POINTS.values(), "^o"):
            ax.plot(wq, cq, mk, mfc="none", mec="k", ms=5)
        tag = r["case_id"].replace("canada_", "").replace("_figA2_6", "").replace("blackhawk_bambi660", "BH").replace("_", " ")
        ax.set_title(f"{tag}\n{(r['speed_m_s'] or 0) * 1.944:.0f} kt, {(r['height_m'] or 0):.0f} m, {v}", fontsize=7)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.supxlabel("required width W (m)")
    fig.supylabel("required coverage C_req (kg/m²)")
    fig.colorbar(im, ax=axes, shrink=0.5, label="line yield η")
    fig.savefig(OUT / "measured_planes.png", dpi=150, bbox_inches="tight")
    hs = hotspot_table(cases)
    with open(OUT / "hotspot_scenarios.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(hs[0]))
        w.writeheader()
        w.writerows(hs)
    sc = scenario_table(cases)
    three_panel(sc)
    with open(OUT / "tagging_scenarios.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(sc[0]))
        w.writeheader()
        w.writerows(sc)
    print(f"{n} water maps -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
