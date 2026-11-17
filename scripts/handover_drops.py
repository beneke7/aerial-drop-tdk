#!/usr/bin/env python3
"""VOF -> drops handover with secondary breakup, carried to the ground (thesis Part I, Ch. 4).

Full-scale version of the Shi et al. (2026) VOF-enhanced-DPM idea, with an exact ledger:

1. Connected liquid regions in one VOF snapshot (cells with alpha >= amin, exported by
   duplicated processor-boundary points are merged first).
   Per region: liquid volume sum(alpha*V), centroid, momentum (each parcel carries exactly these),
   equivalent diameter d_eq = (6V/pi)^(1/3) and asphericity xi = 1 - sqrt(l_min/l_max) of the
   liquid-weighted position covariance (0 = sphere, -> 1 = elongated). The region at the nozzle is
   the attached core and stays in VOF (`--core slices` hands it over as 0.25 m slices along x, a
   sensitivity case: "core assumed to detach as its current slices").
2. Shi et al.'s transfer criterion (d_eq <= d_thr and xi <= xi_c; thresholds unpublished) is
   reported as the mass fraction meeting it on a grid of mesh-scaled thresholds. For this one-way
   snapshot handover every detached region is transferred; the criterion table says how much of
   that mass Shi's rule would still keep in VOF.
3. Secondary breakup (Shi Eq. 10 critical Weber 12; Pilch & Erdman 1987 breakup time T_b): while
   We = rho_g |u_rel|^2 d / sigma stays above 12, the progress P += dt / T_b (current state)
   accumulates, and the parcel's drops split at P = 1 (reset when We drops below 12). This replaced
   the earlier clock, elapsed time >= T_b(current): on the Bambi
   snapshot the largest blobs (~0.27 m) now reach the ground before breaking (d90 2.3 -> 267 mm).
   The drop count n grows so the parcel's mass is kept (Shi Eq. 11, sum d_c^3 = d_p^3); velocity is
   kept, so momentum is exact (Shi perturb child velocities, which does not conserve momentum).
   `--no-breakup` switches it off. `--fragments` picks the child sizes:
   - stable (default): one size, the maximum stable 12 sigma / (rho_g |u_rel|^2);
   - hf: Hsiang & Faeth (1992, Int. J. Multiphase Flow 18(5):635-652, Eq. [14]), checked against
     the paper text: rho_g SMD u^2/sigma = 6.2 (rho_l/rho_g)^(1/4) [mu_l/(rho_l d0 u)]^(1/2) We, i.e.
         SMD / d0 = 6.2 (rho_l/rho_g)^(1/4) (mu_l / (rho_l d0 u))^(1/2),
     d0 the drop diameter before breakup, u the relative speed, mu_l = 1.0e-3 Pa s. Measured range:
     shock-tube air at NTP, Oh < 0.1, We < 1000 (bag, multimode, shear regimes; "catastrophic"
     We > 1e4 is beyond it), rho_l/rho_g 580-12000, Re 300-16000; coefficient s.d. 20%, r = 0.91;
     the authors ask for density-ratio tests. Children follow Simmons' root-normal volume
     distribution with MMD/SMD = 1.2 (HF: bag and multimode breakup only; shear breakup is bimodal
     because of the core drop), as 5 equal-mass quantile classes capped at d0. Mass broken up
     outside We < 1000 / Oh < 0.1 is counted in the JSON. The first split also gets the six-way
     velocity kicks on every class (5 x 6 children), along the +- axes of a random
   orthonormal frame per split (seeded), so no axis is imprinted on the ground pattern.
   A first split replaces the parcel by six children with symmetric velocity kicks.
4. Optional evaporation (`--rh --tair --tdrop`, off by default): Denner (2026) isolated-drop
   model via scripts/evaporation.py `transfer_rates`, with the drag relative speed. Per-parcel
   state: radius, temperature, statistical multiplicity n. Dry-out at r <= min(1e-7 m, 1e-3 r_ref)
   books the rest as evaporated. Ledger M0 = M_evap + M_ground + M_airborne; coverage uses M0.
5. Flight in the ground frame (crossflow speed subtracted, air at rest plus an optional crosswind
   along z): drag after Denner (2026) Eqs. 5-11, gravity, to the ground H below the nozzle. A steady
   release is uniform along track, so the cross-track landing distribution p(z) of the sampled mass
   gives the line's coverage profile c(z) = rho Q / U * p(z) (p normalised by M0, not by landed mass).

    handover_drops.py <liquid.vtu> [--U 70 --Q 1.2566 --H 60 --wind 0 --no-breakup --tag NAME
                                    --fragments hf --rh 0.3 --tair 303.15 --tdrop 303.15 --core slices]

Deterministic: kick frames come from a seeded generator; the input file's SHA-256 and the full argument list are written to
the JSON next to the results.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from statistics import NormalDist

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pyvista as pv  # noqa: E402

RHO_L, RHO_G, MU_G, SIGMA, G = 998.2, 1.18, 1.86e-5, 0.072, 9.81
MU_L = 1.0e-3  # liquid viscosity for the Hsiang-Faeth SMD (Pa s, water, given)
HF_WE_MAX, HF_OH_MAX = 1000.0, 0.1  # Hsiang & Faeth (1992) tested range (rho_l/rho_g 580-12000)
GPC = 0.4075  # kg/m2 (L/m2) per US gal / 100 ft2
CELL = 0.014  # fine-mesh cell size in the jet region (m), for mesh-scaled thresholds
SLICE = 0.25  # core slice width along x (m)
TAU_FLOOR = 2e-3  # semi-implicit drag: parcels relaxing faster than this (s) do not limit the step


def root_normal_classes(n: int = 5, ratio: float = 1.2) -> tuple[np.ndarray, float]:
    """n equal-mass quantile classes of a root-normal volume distribution (sqrt(d) normal), as
    d_k / SMD = ratio (1 + eps z_k)^2 with z_k the normal quantile at (k + 1/2)/n. eps is solved so the
    classes themselves have MMD/SMD = ratio, i.e. mean((1 + eps z)^-2) = ratio (the continuous SMD
    integral is singular at d -> 0). MMD = mu^2 is the median class; the class SMD equals the target."""
    z = np.array([NormalDist().inv_cdf((k + 0.5) / n) for k in range(n)])
    lo, hi = 0.0, 0.99 / abs(z[0])
    for _ in range(80):
        eps = 0.5 * (lo + hi)
        lo, hi = (eps, hi) if np.mean((1 + eps * z) ** -2) < ratio else (lo, eps)
    return ratio * (1 + eps * z) ** 2, eps


RN_D, RN_EPS = root_normal_classes()
_QLO = 0.1  # jittered draws use the quantile range 0.1-0.9, the span of the 5 fixed classes (z = +-1.28), without atoms
_Q400 = np.array([NormalDist().inv_cdf(_QLO + (1 - 2 * _QLO) * (k + 0.5) / 400) for k in range(400)])
_ZCLIP = NormalDist().inv_cdf(1 - _QLO)


def _jit_eps(ratio: float = 1.2) -> float:
    """eps of the clipped continuous root-normal distribution with MMD/SMD = ratio (bisection, as root_normal_classes)."""
    lo, hi = 0.0, 0.99 / _ZCLIP
    for _ in range(80):
        eps = 0.5 * (lo + hi)
        lo, hi = (eps, hi) if np.mean((1 + eps * _Q400) ** -2) < ratio else (lo, eps)
    return 0.5 * (lo + hi)


_JIT_EPS = _jit_eps()


def hf_class_sampler(n: int, jitter: bool):
    """Child-size sampler for fly(hf_sampling=...): f(rng, m) -> (m, n) diameters in units of the HF SMD.
    n equal-mass quantile classes (jitter False; n = 5 is the default RN_D), or, with jitter, the same n strata
    with a random quantile offset per split (stratified draw of the continuous root-normal distribution, so the
    size spectrum has no fixed lines; cost unchanged). The draw is the root-normal distribution truncated to the
    quantile range 0.1-0.9 (z = +-1.28, the span of the 5 fixed classes, 0.5-2.2 SMD; renormalised, no clipped atoms):
    smaller drops would only cost integration time. eps is re-solved for MMD/SMD = 1.2 of the truncated distribution."""
    if not jitter:
        d, _ = root_normal_classes(n)
        return lambda rng, m: np.broadcast_to(d, (m, n))

    def f(rng, m):
        q = _QLO + (1 - 2 * _QLO) * (np.arange(n)[None, :] + rng.random((m, 1))) / n
        z = np.vectorize(NormalDist().inv_cdf)(q)
        return 1.2 * (1 + _JIT_EPS * z) ** 2

    return f


def hf_smd(d: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Hsiang & Faeth (1992) Eq. [14] Sauter mean diameter after secondary breakup (see docstring)."""
    return d * 6.2 * (RHO_L / RHO_G) ** 0.25 * np.sqrt(MU_L / (RHO_L * d * u))


def regions(path: Path, amin: float) -> dict[str, np.ndarray]:
    grid = pv.read(path).clean(tolerance=1e-6)  # merge points duplicated at processor boundaries
    grid = grid.compute_cell_sizes(length=False, area=False)
    a = np.asarray(grid.cell_data["alpha.water"], float)
    grid = grid.extract_cells(np.where(a >= amin)[0])
    lab = np.asarray(grid.connectivity(extraction_mode="all").cell_data["RegionId"])
    a = np.asarray(grid.cell_data["alpha.water"], float)
    v = np.asarray(grid.cell_data["Volume"], float) * a  # liquid volume per cell
    x = np.asarray(grid.cell_centers().points)
    u = np.asarray(grid.cell_data["U"], float)
    n = lab.max() + 1
    vol = np.bincount(lab, v, n)
    cen = np.stack([np.bincount(lab, v * x[:, i], n) for i in range(3)], 1) / vol[:, None]
    vel = np.stack([np.bincount(lab, v * u[:, i], n) for i in range(3)], 1) / vol[:, None]
    dx = x - cen[lab]
    cov = np.zeros((n, 3, 3))
    for i in range(3):
        for j in range(i, 3):
            cov[:, i, j] = cov[:, j, i] = np.bincount(lab, v * dx[:, i] * dx[:, j], n) / vol
    ev = np.linalg.eigvalsh(cov)  # ascending
    xi = 1 - np.sqrt(np.clip(ev[:, 0], 0, None) / np.maximum(ev[:, 2], 1e-30))
    ncell = np.bincount(lab, minlength=n)
    core = lab[np.argmin(np.linalg.norm(x, axis=1))]  # region at the nozzle (origin)
    return {"vol": vol, "cen": cen, "vel": vel, "xi": xi, "ncell": ncell, "core": int(core),
            "total_alpha_volume": float(v.sum()), "cells": {"lab": lab, "v": v, "x": x, "u": u}}


def core_slices(r: dict, width: float = SLICE) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The attached core cut into bins of `width` along x: exact volume, centroid, mean velocity each."""
    c = r["cells"]
    sel = c["lab"] == r["core"]
    v, x, u = c["v"][sel], c["x"][sel], c["u"][sel]
    b = np.floor(x[:, 0] / width).astype(int)
    b -= b.min()
    n = b.max() + 1
    vol = np.bincount(b, v, n)
    cen = np.stack([np.bincount(b, v * x[:, i], n) for i in range(3)], 1)
    vel = np.stack([np.bincount(b, v * u[:, i], n) for i in range(3)], 1)
    k = vol > 0
    return vol[k], cen[k] / vol[k, None], vel[k] / vol[k, None]


DRAG_LAWS = ("hf", "tab")


def drag_coefficient(d: np.ndarray, s: np.ndarray, law: str = "hf") -> tuple[np.ndarray, np.ndarray]:
    """Drag coefficient and Weber number (diameter-based) of drops of diameter d at relative speed s.

    Both laws: rigid-sphere C_D,0 (Clift-Gauvin form) times the distorted-drop factor 1 + 2.632 y of Liu, Mather & Reitz
    (1993, SAE 930072), y the drop distortion (0 sphere, 1 disk). They differ only in y:
    - "hf" (default, every result before 17 Nov 2026): y = 1 - E^-2, E = 1 + 0.19 sqrt(We) the deformed/undeformed radius
      of Hsiang & Faeth's shock-tube drops (maximum deformation under a suddenly applied flow). For steadily falling drops
      it over-flattens them: terminal velocity 1/2/3/4 mm = 3.17/4.92/6.13/7.05 m/s, 21-24 % below Gunn & Kinzer (1949).
    - "tab" (opt-in, 17 Nov 2026): y = min(We/24, 1), the equilibrium distortion of the Taylor-analogy breakup oscillator
      (O'Rourke & Amsden 1987, SAE 872089: y = C_F/(C_k C_b) We_r = We_r/12 with C_F = 1/3, C_k = 8, C_b = 1/2;
      We_r = We/2 radius-based) under a steady load, i.e. the same Liu et al. drag with the TAB steady state in place of the
      shock-tube maximum. No constant is fitted. Terminal velocity 1/2/3/4 mm = 3.78/6.25/7.69/8.54 m/s, within -6 to -3 %
      of Gunn & Kinzer for 0.5-5.5 mm; at large We both laws tend to the same disk limit 3.63 C_D,0, so the breakup of
      large blobs is nearly unchanged.
    """
    re = RHO_G * s * d / MU_G
    we = RHO_G * s**2 * d / SIGMA
    cds = 24 / re * (1 + 0.15 * re**0.687 + 0.0175 * re / (1 + 42500 * re**-1.16))
    if law == "hf":
        rmax = 1 + 0.19 * np.sqrt(we)  # deformed / undeformed radius
        y = 1 - rmax**-2
    elif law == "tab":
        y = np.minimum(we / 24.0, 1.0)
    else:
        raise ValueError(f"unknown drag law {law!r}; one of {DRAG_LAWS}")
    return cds * (1 + 2.632 * y), we


def drag_accel(d: np.ndarray, urel: np.ndarray, law: str = "hf") -> tuple[np.ndarray, np.ndarray]:
    s = np.linalg.norm(urel, axis=1) + 1e-12
    cd, we = drag_coefficient(d, s, law)
    acc = (0.75 * cd * RHO_G / (RHO_L * d) * s)[:, None] * urel
    return acc, we


def terminal_velocity(d: float, law: str = "hf") -> float:
    """Still-air terminal velocity (m/s) of a drop of diameter d (m) under drag law `law` (bisection on drag = gravity)."""
    lo, hi = 1e-3, 40.0
    for _ in range(100):
        v = 0.5 * (lo + hi)
        cd, _ = drag_coefficient(np.array([d]), np.array([v]), law)
        lo, hi = (v, hi) if 0.75 * cd[0] * RHO_G / (RHO_L * d) * v * v < G else (lo, v)
    return 0.5 * (lo + hi)


def breakup_time(d: np.ndarray, s: np.ndarray, we: np.ndarray) -> np.ndarray:
    """Pilch & Erdman (1987) total breakup time, T = t U / d * sqrt(rho_g / rho_l)."""
    w = np.maximum(we, 12.0001)
    t = np.select(
        [w < 18, w < 45, w < 351, w < 2670],
        [6 * (w - 12) ** -0.25, 2.45 * (w - 12) ** 0.25, 14.1 * (w - 12) ** -0.25, 0.766 * (w - 12) ** 0.25],
        5.5,
    )
    return t * d / s * np.sqrt(RHO_L / RHO_G)


def load_evaporation(rh: float, tair: float, tdrop: float) -> dict:
    """Denner (2026) transfer model (scripts/evaporation.py)."""
    old = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        import evaporation as mod
    finally:
        sys.path.pop(0)
        sys.dont_write_bytecode = old
    twb = float(mod.wet_bulb_temperature_k(tair, rh))  # raises outside RH 5-99 %, 0-50 degC
    if twb < 273.15 or not 273.15 <= tdrop <= 373.15:
        raise ValueError(f"wet bulb {twb:.2f} K or drop temperature {tdrop} K outside the Denner liquid-water range")
    return {"mod": mod, "rh": rh, "tair": tair, "tdrop": tdrop, "twb": twb,
            "k_air": float(mod.air_properties(tair, rh).thermal_conductivity_w_m_k)}


def evap_step(ev: dict, r, T, speed, dt: float, rmin, nsub: int = 64):
    """Advance drop radius and temperature over dt at the drag relative speed. Per-parcel sub-steps of
    0.1 min(radius, thermal relaxation time): r^2 explicit (d^2-law, linear in r^2), T semi-implicit in
    the convective term (never more than nsub sub-steps; the forced steps only hit sub-micron drops
    that are about to dry out). Returns r, T and the dry-out flag (r <= rmin)."""
    mod = ev["mod"]
    r2, T = r**2, T.copy()
    rem, dry = np.full(len(r), dt), np.zeros(len(r), bool)
    todo = np.arange(len(r))
    for k in range(nsub):
        rr = np.sqrt(r2[todo])
        q = mod.transfer_rates(rr, T[todo], speed[todo], ev["tair"], ev["rh"])
        rd, td = q.radius_rate_m_s, q.temperature_rate_k_s
        a = 1.5 * q.nusselt_number * ev["k_air"] / (mod.RHO_WATER_KG_M3 * mod.CP_WATER_J_KG_K * rr**2)
        tr = 0.1 * np.minimum(rr / np.maximum(np.abs(rd), 1e-300), 1 / a)  # radius, thermal relaxation
        h = np.minimum(rem[todo], np.maximum(tr, rem[todo] / (nsub - k)))
        r2[todo] += 2 * rr * rd * h
        T[todo] += h * td / (1 + h * a)
        rem[todo] -= h
        gone = r2[todo] <= rmin[todo] ** 2
        dry[todo[gone]] = True
        todo = todo[~gone & (rem[todo] > 1e-12 * dt)]
        if not len(todo):
            break
    return np.sqrt(np.maximum(r2, 0.0)), T, dry


def _psi_m(zeta):
    """Integrated Businger-Dyer stability function psi_m(z/L) (Businger et al. 1971; Dyer 1974, Eq. 8-9 form):
    stable phi_m = 1 + 5 z/L -> psi = -5 z/L; unstable phi_m = (1 - 16 z/L)^(-1/4) -> psi = 2 ln((1+x)/2)
    + ln((1+x^2)/2) - 2 atan x + pi/2, x = (1 - 16 z/L)^(1/4)."""
    zeta = np.asarray(zeta, float)
    x = (1 - 16 * np.minimum(zeta, 0.0)) ** 0.25
    unst = 2 * np.log((1 + x) / 2) + np.log((1 + x**2) / 2) - 2 * np.arctan(x) + np.pi / 2
    return np.where(zeta < 0, unst, -5 * zeta)


def atmosphere(uref, direction_deg, zref, z0=0.03, turb=True, L=None, zi=1000.0, meander_deg=0.0,
               meander_tau=60.0, meander_seed=0, tmax=900.0):
    """Surface layer for fly() (bambi_bucket CASE_RECORD 13.3): wind through uref at zref (direction the wind
    blows TO, degrees from +x). Neutral (L None): log profile; sigma along/vertical/cross = 2.4/1.25/1.9 u*
    (Panofsky & Dutton 1984), T_L = 0.5 z / sigma_w (Hanna 1982). Per-parcel fluctuations without a
    crossing-trajectory correction: an upper bracket of turbulent spreading.
    L (Obukhov length, m; <0 unstable, >0 stable): Monin-Obukhov profile u = u*/k [ln(z/z0) - psi_m(z/L)
    + psi_m(z0/L)] with Businger-Dyer functions; unstable sigma_u = sigma_v = u* (12 - 0.5 zi/L)^(1/3)
    (Panofsky et al. 1977, zi ASSUMED 1000 m), sigma_w = 1.25 u* (1 - 3 z/L)^(1/3) (Panofsky & Dutton 1984);
    stable sigma = (2.4, 1.25, 1.9) u* (Panofsky & Dutton 1984, near-neutral/stable values). T_L stays
    0.5 z / sigma_w with the local sigma_w (ASSUMED; Hanna's 1982 stability forms not verified).
    meander_deg > 0: the mean direction follows a stationary Ornstein-Uhlenbeck angle (sigma meander_deg,
    time scale meander_tau s; ASSUMED values), one realisation shared by all parcels, in global time."""
    th = np.radians(direction_deg)
    ustar = 0.4 * uref / np.log(zref / z0)
    if L is not None:
        ustar = 0.4 * uref / (np.log(zref / z0) - _psi_m(zref / L) + _psi_m(z0 / L))
    basis = np.array([[np.cos(th), 0.0, np.sin(th)], [0.0, 1.0, 0.0], [-np.sin(th), 0.0, np.cos(th)]])
    atm = {"uref": uref, "zref": zref, "z0": z0, "turb": turb, "ustar": ustar, "basis": basis, "L": L,
           "zi": zi, "th": th, "sigma": ustar * np.array([2.4, 1.25, 1.9])}
    if meander_deg > 0:
        r = np.random.default_rng(1000 + meander_seed)
        n = int(tmax / 0.5) + 2
        a = np.exp(-0.5 / meander_tau)
        ang = np.empty(n)
        ang[0] = r.standard_normal()
        for k in range(1, n):
            ang[k] = a * ang[k - 1] + np.sqrt(1 - a * a) * r.standard_normal()
        atm["meander"] = np.radians(meander_deg) * ang  # angle every 0.5 s
    return atm


def wind_speed(atm, z):
    """Mean wind speed at height z (neutral: the original log expression, unchanged)."""
    if atm.get("L") is None:
        return atm["uref"] * np.log(z / atm["z0"]) / np.log(atm["zref"] / atm["z0"])
    L = atm["L"]
    f = lambda zz: np.log(zz / atm["z0"]) - _psi_m(zz / L) + _psi_m(atm["z0"] / L)  # noqa: E731
    return atm["uref"] * f(z) / f(atm["zref"])


def sigmas(atm, z):
    """Turbulence sigma (along, vertical, cross) at heights z, shape (n, 3)."""
    z = np.atleast_1d(z)
    L, us = atm.get("L"), atm["ustar"]
    if L is None or L > 0:
        return np.broadcast_to(atm["sigma"], (len(z), 3))
    su = us * (12 - 0.5 * atm["zi"] / L) ** (1 / 3)
    sw = 1.25 * us * (1 - 3 * z / L) ** (1 / 3)
    return np.stack([np.full(len(z), su), sw, np.full(len(z), su)], axis=1)


def wind_basis(atm, tglob):
    """Per-parcel (n, 3, 3) basis rows along/vertical/cross after the meander rotation at global times tglob."""
    d = np.interp(tglob, 0.5 * np.arange(len(atm["meander"])), atm["meander"]) if "meander" in atm else 0.0
    t = atm["th"] + d
    c, s, o, z = np.cos(t), np.sin(t), np.ones_like(t), np.zeros_like(t)
    return np.stack([np.stack([c, z, s], -1), np.stack([z, o, z], -1), np.stack([-s, z, c], -1)], -2)


def field_at(F: dict, pa: np.ndarray) -> np.ndarray:
    """Trilinear air-velocity perturbation of a gridded field (scripts/air_field_grid.py) at aircraft-frame
    points pa (n, 3); zero outside the grid (beyond the VOF domain the air is the undisturbed wind)."""
    f = (pa - F["origin"]) / F["h"]
    n = np.array(F["du"].shape[:3])
    inside = np.all((f >= 0) & (f <= n - 1), axis=1)
    out = np.zeros_like(pa, dtype=float)
    if inside.any():
        g = f[inside]
        i0 = np.minimum(np.floor(g).astype(int), n - 2)
        r = g - i0
        acc = np.zeros((len(g), 3))
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    wgt = (r[:, 0] if dx else 1 - r[:, 0]) * (r[:, 1] if dy else 1 - r[:, 1]) * (r[:, 2] if dz else 1 - r[:, 2])
                    acc += wgt[:, None] * F["du"][i0[:, 0] + dx, i0[:, 1] + dy, i0[:, 2] + dz]
        out[inside] = acc
    return out


def fly(d, pos, vel, H, breakup: bool, wind: float = 0.0, tmax=600.0, spread: bool = True,
        fragments: str = "stable", evap: dict | None = None, info: dict | None = None, semi: bool = False,
        atm: dict | None = None, record: dict | None = None, field: dict | None = None, hf_sampling=None,
        drag: str = "hf", turb_seed: int = 1):
    """Track parcels to the ground. Returns landing point, landing time, drop size, drops per parcel
    and the parent index of every (child) parcel, so masses can be mapped back (m_child = m_parent*frac).
    `evap` (from load_evaporation) switches evaporation on; a dict passed as `info` is filled with the
    extra per-parcel state (w = remaining mass / mass at handover share, T, alive, replaced-by-children,
    breakup generation, outside-HF-range lineage) and the breakup event log. `semi` (always on with
    `evap`) integrates the drag semi-implicitly: shrinking drops would otherwise set dt -> 1e-5 s.
    `atm` (atmosphere(), optional) replaces the uniform `wind` by a neutral log profile and, with
    turb, an Ornstein-Uhlenbeck fluid velocity seen by each parcel (children inherit it).
    `record` ({"dt": s}) collects every record["dt"] of flight time the alive parcels' parent index, position,
    diameter and mass share (float32) in record["snaps"], for visualization only.
    `field` ({"F": gridded perturbation, "t0": global start time per input parcel, "v_frame": ground velocity
    of the aircraft frame}) adds the time-mean VOF air field (rotor wake, air dragged by the liquid) to the
    wind; it moves with the aircraft, so a parcel at ground point p and global time t0 + t sees it at
    p - v_frame (t0 + t). Superposition on the analytic wind is an assumption.
    `hf_sampling` (hf_class_sampler; default None = the 5 fixed RN_D classes) sets the child size classes.
    `drag` picks the drag law (drag_coefficient: "hf" default, "tab" opt-in free-fall distortion). `turb_seed` seeds the
    turbulence stream only (default 1 = every map before 17 Nov 2026); kick frames and size jitter keep their seeds."""
    if drag not in DRAG_LAWS:
        raise ValueError(f"unknown drag law {drag!r}")
    semi = semi or evap is not None
    air = np.array(wind, float) if np.ndim(wind) else np.array([0.0, 0.0, wind])  # scalar: crosswind along z
    hf = fragments == "hf"
    N = len(d)
    S = {"d": np.array(d, float), "pos": np.array(pos, float), "vel": np.array(vel, float),
         "n": np.ones(N), "frac": np.ones(N), "parent": np.arange(N),  # frac: share of the parent's mass
         "split_once": np.zeros(N, bool), "progress": np.zeros(N), "alive": np.ones(N, bool),
         "land": np.full((N, 3), np.nan), "tland": np.full(N, np.nan), "w": np.ones(N),
         "T": np.full(N, evap["tdrop"] if evap else np.nan), "rmin": np.minimum(1e-7, 5e-4 * np.asarray(d, float)),
         "gen": np.zeros(N, int), "oor": np.zeros(N, bool), "repl": np.zeros(N, bool),
         "up": np.zeros((N, 3))}  # turbulent fluid velocity seen, in (along-wind, vertical, cross-wind)
    kicks6 = np.vstack([np.eye(3), -np.eye(3)])
    rng = np.random.default_rng(0)  # seeded: kick frames are random but reproducible
    srng = np.random.default_rng(2)  # own stream for jittered size classes (kick frames unchanged)
    kicks1 = np.zeros((1, 3))
    log = []
    trng = np.random.default_rng(turb_seed)  # turbulence: its own seeded stream, so the kick frames are unchanged
    if atm is not None and atm["turb"]:  # start in the stationary state: T_L ~ 100 s, so a zero start never spins up
        S["up"] = sigmas(atm, np.maximum(H + S["pos"][:, 1], atm["z0"])) * trng.standard_normal((N, 3))

    def air_at(i):
        if atm is None:
            base = air
        else:
            z = np.maximum(H + S["pos"][i, 1], atm["z0"])  # height above ground
            if "meander" in atm:
                bs = wind_basis(atm, atm["t0"][S["parent"][i]] + t)
                base = wind_speed(atm, z)[:, None] * bs[:, 0] + np.einsum("nk,nkj->nj", S["up"][i], bs)
            else:
                base = np.outer(wind_speed(atm, z), atm["basis"][0]) + S["up"][i] @ atm["basis"]
        if field is None:
            return base
        tg = field["t0"][S["parent"][i]] + t
        return base + field_at(field["F"], S["pos"][i] - np.outer(tg, field["v_frame"]))

    def spawn(js, dcls, sj, kicks):
        """Replace parcels js by len(dcls[0]) size classes x len(kicks) children; mass and momentum exact."""
        nc, nk = dcls.shape[1], len(kicks)
        k = nc * nk
        # symmetric velocity kicks +-delta along x, y, z cancel, so the parent's momentum is exact.
        # delta = |u_rel| sqrt(rho_g/rho_l), the expansion-speed scale of a breakup debris cloud (model parameter)
        delta = sj * np.sqrt(RHO_G / RHO_L)
        if kicks.shape[0] == 6:  # a random orthonormal frame per split: no imprinted axes, momentum exact
            q = np.linalg.qr(rng.standard_normal((len(js), 3, 3)))[0]
            kj = np.einsum("kd,jed->jke", kicks, q)  # (len(js), 6, 3): +-columns of each frame
        else:
            kj = np.broadcast_to(kicks, (len(js),) + kicks.shape)
        cv = S["vel"][js][:, None, None, :] + delta[:, None, None, None] * kj[:, None]
        rep = lambda a: np.repeat(a, k, axis=0)  # noqa: E731
        new = {key: rep(S[key][js]) for key in S}
        new["d"] = dcls.repeat(nk, axis=1).reshape(-1)
        new["vel"] = np.broadcast_to(cv, (len(js), nc, nk, 3)).reshape(-1, 3)
        new["n"] = rep(S["n"][js]) * (rep(S["d"][js]) / new["d"]) ** 3 / k
        new["frac"] = rep(S["frac"][js]) / k
        new["rmin"] = np.minimum(1e-7, 5e-4 * new["d"])
        new["split_once"], new["progress"], new["alive"] = np.ones(len(new["d"]), bool), np.zeros(len(new["d"])), np.ones(len(new["d"]), bool)
        new["repl"], new["land"], new["tland"] = np.zeros(len(new["d"]), bool), np.full((len(new["d"]), 3), np.nan), np.full(len(new["d"]), np.nan)
        S["alive"][js], S["repl"][js] = False, True  # parent replaced by its children
        for key in S:
            S[key] = np.concatenate([S[key], new[key]])

    t = 0.0
    while S["alive"].any() and t < tmax:
        d, pos, vel, alive = S["d"], S["pos"], S["vel"], S["alive"]
        i = np.where(alive)[0]
        ai = air_at(i)
        urel = ai - vel[i]
        acc, we = drag_accel(d[i], urel, drag)
        s0 = np.linalg.norm(urel, axis=1)
        if breakup:
            # Explicit left-endpoint quadrature: P += dt / Tb(start state).
            over_start = we > 12.0
            tb_start = np.ones(len(i), dtype=float)
            tb_start[over_start] = breakup_time(d[i][over_start], s0[over_start], we[over_start])
        # step = a tenth of the fastest drag relaxation time (explicit Euler stays stable)
        tau = s0 / (np.linalg.norm(acc, axis=1) + 1e-12)
        if semi:  # dt resolves how fast |u_rel| changes (0 at terminal velocity); parcels with tau < TAU_FLOOR
            # follow their terminal velocity implicitly and do not set dt
            tr = (s0 / (np.linalg.norm(acc + np.array([0.0, -G, 0.0]), axis=1) + 1e-12))[tau >= TAU_FLOOR]
            dt = float(np.clip(0.1 * tr.min(), 1e-5, 5e-3)) if len(tr) else 5e-3
            k = (np.linalg.norm(acc, axis=1) / (s0 + 1e-12))[:, None]  # = 1/tau; acc = k (air - vel); vel' = (vel + dt (k air + g)) / (1 + dt k)
            vel[i] = (vel[i] + dt * (k * ai + np.array([0.0, -G, 0.0]))) / (1 + dt * k)
        else:
            dt = float(np.clip(0.1 * tau.min(), 1e-5, 5e-3))
            acc[:, 1] -= G
            vel[i] += acc * dt
        pos[i] += vel[i] * dt
        t += dt
        if atm is not None and atm["turb"]:  # OU with T_L(z); sigma is height-independent, so no drift term
            z = np.maximum(H + pos[i, 1], atm["z0"])
            sg = sigmas(atm, z)
            f = np.exp(-dt / (0.5 * z / sg[:, 1]))[:, None]
            S["up"][i] = S["up"][i] * f + sg * np.sqrt(1 - f**2) * trng.standard_normal((len(i), 3))
        if evap:  # drag above used the diameter at the start of the step; same relative speed for Sh and Nu
            r, T, dry = evap_step(evap, d[i] / 2, S["T"][i], s0, dt, S["rmin"][i])
            S["w"][i] = np.where(dry, 0.0, S["w"][i] * (r / (d[i] / 2)) ** 3)
            S["T"][i], d[i] = T, np.where(dry, d[i], 2 * r)
            alive[i[dry]] = False  # dry-out: the rest of the parcel is booked as evaporated
            if breakup:
                over_start, tb_start = over_start[~dry], tb_start[~dry]
            i, we = i[~dry], we[~dry]
            if np.ndim(ai) == 2:
                ai = ai[~dry]
        if breakup and len(i):
            s = np.linalg.norm(ai - vel[i], axis=1)
            over_end = (RHO_G * s**2 * d[i] / SIGMA) > 12.0
            ground_end = pos[i, 1] <= -H
            eligible = over_start & over_end & ~ground_end
            p_next = np.where(eligible, S["progress"][i] + dt / tb_start, 0.0)
            S["progress"][i] = p_next
            split = eligible & (p_next >= 1.0)
            if split.any():
                j, sj = i[split], s[split]
                wej = RHO_G * sj**2 * d[j] / SIGMA
                out = (wej > HF_WE_MAX) | (MU_L / np.sqrt(RHO_L * d[j] * SIGMA) > HF_OH_MAX)
                first = ~S["split_once"][j]
                log.append((S["parent"][j], S["frac"][j] * S["w"][j], first, out, wej))
                S["gen"][j] += 1
                S["oor"][j] |= out
                if hf:  # 5 root-normal classes around the Hsiang-Faeth SMD, capped at the parent size
                    rn = RN_D[None, :] if hf_sampling is None else hf_sampling(srng, len(j))
                    dcls = np.minimum(rn * hf_smd(d[j], sj)[:, None], d[j][:, None])
                else:
                    dcls = np.minimum(12 * SIGMA / (RHO_G * sj**2), d[j])[:, None]
                kick = spread & first  # six-child kicks on the first split only
                new = kick | (dcls.shape[1] > 1)
                S["progress"][j] = 0.0
                stay = j[~new]  # no children: the drops shrink in place, the drop count grows
                if len(stay):
                    dn = dcls[~new, 0]
                    S["n"][stay] *= (d[stay] / dn) ** 3
                    d[stay], S["rmin"][stay] = dn, np.minimum(1e-7, 5e-4 * dn)
                for sel, kk in ((new & kick, kicks6), (new & ~kick, kicks1)):
                    if sel.any():
                        spawn(j[sel], dcls[sel], sj[sel], kk)
        alive, pos = S["alive"], S["pos"]
        hit = np.where(alive & (pos[:, 1] <= -H))[0]
        S["land"][hit], S["tland"][hit], alive[hit] = pos[hit], t, False
        while record is not None and t >= record.setdefault("next", 0.0):
            j = np.where(alive)[0]
            record.setdefault("snaps", []).append((record["next"], S["parent"][j], pos[j].astype(np.float32),
                                                   S["d"][j].astype(np.float32), (S["frac"][j] * S["w"][j]).astype(np.float32)))
            record["next"] += record["dt"]
    if info is not None:
        cat = lambda k: np.concatenate([e[k] for e in log]) if log else np.zeros(0)  # noqa: E731
        info.update({k: S[k] for k in ("w", "T", "alive", "repl", "gen", "oor")}, t_end=t,
                    events={k: cat(c) for c, k in enumerate(("parent", "mass_frac", "first", "out", "we"))})
    return S["land"], S["tland"], S["d"], S["n"], S["frac"], S["parent"]


def book(m, parent, frac, tl, info, group=None) -> dict:
    """Mass ledger (kg) of the leaf parcels: M0 = M_evap + M_ground + M_airborne (+ outside-grid, 0 here)."""
    leaf = ~info["repl"]
    mp = m[parent] * frac
    mw = mp * info["w"]
    landed = leaf & np.isfinite(tl)
    air = leaf & info["alive"]
    out = {"M0_kg": float(m.sum()), "M_evaporated_kg": float((mp - mw)[leaf].sum()),
           "M_ground_kg": float(mw[landed].sum()), "M_ground_outside_grid_kg": 0.0,
           "M_airborne_at_tmax_kg": float(mw[air].sum())}
    out["closure_rel_error"] = abs(out["M_evaporated_kg"] + out["M_ground_kg"] + out["M_airborne_at_tmax_kg"] - out["M0_kg"]) / out["M0_kg"]
    if group is not None:  # per source group (0 detached regions, 1 core slices)
        out["by_source"] = {}
        for g, name in enumerate(("detached", "core_slices")):
            sel, m0 = group[parent] == g, float(m[group == g].sum())
            if m0 > 0:
                e, gr, ar = (float(x[sel].sum()) for x in ((mp - mw) * leaf, mw * landed, mw * air))
                out["by_source"][name] = {"M0_kg": m0, "M_evaporated_kg": e, "M_ground_kg": gr, "M_airborne_at_tmax_kg": ar,
                                          "closure_rel_error": abs(e + gr + ar - m0) / m0}
    return out


def quantiles(d: np.ndarray, m: np.ndarray) -> dict[str, float]:
    o = np.argsort(d)
    cm = np.cumsum(m[o]) / m.sum()
    return {f"d{q}_mm": float(np.interp(q / 100, cm, d[o]) * 1e3) for q in (10, 50, 90)}


def wmedian(x: np.ndarray, w: np.ndarray) -> float:
    o = np.argsort(x)
    return float(np.interp(0.5, np.cumsum(w[o]) / w.sum(), x[o]))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vtu", type=Path)
    ap.add_argument("--amin", type=float, default=0.05)
    ap.add_argument("--U", type=float, default=70.0, help="crossflow (flight) speed, m/s")
    ap.add_argument("--Q", type=float, default=1.2566, help="release rate, m3/s")
    ap.add_argument("--H", type=float, default=60.0, help="nozzle height above ground, m")
    ap.add_argument("--wind", type=float, default=0.0, help="crosswind along z, m/s")
    ap.add_argument("--no-breakup", action="store_true")
    ap.add_argument("--no-spread", action="store_true", help="breakup without the six-child velocity kicks")
    ap.add_argument("--integrator", choices=("euler", "semi"), default="euler",
                    help="drag step: explicit Euler (original) or semi-implicit (many small children, evaporation)")
    ap.add_argument("--fragments", choices=("stable", "hf"), default="stable",
                    help="child sizes: maximum stable size, or Hsiang-Faeth SMD with root-normal classes")
    ap.add_argument("--rh", type=float, default=None, help="relative humidity (fraction); switches evaporation on")
    ap.add_argument("--tair", type=float, default=None, help="air temperature, K (with --rh)")
    ap.add_argument("--tdrop", type=float, default=None, help="initial drop temperature, K (default: --tair)")
    ap.add_argument("--core", choices=("none", "slices"), default="none",
                    help="also hand over the attached core as 0.25 m slices along x (sensitivity case)")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--out", type=Path, default=Path("results/runs/handover"))
    a = ap.parse_args()
    tag = a.tag or ("nobreakup" if a.no_breakup else "breakup")
    ev = None
    if a.rh is not None:
        if a.tair is None:
            ap.error("--rh needs --tair")
        try:
            ev = load_evaporation(a.rh, a.tair, a.tdrop if a.tdrop is not None else a.tair)
        except ValueError as e:
            ap.error(str(e))

    r = regions(a.vtu, a.amin)
    det = np.arange(len(r["vol"])) != r["core"]
    vol, cen, xi = r["vol"][det], r["cen"][det], r["xi"][det]
    vel = r["vel"][det].copy()
    ndet = len(vol)
    if a.core == "slices":
        cv, cc, cu = core_slices(r)
        vol, cen, vel = np.r_[vol, cv], np.vstack([cen, cc]), np.vstack([vel, cu])
    group = np.r_[np.zeros(ndet, int), np.ones(len(vol) - ndet, int)]
    vel[:, 0] -= a.U  # aircraft frame -> ground frame (air at rest)
    d0 = (6 * vol / np.pi) ** (1 / 3)
    m = RHO_L * vol
    p0 = (m[:, None] * vel).sum(0)  # ledger: momentum at handover (ground frame)
    info: dict = {}
    land, tl, d, n, frac, parent = fly(d0.copy(), cen.copy(), vel.copy(), a.H, not a.no_breakup, a.wind,
                                       spread=not a.no_spread, fragments=a.fragments, evap=ev, info=info, semi=a.integrator == "semi")
    mp = m[parent] * frac  # mass of every parcel at its handover share
    mg = mp * info["w"]  # mass left when it lands (evaporation only)
    ok = np.isfinite(tl)
    cen_p = cen[parent]
    led = book(m, parent, frac, tl, info, group if a.core == "slices" else None)

    shi = {}
    for k in (4, 8, 16):
        for xc in (0.3, 0.5, 0.7):
            sel = (d0[:ndet] <= k * CELL) & (xi <= xc)
            shi[f"d_thr={k}cells({k * CELL * 1e3:.0f}mm),xi_c={xc}"] = round(float(m[:ndet][sel].sum() / m[:ndet].sum()), 4)

    # bins cover every landing point (fixed 0.5 m width): a fixed window clipped the crosswind case
    zlo, zhi = np.floor(land[ok, 2].min()) - 1, np.ceil(land[ok, 2].max()) + 1
    edges = np.arange(zlo, zhi + 0.01, 0.5)
    p, _ = np.histogram(land[ok, 2], edges, weights=mg[ok])
    p = p / m.sum() / np.diff(edges)  # per metre cross-track, normalised by M0 (evaporated/airborne mass stays out)
    c = RHO_L * a.Q / a.U * p  # kg/m2 for a steady line
    evs = info["events"]
    hf_range = None
    if len(evs["parent"]):
        me = evs["mass_frac"] * m[evs["parent"].astype(int)]
        fi, out = evs["first"].astype(bool), evs["out"].astype(bool)
        hf_range = {
            "range": f"We < {HF_WE_MAX:g}, Oh < {HF_OH_MAX}, rho_l/rho_g 580-12000 (here {RHO_L / RHO_G:.0f}); counted for every closure",
            "breakup_events": int(len(me)),
            "mass_outside_range_fraction_all_events": float(me[out].sum() / me.sum()),
            "first_breakups": int(fi.sum()),
            "mass_outside_range_fraction_first_breakups": float(me[fi & out].sum() / me[fi].sum()),
            "first_breakup_We_mass_median": wmedian(evs["we"][fi], me[fi]),
            "ground_mass_fraction_with_outside_range_ancestor": float(mg[ok & info["oor"]].sum() / mg[ok].sum()),
            "ground_mass_weighted_mean_breakup_generations": float((mg[ok] * info["gen"][ok]).sum() / mg[ok].sum()),
        }
    res = {
        "command": " ".join(sys.argv),
        "input_sha256": hashlib.sha256(a.vtu.read_bytes()).hexdigest(),
        "ledger": {
            "liquid_in_exported_cells_m3": r["total_alpha_volume"],
            "attached_core_m3": float(r["vol"][r["core"]]),
            "detached_regions": int(det.sum()),
            "detached_m3": float(vol[:ndet].sum()),
            "handover_momentum_ground_frame_kg_m_s": [float(v) for v in p0],
            "landed_fraction_of_detached_mass": float(mg[ok].sum() / m.sum()),
            "parcels_tracked": int(ok.sum()),
        },
        "closure": {"fragments": a.fragments, "drag_integrator": "semi-implicit" if a.integrator == "semi" or ev else "explicit Euler", "breakup": not a.no_breakup, "six_child_kicks": not a.no_spread,
                    **({"hf_mu_l_Pa_s": MU_L, "hf_classes": 5, "hf_MMD_over_SMD": 1.2, "root_normal_eps": float(RN_EPS),
                        "hf_class_d_over_SMD": [float(x) for x in RN_D]} if a.fragments == "hf" else {})},
        "mass_ledger": {**led, "coverage_normalised_by": "M0 (all handed-over mass)",
                        "mass_fractions_of_M0": {k: led[f"M_{k}_kg"] / led["M0_kg"] for k in ("evaporated", "ground", "airborne_at_tmax")}},
        "handover_d_eq_mass_weighted": quantiles(d0, m),
        "handover_mass_fraction_meeting_shi_criterion": shi,
        "ground_drop_size_mass_weighted": quantiles(d[ok], mg[ok]),
        "flight_time_s_mass_median": wmedian(tl[ok], mg[ok]),
        "cross_track_drift_m_mass_median": wmedian(land[ok, 2], mg[ok]),
        "along_track_from_release_m_mass_median": wmedian(land[ok, 0] - cen_p[ok, 0], mg[ok]),
        "cross_track_width_m_at_GPC": {str(g): float((c >= g * GPC).sum() * 0.5) for g in (0.5, 1, 2, 4, 8)},
        "peak_coverage_GPC": float(c.max() / GPC),
        "coverage_integral_check_kg_per_m": [float((c * np.diff(edges)).sum()), float(RHO_L * a.Q / a.U * led["M_ground_kg"] / led["M0_kg"])],
        "assumptions": f"one snapshot as a sample of a steady release; air at rest + crosswind {a.wind} m/s; "
        + ("no evaporation" if ev is None else
           f"evaporation ASSUMED run values RH = {a.rh}, T_air = {a.tair} K, T_drop = {ev['tdrop']} K (not measured)")
        + f"; uniform crossflow {a.U} m/s; H = {a.H} m; Q = {a.Q} m3/s; amin = {a.amin}"
        + ("; core assumed to detach as its current 0.25 m slices (sensitivity case)" if a.core == "slices" else ""),
    }
    if hf_range:
        res["breakup_range_check"] = hf_range
    if ev:
        res["evaporation"] = {
            "model": "Denner (2026) Eqs. 12-33, scripts/evaporation.py transfer_rates (isolated spherical drop, no spray shielding)",
            "status": "run values are ASSUMED, not measured: RH, T_air, T_drop = T_air",
            "relative_humidity": a.rh, "T_air_K": a.tair, "T_drop_initial_K": ev["tdrop"], "wet_bulb_K": ev["twb"],
            "rho_handover_kg_m3": RHO_L, "rho_denner_kg_m3": float(ev["mod"].RHO_WATER_KG_M3), "rho_coverage_kg_m3": RHO_L,
            "dry_out": "r <= min(1e-7 m, 1e-3 r_ref), r_ref = radius when the parcel was created; rest booked as evaporated",
            "substeps": "per parcel 0.1 min(radius, thermal relaxation time), at most 64 per motion step",
            "mass_fraction_of_M0_in_parcels_that_dried_out": float((mp * (~info["repl"] & ~ok & ~info["alive"])).sum() / m.sum()),
            "ground_drop_T_mass_median_K": wmedian(info["T"][ok], mg[ok]),
        }
    if a.core == "slices":
        res["core_slices"] = {"slice_width_m": SLICE, "slices": int((group == 1).sum()), "core_m3": float(vol[ndet:].sum()),
                              "label": "sensitivity case: core assumed to detach as its current slices; 100% of the liquid is handed over"}
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / f"handover_{tag}.json").write_text(json.dumps(res, indent=1))
    np.savetxt(a.out / f"cross_track_{tag}.csv", np.c_[0.5 * (edges[1:] + edges[:-1]), c / GPC],
               delimiter=",", header="z_m,coverage_GPC (steady line)", comments="")
    de = np.geomspace(1e-3, 1e3, 121)  # ground spectrum: cumulative mass below d (mm)
    cum = np.array([mg[ok][d[ok] <= x * 1e-3].sum() for x in de])
    np.savetxt(a.out / f"ground_spectrum_{tag}.csv", np.c_[de, cum / mg[ok].sum(), cum / m.sum()], delimiter=",", fmt="%.6g",
               header="d_mm,cumulative_fraction_of_landed_mass,cumulative_fraction_of_M0", comments="")

    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    shaped = r["ncell"][det] >= 8  # shape is meaningless for regions of a few cells
    ax[0, 0].scatter(d0[:ndet][shaped] * 1e3, xi[shaped], s=4 + 400 * m[:ndet][shaped] / m.max(), c="#2a78d6", alpha=0.6, edgecolors="none")
    ax[0, 0].set(xscale="log", xlabel="handover equivalent diameter (mm)", ylabel="asphericity xi",
                 title="Detached VOF regions with >= 8 cells (size = mass)")
    for k in (4, 8, 16):
        ax[0, 0].axvline(k * CELL * 1e3, color="#8a8984", lw=1, ls=":")
    o = np.argsort(d0)
    ax[0, 1].plot(d0[o] * 1e3, np.cumsum(m[o]) / m.sum(), color="#0b0b0b", label="at handover")
    o = np.argsort(d[ok])
    ax[0, 1].plot(d[ok][o] * 1e3, np.cumsum(mg[ok][o]) / mg[ok].sum(), color="#eb6834", label="at ground")
    ax[0, 1].set(xscale="log", xlabel="drop diameter (mm)", ylabel="cumulative mass fraction",
                 title="Drop size, mass-weighted")
    ax[0, 1].legend(frameon=False)
    zc = 0.5 * (edges[1:] + edges[:-1])
    ax[1, 0].plot(zc, c / GPC, color="#1baf7a")
    ax[1, 0].set(xlabel="cross-track z (m)", ylabel="coverage (GPC)",
                 title=f"Steady line, Q = {a.Q} m3/s at {a.U} m/s, H = {a.H} m")
    ax[1, 1].hist(tl[ok], bins=40, weights=mg[ok] / mg[ok].sum(), color="#eda100")
    ax[1, 1].set(xlabel="flight time to ground (s)", ylabel="mass fraction", title="Flight time")
    fig.suptitle(f"VOF -> drops handover ({tag}): fine A, t = 1.5 s. Computed; assumptions in the JSON.")
    fig.tight_layout()
    fig.savefig(a.out / f"handover_{tag}.png", dpi=130)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
