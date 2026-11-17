"""Hold metric: how long a ground pattern holds a fire, summed over blocks of fireline.

The line runs along x; the fire crosses the treated strip along +y. For each block of length l
along the line the dose profile C(y) (kg/m2, averaged over the block) is crossed by a fire front
that arrives at t_a after the drop:

- wetness s = (C - e t)+ / C*, with C* the industry coverage level for the fuel (assumed to be the
  dose that just reaches extinction moisture) and e the drying rate;
- fuel moisture relative to extinction r = r_dry + (1 - r_dry) s;
- local spread rate R = R0 f(r) / f(r_dry), f(r) = eta_M(r) / Q_ig(r m_ext): Rothermel (1972)
  moisture damping and heat of preignition Q_ig = 250 + 1116 M Btu/lb, relative to the untreated
  fuel. R = 0 for r >= 1: the front waits there until the spot has dried;
- flames reach ahead: burning fuel with spread ratio q = R / R0 reaches fuel up to w q^0.46 ahead
  of it (flame length ~ intensity^0.46, Byram 1959; an assumption), skipping what lies
  between, at its own rate R0 q; the reached fuel then burns at its own rate. The reach w of an untreated front is random, P(w > x) = P_b(x / flame length). A wet band
  therefore holds only while it is wider than the reach of the flames that meet it; slowing and
  blocking share this one width rule, so the delay never falls when the dose rises.

Scores: H_T = sum_blocks l * E[min(D, T)] (metre-seconds, capped at the hold actually needed),
the longest unbroken run of blocks that hold T with probability >= 0.5, and the efficiency
H_T / (T L_ideal), L_ideal = M0 / (W_req (C* + e (t_a + T))) for the perfect strip. With step
damping, w = W_req exactly and no drying this is Restas (2023, Fire 6:351).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def eta_m(r):
    """Rothermel (1972) moisture damping, r = fuel moisture / moisture of extinction."""
    r = np.clip(r, 0.0, 1.0)
    return np.where(r < 1.0, np.maximum(1 - 2.59 * r + 5.11 * r**2 - 3.52 * r**3, 0.0), 0.0)


@dataclass
class Scenario:
    c_star: float  # kg/m2, industry coverage level for the fuel
    flame_length: float  # m
    w_req_factor: float = 2.5  # required wet width / flame length (rule of thumb)
    spread_rate: float = 0.1  # m/s, untreated head-fire spread rate R0
    r_dry: float = 0.3  # untreated fuel moisture / moisture of extinction
    m_ext: float = 0.25  # moisture of extinction (dry-basis fraction), for Q_ig; 0 drops the term
    drying: float = 0.0  # kg/m2/s of applied dose
    arrival: float = 0.0  # s from drop to fire arrival
    damping: str = "rothermel"  # or "step": full speed below C*, blocked at or above (Restas)
    breach: str = "logistic"  # or "step": reach = W_req (Restas); "none": no reach (USFS)
    breach_w50: float = 1.5  # logistic: 50 % bridged at this width / flame length (assumed)
    breach_k: float = 3.9  # logistic slope: ~2 % bridged at 2.5 flame lengths (assumed)
    reach_exponent: float = 0.46  # reach ~ (R / R0)^this (Byram flame length); 0: no feedback

    @property
    def w_req(self) -> float:
        return self.w_req_factor * self.flame_length


def flame_reaches(sc: Scenario, n=32):
    """Reach nodes w_k (m), 0 .. 4 flame lengths, and P(w in [w_k, w_k+1)); the last node carries
    P(w >= w_K) (< 1e-4) and the first P(w <= 0). The delay falls as reach grows, so pairing each
    interval with its two end nodes brackets the expectation (midpoint quantiles
    missed the low tail)."""
    if sc.breach == "none":
        return np.array([0.0]), np.array([1.0])
    if sc.breach == "step":
        return np.array([sc.w_req]), np.array([1.0])
    w = np.linspace(0, 4, n + 1) * sc.flame_length
    cdf = 1 / (1 + np.exp(-sc.breach_k * (w / sc.flame_length - sc.breach_w50)))
    p = np.diff(np.append(cdf, 1.0))
    p[0] += cdf[0]
    return w, p


def _rate(c, t, sc: Scenario):
    """Local spread rate (m/s) of fuel holding dose c (kg/m2, as applied) at time t."""
    s = np.maximum(c - sc.drying * np.maximum(t, 0) if sc.drying else c, 0) / sc.c_star
    if sc.damping == "step":
        return np.where(s >= 1, 0.0, sc.spread_rate)

    def f(r):  # Rothermel moisture terms: damping over heat of preignition
        return eta_m(r) / (250 + 1116 * r * sc.m_ext)

    return sc.spread_rate * f(sc.r_dry + (1 - sc.r_dry) * s) / f(sc.r_dry)


def _march(c, t, dy, sc: Scenario, dt_max=5.0):
    """Time the front leaves a cell of dose c, entering at t, by local spread only."""
    t = t.copy()
    left = np.where(np.isfinite(t), dy, 0.0)
    wait = (c - sc.c_star) / sc.drying + 1e-6 if sc.drying > 0 else np.full_like(c, np.inf)
    while (act := left > 0).any():
        rate = _rate(c, t, sc)
        stuck = act & (rate <= 0)  # wait until the spot dries below extinction, then keep drying
        t = np.where(stuck, np.where(t >= wait, t + dt_max, wait), t)
        go = act & ~stuck
        dt = np.where(go, np.minimum(left / np.where(go, rate, 1), dt_max), 0)
        left = np.where(np.isfinite(t), np.maximum(left - rate * dt, 0), 0)
        t = t + dt
    return t


def block_profiles(C, dx, block):
    n = max(1, round(block / dx))
    nb = C.shape[0] // n
    return C[: nb * n].reshape(nb, n, C.shape[1]).mean(axis=1), n * dx


def delays(C, dx, dy, sc: Scenario, block=1.0, times=False):
    """Delay D (s) of every block for every reach sample, shape (n_reach, n_blocks); block length.
    times=True also returns the time the front leaves each padded cell (rows: reach-major, then blocks),
    the start time t0 of each row, the number of padding cells on each side, the time the front reaches
    each cell and, per cell, the cell whose flames lit it (-1: reached from the cell before)."""
    prof, lb = block_profiles(C, dx, block)
    w, _ = flame_reaches(sc)
    k_max = int(np.ceil(w.max() / dy)) + 1
    prof = np.tile(np.pad(prof, ((0, 0), (k_max, k_max))), (len(w), 1))  # untreated fuel each side
    rows, ny = prof.shape
    # spread ratio of each cell once it can burn, at the fire's arrival: fixed in time so that the
    # score is dose-monotone. ponytail: sources do not dry further; time-dependent q needs a
    # FIFO-safe integrator.
    q = _rate(
        np.minimum(np.maximum(prof - sc.drying * sc.arrival, 0), sc.c_star * (1 - 1e-9)), 0, sc
    )
    q = q / sc.spread_rate
    reach = np.repeat(w, rows // len(w))[:, None] * q**sc.reach_exponent
    t = np.empty((rows, ny))
    # the front reaches the first wet cell at t_a, however much dry margin the map has; the dry
    # fuel before it is crossed at R0 (13 Nov: the clock used to start at the map edge, so a wide
    # margin let the strip dry before the fire arrived)
    wet = prof > 0
    t0 = sc.arrival - np.where(wet.any(axis=1), wet.argmax(axis=1), 0) * dy / sc.spread_rate
    t_in = t0.copy()
    t_enter, src = np.empty((rows, ny)), np.full((rows, ny), -1)
    for j in range(ny):
        i = np.arange(max(0, j - k_max), j - 1)  # sources two or more cells back
        if i.size:  # flames reach this cell's near edge; it then burns at its own rate
            gap = (j - 1 - i) * dy
            ok = (gap <= reach[:, i]) & (q[:, i] > 0)
            hop = gap / (sc.spread_rate * np.where(ok, q[:, i], 1))
            tj = np.where(ok, t[:, i] + hop, np.inf)
            src[:, j] = np.where(tj.min(axis=1) < t_in, i[tj.argmin(axis=1)], -1)
            t_in = np.minimum(t_in, tj.min(axis=1))
        t_enter[:, j] = t_in
        t[:, j] = _march(prof[:, j], t_in, dy, sc)
        t_in = t[:, j]
    d = np.maximum(t[:, -1] - t0 - ny * dy / sc.spread_rate, 0)
    if times:
        return d.reshape(len(w), -1), lb, t, t0, k_max, t_enter, src
    return d.reshape(len(w), -1), lb


def hold_scores(C, dx, dy, sc: Scenario, Ts, payload_kg, block=1.0) -> list[dict]:
    """Capped hold H_T (m s), longest unbroken held run (m) and efficiency, for each T in Ts.
    H_T is the midpoint of the reach-quadrature bracket H_T_bracket_m_s."""
    d, lb = delays(C, dx, dy, sc, block)
    _, p = flame_reaches(sc)
    p = p[:, None]
    out = []
    for T in Ts:
        dc = np.minimum(d, T)
        hi, lo = (p * dc).sum(axis=0), (p * np.vstack([dc[1:], dc[-1:]])).sum(axis=0)
        expected = (hi + lo) / 2
        held = d >= T
        p_hold = ((p * held).sum(axis=0) + (p * np.vstack([held[1:], held[-1:]])).sum(axis=0)) / 2
        h_t = float(np.sum(lb * expected))
        l_ideal = payload_kg / (sc.w_req * (sc.c_star + sc.drying * (sc.arrival + T)))
        out.append(
            {
                "H_T_m_s": h_t,
                "H_T_bracket_m_s": (float(np.sum(lb * lo)), float(np.sum(lb * hi))),
                "held_run_m": float(_longest_runs(p_hold >= 0.5) * lb),
                "efficiency": h_t / (T * l_ideal),
                "L_ideal_m": l_ideal,
                "hold_s": expected,
                # line that holds: each block counts with its chance of holding the full T (no credit for partial delays)
                "p_hold": p_hold,
                "held_m": float(np.sum(lb * p_hold)),
                "held_fraction": float(np.sum(lb * p_hold)) / l_ideal,
            }
        )
    return out


def hold_score(C, dx, dy, sc: Scenario, T, payload_kg, block=1.0) -> dict:
    """hold_scores for a single T. C[ix, iy] in kg/m2."""
    return hold_scores(C, dx, dy, sc, [T], payload_kg, block)[0]


def line_yield(C, dx, dy, c_req, w_req, payload_kg) -> dict:
    """Deposition-only headline (step closure): columns along x whose longest cross-track run with
    C >= c_req is at least w_req wide. Holding for T with uniform drying e is c_req = C* + e (t_a + T).
    eta = line W c_req / payload <= deposited / payload <= 1; Restas (2023) is c_req = C*, w_req = W_req;
    IAB/USFS contour line length is w_req -> 0. C[ix, iy] in kg/m2 at the measurement's resolution."""
    run = _longest_runs(np.asarray(C) >= c_req) * dy
    ok = (run > 0) & (run >= w_req)
    line = float(ok.sum() * dx)
    return {
        "line_m": line,
        "contiguous_m": float(_longest_runs(ok) * dx),
        "eta": line * w_req * c_req / payload_kg,
    }


def _longest_runs(mask):
    """Longest run of True along the last axis, for every row (vectorised)."""
    mask = np.asarray(mask, dtype=bool)
    idx = np.arange(1, mask.shape[-1] + 1)
    last_false = np.maximum.accumulate(np.where(mask, 0, idx), axis=-1)
    return np.where(mask, idx - last_false, 0).max(axis=-1, initial=0)


def _block_chance(prof, dy, sc: Scenario):
    """Chance that the front never crosses each profile (rows), over the flame reaches. Without drying a
    cell at or above C* stops the front for good (infinite delay), so a row is blocked exactly when no
    flame jump carries the fire past its wet core. Reach nodes bracketed as in hold_scores."""
    d, _ = delays(prof, 1.0, dy, sc, block=1.0)
    blocked = np.isinf(d)
    _, p = flame_reaches(sc)
    p = p[:, None]
    return (
        (p * blocked).sum(axis=0) + (p * np.vstack([blocked[1:], blocked[-1:]])).sum(axis=0)
    ) / 2


def blocking_line(
    C,
    dx,
    dy,
    c_star,
    w_req,
    payload_kg,
    block=None,
    need_factor=1.5,
    n_need=10,
    min_p=0.0,
    **scenario,
) -> dict:
    """Blocking line (13 Nov 2026): the flame-jump model of `delays` without a clock (no drying, an
    infinite hold). A column blocks the fire unless the flames jump its wet-enough core; it counts with
    its chance of not being jumped, p_i. The dose the fuel needs is uncertain: log-uniform between
    c_star / need_factor and c_star * need_factor, one value per column, averaged with n_need midpoint
    nodes in log (need_factor 1 or n_need 1: the need is exactly c_star). Below the need a cell only
    slows the front, which is not credited. L_B = sum_i block * p_i, eta_B = L_B w_req c_star / payload,
    normalised with the nominal c_star. Flame length w_req / 2.5. min_p > 0 counts only columns with
    p_i >= min_p (min_p 0.5: columns that block more often than not; 13 Nov test, default off). p_block is
    always the raw chance. C[ix, iy] in kg/m2.

    Spread (14 Nov): the need is one unknown for the whole line (same fuel along it), shared by all columns,
    while flame jumps are independent per column. For each need node c, L(c) = sum_i block * p_i(c)
    (line_by_need). eta_B is the mean of L(c) w_req c_star / payload over the need distribution, eta_B_sd its
    standard deviation; eta_B_sd_jump is the root of the mean within-need jump variance
    sum_i block^2 p_i(c) (1 - p_i(c)), and eta_B_sd_total combines both."""
    C = np.asarray(C, float)
    block = dx if block is None else block
    prof, lb = block_profiles(C, dx, block)
    u = 2 * (np.arange(n_need) + 0.5) / n_need - 1
    needs = c_star * float(need_factor) ** u
    pc = np.zeros((n_need, prof.shape[0]))  # chance by need node (rows) and column
    top = prof.max(axis=1)
    for k, c in enumerate(needs):
        rows = top >= c  # a column with no cell at the need cannot block (finite delay)
        if rows.any():
            sc = Scenario(c_star=float(c), flame_length=w_req / 2.5, drying=0.0, **scenario)
            pc[k, rows] = _block_chance(prof[rows], dy, sc)
    p = pc.mean(axis=0)
    use = pc[:, p >= min_p]
    by_need = lb * use.sum(axis=1)  # L(c): the need is shared by all columns
    # independent flame jumps: variance within one need value, averaged over the needs
    jump_var = float(np.mean(lb**2 * (use * (1 - use)).sum(axis=1)))
    line = float(by_need.mean())
    f = w_req * c_star / payload_kg
    sd = float(by_need.std())
    return {
        "line_m": line,
        "p_block": p,
        "eta_B": line * f,
        "needs_kg_m2": needs,
        "line_by_need": by_need,
        "L_B_sd": sd,
        "eta_B_sd": sd * f,
        "eta_B_sd_jump": np.sqrt(jump_var) * f,
        "eta_B_sd_total": np.sqrt(sd**2 + jump_var) * f,
    }
