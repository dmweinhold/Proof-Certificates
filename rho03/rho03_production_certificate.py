#!/usr/bin/env python3
r"""
rho03_fluid_limit_certificate_picard.py

Validated Picard step-tube certificate for the rho = 0.3 fluid-limit ODE.

The script encloses the hybrid ODE in tail coordinates

    u = 1 - Phi(x),    v = 1 - Phi(c),

using one-sided normal-CDF bounds, verified inverse-normal brackets, IEEE
outward rounding, and a Picard step tube at every integration step.  Each step
constructs a tube T_n satisfying

    B_n + [0,h] F(T_n) \subseteq T_n,

and evaluates all drift, budget, Green-cost and ecological-value integrands on
T_n, not on the endpoint box B_n.

The validated state includes interval enclosures for

    u, v, b_F, b_G, C_G, V_out,

where C_G is Green agricultural cost accumulated from 0 to kappa_low and
V_out is ecological value accumulated only over [m, kappa_low].
"""

from __future__ import annotations

import json
import math
import statistics
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()
OUTPUT_STEM = SCRIPT_PATH.stem
REPORT_KIND = "PRODUCTION" if "production" in OUTPUT_STEM else "FLUID"

norm = statistics.NormalDist()

# ----------------------------------------------------------------------
# Frozen theorem/certificate constants
# ----------------------------------------------------------------------
rho = 0.3
sigma = math.sqrt(1.0 - rho * rho)
ell = 0.1
nu = 10.0
d = nu - ell

m = 0.10
kappa_low = 0.620816
# Coarser threshold supported by the revised paper's theta identity.
theta_03 = 2.36413
budget = 2.525

# Initial diagonal bootstrap and Picard step size.
s0 = 1.0e-5
h = 1.0e-4

# Truncated-normal quadrature for G_A.
GA_T_MIN = -8.0
GA_N = 96

# Picard construction parameters.
PICARD_MAX_ITER = 20
PICARD_INFLATION = 1.25

# ----------------------------------------------------------------------
# Validated Phi / inverse-Phi infrastructure
# ----------------------------------------------------------------------
P = 0.2316419
B1 = 0.319381530
B2 = -0.356563782
B3 = 1.781477937
B4 = -1.821255978
B5 = 1.330274429
PHI_APPROX_ERR = 7.5e-8
ROUND_PAD = 1.0e-12
SQRT_2PI = math.sqrt(2.0 * math.pi)
P_MIN = 1.0e-16
P_MAX = 1.0 - 1.0e-16


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def clamp_prob(p: float) -> float:
    if p <= P_MIN:
        return P_MIN
    if p >= P_MAX:
        return P_MAX
    return p


def phi_pdf(x: float) -> float:
    return math.exp(-0.5 * x * x) / SQRT_2PI


@lru_cache(maxsize=500_000)
def Phi_bounds_float(x: float) -> Tuple[float, float]:
    """Return (lo, hi) such that lo <= Phi(x) <= hi."""
    x = float(x)
    if x >= 0.0:
        t = 1.0 / (1.0 + P * x)
        poly = (((((B5 * t) + B4) * t + B3) * t + B2) * t + B1) * t
        approx = 1.0 - phi_pdf(x) * poly
        lo = max(0.0, approx - PHI_APPROX_ERR - ROUND_PAD)
        hi = min(1.0, approx + PHI_APPROX_ERR + ROUND_PAD)
        return down(lo), up(hi)
    lo_pos, hi_pos = Phi_bounds_float(-x)
    return max(0.0, down(1.0 - hi_pos)), min(1.0, up(1.0 - lo_pos))


@lru_cache(maxsize=500_000)
def invPhi_lower(p: float) -> float:
    """Return x with Phi(x) <= p, verified by Phi_bounds_float."""
    p = clamp_prob(float(p))
    center_p = max(P_MIN, p - (PHI_APPROX_ERR + ROUND_PAD))
    x = norm.inv_cdf(center_p)
    _, hi = Phi_bounds_float(x)
    if hi <= p:
        return down(x)

    fail = x
    step = 1.0e-6
    while True:
        cand = x - step
        _, hi = Phi_bounds_float(cand)
        if hi <= p:
            good = cand
            break
        fail = cand
        step *= 2.0
        if cand < -50.0:
            raise RuntimeError(f"invPhi_lower failed for p={p}")
        x = cand

    lo_x = good
    hi_x = fail
    for _ in range(80):
        mid = (lo_x + hi_x) / 2.0
        _, hi_mid = Phi_bounds_float(mid)
        if hi_mid <= p:
            lo_x = mid
        else:
            hi_x = mid
    return down(lo_x)


@lru_cache(maxsize=500_000)
def invPhi_upper(p: float) -> float:
    """Return x with Phi(x) >= p, verified by Phi_bounds_float."""
    p = clamp_prob(float(p))
    center_p = min(P_MAX, p + (PHI_APPROX_ERR + ROUND_PAD))
    x = norm.inv_cdf(center_p)
    lo, _ = Phi_bounds_float(x)
    if lo >= p:
        return up(x)

    fail = x
    step = 1.0e-6
    while True:
        cand = x + step
        lo, _ = Phi_bounds_float(cand)
        if lo >= p:
            good = cand
            break
        fail = cand
        step *= 2.0
        if cand > 50.0:
            raise RuntimeError(f"invPhi_upper failed for p={p}")
        x = cand

    lo_x = fail
    hi_x = good
    for _ in range(80):
        mid = (lo_x + hi_x) / 2.0
        lo_mid, _ = Phi_bounds_float(mid)
        if lo_mid >= p:
            hi_x = mid
        else:
            lo_x = mid
    return up(hi_x)


def q_interval_from_tail(uL: float, uU: float) -> Tuple[float, float]:
    """For u in [uL,uU], q(u)=Phi^{-1}(1-u) lies in [xL,xU]."""
    p_lo = clamp_prob(1.0 - uU)
    p_hi = clamp_prob(1.0 - uL)
    return invPhi_lower(p_lo), invPhi_upper(p_hi)


# ----------------------------------------------------------------------
# Drift and integrand interval arithmetic
# ----------------------------------------------------------------------

def F1_interval_point(u: float, v: float) -> Tuple[float, float]:
    px = clamp_prob(1.0 - u)
    pc = clamp_prob(1.0 - v)
    x_lo = invPhi_lower(px)
    x_hi = invPhi_upper(px)
    c_lo = invPhi_lower(pc)
    c_hi = invPhi_upper(pc)
    z_lo = down((c_lo - rho * x_hi) / sigma)
    z_hi = up((c_hi - rho * x_lo) / sigma)
    lo_zlo, _ = Phi_bounds_float(z_lo)
    _, hi_zhi = Phi_bounds_float(z_hi)
    if lo_zlo <= 0.0 or hi_zhi <= 0.0:
        raise RuntimeError("Nonpositive Phi bound in F1_interval_point")
    return down(1.0 / hi_zhi), up(1.0 / lo_zlo)


def F2_interval_point(u: float, v: float) -> Tuple[float, float]:
    px = clamp_prob(1.0 - u)
    pc = clamp_prob(1.0 - v)
    x_lo = invPhi_lower(px)
    x_hi = invPhi_upper(px)
    c_lo = invPhi_lower(pc)
    c_hi = invPhi_upper(pc)
    z_lo = down((x_lo - rho * c_hi) / sigma)
    z_hi = up((x_hi - rho * c_lo) / sigma)
    lo_zlo, _ = Phi_bounds_float(z_lo)
    _, hi_zhi = Phi_bounds_float(z_hi)
    if lo_zlo <= 0.0 or hi_zhi <= 0.0:
        raise RuntimeError("Nonpositive Phi bound in F2_interval_point")
    return down(1.0 / hi_zhi), up(1.0 / lo_zlo)


def F1_box(uL: float, uU: float, vL: float, vU: float) -> Tuple[float, float]:
    # F1 decreases in u and increases in v.
    lo, _ = F1_interval_point(uU, vL)
    _, hi = F1_interval_point(uL, vU)
    return lo, hi


def F2_box(uL: float, uU: float, vL: float, vU: float) -> Tuple[float, float]:
    # F2 increases in u and decreases in v.
    lo, _ = F2_interval_point(uL, vU)
    _, hi = F2_interval_point(uU, vL)
    return lo, hi


def phi_max_interval(a: float, b: float) -> float:
    if a > b:
        a, b = b, a
    if a <= 0.0 <= b:
        return phi_pdf(0.0)
    return max(phi_pdf(a), phi_pdf(b))


def abs_z_phi_max_interval(a: float, b: float) -> float:
    if a > b:
        a, b = b, a
    vals = [abs(a) * phi_pdf(a), abs(b) * phi_pdf(b)]
    if a <= 1.0 <= b:
        vals.append(phi_pdf(1.0))
    if a <= -1.0 <= b:
        vals.append(phi_pdf(1.0))
    return max(vals)


def abs_z2_minus1_max_interval(a: float, b: float) -> float:
    if a > b:
        a, b = b, a
    vals = [abs(a * a - 1.0), abs(b * b - 1.0)]
    if a <= 0.0 <= b:
        vals.append(1.0)
    return max(vals)


def ga_numerator_point_bounds(x: float, c: float, *, n: int = GA_N) -> Tuple[float, float]:
    """Bounds numerator int_{-inf}^a Phi(rho*c+sigma*t) phi(t) dt."""
    mu = rho * c
    a = (x - mu) / sigma
    if a <= GA_T_MIN:
        # Crude but safe: 0 <= numerator <= P(T<=a) <= Phi(a).
        _, den_hi = Phi_bounds_float(a)
        return 0.0, den_hi

    width = (a - GA_T_MIN) / float(n)
    lo_total = 0.0
    hi_total = 0.0
    for i in range(n):
        l = GA_T_MIN + width * i
        r = l + width
        mid = 0.5 * (l + r)
        z_mid = mu + sigma * mid
        zL = mu + sigma * l
        zU = mu + sigma * r
        phi_mid = phi_pdf(mid)
        Phi_z_lo, Phi_z_hi = Phi_bounds_float(z_mid)
        f_lo = Phi_z_lo * phi_mid
        f_hi = Phi_z_hi * phi_mid

        # f(t)=A(t)B(t), A=Phi(mu+sigma t), B=phi(t).
        A_max = Phi_bounds_float(zU)[1]
        A1_abs = sigma * phi_max_interval(zL, zU)
        A2_abs = sigma * sigma * abs_z_phi_max_interval(zL, zU)
        B_max = phi_max_interval(l, r)
        B1_abs = abs_z_phi_max_interval(l, r)
        B2_abs = abs_z2_minus1_max_interval(l, r) * B_max
        m2 = A2_abs * B_max + 2.0 * A1_abs * B1_abs + A_max * B2_abs
        rem = m2 * (width ** 3) / 24.0
        lo_total = down(lo_total + width * f_lo - rem)
        hi_total = up(hi_total + width * f_hi + rem)
    if lo_total < 0.0:
        lo_total = 0.0
    # Lower tail below GA_T_MIN.  Upper bound uses Phi(z)<=1.
    _, tail_hi = Phi_bounds_float(GA_T_MIN)
    hi_total = up(hi_total + tail_hi)
    return lo_total, min(1.0, hi_total)


@lru_cache(maxsize=200_000)
def ga_point_bounds_cached(x_key: float, c_key: float) -> Tuple[float, float]:
    # Keys are already rounded endpoint floats; use them directly.
    x = float(x_key)
    c = float(c_key)
    num_lo, num_hi = ga_numerator_point_bounds(x, c)
    a = (x - rho * c) / sigma
    den_lo, den_hi = Phi_bounds_float(a)
    if den_lo <= 0.0 or den_hi <= 0.0:
        raise RuntimeError("Nonpositive denominator in G_A")
    ratio_lo = max(0.0, down(num_lo / den_hi))
    ratio_hi = min(1.0, up(num_hi / den_lo))
    return down(ell + d * ratio_lo), up(ell + d * ratio_hi)


def GA_box(uL: float, uU: float, vL: float, vU: float) -> Tuple[float, float]:
    """Interval for G_A(u,v)=E[A | Y=q(v), X<=q(u)].

    The conditional normal family is monotone in both the truncation point x
    and the conditional mean rho*c, so the lower corner is (xL,cL) and the
    upper corner is (xU,cU), where x=q(u), c=q(v).
    """
    xL, xU = q_interval_from_tail(uL, uU)
    cL, cU = q_interval_from_tail(vL, vU)
    lo, _ = ga_point_bounds_cached(xL, cL)
    _, hi = ga_point_bounds_cached(xU, cU)
    return lo, hi


def value_box(vL: float, vU: float) -> Tuple[float, float]:
    # e(q(v)) = 10 - 9.9 v, decreasing in v.
    return down(nu - d * vU), up(nu - d * vL)


@dataclass
class Box:
    uL: float
    uU: float
    vL: float
    vU: float
    bFL: float
    bFU: float
    bGL: float
    bGU: float
    cGL: float
    cGU: float
    voutL: float
    voutU: float


@dataclass
class DriftBox:
    duL: float
    duU: float
    dvL: float
    dvU: float
    dbFL: float
    dbFU: float
    dbGL: float
    dbGU: float
    dcGL: float
    dcGU: float
    valL: float
    valU: float


@dataclass
class CertificateSummary:
    status: str
    rho: float
    sigma: float
    step_size: float
    s0: float
    kappa_low: float
    theta_03: float
    total_steps: int
    active_certain_steps: int
    uncertain_steps: int
    post_steps: int
    max_picard_iterations: int
    tau_upper_time: Optional[float]
    final_u: Tuple[float, float]
    final_v: Tuple[float, float]
    final_bF: Tuple[float, float]
    final_bG: Tuple[float, float]
    green_cost_interval: Tuple[float, float]
    outside_value_interval: Tuple[float, float]
    cost_slack_lower: float
    value_margin_lower: float


def drift_on_tube(tube: Box, phase: str) -> DriftBox:
    if phase == "active":
        duL, duU = F1_box(tube.uL, tube.uU, tube.vL, tube.vU)
        dbFL = down(-nu + d * tube.uL)
        dbFU = up(-nu + d * tube.uU)
    elif phase == "uncertain":
        _, duU = F1_box(tube.uL, tube.uU, tube.vL, tube.vU)
        duL = 0.0
        dbFL = min(down(-nu + d * tube.uL), 0.0)
        dbFU = up(-nu + d * tube.uU)
    elif phase == "post":
        duL = duU = 0.0
        dbFL = dbFU = 0.0
    else:
        raise ValueError(f"unknown phase {phase}")

    dvL, dvU = F2_box(tube.uL, tube.uU, tube.vL, tube.vU)
    gaL, gaU = GA_box(tube.uL, tube.uU, tube.vL, tube.vU)
    valL, valU = value_box(tube.vL, tube.vU)
    return DriftBox(
        duL=duL,
        duU=duU,
        dvL=dvL,
        dvU=dvU,
        dbFL=dbFL,
        dbFU=dbFU,
        dbGL=down(-gaU),
        dbGU=up(-gaL),
        dcGL=gaL,
        dcGU=gaU,
        valL=valL,
        valU=valU,
    )


def tube_candidate(B: Box, drift: DriftBox, dt: float, value_overlap: float) -> Box:
    return Box(
        uL=B.uL,
        uU=up(B.uU + dt * drift.duU),
        vL=B.vL,
        vU=up(B.vU + dt * drift.dvU),
        bFL=down(B.bFL + dt * drift.dbFL),
        bFU=B.bFU,
        bGL=down(B.bGL + dt * drift.dbGL),
        bGU=B.bGU,
        cGL=B.cGL,
        cGU=up(B.cGU + dt * drift.dcGU),
        voutL=B.voutL,
        voutU=up(B.voutU + value_overlap * drift.valU),
    )


def contains(T: Box, C: Box) -> bool:
    return (
        T.uL <= C.uL and C.uU <= T.uU and
        T.vL <= C.vL and C.vU <= T.vU and
        T.bFL <= C.bFL and C.bFU <= T.bFU and
        T.bGL <= C.bGL and C.bGU <= T.bGU and
        T.cGL <= C.cGL and C.cGU <= T.cGU and
        T.voutL <= C.voutL and C.voutU <= T.voutU
    )


def initial_tube_guess(B: Box, dt: float, phase: str, value_overlap: float) -> Box:
    # A deliberately generous first tube; Picard iteration tightens by verified
    # inclusion, not by trusting these constants.
    fmax = 8.0
    return Box(
        uL=B.uL,
        uU=up(B.uU + PICARD_INFLATION * dt * fmax),
        vL=B.vL,
        vU=up(B.vU + PICARD_INFLATION * dt * fmax),
        bFL=down(B.bFL - PICARD_INFLATION * dt * nu),
        bFU=B.bFU,
        bGL=down(B.bGL - PICARD_INFLATION * dt * nu),
        bGU=B.bGU,
        cGL=B.cGL,
        cGU=up(B.cGU + PICARD_INFLATION * dt * nu),
        voutL=B.voutL,
        voutU=up(B.voutU + PICARD_INFLATION * value_overlap * nu),
    )


def hull(A: Box, B: Box) -> Box:
    return Box(
        uL=min(A.uL, B.uL), uU=max(A.uU, B.uU),
        vL=min(A.vL, B.vL), vU=max(A.vU, B.vU),
        bFL=min(A.bFL, B.bFL), bFU=max(A.bFU, B.bFU),
        bGL=min(A.bGL, B.bGL), bGU=max(A.bGU, B.bGU),
        cGL=min(A.cGL, B.cGL), cGU=max(A.cGU, B.cGU),
        voutL=min(A.voutL, B.voutL), voutU=max(A.voutU, B.voutU),
    )


def picard_step(B: Box, dt: float, phase: str, value_overlap: float) -> Tuple[Box, DriftBox, int]:
    T = initial_tube_guess(B, dt, phase, value_overlap)
    last_drift = None
    for it in range(1, PICARD_MAX_ITER + 1):
        drift = drift_on_tube(T, phase)
        cand = tube_candidate(B, drift, dt, value_overlap)
        if contains(T, cand):
            # Endpoint enclosure, not the whole tube.
            Bnext = Box(
                uL=down(B.uL + dt * drift.duL),
                uU=up(B.uU + dt * drift.duU),
                vL=down(B.vL + dt * drift.dvL),
                vU=up(B.vU + dt * drift.dvU),
                bFL=down(B.bFL + dt * drift.dbFL),
                bFU=up(B.bFU + dt * drift.dbFU),
                bGL=down(B.bGL + dt * drift.dbGL),
                bGU=up(B.bGU + dt * drift.dbGU),
                cGL=down(B.cGL + dt * drift.dcGL),
                cGU=up(B.cGU + dt * drift.dcGU),
                voutL=down(B.voutL + value_overlap * drift.valL),
                voutU=up(B.voutU + value_overlap * drift.valU),
            )
            return Bnext, drift, it
        T = hull(T, cand)
        last_drift = drift
    raise RuntimeError(f"Picard tube inclusion failed after {PICARD_MAX_ITER} iterations in phase={phase}")


# ----------------------------------------------------------------------
# Initial diagonal bootstrap on [0,s0]
# ----------------------------------------------------------------------
U_trial = 2.0 * s0
_, Fdiag_hi = F1_interval_point(U_trial, U_trial)
bootstrap_upper = up(s0 * Fdiag_hi)
if bootstrap_upper > U_trial:
    raise RuntimeError("Initial diagonal bootstrap failed")

uL0 = down(s0)
uU0 = bootstrap_upper
vL0 = down(s0)
vU0 = bootstrap_upper
# Conservative initialization over [0,s0].
bF_L0 = down(budget - nu * s0)
bF_U0 = up(budget + s0 * (-nu + d * uU0))
bG_L0 = down(budget - nu * s0)
bG_U0 = up(budget - ell * s0)
cG_L0 = down(ell * s0)
cG_U0 = up(nu * s0)
vout_L0 = 0.0
vout_U0 = 0.0

B = Box(
    uL=uL0, uU=uU0, vL=vL0, vU=vU0,
    bFL=bF_L0, bFU=bF_U0,
    bGL=bG_L0, bGU=bG_U0,
    cGL=cG_L0, cGU=cG_U0,
    voutL=vout_L0, voutU=vout_U0,
)

# ----------------------------------------------------------------------
# Picard-enclosed hybrid evolution
# ----------------------------------------------------------------------
s = s0
phase = "active"
tau_upper_time: Optional[float] = None
steps = active_steps = uncertain_steps = post_steps = 0
max_iters = 0

while s < kappa_low - 1e-15:
    dt = min(h, kappa_low - s)
    value_overlap = 0.0
    if s + dt > m:
        value_overlap = dt if s >= m else (s + dt - m)

    steps += 1
    if phase == "active":
        active_steps += 1
    elif phase == "uncertain":
        uncertain_steps += 1
    else:
        post_steps += 1

    old = B
    Bnext, drift, iters = picard_step(B, dt, phase, value_overlap)
    max_iters = max(max_iters, iters)

    if phase == "active" and Bnext.bFL <= 0.0:
        # Some admissible path may have exhausted inside this step.  Keep a
        # hybrid enclosure: u may have frozen as early as the pre-step lower
        # endpoint, while its upper endpoint may have followed the active drift.
        Bnext.uL = old.uL
        Bnext.vL = old.vL
        phase = "uncertain"

    if phase in {"active", "uncertain"} and Bnext.bFU <= 0.0:
        # All admissible paths are exhausted by the end of this step.  From the
        # next step onward u is frozen in its current interval.
        if tau_upper_time is None:
            tau_upper_time = s + dt
        phase = "post"

    B = Bnext
    s = up(s + dt)

status = "PASS" if (B.cGU < budget and B.voutL > theta_03 and tau_upper_time is not None) else "FAIL"
summary = CertificateSummary(
    status=status,
    rho=rho,
    sigma=sigma,
    step_size=h,
    s0=s0,
    kappa_low=kappa_low,
    theta_03=theta_03,
    total_steps=steps,
    active_certain_steps=active_steps,
    uncertain_steps=uncertain_steps,
    post_steps=post_steps,
    max_picard_iterations=max_iters,
    tau_upper_time=tau_upper_time,
    final_u=(B.uL, B.uU),
    final_v=(B.vL, B.vU),
    final_bF=(B.bFL, B.bFU),
    final_bG=(B.bGL, B.bGU),
    green_cost_interval=(B.cGL, B.cGU),
    outside_value_interval=(B.voutL, B.voutU),
    cost_slack_lower=budget - B.cGU,
    value_margin_lower=B.voutL - theta_03,
)

lines = []
lines.append("=" * 72)
lines.append(f"rho = 0.3 PICARD STEP-TUBE {REPORT_KIND} CERTIFICATE")
lines.append("=" * 72)
lines.append(f"rho                         = {rho}")
lines.append(f"sigma                       = {sigma:.15f}")
lines.append(f"m                           = {m}")
lines.append(f"kappa_low                   = {kappa_low}")
lines.append(f"threshold Theta_0.3         = {theta_03:.12f}")
lines.append(f"budget                      = {budget:.15f}")
lines.append("")
lines.append("Validated infrastructure")
lines.append(f"  A-S uniform Phi error     = {PHI_APPROX_ERR}")
lines.append(f"  outward round pad         = {ROUND_PAD}")
lines.append(f"  G_A quadrature cells      = {GA_N}")
lines.append("")
lines.append("Initial diagonal bootstrap")
lines.append(f"  s0                        = {s0}")
lines.append(f"  U_trial                   = {U_trial}")
lines.append(f"  F_diag upper at U_trial   = {Fdiag_hi:.15f}")
lines.append(f"  certified u(s0), v(s0) in = [{uL0:.15e}, {uU0:.15e}]")
lines.append(f"  initial bF enclosure      = [{bF_L0:.15f}, {bF_U0:.15f}]")
lines.append(f"  initial bG enclosure      = [{bG_L0:.15f}, {bG_U0:.15f}]")
lines.append("")
lines.append("Picard step-tube settings")
lines.append(f"  step size h               = {h}")
lines.append(f"  total steps               = {steps}")
lines.append(f"  active-certain steps      = {active_steps}")
lines.append(f"  farmer-uncertain steps    = {uncertain_steps}")
lines.append(f"  post-Farmer steps         = {post_steps}")
lines.append(f"  max Picard iterations     = {max_iters}")
lines.append("")
if tau_upper_time is None:
    lines.append("Farmer budget did not certify exhaustion before kappa_low.")
else:
    lines.append(f"Certified upper Farmer-exhaustion time tau <= {tau_upper_time:.9f}")
lines.append(f"Final u interval            = [{B.uL:.15f}, {B.uU:.15f}]")
lines.append(f"Final v interval            = [{B.vL:.15f}, {B.vU:.15f}]")
lines.append("")
lines.append("Validated accumulators")
lines.append(f"  Green cost C_G interval   = [{B.cGL:.12f}, {B.cGU:.12f}]")
lines.append(f"  cost slack lower          = {summary.cost_slack_lower:.12f}")
lines.append(f"  V_out[m,kappa] interval   = [{B.voutL:.12f}, {B.voutU:.12f}]")
lines.append(f"  value margin lower        = {summary.value_margin_lower:.12f}")
lines.append("")
lines.append(f"STATUS: RHO03 PICARD {REPORT_KIND} CERTIFICATE {status}.")
lines.append("This script verifies B_n + [0,h]F(T_n) subset T_n at every step and")
lines.append("evaluates drift, cost, and value integrands on the verified step tube.")
report = "\n".join(lines) + "\n"
print(report, end="")

root = SCRIPT_PATH.parent
(root / f"{OUTPUT_STEM}.txt").write_text(report, encoding="utf-8")
(root / f"{OUTPUT_STEM}.json").write_text(
    json.dumps(asdict(summary), indent=2, sort_keys=True), encoding="utf-8"
)
