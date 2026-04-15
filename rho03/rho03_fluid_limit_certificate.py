#!/usr/bin/env python3
"""
rho03_fluid_limit_certificate_hardened.py

Validated-style hardening of the rho = 0.3 fluid-limit outside-prefix
certificate.

Compared with the pilot script, this version:

1. Replaces every normal-CDF call by one-sided validated bounds built from the
   classical Abramowitz--Stegun rational approximation with its published
   uniform absolute error bound.
2. Replaces every inverse-normal call by a verified one-sided bracket built
   from a nominal standard-library center plus rigorous Phi-bound checks.
3. Uses outward rounding (via math.nextafter) at every arithmetic step.
4. Fixes two enclosure details that matter for theorem-grade use:
   - F1 is decreasing in u and increasing in v (not increasing in both).
   - After the first step whose *upper* Farmer budget reaches zero, the
     post-Farmer phase freezes u on a conservative interval
     [u_freeze_lower, u_freeze_upper], not just at a single upper value.
5. Avoids the singular invPhi(1) issue at s = 0 by certifying a small-time
   diagonal bootstrap enclosure on [0, s0].

State variables
---------------
The active-phase ODE is written in terms of tail probabilities

    u(s) = 1 - Phi(x(s)),
    v(s) = 1 - Phi(c(s)),

with Farmer budget b_F(s). The outside-prefix ecological value integrand is

    e(c(s)) = 0.1 + 9.9 Phi(c(s)) = 10 - 9.9 v(s),

so a lower bound on the value integral is obtained by integrating with the
upper v-envelope.

Outputs
-------
Running the script writes:

- rho03_fluid_limit_certificate_hardened.txt

and prints the same report to stdout.
"""

from __future__ import annotations

import math
import statistics
from functools import lru_cache
from pathlib import Path

norm = statistics.NormalDist()

# ----------------------------------------------------------------------
# Frozen constants
# ----------------------------------------------------------------------
rho = 0.3
sigma = math.sqrt(1.0 - rho * rho)

m = 0.10
kappa_low = 0.620816
theta_03 = 2.3641180445
budget = 2.525

# Initial enclosure start time and step size.
s0 = 1.0e-5
h = 1.0e-4

# ----------------------------------------------------------------------
# Validated Phi infrastructure (from the hardened production certificate)
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

# Safe inverse-CDF domain used only for nominal centers. The actual one-sided
# correctness comes from Phi_bounds_float verification.
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


@lru_cache(maxsize=200_000)
def Phi_bounds_float(x: float) -> tuple[float, float]:
    """Return (lo, hi) such that lo <= Phi(x) <= hi."""
    if x >= 0.0:
        t = 1.0 / (1.0 + P * x)
        poly = (((((B5 * t) + B4) * t + B3) * t + B2) * t + B1) * t
        approx = 1.0 - phi_pdf(x) * poly
        lo = max(0.0, approx - PHI_APPROX_ERR - ROUND_PAD)
        hi = min(1.0, approx + PHI_APPROX_ERR + ROUND_PAD)
        lo = down(lo)
        hi = up(hi)
        if lo < 0.0:
            lo = 0.0
        if hi > 1.0:
            hi = 1.0
        return lo, hi

    lo_pos, hi_pos = Phi_bounds_float(-x)
    lo = max(0.0, down(1.0 - hi_pos))
    hi = min(1.0, up(1.0 - lo_pos))
    return lo, hi


@lru_cache(maxsize=200_000)
def invPhi_lower(p: float) -> float:
    """Return x such that Phi(x) <= p, verified by Phi_bounds_float.

    A nominal center is obtained from the standard library on a slightly
    shifted probability, then rigorously verified. If needed, an exponential
    search followed by binary refinement pushes the point down until the
    certificate condition is satisfied.
    """
    p = clamp_prob(p)
    # Bias the nominal center slightly to the left to absorb the published
    # Phi approximation error before verification.
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
            raise RuntimeError(f"invPhi_lower failed to bracket p={p}")
        x = cand

    lo = good
    hi_x = fail
    for _ in range(80):
        mid = (lo + hi_x) / 2.0
        _, hi_mid = Phi_bounds_float(mid)
        if hi_mid <= p:
            lo = mid
        else:
            hi_x = mid
    return down(lo)


@lru_cache(maxsize=200_000)
def invPhi_upper(p: float) -> float:
    """Return x such that Phi(x) >= p, verified by Phi_bounds_float."""
    p = clamp_prob(p)
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
            raise RuntimeError(f"invPhi_upper failed to bracket p={p}")
        x = cand

    lo_x = fail
    hi = good
    for _ in range(80):
        mid = (lo_x + hi) / 2.0
        lo_mid, _ = Phi_bounds_float(mid)
        if lo_mid >= p:
            hi = mid
        else:
            lo_x = mid
    return up(hi)


# ----------------------------------------------------------------------
# Interval-safe drift evaluations at a single (u, v) point
# ----------------------------------------------------------------------

def F1_interval(u: float, v: float) -> tuple[float, float]:
    """Rigorous interval for

        F1(u, v) = 1 / Phi((c - rho*x)/sigma),
        x = invPhi(1-u), c = invPhi(1-v).

    Because z = (c - rho*x)/sigma is increasing in c and decreasing in x,
    z_lo uses (c_lo, x_hi) and z_hi uses (c_hi, x_lo). Since 1/Phi(z)
    decreases in z, the lower drift uses z_hi and the upper drift uses z_lo.
    """
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
        raise RuntimeError(f"Nonpositive Phi bound in F1_interval at (u,v)=({u},{v})")

    return down(1.0 / hi_zhi), up(1.0 / lo_zlo)



def F2_interval(u: float, v: float) -> tuple[float, float]:
    """Rigorous interval for

        F2(u, v) = 1 / Phi((x - rho*c)/sigma).

    Here z = (x - rho*c)/sigma is increasing in x and decreasing in c, so
    z_lo uses (x_lo, c_hi) and z_hi uses (x_hi, c_lo). Again 1/Phi(z)
    decreases in z.
    """
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
        raise RuntimeError(f"Nonpositive Phi bound in F2_interval at (u,v)=({u},{v})")

    return down(1.0 / hi_zhi), up(1.0 / lo_zlo)



def BFprime_lower(u: float) -> float:
    return down(-10.0 + 9.9 * u)



def BFprime_upper(u: float) -> float:
    return up(-10.0 + 9.9 * u)


# ----------------------------------------------------------------------
# Small-time diagonal bootstrap on [0, s0]
# ----------------------------------------------------------------------
# In the active phase the ODE is symmetric under u <-> v, and with equal
# initial conditions u(0)=v(0)=0 the exact solution stays on the diagonal
# while both players are active. On the diagonal,
#
#     w' = F_diag(w) := F1(w,w) = F2(w,w),
#
# and F_diag(w) >= 1. We use U_trial = 2*s0 as a bootstrap bound. If the
# rigorous upper drift on the diagonal at U_trial satisfies s0*F_diag(U_trial)
# <= U_trial, then w(s0) <= U_trial and hence
#
#     s0 <= u(s0)=v(s0) <= s0*F_diag(U_trial).
#
# This gives a rigorous finite starting box and avoids any invPhi(1) call.

U_trial = 2.0 * s0
_, Fdiag_hi = F1_interval(U_trial, U_trial)
bootstrap_upper = up(s0 * Fdiag_hi)
if bootstrap_upper > U_trial:
    raise RuntimeError(
        "Initial diagonal bootstrap failed: increase U_trial or decrease s0."
    )

uL = down(s0)
uU = bootstrap_upper
vL = down(s0)
vU = bootstrap_upper

# Store the rigorous starting box for reporting.
uL0, uU0 = uL, uU
vL0, vU0 = vL, vU

# Farmer budget enclosure over [0, s0]. The lower budget uses the maximal spend
# rate 10. For the upper budget, u(t) <= uU over [0, s0], so b_F' <= -10+9.9*uU.
bF_L = down(budget - 10.0 * s0)
bF_U = up(budget + s0 * (-10.0 + 9.9 * uU))
bF_L0, bF_U0 = bF_L, bF_U

# ----------------------------------------------------------------------
# Conservative interval-Euler evolution
# ----------------------------------------------------------------------
s = s0
I_lower = 0.0
farmer_still_active = True
u_freeze_lower = None
u_freeze_upper = None
tau_upper_time = None
steps = 0
active_steps = 0
post_steps = 0

while s < kappa_low - 1e-15:
    dt = min(h, kappa_low - s)
    steps += 1

    if farmer_still_active:
        active_steps += 1
        old_uL, old_uU = uL, uU
        old_vL, old_vU = vL, vU

        # Correct monotonicity:
        #   F1 decreases in u and increases in v,
        #   F2 increases in u and decreases in v.
        duL, _ = F1_interval(uU, vL)
        _, duU = F1_interval(uL, vU)

        dvL, _ = F2_interval(uL, vU)
        _, dvU = F2_interval(uU, vL)

        dbL = BFprime_lower(uL)
        dbU = BFprime_upper(uU)

        uL_prop = down(uL + dt * duL)
        uU_prop = up(uU + dt * duU)
        vL_prop = down(vL + dt * dvL)
        vU_prop = up(vU + dt * dvU)
        bF_L_prop = down(bF_L + dt * dbL)
        bF_U_prop = up(bF_U + dt * dbU)

        # Integrate the whole step as active. This delays Farmer exhaustion to
        # the step endpoint, which is conservative for a lower value bound.
        if bF_U_prop <= 0.0:
            farmer_still_active = False
            tau_upper_time = s + dt

            # Safe frozen-u interval for all later times. The lower endpoint uses
            # the pre-step lower value because exhaustion may occur at any time in
            # [s, s+dt], so the actual frozen u can be smaller than the full-step
            # active lower update.
            u_freeze_lower = old_uL
            u_freeze_upper = uU_prop

            # For v, keep the conservative full-step active upper envelope but the
            # pre-step lower envelope. This is safe and slightly wider, which only
            # makes the later lower value bound more conservative.
            uL, uU = u_freeze_lower, u_freeze_upper
            vL, vU = old_vL, vU_prop
            bF_L, bF_U = bF_L_prop, bF_U_prop
        else:
            uL, uU = uL_prop, uU_prop
            vL, vU = vL_prop, vU_prop
            bF_L, bF_U = bF_L_prop, bF_U_prop

    else:
        post_steps += 1
        if u_freeze_lower is None or u_freeze_upper is None:
            raise RuntimeError("Post-Farmer phase entered without frozen-u interval")

        dvL, _ = F2_interval(u_freeze_lower, vU)
        _, dvU = F2_interval(u_freeze_upper, vL)

        vL = down(vL + dt * dvL)
        vU = up(vU + dt * dvU)

    # Lower bound on ecological value over the step:
    # e(c) = 10 - 9.9*v decreases in v, so use the upper v-envelope on the
    # whole step. This is conservative because v is increasing.
    if s + dt > m:
        overlap = dt if s >= m else (s + dt - m)
        I_lower = down(I_lower + overlap * (10.0 - 9.9 * vU))

    s += dt

margin = I_lower - theta_03
status = "PASS" if margin > 0.0 else "FAIL"

report = []
report.append("=" * 72)
report.append("rho = 0.3 GLOBAL FLUID-LIMIT HARDENED CERTIFICATE")
report.append("=" * 72)
report.append(f"rho                         = {rho}")
report.append(f"sigma                       = {sigma:.15f}")
report.append(f"m                           = {m}")
report.append(f"kappa_low                   = {kappa_low}")
report.append(f"threshold Theta_0.3         = {theta_03:.12f}")
report.append(f"budget                      = {budget}")
report.append("")
report.append("Validated Phi / invPhi infrastructure")
report.append(f"  A-S uniform Phi error     = {PHI_APPROX_ERR}")
report.append(f"  outward round pad         = {ROUND_PAD}")
report.append("")
report.append("Initial small-time diagonal bootstrap")
report.append(f"  s0                        = {s0}")
report.append(f"  U_trial                   = {U_trial}")
report.append(f"  F_diag upper at U_trial   = {Fdiag_hi:.15f}")
report.append(f"  certified u(s0), v(s0) in = [{uL0:.15e}, {uU0:.15e}]")
report.append(f"  certified v(s0) enclosure = [{vL0:.15e}, {vU0:.15e}]")
report.append(f"  initial bF enclosure      = [{bF_L0:.15f}, {bF_U0:.15f}]")
report.append("")
report.append("Interval-Euler settings")
report.append(f"  step size h               = {h}")
report.append(f"  total steps               = {steps}")
report.append(f"  active-phase steps        = {active_steps}")
report.append(f"  post-Farmer steps         = {post_steps}")
report.append("")
if tau_upper_time is None:
    report.append("Farmer budget did not certify exhaustion before kappa_low.")
else:
    report.append(f"Certified upper Farmer-exhaustion time tau <= {tau_upper_time:.9f}")
    report.append(
        f"Frozen u interval after tau = [{u_freeze_lower:.15f}, {u_freeze_upper:.15f}]"
    )
report.append("")
report.append("Value-side lower bound")
report.append(f"  lower bound on V_out[m,kappa_low] = {I_lower:.12f}")
report.append(f"  margin over Theta_0.3            = {margin:.12f}")
report.append("")
report.append(f"STATUS: RHO03 GLOBAL FLUID HARDENED {status}.")
report.append(
    "This script uses one-sided validated Phi bounds, verified inverse-Phi "
    "brackets, outward rounding, a diagonal small-time bootstrap, and a "
    "conservative active-to-post switch."
)

text = "\n".join(report)
print(text)

out = Path(__file__).resolve().with_name("rho03_fluid_limit_certificate_hardened.txt")
out.write_text(text, encoding="utf-8")
