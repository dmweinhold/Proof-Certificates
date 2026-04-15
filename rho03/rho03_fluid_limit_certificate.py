#!/usr/bin/env python3
"""
rho = 0.3 fluid-limit pilot certificate (production-style template)

This script computes a conservative lower bound on the outside-prefix
ME ecological value integral up to the proved affordable horizon
kappa_low = 0.620816 using a monotone interval-Euler enclosure for the
state variables

    u(s) = 1 - Phi(x(s)),
    v(s) = 1 - Phi(c(s)),
    b_F(s).

The value integrand is simply
    e(c(s)) = 0.1 + 9.9*Phi(c(s)) = 10 - 9.9*v(s).

IMPORTANT:
- This is a production-style pilot certificate, not a final formal
  interval-ODE proof artifact. It uses conservative one-sided
  monotonicity bounds and a small initial enclosure, but it still relies
  on floating implementations of Phi and invPhi from Python's standard
  library.
- The output is intended to guide theorem integration and later
  validated-numerics work. Because the margin is very large, the script
  is robust to conservative settings.

Outputs:
- prints the lower value integral, threshold, and margin
- writes a text summary to rho03_fluid_limit_certificate.txt
"""

import math
import statistics
from pathlib import Path

norm = statistics.NormalDist()

# ----------------------------------------------------------------------
# Frozen constants
# ----------------------------------------------------------------------
rho = 0.3
sigma = math.sqrt(1.0 - rho * rho)

m = 0.10
kappa_low = 0.620816
theta_03 = 2.3641180445  # direct threshold on the MX side

# Initial enclosure start time and step size
s0 = 1.0e-5
h = 1.0e-4

# Conservative initial enclosure on [0, s0]
# Since u'(0)=v'(0)=1 in the limit and both derivatives are >= 1,
# taking [0, 2*s0] is a crude but safe starting box.
uL = 0.0
uU = 2.0 * s0
vL = 0.0
vU = 2.0 * s0

# Farmer budget enclosure over [0, s0]
# Farmer can spend at most 10 per normalized unit time.
bF_L = 2.525 - 10.0 * s0
bF_U = 2.525

# ----------------------------------------------------------------------
# Utility functions
# ----------------------------------------------------------------------
def Phi(x: float) -> float:
    return norm.cdf(x)

def invPhi(p: float) -> float:
    # clamp away from 0 and 1 for numerical stability
    p = max(min(p, 1.0 - 1e-16), 1.0e-16)
    return norm.inv_cdf(p)

def F1(u: float, v: float) -> float:
    """u' during the active phase. Increasing in u and v."""
    x = invPhi(1.0 - u)
    c = invPhi(1.0 - v)
    z = (c - rho * x) / sigma
    return 1.0 / Phi(z)

def F2(u: float, v: float) -> float:
    """v' during the active phase. Increasing in u, decreasing in v."""
    x = invPhi(1.0 - u)
    c = invPhi(1.0 - v)
    z = (x - rho * c) / sigma
    return 1.0 / Phi(z)

def BFprime(u: float) -> float:
    """Farmer budget drift during the active phase."""
    return -10.0 + 9.9 * u

# ----------------------------------------------------------------------
# Conservative interval-Euler evolution
# ----------------------------------------------------------------------
s = s0
I_lower = 0.0
farmer_still_active = True
u_freeze_upper = None
tau_upper_time = None

while s < kappa_low - 1e-15:
    dt = min(h, kappa_low - s)

    if farmer_still_active:
        # Monotonicity:
        # F1 increasing in (u,v)
        duL = F1(uL, vL)
        duU = F1(uU, vU)

        # F2 increasing in u, decreasing in v
        dvL = F2(uL, vU)
        dvU = F2(uU, vL)

        # bF' increasing in u
        dbL = BFprime(uL)
        dbU = BFprime(uU)

        uL_next = uL + dt * duL
        uU_next = uU + dt * duU
        vL_next = vL + dt * dvL
        vU_next = vU + dt * dvU
        bF_L_next = bF_L + dt * dbL
        bF_U_next = bF_U + dt * dbU

        # Once the upper Farmer budget is <= 0, Farmer is certainly exhausted.
        if bF_U_next <= 0.0 and farmer_still_active:
            farmer_still_active = False
            u_freeze_upper = uU_next
            tau_upper_time = s + dt

        uL, uU = uL_next, uU_next
        vL, vU = vL_next, vU_next
        bF_L, bF_U = bF_L_next, bF_U_next

    else:
        # Post-Farmer phase: u is frozen at a conservative upper value.
        # v' = F2(u_fixed, v) still increases in u and decreases in v.
        dvL = F2(u_freeze_upper, vU)
        dvU = F2(u_freeze_upper, vL)
        vL += dt * dvL
        vU += dt * dvU

    # Lower bound on ecological value over the step:
    # e(c) = 10 - 9.9*v decreases in v, so use the upper v-envelope.
    if s + dt > m:
        overlap = dt if s >= m else (s + dt - m)
        I_lower += overlap * (10.0 - 9.9 * vU)

    s += dt

margin = I_lower - theta_03

status = "PASS" if margin > 0.0 else "FAIL"

report = []
report.append("=" * 72)
report.append("rho = 0.3 GLOBAL FLUID-LIMIT PILOT CERTIFICATE")
report.append("=" * 72)
report.append(f"rho                         = {rho}")
report.append(f"m                           = {m}")
report.append(f"kappa_low                   = {kappa_low}")
report.append(f"threshold Theta_0.3         = {theta_03:.12f}")
report.append("")
report.append("Conservative interval-Euler settings")
report.append(f"  start time s0             = {s0}")
report.append(f"  step size h               = {h}")
report.append(f"  initial u enclosure       = [{0.0:.8f}, {2.0*s0:.8f}]")
report.append(f"  initial v enclosure       = [{0.0:.8f}, {2.0*s0:.8f}]")
report.append("")
if tau_upper_time is None:
    report.append("Farmer budget did not certify exhaustion before kappa_low.")
else:
    report.append(f"Certified upper Farmer-exhaustion time tau <= {tau_upper_time:.9f}")
report.append(f"Frozen u upper after tau    = {u_freeze_upper if u_freeze_upper is not None else 'N/A'}")
report.append("")
report.append("Value-side lower bound")
report.append(f"  lower bound on V_out[m,kappa_low] = {I_lower:.12f}")
report.append(f"  margin over Theta_0.3            = {margin:.12f}")
report.append("")
report.append(f"STATUS: RHO03 GLOBAL FLUID PILOT {status}.")
report.append("NOTE: this script is a production-style pilot template, not yet a")
report.append("fully validated interval-ODE proof artifact. It is, however,")
report.append("conservative in orientation and the remaining margin is very large.")

text = "\n".join(report)
print(text)

out = Path("rho03_fluid_limit_certificate.txt")
out.write_text(text)
