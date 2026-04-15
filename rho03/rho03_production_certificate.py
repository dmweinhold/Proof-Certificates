#!/usr/bin/env python3
"""Production certificate script for the rho = 0.3 dual-host appendix object.

This script combines the validated-style value-side and cost-side pilots into a
single executable certificate for the revised two-host corollary.

What it certifies
-----------------
With frozen theorem inputs

    phi        = 0.295816
    v_low      = 0.325
    V_FA_low   = 1.360955      (taken here as a frozen supporting input)
    threshold  = 2.36411805    (taken here as a frozen supporting input)

it proves the two inequalities required by the revised rho = 0.3 corollary:

1. Value side on the larger host H_val = {A <= Q_A(1-phi)}:

       V_FA_low + Delta M_GO_E,0.3(v_low) > threshold.

2. Cost side on the smaller frontier host
   H_cost = {A <= Q_A(sqrt(1-2 phi))}:

       c_up + C_GO_cost(v_low) < 2.525.

Design notes
------------
- The value side uses Decimal arithmetic and one-sided normal-CDF bounds built
  from the classical Abramowitz--Stegun rational approximation with a published
  uniform absolute error bound.
- The cost side reuses the same one-sided Phi bounds in ordinary double
  precision, which is more than adequate because the explicit Phi error pad is
  orders of magnitude larger than floating roundoff.
- The original floating deterministic script is used only to obtain nominal
  *centers* for band cutoffs and continuation roots. The final inequalities are
  certified by one-sided bracketing and one-sided quadrature, not by trusting
  those floating outputs.

Outputs
-------
Running the script writes three files next to the script:

- rho03_production_certificate.txt   human-readable report
- rho03_production_certificate.csv   certified band table for the cost side
- rho03_production_certificate.json  machine-readable summary

Limitations
-----------
This script treats V_FA_low = 1.360955 and threshold = 2.36411805 as frozen
inputs from the appendix object. Their independent certification is a separate
supporting step.
"""

from __future__ import annotations

import csv
import importlib.machinery
import importlib.util
import json
import math
import sys
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


# =============================================================================
# Locate and load the original floating deterministic script.
# =============================================================================

ROOT = Path(__file__).resolve().parent
SOURCE_CANDIDATES = [
    ROOT / "ORIGINAL_rho03_deterministic_band_reduction.txt",
    ROOT / "rho03_deterministic_band_reduction.txt",
]


def find_source_script() -> Path:
    for candidate in SOURCE_CANDIDATES:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Could not locate the original deterministic band-reduction script. "
        "Expected one of: " + ", ".join(str(p.name) for p in SOURCE_CANDIDATES)
    )


SRC = find_source_script()


def load_band_module():
    name = "rho03_band_src_production_certificate"
    loader = importlib.machinery.SourceFileLoader(name, str(SRC))
    spec = importlib.util.spec_from_loader(name, loader)
    if spec is None:
        raise RuntimeError("Could not create import spec for the source script")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


mod = load_band_module()
rows = mod.build_bands()


# =============================================================================
# Shared frozen theorem constants.
# =============================================================================

PHI_F = float(mod.PHI)
DELTA_F = float(mod.DELTA)
RHO_F = float(mod.RHO)
SIGMA_F = float(mod.SIGMA)
TOTAL_BUDGET_F = float(mod.TOTAL_BUDGET)
V_LOW_F = 0.325
THRESHOLD_F = 2.36411805
V_FA_LOW_F = 1.360955

# Decimal copies for the value-side layer.
getcontext().prec = 80
D = Decimal
PHI_D = D("0.295816")
V_LOW_D = D("0.325")
RHO_D = D("0.3")
SIGMA_D = (D("1") - RHO_D * RHO_D).sqrt()
THRESHOLD_D = D("2.36411805")
V_FA_LOW_D = D("1.360955")
INCREMENT_TARGET_D = THRESHOLD_D - V_FA_LOW_D

ONE_D = D("1")
ZERO_D = D("0")
TWO_D = D("2")
TWENTY4_D = D("24")


# =============================================================================
# A-S Phi approximation constants, both float and Decimal.
# =============================================================================

P_F = 0.2316419
B1_F = 0.319381530
B2_F = -0.356563782
B3_F = 1.781477937
B4_F = -1.821255978
B5_F = 1.330274429
PHI_APPROX_ERR_F = 7.5e-8
ROUND_PAD_F = 1.0e-12

P_D = D("0.2316419")
B1_D = D("0.319381530")
B2_D = D("-0.356563782")
B3_D = D("1.781477937")
B4_D = D("-1.821255978")
B5_D = D("1.330274429")
PHI_APPROX_ERR_D = D("7.5e-8")
ROUND_PAD_D = D("1e-30")

PI_D = D(
    "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899"
)
SQRT_2PI_D = (TWO_D * PI_D).sqrt()
SQRT_2PI_F = math.sqrt(2.0 * math.pi)
BETA_F = RHO_F / SIGMA_F
BETA_D = RHO_D / SIGMA_D


# =============================================================================
# Shared utility helpers.
# =============================================================================


def phi_pdf_float(x: float) -> float:
    return math.exp(-0.5 * x * x) / SQRT_2PI_F


def Phi_bounds_float(x: float) -> Tuple[float, float]:
    """One-sided bounds for the standard normal CDF in float arithmetic."""
    if x >= 0.0:
        t = 1.0 / (1.0 + P_F * x)
        poly = (((((B5_F * t) + B4_F) * t + B3_F) * t + B2_F) * t + B1_F) * t
        approx = 1.0 - phi_pdf_float(x) * poly
        lo = max(0.0, approx - PHI_APPROX_ERR_F - ROUND_PAD_F)
        hi = min(1.0, approx + PHI_APPROX_ERR_F + ROUND_PAD_F)
    else:
        lo_pos, hi_pos = Phi_bounds_float(-x)
        lo = max(0.0, 1.0 - hi_pos)
        hi = min(1.0, 1.0 - lo_pos)
    return lo, hi


def phi_pdf_decimal(x: Decimal) -> Decimal:
    return (-(x * x) / TWO_D).exp() / SQRT_2PI_D


def Phi_bounds_decimal(x: Decimal) -> Tuple[Decimal, Decimal]:
    """One-sided bounds for the standard normal CDF in Decimal arithmetic."""
    if x >= 0:
        t = ONE_D / (ONE_D + P_D * x)
        poly = (((((B5_D * t) + B4_D) * t + B3_D) * t + B2_D) * t + B1_D) * t
        approx = ONE_D - phi_pdf_decimal(x) * poly
        lo = approx - PHI_APPROX_ERR_D - ROUND_PAD_D
        hi = approx + PHI_APPROX_ERR_D + ROUND_PAD_D
    else:
        lo_pos, hi_pos = Phi_bounds_decimal(-x)
        lo = ONE_D - hi_pos
        hi = ONE_D - lo_pos

    if lo < ZERO_D:
        lo = ZERO_D
    if hi > ONE_D:
        hi = ONE_D
    return lo, hi


# =============================================================================
# Value-side production certificate.
# =============================================================================

# Floating centers used only to seed rigorous brackets.
X_VALUE_CENTER_D = D("0.536472553270144")
D_PHI_CENTER_D = D("0.051103921604505")
D_PHI_PLUS_V_CENTER_D = D("-1.3035913758521738")

X_STEP_D = D("1e-7")
ROOT_STEP_D = D("5e-6")
ROOT_YMAX_D = D("8")
ROOT_N_D = 4000
VALUE_N_D = 2000


@dataclass
class ValueCertificateResult:
    x_lo: str
    x_hi: str
    d_phi_lo: str
    d_phi_hi: str
    d_phi_plus_v_lo: str
    d_phi_plus_v_hi: str
    delta_lower: str
    delta_upper: str
    increment_target: str
    increment_margin: str
    total_lower: str
    total_upper: str
    threshold: str
    total_margin: str
    status: str



def max_abs_decimal(a: Decimal, b: Decimal) -> Decimal:
    aa = abs(a)
    bb = abs(b)
    return aa if aa > bb else bb



def phi_max_interval_decimal(a: Decimal, b: Decimal) -> Decimal:
    if a > b:
        a, b = b, a
    if a <= ZERO_D <= b:
        return phi_pdf_decimal(ZERO_D)
    pa = phi_pdf_decimal(a)
    pb = phi_pdf_decimal(b)
    return pa if pa > pb else pb



def abs_z_phi_max_interval_decimal(a: Decimal, b: Decimal) -> Decimal:
    if a > b:
        a, b = b, a
    candidates = [abs(a) * phi_pdf_decimal(a), abs(b) * phi_pdf_decimal(b)]
    if a <= ONE_D <= b:
        candidates.append(phi_pdf_decimal(ONE_D))
    if a <= -ONE_D <= b:
        candidates.append(phi_pdf_decimal(ONE_D))
    return max(candidates)



def abs_z2_minus1_max_interval_decimal(a: Decimal, b: Decimal) -> Decimal:
    if a > b:
        a, b = b, a
    vals = [abs(a * a - ONE_D), abs(b * b - ONE_D)]
    if a <= ZERO_D <= b:
        vals.append(ONE_D)
    return max(vals)



def h_value_decimal(y: Decimal, x_value: Decimal) -> Decimal:
    return (x_value - RHO_D * y) / SIGMA_D



def mass_integrand_midpoint_bounds_value(y: Decimal, x_lo: Decimal, x_hi: Decimal) -> Tuple[Decimal, Decimal]:
    py = phi_pdf_decimal(y)
    lo_h, _ = Phi_bounds_decimal(h_value_decimal(y, x_lo))
    _, hi_h = Phi_bounds_decimal(h_value_decimal(y, x_hi))
    return py * lo_h, py * hi_h



def mass_cell_second_derivative_bound_value(a: Decimal, b: Decimal, x_upper: Decimal) -> Decimal:
    if a > b:
        a, b = b, a

    phi_y_max = phi_max_interval_decimal(a, b)
    y_abs_max = max_abs_decimal(a, b)
    y2m1_max = abs_z2_minus1_max_interval_decimal(a, b)

    h_left = h_value_decimal(a, x_upper)
    h_right = h_value_decimal(b, x_upper)

    _, a_max = Phi_bounds_decimal(h_left)
    a_prime_abs = BETA_D * phi_max_interval_decimal(h_right, h_left)
    a_second_abs = (BETA_D * BETA_D) * abs_z_phi_max_interval_decimal(h_right, h_left)

    return phi_y_max * (
        a_second_abs + TWO_D * y_abs_max * a_prime_abs + a_max * y2m1_max
    )



def midpoint_mass_bounds_value(
    d: Decimal,
    x_lo: Decimal,
    x_hi: Decimal,
    *,
    ymax: Decimal = ROOT_YMAX_D,
    n: int = ROOT_N_D,
) -> Tuple[Decimal, Decimal]:
    """Bounds for m(d) = int_d^inf phi(y) Phi((x-rho y)/sigma) dy."""
    w = (ymax - d) / D(n)
    lo_total = ZERO_D
    hi_total = ZERO_D
    x_upper = x_hi if x_hi > x_lo else x_lo

    for i in range(n):
        left = d + w * D(i)
        right = left + w
        mid = (left + right) / TWO_D

        f_lo_mid, f_hi_mid = mass_integrand_midpoint_bounds_value(mid, x_lo, x_hi)
        m2 = mass_cell_second_derivative_bound_value(left, right, x_upper)
        rem = m2 * (w ** 3) / TWENTY4_D

        lo_total += w * f_lo_mid - rem
        hi_total += w * f_hi_mid + rem

    Phi_ymax_lo, _ = Phi_bounds_decimal(ymax)
    hi_total += ONE_D - Phi_ymax_lo
    if lo_total < ZERO_D:
        lo_total = ZERO_D
    return lo_total, hi_total



def value_integrand_midpoint_bounds(y: Decimal, x_lo: Decimal, x_hi: Decimal) -> Tuple[Decimal, Decimal]:
    py = phi_pdf_decimal(y)

    Ay_lo, Ay_hi = Phi_bounds_decimal(y)
    A_lo = D("0.1") + D("9.9") * Ay_lo
    A_hi = D("0.1") + D("9.9") * Ay_hi

    h_lo, _ = Phi_bounds_decimal(h_value_decimal(y, x_lo))
    _, h_hi = Phi_bounds_decimal(h_value_decimal(y, x_hi))

    return A_lo * py * h_lo, A_hi * py * h_hi



def value_cell_second_derivative_bound(a: Decimal, b: Decimal, x_upper: Decimal) -> Decimal:
    if a > b:
        a, b = b, a

    phi_y_max = phi_max_interval_decimal(a, b)
    y_abs_max = max_abs_decimal(a, b)
    y2m1_max = abs_z2_minus1_max_interval_decimal(a, b)

    _, Phi_b_hi = Phi_bounds_decimal(b)
    A_max = D("0.1") + D("9.9") * Phi_b_hi
    A_prime_abs = D("9.9") * phi_y_max
    A_second_abs = D("9.9") * abs_z_phi_max_interval_decimal(a, b)

    h_left = h_value_decimal(a, x_upper)
    h_right = h_value_decimal(b, x_upper)
    _, small_a_max = Phi_bounds_decimal(h_left)
    small_a_prime_abs = BETA_D * phi_max_interval_decimal(h_right, h_left)
    small_a_second_abs = (BETA_D * BETA_D) * abs_z_phi_max_interval_decimal(h_right, h_left)

    B_max = A_max * small_a_max
    B_prime_abs = A_prime_abs * small_a_max + A_max * small_a_prime_abs
    B_second_abs = (
        A_second_abs * small_a_max
        + TWO_D * A_prime_abs * small_a_prime_abs
        + A_max * small_a_second_abs
    )

    return phi_y_max * (
        B_second_abs + TWO_D * y_abs_max * B_prime_abs + B_max * y2m1_max
    )



def midpoint_value_bounds(left: Decimal, right: Decimal, x_lo: Decimal, x_hi: Decimal, *, n: int = VALUE_N_D) -> Tuple[Decimal, Decimal]:
    if right <= left:
        return ZERO_D, ZERO_D

    w = (right - left) / D(n)
    lo_total = ZERO_D
    hi_total = ZERO_D
    x_upper = x_hi if x_hi > x_lo else x_lo

    for i in range(n):
        a = left + w * D(i)
        b = a + w
        mid = (a + b) / TWO_D

        f_lo_mid, f_hi_mid = value_integrand_midpoint_bounds(mid, x_lo, x_hi)
        m2 = value_cell_second_derivative_bound(a, b, x_upper)
        rem = m2 * (w ** 3) / TWENTY4_D

        lo_total += w * f_lo_mid - rem
        hi_total += w * f_hi_mid + rem

    if lo_total < ZERO_D:
        lo_total = ZERO_D
    return lo_total, hi_total



def bracket_x_value(center: Decimal = X_VALUE_CENTER_D, step: Decimal = X_STEP_D) -> Tuple[Decimal, Decimal]:
    target = ONE_D - PHI_D

    x_lo = center
    while True:
        _, hi = Phi_bounds_decimal(x_lo)
        if hi <= target:
            break
        x_lo -= step

    x_hi = center
    while True:
        lo, _ = Phi_bounds_decimal(x_hi)
        if lo >= target:
            break
        x_hi += step

    return x_lo, x_hi



def bracket_mass_root_value(
    target_mass: Decimal,
    center: Decimal,
    x_lo: Decimal,
    x_hi: Decimal,
    *,
    step: Decimal = ROOT_STEP_D,
    ymax: Decimal = ROOT_YMAX_D,
    n: int = ROOT_N_D,
) -> Tuple[Decimal, Decimal]:
    d_lo = center
    while True:
        lo_val, _ = midpoint_mass_bounds_value(d_lo, x_lo, x_hi, ymax=ymax, n=n)
        if lo_val >= target_mass:
            break
        d_lo -= step

    d_hi = center
    while True:
        _, hi_val = midpoint_mass_bounds_value(d_hi, x_lo, x_hi, ymax=ymax, n=n)
        if hi_val <= target_mass:
            break
        d_hi += step

    return d_lo, d_hi



def run_value_certificate() -> Tuple[ValueCertificateResult, str]:
    x_lo, x_hi = bracket_x_value()
    d0_lo, d0_hi = bracket_mass_root_value(PHI_D, D_PHI_CENTER_D, x_lo, x_hi)
    d1_lo, d1_hi = bracket_mass_root_value(PHI_D + V_LOW_D, D_PHI_PLUS_V_CENTER_D, x_lo, x_hi)

    delta_lo, delta_hi = midpoint_value_bounds(d1_hi, d0_lo, x_lo, x_hi)
    total_lo = V_FA_LOW_D + delta_lo
    total_hi = V_FA_LOW_D + delta_hi

    result = ValueCertificateResult(
        x_lo=str(x_lo),
        x_hi=str(x_hi),
        d_phi_lo=str(d0_lo),
        d_phi_hi=str(d0_hi),
        d_phi_plus_v_lo=str(d1_lo),
        d_phi_plus_v_hi=str(d1_hi),
        delta_lower=str(delta_lo),
        delta_upper=str(delta_hi),
        increment_target=str(INCREMENT_TARGET_D),
        increment_margin=str(delta_lo - INCREMENT_TARGET_D),
        total_lower=str(total_lo),
        total_upper=str(total_hi),
        threshold=str(THRESHOLD_D),
        total_margin=str(total_lo - THRESHOLD_D),
        status="PASS" if delta_lo > INCREMENT_TARGET_D else "FAIL",
    )

    lines: List[str] = []
    lines.append("=" * 72)
    lines.append("VALUE-SIDE CERTIFICATE")
    lines.append("=" * 72)
    lines.append(f"phi                        = {PHI_D}")
    lines.append(f"v_low                      = {V_LOW_D}")
    lines.append(f"rho                        = {RHO_D}")
    lines.append(f"sigma                      = {SIGMA_D}")
    lines.append("")
    lines.append("Brackets for x_value = Phi^{-1}(1-phi)")
    lines.append(f"  x_lo                     = {x_lo}")
    lines.append(f"  x_hi                     = {x_hi}")
    lines.append("")
    lines.append("Brackets for value-host mass cutoffs")
    lines.append(f"  d_val(phi)_lo            = {d0_lo}")
    lines.append(f"  d_val(phi)_hi            = {d0_hi}")
    lines.append(f"  d_val(phi+v)_lo          = {d1_lo}")
    lines.append(f"  d_val(phi+v)_hi          = {d1_hi}")
    lines.append("")
    lines.append("Mass checks at the certified brackets")
    m0_lo, _ = midpoint_mass_bounds_value(d0_lo, x_lo, x_hi)
    _, m0_hi = midpoint_mass_bounds_value(d0_hi, x_lo, x_hi)
    m1_lo, _ = midpoint_mass_bounds_value(d1_lo, x_lo, x_hi)
    _, m1_hi = midpoint_mass_bounds_value(d1_hi, x_lo, x_hi)
    lines.append(f"  mass_lo(d_val(phi)_lo)   = {m0_lo}")
    lines.append(f"  target phi               = {PHI_D}")
    lines.append(f"  mass_hi(d_val(phi)_hi)   = {m0_hi}")
    lines.append(f"  mass_lo(d_val(phi+v)_lo) = {m1_lo}")
    lines.append(f"  target phi+v             = {PHI_D + V_LOW_D}")
    lines.append(f"  mass_hi(d_val(phi+v)_hi) = {m1_hi}")
    lines.append("")
    lines.append("Value increment bounds")
    lines.append(f"  Delta M lower            = {delta_lo}")
    lines.append(f"  Delta M upper            = {delta_hi}")
    lines.append(f"  increment target         = {INCREMENT_TARGET_D}")
    lines.append(f"  lower margin             = {delta_lo - INCREMENT_TARGET_D}")
    lines.append("")
    lines.append("Total value bounds")
    lines.append(f"  V_FA_low                 = {V_FA_LOW_D}")
    lines.append(f"  total lower              = {total_lo}")
    lines.append(f"  total upper              = {total_hi}")
    lines.append(f"  threshold                = {THRESHOLD_D}")
    lines.append(f"  lower margin             = {total_lo - THRESHOLD_D}")
    lines.append("")
    lines.append(
        "STATUS: value-side certificate PASSES."
        if result.status == "PASS"
        else "STATUS: value-side certificate DOES NOT PASS."
    )
    lines.append(
        "NOTE: V_FA_low and the threshold are treated here as frozen supporting inputs; "
        "their independent certification is separate."
    )

    return result, "\n".join(lines) + "\n"


# =============================================================================
# Cost-side production certificate.
# =============================================================================

X_STEP_F = 1.0e-6
C_STEP_F = 2.0e-6
D_STEP_F = 2.0e-6
T_MIN_F = -8.0
T_MAX_F = 8.0
MASS_N_F = 300
COST_N_BAND1_F = 2000
COST_N_BANDS_F = 1200
COST_N_CONT_F = 1800


@dataclass
class CostBandCertificateRow:
    j: int
    x_lo: float
    x_hi: float
    c_prev_lo: float
    c_prev_hi: float
    c_lo: float
    c_hi: float
    nominal_cost: float
    cost_upper: float
    cost_gap: float


@dataclass
class CostCertificateResult:
    nominal_c_up: float
    nominal_c_last: float
    nominal_x_cont: float
    nominal_d_cont: float
    nominal_cont_cost: float
    max_cutoff_width: float
    c_last_lo: float
    c_last_hi: float
    d_lo: float
    d_hi: float
    min_band_lower_margin: float
    min_band_upper_margin: float
    cont_lower_margin: float
    cont_upper_margin: float
    c_up_upper: float
    c_up_gap: float
    continuation_cost_upper: float
    continuation_gap: float
    total_upper: float
    total_budget: float
    slack_lower: float
    worst_band_index: int
    worst_band_gap: float
    status: str



def phi_max_interval_float(a: float, b: float) -> float:
    if a > b:
        a, b = b, a
    if a <= 0.0 <= b:
        return phi_pdf_float(0.0)
    pa = phi_pdf_float(a)
    pb = phi_pdf_float(b)
    return pa if pa > pb else pb



def abs_z_phi_max_interval_float(a: float, b: float) -> float:
    if a > b:
        a, b = b, a
    vals = [abs(a) * phi_pdf_float(a), abs(b) * phi_pdf_float(b)]
    if a <= 1.0 <= b:
        vals.append(phi_pdf_float(1.0))
    if a <= -1.0 <= b:
        vals.append(phi_pdf_float(1.0))
    return max(vals)



def abs_z2_minus1_max_interval_float(a: float, b: float) -> float:
    if a > b:
        a, b = b, a
    vals = [abs(a * a - 1.0), abs(b * b - 1.0)]
    if a <= 0.0 <= b:
        vals.append(1.0)
    return max(vals)



def h_mass_float(y: float, x: float) -> float:
    return (x - RHO_F * y) / SIGMA_F



def mass_integrand_midpoint_bounds_cost(y: float, x_lo: float, x_hi: float) -> Tuple[float, float]:
    py = phi_pdf_float(y)
    lo_h, _ = Phi_bounds_float(h_mass_float(y, x_lo))
    _, hi_h = Phi_bounds_float(h_mass_float(y, x_hi))
    return py * lo_h, py * hi_h



def mass_cell_second_derivative_bound_cost(a: float, b: float, x_upper: float) -> float:
    if a > b:
        a, b = b, a

    phi_y_max = phi_max_interval_float(a, b)
    y_abs_max = max(abs(a), abs(b))
    y2m1_max = abs_z2_minus1_max_interval_float(a, b)

    h_left = h_mass_float(a, x_upper)
    h_right = h_mass_float(b, x_upper)
    _, a_max = Phi_bounds_float(h_left)
    a_prime_abs = BETA_F * phi_max_interval_float(h_right, h_left)
    a_second_abs = (BETA_F * BETA_F) * abs_z_phi_max_interval_float(h_right, h_left)

    return phi_y_max * (
        a_second_abs + 2.0 * y_abs_max * a_prime_abs + a_max * y2m1_max
    )



def midpoint_mass_bounds_cost(
    left: float,
    right: float,
    x_lo: float,
    x_hi: float,
    *,
    n: int = MASS_N_F,
) -> Tuple[float, float]:
    if right <= left:
        return 0.0, 0.0

    w = (right - left) / float(n)
    lo_total = 0.0
    hi_total = 0.0
    x_upper = x_hi if x_hi > x_lo else x_lo

    for i in range(n):
        a = left + w * i
        b = a + w
        mid = 0.5 * (a + b)

        f_lo_mid, f_hi_mid = mass_integrand_midpoint_bounds_cost(mid, x_lo, x_hi)
        m2 = mass_cell_second_derivative_bound_cost(a, b, x_upper)
        rem = m2 * (w ** 3) / 24.0

        lo_total += w * f_lo_mid - rem
        hi_total += w * f_hi_mid + rem

    if lo_total < 0.0:
        lo_total = 0.0
    return lo_total, hi_total



def cost_integrand_upper(mid: float, a_par: float, b_par: float) -> float:
    py = phi_pdf_float(mid)
    _, Phi_mid_hi = Phi_bounds_float(mid)
    A_hi = 0.1 + 9.9 * Phi_mid_hi

    ua = (a_par - RHO_F * mid) / SIGMA_F
    if math.isinf(b_par):
        D_hi = 1.0 - Phi_bounds_float(ua)[0]
    else:
        ub = (b_par - RHO_F * mid) / SIGMA_F
        D_hi = Phi_bounds_float(ub)[1] - Phi_bounds_float(ua)[0]
        if D_hi < 0.0:
            D_hi = 0.0
        elif D_hi > 1.0:
            D_hi = 1.0

    return A_hi * py * D_hi



def cost_cell_second_derivative_bound(l: float, r: float, a_par: float, b_par: float) -> float:
    phi_y_max = phi_max_interval_float(l, r)
    y_abs_max = max(abs(l), abs(r))
    y2m1_max = abs_z2_minus1_max_interval_float(l, r)

    _, Phi_r_hi = Phi_bounds_float(r)
    A_max = 0.1 + 9.9 * Phi_r_hi
    A_prime_abs = 9.9 * phi_y_max
    A_second_abs = 9.9 * abs_z_phi_max_interval_float(l, r)

    ua_low = (a_par - RHO_F * r) / SIGMA_F
    ua_high = (a_par - RHO_F * l) / SIGMA_F

    if math.isinf(b_par):
        D_max = 1.0 - Phi_bounds_float(ua_low)[0]
        Dprime_abs = BETA_F * phi_max_interval_float(ua_low, ua_high)
        Dsecond_abs = (BETA_F * BETA_F) * abs_z_phi_max_interval_float(ua_low, ua_high)
    else:
        ub_low = (b_par - RHO_F * r) / SIGMA_F
        ub_high = (b_par - RHO_F * l) / SIGMA_F
        u_low = ua_low
        u_high = ub_high
        width_u = (b_par - a_par) / SIGMA_F
        phi_u_max = phi_max_interval_float(u_low, u_high)
        absu_phi_max = abs_z_phi_max_interval_float(u_low, u_high)

        D_max = min(1.0, width_u * phi_u_max)
        Dprime_abs = BETA_F * width_u * absu_phi_max
        Dsecond_abs = (BETA_F * BETA_F) * 2.0 * absu_phi_max

    B_max = A_max * D_max
    B_prime_abs = A_prime_abs * D_max + A_max * Dprime_abs
    B_second_abs = (
        A_second_abs * D_max
        + 2.0 * A_prime_abs * Dprime_abs
        + A_max * Dsecond_abs
    )

    return phi_y_max * (
        B_second_abs + 2.0 * y_abs_max * B_prime_abs + y2m1_max * B_max
    )



def midpoint_cost_upper(left: float, right: float, a_par: float, b_par: float, *, n: int) -> float:
    if right <= left:
        return 0.0

    w = (right - left) / float(n)
    total = 0.0

    for i in range(n):
        l = left + w * i
        r = l + w
        mid = 0.5 * (l + r)

        g_hi = cost_integrand_upper(mid, a_par, b_par)
        m2 = cost_cell_second_derivative_bound(l, r, a_par, b_par)
        rem = m2 * (w ** 3) / 24.0

        total += w * g_hi + rem

    return total



def cost_upper_finite_x(x_upper: float, a_par: float, b_par: float, *, n: int) -> float:
    total = midpoint_cost_upper(T_MIN_F, x_upper, a_par, b_par, n=n)
    tail_lo = 10.0 * Phi_bounds_float(T_MIN_F)[1]
    return total + tail_lo



def cost_upper_infinite_x(a_par: float, *, n: int) -> float:
    total = midpoint_cost_upper(T_MIN_F, T_MAX_F, a_par, math.inf, n=n)
    tail_lo = 10.0 * Phi_bounds_float(T_MIN_F)[1]
    tail_hi = 10.0 * (1.0 - Phi_bounds_float(T_MAX_F)[0])
    return total + tail_lo + tail_hi



def bracket_quantile_float(target_p: float, center: float, *, step: float) -> Tuple[float, float]:
    x_lo = float(center)
    while True:
        _, hi = Phi_bounds_float(x_lo)
        if hi <= target_p:
            break
        x_lo -= step

    x_hi = float(center)
    while True:
        lo, _ = Phi_bounds_float(x_hi)
        if lo >= target_p:
            break
        x_hi += step

    return x_lo, x_hi



def bracket_band_cutoff(center: float, x_lo: float, x_hi: float, c_prev_lo: float, c_prev_hi: float) -> Tuple[float, float]:
    c_lo = float(center)
    while True:
        lo_val, _ = midpoint_mass_bounds_cost(c_lo, c_prev_lo, x_lo, x_lo)
        if lo_val >= DELTA_F:
            break
        c_lo -= C_STEP_F

    c_hi = float(center)
    while True:
        _, hi_val = midpoint_mass_bounds_cost(c_hi, c_prev_hi, x_hi, x_hi)
        if hi_val <= DELTA_F:
            break
        c_hi += C_STEP_F

    return c_lo, c_hi



def bracket_continuation_cutoff(center: float, x_lo: float, x_hi: float, c_last_lo: float, c_last_hi: float) -> Tuple[float, float]:
    d_lo = float(center)
    while True:
        lo_val, _ = midpoint_mass_bounds_cost(d_lo, c_last_lo, x_lo, x_lo)
        if lo_val >= V_LOW_F:
            break
        d_lo -= D_STEP_F

    d_hi = float(center)
    while True:
        _, hi_val = midpoint_mass_bounds_cost(d_hi, c_last_hi, x_hi, x_hi)
        if hi_val <= V_LOW_F:
            break
        d_hi += D_STEP_F

    return d_lo, d_hi



def run_cost_certificate() -> Tuple[CostCertificateResult, List[CostBandCertificateRow], str]:
    nominal_c_up = math.fsum(r.cost for r in rows)
    c_last_nominal = float(rows[-1].c_j)
    x_cont_nominal = float(mod.x_frontier(mod.PHI))
    d_cont_nominal = float(mod.continuation_cutoff(V_LOW_F, c_last_nominal, x_cont_nominal))
    nominal_cont_cost = float(mod.continuation_cost(V_LOW_F, c_last_nominal, x_cont_nominal))

    x_brackets: Dict[int | str, Tuple[float, float]] = {1: (math.inf, math.inf)}
    for j, r in enumerate(rows, start=1):
        if j == 1:
            continue
        target_p = math.sqrt(1.0 - 2.0 * float(r.s_prev))
        x_brackets[j] = bracket_quantile_float(target_p, float(r.x_high), step=X_STEP_F)

    x_brackets["cont"] = bracket_quantile_float(
        math.sqrt(1.0 - 2.0 * PHI_F),
        x_cont_nominal,
        step=X_STEP_F,
    )

    c_brackets: Dict[int, Tuple[float, float]] = {}
    c1_center = float(rows[0].c_j)
    c_brackets[1] = bracket_quantile_float(1.0 - DELTA_F, c1_center, step=1.0e-6)

    for j in range(2, len(rows) + 1):
        x_lo, x_hi = x_brackets[j]
        c_prev_lo, c_prev_hi = c_brackets[j - 1]
        c_brackets[j] = bracket_band_cutoff(float(rows[j - 1].c_j), x_lo, x_hi, c_prev_lo, c_prev_hi)

    d_lo, d_hi = bracket_continuation_cutoff(
        d_cont_nominal,
        x_brackets["cont"][0],
        x_brackets["cont"][1],
        c_brackets[len(rows)][0],
        c_brackets[len(rows)][1],
    )

    table_rows: List[CostBandCertificateRow] = []
    cost_uppers: List[float] = []

    for j, r in enumerate(rows, start=1):
        if j == 1:
            cost_upper = cost_upper_infinite_x(c_brackets[1][0], n=COST_N_BAND1_F)
            x_lo = math.inf
            x_hi = math.inf
            c_prev_lo = math.inf
            c_prev_hi = math.inf
        else:
            x_lo, x_hi = x_brackets[j]
            c_prev_lo, c_prev_hi = c_brackets[j - 1]
            cost_upper = cost_upper_finite_x(x_hi, c_brackets[j][0], c_prev_hi, n=COST_N_BANDS_F)

        cost_uppers.append(cost_upper)
        table_rows.append(
            CostBandCertificateRow(
                j=j,
                x_lo=x_lo,
                x_hi=x_hi,
                c_prev_lo=c_prev_lo,
                c_prev_hi=c_prev_hi,
                c_lo=c_brackets[j][0],
                c_hi=c_brackets[j][1],
                nominal_cost=float(r.cost),
                cost_upper=cost_upper,
                cost_gap=cost_upper - float(r.cost),
            )
        )

    c_up_upper = math.fsum(cost_uppers)
    continuation_cost_upper = cost_upper_finite_x(
        x_brackets["cont"][1],
        d_lo,
        c_brackets[len(rows)][1],
        n=COST_N_CONT_F,
    )
    total_upper = c_up_upper + continuation_cost_upper
    slack_lower = TOTAL_BUDGET_F - total_upper

    band_lo_margins: List[float] = []
    band_hi_margins: List[float] = []

    lo_tail_mass = 1.0 - Phi_bounds_float(c_brackets[1][0])[1]
    hi_tail_mass = 1.0 - Phi_bounds_float(c_brackets[1][1])[0]
    band_lo_margins.append(lo_tail_mass - DELTA_F)
    band_hi_margins.append(DELTA_F - hi_tail_mass)

    for j in range(2, len(rows) + 1):
        x_lo, x_hi = x_brackets[j]
        c_prev_lo, c_prev_hi = c_brackets[j - 1]
        c_lo, c_hi = c_brackets[j]
        lo_val, _ = midpoint_mass_bounds_cost(c_lo, c_prev_lo, x_lo, x_lo)
        _, hi_val = midpoint_mass_bounds_cost(c_hi, c_prev_hi, x_hi, x_hi)
        band_lo_margins.append(lo_val - DELTA_F)
        band_hi_margins.append(DELTA_F - hi_val)

    cont_lo_val, _ = midpoint_mass_bounds_cost(d_lo, c_brackets[len(rows)][0], x_brackets["cont"][0], x_brackets["cont"][0])
    _, cont_hi_val = midpoint_mass_bounds_cost(d_hi, c_brackets[len(rows)][1], x_brackets["cont"][1], x_brackets["cont"][1])

    worst_gap = max(table_rows, key=lambda row: row.cost_gap)

    result = CostCertificateResult(
        nominal_c_up=nominal_c_up,
        nominal_c_last=c_last_nominal,
        nominal_x_cont=x_cont_nominal,
        nominal_d_cont=d_cont_nominal,
        nominal_cont_cost=nominal_cont_cost,
        max_cutoff_width=max(c_hi - c_lo for c_lo, c_hi in c_brackets.values()),
        c_last_lo=c_brackets[len(rows)][0],
        c_last_hi=c_brackets[len(rows)][1],
        d_lo=d_lo,
        d_hi=d_hi,
        min_band_lower_margin=min(band_lo_margins),
        min_band_upper_margin=min(band_hi_margins),
        cont_lower_margin=cont_lo_val - V_LOW_F,
        cont_upper_margin=V_LOW_F - cont_hi_val,
        c_up_upper=c_up_upper,
        c_up_gap=c_up_upper - nominal_c_up,
        continuation_cost_upper=continuation_cost_upper,
        continuation_gap=continuation_cost_upper - nominal_cont_cost,
        total_upper=total_upper,
        total_budget=TOTAL_BUDGET_F,
        slack_lower=slack_lower,
        worst_band_index=worst_gap.j,
        worst_band_gap=worst_gap.cost_gap,
        status="PASS" if total_upper < TOTAL_BUDGET_F else "FAIL",
    )

    lines: List[str] = []
    lines.append("=" * 72)
    lines.append("COST-SIDE CERTIFICATE")
    lines.append("=" * 72)
    lines.append(f"phi                             = {PHI_F:.12f}")
    lines.append(f"delta                           = {DELTA_F:.12f}")
    lines.append(f"v_low                           = {V_LOW_F:.12f}")
    lines.append(f"rho                             = {RHO_F:.12f}")
    lines.append(f"sigma                           = {SIGMA_F:.15f}")
    lines.append("")
    lines.append("Nominal floating centers")
    lines.append(f"  nominal c_up                  = {nominal_c_up:.15f}")
    lines.append(f"  nominal c_J                   = {c_last_nominal:.15f}")
    lines.append(f"  nominal x_cont                = {x_cont_nominal:.15f}")
    lines.append(f"  nominal d_cont(v_low)         = {d_cont_nominal:.15f}")
    lines.append(f"  nominal continuation cost     = {nominal_cont_cost:.15f}")
    lines.append("")
    lines.append("Certified cutoff brackets")
    lines.append(f"  max band cutoff width         = {result.max_cutoff_width:.15e}")
    lines.append(f"  final cutoff c_J_lo           = {result.c_last_lo:.15f}")
    lines.append(f"  final cutoff c_J_hi           = {result.c_last_hi:.15f}")
    lines.append(f"  continuation d_lo             = {result.d_lo:.15f}")
    lines.append(f"  continuation d_hi             = {result.d_hi:.15f}")
    lines.append("")
    lines.append("Mass checks")
    lines.append(f"  min band lower margin         = {result.min_band_lower_margin:.15e}")
    lines.append(f"  min band upper margin         = {result.min_band_upper_margin:.15e}")
    lines.append(f"  cont mass_lo(d_lo) - v_low    = {result.cont_lower_margin:.15e}")
    lines.append(f"  v_low - cont mass_hi(d_hi)    = {result.cont_upper_margin:.15e}")
    lines.append("")
    lines.append("Cost bounds")
    lines.append(f"  c_up upper                    = {result.c_up_upper:.15f}")
    lines.append(f"  c_up over nominal             = {result.c_up_gap:.15f}")
    lines.append(f"  continuation cost upper       = {result.continuation_cost_upper:.15f}")
    lines.append(f"  cont cost over nominal        = {result.continuation_gap:.15f}")
    lines.append(f"  total spend upper             = {result.total_upper:.15f}")
    lines.append(f"  total budget                  = {TOTAL_BUDGET_F:.15f}")
    lines.append(f"  certified slack               = {result.slack_lower:.15f}")
    lines.append("")
    lines.append("Band-level diagnostics")
    lines.append(f"  worst band gap index          = {result.worst_band_index}")
    lines.append(f"  worst band gap                = {result.worst_band_gap:.15f}")
    lines.append(f"  band 1 upper                  = {table_rows[0].cost_upper:.15f}")
    lines.append(f"  band 50 upper                 = {table_rows[-1].cost_upper:.15f}")
    lines.append("")
    lines.append(
        "STATUS: cost-side certificate PASSES."
        if result.status == "PASS"
        else "STATUS: cost-side certificate DOES NOT PASS."
    )

    return result, table_rows, "\n".join(lines) + "\n"


# =============================================================================
# Combined report / writers.
# =============================================================================


@dataclass
class CombinedCertificateSummary:
    status: str
    source_script: str
    phi: str
    v_low: str
    threshold: str
    total_budget: str
    value_status: str
    cost_status: str
    value_total_margin: str
    value_increment_margin: str
    cost_slack: str
    note: str



def write_band_csv(path: Path, table_rows: Iterable[CostBandCertificateRow]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "j",
            "x_lo",
            "x_hi",
            "c_prev_lo",
            "c_prev_hi",
            "c_lo",
            "c_hi",
            "nominal_cost",
            "cost_upper",
            "cost_gap",
        ])
        for row in table_rows:
            writer.writerow([
                row.j,
                "inf" if math.isinf(row.x_lo) else f"{row.x_lo:.15f}",
                "inf" if math.isinf(row.x_hi) else f"{row.x_hi:.15f}",
                "inf" if math.isinf(row.c_prev_lo) else f"{row.c_prev_lo:.15f}",
                "inf" if math.isinf(row.c_prev_hi) else f"{row.c_prev_hi:.15f}",
                f"{row.c_lo:.15f}",
                f"{row.c_hi:.15f}",
                f"{row.nominal_cost:.15f}",
                f"{row.cost_upper:.15f}",
                f"{row.cost_gap:.15f}",
            ])



def write_json(path: Path, summary: CombinedCertificateSummary, value_result: ValueCertificateResult, cost_result: CostCertificateResult) -> None:
    payload = {
        "summary": asdict(summary),
        "value": asdict(value_result),
        "cost": asdict(cost_result),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")



def build_combined_report(value_result: ValueCertificateResult, value_report: str, cost_result: CostCertificateResult, cost_report: str) -> Tuple[CombinedCertificateSummary, str]:
    value_margin = D(value_result.total_margin)
    increment_margin = D(value_result.increment_margin)

    final_status = "PASS" if value_result.status == "PASS" and cost_result.status == "PASS" else "FAIL"
    summary = CombinedCertificateSummary(
        status=final_status,
        source_script=str(SRC.name),
        phi=str(PHI_D),
        v_low=str(V_LOW_D),
        threshold=str(THRESHOLD_D),
        total_budget=f"{TOTAL_BUDGET_F:.15f}",
        value_status=value_result.status,
        cost_status=cost_result.status,
        value_total_margin=str(value_margin),
        value_increment_margin=str(increment_margin),
        cost_slack=f"{cost_result.slack_lower:.15f}",
        note=(
            "This production script certifies the two inequalities in the revised "
            "dual-host corollary while treating V_FA_low and the threshold as frozen "
            "supporting inputs."
        ),
    )

    lines: List[str] = []
    lines.append("=" * 72)
    lines.append("rho = 0.3 PRODUCTION CERTIFICATE")
    lines.append("=" * 72)
    lines.append(f"source centers               = {SRC.name}")
    lines.append(f"phi                          = {PHI_D}")
    lines.append(f"v_low                        = {V_LOW_D}")
    lines.append(f"value host                   = A <= Q_A(1-phi)")
    lines.append(f"cost host                    = A <= Q_A(sqrt(1-2phi))")
    lines.append(f"threshold                    = {THRESHOLD_D}")
    lines.append(f"total budget                 = {TOTAL_BUDGET_F:.15f}")
    lines.append("")
    lines.append("Final theorem-object verdict")
    lines.append(f"  value side status          = {value_result.status}")
    lines.append(f"  cost side status           = {cost_result.status}")
    lines.append(f"  value total margin         = {value_margin}")
    lines.append(f"  value increment margin     = {increment_margin}")
    lines.append(f"  cost slack                 = {cost_result.slack_lower:.15f}")
    lines.append("")
    if final_status == "PASS":
        lines.append("STATUS: BOTH CERTIFICATE INEQUALITIES PASS.")
        lines.append(
            "Within the frozen appendix object, the revised rho = 0.3 corollary "
            "is numerically certified at v_low = 0.325."
        )
    else:
        lines.append("STATUS: THE COMBINED CERTIFICATE DOES NOT PASS.")
    lines.append("")
    lines.append(
        "CAVEAT: This script does not yet independently certify the supporting "
        "constants V_FA_low = 1.360955 and threshold = 2.36411805. It treats them "
        "as frozen corollary inputs."
    )
    lines.append("")
    lines.append(value_report.rstrip())
    lines.append("")
    lines.append(cost_report.rstrip())
    lines.append("")
    return summary, "\n".join(lines) + "\n"



def main() -> None:
    value_result, value_report = run_value_certificate()
    cost_result, table_rows, cost_report = run_cost_certificate()
    summary, report = build_combined_report(value_result, value_report, cost_result, cost_report)

    print(report, end="")

    out_txt = Path(__file__).with_suffix(".txt")
    out_csv = Path(__file__).with_suffix(".csv")
    out_json = Path(__file__).with_suffix(".json")

    out_txt.write_text(report, encoding="utf-8")
    write_band_csv(out_csv, table_rows)
    write_json(out_json, summary, value_result, cost_result)


if __name__ == "__main__":
    main()
