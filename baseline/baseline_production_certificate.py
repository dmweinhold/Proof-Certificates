
#!/usr/bin/env python3
"""
baseline_production_certificate.py

Validated closed-form baseline certificate for the repaired Bronze theorem.

This script uses only Python's standard library (Decimal arithmetic with
directed rounding).  It certifies four objects:

1. A bracket for \bar{\tau}, the deterministic Farmer exhaustion share:
      0.30471 < \bar{\tau} < 0.30472
   by checking S_F(0.30471) < 2.525 < S_F(0.30472).

2. Monotonicity of \tau -> L_ME(\tau, 1-\tau) on [0.30, 0.31] by proving
   an explicit upper bound on the derivative that is strictly negative.

3. A certified lower bound on
      L_ME(0.30472, 1-0.30472).

4. The exact threshold
      Theta_0 = E[E 1{A <= q_F}] - 0.901 = 2.65512920
   and the resulting certified margin.

The script writes a human-readable .txt report and a small .json summary.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, localcontext, getcontext, ROUND_FLOOR, ROUND_CEILING

# Global precision
getcontext().prec = 120

# Fixed constants (exact decimals)
TAU_LOW = Decimal("0.30471")
TAU_HIGH = Decimal("0.30472")
PHI = Decimal("0.295816")
TARGET_BUDGET_HALF = Decimal("2.525")
MEAN_E = Decimal("5.05")
PREFIX_E_LOWER = Decimal("0.901")

ONE = Decimal(1)
ZERO = Decimal(0)

def d(x: str) -> Decimal:
    return Decimal(x)

def sqrt_dir(x: Decimal, rounding: str, prec: int = 100) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = prec
        ctx.rounding = rounding
        return x.sqrt()

def S_F_upper(s: Decimal) -> Decimal:
    """
    Upper bound on S_F(s) = 0.1 s + 3.3 (1 - (1-2s)^(3/2)).
    Since the term being subtracted is positive, an upper bound is obtained
    by using a lower bound on sqrt(1-2s).
    """
    t = ONE - d("2.0") * s
    sqrt_low = sqrt_dir(t, ROUND_FLOOR)
    with localcontext() as ctx:
        ctx.prec = 100
        ctx.rounding = ROUND_CEILING
        return d("0.1") * s + d("3.3") * (ONE - t * sqrt_low)

def S_F_lower(s: Decimal) -> Decimal:
    """
    Lower bound on S_F(s) using an upper bound on sqrt(1-2s).
    """
    t = ONE - d("2.0") * s
    sqrt_up = sqrt_dir(t, ROUND_CEILING)
    with localcontext() as ctx:
        ctx.prec = 100
        ctx.rounding = ROUND_FLOOR
        return d("0.1") * s + d("3.3") * (ONE - t * sqrt_up)

def V_FA_lower(tau: Decimal) -> Decimal:
    """
    Lower bound on
      V_FA(tau) = int_{0.10}^{tau} [0.1 + 9.9 sqrt(1-2s)] ds
                = 0.1 (tau-0.1) + 3.3[(0.8)^(3/2) - (1-2tau)^(3/2)].
    """
    a = d("0.8")
    b = ONE - d("2.0") * tau
    a_sqrt_low = sqrt_dir(a, ROUND_FLOOR)
    b_sqrt_up = sqrt_dir(b, ROUND_CEILING)
    with localcontext() as ctx:
        ctx.prec = 100
        ctx.rounding = ROUND_FLOOR
        return d("0.1") * (tau - d("0.1")) + d("3.3") * (a * a_sqrt_low - b * b_sqrt_up)

def V_GO_lower(tau: Decimal) -> Decimal:
    r"""
    Lower bound on
      V_GO(tau,1-tau)
        = (1-\tau) \int_0^{(1-2\tau)/(1-\tau)} (0.1 + 9.9 u) du
        = 0.1(1-2\tau) + 4.95 (1-2\tau)^2/(1-\tau).
    """
    om2 = ONE - d("2.0") * tau
    with localcontext() as ctx:
        ctx.prec = 100
        ctx.rounding = ROUND_FLOOR
        return d("0.1") * om2 + d("4.95") * (om2 * om2) / (ONE - tau)

def L_ME_lower(tau: Decimal) -> Decimal:
    return V_FA_lower(tau) + V_GO_lower(tau)

def theta0_exact() -> Decimal:
    """
    Exact threshold:
      Theta_0 = E[E 1{A <= q_F}] - 0.901
              = 5.05 (1-phi) - 0.901
    because E is independent of A at rho=0.
    """
    return MEAN_E * (ONE - PHI) - PREFIX_E_LOWER

def derivative_upper_on_interval() -> Decimal:
    r"""
    Certified upper bound on d/dtau L_ME(tau,1-tau) for tau in [0.30,0.31].

    d/dtau V_FA = 0.1 + 9.9 sqrt(1-2tau)

    d/dtau V_GO = - (400 tau^2 - 800 tau + 301) / [20 (1-tau)^2]

    On [0.30,0.31]:
      sqrt(1-2tau) <= sqrt(0.4)
      numerator >= 400*(0.31)^2 - 800*(0.31) + 301 = 91.44
      denominator <= 20*(0.70)^2 = 9.8
    so
      d/dtau V_GO <= -91.44/9.8.
    """
    sqrt04_up = sqrt_dir(d("0.4"), ROUND_CEILING)
    first_upper = d("0.1") + d("9.9") * sqrt04_up
    second_lower = d("91.44") / d("9.8")
    with localcontext() as ctx:
        ctx.prec = 100
        ctx.rounding = ROUND_CEILING
        return first_upper - second_lower

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", default="baseline_production_certificate.txt")
    parser.add_argument("--json", default="baseline_production_certificate.json")
    args = parser.parse_args()

    sf_low_upper = S_F_upper(TAU_LOW)
    sf_high_lower = S_F_lower(TAU_HIGH)

    tau_lower_pass = sf_low_upper < TARGET_BUDGET_HALF
    tau_upper_pass = sf_high_lower > TARGET_BUDGET_HALF

    deriv_upper = derivative_upper_on_interval()
    monotonicity_pass = deriv_upper < ZERO

    l_lower = L_ME_lower(TAU_HIGH)
    theta0 = theta0_exact()
    margin_lower = l_lower - theta0
    theorem_pass = tau_lower_pass and tau_upper_pass and monotonicity_pass and (margin_lower > d("0.04620"))

    # Human-readable report
    lines = []
    lines.append("=" * 72)
    lines.append("BASELINE REPAIRED BRONZE THEOREM CERTIFICATE")
    lines.append("=" * 72)
    lines.append(f"phi                           = {PHI}")
    lines.append(f"target half-budget            = {TARGET_BUDGET_HALF}")
    lines.append(f"tau bracket                   = ({TAU_LOW}, {TAU_HIGH})")
    lines.append("")
    lines.append("Step 1. Certified bracket for tau-bar")
    lines.append(f"  S_F({TAU_LOW}) upper        = {sf_low_upper}")
    lines.append(f"  S_F({TAU_HIGH}) lower       = {sf_high_lower}")
    lines.append(f"  lower endpoint check        = {'PASS' if tau_lower_pass else 'FAIL'}")
    lines.append(f"  upper endpoint check        = {'PASS' if tau_upper_pass else 'FAIL'}")
    lines.append("")
    lines.append("Step 2. Certified monotonicity on [0.30, 0.31]")
    lines.append("  derivative formula          = 0.1 + 9.9*sqrt(1-2*tau) "
                 "- (400*tau^2 - 800*tau + 301)/(20*(1-tau)^2)")
    lines.append(f"  certified derivative upper  = {deriv_upper}")
    lines.append(f"  monotonicity status         = {'PASS' if monotonicity_pass else 'FAIL'}")
    lines.append("")
    lines.append("Step 3. Baseline direct reduction")
    lines.append(f"  L_ME({TAU_HIGH}, 1-{TAU_HIGH}) lower = {l_lower}")
    lines.append(f"  Theta_0 exact               = {theta0}")
    lines.append(f"  value margin lower          = {margin_lower}")
    lines.append("")
    lines.append("Certified theorem inequality")
    lines.append(f"  L_ME({TAU_HIGH}, 1-{TAU_HIGH}) > 2.7013307  : {'PASS' if l_lower > d('2.7013307') else 'FAIL'}")
    lines.append(f"  Theta_0 = 2.65512920        : {'PASS' if theta0 == Decimal('2.65512920') else 'FAIL'}")
    lines.append(f"  margin > 0.04620            : {'PASS' if margin_lower > d('0.04620') else 'FAIL'}")
    lines.append("")
    if theorem_pass:
        lines.append("STATUS: BASELINE CERTIFICATE PASSES.")
    else:
        lines.append("STATUS: BASELINE CERTIFICATE FAILS.")
    report = "\n".join(lines) + "\n"

    with open(args.report, "w", encoding="utf-8") as f:
        f.write(report)

    summary = {
        "phi": str(PHI),
        "tau_low": str(TAU_LOW),
        "tau_high": str(TAU_HIGH),
        "S_F_tau_low_upper": str(sf_low_upper),
        "S_F_tau_high_lower": str(sf_high_lower),
        "derivative_upper_on_[0.30,0.31]": str(deriv_upper),
        "L_ME_lower_at_tau_high": str(l_lower),
        "Theta_0_exact": str(theta0),
        "value_margin_lower": str(margin_lower),
        "status": "PASS" if theorem_pass else "FAIL",
    }

    with open(args.json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(report, end="")

if __name__ == "__main__":
    main()
