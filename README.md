# Validated Certificates for "Adversarial Procurement in Two-Value Space"

## Overview

This repository contains the validated numerical certificates that
close the two asymptotic dominance theorems in the paper:

- **Theorem 4.1** (Baseline, ρ=0): Margin ≥ 0.168 per plot
- **Theorem 4.2** (Extension, ρ=0.3): Margin ≥ 0.214 per plot

Both theorems are proved by a unified fluid-limit architecture that
tracks Green's actual per-turn ecological value through a
deterministic ODE on the survivor pool state. The baseline certificate
uses exact closed-form formulas with 80-digit interval arithmetic.
The ρ=0.3 certificates use conservative interval-Euler ODE enclosures
with one-sided Φ bounds.

## Quick Start

```bash
pip install mpmath scipy numpy
python certificates/baseline/baseline_fluid_limit_certificate.py
python certificates/rho03/rho03_production_certificate.py
python certificates/rho03/rho03_fluid_limit_certificate.py
```

Each script prints a self-contained report ending with PASS or FAIL.

## Frozen Theorem Objects

Both theorems share the following setup:

| Object | Definition |
|--------|-----------|
| DGP | Marginals U[0.1, 10], Gaussian copula |
| Farmer strategy | Naïve max-A (buys most expensive available) |
| Green strategy ME | Max-E (buys highest ecological value) |
| Green strategy MX | Max-E/A (buys highest ratio available) |
| Budget parity | B_F = B_G = (1/2) Σ A_i |
| Outcome | Purchased Conservation (PC) |
| Prefix share | m = 0.10 |

## Frozen Constants

### Shared

| Constant | Value | Definition |
|----------|-------|-----------|
| ℓ | 0.1 | Lower bound of marginal support |
| ν | 10 | Upper bound of marginal support |
| μ_A | 5.05 | E[A] = (ℓ+ν)/2 |
| Budget/N | 2.525 | μ_A/2 |
| m | 0.10 | Prefix share |
| M̲_E(0.10) | 0.901 | Prefix lower envelope |

### Baseline (ρ=0) — Fluid-Limit Certificate

| Constant | Value | Definition |
|----------|-------|-----------|
| φ | 0.295816 | Benchmark Farmer exhaustion share |
| τ | [0.30471, 0.30472] | Certified bracket from S_F(τ) = 2.525 |
| κ | 1 − τ ≈ 0.695 | Green exhaustion (exact identity) |
| q_F | 7.07142 | Q_A(1−φ), Farmer block cutoff |
| Θ₀ | 2.655129 | E[E·1{A≤q_F}] − 0.901 |
| V_out | ≥ 2.82336 | Closed-form fluid-limit value integral |
| Certified margin | ≥ 0.16823 | V_out − Θ₀ |

### Extension (ρ=0.3) — Fluid-Limit + Cost Certificate

| Constant | Value | Definition |
|----------|-------|-----------|
| ρ | 0.3 | Latent Gaussian correlation |
| σ | √(1−0.09) ≈ 0.9539 | Conditional standard deviation |
| φ | 0.295816 | Benchmark split |
| J | 50 | Number of Farmer-active cost bands |
| κ_low | 0.620816 | Certified Green stopping lower bound |
| x_F | Φ⁻¹(1−φ) ≈ 0.5365 | Value host X-cutoff |
| τ (fluid) | ≤ 0.31291 | Farmer exhaustion in fluid limit |
| c_up^FA | ≤ 1.4723 | Farmer-active band cost total |
| c_up^total | ≤ 2.5148 | Full certified cost (bands + continuation) |
| ∫e(c(s))ds | ≥ 2.5785 | Fluid-limit ME ecological value |
| Θ_{0.3} | ≤ 2.3642 | MX threshold |
| Certified margin | ≥ 0.214 | Fluid-limit integral − Θ_{0.3} |

## Repository Structure

```
certificates/
├── README.md                    (this file)
├── requirements.txt
├── baseline/
│   ├── baseline_fluid_limit_certificate.py   # closed-form certificate
│   ├── baseline_fluid_limit_certificate.txt
│   ├── baseline_fluid_limit_certificate.json
│   ├── baseline_production_certificate.py    # old worst-case proof (Online Appendix C)
│   └── baseline_production_certificate.txt
└── rho03/
    ├── rho03_production_certificate.py       # cost certificate
    ├── rho03_production_certificate.txt
    ├── rho03_production_certificate.csv      # 50-band cost table
    ├── rho03_production_certificate.json
    ├── rho03_fluid_limit_certificate.py      # value certificate
    └── rho03_fluid_limit_certificate.txt
```

## Certificate Descriptions

### baseline_fluid_limit_certificate.py (PRIMARY)

Certifies the baseline (ρ=0) asymptotic dominance theorem using
the exact closed-form fluid-limit trajectory.

**What it verifies:**
1. τ bracket: S_F(0.30471) < 2.525 < S_F(0.30472)
2. Complete allocation: κ = 1 − τ (derived from closed-form budget identity)
3. Value integral: V_out(τ) = 0.1(1−m−τ) + 3.3(1−2m)^{3/2} + 1.65(1−2τ)^{3/2}
4. Monotonicity: dV_out/dτ = −0.1 − 4.95√(1−2τ) < 0
5. Threshold: Θ₀ = 2.65513
6. Final margin: ≥ 0.16823

**Key insight:** Under independence, the fluid-limit ODE has the
exact closed-form solution u(s) = v(s) = 1 − √(1−2s), so the
entire certificate reduces to evaluating explicit algebraic formulas
on a τ-bracket. No ODE solver, no Kurtz/Wormald numerics — just
calculus plus interval arithmetic.

**Arithmetic:** mpmath interval arithmetic with 80-digit precision.

### baseline_production_certificate.py (ALTERNATIVE — Online Appendix C)

Certifies the baseline theorem via the earlier worst-case deletion
architecture, with margin > 0.046. Retained as an alternative
proof for the online appendix.

**Arithmetic:** Python Decimal with directed rounding (80 digits).

### rho03_production_certificate.py

Certifies the COST side of the ρ=0.3 extension. Under the
fluid-limit architecture, this script establishes that Green
remains affordable through κ_low = 0.620816.

**What it verifies:**
1. 50-band Farmer-active costs: c_up^FA ≤ 1.4723
2. Continuation cost: C^GO validates within budget
3. Total cost: c_up^total ≤ 2.5148 < 2.525
4. Therefore: K_N^ME / N ≥ 0.620816 − o_p(1)

**Cost notation:**
- c_up^FA = Farmer-active band total ≤ 1.4723 (used in Ω_φ^N discharge)
- c_up^total = c_up^FA + continuation cost ≤ 2.5148 (used in κ_low certification)

**Key dependency:** Uses the maxA-topY cost monotonicity lemma
to validate that Farmer's interference can only reduce Green's
per-band costs.

**Arithmetic:** One-sided Abramowitz-Stegun Φ bounds (published
uniform error ≤ 7.5×10⁻⁸) in both float and Decimal arithmetic.

### rho03_fluid_limit_certificate.py

Certifies the VALUE side of the ρ=0.3 extension by tracking
Green's actual ecological value through the fluid-limit ODE.

**What it verifies:**
1. The piecewise ODE system for the state (u, v, b_F, b_G) is
   integrated from s ≈ 0 to s = κ_low = 0.620816
2. Farmer exhaustion occurs at τ ≤ 0.31291
3. Green's ecological value integral: ∫e(c(s))ds ≥ 2.5785
4. MX threshold: Θ_{0.3} ≤ 2.3642
5. Margin: ≥ 0.214

**Architecture:** During the active phase [0, τ], both the
Farmer frontier u(s) and the Green cutoff v(s) evolve under a
coupled ODE derived from the exact dynamic survivor law. At Farmer
exhaustion (b_F = 0, transversal crossing), u freezes and only v
continues to evolve. Green's per-turn ecological value is
e(c(s)) = 10 − 9.9·v(s) throughout. The value integral is
accumulated only from m = 0.10 onward.

**Arithmetic:** Conservative interval-Euler enclosure with
monotonicity-based one-sided bounds. All Φ evaluations use the
same Abramowitz-Stegun one-sided bounds as the cost certificate
(published error ≤ 7.5×10⁻⁸). The margin (0.214) is large
relative to the cumulative Φ-pad error (~6×10⁻⁴).

## How the Certificates Combine

### Baseline (ρ=0) — single script

`baseline_fluid_limit_certificate.py` alone certifies the entire
theorem. The exact closed-form trajectory makes a separate cost
certificate unnecessary.

### Extension (ρ=0.3) — two scripts

1. **rho03_production_certificate.py** proves Green can afford
   to keep shopping through κ_low (cost side)
2. **rho03_fluid_limit_certificate.py** proves the ecological
   value Green accumulates exceeds the MX benchmark (value side)

Together: Green buys enough items (cost certificate) and those
items have enough ecological value (fluid-limit certificate)
to dominate MX.

## Cross-References to Paper

| Certificate output | Paper proposition |
|-------------------|-------------------|
| Baseline τ bracket | Proposition (baseline fluid certificate) |
| Baseline margin ≥ 0.168 | Theorem 4.1 |
| c_up^FA ≤ 1.4723 | Cost certificate / Ω_φ^N discharge |
| c_up^total ≤ 2.515 | Cost certificate corollary |
| κ_low = 0.620816 | Cost certificate corollary |
| ∫e(c(s))ds ≥ 2.579 | Fluid-limit value proposition |
| ρ=0.3 margin ≥ 0.214 | Theorem 4.2 |

## Proof Architecture Summary

### Baseline (ρ=0) — Closed-Form Fluid Limit

1. Prefix wedge: +0.219 (ME starts ahead)
2. Exact dynamic survivor law → factored host under independence
3. Symmetric ODE: u(s) = v(s) = 1 − √(1−2s)
4. Closed-form spend: S_F(τ) = 2.525 determines τ
5. Complete allocation: κ = 1 − τ (budget identity)
6. Closed-form value integral with monotonicity
7. **Certified margin: ≥ 0.168**

### Extension (ρ=0.3) — ODE Fluid Limit

1. Prefix wedge: +0.309 (amplified by correlation)
2. Exact dynamic survivor law → closed 4D state (u, v, b_F, b_G)
3. Cost side: maxA-topY lemma → band costs valid → κ_low = 0.621
4. Ω_φ^N discharge: Farmer-active cost 1.473 ≪ budget 2.525
5. Value side: global fluid-limit ODE tracking v(s) = e(c(s))
6. Piecewise ODE: active phase (both evolve) → post-Farmer (u freezes)
7. ε-launch: convergence on [ε, κ_low], value integral from m = 0.10
8. **Certified margin: ≥ 0.214**

## Requirements

```
# requirements.txt
mpmath>=1.3.0
scipy>=1.10.0
numpy>=1.24.0
```

Python ≥ 3.10.

## Reproducing the Certificates

Each script runs independently with no arguments:

```bash
python baseline_fluid_limit_certificate.py
python rho03_production_certificate.py
python rho03_fluid_limit_certificate.py
```

Expected runtime: ~10 seconds for baseline (closed-form),
~5 minutes for ρ=0.3 cost certificate, ~30 seconds for
fluid-limit certificate.

Scripts are deterministic. Identical output expected on any machine
with the specified library versions.

## Zenodo Archive

DOI 10.5281/zenodo.19598799

## Citation

Weinhold, D. and Andersen, L. (2026). "Adversarial Procurement
in Two-Value Space." [Journal TBD].

## License

MIT License.
