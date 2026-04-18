# Validated Certificates for "Adversarial Procurement in Two-Value Space"

## Overview

This repository contains the validated numerical certificates that
close the two asymptotic dominance theorems in the paper:

- **Theorem 4.1** (Baseline, ρ=0): Asymptotic gap = 0.216135 (exact identity)
- **Theorem 4.2** (Extension, ρ=0.3): Margin ≥ 0.214 per plot

The baseline theorem is proved via an exact E↔A symmetry of the
fluid-limit dynamics under independence, which yields the asymptotic
gap μ_E·φ − S_G(τ) in closed form. The ρ=0.3 extension uses the
general fluid-limit architecture that tracks Green's actual per-turn
ecological value through a deterministic ODE on the survivor pool
state, closed by validated interval-Euler integration.

## Quick Start

```bash
pip install -r requirements.txt
python certificates/baseline/baseline_closed_form_certificate.py
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
| μ_E = μ_A | 5.05 | Marginal means (identical) |
| Budget/N | 2.525 | μ_A/2 |
| m | 0.10 | Prefix share |

### Baseline (ρ=0) — E↔A Symmetry Identity

| Constant | Value | Definition |
|----------|-------|-----------|
| φ | 0.295816024604... | MX Farmer-exhaustion share; root of 9.9φ² − 20φ + 5.05 = 0 |
| τ | 0.304718564833... | ME Farmer-exhaustion share; root of S_F(τ) = 2.525 |
| μ_E·φ | 1.493871 | MX-side total ecological value |
| S_G(τ) | 1.277736 | Farmer's ecological capture under ME (= Green's A-spend, by symmetry) |
| **Certified gap** | **[0.216134, 0.216136]** | **μ_E·φ − S_G(τ) (exact identity)** |

### Extension (ρ=0.3) — Fluid-Limit + Cost Certificate

| Constant | Value | Definition |
|----------|-------|-----------|
| ρ | 0.3 | Latent Gaussian correlation |
| σ | √(1−0.09) ≈ 0.9539 | Conditional standard deviation |
| φ | 0.295816 | Benchmark split |
| J | 50 | Number of Farmer-active cost bands |
| κ_low | 0.620816 | Certified Green affordability horizon |
| τ (fluid) | ≤ 0.31291 | Farmer exhaustion in fluid limit |
| c_up^FA | ≤ 1.4723 | Farmer-active band cost total |
| c_up^total | ≤ 2.5148 | Full certified cost (bands + continuation) |
| ∫e(c(s))ds | ≥ 2.5785 | Fluid-limit ME ecological value |
| Θ_{0.3} | ≤ 2.3642 | MX threshold |
| **Certified margin** | **≥ 0.214** | **Fluid-limit integral − Θ_{0.3}** |

## Repository Structure

certificates/
├── README.md                    (this file)
├── requirements.txt
├── baseline/
│   ├── baseline_closed_form_certificate.py          # E↔A symmetry identity
│   ├── baseline_closed_form_certificate.txt
│   └── baseline_closed_form_certificate_summary.json
└── rho03/
├── rho03_production_certificate.py              # cost certificate
├── rho03_production_certificate.txt
├── rho03_production_certificate.csv             # 50-band cost table
├── rho03_production_certificate.json
├── rho03_fluid_limit_certificate.py             # value certificate
└── rho03_fluid_limit_certificate.txt

## Certificate Descriptions

### baseline_closed_form_certificate.py (PRIMARY for ρ=0)

Certifies the baseline asymptotic dominance theorem using the exact
E↔A symmetry identity that holds under independence and identical
marginals.

**What it verifies:**
1. φ = (10 − √50.005) / 9.9 to 80-digit precision
2. τ bracket: S_F(0.304718) < 2.525 < S_F(0.304719)
3. S_G(τ) bracket via monotonicity of S_G
4. μ_E·φ exactly (μ_E and φ both computed exactly)
5. Gap = μ_E·φ − S_G(τ) with certified enclosure [0.216134, 0.216136]
6. Sanity checks: quadratic residual for φ, S_F residual for τ, positivity of gap

**Key insight:** Under ρ=0 with identical marginals, the fluid-limit
dynamics are symmetric in the two coordinates (u(s) = v(s)). This
symmetry implies that Farmer's per-turn ecological capture equals
Green's per-turn agricultural cost. Integrating over the active
phase gives the closed-form identity μ_E·φ − S_G(τ) for the
asymptotic gap.

**Arithmetic:** mpmath interval arithmetic with 80-digit precision.
All operations use monotonicity of S_G on [0, 1/2] to produce
rigorous enclosures. No ODE solver.

### rho03_production_certificate.py (cost side)

Certifies the COST side of the ρ=0.3 extension: Green remains
affordable through κ_low = 0.620816.

**What it verifies:**
1. 50-band Farmer-active costs: c_up^FA ≤ 1.4723
2. Continuation cost validates within budget
3. Total cost: c_up^total ≤ 2.5148 < 2.525
4. Therefore: K_N^ME / N ≥ 0.620816 − o_p(1)

**Cost notation:**
- c_up^FA = Farmer-active band total ≤ 1.4723 (used in Ω_φ^N discharge)
- c_up^total = c_up^FA + continuation cost ≤ 2.5148 (used in κ_low certification)

**Key dependency:** Uses the maxA-topY cost monotonicity lemma to
validate that Farmer's interference can only reduce Green's
per-band costs.

**Arithmetic:** One-sided Abramowitz-Stegun Φ bounds (published
uniform error ≤ 7.5×10⁻⁸) in both float and Decimal arithmetic.

### rho03_fluid_limit_certificate.py (value side)

Certifies the VALUE side of the ρ=0.3 extension by tracking Green's
actual ecological value through the fluid-limit ODE.

**What it verifies:**
1. The piecewise ODE system for the state (u, v, b_F, b_G) is
   integrated from s = ε to s = κ_low = 0.620816
2. Farmer exhaustion occurs at τ ≤ 0.31291
3. Green's ecological value integral: ∫e(c(s))ds ≥ 2.5785
4. MX threshold: Θ_{0.3} ≤ 2.3642
5. Margin: ≥ 0.214

**Architecture:** During the active phase [ε, τ], both the Farmer
frontier u(s) and the Green cutoff v(s) evolve under a coupled ODE
derived from the exact dynamic survivor law (Lemma lem:dynamic-survivor-law).
At Farmer exhaustion (transversal crossing of b_F = 0), u freezes
and only v continues to evolve. The launch share ε > 0 sidesteps
the non-Lipschitz corner at (u, v) = (0, 0); the trajectory is
glued from this launch via the launch lemma (Proposition
prop:rho03-initial-eps in the paper). Green's per-turn ecological
value is e(c(s)) = 10 − 9.9·v(s) throughout. The value integral
is accumulated from m = 0.10 onward.

**Arithmetic:** Conservative interval-Euler enclosure with
monotonicity-based one-sided bounds. All Φ evaluations use
Abramowitz-Stegun one-sided bounds (published error ≤ 7.5×10⁻⁸).
The margin (0.214) is large relative to the cumulative Φ-pad error
(~6×10⁻⁴).

## How the Certificates Combine

### Baseline (ρ=0) — single script

`baseline_closed_form_certificate.py` alone certifies the theorem
via the closed-form identity. The exact E↔A symmetry makes
separate cost and value certificates unnecessary.

### Extension (ρ=0.3) — two scripts

1. **rho03_production_certificate.py** proves Green can afford to
   keep shopping through κ_low (cost side)
2. **rho03_fluid_limit_certificate.py** proves the ecological
   value Green accumulates exceeds the MX benchmark (value side)

Together: Green buys enough items (cost certificate) and those
items have enough ecological value (fluid-limit certificate) to
dominate MX.

## Cross-References to Paper

| Certificate output | Paper proposition |
|-------------------|-------------------|
| φ to 80-digit precision | Proposition (MX Farmer-exhaustion share) |
| τ bracket | Proof of Theorem 4.1 (numerical evaluation) |
| μ_E·φ − S_G(τ) ∈ [0.216134, 0.216136] | Theorem 4.1 (exact gap) |
| c_up^FA ≤ 1.4723 | Cost certificate / Ω_φ^N discharge |
| c_up^total ≤ 2.515 | Cost certificate corollary (κ_low) |
| κ_low = 0.620816 | Affordability bootstrap (Lemma lem:baseline-affordability-compact, ρ=0.3 version) |
| ∫e(c(s))ds ≥ 2.579 | Proposition prop:rho03-global-fluid-convergence |
| ρ=0.3 margin ≥ 0.214 | Corollary cor:rho03-certificate |

## Proof Architecture Summary

### Baseline (ρ=0) — E↔A Symmetry Identity

1. Dynamic survivor law → factored host under independence
2. Symmetric ODE: u(s) = v(s) = 1 − √(1−2s)
3. E↔A symmetry: Farmer's per-turn E equals Green's per-turn A
4. Closed-form identity: asymptotic gap = μ_E·φ − S_G(τ)
5. **Certified gap: 0.216135 (exact)**

### Extension (ρ=0.3) — ODE Fluid Limit

1. Prefix wedge: +0.309 (amplified by correlation)
2. Dynamic survivor law → closed 4D state (u, v, b_F, b_G)
3. Launch lemma → joint-mass constraint identifies Z(ε)
4. Cost side: maxA-topY lemma → band costs valid → κ_low = 0.621
5. Ω_φ^N discharge: Farmer-active cost 1.473 ≪ budget 2.525
6. Affordability bootstrap: Green plays exact argmax-Y on [m, κ_low]
7. Value side: fluid-limit ODE tracking v(s) = e(c(s))
8. Piecewise ODE: active phase (both evolve) → post-Farmer (u freezes)
9. **Certified margin: ≥ 0.214**

## Requirements

mpmath>=1.3.0
scipy>=1.10.0
numpy>=1.24.0
Python ≥ 3.10.

## Reproducing the Certificates

Each script runs independently with no arguments:

```bash
python baseline_closed_form_certificate.py
python rho03_production_certificate.py
python rho03_fluid_limit_certificate.py
```

Expected runtime: ~2 seconds for baseline (closed-form), ~5 minutes
for ρ=0.3 cost certificate, ~30 seconds for fluid-limit certificate.

Scripts are deterministic. Identical output expected on any machine
with the specified library versions.

## Zenodo Archive

DOI 10.5281/zenodo.19598799

## Citation

Weinhold, D. and Andersen, L. (2026). "Adversarial Procurement in
Two-Value Space." [Journal TBD].

## License

MIT License.







