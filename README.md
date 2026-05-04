# Validated Certificates for "Adversarial Procurement in Two-Value Space"

## Overview

This repository contains the validated numerical certificates that
support the theoretical results in the paper. The proof architecture
distinguishes between **theorem-critical** certificates (whose intervals
close a step in the formal proof chain) and **supplementary** certificates
(which document the quantitative magnitude of effects established
analytically elsewhere).

- **Theorem 5.1** (`thm:claims-tilt`, Claims World, finite-board priority-tilt
  monotonicity): Proved pathwise, distribution-free, and analytically
  in the paper's in-paper appendix. **No numerical certificate required.**
  A supplementary iid-uniform fluid-limit benchmark is archived here for
  scale calibration only.
- **Theorem 5.2** (`thm:bronze`, Budget World, ρ=0): Asymptotic gap
  ∈ [0.216134, 0.216136] (closed-form E↔A symmetry identity, validated
  to 80-digit precision).
- **Theorem 5.3** (`thm:bronze-rho`, Budget World, ρ=0.3): Certified
  margin ≥ 0.198653 per plot (Picard step-tube certificate).

The Budget-World baseline theorem is proved via an exact E↔A symmetry of the
fluid-limit dynamics under independence, which yields the asymptotic gap
μ_E·φ − S_G(τ) in closed form. The ρ=0.3 extension uses the general
fluid-limit architecture that tracks Green's actual per-turn ecological
value and agricultural cost through a deterministic ODE on the survivor
pool state, closed by a validated Picard step-tube enclosure with verified
inclusion at every step.

## Quick Start

```bash
pip install -r requirements.txt

# Theorem-critical (Budget World)
python certificates/baseline/baseline_closed_form_certificate.py
python certificates/rho03/rho03_production_certificate.py
python certificates/rho03/rho03_fluid_limit_certificate.py

# Supplementary (Claims World iid-uniform benchmark)
python certificates/claims/claims_uniform_ratio_frontier_interval_certificate.py
```

Each script prints a self-contained report ending with PASS or FAIL.

## Theorem-Critical Certificates

### Budget World, ρ=0 — `thm:bronze` (Theorem 5.2)

#### baseline_closed_form_certificate.py — PRIMARY for ρ=0

Certifies the baseline asymptotic dominance theorem using the exact
E↔A symmetry identity that holds under independence and identical
marginals.

**What it verifies:**
1. φ = (10 − √50.005) / 9.9 to 80-digit precision
2. τ bracket: S_F(0.304718) < 2.525 < S_F(0.304719)
3. S_G(τ) bracket via monotonicity of S_G
4. μ_E·φ exactly (μ_E and φ both computed exactly)
5. Gap = μ_E·φ − S_G(τ) with certified enclosure [0.216134, 0.216136]
6. Sanity checks: quadratic residual for φ, S_F residual for τ,
   positivity of gap

**Key insight:** Under ρ=0 with identical marginals, the fluid-limit
dynamics are symmetric in the two coordinates (u(s) = v(s)). This
symmetry implies that Farmer's per-turn ecological capture equals
Green's per-turn agricultural cost. Integrating over the active
phase gives the closed-form identity μ_E·φ − S_G(τ) for the
asymptotic gap.

**Arithmetic:** mpmath interval arithmetic with 80-digit precision.
All operations use monotonicity of S_G on [0, 1/2] to produce
rigorous enclosures. No ODE solver.

**Outputs:**
- `baseline_closed_form_certificate.txt` (human-readable report)
- `baseline_closed_form_certificate_summary.json` (machine-readable)

### Budget World, ρ=0.3 — `thm:bronze-rho` (Theorem 5.3)

#### rho03_production_certificate.py (cost side, Picard step-tube)

Certifies the **cost side** of the ρ=0.3 extension through a validated
Picard step-tube enclosure: Green's agricultural cost accumulator
through κ_low = 0.620816 stays strictly below the budget 2.525.

**What it verifies:**
At each of the 6209 integration steps (h = 10⁻⁴):
1. Picard inclusion `B_n + [0,h]·F(T_n) ⊆ T_n` is verified by direct
   coordinate-by-coordinate interval containment check
2. G_A interval quadrature (96 cells) returns intervals with strictly
   positive denominator lower endpoint
3. Diagonal small-time bootstrap `s_0·F_diag(U_trial) ≤ U_trial`
4. Three-phase state machine {active, uncertain, post} across the
   Farmer-exhaustion transition (3127 active + 10 uncertain + 3072 post)

**Outputs:**
- `c_G(κ_low) ∈ [2.353234, 2.371402]` (Green cost accumulator)
- Cost slack lower bound ≥ 0.153597
- Farmer exhaustion upper time τ ≤ 0.31371
- Frozen u-interval at κ_low: [0.426937, 0.430911]

**Arithmetic:** IEEE 754 double with `math.nextafter` outward
rounding. One-sided Abramowitz-Stegun Φ bounds (uniform error
≤ 7.5×10⁻⁸) for all Φ evaluations. Verified inverse-normal brackets
via binary search against the A-S bounds.

#### rho03_fluid_limit_certificate.py (value side, Picard step-tube)

Certifies the **value side** of the ρ=0.3 extension through the same
Picard step-tube integration, with an augmented validated state that
tracks Green's outside-prefix ecological value accumulator.

**What it verifies:**
On the same grid as the cost certificate:
1. Farmer exhaustion: τ ≤ 0.31371
2. Value accumulator: `v_out(κ_low) ∈ [2.562783, 2.585467]`
3. MX threshold: Θ_{0.3} ≤ 2.36413
4. Certified margin: ≥ 0.198653

**Architecture:** During the active phase [s_0, τ], both the Farmer
frontier u(s) and the Green cutoff v(s) evolve under the coupled
fluid-limit ODE with budget and accumulator drifts. At Farmer
exhaustion (transversal b_F = 0 crossing), u freezes and only
(v, b_G, c_G, v_out) continue to evolve. The launch share s_0 = 10⁻⁵
sidesteps the non-Lipschitz corner at (u, v) = (0, 0); the trajectory
is glued from this launch via a validated diagonal bootstrap.

**Green's per-step cost** is `G_A(u, v) = E[A | Y = q(v), X ≤ q(u)]`,
the conditional agricultural cost given Green's max-Y selection over
the Farmer-truncated host. Green's per-step ecological value is
`e(c(s)) = 10 − 9.9·v(s)`. The value integral is accumulated from
m = 0.10 onward.

**Arithmetic:** Same as the cost certificate, plus Picard tube
inclusion verified at every step.

#### rho03_static_quadrature_intervals.json (supplementary scalars)

Static scalar enclosures that enter the ρ=0.3 theorem chain but are
not produced by the Picard step integration: the ratio-cut thresholds
r_m(0.3) and r_φ(0.3), the MX prefix value M_X(0.10, 0.3), and the
MX-side Green affordability margin C_R^{(0.3)}(r_φ(0.3)).

**Keys:**
- `ratio_cuts_rho03` — `r_m_interval`, `r_phi_interval`
- `mx_prefix_value_rho03` — `M_X_interval = [0.59221770, 0.59221775]`
- `mx_affordability_rho03` — `C_R_interval = [0.679, 0.685]`, margin
  ≥ 1.83 below μ_A/2 = 2.525

These intervals are produced by the same outward-rounded quadrature
primitives used in the Picard step certificates (Abramowitz-Stegun
Φ bounds, IEEE nextafter rounding, verified inverse-normal brackets).

## Supplementary Benchmark Certificates

### Claims World iid-uniform fluid benchmark

**Important framing.** Theorem 5.1 (`thm:claims-tilt`) is proved
analytically in the paper's in-paper appendix as a finite-board,
pathwise, distribution-free weak-dominance result under the max-A
Farmer rule. **The certificate below is not part of that proof.** It
is a validated iid-uniform fluid-limit benchmark showing that the
theorem's weak dominance is quantitatively large in the symmetric
baseline, not merely nonnegative. Its role is scale calibration —
documenting that the iid-uniform asymptotic gap is on the order of
half a unit per plot, which connects the analytic theorem to the
magnitude of effects observed in simulations.

#### claims_uniform_ratio_frontier_interval_certificate.py

Certifies the asymptotic Claims World iid-uniform fluid-limit gap at
budget share α = 1/2 (equal claims).

**What it verifies:**
1. V_R(1/2) ∈ [2.8783742, 2.8789214] (MX ratio-frontier ecological value)
2. V_E(1/2) = 3.35 exactly (Max Environment ecological value)
3. Gap = V_E(1/2) − V_R(1/2) ∈ [0.4710786, 0.4716258]
4. Certificate tests: V_R upper bound < 2.88, gap lower bound > 0.47

**Key construction:** The certificate decomposes V_R into three
analytic integrals (I1, I2, I3) over the three regimes of the
ratio-frontier dynamics, computes interval enclosures using analytic
monotonicity of the Regime II and Regime III integrands plus
left/right Riemann sums, and verifies the four transition points
(r12, a12, a1, rF) by interval root-bracketing.

**Arithmetic:** mpmath interval arithmetic at 90-digit precision, with
2000 subdivisions per monotone integral. No simulation, no random
seed. The certificate is deliberately loose — the bounds (V_R < 2.88,
gap > 0.47) are well inside the certified enclosures, providing ample
slack for the calibration claim.

**Outputs:**
- `claims_uniform_ratio_frontier_interval_certificate.txt`
- `claims_uniform_ratio_frontier_interval_certificate.json`
- `claims_uniform_ratio_frontier_interval_certificate.md`

#### claims_uniform_ratio_frontier_certificate.py (diagnostic)

A high-precision analytic / non-interval diagnostic companion to the
interval wrapper above. Useful for inspecting the regime structure
and verifying intermediate quantities; not the validated certificate.

## Frozen Theorem Objects

The Budget-World theorems share the following setup:

| Object | Definition |
|--------|-----------|
| DGP | Marginals U[0.1, 10], Gaussian copula |
| Farmer strategy | Naïve max-A (buys most expensive available) |
| Green strategy ME | Max-E (buys highest ecological value) |
| Green strategy MX | Max-E/A (buys highest ratio available) |
| Budget parity | B_F = B_G = (1/2) Σ A_i |
| Outcome | Purchased Conservation (PC) |
| Prefix share | m = 0.10 |

The Claims World benchmark uses the same DGP and Farmer rule but
operates in the equal-claims, full-leakage regime without prices or
budgets.

## Frozen Constants

### Shared

| Constant | Value | Definition |
|----------|-------|-----------|
| ℓ | 0.1 | Lower bound of marginal support |
| ν | 10 | Upper bound of marginal support |
| μ_E = μ_A | 5.05 | Marginal means (identical) |
| Budget/N | 2.525 | μ_A/2 |
| m | 0.10 | Prefix share |

### Budget World baseline (ρ=0) — E↔A Symmetry Identity

| Constant | Value | Definition |
|----------|-------|-----------|
| φ | 0.295816024604... | MX Farmer-exhaustion share; root of 9.9φ² − 20φ + 5.05 = 0 |
| τ | 0.304718564833... | ME Farmer-exhaustion share; root of S_F(τ) = 2.525 |
| μ_E·φ | 1.493871 | MX-side total ecological value |
| S_G(τ) | 1.277736 | Farmer's ecological capture under ME (= Green's A-spend, by symmetry) |
| **Certified gap** | **[0.216134, 0.216136]** | **μ_E·φ − S_G(τ) (exact identity)** |

### Budget World extension (ρ=0.3) — Picard Step-Tube Certificate

| Constant | Value | Source |
|----------|-------|--------|
| ρ | 0.3 | Latent Gaussian correlation |
| σ | √(1−0.09) ≈ 0.9539 | Conditional standard deviation |
| φ | 0.295816 | Benchmark split |
| κ_low | 0.620816 | Certified Green affordability horizon |
| τ (fluid) | ≤ 0.31371 | Farmer exhaustion in fluid limit (Picard) |
| c_G(κ_low) | ∈ [2.353234, 2.371402] | Green cost accumulator (Picard) |
| Cost slack | ≥ 0.153597 | 2.525 − c_G upper |
| v_out(κ_low) | ∈ [2.562783, 2.585467] | Outside-prefix value accumulator (Picard) |
| Θ_{0.3} | ≤ 2.36413 | MX threshold (static quadrature) |
| **Certified margin** | **≥ 0.198653** | v_out lower − Θ_{0.3} |

### Claims World iid-uniform benchmark (supplementary)

| Constant | Value | Definition |
|----------|-------|-----------|
| V_E(1/2) | 3.35 (exact) | ME asymptotic ecological value at α = 1/2 |
| V_R(1/2) | ∈ [2.8783742, 2.8789214] | MX ratio-frontier asymptotic ecological value at α = 1/2 |
| **Certified gap** | **∈ [0.4710786, 0.4716258]** | **V_E − V_R (validated interval)** |

## Repository Structure

```
certificates/
├── README.md                                              (this file)
├── requirements.txt
├── baseline/
│   ├── baseline_closed_form_certificate.py                # E↔A symmetry identity
│   ├── baseline_closed_form_certificate.txt
│   └── baseline_closed_form_certificate_summary.json
├── rho03/
│   ├── rho03_production_certificate.py                    # Picard cost-side report
│   ├── rho03_production_certificate.txt
│   ├── rho03_production_certificate.json
│   ├── rho03_fluid_limit_certificate.py                   # Picard value-side report
│   ├── rho03_fluid_limit_certificate.txt
│   ├── rho03_fluid_limit_certificate.json
│   └── rho03_static_quadrature_intervals.json             # Supplementary static-quadrature intervals
└── claims/
    ├── claims_uniform_ratio_frontier_certificate.py             # Diagnostic (non-interval)
    ├── claims_uniform_ratio_frontier_interval_certificate.py    # Validated certificate
    ├── claims_uniform_ratio_frontier_interval_certificate.txt
    ├── claims_uniform_ratio_frontier_interval_certificate.json
    └── claims_uniform_ratio_frontier_interval_certificate.md
```

## Cross-References to Paper

| Certificate output | Paper reference | Role |
|-------------------|-----------------|------|
| Pathwise Claims theorem | Theorem 5.1 (`thm:claims-tilt`) | No certificate; analytic proof in `app:claims_proof` |
| φ to 80-digit precision | Proposition (MX Farmer-exhaustion share) | Theorem-critical |
| τ (baseline) bracket | Proof of Theorem 5.2 (`thm:bronze`) | Theorem-critical |
| μ_E·φ − S_G(τ) ∈ [0.216134, 0.216136] | Theorem 5.2 (exact gap) | Theorem-critical |
| c_G(κ_low) ∈ [2.353234, 2.371402] | Lemma `lem:rho03-ode-cost-certificate` | Theorem-critical |
| Cost slack ≥ 0.153597 | Proposition `prop:rho03-kappa-band-appendix` | Theorem-critical |
| κ_low = 0.620816 | Affordability bootstrap | Theorem-critical |
| v_out(κ_low) ∈ [2.562783, 2.585467] | Lemma `lem:rho03-ode-value-certificate` | Theorem-critical |
| ρ=0.3 margin ≥ 0.198653 | Theorem 5.3 (`thm:bronze-rho`), Corollary `cor:rho03-certificate` | Theorem-critical |
| ratio_cuts_rho03 | Equation `eq:ratio-cuts-rho03` | Theorem-critical |
| mx_prefix_value_rho03 | Proposition `prop:rho03-prefix-appendix` | Theorem-critical |
| mx_affordability_rho03 | Lemma `lem:mx-green-affordability-rho03` | Theorem-critical |
| V_R(1/2), V_E(1/2), gap | Claims World iid-uniform fluid benchmark | Supplementary (calibration only) |

## Proof Architecture Summary

### Claims World, finite-board (Theorem 5.1, `thm:claims-tilt`)

Proved analytically and pathwise in the paper's in-paper appendix
(`app:claims_proof`). The result holds on every generic finite board
under the deterministic max-A Farmer rule, distribution-free in
(e, a). No numerical certificate is required. The supplementary
iid-uniform benchmark certificate above documents the asymptotic
magnitude of the gap in the symmetric baseline.

### Budget World, baseline ρ=0 (Theorem 5.2, `thm:bronze`) — E↔A Symmetry Identity

1. Dynamic survivor law → factored host under independence
2. Symmetric ODE: u(s) = v(s) = 1 − √(1−2s)
3. E↔A symmetry: Farmer's per-turn E equals Green's per-turn A
4. Closed-form identity: asymptotic gap = μ_E·φ − S_G(τ)
5. **Certified gap: [0.216134, 0.216136]**

### Budget World, extension ρ=0.3 (Theorem 5.3, `thm:bronze-rho`) — Picard Step-Tube Certificate

1. Fluid-limit ODE in tail coordinates with correct Green drift
   `ḃ_G = −G_A(u, v)` (conditional agricultural cost, not ecological value)
2. Validated state augmented with cost and value accumulators
3. Picard step-tube inclusion `B_n + [0,h]·F(T_n) ⊆ T_n` verified
   at every step by direct interval containment
4. 96-cell interval quadrature encloses G_A on every tube
5. Three-phase state machine {active, uncertain, post} handles the
   Farmer-exhaustion transition rigorously
6. Cost accumulator: c_G(κ_low) ≤ 2.371402 ⟹ Green stays solvent
7. Value accumulator: v_out(κ_low) ≥ 2.562783
8. MX threshold Θ_{0.3} ≤ 2.36413 from static quadrature
9. **Certified margin: ≥ 0.198653**

## Requirements

```
mpmath>=1.3.0
scipy>=1.10.0
numpy>=1.24.0
```

Python ≥ 3.10.

## Reproducing the Certificates

Each script runs independently with no arguments:

```bash
python baseline_closed_form_certificate.py
python rho03_production_certificate.py
python rho03_fluid_limit_certificate.py
python claims_uniform_ratio_frontier_interval_certificate.py
```

Scripts are deterministic. Identical output expected on any machine
with the specified library versions.

## Archive DOIs

This repository hosts **proof certificates** at:
- DOI: **10.5281/zenodo.19598799**

The broader **replication / simulator package** (the dynamic
contest simulator and its outputs) is archived separately at:
- DOI: **10.5281/zenodo.19613421**

The two archives are companions: the proof-certificate archive
backs the formal theorems (Section 5 of the paper); the replication
archive reproduces the simulation framework, figures, and the Bolivia
empirical exercise.

## Citation

Weinhold, D. and Andersen, L. (2026). "Adversarial Procurement in
Two-Value Space: Insights and Evidence for Conservation Siting."
SSRN Working Paper, https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6705959.

## License

MIT License.
