# Certificate Inventory for GitHub/Zenodo Archive (v5)
## Picard step-tube architecture with E↔A symmetry closure at ρ=0

---

## BASELINE THEOREM (ρ=0) — E↔A SYMMETRY IDENTITY

### 1. baseline_closed_form_certificate.py — PRIMARY
**Status:** Complete. PASSES.
**Certifies:**
- φ = (10 − √50.005) / 9.9 ≈ 0.295816024604... (80-digit enclosure)
- τ bracket: 0.304718 < τ < 0.304719 (from S_F(τ) = 2.525)
- S_G(τ) bracket: 1.277734 < S_G(τ) < 1.277738
- μ_E·φ bracket: 1.493871 (exact at 80 digits, μ_E = 5.05 exact)
- **Certified asymptotic gap: [0.216134, 0.216136]**

**Arithmetic:** mpmath interval arithmetic, 80-digit precision.
Pure closed-form evaluation of the E↔A symmetry identity
μ_E·φ − S_G(τ). No ODE solver, no prefix wedge.

**Outputs:**
- `baseline_closed_form_certificate.txt` (report)
- `baseline_closed_form_certificate_summary.json` (machine-readable)

---

## ρ=0.3 EXTENSION — PICARD STEP-TUBE CERTIFICATE

### 2. rho03_production_certificate.py — COST SIDE (Picard)
**Status:** Complete. PASSES.

**Certifies at each of 6209 integration steps (h = 10⁻⁴):**
- Picard inclusion `B_n + [0,h]·F(T_n) ⊆ T_n` by direct
  coordinate-by-coordinate containment
- 96-cell G_A interval quadrature with strictly positive denominator
- Diagonal small-time bootstrap `s_0·F_diag(U_trial) ≤ U_trial`
- Three-phase state machine: 3127 active + 10 uncertain + 3072 post

**Outputs:**
- `c_G(κ_low) ∈ [2.353234170771, 2.371402296351]`
- Cost slack lower bound: ≥ 0.153597703649
- τ ≤ 0.31371 (Farmer exhaustion upper time)
- Frozen u-interval at κ_low: [0.426937, 0.430911]

**Arithmetic:** IEEE 754 double with `math.nextafter` outward
rounding. One-sided Abramowitz-Stegun Φ bounds (uniform error
≤ 7.5×10⁻⁸) for all Φ evaluations; verified inverse-normal brackets
via binary search against A-S bounds; outward pad 10⁻¹².

**Load-bearing for:**
- Lemma `lem:rho03-ode-cost-certificate` in `app:rho03-validated-certificates`
- Proposition `prop:rho03-kappa-band-appendix` (affordability bootstrap)

### 3. rho03_fluid_limit_certificate.py — VALUE SIDE (Picard)
**Status:** Complete. PASSES.

**Certifies on the same grid as the cost certificate:**
- Farmer exhaustion: τ ≤ 0.31371
- Value accumulator: `v_out(κ_low) ∈ [2.562783061818, 2.585467093762]`
- MX threshold: Θ_{0.3} ≤ 2.36413
- **Certified margin: ≥ 0.198653**

**Architecture:** Picard step-tube integration of the augmented
fluid-limit ODE with validated accumulators `c_G` (Green cost) and
`v_out` (outside-prefix ecological value accumulator, starting at 0
at s = m = 0.10). The active phase uses coupled `F_1`, `F_2` frontier
drifts with budget drifts `ḃ_F = −(ν − d·u)` and `ḃ_G = −G_A(u, v)`;
the post-Farmer phase freezes u and continues v. The launch share
s_0 = 10⁻⁵ sidesteps the non-Lipschitz corner at the origin.

**Arithmetic:** Same as cost certificate. The Picard tube inclusion
removes the endpoint-Euler enclosure gap; the margin (0.198653) is
well above all accumulated numerical error.

**Load-bearing for:**
- Lemma `lem:rho03-ode-value-certificate`
- Corollary `cor:rho03-certificate`
- Theorem `thm:bronze-rho`

### 4. rho03_production_certificate and rho03_fluid_limit_certificate outputs
- `.txt`: human-readable reports
- `.json`: machine-readable summaries (both contain the unified
  accumulator intervals from the shared Picard grid)

### 5. rho03_static_quadrature_intervals.json — SUPPLEMENTARY
**Status:** Complete. PASSES.

Static scalar enclosures produced by outward-rounded interval
quadrature under the same Abramowitz-Stegun Φ-bounds and IEEE
nextafter rounding discipline as the Picard scripts.

**Keys:**
- `ratio_cuts_rho03.r_prefix_interval = [3.62406450, 3.62406459]`
- `ratio_cuts_rho03.r_phi_interval = [1.50517830, 1.50517845]`
- `mx_prefix_value_rho03.M_X_interval = [0.59221770, 0.59221775]`
- `mx_affordability_rho03.C_R_interval = [0.679, 0.685]`
  (margin ≥ 1.83 below μ_A/2 = 2.525)

**Load-bearing for:**
- Equation `eq:ratio-cuts-rho03`
- Proposition `prop:rho03-prefix-appendix` (M_X prefix value)
- Lemma `lem:mx-green-affordability-rho03` (C_R affordability margin)

---

## RETIRED (removed from repo)

### baseline_fluid_limit_certificate.py (old, margin 0.168)
Superseded by `baseline_closed_form_certificate.py`. The old
certificate used a prefix-wedge + outside-integral decomposition
that no longer appears in the paper; the current proof uses the
E↔A symmetry identity μ_E·φ − S_G(τ), yielding the exact
asymptotic gap 0.216135 rather than a lower bound of 0.168.

### baseline_production_certificate.py (old worst-case, margin 0.046)
Superseded. The worst-case deletion architecture it certifies is
not present in the current paper.

### rho03_production_certificate.csv (old 50-band cost table)
Superseded by the Picard architecture. The old per-band cost table
applied to the endpoint-Euler integration; the Picard certificate
tracks unified accumulators rather than per-band costs, so the
band-level CSV is no longer produced.

### Endpoint-Euler versions of the ρ=0.3 scripts
The previous endpoint-Euler rho03 certificates (which reported
cost slack 0.01022585, τ ≤ 0.31291, value margin 0.21452) are
superseded by the Picard step-tube versions. The Picard
architecture gives genuinely validated enclosures at the cost of a
slightly smaller value margin (0.198653 vs 0.21452) and a larger
cost slack (0.153597 vs 0.01023).

---

## REPOSITORY STRUCTURE

```
certificates/
├── README.md
├── requirements.txt
├── baseline/
│   ├── baseline_closed_form_certificate.py          # PRIMARY
│   ├── baseline_closed_form_certificate.txt
│   └── baseline_closed_form_certificate_summary.json
└── rho03/
    ├── rho03_production_certificate.py              # cost side (Picard)
    ├── rho03_production_certificate.txt
    ├── rho03_production_certificate.json
    ├── rho03_fluid_limit_certificate.py             # value side (Picard)
    ├── rho03_fluid_limit_certificate.txt
    ├── rho03_fluid_limit_certificate.json
    └── rho03_static_quadrature_intervals.json       # supplementary
```

---

## ACTION ITEMS

1. **Upload to Zenodo** (DOI 10.5281/zenodo.19598799). Replace the
   existing archive contents with the current versions. Confirm
   both JSON reports contain the unified accumulator intervals
   `green_cost_interval = [2.3532341707712034, 2.3714022963505363]`
   and `outside_value_interval = [2.5627830618184935, 2.5854670937623934]`.

2. **Verify locally** — run all scripts and confirm all PASS:
   - `baseline_closed_form_certificate.py` → gap ∈ [0.216134, 0.216136]
   - `rho03_production_certificate.py` → cost slack ≥ 0.153597
   - `rho03_fluid_limit_certificate.py` → margin ≥ 0.198653

3. **Confirm cross-references** between certificate outputs and
   paper propositions match the final manuscript version.
