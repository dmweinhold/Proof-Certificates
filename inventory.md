# Certificate Inventory for GitHub Archive (v4)
## Unified fluid-limit architecture with E↔A symmetry closure at ρ=0

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
**Output:** baseline_closed_form_certificate.txt (report),
baseline_closed_form_certificate_summary.json (machine-readable).

---

## ρ=0.3 EXTENSION — FLUID LIMIT + COST CERTIFICATE

### 2. rho03_production_certificate.py — COST SIDE
**Status:** Complete. Fully hardened. PASSES.
**Certifies:**
- 50-band Farmer-active costs: c_up^FA ≤ 1.4723
- Continuation cost: C^GO validates within budget
- Total cost: c_up^total ≤ 2.5148 < 2.525
- κ_low = 0.620816 (certified Green affordability horizon)
**Arithmetic:** One-sided Abramowitz-Stegun Φ bounds (published
uniform error ≤ 7.5×10⁻⁸) in float and Decimal.
**Cost notation:**
- c_up^FA = Farmer-active band total ≤ 1.4723 (Ω_φ^N discharge)
- c_up^total = c_up^FA + continuation ≤ 2.5148 (κ_low certification)
**Load-bearing for:** Affordability lemma and launch lemma in
Section app:bronze-rho03 of the paper.

### 3. rho03_fluid_limit_certificate.py — VALUE SIDE
**Status:** Being hardened. Pro is replacing standard-library Φ/Φ⁻¹
with interval-safe A-S enclosures (same infrastructure as the cost
certificate). Margin is 0.214 — hardening will not change the
result meaningfully (cumulative Φ-pad error ~6×10⁻⁴).
**Certifies:**
- Farmer exhaustion: τ ≤ 0.31291
- Fluid-limit value integral: ∫_m^{κ_low} e(c(s)) ds ≥ 2.5785
- MX threshold: Θ_{0.3} = 2.3642
- **Certified margin: ≥ 0.214**
**Arithmetic (after hardening):** Conservative interval-Euler
enclosure with monotonicity-based one-sided bounds. All Φ
evaluations use A-S one-sided bounds (error ≤ 7.5×10⁻⁸).

### 4. rho03_production_certificate outputs
**Status:** Complete. Keep for reference.
- .txt: human-readable report
- .csv: 50-band cost table (certified per-band upper bounds)
- .json: machine-readable summary

---

## RETIRED (removed from repo)

### baseline_fluid_limit_certificate.py (old, margin 0.168)
Superseded by baseline_closed_form_certificate.py. The old
certificate used a prefix-wedge + outside-integral decomposition
that no longer appears in the paper; the current proof uses the
E↔A symmetry identity μ_E·φ − S_G(τ), yielding the exact
asymptotic gap 0.216135 rather than a lower bound of 0.168.

### baseline_production_certificate.py (old worst-case, margin 0.046)
Superseded. The worst-case deletion architecture it certifies is
not present in the current paper. No "Online Appendix C" retains
it; the old proof has been fully deleted, not relegated.

---

## REPOSITORY STRUCTURE

certificates/
├── README.md
├── requirements.txt
├── baseline/
│   ├── baseline_closed_form_certificate.py          # PRIMARY
│   ├── baseline_closed_form_certificate.txt
│   └── baseline_closed_form_certificate_summary.json
└── rho03/
├── rho03_production_certificate.py              # cost side
├── rho03_production_certificate.txt
├── rho03_production_certificate.csv             # 50-band cost table
├── rho03_production_certificate.json
├── rho03_fluid_limit_certificate.py             # value side
└── rho03_fluid_limit_certificate.txt


---

## ACTION ITEMS

1. **Receive hardened rho03_fluid_limit_certificate.py from Pro** —
   replace pilot with hardened version using A-S Φ bounds.

2. **Verify locally** — run all scripts and confirm all PASS:
   - baseline_closed_form_certificate.py → gap ∈ [0.216134, 0.216136]
   - rho03_production_certificate.py → cost slack > 0
   - rho03_fluid_limit_certificate.py → margin ≥ 0.214

3. **Upload to Zenodo** (DOI 10.5281/zenodo.19598799).

4. **Confirm cross-references** between certificate outputs and
   paper propositions match the final manuscript version
   (Adversarial_Procurement_Claude_v2.tex).







