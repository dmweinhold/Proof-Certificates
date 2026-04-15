# Certificate Inventory for GitHub Archive (v3)
## Reflecting the unified fluid-limit architecture

---

## BASELINE THEOREM (ρ=0) — CLOSED-FORM FLUID LIMIT

### 1. baseline_fluid_limit_certificate.py — PRIMARY
**Status:** Complete. Fully hardened. PASSES.
**Certifies:**
- τ bracket: 0.30471 < τ < 0.30472 (from S_F(τ) = 2.525)
- Complete allocation: κ = 1 − τ (exact budget identity)
- V_out ≥ 2.82336 (closed-form value integral on τ-bracket)
- Monotonicity: dV_out/dτ = −0.1 − 4.95√(1−2τ) < 0
- Θ₀ = 2.65513
- Final margin: ≥ 0.16823
**Arithmetic:** mpmath interval arithmetic, 80-digit precision.
Pure closed-form evaluation — no ODE solver.
**Location:** Pro produced this script and verified it.
Three independent calculations confirm the margin:
- Pro (mpmath): [0.16823, 0.16863]
- Claude (interval-Euler): ≥ 0.16814
- Claude (closed-form float): 0.16842

### 2. baseline_production_certificate.py — ALTERNATIVE (Online Appendix C)
**Status:** Complete. PASSES. Retained for the old worst-case proof.
**Certifies:**
- τ̄ bracket: 0.30471 < τ̄ < 0.30472
- Monotonicity of τ → L_ME(τ, 1−τ) on [0.30, 0.31]
- L_ME(0.30472, 1−0.30472) > 2.7013307
- Θ₀ = 2.6551292
- Final margin > 0.04620
**Arithmetic:** Python Decimal with directed rounding (80 digits).
**Note:** This certificate supports the old worst-case baseline
proof architecture, which is retained in Online Appendix C as an
alternative proof. It is NOT the primary baseline certificate.

---

## ρ=0.3 EXTENSION — FLUID LIMIT + COST CERTIFICATE

### 3. rho03_production_certificate.py — COST SIDE
**Status:** Complete. Fully hardened. PASSES.
**Certifies:**
- 50-band Farmer-active costs: c_up^FA ≤ 1.4723
- Continuation cost: C^GO validates within budget
- Total cost: c_up^total ≤ 2.5148 < 2.525
- κ_low = 0.620816 (Green buys at least this fraction)
**Arithmetic:** One-sided Abramowitz-Stegun Φ bounds (published
error ≤ 7.5×10⁻⁸) in float and Decimal.
**Location:** /mnt/user-data/uploads/rho03_production_certificate.py
**Note on cost notation:**
- c_up^FA = Farmer-active band total ≤ 1.4723 (Ω_φ^N discharge)
- c_up^total = c_up^FA + continuation ≤ 2.5148 (κ_low certification)

### 4. rho03_fluid_limit_certificate.py — VALUE SIDE
**Status:** Being hardened. Pro is replacing standard-library Φ/Φ⁻¹
with interval-safe A-S enclosures (same infrastructure as the cost
certificate). Margin is 0.214 — hardening will not change the
result meaningfully (cumulative Φ-pad error ~6×10⁻⁴).
**Certifies:**
- Farmer exhaustion: τ ≤ 0.31291
- Fluid-limit value integral: ∫_m^{κ_low} e(c(s))ds ≥ 2.5785
- MX threshold: Θ_{0.3} = 2.3642
- Final margin: ≥ 0.214
**Arithmetic (after hardening):** Conservative interval-Euler
enclosure with monotonicity-based one-sided bounds. All Φ
evaluations use A-S one-sided bounds (error ≤ 7.5×10⁻⁸).
**Location:** Pilot at /mnt/user-data/outputs/rho03_fluid_limit_certificate.py;
hardened version pending from Pro.

### 5. rho03_production_certificate outputs
**Status:** Complete. Keep for reference.
- .txt: human-readable report
- .csv: 50-band cost table (certified per-band upper bounds)
- .json: machine-readable summary
**Note:** The old value-side numbers in these outputs (V^FA_low,
ΔM^GO, old margin 0.00261) are superseded by the fluid-limit
certificate. The COST numbers remain valid and load-bearing.

---

## RETIRED

### Old baseline margin 0.046
Superseded by the fluid-limit margin 0.168 as the primary proof.
The old certificate (baseline_production_certificate.py) supports
Online Appendix C only.

---

## REPOSITORY STRUCTURE

```
certificates/
├── README.md
├── requirements.txt
├── baseline/
│   ├── baseline_fluid_limit_certificate.py   # PRIMARY
│   ├── baseline_fluid_limit_certificate.txt
│   ├── baseline_fluid_limit_certificate.json
│   ├── baseline_production_certificate.py    # Online Appendix C
│   └── baseline_production_certificate.txt
└── rho03/
    ├── rho03_production_certificate.py       # cost side
    ├── rho03_production_certificate.txt
    ├── rho03_production_certificate.csv      # 50-band cost table
    ├── rho03_production_certificate.json
    ├── rho03_fluid_limit_certificate.py      # value side
    └── rho03_fluid_limit_certificate.txt
```

---

## ACTION ITEMS

1. **Receive hardened rho03_fluid_limit_certificate.py from Pro**
   — replace pilot with hardened version using A-S Φ bounds.

2. **Verify locally** — run all scripts on your machine and
   confirm all PASS:
   - baseline_fluid_limit_certificate.py → margin ≥ 0.168
   - rho03_production_certificate.py → cost slack > 0
   - rho03_fluid_limit_certificate.py → margin ≥ 0.214

3. **Update requirements.txt:**
   ```
   mpmath>=1.3.0
   scipy>=1.10.0
   numpy>=1.24.0
   ```

4. **Upload to Zenodo** (DOI: 10.5281/zenodo.17114490)

5. **Confirm cross-references** between certificate outputs and
   paper propositions match the final manuscript version.
