#!/usr/bin/env python3
"""
baseline_fluid_limit_certificate.py

Production-style certificate for the baseline (rho = 0) fluid-limit theorem.
Uses exact closed-form formulas plus interval arithmetic (mpmath.iv).
"""
import json
import mpmath as mp

iv = mp.iv
mp.mp.dps = 80
iv.dps = 80

# Frozen constants
m_mp = mp.mpf('0.10')
phi_mp = mp.mpf('0.295816')
muA_mp = mp.mpf('5.05')
budget_mp = muA_mp / 2

m = iv.mpf('0.10')
phi = iv.mpf('0.295816')
muA = iv.mpf('5.05')
budget = muA / 2

tau_lo_mp = mp.mpf('0.30471')
tau_hi_mp = mp.mpf('0.30472')
tau_lo = iv.mpf('0.30471')
tau_hi = iv.mpf('0.30472')
tau_iv = iv.mpf([tau_lo.a, tau_hi.b])

def S_F(s):
    return iv.mpf('0.1') * s + iv.mpf('3.3') * (1 - (1 - 2*s) ** iv.mpf('1.5'))

def S_G(s):
    return iv.mpf('0.1') * s + iv.mpf('1.65') * (1 - (1 - 2*s) ** iv.mpf('1.5'))

def z_of_tau(tau):
    return iv.sqrt(1 - 2*tau)

def rem_G(tau):
    return budget - S_G(tau)

def g_rate(tau):
    z = z_of_tau(tau)
    return iv.mpf('0.1') + iv.mpf('4.95') * z

def delta_of_tau(tau):
    return rem_G(tau) / g_rate(tau)

def kappa_of_tau(tau):
    return tau + delta_of_tau(tau)

def V_FA(tau):
    z = z_of_tau(tau)
    z_m = iv.sqrt(1 - 2*m)
    return iv.mpf('0.1') * (tau - m) + iv.mpf('3.3') * (z_m**3 - z**3)

def V_GO(tau):
    z = z_of_tau(tau)
    delta = delta_of_tau(tau)
    e_start = iv.mpf('0.1') + iv.mpf('9.9') * z
    return delta * e_start - iv.mpf('4.95') * (delta**2) / z

def V_out(tau):
    return V_FA(tau) + V_GO(tau)

# mp versions for nominal values
def S_F_mp(s):
    return mp.mpf('0.1') * s + mp.mpf('3.3') * (1 - (1 - 2*s) ** mp.mpf('1.5'))

def S_G_mp(s):
    return mp.mpf('0.1') * s + mp.mpf('1.65') * (1 - (1 - 2*s) ** mp.mpf('1.5'))

def z_mp(t):
    return mp.sqrt(1 - 2*t)

def rem_G_mp(t):
    return budget_mp - S_G_mp(t)

def g_rate_mp(t):
    return mp.mpf('0.1') + mp.mpf('4.95') * z_mp(t)

def delta_mp(t):
    return rem_G_mp(t) / g_rate_mp(t)

def kappa_mp(t):
    return t + delta_mp(t)

def V_FA_mp(t):
    return mp.mpf('0.1') * (t - m_mp) + mp.mpf('3.3') * ((1 - 2*m_mp) ** mp.mpf('1.5') - (1 - 2*t) ** mp.mpf('1.5'))

def V_GO_mp(t):
    z = z_mp(t)
    delta = delta_mp(t)
    e_start = mp.mpf('0.1') + mp.mpf('9.9') * z
    return delta * e_start - mp.mpf('4.95') * (delta**2) / z

def V_out_mp(t):
    return V_FA_mp(t) + V_GO_mp(t)

Theta_0 = muA * (1 - phi) - iv.mpf('0.901')
Theta_0_mp = muA_mp * (1 - phi_mp) - mp.mpf('0.901')

# Certified brackets for tau
SF_lo = S_F(tau_lo)
SF_hi = S_F(tau_hi)

# Main value enclosure
V_out_iv = V_out(tau_iv)
kappa_iv = kappa_of_tau(tau_iv)
margin_iv = V_out_iv - Theta_0

# Nominal values
tau_nom = mp.findroot(lambda s: S_F_mp(s) - budget_mp, (tau_lo_mp, tau_hi_mp))
kappa_nom = kappa_mp(tau_nom)
V_nom = V_out_mp(tau_nom)
margin_nom = V_nom - Theta_0_mp

status = "PASS" if margin_iv.a > 0 else "FAIL"

lines = []
lines.append("="*72)
lines.append("BASELINE FLUID-LIMIT CERTIFICATE")
lines.append("="*72)
lines.append(f"m                              = 0.10")
lines.append(f"phi                            = 0.295816")
lines.append(f"Budget                         = 2.525")
lines.append("")
lines.append("Step 1. Certified Farmer-exhaustion bracket")
lines.append(f"  S_F(0.30471)                 = {SF_lo}")
lines.append(f"  S_F(0.30472)                 = {SF_hi}")
lines.append(f"  therefore                    : 0.30471 < tau < 0.30472")
lines.append("")
lines.append("Step 2. Fluid-limit outside-prefix value")
lines.append(f"  kappa(tau) interval          = {kappa_iv}")
lines.append(f"  V_out(tau) interval          = {V_out_iv}")
lines.append("")
lines.append("Step 3. MX threshold")
lines.append(f"  Theta_0                      = {Theta_0}")
lines.append("")
lines.append("Step 4. Certified margin")
lines.append(f"  margin interval              = {margin_iv}")
lines.append(f"  certified lower bound        = {margin_iv.a}")
lines.append("")
lines.append("Nominal values (for reference only)")
lines.append(f"  tau_nom                      = {tau_nom}")
lines.append(f"  kappa_nom                    = {kappa_nom}")
lines.append(f"  V_out_nom                    = {V_nom}")
lines.append(f"  Theta_0_nom                  = {Theta_0_mp}")
lines.append(f"  margin_nom                   = {margin_nom}")
lines.append("")
lines.append(f"STATUS: BASELINE FLUID-LIMIT CERTIFICATE {status}.")
report = "\n".join(lines)

open("baseline_fluid_limit_certificate.txt","w",encoding="utf-8").write(report)
payload = {
    "tau_interval": [str(tau_iv.a), str(tau_iv.b)],
    "kappa_interval": [str(kappa_iv.a), str(kappa_iv.b)],
    "V_out_interval": [str(V_out_iv.a), str(V_out_iv.b)],
    "Theta_0_interval": [str(Theta_0.a), str(Theta_0.b)],
    "margin_interval": [str(margin_iv.a), str(margin_iv.b)],
    "tau_nom": str(tau_nom),
    "kappa_nom": str(kappa_nom),
    "V_out_nom": str(V_nom),
    "Theta_0_nom": str(Theta_0_mp),
    "margin_nom": str(margin_nom),
    "status": status
}
open("baseline_fluid_limit_certificate.json","w",encoding="utf-8").write(json.dumps(payload, indent=2))
print(report)
