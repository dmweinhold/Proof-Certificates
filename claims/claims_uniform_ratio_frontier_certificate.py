#!/usr/bin/env python3
"""
Claims World iid-uniform ratio-frontier certificate-style computation.

This script computes the fluid-limit value of the ratio rule (MX) in the
full-leakage, equal-claims Claims World game when

    A, E iid Uniform[ell, nu],    Farmer score = A,    Green score = R = E/A.

It uses the corrected dynamic ratio-frontier ODE, not a static top-ratio
threshold.  The main output is

    V_R(1/2) = integral of the Green boundary value along the ratio frontier,
    V_E(1/2) = integral_0^1 p Q(p) dp,
    gap      = V_E(1/2) - V_R(1/2).

For ell=0.1, nu=10, the expected output is approximately

    V_R(1/2)  = 2.8786477650228935672546891300...
    V_E(1/2)  = 3.35
    gap       = 0.4713522349771064327453108700...
    N=100 gap = 47.13522349771064...

The calculation below avoids Euler discretization by using the piecewise
analytic reduction of the uniform-support ratio-frontier path:

  Regime I:   r in [r12, nu/ell],        r a >= nu.
  Regime II:  r in [1, r12],             ell <= r a <= nu, r >= 1.
  Regime III: r in [rF, 1],              r < 1.

It is a high-precision, deterministic, self-checking numerical certificate.
It is not a formal machine-verified interval proof, but it is structured so a
paper-grade interval-arithmetic wrapper can be added around the same formulas.

Only dependency: mpmath.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Iterable, Optional

import mpmath as mp


@dataclass
class CertificateResult:
    ell: mp.mpf
    nu: mp.mpf
    d: mp.mpf
    r0: mp.mpf
    r12: mp.mpf
    a12: mp.mpf
    a1: mp.mpf
    rF: mp.mpf
    aF: mp.mpf
    H12: mp.mpf
    H1: mp.mpf
    I1: mp.mpf
    I2: mp.mpf
    I3: mp.mpf
    VR: mp.mpf
    VE: mp.mpf
    gap: mp.mpf


def _fmt(x: mp.mpf, digits: int = 35) -> str:
    return mp.nstr(x, n=digits, strip_zeros=False)


def _bisect_mp(f, lo: mp.mpf, hi: mp.mpf, *, tol: mp.mpf, max_iter: int = 1000) -> mp.mpf:
    """Simple high-precision bisection for a continuous monotone sign change."""
    flo = f(lo)
    fhi = f(hi)
    if flo == 0:
        return lo
    if fhi == 0:
        return hi
    if flo * fhi > 0:
        raise ValueError(f"Bisection bracket does not straddle zero: f(lo)={flo}, f(hi)={fhi}")

    a, b = lo, hi
    fa, fb = flo, fhi
    for _ in range(max_iter):
        mid = (a + b) / 2
        fm = f(mid)
        if abs(b - a) <= tol or fm == 0:
            return mid
        if fa * fm <= 0:
            b, fb = mid, fm
        else:
            a, fa = mid, fm
    raise RuntimeError("Bisection did not converge")


class UniformRatioFrontier:
    def __init__(self, ell: str = "0.1", nu: str = "10", dps: int = 80):
        mp.mp.dps = dps
        self.dps = dps
        self.ell = mp.mpf(ell)
        self.nu = mp.mpf(nu)
        if self.ell <= 0 or self.nu <= self.ell:
            raise ValueError("Need 0 < ell < nu.")
        self.d = self.nu - self.ell
        self.r0 = self.nu / self.ell
        self._tol = mp.mpf(10) ** (-(dps - 20 if dps > 40 else max(20, dps // 2)))

        # Regime I/II transition: a_regime_I(r12) = nu / r12.
        # This root has the closed form below for the smaller root of the quadratic.
        disc = self.nu**4 - self.ell**2 * (3 * self.nu**2 - 2 * self.nu * self.ell)
        self.r12 = (self.nu**2 - mp.sqrt(disc)) / (self.ell**2)
        self.a12 = self.nu / self.r12

        self._J12 = self.J_regime_II(self.a12)
        self._C12 = self.r12 / (self.a12**2 - self.ell**2)

        # Regime II/III transition: r(a1) = 1.
        # For the baseline this root lies safely between ell and a12.  The lower
        # bracket ell*2 is also safe for ell=0.1, nu=10; for other ell/nu values,
        # widen if needed.
        lo = self.ell * mp.mpf("2")
        hi = self.a12
        # If the lower bracket is not below 1, scan upward from ell.
        if self.r_regime_II_of_a(lo) > 1:
            lo = self.ell * (1 + mp.mpf("1e-20"))
        self.a1 = _bisect_mp(lambda a: self.r_regime_II_of_a(a) - 1, lo, hi, tol=self._tol)

        # Regime III solution: a(r) = C sqrt(r) - ell/(3r), matched at r=1.
        self._C3 = self.a1 + self.ell / 3
        self.rF = (4 * self.ell / (3 * self._C3)) ** (mp.mpf(2) / 3)
        self.aF = self.ell / self.rF

    # ---------- Uniform CDF and path geometry ----------

    def F_uniform(self, x: mp.mpf) -> mp.mpf:
        if x <= self.ell:
            return mp.mpf("0")
        if x >= self.nu:
            return mp.mpf("1")
        return (x - self.ell) / self.d

    def H(self, a: mp.mpf, r: mp.mpf) -> mp.mpf:
        """H(a,r)=Pr(A<=a, E/A<=r) under iid Uniform[ell,nu]."""
        L, U, d = self.ell, self.nu, self.d
        if a <= L:
            return mp.mpf("0")
        a = min(a, U)

        # Integrate F_E(r x) against f_A(x)=1/d analytically.
        x_low_linear = max(L, L / r)
        x_high_linear = min(a, U / r)
        linear = mp.mpf("0")
        if x_high_linear > x_low_linear:
            linear = (r * (x_high_linear**2 - x_low_linear**2) / 2 - L * (x_high_linear - x_low_linear)) / d

        x_low_full = max(L, U / r)
        full = mp.mpf("0")
        if a > x_low_full:
            full = a - x_low_full

        return (linear + full) / d

    def a_regime_I(self, r: mp.mpf) -> mp.mpf:
        """Regime I path: r in [r12, r0], with r*a >= nu."""
        L, U, d, r0 = self.ell, self.nu, self.d, self.r0
        return U + ((-U**2 / r - L**2 * r) - (-U**2 / r0 - L**2 * r0)) / (2 * d)

    def J_regime_II(self, a: mp.mpf) -> mp.mpf:
        """Antiderivative of -2 ell / (a^2-ell^2)^2 for Regime II."""
        L = self.ell
        return a / (L * (a**2 - L**2)) + mp.log((a - L) / (a + L)) / (2 * L**2)

    def r_regime_II_of_a(self, a: mp.mpf) -> mp.mpf:
        """Regime II implicit solution r(a), matched at (a12,r12)."""
        return (a**2 - self.ell**2) * (self._C12 + self.J_regime_II(a) - self._J12)

    def dr_da_regime_II(self, a: mp.mpf) -> mp.mpf:
        r = self.r_regime_II_of_a(a)
        return 2 * (r * a - self.ell) / (a**2 - self.ell**2)

    def a_regime_III(self, r: mp.mpf) -> mp.mpf:
        """Regime III path: r in [rF,1], with r<1."""
        return self._C3 * mp.sqrt(r) - self.ell / (3 * r)

    def a_path(self, r: mp.mpf) -> mp.mpf:
        """Return a(r) on the full ratio-frontier path."""
        if r >= self.r12:
            return self.a_regime_I(r)
        if r >= 1:
            # Invert r_regime_II_of_a(a)=r over [a1,a12].
            return _bisect_mp(lambda a: self.r_regime_II_of_a(a) - r, self.a1, self.a12, tol=self._tol)
        return self.a_regime_III(r)

    def M_r_boundary(self, a: mp.mpf, r: mp.mpf) -> mp.mpf:
        """partial_r M_R(a,r), the ecological-value flux through the ratio boundary."""
        L, U, d = self.ell, self.nu, self.d
        b = max(L, L / r)
        c = min(a, U / r)
        if c <= b:
            return mp.mpf("0")
        return r * (c**3 - b**3) / (3 * d**2)

    # ---------- Integrals ----------

    def integral_regime_I(self) -> mp.mpf:
        """Integral of partial_r M from r12 to r0 in Regime I; closed form."""
        L, U, d = self.ell, self.nu, self.d

        def anti(r: mp.mpf) -> mp.mpf:
            return (-U**3 / r - L**3 * r**2 / 2) / (3 * d**2)

        return anti(self.r0) - anti(self.r12)

    def integral_regime_II(self) -> mp.mpf:
        """Integral over Regime II, changing variables from r to a."""
        f = lambda a: self.M_r_boundary(a, self.r_regime_II_of_a(a)) * self.dr_da_regime_II(a)
        return mp.quad(f, [self.a1, self.a12])

    def integral_regime_III(self) -> mp.mpf:
        """Integral of partial_r M over r in [rF,1]."""
        f = lambda r: self.M_r_boundary(self.a_regime_III(r), r)
        return mp.quad(f, [self.rF, 1])

    def compute(self) -> CertificateResult:
        I1 = self.integral_regime_I()
        I2 = self.integral_regime_II()
        I3 = self.integral_regime_III()
        VR = I1 + I2 + I3

        # For Q(p)=ell + (nu-ell)p, integral_0^1 p Q(p) dp = ell/2 + (nu-ell)/3.
        VE = self.ell / 2 + self.d / 3
        gap = VE - VR

        return CertificateResult(
            ell=self.ell,
            nu=self.nu,
            d=self.d,
            r0=self.r0,
            r12=self.r12,
            a12=self.a12,
            a1=self.a1,
            rF=self.rF,
            aF=self.aF,
            H12=self.H(self.a12, self.r12),
            H1=self.H(self.a1, mp.mpf("1")),
            I1=I1,
            I2=I2,
            I3=I3,
            VR=VR,
            VE=VE,
            gap=gap,
        )

    # ---------- Optional frontier table ----------

    def _r_for_target_H(self, target_H: mp.mpf) -> mp.mpf:
        """Invert H(a_path(r),r)=target_H along the frontier path."""
        if target_H <= 0:
            return self.rF
        if target_H >= 1:
            return self.r0
        return _bisect_mp(lambda r: self.H(self.a_path(r), r) - target_H, self.rF, self.r0, tol=self._tol)

    def frontier_row(self, s: mp.mpf):
        target_H = 1 - 2 * s
        r = self._r_for_target_H(target_H)
        aR = self.a_path(r)
        # ME Farmer frontier under iid identical marginals.
        p = mp.sqrt(1 - 2 * s)
        aE = self.ell + self.d * p
        return s, aE, aR, r, target_H


def compute_certificate(ell: str, nu: str, dps: int) -> CertificateResult:
    return UniformRatioFrontier(ell=ell, nu=nu, dps=dps).compute()


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Claims World iid-uniform ratio-frontier certificate-style computation")
    parser.add_argument("--ell", default="0.1", help="Lower support endpoint; default 0.1")
    parser.add_argument("--nu", default="10", help="Upper support endpoint; default 10")
    parser.add_argument("--dps", type=int, default=80, help="mpmath decimal precision; default 80")
    parser.add_argument("--N", type=int, default=100, help="Scale the asymptotic gap by N; default 100")
    parser.add_argument("--no-refinement-check", action="store_true", help="Skip dps+20 stability check")
    parser.add_argument("--frontier-table", action="store_true", help="Print ME/MX Farmer-frontier table at selected s values")
    args = parser.parse_args(list(argv) if argv is not None else None)

    cert = UniformRatioFrontier(ell=args.ell, nu=args.nu, dps=args.dps)
    res = cert.compute()

    print("\nClaims World iid-uniform ratio-frontier certificate-style computation")
    print("==================================================================")
    print(f"support:             ell={_fmt(res.ell, 20)}, nu={_fmt(res.nu, 20)}")
    print(f"precision:           {args.dps} decimal digits")
    print("\nTransition points")
    print(f"r0 = nu/ell:         {_fmt(res.r0)}")
    print(f"Regime I/II r12:     {_fmt(res.r12)}")
    print(f"Regime I/II a12:     {_fmt(res.a12)}")
    print(f"H(a12,r12):          {_fmt(res.H12)}    s12={( _fmt((1-res.H12)/2, 25) )}")
    print(f"Regime II/III a1:    {_fmt(res.a1)}")
    print(f"H(a1,1):             {_fmt(res.H1)}    s1 ={( _fmt((1-res.H1)/2, 25) )}")
    print(f"terminal rF:         {_fmt(res.rF)}")
    print(f"terminal aF=ell/rF:  {_fmt(res.aF)}")

    print("\nMX value decomposition")
    print(f"I1, Regime I:        {_fmt(res.I1)}")
    print(f"I2, Regime II:       {_fmt(res.I2)}")
    print(f"I3, Regime III:      {_fmt(res.I3)}")
    print(f"V_R(1/2):            {_fmt(res.VR)}")
    print("\nME value and gap")
    print(f"V_E(1/2):            {_fmt(res.VE)}")
    print(f"DCP = V_E - V_R:     {_fmt(res.gap)}")
    print(f"N * DCP, N={args.N}:   {_fmt(res.gap * args.N)}")

    if not args.no_refinement_check:
        hi = compute_certificate(args.ell, args.nu, args.dps + 20)
        vr_diff = abs(hi.VR - res.VR)
        gap_diff = abs(hi.gap - res.gap)
        # A conservative stability bracket, not a formal interval proof.
        err = max(10 * vr_diff, mp.mpf(10) ** (-(args.dps // 2)))
        print("\nStability check against dps+20 computation")
        print(f"V_R at dps+20:       {_fmt(hi.VR)}")
        print(f"|delta V_R|:         {_fmt(vr_diff, 12)}")
        print(f"|delta gap|:         {_fmt(gap_diff, 12)}")
        print(f"stability error used:{_fmt(err, 12)}")
        print(f"stability bracket V_R:  [{_fmt(hi.VR - err, 30)}, {_fmt(hi.VR + err, 30)}]")
        print(f"stability bracket gap:  [{_fmt(hi.gap - err, 30)}, {_fmt(hi.gap + err, 30)}]")
        if hi.gap - err > mp.mpf("0.47"):
            print("PASS: stability-bracket lower bound exceeds 0.47.")
        elif hi.gap - err > 0:
            print("PASS: stability-bracket lower bound is positive.")
        else:
            print("WARNING: stability-bracket lower bound is not positive; increase --dps.")

    # Exact consistency checks at the joining points.
    print("\nResidual checks")
    print(f"a_I(r12) - nu/r12:   {_fmt(cert.a_regime_I(res.r12) - res.nu / res.r12, 12)}")
    print(f"r_II(a12) - r12:     {_fmt(cert.r_regime_II_of_a(res.a12) - res.r12, 12)}")
    print(f"r_II(a1) - 1:        {_fmt(cert.r_regime_II_of_a(res.a1) - 1, 12)}")
    print(f"a_III(rF) - ell/rF:  {_fmt(cert.a_regime_III(res.rF) - res.ell / res.rF, 12)}")
    print(f"H(terminal):         {_fmt(cert.H(res.aF, res.rF), 12)}")

    if args.frontier_table:
        print("\nFarmer-frontier table")
        print("s        a_E(s) under ME        a_R(s) under MX        r_R(s)")
        for s_str in ["0.10", "0.20", "0.30", "0.40", "0.45"]:
            s, aE, aR, r, _H = cert.frontier_row(mp.mpf(s_str))
            print(f"{s_str:<7} {_fmt(aE, 16):>20} {_fmt(aR, 16):>20} {_fmt(r, 16):>20}")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
