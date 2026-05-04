#!/usr/bin/env python3
"""
Validated-interval wrapper for the iid-uniform Claims World ratio-frontier
certificate, with optional GitHub-friendly artifact output.

Goal: prove the loose but theorem-sufficient bound

    V_R(1/2) < 2.88,

for A,E iid Uniform[0.1,10] in equal-claims Claims World, where Green uses
R=E/A and Farmer uses A. Since V_E(1/2)=3.35 exactly, this certifies

    V_E(1/2)-V_R(1/2) > 0.47.

The script uses the three-regime analytic reduction of the ratio-frontier path.
It wraps constants in interval arithmetic and bounds the two remaining smooth
integrals by monotone left/right Riemann sums. The monotonicity used here is
analytic:

Regime II integrand:
    f_II(a)= 2 r(a) ((a^3-l^3)/(a^2-l^2)) (r(a)a-l)/(3d^2),
    r'(a)=2(r(a)a-l)/(a^2-l^2)>0,
so each factor is increasing on [a1,a12].

Regime III integrand:
    f_III(r)= r(a_III(r)^3-(l/r)^3)/(3d^2),
    a_III(r)=C sqrt(r)-l/(3r),
and f_III'(r)>0 on [rF,1].

The default subdivision count is intentionally modest because the target bound
has a large margin. Increase --subdivisions for tighter intervals.

Dependency: mpmath.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Iterable, Any

import mpmath as mp

iv = mp.iv


SCHEMA_VERSION = "claims-uniform-ratio-frontier-certificate/v1"


def _set_dps(dps: int) -> None:
    mp.mp.dps = dps
    try:
        mp.iv.dps = dps
    except Exception:
        pass


def _iv(x: mp.mpf | str) -> mp.iv.mpf:
    s = str(x) if not isinstance(x, str) else x
    return iv.mpf([s, s])


def _iv_interval(lo: mp.mpf | str, hi: mp.mpf | str) -> mp.iv.mpf:
    return iv.mpf([str(lo), str(hi)])


def _lo(x) -> mp.mpf:
    return mp.mpf(x.a)


def _hi(x) -> mp.mpf:
    return mp.mpf(x.b)


def _fmt(x, n: int = 35) -> str:
    return mp.nstr(x, n=n, strip_zeros=False)


def _s(x, n: int = 80) -> str:
    """Stable decimal string for JSON/Markdown artifacts."""
    return mp.nstr(x, n=n, strip_zeros=False)


def _interval_dict(lo: mp.mpf, hi: mp.mpf, digits: int = 80) -> dict[str, str]:
    return {"lower": _s(lo, digits), "upper": _s(hi, digits)}


@dataclass
class IntervalCertificate:
    I1: object
    I2_lower: mp.mpf
    I2_upper: mp.mpf
    I3_lower: mp.mpf
    I3_upper: mp.mpf
    VR_lower: mp.mpf
    VR_upper: mp.mpf
    VE: mp.mpf
    gap_lower: mp.mpf
    gap_upper: mp.mpf
    a1_lo: mp.mpf
    a1_hi: mp.mpf
    a12_lo: mp.mpf
    a12_hi: mp.mpf
    r12_lo: mp.mpf
    r12_hi: mp.mpf
    rF_lo: mp.mpf
    rF_hi: mp.mpf
    n: int


class UniformIntervalCertificate:
    def __init__(self, ell: str = "0.1", nu: str = "10", dps: int = 90, root_pad_exp: int = 70):
        _set_dps(dps)
        self.dps = dps
        self.ell = mp.mpf(ell)
        self.nu = mp.mpf(nu)
        if not (0 < self.ell < self.nu):
            raise ValueError("Require 0 < ell < nu.")
        self.d = self.nu - self.ell
        self.r0 = self.nu / self.ell
        self.root_pad = mp.mpf(10) ** (-root_pad_exp)

        # High-precision nominal constants for bracketing and reporting.
        self.r12_nom = (self.nu**2 - mp.sqrt(self.nu**4 - self.ell**2 * (3*self.nu**2 - 2*self.nu*self.ell))) / self.ell**2
        self.a12_nom = self.nu / self.r12_nom

        self.J_nom = lambda a: a/(self.ell*(a*a-self.ell*self.ell)) + mp.log((a-self.ell)/(a+self.ell))/(2*self.ell*self.ell)
        self.C12_nom = self.r12_nom / (self.a12_nom*self.a12_nom - self.ell*self.ell)
        self.J12_nom = self.J_nom(self.a12_nom)

        def rII_nom(a: mp.mpf) -> mp.mpf:
            return (a*a-self.ell*self.ell) * (self.C12_nom + self.J_nom(a) - self.J12_nom)

        # Bracket a1 solving rII(a1)=1.
        lo = 2*self.ell
        hi = self.a12_nom
        flo = rII_nom(lo) - 1
        fhi = rII_nom(hi) - 1
        if flo * fhi > 0:
            raise RuntimeError("Could not bracket a1 root.")
        for _ in range(4*dps):
            mid = (lo+hi)/2
            fm = rII_nom(mid)-1
            if flo*fm <= 0:
                hi = mid
                fhi = fm
            else:
                lo = mid
                flo = fm
            if hi-lo < self.root_pad:
                break
        self.a1_lo_nom = lo - self.root_pad
        self.a1_hi_nom = hi + self.root_pad

        # Interval constants.
        self.L = _iv(self.ell)
        self.U = _iv(self.nu)
        self.D = self.U - self.L
        self.r0_iv = _iv(self.r0)
        self.r12 = (self.U**2 - iv.sqrt(self.U**4 - self.L**2*(3*self.U**2 - 2*self.U*self.L))) / self.L**2
        self.a12 = self.U / self.r12
        self.a1 = _iv_interval(self.a1_lo_nom, self.a1_hi_nom)

        self.J = lambda a: a/(self.L*(a*a-self.L*self.L)) + iv.log((a-self.L)/(a+self.L))/(2*self.L*self.L)
        self.J12 = self.J(self.a12)
        self.C12 = self.r12 / (self.a12*self.a12 - self.L*self.L)
        self.C3 = self.a1 + self.L/3
        self.rF = (4*self.L/(3*self.C3)) ** (mp.mpf(2)/3)

    def rII_at(self, a_iv):
        return (a_iv*a_iv-self.L*self.L) * (self.C12 + self.J(a_iv) - self.J12)

    def fII_at_point(self, a: mp.mpf):
        a_iv = _iv(a)
        r = self.rII_at(a_iv)
        q = (a_iv**3 - self.L**3) / (a_iv**2 - self.L**2)
        return 2*r*q*(r*a_iv - self.L)/(3*self.D**2)

    def aIII_at(self, r_iv):
        return self.C3 * iv.sqrt(r_iv) - self.L/(3*r_iv)

    def fIII_at_point(self, r: mp.mpf):
        r_iv = _iv(r)
        a = self.aIII_at(r_iv)
        return r_iv * (a**3 - (self.L/r_iv)**3) / (3*self.D**2)

    def I1_interval(self):
        def anti(r):
            return (-self.U**3/r - self.L**3*r**2/2) / (3*self.D**2)
        return anti(self.r0_iv) - anti(self.r12)

    @staticmethod
    def monotone_sum_bounds(f_at_point, lo: mp.mpf, hi: mp.mpf, n: int):
        """For increasing f, return interval lower/upper Riemann bounds."""
        lo = mp.mpf(lo)
        hi = mp.mpf(hi)
        dx = (hi-lo)/n
        dx_iv = _iv(dx)
        lower = iv.mpf(["0", "0"])
        upper = iv.mpf(["0", "0"])
        for k in range(n):
            lower += f_at_point(lo + k*dx) * dx_iv
            upper += f_at_point(lo + (k+1)*dx) * dx_iv
        return _lo(lower), _hi(upper)

    def validate_a1_bracket(self) -> tuple[bool, object, object]:
        lo_val = self.rII_at(_iv(self.a1_lo_nom)) - 1
        hi_val = self.rII_at(_iv(self.a1_hi_nom)) - 1
        ok = (_hi(lo_val) < 0) and (_lo(hi_val) > 0)
        return ok, lo_val, hi_val

    def compute(self, subdivisions: int = 2000) -> IntervalCertificate:
        I1 = self.I1_interval()

        # Regime II: enlarge the domain outward; fII is increasing, so right sums upper-bound.
        I2_lo, I2_hi = self.monotone_sum_bounds(
            self.fII_at_point,
            _lo(self.a1),
            _hi(self.a12),
            subdivisions,
        )

        # Regime III: enlarge the domain outward; fIII is increasing, so right sums upper-bound.
        I3_lo, I3_hi = self.monotone_sum_bounds(
            self.fIII_at_point,
            _lo(self.rF),
            mp.mpf("1"),
            subdivisions,
        )

        VR_lo = _lo(I1) + I2_lo + I3_lo
        VR_hi = _hi(I1) + I2_hi + I3_hi
        VE = self.ell/2 + self.d/3
        gap_lo = VE - VR_hi
        gap_hi = VE - VR_lo

        return IntervalCertificate(
            I1=I1,
            I2_lower=I2_lo,
            I2_upper=I2_hi,
            I3_lower=I3_lo,
            I3_upper=I3_hi,
            VR_lower=VR_lo,
            VR_upper=VR_hi,
            VE=VE,
            gap_lower=gap_lo,
            gap_upper=gap_hi,
            a1_lo=_lo(self.a1),
            a1_hi=_hi(self.a1),
            a12_lo=_lo(self.a12),
            a12_hi=_hi(self.a12),
            r12_lo=_lo(self.r12),
            r12_hi=_hi(self.r12),
            rF_lo=_lo(self.rF),
            rF_hi=_hi(self.rF),
            n=subdivisions,
        )


def build_result_payload(
    cert: UniformIntervalCertificate,
    result: IntervalCertificate,
    root_ok: bool,
    root_lo_val: object,
    root_hi_val: object,
    target_vr: mp.mpf,
    target_gap: mp.mpf,
    digits: int = 80,
) -> dict[str, Any]:
    pass_vr = result.VR_upper < target_vr
    pass_gap = result.gap_lower > target_gap
    pass_root = bool(root_ok)
    return {
        "schema_version": SCHEMA_VERSION,
        "theorem_target": "Claims World iid-uniform baseline: V_R(1/2) < 2.88 and V_E(1/2)-V_R(1/2) > 0.47",
        "game": {
            "world": "Claims World",
            "claims": "equal claims, alpha=1/2",
            "farmer_rule": "max A",
            "green_rule": "MX ratio rule, max R=E/A",
            "dgp": "A,E iid Uniform[ell,nu]",
        },
        "parameters": {
            "ell": _s(cert.ell, digits),
            "nu": _s(cert.nu, digits),
            "d": _s(cert.d, digits),
            "precision_decimal_digits": cert.dps,
            "subdivisions_per_monotone_integral": result.n,
            "target_vr_upper": _s(target_vr, digits),
            "target_gap_lower": _s(target_gap, digits),
        },
        "transition_intervals": {
            "r12": _interval_dict(result.r12_lo, result.r12_hi, digits),
            "a12": _interval_dict(result.a12_lo, result.a12_hi, digits),
            "a1": _interval_dict(result.a1_lo, result.a1_hi, digits),
            "rF": _interval_dict(result.rF_lo, result.rF_hi, digits),
        },
        "root_checks": {
            "a1_bracket_pass": pass_root,
            "rII_a1_lo_minus_1": _interval_dict(_lo(root_lo_val), _hi(root_lo_val), digits),
            "rII_a1_hi_minus_1": _interval_dict(_lo(root_hi_val), _hi(root_hi_val), digits),
        },
        "value_bounds": {
            "I1": _interval_dict(_lo(result.I1), _hi(result.I1), digits),
            "I2": _interval_dict(result.I2_lower, result.I2_upper, digits),
            "I3": _interval_dict(result.I3_lower, result.I3_upper, digits),
            "V_R": _interval_dict(result.VR_lower, result.VR_upper, digits),
            "V_E_exact": _s(result.VE, digits),
            "gap": _interval_dict(result.gap_lower, result.gap_upper, digits),
        },
        "certificate_tests": {
            "pass_vr_upper_less_than_target": bool(pass_vr),
            "pass_gap_lower_greater_than_target": bool(pass_gap),
            "pass_all": bool(pass_root and pass_vr and pass_gap),
        },
        "notes": [
            "The interval bounds use analytic monotonicity of the Regime II and Regime III integrands and left/right Riemann sums.",
            "This certificate is deliberately loose: the theorem only requires V_R < 2.88 and gap > 0.47.",
            "No simulation or random seed is used.",
        ],
    }


def render_text(payload: dict[str, Any]) -> str:
    p = payload
    lines: list[str] = []
    lines.append("Claims World iid-uniform validated interval wrapper")
    lines.append("=================================================")
    lines.append(f"support:          ell={p['parameters']['ell']}, nu={p['parameters']['nu']}")
    lines.append(f"precision:        {p['parameters']['precision_decimal_digits']} decimal digits")
    lines.append(f"subdivisions:     {p['parameters']['subdivisions_per_monotone_integral']} per monotone integral")
    lines.append("")
    lines.append("Root/transition interval checks")
    for name in ["r12", "a12", "a1", "rF"]:
        ivd = p["transition_intervals"][name]
        lines.append(f"{name} interval:      [{ivd['lower']}, {ivd['upper']}]")
    lines.append(f"a1 bracket check:  {p['root_checks']['a1_bracket_pass']}")
    lo = p['root_checks']['rII_a1_lo_minus_1']
    hi = p['root_checks']['rII_a1_hi_minus_1']
    lines.append(f"  rII(a1_lo)-1 in [{lo['lower']}, {lo['upper']}]")
    lines.append(f"  rII(a1_hi)-1 in [{hi['lower']}, {hi['upper']}]")
    lines.append("")
    lines.append("Interval value bounds")
    for name in ["I1", "I2", "I3", "V_R"]:
        ivd = p["value_bounds"][name]
        lines.append(f"{name} interval:       [{ivd['lower']}, {ivd['upper']}]")
    lines.append(f"V_E exact:         {p['value_bounds']['V_E_exact']}")
    ivd = p["value_bounds"]["gap"]
    lines.append(f"gap interval:      [{ivd['lower']}, {ivd['upper']}]")
    lines.append("")
    lines.append("Certificate tests")
    tv = p['parameters']['target_vr_upper']
    tg = p['parameters']['target_gap_lower']
    lines.append(("PASS" if p['certificate_tests']['pass_vr_upper_less_than_target'] else "FAIL") + f": V_R upper bound < {tv}")
    lines.append(("PASS" if p['certificate_tests']['pass_gap_lower_greater_than_target'] else "FAIL") + f": gap lower bound > {tg}")
    lines.append(("PASS" if p['certificate_tests']['pass_all'] else "FAIL") + ": all certificate tests")
    return "\n".join(lines) + "\n"


def render_markdown(payload: dict[str, Any]) -> str:
    p = payload
    def iv_md(section: str, name: str) -> str:
        ivd = p[section][name]
        return f"[`{ivd['lower']}`, `{ivd['upper']}`]"

    lines: list[str] = []
    lines.append("# Claims World iid-uniform interval certificate")
    lines.append("")
    lines.append("This artifact records the validated interval wrapper for the iid-uniform Claims World ratio-frontier baseline.")
    lines.append("")
    lines.append("## Theorem target")
    lines.append("")
    lines.append(f"{p['theorem_target']}")
    lines.append("")
    lines.append("## Parameters")
    lines.append("")
    lines.append("| Parameter | Value |")
    lines.append("|---|---:|")
    for key in ["ell", "nu", "d", "precision_decimal_digits", "subdivisions_per_monotone_integral", "target_vr_upper", "target_gap_lower"]:
        lines.append(f"| `{key}` | `{p['parameters'][key]}` |")
    lines.append("")
    lines.append("## Transition intervals")
    lines.append("")
    lines.append("| Quantity | Certified interval |")
    lines.append("|---|---:|")
    for name in ["r12", "a12", "a1", "rF"]:
        lines.append(f"| `{name}` | {iv_md('transition_intervals', name)} |")
    lines.append("")
    lines.append("## Value bounds")
    lines.append("")
    lines.append("| Quantity | Certified interval/value |")
    lines.append("|---|---:|")
    for name in ["I1", "I2", "I3", "V_R"]:
        lines.append(f"| `{name}` | {iv_md('value_bounds', name)} |")
    lines.append(f"| `V_E_exact` | `{p['value_bounds']['V_E_exact']}` |")
    lines.append(f"| `gap = V_E - V_R` | {iv_md('value_bounds', 'gap')} |")
    lines.append("")
    lines.append("## Root checks")
    lines.append("")
    lines.append("| Check | Result |")
    lines.append("|---|---:|")
    lines.append(f"| `a1_bracket_pass` | `{p['root_checks']['a1_bracket_pass']}` |")
    lines.append(f"| `rII(a1_lo)-1` | {iv_md('root_checks', 'rII_a1_lo_minus_1')} |")
    lines.append(f"| `rII(a1_hi)-1` | {iv_md('root_checks', 'rII_a1_hi_minus_1')} |")
    lines.append("")
    lines.append("## Certificate tests")
    lines.append("")
    lines.append("| Test | Result |")
    lines.append("|---|---:|")
    lines.append(f"| `V_R upper bound < target` | `{p['certificate_tests']['pass_vr_upper_less_than_target']}` |")
    lines.append(f"| `gap lower bound > target` | `{p['certificate_tests']['pass_gap_lower_greater_than_target']}` |")
    lines.append(f"| `pass_all` | `{p['certificate_tests']['pass_all']}` |")
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    for note in p["notes"]:
        lines.append(f"- {note}")
    lines.append("")
    return "\n".join(lines)


def write_artifacts(payload: dict[str, Any], outdir: Path, stem: str) -> tuple[Path, Path, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    json_path = outdir / f"{stem}.json"
    md_path = outdir / f"{stem}.md"
    txt_path = outdir / f"{stem}.txt"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md_path.write_text(render_markdown(payload), encoding="utf-8")
    txt_path.write_text(render_text(payload), encoding="utf-8")
    return json_path, md_path, txt_path


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Validated interval wrapper for iid-uniform Claims ratio-frontier certificate")
    parser.add_argument("--ell", default="0.1")
    parser.add_argument("--nu", default="10")
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--subdivisions", type=int, default=2000)
    parser.add_argument("--target-vr-upper", default="2.88")
    parser.add_argument("--target-gap-lower", default="0.47")
    parser.add_argument("--outdir", default=None, help="Optional directory in which to write JSON/Markdown/text certificate artifacts.")
    parser.add_argument("--artifact-stem", default="claims_uniform_ratio_frontier_interval_certificate", help="Filename stem for artifacts written with --outdir.")
    parser.add_argument("--json", default=None, help="Optional explicit JSON output path.")
    parser.add_argument("--md", default=None, help="Optional explicit Markdown output path.")
    parser.add_argument("--txt", default=None, help="Optional explicit text output path.")
    args = parser.parse_args(list(argv) if argv is not None else None)

    cert = UniformIntervalCertificate(args.ell, args.nu, args.dps)
    result = cert.compute(args.subdivisions)
    ok_root, root_lo_val, root_hi_val = cert.validate_a1_bracket()

    target_vr = mp.mpf(args.target_vr_upper)
    target_gap = mp.mpf(args.target_gap_lower)
    payload = build_result_payload(cert, result, ok_root, root_lo_val, root_hi_val, target_vr, target_gap)

    print(render_text(payload))

    written: list[Path] = []
    if args.outdir:
        written.extend(write_artifacts(payload, Path(args.outdir), args.artifact_stem))
    if args.json:
        path = Path(args.json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        written.append(path)
    if args.md:
        path = Path(args.md)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(render_markdown(payload), encoding="utf-8")
        written.append(path)
    if args.txt:
        path = Path(args.txt)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(render_text(payload), encoding="utf-8")
        written.append(path)

    if written:
        print("Artifacts written:")
        for path in written:
            print(f"  {path}")
        print()

    return 0 if payload["certificate_tests"]["pass_all"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
