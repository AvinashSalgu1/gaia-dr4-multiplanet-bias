"""The analyses of step 2, each applied to one data set.

They follow the astrometric orbit models that Gaia DR4 itself publishes
(draft DR4 data model, June 2026: tables nss_two_body_orbit and
nss_multiple_orbits):

    one    one planet               ~ DR4 "Orbital"                      (12 parameters)
           gaiamock's orbit fit, Lammers & Winn detection criteria
    accel  one planet + acceleration ~ DR4 "OrbitalAlternativeAcceleration" (14)
           kept if it improves chi2 by > 2 ln N (BIC), else the one-planet fit
    jerk   one planet + jerk         ~ DR4 "OrbitalAlternativeJerk"         (16)
           kept if it improves chi2 by > 4 ln N, else the one-planet fit
    two    two planets               ~ DR4 "OrbitalAlternativeMultiple"     (19)
           one planet, a search for a second signal with planet 1 held fixed,
           then a joint fit; the second planet is used only if it is itself
           detected ("if justified")
    best   DR4-like model choice (no new fits): of no planet (5), acceleration
           only (7), jerk only (9), one planet (12), + acceleration (14),
           + jerk (16) and two planets (19), keep the lowest BIC = chi2 + k ln N;
           its planets are reported if they pass the usual detection test.
           How DR4 really chooses between its models is not public yet (their
           papers are not out); BIC is our stand-in.

Optional (not run by default): dr3  the DR3-style cascade of gaiamock
(RUWE < 1.4 -> plain star; else acceleration solutions; only then an orbit).

Every planet also gets the Gaia NSS "significance" a0 / sigma_a0 (DR3
definition): DR4 publishes orbits with significance < 5 separately as
"OrbitalPoorlyConstrained".

Every analysis is run on the data WITH the companion (version "w") and on the
same data WITHOUT it (version "a", target alone, identical noise). Results are
stored in flat columns named <analysis>_<version>_<quantity>.
"""
from __future__ import annotations

import numpy as np

FITTED = ["one", "accel", "jerk", "two"]           # analyses with their own fits (default)
ANALYSES = FITTED + ["best"]                         # "best" is chosen from the fits above
VERSIONS = ["w", "a"]
LABEL = {"one": "one planet", "accel": "one planet + acceleration", "jerk": "one planet + jerk",
         "two": "two planets (if justified)", "best": "best model by BIC (DR4-like)",
         "dr3": "DR3-style pipeline"}
DR4_NAME = {"one": "Orbital", "accel": "OrbitalAlternativeAcceleration", "jerk": "OrbitalAlternativeJerk",
            "two": "OrbitalAlternativeMultiple", "best": "lowest BIC of all models"}
SIG_MIN = 5.0                                        # below: DR4 "OrbitalPoorlyConstrained"

# quantities per analysis; those in REQUIRED must be present for a stored result to count
FIELDS = {
    "one": ["det", "P", "m", "e", "dchi2", "F2", "sigP", "a0", "sa0"],
    "accel": ["det", "P", "m", "dchi2", "acc", "sig",
              "o_det", "o_P", "o_m", "o_dchi2", "o_sig",
              "a_det", "a_P", "a_m", "a_sig", "c5", "c7", "c0", "c2"],
    "jerk": ["det", "P", "m", "dchi2", "acc", "sig", "j_det", "j_P", "j_m", "j_sig", "c9", "cJ"],
    "two": ["gain", "P1", "m1", "d1", "s1", "g1", "P2", "m2", "d2", "s2", "g2", "chi2"],
    "dr3": ["kind", "det", "P", "m", "dchi2"],
}
REQUIRED = {
    "one": ["det", "P", "m", "a0", "sa0"],
    "accel": ["det", "P", "m", "acc", "sig", "o_det", "o_P", "o_m", "o_sig", "a_det", "a_sig",
              "c5", "c7", "c0", "c2"],
    "jerk": ["det", "P", "m", "acc", "sig", "j_det", "j_sig", "c9", "cJ"],
    "two": ["P1", "m1", "d1", "s1", "g1", "P2", "m2", "d2", "s2", "g2", "chi2"],
    "dr3": ["kind", "det", "P", "m"],
}
TEXT_FIELDS = {"dr3_w_kind", "dr3_a_kind"}


def cols(analysis: str, version: str, required_only: bool = False) -> list[str]:
    f = REQUIRED[analysis] if required_only else FIELDS[analysis]
    return [f"{analysis}_{version}_{q}" for q in f]


def all_cols(fitted=FITTED) -> list[str]:
    return [c for a in fitted for v in VERSIONS for c in cols(a, v)]


def _nan(analysis, version):
    out = {c: np.nan for c in cols(analysis, version)}
    if analysis == "dr3":
        out[f"dr3_{version}_kind"] = "failed"
    for k in ("det", "o_det", "a_det", "j_det"):
        if f"{analysis}_{version}_{k}" in out:
            out[f"{analysis}_{version}_{k}"] = 0
    return out


def _trend_choice(f0, ft, extra_params, n):
    """Keep the orbit + trend fit if it improves chi2 by more than extra_params * ln N (BIC)."""
    use = (f0.chi2 - ft.chi2) > extra_params * np.log(n)
    return (ft if use else f0), use


def run(ep, obs, err, m_star: float, version: str, which: set[str]) -> dict:
    """Run the analyses in `which` on one data set; returns flat columns."""
    from ..fitting import (chi2_at, detected_lw, f2_statistic, fit_5par, fit_one_companion,
                           fit_orbit_trend, fit_two_planets, lammers_winn_detected, run_cascade)
    from ..orbits import companion_mass_from_a0

    v = version
    n = len(ep)
    out = {}
    det = lambda f: int(detected_lw(f.delta_chi2, f.m_mjup, f.P_days, f.sigma_P))  # noqa: E731
    if "one" in which:
        try:
            f = fit_one_companion(ep, obs, err, m_star)
            out.update({f"one_{v}_det": int(lammers_winn_detected(f)), f"one_{v}_P": f.P_days,
                        f"one_{v}_m": f.m_mjup, f"one_{v}_e": f.ecc, f"one_{v}_dchi2": f.delta_chi2,
                        f"one_{v}_F2": f2_statistic(f.chi2, n), f"one_{v}_sigP": f.sigma_P,
                        f"one_{v}_a0": f.a0_mas, f"one_{v}_sa0": f.sigma_a0})
        except Exception:
            out.update(_nan("one", v))

    f0 = None
    if which & {"accel", "jerk", "two"}:
        try:
            f0 = fit_orbit_trend(ep, obs, err, m_star, n_trend=0)
        except Exception:
            f0 = None
    if "accel" in which:
        try:
            f2 = fit_orbit_trend(ep, obs, err, m_star, n_trend=2)
            fe, use = _trend_choice(f0, f2, 2, n)
            out.update({f"accel_{v}_det": det(fe), f"accel_{v}_P": fe.P_days, f"accel_{v}_m": fe.m_mjup,
                        f"accel_{v}_dchi2": fe.delta_chi2, f"accel_{v}_acc": int(use),
                        f"accel_{v}_sig": fe.significance,
                        f"accel_{v}_o_det": det(f0), f"accel_{v}_o_P": f0.P_days, f"accel_{v}_o_m": f0.m_mjup,
                        f"accel_{v}_o_dchi2": f0.delta_chi2, f"accel_{v}_o_sig": f0.significance,
                        # one planet + acceleration whether or not the BIC rule keeps it
                        f"accel_{v}_a_det": det(f2), f"accel_{v}_a_P": f2.P_days, f"accel_{v}_a_m": f2.m_mjup,
                        f"accel_{v}_a_sig": f2.significance,
                        # chi2 of no planet / acceleration only / one planet / one planet + accel (for BIC)
                        f"accel_{v}_c5": f0.chi2_null, f"accel_{v}_c7": f2.chi2_null,
                        f"accel_{v}_c0": f0.chi2, f"accel_{v}_c2": f2.chi2})
        except Exception:
            out.update(_nan("accel", v))
    if "jerk" in which:
        try:
            f4 = fit_orbit_trend(ep, obs, err, m_star, n_trend=4)
            fe, use = _trend_choice(f0, f4, 4, n)
            out.update({f"jerk_{v}_det": det(fe), f"jerk_{v}_P": fe.P_days, f"jerk_{v}_m": fe.m_mjup,
                        f"jerk_{v}_dchi2": fe.delta_chi2, f"jerk_{v}_acc": int(use), f"jerk_{v}_sig": fe.significance,
                        f"jerk_{v}_j_det": det(f4), f"jerk_{v}_j_P": f4.P_days, f"jerk_{v}_j_m": f4.m_mjup,
                        f"jerk_{v}_j_sig": f4.significance, f"jerk_{v}_c9": f4.chi2_null, f"jerk_{v}_cJ": f4.chi2})
        except Exception:
            out.update(_nan("jerk", v))
    if "two" in which:
        try:
            t = fit_two_planets(ep, obs, err, m_star, one=f0)
            g = t.significance
            out.update({f"two_{v}_gain": t.gain_over_one,
                        f"two_{v}_P1": t.P[0], f"two_{v}_m1": t.m_mjup[0], f"two_{v}_d1": t.drop_dchi2[0],
                        f"two_{v}_s1": t.sigma_P[0], f"two_{v}_g1": g[0],
                        f"two_{v}_P2": t.P[1], f"two_{v}_m2": t.m_mjup[1], f"two_{v}_d2": t.drop_dchi2[1],
                        f"two_{v}_s2": t.sigma_P[1], f"two_{v}_g2": g[1], f"two_{v}_chi2": t.chi2})
        except Exception:
            out.update(_nan("two", v))
    if "dr3" in which:
        try:
            r = run_cascade(ep, obs, err)
            out.update({f"dr3_{v}_kind": r["kind"], f"dr3_{v}_det": 0, f"dr3_{v}_P": np.nan,
                        f"dr3_{v}_m": np.nan, f"dr3_{v}_dchi2": np.nan})
            if r["kind"] == "orbit" and r["parallax"] > 0:
                chi2_o, _ = chi2_at(ep, obs, err, r["P"], r["phi_p"], r["ecc"])
                dchi2 = fit_5par(ep, obs, err)[0] - chi2_o
                m = companion_mass_from_a0(r["a0"], r["P"], m_star, r["parallax"])
                out.update({f"dr3_{v}_det": int(detected_lw(dchi2, m, r["P"], r["sigma_P"])),
                            f"dr3_{v}_P": r["P"], f"dr3_{v}_m": m, f"dr3_{v}_dchi2": dchi2})
        except Exception:
            out.update(_nan("dr3", v))
    return out
