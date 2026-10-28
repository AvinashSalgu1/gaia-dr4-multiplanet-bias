"""Which planets each analysis reports, and whether they are right.

Rules (same as the original step 2):
* A reported planet is "locked" on the target (or companion) if its period is
  within LOCK_TOL in log (about 10%) of that planet's true period.
* TARGET FOUND CORRECTLY: a reported planet locked on the target with period
  within 10% and mass within 20% of a REFERENCE. The reference is what the same
  analysis reports for the target-alone data (identical noise), if it reported
  a planet there; otherwise the truth. This removes the near-threshold noise
  bias that affects both versions equally, so only the companion's effect is
  measured.
* COMPANION FOUND CORRECTLY: locked on the companion, within 10% / 20% of its
  true period / mass (a real planet, correctly catalogued).
* WRONG: any reported planet that is neither of the above.
* WRONG BUT CONFIDENT: a wrong planet with significance a0/sigma_a0 >= 5,
  i.e. one DR4 would publish as a normal orbit rather than as
  "OrbitalPoorlyConstrained".
* FALSE ALARM: on the target-alone data, a reported planet that is not the
  target (period off by more than ~10%), or a second planet.

Best model by BIC: of no planet (5 parameters), acceleration only (7), jerk
only (9), one planet (12), one planet + acceleration (14), one planet + jerk
(16) and two planets (19), the model with the lowest BIC = chi2 + k ln N
(N = number of transits) is chosen; its
planets are then reported only if they pass the usual detection test (one /
one + acceleration: Lammers & Winn criteria; two planets: each planet as below).

Two-planet analysis, "if justified": the joint two-planet solution is used
only when its second planet is itself detected (removing it raises chi2 by
> 50, m < 13 M_J, (P+sigma_P)/(P-sigma_P) < 1.5); otherwise the one-planet
answer (our numpy fit, with Lammers & Winn's criteria) is reported.
"""
from __future__ import annotations

import numpy as np

from .analyses import ANALYSES, SIG_MIN
from .sim import LOCK_TOL

BIC_MODELS = ["no planet", "acceleration only", "jerk only", "one planet", "one planet + acceleration",
              "one planet + jerk", "two planets"]
BIC_K = np.array([5, 7, 9, 12, 14, 16, 19])
NO_PLANET_MODELS = (0, 1, 2)


def _get(d, key, n):
    v = d.get(key)
    return np.full(n, np.nan) if v is None else np.asarray(v, float)


def bic_choice(d: dict, v: str) -> np.ndarray:
    """Index into BIC_MODELS of the lowest-BIC model, per system (-1 if no fits)."""
    n = len(d["n_transits"])
    chi2 = np.column_stack([_get(d, f"accel_{v}_c5", n), _get(d, f"accel_{v}_c7", n), _get(d, f"jerk_{v}_c9", n),
                            _get(d, f"accel_{v}_c0", n), _get(d, f"accel_{v}_c2", n), _get(d, f"jerk_{v}_cJ", n),
                            _get(d, f"two_{v}_chi2", n)])
    bic = chi2 + BIC_K[None, :] * np.log(np.asarray(d["n_transits"], float))[:, None]
    bic = np.where(np.isfinite(bic), bic, np.inf)
    out = np.argmin(bic, axis=1)
    out[~np.isfinite(bic).any(axis=1)] = -1
    return out


def _passes(P, m, dchi2, sigP):
    with np.errstate(invalid="ignore", divide="ignore"):
        lo = P - sigP
        ok = (dchi2 > 50) & (m < 13) & (m > 0) & (lo > 0) & ((P + sigP) / lo < 1.5)
    return np.nan_to_num(ok, nan=0).astype(bool)


def reported(d: dict, analysis: str, v: str) -> list[tuple]:
    """List of (reported mask, P, m, significance a0/sigma_a0) for one analysis and version."""
    n = len(d["P_t"]) if "P_t" in d else len(d["n_transits"])
    g = lambda q: _get(d, f"{analysis}_{v}_{q}", n)  # noqa: E731
    if analysis == "one":
        with np.errstate(invalid="ignore", divide="ignore"):
            sig = g("a0") / g("sa0")
        return [(g("det") == 1, g("P"), g("m"), sig)]
    if analysis == "dr3":
        return [(g("det") == 1, g("P"), g("m"), np.full(n, np.nan))]
    if analysis in ("accel", "jerk"):
        return [(g("det") == 1, g("P"), g("m"), g("sig"))]
    a = lambda q: _get(d, f"accel_{v}_{q}", n)  # noqa: E731
    t = lambda q: _get(d, f"two_{v}_{q}", n)  # noqa: E731
    j = lambda q: _get(d, f"jerk_{v}_{q}", n)  # noqa: E731
    r1 = _passes(t("P1"), t("m1"), t("d1"), t("s1"))
    r2 = _passes(t("P2"), t("m2"), t("d2"), t("s2"))
    if analysis == "best":
        ch = bic_choice(d, v)
        c = [ch == 3, ch == 4, ch == 5, ch == 6]
        rep = np.select(c, [a("o_det") == 1, a("a_det") == 1, j("j_det") == 1, r1], False)
        P = np.select(c, [a("o_P"), a("a_P"), j("j_P"), t("P1")], np.nan)
        m = np.select(c, [a("o_m"), a("a_m"), j("j_m"), t("m1")], np.nan)
        sig = np.select(c, [a("o_sig"), a("a_sig"), j("j_sig"), t("g1")], np.nan)
        return [(rep, P, m, sig), ((ch == 6) & r2, t("P2"), t("m2"), t("g2"))]
    if analysis == "two":
        first = (np.where(r2, r1, a("o_det") == 1), np.where(r2, t("P1"), a("o_P")),
                 np.where(r2, t("m1"), a("o_m")), np.where(r2, t("g1"), a("o_sig")))
        return [first, (r2, t("P2"), t("m2"), t("g2"))]
    raise ValueError(analysis)


def _reference(d, alone):
    """Per system: (P, m) of the target-alone answer nearest the target, else truth."""
    P_ref, m_ref = d["P_t"].copy(), d["m_t"].copy()
    best = np.full(len(P_ref), np.inf)
    with np.errstate(invalid="ignore", divide="ignore"):
        for rep, P, m, _ in alone:
            dist = np.where(rep, np.abs(np.log(P / d["P_t"])), np.inf)
            take = dist < best
            P_ref[take], m_ref[take], best[take] = P[take], m[take], dist[take]
    return P_ref, m_ref


def judge(d: dict, analysis: str, alone: list | None = None) -> dict:
    """Outcome flags for one analysis (arrays over systems).

    `alone` overrides what is taken as this analysis' answer on the target-alone
    data (the reference); by default the same analysis on that data.
    """
    w = reported(d, analysis, "w")
    a = reported(d, analysis, "a") if alone is None else alone
    P_ref, m_ref = _reference(d, a)
    n = len(d["P_t"])
    ok_t, ok_c = np.zeros(n, bool), np.zeros(n, bool)
    n_rep, n_ok = np.zeros(n, int), np.zeros(n, int)
    wrong_conf, wrong_weak = np.zeros(n, bool), np.zeros(n, bool)
    with np.errstate(invalid="ignore", divide="ignore"):
        for rep, P, m, sig in w:
            dt, dc = np.abs(np.log(P / d["P_t"])), np.abs(np.log(P / d["P_c"]))
            lock_t = (dt < LOCK_TOL) & (dt <= dc)
            lock_c = (dc < LOCK_TOL) & ~lock_t
            t_ok = rep & lock_t & (np.abs(m / m_ref - 1) <= 0.2) & (np.abs(P / P_ref - 1) <= 0.1)
            c_ok = rep & lock_c & (np.abs(m / d["m_c"] - 1) <= 0.2) & (np.abs(P / d["P_c"] - 1) <= 0.1)
            bad = rep & ~(t_ok | c_ok)
            wrong_conf |= bad & (sig >= SIG_MIN)
            wrong_weak |= bad & ~(sig >= SIG_MIN)
            ok_t |= t_ok
            ok_c |= c_ok
            n_rep += rep
            n_ok += t_ok | c_ok
        # target-alone data: is anything reported that is not the target?
        alone_rep = np.zeros(n, int)
        alone_bad = np.zeros(n, bool)
        alone_found = np.zeros(n, bool)
        for rep, P, m, _ in a:
            alone_rep += rep
            near = np.abs(np.log(P / d["P_t"])) < LOCK_TOL
            alone_bad |= rep & ~near
            alone_found |= rep & near
        alone_bad |= alone_rep >= 2
    wrong = (n_rep - n_ok) > 0
    accel_only = np.zeros(n, bool)
    if analysis == "dr3":
        accel_only = np.isin(d["dr3_w_kind"], ["accel7", "accel9"]) & (n_rep == 0)
    if analysis == "best":
        accel_only = np.isin(bic_choice(d, "w"), [1, 2]) & (n_rep == 0)
    return dict(target=ok_t, companion=ok_c, wrong=wrong, wrong_confident=wrong_conf,
                wrong_poorly_constrained=wrong_weak & ~wrong_conf, n_reported=n_rep, accel_only=accel_only,
                none=(n_rep == 0) & ~accel_only, alone_detected=alone_found, false_alarm=alone_bad)


def judge_all(d: dict, analyses=ANALYSES) -> dict:
    return {a: judge(d, a) for a in analyses}


def exclusive_categories(o: dict) -> dict:
    """Mutually exclusive outcome shares for stacked bars."""
    t, c, w = o["target"], o["companion"], o["wrong"]
    return {
        "both planets correct": t & c & ~w,
        "target correct": t & ~c & ~w,
        "companion correct (target missed)": c & ~t & ~w,
        "a correct planet + a wrong one": (t | c) & w,
        "WRONG planet only": w & ~(t | c),
        "acceleration / jerk only": o["accel_only"],
        "nothing reported": o["none"],
    }
