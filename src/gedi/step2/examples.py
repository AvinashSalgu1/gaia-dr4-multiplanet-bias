"""Posterior samples and corner plots for three example stars.

Same star and Gaia epochs as the 'outer' sweep reference system:
    no planet   : the star alone
    one planet  : the target alone (3 M_J at 1.5 yr, SNR1 = 4)
    two planets : the target + a companion at 15 yr with 3x the target's wobble

Each star is analysed with the one-planet, one-planet + acceleration and
two-planet models; the title gives the DR4-like choice (lowest BIC of all models).

Sampling: an affine-invariant ensemble sampler (the "stretch move" of Goodman
& Weare 2010, the algorithm used by emcee), written out here so that no extra
package is needed. Only the orbit-shape parameters are sampled (log P, phase,
e for each planet); the parameters that enter linearly (position, proper
motion, parallax, Thiele-Innes constants) are integrated out analytically
(flat priors): log L = -chi2_min / 2 - log det(M^T W M) / 2. Mass for each
sample uses the best-fitting linear parameters at that orbit shape.
Walkers start in a small ball around the best fit, so the plots show the
posterior around the solution the fit found (a second, distant solution
would not be explored).
"""
from __future__ import annotations

import numpy as np

from .sim import RESULTS

LN_P = (np.log(10.0), np.log(1e4))


# ----------------------------------------------------------------------------- sampler
def stretch_sampler(log_prob, p0, n_steps, rng, a=2.0):
    """Goodman & Weare affine-invariant ensemble sampler. Returns (chain[steps, walkers, dim], acc)."""
    p = np.array(p0, float)
    n_w, dim = p.shape
    lp = np.array([log_prob(x) for x in p])
    chain = np.empty((n_steps, n_w, dim))
    acc = 0
    half = n_w // 2
    for s in range(n_steps):
        for first, second in ((slice(0, half), slice(half, n_w)), (slice(half, n_w), slice(0, half))):
            S, C = p[first], p[second]
            z = ((a - 1) * rng.uniform(size=len(S)) + 1) ** 2 / a
            partners = C[rng.integers(0, len(C), len(S))]
            prop = partners + z[:, None] * (S - partners)
            lp_new = np.array([log_prob(x) for x in prop])
            log_r = (dim - 1) * np.log(z) + lp_new - lp[first]
            ok = np.log(rng.uniform(size=len(S))) < log_r
            idx = np.arange(n_w)[first][ok]
            p[idx], lp[idx] = prop[ok], lp_new[ok]
            acc += ok.sum()
        chain[s] = p
    return chain, acc / (n_steps * n_w)


# ----------------------------------------------------------------------------- models
def _design(model, ep, x):
    from ..fitting import _design_9par, _design_orbit_trend, _design_two
    if model == "one":
        return _design_9par(ep, np.exp(x[0]), x[1] % (2 * np.pi), x[2])
    if model == "accel":
        return _design_orbit_trend(ep, np.exp(x[0]), x[1] % (2 * np.pi), x[2], 2)
    return _design_two(ep, np.exp(x[0]), x[1] % (2 * np.pi), x[2], np.exp(x[3]), x[4] % (2 * np.pi), x[5])


def make_log_prob(model, ep, obs, err):
    w = 1.0 / err
    yw = obs * w
    n_orb = 2 if model == "two" else 1

    def log_prob(x):
        for k in range(n_orb):
            lnP, e = x[3 * k], x[3 * k + 2]
            if not (LN_P[0] < lnP < LN_P[1]) or not (0.0 <= e < 0.99):
                return -np.inf
        Mw = _design(model, ep, x) * w[:, None]
        A = Mw.T @ Mw
        try:
            mu = np.linalg.solve(A, Mw.T @ yw)
        except np.linalg.LinAlgError:
            return -np.inf
        r = yw - Mw @ mu
        sign, logdet = np.linalg.slogdet(A)
        return -0.5 * float(r @ r) - 0.5 * logdet if sign > 0 else -np.inf

    return log_prob


def derived(model, ep, obs, err, x, m_star):
    """Physical parameters for one sample: [P1 yr, m1, e1, (P2 yr, m2, e2)]."""
    from ..fitting import _weighted_lstsq
    from ..orbits import campbell_from_thiele_innes, companion_mass_from_a0
    _, mu = _weighted_lstsq(_design(model, ep, x), obs, err)
    out = []
    for k, c0 in enumerate((5,) if model != "two" else (5, 9)):
        B, G, A_, F = mu[c0:c0 + 4]
        a0, *_ = campbell_from_thiele_innes(A_, B, F, G)
        P = np.exp(x[3 * k])
        m = companion_mass_from_a0(a0, P, m_star, mu[4]) if mu[4] > 0 else np.nan
        out += [P / 365.25, m, x[3 * k + 2]]
    return out


def sample(model, ep, obs, err, m_star, best_x, rng, n_steps=None, n_keep=3000):
    """Two-planet model (6 dimensions, often multi-modal): more walkers, 5x longer chains and a
    smaller stretch scale (a = 1.6) so that more proposals are accepted and the chains mix."""
    dim = len(best_x)
    n_w = 32 if dim == 3 else 64
    n_steps = n_steps or (2500 if dim == 3 else 20000)
    scale = np.tile([1e-3, 1e-2, 1e-2], dim // 3)
    p0 = best_x + scale * rng.standard_normal((n_w, dim))
    p0[:, 2::3] = np.clip(p0[:, 2::3], 1e-4, 0.98)
    chain, acc = stretch_sampler(make_log_prob(model, ep, obs, err), p0, n_steps, rng, a=2.0 if dim == 3 else 1.6)
    flat = chain[n_steps // 2:].reshape(-1, dim)
    flat = flat[rng.choice(len(flat), min(n_keep, len(flat)), replace=False)]
    phys = np.array([derived(model, ep, obs, err, x, m_star) for x in flat])
    return phys, acc


# ----------------------------------------------------------------------------- the three stars
def best_fits(ep, obs, err, m_star):
    from ..fitting import fit_orbit_trend, fit_two_planets
    f0 = fit_orbit_trend(ep, obs, err, m_star, n_trend=0)
    f2 = fit_orbit_trend(ep, obs, err, m_star, n_trend=2)
    t = fit_two_planets(ep, obs, err, m_star, one=f0)
    return {"one": (np.array([np.log(f0.P_days), f0.phi_p, f0.ecc]), f0),
            "accel": (np.array([np.log(f2.P_days), f2.phi_p, f2.ecc]), f2),
            "two": (np.array([np.log(t.P[0]), t.phi_p[0], t.ecc[0], np.log(t.P[1]), t.phi_p[1], t.ecc[1]]), t)}


def run(seed=5, quick=False):
    """Sample all models for the three stars; saves results/step2_examples.npz and plots."""
    from ..fitting import fit_orbit_trend
    from .outcomes import BIC_K, BIC_MODELS
    from .sweeps import BASELINES, build, reference_epochs, reference_parallax, simulate_pair

    rng = np.random.default_rng(seed)
    p = BASELINES["outer"]
    _, _, ep = reference_epochs()
    t, c, plx = build(p, reference_parallax(p))
    obs_t, obs_tc, err = simulate_pair(p, plx, t, c, ep, rng)
    from ..simulate import Star, simulate
    star = Star(parallax_mas=plx, m_star_msun=p["m_star"], G_mag=p["G"], pmra=20.0, pmdec=-10.0)
    obs_0 = simulate(ep, star, [], add_noise=False).obs + (obs_t - simulate(ep, star, [t], add_noise=False).obs)
    stars = {"no planet": obs_0, "one planet": obs_t, "two planets": obs_tc}
    truth = dict(P_t=t.P_days / 365.25, m_t=t.m_mjup, e_t=t.ecc, P_c=c.P_days / 365.25, m_c=c.m_mjup, e_c=c.ecc)
    save = {f"truth_{k}": v for k, v in truth.items()}
    for name, obs in stars.items():
        fits = best_fits(ep, obs, err, p["m_star"])
        f0, f2, tw = fits["one"][1], fits["accel"][1], fits["two"][1]
        f4 = fit_orbit_trend(ep, obs, err, p["m_star"], n_trend=4)
        chi2 = np.array([f0.chi2_null, f2.chi2_null, f4.chi2_null, f0.chi2, f2.chi2, f4.chi2, tw.chi2])
        save[f"{name}|bic_choice"] = BIC_MODELS[int(np.argmin(chi2 + BIC_K * np.log(len(ep))))]
        for model, (x, fit) in fits.items():
            steps = 300 if quick else None
            phys, acc = sample(model, ep, obs, err, p["m_star"], x, rng, n_steps=steps,
                               n_keep=400 if quick else 3000)
            save[f"{name}|{model}"] = phys
            gain = fit.chi2_null - fit.chi2 if model != "two" else fit.gain_over_one
            save[f"{name}|{model}|dchi2"] = gain
            print(f"  {name:<12} {model:<6} acceptance {acc:.2f}  "
                  f"{'delta chi2 vs no planet' if model != 'two' else 'gain over one planet'} {gain:.1f}", flush=True)
    out = RESULTS / "step2_examples.npz"
    RESULTS.mkdir(exist_ok=True)
    np.savez(out, **save)
    print(f"saved {out}")
    plot()


# ----------------------------------------------------------------------------- corner plots
def _log_ticks(axis, lo, hi, max_ticks=3):
    """At most 3 readable ticks on a log axis: decades (1, 10, 100, ...) when the range is wide,
    1-2-5 values when it is narrow, so that labels never overlap in the small corner panels."""
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator
    span = np.log10(hi / lo)
    if span > 1.3:
        dec = [10.0 ** k for k in range(int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi))) + 1)]
        t = [v for v in dec if lo <= v <= hi]
        if len(t) < 2:                                     # e.g. 2-200: add 3x values
            t = sorted(set(t + [v for v in (3 * d for d in dec) if lo <= v <= hi]))
    else:
        nice = [m * 10.0 ** k for k in range(-3, 5) for m in (1, 2, 5)]
        t = [v for v in nice if lo <= v <= hi]
        if len(t) < 2:
            t = [m * 10.0 ** k for k in range(-3, 5) for m in (1, 1.5, 2, 3, 5, 7)]
            t = [v for v in t if lo <= v <= hi]
    t = t[::max(1, int(np.ceil(len(t) / max_ticks)))]

    def fmt(v, _):
        if v >= 1e4:
            return f"1e{int(round(np.log10(v)))}"
        return f"{v:g}"

    axis.set_major_locator(FixedLocator(t))
    axis.set_major_formatter(FuncFormatter(fmt))
    axis.set_minor_locator(NullLocator())
    axis.set_minor_formatter(NullFormatter())


def _corner(subfig, data, labels, logs, truths, color):
    """Minimal corner plot into a matplotlib subfigure. truths: list per param of [(value, style)]."""
    from .. import plotstyle as ps
    K = data.shape[1]
    axs = subfig.subplots(K, K, squeeze=False)
    lims = []
    for k in range(K):
        x = data[:, k]
        x = x[np.isfinite(x) & ((x > 0) if logs[k] else True)]
        lo, hi = np.percentile(x, [0.5, 99.5]) if len(x) else (0.1, 1.0)
        for v, _ in truths[k]:
            if np.isfinite(v) and (v > 0 or not logs[k]):
                lo, hi = min(lo, v), max(hi, v)
        if logs[k]:
            lo, hi = lo / 1.15, hi * 1.15
        else:
            pad = 0.05 * (hi - lo or 1)
            lo, hi = lo - pad, hi + pad
        lims.append((lo, hi))
    bins = [np.geomspace(*lims[k], 35) if logs[k] else np.linspace(*lims[k], 35) for k in range(K)]
    for i in range(K):
        for j in range(K):
            a = axs[i, j]
            if j > i:
                a.axis("off")
                continue
            if i == j:
                a.hist(data[:, i], bins=bins[i], color=color, histtype="stepfilled", alpha=0.75)
                a.set_yticks([])
                a.spines["left"].set_visible(False)
            else:
                h, xe, ye = np.histogram2d(data[:, j], data[:, i], bins=[bins[j], bins[i]])
                a.pcolormesh(xe, ye, h.T, cmap="Blues", rasterized=True)
                for v, sty in truths[i]:
                    if np.isfinite(v):
                        a.axhline(v, **sty)
            for v, sty in truths[j]:
                if np.isfinite(v):
                    a.axvline(v, **sty)
            if logs[j]:
                a.set_xscale("log")
                _log_ticks(a.xaxis, *lims[j])
            a.set_xlim(lims[j])
            if i != j:
                if logs[i]:
                    a.set_yscale("log")
                    _log_ticks(a.yaxis, *lims[i])
                a.set_ylim(lims[i])
            a.tick_params(labelsize=7.5, which="both")
            if i == K - 1:
                a.set_xlabel(labels[j], fontsize=8.5)
            else:
                a.tick_params(labelbottom=False, which="both")
            if j == 0 and i > 0:
                a.set_ylabel(labels[i], fontsize=8.5)
            elif i != j:
                a.tick_params(labelleft=False, which="both")
            a.grid(False)
    return axs


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from .. import plotstyle as ps

    ps.apply()
    f = np.load(RESULTS / "step2_examples.npz", allow_pickle=True)
    tr = {k[6:]: float(f[k]) for k in f.files if k.startswith("truth_")}
    t_sty = dict(color=ps.C_TARGET, lw=1.3)
    c_sty = dict(color=ps.C_WRONG, lw=1.3, ls="--")
    titles = {"one": "one-planet model", "accel": "one planet + acceleration", "two": "two-planet model"}
    for name in ("no planet", "one planet", "two planets"):
        fig = plt.figure(figsize=(17, 6.6))
        subs = fig.subfigures(1, 3, wspace=0.03, width_ratios=[1, 1, 1.25])
        has_t, has_c = name != "no planet", name == "two planets"
        for sub, model in zip(subs, ("one", "accel", "two")):
            data = f[f"{name}|{model}"]
            planets = ([("t", t_sty)] if has_t else []) + ([("c", c_sty)] if has_c else [])

            def nearest(col):
                """Truth lines for the fitted planet in column `col` (period): the true planet closest in period."""
                if not planets:
                    return None
                med = np.nanmedian(data[:, col])
                return min(planets, key=lambda p_: abs(np.log(med / tr["P_" + p_[0]])))

            if model == "two":
                cols = [0, 1, 3, 4]
                labels = ["period, planet 1 (yr)", "mass, planet 1 (M_J)", "period, planet 2 (yr)",
                          "mass, planet 2 (M_J)"]
                logs = [True, True, True, True]
                truths = []
                for pc in (0, 3):
                    nb = nearest(pc)
                    truths += [[(tr["P_" + nb[0]], nb[1])], [(tr["m_" + nb[0]], nb[1])]] if nb else [[], []]
            else:
                cols = [0, 1, 2]
                labels = ["period (yr)", "mass (M_J)", "eccentricity"]
                logs = [True, True, False]
                nb = nearest(0)
                truths = ([[(tr["P_" + nb[0]], nb[1])], [(tr["m_" + nb[0]], nb[1])],
                           [(tr["e_" + nb[0]], nb[1])]] if nb else [[], [], []])
            _corner(sub, data[:, cols], labels, logs, truths, ps.C_TARGET)
            dchi = float(f[f"{name}|{model}|dchi2"])
            what = "gain over one planet" if model == "two" else "Δχ² vs no planet"
            sub.suptitle(f"{titles[model]}\n{what}: {dchi:.0f}", fontsize=11)
        gaia_txt = str(f[f"{name}|bic_choice"]) if f"{name}|bic_choice" in f.files else "?"
        handles = []
        if has_t:
            handles.append(Line2D([], [], label=f"true target ({tr['m_t']:.1f} M_J, {tr['P_t']:.1f} yr)", **t_sty))
        if has_c:
            handles.append(Line2D([], [], label=f"true companion ({tr['m_c']:.1f} M_J, {tr['P_c']:.0f} yr)", **c_sty))
        if handles:
            fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=10, bbox_to_anchor=(0.5, -0.02))
        fig.suptitle(f"Example star with {name.upper()}: what each model believes (posterior samples). "
                     f"Lowest BIC (DR4-like model choice): {gaia_txt}.\n"
                     "Lines mark the true planet closest in period to each fitted planet.", fontsize=12.5, y=1.06)
        out = RESULTS / f"step2_examples_{name.replace(' ', '_')}.png"
        fig.savefig(out, dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"saved {out}")
