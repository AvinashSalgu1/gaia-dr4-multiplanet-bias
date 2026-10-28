"""Step-2 summary text and figures."""
from __future__ import annotations

import numpy as np

from . import factors as Fct
from .analyses import ANALYSES, DR4_NAME, LABEL, SIG_MIN
from .outcomes import exclusive_categories, judge_all
from .sim import F2_CUT, RESULTS, T_DR4_DAYS as T

A_COL = {"one": "#9bbcf0", "accel": "#4f86dc", "jerk": "#7a5195", "two": "#1d3f8f", "best": "#2ca58d",
         "dr3": "#3a3a3a"}
SHORT = {"one": "one planet", "accel": "+ acceleration", "jerk": "+ jerk", "two": "two planets",
         "best": "best (BIC)", "dr3": "DR3-style"}
CAT_COL = {"both planets correct": "#0b2a6b", "target correct": "#2a6fdb",
           "companion correct (target missed)": "#8fb8f0", "a correct planet + a wrong one": "#f0b27a",
           "WRONG planet only": "#e07b39", "acceleration / jerk only": "#c9a227", "nothing reported": "#c4c4c4"}


def regimes(d):
    xc, rho = d["P_c"] / T, d["rho"]
    return [("orbit < 3.8 yr,\nweak wobble", (xc < 0.7) & (rho < 0.3)),
            ("orbit < 3.8 yr,\nstrong wobble", (xc < 0.7) & (rho >= 0.3)),
            ("orbit 3.8-11 yr", (xc >= 0.7) & (xc < 2)),
            ("orbit > 11 yr,\nsmaller wobble", (xc >= 2) & (rho < 1)),
            ("orbit > 11 yr,\nbigger wobble", (xc >= 2) & (rho >= 1))]


SHORT_REG = ["< 3.8 yr\nweak", "< 3.8 yr\nstrong", "3.8-11 yr", "> 11 yr\nsmaller", "> 11 yr\nbigger"]


def _pct(x):
    return f"{100 * x:5.1f}"


# ----------------------------------------------------------------------------- summary
def summary(d, outs, fac_tables=None):
    n = len(d["seed"])
    one = outs["one"]
    D = d["one_a_det"] == 1
    B = d["one_w_det"] == 1
    lines = ["=== STEP 2 SUMMARY ===",
             f"systems: {n}, median Gaia transits {np.median(d['n_transits']):.0f}",
             "",
             "[A] One-planet fit (Lammers & Winn), with vs without the companion:",
             f"  target detected alone {_pct(D.mean())}%  ->  with companion {_pct(B.mean())}%  "
             f"(lost {(D & ~B).sum()}, gained {(~D & B).sum()})",
             f"  of planets reported WITH the companion: target correct {_pct(one['target'][B].mean())}%, "
             f"companion correct {_pct(one['companion'][B].mean())}%, WRONG {_pct(one['wrong'][B].mean())}%",
             f"  WRONG planets that pass Gaia's DR3 quality cut (F2 < {F2_CUT:g}): "
             f"{_pct(np.mean(d['one_w_F2'][one['wrong'] & np.isfinite(d['one_w_F2'])] < F2_CUT))}%",
             "",
             "Analyses follow the orbit models Gaia DR4 publishes: " + "; ".join(
                 f"{LABEL[a]} = {DR4_NAME[a]}" for a in ANALYSES),
             "",
             "[B] All analyses. % of ALL systems: target correct / companion correct / "
             "WRONG planet reported / acceleration or jerk only / nothing   [false alarm on target-alone data]"]
    for name, mask in [("all systems", np.ones(n, bool))] + [(l.replace("\n", " "), m) for l, m in regimes(d)]:
        lines.append(f"  {name}  (n={mask.sum()})")
        for a in ANALYSES:
            o = outs[a]
            lines.append(f"     {LABEL[a]:<30}{_pct(o['target'][mask].mean())} /{_pct(o['companion'][mask].mean())} /"
                         f"{_pct(o['wrong'][mask].mean())} /{_pct(o['accel_only'][mask].mean())} /"
                         f"{_pct(o['none'][mask].mean())}   [{_pct(o['false_alarm'][mask].mean())}]")
    two = outs["two"]
    both = two["target"] & two["companion"]
    lines += ["",
              "[C] Two-planet analysis:",
              f"  both planets found correctly: {_pct(both.mean())}% overall; "
              + "; ".join(f"{l.replace(chr(10), ' ')} {_pct(both[m].mean()).strip()}%" for l, m in regimes(d)),
              f"  reports 0 / 1 / 2 planets: " + " / ".join(_pct(np.mean(two['n_reported'] == k)).strip() + "%"
                                                           for k in (0, 1, 2))]
    for other in ("one", "accel"):
        lines.append(f"  vs {LABEL[other]}: target rescued by two planets in {(two['target'] & ~outs[other]['target']).sum()}"
                     f" systems, lost in {(~two['target'] & outs[other]['target']).sum()}")
    from .outcomes import BIC_MODELS, bic_choice
    ch = bic_choice(d, "w")
    ch_a = bic_choice(d, "a")
    lines += ["", "[C2] Best model by BIC, which model it chooses (% of systems): with companion / target alone"]
    for k, nm in enumerate(BIC_MODELS):
        lines.append(f"     {nm:<28}{_pct(np.mean(ch == k))} / {_pct(np.mean(ch_a == k))}")
    lines.append("     by companion type, chosen model with companion: " + "; ".join(
        f"{l.replace(chr(10), ' ')}: " + ", ".join(f"{BIC_MODELS[k]} {100 * np.mean(ch[m] == k):.0f}%"
                                               for k in range(len(BIC_MODELS)) if np.mean(ch[m] == k) >= 0.05)
        for l, m in regimes(d)))
    lines += ["", f"[D] Would DR4's own label warn users? A WRONG planet with significance a0/sigma_a0 >= {SIG_MIN:g} "
              "is published as a normal orbit; below, as 'OrbitalPoorlyConstrained'.",
              "    % of ALL systems with a wrong planet: published as a normal orbit / as poorly constrained"]
    for a in ANALYSES:
        o = outs[a]
        lines.append(f"     {LABEL[a]:<30}{_pct(o['wrong_confident'].mean())} / {_pct(o['wrong_poorly_constrained'].mean())}"
                     f"   (normal-orbit share of the wrong ones: "
                     f"{_pct(o['wrong_confident'].sum() / max(o['wrong'].sum(), 1)).strip()}%)")
    if fac_tables:
        lines += ["", "[E] What else matters (after removing companion period and wobble ratio):",
                  "  importance = largest spread above chance, percentage points of 'target correct'"]
        for a, (table, order, score) in fac_tables.items():
            lines.append(f"  {LABEL[a]:<30}" + ", ".join(f"{k} {100 * score[k]:.0f}" for k in order[:6]))
    lines.append("(small runs: differences of a few points between analyses or regimes can be noise)")
    return lines


# ----------------------------------------------------------------------------- figures
def _setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .. import plotstyle as ps
    ps.apply()
    return plt, ps


def _map(ax, ps, d, value, mask, cmap, vmin=0, vmax=100, norm=None):
    P_yr = d["P_c"] / 365.25
    xe = np.geomspace(30 / 365.25, 10 * T / 365.25, 13)
    ye = np.geomspace(0.01, 100, 11)
    g = ps.binned_map(P_yr, d["rho"], value.astype(float), mask, xe, ye)
    kw = dict(norm=norm) if norm is not None else dict(vmin=vmin, vmax=vmax)
    im = ax.pcolormesh(xe, ye, 100 * g, cmap=cmap, **kw)
    ps.years_axis(ax, "x", lo=xe[0], hi=xe[-1])
    ps.ratio_axis(ax, "y")
    ps.mark_dr4(ax, text=False)
    ax.axhline(1, color=ps.MUTED, lw=1)
    return im


def _grouped_bars(ax, ps, d, outs, key="target", ylabel="systems where the TARGET is found correctly"):
    reg = regimes(d)
    w = 0.8 / len(ANALYSES)
    for k, a in enumerate(ANALYSES):
        ax.bar(np.arange(len(reg)) + (k - (len(ANALYSES) - 1) / 2) * w, [outs[a][key][m].mean() for _, m in reg], w,
               color=A_COL[a], label=LABEL[a])
    alone = [outs["one"]["alone_detected"][m].mean() for _, m in reg]
    for i, v in enumerate(alone):
        ax.plot([i - 0.42, i + 0.42], [v, v], color=ps.INK, lw=1.2, ls=":",
                label="one planet, no companion" if i == 0 else None)
    ax.set_xticks(range(len(reg)), SHORT_REG, fontsize=8.5)
    ax.set_xlabel("companion's orbit (and how big its wobble is compared with the target's)")
    ax.set_ylim(0, 1)
    ps.percent_axis(ax, "y", ylabel)


def _stacked(ax, ps, outs, masks, names):
    cats = list(CAT_COL)
    for i, (a, m) in enumerate(masks):
        ex = exclusive_categories(outs[a])
        left = 0.0
        for c in cats:
            v = ex[c][m].mean()
            ax.barh(i, v, left=left, color=CAT_COL[c], height=0.7, edgecolor="white", lw=0.5)
            left += v
    ax.set_yticks(range(len(masks)), names, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ps.percent_axis(ax, "x", "share of systems")


def _cat_legend(fig_or_ax, **kw):
    from matplotlib.patches import Patch
    fig_or_ax.legend(handles=[Patch(color=c, label=l) for l, c in CAT_COL.items()], **kw)


def overview(d, outs, fac):
    plt, ps = _setup()
    from matplotlib.colors import TwoSlopeNorm
    fig = plt.figure(figsize=(19, 11.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.05])
    one = outs["one"]
    B = d["one_w_det"] == 1
    D = d["one_a_det"] == 1

    a = fig.add_subplot(gs[0, 0])
    im = _map(a, ps, d, one["wrong"], B, "Oranges")
    fig.colorbar(im, ax=a, label="% of reported planets that are WRONG")
    a.set_title("(a) One-planet fit: how often is the\nreported planet wrong?")

    a = fig.add_subplot(gs[0, 1])
    im = _map(a, ps, d, B.astype(float) - D.astype(float), np.ones(len(B), bool), "RdBu",
              norm=TwoSlopeNorm(0, -50, 50))
    fig.colorbar(im, ax=a, label="change in detection chance (percentage points)")
    a.set_title("(b) One-planet fit: does the companion make\ndetection harder (red) or easier (blue)?")

    a = fig.add_subplot(gs[0, 2])
    _grouped_bars(a, ps, d, outs)
    a.set_title("(c) Target found correctly, by analysis\n(dotted: one planet, companion removed)")
    a.legend(fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)

    a = fig.add_subplot(gs[1, 0])
    masks = [(an, np.ones(len(B), bool)) for an in ANALYSES]
    _stacked(a, ps, outs, masks, [LABEL[an] for an in ANALYSES])
    a.set_title("(d) What each analysis reports (all systems)")
    _cat_legend(a, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=8.5)

    a = fig.add_subplot(gs[1, 1])
    x = np.arange(len(ANALYSES))
    conf = np.array([outs[an]["wrong_confident"].mean() for an in ANALYSES])
    weak = np.array([outs[an]["wrong_poorly_constrained"].mean() for an in ANALYSES])
    a.bar(x, conf, color=ps.C_WRONG, label=f"published as a normal orbit (significance ≥ {SIG_MIN:g})")
    a.bar(x, weak, bottom=conf, color="#f3c9a6", label="flagged 'poorly constrained' (significance < 5)")
    for i in x:
        a.text(i, conf[i] + weak[i] + 0.004, f"{100 * (conf[i] + weak[i]):.1f}%", ha="center", fontsize=9)
    fa = np.array([outs[an]["false_alarm"].mean() for an in ANALYSES])
    two_line = {"one": "one\nplanet", "accel": "+ accel-\neration", "jerk": "+ jerk\n", "two": "two\nplanets",
                "best": "best\n(BIC)", "dr3": "DR3-\nstyle"}
    a.set_xticks(x, [f"{two_line[an]}\n[{100 * f:.1f}%]" for an, f in zip(ANALYSES, fa)], fontsize=8.5)
    ps.percent_axis(a, "y", "systems where a WRONG planet is reported")
    a.set_ylim(0, 1.3 * max(0.05, (conf + weak).max()))
    a.legend(fontsize=8.5, loc="upper left")
    a.set_title("(e) Wrong planets: would DR4's own label warn you?\n"
                "[ ] = false alarms: a planet reported with the companion removed")

    a = fig.add_subplot(gs[1, 2])
    table, order, score = fac["one"]
    top = order[:8]
    labels = {k: lab for k, lab, _, _ in Fct.factor_table(d)}
    a.barh(range(len(top)), [100 * score[k] for k in top], color=ps.C_TARGET)
    a.set_yticks(range(len(top)), [labels[k] for k in top], fontsize=8)
    a.invert_yaxis()
    a.axvline(0, color=ps.MUTED, lw=1)
    a.set_xlabel("effect on 'target correct' beyond chance\n(percentage points, one-planet fit)")
    a.set_title("(f) What else matters, besides the\ncompanion's period and wobble size")

    fig.suptitle(f"Step 2: a hidden second planet and the orbit models of Gaia DR4 "
                 f"({len(B)} simulated systems; target always in Gaia's detectable range)", fontsize=14)
    out = RESULTS / "step2_overview.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def detail_maps(d, outs):
    plt, ps = _setup()
    fig, axs = plt.subplots(2, len(ANALYSES), figsize=(5.4 * len(ANALYSES), 9.5), constrained_layout=True)
    for k, a in enumerate(ANALYSES):
        o = outs[a]
        im0 = _map(axs[0, k], ps, d, o["target"], np.ones(len(o["target"]), bool), "Blues")
        axs[0, k].set_title(f"{LABEL[a]}\ntarget found correctly")
        rep = o["n_reported"] > 0
        im1 = _map(axs[1, k], ps, d, o["wrong"], rep, "Oranges")
        axs[1, k].set_title(f"{LABEL[a]}\nreported planets that are WRONG")
    fig.colorbar(im0, ax=axs[0, :], label="% of systems", shrink=0.85)
    fig.colorbar(im1, ax=axs[1, :], label="% of systems with a reported planet", shrink=0.85)
    fig.suptitle("Step 2 detail: where each analysis succeeds and fails "
                 "(dashed line = Gaia DR4 time span)", fontsize=13)
    out = RESULTS / "step2_maps.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def detail_outcomes(d, outs):
    plt, ps = _setup()
    reg = regimes(d)
    fig, axs = plt.subplots(1, 3, figsize=(20, 9.0), constrained_layout=True,
                            gridspec_kw=dict(width_ratios=[1.25, 1.25, 0.9]))
    for ax, part in ((axs[0], reg[:3]), (axs[1], reg[3:])):
        masks = [(a, m) for _, m in part for a in ANALYSES]
        names = [f"{rn.replace(chr(10), ' ')}  ·  {SHORT[a]}" for rn, _ in part for a in ANALYSES]
        _stacked(ax, ps, outs, masks, names)
        for g in range(1, len(part)):
            ax.axhline(len(ANALYSES) * g - 0.5, color=ps.INK, lw=0.8)
    axs[0].set_title("What each analysis reports, by companion type")
    axs[1].set_title("(continued)")
    _cat_legend(fig, loc="lower center", ncol=4, fontsize=9.5, bbox_to_anchor=(0.36, -0.07))

    a = axs[2]
    wrong = outs["one"]["wrong"]
    for i, (name, m) in enumerate(reg):
        w = wrong & m & np.isfinite(d["one_w_F2"])
        if w.sum() == 0:
            continue
        flagged = np.mean(d["one_w_F2"][w] > F2_CUT)
        a.barh(i, 1 - flagged, color=ps.C_WRONG, height=0.6, label="looks fine: would be published" if i == 0 else None)
        a.barh(i, flagged, left=1 - flagged, color=ps.INK, height=0.6,
               label="flagged by Gaia's quality cut" if i == 0 else None)
        a.text(0.02, i, f"{100 * (1 - flagged):.0f}% pass  (n={w.sum()})", va="center", color="white",
               fontsize=9, fontweight="bold")
    a.set_yticks(range(len(reg)), [r[0] for r in reg], fontsize=9)
    a.invert_yaxis()
    a.set_xlim(0, 1)
    ps.percent_axis(a, "x", "share of WRONG planets (one-planet fit)")
    a.set_title(f"Would Gaia's DR3 quality cut (F2 < {F2_CUT:g})\ncatch the wrong ones?")
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), fontsize=9)
    fig.suptitle("Step 2 detail: outcomes for every companion type and analysis", fontsize=13)
    out = RESULTS / "step2_outcomes.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def analyse_population(d, make_plots=True, n_perm=20):
    outs = judge_all(d)
    fac = {}
    for a in ANALYSES:
        r, table = Fct.importance(d, outs[a]["target"], n_perm=n_perm if a == "one" else max(5, n_perm // 4))
        order, score = Fct.ranking(table, [k for k, *_ in Fct.factor_table(d)])
        fac[a] = (table, order, score)
        if a == "one":
            r_one = r
    lines = summary(d, outs, fac)
    print("\n".join(lines))
    (RESULTS / "step2_summary.txt").write_text("\n".join(lines) + "\n")
    if make_plots:
        overview(d, outs, fac)
        detail_maps(d, outs)
        detail_outcomes(d, outs)
        table, order, _ = fac["one"]
        Fct.plot_grid(d, r_one, table, order, RESULTS / "step2_factors.png")
    return outs, fac
