"""Shared, plain-language plot styling for all figures.

Conventions: periods in YEARS on log axes with readable ticks; the Gaia DR4
time span drawn and labelled; the companion/target signal ratio labelled as
"1/10, equal, 10x" rather than log10 values; percentages, not fractions.
"""
from __future__ import annotations

import numpy as np

INK, MUTED, GRID = "#222222", "#8a8a8a", "#e6e6e6"
DR4_YEARS = (2458868.5 - 2456891.5) / 365.25   # 5.41 yr

# fixed categorical colours (same meaning in every figure)
C_TARGET = "#2a6fdb"      # target planet found correctly
C_COMPANION = "#8fb8f0"   # the other planet found correctly
C_WRONG = "#e07b39"       # a wrong planet published
C_ACCEL = "#c9a227"       # only an acceleration published
C_NONE = "#9e9e9e"        # nothing published


def apply():
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "axes.spines.top": False, "axes.spines.right": False, "font.size": 11,
        "axes.titlesize": 12, "axes.titleweight": "bold", "legend.frameon": False,
        "figure.dpi": 100,
    })


def years_axis(ax, which="x", lo=0.05, hi=60.0, label="companion's orbital period"):
    """Log axis in years with plain tick labels."""
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    ticks = [t for t in (0.1, 0.3, 1, 3, 10, 30) if lo <= t <= hi]
    fmt = FuncFormatter(lambda v, _: f"{v:g}")
    axis = ax.xaxis if which == "x" else ax.yaxis
    (ax.set_xscale if which == "x" else ax.set_yscale)("log")
    (ax.set_xlim if which == "x" else ax.set_ylim)(lo, hi)
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(fmt)
    axis.set_minor_formatter(NullFormatter())
    (ax.set_xlabel if which == "x" else ax.set_ylabel)(f"{label} (years)")


def mark_dr4(ax, text=True, y=0.5):
    """Vertical dashed line at the DR4 time span, with a vertical label on it."""
    ax.axvline(DR4_YEARS, color=INK, lw=1.2, ls="--", zorder=5)
    if text:
        ax.text(DR4_YEARS * 1.04, y, "Gaia DR4 time span (5.4 yr)", transform=ax.get_xaxis_transform(),
                rotation=90, va="center", ha="left", fontsize=8.5, color=INK, zorder=6,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))


def ratio_axis(ax, which="y", lo=0.01, hi=100.0,
               label="companion's wobble ÷ target's wobble"):
    """Log axis for the signal ratio with words instead of log10 numbers."""
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
    ticks = [t for t in (0.01, 0.1, 1, 10, 100) if lo <= t <= hi]
    names = {0.01: "1/100", 0.1: "1/10", 1: "equal", 10: "10×", 100: "100×"}
    axis = ax.yaxis if which == "y" else ax.xaxis
    (ax.set_yscale if which == "y" else ax.set_xscale)("log")
    (ax.set_ylim if which == "y" else ax.set_xlim)(lo, hi)
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FixedFormatter([names[t] for t in ticks]))
    axis.set_minor_formatter(NullFormatter())
    (ax.set_ylabel if which == "y" else ax.set_xlabel)(label)


def percent_axis(ax, which="y", label=None):
    from matplotlib.ticker import PercentFormatter
    axis = ax.yaxis if which == "y" else ax.xaxis
    axis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    if label:
        (ax.set_ylabel if which == "y" else ax.set_xlabel)(label)


def binned_map(x, y, value, mask, x_edges, y_edges, min_n=15):
    """Mean of `value` in each (x, y) cell (log-spaced edges); NaN if < min_n."""
    ix = np.digitize(x, x_edges) - 1
    iy = np.digitize(y, y_edges) - 1
    g = np.full((len(y_edges) - 1, len(x_edges) - 1), np.nan)
    for i in range(len(x_edges) - 1):
        for j in range(len(y_edges) - 1):
            m = mask & (ix == i) & (iy == j)
            if m.sum() >= min_n:
                g[j, i] = np.mean(value[m])
    return g
