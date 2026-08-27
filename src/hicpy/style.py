"""Single publication style block for the whole project.

Import this ONCE at the top of every plotting script:

    from hicpy.style import apply_style, OKABE_ITO, stamp
    apply_style()

Rationale (analysis-notebook-standard): every figure in a paper must share one
style block so that fonts, sizes and colours are consistent, and every curve must
stay distinguishable in grayscale and for colour-blind readers.
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# Okabe-Ito colour-blind-safe cycle.
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
             "#E69F00", "#56B4E9", "#F0E442", "#000000"]

# Marker / linestyle cycles: vary shape too, so grayscale printing still works.
MARKERS = ["o", "s", "^", "v", "D", "P", "X", "*"]
LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), (0, (5, 1))]


def apply_style() -> None:
    """Apply the project-wide matplotlib rcParams (HEP conventions)."""
    mpl.rcParams.update({
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 16,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 12,
        "legend.frameon": False,
        "lines.linewidth": 1.8,
        "lines.markersize": 6,
        "figure.figsize": (5.0, 3.8),
        "figure.dpi": 150,
        "savefig.bbox": "tight",
        "errorbar.capsize": 2,
        # HEP convention: inward ticks on all four sides, minor ticks on.
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "axes.prop_cycle": mpl.cycler(color=OKABE_ITO),
        "axes.formatter.use_mathtext": True,
    })


def stamp(ax, lines, loc="upper right", **kw):
    """Put the dataset stamp (system, energy, cuts) on the plot face.

    Physics figures must carry their own provenance: collision system,
    sqrt(s) or sqrt(s_NN), event class and kinematic cuts.
    """
    xy = {"upper right": (0.97, 0.95), "upper left": (0.03, 0.95),
          "lower right": (0.97, 0.05), "lower left": (0.03, 0.05)}[loc]
    ha = "right" if "right" in loc else "left"
    va = "top" if "upper" in loc else "bottom"
    ax.text(xy[0], xy[1], "\n".join(lines), transform=ax.transAxes,
            ha=ha, va=va, fontsize=11, linespacing=1.4, **kw)
    return ax


__all__ = ["apply_style", "stamp", "OKABE_ITO", "MARKERS", "LINESTYLES", "plt"]
