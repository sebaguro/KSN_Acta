#!/usr/bin/env python3
"""
make_fig3_deltaP.py
-------------------
Figure 3 of the paper: the year-over-year increments Delta P_i = P_i - P_(i-1)
of global energy production, 1966-2024, with the Shapiro-Wilk statistics and
the largest falls labelled.

The drawing follows the Delta P panel of the original analysis figure
(ksn_figure_optA.py in the KSN_paper_I repository), from which the submitted
Figure 3 was taken; this script draws it on its own from the data in data/.

Usage:  python3 make_fig3_deltaP.py
Writes outputs/figures/fig3_deltaP.{png,pdf}.
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from ksn_common import load_energy, OUT

BLUE = "#2166ac"

years, P = load_energy()
dP = np.diff(P)
dyears = years[1:]
W, p = stats.shapiro(dP)
skew = stats.skew(dP)

fig, ax = plt.subplots(figsize=(6.6, 5.3))
ax.scatter(dyears, dP / 1e9, color=BLUE, s=28, zorder=5)
ax.axhline(0, color="k", lw=1)

events = [                       # year, label, label position (year, GW)
    (1980, "1980: Oil crisis", (1971.0, -165)),
    (1981, "1981: Recession", (1981.0, -300)),
    (1982, "1982: Recession", (1991.0, -170)),
    (2009, "2009: GFC", (2013.0, -110)),
    (2020, "2020: COVID-19", (2008.0, -560)),
]
for yr, label, (x_t, y_t) in events:
    y_pt = dP[list(dyears).index(yr)] / 1e9
    ax.annotate(label, xy=(yr, y_pt), xytext=(x_t, y_t), ha="center", va="center",
                fontsize=8.5, color="red", fontweight="bold",
                arrowprops=dict(arrowstyle="->", color="red", lw=1.4,
                                shrinkA=2, shrinkB=3),
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="red",
                          alpha=0.9))

ax.set_xlabel("Year", fontsize=11)
ax.set_ylabel(r"$\Delta P_i$  ($\times10^9$ W)", fontsize=11)
ax.grid(True, alpha=0.2)
ax.text(0.05, 0.03,
        f"$\\Delta P_i$  Shapiro\u2013Wilk:\n"
        f"$W={W:.3f}$,  $p={p:.4f}$\n"
        f"Skewness $= {skew:.3f}$\n"
        "Rejects normality ($\\alpha=0.05$)",
        transform=ax.transAxes, va="bottom", fontsize=9,
        bbox=dict(boxstyle="round", fc="white", alpha=0.85))
fig.tight_layout()

figdir = os.path.join(OUT, "figures")
os.makedirs(figdir, exist_ok=True)
fig.savefig(os.path.join(figdir, "fig3_deltaP.png"), dpi=300)
fig.savefig(os.path.join(figdir, "fig3_deltaP.pdf"),
            metadata={"CreationDate": None})   # no timestamp: same file on every run
print(f"Shapiro-Wilk on Delta P: W = {W:.3f}, p = {p:.4f}, skewness = {skew:.3f}")
print("Largest falls: " + ", ".join(f"{y}: {d / 1e9:+.0f} GW"
                                     for y, d in zip(dyears, dP) if d < 0))
