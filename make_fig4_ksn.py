#!/usr/bin/env python3
"""
make_fig4_ksn.py
----------------
Figure 4 of the paper: B(t) = P(t)/H(t) for 2009-2024, and the annual series
behind it (outputs/ksn_annual_series.csv).

H(t): average hashrate over each calendar year. A block mined at difficulty D
requires on average D*2^32 hash evaluations, so the year's average hashrate is
2^32 times the summed difficulty of the blocks mined in that year, divided by
the length of the year. Coin Metrics publishes this quantity day by day
(column HashRate, TH/s), so H(t) is the mean of the daily values over the year
(hashrate_series.py rebuilds it independently from the difficulty record).

P(t): world primary energy from Our World in Data, converted to watts with
Eq. (1) of the paper, P = E * 1e12 / (365.25 * 24).

Usage:  python3 make_fig4_ksn.py
Writes outputs/figures/fig4_ksn.{png,pdf} and outputs/ksn_annual_series.csv.
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ksn_common import load_energy, annual_hashrate, OUT

YEARS = list(range(2009, 2025))
yrs_h, H_arr, NB_arr = annual_hashrate(YEARS)
H = dict(zip(YEARS, H_arr))                                     # H s^-1
NB = dict(zip(YEARS, NB_arr))                                   # blocks in the year
yrs_e, P_arr = load_energy()
P = {int(y): p for y, p in zip(yrs_e, P_arr) if int(y) in YEARS}  # W
B = {y: P[y] / H[y] for y in YEARS}                             # J per hash

FIGDIR = os.path.join(OUT, "figures")
os.makedirs(FIGDIR, exist_ok=True)
with open(os.path.join(OUT, "ksn_annual_series.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["year", "blocks", "H_hash_per_s", "P_W", "B_J_per_hash"])
    for y in YEARS:
        w.writerow([y, NB[y], "%.4e" % H[y], "%.4e" % P[y], "%.4e" % B[y]])

# --- figure (layout and colours as in the submitted Figure 4) ------------------
plt.rcParams.update({"font.size": 10, "axes.labelsize": 11,
                     "xtick.labelsize": 10, "ytick.labelsize": 10,
                     "legend.fontsize": 9.5})
fig, ax = plt.subplots(figsize=(6.3, 5.2))
eras = [(2009.0, 2010.5, "forestgreen", "CPU era begins (Jan 2009)"),
        (2010.5, 2013.0, "royalblue", r"GPU era begins ($\sim$mid 2010)"),
        (2013.0, 2024.5, "crimson", "ASIC era begins (Jan 2013)")]
for x0, x1, c, lab in eras:
    ax.axvspan(x0, x1, color=c, alpha=0.07, lw=0)
for x0, x1, c, lab in eras:
    ax.axvline(x0, color=c, ls="--", lw=2.0, alpha=0.85, label=lab)
ax.plot(YEARS, [B[y] for y in YEARS], ls="none", marker="+", ms=10, mew=2.5,
        color="#1a9641", zorder=5, label=r"$B(t)=\frac{P(t)}{H(t)}$  (KarNak)")
ax.set_yscale("log")
ax.set_yticks([10.0**k for k in range(-7, 8, 2)])
ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.set_xlim(2008.225, 2025.275)
lo, hi = np.log10(min(B.values())), np.log10(max(B.values()))
m = 0.05 * (hi - lo)
ax.set_ylim(10**(lo - m), 10**(hi + m))
ax.grid(True, which="major", color="0.85", lw=0.8, alpha=0.6)
ax.set_xlabel("Year")
ax.set_ylabel(r"$B(t)=\frac{P(t)}{H(t)}$  (KarNak)")
h, l = ax.get_legend_handles_labels()
order = [3, 0, 1, 2]
ax.legend([h[i] for i in order], [l[i] for i in order], loc="upper right",
          framealpha=0.9)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "fig4_ksn.png"), dpi=300)
fig.savefig(os.path.join(FIGDIR, "fig4_ksn.pdf"),
            metadata={"CreationDate": None})   # no timestamp: same file on every run
print("year  blocks   H (H/s)     P (W)        B (J/hash)")
for y in YEARS:
    print(y, "%6d" % NB[y], "%.3e" % H[y], " %.4e" % P[y], " %.3e" % B[y])
print("span of B: %.1f dex; span of H: %.1f dex   [paper: 14 orders of magnitude]"
      % (np.log10(B[2009] / B[2024]), np.log10(H[2024] / H[2009])))
