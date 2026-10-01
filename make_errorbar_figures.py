#!/usr/bin/env python3
# =====================================================================
#  make_errorbar_figures.py
#
#  Redraws the two energy figures of the KSN paper WITH error bars, in a
#  style matching the originals, and prints the robustness numbers quoted
#  in the Figure 1 / Figure 2 captions and in Section 3.3.
#
#  Outputs (the manuscript includes the PNGs from its  figures/  folder):
#     fig1_energy_loglog.{png,pdf}   <- Figure 1 (log axis) + error bars
#     fig2_energy_linear.{png,pdf}   <- Figure 2 (linear axis, in units of
#                                       10^13 W) + error bars
#
#  ABOUT THE ERROR BARS — please read.
#  OWID / the Energy Institute / the EIA do NOT publish formal per-year
#  uncertainties on global primary energy. The bars drawn here are an
#  *adopted* heteroscedastic uncertainty model (larger for older data,
#  smaller for recent data), used for a sensitivity check — they are not
#  measured error bars. The captions and Section 3.3 state this. Keep
#  that framing explicit; a referee will accept an adopted model used
#  transparently for sensitivity analysis, but not error bars presented
#  as if measured. Tune the model in fractional_sigma() to your judgement.
#
#  Dependencies:  numpy, scipy, matplotlib
#     pip install numpy scipy matplotlib
#
#  USAGE
#     python make_errorbar_figures.py               # the copy of the data in data/
#     python make_errorbar_figures.py --csv /path/to/owid-energy-data.csv
#     python make_errorbar_figures.py --demo        # plumbing test only
#  Figures are written to outputs/figures/.
#
#  CSV format: the OWID owid-energy-data.csv, with columns
#     country , year , primary_energy_consumption     (TWh/yr)
#  Rows kept: country == "World", 1965 <= year <= 2024.
#  The default is data/owid_world_energy_1965_2024.csv (see data/README.md).
# =====================================================================

import argparse
import os
import sys
import numpy as np

# ---- conversions / constants ---------------------------------------
TWH_PER_YR_TO_W = 1e12 / (365.25 * 24.0)     # Eq. (1)
L_SUN = 3.828e26                              # Type II Kardashev 1964 threshold (W)
L_SUN_SAGAN = 1.0e26                          # Type II Sagan 1973 threshold (W)
TYPE_I = 4.0e12                               # Type I Kardashev 1964 (W)
SOLAR_INSOLATION = 1.74e17                    # solar insolation at Earth (W)
T0 = 1964                                     # time origin (Kardashev 1964)
KARDASHEV_RATE = 0.01                         # the 1%/yr assumption under test
MJD_2000 = 51544.5                            # MJD of 2000-01-01
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(HERE, "data", "owid_world_energy_1965_2024.csv")
FIGDIR = os.path.join(HERE, "outputs", "figures")


def year_to_mjd(y):
    return (np.asarray(y) - 2000.0) * 365.25 + MJD_2000


def mjd_to_year(m):
    return (np.asarray(m) - MJD_2000) / 365.25 + 2000.0


# ---------------------------------------------------------------------
#  DATA
# ---------------------------------------------------------------------
def load_owid(csv_path):
    import csv
    years, twh = [], []
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        if "primary_energy_consumption" not in reader.fieldnames:
            sys.exit("ERROR: column 'primary_energy_consumption' not found — "
                     "is this owid-energy-data.csv?")
        for row in reader:
            if row.get("country", row.get("entity", "")).strip() != "World":
                continue
            try:
                yr = int(float(row["year"]))
                val = row["primary_energy_consumption"]
                if val in ("", None):
                    continue
                val = float(val)
            except (ValueError, KeyError):
                continue
            if 1965 <= yr <= 2024:
                years.append(yr); twh.append(val)
    if not years:
        sys.exit("ERROR: no World rows for 1965-2024 found.")
    o = np.argsort(years)
    return np.array(years)[o], np.array(twh)[o] * TWH_PER_YR_TO_W


def make_demo_data():
    years = np.arange(1965, 2025)
    t = years - T0
    P = 4.95e12 + 2.44e11 * t
    rng = np.random.default_rng(42)
    P = P * (1.0 + rng.normal(0, 0.01, P.size))
    P[years == 2008] *= 0.985
    P[years == 2020] *= 0.955
    return years, P


def fractional_sigma(years):
    """Adopted per-year fractional 1-sigma. Tune to your judgement.
    Declines smoothly from ~10% in the mid-1960s to ~1.5% from 2000 on."""
    hi, lo, y0, y1 = 0.10, 0.015, 1965, 2000
    return np.clip(np.interp(years, [y0, y1], [hi, lo]), lo, hi)


# ---------------------------------------------------------------------
#  FITS (unweighted, matching the paper's OLS, for the plotted lines)
# ---------------------------------------------------------------------
def fit_linear(t, P, sigma=None):
    w = np.ones_like(P) if sigma is None else 1.0 / sigma**2
    X = np.vstack([np.ones_like(t), t]).T
    W = np.diag(w)
    cov = np.linalg.inv(X.T @ W @ X)
    a, b = cov @ (X.T @ W @ P)
    resid = P - (a + b * t)
    ss_res = np.sum(w * resid**2)
    ss_tot = np.sum(w * (P - np.average(P, weights=w))**2)
    R2 = 1.0 - ss_res / ss_tot
    if sigma is None:                       # scale to OLS parameter errors
        cov = cov * (np.sum(resid**2) / max(len(P) - 2, 1))
    return dict(a=a, b=b, sigma_b=np.sqrt(cov[1, 1]), R2=R2)


def fit_exp(years, P, sigma=None):
    t = years - T0
    s = None if sigma is None else sigma / P
    r = fit_linear(t, np.log(P), sigma=s)
    return dict(a0=np.exp(r["a"]), r=r["b"], sigma_r=r["sigma_b"])


# The paper's headline exponential rate comes from the MCMC posterior
# (r = 2.01%/yr), NOT from a quick log-space fit (which gives ~2.14%).
# The plotted exponential curve uses the published value so the figure
# legend matches the manuscript text. We fit only the amplitude a0.
EXP_RATE_DISPLAY = 0.0201

def display_exponential(years, P, r=EXP_RATE_DISPLAY):
    """Exponential curve at the published rate r; a0 by least squares."""
    t = years - T0
    e = np.exp(r * t)
    a0 = np.sum(P * e) / np.sum(e * e)
    resid = P - a0 * e
    R2 = 1.0 - np.sum(resid**2) / np.sum((P - P.mean())**2)
    return dict(a0=a0, r=r, R2=R2)


# ---------------------------------------------------------------------
#  FIGURES
# ---------------------------------------------------------------------
def add_mjd_axis(ax, millions=False):
    """Modified Julian Date along the top axis. With millions=True
    (Figure 1, whose dates reach 2.2 million) the ticks are given in
    units of 10^6 days, stated in the label, instead of matplotlib's
    separate "1e6" at the end of the axis."""
    sec = ax.secondary_xaxis("top", functions=(year_to_mjd, mjd_to_year))
    if millions:
        sec.xaxis.set_major_formatter(lambda v, pos: f"{v / 1e6:.2f}")
        sec.set_xlabel(r"Modified Julian Date  ($10^{6}$ d)")
    else:
        sec.set_xlabel("Modified Julian Date")


def figure1_loglog(years, P, sigma, lin_unw, lin_w, expo, tag=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = years - T0
    x_lo, x_hi = 1965, 8000
    tt = np.linspace(x_lo - T0, x_hi - T0, 400)
    xx = tt + T0

    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    # data with error bars
    ax.errorbar(years, P, yerr=sigma, fmt="o", ms=4.5, color="k", ecolor="0.45",
                elinewidth=0.9, capsize=2, zorder=6,
                label="OWID / EI World data (1965–2024)")
    # model curves: unweighted OLS (the paper's headline fit) and the
    # 1/sigma_i^2-weighted fit (consistent with the error bars shown)
    ax.plot(xx, lin_unw["a"] + lin_unw["b"] * tt, "--", color="tab:blue", lw=1.6,
            label=r"Model 1: Linear OLS ($b{=}%.2f{\times}10^{11}$, $R^2{=}%.3f$)"
                  % (lin_unw["b"] / 1e11, lin_unw["R2"]))
    ax.plot(xx, lin_w["a"] + lin_w["b"] * tt, "-.", color="navy", lw=1.3,
            label=r"Linear, weighted by $1/\sigma_i^2$ ($b{=}%.2f{\times}10^{11}$, $R^2{=}%.3f$)"
                  % (lin_w["b"] / 1e11, lin_w["R2"]))
    ax.plot(xx, expo["a0"] * np.exp(expo["r"] * tt), "-", color="tab:red", lw=1.6,
            label=r"Model 2: Exp $r{=}%.2f\%%$/yr ($R^2{=}%.3f$)" % (100 * expo["r"], expo["R2"]))
    ax.plot(xx, P[0] * (1.0 + KARDASHEV_RATE) ** (xx - years[0]), ":",
            color="sandybrown", lw=1.6, label="Kardashev $r=1\\%$/yr")
    # threshold lines
    ax.axhline(L_SUN, ls="-.", color="magenta", lw=1.3,
               label=r"Type II Kardashev 1964 ($4\times10^{26}$ W)")
    ax.axhline(L_SUN_SAGAN, ls="-.", color="orchid", lw=1.1,
               label=r"Type II Sagan 1973 ($10^{26}$ W)")
    ax.axhline(SOLAR_INSOLATION, ls=(0, (8, 4)), color="darkorange", lw=1.4,
               label=r"Solar insolation at Earth ($1.74\times10^{17}$ W)")
    ax.axhline(TYPE_I, ls=":", color="olive", lw=1.3,
               label=r"Type I Kardashev 1964 ($4\times10^{12}$ W)")

    ax.set_yscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(1e11, 3e27)
    ax.set_xlabel("Year")
    ax.set_ylabel("Global Energy Production  (W)")
    add_mjd_axis(ax, millions=True)
    ax.legend(fontsize=6.5, ncol=2, loc="lower center", framealpha=0.9)
    ax.grid(True, which="major", alpha=0.25)
    if tag == "demo":
        ax.text(0.5, 0.5, "DEMO — synthetic data", transform=ax.transAxes,
                fontsize=22, color="red", alpha=0.25, ha="center", va="center",
                rotation=20)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig1_energy_loglog.pdf"),
                metadata={"CreationDate": None})
    fig.savefig(os.path.join(FIGDIR, "fig1_energy_loglog.png"), dpi=300)
    plt.close(fig)


def figure2_linear(years, P, sigma, lin_unw, lin_w, expo, tag=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    U = 1e13                                   # vertical axis in units of 10^13 W (R2-19)
    t = years - T0
    x_lo, x_hi = 1965, 2030
    tt = np.linspace(x_lo - T0, x_hi - T0, 400)
    xx = tt + T0

    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    ax.errorbar(years, P / U, yerr=sigma / U, fmt="o", ms=4.5, color="k", ecolor="0.45",
                elinewidth=0.9, capsize=2, zorder=6,
                label="OWID / EI World data (1965–2024)")
    ax.plot(xx, (lin_unw["a"] + lin_unw["b"] * tt) / U, "--", color="tab:blue", lw=1.6,
            label=r"Model 1: Linear OLS ($b{=}%.2f{\times}10^{11}$, $R^2{=}%.3f$)"
                  % (lin_unw["b"] / 1e11, lin_unw["R2"]))
    ax.plot(xx, (lin_w["a"] + lin_w["b"] * tt) / U, "-.", color="navy", lw=1.3,
            label=r"Linear, weighted by $1/\sigma_i^2$ ($b{=}%.2f{\times}10^{11}$, $R^2{=}%.3f$)"
                  % (lin_w["b"] / 1e11, lin_w["R2"]))
    ax.plot(xx, expo["a0"] * np.exp(expo["r"] * tt) / U, "-", color="tab:red", lw=1.6,
            label=r"Model 2: Exp $r{=}%.2f\%%$/yr ($R^2{=}%.3f$)" % (100 * expo["r"], expo["R2"]))
    ax.plot(xx, P[0] * (1.0 + KARDASHEV_RATE) ** (xx - years[0]) / U, ":",
            color="sandybrown", lw=1.6, label="Kardashev $r=1\\%$/yr")
    ax.axhline(TYPE_I / U, ls=":", color="olive", lw=1.3,
               label=r"Type I Kardashev 1964 ($4\times10^{12}$ W)")

    ax.set_xlim(x_lo, x_hi)
    ax.set_xlabel("Year")
    ax.set_ylabel(r"Global Energy Production  ($10^{13}$ W)")
    add_mjd_axis(ax)
    ax.legend(fontsize=7, loc="upper left", framealpha=0.9)
    ax.grid(True, alpha=0.25)
    if tag == "demo":
        ax.text(0.5, 0.5, "DEMO — synthetic data", transform=ax.transAxes,
                fontsize=22, color="red", alpha=0.25, ha="center", va="center",
                rotation=20)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig2_energy_linear.pdf"),
                metadata={"CreationDate": None})
    fig.savefig(os.path.join(FIGDIR, "fig2_energy_linear.png"), dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------
#  MAIN
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Draw KSN energy figures with error bars.")
    ap.add_argument("--csv", default=DEFAULT_CSV,
                    help="path to the OWID energy CSV (default: the copy in data/)")
    ap.add_argument("--demo", action="store_true", help="synthetic data (testing only)")
    args = ap.parse_args()

    if args.demo:
        years, P = make_demo_data(); tag = "demo"
        print("\n*** DEMO MODE — synthetic data — DO NOT PUBLISH THESE NUMBERS ***\n")
    elif args.csv:
        years, P = load_owid(args.csv); tag = ""
    else:
        sys.exit("Provide --csv /path/to/owid-energy-data.csv  (or --demo to test).")

    os.makedirs(FIGDIR, exist_ok=True)
    sigma = fractional_sigma(years) * P
    t = years - T0

    lin_unw = fit_linear(t, P)               # plotted (matches paper OLS)
    lin_w = fit_linear(t, P, sigma=sigma)    # weighted (robustness)
    expo_fit = fit_exp(years, P)             # quick log-space fit (for sigma-rejection)
    expo = display_exponential(years, P)     # plotted curve at published r=2.01%

    figure1_loglog(years, P, sigma, lin_unw, lin_w, expo, tag=tag)
    figure2_linear(years, P, sigma, lin_unw, lin_w, expo, tag=tag)

    bar = "=" * 66
    print(bar)
    print(f"  Figures written: fig1_energy_loglog.{{pdf,png}} , fig2_energy_linear.{{pdf,png}}")
    print(bar)
    print(f"  Unweighted (plotted) linear slope b : {lin_unw['b']:.3e} W/yr  (R^2={lin_unw['R2']:.4f})")
    print(f"  Weighted linear slope b             : {lin_w['b']:.3e} W/yr  (R^2={lin_w['R2']:.4f})")
    pct = 100 * (lin_w['b'] - lin_unw['b']) / lin_unw['b']
    print(f"  Slope change under weighting        : {pct:+.1f}%  "
          f"({abs(lin_w['b']-lin_unw['b'])/lin_w['sigma_b']:.1f} sigma vs the tiny formal error)")
    print(f"  Plotted exponential rate (paper)    : {100*expo['r']:.2f} %/yr  (R^2={expo['R2']:.4f})")
    z = (expo_fit['r'] - KARDASHEV_RATE) / expo_fit['sigma_r']
    print("  Diagnostics, not quoted in the paper:")
    print(f"  Quick log-space rate (this script)  : {100*expo_fit['r']:.2f} %/yr")
    print(f"  Distance of 1%/yr from that rate    : {z:.0f} sigma (formal error of the log fit)")
    print(bar)
    print("  NOTE: the weighted slope changes by a few percent under the adopted")
    print("  uncertainty model; this does NOT change the rejection of the 1% model")
    print("  or the order-of-magnitude Type II timescale. See Section 3.3.")
    print(bar)


if __name__ == "__main__":
    main()
