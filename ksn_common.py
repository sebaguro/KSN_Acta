"""
ksn_common.py
-------------
Data loaders and constants shared by the scripts of this repository.

Every script reads the copies of the data kept in data/ (see data/README.md
for their sources, versions and licences), so the results do not change when
the upstream data sets are revised.
"""

import csv
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "outputs")
OWID_CSV = os.path.join(DATA, "owid_world_energy_1965_2024.csv")
CM_CSV = os.path.join(DATA, "coinmetrics_btc_2009_2024.csv")
CHECKPOINTS = os.path.join(DATA, "electrum_mainnet_checkpoints.json")

# ---- constants -------------------------------------------------------------
TWH_PER_YR_TO_W = 1e12 / (365.25 * 24.0)   # Eq. (1): TWh per year -> W
L_SUN = 3.828e26                           # W, IAU 2015 nominal solar luminosity
TYPE_I_1964 = 4.0e12                       # W, Kardashev's 1964 benchmark
TYPE_II_1964 = 4.0e26                      # W, Kardashev's 1964 Type II value
T0 = 1964                                  # time origin: t = year - 1964
KARDASHEV_RATE = 0.01                      # the one-percent assumption under test
YEAR_S = 365.25 * 24 * 3600                # Julian year in seconds
KM_PER_MPC = 3.0857e19
H0 = 70.0                                  # km/s/Mpc, round value
HUBBLE_TIME_YR = KM_PER_MPC / H0 / YEAR_S  # 1/H0 = 13.97 Gyr, "about 14 Gyr"
AGE_UNIVERSE_YR = 13.8e9                   # yr
HABITABILITY_CEILING_W = (3e14, 1e15)      # Balbi and Lingam (2025), K_max 0.85-0.9


def ensure_out(*parts):
    """Create outputs/<parts> if needed and return the path."""
    path = os.path.join(OUT, *parts)
    os.makedirs(os.path.dirname(path) if os.path.splitext(path)[1] else path,
                exist_ok=True)
    return path


# ---- energy ---------------------------------------------------------------
def load_energy_rows(path=OWID_CSV):
    """Return {year: {column: float}} for the World rows, 1965-2024."""
    rows = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            y = int(r["year"])
            rows[y] = {k: (float(v) if v not in ("", None) else 0.0)
                       for k, v in r.items() if k not in ("country", "year")}
    return rows


def load_energy(path=OWID_CSV, column="primary_energy_consumption"):
    """Years (int array) and power P in W for 1965-2024 (Eq. 1)."""
    rows = load_energy_rows(path)
    years = np.array(sorted(rows))
    P = np.array([rows[y][column] for y in years]) * TWH_PER_YR_TO_W
    return years, P


def fractional_sigma(years):
    """Adopted per-year fractional 1-sigma (Section 3.3): 10% in 1965,
    declining linearly to 1.5% in 2000 and constant thereafter."""
    return np.clip(np.interp(years, [1965, 2000], [0.10, 0.015]), 0.015, 0.10)


# ---- hashrate -------------------------------------------------------------
def load_coinmetrics(path=CM_CSV):
    """Daily Coin Metrics rows: list of (date, blocks, hashrate in H/s).
    Coin Metrics publishes HashRate in TH/s; empty on days with no block."""
    out = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            n = int(float(r["BlkCnt"])) if r["BlkCnt"] else 0
            h = float(r["HashRate"]) * 1e12 if r["HashRate"] else 0.0
            out.append((r["time"], n, h))
    return out


def annual_hashrate(years=range(2009, 2025)):
    """Average hashrate over each calendar year, H/s: the mean of the daily
    Coin Metrics values (a day without a block counts as zero work). For 2009
    the average runs from the genesis block (3 January)."""
    daily = {}
    for date, n, h in load_coinmetrics():
        y = int(date[:4])
        if y in years:
            daily.setdefault(y, []).append((h, n))
    H = np.array([np.mean([h for h, _ in daily[y]]) for y in years])
    blocks = np.array([sum(n for _, n in daily[y]) for y in years])
    return np.array(list(years)), H, blocks


def difficulty_by_height():
    """Function height -> difficulty from Electrum's epoch checkpoints:
    D = (0xFFFF * 2**208) / target, one target per 2016-block epoch."""
    cp = json.load(open(CHECKPOINTS))
    max_target = 0xFFFF * 2 ** 208

    def D(height):
        k = height // 2016
        return 1.0 if k == 0 else max_target / cp[k - 1][1]
    return D


def fmt_sci(x, digits=2):
    """Format x as 'm.mm x 10^e'."""
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    return f"{x / 10 ** e:.{digits}f}e{e:+d}"
