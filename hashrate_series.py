#!/usr/bin/env python3
"""
hashrate_series.py
------------------
The Bitcoin hashrate series H(t) of Section 2.1 and Figure 4.

A block mined at difficulty D requires on average D * 2^32 evaluations of the
hash function (exactly D * 2^48 / 65535, which is 2^32 to within 0.0015
percent), so the average hashrate over any period is 2^32 times the summed
difficulty of the blocks mined in it, divided by its length.

This script
  1. rebuilds that quantity day by day from the difficulty record: the
     difficulty of every 2016-block epoch from Electrum's checkpoints, and the
     blocks mined each day from Coin Metrics (column BlkCnt);
  2. compares it with the daily hashrate Coin Metrics publishes (HashRate),
     showing that the two are the same quantity;
  3. prints the annual means used in the paper, and the statistical
     uncertainty that the random timing of blocks puts on each of them.

Usage:  python3 hashrate_series.py
Writes outputs/hashrate_daily_rebuild.csv.
"""

import csv
import os

import numpy as np

from ksn_common import load_coinmetrics, difficulty_by_height, annual_hashrate, OUT


def main():
    D = difficulty_by_height()
    height = 0                    # the genesis block (height 0) is not in BlkCnt
    rows = []
    for date, n, h_cm in load_coinmetrics():
        first, last = height + 1, height + n
        work = sum(D(k) for k in range(first, last + 1)) * 2 ** 32 / 86400.0
        rows.append((date, n, first if n else "", last if n else "", work, h_cm))
        height = last
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "hashrate_daily_rebuild.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["date", "blocks", "first_height", "last_height",
                    "rebuilt_H_per_s", "coinmetrics_H_per_s"])
        w.writerows(rows)

    days = [r for r in rows if r[1] > 0]
    rel = np.array([abs(r[4] / r[5] - 1) for r in days if r[5] > 0])
    print(f"Days with blocks, 2009-2024: {len(days)}; final height {height}")
    print(f"Rebuilt daily hashrate vs Coin Metrics: largest relative difference "
          f"{rel.max():.1e}; days agreeing to 1e-6: {int((rel < 1e-6).sum())}")

    years, H, blocks = annual_hashrate()
    print("\nAnnual average hashrate (mean of the daily values; 2009 from 3 January)")
    print("year  blocks   H (H/s)     1/sqrt(blocks)")
    for y, h, n in zip(years, H, blocks):
        print(f"{y}  {n:6d}  {h:.3e}   {100 / np.sqrt(n):.2f}%")
    print(f"\nH(2009) = {H[0] / 1e6:.2f} MH/s, H(2024) = {H[-1]:.2e} H/s, span "
          f"{np.log10(H[-1] / H[0]):.1f} dex   [paper: ~4.5 MH/s, ~6.3e20 H/s, 14 orders]")
    print("The random timing of blocks gives each annual mean a statistical uncertainty"
          " of 1/sqrt(blocks), under one percent   [paper, Figure 4: of order one percent]")


if __name__ == "__main__":
    main()
