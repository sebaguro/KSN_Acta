# Kardashev's Conundrum: data and code

This repository holds the data and the scripts behind the paper

> S. Gurovich, "Kardashev's Conundrum: Statistical Falsification of the
> Standard One-Percent Kardashev Growth Model, the Non-Observation of Type II
> Technosignatures, and the Kardashev-Sagan-Nakamoto Resolution", submitted to
> Acta Astronautica (2026).

Every figure of the paper, and every number in it that is computed from data,
is reproduced here from the copies of the data kept in `data/`.

## Running

The scripts need Python 3 with numpy, scipy and matplotlib; `requirements.txt`
gives the versions used. They also run with NumPy 1.x: with NumPy 1.26.4,
SciPy 1.13.1 and Matplotlib 3.8.4 every printed result is identical.

```
pip install -r requirements.txt
sh run_all.sh
```

`run_all.sh` runs every script, writing the printed results to
`outputs/<script>.txt` and the figures to `outputs/figures/`. Each printed
quantity is followed by the value the paper states, as `[paper: ...]`. The run
takes under a minute. The files in `outputs/` were produced in this way, with
Python 3.11 and the versions in `requirements.txt`.

## Where each result comes from

| Paper | Script |
|---|---|
| Section 2.1: Eq. (1), P(1965) and P(2024) | `energy_fits.py` |
| Section 2.1: the hashrate H(t), checked against the difficulty record | `hashrate_series.py` |
| Section 2.2: the one-percent model | `energy_fits.py` |
| Section 2.3: the linear model, the Shapiro-Wilk test on ΔP, Eq. (5) | `energy_fits.py` |
| Section 2.4: the posterior growth rate (MCMC) and the WAIC comparison with its standard error | `waic_indistinguishability.py` |
| Section 2.4: the least-squares exponential and the habitability ceiling | `energy_fits.py` |
| Section 2.4: growth laws between the linear and the exponential | `intermediate_forms.py` |
| Section 2.5: N_b and B_min | `nb_gate_count.py`, `landauer_estimate.py` |
| Section 2.6: f, ε_hw, the hardware factor and the time to the practical floor | `landauer_estimate.py` |
| Section 2.7: serial dependence of the growth rates | `energy_fits.py` |
| Section 3.1: the habitability ceiling relative to L_sun and to present power | `energy_fits.py` |
| Section 3.2, criterion 3 | `landauer_estimate.py` |
| Section 3.3: the energy-accounting conventions | `convention_refit.py` |
| Section 3.3: the weighted re-analysis | `energy_fits.py`, `waic_indistinguishability.py --weighted` |
| Table 1 | `energy_fits.py`; the KSN row, `landauer_estimate.py` |
| Table 2 | `energy_fits.py` |
| Figures 1 and 2 | `make_errorbar_figures.py` |
| Figure 3 | `make_fig3_deltaP.py` |
| Figure 4, and the annual series of H, P and B | `make_fig4_ksn.py` |

`ksn_common.py` holds the data loaders and the constants shared by the
scripts. `hashrate_series.py` also writes the daily rebuild of the hashrate
(`outputs/hashrate_daily_rebuild.csv`), and `make_fig4_ksn.py` the annual
series (`outputs/ksn_annual_series.csv`).

`nb_gate_count.py` builds the double SHA-256 of Bitcoin mining from two-input
XOR, AND and OR gates and checks its digests against Python's `hashlib` for
16,384 nonces, one of them the nonce of block 125552, whose hash it
reproduces. It counts 240,896 gate operations per hash attempt. Of these,
205,828 have inputs that vary with the nonce, and together they destroy an
estimated 199,668 bits of information per attempt (the sum over those gates
of the entropy of the inputs minus that of the output). This is the N_b of
Section 2.5, which `landauer_estimate.py` uses. `run_all.sh` runs the script
with 16,384 samples, the number used for the paper.

## Data

| File | Source | Licence |
|---|---|---|
| `data/owid_world_energy_1965_2024.csv` | Our World in Data, energy dataset: world primary energy, 1965-2024 | CC BY 4.0 for OWID's work; the underlying Energy Institute, EIA and Ember data keep their providers' terms |
| `data/coinmetrics_btc_2009_2024.csv` | Coin Metrics community network data: daily blocks and hashrate of Bitcoin, 2009-2024 | CC BY-NC 4.0 |
| `data/electrum_mainnet_checkpoints.json` | Electrum Bitcoin wallet: the difficulty target of every 2016-block epoch | MIT |

Each file is an extract of a fixed version of its source; `data/README.md`
gives the versions, the columns, the licences and the credits. The scripts
read only these copies, so the results do not change when the sources are
revised. The licence of the code does not extend to the data files, and the
Coin Metrics extract may not be used commercially.

## Notes

- Time is measured from 1964, t = year - 1964, and power is converted from
  TWh per year with Eq. (1): P = E × 10^12 / (365.25 × 24) W.
- In Section 2.2 the one-percent model starts from the production at t = 0,
  P0 = P(1965) / 1.01, as in Eq. (2). `energy_fits.py` also prints, as a
  diagnostic, the ratio of the residual sums of squares when P0 is fitted
  instead.
- The per-year uncertainties of the weighted re-analysis (Section 3.3), drawn
  as error bars in Figures 1 and 2, are an adopted model, not published
  values: a fractional 1σ of 10 percent in 1965, falling linearly to 1.5
  percent in 2000 and constant thereafter.
- The Hubble time is 1/H0 with H0 = 70 km/s/Mpc, that is 13.97 Gyr.
- The MCMC sampler uses a fixed random seed, so its results repeat exactly.

Earlier versions of this repository held the LaTeX source of the submitted
manuscript; it remains in the history.

## Licence

The code in this repository is released under the BSD 3-Clause licence (see
`LICENSE`). The files in `data/` keep the licences of their sources, given in
`data/README.md`. The files in `outputs/` are produced by the scripts from
those data; where an output contains values from a data set, that data set's
terms apply to them (for example the Coin Metrics hashrate in
`outputs/hashrate_daily_rebuild.csv`, CC BY-NC 4.0).

## Citation

If you use these data or scripts, please cite the paper above and the sources
of the data listed in `data/README.md`.
