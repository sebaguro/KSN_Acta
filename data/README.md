# Data

Three public data sets, as used in the paper. Each file is an extract of its
source: rows and columns are selected, and no value is changed.

| File | Source | Version | Licence |
|---|---|---|---|
| `owid_world_energy_1965_2024.csv` | Our World in Data (OWID) energy dataset, [owid/energy-data](https://github.com/owid/energy-data), file `owid-energy-data.csv`: the rows for World, 1965-2024 | commit `e7897aed55b514b8d16e0a3d8f1518789ac3de31`, 27 April 2026 | CC BY 4.0 for OWID's work; the underlying data keep their providers' terms (see Licences) |
| `coinmetrics_btc_2009_2024.csv` | Coin Metrics community network data, [coinmetrics/data](https://github.com/coinmetrics/data), file `csv/btc.csv`: columns `time`, `BlkCnt` and `HashRate`, 3 January 2009 to 31 December 2024 | commit `f1a36afb962731c387bb03982758ab0103063da5`, 24 May 2026 | CC BY-NC 4.0 |
| `electrum_mainnet_checkpoints.json` | Electrum Bitcoin wallet, [spesmilo/electrum](https://github.com/spesmilo/electrum), file `electrum/chains/mainnet/checkpoints.json` | commit `e11c579b0461bb333b806bb41e13eecd47484f1a`, 11 September 2026 | MIT (`LICENCE-electrum.txt`) |

## Columns

**OWID** (`owid_world_energy_1965_2024.csv`). `primary_energy_consumption` is
world primary energy in TWh per year, in the input-equivalent (substitution)
convention of the Energy Institute: electricity generated from nuclear, hydro,
wind and solar sources is counted as the fossil fuel that would have been
needed to generate it. The other `*_consumption` columns are its components in
the same convention, and the `*_electricity` columns give the electricity
generated, in TWh. OWID compiles these data from the Energy Institute
Statistical Review of World Energy, the US Energy Information Administration
(primary energy) and Ember (electricity generation), as its codebook for this
commit states. For 1965-2024 the World series in this commit is identical to
that of the earlier release of 28 January 2026 (commit `fca8f40`).

**Coin Metrics** (`coinmetrics_btc_2009_2024.csv`). `BlkCnt` is the number of
blocks mined in each UTC day. `HashRate` is the mean hashrate of the day in
TH/s: the number of hash evaluations implied by the difficulty of the day's
blocks, divided by the length of the day. It is empty on days without a block.
When this extract was made (28 September 2026), the file had not changed since
the commit above.

**Electrum** (`electrum_mainnet_checkpoints.json`). One entry per 2016-block
difficulty epoch: the hash of the epoch's last block and the target of the
epoch. The difficulty is D = (0xFFFF × 2^208) / target.

## Rebuilding the OWID extract

The OWID extract can be rebuilt from OWID's own file. Download
`owid-energy-data.csv` at commit `e7897aed55b514b8d16e0a3d8f1518789ac3de31`, for
example from
https://raw.githubusercontent.com/owid/energy-data/e7897aed55b514b8d16e0a3d8f1518789ac3de31/owid-energy-data.csv
(9.2 MB; SHA-256
`266f2e2baad7975351bc9bb4aa061d22b1da9fe4c47d51d2ac6071e01e171f76`), and run,
from the top folder of the repository:

```
python3 make_owid_extract.py path/to/owid-energy-data.csv
```

The script checks the SHA-256 of the file, keeps the rows for World from 1965
to 2024 and the sixteen columns listed above, copies each value exactly as
written, writes `outputs/owid_world_energy_1965_2024_rebuilt.csv`, and checks
that it is identical, byte for byte, to `owid_world_energy_1965_2024.csv`
(SHA-256 `40c03ed590e1a6c0042eff4f63643389a10fb80de80b0a46215005d86be4ff0b`).
OWID keeps only this CSV under version control; the JSON and XLSX versions of
the dataset are served from its download site, which holds the latest release
only, so they are not used here.

## Licences

The licence of the code in this repository does not extend to these files;
each keeps the terms of its source.

- **OWID.** Data and code produced by Our World in Data are licensed under
  Creative Commons Attribution 4.0 (CC BY 4.0,
  https://creativecommons.org/licenses/by/4.0/). OWID states that data
  produced by third parties remain subject to the terms of their original
  providers, and asks users to cite both OWID and the underlying sources,
  credited below: the Energy Institute (EI), the US Energy Information
  Administration (EIA), whose data are in the public domain, and Ember,
  whose data are licensed CC BY 4.0. This file is a small extract, kept only
  so that the results of the paper can be reproduced.
- **Coin Metrics.** Creative Commons Attribution-NonCommercial 4.0 (CC BY-NC
  4.0, https://creativecommons.org/licenses/by-nc/4.0/): the extract may be
  used and shared with attribution, but not for commercial purposes. Coin
  Metrics provides the data without warranty.
- **Electrum.** MIT licence. The copyright and permission notice that the
  licence requires to accompany copies is in `LICENCE-electrum.txt`, as
  published in the Electrum repository at the commit above.

## Credits

- Hannah Ritchie, Pablo Rosado and Max Roser, "Energy", Our World in Data,
  https://ourworldindata.org/energy.
- Energy Institute, Statistical Review of World Energy 2025, Energy
  Institute, London, https://www.energyinst.org/statistical-review/.
- U.S. Energy Information Administration, International Energy Data (2025),
  https://www.eia.gov/opendata/bulkfiles.php.
- Ember, Yearly Electricity Data (2026),
  https://ember-energy.org/data/yearly-electricity-data/.
- Coin Metrics, Coin Metrics Community Network Data,
  https://github.com/coinmetrics/data.
- The Electrum developers and Thomas Voegtlin, Electrum,
  https://github.com/spesmilo/electrum.
