#!/usr/bin/env python3
"""
make_owid_tes_extract.py

Rebuilds data/owid_world_total_energy_supply_1965_2024.csv, the second
energy extract that convention_refit.py reads, from Our World in Data's
(OWID) own file, and checks that the rebuilt extract is identical to the one
in data/.

The source is owid_energy.csv, the OWID energy dataset in its release of
10 September 2026. In that release OWID measures primary energy as total
energy supply and has renamed its data columns: the world total is
total_energy_twh. OWID's notes on the release give two addresses for the
file, one that always holds the latest release and one that keeps this
release:

  https://catalog.ourworldindata.org/energy/owid_energy/owid_energy.csv
  https://catalog.ourworldindata.org/garden/energy/2026-09-10/owid_energy/owid_energy.csv

The copy used here was downloaded from the first address on 4 October 2026,
when the latest release was that of 10 September 2026 (9,081,888 bytes); the
second address was not used. Download the file, then run, from the top
folder of this repository:

    python3 make_owid_tes_extract.py path/to/owid_energy.csv

What the script does:
  1. checks the SHA-256 of the source file against SOURCE_SHA256 below, so
     that it is the file used here;
  2. keeps the rows whose country is "World" and whose year is 1965-2024,
     and the eleven columns listed in COLUMNS, copying every value exactly
     as written in the source (nothing is rounded or converted);
  3. compares the rebuilt extract with
     data/owid_world_total_energy_supply_1965_2024.csv, byte for byte, and
     prints the result.

It reads the source file and the extract in data/, writes nothing unless a
second path is given (the rebuilt extract is then saved there; it never
overwrites the file in data/), and uses no network. Needs only the Python
standard library.
"""

import csv
import hashlib
import io
import os
import sys

SOURCE_RELEASE = "10 September 2026"
SOURCE_SHA256 = "9d8efc4a9b1cc2ffdfa7f7f6f1e171b2922e0be2d97c0edf2982e874357ed6de"

COLUMNS = [
    "country", "year",
    "total_energy_twh",
    "nuclear_energy_twh", "nuclear_electricity_twh",
    "hydro_energy_twh", "hydro_electricity_twh",
    "wind_energy_twh", "wind_electricity_twh",
    "solar_energy_twh", "solar_electricity_twh",
]
FIRST_YEAR, LAST_YEAR = 1965, 2024

HERE = os.path.dirname(os.path.abspath(__file__))
EXTRACT = os.path.join(HERE, "data", "owid_world_total_energy_supply_1965_2024.csv")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    if len(sys.argv) not in (2, 3):
        sys.exit("usage: python3 make_owid_tes_extract.py path/to/owid_energy.csv"
                 " [path/to/save/the/rebuilt/extract.csv]")
    source = sys.argv[1]
    save_to = sys.argv[2] if len(sys.argv) == 3 else None
    if save_to and (os.path.realpath(save_to) == os.path.realpath(EXTRACT)
                    or (os.path.exists(save_to) and os.path.exists(EXTRACT)
                        and os.path.samefile(save_to, EXTRACT))):
        sys.exit("the rebuilt extract is never written over the file in data/")

    # 1. Is this the source file used here?
    digest = sha256_file(source)
    print(f"source: {os.path.basename(source)}")
    print(f"  SHA-256 {digest}")
    if digest == SOURCE_SHA256:
        print(f"  matches owid_energy.csv, OWID release of {SOURCE_RELEASE}, as used here")
    else:
        print(f"  DOES NOT match the owid_energy.csv used here (OWID release of"
              f" {SOURCE_RELEASE}; expected {SOURCE_SHA256}); the rebuilt extract may differ")

    # 2. The World rows for 1965-2024, the eleven columns, values as written.
    rows = []
    with open(source, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        missing = [c for c in COLUMNS if c not in reader.fieldnames]
        if missing:
            sys.exit(f"source lacks columns: {missing}")
        for row in reader:
            if row["country"] == "World" and FIRST_YEAR <= int(row["year"]) <= LAST_YEAR:
                rows.append([row[c] for c in COLUMNS])
    rows.sort(key=lambda r: int(r[1]))
    years = [int(r[1]) for r in rows]
    if years != list(range(FIRST_YEAR, LAST_YEAR + 1)):
        sys.exit(f"expected one row per year {FIRST_YEAR}-{LAST_YEAR}, found {len(rows)}")
    print(f"kept {len(rows)} rows (World, {FIRST_YEAR}-{LAST_YEAR}) and {len(COLUMNS)} columns")

    # The rebuilt extract, as Python's csv module writes it (default settings).
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer)
    writer.writerow(COLUMNS)
    writer.writerows(rows)
    rebuilt = buffer.getvalue().encode("utf-8")
    print(f"rebuilt extract: {len(rebuilt)} bytes")
    print(f"  SHA-256 {hashlib.sha256(rebuilt).hexdigest()}")

    # 3. Compare with the extract that convention_refit.py reads.
    print("data/owid_world_total_energy_supply_1965_2024.csv")
    if os.path.exists(EXTRACT):
        print(f"  SHA-256 {sha256_file(EXTRACT)}")
        with open(EXTRACT, "rb") as f:
            same = f.read() == rebuilt
        print("IDENTICAL: the extract in data/ is exactly the rebuilt one" if same
              else "DIFFERENT: the extract in data/ is not the rebuilt one")
    else:
        print("  not found")

    if save_to:
        with open(save_to, "wb") as f:
            f.write(rebuilt)
        print(f"rebuilt extract saved to {save_to}")


if __name__ == "__main__":
    main()
