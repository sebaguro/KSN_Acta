#!/usr/bin/env python3
"""
make_owid_extract.py

Rebuilds data/owid_world_energy_1965_2024.csv, the energy data every other
script reads, from Our World in Data's (OWID) own file, and checks that the
rebuilt extract is identical to the one in data/.

The source is owid-energy-data.csv in OWID's GitHub repository
owid/energy-data, at commit e7897aed55b514b8d16e0a3d8f1518789ac3de31
(27 April 2026). This CSV is the only form of the dataset that OWID keeps
under version control: the JSON and XLSX versions are served from OWID's
download site, which always holds the latest release, so a copy downloaded
there later may carry revised values. Download the source at that commit,
for example from

  https://raw.githubusercontent.com/owid/energy-data/e7897aed55b514b8d16e0a3d8f1518789ac3de31/owid-energy-data.csv

(or open https://github.com/owid/energy-data/tree/e7897aed55b514b8d16e0a3d8f1518789ac3de31
and download owid-energy-data.csv), then run, from the top folder of this
repository:

    python3 make_owid_extract.py path/to/owid-energy-data.csv

What the script does:
  1. checks the SHA-256 of the source file against SOURCE_SHA256 below, so
     that it is the file of that commit;
  2. keeps the rows whose country is "World" and whose year is 1965-2024,
     and the sixteen columns listed in COLUMNS, copying every value exactly
     as written in the source (nothing is rounded or converted);
  3. writes them to outputs/owid_world_energy_1965_2024_rebuilt.csv (it
     never overwrites the file in data/);
  4. compares the rebuilt file with data/owid_world_energy_1965_2024.csv,
     byte for byte, and prints the result.

It reads the source file and data/owid_world_energy_1965_2024.csv, writes
one file in outputs/, and uses no network. Needs only the Python standard
library.
"""

import csv
import hashlib
import os
import sys

SOURCE_COMMIT = "e7897aed55b514b8d16e0a3d8f1518789ac3de31"
SOURCE_SHA256 = "266f2e2baad7975351bc9bb4aa061d22b1da9fe4c47d51d2ac6071e01e171f76"

COLUMNS = [
    "country", "year",
    "primary_energy_consumption", "fossil_fuel_consumption",
    "nuclear_consumption", "renewables_consumption", "hydro_consumption",
    "wind_consumption", "solar_consumption", "biofuel_consumption",
    "other_renewable_consumption",
    "nuclear_electricity", "hydro_electricity", "wind_electricity",
    "solar_electricity", "other_renewable_electricity",
]
FIRST_YEAR, LAST_YEAR = 1965, 2024

HERE = os.path.dirname(os.path.abspath(__file__))
EXTRACT = os.path.join(HERE, "data", "owid_world_energy_1965_2024.csv")
REBUILT = os.path.join(HERE, "outputs", "owid_world_energy_1965_2024_rebuilt.csv")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: python3 make_owid_extract.py path/to/owid-energy-data.csv")
    source = sys.argv[1]

    # 1. Is this the source file of the pinned commit?
    digest = sha256(source)
    print(f"source: {source}")
    print(f"  SHA-256 {digest}")
    if digest == SOURCE_SHA256:
        print(f"  matches owid-energy-data.csv at OWID commit {SOURCE_COMMIT[:7]}")
    else:
        print(f"  DOES NOT match owid-energy-data.csv at OWID commit {SOURCE_COMMIT[:7]}"
              f" (expected {SOURCE_SHA256}); the rebuilt extract may differ")

    # 2. The World rows for 1965-2024, the sixteen columns, values as written.
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

    # 3. Write the rebuilt extract (Python's csv module, default settings).
    os.makedirs(os.path.dirname(REBUILT), exist_ok=True)
    with open(REBUILT, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(COLUMNS)
        writer.writerows(rows)
    print(f"wrote {os.path.relpath(REBUILT, HERE)}")
    print(f"  SHA-256 {sha256(REBUILT)}")

    # 4. Compare with the extract the other scripts read.
    print(f"data/owid_world_energy_1965_2024.csv")
    print(f"  SHA-256 {sha256(EXTRACT)}")
    with open(REBUILT, "rb") as a, open(EXTRACT, "rb") as b:
        same = a.read() == b.read()
    print("IDENTICAL: the extract in data/ is exactly the rebuilt one" if same
          else "DIFFERENT: the extract in data/ is not the rebuilt one")


if __name__ == "__main__":
    main()
