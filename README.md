# FrustX

A tool for calculation of local frustration levels in protein structures at atomistic
resolution.

Reimplements the protocol of Chen et al., *Nat Commun* 11:5944 (2020),
["Surveying biomolecular frustration at atomic resolution"](https://pubmed.ncbi.nlm.nih.gov/33230150/),
which scores frustration with the all-atom Rosetta REF2015 energy function rather than a
coarse-grained one.

See [`docs/method.md`](docs/method.md) for the method spec and
[`CLAUDE.md`](CLAUDE.md) for how the repo is organised.

**Status: early.** Full pipeline runs end to end; not yet validated against a reference implementation.

## Setup

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[dev]"
# PATH prefix is required -- the installer calls a bare `pip`, not the venv's.
PATH=".venv/bin:$PATH" .venv/bin/python -m pip install pyrosetta-installer
PATH=".venv/bin:$PATH" .venv/bin/python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(serialization=True)"
.venv/bin/python -m pytest tests/ -q
```

## Usage

```bash
frustx structure.pdb -o results/
```

Writes `contacts.csv` (one row per contact, the primary result), `residues.csv`
(aggregated per residue), `frustration.pdb` (indices painted into the B-factor column)
and `run.json` (full provenance — results are only comparable between runs with the same
decoy count and relaxation protocol).

Useful flags:

| flag | what |
|---|---|
| `-n, --decoys` | decoy ensemble size (default 1000, the paper's figure) |
| `--protocol` | `relax` (default, FastRelax) or `min` (~10× faster) — this sets Eq. 1's denominator, so it moves every index |
| `--packing-frustration` | keep `fa_rep`, giving the paper's separate structure-quality measure |
| `--cutoff` | Cα–Cα contact cutoff, default 10 Å |

Colour the output in PyMOL with `spectrum b, blue_white_red, all`. High positive =
minimally frustrated, matching frustratometeR's convention (the paper's Eq. 1 uses the
opposite sign; see `docs/method.md`).

Roughly 12 min per 1000 decoys on a 76-residue protein with `--protocol min`, ~1 h with
`relax`.
