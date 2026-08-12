# FrustX

A tool for calculation of local frustration levels in protein structures at atomistic
resolution.

Reimplements the protocol of Chen et al., *Nat Commun* 11:5944 (2020),
["Surveying biomolecular frustration at atomic resolution"](https://pubmed.ncbi.nlm.nih.gov/33230150/),
which scores frustration with the all-atom Rosetta REF2015 energy function rather than a
coarse-grained one.

See [`docs/method.md`](docs/method.md) for the method spec and
[`CLAUDE.md`](CLAUDE.md) for how the repo is organised.

**Status: early.** Contact-map layer only.

## Setup

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e ".[dev]"
# PATH prefix is required -- the installer calls a bare `pip`, not the venv's.
PATH=".venv/bin:$PATH" .venv/bin/python -m pip install pyrosetta-installer
PATH=".venv/bin:$PATH" .venv/bin/python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(serialization=True)"
.venv/bin/python -m pytest tests/ -q
```
