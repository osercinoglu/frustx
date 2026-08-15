# FrustX

A tool for calculation of local frustration levels in protein structures at atomistic
resolution.

Reimplements the protocol of Chen et al., *Nat Commun* 11:5944 (2020),
["Surveying biomolecular frustration at atomic resolution"](https://pubmed.ncbi.nlm.nih.gov/33230150/),
which scores frustration with the all-atom Rosetta REF2015 energy function rather than a
coarse-grained one.

See [`docs/method.md`](docs/method.md) for the method spec and
[`CLAUDE.md`](CLAUDE.md) for how the repo is organised.

**Status: early, not yet validated per-contact.** Full pipeline runs end to end.
Against frustratometeR (AWSEM) on 1UBQ: residue-level agreement is ρ = +0.41 (p = 2e-4,
n = 76) vs `singleresidue` mode, and per-contact ρ = +0.32 (p = 1e-9, n = 344) vs
`mutational` mode. The per-contact comparison against `configurational` mode has been
withdrawn — that mode uses one global decoy mean/sd for the whole protein, so its index is
not a per-contact Z-score.

The deeper limitation: frustratometeR's indices are largely residue-additive (56%
configurational, 95% mutational), while FrustX's is 23%. **frustratometeR cannot validate
per-contact specificity, which is the property FrustX exists to provide.** Reproducing a
figure from Chen et al. (2020) is the remaining test. See `docs/method.md`.

## Setup on a new machine

Requires Python ≥ 3.11 and ~4 GB of disk (PyRosetta is 3.2 GB).

```bash
git clone https://github.com/osercinoglu/frustx.git
cd frustx

python3 -m venv .venv
.venv/bin/python -m pip install -e ".[dev]"

# PyRosetta is not on PyPI and cannot be a declared dependency.
# The PATH prefix is REQUIRED: pyrosetta_installer shells out to a bare `pip`, so
# without it the 3.2 GB lands in whichever python `pip` resolves to -- not the venv.
PATH=".venv/bin:$PATH" .venv/bin/python -m pip install pyrosetta-installer
PATH=".venv/bin:$PATH" .venv/bin/python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(serialization=True)"

.venv/bin/python -m pytest tests/ -q     # 40 tests; all but tests/test_contacts.py need PyRosetta
```

### Data (DVC + GCS)

No input or output data is in git. Structures, decoy runs and result tables live in
`gs://frustx` and are tracked by `results.dvc`.

DVC is intentionally not in `pyproject.toml` — it is data infrastructure, not something
`frustx` imports, and it would make `pip install -e ".[dev]"` ~200 MB heavier:

```bash
.venv/bin/python -m pip install "dvc[gs]"
```

Authenticate with a GCP service-account key that has `roles/storage.objectAdmin` on the
bucket. **Keep the key outside the repo** and readable only by you:

```bash
mkdir -p ~/.gcp && chmod 700 ~/.gcp
# copy your key to ~/.gcp/<key>.json, then:
chmod 600 ~/.gcp/<key>.json
.venv/bin/dvc remote modify --local storage credentialpath ~/.gcp/<key>.json
```

`--local` is not optional: it writes `.dvc/config.local`, which `.dvc/.gitignore`
excludes. Without it the path goes into the tracked `.dvc/config`.

Then fetch the data:

```bash
.venv/bin/dvc pull          # populates results/
.venv/bin/dvc status -c     # "in sync" = local and bucket agree
```

After producing new results, `dvc add results/ && dvc push`, then commit the updated
`results.dvc`. Alternatively `gcloud auth application-default login` works instead of a
key file, if you have the SDK installed.

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
| `--background-weight` | weight `w` on the Eq. 2 background term, default `0` (direct pair energy only). `w=1` is Eq. 2 as literally written, which makes the index degenerate — see `docs/method.md` |

Colour the output in PyMOL with `spectrum b, blue_white_red, all`. High positive =
minimally frustrated, matching frustratometeR's convention (the paper's Eq. 1 uses the
opposite sign; see `docs/method.md`).

Roughly 12 min per 1000 decoys on a 76-residue protein with `--protocol min`, ~1 h with
`relax`.
