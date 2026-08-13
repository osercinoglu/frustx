# The method FrustX implements

Source: Chen M, Chen X, Schafer NP, Clementi C, Komives EA, Ferreiro DU, Wolynes PG.
**"Surveying biomolecular frustration at atomic resolution."**
*Nat Commun* 11:5944 (2020). PMID 33230150 / PMC7683549 / doi:10.1038/s41467-020-19560-9.

Quotations below are verbatim from the paper's Methods.

## 1. Frustration index (Eq. 1)

```
F°_ij = ( E⁰_ij − E^U_{i'j'} ) / σ( E^U_{i'j'} )
```

`E⁰_ij` is the contact energy of residues i,j in the native sequence; `E^U_{i'j'}` is the
same quantity over the decoy ensemble. A Z-score: the native contact energy measured
against the mean and spread of what that same geometric site would give for random
sequence.

## 2. Contact energy (Eq. 2) — the many-body term

```
E_ij = e_ij + ½ Σ_{k≠j} e_ik + ½ Σ_{l≠i} e_jl
```

Not the bare pair energy. `e_ij` is the direct interaction, plus **half of every other
interaction involving i**, plus **half of every one involving j**. Per the paper:

> "Since the fully atomistic force field is a many-body construct, the pairwise energy
> change assigned to forming a contact between residue i and j, Eij is defined by
> considering all the interaction energies that involve changing any of the two residues
> that are in contact."

> "In addition to the direct interactions between i and j (defined as eij), the addition
> of the background interactions accounts for the many-body effects elicited by local
> side-chain reconfigurational changes."

**This is the step the earlier `FrustX.py` prototype was missing** — it used the total
system potential energy, in which a single side-chain change is far below the noise.

## 3. Energy function

**Rosetta REF2015**, all-atom. *"we employed the REF2015 version of the rosetta energy
function, which has a set of well-tested weights for each energy term."*

**The repulsive Lennard-Jones term is removed** when computing `e_ij`:

> "Even after relaxation, the repulsive Lennard–Jones interactions used in computing eij
> give rise to very strong fluctuating clashes... This gives rise to an artificially large
> variance of the decoy energies. These repulsive clashes would make the algorithm too
> sensitive for detecting functionally relevant frustration... it is both easy and
> appropriate to simply remove the harsh rapidly varying repulsive force term."

Keeping `fa_rep` in gives a distinct quantity the authors name **packing frustration**,
which *"turns out to be good at diagnosing the health of a given high-resolution
structure (i.e., predicted structure)"* — i.e. useful for judging predicted models.

## 4. Decoy generation

> "the decoy ensemble of residue contacts is obtained by randomizing both the residue
> identities and locations of the residues in contact"

> "we randomly shuffle the protein sequence and then repack the resulting sequence onto
> the backbone that is provided without perturbing the backbone coordinates within each
> protein chain to make sure that only the side chains are re-packed. A short Monte-Carlo
> relaxation is then performed to better eliminate many of the possible side-chain clashes
> with the backbone fixed."

Native reference energies: *"obtained in a similar fashion by omitting the shuffling
step."*

Key consequences for implementation:

- **Shuffle, not mutate.** A random permutation of the native sequence automatically
  satisfies *"sequence space is randomly sampled according to the native amino acid
  frequency distribution"* — composition is preserved for free.
- **Backbone stays fixed throughout.** Side chains only.
- **1000 decoys per contact.** One shuffle+repack yields an energy for *every* contact at
  once, so this is ~1000 shuffle-repack cycles for the whole protein, not
  1000 × n_contacts.

## 5. Contact definition

> "Protein contacts are defined by the CαCα distances between residues, and a cutoff of
> 10 Å is used."

Implemented in `frustx/contacts.py`.

## Open questions to settle empirically (Phase 6)

1. **Sign convention.** Eq. 1 as written makes a minimally frustrated contact *negative*
   (native energy below the decoy mean). `frustratometeR` reports minimally frustrated as
   `F > 0.78`. One of the two flips a sign; resolve by running both on the same protein
   rather than by guessing.
2. **Sequence separation.** The Methods state a distance cutoff and nothing else, so we
   default to keeping every pair (`min_seq_sep=1`). Other frustration tools exclude
   near-sequence neighbours, whose contact is dictated by covalent geometry rather than by
   packing. On ubiquitin the very first contact found is MET1–GLN2 at 3.8 Å, which shows
   the question is not academic. Settle against reference output.
3. **Repack radius during decoy generation.** "Repack the resulting sequence onto the
   backbone" reads as a global repack; whether restricting to a shell around each contact
   changes the answer is worth measuring, since it dominates runtime.

## Deliberately out of scope for now

- Protein–ligand / binding-pocket frustration (the paper's EGFR and COX sections).
- OpenMM as an alternative energy backend. It cannot reproduce these numbers — different
  force field, different decoy semantics. It is a separate experiment, not a prerequisite.
- Packaging, Docker, web server.

---

# Implementation notes: extracting `e_ij` from Rosetta

Verified empirically on 1UBQ (ubiquitin, 76 residues) with PyRosetta
2026.32 / REF2015. Probe scripts were run in scratch, not committed; the numbers below
are reproducible from a bare `pose_from_pdb` + `get_score_function()`.

## The naive energy graph is NOT the whole energy

Eq. 2 assumes the energy decomposes into pairwise `e_ij`. Reading only the `EnergyGraph`
edges plus one-body terms **silently loses 49.26 REU** on ubiquitin:

| term | weighted energy off the graph |
|---|---|
| `hbond_lr_bb` | −23.13 REU |
| `hbond_sr_bb` | −18.83 REU |
| `rama_prepro` | −7.30 REU |

Losing the backbone–backbone hydrogen bonds would be fatal: they are exactly the
interactions that stabilise secondary structure, i.e. the thing frustration analysis is
supposed to be measuring.

## Fix 1 — put bb–bb hbonds on the graph

```python
opts = sf.energy_method_options()
opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
sf.set_energy_method_options(opts)
```

This is a score-function option, **not** a command-line flag (`-score:decompose_bb_hb...`
does not exist and raises at init). Control: the weighted total is 32.6778 REU with the
option on *and* off — it redistributes where energy is stored, it does not change the
physics.

## Fix 2 — read `rama_prepro` from its long-range container

`rama_prepro` lives in a `PolymerBondedEnergyContainer`
(`LongRangeEnergyType.ramaprepro_lr`), reachable via
`pose.energies().long_range_container(...)`. On ubiquitin it yields 75 pair records, all
strictly (i, i+1), summing to the pose total exactly (−16.2164, difference 0.0000).

Iterate with `it.pre_increment()` — the PyRosetta binding has no `.next()`.

**`rama_prepro` cannot be ignored as sequence-independent.** With the backbone held fixed,
single mutations at V5 of ubiquitin move it by:

| mutation | Δ rama_prepro |
|---|---|
| V5A | +2.89 |
| V5W | +2.58 |
| V5G | +10.35 |
| V5P | +12.51 |

So it does not cancel in the Z-score under sequence shuffling. It is a genuine (i, i+1)
two-body term and Eq. 2 can carry it.

## Result: the decomposition is exact

`one-body + EnergyGraph edges + long-range containers` reproduces the REF2015 total to
**0.00e+00 REU**. Eq. 2 is exactly implementable with no unaccounted energy.

Sweep long-range containers generically rather than hard-coding `ramaprepro_lr`:
`dslf_fa13` appears the same way for disulfide-bonded proteins, and ubiquitin has no
disulfides so it would not have shown up in this test.

## Caveat on raw input structures

A crystal structure scores **+32.7 REU** (unfavourable) before any relaxation, because
`fa_rep` sees clashes that crystallographic refinement tolerates. Never score a bare input
PDB and treat the number as meaningful — this is why the paper specifies a short
Monte-Carlo relaxation, and it previews the `fa_rep` problem solved by deleting the term.

## Validation of the extraction on a real structure

1UBQ, REF2015 with `fa_rep` removed as the paper specifies:

- **Completeness holds on real input**: one-body + pairwise reproduces the REF2015 total
  to 6e-13 REU.
- **The 10 Å Cα cutoff is well chosen**: −323.27 REU of pair energy falls inside the
  contact set, only −10.32 REU outside it. The cutoff captures ~97% of pairwise energy.
- **Magnitudes are sane**: mean e_ij = −0.72 REU over contacts, range −6.26 to +1.24.
- **The strongest contacts are the right ones.** The five most favourable are
  ILE44–HIS68, VAL5–ILE13, MET1–VAL17, ILE3–LEU15, ARG42–VAL70 — the ubiquitin
  hydrophobic core, and Ile44/His68/Val70 are the canonical Ile44 patch that nearly every
  ubiquitin receptor binds. An independent check that the extraction is physically
  meaningful, not merely self-consistent.

### Open question

98 of 550 contacts (18%) have exactly zero e_ij: their Cα atoms are within 10 Å but no
atom pair falls inside Rosetta's interaction cutoff. Whether these should enter the decoy
statistics as genuine zeros or be excluded is not addressed by the paper. They will have
zero variance across decoys where the residues stay far apart, which risks a divide-by-zero
in Eq. 1 — this needs an explicit decision when the frustration index is implemented.

---

# Eq. 2 simplifies: the direct contact term cancels

Expanding Eq. 2 with `R_i = Σ_k e_ik` (residue i's total interaction energy, diagonal
zero):

```
Σ_{k≠j} e_ik = R_i − e_ij
Σ_{l≠i} e_jl = R_j − e_ij

E_ij = e_ij + ½(R_i − e_ij) + ½(R_j − e_ij) = ½ (R_i + R_j)
```

**The `e_ij` term cancels exactly.** Verified numerically against a literal double-loop
transcription of Eq. 2 on 1UBQ: max difference 7.1e-15.

The cancellation is robust to the ambiguity in the summation ranges. Whether `k` runs
over all residues except `j`, or over all except both `i` and `j`, the result is the same,
because `e_ii = 0`.

## What this means

`E_ij` depends only on the two residues' *total* interaction energies — not on how
strongly they touch each other. The paper's contact-level frustration index is therefore
a **pair average of a residue-level quantity**. The contact set enters only by selecting
which (i, j) pairs get reported, not by weighting the direct interaction between them.

This follows from the paper's own stated intent — *"considering all the interaction
energies that involve changing any of the two residues that are in contact"* — so it
appears deliberate rather than a transcription error. It is worth confirming against
`frustratometeR`, whose AWSEM configurational index does not share this property.

Practical consequences:

- **Fast.** Row sums, not a double loop over pairs.
- **The zero-`e_ij` worry is void.** The 98 of 550 ubiquitin contacts with `e_ij = 0`
  raised a divide-by-zero risk in Eq. 1. Since `E_ij` never uses `e_ij`, they get finite
  energies and finite decoy variance like any other pair. The earlier open question is
  closed.

---

# Sign convention: FrustX negates the paper's Eq. 1

The paper writes `F = (E0 − ⟨E_U⟩)/σ`, which makes a **minimally** frustrated contact
come out **negative** — the native energy sits below the decoy mean. Measured on
ubiquitin, the literal formula gives a range of −4.05 to +0.55, with a mean of −0.84:
nearly everything negative, as expected for a well-folded protein.

Every other tool in the field reports the opposite sign. `frustratometeR` classifies
minimally frustrated as `F > 0.78` and highly frustrated as `F < −1`.

**FrustX negates**, so output matches the field:

| FrustX `frustration_index` | meaning |
|---|---|
| high positive | minimally frustrated — native much better than random |
| near zero | neutral |
| negative | highly frustrated — native no better than random |

`tests/test_frustration.py` guards this: a designed native must beat shuffled decoys and
therefore score positive. Anyone "correcting" the sign back to the paper's literal form
breaks that test.

## Classification thresholds are borrowed and provisional

`MINIMALLY_FRUSTRATED = 0.78` and `HIGHLY_FRUSTRATED = −1.0` are taken from
frustratometeR so output is comparable. **They were calibrated on the AWSEM
coarse-grained energy function, not on atomistic REF2015**, and are almost certainly the
wrong cut points here. Recalibrating them is part of validation.

# Known artifact: glycine positions read as minimally frustrated

Whole-sequence shuffling drops bulky residues onto positions where only glycine fits
sterically. Decoy energies there balloon, and the native glycine looks superb by
comparison.

On ubiquitin the C-terminal LRGG tail lands among the *most* minimally frustrated
contacts:

| contact | F (conventional sign) | E0 | ⟨E_decoy⟩ |
|---|---|---|---|
| GLY75–GLY76 | +3.66 | +1.0 | +8.5 |
| ARG74–GLY76 | +3.45 | −0.4 | +4.4 |
| LEU73–GLY76 | +2.96 | −1.5 | +3.3 |

That tail is flexible, solvent-exposed, and the functional conjugation site. It should
read as *frustrated*. What the index is actually reporting is "glycine fits here and
bulky residues do not" — a composition-and-sterics effect, not folding frustration.

**CORRECTED after validation.** This was originally recorded as "inherent to the paper's
shuffling protocol rather than to this implementation". The frustratometeR comparison
disproves that. On the same structure, over 344 shared contacts:

| | GLY contacts (median index) | non-GLY |
|---|---|---|
| FrustX | **+0.39** | +0.55 |
| frustratometeR | **−0.41** | +0.40 |

frustratometeR classifies glycine contacts as *frustrated*; we classify them as
*minimally frustrated*. The signs are opposite, not merely different in magnitude
(PHE45–GLY47: ours +0.57, theirs −0.50).

So this is a defect on our side, not a property of the published protocol. It is very
likely a symptom of the Eq. 2 degeneracy documented below rather than an independent bug:
an index reducible to per-residue burial will read a small, well-fitting glycine as
optimal regardless of what it actually interacts with.

---

# Validation against frustratometeR (1UBQ)

Setup: identical input PDB (frustratometeR's own prepared `1ubq_A.pdb`, 76 residues),
FrustX at 1000 decoys with `--protocol min --min-seq-sep 2` to match their contact set.
frustratometeR at `Mode="configurational"`, defaults.

Note `SeqDist=12` does **not** mean minimum |i−j| ≥ 12, as its name suggests. Their
output contains contacts at |i−j| = 2 (74 of them). `SeqDist` selects the AWSEM
`Welltype` classification — `short` / `long` / `water-mediated` — via the two bundled
LAMMPS binaries. Their effective minimum sequence separation is 2.

## Result: weak agreement

```
FrustX 475 contacts, frustratometeR 389, shared 344
Spearman +0.159   Pearson +0.133   class agreement 56.7%

                 min      median    max      sd
FrustX          -0.18     +0.53    +2.84    0.43
frustratometeR  -2.75     +0.37    +2.89    1.07
```

FrustX finds **zero** highly frustrated contacts in ubiquitin (minimum −0.18, threshold
−1.0). frustratometeR finds 30. A frustration index that never reports frustration is
not functioning.

## Diagnosis: the Eq. 2 cancellation is the cause

Predicting each contact index from the mean of its two residues' indices:

| | ρ |
|---|---|
| **FrustX** | **+0.865** |
| frustratometeR | +0.464 |

Our contact index is 87% reducible to a residue-level quantity; theirs is not. This is
the algebraic cancellation showing up in data: `E_ij = ½(R_i + R_j)` contains no term
describing the specific i–j interaction, so contact-level structure cannot survive.

Corroborating evidence:

- Our per-residue values correlate **better** with their *single-residue* index
  (ρ = +0.348) than our contact values do with their *configurational* index
  (ρ = +0.159). We built a single-residue measure while intending a contact measure.
- By AWSEM well type: ρ = +0.32 (`long`), +0.27 (`short`), but **−0.05
  (`water-mediated`)** — 46% of contacts. Water-mediated contacts are exactly those whose
  energy depends on the specific pair rather than on burial.

## Caveat

AWSEM is not ground truth, and the paper's premise is that atomistic and coarse-grained
frustration differ. Weak correlation alone would not condemn the implementation. But
"no frustrated contacts at all" and "87% reducible to residue identity" are internal
red flags independent of whether frustratometeR is right.

## Experiment: does the Eq. 2 background explain the disagreement?

Two sweeps, 1000 decoys each, same contacts, `protocol="min"`.

**Sweep 1** — `E_ij = ½(R_i + R_j) + α·e_ij` for α ∈ {0, +1, −½}, the three defensible
readings of Eq. 2's summation ranges. All three were indistinguishable
(ρ = 0.147–0.170, zero frustrated contacts, pair-specificity 0.86). The reason is
magnitude: `e_ij` ≈ −0.7 REU against `½(R_i + R_j)` ≈ −10 REU, so every α is a small
perturbation on a dominant term. **The cancellation was real but never the binding
constraint.**

**Sweep 2** — `E_ij = e_ij + w·½(R_i + R_j)`, sweeping the background weight. `w = 1` is
Eq. 2 as written; `w = 0` is the direct pair energy alone, which is structurally what
AWSEM/frustratometeR uses.

```
variant       rho     min   median    max     sd   n<-1   GLYmed   pair-spec
frustratometeR --   -2.75    0.37    2.89   1.07     30    -0.41     +0.464
w=0.0      +0.128   -3.15    0.04    6.17   0.92     11    -0.29     +0.438
w=0.05     +0.176   -1.62    0.55    5.87   0.68      2    +0.25     +0.637
w=0.1      +0.172   -0.80    0.56    5.41   0.59      0    +0.35     +0.716
w=0.25     +0.160   -0.39    0.56    4.17   0.50      0    +0.41     +0.802
w=0.5      +0.157   -0.16    0.56    3.11   0.46      0    +0.43     +0.839
w=1.0      +0.159   -0.11    0.55    2.82   0.44      0    +0.43     +0.855
```

### Two independent conclusions

**1. The background term is what makes the index degenerate.** At `w = 0` every pathology
disappears and the statistical character matches frustratometeR closely:

| | FrustX `w=0` | frustratometeR |
|---|---|---|
| pair-specificity | +0.438 | +0.464 |
| spread (sd) | 0.92 | 1.07 |
| minimum | −3.15 | −2.75 |
| contacts below −1 | 11 | 30 |
| glycine median | **−0.29** | **−0.41** |

The glycine artifact **reverses sign and disappears** at `w = 0`. It was a symptom of the
background term, exactly as suspected — confirming the correction recorded above.

**2. None of this explains the disagreement with frustratometeR.** ρ stays at 0.13–0.18
across the entire sweep, while pair-specificity moves from 0.44 to 0.86. The correlation
is *insensitive to the whole Eq. 2 question*. Whatever drives the disagreement is not the
contact-energy definition.

So: dropping the background gives an index with the right shape, the right spread, real
frustrated contacts and correct glycine behaviour — but it still ranks individual contacts
differently from AWSEM. Right distribution, different assignment.

### What this does not settle

ρ ≈ 0.15 at n = 344 is statistically significant (p ≈ 0.005) but weak. Remaining
candidates, in the order worth testing:

1. **Decoy locality.** We shuffle the entire sequence, so a decoy energy for contact
   (i,j) depends on every other position too. The AWSEM configurational decoy randomises
   the identities and geometry *of the contacting pair*, holding the rest native — a far
   more local perturbation. The atomistic paper does specify whole-sequence shuffling, but
   the two are not statistically equivalent.
2. **Relaxation.** These runs used `protocol="min"`. The default `"relax"` may change the
   decoy spread substantially, and σ is Eq. 1's denominator.
3. **Genuine method difference.** REF2015 and AWSEM are different force fields, and the
   paper's premise is that atomistic resolution reveals what coarse-graining cannot. Low
   correlation is not by itself proof of a bug.

> **All three were tested. See the next section** — candidates 1 and 2 are ruled out, and
> the disagreement is localised to AWSEM's water-mediated contacts (candidate 3).

---

# Resolution: the disagreement is concentrated in AWSEM's water-mediated contacts

## Experiment: relaxation protocol (candidate 2)

Two matched 500-decoy runs on the frustratometeR-prepared `1ubq_A.pdb`, `w = 0`,
`--min-seq-sep 2`, differing only in `--protocol`:

| | ρ vs frustratometeR | sd | median | GLY median |
|---|---|---|---|---|
| `min` | +0.141 | 0.94 | +0.04 | −0.31 |
| `relax` | +0.146 | 0.93 | +0.01 | −0.44 |

The two runs correlate with **each other** at ρ = +0.951 (n = 457 defined contacts).
Side-chain relaxation of the decoys does not materially change the index.

**Candidate 2 is ruled out.** With the contact-energy definition (the `w` sweep) and decoy
locality (pair-local decoys, ρ +0.265 vs +0.273) already eliminated, no mechanical choice
in our implementation accounts for the disagreement.

## Where the disagreement actually lives

Splitting the shared contacts by AWSEM's own well type is decisive, and reproduces
independently in both runs:

| AWSEM well type | n | ρ (`min`) | ρ (`relax`) | p |
|---|---|---|---|---|
| direct — short + long | 186 | **+0.309** | **+0.312** | < 0.0001 |
| water-mediated | 158 | −0.030 | −0.014 | 0.71 / 0.86 |
| all | 344 | +0.141 | +0.146 | ≈ 0.008 |

The methods agree significantly on contacts where **both model direct atomic interaction**,
and not at all — statistically indistinguishable from zero — on contacts where AWSEM
invokes its water-mediated well. Water-mediated contacts are 46 % of the shared set, so
averaging them in is what drags the overall ρ from 0.31 down to 0.14.

Restricted to direct contacts, our index also separates *their* classes in the right
direction: their `highly` contacts sit at our +0.16, their `minimally` at our +0.86
(separation +0.70; +0.67 for `relax`).

### Why this is the expected result, not a defect

AWSEM's water-mediated term is an **explicit** desolvation well: for pairs in the
~6.5–9.5 Å shell its depth is modulated by local residue density, standing in for a
bridging water molecule. REF2015 has no explicit water at all — desolvation enters
implicitly through `fa_sol`. On those pairs the two functions are not approximating the
same quantity, so there is no reason for their frustration indices to track.

This is the paper's own premise stated quantitatively: atomistic resolution and
coarse-grained resolution differ *specifically* where the coarse-grained model substitutes
an effective term for atoms it does not represent.

### Honest statement of what is and is not established

Agreement on the direct-contact subset (ρ ≈ 0.31, p < 10⁻⁴, consistent across two
independent decoy protocols) is evidence the implementation is not fundamentally wrong.
It is **not** proof of correctness: ρ = 0.31 leaves most of the variance unexplained, and
low correlation is also what an undiscovered bug looks like. The only remaining test that
can distinguish "correct but different" from "subtly wrong" is **reproducing a specific
published figure from Chen et al. (2020)** against our own output. Until that is done,
this section records a consistent hypothesis, not a validated one.

## Side finding: σ = 0 contacts at the cutoff edge

17 of 475 contacts have `decoy_std == 0` and `native_energy == 0` — no REF2015 term fires
for the pair in the native structure or in any of 500 decoys. All lie at Cα–Cα 8.5–9.9 Å
(defined pairs: median 7.5 Å), i.e. the outer edge of the 10 Å cutoff, where a Cα-based
contact criterion admits pairs with no atomistic interaction. Eq. 1 is 0/0 there; they are
emitted with `frustration_class = "undefined"` rather than a fabricated number.

None of the 17 appear in frustratometeR's contact set, so they never entered any
comparison above. Worth noting as an intrinsic mismatch between a Cα contact definition
and an all-atom energy — not a bug, but a reason the raw contact counts differ (475 vs 389).
