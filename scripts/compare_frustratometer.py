"""Compare FrustX output against frustratometeR on the same structure.

frustratometeR computes configurational frustration with the AWSEM coarse-grained
energy function; FrustX uses atomistic Rosetta REF2015. The numbers are not expected to
match. What we are testing is whether they agree on the PATTERN -- which contacts are
minimally frustrated and which are frustrated -- and whether our sign convention and
classification thresholds transfer.

Usage:
    python scripts/compare_frustratometer.py <frustx_out_dir> <frustratometeR_file>

where <frustratometeR_file> is
    <ResultsDir>/<pdb>.done/FrustrationData/<pdb>.pdb_configurational

No scipy dependency: Spearman is computed as Pearson on ranks.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd


def spearman(a, b):
    """Rank correlation without scipy: Pearson on the ranks."""
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    return float(np.corrcoef(ra, rb)[0, 1])


def load_frustratometer(path):
    """Their per-contact table: whitespace-separated with a header line."""
    df = pd.read_csv(path, sep=r"\s+")
    df = df.rename(columns={"Res1": "i", "Res2": "j", "FrstIndex": "their_index",
                            "FrstState": "their_class", "AA1": "aa1", "AA2": "aa2"})
    # Normalise pair order so the join cannot miss on (i,j) vs (j,i).
    lo = np.minimum(df["i"], df["j"])
    hi = np.maximum(df["i"], df["j"])
    df["i"], df["j"] = lo, hi
    return df[["i", "j", "aa1", "aa2", "their_index", "their_class", "Welltype"]]


def load_frustx(path):
    df = pd.read_csv(path)
    df = df.rename(columns={"resnum_i": "i", "resnum_j": "j",
                            "frustration_index": "our_index",
                            "frustration_class": "our_class"})
    lo = np.minimum(df["i"], df["j"])
    hi = np.maximum(df["i"], df["j"])
    df["i"], df["j"] = lo, hi
    return df[["i", "j", "resname_i", "resname_j", "our_index", "our_class",
               "native_energy", "decoy_mean", "decoy_std"]]


def main(frustx_dir, frustra_file):
    ours = load_frustx(Path(frustx_dir) / "contacts.csv")
    theirs = load_frustratometer(frustra_file)

    print(f"FrustX contacts          : {len(ours)}")
    print(f"frustratometeR contacts  : {len(theirs)}")

    merged = ours.merge(theirs, on=["i", "j"], how="inner")
    only_ours = len(ours) - len(merged)
    only_theirs = len(theirs) - len(merged)
    print(f"shared contacts          : {len(merged)}"
          f"   (only FrustX: {only_ours}, only frustratometeR: {only_theirs})")

    if merged.empty:
        print("\nNo shared contacts -- residue numbering probably disagrees.")
        return 1

    m = merged.dropna(subset=["our_index", "their_index"])
    print(f"\n--- correlation on {len(m)} shared contacts ---")
    print(f"Spearman : {spearman(m['our_index'], m['their_index']):+.3f}")
    print(f"Pearson  : {np.corrcoef(m['our_index'], m['their_index'])[0,1]:+.3f}")

    print("\n--- index distributions ---")
    for name, col in (("FrustX", "our_index"), ("frustratometeR", "their_index")):
        v = m[col]
        print(f"{name:<16} min {v.min():+.2f}  median {v.median():+.2f}  "
              f"max {v.max():+.2f}  sd {v.std():.2f}")

    print("\n--- classification agreement (rows FrustX, cols frustratometeR) ---")
    print(pd.crosstab(m["our_class"], m["their_class"]).to_string())
    agree = (m["our_class"] == m["their_class"]).mean()
    print(f"\nexact agreement: {agree:.1%}")

    # Correlation by AWSEM well type: a systematic split would say our contact energy
    # behaves differently for buried vs water-mediated contacts.
    print("\n--- Spearman by AWSEM well type ---")
    for well, grp in m.groupby("Welltype"):
        if len(grp) > 5:
            print(f"  {well:<16} n={len(grp):<4} rho={spearman(grp['our_index'], grp['their_index']):+.3f}")

    # The glycine/terminus artifact: does frustratometeR show it too?
    print("\n--- contacts involving GLY (the suspected shuffling artifact) ---")
    gly = m[(m["resname_i"] == "GLY") | (m["resname_j"] == "GLY")]
    if not gly.empty:
        print(f"n={len(gly)}   FrustX median {gly['our_index'].median():+.2f}   "
              f"frustratometeR median {gly['their_index'].median():+.2f}")
        print(f"non-GLY  FrustX median {m[~m.index.isin(gly.index)]['our_index'].median():+.2f}   "
              f"frustratometeR median {m[~m.index.isin(gly.index)]['their_index'].median():+.2f}")
        print("\n  worst offenders by our index:")
        for _, r in gly.nlargest(5, "our_index").iterrows():
            print(f"    {r.resname_i}{r.i:<4}-{r.resname_j}{r.j:<4} "
                  f"ours {r.our_index:+.2f}  theirs {r.their_index:+.2f}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
