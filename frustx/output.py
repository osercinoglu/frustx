"""Turning a FrustrationResult into things you can look at.

Three views, in increasing order of lossiness:

  * `contact_table`  -- one row per contact. The primary result; nothing is aggregated.
  * `residue_table`  -- one row per residue, aggregated over its contacts. Convenient
                        for plotting a profile along the sequence, but note the paper
                        defines frustration per CONTACT, not per residue: this
                        aggregation is ours, not theirs.
  * `write_bfactor_pdb` -- per-residue values painted into the B-factor column, for
                        colouring in PyMOL/ChimeraX.
"""

import numpy as np
import pandas as pd

from frustx.frustration import classify


def contact_table(result):
    """One row per contact. This is the primary output."""
    rows = []
    for i, j in result.contacts:
        ri, rj = result.residues[i], result.residues[j]
        rows.append(
            {
                "chain_i": ri.chain, "resnum_i": ri.resseq, "resname_i": ri.resname,
                "chain_j": rj.chain, "resnum_j": rj.resseq, "resname_j": rj.resname,
                "native_energy": result.native_energy[i, j],
                "decoy_mean": result.decoy_mean[i, j],
                "decoy_std": result.decoy_std[i, j],
                "frustration_index": result.index[i, j],
            }
        )
    df = pd.DataFrame(rows)
    if not df.empty:
        df["frustration_class"] = classify(df["frustration_index"].to_numpy())
    return df


def residue_table(result):
    """One row per residue, aggregating the contacts it participates in.

    `mean_frustration` is the mean index over that residue's contacts. The counts
    mirror what frustratometeR reports, and are often more informative than the mean:
    a residue with many minimally frustrated contacts and a few highly frustrated ones
    is a different thing from a uniformly neutral residue, and averaging hides that.
    """
    n = len(result.residues)
    per_residue = [[] for _ in range(n)]
    for i, j in result.contacts:
        value = result.index[i, j]
        per_residue[i].append(value)
        per_residue[j].append(value)

    rows = []
    for idx, res in enumerate(result.residues):
        values = np.array(per_residue[idx], dtype=float)
        finite = values[np.isfinite(values)]
        labels = classify(finite) if finite.size else np.array([], dtype=object)
        rows.append(
            {
                "chain": res.chain,
                "resnum": res.resseq,
                "resname": res.resname,
                "n_contacts": len(values),
                "mean_frustration": finite.mean() if finite.size else np.nan,
                "n_minimally_frustrated": int((labels == "minimally").sum()),
                "n_neutral": int((labels == "neutral").sum()),
                "n_highly_frustrated": int((labels == "highly").sum()),
            }
        )
    return pd.DataFrame(rows)


def write_bfactor_pdb(pose, result, path, column="mean_frustration"):
    """Write the pose with per-residue frustration in the B-factor column.

    Colour it in PyMOL with:  spectrum b, blue_white_red, all

    Residues with no finite value (no contacts) get 0.0 rather than NaN, since most
    viewers choke on NaN in the B-factor field.
    """
    table = residue_table(result).set_index(["chain", "resnum"])
    info = pose.pdb_info()
    if info is None:
        raise ValueError("pose has no PDB info; cannot write B-factors")

    for i in range(1, pose.total_residue() + 1):
        key = (info.chain(i), info.number(i))
        value = table[column].get(key, np.nan) if key in table.index else np.nan
        if not np.isfinite(value):
            value = 0.0
        # B-factors are per-atom, so paint every atom of the residue the same value.
        for atom in range(1, pose.residue(i).natoms() + 1):
            info.bfactor(i, atom, float(value))

    pose.dump_pdb(str(path))
    return path
