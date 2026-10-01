"""Small, explicit rules shared by residue displays."""
import numpy as np
import pandas as pd
import numbering


def complete_actin_positions(pos):
    # Preserve existing rows (including insertions); add missing reference positions.
    cols = [numbering.to_canon(p) for p in range(1, numbering.ACTIN_LENGTH + 1)]
    missing = sorted(set(cols) - set(pos["canon"]))
    return pd.concat([pos, pd.DataFrame({"canon": missing})], ignore_index=True)


def surface_masks(positions, rsa):
    values = pd.to_numeric(pd.Series([dict(rsa or {}).get(p) for p in positions]),
                           errors="coerce").to_numpy(dtype=float)
    known = np.isfinite(values) & (values >= 0)
    return known & (values >= .2), known & (values < .2)


def pdb_residue_frequencies(observations):
    valid = observations.dropna(subset=["pdb_id", "residue_name"])
    total = valid["pdb_id"].nunique()
    counts = valid.groupby("residue_name")["pdb_id"].nunique().rename("nb").to_frame()
    counts["pct"] = (100 * counts["nb"] / total).round(1) if total else 0.0
    return counts.reset_index().sort_values("nb", ascending=False), total


def select_interaction_chains(residues, interaction_chains):
    """Match the pair (interaction ID, chain), never independent membership lists."""
    expected = residues["interaction_id"].map(interaction_chains)
    return residues[residues["chain"].eq(expected)].copy()
