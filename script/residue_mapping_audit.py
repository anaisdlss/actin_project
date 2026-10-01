"""Audit of actin residue numbering against UniProt P60709.

For every actin chain, find the constant offset between the PDB-chain sequence
number (table 3, residue_number_sequence) and the P60709 numbering that best
explains the recorded residue letters. Writes a per-chain table and a
per-residue table with an explicit P60709 position. Nothing existing is
overwritten; chains without a convincing offset are flagged, not guessed.

Usage: python script/residue_mapping_audit.py   (from the repository root)
"""
from pathlib import Path
import pandas as pd

DETAILS = Path("data/filtered/details")
REF_FASTA = Path("data/P60709_ref.fasta")
OUT_CHAINS = Path("data/filtered/actin_numbering_audit_chains.csv")
OUT_RESIDUES = Path("data/filtered/actin_residue_mapping_p60709.csv")
OFFSETS = range(-15, 16)
MIN_MATCH = 0.85   # isoforms/species differ at a few positions
MIN_RESIDUES = 8


def load_ref():
    return "".join(l.strip() for l in REF_FASTA.read_text().splitlines()
                   if not l.startswith(">"))


def best_offset(seqpos, letters, ref):
    scores = []
    for k in OFFSETS:
        hits = sum(1 for p, a in zip(seqpos, letters) if 1 <= p + k <= len(ref) and ref[p + k - 1] == a)
        scores.append((hits / len(seqpos), k))
    return max(scores)


def main():
    ref = load_ref()
    res = pd.read_csv(DETAILS / "3.interface_residues.csv")
    prot = pd.read_csv(DETAILS / "2.proteins.csv").rename(columns={"chain_id": "chain"})
    actin = prot[prot["protein_name"].str.match(r"Actin, (alpha|cytoplasmic|gamma|beta)", na=False)]
    res = res.merge(actin[["interaction_id", "chain", "protein_name"]].drop_duplicates(),
                    on=["interaction_id", "chain"])
    res = res.dropna(subset=["residue_number_sequence"]).copy()
    res["residue_number_sequence"] = res["residue_number_sequence"].astype(int)

    rows, offsets = [], {}
    for (iid, chain), g in res.groupby(["interaction_id", "chain"]):
        frac, k = best_offset(g["residue_number_sequence"].tolist(), g["residue_name"].tolist(), ref)
        ok = frac >= MIN_MATCH and len(g) >= MIN_RESIDUES
        offsets[(iid, chain)] = k if ok else None
        rows.append(dict(interaction_id=iid, chain=chain, protein_name=g["protein_name"].iloc[0],
                         n_residues=len(g), offset_to_p60709=k, letter_match=round(frac, 3),
                         status="mapped" if ok else "unresolved"))
    chains = pd.DataFrame(rows)
    chains.to_csv(OUT_CHAINS, index=False)

    res["offset"] = [offsets[(i, c)] for i, c in zip(res["interaction_id"], res["chain"])]
    res["p60709_position"] = res["residue_number_sequence"] + res["offset"]
    res.loc[res["offset"].isna(), "p60709_position"] = pd.NA
    res["p60709_letter"] = [ref[int(p) - 1] if pd.notna(p) and 1 <= p <= len(ref) else None
                            for p in res["p60709_position"]]
    keep = ["interaction_id", "chain", "residue_number_structure", "residue_number_sequence",
            "residue_name", "residue_number_canon_mafft", "p60709_position", "p60709_letter"]
    res[keep].to_csv(OUT_RESIDUES, index=False)

    n = len(chains)
    print(f"{n} actin chains: {(chains.status == 'mapped').sum()} mapped, "
          f"{(chains.status == 'unresolved').sum()} unresolved")
    print("offset distribution (mapped):")
    print(chains[chains.status == "mapped"]["offset_to_p60709"].value_counts().head(8).to_string())


if __name__ == "__main__":
    main()
