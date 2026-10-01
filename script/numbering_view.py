"""Inspect the coordinates used by the app without changing source tables."""
from pathlib import Path
import pandas as pd
import numbering


def chain_coordinates(residues, chain):
    columns = ["chain", "residue_number_structure", "residue_number_sequence",
               "residue_name", "residue_number_canon_mafft"]
    result = residues.loc[residues["chain"].eq(chain), columns].drop_duplicates().copy()
    result["P60709 position"] = numbering.uniprot_series(result["residue_number_canon_mafft"])
    result["Mapping"] = result["P60709 position"].map(
        lambda x: "No P60709 position" if pd.isna(x) else "Mapped by actin alignment")
    return result.sort_values(["residue_number_sequence", "residue_number_structure"])


def render_numbering_lookup(read_csv):
    import streamlit as st
    residues_path = Path("data/filtered/details/3.interface_residues.csv")
    proteins_path = Path("data/filtered/proteins_per_pdb.csv")
    if not (residues_path.exists() and proteins_path.exists()):
        return
    with st.expander("Residue numbering: PDB → alignment → P60709"):
        st.caption("Correspondence for actin interface residues present in the source tables. "
                   "This is not the complete sequence of each PDB chain. P60709 positions "
                   "use the same actin-alignment conversion as the app; insertions and "
                   "unmapped positions remain empty. Chain identifiers are case-sensitive.")
        proteins = read_csv(str(proteins_path))
        residues = read_csv(str(residues_path))
        actin = proteins.loc[proteins["is_actin"].astype(str).str.lower().eq("true"), "chain"]
        chains = sorted(set(actin.dropna()) & set(residues["chain"].dropna()))
        if not chains:
            st.info("No actin interface residues available.")
            return
        pdbs = sorted({c.split("_", 1)[0] for c in chains})
        pdb = st.selectbox("PDB for residue numbering", pdbs, key="numbering_pdb")
        chain = st.selectbox("Actin chain", [c for c in chains if c.split("_", 1)[0] == pdb],
                             key="numbering_chain")
        table = chain_coordinates(residues, chain)
        st.dataframe(table.rename(columns={
            "chain": "PDB chain", "residue_number_structure": "PDB residue number",
            "residue_number_sequence": "Chain sequence position", "residue_name": "Observed residue",
            "residue_number_canon_mafft": "Actin alignment column"}), hide_index=True, width="stretch")
        st.download_button("Download this chain's residue correspondence (CSV)",
                           table.to_csv(index=False).encode("utf-8"),
                           file_name=f"{chain}_residue_numbering.csv", mime="text/csv",
                           key="numbering_download")
