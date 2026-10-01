"""Bidirectional actin contacts from the residue-pair source table."""
from pathlib import Path
import pandas as pd
import streamlit as st
import numbering

SOURCES = (Path("data/filtered/details/4.inter-residue_contacts.csv"),
           Path("data/filtered/proteins_per_pdb.csv"))


def orient_actin_contacts(contacts, proteins):
    annotations = proteins.drop_duplicates("chain").set_index("chain")
    flags = annotations["is_actin"].astype(str).str.lower()
    actin = set(flags[flags.eq("true")].index)
    nonactin = set(flags[flags.eq("false")].index)
    parts = []
    for side, other in (("A", "B"), ("B", "A")):
        rows = contacts[contacts[f"chain_{side}_id"].isin(actin)].copy()
        result = pd.DataFrame(index=rows.index)
        result["interaction_id"] = rows["interaction_id"]
        result["Actin chain"] = rows[f"chain_{side}_id"]
        result["PDB"] = result["Actin chain"].str.split("_").str[0]
        result["Alignment column"] = pd.to_numeric(rows[f"residue_{side}_canon_mafft"], errors="coerce")
        result["P60709 position"] = numbering.uniprot_series(result["Alignment column"])
        result["Actin PDB residue"] = rows[f"residue_{side}_structure"]
        result["Observed actin aa"] = rows[f"residue_{side}_name"]
        partner = rows[f"chain_{other}_id"]
        result["Partner chain"] = partner
        result["Partner"] = partner.map(annotations["protein"]).fillna("Unannotated protein")
        result["Interaction type"] = partner.map(lambda p: "Actin–actin" if p in actin else
                                                  "Actin–ABP" if p in nonactin else "Unclassified")
        result["Partner PDB residue"] = rows[f"residue_{other}_structure"]
        result["Partner aa"] = rows[f"residue_{other}_name"]
        result["Actin buried ASA (%)"] = pd.to_numeric(rows[f"asa_pct_{side}"], errors="coerce")
        result["Partner buried ASA (%)"] = pd.to_numeric(rows[f"asa_pct_{other}"], errors="coerce")
        result["Pair contact area (Å²)"] = pd.to_numeric(rows["contact_area"], errors="coerce")
        result["Source contact type"] = rows["contact_type"].fillna("Unspecified")
        parts.append(result)
    return pd.concat(parts, ignore_index=True).drop_duplicates()


@st.cache_data(show_spinner=False)
def load_contacts(mtimes):
    contacts = pd.read_csv(SOURCES[0], low_memory=False,
                           dtype={"residue_A_structure": str, "residue_B_structure": str})
    return orient_actin_contacts(contacts, pd.read_csv(SOURCES[1]))


def render_contacts(canon):
    with st.expander("All recorded residue contacts — actin and ABP partners"):
        if not all(p.exists() for p in SOURCES):
            st.info("Residue-contact source tables unavailable in this dataset.")
            return
        data = load_contacts(tuple(p.stat().st_mtime_ns for p in SOURCES))
        rows = data[data["Alignment column"].eq(canon)].copy()
        st.caption("Both sides of actin–actin contacts are included. ASA is the buried "
                   "percentage for each residue in its interface, not a percentage assigned "
                   "to this pair and not solvent accessibility (RSA). Partner residue numbers "
                   "are PDB coordinates. Blank source measurements remain unavailable.")
        if rows.empty:
            st.info("No interaction recorded for this position in the residue-pair table. "
                    "This alone does not establish whether the residue is buried or accessible.")
            return
        summary = rows.groupby("Interaction type").agg(
            PDBs=("PDB", "nunique"), Interactions=("interaction_id", "nunique"),
            Contact_records=("Partner chain", "size")).reset_index()
        st.dataframe(summary, hide_index=True, width="stretch")
        st.dataframe(rows.drop(columns=["Alignment column"]), hide_index=True, width="stretch")
        st.download_button("Download residue contacts (CSV)", rows.to_csv(index=False).encode(),
                           file_name=f"actin_P60709_{numbering.label(canon)}_contacts.csv",
                           mime="text/csv", key=f"all_contacts_{canon}")
