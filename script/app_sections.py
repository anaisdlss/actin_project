from display_helpers import plotly_chart
"""Scientific navigation and descriptions, shared by local and public builds.

The same scientific sections are used by both native page navigations.
"""
import ast
from pathlib import Path
import pandas as pd
import streamlit as st
from st_io import read_csv

SECTIONS = (
    ("documentation", "Documentation", "Methods, data sources and data management."),
    ("summary-tables", "Summary tables", "Source tables and exploration of individual PDB structures."),
    ("residue-level", "Actin use at a residue level", "Interaction types and the existing residue explorer."),
    ("actin-actin-interfaces", "Actin-actin interfaces", "Binding sites involved in actin–actin contacts, including mixed sites."),
    ("abp-actin-interfaces", "ABP-actin interfaces", "ABP footprints, interface clusters and structures."),
    ("comparative-binding-sites", "Comparative analyses of binding sites", "Compare binding sites, interaction networks and pairs of ABPs."),
    ("actin-conservation", "Actin mutational sensitivity", "ProteoCast sensitivity of actin residues contacted by a selected ABP."),
    ("human-actin-variants", "Human actin variants", ""),
    ("abp-conservation", "ABP mutational sensitivity", "ABP ProteoCast results and sequence alignments."),
    ("interface-properties", "Physico-chemical properties of the interface", "Contact chemistry and structural comparisons for a selected binding site."),
    ("homolog-search", "Homolog search", "Explore existing FoldDisco results for a selected ABP."),
)


def select_cluster_types(table, kind):
    """Mixed sites belong to BOTH views; never treat them as pure homo/hetero."""
    if kind not in {"all", "homo", "hetero"}:
        raise ValueError(f"Unknown interaction view: {kind}")
    if kind == "all":
        return table.copy()
    types = table["interaction_type"].astype(str).str.strip().str.lower()
    return table[types.isin({kind, "mixed", "mixte"})].copy()


CLUSTER_COLUMN_LABELS = {
    "patch": "Binding site",
    "interaction_type": "Type",
    "n_noeuds": "Nodes",
    "n_aretes": "Edges",
    "n_interactions": "Interactions",
    "s2_partners": "Partner proteins",
    "n_s2_partners": "Partners (n)",
    "s2_seq_clusters": "Partner sequence clusters",
    "n_s2_seq_clusters": "Partner sequence clusters (n)",
    "homo_pdb_ids": "PDBs with homo contacts",
    "n_homo_pdbs": "PDBs with homo contacts (n)",
}


def _as_text(value):
    if isinstance(value, str):
        text = value.strip()
        if text.startswith("[") and text.endswith("]"):
            try:
                value = ast.literal_eval(text)
            except (ValueError, SyntaxError):
                return text
        else:
            return text
    if isinstance(value, (list, tuple, set)):
        return ", ".join(str(v) for v in value)
    return value


def _readable_cluster_table(table):
    """Friendly headers, joined lists, and long PDB lists left to the detail views."""
    view = table.drop(columns=["ids_interactions", "homo_pdb_ids"], errors="ignore").copy()
    for col in view.columns:
        if view[col].dtype == object:
            view[col] = view[col].map(_as_text)
    return view.rename(columns=CLUSTER_COLUMN_LABELS)


def render_cluster_table(kind, show_selector=True):
    path = Path("data/filtered/patches_infos_s1_binding_site.csv")
    if not path.exists():
        st.info("Binding-site tables will be available after data generation.")
        return
    table = select_cluster_types(read_csv(str(path)), kind)
    st.markdown("**Binding-site summary**")
    st.caption("Homo = actin–actin only; hetero = actin–ABP only; mixed = both. "
               "A mixed site is included in both interface sections. Counts are those of the existing dataset.")
    if kind in {"homo", "hetero"}:
        st.dataframe(_readable_cluster_table(table), hide_index=True, width="stretch")
    else:
        with st.expander("Show binding-site table"):
            st.dataframe(_readable_cluster_table(table), hide_index=True, width="stretch")
    if show_selector and kind in {"homo", "hetero"} and not table.empty:
        cluster = st.selectbox("Binding site to explore", sorted(table["patch"].astype(str), key=str.casefold), key=f"{kind}_site_link")
        if st.button("Explore selected cluster", key=f"{kind}_site_open"):
            from app_navigation import request_page
            request_page("comparative-binding-sites", "Binding-site clusters", sel_s1=cluster)
            st.rerun()
    if kind == "homo" and not table.empty and "n_homo_pdbs" in table:
        counts = table.set_index("patch")["n_homo_pdbs"].sort_values(ascending=False)
        import plotly.graph_objects as go
        fig = go.Figure(go.Bar(x=counts.index, y=counts.values, marker_color="#0072B2",
                              hovertemplate="Site %{x}<br>%{y} PDBs with homo contacts<extra></extra>"))
        fig.update_layout(xaxis=dict(title="Binding-site cluster", categoryorder="array", categoryarray=counts.index.tolist()),
                          yaxis_title="PDBs with homo contacts", margin=dict(t=12, b=70), height=360)
        plotly_chart(fig, key="homo_pdb_counts", use_container_width=True)
        st.caption("PDB counts measure representation in this dataset, not physiological prevalence. "
                   "Major/minor interface labels have not been assigned. "
                   "Select a binding site below to inspect its contacts, heatmap and 3D structure.")


TABLE_DESCRIPTIONS = {
    "filtered_summary.csv": "PPI3D search results retained after structural filtering. Includes PDB identifiers, protein partners, interface type, buried surface area and contact counts.",
    "filtered_pdb_entry.csv": "Interaction pairs in retained PDB assemblies. Interactor 1/2 identify the partners; interface type distinguishes homo/hetero contacts; counts describe the interface.",
    "filtered_all_data.csv": "Merged structural and clustering data. subunit_1/2 are chain identifiers; s1/s2_actine identify actin; sequence clusters group sequences; binding-site clusters group sites; cluster_data groups interfaces at the stated threshold.",
    "1.interactions.csv": "One record per PPI3D interaction: PDB and chain identifiers, structure metadata, interface area and contact counts. interaction_id joins the detailed tables.",
    "2.proteins.csv": "Protein chains participating in each interaction, with their source name, organism and sequence length. A protein can appear in multiple interactions.",
    "3.interface_residues.csv": "Interface residues with PDB structure numbering, sequence numbering, residue identity and buried ASA. residue_number_canon_mafft is the current MAFFT alignment coordinate; its equivalence to UniProt numbering requires validation.",
    "4.inter-residue_contacts.csv": "Residue-to-residue contacts for chains A/B: PDB and sequence positions, residue identities, contact area and contact type. Canonical MAFFT columns are alignment coordinates; asa_pct_A/B are per-residue buried ASA percentages, not a per-pair decomposition.",
    "5.ligands.csv": "Ligand records associated with the retained interactions, as returned by PPI3D. interaction_id links each record to the interaction table.",
    "6.meta_alignement.csv": "PPI3D query/template alignment quality: scores, identity, similarity, gaps and interface coverage. These are upstream alignment metadata, not the locally generated MAFFT mapping.",
    "7.alignment_sequences.csv": "PPI3D query/template aligned sequences and interface masks. This is not the local MAFFT numbering table. For PDB/sequence/MAFFT correspondences at interface residues, use table 3; for contact pairs, use table 4.",
    "8.structures.csv": "Structure and assembly files associated with each interaction. File paths point to local coordinates or visualization scripts; some files are omitted from the public slim build.",
}

COLUMN_DESCRIPTIONS = {
    "interaction_id": "PPI3D interaction identifier; join key between detail tables.",
    "pdb_id": "PDB structure identifier.",
    "chain": "Chain containing the residue.",
    "residue_number_structure": "Residue identifier in the PDB structure.",
    "residue_number_sequence": "Position in the chain sequence.",
    "residue_number_canon_mafft": "Position in the current MAFFT alignment, not yet validated as a UniProt residue number.",
    "residue_name": "Amino-acid identity recorded for this residue.",
    "buried_ASA_Å²": "Solvent-accessible surface buried at the interface, in Å².",
    "buried_ASA_percent": "Percentage of the residue's accessible surface buried at this interface.",
    "contact_area": "Area attributed to this contact in the source table.",
    "asa_pct_A": "Buried ASA percentage of residue A (shared across its contact rows).",
    "asa_pct_B": "Buried ASA percentage of residue B (shared across its contact rows).",
    "query_sequence": "Aligned query sequence reported by PPI3D.",
    "template_sequence": "Aligned PDB chain sequence reported by PPI3D.",
    "interface_positions": "Mask identifying interface positions in the source alignment.",
}


def describe_table(name, columns):
    st.caption(f"Source file: {name}")
    st.info(TABLE_DESCRIPTIONS.get(name, "Source table from the current dataset."))
    known = [(c, COLUMN_DESCRIPTIONS[c]) for c in columns if c in COLUMN_DESCRIPTIONS]
    if known:
        with st.expander("Column guide"):
            st.dataframe(pd.DataFrame(known, columns=["Column", "Meaning"]),
                         hide_index=True, width="stretch")


TABLE_LABELS = {
    "1.interactions.csv": "Interactions",
    "2.proteins.csv": "Protein chains",
    "3.interface_residues.csv": "Interface residues and numbering",
    "4.inter-residue_contacts.csv": "Residue-to-residue contacts",
    "5.ligands.csv": "Ligands",
    "6.meta_alignement.csv": "Query alignment quality",
    "7.alignment_sequences.csv": "Query and template sequences",
    "8.structures.csv": "Structure files",
    "filtered_summary.csv": "Filtered search results",
    "filtered_pdb_entry.csv": "Interactions in retained structures",
    "filtered_all_data.csv": "Structural and clustering data",
}


def table_label(name):
    return TABLE_LABELS.get(name, name)
