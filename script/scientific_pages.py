"""One scientific question per page, with only the chosen analysis loaded."""
from pathlib import Path
import streamlit as st
from app_sections import SECTIONS, render_cluster_table
from app_navigation import page_view


def render(section):
    title, description = next((title, desc) for key, title, desc in SECTIONS if key == section)
    st.header(title, anchor=section)
    if description:
        st.caption(description)
    # Import the renderers after navigation. Importing does not read datasets.
    import page_views as views

    if section == "documentation":
        st.markdown("Explore actin residues, interfaces and their binding partners using the sections in the sidebar. "
                    "Selected residues, proteins and analysis settings are retained when you return.")
        with st.expander("Data management — reload, download and calculations"):
            st.caption("Re-read the files already installed. This does not download new data or start calculations.")
            if st.button("Clear cache and reload", key="reload_local_cache"):
                st.cache_data.clear()
                st.cache_resource.clear()
                for key in ("viewer_key", "viewer_html"):
                    st.session_state.pop(key, None)
                st.rerun()
            if views.DEPLOY_MODE:
                st.info("This public app uses a prepared dataset. Data updates and calculations are managed in the full project.")
                from local_rebuild_ui import render as render_local_rebuild
                render_local_rebuild()
            else:
                views.render_data_tools()
        guide = Path("GUIDE.md")
        if guide.exists():
            content = guide.read_text()
            if content.startswith("# "):
                content = content.partition("\n")[2]
            st.markdown(content)
    elif section == "summary-tables":
        view = page_view(section, ["Structures", "Source tables", "Residue numbering", "Dataset checks"])
        if view == "Structures": views.render_structures()
        elif view == "Source tables": views.render_source_tables()
        elif view == "Residue numbering": views.render_numbering_lookup(views.read_csv)
        else:
            views.render_structure_annotations()
            views.render_threshold_comparison(views.read_csv)
    elif section == "residue-level":
        view = page_view(section, ["Residue explorer", "Binding-site heatmap"])
        if view == "Residue explorer": views.render_residues()
        else:
            views.render_s1_heatmap()
            render_cluster_table("all")
    elif section == "actin-actin-interfaces":
        view = page_view(section, ["Binding sites", "Compare footprints", "Structural evidence"])
        if view == "Binding sites": views.render_s1_cluster(kind="homo")
        elif view == "Compare footprints": views.render_footprint_comparison()
        else: views.render_interface_evidence()
    elif section == "abp-actin-interfaces":
        view = page_view(section, ["By protein", "Overview and heatmap", "Binding sites"])
        if view == "By protein": views.render_abp_details()
        elif view == "Overview and heatmap": views.render_abp_overview()
        else: views.render_s1_cluster(kind="hetero")
    elif section == "comparative-binding-sites":
        view = page_view(section, ["Binding-site clusters", "Interaction clusters", "ABP networks", "ABP pairs and sequences"])
        if view == "Binding-site clusters": views.render_s1_cluster()
        elif view == "Interaction clusters": views.render_c70_cluster()
        elif view == "ABP networks": views.render_abp_networks()
        else: views.render_sequence_comparison()
    elif section == "actin-conservation":
        view = page_view(section, ["Overview", "By ABP footprint", "Solvent accessibility"])
        if view == "Overview": views.render_conservation()
        elif view == "By ABP footprint": views.render_actin_abp_conservation()
        else: views.render_filament_accessibility()
    elif section == "human-actin-variants":
        views.render_variants()
    elif section == "abp-conservation":
        view = page_view(section, ["Profiles", "3D structure", "Sequence alignments"])
        if view != "Sequence alignments": views.render_abp_proteocast(view)
        else: views.render_alignments()
    elif section == "interface-properties":
        view = page_view(section, ["Actin surface chemistry", "By binding site"])
        if view == "Actin surface chemistry": views.render_interface_properties()
        else: views.render_cluster_chemistry()
    elif section == "homolog-search":
        views.render_homologs()
