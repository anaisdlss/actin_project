"""Native page navigation and persistent scientific choices."""
from pathlib import Path
import streamlit as st
from app_sections import SECTIONS


def request_page(section, view=None, **choices):
    """Queue navigation from an event callback without rendering another page."""
    if section not in {item[0] for item in SECTIONS}:
        raise ValueError(f"Unknown scientific section: {section}")
    st.session_state.update(choices)
    if view is not None:
        st.session_state[f"page_view_{section}"] = view
    st.session_state["_requested_page"] = section


def page_view(section, options):
    key = f"page_view_{section}"
    if st.session_state.get(key) not in options:
        st.session_state[key] = options[0]
    return st.radio("View", options, horizontal=True, key=key)


def run_navigation():
    # Assigning these values before widgets are created detaches them from the
    # previous page, preventing Streamlit's inactive-widget cleanup. Buttons,
    # Plotly selections and calculation triggers are deliberately not retained.
    for key in list(st.session_state):
        if key in PERSISTENT_SELECTIONS or key.startswith(PERSISTENT_PREFIXES):
            st.session_state[key] = st.session_state[key]
    pages = [st.Page(f"views/{key}.py", title=title, url_path=key,
                     default=(key == "documentation")) for key, title, _ in SECTIONS]
    selected = st.navigation(pages, expanded=True)
    with st.sidebar:
        st.caption("Actin–ABP · PPI3D")
        st.caption("Public dataset" if Path("data/.slim_deploy").exists() else "Local research project")
    requested = st.session_state.pop("_requested_page", None)
    if requested:
        st.switch_page(f"views/{requested}.py")
    selected.run()


# Explicitly keyed scientific controls (generated from the existing views).
PERSISTENT_SELECTIONS = {
    'actin_conservation_abp',
    'actin_ov_selbox',
    'actin_surface_mode',
    'all_site_link',
    'chemistry_s1',
    'cons_asa',
    'cons_cluster_abp',
    'cons_cluster_site',
    'conservation_abp',
    'explo_abp',
    'explo_c70',
    'explo_mode',
    'explo_pair_a',
    'explo_pair_b',
    'explo_pair_cl_a',
    'explo_pair_cl_b',
    'explo_res_selbox',
    'explo_seq_acc',
    'explo_seq_in',
    'filament_rsa_cutoff',
    'fp_aligned_abps',
    'fp_all_matrix',
    'fp_comparison',
    'fp_cutoff',
    'fp_reference',
    'global_graph_superclusters',
    'heatmap_norm_mode',
    'hetero_site_link',
    'hide_constant_columns',
    'homo_site_link',
    'homolog_abp',
    'hv_category',
    'hv_conflict_abp',
    'hv_cross_pair',
    'hv_disease',
    'hv_gene',
    'jac_thresh',
    'net_view',
    'numbering_chain',
    'numbering_pdb',
    'pdb_selector',
    'sel_abp_detail',
    'sel_c70',
    'sel_s1',
    'source_table',
    'structure_audit_filter',
    'surface_chem_3d',
    'surface_chem_cutoff',
    'surface_chem_site',
    'surface_chem_weight',
}
PERSISTENT_PREFIXES = ('abp3d_v2_', 'd_asa_thr_rel_', 'disco_cl_', 'fd_query_positions_', 'fd_combination_', 'page_view_', 'pc3d_whole_', 'pc_surface_', 'pc_view_', 'pc_whole_', 'restype_toggle_', 's1posdet_')
