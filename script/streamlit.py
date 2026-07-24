# streamlit run script/streamlit.py
from abp_3d import (
    _s1_abp_3d_options, _build_abp_actin_3d, _build_all_abp_pdb, _build_all_abp_3d, _ABP_MULTI_COLORS)
from proteocast_abp import (
    _abp_actin_focus, _crop_pdb_residues, _render_abp_proteocast, _render_abp_actin_conservation)
from cluster_struct import (
    _render_cluster_struct)
from s1_heatmaps import (
    _build_s1_global_heatmap, _render_s1_global_plotly, _build_s1_patch_detail, _render_s1_patch_plotly, _s1_position_detail, _render_s1_position_detail, _S1_GLOBAL_FILES)
import pipeline_ui
from st_io import (_load_pdb_file, read_csv)
import warnings as _warnings
import logging as _logging
import re as _re_sd
from msa_analysis import (
    norm_chain_id, _cluster_protein_index,
    _msa_section_full, _msa_section_s2_clusters, _msa_section_s1_clusters,
    _MSA_ALN_DIR, _MSA_RIGOR_PDBS,
    _s1_get_ch_maps, _msa_contact_analysis,
    _msa_blast_celegans, _msa_s1_cluster_inline,
)
from data_analysis import (
    _build_abp_heatmap_data, _build_c70_jaccard_edges,
    _build_s1_superclusters, _S1_SUPER_FILES, _ABP_HM_FILES,
)
from residue_passport import build_passport, pp_mtimes
import residue_explorer
import folddisco_view
import proteocast_view
from network_viz import (
    _bip_mtimes, _load_bipartite_base, _build_bipartite_html,
    _load_res4, _build_bipartite_c70_html, _build_s1_3d_html,
    _build_tripartite_graph_html, _build_global_graph_html,
    _global_bfac_max, _BIP_CACHE_VERSION, _BIPARTITE_FILES,
    _AA_RESTYPE_HEX,
)
import matplotlib.pyplot as plt
import streamlit as st
import pandas as pd
import shutil
import ast
import csv
import os
import re
import subprocess
import sys
import numpy as np
from collections import Counter as _Counter
from pathlib import Path as _Path
import matplotlib
matplotlib.use("Agg")

# ── Faire taire les avertissements au runtime ───────────────────────────────
# FutureWarning pandas + dépréciations Streamlit (use_container_width,
# st.components.v1.html) qui polluaient la console à chaque run.
_warnings.filterwarnings("ignore", category=FutureWarning)
_warnings.filterwarnings("ignore", category=DeprecationWarning)
for _ln in ("streamlit", "streamlit.deprecation_util",
            "streamlit.elements.lib.policies", "streamlit.elements"):
    _logging.getLogger(_ln).setLevel(_logging.ERROR)


_GRAPH_COMPONENT = st.components.v1.declare_component(
    "actin_graph",
    path=str(_Path(__file__).parent / "graph_component"),
)

st.set_page_config(layout="wide", page_title="Actin-ABP analysis - PPI3D")

# Style global : chaque GROS titre de section (st.header = h2) est SURLIGNÉ
# (bandeau de fond + accent rouge à gauche) → on voit clairement le début de
# chaque nouvelle partie, sans avoir besoin de traits séparateurs.
st.markdown(
    "<style>"
    "h2{background:#eef0f2;border-left:6px solid #6b7280;"
    "padding:8px 16px!important;border-radius:6px;margin-top:8px}"
    # pas de surlignage dans la barre latérale (titre « Navigation »)
    "[data-testid='stSidebar'] h2{background:none;border-left:none;"
    "padding:0!important;border-radius:0}"
    "</style>",
    unsafe_allow_html=True,
)

# Mode rapide : `FAST=1 streamlit run script/streamlit.py`
# -> saute le rendu lourd (heatmap ABP, Détail par ABP, Filament PyMOL, MSA)
#    et ne garde que les RÉSEAUX (compétition + coopération).
_FAST = os.environ.get("FAST", "") == "1"

st.title("Actin-actin and actin-ABP interaction analysis - PPI3D")

# ── Documentation : tout le sens (méthodes, glossaire, ce que veut dire l'ASA,
#    d'où viennent les données) est dans GUIDE.md, affiché ici. ───────────────
st.header("Documentation", anchor="documentation")
_guide = _Path("GUIDE.md")
with st.expander("How to read this app — methods, glossary, data sources",
                 expanded=False):
    if _guide.exists():
        st.markdown(_guide.read_text(), unsafe_allow_html=False)
    else:
        pass

if _FAST:
    st.info("FAST mode: only the networks are loaded (FAST=1).")

with st.sidebar:
    st.markdown("## Contents")
    st.markdown("""
- [Documentation](#documentation)
- [Data download](#telechargement-des-donnees)
- [Filtered data](#donnees-filtrees-s1-actin)
- [Valid PDB structures](#structures-pdb-valides)
- [ABP footprint on actin](#empreinte-abp)
- [Interaction clusters](#clusters-d-interactions)
- [ABP](#abp)
- [Interactive explorer](#explorateur)
- [MSA — Interface proteins](#msa-proteines)
""")
    st.divider()
    # Remède aux « anciennes données » : vide tous les caches (@st.cache_data /
    # @st.cache_resource) et relit les fichiers depuis le disque.
    if st.button("Clear cache and reload", width="stretch",
                 help="Use after regenerating the data (pipeline) "
                      "to force a re-read and avoid showing stale values."):
        st.cache_data.clear()
        st.cache_resource.clear()
        for _k in ("viewer_key", "viewer_html"):
            st.session_state.pop(_k, None)
        st.rerun()

# ---------------------------------------------------------------------------
# Section téléchargement
# ---------------------------------------------------------------------------

# Mode « déployé » (Streamlit Cloud, données slim pré-embarquées) : on masque le
# pipeline « Run / update » et le calcul ProteoCast, impossibles sur Cloud.
DEPLOY_MODE = os.path.exists("data/.slim_deploy")

if not DEPLOY_MODE:
    pipeline_ui.render()

# ── ProteoCast (opt-in) — campagne de calcul pour TOUS les ABP ─────────────
# Placé juste sous le téléchargement principal : c'est une campagne séparée du
# Run/update (elle dure des heures ; 1 job à la fois ; reprenable). Les
# résultats par ABP s'affichent ensuite dans la section « Per-ABP detail ».
_pc_mt = os.path.getmtime("data/proteocast/abp_inputs/manifest.csv") \
    if os.path.exists("data/proteocast/abp_inputs/manifest.csv") else 0.0
_pc_status = proteocast_view.load_status(_pc_mt)
# Compteur LIVE (lu directement sur le disque à chaque run — pas le cache, qui
# resterait figé tant que le manifest ne change pas) : « remaining » = ABP sans
# résultat ET pas en échec définitif. Les ABP que ProteoCast ne SAIT PAS calculer
# (trop grosses : Myosin/β-myosin ~1900 aa ; fusions ; MSA trop pauvre) sont
# consignés dans _failed_slugs.txt → on les sort du « remaining » (sinon il ne
# tomberait jamais à 0) et on les annonce à part.
import glob as _glob
_pc_all = set(_pc_status["slug"].astype(str)) if _pc_status is not None else set()
_pc_done = {os.path.basename(os.path.dirname(_p))
            for _p in _glob.glob("data/proteocast/abp/*/4.query_ProteoCast.csv")}
_pc_failf = "data/proteocast/abp/_failed_slugs.txt"
_pc_permfail = ({l.strip() for l in open(_pc_failf) if l.strip()}
                if os.path.exists(_pc_failf) else set())
# « remaining » = tous les ABP sans résultat (lu LIVE sur le disque, pas le cache).
_pc_missing = _pc_all - _pc_done
_pc_miss = len(_pc_missing)
# Parmi eux, ceux que ProteoCast ne SAIT PAS calculer (déjà en échec définitif :
# trop grosses ~1900 aa, fusions, MSA trop pauvre) → signalés à part.
_pc_uncomputable = len(_pc_missing & _pc_permfail)
if not DEPLOY_MODE:
    st.markdown("**ABP ProteoCast — compute the mutational landscape for all ABPs**")
    _pc_cap = ("Opt-in, separate from `Run / update` — computing **all** ABPs "
               "can take **several hours** (one ABP at a time). "
               "Resumable: skips those already done.")
    if _pc_uncomputable:
        _pc_cap += (f"  \n_Of these, {_pc_uncomputable} can't be computed by ProteoCast "
                    "(protein too large / fusion / weak MSA) — expected, they will "
                    "keep failing._")
    _pc_label = (f"Compute all missing ProteoCast — {_pc_miss} remaining" if _pc_miss
                 else "Update ProteoCast (refresh + retry)")
    # bouton rouge (comme le téléchargement) mais un peu plus transparent
    st.markdown(
        "<style>[class*='st-key-pc_run_all_missing'] button{"
        "background:rgba(224,82,82,0.75)!important;border-color:rgba(224,82,82,0.4)!important;"
        "color:#fff!important}</style>",
        unsafe_allow_html=True,
    )
    _pc_go = st.button(_pc_label, key="pc_run_all_missing", type="primary")
else:
    _pc_go = False
if _pc_go:
    with st.status("Computing ProteoCast via proteocast.ijm.fr "
                   "(one ABP at a time)…", expanded=True) as _pcall:
        _pclog = st.empty()
        _pclog.code("Preparing the ABP manifest…")
        # 1) régénérer le manifest depuis abp_master (couvre les ABP actuels)
        subprocess.run([sys.executable, "-m", "script.proteocast_prep_manifest"])
        # 2) soumettre les manquants (1 à la fois, reprenable)
        #    bufsize=1 + readline : sortie ligne-par-ligne EN DIRECT (sinon le
        #    buffer de lecture anticipée retient tout par blocs → rien à l'écran).
        # --retry-failed : on re-tente aussi les ABP déjà marqués en échec (le plus
        # souvent des échecs transitoires : timeout serveur, MSA…) → sinon un clic
        # sur « X remaining » ne ferait rien si ces X sont tous en échec.
        _proc = subprocess.Popen(
            [sys.executable, "-u", "script/proteocast_submit_abp.py",
             "--retry-failed", "--jobs", "1"],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)

        # Tableau de bord VIVANT : on n'affiche pas le log brut qui s'accumule,
        # mais seulement l'étape ACTUELLE de chaque job (+ compteurs faits/échoués).
        import re as _re
        _board = {}            # titre ABP -> message d'étape en cours
        _done, _failed = [], []
        _total = [None]

        def _pretty(_s):
            _s = _s.strip().rstrip(".").strip()
            return (_s[:1].upper() + _s[1:] + "…") if _s else "…"

        def _render():
            _rows = [f"**{len(_done)} done** · **{len(_failed)} failed**"
                     + (f" · {_total[0]} to compute" if _total[0] else "")]
            if _board:
                _rows.append(f"\n**In progress ({len(_board)})**")
                _rows += [f"- **{_t}** — {_m}" for _t, _m in _board.items()]
            elif _done or _failed:
                _rows.append("\nWaiting for the next submissions…")
            _pclog.markdown("\n".join(_rows))

        _render()
        for _l in iter(_proc.stdout.readline, ""):
            _l = _l.rstrip()
            if not _l:
                continue
            m = _re.match(r"^(\d+) ABP à soumettre", _l)
            if m:
                _total[0] = int(m.group(1)); _render(); continue
            m = _re.match(r"^\[soumis .*?\]\s+(.+?)\s+\([^)]*\)\s+—\s+job", _l)
            if m:
                _board[m.group(1)] = "Submitted — waiting for first status…"
                _render(); continue
            m = _re.match(r"^\[\d+/\d+\] OK — (.+?)\s+->", _l)
            if m:
                _board.pop(m.group(1), None); _done.append(m.group(1))
                _render(); continue
            if _l.endswith("marqué en échec"):
                m = _re.match(r"^\s+(.+?)\s+:\s+.*—\s+marqué en échec$", _l)
                if m:
                    _board.pop(m.group(1), None); _failed.append(m.group(1))
                    _render()
                continue
            m = _re.match(r"^\s{4,}(.+?)\s+:\s+(.+)$", _l)
            if m and m.group(1) in _board:
                _board[m.group(1)] = _pretty(m.group(2)); _render(); continue
        _proc.wait()
        _pclog.markdown(
            f"**Finished** — {len(_done)} computed · {len(_failed)} failed"
            + (f" (out of {_total[0]})" if _total[0] else ""))
        _pcall.update(label="ProteoCast campaign finished.", state="complete")
    st.cache_data.clear()
    st.rerun()

st.divider()
st.header("Filtered data", anchor="donnees-filtrees-s1-actin")


TABLES = {
    "filtered_summary.csv": "data/filtered/filtered_summary.csv",
    "filtered_pdb_entry.csv": "data/filtered/filtered_pdb_entry.csv",
    "filtered_all_data.csv": "data/filtered/filtered_all_data.csv",
    "1.interactions.csv": "data/filtered/details/1.interactions.csv",
    "2.proteins.csv": "data/filtered/details/2.proteins.csv",
    "3.interface_residues.csv": "data/filtered/details/3.interface_residues.csv",
    "4.inter-residue_contacts.csv": "data/filtered/details/4.inter-residue_contacts.csv",
    "5.ligands.csv": "data/filtered/details/5.ligands.csv",
    "6.meta_alignement.csv": "data/filtered/details/6.meta_alignement.csv",
    "7.alignment_sequences.csv": "data/filtered/details/7.alignment_sequences.csv",
    "8.structures.csv": "data/filtered/details/8.structures.csv",
}

available_tables = {t: p for t, p in TABLES.items() if os.path.exists(p)}

if available_tables:
    selected = st.selectbox("Choose a table",
                            list(available_tables.keys()))
    df = read_csv(available_tables[selected])
    hide_constant = st.checkbox(
        "Hide columns without variation", value=False)
    if hide_constant:
        cols_to_show = [c for c in df.columns if df[c].nunique() > 1]
    else:
        cols_to_show = list(df.columns)
    st.dataframe(df[cols_to_show], width="stretch")

# Métrique PDB
METRICS_FILES = {
    "pdb": "data/filtered/filtered_pdb_entry.csv",
    "interactions": "data/filtered/details/1.interactions.csv",
    "proteins": "data/filtered/details/2.proteins.csv",
    "residues": "data/filtered/details/3.interface_residues.csv",
    "ligands": "data/filtered/details/5.ligands.csv",
}

if all(os.path.exists(p) for p in METRICS_FILES.values()):
    df_pdb = read_csv(METRICS_FILES["pdb"])
    df_int = read_csv(METRICS_FILES["interactions"])
    df_prot = read_csv(METRICS_FILES["proteins"])
    df_res = read_csv(METRICS_FILES["residues"])
    df_lig = read_csv(METRICS_FILES["ligands"])

    nb_homo = int((df_pdb["Interface type"] == "homo").sum()
                  ) if "Interface type" in df_pdb.columns else 0
    nb_hetero = int((df_pdb["Interface type"] == "hetero").sum(
    )) if "Interface type" in df_pdb.columns else 0

    col1, col2, col3 = st.columns(3)
    col1.metric("PDB structures",
                f"{df_pdb['pdb_id'].nunique():,}".replace(',', ' '))
    col2.metric("Protein-protein interactions",
                f"{len(df_int):,}".replace(',', ' '))
    col3.metric("Unique partner proteins (ABPs)",
                f"{df_prot[~df_prot['protein_name'].str.lower().str.contains('actin', na=False)]['protein_name'].nunique():,}".replace(',', ' '))

    # Transparency funnel: how many structures PPI3D returned vs how many are kept
    # after the "≥ 5 connected actin subunits" filter (makes the counts explicit,
    # and explains why the total can go up or down between PPI3D snapshots).
    _raw_summary = "data/raw/ppi3d_actin_summary.csv"
    if os.path.exists(_raw_summary):
        try:
            _df_raw = read_csv(_raw_summary)
            _raw_pdb = int(_df_raw["PDB ID"].astype(
                str).str.upper().str.strip().nunique())
            _kept_pdb = int(df_pdb["pdb_id"].nunique())
            _n_int = f"{len(df_int):,}".replace(",", " ")
            # Visible completeness check: a partial/interrupted details download
            # leaves fewer interactions than the summary — warn instead of failing
            # silently (the pipeline now auto-completes on the next run).
            _fsum = "data/filtered/filtered_summary.csv"
            if os.path.exists(_fsum):
                _n_expected = int(read_csv(_fsum)["interaction_id"].nunique())
                if len(df_int) < _n_expected:
                    _got = f"{len(df_int):,}".replace(",", " ")
                    _exp = f"{_n_expected:,}".replace(",", " ")
                    st.warning(
                        f"Interface details look incomplete: only {_got} of {_exp} "
                        "interactions were downloaded (interrupted run?). Click "
                        "**Run / update** — the pipeline will fetch the missing ones.")
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Section PDB valides — explorateur
# ---------------------------------------------------------------------------

pdb_filt_path = "data/filtered/filtered_pdb_entry.csv"
if os.path.exists(pdb_filt_path):
    st.divider()
    st.header("Valid PDB structures", anchor="structures-pdb-valides")

    df_entry = read_csv(pdb_filt_path)
    pdb_ids = sorted(df_entry["pdb_id"].str.upper().unique())

    summary_path = "data/filtered/filtered_summary.csv"
    if os.path.exists(summary_path):
        df_sum = read_csv(summary_path)
        unique = df_sum[["PDB ID", "Structure title"]
                        ].drop_duplicates("PDB ID")
        title_map = dict(
            zip(unique["PDB ID"].str.upper(), unique["Structure title"]))
    else:
        title_map = {}

    col_sel, _ = st.columns([1, 1])
    with col_sel:
        selected_pdb = st.selectbox(
            f"Choose a structure ({len(pdb_ids)} valid PDBs)",
            pdb_ids,
            format_func=lambda x: f"{x} — {title_map[x]}" if x in title_map else x,
            key="pdb_selector",
        )

    sub = df_entry[df_entry["pdb_id"].str.upper() == selected_pdb]

    # Reset de l'interaction sélectionnée quand on change de PDB
    if st.session_state.get("last_pdb") != selected_pdb:
        st.session_state["sel_inter"] = None
        st.session_state["sel_node"] = None
        st.session_state["last_pdb"] = selected_pdb

    if "sel_inter" not in st.session_state:
        st.session_state["sel_inter"] = None
    if "sel_node" not in st.session_state:
        st.session_state["sel_node"] = None

    int1_path = "data/filtered/details/1.interactions.csv"
    pp_path = "data/filtered/proteins_per_pdb.csv"
    struct_csv = "data/filtered/details/8.structures.csv"
    cont_csv = "data/filtered/details/4.inter-residue_contacts.csv"

    # Données communes aux deux colonnes
    df_int1_g = read_csv(int1_path) if os.path.exists(int1_path) else None
    sub_int1_g = df_int1_g[df_int1_g["pdb_id"].str.upper(
    ) == selected_pdb] if df_int1_g is not None else None

    raw_path = "data/raw/pdb_entry_results.csv"
    sub_raw_g = None
    if os.path.exists(raw_path):
        df_raw_g = read_csv(raw_path)
        sub_raw_g = df_raw_g[df_raw_g["pdb_id"].str.upper() == selected_pdb]

    col_net, col_3d = st.columns(2)

    # ── Colonne gauche : réseau interactif ────────────────────────────────
    with col_net:
        st.markdown(f"**Interaction network — {selected_pdb}**")
        if os.path.exists(pp_path) and sub_int1_g is not None:
            df_pp_g = read_csv(pp_path)
            sub_pp_g = df_pp_g[df_pp_g["pdb_id"].str.upper()
                               == selected_pdb]

            # Nœuds
            g_nodes = []
            for _, nr in sub_pp_g.drop_duplicates("chain").iterrows():
                is_act = bool(nr["is_actin"])
                chain_norm = norm_chain_id(str(nr["chain"]))
                g_nodes.append({
                    "id": chain_norm,
                    "label": chain_norm.split("_")[-1].upper(),
                    "title": str(nr["protein"]),
                    "color": "#E67E22" if is_act else "#2ECC71",
                    "size": 15,
                })

            # Mapping pair → interaction_id
            pair_to_iid = {}
            if sub_int1_g is not None:
                for _, er in sub_int1_g.iterrows():
                    ca = norm_chain_id(str(er["chain_A_id"]))
                    cb = norm_chain_id(str(er["chain_B_id"]))
                    pair_to_iid[frozenset([ca, cb])] = str(
                        er["interaction_id"])

            # Arêtes depuis pdb_entry_results
            g_edges = []
            edge_ids_set = set()
            sel_inter_str = str(st.session_state.get("sel_inter", ""))
            src = sub_raw_g if sub_raw_g is not None else (
                sub_int1_g.rename(
                    columns={"chain_A_id": "Interactor 1", "chain_B_id": "Interactor 2"})
                if sub_int1_g is not None else None
            )
            if src is not None:
                for _, er in src.iterrows():
                    ca = norm_chain_id(str(er["Interactor 1"]))
                    cb = norm_chain_id(str(er["Interactor 2"]))
                    iid = pair_to_iid.get(frozenset([ca, cb]))
                    eid = iid if iid else f"pair_{ca}_{cb}"
                    if iid:
                        edge_ids_set.add(eid)
                    is_sel = (eid == sel_inter_str)
                    g_edges.append({
                        "id": eid, "source": ca, "target": cb,
                        "width": 4 if is_sel else 1,
                        "color": "#FFD700" if is_sel else "#aaaaaa",
                    })

            if "sel_node" not in st.session_state:
                st.session_state["sel_node"] = None

            sel = _GRAPH_COMPONENT(
                nodes=g_nodes, edges=g_edges, height=450, key="actin_graph"
            )

            if sel is not None:
                if sel.get("type") == "edge" and sel.get("id") in edge_ids_set:
                    new_iid = int(sel["id"])
                    if new_iid != st.session_state.get("sel_inter"):
                        st.session_state["sel_inter"] = new_iid
                        st.session_state["sel_node"] = None
                        st.rerun()
                elif sel.get("type") == "node":
                    new_node = sel["id"]
                    if new_node != st.session_state.get("sel_node"):
                        st.session_state["sel_node"] = new_node
                        st.session_state["sel_inter"] = None
                        st.rerun()
            elif sel is None:
                if st.session_state.get("sel_node") or st.session_state.get("sel_inter"):
                    st.session_state["sel_node"] = None
                    st.session_state["sel_inter"] = None
                    st.rerun()

    # ── Colonne droite : viewer 3D ────────────────────────────────────────
    with col_3d:
        st.markdown("**3D visualisation — interface contacts**")
        if os.path.exists(struct_csv) and os.path.exists(cont_csv):
            df_struct8 = read_csv(struct_csv)
            sub_s8 = df_struct8[df_struct8["pdb_id"].str.upper()
                                == selected_pdb]

            pdb_file = None
            for _, r in sub_s8.iterrows():
                p = str(r.get("biological_assembly_pdb_file", "")).strip()
                if p and p.lower() != "nan" and os.path.exists(p):
                    pdb_file = p
                    break
            if pdb_file is None:
                for _, r in sub_s8.iterrows():
                    p = str(r.get("biological_assembly_cif_file", "")).strip()
                    if p and p.lower() != "nan" and os.path.exists(p):
                        pdb_file = p
                        break

            if pdb_file:
                int_ids = sub_s8["interaction_id"].tolist()
                df_cont4 = read_csv(cont_csv)

                sel_inter = st.session_state.get("sel_inter")
                sel_node = st.session_state.get("sel_node")

                # Mapping chaîne → type pour colorier orange/vert dans le viewer
                _chain_is_actin: dict[str, bool] = {}
                if os.path.exists(pp_path):
                    _df_pp_3d = read_csv(pp_path)
                    for _, _nr in _df_pp_3d[
                            _df_pp_3d["pdb_id"].str.upper() == selected_pdb
                    ].iterrows():
                        _letter = norm_chain_id(
                            str(_nr["chain"])).split("_")[-1]
                        _chain_is_actin[_letter] = bool(_nr["is_actin"])

                # chain_colors : par défaut orange=actin / vert=ABP
                # (appliqué même quand rien n'est sélectionné)
                chain_colors = {
                    _letter: ("#E67E22" if _is_actin else "#2ECC71")
                    for _letter, _is_actin in _chain_is_actin.items()
                }

                if sel_inter and sel_inter in int_ids:
                    # Cas 1 : arête sélectionnée — A jaune, B bleu (comme les
                    # séquences d'interface et la vue pairwise).
                    sel_row = df_int1_g[df_int1_g["interaction_id"]
                                        == sel_inter].iloc[0]
                    chain_colors = {
                        str(sel_row["chain_A_id"]).split("_")[-1]: "#FFD700",
                        str(sel_row["chain_B_id"]).split("_")[-1]: "#2E86FF",
                    }
                elif sel_node:
                    # Cas 2 : nœud sélectionné — sélection en jaune,
                    # autres : orange=actin, vert=ABP
                    sel_letter = sel_node.split("_")[-1]
                    chain_colors = {}
                    for _letter, _is_actin in _chain_is_actin.items():
                        if _letter == sel_letter:
                            chain_colors[_letter] = "#FFD700"
                        else:
                            chain_colors[_letter] = (
                                "#E67E22" if _is_actin else "#2ECC71")

                # Résidus des chaînes voisines TOUCHÉS par la chaîne sélectionnée
                # (jaune) → colorés en BLEU sur la surface.
                _BLUE = "#2E86FF"
                _touch: dict[str, set] = {}
                _sel_letter = None
                _rows_c = None
                if sel_inter and sel_inter in int_ids:
                    _sel_letter = str(sel_row["chain_A_id"]).split("_")[-1]
                    _rows_c = df_cont4[df_cont4["interaction_id"] == sel_inter]
                elif sel_node:
                    _sel_letter = sel_node.split("_")[-1]
                    _rows_c = df_cont4[df_cont4["interaction_id"].isin(
                        int_ids)]
                if _rows_c is not None and _sel_letter is not None:
                    for _, _cr in _rows_c.iterrows():
                        _la = str(_cr["chain_A_id"]).split("_")[-1]
                        _lb = str(_cr["chain_B_id"]).split("_")[-1]
                        try:
                            if _la == _sel_letter and _lb != _sel_letter:
                                _touch.setdefault(_lb, set()).add(
                                    int(float(_cr["residue_B_structure"])))
                            elif _lb == _sel_letter and _la != _sel_letter:
                                _touch.setdefault(_la, set()).add(
                                    int(float(_cr["residue_A_structure"])))
                        except (ValueError, TypeError):
                            pass

                fmt = "cif" if pdb_file.endswith(".cif") else "pdb"
                pdb_data = _load_pdb_file(
                    pdb_file, mtime=os.path.getmtime(pdb_file))

                viewer_key = (selected_pdb, sel_inter, sel_node)
                if st.session_state.get("viewer_key") != viewer_key:
                    import py3Dmol
                    view = py3Dmol.view(width=580, height=450)
                    view.addModel(pdb_data, fmt)
                    view.setStyle({}, {})
                    # Tout est affiché en SURFACE. Le SES (lisse) de tout un gros
                    # assemblage fige le navigateur → on l'utilise pour le contexte
                    # / le défaut seulement si l'assemblage n'est pas trop gros
                    # (≤ 15k atomes) ; au-delà, VDW (granuleux mais qui s'affiche).
                    # Les chaînes en surbrillance (peu nombreuses) restent en SES.
                    _natoms = pdb_data.count("\nATOM") + pdb_data.count("\nHETATM")
                    _ctx_surf = py3Dmol.SES if _natoms <= 15000 else py3Dmol.VDW
                    # Simple : couleurs pleines (jaune = chaîne sélectionnée,
                    # bleu = voisines), comme les séquences d'interface.
                    _surf_chains = ({_sel_letter} | set(_touch)) if _sel_letter else set()
                    if _surf_chains:
                        for chain in _surf_chains:
                            col = ("#FFD700" if chain == _sel_letter else _BLUE)
                            view.addSurface(py3Dmol.SES,
                                            {"opacity": 1.0, "color": col},
                                            {"chain": chain})
                        # contexte (chaînes non concernées) → surface blanche légère
                        for chain in [c for c in chain_colors if c not in _surf_chains]:
                            view.addSurface(_ctx_surf,
                                            {"opacity": 1.0, "color": "#e6e6e6"},
                                            {"chain": chain})
                    elif chain_colors:
                        for chain, col in chain_colors.items():
                            view.addSurface(_ctx_surf,
                                            {"opacity": 1.0, "color": col},
                                            {"chain": chain})
                    else:
                        view.addSurface(_ctx_surf,
                                        {"opacity": 0.9, "color": "#cccccc"}, {})
                    view.setBackgroundColor("#ffffff")
                    # Centrer sur les chaînes surfacées (sélection), sinon tout
                    # l'assemblage.
                    if _surf_chains:
                        view.zoomTo({"chain": list(_surf_chains)})
                    else:
                        view.zoomTo()
                    st.session_state["viewer_html"] = view._make_html()
                    st.session_state["viewer_key"] = viewer_key

                st.components.v1.html(
                    st.session_state["viewer_html"], height=470, scrolling=False)
            else:
                st.info("PDB file not available.")

    # ── Séquences avec résidus d'interface colorés ───────────────────────

    all_data_path = "data/filtered/filtered_all_data.csv"
    res3_path_seq = "data/filtered/details/3.interface_residues.csv"
    seq_inter = st.session_state.get("sel_inter")
    seq_node = st.session_state.get("sel_node")

    if (df_int1_g is not None
            and os.path.exists(all_data_path) and os.path.exists(res3_path_seq)
            and (seq_inter or seq_node)):

        def seq_html(seq, iface, col):
            W = 60
            rows = []
            for s in range(0, len(seq), W):
                chunk = seq[s:s + W]
                cells = ""
                for k, aa in enumerate(chunk):
                    _pos = s + k + 1
                    if _pos in iface:
                        _cn = iface[_pos]
                        # tooltip = UNIQUEMENT la position canonique
                        _tt = f"canonical {_cn}" if _cn is not None else ""
                        cells += (
                            f'<span title="{_tt}" style="color:{col};'
                            f'font-weight:bold'
                            f';text-shadow:0 0 4px {col}88">{aa}</span>'
                        )
                    else:
                        cells += f'<span style="color:#aaa">{aa}</span>'
                num = str(s + 1)
                pad = "&nbsp;" * (5 - len(num))
                rows.append(
                    f'<span style="color:#555;user-select:none">'
                    f'{pad}{num}&nbsp;</span>{cells}'
                )
            return (
                '<div style="font-family:\'Courier New\',monospace;'
                'font-size:13px;line-height:2;background:#161b22;'
                'padding:12px 16px;border-radius:8px;overflow-x:auto">'
                + "<br>".join(rows) + "</div>"
            )

        def render_pair(ca, cb, iid, df_all_g, df_r3,
                        col_a="#FFD700", col_b="#2E86FF", primary=None):
            ca_l, cb_l = ca.lower(), cb.lower()
            rseq = df_all_g[
                (df_all_g["subunit_1"].str.lower() == ca_l) &
                (df_all_g["subunit_2"].str.lower() == cb_l)
            ]
            if rseq.empty:
                st.info("Sequences not available for this interaction.")
                return
            s1 = str(rseq.iloc[0].get("s1_sequence", ""))
            s2 = str(rseq.iloc[0].get("s2_sequence", ""))
            t1 = str(rseq.iloc[0].get("subunit_1_title", ca))
            t2 = str(rseq.iloc[0].get("subunit_2_title", cb))
            def _ifmap(sub):
                # {position séquence : position canonique (mafft)} pour le hover
                m = {}
                for _, _rr in sub.iterrows():
                    _sp = _rr.get("residue_number_sequence")
                    if pd.notna(_sp):
                        _cn = _rr.get("residue_number_canon_mafft")
                        m[int(_sp)] = int(_cn) if pd.notna(_cn) else None
                return m
            if iid is not None:
                sub_r3 = df_r3[df_r3["interaction_id"] == iid]
                r3_chain = sub_r3["chain"].str.lower()
                iface_a = _ifmap(sub_r3[r3_chain == ca_l])
                iface_b = _ifmap(sub_r3[r3_chain == cb_l])
            else:
                iface_a, iface_b = {}, {}
            # Si primary == cb, mettre cb à gauche
            swap = primary is not None and primary.lower() == cb.lower()
            if swap:
                pairs = [(cb, t2, s2, iface_b, col_a),
                         (ca, t1, s1, iface_a, col_b)]
            else:
                pairs = [(ca, t1, s1, iface_a, col_a),
                         (cb, t2, s2, iface_b, col_b)]
            sc1, sc2 = st.columns(2)
            for widget, (cid, title, seq, iface, col) in zip([sc1, sc2], pairs):
                with widget:
                    st.markdown(
                        f"**{title}** · `{cid}` · {len(seq)} aa · {len(iface)} interface residues"
                    )
                    if seq and seq != "nan":
                        st.markdown(seq_html(seq, iface, col),
                                    unsafe_allow_html=True)

        df_all_g = read_csv(all_data_path)
        df_r3 = read_csv(res3_path_seq)

        # ── Arête sélectionnée : une paire ───────────────────────────────
        if seq_inter:
            row_i = df_int1_g[df_int1_g["interaction_id"] == seq_inter]
            if not row_i.empty:
                ca = str(row_i.iloc[0]["chain_A_id"])
                cb = str(row_i.iloc[0]["chain_B_id"])
                render_pair(ca, cb, seq_inter, df_all_g, df_r3)

        # ── Nœud sélectionné : toutes ses interactions ───────────────────
        elif seq_node:
            # Source : filtered_all_data pour tout avoir (même sans données d'interface)
            node_rows = df_all_g[
                (df_all_g["subunit_1"].str.lower() == seq_node) |
                (df_all_g["subunit_2"].str.lower() == seq_node)
            ]
            if not node_rows.empty:
                n = len(node_rows)
                st.markdown(
                    f"#### Sequences — `{seq_node.split('_')[-1].upper()}` "
                    f"({n} interaction{'s' if n > 1 else ''})"
                )
                for _, row in node_rows.iterrows():
                    ca = str(row["subunit_1"])
                    cb = str(row["subunit_2"])
                    # Chercher l'interaction_id dans 1.interactions.csv (peut être absent)
                    match = df_int1_g[
                        ((df_int1_g["chain_A_id"].str.lower() == ca.lower()) & (df_int1_g["chain_B_id"].str.lower() == cb.lower())) |
                        ((df_int1_g["chain_A_id"].str.lower() == cb.lower()) & (
                            df_int1_g["chain_B_id"].str.lower() == ca.lower()))
                    ]
                    iid = int(match.iloc[0]["interaction_id"]
                              ) if not match.empty else None
                    render_pair(ca, cb, iid, df_all_g, df_r3,
                                col_a="#FFD700", col_b="#2E86FF",
                                primary=seq_node)

    # ── Protéines impliquées (masqué quand une paire pairwise est choisie :
    #    ça ne fait pas partie des détails de la paire) ──────────────────────
    proteins_path = "data/filtered/proteins_per_pdb.csv"
    if not os.path.exists(proteins_path):
        st.info("Run the pipeline (`graphe_filter.py`) to generate the data.")
    elif not st.session_state.get("sel_inter"):
        st.markdown("**Proteins involved**")
        df_pp = read_csv(proteins_path)
        sub_pp = df_pp[df_pp["pdb_id"].str.upper() == selected_pdb]

        nb_actin = int(sub_pp["is_actin"].sum())
        nb_abp = int((~sub_pp["is_actin"]).sum())
        c1, c2, _ = st.columns([1, 1, 3])
        c1.metric("Actins", nb_actin)
        c2.metric("ABP", nb_abp)

        counts = (
            sub_pp.groupby(["protein"])
            .agg(
                Count=("chain", "count"),
                Chains=("chain", lambda x: ", ".join(
                    sorted(s.split("_")[-1].upper() for s in x)
                )),
            )
            .reset_index()
            .rename(columns={"protein": "Protein"})
            .sort_values("Count", ascending=False)
            .reset_index(drop=True)
        )
        max_len = counts["Protein"].str.len().max() if len(counts) else 0
        col_width = "large" if max_len > 40 else "medium" if max_len > 20 else "small"
        st.dataframe(counts, width="stretch", hide_index=True, column_config={
            "Protein": st.column_config.TextColumn(width=col_width),
            "Chains": st.column_config.TextColumn(width="medium"),
            "Count": st.column_config.NumberColumn(width="small"),
        })

    # ── Tableau interactions du PDB ───────────────────────────────────────
    all_data_path_pdb = "data/filtered/filtered_all_data.csv"
    if os.path.exists(all_data_path_pdb):
        df_all_pdb = read_csv(all_data_path_pdb)
        sub_all_pdb = df_all_pdb[df_all_pdb["pdb_id"].str.upper(
        ) == selected_pdb]
        # Si une paire (arête) ou un nœud est sélectionné → ne montrer que la/les
        # ligne(s) de ces contacts, pas toutes les paires du PDB.
        _tbl_lbl = ""
        _s1l = sub_all_pdb["subunit_1"].str.lower()
        _s2l = sub_all_pdb["subunit_2"].str.lower()
        _sel_i = st.session_state.get("sel_inter")
        _sel_n = st.session_state.get("sel_node")
        if _sel_i and df_int1_g is not None:
            _ri = df_int1_g[df_int1_g["interaction_id"] == _sel_i]
            if not _ri.empty:
                _pa = str(_ri.iloc[0]["chain_A_id"]).lower()
                _pb = str(_ri.iloc[0]["chain_B_id"]).lower()
                sub_all_pdb = sub_all_pdb[
                    ((_s1l == _pa) & (_s2l == _pb)) |
                    ((_s1l == _pb) & (_s2l == _pa))]
                _tbl_lbl = (f" — selected pair `{_pa.split('_')[-1].upper()}`"
                            f"–`{_pb.split('_')[-1].upper()}`")
        elif _sel_n:
            sub_all_pdb = sub_all_pdb[(_s1l == _sel_n) | (_s2l == _sel_n)]
            _tbl_lbl = f" — chain `{_sel_n.split('_')[-1].upper()}`"
        cols_to_show = [c for c in [
            "subunit_1", "subunit_2",
            "cluster_data_70",
            "s1_binding_site_cluster_data_70",
            "s2_binding_site_cluster_data_70",
        ] if c in sub_all_pdb.columns]
        if not sub_all_pdb.empty and cols_to_show:
            # Paire sélectionnée → une seule ligne (celle de la paire) ; nœud →
            # ses lignes ; rien → toutes les paires du PDB.
            st.markdown("**Interaction clusters per pair**" + _tbl_lbl)
            st.dataframe(
                sub_all_pdb[cols_to_show].reset_index(drop=True),
                hide_index=True,
                width="stretch",
            )
            st.divider()

# ── Récap structural par cluster (Foldseek / TM-align / InterPro) ──────────────


# ── ProteoCast par ABP (paysage mutationnel + empreinte d'interface) ───────────


st.divider()

# ── Empreinte ABP sur l'actin : grand actin coloré par nb d'ABP + clic résidu ──
if not _FAST:
    st.header("ABP footprint on actin", anchor="empreinte-abp")
    _pp_ov = build_passport(pp_mtimes())
    if _pp_ov is None:
        st.info("Residue-passport data unavailable "
                "(run the pipeline to generate them).")
    else:
        residue_explorer.render_actin_overview(_pp_ov)

st.divider()
st.header("Interaction clusters", anchor="clusters-d-interactions")
# Défilement auto vers cette section après un clic dans le mini-graphe (explorateur).
if st.session_state.pop("_scroll_clusters", False):
    st.components.v1.html(
        "<script>window.parent.location.hash = 'clusters-d-interactions';"
        "</script>", height=0)

PATCHES_S1_CSV = "data/filtered/patches_infos_s1_binding_site.csv"
PATCHES_C70_CSV = "data/filtered/patches_infos_cluster_data_70.csv"
GLOBAL_GRAPH_IMG = "data/visualisations/patch_graphs/global.png"
GLOBAL_HEATMAP_IMG = "data/visualisations/patch_graphs/heatmap_binary.png"
S1_GRAPHS_DIR = "data/visualisations/patch_graphs/s1_binding_site"
C70_GRAPHS_DIR = "data/visualisations/patch_graphs/cluster_data_70"
S1_HEATMAPS_F4_DIR = "data/visualisations/patch_heatmaps_s1_contacts"
C70_HEATMAPS_F4_DIR = "data/visualisations/patch_heatmaps_c70_contacts"
C70_CONTACTS_DIR = "data/visualisations/actin_c70_contacts"
C70_CONTACTS_SURFACE_AREA_DIR = "data/visualisations/actin_c70_contacts_surface_area"
C70_HEATMAP_SURFACE_AREA = "data/visualisations/actin_c70_heatmap_surface_area.png"
S1_HEATMAP_RAW = "data/visualisations/actin_s1_homo_used_heatmap.png"
S1_HEATMAP_EQUITABLE = "data/visualisations/actin_s1_all_equitable_heatmap.png"
S1_HEATMAP_ABSOLUTE = "data/visualisations/actin_s1_heatmap_absolute.png"
S1_INTERFACE_FREQ_CSV = "data/filtered/actin_s1_canon_area_by_cluster.csv"
S1_INTERFACE_CLUSTER_DIR = "data/visualisations/actin_s1_clusters"
S1_INTERFACE_CLUSTER_BY_C70_DIR = "data/visualisations/actin_s1_clusters_by_c70"
BFACTOR_CLUSTER_DIR = "data/filtered/details/structures_files/bfactor_cluster"


# Longueur de l'actin canonical : les positions au-delà sont des artefacts
# d'alignement (insertions MSA) et ne sont pas affichées sur l'axe des positions.


if not os.path.exists(PATCHES_C70_CSV):
    st.info("Run the pipeline to generate the cluster analysis.")
else:
    _all_data_path = "data/filtered/filtered_all_data.csv"
    _summary_path = "data/filtered/filtered_summary.csv"
    tab_s1, tab_c70 = st.tabs(
        ["Binding Site Cluster Data 70", "Interaction Cluster Data 70"])

    # --- S1 Binding Site ---
    with tab_s1:
        if os.path.exists(_all_data_path) and os.path.exists(_summary_path):
            _use_super = st.toggle(
                "Group actin sites into super-clusters",
                value=False, key="global_graph_superclusters",
                help="Merges actin binding sites whose residues "
                     "canonical se chevauchent (colonnes s1/s2_supercluster) "
                     "into a single representative node.",
            )
            st.components.v1.html(_build_global_graph_html(
                _all_data_path, _summary_path, _use_super), height=780)
            st.divider()
        if os.path.exists(PATCHES_S1_CSV):

            df_s1 = read_csv(PATCHES_S1_CSV)

            df_s1_display = df_s1.drop(
                columns=["ids_interactions"], errors="ignore").copy()

            # Colonnes C70 : clusters d'interaction contenant les interactions du patch S1
            if os.path.exists(PATCHES_C70_CSV):
                import re as _re

                def parse_ids(s):
                    s = str(s)
                    s = _re.sub(r'np\.int64\((\d+)\)', r'\1', s)
                    return set(int(x) for x in _re.findall(r'\d+', s))

                df_c70_tmp = read_csv(PATCHES_C70_CSV)
                # Mapping interaction_id → C70 patch
                id_to_c70 = {}
                for _, r in df_c70_tmp.iterrows():
                    for iid in parse_ids(r['ids_interactions']):
                        id_to_c70[iid] = str(r['patch'])

                c70_clusters_col, n_c70_col = [], []
                for _, row in df_s1.iterrows():
                    iids = parse_ids(row['ids_interactions'])
                    c70s = sorted(set(id_to_c70[i]
                                      for i in iids if i in id_to_c70))
                    c70_clusters_col.append(", ".join(c70s))
                    n_c70_col.append(len(c70s))

                df_s1_display['n_c70_clusters'] = n_c70_col
                df_s1_display['c70_clusters'] = c70_clusters_col

            df_s1_display.index = range(1, len(df_s1_display) + 1)
            st.dataframe(df_s1_display, width="stretch")

            # --- Heatmap global homo + hétéro (repliée : prend de la place) ---
            with st.expander(
                    "Heatmap — canonical actin positions at the S1 interface",
                    expanded=False):
                _s1g_data = _build_s1_global_heatmap(
                    tuple(os.path.getmtime(f) if os.path.exists(f) else 0.0
                          for f in _S1_GLOBAL_FILES))
                if _s1g_data is None:
                    st.info(
                        "Data unavailable (run the pipeline / "
                        "`notebooks/interface_analysis_s1.py`).")
                else:
                    heatmap_mode = st.selectbox(
                        "Normalisation",
                        ["Relative (max cluster = 1)",
                         "Absolute (max = 100 %)"],
                        index=1,   # défaut : valeurs absolues (% ASA)
                        key="heatmap_norm_mode",
                    )
                    _rel = heatmap_mode == "Relative (max cluster = 1)"
                    _render_s1_global_plotly(
                        _s1g_data, _rel,
                        valid_clusters=set(df_s1["patch"].astype(str)))

            st.divider()

            # --- Cluster sélectionné ---
            all_s1 = df_s1.sort_values("n_interactions", ascending=False)[
                "patch"].astype(str).tolist()
            # Navigation depuis le mini-graphe de l'explorateur (si valide).
            _pend_s1 = st.session_state.pop("_pending_s1", None)
            if _pend_s1 in all_s1:
                st.session_state["sel_s1"] = _pend_s1

            # Navigation depuis le graphe : un clic sur un nœud rouge clique (en
            # JS, via window.top) le bouton caché du cluster ci-dessous, qui
            # présélectionne ce cluster. On masque ces boutons (classe st-key-__pt_*).
            st.markdown(
                "<div id='patchsel'></div>"
                "<style>[class*='st-key-__pt_']{position:absolute!important;"
                "left:-9999px!important;top:auto!important;margin:0!important}</style>",
                unsafe_allow_html=True,
            )
            # Défilement auto vers le détail après un clic sur la heatmap globale
            if st.session_state.pop("_scroll_to_s1", False):
                st.components.v1.html(
                    "<script>const d=window.parent.document;"
                    "const el=d.getElementById('patchsel');"
                    "if(el){el.scrollIntoView({behavior:'smooth',block:'start'});}"
                    "</script>", height=0)
            for _cid in all_s1:
                if st.button(_cid, key=f"__pt_{_cid}"):
                    st.session_state["sel_s1"] = _cid

            sel_s1 = st.selectbox("Patch S1 binding site", all_s1, key="sel_s1",
                                  format_func=lambda p: f"{p} — {int(df_s1[df_s1['patch'].astype(str) == p]['n_interactions'].values[0])} interactions")

            row_s1 = df_s1[df_s1["patch"].astype(str) == sel_s1].iloc[0]

            # ── Réseau bipartite interactif + Surface 3D ─────────────────
            col_net_s1, col_3d_s1 = st.columns([3, 3])

            with col_net_s1:
                st.markdown(
                    "**Interactive network — actin residues ↔ partners**")
                _bip_ok = all(os.path.exists(f) for f in _BIPARTITE_FILES)
                if _bip_ok:
                    _html_bip, _n_r, _n_p, _n_t = _build_bipartite_html(
                        sel_s1, _BIP_CACHE_VERSION, *_bip_mtimes())
                    if _html_bip:
                        st.components.v1.html(
                            _html_bip, height=540, scrolling=False)
                    else:
                        st.info("Network not available for this cluster.")
                else:
                    st.info("Run the pipeline to generate the network data.")

            with col_3d_s1:
                # ── Interface 3D interactive : choisir un ABP, ou « tout ensemble » ──
                st.markdown("**Interface 3D — actin ↔ ABP**")
                _abp3d = _s1_abp_3d_options(
                    sel_s1, os.path.getmtime("data/filtered/filtered_all_data.csv"))
                _opts = [{"label": "cluster sur actin", "mode": "none"}]
                if _abp3d:
                    _opts.append({"label": "— All together —", "mode": "all"})
                    _opts += [dict(o, mode="one") for o in _abp3d]
                _sel3d = st.selectbox("Afficher", _opts,
                                      format_func=lambda o: o["label"], key=f"abp3d_v2_{sel_s1}")
                _mode3d = _sel3d.get("mode", "none")
                if _mode3d == "none":
                    _s1solo = _build_s1_3d_html(
                        sel_s1, _BIP_CACHE_VERSION, *_bip_mtimes())
                    if _s1solo:
                        _hs, _smax = _s1solo[0], _s1solo[1]
                        st.components.v1.html(_hs, height=490, scrolling=False)
                    else:
                        pass
                elif _mode3d == "all":
                    _paths = tuple(o["pdb"] for o in _abp3d)
                    _labs = tuple(o["label"] for o in _abp3d)
                    _pdb_c, _ac, _cmap = _build_all_abp_pdb(_paths, _labs)
                    if _cmap:
                        _h3d, _bm = _build_all_abp_3d(_pdb_c, _ac, _cmap)
                        st.components.v1.html(
                            _h3d, height=490, scrolling=False)
                        _leg = " · ".join(
                            f"<span style='color:{_ABP_MULTI_COLORS[i % len(_ABP_MULTI_COLORS)]};"
                            f"font-size:16px'>■</span> {l}"
                            for i, (c, l) in enumerate(_cmap.items()))
                        st.markdown(_leg, unsafe_allow_html=True)
                    else:
                        st.info("Superposition not possible for this site.")
                else:
                    _r3d = _build_abp_actin_3d(_sel3d["pdb"])
                    if _r3d:
                        _h3d, _amax, _bmax = _r3d
                        st.components.v1.html(
                            _h3d, height=490, scrolling=False)
                    else:
                        st.info("3D structure unreadable for this ABP.")

                # Téléchargement du script PyMOL — placé sous la représentation 3D
                _pml_path = os.path.join(
                    "data/filtered/details/structures_files/bfactor_c70_interface/by_s1_gradient",
                    f"{sel_s1}.pml")
                if os.path.exists(_pml_path):
                    with open(_pml_path, "rb") as _pml_f:
                        st.download_button(
                            label="Download the PyMOL script for this cluster",
                            data=_pml_f,
                            file_name=f"{sel_s1}.pml",
                            mime="text/plain",
                            help="Open in PyMOL: File > Run Script…  or  @/path/to/the/file.pml",
                        )

            # ── Détail par position : aa d'actin par organisme + aa d'ABP ──
            _posdet = _s1_position_detail(
                sel_s1, tuple(os.path.getmtime(f) if os.path.exists(f) else 0.0
                              for f in _S1_GLOBAL_FILES))
            if _posdet is not None:
                _render_s1_position_detail(_posdet, sel_s1)

            # Profils interactifs du cluster (positions à leur place + %ASA au
            # survol). Recalculés depuis les données -> la décomposition inclut
            # TOUS les sous-clusters C70 et se met à jour au re-téléchargement.
            _pdetail = _build_s1_patch_detail(
                sel_s1, tuple(os.path.getmtime(f) if os.path.exists(f) else 0.0
                              for f in _S1_GLOBAL_FILES))
            if _pdetail is None:
                st.info("Profile unavailable for this cluster.")
            else:
                _render_s1_patch_plotly(_pdetail, sel_s1)

            # ── Analyse des contacts ABP–actin ──────────────────────────────
            # « MSA du cluster » et « Comparaison structurale » sont des onglets
            # supplémentaires insérés juste avant « Spécificité ».
            st.markdown("**ABP–actin contact analysis**")
            _ch2seq_s1, _ch2title_s1 = _s1_get_ch_maps(sel_s1)
            if _ch2seq_s1:
                def _tab_msa_cluster(_c=sel_s1):
                    st.markdown(
                        "**Cluster MSA** — alignment of this cluster's sequences")
                    _msa_s1_cluster_inline(_c, include_contacts=False)

                def _tab_struct(_c=sel_s1):
                    st.markdown(
                        "**Structural comparison of this site's ABPs**")
                    _render_cluster_struct(_c)

                _msa_contact_analysis(
                    None, f"s1c_{sel_s1}",
                    _ch2seq=_ch2seq_s1, _ch2title=_ch2title_s1,
                    tabs=["C", "D", "E"],
                    extra_tabs=[
                        ("MSA du cluster", _tab_msa_cluster),
                        ("Comparaison structurale", _tab_struct),
                    ],
                )
            else:
                st.info("No contact data available for this cluster.")

    # --- Cluster Data 70 ---
    with tab_c70:
        if os.path.exists(_all_data_path):
            st.components.v1.html(_build_tripartite_graph_html(
                _all_data_path), height=900)
            st.divider()
        if os.path.exists(PATCHES_C70_CSV):
            df_c70 = read_csv(PATCHES_C70_CSV)

            df_c70_display = df_c70.drop(
                columns=["ids_interactions"], errors="ignore").copy()

            df_c70_display.index = range(1, len(df_c70_display) + 1)
            st.dataframe(df_c70_display, width="stretch")

            st.divider()
            all_c70 = df_c70.sort_values("n_interactions", ascending=False)[
                "patch"].astype(str).tolist()
            # Navigation depuis le mini-graphe de l'explorateur (si valide).
            _pend_c70 = st.session_state.pop("_pending_c70", None)
            if _pend_c70 in all_c70:
                st.session_state["sel_c70"] = _pend_c70
            sel_c70 = st.selectbox("Patch cluster_data_70", all_c70, key="sel_c70",
                                   format_func=lambda p: f"{p} — {int(df_c70[df_c70['patch'].astype(str) == p]['n_interactions'].values[0])} interactions")

            row_c70 = df_c70[df_c70["patch"].astype(str) == sel_c70].iloc[0]

            # Clusters S2 depuis df_all — source de vérité pour tout ce qui suit
            df_all_c70_sel = pd.DataFrame()
            if os.path.exists("data/filtered/filtered_all_data.csv"):
                df_all_c70_sel = read_csv(
                    "data/filtered/filtered_all_data.csv")
                df_all_c70_sel = df_all_c70_sel[
                    df_all_c70_sel["cluster_data_70"].astype(
                        str) == str(sel_c70)
                ]
            s2_counts = (
                df_all_c70_sel["s2_sequence_cluster_70"].dropna().astype(int)
                .value_counts().sort_values(ascending=False)
                if not df_all_c70_sel.empty else pd.Series(dtype=int)
            )

            st.markdown(
                "**Interactive network — actin residues ↔ ABP residues**")
            _bip_ok = all(os.path.exists(f) for f in _BIPARTITE_FILES)
            if _bip_ok:
                _color_mode = "restype" if st.toggle(
                    "Physicochemical colouring (hydrophobic/polar/charged…)",
                    value=False, key=f"restype_toggle_{sel_c70}"
                ) else "bfactor"
                (_html_c70, _n_s1, _n_s2, _n_tot, _html_3d_c70,
                 _cmat_c70) = _build_bipartite_c70_html(
                    sel_c70, True, _color_mode, _BIP_CACHE_VERSION, *_bip_mtimes())
                if _html_c70:
                    _net_height = 650
                    if _html_3d_c70:
                        # Biparti + 3D côte à côte (biparti plus large).
                        _cbip, _c3d = st.columns([3, 2])
                        with _cbip:
                            st.components.v1.html(
                                _html_c70, height=_net_height, scrolling=False)
                        with _c3d:
                            st.markdown("**3D interface — representative pair**")
                            st.components.v1.html(
                                _html_3d_c70, height=_net_height - 40,
                                scrolling=False)
                    else:
                        st.components.v1.html(
                            _html_c70, height=_net_height, scrolling=False)
                    # Les 2 séquences de la paire représentative colorées par
                    # % ASA enfouie (au lieu d'une matrice de contact).
                    if _cmat_c70 and _cmat_c70.get("s1_seq"):
                        def _seq_asa_html(seq, asa):
                            _pal = ["#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C",
                                    "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026",
                                    "#800026"]
                            _mx = max(asa.values(), default=1.0) or 1.0
                            W, rows = 60, []
                            for s in range(0, len(seq), W):
                                cells = ""
                                for k, aa in enumerate(seq[s:s + W]):
                                    pos = s + k + 1
                                    if pos in asa:
                                        v = asa[pos]
                                        _c = _pal[min(len(_pal) - 1,
                                                      int(v / _mx * len(_pal)))]
                                        cells += (f'<span title="pos {pos} · ASA '
                                                  f'{v:.0f}%" style="background:{_c};'
                                                  f'color:#111;font-weight:bold">'
                                                  f'{aa}</span>')
                                    else:
                                        cells += f'<span style="color:#999">{aa}</span>'
                                num = str(s + 1)
                                pad = "&nbsp;" * (5 - len(num))
                                rows.append(f'<span style="color:#888">{pad}{num}'
                                            f'&nbsp;</span>{cells}')
                            return ('<div style="font-family:monospace;font-size:13px;'
                                    'line-height:1.9;background:#fff;padding:10px;'
                                    'border:1px solid #eee;border-radius:6px;'
                                    'overflow-x:auto">' + "<br>".join(rows) + "</div>")

                        st.markdown("**Interface sequences — buried %ASA "
                                    "(representative pair)**")
                        _qa1, _qa2 = st.columns(2)
                        with _qa1:
                            st.markdown("Actin (S1)")
                            st.markdown(_seq_asa_html(_cmat_c70["s1_seq"],
                                                      _cmat_c70["s1_asa"]),
                                        unsafe_allow_html=True)
                        with _qa2:
                            st.markdown("Partner (S2)")
                            st.markdown(_seq_asa_html(_cmat_c70["s2_seq"],
                                                      _cmat_c70["s2_asa"]),
                                        unsafe_allow_html=True)
                else:
                    st.info("Network not available for this cluster.")
            else:
                st.info("Run the pipeline to generate the network data.")

st.divider()
st.header("ABP")

# Vue globale des protéines non-actins
proteins_path = "data/filtered/proteins_per_pdb.csv"
_all_data_path = "data/filtered/filtered_all_data.csv"
if os.path.exists(proteins_path):
    st.subheader("Actin-binding proteins (ABPs) — overview")
    df_pp_all = read_csv(proteins_path)
    df_abp = df_pp_all[~df_pp_all["is_actin"]]

    def _fmt_ids(series):
        vals = sorted(set(str(v).strip()
                      for v in series.dropna() if str(v) not in ("nan", "")))
        return ", ".join(vals)

    # Labels de position filament pour les actins
    _fil_pos_path = "data/filtered/actin_filament_positions.csv"
    _fil_label_map = {}
    if os.path.exists(_fil_pos_path):
        _fp = read_csv(_fil_pos_path)[
            ["pdb_id", "chain", "label", "component_size"]]
        for _, _r in _fp.iterrows():
            if _r["component_size"] > 3:
                _fil_label_map[(_r["pdb_id"], _r["chain"])] = _r["label"]

    _label_order = {"-": 0, "-2": 1, "-3": 2,
                    "side": 3, "+3": 4, "+2": 5, "+": 6}

    def _simplify_label(l):
        # seuls les bouts terminaux stricts comptent comme +/- ;
        # +2,+3,-2,-3 et side sont regroupés en latéral ("side")
        if l == "+":
            return "+"
        if l == "-":
            return "-"
        return "side"

    def _fmt_fil_labels(pairs):
        """pairs : série de (pdb_id, actin_chain) tuples — retourne les combinaisons par PDB"""
        pdb_lbls: dict = {}
        for _pk in pairs:
            _l = _fil_label_map.get(_pk)
            if _l:
                pdb_lbls.setdefault(_pk[0], set()).add(_l)
        if not pdb_lbls:
            return None
        combos: set = set()
        for _s in pdb_lbls.values():
            _c = ", ".join(sorted(_s, key=lambda x: _label_order.get(x, 99)))
            if _c:
                combos.add(_c)
        return " ".join(f"({c})" for c in sorted(combos, key=lambda c: [_label_order.get(x.strip(), 99) for x in c.split(",")]))

    def _fmt_fil_summary(pairs):
        """Résumé simplifié : labels uniques +/-/side toutes PDBs confondues"""
        _ls = set()
        for _pk in pairs:
            _l = _fil_label_map.get(_pk)
            if _l:
                _ls.add(_simplify_label(_l))
        if not _ls:
            return None
        return ", ".join(sorted(_ls, key=lambda x: {"-": 0, "side": 1, "+": 2}.get(x, 99)))

    if os.path.exists(_all_data_path):
        df_all_g = read_csv(_all_data_path)
        df_all_g["s1_actine"] = df_all_g["s1_actine"].fillna(
            False).astype(bool)
        df_all_g["s2_actine"] = df_all_g["s2_actine"].fillna(
            False).astype(bool)
        # Direction 1 : ABP = subunit_2, actin = subunit_1
        _h1 = df_all_g[df_all_g["s1_actine"] & ~df_all_g["s2_actine"]][[
            "pdb_id", "subunit_1", "subunit_2",
            "s1_binding_site_cluster_data_70", "cluster_data_70",
        ]].rename(columns={"subunit_2": "_abp_chain",
                           "subunit_1": "_actin_chain",
                           "s1_binding_site_cluster_data_70": "_s1_bs"})
        # Direction 2 : ABP = subunit_1, actin = subunit_2
        _h2 = df_all_g[~df_all_g["s1_actine"] & df_all_g["s2_actine"]][[
            "pdb_id", "subunit_2", "subunit_1",
            "s2_binding_site_cluster_data_70", "cluster_data_70",
        ]].rename(columns={"subunit_1": "_abp_chain",
                           "subunit_2": "_actin_chain",
                           "s2_binding_site_cluster_data_70": "_s1_bs"})
        hetero = pd.concat([_h1, _h2], ignore_index=True)
        hetero["_abp_chain"] = hetero["_abp_chain"].str.lower()
        hetero["_actin_key"] = list(
            zip(hetero["pdb_id"], hetero["_actin_chain"]))
        hetero = hetero.drop(columns=["pdb_id", "_actin_chain"])
        df_abp = df_abp.copy()
        df_abp["chain_low"] = df_abp["chain"].str.lower()
        merged = df_abp.merge(hetero, left_on="chain_low",
                              right_on="_abp_chain", how="inner")
        abp_global = (
            merged.groupby("protein")
            .agg(
                Nb_noeuds=("chain", "nunique"),
                PDB=("pdb_id", lambda x: ", ".join(sorted(x.unique()))),
                **{"Binding site S1": ("_s1_bs", _fmt_ids)},
                **{"Cluster C70": ("cluster_data_70", _fmt_ids)},
                **{"Position": ("_actin_key", _fmt_fil_summary)},
                **{"Localisation filament": ("_actin_key", _fmt_fil_labels)},
            )
            .reset_index()
            .rename(columns={"protein": "Protein", "Nb_noeuds": "# nodes"})
            .sort_values("# nodes", ascending=False)
            .reset_index(drop=True)
        )
    else:
        abp_global = (
            df_abp.groupby("protein")
            .agg(
                Nb_noeuds=("chain", "nunique"),
                PDB=("pdb_id", lambda x: ", ".join(sorted(x.unique()))),
            )
            .reset_index()
            .rename(columns={"protein": "Protein", "Nb_noeuds": "# nodes"})
            .sort_values("# nodes", ascending=False)
            .reset_index(drop=True)
        )

    max_len = abp_global["Protein"].str.len().max() if len(abp_global) else 0
    col_width = "large" if max_len > 40 else "medium" if max_len > 20 else "small"
    st.dataframe(abp_global, width="stretch", hide_index=True, column_config={
        "Protein": st.column_config.TextColumn(width=col_width),
        "# nodes": st.column_config.NumberColumn(width="small"),
        "PDB": st.column_config.TextColumn(width="large"),
        "Binding site S1": st.column_config.TextColumn(width="medium"),
        "Cluster C70": st.column_config.TextColumn(width="medium"),
        "Localisation filament": st.column_config.TextColumn(width="medium"),
    })

# ── Heatmap ABP × résidus actin ──────────────────────────────────────────────
st.subheader("Heatmap — actin residues contacted by ABPs")


if not _FAST and all(os.path.exists(f) for f in _ABP_HM_FILES):
    _hm_mtimes = tuple(
        os.path.getmtime(f) if os.path.exists(f) else 0.0
        for f in _ABP_HM_FILES
    )
    _hm_res = _build_abp_heatmap_data(*_hm_mtimes)
    _hm_ok = (_hm_res is not None and _hm_res[0] is not None
              and not _hm_res[0].empty and _hm_res[0].shape[1] > 0)
    if _hm_res is not None and not _hm_ok:
        st.info("No ABP–actin contact data to display — the derived analyses are "
                "out of date vs the current interactions. Re-run **Run / update** "
                "to regenerate them.")
    if _hm_ok:
        _pivot_full_hm, _abp_freq_hm, _ = _hm_res

        _pivot_hm = _pivot_full_hm.copy()

        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        _n_rows = len(_pivot_hm)

        # Clustering hiérarchique sur les lignes (ABP) — seulement si ≥ 2 lignes,
        # sinon pdist plante ("empty distance matrix"). En cas de jeu vide, on
        # affiche un message plutôt que de crasher.
        _mat = _pivot_hm.fillna(0).values.astype(float)
        if _n_rows >= 2:
            _row_order = leaves_list(
                linkage(pdist(_mat, metric="cosine"), method="average"))
            _pivot_hm = _pivot_hm.iloc[_row_order]

        # Colonnes à leur VRAIE position canonical : on comble les trous
        # (positions non contactées) pour que l'axe X reflète la séquence d'actin.
        _cols_int = [int(c) for c in _pivot_hm.columns]
        _full_cols = list(range(min(_cols_int), max(_cols_int) + 1))
        _piv_ord = _pivot_hm.copy()
        _piv_ord.columns = _cols_int
        _piv_ord = _piv_ord.reindex(columns=_full_cols)  # NaN = pas de contact

        # Bande actin : nb d'ABP qui touchent chaque position (empreinte globale)
        _nabp_pos = (_pivot_hm.notna().sum(axis=0))
        _nabp_pos.index = _cols_int
        _nabp_pos = _nabp_pos.reindex(_full_cols, fill_value=0)

        # Labels Y avec n interactions
        _ylabels = [f"{r}  ({_abp_freq_hm.get(r, 0)})" for r in _piv_ord.index]

        _fig_hm = make_subplots(
            rows=2, cols=1, shared_xaxes=True,
            row_heights=[0.86, 0.14], vertical_spacing=0.02,
        )
        _fig_hm.add_trace(go.Heatmap(
            z=_piv_ord.values, x=_full_cols, y=_ylabels,
            colorscale="YlOrRd", zmin=0, zmax=100,
            colorbar=dict(title="% ASA moy", thickness=12, len=0.86, y=0.57),
            hovertemplate=("ABP : %{y}<br>Position : %{x}"
                           "<br>%ASA moy : %{z:.2f}<extra></extra>"),
            hoverongaps=False,
        ), row=1, col=1)
        _fig_hm.add_trace(go.Bar(
            x=_full_cols, y=_nabp_pos.values,
            marker=dict(color=_nabp_pos.values, colorscale="Blues",
                        cmin=0, showscale=False),
            hovertemplate=("Position: %{x}<br>ABPs touching it: %{y}"
                           "<extra></extra>"),
        ), row=2, col=1)
        _fig_hm.update_layout(
            height=max(520, _n_rows * 15 + 200),
            margin=dict(l=4, r=4, t=28, b=44), bargap=0,
            title=dict(text=f"mean buried % ASA — actin residues × ABP  "
                       f"({_n_rows} ABPs, {len(_cols_int)} contacted positions)",
                       font=dict(size=13)),
        )
        _fig_hm.update_yaxes(autorange="reversed", tickfont=dict(size=9),
                             row=1, col=1)
        _fig_hm.update_yaxes(title_text="n ABP", title_font=dict(size=10),
                             row=2, col=1)
        _fig_hm.update_xaxes(dtick=25, row=2, col=1,
                             title_text="Position canonical de l'actin (MAFFT)",
                             title_font=dict(size=11))
        st.plotly_chart(_fig_hm, use_container_width=True)

# ── Réseaux ABP : compétition OU coopération (deux réseaux indépendants, au choix) ──
_net_view = None
_coop_pairs = {}
_coop_set = set()
if (os.path.exists(proteins_path) and os.path.exists(_all_data_path)
        and len(abp_global) > 0 and "Binding site S1" in abp_global.columns):
    st.subheader("ABP networks")
    _net_view = st.radio(
        "Network to display", ["Competition", "Cooperation"],
        horizontal=True, key="net_view",
        help="Competition = ABPs fighting over the same site (overlapping footprints). "
             "Cooperation = ABPs co-present in the same PDB (coexist on actin).",
    )
    # Paires coopérantes (co-présence dans un même PDB) — servent aux deux vues
    from itertools import combinations as _comb_coop
    _abp_pp_co = df_pp_all[~df_pp_all["is_actin"]]
    _pdb_to_abps = _abp_pp_co.groupby("pdb_id")["protein"].apply(
        lambda s: sorted(set(s.dropna())))
    for _prots_co in _pdb_to_abps:
        for _a_co, _b_co in _comb_coop(_prots_co, 2):
            _kk = frozenset((_a_co, _b_co))
            _coop_pairs[_kk] = _coop_pairs.get(_kk, 0) + 1
    _coop_set = set(_coop_pairs)


# ── Réseau de compétition ABP ─────────────────────────────────────────────────
if (_net_view == "Competition" and os.path.exists(proteins_path)
        and os.path.exists(_all_data_path)
        and len(abp_global) > 0 and "Binding site S1" in abp_global.columns):
    _jac_thresh = st.select_slider(
        "Seuil de recouvrement C70",
        options=[0.25, 0.30, 0.40, 0.50, 0.60,
                 0.70, 0.75, 0.80, 0.90, 1.00],
        value=0.50, key="jac_thresh",
        format_func=lambda x: f"{int(x*100)}%",
    )
    _min_edge_w = 1
    _show_coop = False   # la coopération est une vue séparée (choix ci-dessus)
    # Dataset étendu : actin peut être S1 ou S2 (toutes données hétéro)
    _df_s1_act = df_all_g[
        df_all_g["s1_actine"].fillna(False).astype(bool) &
        ~df_all_g["s2_actine"].fillna(False).astype(bool)
    ][["subunit_2", "s1_binding_site_cluster_data_70", "s2_sequence_cluster_70",
       "cluster_data_70"]].rename(
        columns={"subunit_2": "_ch", "s1_binding_site_cluster_data_70": "_bs",
                 "s2_sequence_cluster_70": "_seqcl", "cluster_data_70": "_c70"})
    _df_s2_act = df_all_g[
        ~df_all_g["s1_actine"].fillna(False).astype(bool) &
        df_all_g["s2_actine"].fillna(False).astype(bool)
    ][["subunit_1", "s2_binding_site_cluster_data_70", "s1_sequence_cluster_70",
       "cluster_data_70"]].rename(
        columns={"subunit_1": "_ch", "s2_binding_site_cluster_data_70": "_bs",
                 "s1_sequence_cluster_70": "_seqcl", "cluster_data_70": "_c70"})
    _hetero_net = pd.concat([_df_s1_act, _df_s2_act], ignore_index=True)
    _hetero_net["_ch"] = _hetero_net["_ch"].str.lower()
    _merged_net = _hetero_net.merge(
        df_abp[["chain_low", "protein"]], left_on="_ch", right_on="chain_low", how="inner"
    )
    # Nombre de clusters C70 uniques par ABP
    _abp_c70_count = (
        _merged_net.dropna(subset=["_c70"])
        .groupby("protein")["_c70"].nunique()
        .to_dict()
    )

    _sc_mtimes = tuple(
        os.path.getmtime(f) if os.path.exists(f) else 0.0
        for f in _S1_SUPER_FILES
    )
    _edge_wts, _abp_c70_comp_count, _c70_footprint, _abp_to_c70s_da = _build_c70_jaccard_edges(
        _jac_thresh, *_sc_mtimes)

    # Filtre localisation filament au niveau C70 :
    # une arête (A,B) est gardée seulement si au moins une paire (ca∈A, cb∈B)
    # a simultanément Jaccard ≥ seuil ET chevauchement de position filament.
    # Plus précis que le filtre ABP-level : évite de garder une arête à cause
    # d'un autre C70 de l'ABP qui, lui, chevauche en position.
    _c70_fil_pos: dict = {}
    for _, _r in merged.dropna(subset=["cluster_data_70"]).iterrows():
        _c70k = str(_r["cluster_data_70"])
        _lbl = _fil_label_map.get(_r["_actin_key"])
        if _lbl:
            _c70_fil_pos.setdefault(_c70k, set()).add(_lbl)

    _n_before = len(_edge_wts)
    _edge_wts_filtered: dict = {}
    _abp_c70_comp_count_filtered: dict = {}  # recalculé après filtre position
    # ABP → set de C70 valides (Jaccard + position)
    _abp_c70_valid: dict = {}
    for (_a, _b) in _edge_wts:
        # Recompter uniquement les paires C70 qui passent Jaccard ET position
        # Paire coopérante (même PDB) : on garde l'arête (dessinée en vert) mais on
        # ne compte PAS ses C70 comme "contestés" -> couleur du nœud = vraie compétition.
        _pair_coop = _show_coop and frozenset((_a, _b)) in _coop_set
        _n_valid = 0
        for _ca in _abp_to_c70s_da.get(_a, set()):
            _fa = _c70_footprint.get(_ca)
            if not _fa:
                continue
            _pos_a = _c70_fil_pos.get(_ca, set())
            for _cb in _abp_to_c70s_da.get(_b, set()):
                _fb = _c70_footprint.get(_cb)
                if not _fb:
                    continue
                _mn = min(len(_fa), len(_fb))
                if _mn > 0 and len(_fa & _fb) / _mn >= _jac_thresh:
                    if _pos_a & _c70_fil_pos.get(_cb, set()):
                        _n_valid += 1
                        if not _pair_coop:
                            _abp_c70_valid.setdefault(_a, set()).add(_ca)
                            _abp_c70_valid.setdefault(_b, set()).add(_cb)
        if _n_valid > 0:
            _edge_wts_filtered[(_a, _b)] = _n_valid
    _n_filtered = _n_before - len(_edge_wts_filtered)
    _edge_wts = _edge_wts_filtered
    # Remplacer _abp_c70_comp_count par le compteur filtré (pour la couleur des noeuds)
    _abp_c70_comp_count = {k: len(v) for k, v in _abp_c70_valid.items()}

    # Position filament simplifiée (+/-/side) par ABP — calculée depuis merged
    _abp_fil_pos_net: dict = {}
    for _, _r in merged.iterrows():
        _lbl = _fil_label_map.get(_r["_actin_key"])
        if _lbl:
            _abp_fil_pos_net.setdefault(
                _r["protein"], set()).add(_simplify_label(_lbl))

    def _pos_border_color(pos_set: set) -> str:
        """Orange = barbée (+ avec ou sans side) ; bleu = pointée (− avec ou sans side) ;
        gris = side seul, mixte +/−, ou inconnu."""
        has_p = "+" in pos_set
        has_m = "-" in pos_set
        if has_p and not has_m:
            return "#e67300"   # orange — barbée (seul ou + side)
        if has_m and not has_p:
            return "#0077cc"   # bleu   — pointée (seul ou + side)
        return "#888888"       # gris   — side seul, barbée+pointée, ou inconnu


    # Tous les ABPs connus (y compris isolés)
    _all_abp_net = set(_merged_net["protein"].dropna().unique())

    if True:  # toujours entrer (même sans arêtes, afficher les nœuds isolés)
        # Nombre de compétiteurs uniques (voisins distincts) par ABP.
        # Si le mode coopération est actif, les paires co-présentes dans un même PDB
        # ne sont PAS des compétiteurs -> exclues du décompte (donc nœud + petit).
        _abp_neighbors_c: dict = {}
        for (_a_c, _b_c) in _edge_wts:
            if _show_coop and frozenset((_a_c, _b_c)) in _coop_set:
                continue
            _abp_neighbors_c.setdefault(_a_c, set()).add(_b_c)
            _abp_neighbors_c.setdefault(_b_c, set()).add(_a_c)
        _abp_n_competitors: dict = {k: len(v)
                                    for k, v in _abp_neighbors_c.items()}

        _max_comp_c = max(_abp_n_competitors.values()
                          ) if _abp_n_competitors else 1

        # Cluster de séquence par ABP → couleur (deux directions)
        _df_s2lnk = _merged_net.dropna(subset=["_seqcl"]).copy()
        _prot_s2c = (
            _df_s2lnk.groupby("protein")["_seqcl"]
            .agg(lambda x: x.value_counts().index[0])
            .to_dict()
        )
        # Palette golden-ratio : chaque couleur saute ~222° → distance perceptuelle maximale
        import colorsys as _cs
        _s2c_all = sorted({v for v in _prot_s2c.values() if pd.notna(v)})
        _n_cls = max(len(_s2c_all), 1)

        def _golden_pal(n: int) -> list:
            _phi = 0.618033988749895
            _h, out = 0.05, []
            for i in range(n):
                _h = (_h + _phi) % 1.0
                _s = 0.88 if i % 2 == 0 else 0.62
                _v = 0.95 if i % 3 != 1 else 0.72
                r, g, b = _cs.hsv_to_rgb(_h, _s, _v)
                out.append("#{:02x}{:02x}{:02x}".format(
                    int(r * 255), int(g * 255), int(b * 255)))
            return out
        _s2c_color = dict(zip(_s2c_all, _golden_pal(_n_cls)))

        # ── Famille biologique par ABP (même définition que le barplot) ──────────
        import re as _re_fam
        _FAMILIES_NET = [
            ('Cofilin', r'cofilin|destrin'),
            ('Tropomodulin', r'tropomodulin|leiomodin'), ('Tropomyosin', r'tropomyosin'),
            ('Myosin', r'myosin'),
            ('CAP', r'cyclase-associated'), ('Cortactin', r'cortactin'),
            ('Capping (CapZ)', r'capping protein|f-actin-capping|capz'),
            ('Cross-linkers', r'spectrin|filamin|plastin|fimbrin|actinin|epididymis secretory protein li 37'),
            ('Utrophin/Dystrophin', r'utrophin|dystrophin'),
            ('Coronin', r'coronin'), ('Vinculin/Catenin/Talin',
                                      r'vinculin|catenin|talin'),
            ('Bacterial toxins', r'sipa|vop|exoy|tccc'), ('Profilin', r'profilin'),
            ('Gelsolin/Villin', r'gelsolin|villin|severin'), ('Formin',
                                                              r'formin|mdia|daam|diaphanous'),
            ('Troponin', r'troponin'), ('Adducin', r'adducin'),
            ('Arp2/3 (ARPC)', r'arpc|complex subunit'),
            ('Arp (Arp2/Arp3)', r'actin-related protein 2|actin-related protein 3|arp2|arp3'),
            ('Dematin', r'dematin'), ('AIP1/WD-repeat', r'wd repeat|aip1'),
            ('Kinase', r'kinase'), ('Afadin',
                                    r'afadin'), ('Thymosin/WH2', r'thymosin|wh2|ciboulot'),
        ]

        def _fam_net(t: str) -> str:
            t = str(t).lower()
            for _n, _pat in _FAMILIES_NET:
                if _re_fam.search(_pat, t):
                    return _n
            return 'Other'
        _fam_of = {_n: _fam_net(_n) for _n in _all_abp_net}
        _fams_present = sorted(set(_fam_of.values()) - {'Other'})
        # palette qualitative à fort contraste (couleurs maximalement distinctes)
        _DISTINCT = [
            "#e6194B", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#42d4f4",
            "#f032e6", "#bfef45", "#fabed4", "#469990", "#dcbeff", "#9A6324",
            "#800000", "#aaffc3", "#808000", "#000075", "#ffd8b1", "#000000",
            "#a9a9a9", "#ffe119", "#e6beff", "#fffac8", "#aed4e6", "#ff7f50",
        ]
        _fam_pal = {_f: _DISTINCT[_i % len(_DISTINCT)]
                    for _i, _f in enumerate(_fams_present)}
        _fam_pal['Autre'] = '#cccccc'

        # Noms courts affichés (le nom complet reste l'id du nœud -> survol/recherche)
        _DISP_NAMES = {
            'Actin-related protein 2': 'Arp2',
            'Actin-related protein 3': 'Arp3',
            'Actin-related protein 2/3 complex subunit 1': 'ARPC1',
            'Actin-related protein 2/3 complex subunit 1B': 'ARPC1B',
            'Actin-related protein 2/3 complex subunit 2': 'ARPC2',
            'Actin-related protein 2/3 complex subunit 3': 'ARPC3',
            'Actin-related protein 2/3 complex subunit 4': 'ARPC4',
            'Actin-related protein 2/3 complex subunit 5': 'ARPC5',
            'Actin-related protein 2/3 complex subunit 5-like protein': 'ARPC5L',
            'Coronin-1B,Methylated-DNA--protein-cysteine methyltransferase': 'Coronin-1B',
            'WD repeat-containing protein 1,Methylated-DNA--protein-cysteine methyltransferase': 'WDR1 (AIP1)',
            'Maltose/maltodextrin-binding periplasmic protein,Adenylate cyclase ExoY': 'ExoY',
            'Maltose/maltodextrin-binding periplasmic protein,TccC3': 'TccC3',
            'Cell invasion protein SipA': 'SipA',
            'Vibrio VopV': 'VopV (Vibrio)',
            'Adenylyl cyclase-associated protein 1': 'CAP1',
            'Epididymis secretory protein Li 37': 'L-plastin',
            'Alpha-catenin-like protein hmp-1': 'α-catenin (hmp-1)',
            'Catenin alpha-1': 'α-catenin-1',
            'Dematin actin binding protein': 'Dematin',
            'Src substrate cortactin': 'Cortactin',
            'Inositol-trisphosphate 3-kinase A': 'IP3 3-kinase A',
            'Protein diaphanous homolog 1': 'Diaphanous-1',
            'Isoform A2 of Troponin T, cardiac muscle': 'Troponin T (A2)',
            'Troponin T2, cardiac type': 'Troponin T2',
            'F-actin-capping protein subunit alpha-1': 'Capping α1',
            'F-actin-capping protein subunit beta': 'Capping β',
            'Isoform 2 of F-actin-capping protein subunit beta': 'Capping β (iso2)',
            'Myosin-14,Alpha-actinin A': 'Myosin-14',
            'Unconventional myosin-Ib': 'Myosin-Ib',
            'beta-cardiac myosin II': 'β-cardiac myosin',
            'Myosin heavy chain 4': 'MyHC-4',
            'Spectrin beta chain': 'Spectrin β',
            'Tropomyosin alpha-1 chain': 'Tropomyosin α1',
            'Tropomyosin beta chain': 'Tropomyosin β',
        }

        def _disp_name(_x):
            return _DISP_NAMES.get(_x, _x)

        def _wrap_lbl(name: str, max_ch: int = 18) -> str:
            words, lines, cur = name.split(), [], ""
            for w in words:
                if cur and len(cur) + 1 + len(w) > max_ch:
                    lines.append(cur)
                    cur = w
                else:
                    cur = (cur + " " + w).strip()
            if cur:
                lines.append(cur)
            return "\n".join(lines)

        # Physique forceAtlas2 : stabilisation auto → groupes visibles, pas de flottement
        from pyvis.network import Network
        # NB : pas de font_color= au constructeur -> il ECRASE le font par-nœud
        # (taille, gras multi:html, halo strokeWidth). Couleur par défaut via set_options.
        _net_c = Network(height="980px", width="100%", bgcolor="#ffffff")
        _net_c.set_options(
            '{"physics": {"enabled": false},'
            ' "nodes": {"shape": "dot", "font": {"size": 14, "color": "#111", "multi": "html"}},'
            ' "edges": {"smooth": false},'
            ' "interaction": {"hover": true, "tooltipDelay": 100, "hideEdgesOnDrag": false}}'
        )
        # Seuil : label visible en permanence seulement si ≥ 3 compétiteurs

        _pct_raw = {
            _nc: (_abp_c70_comp_count.get(_nc, 0) / _abp_c70_count[_nc])
            for _nc in _all_abp_net if _abp_c70_count.get(_nc, 0) > 0
        }
        _pct_min = min(_pct_raw.values()) if _pct_raw else 0.0
        _pct_max = max(_pct_raw.values()) if _pct_raw else 1.0

        def _pct_color(n_comp: int, n_total: int) -> str:
            """Jaune → rouge selon % de C70 en compétition (normalisé min/max réel)."""
            if n_total == 0:
                return "#cccccc"
            pct = n_comp / n_total
            t = (pct - _pct_min) / (_pct_max -
                                    _pct_min) if _pct_max > _pct_min else 1.0
            # jaune (#ffff00) → rouge (#cc0000)
            r = int(255 + (204 - 255) * t)
            g = int(255 * (1 - t))
            b = 0
            return f"#{r:02x}{g:02x}{b:02x}"

        # ── Layout spring regroupé par FAMILLE (arêtes virtuelles intra-famille) ──
        import networkx as _nxg
        _Glay = _nxg.Graph()
        _Glay.add_nodes_from(_all_abp_net)
        for (_ea, _eb), _ew in _edge_wts.items():
            if _ea in _all_abp_net and _eb in _all_abp_net:
                _Glay.add_edge(_ea, _eb, weight=min(_ew, 3))
        _cl_mem_lay: dict = {}
        for _nc in _all_abp_net:
            _cl = _fam_of.get(_nc)
            if _cl and _cl != 'Autre':
                _cl_mem_lay.setdefault(_cl, []).append(_nc)
        for _cl, _mem in _cl_mem_lay.items():
            for _i in range(len(_mem)):
                for _j in range(_i + 1, len(_mem)):
                    if _Glay.has_edge(_mem[_i], _mem[_j]):
                        _Glay[_mem[_i]][_mem[_j]]["weight"] += 20.0
                    else:
                        _Glay.add_edge(_mem[_i], _mem[_j], weight=20.0)
        # TOXINES épinglées au centre -> le réseau se construit autour
        _tox_lay = sorted({_n for _n in _all_abp_net
                           if _fam_of.get(_n) == 'Bacterial toxins'})
        _init_pos = {}
        for _ti, _tn in enumerate(_tox_lay):
            _a = 2 * 3.14159265 * _ti / max(len(_tox_lay), 1)
            _init_pos[_tn] = (0.13 * np.cos(_a), 0.13 * np.sin(_a))
        if _tox_lay:
            _layout_pos = _nxg.spring_layout(_Glay, weight="weight", k=0.62,
                                             iterations=600, seed=3,
                                             pos=_init_pos, fixed=_tox_lay)
        else:
            _layout_pos = _nxg.spring_layout(_Glay, weight="weight", k=0.62,
                                             iterations=600, seed=3)
        # Normalisation robuste : recentre + ramène les nœuds isolés (ex. Profilin
        # à seuil élevé) vers le cœur, pour que vis.js ne dézoome pas et que le
        # réseau reste lisible à TOUS les seuils (sans changer le look organique).
        _parr = np.array(list(_layout_pos.values()), dtype=float)
        _cxl, _cyl = _parr[:, 0].mean(), _parr[:, 1].mean()
        _dl = np.hypot(_parr[:, 0] - _cxl, _parr[:, 1] - _cyl)
        # rayon du cœur (ignore les isolés)
        _rref = float(np.percentile(_dl, 82)) or 1.0
        # rayon max autorisé (rapproche fort les isolés)
        _RMAX = _rref * 1.05

        def _normp(_p):
            _dx, _dy = _p[0] - _cxl, _p[1] - _cyl
            _r = (_dx ** 2 + _dy ** 2) ** 0.5
            if _r > _RMAX:
                _f = _RMAX / _r
                _dx, _dy = _dx * _f, _dy * _f
            return (_dx / _rref, _dy / _rref)
        _layout_pos = {_k: _normp(_v) for _k, _v in _layout_pos.items()}
        _SCALE = 600.0   # rayon du cœur ~1 après normalisation

        _label_thresh = 3
        for _nc in sorted(_all_abp_net):
            _n_c70_total = _abp_c70_count.get(_nc, 0)
            _n_c70_comp = _abp_c70_comp_count.get(_nc, 0)
            _n_comp_nc = _abp_n_competitors.get(_nc, 0)
            _pct_val = round(_n_c70_comp / _n_c70_total *
                             100) if _n_c70_total else 0
            # TAILLE = proportion de sites de liaison en conflit (% C70 contestés)
            _sz = 4 + 0.16 * _pct_val
            # COULEUR = famille biologique de l'ABP
            _fam_nc = _fam_of.get(_nc, 'Autre')
            _col = _fam_pal.get(_fam_nc, '#bbbbbb')
            _show_lbl = _n_comp_nc >= _label_thresh
            _c70_line = (f"Competing C70: {_n_c70_comp} / {_n_c70_total} ({_pct_val}%)"
                         if _n_c70_comp else f"Total C70: {_n_c70_total} (0%)")
            _pos_set = _abp_fil_pos_net.get(_nc, set())
            _pos_str = ", ".join(
                sorted(_pos_set, key=lambda x: {
                       "-": 0, "side": 1, "+": 2}.get(x, 99))
            ) if _pos_set else "?"
            # extrémité du filament : "+" / "−" seulement si l'ABP n'est PAS aussi à
            # l'autre extrémité ; rien si les deux (±), latéral, ou inconnu.
            _has_p, _has_m = ("+" in _pos_set), ("-" in _pos_set)
            _sign = ("+" if (_has_p and not _has_m)
                     else "−" if (_has_m and not _has_p) else "")
            # contour neutre (plus de barbé/pointé/latéral)
            _border_col = "#555555"
            _xy = _layout_pos.get(_nc, (0.0, 0.0))
            _net_c.add_node(
                _nc,
                label=("<b>" + _wrap_lbl(_disp_name(_nc)).replace("\n", "<br>")
                       + (f" ({_sign})" if _sign else "") + "</b>"
                       if _show_lbl else ""),
                title=(f"{_nc}\n"
                       f"Famille : {_fam_nc}\n"
                       f"Position filament : {_pos_str}\n"
                       f"{_c70_line}\n"
                       f"Partenaires ABP uniques en concurrence : {_n_comp_nc}"),
                size=_sz,
                color={"background": _col, "border": _border_col},
                borderWidth=1,
                font={"size": 26, "face": "Arial",
                      "color": "#111111", "multi": "html"},
                x=float(_xy[0]) * _SCALE, y=float(_xy[1]) * _SCALE,
                physics=False,
            )
            # (Le signe d'extrémité +/− est désormais mis dans le NOM entre
            # parenthèses, plus de nœud-texte superposé au centre.)
        _filtered_edges = {k: v for k,
                           v in _edge_wts.items() if v >= _min_edge_w}
        _n_coop_shown = 0
        # nœuds de la famille "Toxines bactériennes" -> leurs arêtes en NOIR (fines)
        _tox_nodes = {_n for _n in _all_abp_net
                      if _fam_of.get(_n) == 'Bacterial toxins'}
        for (_a_c, _b_c), _w_c in sorted(_filtered_edges.items(),
                                         key=lambda x: x[1], reverse=True):
            _is_coop = _show_coop and frozenset((_a_c, _b_c)) in _coop_set
            _is_tox = (_a_c in _tox_nodes) or (_b_c in _tox_nodes)
            if _is_coop:
                _n_coop_shown += 1
                _npdb = _coop_pairs[frozenset((_a_c, _b_c))]
                _net_c.add_edge(
                    _a_c, _b_c, width=2,
                    title=f"Cooperate — co-present in {_npdb} PDBs",
                    color={"color": "#1f9e3a", "opacity": 0.65},
                )
            elif _is_tox:
                _net_c.add_edge(
                    _a_c, _b_c, width=0.8,
                    title=f"Toxin — {_w_c} shared C70 pair(s)",
                    color={"color": "#000000", "opacity": 0.7},
                )
            else:
                _net_c.add_edge(
                    _a_c, _b_c, width=1.8,
                    title=f"{_w_c} shared C70 pair(s)",
                    color={"color": "#000000", "opacity": 0.25},
                )

        if _show_coop and _n_coop_shown:
            pass
        _net_html = _net_c.generate_html()
        # Désactiver la physique dès que la stabilisation est terminée
        for _pat_c in [
            "network = new vis.Network(container, data, options);",
            "var network = new vis.Network(container, data, options);",
        ]:
            if _pat_c in _net_html:
                _net_html = _net_html.replace(
                    _pat_c,
                    _pat_c +
                    "\n  network.once('stabilizationIterationsDone', function(){"
                    "\n    network.setOptions({physics:{enabled:false}});"
                    "\n  });",
                )
                break
        # Barre de recherche flottante à l'intérieur du graphe
        import json as _json_abp
        _search_map_abp = _json_abp.dumps(
            {nc: nc for nc in _all_abp_net}, ensure_ascii=False)
        _inject_search_abp = (
            "(function(){\n"
            f"  var SEARCH_MAP = {_search_map_abp};\n"
            "  var cont = document.getElementById('mynetwork').parentElement;\n"
            "  cont.style.position = 'relative';\n"
            "  var box = document.createElement('input');\n"
            "  box.type = 'text';\n"
            "  box.placeholder = 'Search an ABP…';\n"
            "  box.style.cssText = 'position:absolute;top:10px;right:10px;z-index:9999;'\n"
            "    + 'padding:6px 12px;border:2px solid #aaa;border-radius:20px;'\n"
            "    + 'font-size:13px;width:240px;outline:none;background:white;';\n"
            "  cont.appendChild(box);\n"
            "  var origColors = {};\n"
            "  network.body.data.nodes.get().forEach(function(n){ origColors[n.id] = n.color; });\n"
            "  function _matches(n, q){\n"
            "    if(!q) return false;\n"
            "    var extra = SEARCH_MAP[String(n.id)] || '';\n"
            "    return (extra + ' ' + String(n.label||'')).toLowerCase().indexOf(q) !== -1;\n"
            "  }\n"
            "  box.addEventListener('input', function(){\n"
            "    var q = this.value.trim().toLowerCase();\n"
            "    var all = network.body.data.nodes.get();\n"
            "    var updates = all.map(function(n){\n"
            "      if(!q) return {id:n.id, color:origColors[n.id], font:{color:'#111'}, borderWidth:2};\n"
            "      var m = _matches(n, q);\n"
            "      return {id:n.id,\n"
            "        color: m ? {background:(origColors[n.id]&&origColors[n.id].background)||origColors[n.id], border:'#ff2222'}\n"
            "                 : {background:'#e4e4e4', border:'#cccccc'},\n"
            "        font:{color: m ? '#111' : '#ccc'}, borderWidth: m ? 5 : 1};\n"
            "    });\n"
            "    network.body.data.nodes.update(updates);\n"
            "    if(q){ var hit=all.filter(function(n){return _matches(n,q);});\n"
            "      if(hit.length===1) network.focus(hit[0].id,{scale:1.8,animation:{duration:400}}); }\n"
            "  });\n"
            "})();\n"
        )
        _net_html = _net_html.replace(
            "</body>", "<script>" + _inject_search_abp + "</script>\n</body>")
        with st.expander("Show the network graph", expanded=False):
            st.components.v1.html(_net_html, height=880, scrolling=False)
        # Légende COULEUR = famille (uniquement les familles présentes dans le réseau)
        _fams_in_net = sorted({_fam_of[_n] for _n in _all_abp_net})
        _fam_items = "".join(
            "<span style='display:flex;align-items:center;gap:5px;'>"
            f"<span style='width:13px;height:13px;border-radius:50%;background:{_fam_pal.get(_f, '#bbb')};"
            "display:inline-block;border:1px solid #999'></span> " + _f + "</span>"
            for _f in _fams_in_net
        )
        _legend_fam = (
            "<div style='display:flex;flex-wrap:wrap;gap:10px 14px;align-items:center;"
            "font-size:11px;margin-top:4px;'><b>Couleur = famille :</b>" + _fam_items + "</div>"
        )
        st.markdown(_legend_fam, unsafe_allow_html=True)

    st.divider()

# ── Réseau de coopération ABP ─────────────────────────────────────────────────
if (_net_view == "Cooperation" and os.path.exists(proteins_path)
        and os.path.exists(_all_data_path)
        and len(abp_global) > 0 and "Binding site S1" in abp_global.columns
        and _coop_pairs):
    st.subheader("ABP cooperation network")
    import networkx as _nx_coop
    from pyvis.network import Network as _Net_coop
    _Gco = _nx_coop.Graph()
    for _pr_co, _np_co in _coop_pairs.items():
        _a_n, _b_n = tuple(_pr_co)
        _Gco.add_edge(_a_n, _b_n, weight=_np_co)
    _coop_nodes = list(_Gco.nodes())
    _deg_co = dict(_Gco.degree())
    _maxdeg_co = max(_deg_co.values()) if _deg_co else 1
    # Disposition EN TUILES : chaque module (composant connexe) dans sa propre case,
    # étalé et normalisé -> pas de chevauchement (le spring global tassait tout).
    import math as _math_co
    _comps_co = sorted(_nx_coop.connected_components(_Gco),
                       key=len, reverse=True)
    _node_compsize = {}
    for _comp in _comps_co:
        for _n in _comp:
            _node_compsize[_n] = len(_comp)
    _cols_co = max(1, int(_math_co.ceil(_math_co.sqrt(len(_comps_co)))))
    _CELL_CO = 2.8
    _pos_co = {}
    for _i_co, _comp in enumerate(_comps_co):
        _sub = _Gco.subgraph(_comp)
        if len(_comp) == 1:
            _sp = {next(iter(_comp)): (0.0, 0.0)}
        else:
            _sp = _nx_coop.spring_layout(_sub, k=1.3, iterations=400, seed=3)
        _xs = [p[0] for p in _sp.values()]
        _ys = [p[1] for p in _sp.values()]
        _cx = (max(_xs) + min(_xs)) / 2
        _cy = (max(_ys) + min(_ys)) / 2
        _rx = (max(_xs) - min(_xs)) or 1.0
        _ry = (max(_ys) - min(_ys)) or 1.0
        _rad = 0.45 + 0.12 * len(_comp)
        _gx = (_i_co % _cols_co) * _CELL_CO
        _gy = -(_i_co // _cols_co) * _CELL_CO
        for _n, (x, y) in _sp.items():
            _pos_co[_n] = (_gx + (x - _cx) / _rx * _rad,
                           _gy + (y - _cy) / _ry * _rad)
    _SCALE_CO = 760.0
    _net_co = _Net_coop(height="820px", width="100%", bgcolor="#ffffff",
                        font_color="#111")
    _net_co.set_options('{"physics": {"enabled": false}}')
    for _nc_co in _coop_nodes:
        _dg = _deg_co.get(_nc_co, 0)
        _sz_co = 6 + 20 * (_dg / _maxdeg_co) ** 0.7
        _pset_co = _abp_fil_pos_net.get(_nc_co, set())
        _bcol_co = _pos_border_color(_pset_co)
        _xy_co = _pos_co.get(_nc_co, (0.0, 0.0))
        # étiquette si nœud connecté (>=2 partenaires) ou petit module (<=4 membres)
        _show_lbl_co = (_dg >= 2) or (_node_compsize.get(_nc_co, 1) <= 4)
        _net_co.add_node(
            _nc_co, label=_wrap_lbl(_nc_co) if _show_lbl_co else "",
            title=f"{_nc_co}\nCooperating partners: {_dg}",
            size=_sz_co, color={"background": "#7ec98f", "border": _bcol_co},
            borderWidth=2, font={"size": 4, "face": "Arial", "color": "#111"},
            x=float(_xy_co[0]) * _SCALE_CO, y=float(_xy_co[1]) * _SCALE_CO,
            physics=False,
        )
    for _a_e, _b_e, _d_e in _Gco.edges(data=True):
        _w_e = _d_e.get("weight", 1)
        _net_co.add_edge(
            _a_e, _b_e, width=1 + 1.5 * _w_e,
            title=f"co-present in {_w_e} PDBs",
            color={"color": "#1f9e3a", "opacity": 0.5},
        )
    with st.expander("Show the network graph", expanded=False):
        st.components.v1.html(_net_co.generate_html(), height=780, scrolling=False)

# Mode FAST : on s'arrête après les réseaux (saute Détail/Filament PyMOL/MSA)
if _FAST:
    st.stop()

# ── Détail par ABP ─────────────────────────────────────────────────────────────
if (os.path.exists(proteins_path) and os.path.exists(_all_data_path)
        and len(abp_global) > 0 and "Binding site S1" in abp_global.columns):
    st.subheader("Per-ABP detail")

    _NO_ABP_LABEL = "— PDBs without ABP —"
    abp_names = [_NO_ABP_LABEL] + abp_global["Protein"].tolist()

    # Navigation depuis le graphe : un clic sur un nœud vert (ABP) clique (en JS,
    # via window.top) le bouton caché du même nom → sélectionne cet ABP ici.
    st.markdown(
        "<div id='abpsel'></div>"
        "<style>[class*='st-key-__ab_']{position:absolute!important;"
        "left:-9999px!important;top:auto!important;margin:0!important}</style>",
        unsafe_allow_html=True,
    )
    for _i, _abpn in enumerate(abp_global["Protein"].tolist()):
        if st.button(_abpn, key=f"__ab_{_i}"):
            st.session_state["sel_abp_detail"] = _abpn

    sel_abp = st.selectbox("Select an ABP",
                           abp_names, key="sel_abp_detail")
    _is_no_abp = (sel_abp == _NO_ABP_LABEL)

    # PDB IDs selon la sélection
    if _is_no_abp:
        _pdbs_with_abp_s = set(df_pp_all[~df_pp_all["is_actin"]]["pdb_id"])
        _all_pdbs_s = set(df_all_g["subunit_1"].str.split("_").str[0])
        abp_pdbs = _all_pdbs_s - _pdbs_with_abp_s
    else:
        abp_pdbs = set(merged[merged["protein"] == sel_abp]["pdb_id"])

    if not _is_no_abp:
        # ── Clusters actin–ABP (interactions hétéro) ─────────────────────────
        st.markdown("#### Clusters d'interaction actin–ABP")
        _abp_subunits = set(merged[merged["protein"] == sel_abp]["_abp_chain"])
        _df_hetero = df_all_g.copy()
        _df_hetero["_pdb_h"] = _df_hetero["subunit_1"].str.split("_").str[0]
        abp_int = _df_hetero[
            _df_hetero["cluster_data_70"].notna() & (
                (
                    _df_hetero["s1_actine"].fillna(False).astype(bool) &
                    ~_df_hetero["s2_actine"].fillna(False).astype(bool) &
                    _df_hetero["subunit_2"].str.lower().isin(_abp_subunits)
                ) | (
                    ~_df_hetero["s1_actine"].fillna(False).astype(bool) &
                    _df_hetero["s2_actine"].fillna(False).astype(bool) &
                    _df_hetero["subunit_1"].str.lower().isin(_abp_subunits)
                )
            )
        ].copy()

        if abp_int.empty:
            st.info("No hetero interaction found for this ABP.")
        else:
            def _most_freq_pair(grp):
                """Paire (subunit_1, subunit_2) la plus fréquente dans le groupe."""
                pairs = grp.groupby(["subunit_1", "subunit_2"]).size()
                if pairs.empty:
                    return pd.Series({"rep_actin": "", "rep_abp": ""})
                best = pairs.idxmax()
                return pd.Series({"rep_actin": best[0], "rep_abp": best[1]})

            _grp_cols = [
                "cluster_data_70",
                "s1_binding_site_cluster_data_70",
                "s2_binding_site_cluster_data_70",
            ]
            _stats = (
                abp_int.groupby(_grp_cols, dropna=False)
                .agg(nb_inter=("subunit_1", "count"), nb_pdb=("_pdb_h", "nunique"))
                .reset_index()
            )
            _pairs = (
                abp_int.groupby(_grp_cols, dropna=False)
                .apply(_most_freq_pair)
                .reset_index()
            )
            abp_clust = _stats.merge(_pairs, on=_grp_cols).sort_values(
                "nb_inter", ascending=False
            )
            _total_abp_inter = max(len(abp_int), 1)
            abp_clust["% PDB"] = (
                abp_clust["nb_pdb"] / max(len(abp_pdbs), 1) * 100
            ).round(1)
            abp_clust["% ABP-actin interactions"] = (
                abp_clust["nb_inter"] / _total_abp_inter * 100
            ).round(1)
            abp_clust["Binding site"] = (
                abp_clust["s1_binding_site_cluster_data_70"].astype(str)
                + " × "
                + abp_clust["s2_binding_site_cluster_data_70"].astype(str)
            )
            st.dataframe(
                abp_clust.rename(columns={
                    "cluster_data_70": "Cluster C70",
                    "nb_inter": "# interactions",
                    "nb_pdb": "# PDBs",
                })[["Cluster C70", "Binding site",
                    "# interactions", "% ABP-actin interactions", "# PDBs", "% PDB"]],
                hide_index=True,
                use_container_width=True,
            )

            # ── Viewer 3D inline ──────────────────────────────────────────────
            _bfac_dir_abp = _Path(
                "data/filtered/details/structures_files/bfactor_c70_interface")
            _valid_rows = [
                row for _, row in abp_clust.iterrows()
                if (_bfac_dir_abp / f"{row['cluster_data_70']}.pdb").exists()
            ]
            if not _valid_rows:
                pass
            else:
                _ylord3d = ["#FFFFCC", "#FFF0A9", "#FEE186", "#FECA65", "#FDAA48",
                            "#FC8C3B", "#FC5A2D", "#EC2D21", "#D30F20", "#800026"]
                _grn3d = ["#FFFFCC", "#D9F0A3", "#ADDD8E", "#78C679",
                          "#41AB5D", "#238443", "#006837"]
                for _i3d in range(0, len(_valid_rows), 2):
                    _pair = _valid_rows[_i3d:_i3d + 2]
                    _cols3d = st.columns(len(_pair))
                    for _ci3d, _crow3d in enumerate(_pair):
                        _c70_3d = str(_crow3d["cluster_data_70"])
                        _pdb_3d_txt = (_bfac_dir_abp /
                                       f"{_c70_3d}.pdb").read_text()
                        _bmax_a3d, _bmax_b3d = 1.0, 1.0
                        for _ln in _pdb_3d_txt.splitlines():
                            if _ln.startswith("ATOM") and len(_ln) > 66:
                                try:
                                    _bv3 = float(_ln[60:66].strip())
                                    if _ln[21] == "A":
                                        _bmax_a3d = max(_bmax_a3d, _bv3)
                                    elif _ln[21] == "B":
                                        _bmax_b3d = max(_bmax_b3d, _bv3)
                                except ValueError:
                                    pass
                        import py3Dmol
                        _v3d = py3Dmol.view(width="100%", height=380)
                        _v3d.addModel(_pdb_3d_txt, "pdb")
                        _v3d.setStyle({}, {})
                        _v3d.addSurface(py3Dmol.SES,
                                        {"opacity": 1, "colorscheme": {
                                            "prop": "b", "gradient": "linear",
                                            "colors": _ylord3d, "min": 0, "max": _bmax_a3d,
                                        }}, {"chain": "A"})
                        _v3d.addSurface(py3Dmol.SES,
                                        {"opacity": 1, "colorscheme": {
                                            "prop": "b", "gradient": "linear",
                                            "colors": _grn3d, "min": 0, "max": _bmax_b3d,
                                        }}, {"chain": "B"})
                        _v3d.setBackgroundColor("white")
                        _v3d.zoomTo()
                        with _cols3d[_ci3d]:
                            st.markdown(f"**{_c70_3d}**")
                            st.components.v1.html(
                                _v3d._make_html(), height=390, scrolling=False)

        st.divider()

    # ── Interactions actin-actin (homo) ─────────────────────────────────────
    _homo_label = ("PDBs without ABP" if _is_no_abp
                   else f"the same PDBs as {sel_abp[:30]}")
    st.markdown(f"#### Interactions actin-actin (homo) — {_homo_label}")
    _df_homo = df_all_g.copy()
    _df_homo["_pdb"] = _df_homo["subunit_1"].str.split("_").str[0]
    homo_cooc = _df_homo[
        _df_homo["s1_actine"].fillna(False).astype(bool) &
        _df_homo["s2_actine"].fillna(False).astype(bool) &
        _df_homo["_pdb"].isin(abp_pdbs) &
        _df_homo["cluster_data_70"].notna()
    ][["_pdb", "cluster_data_70",
       "s1_binding_site_cluster_data_70", "s2_binding_site_cluster_data_70"]].copy()

    # Fusionner S1 et S2 en une paire canonical (triée, S1/S2 non distingués)
    homo_cooc["Binding sites"] = homo_cooc.apply(
        lambda r: " × ".join(sorted([
            str(r["s1_binding_site_cluster_data_70"]),
            str(r["s2_binding_site_cluster_data_70"])
        ])), axis=1
    )

    total_homo = len(homo_cooc)
    if homo_cooc.empty:
        st.info("No homo actin-actin interaction in these PDBs.")
    else:
        homo_summary = (
            homo_cooc.groupby(
                ["cluster_data_70", "Binding sites"], dropna=False)
            .agg(nb_pdb=("_pdb", "nunique"), nb_inter=("_pdb", "count"))
            .reset_index()
            .sort_values("nb_pdb", ascending=False)
        )
        homo_summary["% PDB"] = (
            homo_summary["nb_pdb"] / len(abp_pdbs) * 100
        ).round(1)
        homo_summary["% homo interactions"] = (
            homo_summary["nb_inter"] / total_homo * 100
        ).round(1)
        homo_summary = homo_summary.rename(columns={
            "cluster_data_70": "Cluster C70",
            "nb_pdb": "# PDBs",
            "nb_inter": "# homo interactions",
        })[["Cluster C70", "Binding sites",
            "# PDBs", "% PDB", "# homo interactions", "% homo interactions"]]
        st.dataframe(homo_summary, hide_index=True, use_container_width=True)

    # ── Découverte PDB : motif d'interface cherché dans tout le PDB ───────────
    if not _is_no_abp:
        folddisco_view.render_discovery(sel_abp)

    # ── ProteoCast de l'ABP : paysage mutationnel + empreinte + 3D ────────────
    # (Le bouton de calcul global est remonté sous la section téléchargement.)
    st.markdown("#### ABP ProteoCast — mutational landscape + 3D structure")
    if _is_no_abp:
        pass
    else:
        _pc_row = None
        if _pc_status is not None:
            _m = _pc_status[_pc_status["abp_title"] == sel_abp]
            if len(_m):
                _pc_row = _m.iloc[0]
        _pc_slug = (_pc_row["slug"] if _pc_row is not None
                    else _re_sd.sub(r"[^a-zA-Z0-9]+", "_", str(sel_abp)).strip("_")[:60])
        _pc_uni = _pc_row["uniprot"] if _pc_row is not None else None
        # (« Non calculé / échec » + bouton Compute/Retry sont gérés dans l'onglet
        #  « Paysage mutationnel » par _render_abp_proteocast, au plus près du message.)
        _pc_tab1, _pc_tab2 = st.tabs(
            ["Paysage mutationnel", "Visualisation 3D"])
        with _pc_tab1:
            _render_abp_proteocast(sel_abp, _abp_subunits)
            st.markdown(
                "**Actin side — conservation of the residues targeted by this ABP**")
            _render_abp_actin_conservation(sel_abp)
            _imgs = proteocast_view.result_images(_pc_slug)
            if _imgs:
                with st.expander("ProteoCast figures (results folder)"):
                    for _im in _imgs:
                        st.image(_im, use_container_width=True)
            _pc_zip = proteocast_view.folder_zip(_pc_slug)
            if _pc_zip is not None:
                _pc_zip_name, _pc_zip_bytes = _pc_zip
                st.download_button(
                    "Download the full ProteoCast folder (.zip)",
                    data=_pc_zip_bytes, file_name=_pc_zip_name,
                    mime="application/zip", key=f"pc_dl_folder_{_pc_slug}",
                    help="Everything proteocast.ijm.fr returned for this ABP: "
                         "MSA, structures, images and sub-folders.")
        with _pc_tab2:
            _pc_style = st.radio("Display", ["Surface", "Cartoon"],
                                 horizontal=True, key=f"pc_view_{_pc_slug}")
            # gros ABP : restreindre la 3D à la zone qui touche l'actin (ABD)
            _pc_focus = None
            _pc_fz = _abp_actin_focus(sel_abp, _pc_slug)
            if _pc_fz:
                _pflo, _pfhi, _pql = _pc_fz
                if st.toggle(
                        f"Zoom on the actin-binding region (~{_pflo}–{_pfhi}) — "
                        f"otherwise the whole protein ({_pql} aa)",
                        value=False, key=f"pc3d_whole_{_pc_slug}",
                        help="Zoom on the domain that contacts actin; "
                             "off = whole protein."):
                    _pc_focus = (_pflo, _pfhi)
            _pc_struct = proteocast_view.structure_for(
                _pc_slug, _pc_uni, _pc_mt)
            if _pc_struct is None:
                st.info("Structure unavailable (AlphaFold not found for this "
                        "UniProt, or no result provided).")
            else:
                import py3Dmol
                _pc_pdb = (_crop_pdb_residues(_pc_struct["pdb"], *_pc_focus)
                           if _pc_focus else _pc_struct["pdb"])
                _grad = (["#2166ac", "#f7f7f7", "#b2182b"] if _pc_struct["metric"]
                         else ["#ff7f00", "#ffff99", "#1f78b4"])
                _v = py3Dmol.view(width="100%", height=500)
                _v.addModel(_pc_pdb, "pdb")
                _scheme = {"prop": "b", "gradient": "linear", "colors": _grad,
                           "min": _pc_struct["vmin"], "max": _pc_struct["vmax"]}
                if _pc_style == "Surface":
                    _v.setStyle({}, {})
                    _v.addSurface(
                        py3Dmol.SES, {"opacity": 1.0, "colorscheme": _scheme})
                else:
                    _v.setStyle({}, {"cartoon": {"colorscheme": _scheme}})
                _v.zoomTo()
                _v.setBackgroundColor("white")
                st.components.v1.html(
                    _v._make_html(), height=510, scrolling=False)

# ── Explorateur interactif : séquence query / cluster / ABP / paire d'ABP ─────
if not _FAST:
    st.header("Interactive explorer", anchor="explorateur")
    _pp = build_passport(pp_mtimes())
    if _pp is None:
        st.info("Residue-passport data unavailable "
                "(run the pipeline to generate them).")
    else:
        residue_explorer.render_explorer(_pp)

# ── MSA — Protéines d'interface ───────────────────────────────────────────────

st.header("MSA — Interface proteins", anchor="msa-proteines")


if not _Path("data/filtered/filtered_all_data.csv").exists():
    st.warning(
        "filtered_all_data.csv not found — run the filtering steps first.")
else:
    _MSA_ALN_DIR.mkdir(parents=True, exist_ok=True)


    _msa_section_full(
        "Myosins", "myosin",
        lambda t: "myosin" in t.lower() and "tropomyosin" not in t.lower(),
        note="All myosin structures (rigor + ADP), including beta-cardiac.",
    )
    _msa_section_full(
        "Tropomyosins", "tropomyosin",
        lambda t: "tropomyosin" in t.lower(),
    )
    _msa_section_full(
        "Plastin and related", "plastin",
        lambda t: any(x in t.lower()
                      for x in ["plastin", "spectrin beta", "filamin", "utrophin"]),
        note="Plastin-3, Spectrin beta chain, Filamin-A, Utrophin.",
    )
    _msa_section_s2_clusters()
