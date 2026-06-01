# streamlit run script/streamlit.py
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
import matplotlib.pyplot as plt

from network_viz import (
    _bip_mtimes, _load_bipartite_base, _build_bipartite_html,
    _load_res4, _build_bipartite_c70_html, _build_s1_3d_html,
    _build_tripartite_graph_html, _build_global_graph_html,
    _global_bfac_max, _BIP_CACHE_VERSION, _BIPARTITE_FILES,
    _AA_RESTYPE_HEX,
)
from data_analysis import (
    _build_abp_heatmap_data, _build_c70_jaccard_edges,
    _build_s1_superclusters, _S1_SUPER_FILES, _ABP_HM_FILES,
)
from msa_analysis import (
    norm_chain_id, _cluster_protein_index,
    _msa_section_full, _msa_section_s2_clusters, _msa_section_s1_clusters,
    _MSA_ALN_DIR, _MSA_RIGOR_PDBS,
    _s1_get_ch_maps, _msa_contact_analysis,
)

_GRAPH_COMPONENT = st.components.v1.declare_component(
    "actin_graph",
    path=str(_Path(__file__).parent / "graph_component"),
)

st.set_page_config(layout="wide", page_title="Analyse actine-ABP - PPI3D")

st.title("Analyse des interactions actine-actine et actine-ABP - PPI3D")

with st.sidebar:
    st.markdown("## Sommaire")
    st.markdown("""
- [Téléchargement des données](#telechargement-des-donnees)
- [Données filtrées](#donnees-filtrees-s1-actine)
- [Structures PDB valides](#structures-pdb-valides)
- [Clusters d'interactions](#clusters-d-interactions)
- [ABP](#abp)
- [MSA — Protéines d'interface](#msa-proteines)
""")

# ---------------------------------------------------------------------------
# Section téléchargement
# ---------------------------------------------------------------------------

st.header("Téléchargement des données")

STEPS = {
    "1/15":  "Téléchargement du summary PPI3D (BLAST)",
    "2/15":  "Téléchargement des entrées PDB",
    "3/15":  "Téléchargement de toutes les données (cluster table)",
    "4/15":  "Filtrage des structures (≥ 4 actines connectées) - notebook",
    "5/15":  "Téléchargement des interactions d'interface",
    "6/15":  "Alignement MAFFT par cluster de séquences",
    "7/15":  "Analyse des clusters d'interaction C70 - notebook",
    "8/15":  "Calcul B-factors interface C70 par cluster",
    "9/15":  "Génération script PyMOL surface complète C70",
    "10/15": "Génération scripts PyMOL par site S1",
    "11/15": "Analyse interface par cluster C70 - notebook",
    "12/15": "Heatmap S1 binding site et références clusters - notebook",
    "13/15": "Calcul B-factors S1 par cluster",
    "14/15": "Analyse ABP — compétition et interfaces - notebook",
    "15/15": "Sessions PyMOL filament par ABP",
}

# Fichier de sortie attendu pour chaque étape
STEP_OUTPUT_FILES = {
    "1/15":  "data/raw/ppi3d_actin_summary.csv",
    "2/15":  "data/raw/pdb_entry_results.csv",
    "3/15":  "data/raw/all_data.csv",
    "4/15":  "data/filtered/filtered_pdb_entry.csv",
    "5/15":  "data/filtered/details/1.interactions.csv",
    "6/15":  "data/alignments/.done",
    "7/15":  "data/filtered/patches_infos_cluster_data_70.csv",
    "8/15":  "data/filtered/details/structures_files/bfactor_c70_interface",
    "9/15":  "data/filtered/details/structures_files/bfactor_c70_interface/view_full_surface.pml",
    "10/15": "data/filtered/details/structures_files/bfactor_c70_interface/by_s1_cluster",
    "11/15": "data/visualisations/actin_c70_contacts",
    "12/15": "data/visualisations/actin_s1_all_equitable_heatmap.png",
    "13/15": "data/filtered/details/structures_files/bfactor_cluster",
    "14/15": "data/visualisations/abp_analysis_done.flag",
    "15/15": "data/filtered/details/structures_files/filament/by_abp",
}

STEP_KEYS = list(STEPS.keys())
TOTAL = len(STEPS)

SKIP_KEYWORDS = ["Nothing to do", "unchanged",
                 "Using existing", "already up", "Déjà à jour",
                 "No changes detected"]


def step_md(key, label, state):
    if state == "pending":
        return f"⬜ **{key}** — {label}"
    if state == "running":
        return f"🔄 **{key}** — {label} *(en cours...)*"
    if state == "done":
        return f"✅ **{key}** — {label} — *données téléchargées*"
    if state == "skipped":
        return f"✅ **{key}** — {label} — *déjà à jour*"
    if state == "error":
        return f"❌ **{key}** — {label} — *erreur*"
    return f"⬜ **{key}** — {label}"


def initial_state(key):
    path = STEP_OUTPUT_FILES.get(key)
    if path and os.path.exists(path):
        return "skipped"
    return "pending"


def parse_sub_progress(line):
    """
    Retourne (current, total, label) si la ligne contient une progression,
    sinon None.
    Patterns détectés :
      - Étape 2 : "1/148 4A7F"          → PDB X/total
      - Étape 5 : "Downloading 1/3"     → pass X/total
      - Étape 6 : "200,000 lignes traitées (471 Mo en RAM)..."
    """
    # Étape 2 : "X/Y PDBID"
    m = re.match(r"^(\d+)/(\d+)\s+\w+", line.strip())
    if m:
        return int(m.group(1)), int(m.group(2)), f"PDB {line.strip()}"

    # Étape 5 : "Downloading X/Y"
    m = re.search(r"Downloading\s+(\d+)/(\d+)", line)
    if m:
        return int(m.group(1)), int(m.group(2)), line.strip()

    return None


# Bouton en haut, étapes en dessous

clicked = st.button("Lancer le téléchargement")

progress_bar = st.empty()
log_placeholder = st.empty()

placeholders = {}
for key, label in STEPS.items():
    placeholders[key] = st.empty()
    placeholders[key].markdown(step_md(key, label, initial_state(key)))

sub_progress_bar = st.empty()
sub_progress_text = st.empty()

if "just_downloaded" not in st.session_state:
    st.session_state.just_downloaded = False

all_done = all(
    initial_state(k) == "skipped"
    for k in STEP_KEYS
    if STEP_OUTPUT_FILES.get(k)  # ignorer les étapes sans fichier indicateur
)
if all_done and not st.session_state.just_downloaded:
    st.write("✅ Données déjà téléchargées")


if clicked:
    st.session_state.just_downloaded = True

    for key, label in STEPS.items():
        init = initial_state(key)
        placeholders[key].markdown(
            step_md(key, label, init if init == "skipped" else "pending"))

    pb = progress_bar.progress(0, text="Démarrage...")

    current_idx = -1
    step_had_skip = {}
    log_lines = []
    step_log_lines = []

    cmd = [sys.executable, "-u", "-m", "script.data_extract.pipeline_data"]
    if sys.platform == "darwin":
        cmd = ["caffeinate", "-i"] + cmd
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    for line in proc.stdout:
        line = line.rstrip()

        # Détecter le début d'une nouvelle étape
        for i, key in enumerate(STEP_KEYS):
            if f"ETAPE : {key}" in line:
                if current_idx >= 0:
                    prev_key = STEP_KEYS[current_idx]
                    state = "skipped" if step_had_skip.get(
                        prev_key) else "done"
                    placeholders[prev_key].markdown(
                        step_md(prev_key, STEPS[prev_key], state))

                current_idx = i
                step_log_lines = []
                placeholders[key].markdown(step_md(key, STEPS[key], "running"))
                pb.progress(i / TOTAL, text=f"Étape {key} — {STEPS[key]}")
                sub_progress_bar.empty()
                sub_progress_text.empty()
                log_placeholder.empty()
                break

        # Sous-progression de l'étape courante
        result = parse_sub_progress(line)
        if result is not None:
            current, total, label = result
            if current is not None and total is not None and total > 0:
                sub_progress_bar.progress(current / total)
                sub_progress_text.caption(label)
            else:
                sub_progress_text.caption(label)

        # Suivre les mots-clés "skip"
        if current_idx >= 0:
            current_key = STEP_KEYS[current_idx]
            if any(kw in line for kw in SKIP_KEYWORDS):
                step_had_skip[current_key] = True

        if line:
            log_lines.append(line)
            step_log_lines.append(line)
            log_placeholder.text("\n".join(step_log_lines[-15:]))

    proc.wait()

    sub_progress_bar.empty()
    sub_progress_text.empty()
    log_placeholder.empty()

    if current_idx >= 0:
        last_key = STEP_KEYS[current_idx]
        state = "skipped" if step_had_skip.get(last_key) else "done"
        if proc.returncode != 0:
            state = "error"
        placeholders[last_key].markdown(
            step_md(last_key, STEPS[last_key], state))

    if proc.returncode == 0:
        pb.progress(1.0, text="Téléchargement terminé ✓")
    else:
        pb.progress(1.0, text="Erreur — pipeline arrêté")
        st.error("Pipeline arrêté — dernières lignes de log :")
        st.code("\n".join(log_lines[-40:]), language=None)


# ---------------------------------------------------------------------------
# Section données filtrées
# ---------------------------------------------------------------------------

st.divider()
st.header("Données filtrées (s1 actine)")


@st.cache_data
def _load_pdb_file(path, mtime=None):
    with open(path) as f:
        return f.read()


def extract_chain(pdb_data: str, chain_letter: str) -> str:
    lines = [l for l in pdb_data.splitlines()
             if l.startswith('ATOM')
             and len(l) > 21 and l[21] == chain_letter]
    return "\n".join(lines)


@st.cache_data
def load_csv(path, mtime=None):
    with open(path, newline="") as f:
        sample = f.read(10000)
    try:
        sep = csv.Sniffer().sniff(sample, delimiters=";,\t").delimiter
    except csv.Error:
        sep = ","
    return pd.read_csv(path, sep=sep)


def read_csv(path):
    """Charge un CSV en invalidant le cache si le fichier a changé."""
    mtime = os.path.getmtime(path) if os.path.exists(path) else None
    return load_csv(path, mtime=mtime)


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
    selected = st.selectbox("Choisir un tableau",
                            list(available_tables.keys()))
    df = read_csv(available_tables[selected])
    hide_constant = st.checkbox(
        "Masquer les colonnes sans variation", value=False)
    if hide_constant:
        cols_to_show = [c for c in df.columns if df[c].nunique() > 1]
    else:
        cols_to_show = list(df.columns)
    st.caption(
        f"{len(df):,} lignes · {len(cols_to_show)} colonnes affichées (sur {len(df.columns)} total)")
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
    col1.metric("Structures PDB",
                f"{df_pdb['pdb_id'].nunique():,}".replace(',', ' '))
    col2.metric("Interactions protéine-protéine",
                f"{len(df_int):,}".replace(',', ' '))
    col3.metric("Protéines partenaires uniques (ABP)",
                f"{df_prot[~df_prot['protein_name'].str.lower().str.contains('actin', na=False)]['protein_name'].nunique():,}".replace(',', ' '))


# ---------------------------------------------------------------------------
# Section PDB valides — explorateur
# ---------------------------------------------------------------------------

pdb_filt_path = "data/filtered/filtered_pdb_entry.csv"
if os.path.exists(pdb_filt_path):
    st.divider()
    st.header("Structures PDB valides")

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
            f"Choisir une structure ({len(pdb_ids)} PDB valides)",
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
        st.markdown(f"**Réseau d'interactions — {selected_pdb}**")
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
        st.markdown("**Visualisation 3D — contacts d'interface**")
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

                # chain_colors : dict chain_letter → couleur (None = pas de sélection)
                chain_colors = None

                if sel_inter and sel_inter in int_ids:
                    # Cas 1 : arête sélectionnée — A jaune, B vert
                    sel_row = df_int1_g[df_int1_g["interaction_id"]
                                        == sel_inter].iloc[0]
                    chain_colors = {
                        str(sel_row["chain_A_id"]).split("_")[-1]: "#FFD700",
                        str(sel_row["chain_B_id"]).split("_")[-1]: "#39FF14",
                    }
                elif sel_node:
                    # Cas 2 : nœud sélectionné — sélection en jaune,
                    # autres : orange=actine, vert=ABP
                    sel_letter = sel_node.split("_")[-1]
                    chain_colors = {}
                    for _letter, _is_actin in _chain_is_actin.items():
                        if _letter == sel_letter:
                            chain_colors[_letter] = "#FFD700"
                        else:
                            chain_colors[_letter] = (
                                "#E67E22" if _is_actin else "#2ECC71")

                fmt = "cif" if pdb_file.endswith(".cif") else "pdb"
                pdb_data = _load_pdb_file(
                    pdb_file, mtime=os.path.getmtime(pdb_file))

                viewer_key = (selected_pdb, sel_inter, sel_node)
                if st.session_state.get("viewer_key") != viewer_key:
                    import py3Dmol
                    view = py3Dmol.view(width=580, height=450)
                    view.addModel(pdb_data, fmt)
                    view.setStyle({}, {})
                    if chain_colors:
                        for chain, col in chain_colors.items():
                            view.addSurface(py3Dmol.SES,
                                            {"opacity": 1.0, "color": col},
                                            {"chain": chain})
                    else:
                        view.addSurface(py3Dmol.SES,
                                        {"opacity": 0.9, "color": "white"}, {})
                    view.setBackgroundColor("#1e1e1e")
                    view.zoomTo()
                    st.session_state["viewer_html"] = view._make_html()
                    st.session_state["viewer_key"] = viewer_key

                st.components.v1.html(
                    st.session_state["viewer_html"], height=470, scrolling=False)
            else:
                st.info("Fichier PDB non disponible.")

    # ── Tableau interactions du PDB ───────────────────────────────────────
    all_data_path_pdb = "data/filtered/filtered_all_data.csv"
    if os.path.exists(all_data_path_pdb):
        df_all_pdb = read_csv(all_data_path_pdb)
        sub_all_pdb = df_all_pdb[df_all_pdb["pdb_id"].str.upper(
        ) == selected_pdb]
        cols_to_show = [c for c in [
            "subunit_1", "subunit_2",
            "cluster_data_70",
            "s1_binding_site_cluster_data_70",
            "s2_binding_site_cluster_data_70",
        ] if c in sub_all_pdb.columns]
        if not sub_all_pdb.empty and cols_to_show:
            st.markdown("**Clusters d'interaction par paire**")
            st.dataframe(
                sub_all_pdb[cols_to_show].reset_index(drop=True),
                hide_index=True,
                width="stretch",
            )

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
                    if (s + k + 1) in iface:
                        cells += (
                            f'<span style="color:{col};font-weight:bold'
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
                        col_a="#FFD700", col_b="#39FF14", primary=None):
            ca_l, cb_l = ca.lower(), cb.lower()
            rseq = df_all_g[
                (df_all_g["subunit_1"].str.lower() == ca_l) &
                (df_all_g["subunit_2"].str.lower() == cb_l)
            ]
            if rseq.empty:
                st.info("Séquences non disponibles pour cette interaction.")
                return
            s1 = str(rseq.iloc[0].get("s1_sequence", ""))
            s2 = str(rseq.iloc[0].get("s2_sequence", ""))
            t1 = str(rseq.iloc[0].get("subunit_1_title", ca))
            t2 = str(rseq.iloc[0].get("subunit_2_title", cb))
            if iid is not None:
                sub_r3 = df_r3[df_r3["interaction_id"] == iid]
                r3_chain = sub_r3["chain"].str.lower()
                iface_a = set(
                    sub_r3[r3_chain == ca_l]["residue_number_sequence"]
                    .dropna().astype(int)
                )
                iface_b = set(
                    sub_r3[r3_chain == cb_l]["residue_number_sequence"]
                    .dropna().astype(int)
                )
            else:
                iface_a, iface_b = set(), set()
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
                        f"**{title}** · `{cid}` · {len(seq)} aa · {len(iface)} résidus d'interface"
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
                st.markdown("#### Séquences — résidus d'interface")
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
                    f"#### Séquences — `{seq_node.split('_')[-1].upper()}` "
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
                                col_a="#FFD700", col_b="#00E5FF",
                                primary=seq_node)

    # ── Protéines impliquées ──────────────────────────────────────────────
    proteins_path = "data/filtered/proteins_per_pdb.csv"
    if os.path.exists(proteins_path):
        st.markdown("**Protéines impliquées**")
        df_pp = read_csv(proteins_path)
        sub_pp = df_pp[df_pp["pdb_id"].str.upper() == selected_pdb]

        nb_actin = int(sub_pp["is_actin"].sum())
        nb_abp = int((~sub_pp["is_actin"]).sum())
        c1, c2, _ = st.columns([1, 1, 3])
        c1.metric("Actines", nb_actin)
        c2.metric("ABP", nb_abp)

        counts = (
            sub_pp.groupby(["protein"])
            .agg(
                Nombre=("chain", "count"),
                Chaînes=("chain", lambda x: ", ".join(
                    sorted(s.split("_")[-1].upper() for s in x)
                )),
            )
            .reset_index()
            .rename(columns={"protein": "Protéine"})
            .sort_values("Nombre", ascending=False)
            .reset_index(drop=True)
        )
        max_len = counts["Protéine"].str.len().max() if len(counts) else 0
        col_width = "large" if max_len > 40 else "medium" if max_len > 20 else "small"
        st.dataframe(counts, width="stretch", hide_index=True, column_config={
            "Protéine": st.column_config.TextColumn(width=col_width),
            "Chaînes": st.column_config.TextColumn(width="medium"),
            "Nombre": st.column_config.NumberColumn(width="small"),
        })
    else:
        st.info("Relancer le notebook `graphe_filter.ipynb` pour générer les données.")




st.divider()
st.header("Clusters d'intéractions")

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


if not os.path.exists(PATCHES_C70_CSV):
    st.info("Lancer le pipeline pour générer l'analyse des clusters.")
else:
    _all_data_path = "data/filtered/filtered_all_data.csv"
    _summary_path = "data/filtered/filtered_summary.csv"
    tab_s1, tab_c70 = st.tabs(
        ["Binding Site Cluster Data 70", "Interaction Cluster Data 70"])

    # --- S1 Binding Site ---
    with tab_s1:
        if os.path.exists(_all_data_path) and os.path.exists(_summary_path):
            st.components.v1.html(_build_global_graph_html(
                _all_data_path, _summary_path), height=780)
            st.divider()
        if os.path.exists(PATCHES_S1_CSV):

            df_s1 = read_csv(PATCHES_S1_CSV)
            st.caption(f"{len(df_s1)} patchs")

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

            # --- Heatmap global homo + hétéro ---
            st.subheader(
                "Heatmap — positions canoniques actin à l'interface S1")
            if not os.path.exists(S1_INTERFACE_FREQ_CSV):
                st.info(
                    "Lancer le notebook `notebooks/interface_analysis_s1.ipynb` pour générer les données.")
            else:
                heatmap_mode = st.selectbox(
                    "Normalisation",
                    ["Relative (max cluster = 1)", "Absolue (max = 100 %)"],
                    key="heatmap_norm_mode",
                )
                if heatmap_mode == "Relative (max cluster = 1)":
                    st.caption(
                        "Chaque ligne normalisée par son propre max — compare les motifs d'interface entre clusters.")
                    if os.path.exists(S1_HEATMAP_EQUITABLE):
                        st.image(S1_HEATMAP_EQUITABLE, width='content')
                    else:
                        st.info(
                            "Relancer le notebook `interface_analysis_s1.ipynb`.")
                else:
                    st.caption(
                        "Valeurs absolues en % buried ASA — compare l'intensité d'enfouissement entre clusters.")
                    if os.path.exists(S1_HEATMAP_ABSOLUTE):
                        st.image(S1_HEATMAP_ABSOLUTE, width='content')
                    else:
                        st.info(
                            "Relancer le notebook `interface_analysis_s1.ipynb` (cellule 6 — heatmap absolu).")

            st.divider()

            # --- Cluster sélectionné ---
            all_s1 = df_s1.sort_values("n_interactions", ascending=False)[
                "patch"].astype(str).tolist()
            sel_s1 = st.selectbox("Patch S1 binding site", all_s1, key="sel_s1",
                                  format_func=lambda p: f"{p} — {int(df_s1[df_s1['patch'].astype(str) == p]['n_interactions'].values[0])} interactions")

            row_s1 = df_s1[df_s1["patch"].astype(str) == sel_s1].iloc[0]

            # ── Réseau bipartite interactif + Surface 3D ─────────────────
            col_net_s1, col_3d_s1 = st.columns([3, 3])

            with col_net_s1:
                st.markdown(
                    "**Réseau interactif — résidus actine ↔ partenaires**")
                _bip_ok = all(os.path.exists(f) for f in _BIPARTITE_FILES)
                if _bip_ok:
                    _html_bip, _n_r, _n_p, _n_t = _build_bipartite_html(
                        sel_s1, _BIP_CACHE_VERSION, *_bip_mtimes())
                    if _html_bip:
                        st.caption(
                            f"{_n_r} résidus actine · {_n_p} protéines partenaires"
                            f" · n={_n_t} interactions")
                        st.components.v1.html(
                            _html_bip, height=540, scrolling=False)
                    else:
                        st.info("Réseau non disponible pour ce cluster.")
                else:
                    st.info("Lancer le pipeline pour générer les données réseau.")

            with col_3d_s1:
                st.markdown("**Interface 3D — couple représentatif**")
                if _bip_ok:
                    _s1_3d = _build_s1_3d_html(
                        sel_s1, _BIP_CACHE_VERSION, *_bip_mtimes())
                    if _s1_3d:
                        _html_3d_s1, _s1_max, _s2_max, _s1_is_homo = _s1_3d
                        st.components.v1.html(
                            _html_3d_s1, height=470, scrolling=False)
                        if _s1_is_homo:
                            st.caption(
                                f"S1 (actine) : beige → rouge · % ASA buried (0 → 100 %) · max observé : {_s1_max:.1f} %"
                            )
                        else:
                            st.caption(
                                f"S1 (actine) : beige → rouge · % ASA buried (0 → 100 %) · max : {_s1_max:.1f} % "
                                f"· S2 (ABP) : vert · % ASA, max {_s2_max:.1f} %"
                            )
                    else:
                        st.info("Vue 3D non disponible pour ce patch.")
                else:
                    st.info("Lancer le pipeline pour générer les données réseau.")

            # Heatmap 1 ligne pour ce cluster
            s1_cluster_img = os.path.join(
                S1_INTERFACE_CLUSTER_DIR, f"{sel_s1}.png")
            if os.path.exists(s1_cluster_img):
                st.image(s1_cluster_img, width='content')
            elif os.path.exists(S1_INTERFACE_FREQ_CSV):
                st.info("Heatmap non généré pour ce cluster (seuil non atteint).")

            # Heatmap séparé par sous-cluster C70
            s1_by_c70_img = os.path.join(
                S1_INTERFACE_CLUSTER_BY_C70_DIR, f"{sel_s1}.png")
            if os.path.exists(s1_by_c70_img):
                st.caption(
                    "Profil par sous-cluster d'interaction C70 (chaque ligne = 1 cluster d'interface)")
                st.image(s1_by_c70_img, width='content')

            # Téléchargement du script PyMOL pour ce cluster S1
            _pml_path = os.path.join(
                "data/filtered/details/structures_files/bfactor_c70_interface/by_s1_gradient",
                f"{sel_s1}.pml")
            if os.path.exists(_pml_path):
                with open(_pml_path, "rb") as _pml_f:
                    st.download_button(
                        label="⬇ Télécharger le script PyMOL de ce cluster",
                        data=_pml_f,
                        file_name=f"{sel_s1}.pml",
                        mime="text/plain",
                        help="Ouvrir dans PyMOL : File > Run Script…  ou  @/chemin/vers/le/fichier.pml",
                    )

            # ── Analyse des contacts ABP–actine (B–E) ────────────────────────
            st.divider()
            st.markdown("##### 📊 Analyse des contacts ABP–actine")
            _ch2seq_s1, _ch2title_s1 = _s1_get_ch_maps(sel_s1)
            if _ch2seq_s1:
                _msa_contact_analysis(
                    None, f"s1c_{sel_s1}",
                    _ch2seq=_ch2seq_s1, _ch2title=_ch2title_s1,
                    tabs=["B", "C", "D", "E"],
                )
            else:
                st.info("Aucune donnée de contact disponible pour ce cluster.")

    # --- Cluster Data 70 ---
    with tab_c70:
        if os.path.exists(_all_data_path):
            st.caption(
                "Rouge : binding sites actine (S1 + S2 homo) · Violet : clusters C70 · "
                "Vert : binding sites partenaire (hétéro) · "
                "Pointillé gris : interaction actine↔partenaire · Pointillé orange : interaction actine↔actine"
            )
            st.components.v1.html(_build_tripartite_graph_html(
                _all_data_path), height=900)
            st.divider()
        if os.path.exists(PATCHES_C70_CSV):
            df_c70 = read_csv(PATCHES_C70_CSV)
            st.caption(f"{len(df_c70)} patchs")

            df_c70_display = df_c70.drop(
                columns=["ids_interactions"], errors="ignore").copy()

            df_c70_display.index = range(1, len(df_c70_display) + 1)
            st.dataframe(df_c70_display, width="stretch")

            st.divider()
            all_c70 = df_c70.sort_values("n_interactions", ascending=False)[
                "patch"].astype(str).tolist()
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
                "**Réseau interactif — résidus actine ↔ résidus ABP**")
            _bip_ok = all(os.path.exists(f) for f in _BIPARTITE_FILES)
            if _bip_ok:
                _color_mode = "restype" if st.toggle(
                    "Coloration physicochimique (hydrophobe/polaire/chargé…)",
                    value=False, key=f"restype_toggle_{sel_c70}"
                ) else "bfactor"
                _html_c70, _n_s1, _n_s2, _n_tot, _html_3d_c70 = _build_bipartite_c70_html(
                    sel_c70, True, _color_mode, _BIP_CACHE_VERSION, *_bip_mtimes())
                if _html_c70:
                    st.caption(
                        f"{_n_s1} résidus actine (S1) · {_n_s2} résidus partenaire (S2)"
                        f" · n={_n_tot} interactions")
                    _net_height = 650
                    st.components.v1.html(
                        _html_c70, height=_net_height, scrolling=False)
                    if _html_3d_c70:
                        st.markdown("**Interface 3D — couple représentatif**")
                        st.components.v1.html(
                            _html_3d_c70, height=450, scrolling=False)
                else:
                    st.info("Réseau non disponible pour ce cluster.")
            else:
                st.info("Lancer le pipeline pour générer les données réseau.")

st.divider()
st.header("ABP")

# Vue globale des protéines non-actines
proteins_path = "data/filtered/proteins_per_pdb.csv"
_all_data_path = "data/filtered/filtered_all_data.csv"
if os.path.exists(proteins_path):
    st.subheader("Protéines non-actines (ABP) — global")
    df_pp_all = read_csv(proteins_path)
    df_abp = df_pp_all[~df_pp_all["is_actin"]]

    def _fmt_ids(series):
        vals = sorted(set(str(v).strip()
                      for v in series.dropna() if str(v) not in ("nan", "")))
        return ", ".join(vals)

    # Labels de position filament pour les actines
    _fil_pos_path = "data/filtered/actin_filament_positions.csv"
    _fil_label_map = {}
    if os.path.exists(_fil_pos_path):
        _fp = read_csv(_fil_pos_path)[["pdb_id", "chain", "label", "component_size"]]
        for _, _r in _fp.iterrows():
            if _r["component_size"] > 3:
                _fil_label_map[(_r["pdb_id"], _r["chain"])] = _r["label"]

    _label_order = {"-": 0, "-2": 1, "-3": 2, "side": 3, "+3": 4, "+2": 5, "+": 6}

    def _simplify_label(l):
        if l in ("+", "+2", "+3"):
            return "+"
        if l in ("-", "-2", "-3"):
            return "-"
        return l

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
            .rename(columns={"protein": "Protéine", "Nb_noeuds": "Nb nœuds"})
            .sort_values("Nb nœuds", ascending=False)
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
            .rename(columns={"protein": "Protéine", "Nb_noeuds": "Nb nœuds"})
            .sort_values("Nb nœuds", ascending=False)
            .reset_index(drop=True)
        )

    max_len = abp_global["Protéine"].str.len().max() if len(abp_global) else 0
    col_width = "large" if max_len > 40 else "medium" if max_len > 20 else "small"
    st.dataframe(abp_global, width="stretch", hide_index=True, column_config={
        "Protéine": st.column_config.TextColumn(width=col_width),
        "Nb nœuds": st.column_config.NumberColumn(width="small"),
        "PDB": st.column_config.TextColumn(width="large"),
        "Binding site S1": st.column_config.TextColumn(width="medium"),
        "Cluster C70": st.column_config.TextColumn(width="medium"),
        "Localisation filament": st.column_config.TextColumn(width="medium"),
    })

# ── Heatmap ABP × résidus actin ──────────────────────────────────────────────
st.subheader("Heatmap — résidus actin contactés par les ABP")
st.caption(
    "Chaque cellule = % ASA buried moyen équitable du résidu actin (position canonique). "
    "Moyenne en deux étapes : d'abord par cluster C70 (0% si pas de contact), "
    "puis moyenne équitable entre clusters C70."
)


if all(os.path.exists(f) for f in _ABP_HM_FILES):
    _hm_mtimes = tuple(
        os.path.getmtime(f) if os.path.exists(f) else 0.0
        for f in _ABP_HM_FILES
    )
    _hm_res = _build_abp_heatmap_data(*_hm_mtimes)
    if _hm_res is not None:
        _pivot_full_hm, _abp_freq_hm, _ = _hm_res

        _pivot_hm = _pivot_full_hm.copy()

        import seaborn as sns
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist

        _n_rows = len(_pivot_hm)
        _n_cols = len(_pivot_hm.columns)
        _fig_w = max(14, _n_cols * 0.35 + 4)
        _fig_h = max(6, _n_rows * 0.45 + 2)

        # Clustering hiérarchique sur les lignes (ABP) uniquement
        _mat = _pivot_hm.fillna(0).values.astype(float)
        _row_order = leaves_list(
            linkage(pdist(_mat, metric="cosine"), method="average"))
        _pivot_clustered = _pivot_hm.iloc[_row_order].fillna(0)

        # Labels Y avec n interactions (après réordonnancement)
        _ylabels = [
            f"{r}  ({_abp_freq_hm.get(r, 0)})"
            for r in _pivot_clustered.index
        ]

        _fig_hm, _ax_hm = plt.subplots(figsize=(_fig_w, _fig_h))
        sns.heatmap(
            _pivot_clustered,
            ax=_ax_hm,
            cmap="YlOrRd",
            vmin=0, vmax=100,
            xticklabels=True,
            yticklabels=_ylabels,
            linewidths=0,
            cbar_kws={
                "label": "% ASA buried moy (0 si absent)", "shrink": 0.6},
        )
        _ax_hm.set_xlabel("Position canonique (MAFFT)",
                          fontsize=11, labelpad=6)
        _ax_hm.set_ylabel("")
        _ax_hm.set_xticklabels(_ax_hm.get_xticklabels(),
                               rotation=90, fontsize=8)
        _ax_hm.set_yticklabels(_ax_hm.get_yticklabels(),
                               rotation=0, fontsize=22)
        _ax_hm.set_title(
            f"% ASA buried moyen — résidus actin × ABP  ({_n_rows} ABP, {_n_cols} résidus)",
            fontsize=11, pad=8,
        )
        plt.tight_layout()
        st.pyplot(_fig_hm)
        plt.close(_fig_hm)

# ── Réseau de compétition ABP ─────────────────────────────────────────────────
if (os.path.exists(proteins_path) and os.path.exists(_all_data_path)
        and len(abp_global) > 0 and "Binding site S1" in abp_global.columns):
    st.subheader("Réseau de compétition ABP")
    _jac_thresh = st.select_slider(
        "Seuil de recouvrement C70",
        options=[0.25, 0.30, 0.40, 0.50, 0.60,
                 0.70, 0.75, 0.80, 0.90, 1.00],
        value=0.50, key="jac_thresh",
        format_func=lambda x: f"{int(x*100)}%",
    )
    _min_edge_w = 1
    # Dataset étendu : actine peut être S1 ou S2 (toutes données hétéro)
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
    _abp_c70_valid: dict = {}                # ABP → set de C70 valides (Jaccard + position)
    for (_a, _b) in _edge_wts:
        # Recompter uniquement les paires C70 qui passent Jaccard ET position
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
            _abp_fil_pos_net.setdefault(_r["protein"], set()).add(_simplify_label(_lbl))

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

    st.caption(
        f"Seuil de recouvrement {int(_jac_thresh*100)}% — deux ABPs sont liés si au moins "
        f"un de leurs clusters C70 couvre ≥{int(_jac_thresh*100)}% de la plus petite interface. "
        f"{len(_edge_wts)} paires en compétition "
        f"({_n_filtered} exclues car localisations filament disjointes)."
    )

    # Tous les ABPs connus (y compris isolés)
    _all_abp_net = set(_merged_net["protein"].dropna().unique())

    if True:  # toujours entrer (même sans arêtes, afficher les nœuds isolés)
        # Nombre de compétiteurs uniques (voisins distincts) par ABP
        _abp_neighbors_c: dict = {}
        for (_a_c, _b_c) in _edge_wts:
            _abp_neighbors_c.setdefault(_a_c, set()).add(_b_c)
            _abp_neighbors_c.setdefault(_b_c, set()).add(_a_c)
        _abp_n_competitors: dict = {k: len(v)
                                    for k, v in _abp_neighbors_c.items()}

        _max_comp_c = max(_abp_n_competitors.values()) if _abp_n_competitors else 1

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
        _net_c = Network(height="980px", width="100%", bgcolor="#ffffff",
                         font_color="#111")
        _net_c.set_options(
            '{"physics": {"enabled": true, "solver": "barnesHut",'
            '  "barnesHut": {"gravitationalConstant": -8000,'
            '    "centralGravity": 1, "springLength": 100,'
            '    "springConstant": 0.05, "damping": 0.2, "avoidOverlap": 1.0},'
            '  "stabilization": {"enabled": true, "iterations": 2000, "fit": true}},'
            ' "nodes": {"shape": "dot", "font": {"size": 7}},'
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
            t = (pct - _pct_min) / (_pct_max - _pct_min) if _pct_max > _pct_min else 1.0
            # jaune (#ffff00) → rouge (#cc0000)
            r = int(255 + (204 - 255) * t)
            g = int(255 * (1 - t))
            b = 0
            return f"#{r:02x}{g:02x}{b:02x}"

        _label_thresh = 3
        for _nc in sorted(_all_abp_net):
            _n_c70_total = _abp_c70_count.get(_nc, 0)
            _n_c70_comp = _abp_c70_comp_count.get(_nc, 0)
            _n_comp_nc = _abp_n_competitors.get(_nc, 0)
            _sz = 5 + 18 * (_n_comp_nc / _max_comp_c) ** 1.5 if _n_comp_nc > 0 else 5
            _pct_val = round(_n_c70_comp / _n_c70_total *
                             100) if _n_c70_total else 0
            _col = _pct_color(_n_c70_comp, _n_c70_total)
            _show_lbl = _n_comp_nc >= _label_thresh
            _c70_line = (f"C70 en compétition : {_n_c70_comp} / {_n_c70_total} ({_pct_val}%)"
                         if _n_c70_comp else f"C70 total : {_n_c70_total} (0%)")
            _pos_set = _abp_fil_pos_net.get(_nc, set())
            _pos_str = ", ".join(
                sorted(_pos_set, key=lambda x: {"-": 0, "side": 1, "+": 2}.get(x, 99))
            ) if _pos_set else "?"
            _border_col = _pos_border_color(_pos_set)
            _net_c.add_node(
                _nc,
                label=_wrap_lbl(_nc) if _show_lbl else "",
                title=(f"{_nc}\n"
                       f"Position filament : {_pos_str}\n"
                       f"{_c70_line}\n"
                       f"Partenaires ABP uniques en concurrence : {_n_comp_nc}"),
                size=_sz,
                color={"background": _col, "border": _border_col},
                borderWidth=1,
                font={"size": 3, "face": "Arial", "color": "#111111"},
            )
        _filtered_edges = {k: v for k,
                           v in _edge_wts.items() if v >= _min_edge_w}
        for (_a_c, _b_c), _w_c in sorted(_filtered_edges.items(),
                                         key=lambda x: x[1], reverse=True):
            _net_c.add_edge(
                _a_c, _b_c,
                width=1,
                title=f"{_w_c} paire(s) C70 partagée(s)",
                color={"color": "#505050", "opacity": 0.2},
            )
        _net_html = _net_c.generate_html()
        # Désactiver la physique dès que la stabilisation est terminée
        for _pat_c in [
            "network = new vis.Network(container, data, options);",
            "var network = new vis.Network(container, data, options);",
        ]:
            if _pat_c in _net_html:
                _net_html = _net_html.replace(
                    _pat_c,
                    _pat_c + "\n  network.once('stabilizationIterationsDone', function(){"
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
            "  box.placeholder = 'Rechercher un ABP…';\n"
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
        _net_html = _net_html.replace("</body>", "<script>" + _inject_search_abp + "</script>\n</body>")
        st.components.v1.html(_net_html, height=880, scrolling=False)
        st.caption(
            "Remplissage : jaune = % min de clusters C70 en compétition → rouge = % max   ·   "
            "Taille : nombre de partenaires en compétition   ·   "
            "Survol : position filament de l'ABP"
        )
        # Légende des couleurs de contour
        _legend_html = (
            "<div style='display:flex;flex-wrap:wrap;gap:14px;align-items:center;"
            "font-size:12px;margin-top:4px;'>"
            "<b>Contour = position filament :</b>"
            "<span style='display:flex;align-items:center;gap:5px;'>"
            "<span style='width:14px;height:14px;border-radius:50%;border:2px solid #e67300;"
            "background:#eee;display:inline-block'></span> barbée (+ ou + side)</span>"
            "<span style='display:flex;align-items:center;gap:5px;'>"
            "<span style='width:14px;height:14px;border-radius:50%;border:2px solid #0077cc;"
            "background:#eee;display:inline-block'></span> pointée (− ou − side)</span>"
            "<span style='display:flex;align-items:center;gap:5px;'>"
            "<span style='width:14px;height:14px;border-radius:50%;border:2px solid #888888;"
            "background:#eee;display:inline-block'></span> latéral / mixte</span>"
            "</div>"
        )
        st.markdown(_legend_html, unsafe_allow_html=True)

    st.divider()

# ── Détail par ABP ─────────────────────────────────────────────────────────────
if (os.path.exists(proteins_path) and os.path.exists(_all_data_path)
        and len(abp_global) > 0 and "Binding site S1" in abp_global.columns):
    st.subheader("Détail par ABP")

    _NO_ABP_LABEL = "— PDB sans ABP —"
    abp_names = [_NO_ABP_LABEL] + abp_global["Protéine"].tolist()
    sel_abp = st.selectbox("Sélectionner un ABP",
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
        st.markdown("#### Clusters d'interaction actine–ABP")
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
            st.info("Aucune interaction hétéro trouvée pour cet ABP.")
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
            abp_clust["% interactions ABP-actine"] = (
                abp_clust["nb_inter"] / _total_abp_inter * 100
            ).round(1)
            st.caption(
                f"{len(abp_pdbs)} PDB contenant cet ABP · "
                f"{len(abp_int)} interactions hétéro actine–ABP · "
                f"{abp_clust.shape[0]} clusters C70"
            )
            abp_clust["Site liaison"] = (
                abp_clust["s1_binding_site_cluster_data_70"].astype(str)
                + " × "
                + abp_clust["s2_binding_site_cluster_data_70"].astype(str)
            )
            st.dataframe(
                abp_clust.rename(columns={
                    "cluster_data_70": "Cluster C70",
                    "nb_inter": "Nb interactions",
                    "nb_pdb": "Nb PDB",
                })[["Cluster C70", "Site liaison",
                    "Nb interactions", "% interactions ABP-actine", "Nb PDB", "% PDB"]],
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
                st.caption(
                    "Aucun fichier PDB C70 disponible pour ces clusters.")
            else:
                _ylord3d = ["#FFFFCC", "#FFF0A9", "#FEE186", "#FECA65", "#FDAA48",
                            "#FC8C3B", "#FC5A2D", "#EC2D21", "#D30F20", "#800026"]
                _grn3d = ["#FFFFCC", "#D9F0A3", "#ADDD8E", "#78C679",
                          "#41AB5D", "#238443", "#006837"]
                st.caption(
                    "Jaune→rouge = actine (% ASA enfouie) · "
                    "Jaune→vert = ABP · blanc = hors interface"
                )
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

    # ── Interactions actine-actine (homo) ─────────────────────────────────────
    _homo_label = ("PDB sans ABP" if _is_no_abp
                   else f"les mêmes PDB que {sel_abp[:30]}")
    st.markdown(f"#### Interactions actine-actine (homo) — {_homo_label}")
    _df_homo = df_all_g.copy()
    _df_homo["_pdb"] = _df_homo["subunit_1"].str.split("_").str[0]
    homo_cooc = _df_homo[
        _df_homo["s1_actine"].fillna(False).astype(bool) &
        _df_homo["s2_actine"].fillna(False).astype(bool) &
        _df_homo["_pdb"].isin(abp_pdbs) &
        _df_homo["cluster_data_70"].notna()
    ][["_pdb", "cluster_data_70",
       "s1_binding_site_cluster_data_70", "s2_binding_site_cluster_data_70"]].copy()

    # Fusionner S1 et S2 en une paire canonique (triée, S1/S2 non distingués)
    homo_cooc["Binding sites"] = homo_cooc.apply(
        lambda r: " × ".join(sorted([
            str(r["s1_binding_site_cluster_data_70"]),
            str(r["s2_binding_site_cluster_data_70"])
        ])), axis=1
    )

    total_homo = len(homo_cooc)
    st.caption(
        f"{len(abp_pdbs)} PDB · "
        f"{total_homo} interactions homo actine-actine"
    )
    if homo_cooc.empty:
        st.info("Aucune interaction homo actine-actine dans ces PDB.")
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
        homo_summary["% interactions homo"] = (
            homo_summary["nb_inter"] / total_homo * 100
        ).round(1)
        homo_summary = homo_summary.rename(columns={
            "cluster_data_70": "Cluster C70",
            "nb_pdb": "Nb PDB",
            "nb_inter": "Nb interactions homo",
        })[["Cluster C70", "Binding sites",
            "Nb PDB", "% PDB", "Nb interactions homo", "% interactions homo"]]
        st.dataframe(homo_summary, hide_index=True, use_container_width=True)

    # ── Filament actine 10 sous-unités ────────────────────────────────────
    st.subheader("Filament actine — 10 sous-unités (PyMOL)")

    _by_abp_dir = _Path(
        "data/filtered/details/structures_files/filament/by_abp"
    )
    _global_base_pdb = _Path(
        "data/filtered/details/structures_files/filament/filament_global_base.pdb"
    )

    if not _is_no_abp and _by_abp_dir.exists():
        import re as _re
        _abp_sname = _re.sub(r"[^a-zA-Z0-9]+", "_", sel_abp).strip("_")[:60]
        _abp_pml = _by_abp_dir / f"{_abp_sname}.pml"
        _abp_pdb = _by_abp_dir / f"{_abp_sname}_abp.pdb"

        # Lire patches_by_abp.csv pour savoir si même filament que base
        _patches_csv = _by_abp_dir / "patches_by_abp.csv"
        _same_as_base = False
        _abp_patch1 = _abp_patch2 = None
        if _patches_csv.exists():
            _df_patches_abp = pd.read_csv(_patches_csv)
            _row_abp = _df_patches_abp[_df_patches_abp["abp_sname"] == _abp_sname]
            if not _row_abp.empty:
                _same_as_base = bool(_row_abp.iloc[0]["same_as_base"])
                _abp_patch1 = _row_abp.iloc[0]["patch_1"]
                _abp_patch2 = _row_abp.iloc[0]["patch_2"]

        if _same_as_base:
            st.info(
                f"Le filament de {sel_abp} utilise les **memes clusters d'interface que le filament de reference** "
                f"({_abp_patch1} + {_abp_patch2}). "
                "La session PyMOL contient uniquement le filament global (identique)."
            )
        else:
            _p_str = ""
            if _abp_patch1:
                _p_str = f"{_abp_patch1}" + (f" + {_abp_patch2}" if pd.notna(_abp_patch2) else "")
            st.caption(
                "Session PyMOL avec 2 filaments : "
                "**filament_base** (clusters globaux 0_7797_0 + 0_7797_1, blanc→rouge) et "
                f"**filament_abp** (clusters {_p_str}, blanc→orange). "
                "B-factor = somme % ASA buried des sites de liaison S1."
            )

        _fc1, _fc2, _fc3 = st.columns(3)
        if _abp_pml.exists():
            # Lire le PML et remplacer les chemins par des chemins absolus résolus sur cette machine
            _pml_text = _abp_pml.read_text()
            _fc1.download_button(
                "Session PyMOL (.pml)",
                _pml_text.encode(),
                file_name=f"{_abp_sname}.pml",
                mime="text/plain",
                key=f"dl_pml_{_abp_sname}",
            )
        else:
            _fc1.info("Session PML non disponible.")
        if _global_base_pdb.exists():
            with open(_global_base_pdb, "rb") as _f:
                _fc2.download_button(
                    "PDB filament base (.pdb)",
                    _f,
                    file_name="filament_global_base.pdb",
                    mime="chemical/x-pdb",
                    key=f"dl_base_{_abp_sname}",
                )
        if not _same_as_base and _abp_pdb.exists():
            with open(_abp_pdb, "rb") as _f:
                _fc3.download_button(
                    f"PDB filament {sel_abp[:20]} (.pdb)",
                    _f,
                    file_name=f"{_abp_sname}_abp.pdb",
                    mime="chemical/x-pdb",
                    key=f"dl_abp_{_abp_sname}",
                )
    else:
        # Filament générique (PDB sans ABP ou fichiers ABP non disponibles)
        st.caption(
            "Filament reconstruit depuis les structures pairwise 8iah (chaînes A→J). "
            "B-factor = somme % ASA buried des 4 clusters homo principaux (6685_1+2+3+4). "
            "Spectre blanc → rouge : zones les plus contactées dans les interactions homo actine-actine."
        )

        # Tableau récapitulatif : quels ABPs ont le même filament que la base ?
        _patches_csv_noabp = _Path(
            "data/filtered/details/structures_files/filament/by_abp/patches_by_abp.csv"
        )
        if _patches_csv_noabp.exists():
            _df_abp_patches = pd.read_csv(_patches_csv_noabp)
            _df_abp_patches = _df_abp_patches[_df_abp_patches["patch_1"].notna()].copy()
            _df_abp_patches["Filament"] = _df_abp_patches["same_as_base"].map(
                {True: "Identique a la reference", False: "Specifique a l'ABP"}
            )
            _df_abp_patches["Clusters utilises"] = _df_abp_patches.apply(
                lambda r: str(r["patch_1"]) + (
                    f" + {r['patch_2']}" if pd.notna(r["patch_2"]) else ""), axis=1
            )
            # Mapping C70 → S1 supercluster
            _filt_path = _Path("data/filtered/filtered_all_data.csv")
            if _filt_path.exists():
                _df_filt_sc = pd.read_csv(_filt_path, low_memory=False,
                                          usecols=["cluster_data_70", "s1_supercluster"])
                _c70_to_sc = (
                    _df_filt_sc.dropna(subset=["s1_supercluster"])
                    .groupby("cluster_data_70")["s1_supercluster"]
                    .agg(lambda x: x.mode()[0] if len(x) > 0 else "")
                    .to_dict()
                )
                def _get_superfamille(r):
                    sc1 = _c70_to_sc.get(str(r["patch_1"]), "")
                    sc2 = _c70_to_sc.get(str(r["patch_2"]), "") if pd.notna(r.get("patch_2")) else ""
                    parts = sorted({s for s in [sc1, sc2] if s})
                    return " + ".join(parts) if parts else "—"
                _df_abp_patches["Superfamille S1"] = _df_abp_patches.apply(_get_superfamille, axis=1)
                _cols = ["abp_title", "Clusters utilises", "Superfamille S1", "Filament"]
            else:
                _cols = ["abp_title", "Clusters utilises", "Filament"]
            _df_display = _df_abp_patches[_cols].rename(columns={"abp_title": "ABP"})
            st.markdown("**Filament par ABP — comparaison avec la reference sans ABP**")
            st.dataframe(
                _df_display.sort_values("Filament"),
                hide_index=True,
                use_container_width=True,
                column_config={
                    "Filament": st.column_config.TextColumn("Filament", width="medium"),
                    "Superfamille S1": st.column_config.TextColumn("Superfamille S1", width="medium"),
                },
            )
        # Session PyMOL avec tous les filaments spécifiques — générée à la volée
        _by_abp_dir2 = _Path("data/filtered/details/structures_files/filament/by_abp")
        _global_base2 = _Path("data/filtered/details/structures_files/filament/filament_global_base.pdb")
        _patches_csv2 = _by_abp_dir2 / "patches_by_abp.csv"
        if _patches_csv2.exists() and _global_base2.exists():
            import re as _re2
            _df_p2 = pd.read_csv(_patches_csv2)
            _diff2 = _df_p2[_df_p2["same_as_base"] == False].dropna(subset=["patch_1"])
            _pml_lines2 = [
                "# PyMOL — Tous les filaments actine uniques (base + ABPs specifiques)",
                "# ABPs charges mais caches : faire 'show surface, NomObjet' pour les afficher",
                "",
                f"load {_global_base2.resolve().as_posix()}, filament_base",
                "hide everything, filament_base",
                "show surface, filament_base",
                "spectrum b, white_red, filament_base, minimum=0",
                "",
            ]
            for _, _rec2 in _diff2.iterrows():
                _pdb2 = _by_abp_dir2 / f"{_rec2['abp_sname']}_abp.pdb"
                if not _pdb2.exists():
                    continue
                _obj2 = _re2.sub(r"[^a-zA-Z0-9]+", "_", str(_rec2["abp_title"]))[:30].strip("_")
                _pml_lines2 += [
                    f"# {_rec2['abp_title']} — faire 'show surface, {_obj2}' pour afficher",
                    f"load {_pdb2.resolve().as_posix()}, {_obj2}",
                    f"hide everything, {_obj2}",
                    f"spectrum b, white_red, {_obj2}, minimum=0",
                    "",
                ]
            _pml_lines2 += ["set surface_quality, 0", "bg_color white", "zoom filament_base"]
            _pml_content2 = "\n".join(_pml_lines2).encode()
            st.download_button(
                "Session PyMOL — tous les filaments specifiques (.pml)",
                _pml_content2,
                file_name="all_specific_filaments.pml",
                mime="text/plain",
                key="dl_all_filaments",
            )
            st.caption(
                "Contient le filament de reference + les ABPs dont le filament "
                "differe de la base. Chemins absolus generes depuis cette machine."
            )

        _filament_pml = _Path(
            "data/filtered/details/structures_files/filament/actin_filament_10.pml"
        )
        _filament_pdb = _Path(
            "data/filtered/details/structures_files/filament/actin_filament_10.pdb"
        )
        _col_pml, _col_pdb = st.columns(2)
        if _filament_pml.exists():
            with open(_filament_pml, "rb") as _f:
                _col_pml.download_button(
                    "Télécharger le script PyMOL (.pml)",
                    _f,
                    file_name="actin_filament_10.pml",
                    mime="text/plain",
                )
        else:
            _col_pml.info("Script PyMOL non disponible.")
        if _filament_pdb.exists():
            with open(_filament_pdb, "rb") as _f:
                _col_pdb.download_button(
                    "Télécharger le PDB filament (.pdb)",
                    _f,
                    file_name="actin_filament_10.pdb",
                    mime="chemical/x-pdb",
                )

# ── MSA — Protéines d'interface ───────────────────────────────────────────────

st.header("MSA — Protéines d'interface", anchor="msa-proteines")


if not _Path("data/filtered/filtered_all_data.csv").exists():
    st.warning("filtered_all_data.csv non trouvé — lancer d'abord les étapes de filtrage.")
else:
    _MSA_ALN_DIR.mkdir(parents=True, exist_ok=True)

    st.caption(
        "3 familles + clusters : MSA séquence complète ABP (S2) et Actine (S1) — "
        "🟢 vert = majorité en interaction (aa conservé) · "
        "🟣 violet = majorité en interaction (aa variable) · "
        "🔴 rouge = minorité en interaction"
    )

    _msa_section_full(
        "Myosines (rigor — sans ADP)", "myosin",
        lambda t: "myosin" in t.lower() and "tropomyosin" not in t.lower(),
        rigor_pdbs=_MSA_RIGOR_PDBS,
        note="Restreint aux 9 structures rigor (sans ADP sur la chaîne myosine).",
    )
    _msa_section_full(
        "Tropomyosines", "tropomyosin",
        lambda t: "tropomyosin" in t.lower(),
    )
    _msa_section_full(
        "Plastin et apparentés", "plastin",
        lambda t: any(x in t.lower() for x in ["plastin", "spectrin beta", "filamin", "utrophin"]),
        note="Plastin-3, Spectrin beta chain, Filamin-A, Utrophin.",
    )
    _msa_section_s2_clusters()
    _msa_section_s1_clusters()
