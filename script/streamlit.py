# streamlit run script/streamlit.py
import csv
import os
import re
import subprocess
import sys
import pandas as pd
import streamlit as st

st.set_page_config(layout="wide", page_title="Analyse actine-ABP - PPI3D")

st.title("Analyse des interactions actine-actine et actine-ABP - PPI3D")

with st.sidebar:
    st.markdown("## Sommaire")
    st.markdown("""
- [Téléchargement des données](#telechargement-des-donnees)
- [Données filtrées](#donnees-filtrees-s1-actine)
- [Structures PDB valides](#structures-pdb-valides)
- [Clusters d'interactions](#clusters-d-interactions)
""")

# ---------------------------------------------------------------------------
# Section téléchargement
# ---------------------------------------------------------------------------

st.header("Téléchargement des données")

STEPS = {
    "1/6": "Téléchargement du summary PPI3D (BLAST)",
    "2/6": "Téléchargement des PDB issues du summary",
    "3/6": "Téléchargement de toute les données (no clustering)",
    "4/6": "Filtrage des structures (≥ 5 actines) - notebook",
    "5/6": "Filtrage summary + téléchargement des détails d'interface",
    "6/6": "Filtrage all_data",
}

# Fichier de sortie attendu pour chaque étape
STEP_OUTPUT_FILES = {
    "1/6": "data/raw/ppi3d_actin_summary.csv",
    "2/6": "data/raw/pdb_entry_results.csv",
    "3/6": "data/raw/all_data.csv",
    "4/6": "data/filtered/filtered_pdb_entry.csv",
    "5/6": "data/filtered/details/1.interactions.csv",
    "6/6": "data/filtered/filtered_all_data.csv",
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

    # Étape 6 : "X lignes traitées"
    m = re.search(r"([\d,]+)\s+lignes traitées", line)
    if m:
        return None, None, line.strip()

    return None


# Bouton en haut, étapes en dessous

clicked = st.button("Lancer le téléchargement")

progress_bar = st.empty()

placeholders = {}
for key, label in STEPS.items():
    placeholders[key] = st.empty()
    placeholders[key].markdown(step_md(key, label, initial_state(key)))

sub_progress_bar = st.empty()
sub_progress_text = st.empty()

if "just_downloaded" not in st.session_state:
    st.session_state.just_downloaded = False

all_done = all(initial_state(k) == "skipped" for k in STEP_KEYS)
if all_done and not st.session_state.just_downloaded:
    st.write("✅ Données déjà téléchargées")


if clicked:
    st.session_state.just_downloaded = True

    for key, label in STEPS.items():
        placeholders[key].markdown(step_md(key, label, "pending"))

    pb = progress_bar.progress(0, text="Démarrage...")

    current_idx = -1
    step_had_skip = {}

    proc = subprocess.Popen(
        [sys.executable, "-u", "-m", "script.data_extract.pipeline_data"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    for line in proc.stdout:
        line = line.rstrip()

        # Détecter le début d'une nouvelle étape
        for i, key in enumerate(STEP_KEYS):
            if "ETAPE" in line and key in line:
                if current_idx >= 0:
                    prev_key = STEP_KEYS[current_idx]
                    state = "skipped" if step_had_skip.get(
                        prev_key) else "done"
                    placeholders[prev_key].markdown(
                        step_md(prev_key, STEPS[prev_key], state))

                current_idx = i
                placeholders[key].markdown(step_md(key, STEPS[key], "running"))
                pb.progress(i / TOTAL, text=f"Étape {key} — {STEPS[key]}")
                sub_progress_bar.empty()
                sub_progress_text.empty()
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

    proc.wait()

    sub_progress_bar.empty()
    sub_progress_text.empty()

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


# ---------------------------------------------------------------------------
# Section données filtrées
# ---------------------------------------------------------------------------

st.divider()
st.header("Données filtrées (s1 actine)")


@st.cache_data
def load_csv(path):
    with open(path, newline="") as f:
        sep = csv.Sniffer().sniff(f.read(10000), delimiters=";,\t").delimiter
    return pd.read_csv(path, sep=sep)


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
    df = load_csv(available_tables[selected])
    hide_constant = st.checkbox(
        "Masquer les colonnes sans variation", value=False)
    if hide_constant:
        cols_to_show = [c for c in df.columns if df[c].nunique() > 1]
    else:
        cols_to_show = list(df.columns)
    st.caption(
        f"{len(df):,} lignes · {len(cols_to_show)} colonnes affichées (sur {len(df.columns)} total)")
    st.dataframe(df[cols_to_show], use_container_width=True)

# Métrique PDB
METRICS_FILES = {
    "pdb": "data/filtered/filtered_pdb_entry.csv",
    "interactions": "data/filtered/details/1.interactions.csv",
    "proteins": "data/filtered/details/2.proteins.csv",
    "residues": "data/filtered/details/3.interface_residues.csv",
    "ligands": "data/filtered/details/5.ligands.csv",
}

if all(os.path.exists(p) for p in METRICS_FILES.values()):
    df_pdb = load_csv(METRICS_FILES["pdb"])
    df_int = load_csv(METRICS_FILES["interactions"])
    df_prot = load_csv(METRICS_FILES["proteins"])
    df_res = load_csv(METRICS_FILES["residues"])
    df_lig = load_csv(METRICS_FILES["ligands"])

    nb_homo = int((df_pdb["Interface type"] == "homo").sum()
                  ) if "Interface type" in df_pdb.columns else 0
    nb_hetero = int((df_pdb["Interface type"] == "hetero").sum(
    )) if "Interface type" in df_pdb.columns else 0

    col1, col2, col3 = st.columns(3)
    col1.metric("Structures PDB",
                f"{df_pdb['pdb_id'].nunique():,}".replace(',', ' '))
    col2.metric("Interactions protéine-protéine",
                f"{len(df_int):,}".replace(',', ' '))
    col3.metric("Protéines uniques (actines et ABP)",
                f"{df_prot['protein_name'].nunique():,}".replace(',', ' '))

# ---------------------------------------------------------------------------
# Section PDB valides — explorateur
# ---------------------------------------------------------------------------

pdb_filt_path = "data/filtered/filtered_pdb_entry.csv"
if os.path.exists(pdb_filt_path):
    st.divider()
    st.header("Structures PDB valides")

    df_entry = load_csv(pdb_filt_path)
    pdb_ids = sorted(df_entry["pdb_id"].str.upper().unique())

    summary_path = "data/filtered/filtered_summary.csv"
    if os.path.exists(summary_path):
        df_sum = load_csv(summary_path)
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

    col_net, col_prot = st.columns([3, 2])

    with col_net:
        st.markdown(
            f"**Réseau d'interactions — {selected_pdb} — {title_map[selected_pdb]}**")

        img_path = f"visualisations/pdb_graphs/{selected_pdb}.png"
        if os.path.exists(img_path):
            st.image(img_path, width=650)
        else:
            st.info("Image non générée — relancer le téléchagement.")

    with col_prot:
        st.markdown("**Protéines impliquées**")
        proteins_path = "data/filtered/proteins_per_pdb.csv"
        if os.path.exists(proteins_path):
            df_pp = load_csv(proteins_path)
            sub_pp = df_pp[df_pp["pdb_id"] == selected_pdb]

            nb_actin = int(sub_pp["is_actin"].sum())
            nb_abp = int((~sub_pp["is_actin"]).sum())

            c1, c2 = st.columns(2)
            c1.metric("Actines", nb_actin)
            c2.metric("ABP", nb_abp)

            counts = (
                sub_pp.groupby(["protein"])
                .size()
                .reset_index(name="Nombre")
                .rename(columns={"protein": "Protéine"})
                .reset_index(drop=True)
            )
            max_len = counts["Protéine"].str.len().max() if len(counts) else 0
            col_width = "large" if max_len > 40 else "medium" if max_len > 20 else "small"
            st.dataframe(counts, use_container_width=True, hide_index=True, column_config={
                "Protéine": st.column_config.TextColumn(width=col_width),
                "Nombre": st.column_config.NumberColumn(width="small"),
            })
        else:
            st.info(
                "Relancer le notebook `graphe_filter.ipynb` pour générer les données.")

    # Vue globale des protéines non-actines
    proteins_path = "data/filtered/proteins_per_pdb.csv"
    if os.path.exists(proteins_path):
        st.subheader("Protéines non-actines (ABP) — global")
        df_pp_all = load_csv(proteins_path)
        df_abp = df_pp_all[~df_pp_all["is_actin"]]

        abp_global = (
            df_abp.groupby("protein")
            .agg(
                Nb_noeuds=("chain", "count"),
                PDB=("pdb_id", lambda x: ", ".join(sorted(x.unique()))),
            )
            .reset_index()
            .rename(columns={"protein": "Protéine", "Nb_noeuds": "Nb nœuds"})
            .sort_values("Nb nœuds", ascending=False)
            .reset_index(drop=True)
        )

        max_len = abp_global["Protéine"].str.len(
        ).max() if len(abp_global) else 0
        col_width = "large" if max_len > 40 else "medium" if max_len > 20 else "small"
        st.dataframe(abp_global, use_container_width=True, hide_index=True, column_config={
            "Protéine": st.column_config.TextColumn(width=col_width),
            "Nb nœuds": st.column_config.NumberColumn(width="small"),
            "PDB": st.column_config.TextColumn(width="large"),
        })


# ---------------------------------------------------------------------------
# Section clusters
# ---------------------------------------------------------------------------

st.divider()
st.header("Clusters d'intéractions")
