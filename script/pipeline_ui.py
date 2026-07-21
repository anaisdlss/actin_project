"""Data-download / pipeline section, extracted from streamlit.py (keeps it light)."""
import os
import re
import subprocess
import sys
import pandas as pd
import streamlit as st


def render():
    st.header("Data download", anchor="telechargement-des-donnees")
    st.caption("— Technical section (reproducibility). Fetches and updates the "
               "PPI3D data. A biologist can skip ahead to “ABP footprint on "
               "actin”.")

    STEPS = {
        "1/10":  "PPI3D download (summary + PDB + data)",
        "2/10":  "Structure filtering + interface interactions",
        "3/10":  "MAFFT alignment per sequence cluster",
        "4/10":  "C70 interaction cluster analysis",
        "5/10":  "Interface B-factor computation (C70 + S1)",
        "6/10":  "PyMOL script generation per S1 cluster",
        "7/10":  "Heatmaps S1 binding site",
        "8/10":  "ABP analysis — competition and interfaces",
        "9/10":  "ABP structural analyses (Foldseek / InterPro / TM / footprint / chemistry)",
        "10/10": "Per-ABP ProteoCast (sequential submission to proteocast.ijm.fr)",
    }

    # Fichier de sortie attendu pour chaque étape (indicateur de complétion)
    STEP_OUTPUT_FILES = {
        "1/10":  "data/raw/all_data.csv",
        "2/10":  "data/filtered/details/1.interactions.csv",
        "3/10":  "data/alignments/.done",
        "4/10":  "data/filtered/patches_infos_cluster_data_70.csv",
        "5/10":  "data/filtered/details/structures_files/bfactor_cluster",
        "6/10":  "data/filtered/details/structures_files/bfactor_c70_interface/by_s1_gradient",
        "7/10":  "data/visualisations/actin_s1_all_equitable_heatmap.png",
        "8/10":  "data/visualisations/abp_analysis_done.flag",
        "9/10":  "data/exports/abp_site_domain/familles.csv",
        "10/10": "data/proteocast/abp_inputs/manifest.csv",
    }

    STEP_KEYS = list(STEPS.keys())
    TOTAL = len(STEPS)

    SKIP_KEYWORDS = ["Nothing to do", "unchanged",
                     "Using existing", "already up", "Déjà à jour",
                     "No changes detected"]


    def step_md(key, label, state):
        num = key.split("/")[0]
        note = {"running": " — *running…*", "done": " — *done*",
                "skipped": " — *up to date*", "error": " — **error**"}.get(state, "")
        return f"`{num:>2}` &nbsp;{label}{note}"


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


    # ── Section téléchargement / pipeline ─────────────────────────────────────────
    # Date à laquelle le jeu de données a été téléchargé depuis PPI3D.
    _alldata_p = "data/raw/all_data.csv"
    try:
        _dl_date = pd.Timestamp(os.path.getmtime(_alldata_p), unit="s").date()
        st.caption(f"PPI3D dataset downloaded on **{_dl_date}**.")
    except Exception:
        pass

    # Message affiché après un pipeline terminé (via st.rerun) — la date ci-dessus
    # reflète alors déjà les nouvelles données.
    if st.session_state.pop("_just_updated", False):
        st.success("Pipeline finished — data updated.")

    # Regroupement logique des étapes (affichage)
    _STEP_CATEGORIES = {
        "Data": ["1/10", "2/10"],
        "Processing & clustering": ["3/10", "4/10"],
        "Structures & PyMOL scripts": ["5/10", "6/10"],
        "Visualisation": ["7/10"],
        "ABP analyses": ["8/10", "9/10", "10/10"],
    }

    if "just_downloaded" not in st.session_state:
        st.session_state.just_downloaded = False

    _n_done = sum(1 for k in STEP_KEYS if initial_state(k) == "skipped")
    all_done = _n_done == TOTAL

    _c_btn, _c_prog = st.columns([1, 3], vertical_alignment="center")
    with _c_btn:
        clicked = st.button("Run / update", type="primary",
                            width="stretch")
    with _c_prog:
        _ptxt = (f"Pipeline up to date — {_n_done}/{TOTAL} steps"
                 if all_done else
                 f"{_n_done}/{TOTAL} steps up to date — click to complete")
        st.progress(_n_done / TOTAL, text=_ptxt)

    progress_bar = st.empty()

    # Étapes détaillées, groupées par catégorie, repliées si tout est à jour
    placeholders = {}
    with st.expander(f"Details of the {TOTAL} steps", expanded=not all_done):
        for _cat, _keys in _STEP_CATEGORIES.items():
            st.markdown(f"**{_cat}**")
            for key in _keys:
                placeholders[key] = st.empty()
                placeholders[key].markdown(
                    step_md(key, STEPS[key], initial_state(key)))

    sub_progress_bar = st.empty()
    sub_progress_text = st.empty()


    if clicked:
        st.session_state.just_downloaded = True

        for key, label in STEPS.items():
            init = initial_state(key)
            placeholders[key].markdown(
                step_md(key, label, init if init == "skipped" else "pending"))

        pb = progress_bar.progress(0, text="Starting...")

        current_idx = -1
        step_had_skip = {}
        log_lines = []

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
                    placeholders[key].markdown(step_md(key, STEPS[key], "running"))
                    pb.progress(i / TOTAL, text=f"Step {key} — {STEPS[key]}")
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

            if line:
                # conservé uniquement pour l'affichage en cas d'erreur
                log_lines.append(line)

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
            progress_bar.empty()
            # On relance le script pour que TOUT (date affichée, tableaux, vues) se
            # recharge avec les nouvelles données — sinon la date reste celle d'avant
            # le téléchargement (calculée en haut du script, avant ce bloc).
            st.session_state["_just_updated"] = True
            st.rerun()
        else:
            progress_bar.empty()
            st.error("Pipeline stopped on an error — last log lines:")
            st.code("\n".join(log_lines[-40:]), language=None)


    # ---------------------------------------------------------------------------
    # Section données filtrées
    # ---------------------------------------------------------------------------

