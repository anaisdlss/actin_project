import sys
import os
import subprocess
from pathlib import Path

# execution : caffeinate -i python -m script.data_extract.pipeline_data

PROJECT_ROOT = Path(__file__).resolve().parents[2]
# Étapes d'analyse : anciens notebooks convertis en scripts .py (lancés
# directement par le pipeline, plus légers et sans dépendance jupyter au runtime).
FILTER_NOTEBOOK   = PROJECT_ROOT / "notebooks" / "graphe_filter.py"
SUPERCLUSTER_NOTEBOOK = PROJECT_ROOT / "notebooks" / "binding_site_superclusters.py"
CLUSTER_NOTEBOOK  = PROJECT_ROOT / "notebooks" / "cluster_interaction_analysis.py"
C70_NOTEBOOK      = PROJECT_ROOT / "notebooks" / "interface_analysis_c70.py"
S1_NOTEBOOK       = PROJECT_ROOT / "notebooks" / "interface_analysis_s1.py"
ABP_NOTEBOOK      = PROJECT_ROOT / "notebooks" / "abp_analysis.py"

DATA          = PROJECT_ROOT / "data"
RAW           = DATA / "raw"
FILTERED      = DATA / "filtered"
DETAILS       = FILTERED / "details"
ALIGNMENTS    = DATA / "alignments"
VISUALISATIONS = PROJECT_ROOT / "data" / "visualisations"


def _ppi3d_update() -> str | None:
    """Date de dernière mise à jour PPI3D stockée dans data/raw/metadata.json."""
    import json
    try:
        return json.loads((RAW / "metadata.json").read_text()).get("ppi3d_last_update")
    except Exception:
        return None


def _has_supercluster_cols() -> bool:
    """filtered_all_data contient-il déjà les colonnes de super-clusters ?"""
    try:
        import pandas as pd
        cols = pd.read_csv(FILTERED / "filtered_all_data.csv", nrows=0).columns
        return "s1_supercluster" in cols
    except Exception:
        return False


def is_up_to_date(output: Path, *inputs: Path) -> bool:
    """True si output existe et est plus récent que tous les inputs."""
    if not output.exists():
        return False
    out_mtime = output.stat().st_mtime
    return all(inp.exists() and out_mtime >= inp.stat().st_mtime for inp in inputs)


def _details_complete() -> bool:
    """True si 1.interactions.csv couvre TOUTES les interactions du filtered_summary.

    Garde-fou contre un téléchargement des détails interrompu/incomplet : dans ce
    cas on ne doit PAS considérer le pipeline comme « à jour ». Renvoie False (=
    incomplet, on relance) au moindre doute (fichier manquant, erreur de lecture)."""
    try:
        import pandas as pd
        summary = FILTERED / "filtered_summary.csv"
        details = DETAILS / "1.interactions.csv"
        if not summary.exists() or not details.exists():
            return False
        sids = set(pd.read_csv(summary)["interaction_id"]
                   .dropna().astype(int).astype(str))
        dids = set(pd.read_csv(details, dtype={"interaction_id": str})
                   ["interaction_id"].dropna().astype(str))
        # Interactions définitivement introuvables sur PPI3D (après plusieurs
        # tentatives) : consignées dans failed.csv. On ne les compte pas comme
        # bloquantes, sinon le pipeline ne se terminerait jamais.
        fids = set()
        failed = DETAILS / "failed.csv"
        if failed.exists():
            fdf = pd.read_csv(failed, dtype={"interaction_id": str})
            if "interaction_id" in fdf.columns:
                fids = set(fdf["interaction_id"].dropna().astype(str))
        return bool(sids) and sids.issubset(dids | fids)
    except Exception:
        return False



def run_step(label, command, input_text=None, cwd=None):
    print("\n" + "=" * 60)
    print(f"ETAPE : {label}")
    print("=" * 60)
    subprocess.run(command, input=input_text, text=True, check=True, cwd=cwd)


def skip_step(label):
    print("\n" + "=" * 60)
    print(f"ETAPE : {label}")
    print("=" * 60)
    print("Already up to date — Nothing to do")


def run_notebook(label, notebook_path):
    run_step(label, [sys.executable, str(notebook_path)], cwd=PROJECT_ROOT)


def _exec(command, input_text=None):
    subprocess.run(command, input=input_text, text=True, check=True, cwd=PROJECT_ROOT)


def _nb(notebook):
    # Les « notebooks » sont désormais des scripts .py : on les exécute
    # directement (cwd = racine projet, donc les chemins relatifs data/… marchent).
    _exec([sys.executable, str(notebook)])


def run_group(key_title, substeps):
    """Une étape groupée = 1 marqueur ETAPE + plusieurs sous-étapes (cache par sous-étape).
    substeps : liste de (sous-label, deja_a_jour: bool, run: callable, touch: Path|None)."""
    print(f"ETAPE : {key_title}")
    if all(utd for _, utd, _, _ in substeps):
        print("Déjà à jour")
        return
    for sublabel, utd, run, touch in substeps:
        if utd:
            print(f"  - {sublabel} : inchangé")
        else:
            print(f"  - {sublabel}…")
            run()
            if touch is not None:
                touch.touch()


def main():
    py = sys.executable
    os.environ["MPLBACKEND"] = "Agg"
    SF = FILTERED / "details" / "structures_files"
    RG = PROJECT_ROOT / "script" / "data_extract"

    for d in [RAW, FILTERED, DETAILS, ALIGNMENTS]:
        d.mkdir(parents=True, exist_ok=True)

    try:
        # ══ 1/9 — Téléchargement PPI3D (summary + entrées PDB + données) ══════
        # On lance d'abord le summary : il interroge PPI3D et met à jour la date
        # de dernière mise à jour. Si cette date est INCHANGÉE et que le pipeline
        # est déjà complet, il n'y a rien à refaire → on arrête tout ici.
        print("ETAPE : 1/9 — Téléchargement PPI3D (summary + PDB + données)")
        _prev_update = _ppi3d_update()
        print("  - Summary PPI3D (BLAST)…")
        _exec([py, "-m", "script.data_extract.get_summary_results"])
        _now_update = _ppi3d_update()
        _last_output = PROJECT_ROOT / "data/exports/abp_site_domain/familles.csv"
        # On n'arrête tôt QUE si : PPI3D inchangé, ET le jeu de données complet
        # (all_data.csv) est cohérent avec le summary courant (pas plus vieux),
        # ET la dernière sortie du pipeline existe. Ainsi un summary rafraîchi
        # sans re-téléchargement du jeu ne verrouille PAS des données périmées.
        _data_coherent = is_up_to_date(RAW / "all_data.csv", RAW / "ppi3d_actin_summary.csv")
        # Robustesse : on n'arrête tôt que si les détails d'interface sont COMPLETS
        # (pas de download interrompu) ET que la sortie finale (familles.csv) est
        # plus récente que les détails/données filtrées (sinon les analyses en aval
        # sont périmées et doivent se recalculer). Ainsi un utilisateur ne peut pas
        # rester bloqué avec un jeu silencieusement incomplet.
        _details_ok = _details_complete()
        _downstream_fresh = is_up_to_date(
            _last_output, DETAILS / "1.interactions.csv",
            FILTERED / "filtered_all_data.csv")
        if (_prev_update and _prev_update not in ("unknown", None)
                and _prev_update == _now_update and _last_output.exists()
                and _data_coherent and _details_ok and _downstream_fresh):
            print(f"  PPI3D inchangé (dernière mise à jour : {_now_update}), "
                  "jeu de données cohérent, détails complets et analyses à jour.")
            print("Aucune nouvelle donnée — rien à refaire.")
            print("\nPipeline terminé avec succès.")
            return
        if not _details_ok:
            print("  Détails d'interface INCOMPLETS (download précédent interrompu ?) "
                  "→ le pipeline va les compléter.")
        elif not _downstream_fresh:
            print("  Détails mis à jour depuis la dernière analyse "
                  "→ recalcul des étapes en aval.")
        print("  - Entrées PDB…")
        if is_up_to_date(RAW / "pdb_entry_results.csv", RAW / "ppi3d_actin_summary.csv"):
            print("    inchangé")
        else:
            _exec([py, "-m", "script.data_extract.get_pdb_entries"])
        print("  - Toutes les données (cluster table)…")
        _exec([py, "-m", "script.data_extract.get_cluster_table"])

        # ══ 2/9 — Filtrage + interactions d'interface ════════════════════════
        run_group("2/9 — Filtrage des structures + interactions d'interface", [
            ("Filtrage (≥ 4 actines) [notebook]",
             is_up_to_date(FILTERED / "filtered_all_data.csv", RAW / "all_data.csv",
                           RAW / "pdb_entry_results.csv", FILTER_NOTEBOOK),
             lambda: _nb(FILTER_NOTEBOOK), FILTERED / "filtered_all_data.csv"),
            ("Interactions d'interface", False,
             lambda: _exec([py, "-m", "script.data_extract.get_interaction_details"],
                           input_text="f\n"), None),
        ])

        # ══ 3/9 — Alignement MAFFT ═══════════════════════════════════════════
        run_group("3/9 — Alignement MAFFT par cluster de séquences", [
            ("MAFFT + asa_pct",
             is_up_to_date(ALIGNMENTS / ".done", FILTERED / "filtered_all_data.csv",
                           FILTERED / "filtered_summary.csv", DETAILS / "1.interactions.csv",
                           DETAILS / "3.interface_residues.csv", DETAILS / "4.inter-residue_contacts.csv"),
             lambda: _exec([py, "-m", "script.mafft_pipeline"]), None),
        ])

        # ══ 4/9 — Super-clusters + analyse clusters d'interaction C70 ════════
        _sc_flag = FILTERED / ".superclusters.done"
        run_group("4/9 — Analyse des clusters d'interaction C70", [
            ("Super-clusters de sites (ajout colonnes s1/s2_supercluster) [notebook]",
             _has_supercluster_cols() and is_up_to_date(
                 _sc_flag, ALIGNMENTS / ".done", SUPERCLUSTER_NOTEBOOK),
             lambda: _nb(SUPERCLUSTER_NOTEBOOK), _sc_flag),
            ("Clusters C70 [notebook]",
             is_up_to_date(FILTERED / "patches_infos_cluster_data_70.csv",
                           DETAILS / "1.interactions.csv", FILTERED / "filtered_all_data.csv", CLUSTER_NOTEBOOK),
             lambda: _nb(CLUSTER_NOTEBOOK), FILTERED / "patches_infos_cluster_data_70.csv"),
            ("Interface par cluster C70 [notebook]",
             is_up_to_date(VISUALISATIONS / "actin_c70_contacts", FILTERED / "patches_infos_cluster_data_70.csv",
                           DETAILS / "3.interface_residues.csv", DETAILS / "4.inter-residue_contacts.csv", C70_NOTEBOOK),
             lambda: _nb(C70_NOTEBOOK), VISUALISATIONS / "actin_c70_contacts"),
        ])

        # ══ 5/9 — B-factors interface (C70 + S1) ═════════════════════════════
        run_group("5/9 — Calcul des B-factors interface (C70 + S1)", [
            ("B-factors interface C70",
             is_up_to_date(SF / "bfactor_c70_interface", FILTERED / "patches_infos_cluster_data_70.csv",
                           DETAILS / "4.inter-residue_contacts.csv", DETAILS / "8.structures.csv"),
             lambda: _exec([py, "-m", "script.bfactor_c70_interface"]), None),
            ("B-factors S1 par cluster",
             is_up_to_date(SF / "bfactor_cluster", FILTERED / "s1_cluster_reference.csv",
                           FILTERED / "patches_infos_cluster_data_70.csv", DETAILS / "3.interface_residues.csv"),
             lambda: _exec([py, "-m", "script.bfactor"]), None),
        ])

        # ══ 6/9 — Scripts PyMOL par cluster S1 (gradient) ════════════════════
        run_group("6/9 — Génération des scripts PyMOL par cluster S1", [
            ("Scripts PyMOL par site S1",
             is_up_to_date(SF / "bfactor_c70_interface" / "by_s1_cluster",
                           FILTERED / "patches_infos_cluster_data_70.csv",
                           FILTERED / "filtered_all_data.csv", SF / "bfactor_c70_interface"),
             lambda: _exec([py, "-m", "script.bfactor_c70_pymol_by_s1"]), None),
            ("Scripts PyMOL gradient (actine opaque)",
             is_up_to_date(SF / "bfactor_c70_interface" / "by_s1_gradient",
                           FILTERED / "patches_infos_cluster_data_70.csv", SF / "bfactor_c70_interface"),
             lambda: _exec([py, "-m", "script.bfactor_c70_pymol_by_s1_gradient"]), None),
        ])

        # ══ 7/9 — Heatmaps S1 binding site ═══════════════════════════════════
        run_group("7/9 — Heatmaps S1 binding site", [
            ("Profils + CSV équitable C70 [notebook]",
             is_up_to_date(FILTERED / "actin_s1_canon_area_by_cluster.csv",
                           FILTERED / "patches_infos_cluster_data_70.csv",
                           DETAILS / "3.interface_residues.csv", DETAILS / "4.inter-residue_contacts.csv", S1_NOTEBOOK),
             lambda: _nb(S1_NOTEBOOK), FILTERED / "actin_s1_canon_area_by_cluster.csv"),
            ("Heatmaps par cluster + globale (labels lisibles)",
             is_up_to_date(VISUALISATIONS / "actin_s1_all_equitable_heatmap.png",
                           FILTERED / "actin_s1_canon_area_by_cluster.csv"),
             lambda: (_exec([py, str(RG / "regenerate_s1_heatmaps.py")]),
                      _exec([py, str(RG / "regenerate_s1_global_heatmap.py")])), None),
        ])

        # ══ 8/9 — Analyse ABP (compétition + interfaces) [notebook] ══════════
        _flag = VISUALISATIONS / "abp_analysis_done.flag"
        _flag.parent.mkdir(parents=True, exist_ok=True)
        run_group("8/9 — Analyse ABP — compétition et interfaces", [
            ("Compétition, interfaces, PDB sans ABP [notebook]",
             is_up_to_date(_flag, FILTERED / "filtered_all_data.csv",
                           FILTERED / "patches_infos_cluster_data_70.csv",
                           FILTERED / "patches_infos_s1_binding_site.csv", ABP_NOTEBOOK),
             lambda: _nb(ABP_NOTEBOOK), _flag),
        ])

        # ══ 9/9 — Analyses structurales ABP (Foldseek / InterPro / TM / …) ════
        run_group("9/9 — Analyses structurales ABP (convergence)", [
            ("Foldseek + InterPro + TM + empreinte + SS + chimie",
             is_up_to_date(PROJECT_ROOT / "data/exports/abp_site_domain/familles.csv",
                           FILTERED / "filtered_all_data.csv", DETAILS / "3.interface_residues.csv",
                           DETAILS / "1.interactions.csv"),
             lambda: _exec([py, "-m", "script.abp_site_domain.run_all"]), None),
        ])

        print("\nPipeline terminé avec succès.")

    except subprocess.CalledProcessError as e:
        print("\nPipeline interrompu : une étape a échoué.")
        print("Commande :", e.cmd)
        print("Code retour :", e.returncode)
        sys.exit(e.returncode)



if __name__ == "__main__":
    main()
