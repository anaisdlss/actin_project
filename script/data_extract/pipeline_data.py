import sys
import os
import subprocess
from pathlib import Path
import pandas as pd

# execution : caffeinate -i python -m script.data_extract.pipeline_data

PROJECT_ROOT = Path(__file__).resolve().parents[2]
FILTER_NOTEBOOK   = PROJECT_ROOT / "notebooks" / "graphe_filter.ipynb"
CLUSTER_NOTEBOOK  = PROJECT_ROOT / "notebooks" / "cluster_interaction_analysis.ipynb"
C70_NOTEBOOK      = PROJECT_ROOT / "notebooks" / "interface_analysis_c70.ipynb"
S1_NOTEBOOK       = PROJECT_ROOT / "notebooks" / "interface_analysis_s1.ipynb"
ABP_NOTEBOOK      = PROJECT_ROOT / "notebooks" / "abp_analysis.ipynb"

DATA      = PROJECT_ROOT / "data"
RAW       = DATA / "raw"
FILTERED  = DATA / "filtered"
DETAILS   = FILTERED / "details"
ALIGNMENTS = DATA / "alignments"


def is_up_to_date(output: Path, *inputs: Path) -> bool:
    """True si output existe et est plus récent que tous les inputs."""
    if not output.exists():
        return False
    out_mtime = output.stat().st_mtime
    return all(inp.exists() and out_mtime >= inp.stat().st_mtime for inp in inputs)


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


def cleanup_orphan_chains(filtered_dir: Path, details_dir: Path):
    """Supprime de 3.interface_residues.csv les chaînes absentes de filtered_all_data
    (pas de séquence connue → impossible d'aligner avec MAFFT)."""
    df_all = pd.read_csv(filtered_dir / "filtered_all_data.csv")
    valid_chains = set(df_all["subunit_1"]) | set(df_all["subunit_2"])

    df3_path = details_dir / "3.interface_residues.csv"
    df3 = pd.read_csv(df3_path)
    before = len(df3)
    df3_clean = df3[df3["chain"].isin(valid_chains)]
    removed = before - len(df3_clean)
    if removed > 0:
        df3_clean.to_csv(df3_path, index=False)
        print(f"  Chaînes orphelines supprimées : {removed} résidus retirés de df3")
    else:
        print("  Aucune chaîne orpheline détectée dans df3")


def run_notebook(label, notebook_path):
    run_step(label, [
        sys.executable, "-m", "jupyter", "nbconvert",
        "--to", "notebook", "--execute", "--inplace",
        str(notebook_path),
    ], cwd=PROJECT_ROOT)


def main():
    python_exec = sys.executable
    os.environ["MPLBACKEND"] = "Agg"

    for d in [RAW, FILTERED, DETAILS, ALIGNMENTS]:
        d.mkdir(parents=True, exist_ok=True)

    try:
        # 1 — Toujours exécuté : le script détecte lui-même si PPI3D a été mis à jour
        run_step(
            "1/14 — Téléchargement du summary des interactions actine (PPI3D BLAST)",
            [python_exec, "-m", "script.data_extract.get_summary_results"],
        )

        # 2 — Re-télécharge les entrées PDB si le summary a changé
        if is_up_to_date(RAW / "pdb_entry_results.csv", RAW / "ppi3d_actin_summary.csv"):
            skip_step("2/14 — Téléchargement des entrées PDB pour chaque structure")
        else:
            run_step(
                "2/14 — Téléchargement des entrées PDB pour chaque structure",
                [python_exec, "-m", "script.data_extract.get_pdb_entries"],
            )

        # 3 — Toujours exécuté : le script détecte lui-même si PPI3D a été mis à jour
        run_step(
            "3/14 — Téléchargement de toutes les données (cluster table PPI3D)",
            [python_exec, "-m", "script.data_extract.get_cluster_table"],
        )

        # 4 — Filtrage des structures (notebook)
        run_notebook(
            "4/14 — Filtrage des structures (≥ 4 actines connectées)",
            FILTER_NOTEBOOK,
        )

        # 5 — get_interaction_details gère lui-même le cache (hash du summary + progress.csv)
        run_step(
            "5/14 — Téléchargement des interactions d'interface",
            [python_exec, "-m", "script.data_extract.get_interaction_details"],
            input_text="f\n",
        )

        # 5.5 — Nettoyage des chaînes orphelines (absentes de filtered_all_data)
        print("\n" + "=" * 60)
        print("ETAPE : 5.5/14 — Nettoyage des chaînes orphelines")
        print("=" * 60)
        cleanup_orphan_chains(FILTERED, DETAILS)

        # 5.7 — Enrichissement table 4 avec % ASA buried
        # (cellules 15-16 de graphe_filter.ipynb, différées car elles lisent
        #  3.interface_residues.csv qui n'existe qu'après l'étape 5)
        run_notebook(
            "5.7/14 — Enrichissement table 4 avec % ASA buried (graphe_filter.ipynb — phase 2)",
            FILTER_NOTEBOOK,
        )

        # 6 — Alignement MAFFT
        run_step(
            "6/14 — Alignement MAFFT par cluster de séquences",
            [python_exec, "-m", "script.mafft_pipeline"],
            cwd=PROJECT_ROOT,
        )

        # 7 — Analyse clusters d'interaction C70 (notebook)
        run_notebook(
            "7/14 — Analyse des clusters d'interaction C70",
            CLUSTER_NOTEBOOK,
        )

        # 8 — Calcul B-factors interface C70 (doit précéder le notebook C70)
        if is_up_to_date(
            FILTERED / "details" / "structures_files" / "bfactor_c70_interface",
            FILTERED / "patches_infos_cluster_data_70.csv",
            DETAILS / "4.inter-residue_contacts.csv",
            DETAILS / "8.structures.csv",
        ):
            skip_step("8/14 — Calcul B-factors interface C70 par cluster")
        else:
            run_step(
                "8/14 — Calcul B-factors interface C70 par cluster (bfactor_c70_interface.py)",
                [python_exec, "-m", "script.bfactor_c70_interface"],
                cwd=PROJECT_ROOT,
            )

        # 9 — Génération script PyMOL surface complète (dépend des PDB bfactor C70)
        _pml_full = FILTERED / "details" / "structures_files" / "bfactor_c70_interface" / "view_full_surface.pml"
        if is_up_to_date(
            _pml_full,
            FILTERED / "patches_infos_cluster_data_70.csv",
            FILTERED / "details" / "structures_files" / "bfactor_c70_interface",
        ):
            skip_step("9/14 — Génération script PyMOL surface complète C70")
        else:
            run_step(
                "9/14 — Génération script PyMOL surface complète C70 (bfactor_c70_pymol_full_surface.py)",
                [python_exec, "-m", "script.bfactor_c70_pymol_full_surface"],
                cwd=PROJECT_ROOT,
            )

        # 10 — Génération scripts PyMOL par site S1 (dépend des PDB bfactor C70)
        _pml_by_s1_dir = FILTERED / "details" / "structures_files" / "bfactor_c70_interface" / "by_s1_cluster"
        if is_up_to_date(
            _pml_by_s1_dir,
            FILTERED / "patches_infos_cluster_data_70.csv",
            FILTERED / "filtered_all_data.csv",
            FILTERED / "details" / "structures_files" / "bfactor_c70_interface",
        ):
            skip_step("10/14 — Génération scripts PyMOL par site S1")
        else:
            run_step(
                "10/14 — Génération scripts PyMOL par site S1 (bfactor_c70_pymol_by_s1.py)",
                [python_exec, "-m", "script.bfactor_c70_pymol_by_s1"],
                cwd=PROJECT_ROOT,
            )

        # 11 — Analyse interface C70 détaillée (notebook)
        run_notebook(
            "11/14 — Analyse interface par cluster C70",
            C70_NOTEBOOK,
        )

        # 12 — Heatmap S1 binding site + mise à jour références (notebook)
        run_notebook(
            "12/14 — Heatmap S1 binding site et références clusters",
            S1_NOTEBOOK,
        )

        # 13 — Calcul B-factors S1 par cluster (pour PyMOL / Streamlit)
        if is_up_to_date(
            FILTERED / "details" / "structures_files" / "bfactor_cluster",
            FILTERED / "s1_cluster_reference.csv",
            FILTERED / "patches_infos_cluster_data_70.csv",
            DETAILS / "3.interface_residues.csv",
        ):
            skip_step("13/14 — Calcul B-factors S1 par cluster")
        else:
            run_step(
                "13/14 — Calcul B-factors S1 par cluster (bfactor.py)",
                [python_exec, "-m", "script.bfactor"],
                cwd=PROJECT_ROOT,
            )

        # 14 — Analyse ABP : compétition, interfaces, PDB sans ABP (notebook)
        run_notebook(
            "14/14 — Analyse ABP — compétition et interfaces",
            ABP_NOTEBOOK,
        )
        # Flag de sortie pour le suivi Streamlit
        _flag = PROJECT_ROOT / "visualisations" / "abp_analysis_done.flag"
        _flag.parent.mkdir(parents=True, exist_ok=True)
        _flag.touch()

        # 15 — Sessions PyMOL filament par ABP (dépend bfactor_cluster + filtered_all_data)
        _by_abp_dir = FILTERED / "details" / "structures_files" / "filament" / "by_abp"
        if is_up_to_date(
            _by_abp_dir,
            FILTERED / "patches_infos_cluster_data_70.csv",
            FILTERED / "filtered_all_data.csv",
            FILTERED / "details" / "structures_files" / "bfactor_cluster",
        ):
            skip_step("15/15 — Sessions PyMOL filament par ABP")
        else:
            run_step(
                "15/15 — Sessions PyMOL filament par ABP (filament_by_abp.py)",
                [python_exec, "-m", "script.filament_by_abp"],
                cwd=PROJECT_ROOT,
            )

        print("\nPipeline terminé avec succès.")

    except subprocess.CalledProcessError as e:
        print("\nPipeline interrompu : une étape a échoué.")
        print("Commande :", e.cmd)
        print("Code retour :", e.returncode)
        sys.exit(e.returncode)


if __name__ == "__main__":
    main()
