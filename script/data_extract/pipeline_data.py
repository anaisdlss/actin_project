import sys
import os
import subprocess
from pathlib import Path

# execution : caffeinate -i python -m script.data_extract.pipeline_data

PROJECT_ROOT = Path(__file__).resolve().parents[2]
FILTER_NOTEBOOK   = PROJECT_ROOT / "notebooks" / "graphe_filter.ipynb"
CLUSTER_NOTEBOOK  = PROJECT_ROOT / "notebooks" / "cluster_interaction_analysis.ipynb"
C70_NOTEBOOK      = PROJECT_ROOT / "notebooks" / "interface_analysis_c70.ipynb"
S1_NOTEBOOK       = PROJECT_ROOT / "notebooks" / "interface_analysis_s1.ipynb"
ABP_NOTEBOOK      = PROJECT_ROOT / "notebooks" / "abp_analysis.ipynb"

DATA          = PROJECT_ROOT / "data"
RAW           = DATA / "raw"
FILTERED      = DATA / "filtered"
DETAILS       = FILTERED / "details"
ALIGNMENTS    = DATA / "alignments"
VISUALISATIONS = PROJECT_ROOT / "visualisations"


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
            "1/15 — Téléchargement du summary des interactions actine (PPI3D BLAST)",
            [python_exec, "-m", "script.data_extract.get_summary_results"],
        )

        # 2 — Re-télécharge les entrées PDB si le summary a changé
        if is_up_to_date(RAW / "pdb_entry_results.csv", RAW / "ppi3d_actin_summary.csv"):
            skip_step("2/15 — Téléchargement des entrées PDB pour chaque structure")
        else:
            run_step(
                "2/15 — Téléchargement des entrées PDB pour chaque structure",
                [python_exec, "-m", "script.data_extract.get_pdb_entries"],
            )

        # 3 — Toujours exécuté : le script détecte lui-même si PPI3D a été mis à jour
        run_step(
            "3/15 — Téléchargement de toutes les données (cluster table PPI3D)",
            [python_exec, "-m", "script.data_extract.get_cluster_table"],
        )

        # 4 — Filtrage des structures (notebook)
        if is_up_to_date(
            FILTERED / "filtered_all_data.csv",
            RAW / "all_data.csv",
            RAW / "pdb_entry_results.csv",
            FILTER_NOTEBOOK,
        ):
            skip_step("4/15 — Filtrage des structures (≥ 4 actines connectées)")
        else:
            run_notebook(
                "4/15 — Filtrage des structures (≥ 4 actines connectées)",
                FILTER_NOTEBOOK,
            )
            (FILTERED / "filtered_all_data.csv").touch()

        # 5 — get_interaction_details gère lui-même le cache (hash du summary + progress.csv)
        run_step(
            "5/15 — Téléchargement des interactions d'interface",
            [python_exec, "-m", "script.data_extract.get_interaction_details"],
            input_text="f\n",
        )

        # 6 — Alignement MAFFT + enrichissement asa_pct dans table 4
        if is_up_to_date(
            ALIGNMENTS / ".done",
            FILTERED / "filtered_all_data.csv",
            FILTERED / "filtered_summary.csv",
            DETAILS / "1.interactions.csv",
            DETAILS / "3.interface_residues.csv",
            DETAILS / "4.inter-residue_contacts.csv",
        ):
            skip_step("6/15 — Alignement MAFFT par cluster de séquences")
        else:
            run_step(
                "6/15 — Alignement MAFFT par cluster de séquences",
                [python_exec, "-m", "script.mafft_pipeline"],
                cwd=PROJECT_ROOT,
            )

        # 7 — Analyse clusters d'interaction C70 (notebook)
        if is_up_to_date(
            FILTERED / "patches_infos_cluster_data_70.csv",
            DETAILS / "1.interactions.csv",
            FILTERED / "filtered_all_data.csv",
            CLUSTER_NOTEBOOK,
        ):
            skip_step("7/15 — Analyse des clusters d'interaction C70")
        else:
            run_notebook(
                "7/15 — Analyse des clusters d'interaction C70",
                CLUSTER_NOTEBOOK,
            )
            (FILTERED / "patches_infos_cluster_data_70.csv").touch()

        # 8 — Calcul B-factors interface C70 (doit précéder le notebook C70)
        if is_up_to_date(
            FILTERED / "details" / "structures_files" / "bfactor_c70_interface",
            FILTERED / "patches_infos_cluster_data_70.csv",
            DETAILS / "4.inter-residue_contacts.csv",
            DETAILS / "8.structures.csv",
        ):
            skip_step("8/15 — Calcul B-factors interface C70 par cluster")
        else:
            run_step(
                "8/15 — Calcul B-factors interface C70 par cluster (bfactor_c70_interface.py)",
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
            skip_step("9/15 — Génération script PyMOL surface complète C70")
        else:
            run_step(
                "9/15 — Génération script PyMOL surface complète C70 (bfactor_c70_pymol_full_surface.py)",
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
            skip_step("10/15 — Génération scripts PyMOL par site S1")
        else:
            run_step(
                "10/15 — Génération scripts PyMOL par site S1 (bfactor_c70_pymol_by_s1.py)",
                [python_exec, "-m", "script.bfactor_c70_pymol_by_s1"],
                cwd=PROJECT_ROOT,
            )

        # 11 — Analyse interface C70 détaillée (notebook)
        if is_up_to_date(
            VISUALISATIONS / "actin_c70_contacts",
            FILTERED / "patches_infos_cluster_data_70.csv",
            DETAILS / "3.interface_residues.csv",
            DETAILS / "4.inter-residue_contacts.csv",
            C70_NOTEBOOK,
        ):
            skip_step("11/15 — Analyse interface par cluster C70")
        else:
            run_notebook(
                "11/15 — Analyse interface par cluster C70",
                C70_NOTEBOOK,
            )
            (VISUALISATIONS / "actin_c70_contacts").touch()

        # 12 — Heatmap S1 binding site + mise à jour références (notebook)
        if is_up_to_date(
            VISUALISATIONS / "actin_s1_all_equitable_heatmap.png",
            FILTERED / "patches_infos_cluster_data_70.csv",
            DETAILS / "3.interface_residues.csv",
            DETAILS / "4.inter-residue_contacts.csv",
            S1_NOTEBOOK,
        ):
            skip_step("12/15 — Heatmap S1 binding site et références clusters")
        else:
            run_notebook(
                "12/15 — Heatmap S1 binding site et références clusters",
                S1_NOTEBOOK,
            )
            (VISUALISATIONS / "actin_s1_all_equitable_heatmap.png").touch()

        # 13 — Calcul B-factors S1 par cluster (pour PyMOL / Streamlit)
        if is_up_to_date(
            FILTERED / "details" / "structures_files" / "bfactor_cluster",
            FILTERED / "s1_cluster_reference.csv",
            FILTERED / "patches_infos_cluster_data_70.csv",
            DETAILS / "3.interface_residues.csv",
        ):
            skip_step("13/15 — Calcul B-factors S1 par cluster")
        else:
            run_step(
                "13/15 — Calcul B-factors S1 par cluster (bfactor.py)",
                [python_exec, "-m", "script.bfactor"],
                cwd=PROJECT_ROOT,
            )

        # 14 — Analyse ABP : compétition, interfaces, PDB sans ABP (notebook)
        _flag = VISUALISATIONS / "abp_analysis_done.flag"
        _flag.parent.mkdir(parents=True, exist_ok=True)
        if is_up_to_date(
            _flag,
            FILTERED / "filtered_all_data.csv",
            FILTERED / "patches_infos_cluster_data_70.csv",
            FILTERED / "patches_infos_s1_binding_site.csv",
            ABP_NOTEBOOK,
        ):
            skip_step("14/15 — Analyse ABP — compétition et interfaces")
        else:
            run_notebook(
                "14/15 — Analyse ABP — compétition et interfaces",
                ABP_NOTEBOOK,
            )
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
