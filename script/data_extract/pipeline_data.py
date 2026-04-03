import sys
import os
import subprocess
from pathlib import Path

# execution : caffeinate -i python script/data_extract/pipeline_data.py

PROJECT_ROOT = Path(__file__).resolve().parents[2]
FILTER_NOTEBOOK = PROJECT_ROOT / "notebooks" / "graphe_filter.ipynb"


def run_step(label, command, input_text=None, cwd=None):
    print("\n" + "=" * 60)
    print(f"ETAPE : {label}")
    print("=" * 60)

    result = subprocess.run(
        command,
        input=input_text,
        text=True,
        check=True,
        cwd=cwd
    )
    return result


def main():
    python_exec = sys.executable
    os.environ["MPLBACKEND"] = "Agg"

    try:
        # 1. summary principal
        run_step(
            "1/6 — Téléchargement du summary des interactions actine (PPI3D BLAST)",
            [python_exec, "-m", "script.data_extract.get_summary_results"]
        )

        # 2. pdb entries
        run_step(
            "2/6 — Téléchargement des entrées PDB pour chaque structure",
            [python_exec, "-m", "script.data_extract.get_pdb_entries"]
        )

        # 3. table clusters globale PPI3D
        run_step(
            "3/6 — Téléchargement de la table globale des clusters PPI3D",
            [python_exec, "-m", "script.data_extract.get_cluster_table"]
        )

        # 4. notebook de filtrage
        summary_path = PROJECT_ROOT / "data" / "raw" / "ppi3d_actin_summary.csv"
        filtered_path = PROJECT_ROOT / "data" / "filtered" / "filtered_pdb_entry.csv"

        if (
            filtered_path.exists()
            and summary_path.exists()
            and filtered_path.stat().st_mtime >= summary_path.stat().st_mtime
        ):
            print("\n" + "=" * 60)
            print("ETAPE : 4/6 — Filtrage des structures (garder celles avec >= 5 actines)")
            print("=" * 60)
            print("Filtered data already up to date — Nothing to do")
        else:
            run_step(
                "4/6 — Filtrage des structures (garder celles avec >= 5 actines)",
                [
                    python_exec, "-m", "jupyter", "nbconvert",
                    "--to", "notebook", "--execute", "--inplace",
                    str(FILTER_NOTEBOOK),
                ],
                cwd=PROJECT_ROOT
            )

        # 5. details depuis le summary filtré
        run_step(
            "5/6 — Téléchargement des détails d'interface (résidus de contact)",
            [python_exec, "-m", "script.data_extract.get_interaction_details"],
            input_text="f\n"
        )

        # 6. filtrer all_data
        all_data_path = PROJECT_ROOT / "data" / "raw" / "all_data.csv"
        filtered_all_data_path = PROJECT_ROOT / "data" / "filtered" / "filtered_all_data.csv"

        if (
            filtered_all_data_path.exists()
            and all_data_path.exists()
            and filtered_path.exists()
            and filtered_all_data_path.stat().st_mtime >= all_data_path.stat().st_mtime
            and filtered_all_data_path.stat().st_mtime >= filtered_path.stat().st_mtime
        ):
            print("\n" + "=" * 60)
            print("ETAPE : 6/6 — Filtrage de la table all_data selon les PDB retenus")
            print("=" * 60)
            print("Filtered all_data already up to date — Nothing to do")
        else:
            run_step(
                "6/6 — Filtrage de la table all_data selon les PDB retenus",
                [python_exec, "-m", "script.filter.filter_all_data_by_filtered_pdb"]
            )

        print("\nPipeline finished successfully.")

    except subprocess.CalledProcessError as e:
        print("\nPipeline stopped because one step failed.")
        print("Failed command:", e.cmd)
        print("Return code:", e.returncode)

        sys.exit(e.returncode)


if __name__ == "__main__":
    main()
