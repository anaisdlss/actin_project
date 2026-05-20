# filament_by_abp.py
# Pour chaque ABP, génère une session PyMOL avec 2 filaments actine 10 sous-unités :
#   filament_base : top 2 clusters C70 homo globaux (0_7797_0 + 0_7797_1)
#   filament_abp  : top 2 clusters C70 homo spécifiques à l'ABP
# B-factors = somme des % ASA buried des sites de liaison S1 des clusters C70
#
# Exécution : python -m script.filament_by_abp

import ast
import re
import os
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FILTERED = PROJECT_ROOT / "data" / "filtered"
BFACTOR_DIR = FILTERED / "details" / "structures_files" / "bfactor_cluster"
TEMPLATE_PDB = FILTERED / "details" / "structures_files" / "filament" / "actin_filament_10.pdb"
OUT_DIR = FILTERED / "details" / "structures_files" / "filament" / "by_abp"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GLOBAL_BASE_PDB = FILTERED / "details" / "structures_files" / "filament" / "filament_global_base.pdb"

GLOBAL_PATCHES = ["0_7797_0", "0_7797_1"]


# ── helpers ──────────────────────────────────────────────────────────────────

def safe_name(title: str) -> str:
    return re.sub(r"[^a-zA-Z0-9]+", "_", title).strip("_")[:60]


def read_bfactors_sum(s1_sites: list[str]) -> dict[int, float]:
    """Somme les B-factors par numéro de résidu sur tous les fichiers S1."""
    totals: dict[int, float] = {}
    for site in s1_sites:
        path = BFACTOR_DIR / f"{site}.pdb"
        if not path.exists():
            print(f"  [warn] fichier manquant : {path.name}")
            continue
        with open(path) as fh:
            for line in fh:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                resi = int(line[22:26])
                b = float(line[60:66])
                totals[resi] = totals.get(resi, 0.0) + b
    return totals


def apply_bfactors(template_path: Path, bfactors: dict[int, float], out_path: Path) -> float:
    """Recopie le template PDB en remplaçant les B-factors ; retourne le max."""
    out_lines = []
    max_b = 0.0
    with open(template_path) as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resi = int(line[22:26])
                b = bfactors.get(resi, 0.0)
                if b > max_b:
                    max_b = b
                line = line[:60] + f"{b:6.2f}" + line[66:]
            out_lines.append(line)
    with open(out_path, "w") as fh:
        fh.writelines(out_lines)
    return max_b


def get_s1_sites(patches: list[str], df_patches: pd.DataFrame) -> list[str]:
    """Retourne la liste dédupliquée des sites S1 pour les patches donnés."""
    sites = []
    seen = set()
    for patch in patches:
        row = df_patches[df_patches["patch"] == patch]
        if row.empty:
            print(f"  [warn] patch inconnu : {patch}")
            continue
        bsc = ast.literal_eval(row.iloc[0]["binding_site_clusters"])
        for s in bsc:
            if s not in seen:
                sites.append(s)
                seen.add(s)
    return sites


def top2_patches(abp_pdbs: set, df_homo: pd.DataFrame) -> list[str]:
    """Top 2 patches C70 homo par nombre de PDBs distincts dans l'ensemble de l'ABP."""
    sub = df_homo[df_homo["pdb_id"].isin(abp_pdbs)]
    if sub.empty:
        return []
    coverage = (
        sub.groupby("cluster_data_70")["pdb_id"]
        .nunique()
        .sort_values(ascending=False)
    )
    return list(coverage.index[:2])


def write_pml(pml_path: Path, base_pdb: Path, abp_pdb: Path,
              max_base: float, max_abp: float, abp_title: str,
              base_patches: list[str], abp_patches: list[str]) -> None:
    base_rel = base_pdb.as_posix()
    abp_rel = abp_pdb.as_posix()
    base_patch_str = " + ".join(base_patches)
    abp_patch_str = " + ".join(abp_patches) if abp_patches else "aucun"

    content = f"""\
# PyMOL — Filament actine 10 s-u — {abp_title}
# filament_base : clusters globaux ({base_patch_str}), white_red, max={max_base:.2f}
# filament_abp  : clusters spécifiques ({abp_patch_str}), white_orange, max={max_abp:.2f}

# ── filament_base ──
load {base_rel}, filament_base
hide everything, filament_base
show surface, filament_base
color white, filament_base
spectrum b, white_red, filament_base, minimum=0, maximum={max_base:.2f}

# ── filament_abp ──
load {abp_rel}, filament_abp
hide everything, filament_abp
show surface, filament_abp
color white, filament_abp
align filament_abp, filament_base
spectrum b, white_orange, filament_abp, minimum=0, maximum={max_abp:.2f}

set surface_quality, 1
bg_color white
zoom filament_base
"""
    with open(pml_path, "w") as fh:
        fh.write(content)


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    df_all = pd.read_csv(FILTERED / "filtered_all_data.csv")
    df_patches = pd.read_csv(FILTERED / "patches_infos_cluster_data_70.csv")

    df_homo = df_all[df_all["homo"] == True].copy()
    df_hetero = df_all[df_all["homo"] == False].copy()

    # ── filament global base (même pour tous les ABPs) ──────────────────────
    global_sites = get_s1_sites(GLOBAL_PATCHES, df_patches)
    print(f"Sites globaux ({' + '.join(GLOBAL_PATCHES)}) : {global_sites}")

    if not GLOBAL_BASE_PDB.exists():
        print("Génération filament_global_base.pdb …")
        global_bf = read_bfactors_sum(global_sites)
        max_global = apply_bfactors(TEMPLATE_PDB, global_bf, GLOBAL_BASE_PDB)
        print(f"  → max B-factor global : {max_global:.2f}")
    else:
        # Lire le max depuis le fichier existant
        max_global = 0.0
        with open(GLOBAL_BASE_PDB) as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    b = float(line[60:66])
                    if b > max_global:
                        max_global = b
        print(f"filament_global_base.pdb existant, max B = {max_global:.2f}")

    # ── par ABP ─────────────────────────────────────────────────────────────
    abp_groups = df_hetero.groupby("subunit_2_title")
    n = len(abp_groups)
    print(f"\n{n} ABPs à traiter\n")

    for i, (abp_title, group) in enumerate(abp_groups, 1):
        sname = safe_name(abp_title)
        pml_path = OUT_DIR / f"{sname}.pml"
        abp_pdb_path = OUT_DIR / f"{sname}_abp.pdb"

        print(f"[{i:02d}/{n}] {abp_title[:50]}")

        abp_pdbs = set(group["pdb_id"])
        patches = top2_patches(abp_pdbs, df_homo)
        print(f"       top2 patches : {patches}")

        if not patches:
            print("       [skip] aucun patch homo trouvé")
            continue

        abp_sites = get_s1_sites(patches, df_patches)
        print(f"       sites S1 : {abp_sites}")

        # PDB ABP-spécifique
        if not abp_pdb_path.exists():
            abp_bf = read_bfactors_sum(abp_sites)
            max_abp = apply_bfactors(TEMPLATE_PDB, abp_bf, abp_pdb_path)
        else:
            max_abp = 0.0
            with open(abp_pdb_path) as fh:
                for line in fh:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        b = float(line[60:66])
                        if b > max_abp:
                            max_abp = b

        print(f"       max B ABP : {max_abp:.2f}")

        write_pml(pml_path, GLOBAL_BASE_PDB, abp_pdb_path,
                  max_global, max_abp, abp_title,
                  GLOBAL_PATCHES, patches)

    print("\nTerminé.")


if __name__ == "__main__":
    main()
