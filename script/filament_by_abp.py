# filament_by_abp.py
# Pour chaque ABP, génère une session PyMOL avec 2 filaments actine 10 sous-unités :
#   filament_base : top 2 clusters C70 homo GLOBAUX (0_7797_0 + 0_7797_1),
#                   construit depuis le template 8iah, white→red
#   filament_abp  : top 2 clusters C70 homo SPÉCIFIQUES à l'ABP,
#                   construit depuis une structure pairwise du cluster dominant, white→orange
# B-factors = somme des % ASA buried des sites de liaison S1 des clusters C70
#
# Exécution : python -m script.filament_by_abp

import ast
import re
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
FILTERED     = PROJECT_ROOT / "data" / "filtered"
BFACTOR_DIR  = FILTERED / "details" / "structures_files" / "bfactor_cluster"
PAIRWISE_DIR = FILTERED / "details" / "structures_files" / "pairwise"
TEMPLATE_PDB = FILTERED / "details" / "structures_files" / "filament" / "actin_filament_10.pdb"
OUT_DIR      = FILTERED / "details" / "structures_files" / "filament" / "by_abp"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GLOBAL_BASE_PDB = FILTERED / "details" / "structures_files" / "filament" / "filament_global_base.pdb"
GLOBAL_PATCHES  = ["0_7797_0", "0_7797_1"]
CHAIN_LABELS    = "ABCDEFGHIJ"
N_SUBUNITS      = 10


# ── utilitaires ──────────────────────────────────────────────────────────────

def safe_name(title: str) -> str:
    return re.sub(r"[^a-zA-Z0-9]+", "_", title).strip("_")[:60]


def parse_ids(s: str) -> list[int]:
    s = re.sub(r"np\.int64\((\d+)\)", r"\1", str(s))
    return [int(x) for x in re.findall(r"\d+", s)]


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


def apply_bfactors_to_template(bfactors: dict[int, float], out_path: Path) -> float:
    """Applique les B-factors au template 8iah ; retourne le max."""
    out_lines = []
    max_b = 0.0
    with open(TEMPLATE_PDB) as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resi = int(line[22:26])
                b = bfactors.get(resi, 0.0)
                max_b = max(max_b, b)
                line = line[:60] + f"{b:6.2f}" + line[66:]
            out_lines.append(line)
    with open(out_path, "w") as fh:
        fh.writelines(out_lines)
    return max_b


def get_s1_sites(patches: list[str], df_patches: pd.DataFrame) -> list[str]:
    sites, seen = [], set()
    for patch in patches:
        row = df_patches[df_patches["patch"] == patch]
        if row.empty:
            continue
        for s in ast.literal_eval(row.iloc[0]["binding_site_clusters"]):
            if s not in seen:
                sites.append(s)
                seen.add(s)
    return sites


def top2_patches(abp_pdbs: set, df_homo: pd.DataFrame) -> list[str]:
    sub = df_homo[df_homo["pdb_id"].isin(abp_pdbs)]
    if sub.empty:
        return []
    cov = sub.groupby("cluster_data_70")["pdb_id"].nunique().sort_values(ascending=False)
    return list(cov.index[:2])


# ── construction de filament depuis pairwise ─────────────────────────────────

def _load_pairwise_atoms(path: Path) -> tuple[list[str], list[str]]:
    """Retourne (atoms_chain_A, atoms_chain_B) — toutes les lignes ATOM/HETATM."""
    a_lines, b_lines = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[21] == "A":
                    a_lines.append(line)
                elif line[21] == "B":
                    b_lines.append(line)
    return a_lines, b_lines


def _get_ca(atom_lines: list[str]) -> tuple[np.ndarray, list[int]]:
    coords, resis = [], []
    for line in atom_lines:
        if line[12:16].strip() == "CA":
            resi = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            coords.append([x, y, z])
            resis.append(resi)
    return np.array(coords), resis


def _kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = cQ - R @ cP
    return R, t


def _apply_helical(coords: np.ndarray, R: np.ndarray, t: np.ndarray,
                   k: int) -> np.ndarray:
    """Applique le transform hélicoïdal k fois à coords."""
    result = coords.copy()
    for _ in range(k):
        result = (R @ result.T).T + t
    return result


def _atom_xyz(line: str) -> np.ndarray:
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def _set_atom_xyz(line: str, xyz: np.ndarray) -> str:
    return (line[:30]
            + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            + line[54:])


def build_filament_from_pairwise(
    pairwise_path: Path,
    bfactors: dict[int, float],
    out_path: Path,
    n_subunits: int = N_SUBUNITS,
) -> float:
    """
    Construit un filament à n_subunits sous-unités en propageant le transform
    hélicoïdal A→B de la structure pairwise.
    Applique bfactors (par résidu) et écrit un PDB multi-chaîne.
    Retourne le max B-factor observé.
    """
    a_lines, b_lines = _load_pairwise_atoms(pairwise_path)
    if not a_lines or not b_lines:
        raise ValueError(f"Structure pairwise invalide : {pairwise_path.name}")

    ca_A, rA = _get_ca(a_lines)
    ca_B, rB = _get_ca(b_lines)
    common = sorted(set(rA) & set(rB))
    if len(common) < 50:
        raise ValueError(f"Trop peu de Cα communs ({len(common)}) dans {pairwise_path.name}")

    rA_idx = {r: i for i, r in enumerate(rA)}
    rB_idx = {r: i for i, r in enumerate(rB)}
    P = ca_A[[rA_idx[r] for r in common]]
    Q = ca_B[[rB_idx[r] for r in common]]
    R, t = _kabsch(P, Q)

    # Extraire les coordonnées XYZ de toutes les lignes de la chaîne A
    all_xyz_A = np.array([_atom_xyz(l) for l in a_lines])

    max_b = 0.0
    out_lines: list[str] = []
    serial = 1

    for k in range(n_subunits):
        chain_label = CHAIN_LABELS[k]
        # Coordonnées de ce subunit = transform hélicoïdal k fois sur chain A
        xyz_k = _apply_helical(all_xyz_A, R, t, k)

        for i, line in enumerate(a_lines):
            resi = int(line[22:26])
            b = bfactors.get(resi, 0.0)
            max_b = max(max_b, b)
            new_line = (
                line[:6]
                + f"{serial:5d}"
                + line[11:21]
                + chain_label
                + line[22:60]
                + f"{b:6.2f}"
                + line[66:]
            )
            new_line = _set_atom_xyz(new_line, xyz_k[i])
            out_lines.append(new_line)
            serial += 1

        out_lines.append(f"TER   {serial:5d}      {chain_label}\n")
        serial += 1

    out_lines.append("END\n")
    with open(out_path, "w") as fh:
        fh.writelines(out_lines)
    return max_b


C70_INTERFACE_DIR = FILTERED / "details" / "structures_files" / "bfactor_c70_interface"


def _helix_angle_and_longitudinal(path: Path) -> float | None:
    """Retourne l'angle hélicoïdal A→B si c'est un contact longitudinal, sinon None."""
    ca_A, ca_B = {}, {}
    try:
        with open(path) as fh:
            for line in fh:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                ch = line[21]
                if line[12:16].strip() != "CA":
                    continue
                resi = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                if ch == "A":
                    ca_A[resi] = xyz
                elif ch == "B":
                    ca_B[resi] = xyz
    except Exception:
        return None
    common = sorted(set(ca_A) & set(ca_B))
    if len(common) < 50:
        return None
    cA = np.array(list(ca_A.values())).mean(0)
    cB = np.array(list(ca_B.values())).mean(0)
    if not (25 < np.linalg.norm(cB - cA) < 55):
        return None  # contact non-longitudinal
    P = np.array([ca_A[r] for r in common])
    Q = np.array([ca_B[r] for r in common])
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))


def find_representative_pdb(patches: list[str]) -> Path | None:
    """
    Cherche le représentant C70 longitudinal (bfactor_c70_interface/) pour
    les patches donnés, en priorité le premier patch ayant un contact longitudinal.
    Fallback sur le premier fichier disponible.
    """
    fallback: Path | None = None
    for patch in patches:
        path = C70_INTERFACE_DIR / f"{patch}.pdb"
        if not path.exists():
            continue
        angle = _helix_angle_and_longitudinal(path)
        if angle is not None:
            return path  # contact longitudinal trouvé
        if fallback is None:
            fallback = path
    return fallback


# ── PML ──────────────────────────────────────────────────────────────────────

def write_pml(pml_path: Path, base_pdb: Path, abp_pdb: Path,
              max_base: float, max_abp: float, abp_title: str,
              base_patches: list[str], abp_patches: list[str]) -> None:
    content = f"""\
# PyMOL — Filament actine 10 s-u — {abp_title}
# filament_base : clusters globaux ({' + '.join(base_patches)}), white_red, max={max_base:.2f}
# filament_abp  : clusters spécifiques ({' + '.join(abp_patches) if abp_patches else 'aucun'}), white_orange, max={max_abp:.2f}

# ── filament_base (structure 8iah — B-factors clusters globaux) ──
load {base_pdb.as_posix()}, filament_base
hide everything, filament_base
show surface, filament_base
color white, filament_base
spectrum b, white_red, filament_base, minimum=0, maximum={max_base:.2f}

# ── filament_abp (structure pairwise ABP-spécifique — B-factors clusters ABP) ──
load {abp_pdb.as_posix()}, filament_abp
hide everything, filament_abp
show surface, filament_abp
color white, filament_abp
spectrum b, white_orange, filament_abp, minimum=0, maximum={max_abp:.2f}

set surface_quality, 1
bg_color white
zoom filament_base
"""
    with open(pml_path, "w") as fh:
        fh.write(content)


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    df_all     = pd.read_csv(FILTERED / "filtered_all_data.csv")
    df_patches = pd.read_csv(FILTERED / "patches_infos_cluster_data_70.csv")

    df_homo   = df_all[df_all["homo"] == True].copy()
    df_hetero = df_all[df_all["homo"] == False].copy()

    # Mapping patch → IDs d'interactions
    patch_to_ids: dict[str, list[int]] = {}
    for _, row in df_patches.iterrows():
        patch_to_ids[row["patch"]] = parse_ids(row["ids_interactions"])

    # ── filament global base (même pour tous les ABPs) ──────────────────────
    global_sites = get_s1_sites(GLOBAL_PATCHES, df_patches)
    print(f"Sites globaux ({' + '.join(GLOBAL_PATCHES)}) : {global_sites}")

    if not GLOBAL_BASE_PDB.exists():
        print("Génération filament_global_base.pdb …")
        global_rep = find_representative_pdb(GLOBAL_PATCHES)
        global_bf = read_bfactors_sum(global_sites)
        if global_rep is not None:
            print(f"  représentant C70 global : {global_rep.name}")
            max_global = build_filament_from_pairwise(global_rep, global_bf, GLOBAL_BASE_PDB)
        else:
            max_global = apply_bfactors_to_template(global_bf, GLOBAL_BASE_PDB)
        print(f"  → max B-factor global : {max_global:.2f}")
    else:
        max_global = 0.0
        with open(GLOBAL_BASE_PDB) as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        max_global = max(max_global, float(line[60:66]))
                    except ValueError:
                        pass
        print(f"filament_global_base.pdb existant, max B = {max_global:.2f}")

    # ── par ABP ─────────────────────────────────────────────────────────────
    abp_groups = df_hetero.groupby("subunit_2_title")
    n = len(abp_groups)
    print(f"\n{n} ABPs à traiter\n")

    for i, (abp_title, group) in enumerate(abp_groups, 1):
        sname    = safe_name(abp_title)
        pml_path = OUT_DIR / f"{sname}.pml"
        abp_pdb  = OUT_DIR / f"{sname}_abp.pdb"

        print(f"[{i:02d}/{n}] {abp_title[:55]}")

        abp_pdbs = set(group["pdb_id"])
        patches  = top2_patches(abp_pdbs, df_homo)
        print(f"       top2 patches : {patches}")

        if not patches:
            print("       [skip] aucun patch homo trouvé")
            continue

        abp_sites = get_s1_sites(patches, df_patches)
        print(f"       sites S1 : {abp_sites}")

        if not abp_pdb.exists():
            # Utiliser le représentant C70 officiel (bfactor_c70_interface/)
            rep_pdb = find_representative_pdb(patches)

            abp_bf = read_bfactors_sum(abp_sites)
            if rep_pdb is None:
                print(f"       [warn] aucun représentant C70 trouvé, utilisation du template 8iah")
                max_abp = apply_bfactors_to_template(abp_bf, abp_pdb)
            else:
                print(f"       représentant C70 : {rep_pdb.name}")
                try:
                    max_abp = build_filament_from_pairwise(rep_pdb, abp_bf, abp_pdb)
                except Exception as e:
                    print(f"       [warn] construction échouée ({e}), fallback template")
                    max_abp = apply_bfactors_to_template(abp_bf, abp_pdb)
        else:
            max_abp = 0.0
            with open(abp_pdb) as fh:
                for line in fh:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            max_abp = max(max_abp, float(line[60:66]))
                        except ValueError:
                            pass

        print(f"       max B ABP : {max_abp:.2f}")
        write_pml(pml_path, GLOBAL_BASE_PDB, abp_pdb,
                  max_global, max_abp, abp_title,
                  GLOBAL_PATCHES, patches)

    print("\nTerminé.")


if __name__ == "__main__":
    main()
