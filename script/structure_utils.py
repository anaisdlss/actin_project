"""Reconstruction d'un PDB « pairwise » (2 chaînes) à partir de l'assembly.

Les fichiers pairwise téléchargés depuis PPI3D ne sont qu'un extrait de
l'assembly biologique : 2 chaînes, renommées A/B, la chaîne d'origine étant
recopiée à l'index 71 (colonne 72). Coordonnées, numéros de série, occupancy et
B-factor sont IDENTIQUES à l'assembly (vérifié atome par atome).

On ne télécharge donc plus les pairwise : on les reconstruit à la volée depuis
`biological_assembly_pdb_file` quand un consommateur en a besoin.
"""
from pathlib import Path


def pairwise_text_from_assembly(assembly_path, chain_a, chain_b):
    """Reconstruit le texte PDB pairwise (chaîne A puis B) depuis l'assembly.

    chain_a / chain_b = lettres de chaîne PHYSIQUES dans l'assembly (ex. 'G',
    'H'). Renvoie une chaîne PDB (ou None si l'assembly est introuvable).
    Format identique à PPI3D : A/B en col 22, chaîne d'origine en col 72.
    """
    assembly_path = Path(assembly_path)
    if not assembly_path.exists():
        return None
    ren = {chain_a: "A", chain_b: "B"}
    # On regroupe par nouvelle chaîne (A d'abord, puis B) comme le format PPI3D.
    buckets = {"A": [], "B": []}
    for line in assembly_path.open():
        if line.startswith(("ATOM", "HETATM")) and len(line) > 21 and line[21] in ren:
            orig = line[21]
            new = ren[orig]
            base = line.rstrip("\n").ljust(72)
            # index 21 = nouvelle chaîne ; index 71 (col 72) = chaîne d'origine
            newl = base[:21] + new + base[22:71] + orig + base[72:]
            buckets[new].append(newl)
    if not buckets["A"] and not buckets["B"]:
        return None
    out = buckets["A"] + ["TER"] + buckets["B"] + ["TER", "END"]
    return "\n".join(out) + "\n"


def pairwise_filename_parts(name):
    """Décompose un nom de fichier pairwise `<iid>_<pdb>_<chA>_<chB>.pdb`.
    Renvoie (pdb_id, chain_a, chain_b) ou None."""
    stem = Path(name).stem
    parts = stem.split("_")
    if len(parts) >= 4:
        return parts[1], parts[-2], parts[-1]
    return None


PAIRWISE_DIR = Path("data/filtered/details/structures_files/pairwise")


def expected_pairwise_path(iid, pdb_id, chain_a, chain_b, pairwise_dir=None):
    """Chemin attendu du fichier pairwise `<iid>_<pdb>_<chA>_<chB>.pdb`."""
    d = Path(pairwise_dir) if pairwise_dir else PAIRWISE_DIR
    return d / f"{iid}_{str(pdb_id).lower()}_{chain_a}_{chain_b}.pdb"


def pairwise_text(iid, pdb_id, chain_a, chain_b, assembly_pdb=None,
                  pairwise_dir=None):
    """Texte PDB pairwise pour une interaction.

    1) si un assembly PDB existe → on le reconstruit (chaînes A/B) ;
    2) sinon (assembly CIF-only) → on lit le fichier pairwise conservé.
    Renvoie None si rien de disponible.
    """
    if assembly_pdb and Path(assembly_pdb).exists():
        t = pairwise_text_from_assembly(assembly_pdb, chain_a, chain_b)
        if t:
            return t
    fb = expected_pairwise_path(iid, pdb_id, chain_a, chain_b, pairwise_dir)
    if fb.exists():
        return fb.read_text()
    return None


def prune_pairwise(df8, pairwise_dir=None, used_iids=None):
    """Purge les fichiers pairwise inutiles. On SUPPRIME :
      - ceux dérivables d'un assembly .pdb (reconstructibles à la volée) ;
      - ceux d'interactions non filtrées (si `used_iids` est fourni).
    On ne GARDE que les pairwise d'interactions filtrées dont l'assembly n'est
    qu'en CIF (non reconstructible). Renvoie (supprimés, gardés).
    `df8` = 8.structures (pdb_id, biological_assembly_pdb_file)."""
    import os
    d = Path(pairwise_dir) if pairwise_dir else PAIRWISE_DIR
    if not d.exists():
        return 0, 0
    has_pdb_asm = set()
    for _, r in df8.iterrows():
        p = str(r.get("biological_assembly_pdb_file", ""))
        if p and p.lower() != "nan" and os.path.exists(p):
            has_pdb_asm.add(str(r.get("pdb_id", "")).lower())
    used = None if used_iids is None else {str(i) for i in used_iids}
    removed = kept = 0
    for f in d.glob("*.pdb"):
        iid = f.name.split("_")[0]
        parts = pairwise_filename_parts(f.name)
        pdb = parts[0].lower() if parts else ""
        drop = (pdb in has_pdb_asm) or (used is not None and iid not in used)
        if drop:
            f.unlink()
            removed += 1
        else:
            kept += 1
    return removed, kept
