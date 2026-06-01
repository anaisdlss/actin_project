"""
Génère des PDB et scripts PyMOL où le B-factor encode la classe physicochimique.

B-factor par résidu d'interface (b_original > 0 dans le PDB source) :
  1 = hydrophobe         (ALA GLY ILE LEU MET VAL)
  2 = aromatique         (PHE TRP TYR)
  3 = polaire non chargé (HIS ASN GLN SER THR)
  4 = positif            (LYS ARG)
  5 = négatif            (ASP GLU)
  6 = cystéine           (CYS)
  7 = proline            (PRO)
  0 = non-interface (ou acide aminé inconnu)

Usage :
  python -m script.bfactor_physicochemical_pdb

Sortie :
  data/filtered/details/structures_files/bfactor_physicochemical/
    {site}.pdb   ← B-factor = classe physicochimique
    {site}.pml   ← script PyMOL prêt à l'emploi
"""

import re
from pathlib import Path

SRC_DIR = Path("data/filtered/details/structures_files/bfactor_cluster")
OUT_DIR = Path("data/filtered/details/structures_files/bfactor_physicochemical")
OUT_DIR.mkdir(parents=True, exist_ok=True)

_PHYSICO: dict[str, int] = {}
for aa in ["ALA", "GLY", "ILE", "LEU", "MET", "VAL"]:
    _PHYSICO[aa] = 1
for aa in ["PHE", "TRP", "TYR"]:
    _PHYSICO[aa] = 2
for aa in ["HIS", "ASN", "GLN", "SER", "THR"]:
    _PHYSICO[aa] = 3
for aa in ["LYS", "ARG"]:
    _PHYSICO[aa] = 4
for aa in ["ASP", "GLU"]:
    _PHYSICO[aa] = 5
_PHYSICO["CYS"] = 6
_PHYSICO["PRO"] = 7

COLOR_MAP = {
    0: ("orange",    "non-interface"),
    1: ("grey60",    "hydrophobe (Ala Gly Ile Leu Met Val)"),
    2: ("hotpink",   "aromatique (Phe Trp Tyr)"),
    3: ("cyan",      "polaire (His Asn Gln Ser Thr)"),
    4: ("blue",      "positif (Lys Arg)"),
    5: ("red",       "négatif (Asp Glu)"),
    6: ("yellow",    "cystéine (Cys)"),
    7: ("limegreen", "proline (Pro)"),
}


def _recode_pdb(src: Path, dst: Path) -> int:
    """Réécrit le PDB en remplaçant le B-factor par la classe physicochimique.
    Retourne le nombre de résidus d'interface recodés."""
    n_interface = 0
    lines_out = []
    with open(src) as f:
        for line in f:
            if line.startswith("ATOM") and len(line) > 60:
                resn = line[17:20].strip()
                try:
                    b_orig = float(line[60:66])
                except ValueError:
                    lines_out.append(line)
                    continue
                if b_orig > 0:
                    code = _PHYSICO.get(resn, 0)
                    n_interface += 1
                else:
                    code = 0
                new_b = f"{code:6.2f}"
                line = line[:60] + new_b + line[66:]
                lines_out.append(line)
            else:
                lines_out.append(line)
    dst.write_text("".join(lines_out))
    return n_interface


def _write_pml(pdb_abs: Path, pml_path: Path, site: str) -> None:
    obj = re.sub(r"[^A-Za-z0-9_]", "_", site)
    cmds = [
        f"# PyMOL — classes physicochimiques — site actine : {site}",
        f"# B-factor : 0=non-interface · 1=hydrophobe · 2=aromatique · 3=polaire",
        f"#            4=positif · 5=négatif · 6=Cys · 7=Pro",
        "",
        f"load {pdb_abs}, {obj}",
        f"hide everything, {obj}",
        f"show surface, {obj}",
        "",
        "# ── Couleurs physicochimiques ───────────────────────────────────────",
        f"color orange,    {obj} and b = 0   # non-interface",
        f"color grey60,    {obj} and b = 1   # hydrophobe",
        f"color hotpink,   {obj} and b = 2   # aromatique",
        f"color cyan,      {obj} and b = 3   # polaire",
        f"color blue,      {obj} and b = 4   # positif",
        f"color red,       {obj} and b = 5   # négatif",
        f"color yellow,    {obj} and b = 6   # cystéine",
        f"color limegreen, {obj} and b = 7   # proline",
        "",
        "set surface_quality, 1",
        "bg_color white",
        f"zoom {obj}",
    ]
    pml_path.write_text("\n".join(cmds))


n_ok = 0
for src in sorted(SRC_DIR.glob("*.pdb")):
    site = src.stem
    dst_pdb = OUT_DIR / src.name
    dst_pml = OUT_DIR / f"{site}.pml"
    n = _recode_pdb(src, dst_pdb)
    _write_pml(dst_pdb.resolve(), dst_pml, site)
    n_ok += 1

print(f"{n_ok} PDB + PML générés dans {OUT_DIR}/")
print()
print("Usage dans PyMOL :")
print("  1. Ouvrir le .pml avec 'File > Run Script...'  (ou '@chemin/site.pml' dans la console)")
print("  2. Les résidus d'interface sont colorés par classe physicochimique")
print("  3. La surface orange = résidus non-interface de l'actine")
print()
print("Légende B-factor :")
for code, (color, label) in COLOR_MAP.items():
    print(f"  b = {code} → {color:12s} ({label})")
