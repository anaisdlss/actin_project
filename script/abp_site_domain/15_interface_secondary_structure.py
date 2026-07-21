"""
Structure secondaire de la ZONE qui touche l'actine, pour chaque ABP.

Pour chaque ABP représentatif : on assigne la structure secondaire (PyMOL dss)
sur sa chaîne, puis on regarde la SS des résidus qui contactent l'actine.
=> dit si l'ABP agrippe l'actine via une hélice (H), un brin β (S) ou une boucle (L).

Lancer : pixi run pymol -cq script/abp_site_domain/15_interface_secondary_structure.py
Sortie : data/exports/abp_site_domain/interface_secondary_structure.csv
"""
import csv
import re
from pathlib import Path
from collections import Counter
from pymol import cmd

ROOT = Path("/Users/user/Desktop/stage/actin_project")
OUT = ROOT / "data/exports/abp_site_domain"
CHAINS = OUT / "abp_chains"
REP = OUT / "abp_representatives.csv"
RES = ROOT / "data/filtered/details/3.interface_residues.csv"


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")[:50]


# représentants : abp_title -> (iid, pdb, abp_chain)
reps = {}
with open(REP) as f:
    for r in csv.DictReader(f):
        reps[r["abp_title"]] = (int(float(r["interaction_id"])), r["pdb"], r["abp_chain"])

# résidus d'interface côté ABP : (iid, "pdb_chain") -> set(residue_number_structure)
contacts = {}
with open(RES) as f:
    for r in csv.DictReader(f):
        try:
            iid = int(float(r["interaction_id"]))
            resi = int(float(r["residue_number_structure"]))
        except (ValueError, KeyError):
            continue
        contacts.setdefault((iid, r["chain"]), set()).add(resi)

rows = []
for title, (iid, pdb, ch) in reps.items():
    pdbfile = CHAINS / f"{slug(title)}__{pdb}_{ch}.pdb"
    cset = contacts.get((iid, f"{pdb}_{ch}"), set())
    if not pdbfile.exists() or not cset:
        continue
    cmd.delete("all")
    cmd.load(str(pdbfile), "m")
    cmd.dss("m")                         # assignation structure secondaire
    ss_by_resi = {}
    cmd.iterate("m and name CA",
                "ss_by_resi[int(resi)] = (ss or 'L')",
                space={"ss_by_resi": ss_by_resi})
    # SS des résidus de contact
    cc = Counter(ss_by_resi.get(r, "L") for r in cset if r in ss_by_resi)
    n = sum(cc.values())
    if n == 0:
        continue
    h, s, l = cc.get("H", 0), cc.get("S", 0), cc.get("L", 0)
    dom = max([("hélice", h), ("brin β", s), ("boucle", l)], key=lambda x: x[1])[0]
    rows.append({
        "abp_title": title, "pdb": pdb, "chain": ch,
        "n_contacts_SS": n,
        "pct_helice": round(100 * h / n),
        "pct_brin": round(100 * s / n),
        "pct_boucle": round(100 * l / n),
        "SS_dominante": dom,
    })

rows.sort(key=lambda r: r["abp_title"])
with open(OUT / "interface_secondary_structure.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print(f"{len(rows)} ABP annotés -> interface_secondary_structure.csv\n")
print(f"{'ABP':40s} {'n':>3} {'%H':>4} {'%S':>4} {'%L':>4}  dominante")
for r in rows:
    print(f"{r['abp_title'][:40]:40s} {r['n_contacts_SS']:>3} "
          f"{r['pct_helice']:>4} {r['pct_brin']:>4} {r['pct_boucle']:>4}  {r['SS_dominante']}")
