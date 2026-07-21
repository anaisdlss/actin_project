#!/usr/bin/env python
"""
Balayage TM-score de la ZONE qui touche l'actine, à contexte croissant.

Méthode : TM-align pairwise DIRECT (tmtools) — pas de préfiltre, fonctionne sur
les fragments d'interface. Pour PAD = 0, 6, 12, 20, 35 résidus de contexte autour
de chaque résidu de contact, + la protéine ENTIÈRE.

On compare, pour les paires d'ABP du MÊME cluster de site :
  - "même famille"        (≥1 domaine Pfam commun ou même UniProt)
  - "familles différentes"

Sorties :
  data/exports/abp_site_domain/interface_tm_sweep.csv
  data/exports/abp_site_domain/figure_interface_tm_sweep.png
"""
import re
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.SeqUtils import seq1
from tmtools import tm_align

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
ASM = ROOT / "data/filtered/details/structures_files/assembly"
CHAINS = OUT / "abp_chains"
TAG_DROP = {"PF01547", "PF01035", "PF02870"}
PADS = [0, 6, 12, 20, 35]

rep = pd.read_csv(OUT / "abp_representatives.csv")
res = pd.read_csv(ROOT / "data/filtered/details/3.interface_residues.csv")
ip = pd.read_csv(OUT / "abp_interpro.csv")
t = pd.read_csv(OUT / "abp_site_table.csv")


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")[:50]


pfam = {r.abp_title: ((set() if pd.isna(r.pfam_acc) else set(str(r.pfam_acc).split(";"))) - TAG_DROP) for _, r in ip.iterrows()}
uni = {r.abp_title: (set() if pd.isna(r.uniprot) else set(str(r.uniprot).split(";"))) for _, r in ip.iterrows()}
stem = {r.abp_title: f"{slug(r.abp_title)}__{r.pdb}_{r.abp_chain}" for _, r in rep.iterrows()}

# paires intra-cluster
pt = t.dropna(subset=["actin_site_cluster"]).drop_duplicates(["actin_site_cluster", "abp_title"])
multi = pt.groupby("actin_site_cluster").filter(lambda g: g.abp_title.nunique() >= 2)
intra = {frozenset((a, b)) for site, g in multi.groupby("actin_site_cluster")
         for a, b in combinations(sorted(g.abp_title.unique()), 2)}


def same_family(a, b):
    if uni.get(a) and uni.get(b) and (uni[a] & uni[b]):
        return True
    return bool(pfam.get(a) and pfam.get(b) and (pfam[a] & pfam[b]))


class ResSel(Select):
    def __init__(self, ch, ids):
        self.ch, self.ids = ch, ids
    def accept_chain(self, c):
        return c.id == self.ch
    def accept_residue(self, r):
        return r.id[0] == " " and r.id[1] in self.ids


io = PDBIO()


def extract(pad):
    dst = OUT / f"_sweep_pad{pad}"
    dst.mkdir(exist_ok=True)
    for f in dst.glob("*.pdb"):
        f.unlink()
    for _, r in rep.iterrows():
        sub = res[(res.interaction_id == int(r.interaction_id)) & (res.chain == f"{r.pdb}_{r.abp_chain}")]
        contacts = sorted(pd.to_numeric(sub.residue_number_structure, errors="coerce").dropna().astype(int))
        if not contacts:
            continue
        keep = {c + d for c in contacts for d in range(-pad, pad + 1)}
        base = f"{r.pdb}"
        f = ASM / (base + ".pdb"); parser = PDBParser(QUIET=True)
        if not f.exists():
            f = ASM / (base + ".cif"); parser = MMCIFParser(QUIET=True)
        st = parser.get_structure(base, str(f))
        io.set_structure(st)
        io.save(str(dst / (stem[r.abp_title] + ".pdb")), ResSel(r.abp_chain, keep))
    return dst


def load_ca(pdb_path):
    """Renvoie (coords Nx3, seq) des CA d'un fichier PDB mono-chaîne."""
    st = PDBParser(QUIET=True).get_structure("x", str(pdb_path))
    coords, seq = [], []
    for res_ in st.get_residues():
        if "CA" in res_:
            coords.append(res_["CA"].coord)
            seq.append(seq1(res_.resname, custom_map={}) or "X")
    return np.array(coords, dtype=float), "".join(seq)


def tm_pair(p1, p2):
    if not p1.exists() or not p2.exists():
        return np.nan
    c1, s1 = load_ca(p1)
    c2, s2 = load_ca(p2)
    if len(c1) < 3 or len(c2) < 3:
        return np.nan
    r = tm_align(c1, c2, s1, s2)
    return max(r.tm_norm_chain1, r.tm_norm_chain2)


rows = []
levels = []
for pad in PADS:
    d = extract(pad)
    lab = f"±{pad}"
    levels.append(lab)
    n_ok = 0
    for pr in intra:
        a, b = tuple(pr)
        tm = tm_pair(d / (stem[a] + ".pdb"), d / (stem[b] + ".pdb"))
        if not np.isnan(tm):
            n_ok += 1
        rows.append(dict(level=lab, a=a, b=b, tm=tm, fam=same_family(a, b)))
    print(f"PAD={pad}: {n_ok}/{len(intra)} paires alignées")

# protéine entière (même méthode tmtools)
levels.append("entière")
for pr in intra:
    a, b = tuple(pr)
    tm = tm_pair(CHAINS / (stem[a] + ".pdb"), CHAINS / (stem[b] + ".pdb"))
    rows.append(dict(level="entière", a=a, b=b, tm=tm, fam=same_family(a, b)))
print(f"protéine entière : calculée")

df = pd.DataFrame(rows)
df.to_csv(OUT / "interface_tm_sweep.csv", index=False)

print("\n=== TM médian des paires intra-cluster ===")
print(f"{'niveau':>10} | {'même fam. (méd)':>16} | {'fam. diff. (méd)':>16} | {'fam.diff TM≥0.5':>15}")
for lab in levels:
    sub = df[df.level == lab]
    mf = sub[sub.fam].tm.median()
    fd = sub[~sub.fam].tm
    print(f"{lab:>10} | {mf:>16.2f} | {fd.median():>16.2f} | {(fd>=0.5).sum():>3}/{fd.notna().sum()}")

# figure
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(11, 6))
xs = range(len(levels))
for fam, color, name in [(True, "#2a9d8f", "même famille"), (False, "#e76f51", "familles différentes")]:
    med = [df[(df.level == l) & (df.fam == fam)].tm.median() for l in levels]
    ax.plot(list(xs), med, "o-", color=color, lw=2, label=name)
    for xi, l in zip(xs, levels):
        v = df[(df.level == l) & (df.fam == fam)].tm.dropna()
        ax.scatter([xi + (0.08 if fam else -0.08)] * len(v), v, color=color, s=16, alpha=0.5)
ax.axhline(0.5, ls="--", color="grey")
ax.text(0, 0.52, "seuil 'même fold' = 0.5", color="grey", fontsize=9)
ax.set_xticks(list(xs))
ax.set_xticklabels([l if l != "entière" else "protéine\nentière" for l in levels])
ax.set_xlabel("Étendue comparée : ±N résidus de contexte autour des contacts → protéine entière")
ax.set_ylabel("TM-score (paires d'ABP du même site)")
ax.set_title("La zone qui touche l'actine est-elle structurellement comparable ? (TM-align direct)")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "figure_interface_tm_sweep.png", dpi=150, bbox_inches="tight")
print(f"\nfigure : {OUT/'figure_interface_tm_sweep.png'}")
