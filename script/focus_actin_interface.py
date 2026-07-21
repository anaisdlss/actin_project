#!/usr/bin/env python3
"""Focus « interface actine » — figures pour le rapport.

Produit 4 livrables autour de la vue fusionnée (résidus ABP projetés sur l'actine) :
  1. PyMOL : surface actine colorée par nb d'ABP contactant (gradient hotspots)  -> .pml
  2. Conservation physico-chimique : identité vs propriété du résidu ABP par position
  3. Spécificité / entropie du résidu ABP par position actine
  4. Hubs vs conservation actine (MSA) + enfouissement (ASA)

Sorties : data/visualisations/focus/*.png  et  *_hotspots.pml

Usage : python -m script.focus_actin_interface
"""
import collections
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner

CONTACTS = Path("data/filtered/details/4.inter-residue_contacts.csv")
FILTERED = Path("data/filtered/filtered_all_data.csv")
PDB      = Path("actin_human_pdb.pdb")                       # P60709, repère pour GEMME/MSA
REF_STRUCT = Path("actin_ref_alpha_cardiac_8zbn.pdb")        # α-cardiaque 8zbn (chaîne B) — repère 3D
REF_CHAIN  = "B"
A3M      = Path("data/AF-P60709-F1-msa_v6.a3m")
GEMME    = Path("proteocast_evalue0_actine/14.query_GEMME_pLDDT.csv")
OUTDIR   = Path("data/visualisations/focus")

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "HIC": "H"}

# classes physico-chimiques
CLASS = {}
for a in "KRH": CLASS[a] = "Basique (+)"
for a in "DE":  CLASS[a] = "Acide (−)"
for a in "STNQCY": CLASS[a] = "Polaire"
for a in "AVLIMFW": CLASS[a] = "Hydrophobe"
for a in "GP":  CLASS[a] = "Gly/Pro"
CLASS_COL = {"Basique (+)": "#2166ac", "Acide (−)": "#b2182b",
             "Polaire": "#1b9e77", "Hydrophobe": "#e6ab02", "Gly/Pro": "#7570b3"}


def to1(name):
    name = str(name).strip().upper()
    return name if len(name) == 1 else AA3.get(name, "X")


def is_actin(t):
    t = str(t).lower()
    return t.startswith("actin") or t in ("actin", "actin-1")


# ── 1. chargement + assemblage ────────────────────────────────────────────────
def load_contacts():
    ct = pd.read_csv(CONTACTS)
    fa = pd.read_csv(FILTERED, low_memory=False)
    s2 = fa[["subunit_2", "subunit_2_title", "s2_sequence_cluster_70", "pdb_id"]].drop_duplicates("subunit_2")
    title = dict(zip(s2["subunit_2"], s2["subunit_2_title"]))
    pdbof = dict(zip(s2["subunit_2"], s2["pdb_id"]))
    clof  = dict(zip(s2["subunit_2"], s2["s2_sequence_cluster_70"]))
    ct["abp"]  = ct["chain_B_id"].map(title)
    ct["pdb"]  = ct["chain_B_id"].map(pdbof)
    ct["clu"]  = ct["chain_B_id"].map(clof)
    ct["ca"]   = pd.to_numeric(ct["residue_A_canon_mafft"], errors="coerce")
    ct["area"] = pd.to_numeric(ct["contact_area"].astype(str).str.replace("<", "", regex=False),
                               errors="coerce").fillna(0.05)
    ct["asa"]  = pd.to_numeric(ct["asa_pct_A"], errors="coerce")
    ct["aa_b"] = ct["residue_B_name"].map(to1)
    ct["aa_a"] = ct["residue_A_name"].map(to1)
    ct = ct.dropna(subset=["ca", "abp"])
    ct = ct[~ct["abp"].map(is_actin)]
    ct = ct[ct["pdb"].astype(str).str.lower() != "4b1z"]      # exclusion 4b1z (PHACTR1)
    ct["ca"] = ct["ca"].astype(int)
    return ct


# ── 2. mapping position canonique -> resid PDB (= colonne a3m) ─────────────────
def canon_to_pdb(ct, pdb_path=PDB, chain="A"):
    canon_aa = {p: collections.Counter(g).most_common(1)[0][0]
                for p, g in ct.groupby("ca")["aa_a"]}
    canon_list = sorted(canon_aa)
    canon_seq  = "".join(canon_aa[c] for c in canon_list)
    res = {}
    for ln in pdb_path.read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA" and ln[21] == chain:
            res[int(ln[22:26])] = AA3.get(ln[17:20].strip(), "X")
    pdb_ids = sorted(res)
    pdb_seq = "".join(res[i] for i in pdb_ids)
    al = PairwiseAligner()
    al.mode = "global"; al.match_score = 2; al.mismatch_score = -1
    al.open_gap_score = -10; al.extend_gap_score = -0.5
    aln = al.align(canon_seq, pdb_seq)[0]
    canon2pdb = {}
    for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(a1 - a0):
            canon2pdb[canon_list[a0 + k]] = pdb_ids[b0 + k]
    return canon2pdb


# ── 3. conservation actine depuis la MSA a3m ──────────────────────────────────
def actin_conservation():
    """Retourne dict resid_pdb(1..375) -> conservation (0..1, fraction AA majoritaire)."""
    seqs = []
    cur = None
    for l in A3M.read_text().splitlines():
        if l.startswith(">"):
            if cur is not None:
                seqs.append(cur)
            cur = ""
        else:
            cur += l
    if cur:
        seqs.append(cur)
    # garder seulement colonnes majuscules/gap (positions de la query), retirer insertions
    keep = [c.isupper() or c == "-" for c in seqs[0].replace(".", "-")]
    cons = {}
    L = sum(keep)
    cols = [[] for _ in range(L)]
    for s in seqs:
        s = s.replace(".", "-")
        j = 0
        for c, k in zip(s, keep):
            if k:
                cols[j].append(c.upper())
                j += 1
    for i, col in enumerate(cols):
        vals = [c for c in col if c != "-"]
        if not vals:
            cons[i + 1] = np.nan
            continue
        cnt = collections.Counter(vals)
        cons[i + 1] = cnt.most_common(1)[0][1] / len(vals)
    return cons


# ── 4. métriques par position canonique ───────────────────────────────────────
def per_position(ct):
    # résidu ABP dominant (area max) par (abp, canon)
    dom = (ct.groupby(["abp", "ca", "aa_b"])["area"].sum().reset_index()
           .sort_values("area", ascending=False).drop_duplicates(["abp", "ca"]))
    rows = []
    for ca, g in dom.groupby("ca"):
        aas = list(g["aa_b"])
        n = len(aas)
        cnt = collections.Counter(aas)
        id_cons = cnt.most_common(1)[0][1] / n                       # conservation d'identité
        classes = [CLASS.get(a, "Autre") for a in aas]
        ccnt = collections.Counter(classes)
        prop_cons = ccnt.most_common(1)[0][1] / n                    # conservation de propriété
        p = np.array([v / n for v in cnt.values()])
        entropy = -(p * np.log2(p)).sum()
        rows.append({"ca": ca, "n_abp": n, "id_cons": id_cons,
                     "prop_cons": prop_cons, "entropy": entropy,
                     "dom_class": ccnt.most_common(1)[0][0]})
    df = pd.DataFrame(rows)
    asa = ct.groupby("ca")["asa"].median()
    aa_a = ct.groupby("ca")["aa_a"].agg(lambda x: collections.Counter(x).most_common(1)[0][0])
    df["asa"]  = df["ca"].map(asa)
    df["aa_a"] = df["ca"].map(aa_a)
    return df


# ── figures ───────────────────────────────────────────────────────────────────
def fig_entropy(df):
    d = df.sort_values("ca")
    fig, ax = plt.subplots(figsize=(13, 3.2))
    ax.bar(range(len(d)), d["entropy"], color="#34495e", width=0.9)
    ax.set_xticks(range(0, len(d), 12))
    ax.set_xticklabels(d["ca"].iloc[::12], fontsize=7, rotation=90)
    ax.set_xlabel("Position canonique actine"); ax.set_ylabel("Entropie de Shannon\n(résidus ABP)")
    ax.set_title("Diversité de la chimie de liaison ABP par position actine "
                 "(haute = chaque ABP utilise un résidu différent)", fontsize=10)
    ax.margins(x=0.005)
    fig.tight_layout(); fig.savefig(OUTDIR / "fig3_entropie_par_position.png", dpi=160)
    plt.close(fig)


def fig_physchem(df):
    fig, ax = plt.subplots(figsize=(6.2, 6))
    for cl, sub in df.groupby("dom_class"):
        ax.scatter(sub["id_cons"], sub["prop_cons"], s=18 + sub["n_abp"] * 1.2,
                   color=CLASS_COL.get(cl, "#888"), alpha=0.7, edgecolor="none", label=cl)
    ax.plot([0, 1], [0, 1], "--", color="#999", lw=1)
    ax.fill_between([0, 1], [0, 1], [1, 1], color="#2ecc71", alpha=0.06)
    ax.set_xlabel("Conservation d'identité du résidu ABP\n(fraction AA majoritaire)")
    ax.set_ylabel("Conservation de propriété physico-chimique\n(fraction classe majoritaire)")
    ax.set_title("Chimie conservée malgré identité variable\n"
                 "(points au-dessus de la diagonale = propriété conservée)", fontsize=10)
    ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.02)
    ax.legend(fontsize=7, loc="lower right", title="Classe dominante", title_fontsize=7)
    fig.tight_layout(); fig.savefig(OUTDIR / "fig2_conservation_physicochimique.png", dpi=160)
    plt.close(fig)


def fig_hubs(df):
    d = df.dropna(subset=["gemme"])
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    sc = ax.scatter(d["n_abp"], d["gemme"], c=d["id_cons"], cmap="plasma_r",
                    s=30, alpha=0.85, edgecolor="white", linewidth=0.3)
    cb = fig.colorbar(sc); cb.set_label("Conservation d'identité du résidu ABP")
    if len(d) > 3:
        z = np.polyfit(d["n_abp"], d["gemme"], 1)
        xs = np.linspace(d["n_abp"].min(), d["n_abp"].max(), 50)
        ax.plot(xs, np.polyval(z, xs), color="#2c3e50", lw=1.4)
        r = d["n_abp"].corr(d["gemme"], method="spearman")
        ax.text(0.04, 0.06, f"Spearman ρ = {r:.2f}", transform=ax.transAxes,
                fontsize=9, color="#2c3e50")
    ax.set_xlabel("Nb d'ABP contactant la position (hub)")
    ax.set_ylabel("Contrainte évolutive actine\n(GEMME, élevé = conservé)")
    ax.set_title("Les hubs d'interaction sont-ils plus conservés ?", fontsize=10)
    fig.tight_layout(); fig.savefig(OUTDIR / "fig4_hubs_vs_conservation.png", dpi=160)
    plt.close(fig)


def write_pml(df, canon2pdb, pdb_path, chain):
    # b-factor = n_abp, spectrum sur la surface
    lines = ["# Hotspots ABP sur l'actine (α-cardiaque 8zbn, représentant le plus fréquent myosine)",
             "# surface colorée par nb d'ABP contactant la position",
             f"load {pdb_path.resolve()}, actin", "hide everything, actin",
             f"show surface, actin and chain {chain}", "color grey90, actin", "alter actin, b=0"]
    for _, r in df.iterrows():
        rp = canon2pdb.get(int(r["ca"]))
        if rp is not None:
            lines.append(f"alter actin and chain {chain} and resi {rp}, b={int(r['n_abp'])}")
    vmax = int(df["n_abp"].max())
    lines += [f"spectrum b, white_purple, actin, minimum=0, maximum={vmax}",
              "bg_color white", "set ray_shadows, 0", "rebuild", "deselect",
              f"# échelle b-factor 0..{vmax} ABP"]
    (OUTDIR / "fig1_hotspots.pml").write_text("\n".join(lines) + "\n")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ct = load_contacts()
    canon2cons   = canon_to_pdb(ct, PDB, "A")                     # -> P60709 (GEMME / MSA)
    canon2struct = canon_to_pdb(ct, REF_STRUCT, REF_CHAIN)        # -> 8zbn (rendu 3D)
    cons = actin_conservation()
    gemme = pd.read_csv(GEMME).set_index("Residue")["GEMME_mean"].to_dict()
    df = per_position(ct)
    df["cons"]  = df["ca"].map(lambda c: cons.get(canon2cons.get(c), np.nan))
    df["gemme"] = df["ca"].map(lambda c: gemme.get(canon2cons.get(c), np.nan))
    df.to_csv(OUTDIR / "per_position_metrics.csv", index=False)

    fig_entropy(df); fig_physchem(df); fig_hubs(df)

    # empreinte actine 3D restreinte aux MYOSINES (lisible, vs 46 ABP « violet partout »)
    _t = ct["abp"].astype(str).str.lower()
    myo_ct = ct[_t.str.contains("myosin") & ~_t.str.contains("tropomyosin")]
    df_myo = myo_ct.groupby("ca")["abp"].nunique().reset_index(name="n_abp")
    write_pml(df_myo, canon2struct, REF_STRUCT, REF_CHAIN)
    print(f"empreinte actine (myosines) : {len(df_myo)} positions, "
          f"max {int(df_myo['n_abp'].max())} myosines/position")

    n_above = ((df["prop_cons"] - df["id_cons"]) > 0.15).sum()
    print(f"{len(df)} positions actine · {ct['abp'].nunique()} ABP distinctes")
    print(f"hotspot max : pos {int(df.loc[df['n_abp'].idxmax(),'ca'])} "
          f"({int(df['n_abp'].max())} ABP)")
    print(f"positions à chimie conservée mais identité variable (Δ>0.15) : {n_above}")
    rho = df.dropna(subset=['gemme'])[['n_abp', 'gemme']].corr(method='spearman').iloc[0, 1]
    print(f"corrélation hubs↔contrainte GEMME (Spearman) : {rho:.2f}")
    rho2 = df[['n_abp', 'id_cons']].corr(method='spearman').iloc[0, 1]
    print(f"corrélation hubs↔conservation identité ABP (Spearman) : {rho2:.2f}")
    print(f"figures + .pml dans {OUTDIR}/")


if __name__ == "__main__":
    main()
