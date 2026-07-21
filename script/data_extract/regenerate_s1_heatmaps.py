#!/usr/bin/env python
"""
Régénère les heatmaps S1 PAR CLUSTER (profil + décomposition C70) — RÉGLAGES D'ORIGINE
(toutes les positions, polices d'origine). Sert à restaurer l'état initial.

Sorties (écrase) :
  data/visualisations/actin_s1_clusters/{patch}.png
  data/visualisations/actin_s1_clusters_by_c70/{patch}.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
VIS = ROOT / "data/visualisations"
OUT1 = VIS / "actin_s1_clusters"
OUT2 = VIS / "actin_s1_clusters_by_c70"
OUT1.mkdir(parents=True, exist_ok=True)
OUT2.mkdir(parents=True, exist_ok=True)

area = pd.read_csv(ROOT / "data/filtered/actin_s1_canon_area_by_cluster.csv", index_col="patch")
positions = [int(p) for p in area.columns]
vmax_global = max(float(area.values.max()), 1.0)


def _xticks(ax, pos):
    step = max(1, len(pos) // 50)
    idx = list(range(0, len(pos), step))
    ax.set_xticks(idx)
    ax.set_xticklabels([pos[i] for i in idx], fontsize=6, rotation=90)


# ── 1) profil mono-ligne par patch (toutes positions) ─────────────────────────
for patch, row in area.iterrows():
    vals = row.values.astype(float).reshape(1, -1)
    fig, ax = plt.subplots(figsize=(max(14, len(positions) * 0.09), 1.5))
    im = ax.imshow(vals, aspect="auto", cmap="YlOrRd", interpolation="none",
                   vmin=0, vmax=vmax_global)
    _xticks(ax, positions)
    ax.set_yticks([])
    ax.set_title(f"Patch {patch} — profil interface S1 (% ASA buried absolu, équitable C70)",
                 fontsize=9)
    ax.set_xlabel("Position canonique MAFFT (actine)", fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.8, label="% ASA buried moyen")
    plt.tight_layout()
    fig.savefig(OUT1 / f"{patch}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
print(f"profils : {OUT1}")

# ── 2) décomposition par sous-cluster C70 (toutes positions) ──────────────────
df = pd.read_csv(ROOT / "data/filtered/filtered_all_data.csv", low_memory=False)
di = pd.read_csv(ROOT / "data/filtered/details/1.interactions.csv")[
    ["interaction_id", "chain_A_id", "chain_B_id"]]
res = pd.read_csv(ROOT / "data/filtered/details/3.interface_residues.csv")
res["canon"] = pd.to_numeric(res["residue_number_canon_mafft"], errors="coerce")
res["basa"] = pd.to_numeric(
    res["buried_ASA_percent"].astype(str).str.replace("%", "", regex=False), errors="coerce")

m = df.merge(di, left_on=["subunit_1", "subunit_2"],
             right_on=["chain_A_id", "chain_B_id"], how="left")
m = m[m["s1_binding_site_cluster_data_70"].notna()]

MIN_IIDS = 2
saved = 0
for patch, g in m.groupby("s1_binding_site_cluster_data_70"):
    c70_iids, iid_chain, c70_prot = {}, {}, {}
    for _, r in g.iterrows():
        if pd.isna(r.interaction_id) or pd.isna(r.cluster_data_70):
            continue
        iid = int(r.interaction_id)
        _c70 = str(r.cluster_data_70)
        c70_iids.setdefault(_c70, set()).add(iid)
        iid_chain[iid] = r.subunit_1
        _pt = r.get("subunit_2_title")
        if isinstance(_pt, str) and _pt.strip():
            c70_prot.setdefault(_c70, set()).add(_pt.strip())
    c70_iids = {c: i for c, i in c70_iids.items() if len(i) >= MIN_IIDS}
    if not c70_iids:
        continue
    rows, labels = [], []
    for c70, iids in sorted(c70_iids.items(), key=lambda x: -len(x[1])):
        chains = {iid_chain[i] for i in iids}
        sub = res[(res.interaction_id.isin(iids)) & (res.chain.isin(chains))]
        if sub.empty:
            continue
        prof = (sub.groupby(["interaction_id", "canon"])["basa"].max()
                .groupby(level="canon").sum() / len(iids))
        prof = prof.reindex(positions, fill_value=0).fillna(0).values.astype(float)
        rows.append(prof)
        # protéines partenaires de ce C70 (tronquées si trop nombreuses)
        _prots = sorted(c70_prot.get(c70, set()))
        _pstr = ", ".join(_prots)
        if len(_pstr) > 55:
            _pstr = _pstr[:54] + "…"
        _lab = f"C70={c70} (n={len(iids)})"
        if _pstr:
            _lab += f"\n{_pstr}"
        labels.append(_lab)
    if not rows:
        continue
    mat = np.array(rows); nr = len(rows)
    vmax = mat.max() if mat.max() > 0 else 1.0
    fig, ax = plt.subplots(figsize=(max(16, len(positions) * 0.09), max(2.5, nr * 0.7 + 1.4)))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", interpolation="none", vmin=0, vmax=vmax)
    _xticks(ax, positions)
    ax.set_yticks(range(nr)); ax.set_yticklabels(labels, fontsize=11)
    ax.set_title(f"Patch {patch} — décomposition par sous-cluster C70", fontsize=9)
    ax.set_xlabel("Position canonique MAFFT (actine)", fontsize=8)
    plt.colorbar(im, ax=ax, shrink=0.8, label="% ASA buried moyen / interaction")
    plt.tight_layout()
    fig.savefig(OUT2 / f"{patch}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    saved += 1
print(f"décompositions C70 : {saved} -> {OUT2}")
