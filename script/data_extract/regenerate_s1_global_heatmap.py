#!/usr/bin/env python
"""
Régénère la HEATMAP GLOBALE S1 (HOMO + HÉTÉRO).

- profils recalculés depuis les données brutes en séparant proprement
  HOMO (s2 = actine) et HÉTÉRO (s2 = ABP) — un patch "mixed" apparaît donc
  correctement dans les DEUX sections ;
- clusters triés NUMÉRIQUEMENT (6685_0, 6685_1, 6685_2 …) pour s'y retrouver ;
- labels de lignes plus gros, pas d'espace blanc au-dessus.

Sorties (écrase) :
  data/visualisations/actin_s1_all_equitable_heatmap.png   (relatif : ligne / son max)
  data/visualisations/actin_s1_heatmap_absolute.png        (absolu : % ASA buried)
"""
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

ROOT = Path(__file__).resolve().parents[2]
VIS = ROOT / "data/visualisations"
FS_ROW, FS_TITLE, FS_XTICK, FS_AX, FS_SUP = 9, 12, 8, 11, 13

# axe des positions = colonnes du CSV équitable (toutes positions canoniques)
area = pd.read_csv(ROOT / "data/filtered/actin_s1_canon_area_by_cluster.csv", index_col="patch")
positions = [int(p) for p in area.columns]
pos_index = {p: i for i, p in enumerate(positions)}

df = pd.read_csv(ROOT / "data/filtered/filtered_all_data.csv", low_memory=False)
di = pd.read_csv(ROOT / "data/filtered/details/1.interactions.csv")[
    ["interaction_id", "chain_A_id", "chain_B_id"]]
res = pd.read_csv(ROOT / "data/filtered/details/3.interface_residues.csv")
res["canon"] = pd.to_numeric(res["residue_number_canon_mafft"], errors="coerce")
res["basa"] = pd.to_numeric(
    res["buried_ASA_percent"].astype(str).str.replace("%", "", regex=False), errors="coerce")

m = df.merge(di, left_on=["subunit_1", "subunit_2"],
             right_on=["chain_A_id", "chain_B_id"], how="left")
m = m[m["s1_binding_site_cluster_data_70"].notna() & m["interaction_id"].notna()]


def _profile(sub_rows):
    """Profil équitable C70 (moyenne par C70 puis moyenne des C70) sur l'axe positions."""
    by_c70 = defaultdict(dict)   # c70 -> {iid: actin_chain}
    for _, r in sub_rows.iterrows():
        if pd.isna(r.cluster_data_70):
            continue
        by_c70[str(r.cluster_data_70)][int(r.interaction_id)] = r.subunit_1
    if not by_c70:
        return None
    c70_profiles = []
    for c70, iid_ch in by_c70.items():
        iids = list(iid_ch)
        chains = set(iid_ch.values())
        s = res[(res.interaction_id.isin(iids)) & (res.chain.isin(chains))]
        if s.empty:
            continue
        prof = (s.groupby(["interaction_id", "canon"])["basa"].max()
                .groupby(level="canon").sum() / len(iids))
        v = np.zeros(len(positions))
        for canon, val in prof.items():
            if canon in pos_index:
                v[pos_index[int(canon)]] = val
        c70_profiles.append(v)
    if not c70_profiles:
        return None
    return np.mean(c70_profiles, axis=0)   # équitable entre C70


homo_prof, hetero_prof = {}, {}
for patch, g in m.groupby("s1_binding_site_cluster_data_70"):
    g_homo = g[g["s2_actine"].fillna(False).astype(bool)]
    g_het = g[~g["s2_actine"].fillna(False).astype(bool)]
    if len(g_homo):
        p = _profile(g_homo)
        if p is not None:
            homo_prof[str(patch)] = p
    if len(g_het):
        p = _profile(g_het)
        if p is not None:
            hetero_prof[str(patch)] = p


def _numkey(patch):
    try:
        return int(patch.split("_")[1])
    except (IndexError, ValueError):
        return 10 ** 9


homo = sorted(homo_prof, key=_numkey)
hetero = sorted(hetero_prof, key=_numkey)
print(f"HOMO : {len(homo)} clusters · HÉTÉRO : {len(hetero)} clusters")
print("6685_3 dans HOMO ?", "6685_3" in homo, "| dans HÉTÉRO ?", "6685_3" in hetero)


def _mat(patches, prof, relative):
    if not patches:
        return np.empty((0, len(positions)))
    mm = np.array([prof[p] for p in patches], dtype=float)
    if relative:
        mx = mm.max(axis=1, keepdims=True)
        mm = np.divide(mm, mx, out=np.zeros_like(mm), where=mx > 0)
    return mm


def _draw(relative, vmax, cbar_label, title, out_path):
    mh = _mat(homo, homo_prof, relative)
    me = _mat(hetero, hetero_prof, relative)
    nh, ne = len(homo), len(hetero)
    row_h = 0.26
    fig_h = (nh + ne) * row_h + 1.8
    fig_w = max(14, len(positions) * 0.09)
    fig = plt.figure(figsize=(fig_w, fig_h))
    _title_frac = min(0.06, 1.2 / fig_h)
    gs = gridspec.GridSpec(2, 2, height_ratios=[max(nh, 1), max(ne, 1)],
                           width_ratios=[1, 0.02], hspace=0.06, wspace=0.02,
                           top=1 - _title_frac, bottom=0.02, left=0.09, right=0.93)
    ax_h = fig.add_subplot(gs[0, 0]); ax_e = fig.add_subplot(gs[1, 0]); ax_cb = fig.add_subplot(gs[:, 1])
    step = max(1, len(positions) // 55)
    xt = list(range(0, len(positions), step))
    im = None
    if nh:
        im = ax_h.imshow(mh, aspect="auto", cmap="YlOrRd", interpolation="none", vmin=0, vmax=vmax)
        ax_h.set_xticks(xt); ax_h.set_xticklabels([])
        ax_h.set_yticks(range(nh)); ax_h.set_yticklabels(homo, fontsize=FS_ROW)
        ax_h.set_title(f"HOMO — actine / actine ({nh} clusters)", fontsize=FS_TITLE, loc="left", pad=3)
    else:
        ax_h.axis("off")
    if ne:
        im = ax_e.imshow(me, aspect="auto", cmap="YlOrRd", interpolation="none", vmin=0, vmax=vmax)
        ax_e.set_xticks(xt)
        ax_e.set_xticklabels([str(positions[i]) for i in xt], fontsize=FS_XTICK, rotation=90)
        ax_e.set_xlabel("Position canonique MAFFT (actine)", fontsize=FS_AX)
        ax_e.set_yticks(range(ne)); ax_e.set_yticklabels(hetero, fontsize=FS_ROW)
        ax_e.set_title(f"HÉTÉRO — actine / ABP ({ne} clusters)", fontsize=FS_TITLE, loc="left", pad=3)
    else:
        ax_e.axis("off")
    if im is not None:
        plt.colorbar(im, cax=ax_cb, label=cbar_label)
    fig.suptitle(title, fontsize=FS_SUP, y=1 - _title_frac * 0.35)
    fig.savefig(out_path, dpi=150, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"écrit : {out_path}")


_draw(True, 1.0, "Score normalisé (0–1, max par cluster)",
      "Heatmap globale — tous clusters S1 binding site (relatif : ligne / son max, équitable C70)",
      VIS / "actin_s1_all_equitable_heatmap.png")
vmax_abs = max(max((p.max() for p in homo_prof.values()), default=0),
               max((p.max() for p in hetero_prof.values()), default=0), 1.0)
_draw(False, vmax_abs, "% ASA buried moyen (équitable C70)",
      "Heatmap globale — tous clusters S1 binding site (absolu : % ASA buried, équitable C70)",
      VIS / "actin_s1_heatmap_absolute.png")
