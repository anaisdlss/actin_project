#!/usr/bin/env python
"""
Chimie de l'interface : que partagent des ABP différentes au même site, au-delà
de la structure ?
  - composition physicochimique des résidus de contact côté ABP (chargés +/-, hydrophobes,
    polaires, aromatiques) + charge nette
  - charge nette du patch d'actine ciblé (complémentarité électrostatique ?)
  - nature physique du contact : aire, nb contacts, liaisons H, ponts salins (1.interactions.csv)

Sortie : data/exports/abp_site_domain/interface_chemistry.csv
"""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"

rep = pd.read_csv(OUT / "abp_representatives.csv")
res = pd.read_csv(ROOT / "data/filtered/details/3.interface_residues.csv")
di  = pd.read_csv(ROOT / "data/filtered/details/1.interactions.csv")
df  = pd.read_csv(ROOT / "data/filtered/filtered_all_data.csv", low_memory=False)

POS = set("KR"); NEG = set("DE"); HYD = set("AVLIMFWPC"); ARO = set("FWY")
POL = set("STNQYHGC")  # polaires/neutres (H compté à part en charge)


def classify(names):
    """names = liste de codes 1-lettre. Renvoie compo % + charge nette."""
    n = len(names)
    if n == 0:
        return None
    pos = sum(c in POS for c in names)
    neg = sum(c in NEG for c in names)
    hyd = sum(c in HYD for c in names)
    aro = sum(c in ARO for c in names)
    return {
        "n_res": n,
        "pct_charge_pos": round(100 * pos / n),
        "pct_charge_neg": round(100 * neg / n),
        "pct_hydrophobe": round(100 * hyd / n),
        "pct_aromatique": round(100 * aro / n),
        "charge_nette": pos - neg,
    }


THREE2ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}
def one(nm):
    nm = str(nm).strip().upper()
    return THREE2ONE.get(nm, nm if len(nm) == 1 else "X")


# stats de contact par interaction
di_idx = di.set_index("interaction_id")[["interface_area", "num_contacts", "num_hbonds", "num_salt_bridges"]]

# actine = subunit_1 par interaction (pour le patch actine)
df_actin = df[(df.s1_actine) & (~df.s2_actine)][["subunit_1", "subunit_2"]]
abp_to_actinchain = {}
for _, r in df_actin.iterrows():
    abp_to_actinchain.setdefault(r.subunit_2, r.subunit_1)

rows = []
for _, r in rep.iterrows():
    iid = int(r.interaction_id)
    abp_full = f"{r.pdb}_{r.abp_chain}"
    # côté ABP
    sa = res[(res.interaction_id == iid) & (res.chain == abp_full)]
    abp_names = [one(x) for x in sa.residue_name]
    chem = classify(abp_names)
    if chem is None:
        continue
    # côté actine (patch)
    actin_chain = abp_to_actinchain.get(r.abp_subunit, None)
    actin_charge = None
    if actin_chain is not None:
        sact = res[(res.interaction_id == iid) & (res.chain == actin_chain)]
        an = [one(x) for x in sact.residue_name]
        if an:
            actin_charge = sum(c in POS for c in an) - sum(c in NEG for c in an)
    # stats contact
    import re as _re
    def _num(v):
        if v is None or pd.isna(v):
            return None
        m = _re.search(r"-?\d+\.?\d*", str(v))
        return float(m.group()) if m else None
    stt = di_idx.loc[iid] if iid in di_idx.index else None
    rows.append({
        "abp_title": r.abp_title,
        **chem,
        "charge_actine_patch": actin_charge,
        "interface_area": (round(_num(stt["interface_area"])) if stt is not None and _num(stt["interface_area"]) is not None else None),
        "num_contacts": (int(_num(stt["num_contacts"])) if stt is not None and _num(stt["num_contacts"]) is not None else None),
        "num_hbonds": (int(_num(stt["num_hbonds"])) if stt is not None and _num(stt["num_hbonds"]) is not None else None),
        "num_salt_bridges": (int(_num(stt["num_salt_bridges"])) if stt is not None and _num(stt["num_salt_bridges"]) is not None else None),
    })

out = pd.DataFrame(rows).sort_values("abp_title")
out.to_csv(OUT / "interface_chemistry.csv", index=False)
print(f"{len(out)} ABP -> interface_chemistry.csv\n")
print(f"{'ABP':38s} {'+%':>3} {'-%':>3} {'hyd%':>4} {'aro%':>4} {'chrg':>4} {'actNet':>6} {'aire':>5} {'Hb':>3} {'salt':>4}")
for _, r in out.iterrows():
    print(f"{str(r.abp_title)[:38]:38s} {r.pct_charge_pos:>3} {r.pct_charge_neg:>3} "
          f"{r.pct_hydrophobe:>4} {r.pct_aromatique:>4} {r.charge_nette:>4} "
          f"{str(r.charge_actine_patch):>6} {str(r.interface_area):>5} "
          f"{str(r.num_hbonds):>3} {str(r.num_salt_bridges):>4}")
