#!/usr/bin/env python
"""
Domaine RÉELLEMENT utilisé au contact de l'actine (par ABP), puis re-test des
familles « au niveau de l'interface ».

Idée : un ABP peut avoir plusieurs domaines Pfam, mais à un site il ne touche
l'actine qu'avec une région → un seul domaine est « utilisé ». On mappe les
résidus de contact de l'ABP (numérotation PDB) vers UniProt (via SIFTS), puis on
regarde dans quel domaine Pfam (bornes via InterPro) ils tombent.

Deux ABP au même site sont « même famille d'interface » s'ils lient l'actine par
le MÊME domaine Pfam (même accession), et non plus s'ils partagent n'importe quel
domaine ailleurs dans la séquence.

Sorties : data/exports/abp_site_domain/interface_domain_by_abp.csv
          data/exports/abp_site_domain/interface_family_refined.csv
"""
import json
import re
import time
from itertools import combinations
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
CACHE = OUT / "_api_cache"
CACHE.mkdir(parents=True, exist_ok=True)
RES = ROOT / "data/filtered/details/3.interface_residues.csv"

SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb}"
IPR_POS = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/"
# repli : entrées InterPro (familles / superfamilles homologues, inter-espèces)
# quand Pfam n'annote que des fragments partiels ne couvrant pas les contacts.
IPR_ALL = "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{acc}/"
KEEP_IPR = {"family", "homologous_superfamily", "domain"}


def get_json(url, key):
    p = CACHE / f"{key}.json"
    if p.exists():
        try:
            return json.loads(p.read_text())
        except Exception:
            pass
    try:
        r = requests.get(url, timeout=60, headers={"Accept": "application/json"})
        if r.status_code == 204:
            p.write_text("{}"); return {}
        r.raise_for_status()
        d = r.json()
        p.write_text(json.dumps(d)); time.sleep(0.2)
        return d
    except Exception as e:
        print(f"  ! {key}: {e}")
        return None


def sifts_offset(pdb, chain, acc):
    """Décalage résidu PDB(author) -> UniProt pour (pdb, chain, acc). Renvoie
    (offset, unp_min, unp_max) ou None."""
    d = get_json(SIFTS.format(pdb=pdb.lower()), f"sifts_{pdb.lower()}")
    if not d:
        return None
    up = d.get(pdb.lower(), {}).get("UniProt", {})
    info = up.get(acc)
    if not info:
        return None

    def _off(m):
        if m.get("unp_start") is None:
            return None
        st = m.get("start", {})
        ps = st.get("author_residue_number")
        if ps is None:
            # numéro d'auteur absent → la structure est presque toujours numérotée
            # directement en UniProt → offset 0 (résidu PDB == résidu UniProt).
            return (0, m["unp_start"], m.get("unp_end"))
        return (m["unp_start"] - ps, m["unp_start"], m.get("unp_end"))

    # chaîne correspondante d'abord, sinon 1er mapping exploitable
    for m in info.get("mappings", []):
        if m.get("chain_id") == chain or m.get("struct_asym_id") == chain:
            o = _off(m)
            if o:
                return o
    for m in info.get("mappings", []):
        o = _off(m)
        if o:
            return o
    return None


def _parse_pos(d, keep_types=None):
    out = []
    if not d:
        return out
    for r in d.get("results", []):
        md = r["metadata"]
        if keep_types and md.get("type") not in keep_types:
            continue
        ranges = []
        for prot in r.get("proteins", []):
            for loc in prot.get("entry_protein_locations", []):
                for fr in loc.get("fragments", []):
                    ranges.append((fr["start"], fr["end"]))
        if ranges:
            out.append((md["accession"], md.get("name", ""), ranges))
    return out


def pfam_domains_pos(acc):
    """Domaines Pfam (start,end). Repli sur les entrées InterPro
    (familles/superfamilles homologues) si besoin — fait via la sélection."""
    return _parse_pos(get_json(IPR_POS.format(acc=acc), f"iprpos_{acc}"))


def interpro_domains_pos(acc):
    return _parse_pos(get_json(IPR_ALL.format(acc=acc), f"iprall_{acc}"), KEEP_IPR)


# ── données ───────────────────────────────────────────────────────────────
rep = pd.read_csv(OUT / "abp_representatives.csv")
ip = pd.read_csv(OUT / "abp_interpro.csv")
uni_of = {r.abp_title: (str(r.uniprot).split(";")[0] if pd.notna(r.uniprot) else None)
          for _, r in ip.iterrows()}
res = pd.read_csv(RES, usecols=["interaction_id", "chain", "residue_number_structure"])
res["rn"] = pd.to_numeric(res["residue_number_structure"], errors="coerce")

rows = []
for _, r in rep.drop_duplicates("abp_title").iterrows():
    acc = uni_of.get(r.abp_title)
    chain = f"{r.pdb}_{r.abp_chain}"
    contacts = res[(res.interaction_id == r.interaction_id) & (res.chain == chain)]["rn"].dropna().astype(int)
    domdom = "—"
    dom_name = "—"
    n_dom_prot = 0
    if acc and len(contacts):
        off = sifts_offset(r.pdb, r.abp_chain, acc)
        doms = pfam_domains_pos(acc)
        n_dom_prot = len(doms)
        if off:
            unp = [c + off[0] for c in contacts]

            def _best(dset):
                b, bn, bname = None, 0, "—"
                for pacc, pname, ranges in dset:
                    n = sum(any(s <= u <= e for s, e in ranges) for u in unp)
                    if n > bn:
                        b, bn, bname = pacc, n, pname
                return b, bn, bname
            # 1) Pfam qui couvre le plus de contacts
            best, bestn, dom_name = _best(doms)
            # 2) repli InterPro (superfamille homologue) si aucun Pfam ne couvre
            if not best:
                ib, ibn, iname = _best(interpro_domains_pos(acc))
                if ib:
                    best, dom_name = ib, iname + " (InterPro)"
            if best:
                domdom = best
    rows.append(dict(abp_title=r.abp_title, uniprot=acc, n_contacts=len(contacts),
                     n_domaines_proteine=n_dom_prot,
                     domaine_interface_acc=domdom, domaine_interface=dom_name))

by = pd.DataFrame(rows)
by.to_csv(OUT / "interface_domain_by_abp.csv", index=False)
print("== Domaine utilisé au contact, par ABP ==")
print(by[["abp_title", "n_domaines_proteine", "domaine_interface"]].to_string(index=False))

# ── re-test des familles AU NIVEAU DE L'INTERFACE, par cluster ──────────────
fo = pd.read_csv(OUT / "actin_footprint_overlap.csv")
dib = dict(zip(by.abp_title, by.domaine_interface_acc))
dname = dict(zip(by.abp_title, by.domaine_interface))
ref = []
for _, r in fo.iterrows():
    da, db = dib.get(r.a, "—"), dib.get(r.b, "—")
    same_iface = (da != "—" and da == db)
    ref.append(dict(site=r.site, a=r.a, b=r.b,
                    dom_a=dname.get(r.a, "—"), dom_b=dname.get(r.b, "—"),
                    meme_domaine_interface="oui" if same_iface else "non",
                    ancienne_meme_famille="oui" if r.same_fam else "non",
                    empreinte=round(r.recouvrement_min, 2)))
rf = pd.DataFrame(ref)
rf.to_csv(OUT / "interface_family_refined.csv", index=False)

print("\n== Comparaison ancienne famille (protéine entière) vs domaine d'interface ==")
chg = rf[(rf.ancienne_meme_famille != rf.meme_domaine_interface)]
print(f"paires où la conclusion CHANGE : {len(chg)} / {len(rf)}")
print(rf.groupby(["ancienne_meme_famille", "meme_domaine_interface"]).size().to_string())
print(f"\nécrit : {OUT/'interface_domain_by_abp.csv'} et interface_family_refined.csv")
