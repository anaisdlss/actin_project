#!/usr/bin/env python
"""
Etape 4 — Annotation des domaines (Pfam / InterPro) de chaque ABP.

Pour chaque ABP representatif : PDB+chaine -> UniProt (SIFTS PDBe) -> domaines (API InterPro).
La chaine actine eventuellement co-mappee est filtree.

Sortie : data/exports/abp_site_domain/abp_interpro.csv
         (abp_title, uniprot, pfam_domains, interpro_domains)
"""
import json
import re
import time
from pathlib import Path
import urllib.request
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/exports/abp_site_domain"
CACHE = OUT / "_api_cache"
CACHE.mkdir(parents=True, exist_ok=True)

SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb}"
INTERPRO = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}/?page_size=200"
ACTIN_NAME = re.compile(r"^ACT[A-Z0-9]?_|ACTIN", re.I)


def get_json(url, key):
    f = CACHE / (re.sub(r"[^A-Za-z0-9]+", "_", key) + ".json")
    if f.exists():
        return json.loads(f.read_text())
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                d = json.load(r)
            f.write_text(json.dumps(d))
            time.sleep(0.3)
            return d
        except Exception as e:
            if attempt == 2:
                print(f"  ! echec {key}: {e}")
                return None
            time.sleep(2)


def uniprot_for_chain(pdb, chain):
    """UniProt accessions (non-actine) mappes a une chaine."""
    d = get_json(SIFTS.format(pdb=pdb), f"sifts_{pdb}")
    if not d or pdb not in d:
        return []
    accs = []
    for acc, v in d[pdb].get("UniProt", {}).items():
        if ACTIN_NAME.search(v.get("name", "")):
            continue
        if any(m.get("chain_id") == chain for m in v.get("mappings", [])):
            accs.append((acc, v.get("name", "")))
    return accs


# on retient les types qui caracterisent le repliement / la famille de l'ABP
KEEP_TYPES = {"domain", "family", "repeat", "homologous_superfamily"}


def domains_for_acc(acc):
    d = get_json(INTERPRO.format(acc=acc), f"interpro_{acc}")
    pfam, ipr = {}, {}
    if not d:
        return pfam, ipr
    for r in d.get("results", []):
        md = r["metadata"]
        if md.get("type") not in KEEP_TYPES:
            continue
        if md["source_database"] == "pfam":
            pfam[md["accession"]] = md.get("name")
        elif md["source_database"] == "interpro":
            ipr[md["accession"]] = md.get("name")
    return pfam, ipr


rep = pd.read_csv(OUT / "abp_representatives.csv")
rows = []
for _, r in rep.iterrows():
    accs = uniprot_for_chain(r.pdb, r.abp_chain)
    pfam, ipr = {}, {}
    for acc, _name in accs:
        p, i = domains_for_acc(acc)
        pfam.update(p)
        ipr.update(i)
    rows.append(dict(
        abp_title=r.abp_title, pdb=r.pdb, abp_chain=r.abp_chain,
        uniprot=";".join(a for a, _ in accs),
        pfam_acc=";".join(sorted(pfam)),
        pfam_domains=" | ".join(sorted(set(pfam.values()))),
        interpro_acc=";".join(sorted(ipr)),
        interpro_domains=" | ".join(sorted(set(ipr.values()))),
    ))
    print(f"{r.abp_title[:40]:40s} {r.pdb}_{r.abp_chain}  pfam={list(pfam.values())}")

pd.DataFrame(rows).to_csv(OUT / "abp_interpro.csv", index=False)
print(f"\necrit : {OUT/'abp_interpro.csv'}")
miss = [row["abp_title"] for row in rows if not row["pfam_domains"] and not row["interpro_domains"]]
if miss:
    print(f"sans domaine ({len(miss)}): {miss}")
