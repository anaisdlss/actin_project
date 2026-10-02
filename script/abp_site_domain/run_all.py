#!/usr/bin/env python
"""
Étape 10 du pipeline — Analyses structurales ABP (reproductible).

Enchaîne tout ce qui alimente le panneau « Comparaison structurale des ABP »
du dashboard : table ABP×site, extraction de chaînes, Foldseek (folds), InterPro,
TM-align interface, familles, empreinte actine, structure secondaire, chimie.

Lancer :  pixi run python -m script.abp_site_domain.run_all
Sorties : data/exports/abp_site_domain/*.csv/tsv (+ figures)
"""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SD = ROOT / "script/abp_site_domain"
OUT = ROOT / "data/exports/abp_site_domain"
CHAINS = OUT / "abp_chains"
PY = sys.executable


_failed = []


def _py(script):
    # fail-soft : un sous-script qui échoue (réseau transitoire, données…) est
    # loggé et on CONTINUE — l'étape ne doit pas tomber pour une brique annexe.
    # La complétude réelle est jugée sur la sortie clé (familles.csv) à la fin.
    print(f"  → {script}", flush=True)
    rc = subprocess.run([PY, str(SD / script)], cwd=ROOT).returncode
    if rc != 0:
        print(f"  ! {script} a échoué (code {rc}) — on continue", flush=True)
        _failed.append(script)
    return rc == 0


def _sh(cmd, desc):
    print(f"  → {desc}", flush=True)
    try:
        rc = subprocess.run(cmd, cwd=ROOT).returncode
    except FileNotFoundError as e:
        print(f"  ! {desc} : binaire introuvable ({e}) — on continue", flush=True)
        _failed.append(desc)
        return False
    if rc != 0:
        print(f"  ! {desc} a échoué (code {rc}) — on continue", flush=True)
        _failed.append(desc)
    return rc == 0


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    marker = OUT / "pipeline_incomplete.json"
    marker.write_text(json.dumps({"state": "running"}))
    print("== Étape 10 : analyses structurales ABP ==")
    _py("01_build_table.py")
    _py("02_extract_chains.py")

    # Foldseek : clustering structural (folds) + recherche all-vs-all (TM protéine entière)
    _sh(["foldseek", "easy-cluster", str(CHAINS), str(OUT / "fold_cluster"),
         str(OUT / "foldseek_tmp"), "--alignment-type", "2", "-c", "0.6",
         "--cov-mode", "0", "--tmscore-threshold", "0.5"],
        "foldseek easy-cluster (folds)")
    _sh(["foldseek", "easy-search", str(CHAINS), str(CHAINS),
         str(OUT / "whole_pairs_all.tsv"), str(OUT / "whole_tmp"),
         "--alignment-type", "1", "--tmalign-fast", "0", "-e", "inf", "-c", "0",
         "--max-seqs", "200", "--format-output",
         "query,target,qtmscore,ttmscore,rmsd,fident,alnlen"],
        "foldseek easy-search (TM protéine entière)")

    _py("03_interpro.py")          # domaines (API InterPro, cache dans _api_cache/)
    _py("04_crosstab.py")          # abp_master.csv + croisement
    _py("05_extract_interface.py")  # interfaces (PAD=0) pour le sweep TM
    _py("08_interface_tm_sweep.py")  # TM interface→entière (tmtools)
    _py("13_folddisco_interface.py")  # motifs d'interface discontinus (FoldDisco)
    _py("14_folddisco_discovery.py")  # découverte PDB+AlphaFold (API FoldDisco, réseau, lent ; reprise auto)
    _py("14b_resolve_names.py")       # noms de protéines des hits (UniProt + RCSB, cache)
    _py("11_family_regroup.py")    # familles.csv + convergences
    _py("13_actin_footprint_overlap.py")  # empreinte actine
    _py("16_interface_chemistry.py")      # chimie de l'interface
    _py("18_site_determinants.py")        # déterminants (RSA/conservation)

    # Structure secondaire de l'interface (PyMOL dss)
    _sh(["pymol", "-cq", str(SD / "15_interface_secondary_structure.py")],
        "structure secondaire de l'interface (PyMOL)")

    _py("19_interface_domain_family.py")  # domaine Pfam RÉEL au contact + familles d'interface

    # nettoyage des tmp foldseek
    for d in (OUT / "foldseek_tmp", OUT / "whole_tmp"):
        subprocess.run(["rm", "-rf", str(d)], check=False)

    if _failed:
        marker.write_text(json.dumps({"state": "failed", "steps": _failed}, indent=2))
        print(f"Structural analysis incomplete: {', '.join(_failed)}. Previous files do not establish success.", flush=True)
        sys.exit(1)
    familles = OUT / "familles.csv"
    if not familles.exists():
        print("== Étape 10 INCOMPLÈTE : familles.csv (sortie clé) absent — "
              "relance Run/update pour re-tenter. ==", flush=True)
        sys.exit(1)
    marker.unlink(missing_ok=True)
    print("== Étape 10 terminée ==")


if __name__ == "__main__":
    main()
