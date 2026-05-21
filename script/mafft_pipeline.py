"""
Pipeline MAFFT — alignement des séquences par cluster et mise à jour des tables 3 et 4.

Étapes :
  ETAPE 1/4 — Chargement des données
  ETAPE 2/4 — Création des fichiers FASTA
  ETAPE 3/4 — Alignements MAFFT
  ETAPE 4/4 — Mapping canon_mafft → tables 3 et 4

Usage :
  python -m script.mafft_pipeline
"""

import hashlib
import shutil
import subprocess
import sys
import urllib.request
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def _fasta_hash(path: Path) -> str:
    """SHA-1 du contenu d'un fichier FASTA (détecte tout changement de séquences)."""
    return hashlib.sha1(path.read_bytes()).hexdigest()


_P60709_CACHE = Path("data/P60709_ref.fasta")
REF_ID = "P60709_ref"

def _load_ref_seq() -> SeqRecord | None:
    """Charge P60709 (beta-actine humaine) depuis le cache local ou UniProt."""
    if _P60709_CACHE.exists():
        records = list(SeqIO.parse(str(_P60709_CACHE), "fasta"))
        if records:
            return SeqRecord(Seq(str(records[0].seq)), id=REF_ID, description="")
    try:
        url = "https://rest.uniprot.org/uniprotkb/P60709.fasta"
        with urllib.request.urlopen(url, timeout=15) as r:
            fasta_data = r.read().decode()
        _P60709_CACHE.parent.mkdir(parents=True, exist_ok=True)
        _P60709_CACHE.write_text(fasta_data)
        records = list(SeqIO.parse(StringIO(fasta_data), "fasta"))
        if records:
            print(f"  P60709 téléchargé depuis UniProt ({len(records[0].seq)} aa) → {_P60709_CACHE}")
            return SeqRecord(Seq(str(records[0].seq)), id=REF_ID, description="")
    except Exception as e:
        print(f"  Impossible de charger P60709 : {e} — alignements sans référence fixe")
    return None


def _hash_path(aln_path: Path) -> Path:
    return aln_path.with_suffix(".aln.hash")

# ── Chemins ───────────────────────────────────────────────────────────────
DATA_DIR = Path("data")
DETAILS_PATH = DATA_DIR / "filtered" / "details"
EXPORT_DIR = DATA_DIR / "exports"
FASTA_DIR = EXPORT_DIR / "fasta"
ALN_DIR = DATA_DIR / "alignments"


def run_mafft(fasta_path: Path, aln_path: Path) -> bool:
    """Aligne un FASTA avec MAFFT via le binaire disponible dans le PATH."""
    mafft_bin = shutil.which("mafft")

    if mafft_bin is None:
        print(
            "MAFFT introuvable dans le PATH. Vérifie avec `mafft --version`.",
            file=sys.stderr,
        )
        return False

    try:
        with open(aln_path, "w") as out:
            result = subprocess.run(
                [mafft_bin, "--auto", str(fasta_path)],
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                timeout=120,
            )

        if result.returncode != 0:
            if result.stderr:
                print(result.stderr, file=sys.stderr)
            return False

        return aln_path.exists() and aln_path.stat().st_size > 0

    except subprocess.TimeoutExpired:
        print(f"Timeout MAFFT pour {fasta_path.name}", file=sys.stderr)
        return False


def build_seq_to_col_mapping(aligned_seq: str) -> dict:
    """Position séquence (1-based) → colonne alignement (1-based)."""
    mapping, seq_pos = {}, 0
    for col, aa in enumerate(aligned_seq, start=1):
        if aa != "-":
            seq_pos += 1
            mapping[seq_pos] = col
    return mapping


def main():
    for d in (FASTA_DIR, ALN_DIR):
        d.mkdir(parents=True, exist_ok=True)

    # Invalide le sentinel .done si les données sources ont changé depuis
    done_sentinel = ALN_DIR / ".done"
    source_files  = [
        DATA_DIR / "filtered/filtered_all_data.csv",
        DATA_DIR / "filtered/filtered_summary.csv",
    ]
    if done_sentinel.exists():
        done_mtime = done_sentinel.stat().st_mtime
        if any(f.exists() and f.stat().st_mtime > done_mtime for f in source_files):
            done_sentinel.unlink()
            print("  Sources modifiées depuis le dernier alignement — recalcul forcé")

    # ── ETAPE 1/4 : Chargement ────────────────────────────────────────────
    print("ETAPE 1/4 — Chargement des données")
    df_all = pd.read_csv(DATA_DIR / "filtered/filtered_all_data.csv")
    df_summary = pd.read_csv(DATA_DIR / "filtered/filtered_summary.csv")
    df_merged = df_all.merge(
        df_summary,
        left_on=["subunit_1", "subunit_2"],
        right_on=["Protein ID", "Interactor ID"],
        how="left",
    )
    print(f"  {len(df_merged)} interactions chargées")

    # ── ETAPE 2/4 : FASTA ─────────────────────────────────────────────────
    print("ETAPE 2/4 — Création des fichiers FASTA")
    ref_record = _load_ref_seq()
    if ref_record:
        print(f"  Référence P60709 ({len(ref_record.seq)} aa) ajoutée à tous les clusters")

    subunit_title = {}
    for _, row in df_merged.iterrows():
        for scol, tcol in [("subunit_1", "subunit_1_title"), ("subunit_2", "subunit_2_title")]:
            t = str(row.get(tcol, "")).strip()
            if t and t != "nan":
                subunit_title[row[scol]] = t

    seq_entries = []
    for _, row in df_merged.iterrows():
        for seq_col, sub_col, clust_col in [
            ("s1_sequence", "subunit_1", "s1_sequence_cluster_70"),
            ("s2_sequence", "subunit_2", "s2_sequence_cluster_70"),
        ]:
            seq = str(row.get(seq_col, "")).strip()
            cluster = row.get(clust_col)
            subunit = row[sub_col]
            if seq and seq != "nan" and pd.notna(cluster):
                seq_entries.append((int(cluster), seq, subunit))

    cluster_data = defaultdict(lambda: defaultdict(list))
    for cluster_id, seq, subunit in seq_entries:
        if subunit not in cluster_data[cluster_id][seq]:
            cluster_data[cluster_id][seq].append(subunit)

    # seq_key_to_subunits[(cluster_id, first_subunit)] = [all subunits with same seq]
    # Utilisé après alignement pour propager le mapping à tous les subunits identiques.
    seq_key_to_subunits: dict = {}
    ref_seq_str = str(ref_record.seq) if ref_record else None
    n_files = 0
    for cluster_id, seq_to_subunits in sorted(cluster_data.items()):
        records = []
        for seq, subunits in seq_to_subunits.items():
            first = subunits[0]
            seq_key_to_subunits[(cluster_id, first)] = subunits
            records.append(SeqRecord(Seq(seq), id=first, description=""))

        if len(records) >= 2:
            # P60709 en tête comme ancre de référence globale (si non déjà présente)
            # puis les autres triées par longueur décroissante pour couvrir les positions >375
            cluster_seqs = {str(r.seq) for r in records}
            if ref_record and ref_seq_str not in cluster_seqs:
                records = [ref_record] + sorted(records, key=lambda r: len(r.seq), reverse=True)
            else:
                records.sort(key=lambda r: len(r.seq), reverse=True)
            path = FASTA_DIR / f"cluster_{cluster_id}.fasta"
            SeqIO.write(records, path, "fasta")
            n_files += 1
    print(f"  {n_files} fichiers FASTA → {FASTA_DIR}/")

    # ── ETAPE 3/4 : Alignements ───────────────────────────────────────────
    print("ETAPE 3/4 — Alignements MAFFT")
    fastas = sorted(FASTA_DIR.glob("*.fasta"))
    ok = skipped = 0
    for i, fasta in enumerate(fastas, 1):
        aln = ALN_DIR / fasta.name.replace(".fasta", ".aln")
        hp = _hash_path(aln)
        current_hash = _fasta_hash(fasta)
        if aln.exists() and aln.stat().st_size > 0 and hp.exists() and hp.read_text().strip() == current_hash:
            ok += 1
            skipped += 1
            continue
        print(f"  {i}/{len(fastas)} {fasta.stem}", flush=True)
        if run_mafft(fasta, aln):
            hp.write_text(current_hash)
            ok += 1
    print(f"  Alignements : {ok}/{len(fastas)}  ({skipped} inchangés, ignorés)")

    # ── ETAPE 4/4 : Mapping canon_mafft ───────────────────────────────────
    print("ETAPE 4/4 — Mapping canon_mafft → tables 3 et 4")

    subunit_to_mapping = {}
    for aln_path in sorted(ALN_DIR.glob("*.aln")):
        cluster_id_str = aln_path.stem.replace("cluster_", "")
        try:
            cluster_id_int = int(cluster_id_str)
        except ValueError:
            cluster_id_int = None
        try:
            alignment = AlignIO.read(str(aln_path), "fasta")
        except Exception:
            continue
        for record in alignment:
            mapping = build_seq_to_col_mapping(str(record.seq))
            # record.id = premier subunit (ID court, non tronqué)
            # On propage le mapping à TOUS les subunits partageant la même séquence
            all_subs = [record.id]
            if cluster_id_int is not None:
                all_subs = seq_key_to_subunits.get(
                    (cluster_id_int, record.id), [record.id]
                )
            for subunit in all_subs:
                subunit_to_mapping[subunit] = mapping

    seq_by_subunit = {}
    for _, row in df_merged.iterrows():
        seq_by_subunit[row["subunit_1"]] = str(row["s1_sequence"])
        seq_by_subunit[row["subunit_2"]] = str(row["s2_sequence"])

    n_identity = 0
    for subunit, seq in seq_by_subunit.items():
        if subunit not in subunit_to_mapping and seq and seq.upper() != "NAN":
            subunit_to_mapping[subunit] = {
                i: i for i in range(1, len(seq) + 1)}
            n_identity += 1

    print(
        f"  Mappings alignements : {len(subunit_to_mapping) - n_identity} subunits")
    print(f"  Mappings identité    : {n_identity} subunits")

    df_res3 = pd.read_csv(DETAILS_PATH / "3.interface_residues.csv")
    df_res4 = pd.read_csv(DETAILS_PATH / "4.inter-residue_contacts.csv")
    # Supprimer les éventuelles colonnes canon_mafft d'un run précédent
    df_res4.drop(columns=[c for c in df_res4.columns if "canon" in c.lower()], inplace=True)
    df_res3.drop(columns=[c for c in df_res3.columns if "canon" in c.lower()], inplace=True)
    df_int1 = pd.read_csv(DETAILS_PATH / "1.interactions.csv")
    chain_a_map = dict(zip(df_int1["interaction_id"], df_int1["chain_A_id"]))
    chain_b_map = dict(zip(df_int1["interaction_id"], df_int1["chain_B_id"]))

    def map_t4(row):
        subunit = chain_a_map.get(int(row["interaction_id"]))
        mapping = subunit_to_mapping.get(subunit, {})
        seq_pos = row["residue_A_sequence"]
        if pd.isna(seq_pos):
            return pd.NA
        return mapping.get(int(seq_pos), pd.NA)

    def map_t4_b(row):
        subunit = chain_b_map.get(int(row["interaction_id"]))
        mapping = subunit_to_mapping.get(subunit, {})
        seq_pos = row["residue_B_sequence"]
        if pd.isna(seq_pos):
            return pd.NA
        return mapping.get(int(seq_pos), pd.NA)

    def map_t3(row):
        iid = int(row["interaction_id"])
        chain = row["chain"]
        # Accepte chain_A ET chain_B (utile pour les homodimères actine-actine)
        if chain == chain_a_map.get(iid) or chain == chain_b_map.get(iid):
            subunit = chain
        else:
            return pd.NA
        mapping = subunit_to_mapping.get(subunit, {})
        seq_pos = row["residue_number_sequence"]
        if pd.isna(seq_pos):
            return pd.NA
        return mapping.get(int(seq_pos), pd.NA)

    df_res4["residue_A_canon_mafft"] = df_res4.apply(map_t4, axis=1)
    df_res4["residue_B_canon_mafft"] = df_res4.apply(map_t4_b, axis=1)
    df_res3["residue_number_canon_mafft"] = df_res3.apply(map_t3, axis=1)

    n4a = df_res4["residue_A_canon_mafft"].notna().sum()
    n4b = df_res4["residue_B_canon_mafft"].notna().sum()
    n3 = df_res3["residue_number_canon_mafft"].notna().sum()
    print(f"  Table 4 : {n4a}/{len(df_res4)} contacts avec residue_A_canon_mafft  |  {n4b}/{len(df_res4)} avec residue_B_canon_mafft")
    print(f"  Table 3 : {n3}/{len(df_res3)} résidus avec residue_number_canon_mafft")

    # Enrichissement asa_pct_A/B : jointure table 3 → table 4
    # (table 4 brute téléchargée par get_interaction_details ne contient pas ces colonnes)
    _df3_asa = df_res3.copy()
    _df3_asa["buried_ASA_percent"] = pd.to_numeric(
        _df3_asa["buried_ASA_percent"].astype(str).str.replace("%", "", regex=False),
        errors="coerce",
    )
    _df3_dedup = (
        _df3_asa.groupby(["interaction_id", "chain", "residue_number_sequence"], sort=False)
        ["buried_ASA_percent"].max().reset_index()
    )
    df_res4 = df_res4.drop(columns=["asa_pct_A", "asa_pct_B"], errors="ignore")
    df_res4 = df_res4.merge(
        _df3_dedup.rename(columns={"chain": "chain_A_id",
                                   "residue_number_sequence": "residue_A_sequence",
                                   "buried_ASA_percent": "asa_pct_A"}),
        on=["interaction_id", "chain_A_id", "residue_A_sequence"], how="left",
    )
    df_res4 = df_res4.merge(
        _df3_dedup.rename(columns={"chain": "chain_B_id",
                                   "residue_number_sequence": "residue_B_sequence",
                                   "buried_ASA_percent": "asa_pct_B"}),
        on=["interaction_id", "chain_B_id", "residue_B_sequence"], how="left",
    )
    ok_A = df_res4["asa_pct_A"].notna().sum()
    ok_B = df_res4["asa_pct_B"].notna().sum()
    print(f"  asa_pct_A/B : {ok_A}/{len(df_res4)} renseignés")

    df_res4.to_csv(DETAILS_PATH / "4.inter-residue_contacts.csv", index=False)
    df_res3.to_csv(DETAILS_PATH / "3.interface_residues.csv", index=False)
    print("  Tables 3 et 4 mises à jour ✓")

    (ALN_DIR / ".done").touch()
    print("  Sentinel .done créé ✓")


if __name__ == "__main__":
    main()
