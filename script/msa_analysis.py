import streamlit as st
import pandas as pd
import numpy as np
import os
import re
import sys
import shutil
import subprocess
from pathlib import Path as _Path
from collections import Counter as _Counter, defaultdict
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


_CHAIN_ID_RE = re.compile(r'^[0-9][a-z0-9]{3}_[A-Za-z]')


def norm_chain_id(s):
    """PDB ID en minuscules, lettre de chaîne préservée. Ex: '6VEC_A' → '6vec_A'."""
    parts = str(s).split("_", 1)
    return parts[0].lower() + "_" + parts[1] if len(parts) == 2 else str(s).lower()


@st.cache_data
def _cluster_protein_index(aln_dir: str) -> dict:
    """Pour chaque cluster : liste de {name, ids} avec fallback sur le nom le plus commun."""
    result = {}
    for f in sorted(_Path(aln_dir).glob("*.aln")):
        try:
            aln = AlignIO.read(str(f), "fasta")
            entries = []
            for r in aln:
                subunits = r.id.split("|")
                desc = r.description
                if desc.startswith(r.id):
                    desc = desc[len(r.id):].strip()
                names = [p.strip() for p in desc.split(" | ") if p.strip()]
                entries.append({"ids": subunits, "names": names})

            # Noms "réels" (pas des IDs de chaîne PDB) présents dans l'alignement
            real_names = [
                n for e in entries for n in e["names"]
                if not _CHAIN_ID_RE.match(n)
            ]
            fallback = _Counter(real_names).most_common(1)[
                0][0] if real_names else ""

            result[f.stem] = [
                {
                    "name": (e["names"][0] if e["names"] else "") or fallback or e["ids"][0],
                    "ids": e["ids"],
                }
                for e in entries
            ]
        except Exception:
            result[f.stem] = []
    return result


_MSA_ALN_DIR    = _Path("data/alignments")
_MSA_RIGOR_PDBS = {"5jlh", "6c1h", "7aln", "7jh7", "8efi", "8enc", "8zb7", "8zbm", "8zbn"}

_MSA_TAX = {
    9606: "H.sapiens", 10090: "M.musculus", 9986: "O.cuniculus",
    10116: "R.norvegicus", 9823: "S.scrofa", 9031: "G.gallus",
    9913: "B.taurus", 44689: "D.discoideum", 5811: "T.gondii",
    36329: "P.falciparum", 5843: "P.knowlesi", 5833: "P.vivax",
    5791: "Plasmodium", 5850: "P.falciparum3D7", 7227: "D.melanogaster",
    4932: "S.cerevisiae", 559292: "S.cerevisiae", 4896: "S.pombe",
    6577: "A.irradians", 31199: "A.irradians", 6637: "Mollusc",
    6100: "other", 229533: "Acanthamoeba", 5754: "Acanthamoeba",
    5759: "Acanthamoeba", 1051067: "bacteria", 1579: "bacteria",
    32630: "synthetic", 32644: "unknown", 0: "chimera",
}


def _msa_extract_interface_seqs(filter_fn, rigor_pdbs=None):
    """
    Pour chaque séquence unique du groupe, extrait les résidus d'interface
    (union des positions sur toutes les chaînes de cette séquence) et retourne
    une DataFrame : seq_id, title, organism, length_full, n_interface, interface_seq.
    Clé de coloration : interface_seq.lower() → set(range(1, n+1)) (toutes positions = interface).
    """
    from collections import defaultdict
    filt_path = _Path("data/filtered/filtered_all_data.csv")
    int3_path = _Path("data/filtered/details/3.interface_residues.csv")
    if not filt_path.exists() or not int3_path.exists():
        return pd.DataFrame()

    df_filt = pd.read_csv(filt_path, low_memory=False)
    if rigor_pdbs:
        df_filt = df_filt[df_filt["pdb_id"].isin(rigor_pdbs)]
    df3 = pd.read_csv(int3_path)

    mask = df_filt["subunit_2_title"].apply(lambda t: filter_fn(str(t)))
    rows = df_filt[mask][
        ["subunit_2", "subunit_2_title", "s2_sequence", "s2_taxonomy_id"]
    ].drop_duplicates(subset=["subunit_2"])

    chain_to_fullseq = {r["subunit_2"]: str(r["s2_sequence"]).strip() for _, r in rows.iterrows()}
    chain_to_seqlow  = {ch: s.lower() for ch, s in chain_to_fullseq.items()}
    chain_to_title   = {r["subunit_2"]: str(r["subunit_2_title"]) for _, r in rows.iterrows()}
    chain_to_taxid   = {r["subunit_2"]: r["s2_taxonomy_id"] for _, r in rows.iterrows()}

    res = df3[df3["chain"].isin(chain_to_fullseq.keys())][
        ["chain", "residue_number_sequence"]
    ].dropna()

    # seqlow → {chain → set(pos)}
    seq_chains_pos: dict = defaultdict(dict)
    for _, row in res.iterrows():
        ch = row["chain"]
        seqlow = chain_to_seqlow.get(ch)
        if seqlow is None:
            continue
        seq_chains_pos[seqlow].setdefault(ch, set()).add(int(row["residue_number_sequence"]))

    # Premier chain représentatif par séquence unique
    seqlow_to_chain: dict = {}
    for ch, seqlow in chain_to_seqlow.items():
        seqlow_to_chain.setdefault(seqlow, ch)

    records = []
    for seqlow, chains_pos in seq_chains_pos.items():
        union_pos = set()
        for s in chains_pos.values():
            union_pos |= s
        if not union_pos:
            continue
        first_ch = seqlow_to_chain[seqlow]
        full_seq  = chain_to_fullseq[first_ch]
        iface_seq = "".join(
            full_seq[p - 1] if 0 < p <= len(full_seq) else "X"
            for p in sorted(union_pos)
        )
        taxid    = chain_to_taxid[first_ch]
        organism = _MSA_TAX.get(int(taxid), f"taxid:{int(taxid)}") if pd.notna(taxid) else "unknown"
        records.append({
            "title":        chain_to_title[first_ch],
            "organism":     organism,
            "length_full":  len(full_seq),
            "n_interface":  len(union_pos),
            "interface_seq": iface_seq,
        })

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)
    base_ids = (
        df["title"].str[:24].str.replace(r"[^\w]", "_", regex=True)
        + "_" + df["organism"]
        + "_n" + df["n_interface"].astype(str)
    )
    seen: dict = {}
    uid: list = []
    for bid in base_ids:
        cnt = seen.get(bid, 0)
        uid.append(bid if cnt == 0 else f"{bid}_{cnt}")
        seen[bid] = cnt + 1
    df["seq_id"] = uid
    return df[["seq_id", "title", "organism", "length_full", "n_interface", "interface_seq"]].reset_index(drop=True)


def _msa_run_mafft(fasta_path: _Path, aln_path: _Path):
    mafft_bin = shutil.which("mafft")
    if mafft_bin is None:
        return False, "mafft introuvable dans le PATH"
    try:
        with open(aln_path, "w") as out:
            res = subprocess.run(
                [mafft_bin, "--auto", "--reorder", str(fasta_path)],
                stdout=out, stderr=subprocess.PIPE, text=True, timeout=600,
            )
        if res.returncode != 0:
            return False, res.stderr[:500]
        return aln_path.exists() and aln_path.stat().st_size > 0, res.stderr
    except subprocess.TimeoutExpired:
        return False, "Timeout (>10 min)"


def _msa_col_colors(alignment, core_by_seqlow: dict, var_by_seqlow: dict) -> dict:
    """
    Précalcule pour chaque (record.id, colonne_alignement) la couleur :
      'orange' : majorité (≥50 %) des séquences en interaction à cette colonne ET aa conservé.
      'red'    : majorité (≥50 %) des séquences en interaction ET aa variable.
      'yellow' : minorité (<50 %) des séquences en interaction à cette colonne.
    Retourne dict{ record.id → dict{ aln_col → 'orange'|'red'|'yellow' } }
    """
    from collections import Counter

    n_seqs = len(alignment)

    # Étape 1 : seq_pos → colonne alignement pour chaque séquence
    seq_to_aln: dict = {}
    for record in alignment:
        m: dict = {}
        sp = 0
        for col, aa in enumerate(str(record.seq)):
            if aa != "-":
                sp += 1
                m[sp] = col
        seq_to_aln[str(record.id)] = m

    # Étape 2 : nb séquences non-gap par colonne
    non_gap_count: dict = {}
    for record in alignment:
        for col, aa in enumerate(str(record.seq)):
            if aa != "-":
                non_gap_count[col] = non_gap_count.get(col, 0) + 1

    # Étape 3 : pour chaque colonne, compter les séquences avec résidu d'interface + Counter AA
    col_counter: dict = {}        # aln_col → Counter{aa → count}
    col_interface_count: dict = {}  # aln_col → nb séquences avec interface ici
    for record in alignment:
        seq_low  = str(record.seq).replace("-", "").lower()
        core_pos = core_by_seqlow.get(seq_low, set())
        var_pos  = var_by_seqlow.get(seq_low,  set())
        seq_str  = str(record.seq)
        for sp, col in seq_to_aln[str(record.id)].items():
            if sp in core_pos or sp in var_pos:
                aa = seq_str[col].upper()
                if col not in col_counter:
                    col_counter[col] = Counter()
                col_counter[col][aa] += 1
                col_interface_count[col] = col_interface_count.get(col, 0) + 1

    # Étape 4 : couleur par (record, colonne)
    result: dict = {}
    for record in alignment:
        seq_low  = str(record.seq).replace("-", "").lower()
        core_pos = core_by_seqlow.get(seq_low, set())
        var_pos  = var_by_seqlow.get(seq_low,  set())
        seq_str  = str(record.seq)
        colors: dict = {}
        for sp, col in seq_to_aln[str(record.id)].items():
            if sp in core_pos or sp in var_pos:
                n_iface  = col_interface_count.get(col, 0)
                n_nongap = non_gap_count.get(col, n_seqs)
                if n_nongap > 0 and n_iface / n_nongap >= 0.5:
                    # Majorité en interaction → vert (conservé) ou violet (variable)
                    counter   = col_counter.get(col, Counter())
                    max_count = max(counter.values()) if counter else 1
                    aa        = seq_str[col].upper()
                    colors[col] = "orange" if counter[aa] == max_count else "red"
                else:
                    # Minorité en interaction → rouge
                    colors[col] = "yellow"
        result[str(record.id)] = colors
    return result


def _msa_render_full(alignment, core_by_seqlow: dict, var_by_seqlow: dict,
                     cols_per_line: int = 60, label_map: dict = None) -> str:
    """
    HTML de l'alignement complet.
    Orange : cœur d'interface ET aa conservé à cette colonne MSA.
    Rouge  : interface mais aa variable ou position non partagée.
    Gris   : non-interface.
    """
    col_colors = _msa_col_colors(alignment, core_by_seqlow, var_by_seqlow)
    aln_len = alignment.get_alignment_length()
    n_lines = (aln_len + cols_per_line - 1) // cols_per_line
    _lm = label_map or {}
    label_w = min(max(len(_lm.get(str(r.id), str(r.id))) for r in alignment), 36)
    parts = [
        '<div style="font-family:\'Courier New\',Courier,monospace;font-size:11px;'
        'line-height:1.6;background:#fafafa;padding:10px;'
        'overflow-x:auto;'
        'user-select:none;-webkit-user-select:none;cursor:default;">'
    ]
    for li in range(n_lines):
        c0 = li * cols_per_line
        c1 = min(c0 + cols_per_line, aln_len)
        parts.append(
            f'<div style="color:#aaa;font-size:10px;margin:{("14px" if li else "0")} 0 2px 0">'
            f'Positions {c0 + 1}–{c1}</div>'
        )
        parts.append('<table style="border-collapse:collapse;">')
        for record in alignment:
            rec_colors = col_colors.get(str(record.id), {})
            _disp      = _lm.get(str(record.id), str(record.id))
            label      = _disp[:label_w].ljust(label_w).replace(" ", "&nbsp;")
            block      = str(record.seq)[c0:c1]
            seg_parts  = []
            buf_grey   = []
            for i, aa in enumerate(block):
                aln_col = c0 + i
                if aa == "-":
                    if buf_grey:
                        seg_parts.append(f'<span style="color:#ccc">{"".join(buf_grey)}</span>')
                        buf_grey = []
                    seg_parts.append('<span style="color:#e0e0e0">-</span>')
                    continue
                aa_u  = aa.upper()
                color = rec_colors.get(aln_col)
                if color == "orange":
                    if buf_grey:
                        seg_parts.append(f'<span style="color:#ccc">{"".join(buf_grey)}</span>')
                        buf_grey = []
                    seg_parts.append(f'<span style="background:#27ae60;color:#fff">{aa_u}</span>')
                elif color == "red":
                    if buf_grey:
                        seg_parts.append(f'<span style="color:#ccc">{"".join(buf_grey)}</span>')
                        buf_grey = []
                    seg_parts.append(f'<span style="background:#8e44ad;color:#fff">{aa_u}</span>')
                elif color == "yellow":
                    if buf_grey:
                        seg_parts.append(f'<span style="color:#ccc">{"".join(buf_grey)}</span>')
                        buf_grey = []
                    seg_parts.append(f'<span style="background:#e74c3c;color:#fff">{aa_u}</span>')
                else:
                    buf_grey.append(aa_u)
            if buf_grey:
                seg_parts.append(f'<span style="color:#ccc">{"".join(buf_grey)}</span>')
            parts.append(
                f'<tr>'
                f'<td style="padding-right:12px;color:#555;font-size:10px;'
                f'white-space:pre;vertical-align:top">{label}</td>'
                f'<td style="white-space:pre">{"".join(seg_parts)}</td>'
                f'</tr>'
            )
        parts.append("</table>")
    # Légende
    parts.append(
        '<div style="margin-top:14px;font-size:10px;color:#666">'
        '<span style="background:#27ae60;color:#fff;padding:1px 6px;border-radius:2px;margin-right:6px">'
        '&#9632; Majorité en interaction (aa conservé)</span>'
        '<span style="background:#8e44ad;color:#fff;padding:1px 6px;border-radius:2px;margin-right:6px">'
        '&#9632; Majorité en interaction (aa variable)</span>'
        '<span style="background:#e74c3c;color:#fff;padding:1px 6px;border-radius:2px;margin-right:6px">'
        '&#9632; Minorité en interaction</span>'
        '<span style="color:#ccc;margin-right:6px">&#9632; Non-interface</span>'
        '</div></div>'
    )
    return "".join(parts)


# ── Helpers MSA ───────────────────────────────────────────────────────────────

def _msa_load_seqs_g(filter_fn, rigor_pdbs=None):
    """Charge les séquences COMPLÈTES pour un groupe (dédupliquées par séquence lowercase)."""
    filt_path = _Path("data/filtered/filtered_all_data.csv")
    if not filt_path.exists():
        return pd.DataFrame()
    df = pd.read_csv(filt_path, low_memory=False)
    if rigor_pdbs:
        df = df[df["pdb_id"].isin(rigor_pdbs)]
    mask2 = df["subunit_2_title"].apply(lambda t: filter_fn(str(t)))
    s2 = df[mask2][["subunit_2_title", "s2_sequence", "s2_taxonomy_id"]].rename(
        columns={"subunit_2_title": "title", "s2_sequence": "seq", "s2_taxonomy_id": "taxid"}
    )
    combined = s2.dropna(subset=["seq"]).copy()
    combined = combined[combined["seq"].str.strip().str.lower() != "nan"]
    combined["seq"] = combined["seq"].str.strip()
    combined["_key"] = combined["seq"].str.lower()
    combined = combined.drop_duplicates(subset=["_key"]).drop(columns=["_key"]).reset_index(drop=True)
    combined["length"] = combined["seq"].str.len()
    combined["organism"] = combined["taxid"].apply(
        lambda t: _MSA_TAX.get(int(t), f"taxid:{int(t)}") if pd.notna(t) else "unknown"
    )
    base_ids = (
        combined["title"].str[:28].str.replace(r"[^\w]", "_", regex=True)
        + "_" + combined["organism"]
        + "_" + combined["length"].astype(str)
    )
    seen: dict = {}
    uid: list = []
    for bid in base_ids:
        cnt = seen.get(bid, 0)
        uid.append(bid if cnt == 0 else f"{bid}_{cnt}")
        seen[bid] = cnt + 1
    combined["seq_id"] = uid
    return combined[["seq_id", "title", "organism", "length", "seq"]].reset_index(drop=True)


def _msa_load_interface_g(filter_fn, rigor_pdbs=None):
    """Retourne (core, variable) : seqlow → set(int) pour un groupe."""
    from collections import defaultdict
    filt_path = _Path("data/filtered/filtered_all_data.csv")
    int3_path = _Path("data/filtered/details/3.interface_residues.csv")
    if not filt_path.exists() or not int3_path.exists():
        return {}, {}
    df_filt = pd.read_csv(filt_path, low_memory=False)
    if rigor_pdbs:
        df_filt = df_filt[df_filt["pdb_id"].isin(rigor_pdbs)]
    df3 = pd.read_csv(int3_path)
    rows = df_filt[df_filt["subunit_2_title"].apply(lambda t: filter_fn(str(t)))][
        ["subunit_2", "s2_sequence"]
    ].drop_duplicates(subset=["subunit_2"])
    chain_to_seqlow = {r["subunit_2"]: str(r["s2_sequence"]).strip().lower() for _, r in rows.iterrows()}
    res = df3[df3["chain"].isin(chain_to_seqlow)][["chain", "residue_number_sequence"]].dropna()
    per_chain: dict = defaultdict(dict)
    for _, row in res.iterrows():
        seqlow = chain_to_seqlow.get(row["chain"])
        if seqlow:
            per_chain[seqlow].setdefault(row["chain"], set()).add(int(row["residue_number_sequence"]))
    core: dict = {}
    variable: dict = {}
    for seqlow, chains_pos in per_chain.items():
        all_sets = list(chains_pos.values())
        if len(all_sets) == 1:
            core[seqlow], variable[seqlow] = all_sets[0].copy(), set()
        else:
            inter = all_sets[0].copy()
            union = all_sets[0].copy()
            for s in all_sets[1:]:
                inter &= s; union |= s
            core[seqlow], variable[seqlow] = inter, union - inter
    return core, variable


def _msa_actin_contacts_from_pairs(df_pairs, df1, df3):
    """
    Calcule, pour chaque séquence ABP unique (S2), les positions canoniques actine (S1) contactées.

    df_pairs : DataFrame avec colonnes subunit_1, subunit_2, subunit_2_title, s2_sequence, s2_taxonomy_id
    df1      : 1.interactions.csv (déjà chargé)
    df3      : 3.interface_residues.csv (déjà chargé, avec chain_lower et residue_number_canon_mafft)

    Retourne (abp_rows, aa_at_canon, col_color) :
      abp_rows   : list de {seq_id, title, organism, contacts: set(int)}
      aa_at_canon: dict{canon_pos → actin_AA_1_lettre}
      col_color  : dict{canon_pos → 'orange'|'red'|'yellow'}
    """
    from collections import defaultdict, Counter as _Ctr

    if df_pairs.empty:
        return [], {}, {}

    df_joined = df_pairs.merge(
        df1[["interaction_id", "chain_A_id", "chain_B_id"]],
        left_on=["subunit_1", "subunit_2"],
        right_on=["chain_A_id", "chain_B_id"],
        how="left",
    ).dropna(subset=["interaction_id"])
    df_joined["interaction_id"] = df_joined["interaction_id"].astype(int)
    if df_joined.empty:
        return [], {}, {}

    valid_iids = set(df_joined["interaction_id"])
    iid_to_s1 = (
        df_joined.drop_duplicates("interaction_id")
        .set_index("interaction_id")["subunit_1"]
        .str.lower()
    )

    # Résidus S1 (actine) avec position canonique
    df3_rel = df3[
        df3["interaction_id"].isin(valid_iids) &
        df3["residue_number_canon_mafft"].notna()
    ].copy()
    if df3_rel.empty:
        return [], {}, {}
    df3_rel["canon"] = df3_rel["residue_number_canon_mafft"].astype(int)
    df3_rel["expected_s1"] = df3_rel["interaction_id"].map(iid_to_s1)
    df3_s1 = df3_rel[df3_rel["chain_lower"] == df3_rel["expected_s1"]].copy()
    if df3_s1.empty:
        return [], {}, {}

    # iid → set(canon_pos)
    iid_to_canon: dict = df3_s1.groupby("interaction_id")["canon"].apply(set).to_dict()

    # AA actine à chaque position canonique (groupe filtré — pour col_color)
    _VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
    df3_s1["aa1"] = df3_s1["residue_name"].str.strip().str.upper()
    df3_s1_v = df3_s1[df3_s1["aa1"].isin(_VALID_AA)]
    aa_counter: dict = defaultdict(_Ctr)
    for _, row in df3_s1_v.iterrows():
        aa_counter[row["canon"]][row["aa1"]] += 1
    aa_at_canon = {pos: cnt.most_common(1)[0][0] for pos, cnt in aa_counter.items() if cnt}

    # Séquence canonique complète : toutes les interactions S1 (pas juste le groupe filtré)
    _iid_to_s1_all = df1.set_index("interaction_id")["chain_A_id"].str.lower()
    _df3_full = df3[pd.to_numeric(df3["residue_number_canon_mafft"], errors="coerce").notna()].copy()
    _df3_full["canon_f"] = pd.to_numeric(_df3_full["residue_number_canon_mafft"], errors="coerce").astype(int)
    _df3_full["aa_f"]    = _df3_full["residue_name"].str.strip().str.upper()
    if "chain_lower" not in _df3_full.columns:
        _df3_full["chain_lower"] = _df3_full["chain"].str.lower()
    _df3_full["exp_s1_f"] = _df3_full["interaction_id"].map(_iid_to_s1_all)
    _df3_full_s1 = _df3_full[
        (_df3_full["chain_lower"] == _df3_full["exp_s1_f"]) &
        _df3_full["aa_f"].isin(_VALID_AA)
    ]
    aa_at_canon_full: dict = {
        int(pos): grp.value_counts().index[0]
        for pos, grp in _df3_full_s1.groupby("canon_f")["aa_f"]
        if not grp.empty
    }

    # Une ligne par séquence ABP unique (S2)
    s2seq_contacts: dict = defaultdict(set)
    s2seq_meta: dict = {}
    for _, row in df_joined.iterrows():
        iid = int(row["interaction_id"])
        s2seq = str(row.get("s2_sequence", "")).strip().lower()
        if not s2seq or s2seq == "nan":
            continue
        s2seq_contacts[s2seq] |= iid_to_canon.get(iid, set())
        if s2seq not in s2seq_meta:
            taxid = row.get("s2_taxonomy_id")
            org = _MSA_TAX.get(int(taxid), f"taxid:{int(taxid)}") if pd.notna(taxid) else "unknown"
            s2seq_meta[s2seq] = {"title": str(row.get("subunit_2_title", "")), "organism": org}

    abp_rows: list = []
    seen_ids: dict = {}
    for s2seq, meta in s2seq_meta.items():
        title, org = meta["title"], meta["organism"]
        base_id = (title[:24].replace(" ", "_").replace(",", "") + "_" + org)[:40]
        cnt = seen_ids.get(base_id, 0)
        seq_id = base_id if cnt == 0 else f"{base_id}_{cnt}"
        seen_ids[base_id] = cnt + 1
        abp_rows.append({"seq_id": seq_id, "title": title, "organism": org,
                         "contacts": s2seq_contacts[s2seq], "s2seq": s2seq})

    if not abp_rows:
        return [], {}, {}

    # Couleur par position canonique
    n_abp = len(abp_rows)
    col_color: dict = {}
    for pos in aa_at_canon:
        n_contact = sum(1 for r in abp_rows if pos in r["contacts"])
        if n_abp > 0 and n_contact / n_abp >= 0.5:
            cnt = aa_counter.get(pos, _Ctr())
            if cnt:
                max_c = max(cnt.values()); tot = sum(cnt.values())
                col_color[pos] = "orange" if max_c / tot >= 0.8 else "red"
            else:
                col_color[pos] = "orange"
        else:
            col_color[pos] = "yellow"

    return abp_rows, aa_at_canon_full, col_color


def _msa_actin_contacts_per_abp(filter_fn, rigor_pdbs=None):
    """Wrapper : charge les données puis appelle _msa_actin_contacts_from_pairs."""
    filt_path = _Path("data/filtered/filtered_all_data.csv")
    int1_path = _Path("data/filtered/details/1.interactions.csv")
    int3_path = _Path("data/filtered/details/3.interface_residues.csv")
    if not filt_path.exists() or not int1_path.exists() or not int3_path.exists():
        return [], {}, {}
    df_filt = pd.read_csv(filt_path, low_memory=False)
    if rigor_pdbs:
        df_filt = df_filt[df_filt["pdb_id"].isin(rigor_pdbs)]
    mask = df_filt["subunit_2_title"].apply(lambda t: filter_fn(str(t)))
    df_pairs = df_filt[mask][
        ["subunit_1", "subunit_2", "subunit_2_title", "s2_sequence", "s2_taxonomy_id"]
    ].copy()
    df1 = pd.read_csv(int1_path)
    df3 = pd.read_csv(int3_path)
    df3["residue_number_canon_mafft"] = pd.to_numeric(df3["residue_number_canon_mafft"], errors="coerce")
    df3["chain_lower"] = df3["chain"].str.lower()
    return _msa_actin_contacts_from_pairs(df_pairs, df1, df3)


def _msa_render_actin_contacts(abp_rows, aa_at_canon, col_color, cols_per_line=60, label_map=None):
    """
    Affiche les contacts actine-ABP : lignes = séquences ABP, colonnes = positions canoniques actine.
    L'AA actine est affiché à chaque position, coloré si cette séquence ABP la contacte.
    """
    all_positions = sorted(aa_at_canon.keys())
    if not all_positions or not abp_rows:
        return ""
    _lm = label_map or {}
    n_pos = len(all_positions)
    n_lines = (n_pos + cols_per_line - 1) // cols_per_line
    label_w = min(max(len(_lm.get(r["seq_id"], r["seq_id"])) for r in abp_rows), 36)

    parts = [
        '<div style="font-family:\'Courier New\',Courier,monospace;font-size:11px;'
        'line-height:1.6;background:#fafafa;padding:10px;'
        'overflow-x:auto;'
        'user-select:none;-webkit-user-select:none;cursor:default;">'
    ]
    for li in range(n_lines):
        c0 = li * cols_per_line
        c1 = min(c0 + cols_per_line, n_pos)
        pos_slice = all_positions[c0:c1]
        parts.append(
            f'<div style="color:#aaa;font-size:10px;margin:{("14px" if li else "0")} 0 2px 0">'
            f'Positions actine canoniques {pos_slice[0]}–{pos_slice[-1]}</div>'
        )
        parts.append('<table style="border-collapse:collapse;">')
        for row in abp_rows:
            _disp = _lm.get(row["seq_id"], row["seq_id"])
            label = _disp[:label_w].ljust(label_w).replace(" ", "&nbsp;")
            contacts = row["contacts"]
            seg_parts: list = []
            for pos in pos_slice:
                aa = aa_at_canon.get(pos, "-")
                if pos in contacts:
                    color = col_color.get(pos, "yellow")
                    bg = "#27ae60" if color == "orange" else ("#8e44ad" if color == "red" else "#e74c3c")
                    seg_parts.append(f'<span style="background:{bg};color:#fff">{aa}</span>')
                else:
                    seg_parts.append(f'<span style="color:#ccc">{aa}</span>')
            parts.append(
                f'<tr>'
                f'<td style="padding-right:12px;color:#555;font-size:10px;white-space:pre;vertical-align:top">{label}</td>'
                f'<td style="white-space:pre">{"".join(seg_parts)}</td>'
                f'</tr>'
            )
        parts.append("</table>")
    parts.append(
        '<div style="margin-top:14px;font-size:10px;color:#666">'
        '<span style="background:#27ae60;color:#fff;padding:1px 6px;border-radius:2px;margin-right:6px">&#9632; Majorité en interaction (aa actine conservé)</span>'
        '<span style="background:#8e44ad;color:#fff;padding:1px 6px;border-radius:2px;margin-right:6px">&#9632; Majorité en interaction (aa actine variable)</span>'
        '<span style="background:#e74c3c;color:#fff;padding:1px 6px;border-radius:2px;margin-right:6px">&#9632; Minorité en interaction</span>'
        '<span style="color:#ccc;margin-right:6px">AA Non-contacté</span>'
        '</div></div>'
    )
    return "".join(parts)


def _msa_render_projected(full_seqs_by_id, aln_iface, iface_pos_by_id, cols_per_line=100):
    """
    Affiche les séquences COMPLÈTES avec les résidus d'interface colorés
    selon un alignement MAFFT fait sur les sous-séquences d'interface.
    """
    # Couleurs depuis l'alignement interface (toutes positions = interface)
    core_iface: dict = {}
    for record in aln_iface:
        iseq_low = str(record.seq).replace("-", "").lower()
        core_iface[iseq_low] = set(range(1, len(iseq_low) + 1))
    col_colors_iface = _msa_col_colors(aln_iface, core_iface, {})

    # Projection : aln_col → position dans la séquence complète → couleur
    full_seq_colors: dict = {}
    for record in aln_iface:
        sid = str(record.id)
        iface_positions = iface_pos_by_id.get(sid, [])
        rec_colors = col_colors_iface.get(sid, {})
        pos_map: dict = {}
        iface_idx = 0
        for aln_col, char in enumerate(str(record.seq)):
            if char != "-":
                if iface_idx < len(iface_positions):
                    full_pos = iface_positions[iface_idx]
                    if aln_col in rec_colors:
                        pos_map[full_pos] = rec_colors[aln_col]
                    iface_idx += 1
        full_seq_colors[sid] = pos_map

    seq_ids_ord = [str(r.id) for r in aln_iface]
    max_len  = max((len(full_seqs_by_id.get(s, "")) for s in seq_ids_ord), default=0)
    label_w  = min(max((len(s) for s in seq_ids_ord), default=10), 36)
    n_lines  = (max_len + cols_per_line - 1) // cols_per_line

    parts = [
        '<div style="font-family:\'Courier New\',Courier,monospace;font-size:11px;'
        'line-height:1.6;background:#fafafa;padding:10px;'
        'overflow-x:auto;'
        'user-select:none;-webkit-user-select:none;cursor:default;">'
    ]
    for li in range(n_lines):
        c0 = li * cols_per_line
        c1 = min(c0 + cols_per_line, max_len)
        parts.append(
            f'<div style="color:#aaa;font-size:10px;margin:{("14px" if li else "0")} 0 2px 0">'
            f'Positions {c0+1}–{c1}</div>'
        )
        parts.append('<table style="border-collapse:collapse;">')
        for sid in seq_ids_ord:
            full_seq  = full_seqs_by_id.get(sid, "")
            pos_colors = full_seq_colors.get(sid, {})
            label = sid[:label_w].ljust(label_w).replace(" ", "&nbsp;")
            block = full_seq[c0:c1]
            seg_parts: list = []
            buf_grey: list  = []
            for i, aa in enumerate(block):
                full_pos = c0 + i + 1
                aa_u  = aa.upper()
                color = pos_colors.get(full_pos)
                if color in ("orange", "red", "yellow"):
                    if buf_grey:
                        seg_parts.append(f'<span style="color:#ccc">{"".join(buf_grey)}</span>')
                        buf_grey = []
                    bg = "#27ae60" if color == "orange" else ("#8e44ad" if color == "red" else "#e74c3c")
                    seg_parts.append(f'<span style="background:{bg};color:#fff">{aa_u}</span>')
                else:
                    buf_grey.append(aa_u)
            if buf_grey:
                seg_parts.append(f'<span style="color:#ccc">{"".join(buf_grey)}</span>')
            parts.append(
                f'<tr>'
                f'<td style="padding-right:12px;color:#555;font-size:10px;white-space:pre;vertical-align:top">{label}</td>'
                f'<td style="white-space:pre">{"".join(seg_parts)}</td>'
                f'</tr>'
            )
        parts.append("</table>")
    parts.append(
        '<div style="margin-top:14px;font-size:10px;color:#666">'
        '<span style="background:#27ae60;color:#fff;padding:1px 6px;border-radius:2px;margin-right:6px">&#9632; Majorité en interaction (aa conservé)</span>'
        '<span style="background:#8e44ad;color:#fff;padding:1px 6px;border-radius:2px;margin-right:6px">&#9632; Majorité en interaction (aa variable)</span>'
        '<span style="background:#e74c3c;color:#fff;padding:1px 6px;border-radius:2px;margin-right:6px">&#9632; Minorité en interaction</span>'
        '<span style="color:#ccc;margin-right:6px">&#9632; Non-interface</span>'
        '</div></div>'
    )
    return "".join(parts)


def _msa_contact_analysis(filter_fn, group_key, rigor_pdbs=None,
                           _ch2seq=None, _ch2title=None, tabs=None):
    """
    Heatmaps interactives (survol = tooltip) + classification AA des contacts ABP–actine.
      A – Heatmap ABP    : séquences × positions canoniques ABP, couleur = aire de contact
      B – Heatmap Actine : séquences × positions canoniques actine, couleur = aire de contact
      C – Classification : répartition des AA de contact par famille physicochimique

    _ch2seq / _ch2title : mappings chain_id → seq_low / titre, fournis directement
    (évite de relire filtered_all_data.csv pour les clusters S1).
    tabs : liste de clés à afficher, ex. ["A"] ou ["B","C","D","E"]. None = tous.
    """
    import matplotlib.pyplot as plt

    int4_path = _Path("data/filtered/details/4.inter-residue_contacts.csv")
    if not int4_path.exists():
        st.warning("Fichier 4.inter-residue_contacts.csv manquant.")
        return

    if _ch2seq is not None and _ch2title is not None:
        ch2seq   = _ch2seq
        ch2title = _ch2title
    else:
        filt_path = _Path("data/filtered/filtered_all_data.csv")
        if not filt_path.exists():
            st.warning("filtered_all_data.csv manquant.")
            return
        df_filt = pd.read_csv(filt_path, low_memory=False)
        if rigor_pdbs:
            df_filt = df_filt[df_filt["pdb_id"].isin(rigor_pdbs)]
        mask = df_filt["subunit_2_title"].apply(lambda t: filter_fn(str(t)))
        df_grp = df_filt[mask][["subunit_2", "subunit_2_title", "s2_sequence"]].drop_duplicates("subunit_2")
        if df_grp.empty:
            st.info("Aucune séquence correspondante.")
            return
        ch2seq   = {r["subunit_2"]: str(r["s2_sequence"]).strip().lower() for _, r in df_grp.iterrows()}
        ch2title = {r["subunit_2"]: str(r["subunit_2_title"]) for _, r in df_grp.iterrows()}

    df4 = pd.read_csv(int4_path)
    df4 = df4[df4["chain_B_id"].isin(ch2seq)].copy()
    if df4.empty:
        st.info("Aucun contact trouvé pour ce groupe.")
        return

    df4["area_f"]  = pd.to_numeric(
        df4["contact_area"].astype(str).str.replace("<", "", regex=False),
        errors="coerce",
    ).fillna(0.05)
    df4["canon_b"] = pd.to_numeric(df4["residue_B_canon_mafft"], errors="coerce")
    df4["canon_a"] = pd.to_numeric(df4["residue_A_canon_mafft"], errors="coerce")
    df4["aa_b"]    = df4["residue_B_name"].str.strip().str.upper()   # résidu ABP
    df4["aa_a"]    = df4["residue_A_name"].str.strip().str.upper()   # résidu actine
    df4["seq_low"] = df4["chain_B_id"].map(ch2seq)
    df4["title"]   = df4["chain_B_id"].map(ch2title)
    df4 = df4.dropna(subset=["canon_b", "canon_a", "seq_low"])
    df4["canon_b"] = df4["canon_b"].astype(int)
    df4["canon_a"] = df4["canon_a"].astype(int)

    # Labels uniques (dédupliqués si plusieurs chaînes avec le même titre)
    seen_cnt: dict = {}
    title_to_label: dict = {}
    for sl in df4["seq_low"].unique():
        t = df4[df4["seq_low"] == sl]["title"].iloc[0]
        cnt = seen_cnt.get(t, 0)
        title_to_label[sl] = t if cnt == 0 else f"{t} ({cnt})"
        seen_cnt[t] = cnt + 1
    df4["label"] = df4["seq_low"].map(title_to_label)

    # Ordre des lignes : identique à la MSA S2 (fichier .aln MAFFT, --reorder)
    _aln_path_hm = _MSA_ALN_DIR / f"{group_key}_msa.aln"
    _label_rank: dict = {}
    if _aln_path_hm.exists():
        try:
            _aln_hm = AlignIO.read(str(_aln_path_hm), "fasta")
            _seqlow_rank = {
                str(rec.seq).replace("-", "").lower(): i
                for i, rec in enumerate(_aln_hm)
            }
            _label_rank = {
                lbl: _seqlow_rank.get(sl, 9999)
                for sl, lbl in title_to_label.items()
            }
        except Exception:
            pass
    if not _label_rank:
        _label_rank = {lbl: i for i, lbl in enumerate(sorted(title_to_label.values()))}

    # ── Agrégations corrigées ─────────────────────────────────────────────────
    # Moyenne corrigée : somme par structure PDB, puis division par n_structures TOTAL
    # (inclut les structures sans contact à cette position = zéros implicites)
    n_struct_map: dict = df4.groupby("seq_low")["chain_B_id"].nunique().to_dict()

    # ABP : somme par (seq, chain, canon_b) → somme sur toutes les chaînes / n_struct
    _ps_b = df4.groupby(["seq_low", "label", "chain_B_id", "canon_b"])["area_f"].sum().reset_index()
    corr_b = _ps_b.groupby(["seq_low", "label", "canon_b"])["area_f"].sum().reset_index()
    corr_b["corr_mean"] = corr_b["area_f"] / corr_b["seq_low"].map(n_struct_map)
    _tot_b = corr_b.groupby("seq_low")["corr_mean"].sum().to_dict()
    corr_b["pct"] = corr_b.apply(
        lambda r: r["corr_mean"] / _tot_b.get(r["seq_low"], 1) * 100, axis=1
    )

    # Actine : idem pour canon_a
    _ps_a = df4.groupby(["seq_low", "label", "chain_B_id", "canon_a"])["area_f"].sum().reset_index()
    corr_a = _ps_a.groupby(["seq_low", "label", "canon_a"])["area_f"].sum().reset_index()
    corr_a["corr_mean"] = corr_a["area_f"] / corr_a["seq_low"].map(n_struct_map)
    _tot_a = corr_a.groupby("seq_low")["corr_mean"].sum().to_dict()
    corr_a["pct"] = corr_a.apply(
        lambda r: r["corr_mean"] / _tot_a.get(r["seq_low"], 1) * 100, axis=1
    )

    # AA dominant à chaque (seq, position) — pour les tooltips et la classification
    agg_b_aa = (
        df4.groupby(["seq_low", "label", "canon_b", "aa_b"])["area_f"]
        .sum().reset_index()
        .sort_values("area_f", ascending=False)
        .drop_duplicates(["seq_low", "canon_b"])
    )
    agg_a_aa = (
        df4.groupby(["seq_low", "label", "canon_a", "aa_a"])["area_f"]
        .sum().reset_index()
        .sort_values("area_f", ascending=False)
        .drop_duplicates(["seq_low", "canon_a"])
    )
    # garder aussi agg_b / agg_a pour le tab C classification (besoin du n_contacts)
    agg_b = agg_b_aa.rename(columns={"area_f": "mean_area"}).copy()
    agg_a = agg_a_aa.rename(columns={"area_f": "mean_area"}).copy()

    # AA actine majoritaire à chaque position canonique ABP (tous ABPs confondus)
    aa_maj_actin_for_abp: dict = (
        df4.groupby("canon_b")["aa_a"]
        .agg(lambda x: x.value_counts().index[0] if len(x) > 0 else "?")
        .to_dict()
    )
    # Pour chaque (seq_low, canon_a) → AA ABP spécifique le plus contacté
    aa_b_specific: dict = (
        df4.groupby(["seq_low", "canon_a"])["aa_b"]
        .agg(lambda x: x.value_counts().index[0] if len(x) > 0 else "?")
        .to_dict()
    )

    # ── Générateur HTML heatmap interactive ──────────────────────────────────
    def _html_interactive_heatmap(pivot_df, row_labels, col_labels, max_area,
                                   tooltip_fn, title_str, subtitle_str, cell_px=9):
        def _bg(v):
            if pd.isna(v) or v == 0:
                return "background:#f0f0f0", "#bbb"
            t = min(1.0, v / max_area)
            r  = int(255 - (255 - 139) * t)
            gc = int(255 * (1 - t))
            bc = int(255 * (1 - t))
            fg = "#fff" if t > 0.45 else "#222"
            return f"background:rgb({max(0,r)},{max(0,gc)},{max(0,bc)})", fg

        col_step = max(1, len(col_labels) // 25)
        LABEL_W  = 200

        js = """
<div id="htt" style="position:fixed;background:rgba(20,20,20,0.88);color:#fff;
padding:7px 11px;border-radius:7px;font-size:11px;line-height:1.6;
pointer-events:none;display:none;z-index:99999;max-width:300px;
font-family:'Helvetica Neue',Arial,sans-serif;white-space:pre-wrap"></div>
<script>
(function(){
var t=document.getElementById('htt');
function show(el,txt){t.innerHTML=txt;t.style.display='block';}
function move(e){t.style.left=(e.clientX+16)+'px';t.style.top=(e.clientY+16)+'px';}
function hide(){t.style.display='none';}
document.querySelectorAll('[data-tt]').forEach(function(el){
  el.addEventListener('mouseenter',function(e){show(el,el.getAttribute('data-tt'));move(e);});
  el.addEventListener('mousemove',move);
  el.addEventListener('mouseleave',hide);
});
})();
</script>"""

        grad_cells = "".join(
            f'<span style="background:rgb({max(0,int(255-(255-139)*(i/11)))},{max(0,int(255*((11-i)/11)))},{max(0,int(255*((11-i)/11)))});'
            f'width:18px;height:12px;display:inline-block;border:1px solid #ddd"></span>'
            for i in range(12)
        )

        parts = [
            '<div style="font-family:\'Helvetica Neue\',Arial,sans-serif;font-size:11px;'
            'background:#fff;padding:12px;overflow-x:auto;">',
            f'<div style="font-weight:600;font-size:12px;margin-bottom:3px">{title_str}</div>',
            f'<div style="font-size:10px;color:#888;margin-bottom:8px">{subtitle_str}</div>',
            '<div style="display:flex;align-items:center;gap:3px;margin-bottom:10px;font-size:10px;color:#666">',
            '<span style="margin-right:3px">0 Å²</span>',
            grad_cells,
            f'<span style="margin-left:3px">{max_area:.1f} Å²</span></div>',
            js,
            '<table style="border-collapse:collapse;table-layout:fixed"><thead><tr>',
            f'<th style="width:{LABEL_W}px;min-width:{LABEL_W}px"></th>',
        ]
        for j, c in enumerate(col_labels):
            lbl = str(c) if j % col_step == 0 else ""
            parts.append(
                f'<th style="width:{cell_px}px;min-width:{cell_px}px;font-size:7px;font-weight:normal;'
                f'color:#aaa;writing-mode:vertical-rl;text-align:left;padding:0;'
                f'height:26px;vertical-align:bottom">{lbl}</th>'
            )
        parts.append("</tr></thead><tbody>")

        data_arr = pivot_df.values
        for i, rl in enumerate(row_labels):
            short = (rl[:30] + "…") if len(rl) > 31 else rl
            parts.append("<tr>")
            parts.append(
                f'<td style="width:{LABEL_W}px;max-width:{LABEL_W}px;overflow:hidden;'
                f'white-space:nowrap;text-overflow:ellipsis;font-size:9px;color:#444;'
                f'padding:1px 6px 1px 0;text-align:right" title="{rl}">{short}</td>'
            )
            for j, v in enumerate(data_arr[i]):
                bg_style, _fg = _bg(v)
                tt = tooltip_fn(i, j, v)
                tt_esc = tt.replace("&", "&amp;").replace('"', "&quot;")
                parts.append(
                    f'<td style="width:{cell_px}px;min-width:{cell_px}px;height:14px;'
                    f'{bg_style};cursor:default;padding:0" data-tt="{tt_esc}"></td>'
                )
            parts.append("</tr>")
        parts.append("</tbody></table></div>")
        return "".join(parts)

    _ALL_TAB_KEYS  = ["A", "B", "C", "D", "E"]
    _ALL_TAB_NAMES = [
        "A — Heatmap ABP",
        "B — Heatmap Actine",
        "C — Classification AA",
        "D — Paires AA",
        "E — Spécificité",
    ]
    _show_tabs  = set(tabs) if tabs else set(_ALL_TAB_KEYS)
    _tab_labels = [n for k, n in zip(_ALL_TAB_KEYS, _ALL_TAB_NAMES) if k in _show_tabs]
    _tab_objs   = st.tabs(_tab_labels)
    _tab_map: dict = {}
    _idx = 0
    for _k in _ALL_TAB_KEYS:
        if _k in _show_tabs:
            _tab_map[_k] = _tab_objs[_idx]
            _idx += 1

    # ── TAB A : Heatmap ABP (positions canoniques ABP) ────────────────────────
    if "A" in _tab_map:
        with _tab_map["A"]:
            st.caption(
                "Lignes = séquences ABP · Colonnes = positions canoniques ABP (MAFFT) · "
                "Couleur = **% de l'interface totale** de cette séquence · "
                "Moyenne corrigée (zéros inclus pour les structures sans contact) · "
                "**Survolez** une case pour les détails."
            )
            _corr_b_aa = corr_b.merge(
                agg_b_aa[["seq_low", "canon_b", "aa_b"]], on=["seq_low", "canon_b"], how="left"
            )
            pivot_b = _corr_b_aa.pivot_table(
                index="label", columns="canon_b", values="pct", aggfunc="mean"
            )
            pivot_b = pivot_b.loc[sorted(pivot_b.index, key=lambda lb: _label_rank.get(lb, 9999))]
            row_labels_b = list(pivot_b.index)
            col_labels_b = list(pivot_b.columns)
            max_pct_b    = float(corr_b["pct"].max()) if not corr_b.empty else 1.0
            aa_b_lut     = {(r["label"], int(r["canon_b"])): r["aa_b"]
                            for _, r in _corr_b_aa.iterrows()}
            cm_b_lut     = {(r["label"], int(r["canon_b"])): r["corr_mean"]
                            for _, r in _corr_b_aa.iterrows()}

            def _tt_b(i, j, v):
                rl  = row_labels_b[i]
                cb  = col_labels_b[j]
                aa  = aa_b_lut.get((rl, cb), "?")
                aaa = aa_maj_actin_for_abp.get(cb, "?")
                if pd.isna(v) or v == 0:
                    return f"Pos ABP: {cb}\nAucun contact"
                cm = cm_b_lut.get((rl, cb), 0.0)
                return (
                    f"Pos ABP : {cb}\n"
                    f"Résidu ABP ({rl[:22]}) : {aa}\n"
                    f"Résidu actine fréq.  : {aaa}\n"
                    f"% interface          : {v:.2f}%\n"
                    f"Aire moy. corrigée   : {cm:.1f} Å²"
                )

            _h_a = _html_interactive_heatmap(
                pivot_b, row_labels_b, col_labels_b, max_pct_b, _tt_b,
                "Contacts ABP–actine — côté ABP",
                f"{len(row_labels_b)} séquences · {len(col_labels_b)} positions · "
                f"valeur = % de l'interface totale par séquence · max {max_pct_b:.2f}%",
            )
            st.components.v1.html(_h_a, height=max(len(row_labels_b) * 16 + 180, 300), scrolling=True)

    # ── TAB B : Heatmap Actine (positions canoniques actine) ──────────────────
    if "B" in _tab_map:
        with _tab_map["B"]:
            st.caption(
                "Lignes = séquences ABP · Colonnes = positions canoniques actine (MAFFT) · "
                "Couleur = **% de l'interface totale** de cette séquence · "
                "Moyenne corrigée (zéros inclus) · **Survolez** une case pour les détails."
            )
            _corr_a_aa = corr_a.merge(
                agg_a_aa[["seq_low", "canon_a", "aa_a"]], on=["seq_low", "canon_a"], how="left"
            )
            pivot_a = _corr_a_aa.pivot_table(
                index="label", columns="canon_a", values="pct", aggfunc="mean"
            )
            pivot_a = pivot_a.loc[sorted(pivot_a.index, key=lambda lb: _label_rank.get(lb, 9999))]
            row_labels_a = list(pivot_a.index)
            col_labels_a = list(pivot_a.columns)
            max_pct_a    = float(corr_a["pct"].max()) if not corr_a.empty else 1.0
            aa_a_lut     = {(r["label"], int(r["canon_a"])): r["aa_a"]
                            for _, r in _corr_a_aa.iterrows()}
            cm_a_lut     = {(r["label"], int(r["canon_a"])): r["corr_mean"]
                            for _, r in _corr_a_aa.iterrows()}

            def _tt_a(i, j, v):
                rl        = row_labels_a[i]
                ca        = col_labels_a[j]
                aa_a_val  = aa_a_lut.get((rl, ca), "?")
                sl        = next((s for s, lb in title_to_label.items() if lb == rl), None)
                aa_b_spec = aa_b_specific.get((sl, ca), "?") if sl else "?"
                if pd.isna(v) or v == 0:
                    return f"Pos actine: {ca}\nAucun contact"
                cm = cm_a_lut.get((rl, ca), 0.0)
                return (
                    f"Pos actine : {ca}\n"
                    f"Résidu actine        : {aa_a_val}\n"
                    f"Résidu ABP ({rl[:18]}) : {aa_b_spec}\n"
                    f"% interface          : {v:.2f}%\n"
                    f"Aire moy. corrigée   : {cm:.1f} Å²"
                )

            _h_b = _html_interactive_heatmap(
                pivot_a, row_labels_a, col_labels_a, max_pct_a, _tt_a,
                "Contacts ABP–actine — côté actine",
                f"{len(row_labels_a)} séquences · {len(col_labels_a)} positions · "
                f"valeur = % de l'interface totale par séquence · max {max_pct_a:.2f}%",
            )
            st.components.v1.html(_h_b, height=max(len(row_labels_a) * 16 + 180, 300), scrolling=True)

    # ── TAB C : Classification AA ─────────────────────────────────────────────
    if "C" in _tab_map:
        with _tab_map["C"]:
            _CLASSES = {
                "Hydrophobe":      set("AVILMFWP"),
                "Polaire (neutre)": set("STCYNQ"),
                "Chargé (+)":      set("KRH"),
                "Chargé (−)":      set("DE"),
            }
            _CLASS_COL = {
                "Hydrophobe":       "#66bb6a",
                "Polaire (neutre)": "#42a5f5",
                "Chargé (+)":       "#ef5350",
                "Chargé (−)":       "#ab47bc",
            }
            def _cls(aa):
                for c, s in _CLASSES.items():
                    if aa in s:
                        return c
                return "Autre"

            agg_b["class"] = agg_b["aa_b"].apply(_cls)
            agg_a["class"] = agg_a["aa_a"].apply(_cls)

            col_left, col_right = st.columns(2)

            for side, corr_s, aa_s, aa_col, lbl in [
                ("ABP",    corr_b, agg_b_aa, "aa_b", col_left),
                ("Actine", corr_a, agg_a_aa, "aa_a", col_right),
            ]:
                cls_data = corr_s.merge(aa_s[["seq_low", aa_col.replace("aa_","canon_"), aa_col]],
                                        on=["seq_low", aa_col.replace("aa_","canon_")], how="left")
                cls_data["class"] = cls_data[aa_col].apply(_cls)
                tot_corr = cls_data["corr_mean"].sum() or 1.0
                cls_sum = (
                    cls_data.groupby("class")["corr_mean"]
                    .agg(aire_totale="sum", n_contacts="count")
                    .reset_index()
                )
                cls_sum["pct"] = cls_sum["aire_totale"] / tot_corr * 100
                cls_sum = cls_sum.sort_values("pct", ascending=False)
                aa_sum = (
                    cls_data.groupby([aa_col, "class"])["corr_mean"]
                    .sum().reset_index()
                )
                aa_sum["pct"] = aa_sum["corr_mean"] / tot_corr * 100
                aa_sum = aa_sum.sort_values("pct", ascending=True)
                with lbl:
                    st.markdown(f"**Résidus {side} en contact**")
                    fig_c, axes = plt.subplots(1, 2, figsize=(6, 3))
                    colors_c = [_CLASS_COL.get(c, "#aaa") for c in cls_sum["class"]]
                    axes[0].barh(cls_sum["class"], cls_sum["pct"], color=colors_c)
                    axes[0].set_xlabel("% de l'interface totale", fontsize=8)
                    axes[0].tick_params(labelsize=8)
                    axes[0].set_title("Par famille", fontsize=9)
                    colors_aa = [_CLASS_COL.get(c, "#aaa") for c in aa_sum["class"]]
                    axes[1].barh(aa_sum[aa_col], aa_sum["pct"], color=colors_aa)
                    axes[1].set_xlabel("% de l'interface totale", fontsize=8)
                    axes[1].tick_params(labelsize=8)
                    axes[1].set_title("Par AA", fontsize=9)
                    from matplotlib.patches import Patch
                    legend_els = [Patch(facecolor=c, label=n) for n, c in _CLASS_COL.items()]
                    axes[1].legend(handles=legend_els, fontsize=6, loc="lower right")
                    plt.tight_layout()
                    st.pyplot(fig_c, use_container_width=True)
                    plt.close(fig_c)
                    st.dataframe(
                        cls_sum[["class","pct","n_contacts"]].rename(
                            columns={"class":"Famille","pct":"% interface","n_contacts":"Contacts"}
                        ).assign(**{"% interface": lambda d: d["% interface"].round(1)}),
                        use_container_width=True, hide_index=True,
                    )

    # ── TAB D : Paires AA spécifiques ABP ↔ Actine ───────────────────────────
    if "D" in _tab_map:
        with _tab_map["D"]:
            st.caption(
                "**Matrice AA-AA** : quels acides aminés ABP interagissent avec quels AA actine. "
                "**Table par séquence** : choisir un ABP pour voir tous ses contacts spécifiques."
            )

            _AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")
            pair_agg = (
                df4.groupby(["aa_b", "aa_a"])
                .agg(aire_tot=("area_f", "sum"), n=("area_f", "count"))
                .reset_index()
            )
            pair_agg = pair_agg[
                pair_agg["aa_b"].isin(_AA_ORDER) & pair_agg["aa_a"].isin(_AA_ORDER)
            ]

            pivot_aa_raw = pair_agg.pivot_table(
                index="aa_a", columns="aa_b", values="aire_tot", aggfunc="sum", fill_value=0
            ).reindex(index=_AA_ORDER, columns=_AA_ORDER, fill_value=0)

            _total_aa = float(pivot_aa_raw.values.sum()) or 1.0
            pivot_aa  = pivot_aa_raw / _total_aa * 100

            max_aa = float(pivot_aa.values.max()) if pivot_aa.values.max() > 0 else 1.0

            pivot_n = pair_agg.pivot_table(
                index="aa_a", columns="aa_b", values="n", aggfunc="sum", fill_value=0
            ).reindex(index=_AA_ORDER, columns=_AA_ORDER, fill_value=0)

            def _bg_aa(v):
                if v == 0:
                    return "background:#f5f5f5"
                t = min(1.0, v / max_aa)
                r  = int(255 - (255 - 139) * t)
                gc = int(255 * (1 - t))
                bc = int(255 * (1 - t))
                return f"background:rgb({max(0,r)},{max(0,gc)},{max(0,bc)})"

            grad_aa = "".join(
                f'<span style="background:rgb({max(0,int(255-(255-139)*(i/11)))},{max(0,int(255*((11-i)/11)))},{max(0,int(255*((11-i)/11)))});'
                f'width:16px;height:11px;display:inline-block;border:1px solid #eee"></span>'
                for i in range(12)
            )
            js_aa = """
<div id="htt2" style="position:fixed;background:rgba(20,20,20,0.88);color:#fff;
padding:7px 11px;border-radius:7px;font-size:11px;line-height:1.6;
pointer-events:none;display:none;z-index:99999;max-width:260px;
font-family:'Helvetica Neue',Arial,sans-serif;white-space:pre-wrap"></div>
<script>
(function(){
var t=document.getElementById('htt2');
document.querySelectorAll('[data-tt2]').forEach(function(el){
  el.addEventListener('mouseenter',function(e){
    t.innerHTML=el.getAttribute('data-tt2');t.style.display='block';
    t.style.left=(e.clientX+16)+'px';t.style.top=(e.clientY+16)+'px';
  });
  el.addEventListener('mousemove',function(e){
    t.style.left=(e.clientX+16)+'px';t.style.top=(e.clientY+16)+'px';
  });
  el.addEventListener('mouseleave',function(){t.style.display='none';});
});
})();
</script>"""

            CELL = 22
            parts_aa = [
                '<div style="font-family:\'Helvetica Neue\',Arial,sans-serif;font-size:11px;'
                'background:#fff;padding:12px;overflow-x:auto;">',
                '<div style="font-weight:600;font-size:12px;margin-bottom:3px">'
                'Matrice AA ABP × AA Actine — % de l\'aire de contact totale</div>',
                '<div style="font-size:10px;color:#888;margin-bottom:6px">'
                'Lignes = AA actine · Colonnes = AA ABP · Valeur = % de l\'interface totale · Survol = détail</div>',
                '<div style="display:flex;align-items:center;gap:3px;margin-bottom:10px;font-size:10px;color:#666">',
                '<span style="margin-right:3px">0%</span>', grad_aa,
                f'<span style="margin-left:3px">{max_aa:.1f}%</span></div>',
                js_aa,
                '<table style="border-collapse:collapse"><thead><tr>',
                f'<th style="width:28px;font-size:9px;color:#999;text-align:right;padding-right:4px">Actine↓ ABP→</th>',
            ]
            for ab in _AA_ORDER:
                parts_aa.append(
                    f'<th style="width:{CELL}px;font-size:9px;color:#555;text-align:center;padding:1px">{ab}</th>'
                )
            parts_aa.append('</tr></thead><tbody>')

            for aa_row in _AA_ORDER:
                parts_aa.append(f'<tr><td style="font-size:9px;color:#555;text-align:right;padding-right:4px;font-weight:600">{aa_row}</td>')
                for aa_col_aa in _AA_ORDER:
                    v = pivot_aa.at[aa_row, aa_col_aa] if (aa_row in pivot_aa.index and aa_col_aa in pivot_aa.columns) else 0
                    n = int(pivot_n.at[aa_row, aa_col_aa]) if (aa_row in pivot_n.index and aa_col_aa in pivot_n.columns) else 0
                    bg = _bg_aa(v)
                    fg = "#fff" if v / max_aa > 0.45 else "#333"
                    if v > 0:
                        raw_area = pivot_aa_raw.at[aa_row, aa_col_aa] if (aa_row in pivot_aa_raw.index and aa_col_aa in pivot_aa_raw.columns) else 0
                        tt = f"ABP: {aa_col_aa} ↔ Actine: {aa_row}\\n% interface: {v:.2f}%\\nAire totale: {raw_area:.0f} Å²\\nNb contacts: {n}"
                        tt_esc = tt.replace('"', '&quot;')
                        parts_aa.append(
                            f'<td style="width:{CELL}px;height:{CELL}px;{bg};color:{fg};'
                            f'text-align:center;font-size:8px;cursor:default;font-weight:600" '
                            f'data-tt2="{tt_esc}">{v:.1f}%</td>'
                        )
                    else:
                        parts_aa.append(
                            f'<td style="width:{CELL}px;height:{CELL}px;background:#f8f8f8;'
                            f'color:#ddd;text-align:center;font-size:8px">·</td>'
                        )
                parts_aa.append('</tr>')
            parts_aa.append('</tbody></table></div>')

            st.components.v1.html("".join(parts_aa), height=20 * CELL + 180, scrolling=False)

            st.markdown("---")
            st.markdown("**Contacts spécifiques par séquence ABP**")
            st.caption("Chaque ligne = une paire de résidus en contact direct (ABP ↔ actine).")

            seq_options = sorted(title_to_label.values())
            if not seq_options:
                st.info("Aucune séquence disponible.")
            else:
                sel_label = st.selectbox(
                    "Séquence ABP :", seq_options,
                    key=f"d_seq_sel_{group_key}",
                )
                sel_sl = next((sl for sl, lb in title_to_label.items() if lb == sel_label), None)
                if sel_sl:
                    df4_sel = df4[df4["seq_low"] == sel_sl].copy()
                    pairs_sel = (
                        df4_sel.groupby(["canon_b", "aa_b", "canon_a", "aa_a"])
                        .agg(aire_moy=("area_f", "mean"), n=("area_f", "count"))
                        .reset_index()
                        .sort_values("aire_moy", ascending=False)
                    )
                    pairs_sel.columns = [
                        "Pos ABP", "AA ABP", "Pos Actine", "AA Actine",
                        "Aire moy. (Å²)", "N contacts",
                    ]
                    pairs_sel["Aire moy. (Å²)"] = pairs_sel["Aire moy. (Å²)"].round(1)
                    pairs_sel["Paire"] = (
                        pairs_sel["AA ABP"] + pairs_sel["Pos ABP"].astype(str)
                        + " ↔ "
                        + pairs_sel["AA Actine"] + pairs_sel["Pos Actine"].astype(str)
                    )
                    cols_order = ["Paire", "Pos ABP", "AA ABP", "Pos Actine", "AA Actine",
                                  "Aire moy. (Å²)", "N contacts"]
                    st.dataframe(
                        pairs_sel[cols_order].reset_index(drop=True),
                        use_container_width=True, hide_index=True,
                    )
                    st.caption(
                        f"{len(pairs_sel)} paires de résidus · "
                        f"{pairs_sel['Pos ABP'].nunique()} positions ABP · "
                        f"{pairs_sel['Pos Actine'].nunique()} positions actine"
                    )

    # ── TAB E : Spécificité ───────────────────────────────────────────────────
    if "E" in _tab_map:
        with _tab_map["E"]:
            st.caption(
                "Pour chaque position canonique ABP en contact, on compare l'AA de chaque séquence. "
                "Un résidu est **unique** s'il n'apparaît que dans une seule séquence à cette position. "
                "Ces résidus sont candidats pour expliquer la spécificité d'interaction avec l'actine."
            )

            pos_aa = (
                df4.groupby(["label", "canon_b", "aa_b"])["area_f"]
                .sum().reset_index()
                .sort_values("area_f", ascending=False)
                .drop_duplicates(["label", "canon_b"])
                [["label", "canon_b", "aa_b"]]
            )
            pivot_sp = pos_aa.pivot(index="label", columns="canon_b", values="aa_b")
            pivot_sp = pivot_sp.loc[sorted(pivot_sp.index, key=lambda lb: _label_rank.get(lb, 9999))]

            all_labels = list(pivot_sp.index)
            n_seqs     = len(all_labels)
            all_pos    = list(pivot_sp.columns)

            pos_to_actin: dict = (
                df4.groupby(["label", "canon_b"])
                .apply(lambda g: sorted(g["canon_a"].unique().astype(int).tolist()))
                .to_dict()
            )
            pos_to_actin_aa: dict = (
                df4.groupby(["label", "canon_b"])
                .apply(lambda g: g.groupby("canon_a")["aa_a"].agg(
                    lambda x: x.value_counts().index[0]
                ).to_dict())
                .to_dict()
            )

            def _uniq_score(col_vals):
                aa_counts: dict = {}
                for lbl, aa in col_vals:
                    if pd.notna(aa):
                        aa_counts[lbl] = aa
                freq: dict = {}
                for aa in aa_counts.values():
                    freq[aa] = freq.get(aa, 0) + 1
                return aa_counts, freq

            def _spec_bg(aa, freq, n_total):
                if aa is None or pd.isna(aa):
                    return "background:#f5f5f5", "#ccc"
                c = freq.get(aa, 0)
                if c == 1:
                    return "background:#e74c3c", "#fff"
                if c == 2:
                    return "background:#e67e22", "#fff"
                if c <= n_total // 2:
                    return "background:#f1c40f", "#333"
                return "background:#95a5a6", "#fff"

            col_step_sp = max(1, len(all_pos) // 20)
            CELL_SP     = 18
            LABEL_SP    = 200

            js_sp = """
<div id="htt3" style="position:fixed;background:rgba(20,20,20,0.9);color:#fff;
padding:8px 12px;border-radius:8px;font-size:11px;line-height:1.7;
pointer-events:none;display:none;z-index:99999;max-width:320px;
font-family:'Helvetica Neue',Arial,sans-serif;white-space:pre-wrap"></div>
<script>
(function(){
var t=document.getElementById('htt3');
document.querySelectorAll('[data-sp]').forEach(function(el){
  el.addEventListener('mouseenter',function(e){
    t.innerHTML=el.getAttribute('data-sp');t.style.display='block';
    t.style.left=(e.clientX+16)+'px';t.style.top=(e.clientY+16)+'px';
  });
  el.addEventListener('mousemove',function(e){
    t.style.left=(e.clientX+16)+'px';t.style.top=(e.clientY+16)+'px';
  });
  el.addEventListener('mouseleave',function(){t.style.display='none';});
});
})();
</script>"""

            parts_sp = [
                '<div style="font-family:\'Helvetica Neue\',Arial,sans-serif;font-size:11px;'
                'background:#fff;padding:12px;overflow-x:auto;">',
                '<div style="font-weight:600;font-size:12px;margin-bottom:3px">'
                'Spécificité des résidus ABP aux positions de contact</div>',
                '<div style="font-size:10px;color:#888;margin-bottom:6px">'
                'Lignes = séquences · Colonnes = positions canoniques ABP en contact avec l\'actine · '
                'Survol = acide aminé + positions actine contactées</div>',
                '<div style="display:flex;align-items:center;gap:8px;margin-bottom:10px;font-size:10px">',
                '<span style="background:#e74c3c;color:#fff;padding:1px 7px;border-radius:3px">Unique (1/n)</span>',
                '<span style="background:#e67e22;color:#fff;padding:1px 7px;border-radius:3px">Rare (2/n)</span>',
                '<span style="background:#f1c40f;color:#333;padding:1px 7px;border-radius:3px">Minorité</span>',
                '<span style="background:#95a5a6;color:#fff;padding:1px 7px;border-radius:3px">Majoritaire</span>',
                '<span style="background:#f5f5f5;color:#ccc;padding:1px 7px;border-radius:3px">Absent</span>',
                '</div>',
                js_sp,
                '<table style="border-collapse:collapse;table-layout:fixed"><thead><tr>',
                f'<th style="width:{LABEL_SP}px;min-width:{LABEL_SP}px"></th>',
            ]
            for j, p in enumerate(all_pos):
                lbl_p = str(p) if j % col_step_sp == 0 else ""
                parts_sp.append(
                    f'<th style="width:{CELL_SP}px;min-width:{CELL_SP}px;font-size:7px;font-weight:normal;'
                    f'color:#aaa;writing-mode:vertical-rl;text-align:left;padding:0;'
                    f'height:26px;vertical-align:bottom">{lbl_p}</th>'
                )
            parts_sp.append('</tr></thead><tbody>')

            col_aa_counts: dict = {}
            col_freq: dict = {}
            for p in all_pos:
                col_vals = [(lb, pivot_sp.at[lb, p] if lb in pivot_sp.index and p in pivot_sp.columns else None)
                            for lb in all_labels]
                aa_c, fr = _uniq_score(col_vals)
                col_aa_counts[p] = aa_c
                col_freq[p] = fr

            for lb in all_labels:
                short = (lb[:28] + "…") if len(lb) > 29 else lb
                parts_sp.append('<tr>')
                parts_sp.append(
                    f'<td style="width:{LABEL_SP}px;max-width:{LABEL_SP}px;overflow:hidden;'
                    f'white-space:nowrap;text-overflow:ellipsis;font-size:9px;color:#444;'
                    f'padding:1px 6px 1px 0;text-align:right" title="{lb}">{short}</td>'
                )
                for p in all_pos:
                    aa   = col_aa_counts[p].get(lb)
                    freq = col_freq[p]
                    bg, fg = _spec_bg(aa, freq, n_seqs)
                    if aa and not pd.isna(aa):
                        c_aa = freq.get(aa, 0)
                        actin_pos_list = pos_to_actin.get((lb, p), [])
                        actin_aa_map   = pos_to_actin_aa.get((lb, p), {})
                        actin_str = ", ".join(
                            f"{actin_aa_map.get(ap,'?')}{ap}" for ap in actin_pos_list[:6]
                        )
                        if len(actin_pos_list) > 6:
                            actin_str += f" (+{len(actin_pos_list)-6})"
                        uniq_lbl = (
                            "UNIQUE" if c_aa == 1
                            else f"rare ({c_aa}/{n_seqs})" if c_aa <= 2
                            else f"minorité ({c_aa}/{n_seqs})" if c_aa <= n_seqs // 2
                            else f"majoritaire ({c_aa}/{n_seqs})"
                        )
                        tt = (
                            f"Pos ABP: {p} | {lb[:28]}\n"
                            f"AA ABP: {aa} — {uniq_lbl}\n"
                            f"Contact actine: {actin_str if actin_str else '?'}"
                        )
                        others = [f"{al}:{col_aa_counts[p].get(al,'–')}"
                                  for al in all_labels if al != lb and col_aa_counts[p].get(al)]
                        if others:
                            tt += "\nAutres: " + ", ".join(others[:4])
                            if len(others) > 4:
                                tt += f" (+{len(others)-4})"
                        tt_esc = tt.replace('"', '&quot;')
                        parts_sp.append(
                            f'<td style="width:{CELL_SP}px;min-width:{CELL_SP}px;height:{CELL_SP}px;'
                            f'{bg};color:{fg};text-align:center;font-size:8px;font-weight:600;'
                            f'cursor:default;padding:0" data-sp="{tt_esc}">{aa}</td>'
                        )
                    else:
                        parts_sp.append(
                            f'<td style="width:{CELL_SP}px;min-width:{CELL_SP}px;height:{CELL_SP}px;'
                            f'background:#f5f5f5;color:#ddd;text-align:center;font-size:8px;'
                            f'padding:0">·</td>'
                        )
                parts_sp.append('</tr>')
            parts_sp.append('</tbody></table>')

            uniq_rows = []
            for p in all_pos:
                for lb, aa in col_aa_counts[p].items():
                    if aa and not pd.isna(aa) and col_freq[p].get(aa, 0) == 1:
                        ap_list = pos_to_actin.get((lb, p), [])
                        ap_aa   = pos_to_actin_aa.get((lb, p), {})
                        actin_s = ", ".join(f"{ap_aa.get(a,'?')}{a}" for a in ap_list[:5])
                        uniq_rows.append({
                            "Séquence ABP": lb,
                            "Pos ABP": p,
                            "AA ABP unique": aa,
                            "Contacte actine (AA+pos)": actin_s,
                        })

            parts_sp.append('</div>')
            _ht_sp = max(n_seqs * CELL_SP + 120, 250)
            st.components.v1.html("".join(parts_sp), height=_ht_sp, scrolling=True)

            if uniq_rows:
                st.markdown(f"**{len(uniq_rows)} résidus ABP uniques** (AA absent de toutes les autres séquences à cette position)")
                df_uniq = pd.DataFrame(uniq_rows).sort_values(["Séquence ABP", "Pos ABP"])
                st.dataframe(df_uniq.reset_index(drop=True), use_container_width=True, hide_index=True)
            else:
                st.info("Aucun résidu ABP strictement unique — tous les AA de contact sont partagés par ≥2 séquences.")

            st.divider()
            st.markdown("#### Côté actine — positions canoniques contactées spécifiquement")
            st.caption(
                "Lignes = séquences ABP · Colonnes = positions canoniques actine · "
                "Lettre = AA actine à cette position · "
                "Couleur = nombre de séquences qui contactent cette position : "
                "**rouge** = unique (1 seule myosine), orange = rare, jaune = minorité, gris = partagé."
            )

            actin_aa_specific: dict = (
                df4.groupby(["label", "canon_a"])["aa_a"]
                .agg(lambda x: x.value_counts().index[0] if len(x) > 0 else "?")
                .to_dict()
            )
            pres_a = (
                df4.groupby(["label", "canon_a"])["area_f"]
                .sum()
                .gt(0)
                .reset_index()
                .rename(columns={"area_f": "has_contact"})
            )
            contact_count_a: dict = pres_a[pres_a["has_contact"]].groupby("canon_a")["label"].nunique().to_dict()

            abp_per_actin: dict = (
                pres_a[pres_a["has_contact"]].groupby("canon_a")["label"]
                .apply(list).to_dict()
            )

            all_pos_a = sorted(contact_count_a.keys())

            def _spec_bg_a(n_contact, n_total):
                if n_contact == 0:
                    return "background:#f5f5f5", "#ccc"
                if n_contact == 1:
                    return "background:#e74c3c", "#fff"
                if n_contact == 2:
                    return "background:#e67e22", "#fff"
                if n_contact <= n_total // 2:
                    return "background:#f1c40f", "#333"
                return "background:#95a5a6", "#fff"

            col_step_a = max(1, len(all_pos_a) // 20)
            js_spa = js_sp.replace("htt3", "htt4").replace("data-sp", "data-spa")

            parts_spa = [
                '<div style="font-family:\'Helvetica Neue\',Arial,sans-serif;font-size:11px;'
                'background:#fff;padding:12px;overflow-x:auto;">',
                js_spa,
                '<table style="border-collapse:collapse;table-layout:fixed"><thead><tr>',
                f'<th style="width:{LABEL_SP}px;min-width:{LABEL_SP}px"></th>',
            ]
            for j, p in enumerate(all_pos_a):
                lbl_p = str(p) if j % col_step_a == 0 else ""
                parts_spa.append(
                    f'<th style="width:{CELL_SP}px;min-width:{CELL_SP}px;font-size:7px;font-weight:normal;'
                    f'color:#aaa;writing-mode:vertical-rl;text-align:left;padding:0;'
                    f'height:26px;vertical-align:bottom">{lbl_p}</th>'
                )
            parts_spa.append('</tr></thead><tbody>')

            pres_dict: dict = {
                (r["label"], int(r["canon_a"])): r["has_contact"]
                for _, r in pres_a.iterrows()
            }

            for lb in all_labels:
                short = (lb[:28] + "…") if len(lb) > 29 else lb
                parts_spa.append('<tr>')
                parts_spa.append(
                    f'<td style="width:{LABEL_SP}px;max-width:{LABEL_SP}px;overflow:hidden;'
                    f'white-space:nowrap;text-overflow:ellipsis;font-size:9px;color:#444;'
                    f'padding:1px 6px 1px 0;text-align:right" title="{lb}">{short}</td>'
                )
                for p in all_pos_a:
                    has = pres_dict.get((lb, int(p)), False)
                    n_c = contact_count_a.get(p, 0)
                    aa_act = actin_aa_specific.get((lb, p), "?")
                    bg, fg = _spec_bg_a(n_c if has else 0, n_seqs)
                    if has:
                        others = [lb2 for lb2 in abp_per_actin.get(p, []) if lb2 != lb]
                        uniq_lbl = (
                            "UNIQUE" if n_c == 1
                            else f"rare ({n_c}/{n_seqs})" if n_c <= 2
                            else f"minorité ({n_c}/{n_seqs})" if n_c <= n_seqs // 2
                            else f"partagé ({n_c}/{n_seqs})"
                        )
                        tt = (
                            f"Pos actine : {p}\n"
                            f"AA actine  : {aa_act}\n"
                            f"Contact    : {lb[:28]} — {uniq_lbl}"
                        )
                        if others:
                            tt += "\nAussi : " + ", ".join(o[:20] for o in others[:4])
                        tt_esc = tt.replace('"', '&quot;')
                        parts_spa.append(
                            f'<td style="width:{CELL_SP}px;min-width:{CELL_SP}px;height:{CELL_SP}px;'
                            f'{bg};color:{fg};text-align:center;font-size:8px;font-weight:600;'
                            f'cursor:default;padding:0" data-spa="{tt_esc}">{aa_act}</td>'
                        )
                    else:
                        parts_spa.append(
                            f'<td style="width:{CELL_SP}px;min-width:{CELL_SP}px;height:{CELL_SP}px;'
                            f'background:#f5f5f5;color:#ddd;text-align:center;font-size:8px;'
                            f'padding:0">·</td>'
                        )
                parts_spa.append('</tr>')
            parts_spa.append('</tbody></table></div>')

            _html_spa = "".join(parts_spa).replace("data-sp\"", "data-spa\"")
            st.components.v1.html(_html_spa, height=max(n_seqs * CELL_SP + 120, 250), scrolling=True)

            uniq_rows_a = []
            for p in all_pos_a:
                if contact_count_a.get(p, 0) == 1:
                    lb = abp_per_actin[p][0]
                    aa_act = actin_aa_specific.get((lb, p), "?")
                    aa_abp_here = df4[(df4["label"] == lb) & (df4["canon_a"] == p)]["aa_b"].value_counts()
                    aa_abp_str = ", ".join(f"{aa}({cnt})" for aa, cnt in list(aa_abp_here.items())[:3])
                    uniq_rows_a.append({
                        "Pos actine": p,
                        "AA actine": aa_act,
                        "Contactée uniquement par": lb,
                        "AA ABP correspondants": aa_abp_str,
                    })

            if uniq_rows_a:
                st.markdown(f"**{len(uniq_rows_a)} positions actine contactées par 1 seule séquence ABP**")
                df_uniq_a = pd.DataFrame(uniq_rows_a).sort_values("Pos actine")
                st.dataframe(df_uniq_a.reset_index(drop=True), use_container_width=True, hide_index=True)
            else:
                st.info("Toutes les positions actine sont contactées par ≥2 séquences ABP.")


# ── Comparaison C. elegans myosines ───────────────────────────────────────────

def _s1_get_ch_maps(cid: str):
    """Retourne (ch2seqlow, ch2title) pour un cluster S1 donné.
    Utilisé par streamlit.py pour afficher les onglets B–E dans la section Binding Site.
    Inclut tous les partenaires subunit_2 sans filtre (actin pour clusters homo,
    ABP pour clusters hetero, les deux pour clusters mixed).
    """
    filt_path = _Path("data/filtered/filtered_all_data.csv")
    if not filt_path.exists():
        return {}, {}
    df_filt = pd.read_csv(filt_path, low_memory=False)

    df_pairs = df_filt[df_filt["s1_binding_site_cluster_data_70"].astype(str) == str(cid)][
        ["subunit_2", "subunit_2_title", "s2_sequence"]
    ].drop_duplicates("subunit_2")
    if df_pairs.empty:
        return {}, {}
    ch2seqlow = {r["subunit_2"]: str(r["s2_sequence"]).strip().lower() for _, r in df_pairs.iterrows()}
    ch2title  = {r["subunit_2"]: str(r["subunit_2_title"])             for _, r in df_pairs.iterrows()}
    return ch2seqlow, ch2title


def _msa_celegans_comparison():
    """Affiche la comparaison des myosines C. elegans (myo-1/2/3) avec les séquences du dataset."""
    ALN_PATH = _MSA_ALN_DIR / "myosin_celegans_msa.aln"
    if not ALN_PATH.exists():
        st.info("Fichier myosin_celegans_msa.aln manquant — relancer MAFFT avec les séquences C. elegans.")
        return

    try:
        aln = AlignIO.read(str(ALN_PATH), "fasta")
    except Exception as _e:
        st.error(f"Erreur lecture alignement C. elegans : {_e}")
        return

    # Colonnes de liaison actine calculées à partir de MYO1_CAEEL (pos UniProt 660-682 et 764-778)
    # ↔ domaine moteur résidus 577-599 et 681-695 → colonnes alignement 735-757 et 841-858
    _AB_COLS: set = set(range(735, 758)) | set(range(841, 859))

    celegans_recs = [r for r in aln if "CAEEL" in r.id]
    other_recs    = [r for r in aln if "CAEEL" not in r.id]
    all_recs      = celegans_recs + other_recs

    # ── Tableau de similarité ─────────────────────────────────────────────────
    st.markdown("**Identité de séquence — domaine moteur (%) :**")
    _CLEAN = lambda s: (s.replace("_motor_84-785", "")
                         .replace("_CAEEL", "")
                         .replace("_", " ")
                         .replace("  ", " "))
    sim_rows = []
    for crec in celegans_recs:
        row: dict = {"C. elegans": _CLEAN(crec.id)}
        best_val, best_name = 0.0, ""
        for orec in other_recs:
            s1, s2 = str(crec.seq), str(orec.seq)
            matches = sum(1 for a, b in zip(s1, s2) if a not in "-" and b not in "-" and a == b)
            aligned = sum(1 for a, b in zip(s1, s2) if a not in "-" and b not in "-")
            pct = matches / aligned * 100 if aligned > 0 else 0.0
            row[_CLEAN(orec.id)] = f"{pct:.1f}%"
            if pct > best_val:
                best_val, best_name = pct, _CLEAN(orec.id)
        row["→ Plus proche"] = f"{best_name} ({best_val:.1f}%)"
        sim_rows.append(row)
    df_sim = pd.DataFrame(sim_rows).set_index("C. elegans")
    st.dataframe(df_sim, use_container_width=True)

    # ── MSA HTML ─────────────────────────────────────────────────────────────
    st.markdown(
        "**Alignement — zones de liaison actine en gris** "
        "(UniProt myo-1 : résidus 660–682 et 764–778) :"
    )
    aln_len = aln.get_alignment_length()
    label_w = min(max(len(_CLEAN(r.id)) for r in all_recs), 38)

    STYLE: dict = {
        "ab_ce":    "background:#95a5a6;color:#fff;font-weight:700",
        "ab_ot":    "background:#dfe6e9;color:#636e72",
        "ce":       "color:#2980b9",
        "ot":       "color:#b2bec3",
        "gap":      "color:#ececec",
        "ab_gap":   "background:#dfe6e9;color:#b2bec3",
    }

    def _span(style_key: str, text: str) -> str:
        return f'<span style="{STYLE[style_key]}">{text}</span>'

    parts = [
        '<div style="font-family:\'Courier New\',Courier,monospace;font-size:10px;'
        'line-height:1.7;background:#fafafa;padding:10px;overflow-x:auto;white-space:nowrap;">'
    ]

    # Barre d'annotation (en-tête)
    bar = "&nbsp;" * (label_w + 2)
    for col in range(aln_len):
        if col in _AB_COLS:
            bar += '<span style="background:#7f8c8d;color:#fff">▪</span>'
        else:
            bar += '<span style="color:#eee">·</span>'
    parts.append(f"<div>{bar}</div>")

    parts.append('<table style="border-collapse:collapse;margin-top:4px">')
    for rec in all_recs:
        is_ce = "CAEEL" in rec.id
        label = _CLEAN(rec.id)[:label_w].ljust(label_w).replace(" ", "&nbsp;")
        if is_ce:
            lbl_html = f'<b style="color:#2980b9">{label}</b>'
        else:
            lbl_html = f'<span style="color:#636e72">{label}</span>'

        parts.append(f'<tr><td style="padding-right:6px;white-space:nowrap">{lbl_html}</td><td>')

        seq_str = str(rec.seq)
        segs: list = []
        cur_key: str | None = None
        buf: list = []

        def _flush():
            if buf and cur_key:
                segs.append(_span(cur_key, "".join(buf)))

        for col, aa in enumerate(seq_str):
            in_ab = col in _AB_COLS
            if aa == "-":
                new_key = "ab_gap" if in_ab else "gap"
            else:
                if in_ab:
                    new_key = "ab_ce" if is_ce else "ab_ot"
                else:
                    new_key = "ce" if is_ce else "ot"

            if new_key != cur_key:
                _flush()
                buf = []
                cur_key = new_key
            buf.append(aa.upper() if aa != "-" else "-")

        _flush()
        parts.append("".join(segs))
        parts.append("</td></tr>")

    parts.append("</table></div>")

    height = max(len(all_recs) * 20 + 120, 300)
    st.components.v1.html("".join(parts), height=height, scrolling=True)


# ── UI MSA ────────────────────────────────────────────────────────────────────

def _msa_section_full(group_label, group_key, filter_fn, rigor_pdbs=None, note=None):
    """Section MSA séquence COMPLÈTE — ABP (S2) et Actine (S1) en colonnes parallèles."""
    with st.expander(f"**{group_label}**", expanded=(group_key == "myosin")):
        if note:
            st.caption(note)

        # Calcul des contacts actine et chargement des séquences ABP en avance.
        # Ordre partagé : S2 trié par titre, S1 trié par rang de la séquence dans S2.
        _abp_rows_a, _aa_at_a, _col_col_a = _msa_actin_contacts_per_abp(filter_fn, rigor_pdbs)

        df_seqs = _msa_load_seqs_g(filter_fn, rigor_pdbs)
        if not df_seqs.empty:
            df_seqs = df_seqs.sort_values("title").reset_index(drop=True)

        # Construire label_map : seq_id → titre propre (avec suffixe si doublons)
        _label_map: dict = {}
        if not df_seqs.empty:
            from collections import Counter as _Ctr2
            _tcounts = _Ctr2(df_seqs["title"])
            _tseen: dict = {}
            for _, _row in df_seqs.iterrows():
                _t = _row["title"]
                if _tcounts[_t] > 1:
                    _n = _tseen.get(_t, 1)
                    _label_map[_row["seq_id"]] = f"{_t} ({_n})"
                    _tseen[_t] = _n + 1
                else:
                    _label_map[_row["seq_id"]] = _t
        if group_key == "myosin":
            _label_map["MYO1_CAEEL_motor_84-785"] = "myo-1 (C. elegans)"
            _label_map["MYO2_CAEEL_motor_84-785"] = "myo-2 (C. elegans)"
            _label_map["MYO3_CAEEL_motor_84-785"] = "myo-3 (C. elegans)"

        # Ordre S1 : lire l'ordre réel du fichier .aln (MAFFT --reorder change l'ordre d'entrée).
        # Mapping : seq_id[:50] (ID MAFFT) → séquence lower → rang dans l'alignement.
        _aln_path_check = _MSA_ALN_DIR / f"{group_key}_msa.aln"
        _sid_to_seq = (
            {r["seq_id"][:50]: r["seq"].strip().lower() for _, r in df_seqs.iterrows()}
            if not df_seqs.empty else {}
        )
        if _aln_path_check.exists() and _sid_to_seq:
            try:
                _aln_ord = AlignIO.read(str(_aln_path_check), "fasta")
                _s2seq_rank = {
                    _sid_to_seq[rec.id]: i
                    for i, rec in enumerate(_aln_ord)
                    if rec.id in _sid_to_seq
                }
                _abp_rows_a = sorted(_abp_rows_a,
                                     key=lambda r: _s2seq_rank.get(r.get("s2seq", ""), 9999))
            except Exception:
                _abp_rows_a = sorted(_abp_rows_a, key=lambda r: r["title"])
        else:
            _abp_rows_a = sorted(_abp_rows_a, key=lambda r: r["title"])

        # ── S2 : ABP — séquence complète ──────────────────────────────────────
        st.markdown("##### 🔵 ABP — séquence complète (S2)")
        if df_seqs.empty:
            st.info("Aucune séquence trouvée.")
        else:
            st.caption(
                f"**{len(df_seqs)} séquences uniques** "
                f"({df_seqs['length'].min()}–{df_seqs['length'].max()} aa)"
            )
            if len(df_seqs) < 2:
                st.info("Moins de 2 séquences — alignement impossible.")
            else:
                fasta_path = _MSA_ALN_DIR / f"{group_key}_msa.fasta"
                aln_path   = _MSA_ALN_DIR / f"{group_key}_msa.aln"
                _btn_c, _force_c = st.columns([1, 2])
                _run   = _btn_c.button("Lancer MAFFT", key=f"msa_{group_key}_btn")
                _force = _force_c.checkbox("Forcer recalcul", key=f"msa_{group_key}_force")

                if _run or aln_path.exists():
                    if _run or not aln_path.exists() or _force:
                        _MSA_ALN_DIR.mkdir(parents=True, exist_ok=True)
                        SeqIO.write(
                            [SeqRecord(Seq(r["seq"]), id=r["seq_id"][:50], description="")
                             for _, r in df_seqs.iterrows()],
                            fasta_path, "fasta",
                        )
                        with st.spinner(f"MAFFT sur {len(df_seqs)} séquences…"):
                            _ok, _err = _msa_run_mafft(fasta_path, aln_path)
                        if not _ok:
                            st.error(f"Erreur MAFFT : {_err}")

                    if aln_path.exists():
                        try:
                            _aln = AlignIO.read(str(aln_path), "fasta")
                            # Pour myosine : utiliser l'alignement avec C. elegans si disponible
                            if group_key == "myosin":
                                _ce_path = _MSA_ALN_DIR / "myosin_celegans_msa.aln"
                                if _ce_path.exists():
                                    try:
                                        from Bio.Align import MultipleSeqAlignment
                                        _aln_ce = AlignIO.read(str(_ce_path), "fasta")
                                        _ce_recs  = [r for r in _aln_ce if "CAEEL" in r.id]
                                        _oth_recs = [r for r in _aln_ce if "CAEEL" not in r.id]
                                        # Tri par similarité d'interface (mêmes colonnes que le tableau)
                                        _core_pre, _var_pre = _msa_load_interface_g(
                                            filter_fn, rigor_pdbs)
                                        def _s2aln_map(seq_str):
                                            _m, _sp = {}, 0
                                            for _c, _a in enumerate(seq_str):
                                                if _a != "-":
                                                    _sp += 1
                                                    _m[_sp] = _c
                                            return _m
                                        def _iface_pid(rec_s, ce_s):
                                            _ug = rec_s.replace("-", "").lower()
                                            _ipos = (
                                                _core_pre.get(_ug, set())
                                                | _var_pre.get(_ug, set())
                                            )
                                            if not _ipos:
                                                _ps = [(a, b) for a, b in zip(rec_s, ce_s)
                                                       if a != "-" and b != "-"]
                                                return sum(a == b for a, b in _ps) / max(len(_ps), 1)
                                            _icols = {
                                                _s2aln_map(rec_s)[sp]
                                                for sp in _ipos
                                                if sp in _s2aln_map(rec_s)
                                            }
                                            _m = _p = 0
                                            for _col in _icols:
                                                _a = rec_s[_col] if _col < len(rec_s) else "-"
                                                _b = ce_s[_col] if _col < len(ce_s) else "-"
                                                if _a != "-" and _b != "-":
                                                    _p += 1
                                                    if _a.upper() == _b.upper():
                                                        _m += 1
                                            return _m / max(_p, 1)
                                        # Pour chaque CAEEL → indice du dataset le plus proche
                                        _ce_to_oth = {}
                                        for _cr in _ce_recs:
                                            _best_i, _best_v = 0, -1.0
                                            for _oi, _or in enumerate(_oth_recs):
                                                _v = _iface_pid(str(_or.seq), str(_cr.seq))
                                                if _v > _best_v:
                                                    _best_v, _best_i = _v, _oi
                                            _ce_to_oth[_cr.id] = _best_i
                                        _after: dict = {}
                                        for _cr in _ce_recs:
                                            _after.setdefault(_ce_to_oth[_cr.id], []).append(_cr)
                                        _ordered = []
                                        for _oi, _or in enumerate(_oth_recs):
                                            _ordered.append(_or)
                                            _ordered.extend(_after.get(_oi, []))
                                        _aln = MultipleSeqAlignment(_ordered)
                                    except Exception:
                                        pass
                        except Exception as _e:
                            st.error(f"Impossible de lire l'alignement : {_e}")
                        else:
                            _nseqs = len(_aln)
                            _alen  = _aln.get_alignment_length()
                            st.success(f"**{_nseqs} séq. × {_alen} col.**")
                            # Re-trier S1 selon l'ordre réel de cet alignement
                            _s2seq_rank_fresh = {
                                _sid_to_seq[rec.id]: i
                                for i, rec in enumerate(_aln)
                                if rec.id in _sid_to_seq
                            }
                            if _s2seq_rank_fresh:
                                _abp_rows_a = sorted(
                                    _abp_rows_a,
                                    key=lambda r: _s2seq_rank_fresh.get(r.get("s2seq", ""), 9999),
                                )
                            _core, _var = _msa_load_interface_g(filter_fn, rigor_pdbs)
                            _n_core = sum(len(v) for v in _core.values())
                            _n_var  = sum(len(v) for v in _var.values())
                            st.caption(f"{_n_core} résidus core · {_n_var} variables")
                            _html   = _msa_render_full(_aln, _core, _var, 9999, _label_map)
                            _height = min(_nseqs * 18 + 80, 6000)
                            st.components.v1.html(_html, height=_height, scrolling=True)
                            with open(aln_path, "rb") as _f:
                                st.download_button(
                                    f"Télécharger {group_key}_msa.aln", _f,
                                    file_name=f"{group_key}_msa.aln", mime="text/plain",
                                    key=f"msa_{group_key}_dl",
                                )
                            if group_key == "myosin":
                                _caeel_recs = {
                                    r.id: str(r.seq)
                                    for r in _aln if "CAEEL" in r.id
                                }
                                if _caeel_recs:
                                    def _build_s2aln(seq_str):
                                        """seq_pos (1-based) → aln_col"""
                                        _m, _sp = {}, 0
                                        for _c, _a in enumerate(seq_str):
                                            if _a != "-":
                                                _sp += 1
                                                _m[_sp] = _c
                                        return _m
                                    _pid_rows = []
                                    for _rec in _aln:
                                        if "CAEEL" in _rec.id:
                                            continue
                                        _seq_s = str(_rec.seq)
                                        _ungap_low = _seq_s.replace("-", "").lower()
                                        _iface_pos = (
                                            _core.get(_ungap_low, set())
                                            | _var.get(_ungap_low, set())
                                        )
                                        if not _iface_pos:
                                            continue
                                        # positions séquence → colonnes alignement
                                        _s2aln = _build_s2aln(_seq_s)
                                        _iface_cols = {
                                            _s2aln[sp] for sp in _iface_pos if sp in _s2aln
                                        }
                                        if not _iface_cols:
                                            continue
                                        _row = {"Séquence": _label_map.get(_rec.id, _rec.id)}
                                        for _ce_id, _ce_seq in _caeel_recs.items():
                                            _matches = _pairs = 0
                                            for _col in _iface_cols:
                                                _aa_s = _seq_s[_col] if _col < len(_seq_s) else "-"
                                                _aa_c = _ce_seq[_col] if _col < len(_ce_seq) else "-"
                                                if _aa_s != "-" and _aa_c != "-":
                                                    _pairs += 1
                                                    if _aa_s.upper() == _aa_c.upper():
                                                        _matches += 1
                                            _pct = f"{_matches / _pairs * 100:.1f}%" if _pairs > 0 else "N/A"
                                            _row[_label_map.get(_ce_id, _ce_id)] = _pct
                                        _pid_rows.append(_row)
                                    if _pid_rows:
                                        st.divider()
                                        st.markdown(
                                            "##### 🪱 Identité sur résidus d'interface — vs *C. elegans*"
                                        )
                                        st.caption(
                                            "% d'acides aminés identiques uniquement sur les colonnes "
                                            "colorées (résidus en contact avec l'actine)."
                                        )
                                        st.dataframe(
                                            pd.DataFrame(_pid_rows).set_index("Séquence"),
                                            use_container_width=True,
                                        )

        # ── S1 : Actine — positions canoniques ────────────────────────────────
        st.divider()
        st.markdown("##### 🔴 Actine — positions canoniques (S1)")
        st.caption(
            "Même ordre de lignes que S2. "
            "Colonnes = positions canoniques actine. "
            "Scroll horizontal pour voir toute la séquence."
        )
        if not _abp_rows_a:
            st.info("Aucune donnée d'interface actine trouvée.")
        else:
            st.caption(
                f"**{len(_abp_rows_a)} séquences ABP** · "
                f"{len(_aa_at_a)} positions canoniques"
            )
            # Label map S1 : même noms propres que S2
            from collections import Counter as _Ctr3
            _s1_tcounts = _Ctr3(r["title"] for r in _abp_rows_a)
            _s1_tseen: dict = {}
            _label_map_s1: dict = {}
            for _r in _abp_rows_a:
                _t = _r["title"]
                if _s1_tcounts[_t] > 1:
                    _n = _s1_tseen.get(_t, 1)
                    _label_map_s1[_r["seq_id"]] = f"{_t} ({_n})"
                    _s1_tseen[_t] = _n + 1
                else:
                    _label_map_s1[_r["seq_id"]] = _t
            _html_a   = _msa_render_actin_contacts(_abp_rows_a, _aa_at_a, _col_col_a, 9999, _label_map_s1)
            _height_a = min(len(_abp_rows_a) * 18 + 80, 6000)
            st.components.v1.html(_html_a, height=_height_a, scrolling=True)

        # ── Analyse des contacts ABP–actine (onglets A + B C D E) ──
        st.divider()
        st.markdown("##### 📊 Analyse des contacts ABP–actine")
        _msa_contact_analysis(filter_fn, group_key, rigor_pdbs, tabs=["A", "B", "C", "D", "E"])


def _msa_section_s2_clusters():
    """
    Clusters S2 (s2_sequence_cluster_70) hors 3 familles nommées.
    MAFFT sur séquence COMPLÈTE + interface colorée — identique aux 3 familles.
    """
    _EXCLUDE_FN = lambda t: (
        ("myosin" in t.lower() and "tropomyosin" not in t.lower())
        or "tropomyosin" in t.lower()
        or any(x in t.lower() for x in ["plastin", "spectrin beta", "filamin", "utrophin"])
        or t.lower().startswith("actin,")
        or t.lower() in ("actin", "actin-1")
        or (t.lower().startswith("actin") and "actin-related" not in t.lower()
            and "actin-depolymerizing" not in t.lower())
    )

    with st.expander("**Autres clusters de protéines — séquence complète**"):
        st.caption("MAFFT sur séquence complète, résidus d'interface colorés (s2_sequence_cluster_70).")

        filt_path = _Path("data/filtered/filtered_all_data.csv")
        int3_path = _Path("data/filtered/details/3.interface_residues.csv")
        if not filt_path.exists() or not int3_path.exists():
            st.warning("Données manquantes.")
            return

        from collections import defaultdict
        _int1_path_s2cl = _Path("data/filtered/details/1.interactions.csv")
        df_filt = pd.read_csv(filt_path, low_memory=False)
        df3     = pd.read_csv(int3_path)
        _df1_s2cl = pd.read_csv(_int1_path_s2cl) if _int1_path_s2cl.exists() else pd.DataFrame()
        _df3_s2cl = df3.copy()
        _df3_s2cl["residue_number_canon_mafft"] = pd.to_numeric(
            _df3_s2cl["residue_number_canon_mafft"], errors="coerce")
        _df3_s2cl["chain_lower"] = _df3_s2cl["chain"].str.lower()

        mask    = ~df_filt["subunit_2_title"].apply(lambda t: _EXCLUDE_FN(str(t)))
        df_other = df_filt[mask].copy()

        clusters = (
            df_other[["s2_sequence_cluster_70", "subunit_2_title"]]
            .dropna(subset=["s2_sequence_cluster_70"]).drop_duplicates()
            .groupby("s2_sequence_cluster_70")["subunit_2_title"]
            .apply(lambda x: ", ".join(sorted(set(x))[:3])).reset_index()
        )
        clusters.columns = ["cluster_id", "titles"]

        for _, crow in clusters.sort_values("cluster_id").iterrows():
            cid   = int(crow["cluster_id"])
            ckey  = f"s2c_{cid}"
            rows_c = df_other[df_other["s2_sequence_cluster_70"] == cid][
                ["subunit_2", "subunit_2_title", "s2_sequence", "s2_taxonomy_id"]
            ].drop_duplicates(subset=["subunit_2"])

            chain_to_full   = {r["subunit_2"]: str(r["s2_sequence"]).strip() for _, r in rows_c.iterrows()}
            chain_to_seqlow = {ch: s.lower() for ch, s in chain_to_full.items()}
            chain_to_title  = {r["subunit_2"]: str(r["subunit_2_title"]) for _, r in rows_c.iterrows()}
            chain_to_taxid  = {r["subunit_2"]: r["s2_taxonomy_id"] for _, r in rows_c.iterrows()}

            # Interface positions per unique sequence
            res_c = df3[df3["chain"].isin(chain_to_full)][["chain","residue_number_sequence"]].dropna()
            seq_iface_chains: dict = defaultdict(dict)
            for _, rrow in res_c.iterrows():
                seqlow = chain_to_seqlow.get(rrow["chain"])
                if seqlow:
                    seq_iface_chains[seqlow].setdefault(rrow["chain"], set()).add(int(rrow["residue_number_sequence"]))

            # Build core/variable for _msa_render_full
            core_c: dict = {}
            var_c:  dict = {}
            for seqlow, cpos in seq_iface_chains.items():
                sets = list(cpos.values())
                if len(sets) == 1:
                    core_c[seqlow], var_c[seqlow] = sets[0].copy(), set()
                else:
                    inter = sets[0].copy(); union = sets[0].copy()
                    for s in sets[1:]: inter &= s; union |= s
                    core_c[seqlow], var_c[seqlow] = inter, union - inter

            # Unique sequences (full) with interface data
            seen: set = set()
            uniq: list = []
            for ch, seqlow in chain_to_seqlow.items():
                if seqlow not in seen and seqlow in seq_iface_chains:
                    seen.add(seqlow)
                    taxid = chain_to_taxid[ch]
                    org   = _MSA_TAX.get(int(taxid), f"taxid:{int(taxid)}") if pd.notna(taxid) else "unk"
                    title = chain_to_title[ch]
                    bid   = (title[:20].replace(" ","_").replace(",","") + "_" + org)[:40]
                    uniq.append({"seq_id": bid, "seq": chain_to_full[ch]})

            if len(uniq) < 2:
                continue

            # Actine S1 pour ce cluster
            rows_act = df_other[df_other["s2_sequence_cluster_70"] == cid][
                ["subunit_1", "subunit_1_title", "s1_sequence", "s1_taxonomy_id"]
            ].drop_duplicates(subset=["subunit_1"])
            chain_to_full_act   = {r["subunit_1"]: str(r["s1_sequence"]).strip() for _, r in rows_act.iterrows()}
            chain_to_seqlow_act = {ch: s.lower() for ch, s in chain_to_full_act.items()}
            chain_to_taxid_act  = {r["subunit_1"]: r["s1_taxonomy_id"] for _, r in rows_act.iterrows()}
            chain_to_title_act  = {r["subunit_1"]: str(r["subunit_1_title"]) for _, r in rows_act.iterrows()}
            res_act = df3[df3["chain"].isin(chain_to_full_act)][["chain", "residue_number_sequence"]].dropna()
            seq_iface_act: dict = defaultdict(dict)
            for _, rrow in res_act.iterrows():
                sl = chain_to_seqlow_act.get(rrow["chain"])
                if sl:
                    seq_iface_act[sl].setdefault(rrow["chain"], set()).add(int(rrow["residue_number_sequence"]))
            core_act: dict = {}; var_act: dict = {}
            for sl, cpos in seq_iface_act.items():
                sets = list(cpos.values())
                if len(sets) == 1:
                    core_act[sl], var_act[sl] = sets[0].copy(), set()
                else:
                    inter = sets[0].copy(); union = sets[0].copy()
                    for s in sets[1:]: inter &= s; union |= s
                    core_act[sl], var_act[sl] = inter, union - inter
            seen_act: set = set(); uniq_act: list = []
            for ch, sl in chain_to_seqlow_act.items():
                if sl not in seen_act and sl in seq_iface_act:
                    seen_act.add(sl)
                    taxid = chain_to_taxid_act[ch]
                    org   = _MSA_TAX.get(int(taxid), f"taxid:{int(taxid)}") if pd.notna(taxid) else "unk"
                    bid   = (chain_to_title_act[ch][:20].replace(" ", "_").replace(",", "") + "_" + org)[:40]
                    uniq_act.append({"seq_id": bid, "seq": chain_to_full_act[ch]})

            # Calcul contacts actine pour ce cluster (avant les colonnes, pour trier)
            _df_pairs_c2 = df_other[df_other["s2_sequence_cluster_70"] == cid][
                ["subunit_1", "subunit_2", "subunit_2_title", "s2_sequence", "s2_taxonomy_id"]
            ].copy()
            _abp_rows_c2, _aa_at_c2, _col_col_c2 = _msa_actin_contacts_from_pairs(
                _df_pairs_c2, _df1_s2cl, _df3_s2cl
            )
            uniq = sorted(uniq, key=lambda s: s["seq_id"])
            _seq_rank_c2 = {s["seq"].strip().lower(): i for i, s in enumerate(uniq)}
            _abp_rows_c2 = sorted(_abp_rows_c2,
                                  key=lambda r: _seq_rank_c2.get(r.get("s2seq", ""), 9999))

            clabel = crow["titles"][:60]
            with st.expander(f"Cluster {cid} — {clabel} ({len(uniq)} séq.)"):
                st.markdown("##### 🔵 ABP — séquence complète (S2)")
                fasta_path = _MSA_ALN_DIR / f"{ckey}_msa.fasta"
                aln_path   = _MSA_ALN_DIR / f"{ckey}_msa.aln"
                _bc, _fc = st.columns([1, 2])
                _run2   = _bc.button("Lancer MAFFT", key=f"msa_{ckey}_btn")
                _force2 = _fc.checkbox("Forcer recalcul", key=f"msa_{ckey}_force")

                if _run2 or aln_path.exists():
                    if _run2 or not aln_path.exists() or _force2:
                        _MSA_ALN_DIR.mkdir(parents=True, exist_ok=True)
                        SeqIO.write(
                            [SeqRecord(Seq(s["seq"]), id=s["seq_id"][:50], description="") for s in uniq],
                            fasta_path, "fasta",
                        )
                        with st.spinner(f"MAFFT ({len(uniq)} séq.)…"):
                            _ok2, _err2 = _msa_run_mafft(fasta_path, aln_path)
                        if not _ok2:
                            st.error(f"Erreur MAFFT : {_err2}")

                    if aln_path.exists():
                        try:
                            _aln2 = AlignIO.read(str(aln_path), "fasta")
                        except Exception as _e2:
                            st.error(f"Impossible de lire l'alignement : {_e2}")
                        else:
                            _ns2 = len(_aln2); _al2 = _aln2.get_alignment_length()
                            st.success(f"**{_ns2} séq. × {_al2} col.**")
                            _html2 = _msa_render_full(_aln2, core_c, var_c, 9999)
                            _height2 = min(_ns2 * 18 + 80, 6000)
                            st.components.v1.html(_html2, height=_height2, scrolling=True)

                st.divider()
                st.markdown("##### 🔴 Actine — positions canoniques (S1)")
                st.caption("Même ordre de lignes que S2 (tri par titre). Scroll horizontal.")
                if not _abp_rows_c2:
                    st.info("Aucune donnée d'interface actine pour ce cluster.")
                else:
                    st.caption(f"**{len(_abp_rows_c2)} séq. ABP** · {len(_aa_at_c2)} positions canoniques")
                    _html_a2   = _msa_render_actin_contacts(_abp_rows_c2, _aa_at_c2, _col_col_c2, 9999)
                    _height_a2 = min(len(_abp_rows_c2) * 18 + 80, 6000)
                    st.components.v1.html(_html_a2, height=_height_a2, scrolling=True)


def _msa_section_s1_clusters():
    """
    Clusters S1 binding site — lignes = séquences ABP, colonnes = positions canoniques actine.
    """
    def _is_actin(t):
        t = t.lower()
        return t.startswith("actin") and "actin-related" not in t and "actin-depolymerizing" not in t

    with st.expander("**Clusters S1 — site de liaison actine**"):
        st.caption(
            "Chaque cluster regroupe les interactions partageant le même site de liaison sur l'actine. "
            "Chaque ligne = une séquence ABP. "
            "Les colonnes sont les positions canoniques actine contactées par au moins un ABP du cluster."
        )

        filt_path = _Path("data/filtered/filtered_all_data.csv")
        int1_path = _Path("data/filtered/details/1.interactions.csv")
        int3_path = _Path("data/filtered/details/3.interface_residues.csv")
        if not filt_path.exists() or not int1_path.exists() or not int3_path.exists():
            st.warning("Données manquantes.")
            return

        df_filt   = pd.read_csv(filt_path, low_memory=False)
        _df1_s1   = pd.read_csv(int1_path)
        _df3_s1   = pd.read_csv(int3_path)
        _df3_s1["residue_number_canon_mafft"] = pd.to_numeric(
            _df3_s1["residue_number_canon_mafft"], errors="coerce")
        _df3_s1["chain_lower"] = _df3_s1["chain"].str.lower()

        df_h = df_filt[~df_filt["subunit_2_title"].apply(lambda t: _is_actin(str(t)))]

        clusters_s1 = (
            df_h[["s1_binding_site_cluster_data_70", "subunit_2_title"]]
            .dropna(subset=["s1_binding_site_cluster_data_70"]).drop_duplicates()
            .groupby("s1_binding_site_cluster_data_70")["subunit_2_title"]
            .apply(lambda x: ", ".join(sorted(set(x))[:4])).reset_index()
        )
        clusters_s1.columns = ["cluster_id", "partners"]
        clusters_s1 = clusters_s1.sort_values("cluster_id")

        from collections import defaultdict as _ddc

        for _, crow in clusters_s1.iterrows():
            cid  = str(crow["cluster_id"])
            ckey = f"s1c_{cid}"

            df_pairs_s1 = df_h[df_h["s1_binding_site_cluster_data_70"] == cid][
                ["subunit_1", "subunit_2", "subunit_2_title", "s2_sequence", "s2_taxonomy_id"]
            ].copy()

            abp_rows_s1, aa_at_s1, col_col_s1 = _msa_actin_contacts_from_pairs(
                df_pairs_s1, _df1_s1, _df3_s1
            )
            if not abp_rows_s1:
                continue

            # ── Séquences ABP (S2) + positions d'interface pour MAFFT ──────────
            rows_s2 = df_pairs_s1[
                ["subunit_2", "subunit_2_title", "s2_sequence", "s2_taxonomy_id"]
            ].drop_duplicates("subunit_2")
            ch2full   = {r["subunit_2"]: str(r["s2_sequence"]).strip()  for _, r in rows_s2.iterrows()}
            ch2seqlow = {ch: s.lower() for ch, s in ch2full.items()}
            ch2title  = {r["subunit_2"]: str(r["subunit_2_title"])       for _, r in rows_s2.iterrows()}
            ch2taxid  = {r["subunit_2"]: r["s2_taxonomy_id"]             for _, r in rows_s2.iterrows()}

            res_s2 = _df3_s1[_df3_s1["chain"].isin(ch2full)][
                ["chain", "residue_number_sequence"]
            ].dropna()
            seq_iface: dict = _ddc(dict)
            for _, rrow in res_s2.iterrows():
                sl = ch2seqlow.get(rrow["chain"])
                if sl:
                    seq_iface[sl].setdefault(rrow["chain"], set()).add(
                        int(rrow["residue_number_sequence"])
                    )

            # Sous-séquence d'interface : seulement les résidus qui contactent l'actine
            seen2: set = set(); uniq_s2: list = []
            core_iface: dict = {}  # iface_seq.lower() → toutes positions (pour _msa_render_full)
            for ch, sl in ch2seqlow.items():
                if sl not in seen2 and sl in seq_iface:
                    seen2.add(sl)
                    taxid = ch2taxid[ch]
                    org   = _MSA_TAX.get(int(taxid), f"taxid:{int(taxid)}") if pd.notna(taxid) else "unk"
                    title = ch2title[ch]
                    bid   = (title[:20].replace(" ", "_").replace(",", "") + "_" + org)[:40]
                    # Union de toutes les positions d'interface pour cette séquence
                    union_pos: set = set()
                    for pos_set in seq_iface[sl].values():
                        union_pos |= pos_set
                    full_seq = ch2full[ch]
                    iface_seq = "".join(
                        full_seq[p - 1] if 0 < p <= len(full_seq) else "X"
                        for p in sorted(union_pos)
                    )
                    uniq_s2.append({"seq_id": bid, "seq": iface_seq, "seqlow": sl,
                                    "n_iface": len(union_pos)})
                    core_iface[iface_seq.lower()] = set(range(1, len(iface_seq) + 1))

            # Trier S1 selon l'ordre du .aln s'il existe déjà
            _aln_path_s1c = _MSA_ALN_DIR / f"{ckey}_msa.aln"
            _sid2seqlow   = {u["seq_id"][:50]: u["seqlow"] for u in uniq_s2}
            if _aln_path_s1c.exists() and _sid2seqlow:
                try:
                    _aln_ord_s1c = AlignIO.read(str(_aln_path_s1c), "fasta")
                    _rank_s1c    = {
                        _sid2seqlow[rec.id]: i
                        for i, rec in enumerate(_aln_ord_s1c)
                        if rec.id in _sid2seqlow
                    }
                    abp_rows_s1 = sorted(abp_rows_s1,
                                         key=lambda r: _rank_s1c.get(r.get("s2seq", ""), 9999))
                except Exception:
                    pass

            partners = crow["partners"][:70]
            with st.expander(f"**{cid}** — {partners} ({len(abp_rows_s1)} séq. ABP)"):

                # ── S2 : MAFFT interface ABP ───────────────────────────────────
                st.markdown("##### 🔵 ABP — interface avec l'actine (S2)")
                if len(uniq_s2) < 2:
                    st.info("Moins de 2 séquences ABP — alignement impossible.")
                else:
                    fasta_path_s1c = _MSA_ALN_DIR / f"{ckey}_msa.fasta"
                    _bc1, _fc1 = st.columns([1, 2])
                    _run_s1c   = _bc1.button("Lancer MAFFT", key=f"msa_{ckey}_btn")
                    _force_s1c = _fc1.checkbox("Forcer recalcul", key=f"msa_{ckey}_force")

                    if _run_s1c or _aln_path_s1c.exists():
                        if _run_s1c or not _aln_path_s1c.exists() or _force_s1c:
                            _MSA_ALN_DIR.mkdir(parents=True, exist_ok=True)
                            SeqIO.write(
                                [SeqRecord(Seq(u["seq"]), id=u["seq_id"][:50], description="")
                                 for u in uniq_s2],
                                fasta_path_s1c, "fasta",
                            )
                            with st.spinner(f"MAFFT ({len(uniq_s2)} séq.)…"):
                                _ok_s1c, _err_s1c = _msa_run_mafft(fasta_path_s1c, _aln_path_s1c)
                            if not _ok_s1c:
                                st.error(f"Erreur MAFFT : {_err_s1c}")

                        if _aln_path_s1c.exists():
                            try:
                                _aln_s1c = AlignIO.read(str(_aln_path_s1c), "fasta")
                            except Exception as _es1c:
                                st.error(f"Impossible de lire l'alignement : {_es1c}")
                            else:
                                _ns1c = len(_aln_s1c); _als1c = _aln_s1c.get_alignment_length()
                                st.success(f"**{_ns1c} séq. × {_als1c} col.**")
                                # Re-trier S1 selon l'ordre réel du nouvel alignement
                                _rank_fresh = {
                                    _sid2seqlow[rec.id]: i
                                    for i, rec in enumerate(_aln_s1c)
                                    if rec.id in _sid2seqlow
                                }
                                if _rank_fresh:
                                    abp_rows_s1 = sorted(
                                        abp_rows_s1,
                                        key=lambda r: _rank_fresh.get(r.get("s2seq", ""), 9999),
                                    )
                                n_iface_total = sum(u["n_iface"] for u in uniq_s2)
                                st.caption(
                                    f"Sous-séquences d'interface — "
                                    f"moy. {n_iface_total // max(len(uniq_s2), 1)} résidus/séq."
                                )
                                _html_s2c = _msa_render_full(_aln_s1c, core_iface, {}, 9999)
                                _h_s2c    = min(_ns1c * 18 + 80, 6000)
                                st.components.v1.html(_html_s2c, height=_h_s2c, scrolling=True)

                # ── S1 : positions canoniques actine ──────────────────────────
                st.divider()
                st.markdown("##### 🔴 Actine — positions canoniques (S1)")
                st.caption(
                    f"**{len(abp_rows_s1)} séquences ABP** — "
                    f"{len(aa_at_s1)} positions canoniques actine"
                )
                _html_s1   = _msa_render_actin_contacts(abp_rows_s1, aa_at_s1, col_col_s1, 9999)
                _height_s1 = min(len(abp_rows_s1) * 18 + 80, 6000)
                st.components.v1.html(_html_s1, height=_height_s1, scrolling=True)

                # ── Analyse des contacts (B–E déplacés dans la section Binding Site) ──
                st.divider()
                st.markdown("##### 📊 Heatmap ABP — contacts avec l'actine")
                _msa_contact_analysis(
                    None, ckey,
                    _ch2seq=ch2seqlow,
                    _ch2title=ch2title,
                    tabs=["A"],
                )

