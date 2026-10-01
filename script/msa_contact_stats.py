"""Exact contact scope and conditional chain/PDB means for legacy MSA panels."""
import numpy as np
import pandas as pd


def normalize_chain(value):
    parts = str(value).split("_", 1)
    return parts[0].lower() + "_" + parts[1] if len(parts) == 2 else str(value).lower()


def oriented_pairs(all_data, interactions):
    """Match retained chain pairs to interaction IDs, then orient each actin side."""
    metadata = all_data.copy()
    ids = interactions[["interaction_id", "chain_A_id", "chain_B_id"]].copy()
    for col in ("subunit_1", "subunit_2"):
        metadata[col] = metadata[col].map(normalize_chain)
    for col in ("chain_A_id", "chain_B_id"):
        ids[col] = ids[col].map(normalize_chain)
    ids = pd.concat([ids, ids.rename(columns={"chain_A_id": "chain_B_id", "chain_B_id": "chain_A_id"})]).drop_duplicates()
    joined = metadata.merge(ids, left_on=["subunit_1", "subunit_2"],
                            right_on=["chain_A_id", "chain_B_id"], how="inner")
    parts = []
    for side, other in ((1, 2), (2, 1)):
        rows = joined[joined[f"s{side}_actine"].astype(str).str.lower().eq("true")]
        parts.append(pd.DataFrame({
            "interaction_id": rows.interaction_id,
            "actin_chain": rows[f"subunit_{side}"],
            "partner_chain": rows[f"subunit_{other}"],
            "site": rows[f"s{side}_binding_site_cluster_data_70"],
            "seq_low": rows[f"s{other}_sequence"].astype("string").str.strip().str.lower(),
            "title": rows[f"subunit_{other}_title"],
            "partner_is_actin": rows[f"s{other}_actine"].astype(str).str.lower().eq("true"),
        }))
    result = pd.concat(parts, ignore_index=True).drop_duplicates()
    return result.dropna(subset=["seq_low"])


def scoped_contacts(contacts, pairs, site=None, filter_fn=None, chains=None, rigor_pdbs=None):
    """Preserve chain case/insertion codes, and filter on exact oriented triples."""
    selected = pairs.copy()
    if site is not None:
        selected = selected[selected.site.astype(str).eq(str(site))]
    elif chains is not None:
        selected = selected[selected.partner_chain.isin({normalize_chain(c) for c in chains})]
    elif filter_fn is not None:
        selected = selected[selected.title.map(lambda value: filter_fn(str(value)))]
    if rigor_pdbs:
        selected = selected[selected.actin_chain.str.split("_").str[0].isin({str(p).lower() for p in rigor_pdbs})]
    source = contacts.copy()
    for col in ("chain_A_id", "chain_B_id"):
        source[col] = source[col].map(normalize_chain)
    renames = {}
    for col in source:
        if "_A_" in col:
            renames[col] = col.replace("_A_", "_B_")
        elif "_B_" in col:
            renames[col] = col.replace("_B_", "_A_")
        elif col.endswith("_A"):
            renames[col] = col[:-1] + "B"
        elif col.endswith("_B"):
            renames[col] = col[:-1] + "A"
    parts = []
    for oriented in (source, source.rename(columns=renames)):
        parts.append(oriented.merge(selected, left_on=["interaction_id", "chain_A_id", "chain_B_id"],
                                    right_on=["interaction_id", "actin_chain", "partner_chain"], how="inner"))
    result = pd.concat(parts, ignore_index=True).drop_duplicates()
    result["pdb"] = result.actin_chain.str.split("_").str[0]
    # Censored strings such as '<0.1' remain NaN; no artificial contact area.
    numeric = pd.to_numeric(result.contact_area, errors="coerce")
    result["area_f"] = numeric.where(np.isfinite(numeric) & numeric.gt(0))
    for side in ("a", "b"):
        result[f"canon_{side}"] = pd.to_numeric(result[f"residue_{side.upper()}_canon_mafft"], errors="coerce")
        result[f"aa_{side}"] = result[f"residue_{side.upper()}_name"].astype(str).str.strip().str.upper()
    return result


def conditional_area_profile(contacts, side):
    """Return position means and additive AA contributions with equal PDB weights.

    Areas are summed per observed chain/position (for the displayed side),
    averaged across observed chains within each PDB, then across observed PDBs.
    An absent position contributes neither a zero nor a denominator. At an
    observed position an absent AA contributes zero to that position's mixture.
    """
    if side not in ("a", "b"):
        raise ValueError("Side must be 'a' (actin) or 'b' (partner).")
    position, aa, chain = f"canon_{side}", f"aa_{side}", f"chain_{side.upper()}_id"
    keys = ["seq_low", "label", position]
    rows = contacts[np.isfinite(contacts.area_f) & contacts.area_f.gt(0)].dropna(subset=[position]).copy()
    output_columns = keys + ["corr_mean", "pct", "observed_pdbs", "observed_chains"]
    if rows.empty:
        return pd.DataFrame(columns=output_columns), pd.DataFrame(columns=keys + [aa, "corr_mean"])
    per_chain = rows.groupby(keys + ["pdb", chain, aa], dropna=False).area_f.sum().rename("area")
    chain_denominator = rows.groupby(keys + ["pdb"])[chain].nunique().rename("chain_count")
    pdb_aa = per_chain.groupby(level=keys + ["pdb", aa]).sum().reset_index()
    pdb_aa = pdb_aa.merge(chain_denominator.reset_index(), on=keys + ["pdb"])
    pdb_aa["pdb_mean"] = pdb_aa.area / pdb_aa.chain_count
    pdb_denominator = rows.groupby(keys).pdb.nunique().rename("observed_pdbs")
    aa_means = pdb_aa.groupby(keys + [aa]).pdb_mean.sum().reset_index()
    aa_means = aa_means.merge(pdb_denominator.reset_index(), on=keys)
    aa_means["corr_mean"] = aa_means.pdb_mean / aa_means.observed_pdbs
    profile = aa_means.groupby(keys).corr_mean.sum().reset_index()
    profile = profile.merge(pdb_denominator.reset_index(), on=keys)
    chain_counts = rows.groupby(keys)[chain].nunique().rename("observed_chains")
    profile = profile.merge(chain_counts.reset_index(), on=keys)
    totals = profile.groupby("seq_low").corr_mean.transform("sum")
    profile["pct"] = profile.corr_mean / totals * 100
    return profile[output_columns], aa_means[keys + [aa, "corr_mean"]]
