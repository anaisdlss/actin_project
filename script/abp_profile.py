"""Audited residue profiles from existing ABP ProteoCast score files.

The submitted query, when present, is authoritative. Without that file the
reference sequence can be reconstructed from a complete, internally consistent
20-substitution score grid; input PDB-chain FASTAs are not query substitutes.
"""
from pathlib import Path

import numpy as np
import pandas as pd


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def query_file(csv_path):
    """Locate the query for either supported nested or flat result layout."""
    csv_path = Path(csv_path)
    folder = csv_path.parent if csv_path.name == "4.query_ProteoCast.csv" else csv_path.parent / csv_path.stem
    return folder / "1.query.fasta"


def read_query(path):
    """Read one nonempty FASTA record; reject ambiguous/multiple queries."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return None
    lines = path.read_text().splitlines()
    if sum(line.startswith(">") for line in lines) != 1:
        raise ValueError("The ProteoCast query must contain exactly one FASTA record.")
    sequence = "".join(line.strip() for line in lines if not line.startswith(">")).upper()
    if not sequence or set(sequence) - set(AMINO_ACIDS):
        raise ValueError("The ProteoCast query has missing or noncanonical amino acids.")
    return sequence


def load_abp_scores(csv_path, query_path=None):
    """Return raw scores and a full-length profile without filling missing scores.

    A mean requires all 20 distinct canonical substitutions with finite scores,
    including the unchanged amino acid. Incomplete positions remain NaN.
    Ambiguous mutation identifiers, duplicate substitutions, inconsistent
    reference letters and query contradictions are rejected, never averaged.
    """
    csv_path = Path(csv_path)
    raw = pd.read_csv(csv_path)
    required = {"Mutation", "Residue", "Variant_score"}
    if raw.empty or not required.issubset(raw.columns):
        raise ValueError("ProteoCast scores require Mutation, Residue and Variant_score.")
    mutation = raw.Mutation.astype(str).str.extract(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")
    if mutation.isna().any().any():
        raise ValueError("Invalid ProteoCast mutation identifier.")
    raw = raw.copy()
    raw["position"] = mutation[1].astype(int)
    raw["reference_aa"] = mutation[0]
    raw["alternate_aa"] = mutation[2]
    stated_positions = pd.to_numeric(raw.Residue, errors="coerce")
    if not stated_positions.eq(raw.position).all() or raw.position.lt(1).any():
        raise ValueError("Residue positions contradict the ProteoCast mutation identifiers.")
    if raw.duplicated(["position", "alternate_aa"]).any():
        raise ValueError("Duplicate substitution scores in ProteoCast results.")
    references = raw.groupby("position").reference_aa
    if references.nunique().ne(1).any():
        raise ValueError("Contradictory reference amino acids within a ProteoCast position.")
    query_path = Path(query_path) if query_path else query_file(csv_path)
    query = read_query(query_path)
    refs = references.first()
    length = len(query) if query else int(raw.position.max())
    if query and (raw.position.gt(length).any() or any(query[int(p) - 1] != aa for p, aa in refs.items())):
        raise ValueError("ProteoCast scores contradict the submitted query sequence.")
    raw["Variant_score"] = pd.to_numeric(raw.Variant_score, errors="coerce")
    raw.loc[~np.isfinite(raw.Variant_score), "Variant_score"] = np.nan
    positions = pd.RangeIndex(1, length + 1, name="position")
    grid = raw.pivot(index="position", columns="alternate_aa", values="Variant_score")
    grid = grid.reindex(index=positions, columns=list(AMINO_ACIDS))
    profile = pd.DataFrame(index=positions)
    profile["reference_aa"] = list(query) if query else refs.reindex(positions)
    profile["finite_score_count"] = grid.notna().sum(axis=1)
    profile["complete_score_grid"] = profile.finite_score_count.eq(20)
    profile["mean_variant_score"] = grid.mean(axis=1).where(profile.complete_score_grid)
    profile["mutational_sensitivity"] = -profile.mean_variant_score
    profile["score_status"] = np.where(profile.complete_score_grid, "complete (20 finite scores)",
                                        "incomplete; sensitivity unavailable")
    provenance = ("submitted ProteoCast query FASTA" if query else
                  "reference reconstructed from ProteoCast mutation labels")
    profile["sequence_source"] = provenance
    return raw, profile.reset_index()


def score_sequence(csv_path, query_path=None):
    """Safe sequence fallback for contact mapping; no incomplete grid accepted."""
    _, profile = load_abp_scores(csv_path, query_path=query_path)
    if not profile.complete_score_grid.all() or profile.reference_aa.isna().any():
        raise ValueError("A complete score grid is required to reconstruct the ProteoCast query.")
    return "".join(profile.reference_aa)


def annotate_profile(profile, iface_asa=None, domains=None, rsa=None, surface_only=False):
    """Attach observed measurements; unavailable ASA/RSA remain NaN, not zero."""
    frame = profile.copy()
    for column, values in (("buried_ASA_percent_max", iface_asa), ("rsa", rsa)):
        numeric = pd.to_numeric(frame.position.map(values or {}), errors="coerce")
        frame[column] = numeric.where(np.isfinite(numeric))
    domain_names = []
    for position in frame.position:
        names = [f"{d['name']} ({d['db']})" for d in (domains or [])
                 if any(start <= position <= end for start, end in d.get("spans", []))]
        domain_names.append("; ".join(dict.fromkeys(names)))
    frame["domains"] = domain_names
    frame["displayed"] = frame.rsa.ge(0.2) if surface_only else True
    return frame


def profile_figure(profile, domains=None, title="", focus=None):
    """Continuous sensitivity and observed ASA tracks with shared position hover."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from plot_interaction import position_hover

    domains = list(domains or [])
    count = len(domains)
    fig = make_subplots(rows=3 if count else 2, cols=1, shared_xaxes=True,
                        vertical_spacing=0.06,
                        row_heights=[0.52, 0.28, 0.20] if count else [0.65, 0.35],
                        subplot_titles=["Mutational sensitivity", "Observed buried ASA at actin contacts"]
                        + (["Annotated domains"] if count else []))
    shown = profile.displayed
    domain_text = profile.domains.replace("", "No domain annotation at this position")
    custom = np.column_stack([profile.reference_aa.fillna("?"), domain_text, profile.score_status])
    fig.add_trace(go.Scatter(
        x=profile.position, y=profile.mutational_sensitivity.where(shown), mode="lines",
        line=dict(color="#0072B2", width=1.7), connectgaps=False,
        name="Mutational sensitivity", customdata=custom,
        hovertemplate="%{customdata[0]}%{x}<br>Sensitivity: %{y:.3f}"
                      "<br>%{customdata[1]}<extra></extra>"), row=1, col=1)
    fig.add_trace(go.Scatter(
        x=profile.position, y=profile.buried_ASA_percent_max.where(shown), mode="lines+markers",
        line=dict(color="#E69F00", width=1.5), marker=dict(size=3), connectgaps=False,
        name="Observed contact ASA", hovertemplate="Buried ASA: %{y:.2f}%<extra></extra>"), row=2, col=1)
    palette = ["#0072B2", "#E69F00", "#CC79A7", "#56B4E9", "#6A51A3", "#777777"]
    for i, domain in enumerate(domains):
        label = f"{domain['name']} ({domain['db']})"
        covered = [any(start <= p <= end for start, end in domain.get("spans", []))
                   for p in profile.position]
        # Dense domain positions provide an accurate cursor inside spans, too.
        fig.add_trace(go.Scatter(x=profile.position,
                                y=[label if yes else None for yes in covered],
                                mode="lines", line=dict(color=palette[i % len(palette)], width=8),
                                connectgaps=False, showlegend=False, hoverinfo="skip"), row=3, col=1)
    fig.update_layout(height=560 + min(count * 22, 440), showlegend=False,
                      margin=dict(l=8, r=8, t=65, b=40), title=title,
                      plot_bgcolor="white")
    fig.update_yaxes(title_text="−mean variant score", row=1, col=1)
    fig.update_yaxes(title_text="Buried ASA (%)", range=[0, 100], row=2, col=1)
    if count:
        fig.update_yaxes(autorange="reversed", tickfont=dict(size=9), row=3, col=1)
    fig.update_xaxes(title_text="ABP residue (ProteoCast query)", row=3 if count else 2, col=1)
    lower, upper = focus if focus else (1, int(profile.position.max()))
    fig.update_xaxes(range=[lower - 0.5, upper + 0.5])
    return position_hover(fig)


def nonempty_alignment(directory):
    """Offer only actual, nonempty multi-sequence ProteoCast FASTA alignments."""
    for path in sorted(Path(directory).glob("2.*.fasta")):
        if path.stat().st_size:
            text = path.read_text()
            records = ["".join(record.splitlines()[1:]).strip() for record in text.split(">")[1:]]
            if len(records) >= 2 and all(records):
                return path
    return None
