import pandas as pd
from pathlib import Path


def load_data(data_dir: str | Path = "data") -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Charge et fusionne les données filtrées de façon cohérente avec le pipeline.

    Returns
    -------
    df_all     : interactions actine (filtered_all_data.csv)
    df_summary : summary filtré avec interaction_id séquentiel (= position dans le fichier,
                 cohérent avec les IDs des tables de détails téléchargées par get_interaction_details)
    df_merged  : df_all LEFT JOIN df_summary sur (subunit_1/2 ↔ Protein ID/Interactor ID)
    """
    data_dir = Path(data_dir)

    df_all = pd.read_csv(data_dir / "filtered/filtered_all_data.csv")

    df_summary = pd.read_csv(data_dir / "filtered/filtered_summary.csv")
    df_summary["interaction_id"] = range(1, len(df_summary) + 1)

    df_merged = df_all.merge(
        df_summary,
        left_on=["subunit_1", "subunit_2"],
        right_on=["Protein ID", "Interactor ID"],
        how="left",
    )

    return df_all, df_summary, df_merged
