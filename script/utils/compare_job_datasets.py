import argparse
from pathlib import Path

import pandas as pd


def load_table(path, sep):
    return pd.read_csv(path, sep=sep, dtype=str).fillna("")


def normalize_dataframe(df):
    normalized = df.copy()
    for column in normalized.columns:
        normalized[column] = normalized[column].astype(str).str.strip()
    return normalized


def row_signatures(df, columns):
    if not columns:
        return set()

    records = df[columns].astype(str).to_dict("records")
    return {
        tuple((column, row[column]) for column in columns)
        for row in records
    }


def print_list(title, values):
    print(f"\n{title} ({len(values)})")
    if not values:
        print("none")
        return

    print(list(values))


def compare_by_key(dfs, labels, key_column, shared_columns, max_diffs):
    print(f"\nComparison by key column: {key_column}")

    key_sets = [
        set(df[key_column].astype(str).str.strip()) - {""}
        for df in dfs
    ]

    common_keys = set.intersection(*key_sets)
    print(f"Common keys in all 3 files: {len(common_keys)}")

    for label, key_set in zip(labels, key_sets):
        print(f"{label}: {len(key_set)} unique keys")

    print_list(f"Only in {labels[0]}", sorted(
        key_sets[0] - key_sets[1] - key_sets[2]))
    print_list(f"Only in {labels[1]}", sorted(
        key_sets[1] - key_sets[0] - key_sets[2]))
    print_list(f"Only in {labels[2]}", sorted(
        key_sets[2] - key_sets[0] - key_sets[1]))

    indexed = []
    for df in dfs:
        dedup = df.drop_duplicates(subset=[key_column], keep="first")
        indexed.append(dedup.set_index(key_column))

    differing_keys = []
    comparable_columns = [col for col in shared_columns if col != key_column]

    for key in sorted(common_keys):
        row_values = []
        for idx in indexed:
            row = idx.loc[key, comparable_columns]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            row_values.append(row.to_dict())

        changed_columns = []
        for column in comparable_columns:
            values = [row[column] for row in row_values]
            if len(set(values)) > 1:
                changed_columns.append((column, values))

        if changed_columns:
            differing_keys.append((key, changed_columns))

    print(
        f"\nKeys with at least one differing shared column: {len(differing_keys)}")

    if not differing_keys:
        return

    print(
        f"Showing first {min(max_diffs, len(differing_keys))} differing keys:")
    for key, changed_columns in differing_keys[:max_diffs]:
        print(f"\n[{key}]")
        for column, values in changed_columns:
            print(f"  {column}")
            print(f"    {labels[0]}: {values[0]}")
            print(f"    {labels[1]}: {values[1]}")
            print(f"    {labels[2]}: {values[2]}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare 3 datasets from different jobs and report what changes."
    )
    parser.add_argument("file1", help="First dataset file")
    parser.add_argument("file2", help="Second dataset file")
    parser.add_argument("file3", help="Third dataset file")
    parser.add_argument(
        "--sep",
        default=",",
        help="Field separator used by all files (default: tab)"
    )
    parser.add_argument(
        "--key-column",
        help="Optional key column for record-by-record comparison"
    )
    parser.add_argument(
        "--max-diffs",
        type=int,
        default=20,
        help="Maximum number of differing keys to print in detail"
    )

    args = parser.parse_args()

    paths = [Path(args.file1), Path(args.file2), Path(args.file3)]
    labels = [path.name for path in paths]

    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")

    dfs = [normalize_dataframe(load_table(path, args.sep)) for path in paths]

    print("Files compared:")
    for path, df in zip(paths, dfs):
        print(f"- {path}: {len(df)} rows, {len(df.columns)} columns")

    column_sets = [set(df.columns) for df in dfs]
    shared_columns = sorted(set.intersection(*column_sets))

    print_list("Shared columns", shared_columns)
    print_list(f"Columns only in {labels[0]}", sorted(
        column_sets[0] - column_sets[1] - column_sets[2]))
    print_list(f"Columns only in {labels[1]}", sorted(
        column_sets[1] - column_sets[0] - column_sets[2]))
    print_list(f"Columns only in {labels[2]}", sorted(
        column_sets[2] - column_sets[0] - column_sets[1]))

    shared_row_sets = [row_signatures(df, shared_columns) for df in dfs]
    common_rows = set.intersection(
        *shared_row_sets) if shared_columns else set()

    print(f"\nExact row comparison on shared columns:")
    print(f"- common rows in all 3 files: {len(common_rows)}")
    print(
        f"- rows only in {labels[0]}: {len(shared_row_sets[0] - shared_row_sets[1] - shared_row_sets[2])}")
    print(
        f"- rows only in {labels[1]}: {len(shared_row_sets[1] - shared_row_sets[0] - shared_row_sets[2])}")
    print(
        f"- rows only in {labels[2]}: {len(shared_row_sets[2] - shared_row_sets[0] - shared_row_sets[1])}")

    if args.key_column:
        if any(args.key_column not in df.columns for df in dfs):
            raise ValueError(
                f"Key column '{args.key_column}' must exist in all 3 files"
            )
        compare_by_key(
            dfs,
            labels,
            args.key_column,
            shared_columns,
            args.max_diffs
        )


if __name__ == "__main__":
    main()
