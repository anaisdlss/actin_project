import argparse
from pathlib import Path

import pandas as pd


def load_unique_values(path, column, sep):
    df = pd.read_csv(path, sep=sep, dtype=str)

    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in {path}")

    values = (
        df[column]
        .dropna()
        .astype(str)
        .str.strip()
    )

    values = values[values != ""]
    return set(values.unique())


def print_set(label, values):
    print(f"\n{label} ({len(values)})")
    print(sorted(values))


def print_values_only(values):
    print(sorted(values))


def main():
    parser = argparse.ArgumentParser(
        description="Compare unique values from one column across 3 files."
    )
    parser.add_argument("file1", help="Path to the first file")
    parser.add_argument("file2", help="Path to the second file")
    parser.add_argument("file3", help="Path to the third file")
    parser.add_argument("column", help="Column name to compare")
    parser.add_argument(
        "--sep",
        default=",",
        help="CSV separator used for all files (default: ,)"
    )
    parser.add_argument(
        "--only",
        choices=[
            "common",
            "only1",
            "only2",
            "only3"
        ],
        help="Print only one result set, one value per line"
    )

    args = parser.parse_args()

    paths = [Path(args.file1), Path(args.file2), Path(args.file3)]

    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")

    uniques = [
        load_unique_values(path, args.column, args.sep)
        for path in paths
    ]

    labels = [str(path) for path in paths]
    a, b, c = uniques

    print(f"Column compared: {args.column}")
    for label, values in zip(labels, uniques):
        print(f"{label}: {len(values)} unique values")

    common_all = a & b & c
    only_a = a - b - c
    only_b = b - a - c
    only_c = c - a - b

    result_sets = {
        "common": common_all,
        "only1": only_a,
        "only2": only_b,
        "only3": only_c
    }

    if args.only:
        print_values_only(result_sets[args.only])
        return

    print_set("Common to all 3 files", common_all)
    print_set(f"Only in {labels[0]}", only_a)
    print_set(f"Only in {labels[1]}", only_b)
    print_set(f"Only in {labels[2]}", only_c)


if __name__ == "__main__":
    main()
