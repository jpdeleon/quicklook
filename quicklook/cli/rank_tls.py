#!/usr/bin/env python
import os
import shutil
import argparse
import pandas as pd
from pathlib import Path
import numpy as np


def is_near_integer(series, tolerance):
    """
    Helper function to check if a ratio is close to an integer
    """
    fractional_part = series % 1
    return np.isclose(fractional_part, 0, atol=tolerance) | np.isclose(
        fractional_part, 1, atol=tolerance
    )


def apply_filters(df, tolerance=0.2, depth_min=1, depth_max=100, min_SDE=5):
    """
    tolerance in days
    depth_min, depth_max in ppt
    """
    df2 = df.copy()
    # Condition 1: Absolute difference is not near an integer
    # not_near_period = abs(df['Prot_gls'] - df['Porb_tls']) >= period_tolerance
    depth_range = (df2["depth"] >= depth_min) & (df2["depth"] < depth_max)
    errmsg = f"No candidates satisfy `{depth_min}<depth_range<{depth_max}` ppt."
    assert sum(depth_range) > 0, errmsg

    # Condition 2: Ratio is not close to an integer
    ratio_forward_is_near_int = is_near_integer(df2["Prot_gls"] / df2["Porb_tls"], tolerance)
    ratio_inverse_is_near_int = is_near_integer(df2["Porb_tls"] / df2["Prot_gls"], tolerance)

    # Filter for rows where NEITHER ratio is near an integer.
    # The ~ operator inverts the boolean mask.
    not_near_period_ratio = ~(ratio_forward_is_near_int | ratio_inverse_is_near_int)
    errmsg = "No candidates satisfy `not_near_period_ratio`."
    assert sum(not_near_period_ratio) > 0, errmsg

    # Condition 3: Not EB!
    not_EB = ~df2["simbad_object"].str.lower().isin(["eclipsing binary", "eclbin"])
    errmsg = "No candidates satisfy `not_EB`."
    assert sum(not_EB) > 0, errmsg

    # Condition 4: SDE>5
    strong_signal = df["SDE_tls"] > 5
    errmsg = f"No candidates satisfy `SDE>{min_SDE}`."
    assert sum(strong_signal) > 0, errmsg

    # Condition 5: not a known TOI
    not_TOI = df2["TOI"].isna() | (df2["TOI"] == "")
    idx = depth_range & not_near_period_ratio & not_EB & strong_signal & not_TOI
    return df2[idx]


def rename_and_copy(src_dir, csv_path, dst_dir, column_name=None):
    """ """
    # Create output directory if it doesn't exist
    Path(dst_dir).mkdir(exist_ok=True)

    # Load CSV with pandas; already sorted in decreasing SDE
    df = pd.read_csv(csv_path)

    # Determine the column to use
    if column_name is not None:
        if column_name not in df.columns:
            raise ValueError(f"Column '{column_name}' not found in CSV.")
        df = df.sort_values(by=column_name, ascending=False)

    # import pdb; pdb.set_trace()
    df2 = apply_filters(df)

    filenames = df2.filename.apply(lambda x: x.split("/")[-1].split("_tls.")[0] + ".png").values

    # Determine padding width
    num_digits = len(str(len(filenames)))

    # Copy and rename files
    for idx, filename in enumerate(filenames, start=1):
        src_path = Path(src_dir, filename)
        if not src_path.exists():
            print(f"Warning: File not found: {src_path}")
            continue
        prefix = str(idx).zfill(num_digits)
        new_name = f"{prefix}_{filename}"
        dst_path = Path(dst_dir, new_name)

        shutil.copy2(src_path, dst_path)
        print(f"Copied: {src_path} -> {dst_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Copy and prefix-rank files based on CSV order using pandas."
    )
    parser.add_argument("input_dir", help="Directory where the original files are located.")
    parser.add_argument(
        "--csv_path", help="Path to the CSV file containing the ranked filenames.", default=None
    )
    parser.add_argument(
        "--output_dir", help="Directory where renamed files will be copied.", default=None
    )
    parser.add_argument(
        "--column", help="Column name in csv to rank the files in order.", default=None
    )
    parser.add_argument(
        "--ascending", help="Ranking order. Default is False (descending).", default=False
    )

    args = parser.parse_args()
    if args.csv_path is None:
        csv_path = f"{args.input_dir}_tls.csv"
        if not os.path.exists(csv_path):
            os.system(f"read_tls {args.input_dir}")
    if args.output_dir is None:
        output_dir = f"{args.input_dir}/temp"
    rename_and_copy(args.input_dir, csv_path, output_dir, column_name=args.column)
