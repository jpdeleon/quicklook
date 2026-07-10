#!/usr/bin/env python
import os
import shutil
import argparse
import subprocess
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


def apply_filters(
    df: pd.DataFrame,
    tolerance: float = 0.2,
    min_depth: float = 1,
    max_depth: float = 100,
    min_SDE: float = 5,
) -> pd.DataFrame:
    """
    tolerance in days
    min_depth, max_depth in ppt
    """
    df2 = df.copy()
    # Condition 1: Absolute difference is not near an integer
    # not_near_period = abs(df['Prot_gls'] - df['Porb_tls']) >= period_tolerance
    depth_range = (df2["depth"] >= min_depth) & (df2["depth"] < max_depth)
    errmsg = f"No candidates satisfy `{min_depth}<depth_range<{max_depth}` ppt."
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
    # not_EB = ~df2["simbad_object"].str.lower().isin(["eclipsing binary", "eclbin"])
    not_EB = ~df2["simbad_object"].fillna("").str.lower().isin(["eclipsing binary", "eclbin"])
    errmsg = "No candidates satisfy `not_EB`."
    assert sum(not_EB) > 0, errmsg

    # Condition 4: SDE>5
    strong_signal = df["SDE_tls"] > 5
    errmsg = f"No candidates satisfy `SDE>{min_SDE}`."
    assert sum(strong_signal) > 0, errmsg

    # Condition 5: not a known TOI
    not_TOI = df2["TOI"].isna() | (df2["TOI"] == "")

    # Condition 6: In transit data is not sparse
    exp = pd.to_numeric(
        df2.get("exptime", pd.Series(9999, index=df2.index)), errors="coerce"
    ).fillna(9999)
    count = df2.get("per_transit_count", pd.Series(0, index=df2.index))
    # per_transit_count is a tuple stored as string in CSV; extract the minimum count
    if count.dtype == object:
        import ast

        def _min_count(val):
            if pd.isna(val) or val == "" or val == "0":
                return 0
            try:
                parsed = ast.literal_eval(str(val))
                if isinstance(parsed, (tuple, list)) and len(parsed) > 0:
                    return min(int(x) for x in parsed)
                return int(parsed)
            except (ValueError, SyntaxError):
                return 0

        count = count.apply(_min_count)

    # Define the two valid "not sparse" conditions
    short_exp_valid = (exp <= 120) & (count > 10)
    long_exp_valid = (exp > 120) & (count > 5)

    # Data is not sparse if it meets either valid condition OR if count is missing (0)
    not_sparse = short_exp_valid | long_exp_valid | (count == 0)

    idx = depth_range & not_near_period_ratio & not_EB & strong_signal & not_TOI & not_sparse
    return df2[idx]


def rename_and_copy(
    src_dir: str,
    csv_path: str,
    dst_dir: str,
    column_name: str = None,
    remove_filter: bool = True,
    use_ascending: bool = False,
    min_SDE: float = 5,
    tolerance: float = 0.2,
    min_depth: float = 1,
    max_depth: float = 100,
):
    """ """
    # Create output directory if it doesn't exist
    Path(dst_dir).mkdir(exist_ok=True)

    # Load CSV with pandas; already sorted in decreasing SDE
    df = pd.read_csv(csv_path)

    # Determine the column to use
    if column_name is not None:
        if column_name not in df.columns:
            raise ValueError(f"Column '{column_name}' not found in CSV.")
        df = df.sort_values(by=column_name, ascending=use_ascending)

    # import pdb; pdb.set_trace()
    if remove_filter:
        print("Ranking all files...")
        df2 = df.copy()
    else:
        print("Applying filters to files...")
        df2 = apply_filters(
            df, tolerance=tolerance, min_depth=min_depth, max_depth=max_depth, min_SDE=min_SDE
        )
        err_msg = "No good candidates left after applying the filters. Try --remove_filter."
        assert len(df2) > 0, err_msg

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


def main():
    parser = argparse.ArgumentParser(
        description="Copy and add prefix to png files ordered by `SDE_tls` (default) defined in csv output of `read_tls` script. Only ranks filtered results."
    )
    parser.add_argument("input_dir", help="Directory where the original files are located.")
    parser.add_argument(
        "--csv_path",
        help="Path to the CSV file containing the ranked filenames.",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--output_dir",
        help="Directory where renamed files will be copied. Default is `temp`.",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--column",
        help="Column name in csv to rank the files in order. Default is `SDE_tls`.",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--remove_filter",
        help="Remove default filters. Default is False (use filters).",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--use_ascending",
        help="Ranking order. Default is False (descending).",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--min_sde", help="Minimum SDE to copy. Default is 5.", default=5, type=float
    )
    parser.add_argument(
        "--min_depth",
        help="Minimum transit depth in ppt to copy. Default is 1 ppt.",
        default=1,
        type=float,
    )
    parser.add_argument(
        "--max_depth",
        help="Maximum transit depth in ppt to copy. Default is 100 ppt.",
        default=100,
        type=float,
    )

    args = parser.parse_args()
    if args.csv_path is None:
        csv_path = f"{args.input_dir}_tls.csv"
        if not os.path.exists(csv_path):
            # No shell: input_dir is a user-supplied path and would otherwise
            # be word-split and glob/metachar-expanded by /bin/sh. check=True
            # so a failed extraction surfaces here rather than as a confusing
            # "file not found" from the read_csv below.
            subprocess.run(["read_tls", args.input_dir], check=True)
    else:
        csv_path = args.csv_path
    if args.output_dir is None:
        output_dir = f"{args.input_dir}/temp"
    else:
        output_dir = f"{args.input_dir}/{args.output_dir}"
    rename_and_copy(
        args.input_dir,
        csv_path,
        output_dir,
        remove_filter=args.remove_filter,
        # use_ascending=args.use_ascending,
        min_SDE=args.min_sde,
        min_depth=args.min_depth,
        max_depth=args.max_depth,
        column_name=args.column,
    )


if __name__ == "__main__":
    main()
