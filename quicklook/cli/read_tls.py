#!/usr/bin/env python

import flammkuchen as fk
from glob import glob
import pandas as pd
from tqdm import tqdm
import argparse


def main():
    parser = argparse.ArgumentParser(description="Summarize TLS results.")
    parser.add_argument("input_dir", type=str, help="Input directory containing *.h5 files")
    parser.add_argument(
        "-param",
        choices=["TOI", "power_gls", "Prot_gls", "amp_gls", "SDE_tls", "Porb_tls"],
        default="SDE_tls",
        type=str,
        help="Parameter to sort by",
    )
    args = parser.parse_args()

    input_dir = args.input_dir
    param = args.param

    files = glob(input_dir + "/*.h5")
    assert len(files) > 0, "No *.h5 files found!"

    ss = []
    errors = []

    print("Reading *.tls files...")
    for file in tqdm(files):
        try:
            data = fk.load(file)
            d = {
                "TOI": data.get("toiid"),
                "TIC": data.get("ticid"),
                "gaiaid": data.get("gaiaid"),
                "power_gls": data.get("power_gls")[0] if data.get("power_gls") else None,
                "SDE_tls": data.get("SDE"),
                "Porb_tls": data.get("period"),
                "Prot_gls": data.get("Prot_gls")[0] if data.get("Prot_gls") else None,
                "amp_gls": data.get("amp_gls")[0] if data.get("amp_gls") else None,
                "depth": (1 - data.get("depth")) * 1e3 if data.get("depth") is not None else None,
                "simbad_object": data.get("simbad_obj"),
                "filename": file,
            }
            ss.append(pd.Series(d))
        except Exception as e:
            print(e)
            errors.append(data.get("ticid") if data else file)

    print(f"Sorting by {param}...")
    df = (
        pd.DataFrame(ss)
        .sort_values(by=param, ascending=False)
        .reset_index(drop=True)
        .drop_duplicates()
    )

    fp = input_dir + "_tls.csv"
    df.to_csv(fp, index=False)
    print(f"Saved: {fp}")
    if errors:
        print("Errors:", errors)


if __name__ == "__main__":
    main()
