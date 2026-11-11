"""Merge Thermometer scores into the feature dataframe produced by
FeatureCalculatorPcAfGraphs.py.

    python combine_thermo_features.py -D <job_folder> [--thermo-dir <dir>]

Reads  <job_folder>/cmb_df.pkl
       <thermo_dir>/RESULTS/*_Ts_global_score.csv    (default: <job_folder>/thermo)
Writes <job_folder>/cmb_df.pkl                       (9 Ts_* columns added)
       <job_folder>/thermo_Ts_global_score.csv       (the concatenated table)

Key alignment
-------------
cmb_df is indexed by the PDB filename stem (e.g. "AF-O00168-F1-model_v6") --
all three of its sources use os.path.basename(f).split('.')[0]. The PDB column
of *_Ts_global_score.csv now carries the same stem, because Ts_global_computer.R
takes it from the name the caller passes instead of strsplit(filename,"_")[1]
(which used to cut it down to "AF-O00168-F1-model").

Structures that failed in Thermometer simply get NaN in the Ts_* columns rather
than being dropped, so the row count of cmb_df never changes.

This replaces the never-found combine_feature_dfs.py referenced by
calculate_all_features.sh.
"""

import argparse
import os
import sys
from glob import glob

import joblib
import pandas as pd

TS_COLUMNS = [
    "Ts_tot_PN", "Ts_tot_N", "Ts_tot_P",
    "Ts_C_PN", "Ts_C_N", "Ts_C_P",
    "Ts_lj_PN", "Ts_lj_N", "Ts_lj_P",
]


def load_ts_table(results_dir):
    """Concatenate every *_Ts_global_score.csv into one frame indexed by PDB."""
    files = sorted(glob(os.path.join(results_dir, "*_Ts_global_score.csv")))
    if not files:
        return None, files

    frames = []
    for f in files:
        df = pd.read_csv(f, index_col=0)
        missing = [c for c in ["PDB"] + TS_COLUMNS if c not in df.columns]
        if missing:
            print(f"  skipping {os.path.basename(f)}: missing {missing}",
                  file=sys.stderr)
            continue
        frames.append(df[["PDB"] + TS_COLUMNS])

    if not frames:
        return None, files

    ts = pd.concat(frames, ignore_index=True)
    ts = ts.drop_duplicates(subset="PDB", keep="last").set_index("PDB")
    # R writes these quoted, so they arrive as object dtype
    ts[TS_COLUMNS] = ts[TS_COLUMNS].apply(pd.to_numeric, errors="coerce")
    return ts, files


def main():
    parser = argparse.ArgumentParser(
        description="Merge Thermometer Ts_* scores into cmb_df.pkl")
    parser.add_argument("-D", "--dir", required=True, help="job folder")
    parser.add_argument("--thermo-dir", default=None,
                        help="dir containing RESULTS/ (default: <job folder>/thermo)")
    parser.add_argument("--require-all", action="store_true",
                        help="exit non-zero if any structure has no Ts scores")
    args = parser.parse_args()

    folder = os.path.abspath(args.dir)
    thermo_dir = os.path.abspath(args.thermo_dir or os.path.join(folder, "thermo"))
    results_dir = os.path.join(thermo_dir, "RESULTS")

    cmb_path = os.path.join(folder, "cmb_df.pkl")
    if not os.path.isfile(cmb_path):
        print(f"no cmb_df.pkl in {folder}", file=sys.stderr)
        return 1

    cmb = joblib.load(cmb_path)
    ts, files = load_ts_table(results_dir)

    if ts is None:
        print(f"no usable *_Ts_global_score.csv in {results_dir}", file=sys.stderr)
        return 1

    ts.to_csv(os.path.join(folder, "thermo_Ts_global_score.csv"))

    # Drop any stale Ts_* columns so re-running is idempotent.
    cmb = cmb.drop(columns=[c for c in TS_COLUMNS if c in cmb.columns])
    merged = cmb.join(ts[TS_COLUMNS], how="left")

    if len(merged) != len(cmb):
        print(f"join changed the row count ({len(cmb)} -> {len(merged)})",
              file=sys.stderr)
        return 1

    joblib.dump(merged, cmb_path)

    matched = int(merged[TS_COLUMNS[0]].notna().sum())
    print(f"merged {len(files)} Ts files; {matched}/{len(merged)} structures "
          f"have thermostability scores -> {cmb_path}")

    if matched < len(merged):
        unmatched = list(merged.index[merged[TS_COLUMNS[0]].isna()])
        print(f"  no scores for: {', '.join(map(str, unmatched))}", file=sys.stderr)
        if args.require_all:
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
