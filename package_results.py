"""Package a finished job's features into the single file users download.

    python package_results.py -D <job_folder>

Reads  <job_folder>/cmb_df.pkl
Writes <job_folder>/phanes_results.zip containing

    features_scalar.csv        one row per structure, the scalar features
    features_per_residue.pkl   one row per structure, one array per residue
    README.txt                 how to load them and how to rejoin them

Why split, and why compress
---------------------------
The per-residue profiles dominate the size: 96 of the 255 columns are arrays
whose length is the structure's residue count, which measured at ~1.2 KB per
residue. At MAX_PDB = 10 structures of 2700 residues that is a ~32 MB pickle,
and the in-app download used to base64 it into a Dash callback response (a
further 1.33x). Most users only want the scalar table, which is a few KB.

So: scalars go out as a plain CSV that opens in anything, profiles stay a
pickle (a CSV of arrays would be unusable and far larger), and the pickle is
written with joblib compress=3 -- measured 1.8x on real output.

The pickle is added to the zip with ZIP_STORED because joblib has already
compressed it; deflating it a second time only costs CPU. The CSV and README
are deflated normally.
"""

import argparse
import io
import os
import sys
import zipfile

import joblib
import numpy as np
import pandas as pd

ARCHIVE_NAME = "phanes_results.zip"
SCALAR_NAME = "features_scalar.csv"
PROFILE_NAME = "features_per_residue.pkl"

README = """PHANES feature calculation results
==================================

features_scalar.csv
    One row per input structure, indexed by the PDB file name (without the
    .pdb extension). All single-valued features: physico-chemical scale
    averages, structural/AlphaFold descriptors, graph descriptors, and the
    nine Ts_* thermostability scores.

features_per_residue.pkl
    One row per input structure, same index. Every cell is an array with one
    value per residue. Column names end in "_profile".

Loading them:

    import pandas as pd, joblib
    scalar  = pd.read_csv("features_scalar.csv", index_col=0,
                          float_precision="round_trip")
    profile = joblib.load("features_per_residue.pkl")

The CSV stores full float64 precision. float_precision="round_trip" makes
pandas read it back exactly; without it pandas' fast parser loses the last
couple of digits (~1e-14 relative), which is harmless for these features but
means the values will not compare equal to the originals.

To get back the single combined table:

    combined = scalar.join(profile)

features_per_residue.pkl needs pandas, numpy and joblib to read.
"""


def split_columns(df):
    """Columns whose cells are per-residue arrays vs. everything else.

    Checks every row, not just the first: a column is per-residue if any
    structure has an array there. `degrees_profile` is a plain Python list
    rather than an ndarray, so both are accepted.
    """
    profile_cols = [
        c for c in df.columns
        if any(isinstance(v, (np.ndarray, list)) for v in df[c])
    ]
    scalar_cols = [c for c in df.columns if c not in profile_cols]
    return scalar_cols, profile_cols


def main():
    parser = argparse.ArgumentParser(
        description="Package job results into a single downloadable archive")
    parser.add_argument("-D", "--dir", required=True, help="job folder")
    parser.add_argument("--keep-cmb", action="store_true", default=True,
                        help="keep the intermediate cmb_df.pkl (default)")
    args = parser.parse_args()

    folder = os.path.abspath(args.dir)
    cmb_path = os.path.join(folder, "cmb_df.pkl")
    if not os.path.isfile(cmb_path):
        print(f"no cmb_df.pkl in {folder}", file=sys.stderr)
        return 1

    df = joblib.load(cmb_path)
    scalar_cols, profile_cols = split_columns(df)
    print(f"{len(df)} structures, {df.shape[1]} columns "
          f"-> {len(scalar_cols)} scalar, {len(profile_cols)} per-residue")

    scalar_bytes = df[scalar_cols].to_csv().encode("utf-8")

    buf = io.BytesIO()
    joblib.dump(df[profile_cols], buf, compress=3)
    profile_bytes = buf.getvalue()

    archive = os.path.join(folder, ARCHIVE_NAME)
    tmp = archive + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(SCALAR_NAME, scalar_bytes)
        # already compressed by joblib -- do not deflate it again
        z.writestr(PROFILE_NAME, profile_bytes, compress_type=zipfile.ZIP_STORED)
        z.writestr("README.txt", README)
    os.replace(tmp, archive)          # never leave a half-written archive

    raw = os.path.getsize(cmb_path)
    out = os.path.getsize(archive)
    print(f"  {SCALAR_NAME:26} {len(scalar_bytes)/1024:8.1f} KB")
    print(f"  {PROFILE_NAME:26} {len(profile_bytes)/1024:8.1f} KB (joblib compress=3)")
    print(f"{ARCHIVE_NAME} = {out/1024:.1f} KB "
          f"(cmb_df.pkl was {raw/1024:.1f} KB, {raw/max(out,1):.1f}x smaller)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
