#!/bin/bash
#
# Extract the final minimised frame of each NAMD trajectory as a PDB.
#
#   usage: catDCD.sh WORK_DIR RUN_DIR
#
#   WORK_DIR  dir containing NAMD-MIN/  (a per-job work dir, or "." for the
#             historical in-tool-dir layout)
#   RUN_DIR   the NAMD-MIN/<name> subdirectory for this structure
#
#   THERMO_HOME  root of the thermometer tool, for bin/catdcd.
#                Falls back to WORK_DIR so old call sites keep working.
#
# Changes from the original:
#   - was "#bin/bash" (not a shebang at all)
#   - split tool root from work root; the original used $1 for both
#   - set -euo pipefail + nullglob; the original swallowed every failure and
#     still exited 0, which is how a broken chain looked like success
#   - quoted every path, mkdir -p, and the catdcd exit code is now checked

set -euo pipefail
shopt -s nullglob

if [ $# -ne 2 ]; then
    echo "usage: $0 WORK_DIR RUN_DIR" >&2
    exit 2
fi

work_dir=$1
run_dir=$2
tool_home=${THERMO_HOME:-$work_dir}

# The frame to extract. minimize 10000 + dcdfreq 100 (files/MIN_general.conf)
# produce 100 frames, and we want the last one. Keep these in sync if you
# change the minimisation length.
first_step=${CATDCD_FIRST_STEP:-100}
last_step=${CATDCD_LAST_STEP:-100}

catdcd="$tool_home/bin/catdcd"
pdb_input="$work_dir/NAMD-MIN/$run_dir/pdb_dataset"
pdb_output="$work_dir/NAMD-MIN/$run_dir/pdbMin"
dcd_dir="$work_dir/NAMD-MIN/$run_dir/DCD"

mkdir -p "$pdb_output"

found=0
for filename in "$dcd_dir"/*.dcd; do
    found=1
    name=${filename##*/}
    base=${name%.dcd}

    structure="$pdb_input/${base}_autopsf.pdb"
    if [ ! -f "$structure" ]; then
        # Before run_namd_2.pl was fixed, NAMD truncated outputName to 4 chars
        # and this is exactly where the chain broke -- silently.
        echo "catDCD.sh: no structure file for '$base' at $structure" >&2
        exit 1
    fi

    "$catdcd" -o "$pdb_output/${base}_min.pdb" -otype pdb \
              -s "$structure" -stype pdb \
              -first "$first_step" -last "$last_step" \
              "$filename"
done

if [ "$found" -eq 0 ]; then
    echo "catDCD.sh: no .dcd files in $dcd_dir -- NAMD produced no trajectory" >&2
    exit 1
fi
