"""Run the Thermometer pipeline over a directory of PDB structures.

Usage
-----
    python thermometer_new_parallel.py -D <pdb_dir> -O <output_dir>

Produces ``<output_dir>/RESULTS/`` with, per structure:

    <name>_Ts_global_score.csv        <- the 9 Ts_* features the app uses
    <name>_Strength_Global.csv
    <name>_RealRandom_Strength.csv
    <name>_Strength_SingleResidue.csv
    ENERGY/<name>.csv

and ``<output_dir>/thermometer_failures.txt`` listing structures that failed.

Exit status is non-zero if any structure failed.

Changes from the previous version
---------------------------------
* ``get_ca_count`` referenced an undefined ``pdb`` and raised NameError for
  every structure. ``executor.map`` results were never consumed, so
  ProcessPoolExecutor discarded those exceptions and the whole run was a
  silent no-op that still exited 0. The function is gone (the ``max_res``
  value it computed could never fail the assertion it fed) and results are
  now consumed so failures surface.
* Every structure gets its own work directory. Nothing is written into the
  tool directory, so concurrent jobs cannot collide and cleanup is a single
  rmtree that runs on success *and* on failure.
* The tool directory is located via THERMO_HOME, which the R/perl/shell
  scripts now honour too.
* The hardcoded ``mv NAMD-MIN/<dir>/AF*`` rename is gone. It existed to work
  around run_namd_2.pl truncating outputName to 4 characters, and it silently
  did nothing for any structure not named ``AF-*`` -- which is every
  user-uploaded file. run_namd_2.pl now uses the full basename.
* No ``shell=True``. Structure names come from user uploads.
* Every subprocess return code is checked; before, six ``check_*`` variables
  were assigned and never looked at.
"""

import argparse
import concurrent.futures as cf
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from glob import glob


THERMO_HOME = os.environ.get(
    "THERMO_HOME", os.path.dirname(os.path.abspath(__file__))
)
NAMD_PROCS = int(os.environ.get("NAMD_PROCS", "2"))

# check_pdb.R asserts nrow(CA) <= this. The old caller passed
# max(1000, ca_count), which can never fail, so there was effectively no cap.
# 0 means "no cap"; the web app should set a real one.
DEFAULT_MAX_RESIDUES = 0


class StepError(RuntimeError):
    """A pipeline step failed for one structure."""


def run_step(cmd, cwd, name, log):
    """Run one pipeline step, recording output and raising on failure."""
    env = dict(os.environ, THERMO_HOME=THERMO_HOME, NAMD_PROCS=str(NAMD_PROCS))
    proc = subprocess.run(
        cmd, cwd=cwd, env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    log.write(f"\n$ {' '.join(cmd)}\n{proc.stdout or ''}")
    if proc.returncode != 0:
        tail = "\n".join((proc.stdout or "").strip().splitlines()[-15:])
        raise StepError(f"{name}: {cmd[0]} exited {proc.returncode}\n{tail}")
    return proc


def write_tcl(path, clean_pdb, out_pdb, out_psf):
    """psfgen input. All paths absolute, so vmd's cwd does not matter and the
    outputs land where we want instead of needing a follow-up mv."""
    topology = os.path.join(THERMO_HOME, "files", "top_all27_prot_lipid.inp")
    with open(path, "w") as fp:
        fp.write(
            "package require psfgen\n"
            f"topology {topology}\n"
            "pdbalias residue HIS HSD\n"
            "pdbalias atom ILE CD1 CD\n"
            f"segment U {{pdb {clean_pdb}}}\n"
            f"coordpdb {clean_pdb} U\n"
            "guesscoord\n"
            f"writepdb {out_pdb}\n"
            f"writepsf {out_psf}\n"
            "exit\n"
        )


def link_or_copy(src, dst):
    """Symlink the big shared parameter files instead of copying them per
    structure (prova.inp alone is 183 KB)."""
    try:
        os.symlink(src, dst)
    except (OSError, NotImplementedError, AttributeError):
        shutil.copy(src, dst)


def run_thermometer(args):
    """Process one structure. Returns (name, error_or_None)."""
    pdb_path, work_root, results_dir, max_residues = args
    name = os.path.splitext(os.path.basename(pdb_path))[0]

    work = os.path.join(work_root, name)
    run_dir = os.path.join(work, "NAMD-MIN", name)
    pdb_dataset = os.path.join(run_dir, "pdb_dataset")
    psf_dataset = os.path.join(run_dir, "psf_dataset")
    for d in (pdb_dataset, psf_dataset):
        os.makedirs(d, exist_ok=True)

    log_path = os.path.join(work, f"{name}.log")
    print(f"running for {name}", flush=True)

    try:
        with open(log_path, "w") as log:
            local_pdb = os.path.join(work, f"{name}.pdb")
            shutil.copy(pdb_path, local_pdb)

            # 1. validate + renumber -> <name>_clean.pdb, <name>_dict_tmp.csv
            cap = str(max_residues) if max_residues else "100000"
            run_step(
                ["Rscript", os.path.join(THERMO_HOME, "scripts", "check_pdb.R"),
                 local_pdb, cap],
                work, name, log,
            )
            clean_pdb = os.path.join(work, f"{name}_clean.pdb")
            if not os.path.isfile(clean_pdb):
                raise StepError("check_pdb.R produced no _clean.pdb "
                                "(structure rejected -- see log)")

            # 2. psfgen
            tcl = os.path.join(run_dir, "autopsf.tcl")
            write_tcl(
                tcl, clean_pdb,
                os.path.join(pdb_dataset, f"{name}_autopsf.pdb"),
                os.path.join(psf_dataset, f"{name}_autopsf.psf"),
            )
            run_step(["vmd", "-dispdev", "text", "-e", tcl], work, name, log)
            if not os.path.isfile(os.path.join(pdb_dataset, f"{name}_autopsf.pdb")):
                raise StepError("vmd/psfgen produced no _autopsf.pdb")

            # 3. NAMD minimisation
            for f in ("MIN_general.conf", "prova.inp"):
                link_or_copy(os.path.join(THERMO_HOME, "files", f),
                             os.path.join(run_dir, f))
            shutil.copy(os.path.join(THERMO_HOME, "scripts", "run_namd_2.pl"),
                        os.path.join(run_dir, "run_namd_2.pl"))
            run_step(["perl", "run_namd_2.pl", "./pdb_dataset"], run_dir, name, log)

            # 4. last frame -> pdbMin/<name>_min.pdb
            dcd_dir = os.path.join(run_dir, "DCD")
            os.makedirs(dcd_dir, exist_ok=True)
            for f in os.listdir(run_dir):
                if f.endswith(".dcd"):
                    shutil.move(os.path.join(run_dir, f), os.path.join(dcd_dir, f))
            run_step(
                ["bash", os.path.join(THERMO_HOME, "scripts", "catDCD.sh"),
                 work, name],
                work, name, log,
            )

            min_pdb = os.path.join(work, f"{name}_min.pdb")
            shutil.move(os.path.join(run_dir, "pdbMin", f"{name}_min.pdb"), min_pdb)

            # 5. restore the original residue numbering (in place)
            run_step(
                ["Rscript", os.path.join(THERMO_HOME, "scripts", "Renumbering_Res.R"),
                 os.path.join(work, f"{name}_dict_tmp.csv"), min_pdb],
                work, name, log,
            )

            # 6. energies -> <name>.csv
            run_step(
                ["Rscript", os.path.join(THERMO_HOME, "scripts", "Energy_computer.R"),
                 min_pdb, name, "1"],
                work, name, log,
            )
            energy_csv = os.path.join(work, f"{name}.csv")
            if not os.path.isfile(energy_csv):
                raise StepError("Energy_computer.R produced no csv")

            # 7. global + per-residue Ts scores
            run_step(
                ["Rscript", os.path.join(THERMO_HOME, "scripts", "Ts_global_computer.R"),
                 energy_csv, "1", name],
                work, name, log,
            )
            run_step(
                ["Rscript", os.path.join(THERMO_HOME, "scripts", "Ts_local_computer.R"),
                 energy_csv, "1", name],
                work, name, log,
            )

            # 8. keep the results, drop everything else
            os.makedirs(os.path.join(results_dir, "ENERGY"), exist_ok=True)
            shutil.move(energy_csv,
                        os.path.join(results_dir, "ENERGY", f"{name}.csv"))
            for suffix in ("_Ts_global_score.csv", "_Strength_Global.csv",
                           "_RealRandom_Strength.csv", "_Strength_SingleResidue.csv"):
                src = os.path.join(work, f"{name}{suffix}")
                if os.path.isfile(src):
                    shutil.move(src, os.path.join(results_dir, f"{name}{suffix}"))

            ts = os.path.join(results_dir, f"{name}_Ts_global_score.csv")
            if not os.path.isfile(ts):
                raise StepError("no Ts_global_score.csv was produced")

    except Exception as exc:                      # noqa: BLE001 - reported, not raised
        # Keep this structure's work dir so the failure can be inspected.
        print(f"FAILED {name}: {exc}", flush=True)
        return name, str(exc)

    # NAMD-MIN is ~7 MB per structure; nothing downstream needs it.
    shutil.rmtree(work, ignore_errors=True)
    print(f"ending {name}", flush=True)
    return name, None


def report(progress_file, pct, status):
    if not progress_file:
        return
    tmp = f"{progress_file}.tmp"
    with open(tmp, "w") as f:
        json.dump({"progress": round(pct), "status": status}, f)
    os.replace(tmp, progress_file)      # atomic: the poller never sees a partial file


def main():
    parser = argparse.ArgumentParser(
        description="Run Thermometer over a directory of PDB files")
    parser.add_argument("-D", "--dir", required=True,
                        help="directory containing the pdb files")
    parser.add_argument("-O", "--output", required=True,
                        help="directory to receive RESULTS/")
    parser.add_argument("-W", "--work", default=None,
                        help="work directory (default: a temp dir under --output)")
    parser.add_argument("--keep-work", action="store_true",
                        help="do not delete the work directory")
    parser.add_argument("-j", "--jobs", type=int, default=0,
                        help="parallel structures (default: cpus // NAMD_PROCS)")
    parser.add_argument("--max-residues", type=int, default=DEFAULT_MAX_RESIDUES,
                        help="reject structures with more CA atoms than this "
                             "(0 = no cap)")
    parser.add_argument("--progress-file", default=None,
                        help="write progress.json here while running")
    parser.add_argument("--progress-start", type=float, default=0.0)
    parser.add_argument("--progress-end", type=float, default=100.0)
    args = parser.parse_args()

    pdbs = sorted(glob(os.path.join(args.dir.rstrip("/"), "*.pdb")))
    if not pdbs:
        print(f"no *.pdb found in {args.dir}", file=sys.stderr)
        return 1

    output = os.path.abspath(args.output)
    results_dir = os.path.join(output, "RESULTS")
    os.makedirs(results_dir, exist_ok=True)

    if args.work:
        work_root, owns_work = os.path.abspath(args.work), False
        os.makedirs(work_root, exist_ok=True)
    else:
        work_root, owns_work = tempfile.mkdtemp(prefix="thermo_", dir=output), True

    workers = args.jobs or max(1, (os.cpu_count() or 2) // max(1, NAMD_PROCS))
    workers = min(workers, len(pdbs))

    print(f"THERMO_HOME={THERMO_HOME}")
    print(f"work={work_root}  structures={len(pdbs)}  "
          f"workers={workers}  namd_procs={NAMD_PROCS}")

    span = args.progress_end - args.progress_start
    report(args.progress_file, args.progress_start,
           f"Thermostability (0/{len(pdbs)})...")

    tic = time.time()
    failures = []
    try:
        payload = [(p, work_root, results_dir, args.max_residues) for p in pdbs]
        with cf.ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(run_thermometer, x) for x in payload]
            for done, fut in enumerate(cf.as_completed(futures), start=1):
                name, err = fut.result()        # consumed, so nothing is swallowed
                if err:
                    failures.append((name, err))
                report(args.progress_file,
                       args.progress_start + span * done / len(pdbs),
                       f"Thermostability ({done}/{len(pdbs)})...")
    finally:
        # Runs on success, on exception and on KeyboardInterrupt -- unlike the
        # rm lines that used to sit at the end of the wrapper script.
        if owns_work and not args.keep_work and not failures:
            shutil.rmtree(work_root, ignore_errors=True)

    elapsed = (time.time() - tic) / 60
    print(f"elapsed time for {len(pdbs)} structures is {elapsed:.2f} minutes")

    fail_file = os.path.join(output, "thermometer_failures.txt")
    if failures:
        with open(fail_file, "w") as f:
            for name, err in failures:
                f.write(f"{name}\t{err.splitlines()[0]}\n")
        print(f"{len(failures)}/{len(pdbs)} structures FAILED "
              f"(details in {fail_file}, work kept at {work_root})",
              file=sys.stderr)
        return 1

    if os.path.exists(fail_file):
        os.remove(fail_file)
    print(f"{len(pdbs)}/{len(pdbs)} structures OK -> {results_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
