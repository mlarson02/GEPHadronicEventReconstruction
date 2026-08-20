#!/usr/bin/env python3
"""
Generate and submit MET emulation Condor jobs.

Jobs are parallelized over both algorithm configurations AND individual input
ntuples — one Condor job per (algorithm config, input file).  The INPUT_DIRS
dict maps each sample name to the directory containing its HERNTupler output
files; all *.root files found there become separate jobs.

Usage:
  python3 submit_met_emulation.py [--dry-run] [--max-jobs N] [--label LABEL] [--pu 140]

  --dry-run    print the .sub file without submitting
  --max-jobs N limit to first N jobs (useful for testing)
  --label STR  override the auto-generated label
  --pu N       pileup scenario (200 = default, 140 = the ntuples_PU140 mirrors)
"""

import argparse
import itertools
import os
import subprocess
import sys
from pathlib import Path

WRAPPER = Path(__file__).parent / "run_met_emulation_job.sh"

# ---------------------------------------------------------------------------
# Parameter grid — edit these to match metEmulationConfigLocal.sh
# ---------------------------------------------------------------------------
SIGNALS        = [True, False]
SIGNAL_STRINGS = ["ZvvHbb", "ttbar_semilep", "ttbar_dilep", "Zmumu"]   # only used when signal=True
PU_SUPPRESSION = [True, False]
JET_ET_THRESHOLDS        = [15.0]
DO_JET_TOWER_OR          = [False]
TOWER_ET_THRESHOLDS      = [2.0]
ETA_SK_OBJECTS           = [True, False]
# (towerScaleFactor, jetScaleFactor) pairs applied in the totalMET sum.
# Default (1.0, 1.0); test pair (0.4, 1.0) to down-weight tower contribution.
#SCALE_FACTOR_PAIRS       = [(1.0, 1.0), (1.0, 0.5), (0.2, 0.5), (0.4, 1.0)]
#SCALE_FACTOR_PAIRS       = [(0.4, 1.0),(0.3, 1.0),(0.2, 1.0),(0.6, 1.0),(0.1, 1.0),(0.7, 1.0)]
SCALE_FACTOR_PAIRS       = [(1.0, 1.0)]

# --- Produced 08192026: the single point missing from the Z->mumu OR comparison --------------
# EtaSK / jetEt15 / towerEt2 / NoOR / twrSF0p5 / jetSF1, for the Zmumu signal and the dijet
# background, at the default --pu 200. Both outputs exist and are merged, so this grid is done;
# kept for the record and to re-run from if either output ever has to be rebuilt.
#SIGNALS        = [True, False]        # True -> Zmumu signal, False -> BACKGROUND (JZ dijet)
#SIGNAL_STRINGS = ["Zmumu"]            # only used when signal=True
## EtaSK is selected by ETA_SK_OBJECTS=True; PU_SUPPRESSION must be True alongside it, since
## the output tag is EtaSK either way and running both values would give two jobs writing the
## same filename.
#PU_SUPPRESSION = [True]
#JET_ET_THRESHOLDS        = [15.0]
#DO_JET_TOWER_OR          = [False]    # NoOR
#TOWER_ET_THRESHOLDS      = [2.0]
#ETA_SK_OBJECTS           = [True]     # EtaSK
#SCALE_FACTOR_PAIRS       = [(0.5, 1.0)]
# Number of input HERNTupler ntuples processed per Condor job.
# Each job hadds its assigned inputs on the worker (preserving the sorted
# order) and runs metEmulation once on the merged input — this amortises
# the asetup/ROOT/ACLiC startup cost across many files.
#
# Ordering note: chunks are taken from the *sorted* input list, so chunk K
# always contains files with fidx [K*FILES_PER_JOB, (K+1)*FILES_PER_JOB).
# The output is named with the chunk's first fidx, so the global hadd of
# the chunk outputs preserves event order matching a global hadd of the
# original inputs in sorted order.
#
# Worker disk: with FILES_PER_JOB=10 and ~250 MB per input ntuple, the
# merged input on /scratch is ~2.5 GB — request_disk must cover that.
FILES_PER_JOB            = 5

# Directories containing HERNTupler output ntuples for each sample.
# All *.root files found there become separate parallel jobs.
# Key must match a signal string from SIGNAL_STRINGS, or "BACKGROUND".
_NTUPLE_BASE = "/data/larsonma/GEPHadronicEventReconstruction/ntuples"
# PU140 ntuples live in a parallel directory tree with identical sample
# subdirectory names (see submit_all_ntupler.sh --pu 140). Selected with --pu 140;
# every INPUT_DIRS path below is remapped onto this base by ntuple_dir().
_NTUPLE_BASE_BY_PU = {
    200: _NTUPLE_BASE,
    140: "/data/larsonma/GEPHadronicEventReconstruction/ntuples_PU140",
}
# Pileup scenario for this submission; set from --pu in main().
PILEUP = 200
INPUT_DIRS = {
    "ZvvHbb":        f"{_NTUPLE_BASE}/ZvvHbb_v4/",
    "ttbar_semilep": f"{_NTUPLE_BASE}/ttbar_semilep_v4/",
    "ttbar_dilep":   f"{_NTUPLE_BASE}/ttbar_dilep_v4/",
    "Zmumu":         f"{_NTUPLE_BASE}/Zmumu_v4/",
    "BACKGROUND": [
        f"{_NTUPLE_BASE}/QCD_Dijet_JZ0_v4/",
        f"{_NTUPLE_BASE}/QCD_Dijet_JZ1_v4/",
        f"{_NTUPLE_BASE}/QCD_Dijet_JZ2_v4/",
        f"{_NTUPLE_BASE}/QCD_Dijet_JZ3_v4/",
        f"{_NTUPLE_BASE}/QCD_Dijet_JZ4_v4/",
        f"{_NTUPLE_BASE}/QCD_Dijet_JZ5_v4/",
        f"{_NTUPLE_BASE}/QCD_Dijet_JZ6_v4/",
        f"{_NTUPLE_BASE}/QCD_Dijet_JZ7_v4/",
        f"{_NTUPLE_BASE}/QCD_Dijet_JZ8_v4/",
        f"{_NTUPLE_BASE}/QCD_Dijet_JZ9_v4/",
    ],
}

# ---------------------------------------------------------------------------


def fmt_float(v: float) -> str:
    int_v = int(v)
    return str(int_v) if v == int_v else f"{v:.4g}"


def fmt_scale_factor_tag(v: float) -> str:
    """Filesystem-safe scale-factor tag, mirrors makeOutputMETFileName.
    Assumes SF granularity of 0.1: 1.0 -> "1", 0.4 -> "0p4", 1.2 -> "1p2"."""
    int_v = int(v)
    if v == int_v:
        return str(int_v)
    return f"{v:.1f}".replace(".", "p")


def bool_str(b: bool) -> str:
    return "true" if b else "false"


def ntuple_dir(path: str) -> str:
    """Remap an INPUT_DIRS path (written for PU200) onto the selected pileup's
    ntuple base. PU140 mirrors the PU200 tree, so only the base differs."""
    base = _NTUPLE_BASE_BY_PU[PILEUP]
    return path.replace(_NTUPLE_BASE, base, 1) if path.startswith(_NTUPLE_BASE) else path


def sample_dirs(sample: str) -> list[str]:
    """Directories to search for the given sample, for the selected pileup."""
    entry = INPUT_DIRS.get(sample)
    if not entry:
        return []
    dirs = [entry] if isinstance(entry, str) else entry
    return [ntuple_dir(d) for d in dirs]


def find_input_files(sample: str) -> list[str]:
    """Return sorted list of *.root files for the given sample.
    Accepts a single directory or a list of directories (used for JZ slices)."""
    dirs = sample_dirs(sample)
    if not dirs:
        print(f"[warn] No INPUT_DIRS entry for sample '{sample}', skipping file discovery")
        return []
    files = []
    for d in dirs:
        p = Path(d)
        if not p.is_dir():
            print(f"[warn] Input directory not found for '{sample}': {d}")
            continue
        # Exclude merged ntuples (end in GEP.root rather than GEP_NNNNNN.root)
        files.extend(sorted(str(f) for f in p.glob("*.root") if not f.name.endswith("GEP.root")))
    if not files:
        print(f"[warn] No *.root files found for sample '{sample}'")
    return files


def print_file_index_map(samples: list[str]) -> None:
    """Print the fidx → input file mapping for each sample."""
    for sample in samples:
        dirs = sample_dirs(sample)
        print(f"\n=== Sample: {sample} (PU{PILEUP}) ===")
        if not dirs:
            print("  [no INPUT_DIRS entry]")
            continue
        print(f"  Source directories ({len(dirs)}):")
        for d in dirs:
            exists = "OK" if Path(d).is_dir() else "MISSING"
            print(f"    [{exists}] {d}")
        files = find_input_files(sample)
        if not files:
            print("  -> No files found; jobs will use fidx=0 infile='' (makeInputFileName fallback)")
            continue
        print(f"  fidx  basename                                   full path")
        print(f"  ----  -----------------------------------------  ---------")
        for fidx, fpath in enumerate(files):
            print(f"  {fidx:4d}  {Path(fpath).name:<41}  {fpath}")
        print(f"  -> {len(files)} file(s) total")


def make_submit_file(jobs: list[dict], wrapper: str, label: str,
                     itemdata_path: str) -> tuple[str, str]:
    log_dir = jobs[0]["log_dir"]

    submit_lines = [
        f"# Auto-generated by submit_met_emulation.py for: {label}",
        "",
        f"executable            = {wrapper}",
        "universe              = vanilla",
        "request_cpus          = 1",
        "request_memory        = 2048",
        # Disk in KB. Wrapper hadds FILES_PER_JOB inputs on /scratch before
        # running emulation; size with headroom for ~250 MB/input * 2 (original
        # + merged) plus ~200 MB for ROOT/asetup files.
        f"request_disk          = {max(4096, FILES_PER_JOB * 500_000 + 500_000)}",
        # Default queue (no 4 h wall-clock cap). Re-add `+queue = "short"` once
        # jobs reliably finish well under 4 h, to access the larger short-queue slot pool.
        "getenv                = false",
        "",
        # pileup is fixed for the whole submission, so it goes in the header rather
        # than the per-job ItemData.
        f"arguments = $(SIG) $(PUSUP) $(SIGSTR) $(JETET) $(JTOR) $(TOWERET) $(ETASK) $(INFILE) $(FIDX) $(TWRSF) $(JETSF) {PILEUP}",
        "log       = $(LOG)",
        "output    = $(OUT)",
        "error     = $(ERR)",
        "",
        f"queue SIG, PUSUP, SIGSTR, JETET, JTOR, TOWERET, ETASK, INFILE, FIDX, TWRSF, JETSF, LOG, OUT, ERR from {itemdata_path}",
    ]

    data_lines = []
    for j in jobs:
        data_lines.append(
            f"{j['signal']}, {j['pusup']}, {j['sigstr']}, {j['jetet']}, "
            f"{j['jtor']}, {j['toweret']}, {j['etask']}, "
            f"{j['infile']}, {j['fidx']}, "
            f"{j['twrsf']}, {j['jetsf']}, "
            f"{j['log']}, {j['stdout']}, {j['stderr']}"
        )

    return "\n".join(submit_lines), "\n".join(data_lines)


def enumerate_jobs(log_dir: str, label: str) -> list[dict]:
    jobs = []
    for pusup, jetet, jtor, toweret, etask, sfpair in itertools.product(
            PU_SUPPRESSION, JET_ET_THRESHOLDS, DO_JET_TOWER_OR,
            TOWER_ET_THRESHOLDS, ETA_SK_OBJECTS, SCALE_FACTOR_PAIRS):
        twrsf, jetsf = sfpair
        for signal in SIGNALS:
            sig_list = SIGNAL_STRINGS if signal else ["BACKGROUND"]
            for sigstr in sig_list:
                algo_tag = (f"sig_{sigstr}" if signal else "bkg"
                            ) + (f"_{('EtaSK' if etask else 'SK' if pusup else 'NoSK')}"
                                 f"_jetEt{fmt_float(jetet)}"
                                 f"_towerEt{fmt_float(toweret)}"
                                 f"_{'OR' if jtor else 'NoOR'}"
                                 f"_twrSF{fmt_scale_factor_tag(twrsf)}"
                                 f"_jetSF{fmt_scale_factor_tag(jetsf)}")

                input_files = find_input_files(sigstr)
                if not input_files:
                    input_files = [""]

                # Chunk the sorted input list into groups of FILES_PER_JOB.
                # Each chunk becomes one Condor job; the wrapper hadds the
                # chunk's inputs on the worker (preserving this order) and
                # runs metEmulation once on the merged file. The chunk's
                # output is named with its first fidx so the global hadd
                # of all jobs' outputs preserves event order.
                step = max(1, FILES_PER_JOB)
                for chunk_start in range(0, len(input_files), step):
                    chunk_paths = input_files[chunk_start:chunk_start + step]
                    chunk_fidxs = list(range(chunk_start, chunk_start + len(chunk_paths)))
                    first_fidx  = chunk_fidxs[0]
                    last_fidx   = chunk_fidxs[-1]
                    # `;` chosen as separator since Condor's submit-file ItemData
                    # parses on `,`. The wrapper splits on `;` and hadds in order.
                    infile_csv = ";".join(chunk_paths)
                    fidx_csv   = ";".join(str(i) for i in chunk_fidxs)
                    chunk_tag  = (f"file{first_fidx}" if first_fidx == last_fidx
                                  else f"files{first_fidx}to{last_fidx}")
                    safe = f"{label}_{algo_tag}_{chunk_tag}"
                    jobs.append({
                        "signal":  bool_str(signal),
                        "pusup":   bool_str(pusup),
                        "sigstr":  sigstr,
                        "jetet":   fmt_float(jetet),
                        "jtor":    bool_str(jtor),
                        "toweret": fmt_float(toweret),
                        "etask":   bool_str(etask),
                        "infile":  infile_csv,
                        "fidx":    first_fidx,    # used by metEmulation to name the output
                        "fidx_list": fidx_csv,    # full list for the wrapper's log
                        "twrsf":   fmt_float(twrsf),
                        "jetsf":   fmt_float(jetsf),
                        "log_dir": log_dir,
                        "log":     os.path.join(log_dir, f"{safe}.log"),
                        "stdout":  os.path.join(log_dir, f"{safe}.out"),
                        "stderr":  os.path.join(log_dir, f"{safe}.err"),
                    })
    return jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the .sub file but do not submit")
    parser.add_argument("--max-jobs", type=int, default=0, metavar="N",
                        help="Limit to first N jobs (0 = no limit)")
    parser.add_argument("--label", default="met_emulation",
                        help="Label prefix for log files and submit file (default: met_emulation)")
    parser.add_argument("--debug-files", action="store_true",
                        help="Print the fidx→filename mapping for every sample and exit")
    parser.add_argument("--lookup-fidx", type=int, nargs="+", metavar="N",
                        help="Look up which input file(s) correspond to the given output _fileN number(s) and exit")
    parser.add_argument("--pu", type=int, default=200, choices=[140, 200], metavar="PILEUP",
                        help="Pileup scenario (default: 200). 140 reads the ntuples_PU140 "
                             "mirror tree instead of ntuples.")
    args = parser.parse_args()

    global PILEUP
    PILEUP = args.pu

    all_samples = list(SIGNAL_STRINGS) + ["BACKGROUND"]

    if args.debug_files:
        print("File index map (fidx = value passed as $(FIDX) to each job):")
        print_file_index_map(all_samples)
        return

    if args.lookup_fidx:
        targets = set(args.lookup_fidx)
        print(f"Looking up fidx values: {sorted(targets)}\n")
        for sample in all_samples:
            files = find_input_files(sample)
            hits = {i: f for i, f in enumerate(files) if i in targets}
            if hits:
                print(f"Sample: {sample}")
                for fidx in sorted(hits):
                    print(f"  fidx {fidx:5d}  ->  {hits[fidx]}")
        return

    # Tag non-default pileup so PU140 logs/itemdata/sub files don't collide with PU200.
    if args.pu != 200:
        args.label += f"_PU{args.pu}"
        print(f"[pileup] PU{args.pu}: reading inputs from {_NTUPLE_BASE_BY_PU[args.pu]}")
        print(f"[pileup] Outputs are tagged r16129 (PU140) instead of r16130 (PU200), "
              "so they sit alongside the PU200 ones without overwriting them.")

    log_dir = str(Path.home() / "condor_logs" / args.label)
    os.makedirs(log_dir, exist_ok=True)

    jobs = enumerate_jobs(log_dir, args.label)
    print(f"Total jobs (algorithm configs × input files): {len(jobs)}")

    if args.max_jobs > 0:
        jobs = jobs[:args.max_jobs]
        print(f"Limiting to first {args.max_jobs} job(s) (--max-jobs)")

    itemdata_path = os.path.join(log_dir, f"{args.label}.dat")
    sub_path      = os.path.join(log_dir, f"{args.label}.sub")

    submit_text, itemdata_text = make_submit_file(jobs, str(WRAPPER), args.label, itemdata_path)

    with open(sub_path, "w") as fh:
        fh.write(submit_text)
    with open(itemdata_path, "w") as fh:
        fh.write(itemdata_text)
    print(f"Wrote submit file ({len(jobs)} job(s)): {sub_path}")
    print(f"Log directory: {log_dir}")

    if args.dry_run:
        print("\n--- file index map (fidx → input file per sample) ---")
        print_file_index_map(all_samples)
        print("\n--- .sub file ---")
        print(submit_text)
        print("\n--- itemdata (first 3 lines) ---")
        for line in itemdata_text.splitlines()[:3]:
            print(line)
        print("--- end (dry run, not submitted) ---")
        return

    result = subprocess.run(["condor_submit", sub_path], capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print("[error] condor_submit failed:", result.stderr)
        sys.exit(result.returncode)


if __name__ == "__main__":
    main()
