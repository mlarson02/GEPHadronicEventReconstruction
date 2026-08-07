#!/usr/bin/env python3
"""
Generate and submit JetTagger emulation Condor jobs.

Jobs are parallelized over both algorithm configurations AND individual input
ntuples — one Condor job per (algorithm config, input file).  The INPUT_DIRS
dict maps each sample name to the directory containing its HERNTupler output
files; all *.root files found there become separate jobs.

Usage:
  python3 submit_jet_tagger_emulation.py [--dry-run] [--max-jobs N] [--label LABEL]

  --dry-run    print the .sub file without submitting
  --max-jobs N limit to first N jobs (useful for testing)
  --label STR  override the auto-generated label
"""

import argparse
import itertools
import os
import subprocess
import sys
from pathlib import Path

WRAPPER = Path(__file__).parent / "run_jet_tagger_emulation_job.sh"

# ---------------------------------------------------------------------------
# Parameter grid — edit these to match jetTaggerConfigLocal.sh
# ---------------------------------------------------------------------------
ALGO_VERSIONS  = [3]
R_MERGE_CUTS   = [2, 0.001]
R_SQUARED_CUTS = [1.21]
N_IOS          = [128]
N_SEEDS        = [2]
SIGNALS        = [True, False]
SIGNAL_STRINGS = ["ggF_hh_bbbb", "Zprime_ttbar_allhad_flatpT", "ttbar_allhad"]   # only used when signal=True
INPUT_OBJECTS  = ["gepCellsTowers"]
SEED_OBJECTS   = ["gepWTAConeCellsTowersJets"]
PU_SUPPRESSION = [True]
ETA_SK_OBJECTS = [False]   # True = use EtaSK PU-suppressed towers+jets (gepCellsTowers and WTAConeJets only)
ET_WEIGHTED_MIDPOINTS      = [False]
MIN_ET_SEED_POS_OPT        = [True]
MIN_ET_SEED_POS_OPT_CUTS   = [20.0]

SUBJET_ET_BY_SEED = {
    "gFEXSRJ":                   25,
    "jFEXSRJ":                   35,
    "gepWTAConeCellsTowersJets":  25,
}

# Directories containing HERNTupler output ntuples for each sample.
# All *.root files found are submitted as separate parallel jobs.
# Key must match a signal string from SIGNAL_STRINGS, or "BACKGROUND".
# Values may be a single path string or a list of paths (used for the
# 10 JZ dijet slices, which are pooled together under "BACKGROUND").
_NTUPLE_BASE = "/data/larsonma/GEPHadronicEventReconstruction/ntuples"
INPUT_DIRS = {
    "ggF_hh_bbbb":                f"{_NTUPLE_BASE}/ggF_HHbbbb_v4/",
    "ZvvHbb":                     f"{_NTUPLE_BASE}/ZvvHbb_v4/",
    "Zprime_ttbar_allhad_flatpT": f"{_NTUPLE_BASE}/Zprime_ttbar_allhad_flatpT_v4/",
    "ttbar_allhad":               f"{_NTUPLE_BASE}/ttbar_allhad_v4/",
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


def subjet_et(seed_obj: str) -> int:
    return SUBJET_ET_BY_SEED.get(seed_obj, 25)


def fmt_rmerge(v: float) -> str:
    return f"{v:.4g}"


def fmt_r2(v: float) -> str:
    return f"{v:.3g}"


def bool_str(b: bool) -> str:
    return "true" if b else "false"


def find_input_files(sample: str) -> list[str]:
    """Return sorted list of *.root files for the given sample.
    Accepts a single directory or a list of directories (used for JZ slices)."""
    entry = INPUT_DIRS.get(sample)
    if not entry:
        print(f"[warn] No INPUT_DIRS entry for sample '{sample}', skipping file discovery")
        return []
    dirs = [entry] if isinstance(entry, str) else entry
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
    """Print the fidx → input file mapping for each sample.

    This is the authoritative index that jobs receive as $(FIDX) and $(INFILE).
    Use it to verify alignment between input ntuple numbering and output numbering.
    """
    for sample in samples:
        entry = INPUT_DIRS.get(sample)
        dirs = [entry] if isinstance(entry, str) else (entry or [])
        print(f"\n=== Sample: {sample} ===")
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
        f"# Auto-generated by submit_jet_tagger_emulation.py for: {label}",
        "",
        f"executable            = {wrapper}",
        "universe              = vanilla",
        "request_cpus          = 1",
        "request_memory        = 4096",
        "request_disk          = 4096",
        '+queue                = "short"',
        "getenv                = false",
        "",
        "arguments = $(RMRG) $(R2) $(NIOS) $(NSEEDS) $(ALGOV) $(SIG) $(SIGSTR) $(PUSUP) $(INOBJ) $(SEEDOBJ) $(SJETT) $(ETWM) $(MINETO) $(MINETC) $(INFILE) $(FIDX) $(ETASK)",
        "log       = $(LOG)",
        "output    = $(OUT)",
        "error     = $(ERR)",
        "",
        f"queue RMRG, R2, NIOS, NSEEDS, ALGOV, SIG, SIGSTR, PUSUP, INOBJ, SEEDOBJ, SJETT, ETWM, MINETO, MINETC, INFILE, FIDX, ETASK, LOG, OUT, ERR from {itemdata_path}",
    ]

    data_lines = []
    for j in jobs:
        data_lines.append(
            f"{j['rmrg']}, {j['r2']}, {j['nios']}, {j['nseeds']}, {j['algov']}, "
            f"{j['signal']}, {j['sigstr']}, {j['pusup']}, {j['inobj']}, {j['seedobj']}, "
            f"{j['sjett']}, {j['etwm']}, {j['mineto']}, {j['minetc']}, "
            f"{j['infile']}, {j['fidx']}, {j['etask']}, "
            f"{j['log']}, {j['stdout']}, {j['stderr']}"
        )

    return "\n".join(submit_lines), "\n".join(data_lines)


def enumerate_jobs(log_dir: str, label: str) -> list[dict]:
    jobs = []
    for algov, rmrg, r2, nios, nseeds in itertools.product(
            ALGO_VERSIONS, R_MERGE_CUTS, R_SQUARED_CUTS, N_IOS, N_SEEDS):
        if algov == 2 and rmrg != 0.001:
            continue
        for inobj, seedobj in itertools.product(INPUT_OBJECTS, SEED_OBJECTS):
            sjett = subjet_et(seedobj)
            for etwm, mineto, minetc in itertools.product(
                    ET_WEIGHTED_MIDPOINTS, MIN_ET_SEED_POS_OPT, MIN_ET_SEED_POS_OPT_CUTS):
                for signal in SIGNALS:
                    sig_list = SIGNAL_STRINGS if signal else ["BACKGROUND"]
                    for sigstr in sig_list:
                        for pusup in PU_SUPPRESSION:
                            for etask in ETA_SK_OBJECTS:
                                # EtaSK only supported for gepCellsTowers input + WTAConeJets seeds
                                if etask and inobj not in ("gepCellsTowers",) or \
                                   etask and seedobj not in ("gepWTAConeCellsTowersJets",):
                                    continue
                                if etask:
                                    pu_tag = "EtaSK"
                                elif pusup:
                                    pu_tag = "SK"
                                else:
                                    pu_tag = "NoSK"
                                algo_tag = (f"rM{fmt_rmerge(rmrg)}_R2{fmt_r2(r2)}_IOs{nios}"
                                            f"_Sd{nseeds}_v{algov}"
                                            f"_{'sig' if signal else 'bkg'}"
                                            f"_{sigstr if signal else 'bkg'}"
                                            f"_{pu_tag}"
                                            f"_etwm{int(etwm)}_mep{int(mineto)}"
                                            f"_mec{minetc:.4g}")

                                # Discover per-file parallelism
                                input_files = find_input_files(sigstr)
                                if not input_files:
                                    # Fall back: single job with empty path (uses makeInputFileName)
                                    input_files = [""]

                                for fidx, infile_path in enumerate(input_files):
                                    safe = f"{label}_{algo_tag}_file{fidx}"
                                    jobs.append({
                                        "rmrg":    fmt_rmerge(rmrg),
                                        "r2":      fmt_r2(r2),
                                        "nios":    nios,
                                        "nseeds":  nseeds,
                                        "algov":   algov,
                                        "signal":  bool_str(signal),
                                        "sigstr":  sigstr,
                                        "pusup":   bool_str(pusup),
                                        "inobj":   inobj,
                                        "seedobj": seedobj,
                                        "sjett":   sjett,
                                        "etwm":    bool_str(etwm),
                                        "mineto":  bool_str(mineto),
                                        "minetc":  minetc,
                                        "infile":  infile_path,
                                        "fidx":    fidx,
                                        "etask":   bool_str(etask),
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
    parser.add_argument("--label", default="jet_tagger_emulation",
                        help="Label prefix for log files and submit file (default: jet_tagger_emulation)")
    parser.add_argument("--debug-files", action="store_true",
                        help="Print the fidx→filename mapping for every sample and exit")
    parser.add_argument("--lookup-fidx", type=int, nargs="+", metavar="N",
                        help="Look up which input file(s) correspond to the given output _fileN number(s) and exit")
    args = parser.parse_args()

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
