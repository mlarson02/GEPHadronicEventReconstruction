#!/usr/bin/env python3
"""
Generate and submit caloShowerShapeNTupler Condor jobs for all DAOD files in a
sample directory.

This is the shower-shape-study analogue of submit_ntupler.py: each job processes
one DAOD_JETM42 pool.root file paired with its matching GEPOutputReader ntuple
file (paired by the 6-digit job number embedded in both filenames as _NNNNNN).
When only one GEP file exists in the GEP directory (no job number), it is used
for all jobs.

Unlike submit_ntupler.py there is no --algo / --special-jz0: the shower shape
ntupler just reads EtaSK GEP objects + truth BSM. The per-sample knobs are the
signal label (used to tag the output), signal vs background, and — for QCD dijet
background — the JZ slice, which the ntupler needs to compute each event's
cross-section weight. Each slice produces its own caloShowerShape_dijet_JZ<N>
ntuple; the analysis macros chain the ten of them with a glob.

Condor job event logs are written under ~/condor_logs/<label>/ (local
filesystem, required by the UChicago AF policy).

Usage examples:

  # Displaced dark photon signal — with the standard CaloShowerShapeTriggers
  # data layout, --signal is all you need (paths are auto-detected):
  python3 submit_caloShowerShape.py --signal displaced_dark_photon

  # Override any auto-detected path explicitly:
  python3 submit_caloShowerShape.py \\
      --signal displaced_dark_photon \\
      --daod-dir /data/larsonma/CaloShowerShapeTriggers/JETM42_DAODs/DisplacedDarkPhoton/<container> \\
      --gep-dir  /data/larsonma/CaloShowerShapeTriggers/GEPOutputReaderNTuples/DisplacedDarkPhoton/<container> \\
      --output-dir /data/larsonma/CaloShowerShapeTriggers/ntuples

  # QCD dijet background, JZ slice 4 — paths auto-detected from the
  # GEPHadronicEventReconstruction QCD_Dijet layout:
  python3 submit_caloShowerShape.py --background --jz 4

  # Dry run — prints the .sub file without submitting:
  python3 submit_caloShowerShape.py --signal displaced_dark_photon --dry-run
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

WRAPPER = Path(__file__).parent / "run_caloShowerShape_job.sh"

SIGNAL_STRINGS = {
    "displaced_dark_photon",
    "emerging_jets",
}

# Default CaloShowerShapeTriggers data layout, used to auto-fill --daod-dir /
# --gep-dir / --output-dir when they are omitted for a --signal sample.
CALO_BASE = "/data/larsonma/CaloShowerShapeTriggers"
DEFAULT_DAOD_BASE = f"{CALO_BASE}/JETM42_DAODs"
DEFAULT_GEP_BASE = f"{CALO_BASE}/GEPOutputReaderNTuples"
DEFAULT_OUTPUT_DIR = f"{CALO_BASE}/ntuples"

# Signal key -> sample sub-directory name under the *_BASE paths above.
SIGNAL_DIRS = {
    "displaced_dark_photon": "DisplacedDarkPhoton",
    "emerging_jets":         "EmergingJets",
}

# QCD dijet is not part of the CaloShowerShapeTriggers download: its DAODs and
# GEP ntuples live in the main GEPHadronicEventReconstruction layout, one
# sub-directory per JZ slice (the same inputs submit_ntupler.py uses).
GEP_PROJECT_BASE = "/data/larsonma/GEPHadronicEventReconstruction"
DIJET_DAOD_BASE = f"{GEP_PROJECT_BASE}/JETM42_DAODs/QCD_Dijet"
DIJET_GEP_BASE = f"{GEP_PROJECT_BASE}/GEPOutputReaderNTuples/QCD_Dijet"
N_JZ_SLICES = 10

_JOB_NUM_RE = re.compile(r'\.?_([\d]{5,6})\.')


def extract_job_num(path: str) -> str | None:
    m = _JOB_NUM_RE.search(Path(path).name)
    return m.group(1) if m else None


def find_container(base_dir: str, ext: str) -> str:
    """Return the single user.mlarson*_<ext> Rucio container dir under base_dir.

    Used to auto-resolve a sample's DAOD (EXT1) / GEP (EXT0) container directory
    when --daod-dir / --gep-dir are not given explicitly.
    """
    matches = sorted(str(p) for p in Path(base_dir).glob(f"user.mlarson*_{ext}")
                     if p.is_dir())
    if not matches:
        sys.exit(f"[error] no {ext} container directory found under: {base_dir}")
    if len(matches) > 1:
        print(f"[warn] multiple {ext} containers under {base_dir}, using: {matches[0]}")
    return matches[0]


def find_daod_files(daod_dir: str) -> list[tuple[str, str]]:
    """Return sorted [(job_num, filepath)] for DAOD pool.root files."""
    results = []
    for p in Path(daod_dir).iterdir():
        if not p.is_file():
            continue
        name = p.name
        if not ("DAOD_JETM42" in name and name.endswith(".root") and "pool" in name):
            continue
        num = extract_job_num(str(p))
        if num is None:
            print(f"  [warn] no job number in: {name} — skipping")
            continue
        results.append((num, str(p)))
    results.sort(key=lambda x: x[0])
    return results


def find_gep_files(gep_dir: str) -> dict[str, str]:
    """Return {job_num: filepath} for ntuple *.root files in gep_dir.
    Files without a job number are stored under key ''.
    """
    result: dict[str, str] = {}
    for p in Path(gep_dir).iterdir():
        if not p.is_file():
            continue
        if not p.name.endswith(".root"):
            continue
        num = extract_job_num(str(p))
        key = num if num is not None else ""
        result[key] = str(p)
    return result


def make_submit_file(jobs: list[dict], wrapper: str, label: str,
                     itemdata_path: str) -> tuple[str, str]:
    """Return (submit_text, itemdata_text).

    Fixed arguments (same for all jobs: signal flag, tag, output dir, JZ slice,
    pileup) go in the submit header. Only the per-file parts (DAOD, GEP, SUFFIX,
    LOG, OUT, ERR) go in the comma-separated itemdata file.

    Background jobs pass "BACKGROUND" as signal_string since the empty-string
    token would shift all positional bash arguments (the wrapper maps it to "").
    """
    j0 = jobs[0]
    sig_str = j0['signal_string'] if j0['signal_string'] else "BACKGROUND"
    fixed_args = (f"{j0['signal_bool']} {sig_str} {j0['output_dir']} "
                  f"{j0['jz_slice']} {j0['pileup']}")

    submit_lines = [
        f"# Auto-generated by submit_caloShowerShape.py for: {label}",
        "",
        f"executable            = {wrapper}",
        "universe              = vanilla",
        "request_cpus          = 1",
        "request_memory        = 4096",
        "request_disk          = 2048",
        '+queue                = "short"',
        "use_x509userproxy     = true",
        "x509userproxy         = /home/larsonma/x509proxy",
        "getenv                = false",
        "",
        f"arguments = {fixed_args} $(DAOD) $(GEP) $(SUFFIX)",
        "log       = $(LOG)",
        "output    = $(OUT)",
        "error     = $(ERR)",
        "",
        f"queue DAOD, GEP, SUFFIX, LOG, OUT, ERR from {itemdata_path}",
    ]

    data_lines = []
    for j in jobs:
        data_lines.append(
            f"{j['daod_file']}, {j['gep_file']}, {j['file_suffix']}, "
            f"{j['log']}, {j['stdout']}, {j['stderr']}"
        )

    return "\n".join(submit_lines), "\n".join(data_lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--signal", metavar="STR",
                       help="Signal sample key (e.g. displaced_dark_photon)")
    group.add_argument("--background", action="store_true",
                       help="Process as background/dijet (signal flag = 0)")

    parser.add_argument("--label", metavar="STR",
                        help="Override the sample label used for output naming / "
                             "condor logs (e.g. dijet_JZ4). Required to keep multiple "
                             "background samples from colliding; also becomes the "
                             "output-file tag. Defaults to dijet_JZ<N> when --jz is given.")
    parser.add_argument("--jz", type=int, default=None, metavar="N",
                        help=f"QCD dijet JZ slice (0-{N_JZ_SLICES - 1}) for --background. The "
                             "ntupler needs it to weight each event by the slice cross section, "
                             "and it selects the default input directories.")
    parser.add_argument("--pu", type=int, default=200, choices=(140, 200), metavar="N",
                        help="Pileup scenario selecting the sum-of-weights table (default: 200)")
    parser.add_argument("--daod-dir", metavar="DIR",
                        help="Directory containing DAOD_JETM42 *.pool.root files "
                             "(default: auto-detected EXT1 container for --signal under "
                             f"{DEFAULT_DAOD_BASE}/<Sample>)")
    parser.add_argument("--gep-dir", metavar="DIR",
                        help="Directory containing GEPOutputReader ntuple *.root files "
                             "(default: auto-detected EXT0 container for --signal under "
                             f"{DEFAULT_GEP_BASE}/<Sample>)")
    parser.add_argument("--output-dir", metavar="DIR",
                        help="Output directory for per-job ntuple ROOT files "
                             f"(default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("--log-dir", metavar="DIR",
                        help="Condor log/out/err directory (default: ~/condor_logs/<label>)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the .sub file but do not submit")
    parser.add_argument("--max-jobs", type=int, default=0, metavar="N",
                        help="Limit to first N job pairs (0 = no limit, useful for testing)")
    args = parser.parse_args()

    if args.signal and args.signal not in SIGNAL_STRINGS:
        parser.error(f"unknown signal string '{args.signal}'\n  valid: {sorted(SIGNAL_STRINGS)}")
    if args.jz is not None:
        if not args.background:
            parser.error("--jz applies to --background only")
        if not 0 <= args.jz < N_JZ_SLICES:
            parser.error(f"--jz must be in [0, {N_JZ_SLICES - 1}]")
    if args.background and args.jz is None:
        print("[warn] --background without --jz: the ntuple gets sampleJZSlice = -1 and "
              "mcEventWeight as its weight, so it cannot be mixed with other slices")

    signal_bool = 0 if args.background else 1
    # jz_slice reaches the macro as its jzSlice argument (-1 = no slice / signal).
    jz_slice = args.jz if args.jz is not None else -1
    # label: used for output-file tag, condor logs, and submit filenames.
    # For background/dijet samples --label (or --jz) distinguishes JZ slices.
    if args.background:
        label = args.label or (f"dijet_JZ{args.jz}" if args.jz is not None else "background")
    else:
        label = args.label if args.label else args.signal
    # signal_string is the tag the macro stamps on its output ROOT file. For a
    # real signal it is the signal key; for background it is the label, i.e.
    # dijet_JZ<N> for a slice (empty -> the wrapper/macro fall back to the
    # "background" tag).
    if args.background:
        signal_string = "" if label == "background" else label
    else:
        signal_string = args.signal

    # Fill in any omitted paths. Signals come from the CaloShowerShapeTriggers
    # layout; a --jz background slice comes from the QCD_Dijet/JZ<N> directories of
    # the main GEPHadronicEventReconstruction layout. Background without --jz has no
    # known layout and must be pointed at its inputs explicitly.
    if not args.daod_dir or not args.gep_dir:
        if args.background and args.jz is None:
            parser.error("--background without --jz requires explicit --daod-dir and --gep-dir")
        if args.background:
            if not args.daod_dir:
                args.daod_dir = find_container(os.path.join(DIJET_DAOD_BASE, f"JZ{args.jz}"), "EXT1")
            if not args.gep_dir:
                args.gep_dir = find_container(os.path.join(DIJET_GEP_BASE, f"JZ{args.jz}"), "EXT0")
        else:
            sample_dir = SIGNAL_DIRS.get(args.signal)
            if sample_dir is None:
                parser.error(f"no default data layout for signal '{args.signal}'; "
                             "pass --daod-dir/--gep-dir explicitly")
            if not args.daod_dir:
                args.daod_dir = find_container(os.path.join(DEFAULT_DAOD_BASE, sample_dir), "EXT1")
            if not args.gep_dir:
                args.gep_dir = find_container(os.path.join(DEFAULT_GEP_BASE, sample_dir), "EXT0")
    if not args.output_dir:
        args.output_dir = DEFAULT_OUTPUT_DIR

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Condor job event logs MUST be on a local filesystem (not Ceph /data).
    if args.log_dir:
        log_dir = os.path.abspath(args.log_dir)
    else:
        log_dir = str(Path.home() / "condor_logs" / f"caloShowerShape_{label}")
    os.makedirs(log_dir, exist_ok=True)

    print(f"Sample: {label}  (signal={signal_bool}, jzSlice={jz_slice}, PU{args.pu})")
    print(f"  DAOD dir: {args.daod_dir}")
    print(f"  GEP  dir: {args.gep_dir}")

    # Find DAOD files
    daod_files = find_daod_files(args.daod_dir)
    if not daod_files:
        sys.exit(f"[error] no DAOD pool.root files found in: {args.daod_dir}")
    print(f"Found {len(daod_files)} DAOD file(s)")

    # Find GEP files
    gep_map = find_gep_files(args.gep_dir)
    if not gep_map:
        sys.exit(f"[error] no GEP*.root files found in: {args.gep_dir}")

    single_gep: str | None = None
    if len(gep_map) == 1 and "" in gep_map:
        single_gep = gep_map[""]
        print(f"Found 1 merged GEP file (used for all jobs): {single_gep}")
    else:
        print(f"Found {len(gep_map)} numbered GEP file(s) — pairing by job number")

    # Pair DAOD with GEP files
    jobs: list[dict] = []
    skipped = 0
    for num, daod_path in daod_files:
        if single_gep:
            gep_path = single_gep
        else:
            gep_path = gep_map.get(num)
            if gep_path is None:
                print(f"  [warn] no GEP match for job number _{num} ({Path(daod_path).name}) — skipping")
                skipped += 1
                continue

        safe = f"caloShowerShape_{label}_{num}"
        jobs.append({
            "signal_bool":   signal_bool,
            "signal_string": signal_string,
            "output_dir":    output_dir,
            "jz_slice":      jz_slice,
            "pileup":        args.pu,
            "daod_file":     daod_path,
            "gep_file":      gep_path,
            "file_suffix":   f"_{num}",
            "log":    os.path.join(log_dir, f"{safe}.log"),
            "stdout": os.path.join(log_dir, f"{safe}.out"),
            "stderr": os.path.join(log_dir, f"{safe}.err"),
        })

    if skipped:
        print(f"[warn] {skipped} DAOD file(s) skipped — no matching GEP file by job number")
    if not jobs:
        sys.exit("[error] no jobs to submit after pairing")

    if args.max_jobs > 0:
        jobs = jobs[:args.max_jobs]
        print(f"Limiting to first {args.max_jobs} job(s) (--max-jobs)")

    # Write submit file + itemdata file
    itemdata_path = os.path.join(log_dir, f"caloShowerShape_{label}.dat")
    submit_text, itemdata_text = make_submit_file(jobs, str(WRAPPER), label, itemdata_path)

    sub_path = os.path.join(log_dir, f"caloShowerShape_{label}.sub")
    with open(sub_path, "w") as fh:
        fh.write(submit_text)
    with open(itemdata_path, "w") as fh:
        fh.write(itemdata_text)
    print(f"\nWrote submit file ({len(jobs)} job(s)): {sub_path}")
    print(f"Log directory: {log_dir}")

    if args.dry_run:
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
