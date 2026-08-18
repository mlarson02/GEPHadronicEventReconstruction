#!/usr/bin/env python3
"""
Generate and submit HERNTupler Condor jobs for all DAOD files in a sample directory.

Each job processes one DAOD_JETM42 pool.root file paired with its matching
TrigGepPerf GEP ntuple file.  File pairing is done by the 6-digit job number
embedded in both filenames as  _NNNNNN  (e.g. _000001).  When only one GEP file
exists in the GEP directory (no job number), it is used for all jobs.

Condor job event logs are written under ~/condor_logs/<label>/ (local filesystem,
required by the UChicago AF policy — Ceph /data is not allowed for log files).
Job stdout/stderr go to the same directory unless --log-dir is specified.

Usage examples:

  # Signal sample (ggF HH->bbbb resim, algoVersion=3):
  python submit_ntupler.py --signal ggF_hh_bbbb_resim --algo 3 --daod-dir /data/larsonma/GEPHadronicEventReconstruction/JETM42_DAODs/ggF_HHbbbb/user.mlarson.GEPNtupleJETM42.ggF_HHbbbb.v2_TGP_JETM42_ggF_realsubmit_lessmemory_EXT1 --gep-dir /data/larsonma/GEPHadronicEventReconstruction/GEPOutputReaderNTuples/ggF_HHbbbb/user.mlarson.GEPNtupleJETM42.ggF_HHbbbb.v2_TGP_JETM42_ggF_realsubmit_lessmemory_EXT0 --output-dir /data/larsonma/GEPHadronicEventReconstruction/DAOD_TrigGepPerf/ggF_HHbbbb

  # Background JZ3 (algoVersion=3):
  python3 submit_ntupler.py \\
      --background --jz 3  --algo 3 \\
      --daod-dir /data/larsonma/GEPHadronicEventReconstruction/JETM42_DAODs/QCD_Dijet/JZ3/<container>/ \\
      --gep-dir  /data/larsonma/GEPHadronicEventReconstruction/GEPOutputReaderNTuples/QCD_Dijet/JZ3/<container>/ \\
      --output-dir /data/larsonma/GEPHadronicEventReconstruction/ntuples/QCD_Dijet_JZ3_v3/

  # PU140 background JZ3 — point dirs at the _PU140 mirrors and pass --pu 140:
  python3 submit_ntupler.py \\
      --background --jz 3 --algo 3 --pu 140 \\
      --daod-dir /data/larsonma/GEPHadronicEventReconstruction/JETM42_DAODs_PU140/QCD_Dijet/JZ3/<container>/ \\
      --gep-dir  /data/larsonma/GEPHadronicEventReconstruction/GEPOutputReaderNTuples_PU140/QCD_Dijet/JZ3/<container>/ \\
      --output-dir /data/larsonma/GEPHadronicEventReconstruction/ntuples_PU140/QCD_Dijet_JZ3_v3/

  # Dry run — prints the .sub file without submitting:
  python3 submit_ntupler.py ... --dry-run
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

WRAPPER = Path(__file__).parent / "run_ntupler_job.sh"

SIGNAL_STRINGS = {
    "ggF_hh_bbbb",
    "ZvvHbb",
    "VBF_HHbbbb",
    "ttbar_allhad",
    "ttbar_semilep",
    "ttbar_dilep",
    "Zprime_ttbar_allhad_flatpT",
    # Displaced dark photon: mc21_14TeV.543785.MGPy8EG_A14NN23LO_HAHM_ZdZd4e_110_30_10ns.recon.AOD.e8557_s4422_r16130
    #"displaced_dark_photon",
    # Emerging jets: mc21_14TeV.801966.Py8EG_Zprime2EJs_Ld20_rho40_pi10_Zp1500_l50.recon.AOD.e8532_s4422_r16130
    #"emerging_jets",
    # Z -> mu mu: mc21_14TeV.601190.PhPy8EG_AZNLO_Zmumu.recon.AOD.e8557_s4422_r16130
    "Zmumu",
    # Stau-Stau LLP (10 ns): mc21_14TeV.516672.MGPy8EG_A14NNPDF23LO_StauStauLLP_500_0_10ns.recon.AOD.e8557_s4422_r16130
    #"StauStauLLP_500_0_10ns",
}

_JOB_NUM_RE = re.compile(r'\.?_([\d]{5,6})\.')


def extract_job_num(path: str) -> str | None:
    m = _JOB_NUM_RE.search(Path(path).name)
    return m.group(1) if m else None


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

    Fixed arguments (same for all jobs) go in the submit header.
    Only the per-file parts (DAOD, GEP, SUFFIX, LOG, OUT, ERR) go in the
    comma-separated itemdata file — none of those values contain spaces or
    commas, so Condor parses them correctly.

    Background jobs pass "BACKGROUND" as signal_string since the empty-string
    token would shift all positional bash arguments.
    """
    j0 = jobs[0]
    sig_str = j0['signal_string'] if j0['signal_string'] else "BACKGROUND"
    fixed_args = (f"{j0['signal_bool']} {sig_str} {j0['algo_version']} "
                  f"{j0['jz_slice']} {j0['special_jz0']} {j0['output_dir']}")
    # pileup is fixed for the whole submission; it is the trailing ($10) argument.
    pileup = j0['pileup']

    submit_lines = [
        f"# Auto-generated by submit_ntupler.py for: {label}",
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
        f"arguments = {fixed_args} $(DAOD) $(GEP) $(SUFFIX) {pileup}",
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
                       help="Signal sample key (e.g. ggF_hh_bbbb_resim)")
    group.add_argument("--background", action="store_true",
                       help="Process as QCD dijet background")

    parser.add_argument("--jz", type=int, default=3, metavar="N",
                        help="JZ slice (0-9) when using --background (default: 3)")
    parser.add_argument("--special-jz0", action="store_true",
                        help="specialJZ0=true: the noHSTP JZ0 background variant")
    parser.add_argument("--algo", type=int, required=True, choices=[2, 3, 4], metavar="VERSION",
                        help="Algorithm version (2 or 3)")
    parser.add_argument("--pu", type=int, default=200, choices=[140, 200], metavar="PILEUP",
                        help="Pileup scenario; selects the reweighting constants in HERNTupler "
                             "(default: 200). Point --daod-dir/--gep-dir/--output-dir at the "
                             "_PU140 mirrors when using --pu 140.")
    parser.add_argument("--daod-dir", required=True, metavar="DIR",
                        help="Directory containing DAOD_JETM42 *.pool.root files")
    parser.add_argument("--gep-dir", required=True, metavar="DIR",
                        help="Directory containing GEP ntuple *.root files")
    parser.add_argument("--output-dir", required=True, metavar="DIR",
                        help="Output directory for per-job ntuple ROOT files")
    parser.add_argument("--log-dir", metavar="DIR",
                        help="Condor log/out/err directory (default: ~/condor_logs/<label>)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the .sub file but do not submit")
    parser.add_argument("--max-jobs", type=int, default=0, metavar="N",
                        help="Limit to first N job pairs (0 = no limit, useful for testing)")
    args = parser.parse_args()

    if args.signal and args.signal not in SIGNAL_STRINGS:
        parser.error(f"unknown signal string '{args.signal}'\n  valid: {sorted(SIGNAL_STRINGS)}")

    signal_bool   = 0 if args.background else 1
    signal_string = "" if args.background else args.signal
    jz_slice      = args.jz
    special_jz0   = 1 if args.special_jz0 else 0

    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    if args.background:
        label = f"background_JZ{jz_slice}_v{args.algo}"
        if args.special_jz0:
            label += "_noHSTP"
    else:
        label = f"{args.signal}_v{args.algo}"
    # Tag non-default pileup so PU140 logs/itemdata/sub files don't collide with PU200.
    if args.pu != 200:
        label += f"_PU{args.pu}"

    # Condor job event logs MUST be on a local filesystem (not Ceph /data).
    if args.log_dir:
        log_dir = os.path.abspath(args.log_dir)
    else:
        log_dir = str(Path.home() / "condor_logs" / label)
    os.makedirs(log_dir, exist_ok=True)

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

        safe = f"{label}_{num}"
        jobs.append({
            "signal_bool":   signal_bool,
            "signal_string": signal_string,
            "algo_version":  args.algo,
            "jz_slice":      jz_slice,
            "special_jz0":   special_jz0,
            "output_dir":    output_dir,
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
    itemdata_path = os.path.join(log_dir, f"ntupler_{label}.dat")
    submit_text, itemdata_text = make_submit_file(jobs, str(WRAPPER), label, itemdata_path)

    sub_path = os.path.join(log_dir, f"ntupler_{label}.sub")
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
