#!/usr/bin/env python3
"""
Discover, submit, and select largeRJet compute/plot configurations.

DISCOVERY instead of a hand-typed grid: every merged background jet-tagger output in
TAGGER_DIR is one algorithm configuration (the filename encodes it — the same tags
analysisHelperFunctions.h ParseFileName decodes). Each configuration is paired with
every SIGNAL_STRINGS sample whose matching signal tagger output exists. Compute
therefore covers "everything reasonable that exists on disk" and never needs editing
when new emulation outputs appear; plotting later is a --select away.

State files are keyed by the background tagger basename and their contents carry no
pair indices, so:
  * compute jobs are submitted ONCE per distinct background configuration (all signal
    samples share it — signal loops run in the plot stage anyway);
  * any plot-time subset/order of pairs works against previously computed state.

Usage:
  python3 submit_largeRJetCompute.py [--dry-run] [--pairs-only] [--pu 140]
                                     [--select SUBSTR ...] [--n-jobs N] [--label STR]

  --select S    keep only pairs whose description or filename contains S. Repeatable;
                entries are AND-ed, and '|' inside ONE entry gives OR-ed alternatives
                (e.g. --select "_IOs_128_|_IOs_256_"). Same rule as lrj_plot_select_
                in the plot macro. Applies to the written pair list and, without
                --pairs-only, to what is submitted.
  --pairs-only  write ../analysis/largeRJetFilePairs.txt for a plot run, submit
                nothing; warns for selected pairs whose state file is missing.
  --dry-run     print everything, write and submit nothing
  --pu N        200 (r16130, default), 140 (r16129), or "both" = one list holding
                BOTH pileups, which is what a PU140-vs-PU200 overlay needs
  --n-jobs N    background slices per config (default 1; >1 unproven, see splitPlan 5b)
  --label STR   study label: condor log dir + overlay output dir name

Typical flows:
  compute everything available:   python3 submit_largeRJetCompute.py
  plot a subset later:            python3 submit_largeRJetCompute.py --pairs-only \
                                      --select ggF --select rMerge_0.001 --select IOs_128 \
                                      --label seededConeOnly
                                  cd ../analysis && root -b -l -q 'largeRJetPlot.C'
"""

import argparse
import glob
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

WRAPPER   = Path(__file__).parent / "run_largeRJetCompute_job.sh"
PAIRS_TXT = Path(__file__).parent.parent / "analysis" / "largeRJetFilePairs.txt"

# ---------------------------------------------------------------------------
# What to look at — these rarely change.
# ---------------------------------------------------------------------------
# Jet-tagger outputs. ONE directory holds every production: the 2 GeV input-tower E_T cut
# (jetTaggerEmulation.cc `apply_input_tower_et_cut_`) is recorded in the FILENAME as _T2_
# (cut applied) or _T0_ (no cut), not in the path, so both are discovered together and can
# be overlaid. Select one with --select _T0_ / --select _T2_. Override the directory
# per-invocation with --tagger-dir, or by name with --production.
# Each production writes the SAME filenames for the same configuration, so a production is
# only identified by its directory — and state files are keyed on the tagger BASENAME alone.
# Every entry therefore needs its own state directory in STATE_DIRS below, and the two are
# switched together via --production.
TAGGER_DIRS = {
    # Current: the input-tower E_T cut is a per-job parameter and is applied consistently at
    # every PU-suppression branch (see default_input_tower_et_cut_gev_ in jetTaggerEmulation.cc).
    "towerEtCutFixed": "/data/larsonma/LargeRadiusJets/outputNTuplesDev_gjTowerSamples_towerEtCutFixed",
    # Previous: _T2_ meant a 4 GeV cut on the EtaSK arm and 2 GeV on SK/NoSK. Kept so the old
    # plots remain reproducible; do not mix its files with the ones above in a single study.
    "all": "/data/larsonma/LargeRadiusJets/outputNTuplesDev_gjTowerSamples",
}
STATE_DIRS = {
    "towerEtCutFixed": "/data/larsonma/LargeRadiusJets/largeRJetIntermediateState_towerEtCutFixed",
    "all":             "/data/larsonma/LargeRadiusJets/largeRJetIntermediateState",
}
DEFAULT_PRODUCTION = "towerEtCutFixed"
TAGGER_DIR = TAGGER_DIRS[DEFAULT_PRODUCTION]
STATE_DIR  = Path(STATE_DIRS[DEFAULT_PRODUCTION])

# Signal samples to pair with each discovered configuration (only pairs whose signal
# tagger output exists are kept, so listing a sample here is free).
SIGNAL_STRINGS = ["ggF_hh_bbbb", "ttbar_allhad", "Zprime_ttbar_allhad_flatpT"]

# Discovered configurations matching any of these substrings are dropped entirely
# (e.g. one-off debug encodings). Checked against the config suffix.
EXCLUDE_CONFIGS = []

NTUPLE_BASE_BY_PU = {
    200: "/data/larsonma/GEPHadronicEventReconstruction/ntuples",
    140: "/data/larsonma/GEPHadronicEventReconstruction/ntuples_PU140",
}
RTAG_BY_PU = {200: "r16130", 140: "r16129"}

# signal string -> (ntuple subdirectory, mc tag up to and excluding the r-tag)
SAMPLES = {
    "ggF_hh_bbbb":                ("ggF_HHbbbb_v4",                "mc21_14TeV_HHbbbb_HLLHC_e8564_s4422"),
    "VBF_hh_bbbb":                ("VBF_HHbbbb_v4",                "mc21_14TeV_HHbbbb_HLLHC_VBF_e8557_s4422"),
    "ZvvHbb":                     ("ZvvHbb_v4",                    "mc21_14TeV_ZvvH125_bb_e8557_s4422"),
    "ttbar_allhad":               ("ttbar_allhad_v4",              "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422"),
    "Zprime_ttbar_allhad_flatpT": ("Zprime_ttbar_allhad_flatpT_v4","mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422"),
}
BACKGROUND_TAG = "mc21_14TeV_jj_JZ_e8557_s4422"

# ---------------------------------------------------------------------------


def discover_configs(pu: int) -> list[str]:
    """Config suffixes (everything after '<BACKGROUND_TAG>_<rtag>_', .root included)
    of every merged background tagger output at this pileup."""
    rtag = RTAG_BY_PU[pu]
    prefix = f"{BACKGROUND_TAG}_{rtag}_"
    suffixes = []
    for f in sorted(glob.glob(f"{TAGGER_DIR}/{prefix}*.root")):
        name = Path(f).name
        if re.search(r"_file\d+\.root$", name):
            continue                      # per-file emulation chunk, not a merged output
        suffix = name[len(prefix):]
        if any(x in suffix for x in EXCLUDE_CONFIGS):
            continue
        suffixes.append(suffix)
    return suffixes


def describe(suffix: str) -> str:
    """Compact human/select-friendly description parsed from the config suffix."""
    fields = []
    for pat, label in [(r"rMerge_([\d.]+)", "rMerge={}"), (r"IOs_(\d+)", "IOs={}"),
                       (r"_IO_([A-Za-z]+)_Seed", "IO={}"), (r"_Seed_([A-Za-z]+)_", "Seed={}"),
                       (r"_(EtaSK|SK|NoSK)_", "{}"), (r"_T(\d+)_subjetEt", "T{}"),
                       (r"subjetEt(\d+)GeV", "subjetEt={}"),
                       (r"_ewm(\d)_", "ewm{}"), (r"_mep(\d)_", "mep{}"),
                       (r"_mec(\d+)GeV", "mec={}"), (r"_v(\d)\.root", "v{}")]:
        m = re.search(pat, suffix)
        if m:
            fields.append(label.format(m.group(1)))
    return " ".join(fields) if fields else suffix


def build_pairs(pu: int):
    """(desc, sigNtuple, sigTagger, bkgNtuple, bkgTagger) per (config, signal) with an
    existing signal tagger; deterministic order: config suffix, then SIGNAL_STRINGS."""
    base, rtag = NTUPLE_BASE_BY_PU[pu], RTAG_BY_PU[pu]
    pairs, missing_sig = [], []
    for suffix in discover_configs(pu):
        bkg_tagger = f"{TAGGER_DIR}/{BACKGROUND_TAG}_{rtag}_{suffix}"
        bkg_ntuple = (f"{base}/QCD_Dijet_JZ*_v4/"
                      f"mc21_14TeV_jj_JZ*_e8557_s4422_{rtag}_DAOD_NTUPLE_GEP.root")
        for sig in SIGNAL_STRINGS:
            subdir, mc_tag = SAMPLES[sig]
            sig_tagger = f"{TAGGER_DIR}/{mc_tag}_{rtag}_{suffix}"
            if not Path(sig_tagger).exists():
                missing_sig.append(f"{sig}: {Path(sig_tagger).name}")
                continue
            sig_ntuple = f"{base}/{subdir}/{mc_tag}_{rtag}_DAOD_NTUPLE_GEP.root"
            desc = f"sample={sig} PU{pu} {describe(suffix)}"
            pairs.append((desc, sig_ntuple, sig_tagger, bkg_ntuple, bkg_tagger))
    return pairs, missing_sig


def select_pairs(pairs, selects):
    """Entries are AND-ed; alternatives within one entry are OR-ed on '|'. Same rule as
    lrj_plot_select_ in largeRJetAnalysisCore.h — the two selectors must agree, since a
    study is normally carved here and then mirrored in largeRJetPlot.C."""
    if not selects:
        return pairs
    def keep(p):
        haystack = p[0] + " " + p[2] + " " + p[4]
        return all(any(alt and alt in haystack for alt in s.split("|")) for s in selects)
    return [p for p in pairs if keep(p)]


def write_pairs_file(pairs, pu, label, dry_run, dest: Path = None):
    lines = [
        f"# Generated by submit_largeRJetCompute.py on {datetime.now():%Y-%m-%d %H:%M:%S}",
        f"# --pu {pu}  label: {label}",
        "# Read by largeRJetCompute.C / largeRJetPlot.C. Regenerate (with --pairs-only",
        "# and --select filters) rather than hand-editing.",
        f"overlaydir overlayMultipleFiles/largeRJet_{label}/",
    ]
    for i, (desc, sn, st, bn, bt) in enumerate(pairs):
        lines.append(f"# pair {i}: {desc}")
        lines.append(f"pair {sn} {st} {bn} {bt}")
    text = "\n".join(lines) + "\n"
    out = dest or PAIRS_TXT
    if dry_run:
        print(f"--- {out} (dry run, not written) ---\n{text}")
    else:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
        print(f"Wrote {len(pairs)} pair(s) to {out}")


def make_submit_file(compute_jobs, n_jobs, label, log_dir, pairs_path):
    lines = [
        f"# Auto-generated by submit_largeRJetCompute.py for: {label}",
        "",
        f"executable            = {WRAPPER}",
        "universe              = vanilla",
        "request_cpus          = 1",
        "request_memory        = 8192",
        "request_disk          = 2048",
        '+queue                = "short"',
        "getenv                = false",
        "",
        # Per-submission SNAPSHOT, not the shared analysis/largeRJetFilePairs.txt: jobs read
        # this path when they START, so a second submission (or a --pairs-only plot
        # selection) overwriting the shared file would silently repoint jobs already in the
        # queue at a different pair list — which is exactly what happened on 2026-08-21.
        f"arguments = $(PAIR_INDEX) $(JOB_INDEX) {n_jobs} {pairs_path}",
        f"log       = {log_dir}/job_$(PAIR_INDEX)_$(JOB_INDEX).log",
        f"output    = {log_dir}/job_$(PAIR_INDEX)_$(JOB_INDEX).out",
        f"error     = {log_dir}/job_$(PAIR_INDEX)_$(JOB_INDEX).err",
        "",
        "queue PAIR_INDEX, JOB_INDEX from (",
    ]
    for pair_index in compute_jobs:
        for j in range(n_jobs):
            lines.append(f"  {pair_index}, {j}")
    lines.append(")")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run",    action="store_true")
    ap.add_argument("--pairs-only", action="store_true")
    ap.add_argument("--select", action="append", default=[], metavar="SUBSTR")
    # "both" puts PU140 AND PU200 pairs in ONE list, which is what a cross-pileup
    # overlay needs: the plot stage labels mixed runs "<PU> = 140, 200" and appends
    # "<PU> = 140"/"= 200" to each legend entry (see overlayMixesPileup in the core).
    ap.add_argument("--pu", default="200", choices=("200", "140", "both"))
    # Default 1 = one job per configuration, over the full background. Chunking
    # (--n-jobs > 1) was originally suspected of causing the 2026-08-21 bad_alloc
    # failures; that was wrong (the pair-skip UB was — see splitPlan 5b, now fixed).
    # Chunking is simply unvalidated: leave at 1 until a clean run is diffed against
    # the monolith.
    ap.add_argument("--n-jobs", type=int, default=1,
                    help="background slices per config; >1 is UNSAFE, see splitPlan 5b")
    ap.add_argument("--label",  default=None)
    # Which jet-tagger production to discover configurations from. --production names one of
    # the known ones; --tagger-dir takes an arbitrary path and wins if both are given.
    ap.add_argument("--production", choices=sorted(TAGGER_DIRS), default=None,
                    help=f"named tagger production (default: the TAGGER_DIR set in this file)")
    ap.add_argument("--tagger-dir", default=None, metavar="DIR",
                    help="explicit tagger output directory; overrides --production")
    args = ap.parse_args()

    global TAGGER_DIR, STATE_DIR
    if args.tagger_dir:
        TAGGER_DIR = args.tagger_dir.rstrip("/")
        print(f"[state] --tagger-dir given; state directory stays {STATE_DIR}."
              " Pass --production instead if this is a separate production, or its chunks will"
              " land on top of that one's.")
    elif args.production:
        TAGGER_DIR = TAGGER_DIRS[args.production]
        # State is keyed on the tagger basename, which is identical across productions, so the
        # state directory has to follow or this production would read the other one's results.
        STATE_DIR = Path(STATE_DIRS[args.production])
        print(f"[state] using {STATE_DIR}")
    if not Path(TAGGER_DIR).is_dir():
        print(f"[error] tagger directory does not exist: {TAGGER_DIR}")
        sys.exit(1)
    print(f"[tagger] discovering configurations from {TAGGER_DIR}")

    label = args.label or f"{datetime.now():%Y%m%d_%H%M%S}"
    pu_list = [140, 200] if args.pu == "both" else [int(args.pu)]
    pairs, missing_sig = [], []
    for _pu in pu_list:
        _p, _m = build_pairs(_pu)
        pairs += _p; missing_sig += _m
    if missing_sig:
        print(f"[info] {len(missing_sig)} (config, signal) combination(s) skipped —"
              " no signal tagger output on disk. First few:")
        for m in missing_sig[:5]:
            print(f"       {m}")
    n_before = len(pairs)
    pairs = select_pairs(pairs, args.select)
    if args.select:
        print(f"--select {args.select} kept {len(pairs)} of {n_before} pair(s)")
    if not pairs:
        print("No pairs survive discovery + --select filters."); sys.exit(1)

    print(f"=== {len(pairs)} pair(s) discovered, PU{args.pu} ===")
    for i, (desc, *_r) in enumerate(pairs):
        print(f"  pair {i}: {desc}")

    write_pairs_file(pairs, args.pu, label, args.dry_run)

    if args.pairs_only:
        # Plot-run generation: check the state each selected pair needs is in place.
        missing = []
        for desc, _sn, _st, _bn, bt in pairs:
            state = STATE_DIR / (Path(bt).name[:-5] + "__state.root")
            if not state.exists():
                missing.append(str(state))
        if missing:
            print(f"\n[WARN] {len(missing)} selected pair(s) have NO merged state file yet"
                  " — run compute (+ hadd) for them first:")
            for m in sorted(set(missing)):
                print(f"       {m}")
        else:
            print("All selected pairs have merged state files.")
        print("Plot with: cd ../analysis && root -b -l -q 'largeRJetPlot.C'")
        return

    # Compute submission: one job set per DISTINCT background tagger (state files are
    # shared across the signal samples of a configuration).
    seen, compute_jobs = set(), []
    for i, (_d, _sn, _st, _bn, bt) in enumerate(pairs):
        if bt not in seen:
            seen.add(bt)
            compute_jobs.append(i)
    print(f"\n{len(compute_jobs)} distinct background configuration(s) to compute"
          f" ({len(pairs)} pairs share them); {args.n_jobs} slice(s) each"
          f" = {len(compute_jobs) * args.n_jobs} jobs")
    if args.n_jobs > 1:
        print("[WARN] --n-jobs > 1 splits the background across jobs; each job then runs"
              " the rate-vs-eff\n       derivations on a partial background. Unvalidated —"
              " see ../analysis/largeRJetSplitPlan.md 5b.")

    log_dir = Path.home() / "condor_logs" / f"largeRJetCompute_{label}"
    # Snapshot the pair list beside the submit file and point the jobs at THAT copy, so a
    # later submission or plot selection cannot repoint queued jobs (see make_submit_file).
    pairs_snapshot = log_dir / f"largeRJetFilePairs_{label}.txt"
    sub_text = make_submit_file(compute_jobs, args.n_jobs, label, log_dir, pairs_snapshot)
    if args.dry_run:
        print(f"--- pairs snapshot would go to {pairs_snapshot} ---")
        print(f"--- submit file (dry run) ---\n{sub_text}")
        return
    log_dir.mkdir(parents=True, exist_ok=True)
    write_pairs_file(pairs, args.pu, label, False, dest=pairs_snapshot)
    sub_path = log_dir / f"largeRJetCompute_{label}.sub"
    sub_path.write_text(sub_text)
    subprocess.run(["condor_submit", str(sub_path)], check=True)
    print("\nMonitor: condor_watch_q\nThen:    ./hadd_largeRJetState.sh"
          "\nThen:    python3 submit_largeRJetCompute.py --pairs-only --select ..."
          " && cd ../analysis && root -b -l -q 'largeRJetPlot.C'")


if __name__ == "__main__":
    main()
