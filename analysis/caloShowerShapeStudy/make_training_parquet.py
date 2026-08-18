#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.10"
# dependencies = ["uproot>=5", "awkward>=2", "numpy>=1.24", "pyarrow>=14"]
# ///
"""
make_training_parquet.py
========================
Turn caloShowerShape ntuples into a *digitized* per-jet feature table (parquet)
for training an HLS4ML-implementable displaced-jet trigger network.

Ultimate goal of the study: a firmware ML trigger that predicts jet displacement
from the calorimeter's longitudinal (per-layer) shower profile, so a displaced
jet can pass at *lower* E_T than a prompt-jet threshold would allow (lower rate).
This script builds the training set that discriminates prompt QCD (JZ slices,
label 0) from displaced signal (displaced dark photon / emerging jets, label 1).

Digitization matches the GEP tower FW/emulator quantization (see constants below):
  * Et  : LSB = 0.125 GeV, 13-bit field, clamped to [0, 1024] GeV  (et code).
  * eta : 0.1 grid over [-4.85, 4.95], 7-bit field  (eta code, 0..97).
  * phi : 0.1 grid over [-3.2, 3.2],   6-bit field  (phi code, 0..63).

Per jet (one parquet row per jet):
  * Associate GEP EtaSK towers to the jet by dR < Rassoc (WTACone 0.4, LRJ 1.0).
  * DIGITIZE each tower's per-layer Et to integer counts (the value the firmware
    sums), then sum the counts of the in-cone towers per layer -> q_l0..q_l6.
    These are the FW-implementable NN inputs.
  * Also store the raw GeV per-layer sums (reference), the jet eta/phi/Et codes,
    the classification label, and a best-effort truth displacement target.

The full quantization is recorded in the parquet file-level metadata.

QCD background comes as one ntuple per JZ slice (caloShowerShape_dijet_JZ<N>.root,
JZ0-9), so with no arguments the two merged signals AND all ten dijet slices are
read. Each row carries its slice index (jz_slice) and the slice cross-section
weight the ntupler stamped on the event (jz_weight): the slices span thirteen
orders of magnitude in cross section, so any rate-like use of the background must
weight by jz_weight rather than counting rows. Background rows failing the HSTP
filter are dropped by default (--no-hstp-filter keeps them).

Run — with no arguments it processes signals + the JZ background and writes the
parquet to EOS:
  cd analysis/caloShowerShapeStudy
  python make_training_parquet.py                    # wtacone (default)
  python make_training_parquet.py --collection lrj
  python make_training_parquet.py --no-background    # signals only

  # or pass explicit inputs / globs (quote them so python does the expansion):
  python make_training_parquet.py 'caloShowerShape_dijet_JZ[0-9].root'

(Requires uproot + pyarrow in the active python env; `uv run make_training_parquet.py`
also works via the inline script metadata above.)
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import uproot
import pyarrow as pa
import pyarrow.parquet as pq

NLAYERS = 7
NTUPLE_DIR = "/data/larsonma/CaloShowerShapeTriggers/ntuples"   # merged input ntuples
OUT_DIR = "/eos/user/m/mlarson/CaloShowerDisplacement"          # parquet output (EOS)
# Default samples processed when no input paths are given: the two merged signals
# plus the ten per-slice QCD dijet ntuples (label 0).
DEFAULT_SIGNALS = ["displaced_dark_photon", "emerging_jets"]
# [0-9].root, not *.root: the per-job outputs (caloShowerShape_dijet_JZ9_000510.root)
# stay in this flat directory after the hadd, so a bare wildcard would read every
# event twice -- once per-job and once from the merged slice file.
DEFAULT_DIJET_GLOB = os.path.join(NTUPLE_DIR, "caloShowerShape_dijet_JZ[0-9].root")

# ---------------------------------------------------------------------------
# GEP tower FW quantization (mirrors the emulator constexpr constants):
#   et_granularity_=0.125  et_bit_length_=13   et_min_=0  et_max_=1024
#   eta_bit_length_=7  eta_range_=98  eta_min_=-4.85 eta_max_=4.95  eta_gran=0.1
#   phi_bit_length_=6  phi_min_=-3.2  phi_max_=3.2   phi_gran=0.1
#   (a tower packs as eta|et|phi|pad into 64 bits; pad = 64-7-13-6 = 38)
# ---------------------------------------------------------------------------
ET_GRANULARITY  = 0.125          # GeV / Et count
ET_BIT_LENGTH   = 13
ET_MIN_GEV      = 0.0
ET_MAX_GEV      = 1024.0
ETA_BIT_LENGTH  = 7
ETA_RANGE       = 98
ETA_MIN         = -4.85
ETA_MAX         = 4.95
ETA_GRANULARITY = 0.1
PHI_BIT_LENGTH  = 6
PHI_MIN         = -3.2
PHI_MAX         = 3.2
PHI_GRANULARITY = 0.1

ET_MAX_CODE  = min((1 << ET_BIT_LENGTH) - 1, int(round(ET_MAX_GEV / ET_GRANULARITY)))  # 8191
ETA_MAX_CODE = ETA_RANGE - 1                                                            # 97
PHI_MAX_CODE = (1 << PHI_BIT_LENGTH) - 1                                                # 63


def quant_et(et):
    """GeV -> integer Et counts (ADC-like truncation, clamped to [0, ET_MAX_CODE])."""
    c = np.floor(np.clip(et, ET_MIN_GEV, ET_MAX_GEV) / ET_GRANULARITY)
    return np.clip(c, 0, ET_MAX_CODE).astype(np.int64)


def quant_eta(eta):
    """eta -> integer grid code (nearest, clamped to [0, ETA_MAX_CODE])."""
    c = np.round((np.clip(eta, ETA_MIN, ETA_MAX) - ETA_MIN) / ETA_GRANULARITY)
    return np.clip(c, 0, ETA_MAX_CODE).astype(np.int64)


def quant_phi(phi):
    """phi -> integer grid code (nearest, clamped to [0, PHI_MAX_CODE])."""
    c = np.round((np.clip(phi, PHI_MIN, PHI_MAX) - PHI_MIN) / PHI_GRANULARITY)
    return np.clip(c, 0, PHI_MAX_CODE).astype(np.int64)


# Jet collection -> its GEP tree and the branch names / units used for kinematics.
# pt_scale converts the stored jet energy to GeV (WTACone Pt is already GeV; the
# LRJ Et branch is a GEP pass-through in MeV).
COLLECTIONS = {
    "wtacone": dict(tree="wtaConeCellsTowersEtaSKTree", pt="Pt", eta="Eta", phi="Phi",
                    pt_scale=1.0, rassoc=0.4),
    "lrj":     dict(tree="jetTaggerLRJEtaSKTree",       pt="Et", eta="Eta", phi="Phi",
                    pt_scale=1.0e-3, rassoc=1.0),
}
TOW_TREE = "gepCellsTowersEtaSKTree"
EVT_TREE = "eventInfoTree"
SP_TREE = "showerParentTree"   # shower-initiating LLP parents (decay vertex = label)


def wrap_dphi(d):
    """Wrap phi difference into (-pi, pi]."""
    return (d + np.pi) % (2.0 * np.pi) - np.pi


def sample_tag(path):
    """caloShowerShape_<tag>.root -> <tag> (falls back to the bare filename)."""
    b = os.path.basename(path)
    if b.startswith("caloShowerShape_") and b.endswith(".root"):
        return b[len("caloShowerShape_"):-len(".root")]
    return b[:-5] if b.endswith(".root") else b


def natural_key(path):
    """Sort key that orders digit runs numerically, so JZ10 follows JZ9."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", os.path.basename(path))]


def expand_inputs(patterns):
    """Expand globs (and plain paths) into an ordered, de-duplicated file list.

    The QCD background is ten per-slice files, so every input is treated as a
    possible glob; a literal path that matches nothing is reported, not silently
    dropped.
    """
    out = []
    for pat in patterns:
        matches = sorted(glob.glob(pat), key=natural_key)
        if not matches:
            print(f"  [warn] no file matches: {pat}")
            continue
        for m in matches:
            if m not in out:
                out.append(m)
    return out


def process_file(path, coll, et_min_tower, seed_cone, grid_hw=0, n_jets=2,
                 hstp_filter=True):
    """Return a dict of column-name -> python list of per-jet rows for one file.

    If grid_hw > 0, also emit per jet a fixed (side x side x NLAYERS) tower grid
    (tow_grid, side = 2*grid_hw+1): the digitized per-layer Et of the in-cone towers
    scattered into 0.1x0.1 (eta,phi) cells centered on the jet, row-major
    [dEta cell][dPhi cell][layer]. A fixed grid needs no mask -- empty cells are a
    physical zero. grid_hw=4 -> a 9x9 window (+/-0.4) that bounds the R=0.4 cone
    (<=45 of the 81 cells ever non-zero).

    n_jets limits processing to the leading-n jets by Et per event (0 = all);
    jet_index is then the Et rank (0 = leading, 1 = subleading).

    hstp_filter drops background events with passHSTP == False, matching the
    background selection the rate macros apply to the JZ chain.
    """
    cfg = COLLECTIONS[coll]
    tag = sample_tag(path)
    seed_cone = seed_cone if seed_cone > 0 else cfg["rassoc"]   # SEED->jet match cone

    f = uproot.open(path)
    tow = f[TOW_TREE].arrays(["Eta", "Phi"] + [f"Et_l{l}" for l in range(NLAYERS)], library="ak")
    jet = f[cfg["tree"]].arrays([cfg["pt"], cfg["eta"], cfg["phi"]], library="ak")
    # JZ bookkeeping is optional: ntuples made before the JZ-slice support have
    # none of these branches and simply get jz_slice=-1 / jz_weight=mcEventWeight.
    evt_keys = set(k.split(";")[0] for k in f[EVT_TREE].keys())
    ev = f[EVT_TREE].arrays(
        ["eventNumber", "runNumber", "mcChannelNumber", "mcEventWeight", "signal"], library="np")
    jz_slice = (f[EVT_TREE]["sampleJZSlice"].array(library="np")
                if "sampleJZSlice" in evt_keys else np.full(len(ev["eventNumber"]), -1))
    pass_hstp = (f[EVT_TREE]["passHSTP"].array(library="np").astype(bool)
                 if "passHSTP" in evt_keys else np.ones(len(ev["eventNumber"]), dtype=bool))
    if "eventWeights" in evt_keys:
        jzw = f[EVT_TREE]["eventWeights"].array(library="np")
        jz_weight = np.array([w[0] if len(w) else np.nan for w in jzw], dtype=np.float64)
    else:
        jz_weight = np.asarray(ev["mcEventWeight"], dtype=np.float64)
    n_hstp_dropped = 0
    have_sp = SP_TREE in f
    if have_sp:
        sp = f[SP_TREE].arrays(["pt", "eta", "phi", "decayVtx_Lxy", "decayVtx_r3d"], library="ak")
    else:
        print(f"    [warn] {SP_TREE} not in {os.path.basename(path)} -- regression labels will be "
              "NaN (re-run the ntupler to produce showerParentTree)")

    n = len(ev["eventNumber"])
    cols = {k: [] for k in (
        ["sample", "mcChannelNumber", "eventNumber", "runNumber", "event_weight",
         "jz_slice", "jz_weight", "pass_hstp",
         "collection", "jet_index", "jet_et_gev", "jet_eta", "jet_phi",
         "jet_et_code", "jet_eta_code", "jet_phi_code", "n_towers_in_cone"]
        + [f"q_l{l}" for l in range(NLAYERS)]
        + [f"et_l{l}_gev" for l in range(NLAYERS)]
        + ["q_total", "label", "truth_decay_Lxy_mm", "truth_decay_r3d_mm"]
        + (["tow_grid"] if grid_hw > 0 else [])
    )}

    for i in range(n):
        if hstp_filter and not ev["signal"][i] and not pass_hstp[i]:
            n_hstp_dropped += 1
            continue
        jeta = np.asarray(jet[cfg["eta"]][i], dtype=np.float64)
        if jeta.size == 0:
            continue
        jphi = np.asarray(jet[cfg["phi"]][i], dtype=np.float64)
        jpt = np.asarray(jet[cfg["pt"]][i], dtype=np.float64) * cfg["pt_scale"]
        jet_etc = quant_et(jpt)
        jet_etac = quant_eta(jeta)
        jet_phic = quant_phi(jphi)

        teta = np.asarray(tow["Eta"][i], dtype=np.float64)
        tphi = np.asarray(tow["Phi"][i], dtype=np.float64)
        # (ntow, NLAYERS) tower per-layer Et in GeV
        tetl = (np.stack([np.asarray(tow[f"Et_l{l}"][i], dtype=np.float64) for l in range(NLAYERS)], axis=1)
                if teta.size else np.zeros((0, NLAYERS)))
        if et_min_tower > 0.0:
            tetl = np.where(tetl > et_min_tower, tetl, 0.0)
        tow_qetl = quant_et(tetl)   # (ntow, NLAYERS) integer Et counts (FW digitization)

        if teta.size:
            deta = jeta[:, None] - teta[None, :]
            dphi = wrap_dphi(jphi[:, None] - tphi[None, :])
            in_cone = (deta * deta + dphi * dphi) < (cfg["rassoc"] ** 2)   # (njet, ntow) bool
            per_layer_counts = in_cone.astype(np.int64) @ tow_qetl         # (njet, NLAYERS)
            per_layer_gev = in_cone.astype(np.float64) @ np.where(tetl > 0, tetl, 0.0)
            n_in_cone = in_cone.sum(axis=1)
        else:
            per_layer_counts = np.zeros((jeta.size, NLAYERS), dtype=np.int64)
            per_layer_gev = np.zeros((jeta.size, NLAYERS), dtype=np.float64)
            n_in_cone = np.zeros(jeta.size, dtype=np.int64)

        # regression label: decay vertex of the highest-pt shower-parent (SEED)
        # whose momentum direction falls within the jet cone (dR < seed_cone).
        # Keep all radii (no cap). NaN if the jet has no in-cone SEED.
        Lxy = np.full(jeta.size, np.nan)
        R3d = np.full(jeta.size, np.nan)
        if have_sp:
            se = np.asarray(sp["eta"][i], dtype=np.float64)
            if se.size:
                sph = np.asarray(sp["phi"][i], dtype=np.float64)
                spt = np.asarray(sp["pt"][i], dtype=np.float64)
                slxy = np.asarray(sp["decayVtx_Lxy"][i], dtype=np.float64)
                sr3d = np.asarray(sp["decayVtx_r3d"][i], dtype=np.float64)
                dd = jeta[:, None] - se[None, :]
                dp = wrap_dphi(jphi[:, None] - sph[None, :])
                in_cone_sp = (dd * dd + dp * dp) < (seed_cone ** 2)        # (njet, nseed)
                pt_in = np.where(in_cone_sp, spt[None, :], -1.0)
                has = in_cone_sp.any(axis=1)
                jbest = np.argmax(pt_in, axis=1)                          # highest-pt in-cone SEED
                Lxy = np.where(has, slxy[jbest], np.nan)
                R3d = np.where(has, sr3d[jbest], np.nan)

        # process only the leading-n jets by Et (0 = all); jet_index = Et rank
        order = np.argsort(jpt)[::-1]
        if n_jets > 0:
            order = order[:n_jets]
        for rank, j in enumerate(order):
            cols["sample"].append(tag)
            cols["mcChannelNumber"].append(int(ev["mcChannelNumber"][i]))
            cols["eventNumber"].append(int(ev["eventNumber"][i]))
            cols["runNumber"].append(int(ev["runNumber"][i]))
            cols["event_weight"].append(float(ev["mcEventWeight"][i]))
            cols["jz_slice"].append(int(jz_slice[i]))
            cols["jz_weight"].append(float(jz_weight[i]))
            cols["pass_hstp"].append(int(pass_hstp[i]))
            cols["collection"].append(coll)
            cols["jet_index"].append(int(rank))
            cols["jet_et_gev"].append(float(jpt[j]))
            cols["jet_eta"].append(float(jeta[j]))
            cols["jet_phi"].append(float(jphi[j]))
            cols["jet_et_code"].append(int(jet_etc[j]))
            cols["jet_eta_code"].append(int(jet_etac[j]))
            cols["jet_phi_code"].append(int(jet_phic[j]))
            cols["n_towers_in_cone"].append(int(n_in_cone[j]))
            for l in range(NLAYERS):
                cols[f"q_l{l}"].append(int(per_layer_counts[j, l]))
                cols[f"et_l{l}_gev"].append(float(per_layer_gev[j, l]))
            cols["q_total"].append(int(per_layer_counts[j].sum()))
            cols["label"].append(int(ev["signal"][i]))
            cols["truth_decay_Lxy_mm"].append(float(Lxy[j]))
            cols["truth_decay_r3d_mm"].append(float(R3d[j]))

            if grid_hw > 0:
                side = 2 * grid_hw + 1
                grid = np.zeros((side, side, NLAYERS), dtype=np.float32)
                idx = np.nonzero(in_cone[j])[0] if teta.size else np.empty(0, dtype=int)
                if idx.size:
                    ci = np.round((teta[idx] - jeta[j]) / ETA_GRANULARITY).astype(int) + grid_hw
                    cj = np.round(wrap_dphi(tphi[idx] - jphi[j]) / PHI_GRANULARITY).astype(int) + grid_hw
                    ok = (ci >= 0) & (ci < side) & (cj >= 0) & (cj < side)
                    for a, b, q in zip(ci[ok], cj[ok], tow_qetl[idx][ok]):
                        grid[a, b] += q
                cols["tow_grid"].append(grid.reshape(-1).tolist())

        if (i + 1) % 2000 == 0:
            print(f"    ...{i + 1}/{n} events", flush=True)

    if n_hstp_dropped:
        print(f"    HSTP filter dropped {n_hstp_dropped}/{n} background events")
    return cols


def to_arrays(cols):
    """Cast the accumulated python lists to typed numpy arrays for pyarrow."""
    out = {}
    for k, v in cols.items():
        if k == "tow_grid":       # list column, built separately in main
            continue
        if k in ("sample", "collection"):
            out[k] = np.asarray(v, dtype=object)
        elif k == "eventNumber":
            out[k] = np.asarray(v, dtype=np.int64)
        elif (k.startswith("q_") or k.endswith("_code")
              or k in ("mcChannelNumber", "runNumber", "jet_index", "n_towers_in_cone", "label",
                       "jz_slice", "pass_hstp")):
            out[k] = np.asarray(v, dtype=np.int32)
        else:
            out[k] = np.asarray(v, dtype=np.float32)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="*", metavar="ROOT",
                    help="caloShowerShape_*.root files or globs "
                         f"(default: the two merged signals + the JZ dijet slices in {NTUPLE_DIR})")
    ap.add_argument("--no-background", action="store_true",
                    help="skip the default QCD dijet JZ slices (signals only)")
    ap.add_argument("--no-hstp-filter", action="store_true",
                    help="keep background events failing the HSTP filter (default: drop them)")
    ap.add_argument("--collection", choices=sorted(COLLECTIONS), default="wtacone",
                    help="jet collection to build features from (default: wtacone)")
    ap.add_argument("--et-min-tower", type=float, default=0.0, metavar="GEV",
                    help="optional per-layer tower zero-suppression before digitizing "
                         "(default: 0.0 = keep all, matching et_min_=0)")
    ap.add_argument("--rassoc", type=float, default=None, metavar="R",
                    help="override the jet-tower association cone (default: per collection)")
    ap.add_argument("--seed-cone", type=float, default=0.0, metavar="R",
                    help="dR cone to match a shower-parent SEED (showerParentTree) to a jet for "
                         "the regression label; the highest-pt in-cone SEED wins. "
                         "0 = use the collection Rassoc. Default 0.")
    ap.add_argument("--tower-grid-hw", type=int, default=4, metavar="W",
                    help="half-width (in 0.1 cells) of the per-jet tower grid; side=2W+1. "
                         "W=4 -> 9x9 bounding R=0.4 (0 = don't emit tow_grid). Default 4.")
    ap.add_argument("--n-jets", type=int, default=2, metavar="N",
                    help="process only the leading-N jets by Et per event (0 = all). Default 2.")
    ap.add_argument("--out", default=None, metavar="PARQUET",
                    help=f"output parquet (default: {OUT_DIR}/training_<collection>.parquet)")
    args = ap.parse_args()

    if args.rassoc is not None:
        COLLECTIONS[args.collection]["rassoc"] = args.rassoc

    patterns = list(args.inputs)
    if not patterns:
        patterns = [os.path.join(NTUPLE_DIR, f"caloShowerShape_{s}.root") for s in DEFAULT_SIGNALS]
        if not args.no_background:
            patterns.append(DEFAULT_DIJET_GLOB)
    inputs = expand_inputs(patterns)
    if not inputs:
        sys.exit(f"[error] no input ntuples found (looked in {NTUPLE_DIR})")

    out = args.out or os.path.join(OUT_DIR, f"training_{args.collection}.parquet")

    print(f"Collection : {args.collection}  (Rassoc={COLLECTIONS[args.collection]['rassoc']})")
    print(f"Et quant   : LSB={ET_GRANULARITY} GeV, {ET_BIT_LENGTH}-bit (max code {ET_MAX_CODE} "
          f"= {ET_MAX_CODE * ET_GRANULARITY:.1f} GeV), zero-suppress Et_l>{args.et_min_tower} GeV")
    print(f"eta/phi    : {ETA_GRANULARITY}/{PHI_GRANULARITY} grid, {ETA_BIT_LENGTH}/{PHI_BIT_LENGTH}-bit")
    if args.tower_grid_hw > 0:
        side = 2 * args.tower_grid_hw + 1
        print(f"Tower grid : {side}x{side}x{NLAYERS} per jet (hw={args.tower_grid_hw})")
    print(f"Jets       : leading {args.n_jets if args.n_jets > 0 else 'all'} per event")
    print(f"Background : HSTP filter {'off' if args.no_hstp_filter else 'on'}; "
          "JZ slices weighted via the jz_weight column (not by row count)")
    print(f"Inputs     : {len(inputs)} file(s)")

    merged = None
    for path in inputs:
        print(f"  reading {path}", flush=True)
        cols = process_file(path, args.collection, args.et_min_tower, args.seed_cone,
                            grid_hw=args.tower_grid_hw, n_jets=args.n_jets,
                            hstp_filter=not args.no_hstp_filter)
        if merged is None:
            merged = cols
        else:
            for k in merged:
                merged[k].extend(cols[k])

    arrays = to_arrays(merged)
    nrows = len(arrays["label"])
    if nrows == 0:
        sys.exit("[error] no jets found in the inputs — nothing to write")

    table_cols = dict(arrays)
    if args.tower_grid_hw > 0:
        table_cols["tow_grid"] = pa.array(merged["tow_grid"], type=pa.list_(pa.float32()))
    table = pa.table(table_cols)

    meta = {
        "producer": "make_training_parquet.py",
        "collection": args.collection,
        "rassoc": str(COLLECTIONS[args.collection]["rassoc"]),
        "n_jets": str(args.n_jets),
        "et_granularity_gev": str(ET_GRANULARITY),
        "et_bit_length": str(ET_BIT_LENGTH),
        "et_min_gev": str(ET_MIN_GEV),
        "et_max_gev": str(ET_MAX_GEV),
        "eta_grid": f"{ETA_GRANULARITY} over [{ETA_MIN},{ETA_MAX}], {ETA_BIT_LENGTH}-bit",
        "phi_grid": f"{PHI_GRANULARITY} over [{PHI_MIN},{PHI_MAX}], {PHI_BIT_LENGTH}-bit",
        "et_min_tower_gev": str(args.et_min_tower),
        "seed_cone": str(args.seed_cone if args.seed_cone > 0 else COLLECTIONS[args.collection]["rassoc"]),
        "features": ",".join(f"q_l{l}" for l in range(NLAYERS)),
        "label": "signal (1=displaced signal, 0=prompt/QCD)",
        "jz_slice": "QCD JZ slice of the event (-1 for signal)",
        "jz_weight": ("per-event JZ cross-section weight from the ntupler "
                      "(mcEventWeight * sigma * filterEff * L / sumOfWeights); use it for any "
                      "rate-like combination of the slices, whose row counts are not physical"),
        "hstp_filter": "off" if args.no_hstp_filter else "background events with passHSTP==0 dropped",
        "regression_target": ("truth_decay_r3d_mm / truth_decay_Lxy_mm = decay vertex of the "
                              "highest-pt showerParentTree SEED within the jet cone (NaN if none)"),
    }
    if args.tower_grid_hw > 0:
        side = 2 * args.tower_grid_hw + 1
        meta["tower_grid"] = (f"{side}x{side}x{NLAYERS} (dEta_row,dPhi_col,layer), "
                              f"row-major flat len {side * side * NLAYERS}")
        meta["tower_grid_side"] = str(side)
        meta["tower_grid_layers"] = str(NLAYERS)
    table = table.replace_schema_metadata(meta)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    pq.write_table(table, out)

    # summary
    labels = arrays["label"]
    print(f"\nWrote {nrows} jet rows -> {out}")
    for lab in np.unique(labels):
        print(f"  label={lab}: {int((labels == lab).sum())} jets")
    for s in np.unique(arrays["sample"]):
        print(f"  sample={s}: {int((arrays['sample'] == s).sum())} jets")
    jz = arrays["jz_slice"]
    for slice_id in sorted(int(v) for v in np.unique(jz) if v >= 0):
        sel = jz == slice_id
        print(f"  JZ{slice_id}: {int(sel.sum())} jets, sum(jz_weight)="
              f"{float(arrays['jz_weight'][sel].sum()):.6g}")


if __name__ == "__main__":
    main()
