#!/bin/bash
# Derive the GEP jet energy scale (JES) LUTs for the WTACone CellsTowers jets
# (no PU suppression / SK / EtaSK) from the QCD dijet JZ0-JZ9 HERNTupler ntuples
# of one pileup scenario. Runs locally, no condor jobs.
#
# Writes GepJetCalibHist.root (response histograms) and GEP_JES_LUTs.root (the
# LUTs) to gepJetCalib/PU<pu>/, where HERNTupler.C looks for them
# (gep_jes_lut_file_pu200_ / gep_jes_lut_file_pu140_).
set -e
ALGO=4
PU=200

usage() {
    cat <<EOF
Usage: $(basename "$0") [--pu 140|200] [--algo N]

Run once every JZ0-JZ9 ntupling job of that pileup has finished, e.g. after
  submit_all_ntupler.sh --jz-only [--pu 140]

  --pu 140    Use the PU140 ntuples (ntuples_PU140/, r16129).
  --pu 200    Default: PU200 ntuples (ntuples/, r16130).
  --algo N    Algorithm version of the ntuple directories (default: $ALGO).
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --pu)     PU="$2";   shift 2 ;;
        --algo)   ALGO="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *)        echo "[error] unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# The reconstruction tag in the file names encodes the pileup (r16130 = PU200, r16129 = PU140).
case "$PU" in
    200) PU_SUFFIX="";       RTAG=r16130 ;;
    140) PU_SUFFIX="_PU140"; RTAG=r16129 ;;
    *)   echo "[error] --pu must be 140 or 200 (got: $PU)" >&2; exit 1 ;;
esac

CALIBRATOR=/home/larsonma/gepjetcalibrator
NTUPLE_BASE=/data/larsonma/GEPHadronicEventReconstruction/ntuples${PU_SUFFIX}
OUT_DIR=/data/larsonma/GEPHadronicEventReconstruction/gepJetCalib/PU${PU}

# Per-job ntuples only (..._GEP_NNNNNN.root): the hadd'ed ..._GEP.root would count every event twice.
shopt -s nullglob
FILES=()
for JZ in $(seq 0 9); do
    jz_files=( "$NTUPLE_BASE"/QCD_Dijet_JZ${JZ}_v${ALGO}/mc21_14TeV_jj_JZ${JZ}_e8557_s4422_${RTAG}_DAOD_NTUPLE_GEP_[0-9]*.root )
    echo "JZ${JZ}: ${#jz_files[@]} ntuple(s)"
    if [ "${#jz_files[@]}" -eq 0 ]; then
        echo "[error] no JZ${JZ} ntuples in $NTUPLE_BASE/QCD_Dijet_JZ${JZ}_v${ALGO}/" >&2; exit 1
    fi
    FILES+=( "${jz_files[@]}" )
done

# Uses the ROOT of the calling shell, e.g. after lsetup "views LCG_107 x86_64-el9-gcc13-opt"
set -o pipefail

# --- Build (incremental after the first time) ---
for TOOL in HERCalibHistMaker CalibLUTMaker; do
    echo "==> Building $TOOL"
    mkdir -p "$CALIBRATOR/$TOOL/build"
    cd "$CALIBRATOR/$TOOL/build"
    [ -f Makefile ] || cmake ../source
    make -j4
done

# --- Response histograms, then the LUT fits ---
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
echo "==> Filling response histograms from ${#FILES[@]} PU${PU} ntuples"
"$CALIBRATOR/HERCalibHistMaker/build/runHERCalibHistMaker" GepJetCalibHist.root "${FILES[@]}" 2>&1 | tee runHERCalibHistMaker.log

# One Gaussian fit per (eta, pT) bin and jet collection: far too verbose for the terminal.
echo "==> Fitting LUTs (log: $OUT_DIR/runCalibLUTMakerHER.log)"
"$CALIBRATOR/CalibLUTMaker/build/runCalibLUTMakerHER" GepJetCalibHist.root > runCalibLUTMakerHER.log 2>&1

echo "==> Done: $OUT_DIR/GEP_JES_LUTs.root"
