#!/bin/bash
# Condor job wrapper for metEmulation.cc.
# Arguments (passed by submit_met_emulation.py via Condor):
#   $1   signal              (true or false)
#   $2   puSuppression       (true or false)
#   $3   signalString        (e.g. ZvvHbb, or BACKGROUND)
#   $4   jetEtThreshold      (e.g. 20.0)
#   $5   doJetTowerOR        (true or false)
#   $6   towerEtThreshold    (e.g. 1.0)
#   $7   useEtaSKObjects     (true or false)
#   $8   inputFilePath(s)    (single path, OR `;`-separated list to hadd on worker)
#   $9   fileIndex           (integer index appended to output name; for batched
#                             jobs this is the FIRST fidx in the chunk)
#   ${10} towerScaleFactor   (scalar weight on tower MET in totalMET sum, e.g. 1.0)
#   ${11} jetScaleFactor     (scalar weight on jet MET in totalMET sum, e.g. 0.4)
#
# Ordering: when $8 is a `;`-separated list, the wrapper hadds the inputs in
# the EXACT order given. The submit script provides them in sorted-fidx order,
# so the resulting merged input — and therefore the chunked output — preserves
# global event order across jobs.

SIGNAL="$1"
PUSUP="$2"
SIGSTR="$3"
JETET="$4"
JTOR="$5"
TOWERET="$6"
ETASK="$7"
INFILE="$8"
FIDX="$9"
TWRSF="${10}"
JETSF="${11}"

echo "=== MET Emulation Condor job ==="
echo "  signal=$SIGNAL  signalString=$SIGSTR  puSup=$PUSUP  etaSK=$ETASK"
echo "  jetEtThreshold=$JETET  doJetTowerOR=$JTOR  towerEtThreshold=$TOWERET"
echo "  inputFile=$INFILE  fileIndex=$FIDX"
echo "  towerScaleFactor=$TWRSF  jetScaleFactor=$JETSF"
echo ""

# --- ATLAS/ROOT environment setup ---
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
export ALRB_localConfigDir=$HOME/localConfig
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
asetup AnalysisBase,25.2.29

set -e

# --- Per-job temp directory on local worker scratch (avoids NFS concurrency) ---
EMUL_DIR=/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation
JOBDIR=$(mktemp -d /tmp/met_emulation_XXXXXXXX)
trap "rm -rf '$JOBDIR'" EXIT

cp "$EMUL_DIR/metEmulation.cc"                "$JOBDIR/"
cp "$EMUL_DIR/emulationHelperFunctions_MET.h" "$JOBDIR/"
cp -r "$EMUL_DIR/metConstants"                "$JOBDIR/"

cd "$JOBDIR"

# --- Resolve the actual input to feed metEmulation ---
# Split INFILE on `;` (used since Condor's submit-file ItemData uses `,`).
# Order MUST be preserved end-to-end: we hadd in the order given, so the
# resulting metTree events appear in the same global order as a top-level
# hadd of the original inputs in sorted-fidx order.
IFS=';' read -r -a INPUT_LIST <<< "$INFILE"
NUM_INPUTS=${#INPUT_LIST[@]}
echo "Inputs in this chunk ($NUM_INPUTS):"
for i in "${!INPUT_LIST[@]}"; do
    echo "  [$i] ${INPUT_LIST[$i]}"
done

if [[ $NUM_INPUTS -le 1 ]]; then
    # Single-input case — identical to pre-batching behavior.
    EMU_INPUT="${INPUT_LIST[0]}"
else
    # Multi-input case — hadd locally on worker scratch, preserving order.
    EMU_INPUT="$JOBDIR/chunk_merged_input.root"
    echo "Merging $NUM_INPUTS inputs into $EMU_INPUT (preserving listed order)"
    hadd -f "$EMU_INPUT" "${INPUT_LIST[@]}"
    echo "Merged input ready ($(du -h "$EMU_INPUT" | cut -f1))"
fi

# --- Run emulation ---
echo "Running emulation in: $JOBDIR"
root -l -b -q "metEmulation.cc+(${SIGNAL}, ${PUSUP}, \"${SIGSTR}\", ${JETET}, ${JTOR}, ${TOWERET}, ${ETASK}, \"${EMU_INPUT}\", ${FIDX}, ${TWRSF}, ${JETSF})"

echo "=== Job complete ==="
