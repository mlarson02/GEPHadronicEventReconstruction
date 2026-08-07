#!/bin/bash
# Condor job wrapper for jetTaggerEmulation.cc single-config mode.
# Arguments (passed by submit_jet_tagger_emulation.py via Condor):
#   $1  rMerge              (e.g. 2)
#   $2  rSquared            (e.g. 1.21)
#   $3  nIOs                (e.g. 128)
#   $4  nSeeds              (e.g. 2)
#   $5  algoVersion         (e.g. 3)
#   $6  signal              (true or false)
#   $7  signalString        (e.g. ggF_hh_bbbb, or BACKGROUND)
#   $8  puSuppression       (true or false)
#   $9  inputObject         (e.g. gepCellsTowers)
#   $10 seedObject          (e.g. gepWTAConeCellsTowersJets)
#   $11 subjetEtThreshold   (e.g. 25)
#   $12 etWeightedMidpoint  (true or false)
#   $13 minEtSeedPosOpt     (true or false)
#   $14 minEtSeedPosCut     (e.g. 25.0)
#   $15 inputFilePath       (explicit path to input ntuple for per-file parallelism)
#   $16 fileIndex           (integer index appended to output name to avoid collisions)
#   $17 useEtaSKObjects     (true or false)

RMRG="$1"
R2="$2"
NIOS="$3"
NSEEDS="$4"
ALGOV="$5"
SIGNAL="$6"
SIGSTR="$7"
PUSUP="$8"
INOBJ="$9"
SEEDOBJ="${10}"
SJETT="${11}"
ETWM="${12}"
MINETO="${13}"
MINETC="${14}"
INFILE="${15}"
FIDX="${16}"
ETASK="${17}"

echo "=== JetTagger Emulation Condor job ==="
echo "  rMerge=$RMRG  rSquared=$R2  nIOs=$NIOS  nSeeds=$NSEEDS  algo=$ALGOV"
echo "  signal=$SIGNAL  signalString=$SIGSTR  puSup=$PUSUP  etaSK=$ETASK"
echo "  inputObj=$INOBJ  seedObj=$SEEDOBJ  subjetEt=$SJETT"
echo "  etWtMid=$ETWM  minEtOpt=$MINETO  minEtCut=$MINETC"
echo "  inputFile=$INFILE  fileIndex=$FIDX"
echo ""

# --- ATLAS/ROOT environment setup ---
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
export ALRB_localConfigDir=$HOME/localConfig
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
asetup AnalysisBase,25.2.29

set -e

# --- Reconstruct LUT/constants file paths (mirrors local.sh logic) ---
EMUL_DIR=/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation
LUT_DIR=$EMUL_DIR/LUT_Constants_Generation/LUTs
CONST_DIR=$EMUL_DIR/LUT_Constants_Generation/constants

rMergeFmt=$(python3 -c "print(f'{float(\"$RMRG\"):.4g}')")
r2Fmt=$(python3    -c "print(f'{float(\"$R2\"):.3g}')")

CONSTANTS_FILE="$CONST_DIR/constants_rMerge_${rMergeFmt}_R2_${r2Fmt}_IOs_${NIOS}_Seeds_${NSEEDS}_v${ALGOV}.h"
LUT_FILE="$LUT_DIR/LUT_deltaR_8b_rMerge_${rMergeFmt}_R2_${r2Fmt}.h"

echo "  constants: $CONSTANTS_FILE"
echo "  LUT:       $LUT_FILE"

if [ ! -f "$CONSTANTS_FILE" ]; then
    echo "[error] constants file not found: $CONSTANTS_FILE" >&2; exit 1
fi
if [ ! -f "$LUT_FILE" ]; then
    echo "[error] LUT file not found: $LUT_FILE" >&2; exit 1
fi

# --- Per-job temp directory on local worker scratch (avoids NFS concurrency) ---
JOBDIR=$(mktemp -d /tmp/jet_tagger_emulation_XXXXXXXX)
trap "rm -rf '$JOBDIR'" EXIT

cp "$EMUL_DIR/jetTaggerEmulation.cc"      "$JOBDIR/"
cp "$EMUL_DIR/emulationHelperFunctions.h" "$JOBDIR/"
mkdir "$JOBDIR/UsedConstantsLUTs"
cp "$CONSTANTS_FILE" "$JOBDIR/UsedConstantsLUTs/constants.h"
cp "$LUT_FILE"       "$JOBDIR/UsedConstantsLUTs/deltaRLUT_8b.h"

cd "$JOBDIR"

# --- Run emulation ---
echo "Running emulation in: $JOBDIR"
root -l -b -q "jetTaggerEmulation.cc+(${RMRG}, ${NIOS}, ${NSEEDS}, ${R2}, ${SIGNAL}, true, ${PUSUP}, \"${SIGSTR}\", \"${INOBJ}\", \"${SEEDOBJ}\", ${SJETT}, ${ETWM}, ${MINETO}, ${MINETC}, false, \"${INFILE}\", ${FIDX}, ${ETASK})"

echo "=== Job complete ==="
