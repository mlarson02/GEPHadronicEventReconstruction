#!/bin/bash
# Condor job wrapper for HERNTupler.C single-file mode.
# Arguments (passed by the submit script via Condor Arguments):
#   $1  signalBool   (0 or 1)
#   $2  signalString (e.g. "ggF_hh_bbbb_resim", or "" for background)
#   $3  algoVersion  (2 or 3)
#   $4  jzSlice      (0-9; ignored for signal)
#   $5  specialJZ0   (0 or 1; only relevant for JZ0 background)
#   $6  outputDir    (absolute path, created if absent)
#   $7  daodFile     (absolute path to the specific DAOD_JETM42 pool.root file)
#   $8  gepFile      (absolute path to the matching TrigGepPerf ntuple file)
#   $9  fileSuffix   (string appended before .root in output name, e.g. _000001)

SIGNAL_BOOL="$1"
SIGNAL_STRING="$2"
ALGO_VERSION="$3"
JZ_SLICE="$4"
SPECIAL_JZ0="$5"
OUTPUT_DIR="$6"
DAOD_FILE="$7"
GEP_FILE="$8"
FILE_SUFFIX="$9"

echo "=== HERNTupler Condor job ==="
echo "  SIGNAL_BOOL    = $SIGNAL_BOOL"
echo "  SIGNAL_STRING  = $SIGNAL_STRING"
echo "  ALGO_VERSION   = $ALGO_VERSION"
echo "  JZ_SLICE       = $JZ_SLICE"
echo "  SPECIAL_JZ0    = $SPECIAL_JZ0"
echo "  OUTPUT_DIR     = $OUTPUT_DIR"
echo "  DAOD_FILE      = $DAOD_FILE"
echo "  GEP_FILE       = $GEP_FILE"
echo "  FILE_SUFFIX    = $FILE_SUFFIX"
echo ""

# --- ATLAS environment setup ---
# atlasLocalSetup.sh and asetup intentionally return non-zero at points during
# setup, so set -e must not be active while they run.
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
export ALRB_localConfigDir=$HOME/localConfig
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
asetup AnalysisBase,25.2.29

set -e

# --- Grid proxy (pre-created and transferred via Condor) ---
if [ -f "$HOME/x509proxy" ]; then
    export X509_USER_PROXY=$HOME/x509proxy
fi

# --- Create output directory if needed ---
mkdir -p "$OUTPUT_DIR"

# --- Run the ntupler ---
ANALYSIS_DIR=/home/larsonma/GEPHadronicEventReconstruction/analysis
cd "$ANALYSIS_DIR"

root -b -q "HERNTupler.C(${SIGNAL_BOOL}, \"${SIGNAL_STRING}\", ${ALGO_VERSION}, ${JZ_SLICE}, ${SPECIAL_JZ0}, \"${OUTPUT_DIR}/\", \"${DAOD_FILE}\", \"${GEP_FILE}\", \"${FILE_SUFFIX}\")"

echo "=== Job complete ==="
