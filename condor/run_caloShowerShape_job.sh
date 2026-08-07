#!/bin/bash
# Condor job wrapper for caloShowerShapeNTupler.C single-file mode.
# Arguments (passed by the submit script via Condor Arguments):
#   $1  signalBool   (0 or 1)
#   $2  signalString (e.g. "displaced_dark_photon"; "BACKGROUND" -> treated as "")
#   $3  outputDir    (absolute path, created if absent)
#   $4  daodFile     (absolute path to the specific DAOD_JETM42 pool.root file)
#   $5  gepFile      (absolute path to the matching GEPOutputReader ntuple file)
#   $6  fileSuffix   (string appended before .root in output name, e.g. _000001)

SIGNAL_BOOL="$1"
SIGNAL_STRING="$2"
OUTPUT_DIR="$3"
DAOD_FILE="$4"
GEP_FILE="$5"
FILE_SUFFIX="$6"

# The submit script passes "BACKGROUND" instead of an empty string so the
# positional arguments don't shift; translate it back to "" for the macro.
if [ "$SIGNAL_STRING" = "BACKGROUND" ]; then
    SIGNAL_STRING=""
fi

echo "=== caloShowerShapeNTupler Condor job ==="
echo "  SIGNAL_BOOL    = $SIGNAL_BOOL"
echo "  SIGNAL_STRING  = $SIGNAL_STRING"
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

# --- Run the ntupler (interpreted, xAOD types via autoloading) ---
ANALYSIS_DIR=/home/larsonma/GEPHadronicEventReconstruction/analysis
cd "$ANALYSIS_DIR"

root -b -q "caloShowerShapeStudy/caloShowerShapeNTupler.C(${SIGNAL_BOOL}, \"${SIGNAL_STRING}\", \"${OUTPUT_DIR}/\", \"${DAOD_FILE}\", \"${GEP_FILE}\", \"${FILE_SUFFIX}\")"

echo "=== Job complete ==="
