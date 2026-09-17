#!/bin/bash
# Condor job wrapper for the largeRJetAnalysisAndRates compute stage
# (analysis/largeRJetCompute.C — see analysis/largeRJetSplitPlan.md).
# Arguments (passed by submit_largeRJetCompute.py via Condor Arguments):
#   $1  pairIndex  (index into the generated pair list)
#   $2  jobIndex   (which 1/nJobs slice of the background entries)
#   $3  nJobs      (total slices per pair)
#   $4  pairsFile  (absolute path to the generated largeRJetFilePairs.txt)

PAIR_INDEX="$1"
JOB_INDEX="$2"
N_JOBS="$3"
PAIRS_FILE="$4"

echo "=== largeRJetCompute Condor job ==="
echo "  PAIR_INDEX = $PAIR_INDEX"
echo "  JOB_INDEX  = $JOB_INDEX"
echo "  N_JOBS     = $N_JOBS"
echo "  PAIRS_FILE = $PAIRS_FILE"
echo ""

# --- ATLAS environment setup ---
# atlasLocalSetup.sh intentionally returns non-zero at points during setup, so
# set -e must not be active while it runs.
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
export ALRB_localConfigDir=$HOME/localConfig
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh --quiet
# Same ROOT the other analysis-stage jobs use (run_jet_tagger_emulation_job.sh etc.).
asetup AnalysisBase,25.2.29

# The macro uses relative includes and is normally launched from the analysis
# directory; state output goes to /data (absolute), so nothing of value lands in
# the job scratch.
cd /home/larsonma/GEPHadronicEventReconstruction/analysis || exit 1

root -b -l -q "largeRJetCompute.C(${PAIR_INDEX}, ${JOB_INDEX}, ${N_JOBS}, \"${PAIRS_FILE}\")"
STATUS=$?
echo "root exit status: ${STATUS}"
exit ${STATUS}
