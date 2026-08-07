#!/bin/bash
set -e
SCRIPT="$(dirname "$0")/submit_caloShowerShape.py"

DAOD_BASE=/data/larsonma/CaloShowerShapeTriggers/JETM42_DAODs
GEP_BASE=/data/larsonma/CaloShowerShapeTriggers/GEPOutputReaderNTuples
OUT_BASE=/data/larsonma/CaloShowerShapeTriggers/ntuples

# Find the unique user.mlarson*_EXT{0,1} subdirectory under a given base path.
find_container() {
    local base="$1" ext="$2"
    local matches=("$base"/user.mlarson*_${ext})
    if [ ! -d "${matches[0]}" ]; then
        echo "[error] no ${ext} container found in: $base" >&2; exit 1
    fi
    if [ "${#matches[@]}" -gt 1 ]; then
        echo "[warn] multiple ${ext} containers in $base, using: ${matches[0]}" >&2
    fi
    echo "${matches[0]}"
}

submit() {
    local label="$1"; shift
    echo "==> Submitting: $label"
    python3 "$SCRIPT" "$@"
}

# --- Signal: displaced dark photon (HAHM Zd->4e) ---
submit displaced_dark_photon \
    --signal displaced_dark_photon \
    --daod-dir "$(find_container $DAOD_BASE/DisplacedDarkPhoton EXT1)" \
    --gep-dir  "$(find_container $GEP_BASE/DisplacedDarkPhoton  EXT0)" \
    --output-dir $OUT_BASE

# --- Signal: emerging jets (Z' -> dark quarks) ---
submit emerging_jets \
    --signal emerging_jets \
    --daod-dir "$(find_container $DAOD_BASE/EmergingJets EXT1)" \
    --gep-dir  "$(find_container $GEP_BASE/EmergingJets  EXT0)" \
    --output-dir $OUT_BASE

# --- Background: QCD dijet, for the "what does background look like" comparison.
# One or more JZ slices; --label keeps their outputs/logs from colliding.
# Edit DIJET_JZ to pick which slices to run (a couple is usually enough).
# NOTE: QCD dijet is not part of the CaloShowerShapeTriggers download yet — update
# the QCD_Dijet/JZ* sub-paths below to match the new layout once it is staged.
DIJET_JZ="2 4"
for JZ in $DIJET_JZ; do
    submit dijet_JZ${JZ} \
        --background --label dijet_JZ${JZ} \
        --daod-dir "$(find_container $DAOD_BASE/QCD_Dijet/JZ${JZ} EXT1)" \
        --gep-dir  "$(find_container $GEP_BASE/QCD_Dijet/JZ${JZ}  EXT0)" \
        --output-dir $OUT_BASE
done
