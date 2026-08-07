#!/bin/bash
set -e
SCRIPT="$(dirname "$0")/submit_ntupler.py"
ALGO=4

DAOD_BASE=/data/larsonma/GEPHadronicEventReconstruction/JETM42_DAODs
GEP_BASE=/data/larsonma/GEPHadronicEventReconstruction/GEPOutputReaderNTuples
OUT_BASE=/data/larsonma/GEPHadronicEventReconstruction/ntuples

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

# --- Signal: ggF HH->bbbb ---
submit ggF_HHbbbb \
    --signal ggF_hh_bbbb --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/ggF_HHbbbb EXT1)" \
    --gep-dir  "$(find_container $GEP_BASE/ggF_HHbbbb  EXT0)" \
    --output-dir $OUT_BASE/ggF_HHbbbb_v$ALGO

# --- Signal: ZvvHbb ---
submit ZvvHbb \
    --signal ZvvHbb --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/ZvvHbb EXT1)" \
    --gep-dir  "$(find_container $GEP_BASE/ZvvHbb  EXT0)" \
    --output-dir $OUT_BASE/ZvvHbb_v$ALGO

# --- Signal: VBF HH->bbbb ---
submit VBF_HHbbbb \
    --signal VBF_HHbbbb --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/VBF_HHbbbb EXT1)" \
    --gep-dir  "$(find_container $GEP_BASE/VBF_HHbbbb  EXT0)" \
    --output-dir $OUT_BASE/VBF_HHbbbb_v$ALGO

# --- Signal: Z' -> ttbar (all-hadronic, flat pT) ---
submit Zprime_ttbar_allhad_flatpT \
    --signal Zprime_ttbar_allhad_flatpT --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/Zprime_ttbar_allhad_flatpT EXT1)" \
    --gep-dir  "$(find_container $GEP_BASE/Zprime_ttbar_allhad_flatpT  EXT0)" \
    --output-dir $OUT_BASE/Zprime_ttbar_allhad_flatpT_v$ALGO

# --- Signal: ttbar (all-hadronic) ---
submit ttbar_allhad \
    --signal ttbar_allhad --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/ttbar_allhad EXT1)" \
    --gep-dir  "$(find_container $GEP_BASE/ttbar_allhad  EXT0)" \
    --output-dir $OUT_BASE/ttbar_allhad_v$ALGO

# --- Signal: ttbar (semi-leptonic) ---
submit ttbar_semilep \
    --signal ttbar_semilep --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/ttbar_semilep EXT1)" \
    --gep-dir  "$(find_container $GEP_BASE/ttbar_semilep  EXT0)" \
    --output-dir $OUT_BASE/ttbar_semilep_v$ALGO

# --- Signal: ttbar (di-leptonic) ---
submit ttbar_dilep \
    --signal ttbar_dilep --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/ttbar_dilep EXT1)" \
    --gep-dir  "$(find_container $GEP_BASE/ttbar_dilep  EXT0)" \
    --output-dir $OUT_BASE/ttbar_dilep_v$ALGO

# --- Background: QCD Dijet JZ0-JZ9 ---
for JZ in $(seq 0 9); do
    submit QCD_Dijet_JZ${JZ} \
        --background --jz $JZ --algo $ALGO \
        --daod-dir "$(find_container $DAOD_BASE/QCD_Dijet/JZ${JZ} EXT1)" \
        --gep-dir  "$(find_container $GEP_BASE/QCD_Dijet/JZ${JZ}  EXT0)" \
        --output-dir $OUT_BASE/QCD_Dijet_JZ${JZ}_v$ALGO
done
