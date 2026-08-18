#!/bin/bash
set -e
SCRIPT="$(dirname "$0")/submit_ntupler.py"
ALGO=4
PU=200
EXTRA_ARGS=()

usage() {
    cat <<EOF
Usage: $(basename "$0") [--pu 140|200] [--algo N] [extra submit_ntupler.py args...]

  --pu 140    Process the PU140 samples: reads from the *_PU140 directory
              mirrors and passes --pu 140 to submit_ntupler.py (selects the
              PU140 reweighting constants in HERNTupler).
  --pu 200    Default: nominal PU200 directories.
  --algo N    Algorithm version (default: $ALGO).

Any other argument is forwarded verbatim to submit_ntupler.py, e.g.
  $(basename "$0") --pu 140 --dry-run
  $(basename "$0") --pu 140 --max-jobs 2
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --pu)     PU="$2";   shift 2 ;;
        --algo)   ALGO="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *)        EXTRA_ARGS+=("$1"); shift ;;
    esac
done

# PU140 samples live in parallel *_PU140 directory trees.
case "$PU" in
    200) PU_SUFFIX="" ;;
    140) PU_SUFFIX="_PU140" ;;
    *)   echo "[error] --pu must be 140 or 200 (got: $PU)" >&2; exit 1 ;;
esac

DAOD_BASE=/data/larsonma/GEPHadronicEventReconstruction/JETM42_DAODs${PU_SUFFIX}
GEP_BASE=/data/larsonma/GEPHadronicEventReconstruction/GEPOutputReaderNTuples${PU_SUFFIX}
OUT_BASE=/data/larsonma/GEPHadronicEventReconstruction/ntuples${PU_SUFFIX}

# Find the unique user.mlarson*_EXT{0,1} subdirectory under a given base path.
# An optional third argument is a substring the container name must contain —
# used to reject stray containers copied into the wrong sample directory
# (e.g. a JZ7 container sitting in .../QCD_Dijet/JZ6/).
find_container() {
    local base="$1" ext="$2" tag="$3"
    local matches=() d
    for d in "$base"/user.mlarson*_${ext}; do
        [ -d "$d" ] || continue
        if [ -n "$tag" ] && [[ "$(basename "$d")" != *"$tag"* ]]; then
            continue
        fi
        matches+=("$d")
    done
    if [ "${#matches[@]}" -eq 0 ]; then
        echo "[error] no ${ext} container found in: $base" >&2; exit 1
    fi
    if [ "${#matches[@]}" -gt 1 ]; then
        echo "[warn] multiple ${ext} containers in $base, using: ${matches[0]}" >&2
    fi
    echo "${matches[0]}"
}

submit() {
    local label="$1"; shift
    echo "==> Submitting: $label  (PU${PU}, algo v${ALGO})"
    python3 "$SCRIPT" --pu "$PU" "$@" "${EXTRA_ARGS[@]}"
}

# --- Signal: ggF HH->bbbb ---
submit ggF_HHbbbb \
    --signal ggF_hh_bbbb --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/ggF_HHbbbb EXT1 .ggF_HHbbbb.)" \
    --gep-dir  "$(find_container $GEP_BASE/ggF_HHbbbb  EXT0 .ggF_HHbbbb.)" \
    --output-dir $OUT_BASE/ggF_HHbbbb_v$ALGO

# --- Signal: ZvvHbb ---
submit ZvvHbb \
    --signal ZvvHbb --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/ZvvHbb EXT1 .ZvvHbb.)" \
    --gep-dir  "$(find_container $GEP_BASE/ZvvHbb  EXT0 .ZvvHbb.)" \
    --output-dir $OUT_BASE/ZvvHbb_v$ALGO

# --- Signal: VBF HH->bbbb ---
submit VBF_HHbbbb \
    --signal VBF_HHbbbb --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/VBF_HHbbbb EXT1 .VBF_HHbbbb.)" \
    --gep-dir  "$(find_container $GEP_BASE/VBF_HHbbbb  EXT0 .VBF_HHbbbb.)" \
    --output-dir $OUT_BASE/VBF_HHbbbb_v$ALGO

# --- Signal: Z' -> ttbar (all-hadronic, flat pT) ---
submit Zprime_ttbar_allhad_flatpT \
    --signal Zprime_ttbar_allhad_flatpT --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/Zprime_ttbar_allhad_flatpT EXT1 .Zprime_ttbar_allhad_flatpT.)" \
    --gep-dir  "$(find_container $GEP_BASE/Zprime_ttbar_allhad_flatpT  EXT0 .Zprime_ttbar_allhad_flatpT.)" \
    --output-dir $OUT_BASE/Zprime_ttbar_allhad_flatpT_v$ALGO

# --- Signal: ttbar (all-hadronic) ---
submit ttbar_allhad \
    --signal ttbar_allhad --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/ttbar_allhad EXT1 .ttbar_allhad.)" \
    --gep-dir  "$(find_container $GEP_BASE/ttbar_allhad  EXT0 .ttbar_allhad.)" \
    --output-dir $OUT_BASE/ttbar_allhad_v$ALGO

# --- Signal: ttbar (semi-leptonic) ---
submit ttbar_semilep \
    --signal ttbar_semilep --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/ttbar_semilep EXT1 .ttbar_semilep.)" \
    --gep-dir  "$(find_container $GEP_BASE/ttbar_semilep  EXT0 .ttbar_semilep.)" \
    --output-dir $OUT_BASE/ttbar_semilep_v$ALGO

# --- Signal: ttbar (di-leptonic) ---
submit ttbar_dilep \
    --signal ttbar_dilep --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/ttbar_dilep EXT1 .ttbar_dilep.)" \
    --gep-dir  "$(find_container $GEP_BASE/ttbar_dilep  EXT0 .ttbar_dilep.)" \
    --output-dir $OUT_BASE/ttbar_dilep_v$ALGO

# --- Signal: Z -> mu mu ---
submit Zmumu \
    --signal Zmumu --algo $ALGO \
    --daod-dir "$(find_container $DAOD_BASE/Zmumu EXT1 .Zmumu.)" \
    --gep-dir  "$(find_container $GEP_BASE/Zmumu  EXT0 .Zmumu.)" \
    --output-dir $OUT_BASE/Zmumu_v$ALGO

# --- Background: QCD Dijet JZ0-JZ9 ---
for JZ in $(seq 0 9); do
    submit QCD_Dijet_JZ${JZ} \
        --background --jz $JZ --algo $ALGO \
        --daod-dir "$(find_container $DAOD_BASE/QCD_Dijet/JZ${JZ} EXT1 .QCD_Dijet_JZ${JZ}.)" \
        --gep-dir  "$(find_container $GEP_BASE/QCD_Dijet/JZ${JZ}  EXT0 .QCD_Dijet_JZ${JZ}.)" \
        --output-dir $OUT_BASE/QCD_Dijet_JZ${JZ}_v$ALGO
done
