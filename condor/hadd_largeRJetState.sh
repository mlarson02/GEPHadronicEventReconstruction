#!/bin/bash
# Stage 2 of the largeRJetAnalysisAndRates split: merge each pair's compute chunks
#   <base>__job<i>of<N>.root  ->  <base>__state.root
# in the intermediate-state directory. hadd sums histogram bins (raw counts and the
# 1-bin scalar histograms) and appends the per-event cache trees — both exactly the
# right merge. See analysis/largeRJetSplitPlan.md.
#
# Usage: ./hadd_largeRJetState.sh [options] [pattern]
#
#   pattern         only merge bases matching this glob (default: all)
#   --threads N     pass -j N to hadd for multi-threaded merging (default: off).
#                   Same meaning as --threads in hadd_emulator_outputs.sh.
#   --jobs N        run N hadd processes concurrently (default: 4). Each state file is
#                   independent, so this is the bigger win when several configurations
#                   finished at once; --threads helps a single large merge.
#   --force         Overwrite existing __state.root files (this is the DEFAULT)
#   --skip-existing Merge only bases whose __state.root does not exist yet, leaving
#                   already-merged configurations untouched. The fast path when only a
#                   few configurations were just computed.
#                   Skips on existence alone: if new chunks landed after a merge was
#                   made, that merge is stale and STAYS stale — the skip message flags
#                   it, and --force is how you rebuild it. (Same rule as
#                   hadd_emulator_outputs.sh, so the two behave alike.)
#
# Examples:
#   ./hadd_largeRJetState.sh                        # all, 4 concurrent merges
#   ./hadd_largeRJetState.sh --threads 4 --jobs 8   # 8 merges at a time, 4 threads each
#   ./hadd_largeRJetState.sh --skip-existing        # only newly computed configurations
#   ./hadd_largeRJetState.sh '*r16129*'             # PU140 only

set -e

# Must match gLRJState.stateDir (analysis/largeRJetStateRegistry.h) and STATE_DIRS in
# submit_largeRJetCompute.py — one state directory per emulation production.
STATE_DIR="${LRJ_STATE_DIR:-/data/larsonma/LargeRadiusJets/largeRJetIntermediateState_towerEtCutFixed}"
PATTERN="*"
HADD_THREADS=0
PARALLEL_JOBS=4
# Overwrite merged outputs by default; --skip-existing is what turns this off.
FORCE=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)       HADD_THREADS="$2"; shift 2 ;;
        --jobs)          PARALLEL_JOBS="$2"; shift 2 ;;
        --force)         FORCE=1; shift ;;
        --skip-existing) FORCE=0; shift ;;
        -h|--help)       sed -n '2,28p' "$0"; exit 0 ;;
        -*)              echo "Unknown option: $1" >&2; exit 2 ;;
        *)               PATTERN="$1"; shift ;;
    esac
done

cd "${STATE_DIR}"

# Distinct pair bases, from the chunk-file names.
mapfile -t BASES < <(ls ${PATTERN}__job*of*.root 2>/dev/null | sed 's/__job[0-9]*of[0-9]*\.root$//' | sort -u)
if [ ${#BASES[@]} -eq 0 ]; then
    echo "No chunk files matching '${PATTERN}__job*of*.root' in ${STATE_DIR}" >&2
    exit 1
fi

# hadd's own -j flag; empty when --threads was not given, so the command line is
# unchanged from the single-threaded version in that case.
HADD_J=()
[ "${HADD_THREADS}" -gt 0 ] && HADD_J=( -j "${HADD_THREADS}" )

FAILED=0
SKIPPED=0
for BASE in "${BASES[@]}"; do
    # --skip-existing: leave an already-merged configuration alone. Checked before the
    # chunk glob so a skipped base costs one stat rather than a glob plus a count.
    if [ -f "${BASE}__state.root" ] && [ "${FORCE}" -eq 0 ]; then
        # A state file older than one of its chunks is a PARTIAL merge — chunks landed
        # after it was built. Still skipped (that is what was asked for), but said out
        # loud, because silently keeping a stale merge is the one way this option bites.
        NEWER=$(find . -maxdepth 1 -name "$(basename "${BASE}")__job*of*.root" \
                     -newer "${BASE}__state.root" -print -quit 2>/dev/null || true)
        if [ -n "${NEWER}" ]; then
            echo "[skip] exists but CHUNKS ARE NEWER — stale, use --force: ${BASE}__state.root"
        else
            echo "[skip] exists (--force to overwrite): ${BASE}__state.root"
        fi
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    CHUNKS=( ${BASE}__job*of*.root )
    # Every chunk of a pair declares its N in the name; refuse to merge an incomplete set.
    N_EXPECTED=$(echo "${CHUNKS[0]}" | sed 's/.*__job[0-9]*of\([0-9]*\)\.root/\1/')
    if [ ${#CHUNKS[@]} -ne "${N_EXPECTED}" ]; then
        echo "SKIP ${BASE}: found ${#CHUNKS[@]} of ${N_EXPECTED} chunks (jobs still running or failed?)" >&2
        continue
    fi
    echo "hadd: ${BASE}__state.root  <-  ${#CHUNKS[@]} chunks"
    # Each state file is independent, so merges run concurrently up to --jobs. Output is
    # per-merge and interleaves; the summary below is what says whether all of them worked.
    hadd -f "${HADD_J[@]}" "${BASE}__state.root" "${CHUNKS[@]}" \
        > "${BASE}__hadd.log" 2>&1 || { echo "FAILED: ${BASE} (see ${BASE}__hadd.log)" >&2; } &
    # Throttle: wait whenever PARALLEL_JOBS merges are already in flight.
    while [ "$(jobs -rp | wc -l)" -ge "${PARALLEL_JOBS}" ]; do wait -n; done
done
# set -e would abort on a non-zero background job, so collect statuses explicitly.
set +e
wait
set -e

# With merges running concurrently their output interleaves, so verify by result rather
# than by reading the log: every base that had a complete chunk set must now have a state
# file. A missing one means that merge failed — its __hadd.log has the reason.
# A base with no __hadd.log was never attempted this run (skipped, or an incomplete chunk
# set), so only bases we actually tried are judged here.
echo ""
MADE=0
for BASE in "${BASES[@]}"; do
    [ -f "${BASE}__hadd.log" ] || continue
    if [ -f "${BASE}__state.root" ]; then
        MADE=$((MADE + 1))
        rm -f "${BASE}__hadd.log"
    else
        FAILED=$((FAILED + 1))
        echo "MISSING: ${BASE}__state.root — check ${BASE}__hadd.log" >&2
    fi
done
echo "Merged ${MADE} state file(s); ${SKIPPED} skipped; ${FAILED} failure(s)."
[ "${FAILED}" -gt 0 ] && exit 1

echo "Done. Next: cd /home/larsonma/GEPHadronicEventReconstruction/analysis && root -b -l -q 'largeRJetPlot.C'"
