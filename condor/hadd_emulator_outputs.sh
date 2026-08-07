#!/bin/bash
# Merge per-file outputs produced by parallel Condor jobs.
#
# Three modes:
#
#   emulator (default)
#     Input:  {base}_file{N}.root  (emulator per-job outputs)
#     Output: {base}.root
#     DIR defaults to /data/larsonma/LargeRadiusJets/outputNTuplesDev_CondorSubmission_NewSamples/
#
#   --met
#     Input:  {base}_file{N}.root  (MET emulator per-job outputs)
#     Output: {base}.root
#     DIR defaults to /data/larsonma/GEPMET/outputNTuplesDev_METv2/
#
#   --ntuples
#     Input:  {base}_{NNNNNN}.root  (HERNTupler per-job outputs, 6-digit index)
#     Output: {base}.root  written into each sample subdirectory
#     DIR defaults to /data/larsonma/GEPHadronicEventReconstruction/ntuples/
#
#   --caloshowershape
#     Input:  caloShowerShape_{tag}_{NNNNNN}.root  (caloShowerShapeNTupler per-job
#             outputs, 6-digit index; {tag} is the signal/background label)
#     Output: caloShowerShape_{tag}.root  (one merged file per tag, written into
#             the same flat directory — the tag keeps samples from colliding)
#     DIR defaults to /data/larsonma/CaloShowerShapeTriggers/ntuples/
#
# Usage:
#   ./hadd_emulator_outputs.sh [OPTIONS] [DIR]
#
# Options:
#   --met           Merge MET emulator outputs instead of jet emulator outputs
#   --ntuples       Merge HERNTupler ntuples instead of emulator outputs
#   --caloshowershape  Merge caloShowerShapeNTupler per-job outputs (flat dir; one
#                   merged file per signal/background tag)
#   --dry-run       Print what would be done without running hadd
#   --force         Overwrite existing merged output files
#   --cleanup       Delete per-file inputs after a successful hadd
#   --jobs N        Run N hadd processes in parallel (default: 4)
#   --threads N     Pass -j N to hadd for multi-threaded merging (default: off)
#   --merge-jz      After per-sample merges, combine all JZ-slice outputs into
#                   one file per group (e.g. QCD_Dijet_JZ{0..9}_v3 → QCD_Dijet_JZ_v3.root)
#   --jz-only       Skip per-sample merges; only run the JZ slice combine step
#   --sample NAME   --ntuples only: restrict the merge to a single sub-directory
#                   of NTUPLE_DIR (e.g. --sample ZvvHbb_v3). Useful when one
#                   sample has been re-run and the rest do not need re-hadding.
#   --debug         Print sorted input file lists before each hadd to verify ordering
#   --version N     --ntuples only: restrict processing to sample directories ending
#                   in _v{N} (e.g. --version 4 matches ggF_HHbbbb_v4, QCD_Dijet_JZ0_v4).
#                   Also restricts --merge-jz to only combine files from _v{N} dirs.

set -euo pipefail

MODE="emulator"
EMULATOR_DIR="/data/larsonma/LargeRadiusJets/outputNTuplesDev_CondorSubmission_NewSamples/"
MET_DIR="/data/larsonma/GEPMET/outputNTuplesDev_METv2/"
NTUPLE_DIR="/data/larsonma/GEPHadronicEventReconstruction/ntuples/"
CALO_DIR="/data/larsonma/CaloShowerShapeTriggers/ntuples/"
DRY_RUN=0
FORCE=1
CLEANUP=0
PARALLEL_JOBS=4
HADD_THREADS=0
MERGE_JZ=0
JZ_ONLY=0
DEBUG=0
EXPLICIT_DIR=""
SAMPLE_FILTER=""
VERSION_FILTER=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ntuples)  MODE="ntuples" ;;
        --caloshowershape) MODE="caloshowershape" ;;
        --met)      MODE="met" ;;
        --dry-run)  DRY_RUN=1 ;;
        --force)    FORCE=1 ;;
        --cleanup)  CLEANUP=1 ;;
        --jobs)     PARALLEL_JOBS="$2"; shift ;;
        --threads)  HADD_THREADS="$2"; shift ;;
        --merge-jz) MERGE_JZ=1 ;;   # only combine 10 JZ slices → 1 file (skip per-sample merges)
        --jz-only)  JZ_ONLY=1 ;;   # only do per-sample merges for JZ dirs (no combine)
        --sample)   SAMPLE_FILTER="$2"; shift ;;
        --version)  VERSION_FILTER="$2"; shift ;;
        --debug)    DEBUG=1 ;;
        -*)         echo "Unknown option: $1"; exit 1 ;;
        *)          EXPLICIT_DIR="$1" ;;
    esac
    shift
done

if [[ -n "$SAMPLE_FILTER" && "$MODE" != "ntuples" ]]; then
    echo "[error] --sample is only valid in --ntuples mode"; exit 1
fi
if [[ -n "$VERSION_FILTER" && "$MODE" != "ntuples" ]]; then
    echo "[error] --version is only valid in --ntuples mode"; exit 1
fi

# --- hadd helper (runs in a subshell for parallel execution) ---
run_hadd() {
    local base="$1"
    local pattern="$2"   # glob suffix to find inputs, e.g. '_file*.root' or '_[0-9]*.root'
    local force="$3"
    local cleanup="$4"
    local dry="$5"
    local threads="$6"
    local debug="$7"

    local out="${base}.root"
    # Sort inputs with version sort so _file2.root < _file10.root (not lexicographic)
    mapfile -t inputs < <(printf '%s\n' ${base}${pattern} | sort -V)

    if [[ ! -e "${inputs[0]}" ]]; then
        echo "[skip] No input files match: ${base}${pattern}"
        return 0
    fi

    local n="${#inputs[@]}"

    if [[ "$debug" -eq 1 ]]; then
        echo "[debug] Sorted input order for $(basename "$out"):"
        for i in "${!inputs[@]}"; do
            printf "        [%3d] %s\n" "$i" "$(basename "${inputs[$i]}")"
        done
    fi

    if [[ -f "$out" && "$force" -eq 0 ]]; then
        echo "[skip] Already exists (--force to overwrite): $(basename "$out")"
        return 0
    fi

    echo "[hadd] $(basename "$out")  ($n files)"

    if [[ "$dry" -eq 1 ]]; then
        echo "       (dry run — not executing)"
        return 0
    fi

    local hadd_j_flag=()
    [[ "$threads" -gt 0 ]] 2>/dev/null && hadd_j_flag=(-j "$threads")

    if hadd -f "${hadd_j_flag[@]}" "$out" "${inputs[@]}" > /dev/null 2>&1; then
        echo "  OK"
        if [[ "$cleanup" -eq 1 ]]; then
            rm -f "${inputs[@]}"
            echo "  cleaned up $n input files"
        fi
    else
        echo "  [error] hadd failed for $(basename "$out")" >&2
        return 1
    fi
}

export -f run_hadd
export DEBUG

# --- Parallel job pool ---
run_parallel() {
    local pattern="$1"
    shift
    local base_list=("$@")

    local pids=()
    local errors=0

    for base in "${base_list[@]}"; do
        run_hadd "$base" "$pattern" "$FORCE" "$CLEANUP" "$DRY_RUN" "$HADD_THREADS" "$DEBUG" &
        pids+=($!)

        while [[ ${#pids[@]} -ge $PARALLEL_JOBS ]]; do
            local new_pids=()
            for pid in "${pids[@]}"; do
                if kill -0 "$pid" 2>/dev/null; then
                    new_pids+=("$pid")
                else
                    wait "$pid" || (( errors++ )) || true
                fi
            done
            pids=("${new_pids[@]}")
            [[ ${#pids[@]} -ge $PARALLEL_JOBS ]] && sleep 0.5
        done
    done

    for pid in "${pids[@]}"; do
        wait "$pid" || (( errors++ )) || true
    done

    echo ""
    if [[ $errors -eq 0 ]]; then
        echo "Done. All merges completed successfully."
    else
        echo "[warn] $errors merge(s) failed."
        return 1
    fi
}

# ===========================================================================
# Emulator / MET mode: merge {base}_file{N}.root → {base}.root in a flat directory
# ===========================================================================
if [[ "$MODE" == "emulator" || "$MODE" == "met" ]]; then
    if [[ "$MODE" == "met" ]]; then
        DIR="${EXPLICIT_DIR:-$MET_DIR}"
    else
        DIR="${EXPLICIT_DIR:-$EMULATOR_DIR}"
    fi
    DIR="${DIR%/}/"

    if [[ ! -d "$DIR" ]]; then
        echo "[error] Directory not found: $DIR"; exit 1
    fi

    declare -A bases
    while IFS= read -r -d '' f; do
        base=$(basename "$f" | sed 's/_file[0-9][0-9]*\.root$//')
        bases["${DIR}${base}"]=1
    done < <(find "$DIR" -maxdepth 1 -name '*_file[0-9]*.root' -print0 | sort -zV)

    if [[ ${#bases[@]} -eq 0 ]]; then
        echo "No per-file emulator outputs (*_file<N>.root) found in $DIR"; exit 0
    fi

    echo "Merging ${#bases[@]} algorithm configuration(s) in $DIR"
    echo ""

    mapfile -t sorted_bases < <(printf '%s\n' "${!bases[@]}" | sort)
    run_parallel "_file*.root" "${sorted_bases[@]}"

# ===========================================================================
# caloShowerShape mode: merge caloShowerShape_{tag}_{NNNNNN}.root →
# caloShowerShape_{tag}.root in one flat directory. The per-sample tag embedded in
# each filename keeps different samples (signals, background) from colliding.
# ===========================================================================
elif [[ "$MODE" == "caloshowershape" ]]; then
    DIR="${EXPLICIT_DIR:-$CALO_DIR}"
    DIR="${DIR%/}/"

    if [[ ! -d "$DIR" ]]; then
        echo "[error] Directory not found: $DIR"; exit 1
    fi

    declare -A bases
    while IFS= read -r -d '' f; do
        base=$(basename "$f" | sed 's/_[0-9][0-9][0-9][0-9][0-9][0-9]\.root$//')
        bases["${DIR}${base}"]=1
    done < <(find "$DIR" -maxdepth 1 -name '*_[0-9][0-9][0-9][0-9][0-9][0-9].root' -print0 | sort -zV)

    if [[ ${#bases[@]} -eq 0 ]]; then
        echo "No per-file caloShowerShape outputs (*_<NNNNNN>.root) found in $DIR"; exit 0
    fi

    echo "Merging ${#bases[@]} caloShowerShape sample tag(s) in $DIR"
    echo ""

    mapfile -t sorted_bases < <(printf '%s\n' "${!bases[@]}" | sort)
    run_parallel "_[0-9]*.root" "${sorted_bases[@]}"

# ===========================================================================
# Ntuple mode: merge {base}_{NNNNNN}.root → {base}.root per sample directory
# ===========================================================================
else
    DIR="${EXPLICIT_DIR:-$NTUPLE_DIR}"
    DIR="${DIR%/}/"

    if [[ ! -d "$DIR" ]]; then
        echo "[error] Directory not found: $DIR"; exit 1
    fi

    # Find all sample subdirectories (filtered to a single one if --sample given,
    # or to a version suffix if --version given)
    if [[ -n "$SAMPLE_FILTER" ]]; then
        candidate_dir="${DIR}${SAMPLE_FILTER}"
        if [[ ! -d "$candidate_dir" ]]; then
            echo "[error] --sample directory not found: $candidate_dir"; exit 1
        fi
        sample_dirs=("$candidate_dir")
        echo "Restricting --ntuples merge to single sample: $SAMPLE_FILTER"
    elif [[ -n "$VERSION_FILTER" ]]; then
        mapfile -t sample_dirs < <(find "$DIR" -mindepth 1 -maxdepth 1 -type d -name "*_v${VERSION_FILTER}" | sort)
        echo "Restricting --ntuples merge to version: v${VERSION_FILTER}"
    else
        mapfile -t sample_dirs < <(find "$DIR" -mindepth 1 -maxdepth 1 -type d | sort)
    fi

    if [[ ${#sample_dirs[@]} -eq 0 ]]; then
        echo "No subdirectories found in $DIR"; exit 0
    fi

    # When --jz-only: restrict per-sample merge to JZ slice directories only
    if [[ $JZ_ONLY -eq 1 ]]; then
        mapfile -t sample_dirs < <(printf '%s\n' "${sample_dirs[@]}" | grep -E 'JZ[0-9]')
        echo "Restricting per-sample merge to JZ slice directories only"
    fi

    all_errors=0
    if [[ $MERGE_JZ -eq 0 || $JZ_ONLY -eq 1 ]]; then
    for sample_dir in "${sample_dirs[@]}"; do
        sample_dir="${sample_dir%/}/"
        sample_name=$(basename "${sample_dir%/}")

        n_files=$(find "$sample_dir" -maxdepth 1 -name '*_[0-9][0-9][0-9][0-9][0-9][0-9].root' | wc -l)
        if [[ $n_files -eq 0 ]]; then
            echo "[skip] No per-file ntuples found in $sample_dir"
            continue
        fi

        unset bases
        declare -A bases
        while IFS= read -r -d '' f; do
            base=$(basename "$f" | sed 's/_[0-9][0-9][0-9][0-9][0-9][0-9]\.root$//')
            bases["${sample_dir}${base}"]=1
        done < <(find "$sample_dir" -maxdepth 1 -name '*_[0-9][0-9][0-9][0-9][0-9][0-9].root' -print0 | sort -z)

        echo "=== $sample_name  ($n_files files → ${#bases[@]} merged output(s)) ==="

        mapfile -t sorted_bases < <(printf '%s\n' "${!bases[@]}" | sort)
        run_parallel "_[0-9]*.root" "${sorted_bases[@]}" || (( all_errors++ )) || true
    done
    fi  # per-sample merges

    if [[ $all_errors -gt 0 ]]; then
        echo "[warn] $all_errors sample(s) had failures."
        exit 1
    fi

    # -----------------------------------------------------------------------
    # Optional: combine all JZ-slice merged outputs into one file
    # -----------------------------------------------------------------------
    if [[ $MERGE_JZ -eq 1 ]]; then
        echo ""
        echo "=== Merging JZ slices ==="

        jz_path_filter="${DIR}*"
        [[ -n "$VERSION_FILTER" ]] && jz_path_filter="${DIR}*_v${VERSION_FILTER}"
        mapfile -t jz_inputs < <(
            find "$DIR" -mindepth 2 -maxdepth 2 \
                -path "${jz_path_filter}/mc21_14TeV_jj_JZ[0-9]*_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root" | sort
        )

        if [[ ${#jz_inputs[@]} -eq 0 ]]; then
            echo "[skip] No JZ-slice merged outputs found in $DIR"
        else
            out="${DIR}mc21_14TeV_jj_JZ_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root"
            echo "[hadd-jz] $(basename "$out")  (${#jz_inputs[@]} slices)"
            for f in "${jz_inputs[@]}"; do echo "          $f"; done

            if [[ $DRY_RUN -eq 1 ]]; then
                echo "       (dry run — not executing)"
            else
                hadd_j_flag=()
                [[ "$HADD_THREADS" -gt 0 ]] 2>/dev/null && hadd_j_flag=(-j "$HADD_THREADS")

                if hadd -f "${hadd_j_flag[@]}" "$out" "${jz_inputs[@]}" > /dev/null 2>&1; then
                    echo "  OK -> $out"
                else
                    echo "  [error] hadd failed for $(basename "$out")" >&2
                    exit 1
                fi
            fi
        fi
    fi
fi
