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
#     (--pu 140: /data/larsonma/GEPHadronicEventReconstruction/ntuples_PU140/)
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
#   --force         Overwrite existing merged output files (this is the DEFAULT)
#   --skip-existing Merge only the bases whose {base}.root does not exist yet, leaving
#                   already-merged configurations untouched. The fast path when the output
#                   directory holds many configurations and only a few were just produced.
#                   Note it skips on existence alone: if more per-job files landed after a
#                   merge was made, that merge is stale and stays stale — the skip message
#                   flags this case, and --force is how you rebuild it.
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
#   --pu N          Pileup scenario, 140 or 200. The reconstruction tag in the file
#                   name is what distinguishes them: r16130 = PU200, r16129 = PU140.
#                   --ntuples: selects the ntuples_PU140 tree and the r16129 JZ files.
#                   emulator/MET: restricts the merge to that pileup's outputs, which
#                   share one output directory with the other pileup's.
#                   Omit to merge everything found, regardless of pileup (the default,
#                   and correct as long as only one pileup is present).

set -euo pipefail

# This script only orchestrates hadd; it does no ROOT setup of its own. Without a
# ROOT environment every merge fails identically and instantly, which looks like a
# data problem but is not — so check once, up front, with a message that says what
# to do about it.
if ! command -v hadd > /dev/null 2>&1; then
    echo "[error] hadd not found on PATH: this script needs a ROOT environment." >&2
    echo "        export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase" >&2
    echo "        source \${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh" >&2
    echo "        lsetup root" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MODE="emulator"
EMULATOR_DIR="/data/larsonma/LargeRadiusJets/outputNTuplesDev_gjTowerSamples/"
MET_DIR="/data/larsonma/GEPMET/outputNTuplesDev_METv2/"
NTUPLE_DIR="/data/larsonma/GEPHadronicEventReconstruction/ntuples/"
NTUPLE_DIR_PU140="/data/larsonma/GEPHadronicEventReconstruction/ntuples_PU140/"
CALO_DIR="/data/larsonma/CaloShowerShapeTriggers/ntuples/"
DRY_RUN=0
# Overwrite merged outputs by default; --skip-existing is what turns this off.
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
PU=""   # empty = merge every pileup scenario found

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ntuples)  MODE="ntuples" ;;
        --caloshowershape) MODE="caloshowershape" ;;
        --met)      MODE="met" ;;
        --dry-run)  DRY_RUN=1 ;;
        --force)    FORCE=1 ;;
        --skip-existing) FORCE=0 ;;
        --cleanup)  CLEANUP=1 ;;
        --jobs)     PARALLEL_JOBS="$2"; shift ;;
        --threads)  HADD_THREADS="$2"; shift ;;
        --merge-jz) MERGE_JZ=1 ;;   # only combine 10 JZ slices → 1 file (skip per-sample merges)
        --jz-only)  JZ_ONLY=1 ;;   # only do per-sample merges for JZ dirs (no combine)
        --sample)   SAMPLE_FILTER="$2"; shift ;;
        --version)  VERSION_FILTER="$2"; shift ;;
        --pu)       PU="$2"; shift ;;
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
if [[ -n "$PU" && "$PU" != "140" && "$PU" != "200" ]]; then
    echo "[error] --pu must be 140 or 200 (got: $PU)"; exit 1
fi

# The reconstruction tag in the file name encodes the pileup scenario: PU200 samples
# are r16130, PU140 samples r16129 (stamped by HERNTupler and carried through the
# emulation output names). It is the only differentiator for the emulator/MET outputs,
# which land in the same directory for both scenarios.
pu_reco_tag() {
    case "$1" in
        140) echo "r16129" ;;
        *)   echo "r16130" ;;
    esac
}
RECO_TAG="$(pu_reco_tag "${PU:-200}")"

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

    # --skip-existing: leave an already-merged configuration alone. Checked before the input
    # glob so a skipped base costs a single stat rather than a glob plus a version sort —
    # which is the whole point when the directory holds hundreds of finished configurations
    # and only a handful were just produced.
    if [[ -f "$out" && "$force" -eq 0 ]]; then
        # An output older than one of its inputs means per-job files landed after it was
        # merged, so it is a partial merge. Still skipped, since that is what was asked for,
        # but said out loud: silently keeping a stale merge is the one way this option bites.
        local newer
        newer=$(find "$(dirname "$out")" -maxdepth 1 \
                     -name "$(basename "$base")${pattern}" -newer "$out" -print -quit 2>/dev/null || true)
        if [[ -n "$newer" ]]; then
            echo "[skip] Already exists but INPUTS ARE NEWER — stale, use --force: $(basename "$out")"
        else
            echo "[skip] Already exists (--force to overwrite): $(basename "$out")"
        fi
        return 0
    fi

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

    echo "[hadd] $(basename "$out")  ($n files)"

    if [[ "$dry" -eq 1 ]]; then
        echo "       (dry run — not executing)"
        return 0
    fi

    local hadd_j_flag=()
    [[ "$threads" -gt 0 ]] 2>/dev/null && hadd_j_flag=(-j "$threads")

    # Keep hadd's output: a swallowed error message turns every failure mode
    # (missing ROOT, corrupt input, exceeded fgMaxTreeSize) into the same
    # uninformative "hadd failed" line.
    local log="${out}.haddlog"
    if hadd -f "${hadd_j_flag[@]}" "$out" "${inputs[@]}" > "$log" 2>&1; then
        echo "  OK"
        rm -f "$log"
        if [[ "$cleanup" -eq 1 ]]; then
            rm -f "${inputs[@]}"
            echo "  cleaned up $n input files"
        fi
    else
        echo "  [error] hadd failed for $(basename "$out")  (log: $log)" >&2
        tail -n 5 "$log" | sed 's/^/          /' >&2
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

    # PU200 and PU140 outputs share this directory, so restrict by reconstruction
    # tag when --pu is given. Without it every base found is merged, each against
    # its own _file<N> inputs.
    [[ -n "$PU" ]] && echo "Restricting merge to PU${PU} outputs (${RECO_TAG})"

    declare -A bases
    while IFS= read -r -d '' f; do
        base=$(basename "$f" | sed 's/_file[0-9][0-9]*\.root$//')
        [[ -n "$PU" && "$base" != *"_${RECO_TAG}_"* ]] && continue
        bases["${DIR}${base}"]=1
    done < <(find "$DIR" -maxdepth 1 -name '*_file[0-9]*.root' -print0 | sort -zV)

    if [[ ${#bases[@]} -eq 0 ]]; then
        echo "No per-file emulator outputs (*_file<N>.root${PU:+, PU$PU}) found in $DIR"; exit 0
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
    if [[ "$PU" == "140" ]]; then
        DIR="${EXPLICIT_DIR:-$NTUPLE_DIR_PU140}"
    else
        DIR="${EXPLICIT_DIR:-$NTUPLE_DIR}"
    fi
    DIR="${DIR%/}/"

    if [[ ! -d "$DIR" ]]; then
        echo "[error] Directory not found: $DIR"; exit 1
    fi

    # An explicitly given PU140 tree without --pu would otherwise look for r16130
    # JZ slice files and find nothing, so infer the tag from the directory name.
    if [[ -z "$PU" && "$DIR" == *PU140* ]]; then
        RECO_TAG="$(pu_reco_tag 140)"
        echo "[pileup] $DIR looks like a PU140 tree — using $RECO_TAG for the JZ merge"
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
                -path "${jz_path_filter}/mc21_14TeV_jj_JZ[0-9]*_e8557_s4422_${RECO_TAG}_DAOD_NTUPLE_GEP.root" | sort
        )

        if [[ ${#jz_inputs[@]} -eq 0 ]]; then
            echo "[skip] No JZ-slice merged outputs ($RECO_TAG) found in $DIR"
        else
            out="${DIR}mc21_14TeV_jj_JZ_e8557_s4422_${RECO_TAG}_DAOD_NTUPLE_GEP.root"
            echo "[hadd-jz] $(basename "$out")  (${#jz_inputs[@]} slices)"
            for f in "${jz_inputs[@]}"; do echo "          $f"; done

            if [[ $DRY_RUN -eq 1 ]]; then
                echo "       (dry run — not executing)"
            else
                # The combined JZ ntuple exceeds ROOT's default 100 GB TTree::fgMaxTreeSize.
                # At that point plain hadd calls TTree::ChangeFile(), swaps the output for
                # {out}_1.root mid-merge and aborts, leaving a truncated file whose later
                # trees are silently missing. haddBigTree.C raises the limit first. It is
                # single-threaded on purpose: hadd -j forks, and the raised limit would not
                # propagate into the child processes.
                jz_list=$(mktemp "${TMPDIR:-/tmp}/jz_inputs.XXXXXX")
                printf '%s\n' "${jz_inputs[@]}" > "$jz_list"
                rm -f "${out%.root}_1.root"

                if root -b -l -q "${SCRIPT_DIR}/haddBigTree.C(\"${out}\",\"${jz_list}\")"; then
                    echo "  OK -> $out"
                    rm -f "$jz_list"
                else
                    echo "  [error] haddBigTree failed for $(basename "$out")" >&2
                    echo "          input list kept for debugging: $jz_list" >&2
                    exit 1
                fi
            fi
        fi
    fi
fi
