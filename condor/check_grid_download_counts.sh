#!/bin/bash
# Count the downloaded grid outputs per process, so you can see at a glance which
# pathena jobs have landed and which are still coming in. Walks the four production
# directories written by the rucio downloads:
#
#   GEPOutputReaderNTuples        PU200 GEP ntuples  (EXT0)
#   JETM42_DAODs                  PU200 JETM42 DAODs (EXT1)
#   GEPOutputReaderNTuples_PU140  PU140 GEP ntuples  (EXT0)
#   JETM42_DAODs_PU140            PU140 JETM42 DAODs (EXT1)
#
# The EXT0/EXT1 pair comes from the same grid job (pathena_submit.sh with
# SUBMIT_DERIVATION=true), so within one pile-up the ntuple and DAOD counts should
# match. A row where they differ is flagged with '*' — use --diff to list exactly
# which file indices are missing on which side.
#
# With --events it then opens every GEP ntuple and reports the event count (entries
# in the "ntuple" tree) per process, with a QCD dijet JZ0-JZ9 subtotal, so the sample
# summary table can be read straight off. Events are counted from the ntuples only:
# the DAODs are the same events.
#
# Usage: ./check_grid_download_counts.sh [options] [process ...]
#
#   process      only report these processes (default: all). Names are the sample
#                directory names, e.g. Zmumu, ttbar_dilep, or a QCD slice as JZ3.
#   --pu 140|200 report only one pile-up (default: both)
#   --diff       for every row whose NTuple and DAOD counts disagree, list the
#                file indices present on one side but not the other
#   --part       list the individual *.root.part files (downloads still in flight)
#   --events     also print the per-process event-count table. Needs PyROOT
#                (lsetup root) and opens every file, so it takes minutes on a full
#                production — off by default to keep a progress check instant.
#   --jobs N     files opened in parallel by --events (default: 8)
#
# Examples:
#   ./check_grid_download_counts.sh                  # full table
#   ./check_grid_download_counts.sh --pu 140         # PU140 only
#   ./check_grid_download_counts.sh --diff Zmumu JZ2 # chase two specific rows
#   ./check_grid_download_counts.sh JZ0 JZ1          # two slices only
#   ./check_grid_download_counts.sh --events         # file counts + event counts
#   ./check_grid_download_counts.sh --events JZ0 JZ1 # events for two slices only
#
# Counts are of complete *.root files. A trailing "+Np" means N *.root.part files
# are also present, i.e. rucio was still downloading when this ran.

set -e

# Parent of the four production directories. Override to point at a copy.
DATA_BASE="${GEP_DATA_BASE:-/data/larsonma/GEPHadronicEventReconstruction}"

SHOW_PU200=1
SHOW_PU140=1
SHOW_DIFF=0
SHOW_PART=0
COUNT_EVENTS=0
EVENT_JOBS=8
SELECT=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pu)      case "$2" in
                       140) SHOW_PU200=0 ;;
                       200) SHOW_PU140=0 ;;
                       *)   echo "ERROR: unsupported --pu $2 (expected 140 or 200)" >&2; exit 2 ;;
                   esac
                   shift 2 ;;
        --diff)    SHOW_DIFF=1; shift ;;
        --part)    SHOW_PART=1; shift ;;
        --events)  COUNT_EVENTS=1; shift ;;
        --jobs)    EVENT_JOBS="$2"; shift 2 ;;
        -h|--help) sed -n '2,43p' "$0"; exit 0 ;;
        -*)        echo "Unknown option: $1" >&2; exit 2 ;;
        *)         SELECT+=( "$1" ); shift ;;
    esac
done

# Process list and order follow the submission toggles in pathena_submit.sh, so this
# table reads in the same order you turned the jobs on. QCD_Dijet is expanded into its
# ten JZ slices because each slice is a separate dataset with its own job count.
PROCESSES=(
    ggF_HHbbbb
    VBF_HHbbbb
    QCD_Dijet/JZ0 QCD_Dijet/JZ1 QCD_Dijet/JZ2 QCD_Dijet/JZ3 QCD_Dijet/JZ4
    QCD_Dijet/JZ5 QCD_Dijet/JZ6 QCD_Dijet/JZ7 QCD_Dijet/JZ8 QCD_Dijet/JZ9
    ttbar_allhad
    ttbar_semilep
    ttbar_dilep
    Zprime_ttbar_allhad_flatpT
    ZvvHbb
    Zmumu
    StauStau
)

# The dataset directory (user.mlarson.<...>_EXT0) sits one level below the process
# directory, so the .root files are always at depth 2. Fixing the depth is what keeps
# the per-sample nohup.out (depth 1) out of the counts.
count_root() {
    local dir="$1"
    [ -d "${dir}" ] || { echo "-1 -1"; return; }
    local n p
    n=$(find "${dir}" -mindepth 2 -maxdepth 2 -type f -name '*.root'      2>/dev/null | wc -l)
    p=$(find "${dir}" -mindepth 2 -maxdepth 2 -type f -name '*.root.part' 2>/dev/null | wc -l)
    echo "${n} ${p}"
}

# Grid outputs carry a zero-padded six-digit index (user.mlarson.<jobid>.EXT0._000303.*),
# unique per file within a dataset, which is what makes a set difference meaningful.
idx_list() {
    find "$1" -mindepth 2 -maxdepth 2 -type f -name '*.root' -printf '%f\n' 2>/dev/null \
        | grep -oE '_[0-9]{6}\.' | tr -d '_.' | sort -u
}

fmt_cell() {
    local n="$1" p="$2"
    if [ "${n}" -lt 0 ]; then printf -- "--"; return; fi
    if [ "${p}" -gt 0 ]; then printf "%d+%dp" "${n}" "${p}"; else printf "%d" "${n}"; fi
}

# A row is selected if no filter was given, or if any filter term matches its name.
# Matching on the JZ leaf (not the full QCD_Dijet/JZ3 path) is what lets you type "JZ3".
selected() {
    local proc="$1" leaf="${1##*/}" term
    [ ${#SELECT[@]} -eq 0 ] && return 0
    for term in "${SELECT[@]}"; do
        [[ "${proc}" == "${term}" || "${leaf}" == "${term}" ]] && return 0
    done
    return 1
}

# Resolve the filter before anything is printed, so an unmatched name is a clean error
# rather than an empty table with a header.
ROWS=()
for PROC in "${PROCESSES[@]}"; do
    selected "${PROC}" && ROWS+=( "${PROC}" )
done
if [ ${#ROWS[@]} -eq 0 ]; then
    echo "No process matched: ${SELECT[*]}" >&2
    exit 1
fi

# Header: only the requested pile-ups get columns.
HEADER=$(printf "%-28s" "PROCESS")
[ "${SHOW_PU200}" -eq 1 ] && HEADER+=$(printf " %14s %14s" "PU200 NTuple" "PU200 DAOD")
[ "${SHOW_PU140}" -eq 1 ] && HEADER+=$(printf " %14s %14s" "PU140 NTuple" "PU140 DAOD")
echo ""
echo "${DATA_BASE}"
echo ""
echo "${HEADER}"
printf '%*s\n' "${#HEADER}" '' | tr ' ' '-'

TOT_NT200=0; TOT_DA200=0; TOT_NT140=0; TOT_DA140=0
TOT_P200=0;  TOT_P140=0
MISMATCHES=()

for PROC in "${ROWS[@]}"; do
    DIR_NT200="${DATA_BASE}/GEPOutputReaderNTuples/${PROC}"
    DIR_DA200="${DATA_BASE}/JETM42_DAODs/${PROC}"
    DIR_NT140="${DATA_BASE}/GEPOutputReaderNTuples_PU140/${PROC}"
    DIR_DA140="${DATA_BASE}/JETM42_DAODs_PU140/${PROC}"

    read -r N_NT200 P_NT200 < <(count_root "${DIR_NT200}")
    read -r N_DA200 P_DA200 < <(count_root "${DIR_DA200}")
    read -r N_NT140 P_NT140 < <(count_root "${DIR_NT140}")
    read -r N_DA140 P_DA140 < <(count_root "${DIR_DA140}")

    # Only compare when both sides of a pair exist; a process that was never run at one
    # pile-up (StauStau at PU140) is absent, not mismatched.
    FLAG=""
    if [ "${SHOW_PU200}" -eq 1 ] && [ "${N_NT200}" -ge 0 ] && [ "${N_DA200}" -ge 0 ] \
       && [ "${N_NT200}" -ne "${N_DA200}" ]; then
        FLAG="*"; MISMATCHES+=( "200 ${PROC}" )
    fi
    if [ "${SHOW_PU140}" -eq 1 ] && [ "${N_NT140}" -ge 0 ] && [ "${N_DA140}" -ge 0 ] \
       && [ "${N_NT140}" -ne "${N_DA140}" ]; then
        FLAG="*"; MISMATCHES+=( "140 ${PROC}" )
    fi

    ROW=$(printf "%-28s" "${PROC}${FLAG}")
    if [ "${SHOW_PU200}" -eq 1 ]; then
        ROW+=$(printf " %14s %14s" "$(fmt_cell "${N_NT200}" "${P_NT200}")" \
                                   "$(fmt_cell "${N_DA200}" "${P_DA200}")")
        [ "${N_NT200}" -gt 0 ] && TOT_NT200=$((TOT_NT200 + N_NT200))
        [ "${N_DA200}" -gt 0 ] && TOT_DA200=$((TOT_DA200 + N_DA200))
        [ "${P_NT200}" -gt 0 ] && TOT_P200=$((TOT_P200 + P_NT200))
        [ "${P_DA200}" -gt 0 ] && TOT_P200=$((TOT_P200 + P_DA200))
    fi
    if [ "${SHOW_PU140}" -eq 1 ]; then
        ROW+=$(printf " %14s %14s" "$(fmt_cell "${N_NT140}" "${P_NT140}")" \
                                   "$(fmt_cell "${N_DA140}" "${P_DA140}")")
        [ "${N_NT140}" -gt 0 ] && TOT_NT140=$((TOT_NT140 + N_NT140))
        [ "${N_DA140}" -gt 0 ] && TOT_DA140=$((TOT_DA140 + N_DA140))
        [ "${P_NT140}" -gt 0 ] && TOT_P140=$((TOT_P140 + P_NT140))
        [ "${P_DA140}" -gt 0 ] && TOT_P140=$((TOT_P140 + P_DA140))
    fi
    echo "${ROW}"
done

printf '%*s\n' "${#HEADER}" '' | tr ' ' '-'
TOTAL=$(printf "%-28s" "TOTAL")
[ "${SHOW_PU200}" -eq 1 ] && TOTAL+=$(printf " %14s %14s" "${TOT_NT200}" "${TOT_DA200}")
[ "${SHOW_PU140}" -eq 1 ] && TOTAL+=$(printf " %14s %14s" "${TOT_NT140}" "${TOT_DA140}")
echo "${TOTAL}"

echo ""
echo "  N+Mp  = N complete .root files, M .root.part still downloading"
echo "  *     = NTuple and DAOD counts differ for this process (--diff to list them)"

IN_FLIGHT=$((TOT_P200 + TOT_P140))
if [ "${IN_FLIGHT}" -gt 0 ]; then
    echo ""
    echo "${IN_FLIGHT} partial file(s) present — a download was in progress when this ran."
    if [ "${SHOW_PART}" -eq 1 ]; then
        for D in GEPOutputReaderNTuples JETM42_DAODs GEPOutputReaderNTuples_PU140 JETM42_DAODs_PU140; do
            find "${DATA_BASE}/${D}" -type f -name '*.root.part' 2>/dev/null | sort
        done
    fi
fi

# The two-way difference matters: the earlier PU140 JZ2 case had indices missing from
# each side, which a one-way "what is the DAOD missing" check would have half-reported.
if [ "${SHOW_DIFF}" -eq 1 ] && [ ${#MISMATCHES[@]} -gt 0 ]; then
    for ENTRY in "${MISMATCHES[@]}"; do
        PU="${ENTRY%% *}"
        PROC="${ENTRY#* }"
        case "${PU}" in
            200) NT_DIR="${DATA_BASE}/GEPOutputReaderNTuples/${PROC}"
                 DA_DIR="${DATA_BASE}/JETM42_DAODs/${PROC}" ;;
            140) NT_DIR="${DATA_BASE}/GEPOutputReaderNTuples_PU140/${PROC}"
                 DA_DIR="${DATA_BASE}/JETM42_DAODs_PU140/${PROC}" ;;
        esac
        echo ""
        echo "### PU${PU} ${PROC}"
        ONLY_NT=$(comm -23 <(idx_list "${NT_DIR}") <(idx_list "${DA_DIR}") | tr '\n' ' ')
        ONLY_DA=$(comm -13 <(idx_list "${NT_DIR}") <(idx_list "${DA_DIR}") | tr '\n' ' ')
        # '|| true': with an empty difference the && chain returns 1, and as the last
        # command in the loop body that would trip set -e and skip the event table.
        [ -n "${ONLY_NT}" ] && echo "  missing DAOD   (ntuple only): ${ONLY_NT}" || true
        [ -n "${ONLY_DA}" ] && echo "  missing ntuple (DAOD only):   ${ONLY_DA}" || true
    done
fi

# ---------------------------------------------------------------------------
# Event counts (--events)
# ---------------------------------------------------------------------------
# Everything above is a stat() walk over filenames; this part opens every ntuple and
# reads its tree header, so it is the slow half and is opt-in. Counted from the EXT0
# ntuples only: the EXT1 DAODs are derived from the same events, so counting both
# sides would print the same number twice.
#
# Early return rather than wrapping the rest in an if, so the default path ends here.
if [ "${COUNT_EVENTS}" -eq 0 ]; then
    exit 0
fi

# Plain ifs, not `[ ... ] && VAR=...`: under set -e a false test there is a failing
# command and would abort the script whenever a pile-up was deselected.
PU_CSV=""
if [ "${SHOW_PU200}" -eq 1 ]; then PU_CSV="200"; fi
if [ "${SHOW_PU140}" -eq 1 ]; then PU_CSV="${PU_CSV:+${PU_CSV},}140"; fi

# PyROOT comes from the ATLAS release rather than the system python, and which of
# python/python3 has it depends on the release, so probe for an interpreter that can
# actually import ROOT instead of assuming one.
PYBIN=""
for CAND in python python3; do
    command -v "${CAND}" >/dev/null 2>&1 || continue
    if "${CAND}" -c "import ROOT" >/dev/null 2>&1; then PYBIN="${CAND}"; break; fi
done
if [ -z "${PYBIN}" ]; then
    echo ""
    echo "Event counts skipped: PyROOT is not importable. Set up ROOT first, e.g." >&2
    echo "  source \${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh && lsetup root" >&2
    exit 0
fi

"${PYBIN}" - "${DATA_BASE}" "${PU_CSV}" "${EVENT_JOBS}" "${ROWS[@]}" <<'PYEOF'
"""Per-process event counts, read from the GEP ntuples (EXT0)."""
import os
import sys
import glob
from multiprocessing import Pool

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

BASE, PU_CSV, NJOBS = sys.argv[1], sys.argv[2], max(1, int(sys.argv[3]))
PUS = [p for p in PU_CSV.split(",") if p]
PROCS = sys.argv[4:]

# Same mapping as the bash half; only the ntuple (EXT0) side is opened.
NT_DIR = {"200": "GEPOutputReaderNTuples", "140": "GEPOutputReaderNTuples_PU140"}
TREE = "ntuple"          # tree written by GepOutputReader, one entry per event


def count(path):
    """Entries in one file, or -1 if it could not be read.

    Unreadable files are reported separately rather than counted as zero, so a
    truncated download cannot silently deflate a process total into looking complete.
    """
    try:
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            return -1
        t = f.Get(TREE)
        n = int(t.GetEntries()) if isinstance(t, ROOT.TTree) else -1
        f.Close()
        return n
    except Exception:
        return -1


# One work item per file rather than per process: the slices are wildly uneven
# (JZ0 has 500 files, StauStau 17), so per-process chunks would leave workers idle.
work = []
for proc in PROCS:
    for pu in PUS:
        pattern = os.path.join(BASE, NT_DIR[pu], proc, "*", "*.root")
        for path in sorted(glob.glob(pattern)):
            work.append((proc, pu, path))

# (proc, pu) -> [events, files_read, files_bad]
tally = dict(((p, u), [0, 0, 0]) for p in PROCS for u in PUS)

if work:
    sys.stderr.write("\nCounting events in %d ntuple file(s) with %d worker(s)...\n"
                     % (len(work), NJOBS))
    sys.stderr.flush()
    pool = Pool(NJOBS)
    try:
        # imap preserves input order, so results line up with `work` element-wise.
        counts = pool.imap(count, [w[2] for w in work], chunksize=4)
        for (proc, pu, _), n in zip(work, counts):
            rec = tally[(proc, pu)]
            if n < 0:
                rec[2] += 1
            else:
                rec[0] += n
                rec[1] += 1
    finally:
        pool.close()
        pool.join()

WIDTH = 28


def cell(proc, pu):
    events, nfiles, nbad = tally[(proc, pu)]
    if nfiles == 0 and nbad == 0:
        return "--"                      # process not produced at this pile-up
    text = "{:,}".format(events)
    if nbad:
        text += " (%db)" % nbad
    return text


header = "%-*s" % (WIDTH, "PROCESS")
for pu in PUS:
    header += " %18s" % ("PU%s Events" % pu)
rule = "-" * len(header)

print("")
print(header)
print(rule)

qcd = [p for p in PROCS if p.startswith("QCD_Dijet/")]
totals = dict((pu, [0, 0]) for pu in PUS)

for proc in PROCS:
    row = "%-*s" % (WIDTH, proc)
    for pu in PUS:
        row += " %18s" % cell(proc, pu)
        totals[pu][0] += tally[(proc, pu)][0]
        totals[pu][1] += tally[(proc, pu)][2]
    print(row)
    # Subtotal sits directly under the slices it sums, which is the row the sample
    # summary table quotes as a single "QCD dijet (JZ0-JZ9)" line.
    if qcd and proc == qcd[-1]:
        sub = "%-*s" % (WIDTH, "  QCD dijet (JZ0-JZ9)")
        for pu in PUS:
            events = sum(tally[(q, pu)][0] for q in qcd)
            nfiles = sum(tally[(q, pu)][1] for q in qcd)
            sub += " %18s" % ("{:,}".format(events) if nfiles else "--")
        print(sub)

print(rule)
total = "%-*s" % (WIDTH, "TOTAL")
for pu in PUS:
    total += " %18s" % "{:,}".format(totals[pu][0])
print(total)

nbad = sum(totals[pu][1] for pu in PUS)
if nbad:
    print("")
    print("  (Nb) = N file(s) could not be opened or had no '%s' tree; their events "
          "are NOT in the totals." % TREE)
PYEOF
