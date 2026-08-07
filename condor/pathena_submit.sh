#!/bin/bash
# pathena_submit.sh
# Submit GEP ntuple production to the grid for all datasets.
# Launch with: ./pathena_submit.sh

# ============================================================
# Global job configuration
# ============================================================
#VERSION="v1_HHbbbb_QCD_dijet_ZvvHbb"
VERSION="v22_GEPHadronicReconstruction_gjTowers_PU200"
ATHENA_TAG="Athena,main,latest"

# ============================================================
# Pileup scenario: 200 or 140
# The <mu> value is encoded in the reconstruction (r) tag of each dataset.
# Every dataset string below uses ${RECO_TAG} in place of a hard-coded r-tag,
# so switching PILEUP re-points all samples at the matching production.
# The PU level is also baked into the output dataset name so <mu>=140 and
# <mu>=200 productions never collide.
# ============================================================
PILEUP=200   # 140 or 200
case "${PILEUP}" in
    200) RECO_TAG="r16130" ;;
    140) RECO_TAG="r16129" ;;
    *)   echo "ERROR: unsupported PILEUP=${PILEUP} (expected 140 or 200)" >&2; exit 1 ;;
esac

N_FILES_PER_JOB=5   # input files per grid job
N_JOBS=""           # max jobs per dataset    (empty = no limit)
N_FILES=1000         # max input files used    (empty = no limit; set e.g. N_FILES=5 for testing)
MAX_N_FILES=1000    # cap when N_FILES is empty and rucio count succeeds
MAX_EVENTS=""       # max events per job      (empty = all; useful for testing)

# ============================================================
# Toggle which physics processes to submit (true/false)
# ============================================================
SUBMIT_ggF_HHbbbb=false
SUBMIT_VBF_HHbbbb=false
SUBMIT_QCD_Dijet=false
# Which QCD Dijet JZ slices to submit (space-separated list of 0-9; empty = all).
# e.g. "1" to re-submit only JZ1, "0 1 2" for the first three.
#QCD_JZ_SLICES="0 2 3 4 5 6 7 8 9"
QCD_JZ_SLICES="0 8"
SUBMIT_ttbar_allhad=false
SUBMIT_ttbar_semilep=false
SUBMIT_ttbar_dilep=false
SUBMIT_Zprime=false
SUBMIT_ZvvHbb=false
SUBMIT_Zmumu=true
SUBMIT_displaced_dark_photon=false
SUBMIT_emerging_jets=false
SUBMIT_StauStauLLP=false


# Run JETM42 derivation alongside GEP ntuple production.
# Set false to run GEP only (uses run_gep_only.sh, drops DAOD output).
SUBMIT_DERIVATION=true

# ============================================================
# Helper: query number of files in a dataset via rucio.
# Returns empty string on failure (rucio not set up, dataset not found).
#   get_nfiles <full_dataset_name>
# ============================================================
get_nfiles() {
    local ds=$1
    local scope="${ds%%.*}"
    rucio info "${scope}:${ds}" 2>/dev/null | awk '/^length[[:space:]]/ {print $NF}'
}

get_total_events() {
    local ds=$1
    local scope="${ds%%.*}"
    rucio info "${scope}:${ds}" 2>/dev/null | awk '/^events[[:space:]]/ {print $NF}'
}

# ============================================================
# Helper: build optional pathena flags from config above.
# If N_FILES is empty, auto-detect via rucio and cap at MAX_N_FILES.
# Prints events/file and warns if != 100.
#   extra_opts <dataset>
# ============================================================
extra_opts() {
    local ds="${1:-}"
    local opts=""
    [[ -n "${N_JOBS}"     ]] && opts+=" --nJobs ${N_JOBS}"
    if [[ -n "${N_FILES}" ]]; then
        opts+=" --nFiles ${N_FILES}"
    elif [[ -n "${ds}" ]]; then
        local nfiles nevents
        nfiles=$(get_nfiles "${ds}")
        nevents=$(get_total_events "${ds}")
        if [[ -n "${nfiles}" && "${nfiles}" -gt 0 ]]; then
            local limit=$(( nfiles < MAX_N_FILES ? nfiles : MAX_N_FILES ))
            echo "  rucio: ${nfiles} files available, using --nFiles ${limit}" >&2
            opts+=" --nFiles ${limit}"
            if [[ -n "${nevents}" && "${nevents}" -gt 0 ]]; then
                local epf=$(( nevents / nfiles ))
                echo "  rucio: ~${epf} events/file (${nevents} total)" >&2
                if [[ "${epf}" -ne 100 ]]; then
                    echo "  WARNING: events/file = ${epf}, expected 100" >&2
                fi
            fi
        else
            echo "  rucio: could not determine file count for ${ds}, submitting without --nFiles" >&2
        fi
    fi
    [[ -n "${MAX_EVENTS}" ]] && opts+=" --maxEvents ${MAX_EVENTS}"
    echo "${opts}"
}

# Helper: submit one dataset
#   submit <AMI_string> <short_label> [extra_runconfig_opts]
submit() {
    local inds=$1
    local label=$2
    local run_opts="${3:-}"
    local outds="user.mlarson.GEPNtupleJETM42.${label}.PU${PILEUP}.${VERSION}"

    echo "------------------------------------------------------------"
    echo "  label : ${label}"
    echo "  inDS  : ${inds}"
    echo "  outDS : ${outds}"
    [[ -n "${run_opts}" ]] && echo "  opts  : ${run_opts}"
    echo "------------------------------------------------------------"

    # GEP-only mode: skip JETM42 derivation entirely.
    if ! ${SUBMIT_DERIVATION}; then
        local trf="bash run_gep_only.sh %IN${run_opts:+ ${run_opts}}"
        pathena --trf           "${trf}"                                         \
            --inDS              "${inds}"                                        \
            --memory            "4000"                                           \
            --outDS             "${outds}"                                       \
            --extOutFile        "outputGEPNtuple.root"                           \
            --nFilesPerJob      "${N_FILES_PER_JOB}"                             \
            $(extra_opts "${inds}")
        return
    fi

    local trf="bash run_both.sh %IN${run_opts:+ ${run_opts}}"

    pathena --trf           "${trf}"                                         \
        --inDS              "${inds}"                                        \
        --memory            "4000"                                           \
        --outDS             "${outds}"                                       \
        --extOutFile        "outputGEPNtuple.root,DAOD_JETM42.outputDAOD_JETM42.pool.root" \
        --nFilesPerJob      "${N_FILES_PER_JOB}"                             \
        $(extra_opts "${inds}")
}

# ============================================================
# ggF HH->bbbb
# ============================================================
if ${SUBMIT_ggF_HHbbbb}; then
    submit \
        "mc21_14TeV.603277.PhPy8EG_PDF4LHC21_HHbbbb_HLLHC_chhh1p0.recon.AOD.e8564_s4422_${RECO_TAG}" \
        "ggF_HHbbbb" \
        "-higgs"
fi

# ============================================================
# VBF HH->bbbb
# ============================================================
if ${SUBMIT_VBF_HHbbbb}; then
    submit \
        "mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_${RECO_TAG}" \
        "VBF_HHbbbb" \
        "-higgs"
fi

# ============================================================
# QCD Dijet (JZ0-JZ8 share a naming pattern; JZ9 is "incl")
# ============================================================
if ${SUBMIT_QCD_Dijet}; then
    # Default to all slices 0-9 when QCD_JZ_SLICES is empty.
    jz_slices="${QCD_JZ_SLICES:-0 1 2 3 4 5 6 7 8 9}"
    for jz in ${jz_slices}; do
        dsid=$(( 801165 + jz ))
        if [[ "${jz}" -eq 9 ]]; then
            # JZ9incl has a different suffix
            inds="mc21_14TeV.${dsid}.Py8EG_A14NNPDF23LO_jj_JZ9incl.recon.AOD.e8557_s4422_${RECO_TAG}"
        else
            inds="mc21_14TeV.${dsid}.Py8EG_A14NNPDF23LO_jj_JZ${jz}.recon.AOD.e8557_s4422_${RECO_TAG}"
        fi
        submit "${inds}" "QCD_Dijet_JZ${jz}"
    done
fi

# ============================================================
# ttbar
# ============================================================
if ${SUBMIT_ttbar_allhad}; then
    submit \
        "mc21_14TeV.601237.PhPy8EG_A14_ttbar_hdamp258p75_allhad.recon.AOD.e8557_s4422_${RECO_TAG}" \
        "ttbar_allhad"
fi

if ${SUBMIT_ttbar_semilep}; then
    submit \
        "mc21_14TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.recon.AOD.e8514_s4422_${RECO_TAG}" \
        "ttbar_semilep"
fi

if ${SUBMIT_ttbar_dilep}; then
    submit \
        "mc21_14TeV.601230.PhPy8EG_A14_ttbar_hdamp258p75_dil.recon.AOD.e8557_s4422_${RECO_TAG}" \
        "ttbar_dilep"
fi

# ============================================================
# Z' -> ttbar (flat pT)
# ============================================================
if ${SUBMIT_Zprime}; then
    submit \
        "mc21_14TeV.800301.Py8EG_A14NNPDF23LO_flatpT_Zprime_tthad.recon.AOD.e8557_s4422_${RECO_TAG}" \
        "Zprime_ttbar_allhad_flatpT"
fi

# ============================================================
# ZvvH->bb
# ============================================================
if ${SUBMIT_ZvvHbb}; then
    submit \
        "mc21_14TeV.800300.Py8EG_A14NNPDF23LO_ZvvH125_bb.recon.AOD.e8557_s4422_${RECO_TAG}" \
        "ZvvHbb" \
        "-higgs"
fi

# ============================================================
# Z -> mu mu
# ============================================================
if ${SUBMIT_Zmumu}; then
    submit \
        "mc21_14TeV.601190.PhPy8EG_AZNLO_Zmumu.recon.AOD.e8557_s4422_${RECO_TAG}" \
        "Zmumu"
fi

# ============================================================
# Displaced dark photon (HAHM Zd->Zd->4e)
# ============================================================
if ${SUBMIT_displaced_dark_photon}; then
    submit \
        "mc21_14TeV.543785.MGPy8EG_A14NN23LO_HAHM_ZdZd4e_110_30_10ns.recon.AOD.e8557_s4422_${RECO_TAG}" \
        "displaced_dark_photon"
fi

# ============================================================
# Emerging jets (Z' -> dark quarks)
# ============================================================
if ${SUBMIT_emerging_jets}; then
    submit \
        "mc21_14TeV.801966.Py8EG_Zprime2EJs_Ld20_rho40_pi10_Zp1500_l50.recon.AOD.e8532_s4422_${RECO_TAG}" \
        "emerging_jets"
fi

# ============================================================
# Stau-Stau LLP (long-lived, 10ns lifetime)
# r16129 = <mu>=140, r16130 = <mu>=200, selected via ${RECO_TAG}.
# ============================================================
if ${SUBMIT_StauStauLLP}; then
    submit \
        "mc21_14TeV.516672.MGPy8EG_A14NNPDF23LO_StauStauLLP_500_0_10ns.recon.AOD.e8557_s4422_${RECO_TAG}" \
        "StauStauLLP_500_0_10ns"
fi
