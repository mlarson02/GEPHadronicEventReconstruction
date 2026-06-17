#!/usr/bin/env bash
# Loops over all requested MET emulation configurations and runs metEmulation.cc via ROOT.
# Mirrors the argument set used by condor/run_met_emulation_job.sh so the local
# entry-point stays available for simple debugging without batching/hadd logic.
set -euo pipefail

# ---- parameter grids ----
signals=(true)

signalStrings=("ZvvHbb")
#signalStrings=("ZvvHbb")

puSuppression=(true)
#puSuppression=(true)

doJetTowerOverlapRemovalBools=(true)

#jetEtThresholds=(0.0 10.0 20.0 30.0)
jetEtThresholds=(20.0)

towerEtThresholds=(2.0)

# EtaSK PU-suppressed objects (towers + WTAConeJets). The metAnalysisAndRates
# config consumes EtaSK_OR ntuples, so default to true here.
useEtaSKObjectsBools=(true)

# Scalar weights on tower / jet MET in the totalMET sum. The current analysis
# grid uses (1.0, 1.0) and (1.0, 0.4); leave just (1.0, 1.0) for fast debug.
towerScaleFactors=(1.0)
jetScaleFactors=(1.0)

# Optional: leave explicitInputPath empty to let metEmulation resolve the input
# via makeInputFileName(signalBool, signalString). Set a path here to override.
explicitInputPath=""

# Optional: fileIndex < 0 disables the _fileN suffix on the output name.
fileIndex=-1

# ---- paths ----
src_cc="/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation/metEmulation.cc"

# Clean up any stale ACLiC compiled libraries so ROOT always recompiles fresh
rm -rf metEmulation_cc*

# ---- main loops ----
for signal in "${signals[@]}"; do
  for puSupBool in "${puSuppression[@]}"; do
    for jetEtThr in "${jetEtThresholds[@]}"; do
      for doJetTowerOverlapRemoval in "${doJetTowerOverlapRemovalBools[@]}"; do
        for towerEtThr in "${towerEtThresholds[@]}"; do
          for useEtaSK in "${useEtaSKObjectsBools[@]}"; do
            for twrSF in "${towerScaleFactors[@]}"; do
              for jetSF in "${jetScaleFactors[@]}"; do

                if [[ "$signal" == true ]]; then
                  sigList=("${signalStrings[@]}")
                else
                  sigList=("ZvvHbb")  # signalString unused for background, but required by function signature
                fi

                for signalString in "${sigList[@]}"; do
                  echo "Running: signal=$signal useSK=$puSupBool signalString=$signalString jetEtThreshold=$jetEtThr doJetTowerOverlapRemoval=$doJetTowerOverlapRemoval towerEtThreshold=$towerEtThr useEtaSK=$useEtaSK explicitInputPath=\"${explicitInputPath}\" fileIndex=$fileIndex twrSF=$twrSF jetSF=$jetSF"
                  root -b -l -q "${src_cc}+(${signal}, ${puSupBool}, \"${signalString}\", ${jetEtThr}, ${doJetTowerOverlapRemoval}, ${towerEtThr}, ${useEtaSK}, \"${explicitInputPath}\", ${fileIndex}, ${twrSF}, ${jetSF})"
                done

              done
            done
          done
        done
      done
    done
  done
done

rm -rf metEmulation_cc*

echo "All configurations complete."
