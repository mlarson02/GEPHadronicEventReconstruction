#!/usr/bin/env bash
# loops over all requested parameter values (does nothing else)
set -euo pipefail
# ---- parameter grids ----
#rMergeCuts=(0.001 1.5 2.0 2.5 3.5)
#rSquaredCuts=(0.81 1.0 1.21 1.44)
#nIOs=(64 128 256 512)

#signals=(true false)
#rMergeCuts=(0.001 1.5 2.0 2.5) 
#rMergeCuts=(0.001 2) 
algoVersions=(3)
#rMergeCuts=(1.75 1.875 2.0) 
rMergeCuts=(2) 
#rMergeCuts=(1.25 1.75)
rSquaredCuts=(1.21)
nIOs=(128)
nSeeds=(2)
signals=(true false)
puSuppression=(true)
#vbfBools=(true false)
#signalStrings=("Zprime_ttbar")
#signalStrings=("VBF_hh_bbbb_cvv0" "VBF_hh_bbbb_cvv1" "ggF_hh_bbbb" "ZvvHbb" "ttbar_had" "Zprime_ttbar")
signalStrings=("ggF_hh_bbbb")
#signalStrings=("VBF_hh_bbbb_cvv1")

#signalStrings=("ttbar_had")
#signalStrings=("ggF_hh_bbbb")

#inputObjects=("gepWTAConeCellsTowersJets" "gepCellsTowers")
#seedObjects=("gepWTAConeCellsTowersJets" "jFEXSRJ" "gFEXSRJ")
inputObjects=("gepCellsTowers")
#seedObjects=("gepWTAConeCellsTowersJets")
seedObjects=("gepWTAConeCellsTowersJets")
condor=false

# ---- new algorithm behaviour parameters ----
enableEtWeightedMidpoints=(false)
minEtSeedPosOptimizations=(false true)
minEtSeedPosOptimizationCuts=(25.0)

# ---- paths ----
inputLUTFilePath="/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation/LUT_Constants_Generation/LUTs/"
constants_base="/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation/LUT_Constants_Generation/constants/"
dest_dir="/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation/UsedConstantsLUTs/"
src_cc="/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation/jetTaggerEmulation.cc"

# helper function to generate correct LUT names for each algorithm configuration
make_input_LUT_file_name() {
  local rMerge="$1" r2="$2" type="$3"
  local base="${4:-/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation/LUT_Constants_Generation/LUTs/}"
  # match C++: setprecision(4) for rMerge, setprecision(3) for R2 (no std::fixed)
  local rMergeFmt r2Fmt
  rMergeFmt=$(printf "%.4g" "$rMerge")
  r2Fmt=$(printf "%.3g" "$r2")
  printf "%sLUT_%s_rMerge_%s_R2_%s.h" "$base" "$type" "$rMergeFmt" "$r2Fmt"
}

# helper function to generate correct constants file names for each algorithm configuration
make_input_constants_file_name() {
  local rMerge="$1" nIOs="$2" nSeeds="$3" r2="$4" vers="$5"
  local base="${6:-/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation/LUT_Constants_Generation/constants/}"
  # Match C++ setprecision behavior (no fixed): 4 sig figs for rMerge, 3 for R2
  local rMergeFmt r2Fmt
  rMergeFmt=$(printf "%.4g" "$rMerge")
  r2Fmt=$(printf "%.3g" "$r2")
  printf "%sconstants_rMerge_%s_R2_%s_IOs_%s_Seeds_%s_v%s.h" \
         "$base" "$rMergeFmt" "$r2Fmt" "$nIOs" "$nSeeds" "$vers"
}

# ---- main loops ----
for rMerge in "${rMergeCuts[@]}"; do
  for r2 in "${rSquaredCuts[@]}"; do
    for ios in "${nIOs[@]}"; do
      for seeds in "${nSeeds[@]}"; do
        for algoVersion in "${algoVersions[@]}"; do
          if [[ "$algoVersion" -eq 2 && "$rMerge" != "0.001" ]]; then
            continue
          fi
          # Build filenames
          constants_file=$(make_input_constants_file_name "$rMerge" "$ios" "$seeds" "$r2" "$algoVersion" "$constants_base")

          #lut_output_path=$(make_input_LUT_file_name "$rMerge" "$r2" "deltaR2Cut" "$inputLUTFilePath")
          #lutR_output_path=$(make_input_LUT_file_name "$rMerge" "$r2" "deltaR" "$inputLUTFilePath")
          lutR_8b_output_path=$(make_input_LUT_file_name "$rMerge" "$r2" "deltaR_8b" "$inputLUTFilePath")

          echo "Config: rMerge=$rMerge r2=$r2 nIOs=$ios nSeeds=$seeds"
          echo "  constants: $constants_file"
          echo "  LUTs:"
          #echo "    $lut_output_path"
          #echo "    $lutR_output_path"
          echo "    $lutR_8b_output_path"

          echo "copying files"
          cp -f "$constants_file" "$dest_dir/constants.h"

          #cp -f "$lut_output_path" "$dest_dir/deltaR2LUT.h"
          #cp -f "$lutR_output_path" "$dest_dir/deltaRLUT.h"
          cp -f "$lutR_8b_output_path" "$dest_dir/deltaRLUT_8b.h"
  
          # Run ROOT for each (signal, PU suppression, VBF, inputObject, seedObject)
          for inputObj in "${inputObjects[@]}"; do
            for seedObj in "${seedObjects[@]}"; do
              # Subjet Et threshold depends on seed object type
              case "$seedObj" in
                gFEXSRJ)                    subjetEtThreshold=35 ;;
                jFEXSRJ)                    subjetEtThreshold=50 ;;
                gepWTAConeCellsTowersJets)  subjetEtThreshold=25 ;;
                *)                          subjetEtThreshold=25 ;;
              esac
              for etWtMid in "${enableEtWeightedMidpoints[@]}"; do
              for minEtOpt in "${minEtSeedPosOptimizations[@]}"; do
              for minEtCut in "${minEtSeedPosOptimizationCuts[@]}"; do
              for signal in "${signals[@]}"; do
                for puSupBool in "${puSuppression[@]}"; do

                  if [[ "$signal" == true ]]; then
                    sigList=("${signalStrings[@]}")
                  else
                    sigList=("VBF_hh_bbbb")   # run once; just use VBF_hh_bbbb, but it shouldn't matter which signalString is fed into algorithm FIXME allow for a null or empty string to be passed in for backgrounds?
                  fi

                  for signalString in "${sigList[@]}"; do
                    echo "Running ROOT: signal=$signal rMerge=$rMerge r2=$r2 nIOs=$ios nSeeds=$seeds puSup=$puSupBool signalString=$signalString inputObject=$inputObj seedObject=$seedObj subjetEtThreshold=$subjetEtThreshold etWtMid=$etWtMid minEtOpt=$minEtOpt minEtCut=$minEtCut"
                    root -l -b -q "${src_cc}+(${rMerge}, ${ios}, ${seeds}, ${r2}, ${signal}, ${condor}, ${puSupBool}, \"${signalString}\", \"${inputObj}\", \"${seedObj}\", ${subjetEtThreshold}, ${etWtMid}, ${minEtOpt}, ${minEtCut})"
                  done
                done
              done          # signal
              done          # minEtCut
              done          # minEtOpt
              done          # etWtMid
            done            # seedObj
          done              # inputObj
        done
      done
    done
  done
done

rm -rf jetTaggerEmulation_cc*