# Large-R Jet Tagging Algorithm (ATLAS GEP)

This repository contains code for the development and performance
studies of a large-R jet tagging algorithm for the ATLAS Global Event
Processor.

------------------------------------------------------------------------

## Directory Structure

### `algorithm/`

Contains the HLS-based basic (v2) & advanced (v3) algorithms, their
respective testbenches, code for reading and writing test vector files,
and a standalone C++-based emulation which is synchronized to be
bit-wise (input-output) accurate with the HLS algorithm.

### `analysis/`

Contains the NTupling code (reads in TrigGepPerf/GEPOutputReader output
& JETM42 DAOD) to generate algorithm input test vectors and ROOT
NTuples, code for generating event displays, and code for analyzing the
output of the jet tagger algorithm and performing rate estimates. Only
LRJNTupler.C, largeRJetAnalysisAndRates.C, and makeEventDisplays.C are
actively used and under-development.

### `config/`

Contains .py scripts for running/synthesizing the HLS basic & advanced
algorithms.

### `data/`

Is a placeholder directory, though contains scripts for converting from
jet tagger test vector formats --\> GEP-FW test vector format, and
generating random data in jet tagger test vector format.

### `project/`

Contains a deprecated .tcl script for running/synthesizing the HLS
algorithm.

------------------------------------------------------------------------

## Instructions for Different Scripts

NOTE: for all different scripts, most parameters are changed within the
script rather than as arguments in the terminal, and all paths are
configurable, but take default values of what the author uses on either
the U Chicago Tier 3 Analaysis Facility (login.af.uchicago.edu) or the
Miller group workstation (millerlabml01.uchicago.edu).

------------------------------------------------------------------------

## HLS Algorithm

To run the HLS algorithm or synthesis, simply do:

``` bash
python run_hls(_adv).py
```

After configuring within the script which algorithm configurations to
run over and the arguments for jet_tagger_hls(\_adv).py, which actually
calls HLS.

Use: - hls_exec = 0 for C-sim only
- hls_exec = 1 for synthesis only
- hls_exec = 2 for C-sim, synthesis & co-sim
- hls_exec = 3 for C-sim, synthesis, co-sim, & implementation

When running, it will output within config/w/ a directory which contains
the output from running HLS, which is configured based on the algorithm
parameters used for that run.

------------------------------------------------------------------------

## Standalone Emulation

The standalone emulation is run through a configuration script:

``` bash
algorithm/emulation/jetTaggerConfigLocal.sh
```

### Parameters

1.  algorithm version (2, 3)

2.  rMergeCuts (d_search): configures the distance away from the
    original seeds to search for proto-seeds to take the midpoint
    between and shift the final jet seed positions there. Typical values
    used include 0.001 (disabled), 1.5, 1.75 & 2. Values should not
    exceed 2 times the jet radius to prevent loss of E_T from the
    original jet.

3.  rSquaredCuts: jet radius squared, typical values include 1, 1.21,
    1.44

4.  nIOs: number of considered input objects to be used for reclustering
    and E_T computation of each jet. Typical values when using towers as
    input objects: 64, 128, 256. Tower multiplicity when applying
    SoftKiller PU-suppression typically does not exceed 128. The value
    used when using WTA-cone jets as input objects is 8 (10 = max jets
    per event - 2 = fixed number of seeds).

5.  nSeeds: number of seeds and consequently outputted large-R jets per
    event.

6.  signals: use (true false) to process both signals and background, or
    only (true) or (false) to process only signals or backgrounds,
    respectively.

7.  signalStrings: a list of signal strings to process, currently
    supported signal strings are: ("VBF_hh_bbbb_cvv0" "VBF_hh_bbbb_cvv1"
    "ggF_hh_bbbb" "ZvvHbb" "ttbar_had" "Zprime_ttbar")

8.  inputObjects: types of input objects used for reclustering,
    currently supported are: ("gepWTAConeCellsTowersJets"
    "gepCellsTowers"). Note that the basic algorithm can only use
    gepWTAConeCellsTowersJets. Deprecated input objects include
    CaloTopoTowers, GEPBasicClusters (still functional), and 422 Topo
    Clusters.

9.  seedObjects: type of objects used for seeding & subjet finding,
    currently supported are: ("gepWTAConeCellsTowersJets" "jFEXSRJ"
    "gFEXSRJ")

10. condor: whether to run on condor, currently unsupported and
    unnecessary since emulation runs very quickly.
11. subjetEtThreshold: minimum E_T for small-R jets to be considered as subjet candidates, is controlled separately in this script for gFEX, jFEX and WTA-cone jets.
    
12. enableEtWeightedMidpoints: whether the midpoint computation used for seed optimization is weighted by the E_T of each proto-seed used in the midpoint computation.
    
13. minEtSeedPosOptimizations and minEtSeedPosOptimizationsCuts: whether (and the value of) there is a minimum E_T threshold for small-R jets to be considered as proto-seeds for seed         position optimization (nominally is 25 GeV targetted for WTA-cone jets). 

------------------------------------------------------------------------

## LUT & Constants Generation

The emulation also relies on having already run:

``` bash
root -b
.L writeLUTsConstantsEmulation.cc
main()
```

This writes the used constants and Look Up Table (LUT) header files for
the algorithm provided the hard-coded constants for each algorithm
version. Parameters here include the number of bits allocated to
different data-types, total bits for jet tagger TOBs, & the rSquared and
rMerge cuts used.

In order to run a certain configuration of the algorithm, the constants
and LUTs must have been already generated by this script.

------------------------------------------------------------------------

## Emulation Workflow

-   jetTaggerConfigLocal.sh routes the correct constants and LUT header
    files for each algorithm configuration to be included in
    algorithm/emulation/UsedConstantsLUTs/

-   jetTaggerEmulation.cc routes the correct input NTuple based on the
    signal / background options, and sets the correct output file name
    based on the algorithm configuration processed.

------------------------------------------------------------------------

## NTupling Script

Running the standalone emulation requires having generated NTuples from
combining a GEPOutputReader NTuple and JETM42 Derivation DAOD file. In
the future, it should just take GEPOutputReader NTuple as input.

The NTupling script (LRJNTupler.C) uses constants (bit-widths, ranges)
hard-coded in analysis/analyisHelperFunctions.h. These constants must be
synchronized with what are used when running the emulation.

The script can be run through:

``` bash
root -b -l -q 'LRJNTupler.C'
```

Paths used are hard-coded in functions: - makeFileDir
- makeMemPrintFileNames
- makeOutputFileName

------------------------------------------------------------------------

## Analysis Script

Running the analysis script requires having already run the standalone
emulation, and providing the correct path to the signal and background
used.

The script allows processing multiple signals and backgrounds in one
run, and will output object property & performance plots for each signal
and background, and all files processed overlaid.

### Output Directories

-   largeRJetPlots/: plots of object properties and other misc. plots
-   rateVsEff/: performance plots evaluating algorithm by rate estimates
    and ROC curves for different selections applied to jet tagger TOBs
-   overlayJZSliceHistograms/: plots showing reconstructed object E_T
    distributions for different JZ slices (0-9) which can be used to
    evaluate energy resolution and algorithm performance compared to
    truth & offline jets
-   overlayMultipleFiles/: object property & performance plots for all
    different files looped through

Run with:

``` bash
root -b -l -q 'largeRJetAnalysisAndRates.C'
```

------------------------------------------------------------------------

## Workflow Diagram

<img width="1685" height="1028" alt="image" src="https://github.com/user-attachments/assets/a371f77c-457a-4178-8ed9-f9daa164a609" />

------------------------------------------------------------------------

## Contact

Any questions can be directed to the repository author:

mark.steven.larson@cern.ch
