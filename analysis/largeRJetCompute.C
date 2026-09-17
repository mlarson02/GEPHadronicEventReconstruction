// Compute stage of the largeRJetAnalysisAndRates split (largeRJetSplitPlan.md).
//
// Runs ONE file pair's two raw background event loops over a slice of the global entry
// index and writes the filled histograms + per-event cache into
//   /data/larsonma/LargeRadiusJets/largeRJetIntermediateState/<taggerBase>__<chunk>.root
//
// Usage (after the usual ATLAS env + lsetup root):
//   root -b -l -q 'largeRJetCompute.C(<pairIndex>, <jobIndex>, <nJobs>[, "<pairsFile>"])'
//
//   pairIndex  index into the pair list, selecting which pair this job computes. State
//              contents carry NO pair index, so later plot runs may use any subset or
//              ordering of pairs; pairs sharing a background config share one state file
//              (submit_largeRJetCompute.py submits one job set per distinct config).
//   jobIndex   which 1/nJobs slice of the background entries this job processes.
//   nJobs      how many slices the pair is split into. (0, 1) = one whole-sample job.
//   pairsFile  generated pair-list file (condor/submit_largeRJetCompute.py writes it —
//              that grid is where studies are configured). Default:
//              largeRJetFilePairs.txt in this directory. Pass "" to fall back to the
//              hardcoded lists in largeRJetAnalysisCore.h.
//
// After all nJobs chunks of a pair finish:
//   hadd -f <stateDir>/<taggerBase>__state.root <stateDir>/<taggerBase>__job*.root
// (see hadd_largeRJetState.sh) then run largeRJetPlot.C.
#include "largeRJetAnalysisCore.h"

void largeRJetCompute(int pairIndex = 0, int jobIndex = 0, int nJobs = 1,
                      const char* pairsFile = "largeRJetFilePairs.txt") {
    if (nJobs < 1 || jobIndex < 0 || jobIndex >= nJobs) {
        std::cerr << "largeRJetCompute: need 0 <= jobIndex < nJobs (got " << jobIndex
                  << " of " << nJobs << ")\n";
        gSystem->Exit(2);
    }
    if (kMaxEventsPerSlice >= 0 && nJobs > 1) {
        // The per-JZ event caps count within each job, so chunked jobs would apply the cap
        // per chunk instead of per sample — silently different physics. Refuse.
        std::cerr << "largeRJetCompute: kMaxEventsPerSlice >= 0 is a whole-sample debug cap"
                  << " and is not compatible with chunked jobs\n";
        gSystem->Exit(2);
    }
    gLRJState.stage      = LRJStage::kCompute;
    lrj_file_pairs_path_ = pairsFile ? pairsFile : "";
    lrj_only_pair_index_ = pairIndex;
    lrj_job_index_       = jobIndex;
    lrj_n_jobs_          = nJobs;
    gLRJState.chunkTag   = Form("job%dof%d", jobIndex, nJobs);
    largeRJetRun();
    gSystem->Exit(0);
}
