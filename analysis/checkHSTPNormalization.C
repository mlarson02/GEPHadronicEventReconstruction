// checkHSTPNormalization.C
//
// Standalone diagnostic for the background rate normalization. Reads only the eventInfoTree
// of each JZ-slice ntuple, so it runs in seconds and needs no re-ntupling.
//
// Part 1 — measurement. Per slice:
//   * the number of ntupled events,
//   * sum(w) over all events            -> the rate the slice carries before the HSTP filter,
//   * sum(w) over passHSTP events       -> what survives the filter in the analysis macros,
//   * the weighted HSTP pass fraction,
//   * the largest single-event weight   -> how concentrated the tail of a slice is,
//   * sum(mcEventWeight) vs sumOfWeightsForSample -> the fraction of the generated sample that
//     was actually ntupled. If this is well below 1, the absolute rate is scaled down by it.
//
// Part 2 — the global post-HSTP scale used by metAnalysisAndRates.C, which rescales the rate
// histograms so the rate at threshold 0 is the colliding-bunch crossing rate. This reports the
// factor and what it does to each slice's share of the total, against the slice's true
// xsec x filterEff share.
//
// Run with:
//   root -l -b -q 'checkHSTPNormalization.C("/data/larsonma/GEPHadronicEventReconstruction/ntuples")'
//
// Measured on the v4 PU200 ntuples (2026-08-12): sum(w) over all slices = 30 MHz as designed,
// with JZ0 alone carrying 28.6 MHz (95.4%) and JZ1 1.36 MHz (4.5%) — but only 323 kHz survives
// the HSTP filter, because JZ0 keeps 0.064% of its weight and JZ1 21.5%. Note JZ1 passes 63% of
// its events but only 21.5% of its weight, i.e. the filter preferentially removes the heavy
// events, so the surviving sample is not a representative subset of the slice.

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdio>

// ---------------------------------------------------------------------------
// MIRROR of the normalization constants in HERNTupler.C — keep the two in sync.
// ---------------------------------------------------------------------------
namespace hstpRef {
    const unsigned int nJZ = 10;
    const double targetRate      = 30.0e6;   // HERNTupler targetRate (rate mode) [Hz]
    const double crossingRateHz  = 30.9e6;   // colliding-bunch crossing rate, the MET scale target
    const double crossSection[nJZ] = {0.07893, 0.09679, 0.0026805, 0.000029984, 2.972e-7,
                                      5.5384e-09, 3.2616e-10, 2.1734e-11, 9.2995e-13, 3.4519e-14};
    const double filterEff[nJZ]    = {0.9716436, 0.03777559, 0.01136654, 0.01367042, 0.01628158,
                                      0.01905588, 0.01352844, 0.01764909, 0.01887484, 0.02827565};

    // Slice share of targetRate: what sum(w) over ALL events of the slice should be.
    inline double sliceShare(unsigned int jz) {
        double sigmaRef = 0.0;
        for (unsigned int i = 0; i < nJZ; ++i) sigmaRef += crossSection[i] * filterEff[i];
        return crossSection[jz] * filterEff[jz] * (targetRate / sigmaRef);
    }
}

void checkHSTPNormalization(const char* ntupleRoot = "/data/larsonma/GEPHadronicEventReconstruction/ntuples",
                            const char* versionTag = "v4",
                            unsigned int nSlices = 10) {
    if (nSlices > hstpRef::nJZ) nSlices = hstpRef::nJZ;

    // ----------------------------------------------------------------- Part 1
    printf("\n  === Measured normalization ===\n");
    printf("\n  JZ   NEvents   NEvt HSTP    Sum(w) [Hz]   Sum(w|HSTP) [Hz]   HSTP keep(w)   HSTP keep(N)   Sum(mcw)     sumOfWeights   ntupled frac   max w [Hz]\n");
    printf("  ---------------------------------------------------------------------------------------------------------------------------------------------\n");

    double totW = 0.0, totWHSTP = 0.0;
    long   totN = 0,   totNHSTP = 0;
    std::vector<double> keepFractions(nSlices, 0.0);
    std::vector<double> sumWs(nSlices, 0.0), sumWHSTPs(nSlices, 0.0);
    std::vector<bool>   haveSlice(nSlices, false);

    for (unsigned int jz = 0; jz < nSlices; ++jz) {
        TString path = TString::Format("%s/QCD_Dijet_JZ%u_%s/mc21_14TeV_jj_JZ%u_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root",
                                       ntupleRoot, jz, versionTag, jz);
        TFile* f = TFile::Open(path, "READ");
        if (!f || f->IsZombie()) { printf("  %2u   <cannot open %s>\n", jz, path.Data()); continue; }
        TTree* t = (TTree*)f->Get("eventInfoTree");
        if (!t) { printf("  %2u   <no eventInfoTree in %s>\n", jz, path.Data()); f->Close(); continue; }

        std::vector<double>* eventWeights   = nullptr;
        std::vector<double>* mcEventWeights = nullptr;
        double sumOfWeightsForSample = 0.0;
        bool   passHSTP = true;
        int    sampleJZSlice = -1;
        t->SetBranchAddress("eventWeights",          &eventWeights);
        t->SetBranchAddress("mcEventWeight",         &mcEventWeights);
        t->SetBranchAddress("sumOfWeightsForSample", &sumOfWeightsForSample);
        t->SetBranchAddress("passHSTP",              &passHSTP);
        t->SetBranchAddress("sampleJZSlice",         &sampleJZSlice);

        double sumW = 0.0, sumWHSTP = 0.0, sumMcW = 0.0, maxW = 0.0, sowSample = 0.0;
        long   n = 0, nHSTP = 0, nSliceMismatch = 0, nEmptyWeights = 0;
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            // Events with an empty eventWeights vector exist and would throw on at(0).
            if (!eventWeights || eventWeights->empty()) { nEmptyWeights++; continue; }
            double w = eventWeights->at(0);
            sumW += w; n++;
            if (w > maxW) maxW = w;
            if (mcEventWeights && !mcEventWeights->empty()) sumMcW += mcEventWeights->at(0);
            if (passHSTP) { sumWHSTP += w; nHSTP++; }
            if (sampleJZSlice != (int)jz) nSliceMismatch++;
            sowSample = sumOfWeightsForSample;
        }
        double keepW = sumW > 0.0 ? sumWHSTP / sumW : 0.0;
        double keepN = n > 0 ? (double)nHSTP / (double)n : 0.0;
        double ntupledFrac = sowSample > 0.0 ? sumMcW / sowSample : 0.0;
        keepFractions[jz] = keepW;
        sumWs[jz] = sumW; sumWHSTPs[jz] = sumWHSTP; haveSlice[jz] = (n > 0);
        printf("  %2u %9ld %11ld %14.5g %18.5g %14.5f %14.5f %12.5g %14.5g %14.5f %12.5g\n",
               jz, n, nHSTP, sumW, sumWHSTP, keepW, keepN, sumMcW, sowSample, ntupledFrac, maxW);
        if (nSliceMismatch > 0)
            printf("       WARNING: %ld events have sampleJZSlice != %u\n", nSliceMismatch, jz);
        if (nEmptyWeights > 0)
            printf("       WARNING: %ld events have an empty eventWeights vector and were skipped\n", nEmptyWeights);
        totW += sumW; totWHSTP += sumWHSTP; totN += n; totNHSTP += nHSTP;
        f->Close();
    }

    printf("  ---------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("  ALL %8ld %11ld %14.5g %18.5g %14.5f %14.5f\n",
           totN, totNHSTP, totW, totWHSTP,
           totW > 0.0 ? totWHSTP / totW : 0.0,
           totN > 0 ? (double)totNHSTP / (double)totN : 0.0);
    printf("\n  Total rate before HSTP = %.5g MHz   (targetRate in HERNTupler.C is %.4g MHz)\n",
           totW / 1e6, hstpRef::targetRate / 1e6);
    printf("  Total rate after  HSTP = %.5g kHz\n\n", totWHSTP / 1e3);
    printf("  Where HSTP keep(N) and keep(w) differ, the filter is removing events of atypical\n");
    printf("  weight, so the survivors are not a representative subset of the slice.\n\n");

    // ----------------------------------------------------------------- Part 2
    // The global scale metAnalysisAndRates.C applies to its rate histograms: one factor putting
    // the rate at threshold 0 at the crossing rate. Relative weights are untouched, but the
    // slice mixture of the surviving sample becomes the mixture the full rate is handed to.
    const double globalScale = totWHSTP > 0.0 ? hstpRef::crossingRateHz / totWHSTP : 0.0;
    printf("  === Global post-HSTP scale = f_BX / sum(w|HSTP) = %.4g ===\n", globalScale);
    printf("\n  JZ   true xsec share   share after global scale   ratio (after / true)\n");
    printf("  --------------------------------------------------------------------------\n");
    for (unsigned int jz = 0; jz < nSlices; ++jz) {
        if (!haveSlice[jz]) continue;
        const double trueShare = 100.0 * hstpRef::sliceShare(jz) / hstpRef::targetRate;
        const double globShare = 100.0 * sumWHSTPs[jz] * globalScale / hstpRef::crossingRateHz;
        printf("  %2u %16.3f%% %25.3f%% %22.4g\n",
               jz, trueShare, globShare, trueShare > 0.0 ? globShare / trueShare : 0.0);
    }
    printf("  --------------------------------------------------------------------------\n");
    printf("\n  Every threshold moves by the same %.4g, including the high-threshold region fed by\n", globalScale);
    printf("  JZ2 and JZ3, which keep %.3f and %.3f of their weight through the filter and so were\n",
           nSlices > 2 ? keepFractions[2] : 0.0, nSlices > 3 ? keepFractions[3] : 0.0);
    printf("  already normalized correctly. Check JZSlices/*_JZSlices_Rate.pdf for which slices\n");
    printf("  populate your working point before reading the scaled rate off the plot.\n\n");
}
