//
// compareJFexMET.C
//
// Quick overlay of the jFEX MET branch (`jFex_met`, MeV) from two TrigGepPerf
// validation ntuples, matched event-by-event on (runNumber, eventNumber).
//
//   TEST : v26 ..._jFEXDBOverride_PU200_TEST_EXT0/validation.root
//   REF  : v24 ..._jFEXDBOverride_PU200_EXT0/validation_original.root
//
// Produces
//   <outTag>.pdf   overlay of the two MET spectra + TEST/REF ratio panel
//   <outTag>.root  TTree "jFexMETCompare" with both branches (GeV) per matched
//                  event, plus the two TH1Fs
//
// Usage:
//   root -l -b -q 'compareJFexMET.C'
//   root -l -b -q 'compareJFexMET.C("test.root","ref.root","jFexMET_v26_vs_v24",100,0,500,-1)'
//
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>

#include <cstdio>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

constexpr double MEV2GEV = 1.0e-3;

// Petroff kP6: blue / orange
const char* kTestHex = "#5790fc";
const char* kRefHex  = "#f89c20";

ULong64_t makeKey(int run, int event)
{
    return (static_cast<ULong64_t>(static_cast<UInt_t>(run)) << 32) |
            static_cast<ULong64_t>(static_cast<UInt_t>(event));
}

// Open <path>, grab the TrigGepPerf "ntuple" tree and switch off everything
// except the three branches we need.
TTree* openNtuple(const std::string& path, TFile*& file)
{
    file = TFile::Open(path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Could not open: " << path << "\n";
        return nullptr;
    }

    TTree* t = nullptr;
    file->GetObject("ntuple", t);
    if (!t) {
        std::cerr << "No TTree named 'ntuple' in: " << path << "\n";
        return nullptr;
    }

    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("runNumber",   1);
    t->SetBranchStatus("eventNumber", 1);
    t->SetBranchStatus("jFex_met",    1);
    return t;
}

}  // namespace

void compareJFexMET(
    const char* testPath = "/data/larsonma/jFEXValidation/"
        "user.mlarson.GEPNtupleJETM42.QCD_Dijet_JZ4.PU200.v26_GEPHadronicReconstruction_jFEXDBOverride_PU200_TEST_EXT0/"
        "validation.root",
    const char* refPath  = "/data/larsonma/GEPHadronicEventReconstruction/GEPOutputReaderNTuples/QCD_Dijet/JZ4/"
        "user.mlarson.GEPNtupleJETM42.QCD_Dijet_JZ4.PU200.v24_GEPHadronicReconstruction_jFEXDBOverride_PU200_EXT0/"
        "validation_original.root",
    const char* outTag   = "jFexMET_v26_vs_v24",
    int    nBins   = 100,
    double metMin  = 0.0,
    double metMax  = 500.0,
    Long64_t maxEntries = -1)   // -1 = all events
{
    gStyle->SetOptStat(0);

    // ---------------------------------------------------------------- inputs
    TFile* fTest = nullptr;
    TFile* fRef  = nullptr;
    TTree* tTest = openNtuple(testPath, fTest);
    TTree* tRef  = openNtuple(refPath,  fRef);
    if (!tTest || !tRef) return;

    int   runTest = 0, evtTest = 0;
    int   runRef  = 0, evtRef  = 0;
    float metTest = 0.0f, metRef = 0.0f;

    tTest->SetBranchAddress("runNumber",   &runTest);
    tTest->SetBranchAddress("eventNumber", &evtTest);
    tTest->SetBranchAddress("jFex_met",    &metTest);

    tRef->SetBranchAddress("runNumber",   &runRef);
    tRef->SetBranchAddress("eventNumber", &evtRef);
    tRef->SetBranchAddress("jFex_met",    &metRef);

    // ----------------------------------------- pass 1: index the reference file
    const Long64_t nRef = (maxEntries > 0 && maxEntries < tRef->GetEntries())
                            ? maxEntries : tRef->GetEntries();
    std::unordered_map<ULong64_t, float> refMET;
    refMET.reserve(static_cast<size_t>(nRef) * 2);

    std::cout << "Indexing REF  (" << nRef << " events): " << refPath << "\n";
    for (Long64_t i = 0; i < nRef; ++i) {
        tRef->GetEntry(i);
        refMET[makeKey(runRef, evtRef)] = metRef;
        if (i % 100000 == 0 && i > 0) std::cout << "  ref entry " << i << "\n";
    }

    // ------------------------------------------------------------ output setup
    TFile* fOut = TFile::Open(Form("%s.root", outTag), "RECREATE");

    int    outRun = 0, outEvt = 0;
    double outTestMET = 0.0, outRefMET = 0.0, outDiff = 0.0;

    TTree* outTree = new TTree("jFexMETCompare", "jFEX MET, TEST (v26) vs REF (v24)");
    outTree->Branch("runNumber",    &outRun);
    outTree->Branch("eventNumber",  &outEvt);
    outTree->Branch("jFexMET_test", &outTestMET);   // GeV
    outTree->Branch("jFexMET_ref",  &outRefMET);    // GeV
    outTree->Branch("jFexMET_diff", &outDiff);      // GeV, test - ref

    TH1F* hTest = new TH1F("hJFexMET_test", ";jFEX E_{T}^{miss} [GeV];Events",
                           nBins, metMin, metMax);
    TH1F* hRef  = new TH1F("hJFexMET_ref",  ";jFEX E_{T}^{miss} [GeV];Events",
                           nBins, metMin, metMax);
    hTest->Sumw2();
    hRef->Sumw2();

    // ------------------------------------------- pass 2: loop the test file
    const Long64_t nTest = (maxEntries > 0 && maxEntries < tTest->GetEntries())
                             ? maxEntries : tTest->GetEntries();
    Long64_t nMatched = 0, nUnmatched = 0, nDiffer = 0;

    std::cout << "Looping TEST (" << nTest << " events): " << testPath << "\n";
    for (Long64_t i = 0; i < nTest; ++i) {
        tTest->GetEntry(i);

        auto it = refMET.find(makeKey(runTest, evtTest));
        if (it == refMET.end()) { ++nUnmatched; continue; }
        ++nMatched;

        outRun     = runTest;
        outEvt     = evtTest;
        outTestMET = metTest      * MEV2GEV;
        outRefMET  = it->second   * MEV2GEV;
        outDiff    = outTestMET - outRefMET;
        if (outDiff != 0.0) ++nDiffer;

        outTree->Fill();
        hTest->Fill(outTestMET);
        hRef->Fill(outRefMET);

        if (i % 100000 == 0 && i > 0) std::cout << "  test entry " << i << "\n";
    }

    std::cout << "matched: " << nMatched
              << "  unmatched (test w/o ref): " << nUnmatched
              << "  differing MET: " << nDiffer << "\n";

    // ---------------------------------------------------------------- plotting
    const Color_t cTest = TColor::GetColor(kTestHex);
    const Color_t cRef  = TColor::GetColor(kRefHex);

    hTest->SetLineColor(cTest);
    hTest->SetMarkerColor(cTest);
    hTest->SetLineWidth(2);
    hRef->SetLineColor(cRef);
    hRef->SetMarkerColor(cRef);
    hRef->SetLineWidth(2);
    hRef->SetLineStyle(2);

    TCanvas* c = new TCanvas("cJFexMET", "jFEX MET comparison", 800, 800);

    TPad* pTop = new TPad("pTop", "pTop", 0.0, 0.30, 1.0, 1.0);
    TPad* pBot = new TPad("pBot", "pBot", 0.0, 0.00, 1.0, 0.30);
    pTop->SetBottomMargin(0.02);
    pTop->SetLogy();
    pBot->SetTopMargin(0.04);
    pBot->SetBottomMargin(0.32);
    pTop->Draw();
    pBot->Draw();

    pTop->cd();
    hTest->GetXaxis()->SetLabelSize(0);
    hTest->GetYaxis()->SetTitleSize(0.05);
    hTest->GetYaxis()->SetTitleOffset(1.0);
    hTest->SetMaximum(hTest->GetMaximum() * 10.0);
    hTest->Draw("HIST E");
    hRef->Draw("HIST E SAME");

    TLegend* leg = new TLegend(0.52, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hTest, Form("v26 TEST (#mu = %.1f GeV)", hTest->GetMean()), "l");
    leg->AddEntry(hRef,  Form("v24 original (#mu = %.1f GeV)", hRef->GetMean()), "l");
    leg->Draw();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.DrawLatex(0.16, 0.86, "QCD Dijet JZ4, PU200");
    lat.DrawLatex(0.16, 0.81, Form("%lld matched events", nMatched));

    pBot->cd();
    TH1F* hRatio = static_cast<TH1F*>(hTest->Clone("hJFexMET_ratio"));
    hRatio->Divide(hRef);
    hRatio->SetTitle("");
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.7);
    hRatio->GetYaxis()->SetTitle("TEST / ref");
    hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.11);
    hRatio->GetYaxis()->SetTitleOffset(0.45);
    hRatio->GetYaxis()->SetLabelSize(0.09);
    hRatio->GetXaxis()->SetTitleSize(0.13);
    hRatio->GetXaxis()->SetTitleOffset(1.05);
    hRatio->GetXaxis()->SetLabelSize(0.11);
    hRatio->Draw("E1");

    TLine* one = new TLine(metMin, 1.0, metMax, 1.0);
    one->SetLineColor(kGray + 2);
    one->SetLineStyle(2);
    one->Draw();

    c->SaveAs(Form("%s.pdf", outTag));

    // ------------------------------------------------------------------ write
    fOut->cd();
    outTree->Write("", TObject::kOverwrite);
    hTest->Write();
    hRef->Write();
    hRatio->Write();
    c->Write("cJFexMET");
    fOut->Close();

    fTest->Close();
    fRef->Close();

    std::cout << "wrote " << outTag << ".pdf and " << outTag << ".root\n";
}
