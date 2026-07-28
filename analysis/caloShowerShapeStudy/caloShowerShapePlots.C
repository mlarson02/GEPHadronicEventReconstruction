// caloShowerShapePlots.C
// ---------------------------------------------------------------------------
// Aggregate per-layer shower-shape distributions for the displaced-jet study.
// This is the *quantitative* companion to caloShowerEventDisplays.C: for each
// jet (ΔR-matched to a displaced truth-BSM particle in signal, or the two
// leading jets in dijet background) it sums the constituent-tower Et in each
// calorimeter layer l0..l6 and forms the longitudinal energy profile. Signal
// and dijet are overlaid so the layer that discriminates a displaced/late
// shower stands out.
//
//   # both signals, once each, into plots/<signal>/ (no arguments needed):
//   root -b -l -q 'caloShowerShapePlots.C'
//   # or one sample explicitly (optionally vs a dijet file):
//   root -b -q 'caloShowerShapePlots.C("signal.root","dijet.root","plots/")'
//
// One PDF per jet collection:
//   <outDir>/caloShowerShapePlots_WTACone.pdf
//   <outDir>/caloShowerShapePlots_LRJ.pdf
//
// Distributions (signal vs dijet, area-normalized):
//   * Et fraction in each layer l0..l6
//   * longitudinal shower "depth" = Et-weighted mean layer index (0..6)
//   * EM fraction  = (l0+l1+l2+l3)/total   [displaced/hadronic -> lower]
//   * n constituent towers per jet
// ---------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"

static inline double dRp(double e1,double p1,double e2,double p2){
    double dp=std::fabs(p1-p2); while(dp>M_PI) dp=std::fabs(dp-2*M_PI);
    double de=e1-e2; return std::sqrt(de*de+dp*dp);
}

// Fill one sample's shower-shape histograms.
static long fillShowerHists(const std::string& file, bool isDijet, bool isLRJ, double Rassoc,
                            double ptMinBSM, double etMinTower, bool useTruth,
                            TH1F* hfrac[7], TH1F* hdepth, TH1F* hEMfrac, TH1F* hntow) {
    TFile* fin = TFile::Open(file.c_str());
    if (!fin || fin->IsZombie()) { std::cerr << "cannot open " << file << "\n"; return 0; }

    TTree* towTree=nullptr; fin->GetObject("gepCellsTowersEtaSKTree", towTree);
    TTree* wtaTree=nullptr; fin->GetObject("wtaConeCellsTowersEtaSKTree", wtaTree);
    TTree* lrjTree=nullptr; fin->GetObject("jetTaggerLRJEtaSKTree", lrjTree);
    TTree* bsmTree=nullptr; fin->GetObject("truthBSMTree", bsmTree);
    TTree* jetTree = isLRJ ? lrjTree : wtaTree;
    // ---- diagnostics: tree presence + entries ----
    std::cout << "[fillShowerHists] " << file << (isLRJ ? "  [LRJ]" : "  [WTACone]")
              << (useTruth ? "  useTruth" : "  leadingJets") << "\n";
    for (auto pr : {std::make_pair("gepCellsTowersEtaSKTree", towTree),
                    std::make_pair("jetTree(selected)",       jetTree),
                    std::make_pair("truthBSMTree",            bsmTree)}) {
        Long64_t ne = pr.second ? pr.second->GetEntries() : -1;
        std::cout << "    " << pr.first << " entries=" << ne
                  << (!pr.second ? "  [MISSING]" : (ne == 0 ? "  [EMPTY]" : "")) << "\n";
    }
    if (!towTree || !jetTree) { std::cerr << "missing trees in " << file << "\n"; fin->Close(); return 0; }

    std::vector<double> *tow_Eta=nullptr,*tow_Phi=nullptr,*tow_Et_l[7]={nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
    towTree->SetBranchAddress("Eta",&tow_Eta); towTree->SetBranchAddress("Phi",&tow_Phi);
    for(int l=0;l<7;++l) towTree->SetBranchAddress(Form("Et_l%d",l),&tow_Et_l[l]);

    std::vector<double> *wta_pt=nullptr,*wta_eta=nullptr,*wta_phi=nullptr;
    std::vector<float>  *lrj_Et=nullptr,*lrj_eta=nullptr,*lrj_phi=nullptr;
    if (isLRJ){ lrjTree->SetBranchAddress("Et",&lrj_Et); lrjTree->SetBranchAddress("Eta",&lrj_eta); lrjTree->SetBranchAddress("Phi",&lrj_phi); }
    else      { wtaTree->SetBranchAddress("Pt",&wta_pt); wtaTree->SetBranchAddress("Eta",&wta_eta); wtaTree->SetBranchAddress("Phi",&wta_phi); }

    std::vector<int>    *bsm_hasDecayVtx=nullptr;
    std::vector<double> *bsm_pt=nullptr,*bsm_eta=nullptr,*bsm_phi=nullptr;
    if (bsmTree && !isDijet && useTruth){
        bsmTree->SetBranchAddress("pt",&bsm_pt); bsmTree->SetBranchAddress("eta",&bsm_eta);
        bsmTree->SetBranchAddress("phi",&bsm_phi); bsmTree->SetBranchAddress("hasDecayVtx",&bsm_hasDecayVtx);
    }

    long nJets=0;
    Long64_t nEvt=towTree->GetEntries();
    for (Long64_t iEvt=0; iEvt<nEvt; ++iEvt){
        towTree->GetEntry(iEvt); jetTree->GetEntry(iEvt);
        if (bsmTree && !isDijet) bsmTree->GetEntry(iEvt);

        std::vector<std::pair<double,double>> jetEP; // eta,phi
        std::vector<double> jetPt;
        if (isLRJ){ if(lrj_Et) for(size_t j=0;j<lrj_Et->size();++j){ jetEP.push_back({(double)(*lrj_eta)[j],(double)(*lrj_phi)[j]}); jetPt.push_back((*lrj_Et)[j]); } }
        else      { if(wta_pt) for(size_t j=0;j<wta_pt->size();++j){ jetEP.push_back({(*wta_eta)[j],(*wta_phi)[j]}); jetPt.push_back((*wta_pt)[j]); } }
        if (jetEP.empty()) continue;

        std::vector<int> selJets;
        if (!isDijet && useTruth && bsm_pt){
            std::vector<int> boi;
            for(size_t b=0;b<bsm_pt->size();++b)
                if (bsm_hasDecayVtx && (*bsm_hasDecayVtx)[b]==1 && (*bsm_pt)[b]>ptMinBSM) boi.push_back((int)b);
            for(size_t j=0;j<jetEP.size();++j){
                for(int b:boi) if (dRp(jetEP[j].first,jetEP[j].second,(*bsm_eta)[b],(*bsm_phi)[b])<Rassoc){ selJets.push_back((int)j); break; }
            }
        } else {
            std::vector<int> ord(jetEP.size()); for(size_t j=0;j<jetEP.size();++j) ord[j]=(int)j;
            std::sort(ord.begin(),ord.end(),[&](int a,int b){return jetPt[a]>jetPt[b];});
            for(int k=0;k<(int)ord.size() && k<2;++k) selJets.push_back(ord[k]);
        }

        for (int sj : selJets){
            double je=jetEP[sj].first, jp=jetEP[sj].second;
            double sum[7]={0,0,0,0,0,0,0}; double total=0; int ntow=0;
            for (size_t t=0;t<tow_Eta->size();++t){
                if (dRp((*tow_Eta)[t],(*tow_Phi)[t],je,jp) >= Rassoc) continue;
                bool any=false;
                for(int l=0;l<7;++l){ double et = tow_Et_l[l]?(*tow_Et_l[l])[t]:0.0; if(et>etMinTower){ sum[l]+=et; total+=et; any=true; } }
                if (any) ++ntow;
            }
            if (total<=0) continue;
            double depth=0, emf=0;
            for(int l=0;l<7;++l){ double f=sum[l]/total; hfrac[l]->Fill(f); depth+=l*f; if(l<=3) emf+=f; }
            hdepth->Fill(depth); hEMfrac->Fill(emf); hntow->Fill(ntow);
            ++nJets;
        }
    }
    fin->Close();
    return nJets;
}

static void styleOverlay(TH1F* hs, TH1F* hd){
    if(hs){ hs->SetLineColor(kRed+1);  hs->SetLineWidth(2); if(hs->Integral()>0) hs->Scale(1.0/hs->Integral()); }
    if(hd){ hd->SetLineColor(kBlue+1); hd->SetLineWidth(2); hd->SetLineStyle(2); if(hd->Integral()>0) hd->Scale(1.0/hd->Integral()); }
}

// ===========================================================================
void caloShowerShapePlots(std::string signalFile = "",
                          std::string dijetFile = "",
                          std::string outDir    = ".",
                          double      ptMinBSM   = 5.0,
                          double      etMinTower = 0.5,
                          bool        useTruth   = false) {
    // No-arg default (root -b -l -q 'caloShowerShapePlots.C'): run both signals
    // once each, each written into its own plots/<label>/ sub-directory. Pass a
    // signalFile explicitly to process one sample (optionally vs a dijet file).
    if (signalFile.empty()) {
        const std::vector<std::pair<std::string,std::string>> samples = {
            {"/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_displaced_dark_photon.root", "displaced_dark_photon"},
            {"/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_emerging_jets.root",         "emerging_jets"},
        };
        for (const auto& s : samples) {
            std::string od = "plots/" + s.second;
            gSystem->mkdir(od.c_str(), kTRUE);
            std::cout << "[caloShowerShapePlots] === " << s.second << " -> " << od << "/ ===\n";
            caloShowerShapePlots(s.first, "", od, ptMinBSM, etMinTower);
        }
        return;
    }

    gStyle->SetOptStat(0);

    struct JetColl { const char* tag; bool isLRJ; double Rassoc; };
    std::vector<JetColl> colls = { { "WTACone", false, 0.4 }, { "LRJ", true, 1.0 } };

    for (const auto& coll : colls) {
        TH1F* sfrac[7]; TH1F* dfrac[7];
        for(int l=0;l<7;++l){
            sfrac[l]=new TH1F(Form("s_frac%d_%s",l,coll.tag),Form("layer %d E_{T} fraction;E_{T,l%d}/E_{T};jets (norm)",l,l),40,0,1);
            dfrac[l]=new TH1F(Form("d_frac%d_%s",l,coll.tag),"",40,0,1);
            sfrac[l]->SetDirectory(nullptr); dfrac[l]->SetDirectory(nullptr);
        }
        TH1F* sdepth=new TH1F(Form("s_depth_%s",coll.tag),"shower depth;E_{T}-weighted mean layer;jets (norm)",42,0,6.3);
        TH1F* ddepth=new TH1F(Form("d_depth_%s",coll.tag),"",42,0,6.3);
        TH1F* sEM   =new TH1F(Form("s_emf_%s",coll.tag),"EM fraction;(l0+l1+l2+l3)/E_{T};jets (norm)",40,0,1);
        TH1F* dEM   =new TH1F(Form("d_emf_%s",coll.tag),"",40,0,1);
        TH1F* sntow =new TH1F(Form("s_ntow_%s",coll.tag),"n constituent towers;n towers;jets (norm)",40,0,40);
        TH1F* dntow =new TH1F(Form("d_ntow_%s",coll.tag),"",40,0,40);
        for (TH1F* h : {sdepth,ddepth,sEM,dEM,sntow,dntow}) h->SetDirectory(nullptr);

        long nS=0,nD=0;
        if (!signalFile.empty()) nS=fillShowerHists(signalFile,false,coll.isLRJ,coll.Rassoc,ptMinBSM,etMinTower,useTruth,sfrac,sdepth,sEM,sntow);
        if (!dijetFile.empty())  nD=fillShowerHists(dijetFile, true, coll.isLRJ,coll.Rassoc,ptMinBSM,etMinTower,useTruth,dfrac,ddepth,dEM,dntow);
        std::cout << "[caloShowerShapePlots] " << coll.tag << ": signal jets=" << nS << "  dijet jets=" << nD << "\n";

        for(int l=0;l<7;++l) styleOverlay(sfrac[l],dfrac[l]);
        styleOverlay(sdepth,ddepth); styleOverlay(sEM,dEM); styleOverlay(sntow,dntow);

        TString pdf = TString(outDir) + "/caloShowerShapePlots_" + coll.tag + ".pdf";
        TCanvas* c = new TCanvas("cShower","shower shapes",1500,1000);
        c->Print(pdf + "[");

        auto drawPair=[&](TH1F* hs, TH1F* hd, const char* title){
            double mx=0; if(hs) mx=std::max(mx,hs->GetMaximum()); if(hd) mx=std::max(mx,hd->GetMaximum());
            if(hs){ hs->SetTitle(title); hs->SetMaximum(1.3*mx); hs->Draw("hist"); if(hd) hd->Draw("hist same"); }
            else if(hd){ hd->SetTitle(title); hd->SetMaximum(1.3*mx); hd->Draw("hist"); }
            TLegend* lg=new TLegend(0.62,0.75,0.9,0.9); lg->SetBorderSize(0);
            if(hs) lg->AddEntry(hs,Form("signal (%ld)",nS),"l");
            if(hd) lg->AddEntry(hd,Form("dijet (%ld)",nD),"l");
            lg->Draw();
        };

        // page 1: the 7 layer fractions (4x2)
        c->Clear(); c->Divide(4,2);
        for(int l=0;l<7;++l){ c->cd(l+1); drawPair(sfrac[l],dfrac[l],Form("layer %d E_{T} fraction;E_{T,l%d}/E_{T};jets (norm)",l,l)); }
        c->cd(8); { TLatex tl; tl.SetNDC(); tl.SetTextSize(0.06);
            tl.DrawLatex(0.1,0.6,Form("%s jets",coll.tag)); tl.DrawLatex(0.1,0.45,"red = signal (solid)"); tl.DrawLatex(0.1,0.33,"blue = dijet (dashed)"); }
        c->Print(pdf);

        // page 2: summary observables (2x2)
        c->Clear(); c->Divide(2,2);
        c->cd(1); drawPair(sdepth,ddepth,"shower depth;E_{T}-weighted mean layer;jets (norm)");
        c->cd(2); drawPair(sEM,dEM,"EM fraction;(l0+l1+l2+l3)/E_{T};jets (norm)");
        c->cd(3); drawPair(sntow,dntow,"n constituent towers;n towers;jets (norm)");
        c->Print(pdf);

        c->Print(pdf + "]");
        std::cout << "[caloShowerShapePlots] wrote -> " << pdf << "\n";
        delete c;
    }
}
