// checkJZSliceContents.C
// ---------------------------------------------------------------------------
// Per-JZ-slice content audit for the caloShowerShape ntuples.
//
// Motivation: the JZ1 slice contributed 174k WTACone jets but ZERO JetTaggerLRJ
// jets to caloShowerShapePlots.C, and JZ1 is the one slice reprocessed from a
// different production (v22 "PU200"; all others are v17/v18). This prints, per
// slice, how many events actually carry each object collection, so an empty
// collection is obvious rather than silently dropping a slice from the plots.
//
// It also prints the branch CLASS NAME for one branch of each collection: a
// vector<float> -> vector<double> change between productions would bind without
// the ntupler's "branch missing" warning and give empty/garbage reads.
//
//   root -b -l -q 'caloShowerShapeStudy/checkJZSliceContents.C'
//   # or audit the upstream GEPOutputReader ntuples instead (pass a GEP file):
//   root -b -q 'caloShowerShapeStudy/checkJZSliceContents.C("<gepFile>.root")'
// ---------------------------------------------------------------------------

#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

static const char* kNtupleDir = "/data/larsonma/CaloShowerShapeTriggers/ntuples";
static const Long64_t kMaxEvents = 20000;   // per slice; -1 = all

// Class name of a branch ("vector<float>", "vector<double>", ...) or "ABSENT".
static std::string branchType(TTree* t, const char* name) {
    if (!t) return "no tree";
    TBranch* b = t->GetBranch(name);
    if (!b) return "ABSENT";
    std::string cn = b->GetClassName();
    if (!cn.empty()) return cn;
    TLeaf* l = b->GetLeaf(name);
    return l ? l->GetTypeName() : "?";
}

// Fraction of events with a non-empty vector branch, plus the mean multiplicity.
template <typename T>
static void vectorOccupancy(TTree* t, const char* branch, Long64_t nMax,
                            double& fracNonEmpty, double& meanN) {
    fracNonEmpty = -1; meanN = -1;
    if (!t || !t->GetBranch(branch)) return;
    std::vector<T>* v = nullptr;
    t->SetBranchAddress(branch, &v);
    Long64_t n = t->GetEntries();
    if (nMax >= 0 && n > nMax) n = nMax;
    Long64_t nonEmpty = 0, total = 0;
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (v && !v->empty()) ++nonEmpty;
        if (v) total += (Long64_t)v->size();
    }
    fracNonEmpty = n > 0 ? (double)nonEmpty / n : 0.0;
    meanN        = n > 0 ? (double)total / n    : 0.0;
    t->ResetBranchAddresses();
}

static void auditFile(const std::string& path, const std::string& label) {
    TFile* f = TFile::Open(path.c_str());
    if (!f || f->IsZombie()) { std::cout << "  " << label << ": CANNOT OPEN\n"; return; }

    TTree* wta = nullptr; f->GetObject("wtaConeCellsTowersEtaSKTree", wta);
    TTree* lrj = nullptr; f->GetObject("jetTaggerLRJEtaSKTree",       lrj);
    TTree* tow = nullptr; f->GetObject("gepCellsTowersEtaSKTree",     tow);
    TTree* gep = nullptr; f->GetObject("ntuple",                      gep);   // upstream GEP file

    if (gep) {
        // Upstream GEPOutputReader ntuple: check the source branches directly.
        double fw=0,mw=0,fl=0,ml=0;
        vectorOccupancy<float>(gep, "WTAConeGEPCellsTowerEtaSKJets_pt",        kMaxEvents, fw, mw);
        vectorOccupancy<float>(gep, "JetTaggerLRJGEPCellsTowerEtaSKJets_Et",   kMaxEvents, fl, ml);
        std::cout << "  " << label << "  [GEP ntuple, " << gep->GetEntries() << " evt]\n"
                  << "      WTACone pt   type=" << branchType(gep,"WTAConeGEPCellsTowerEtaSKJets_pt")
                  << "  nonEmpty=" << fw << "  <n>=" << mw << "\n"
                  << "      LRJ Et       type=" << branchType(gep,"JetTaggerLRJGEPCellsTowerEtaSKJets_Et")
                  << "  nonEmpty=" << fl << "  <n>=" << ml << "\n"
                  << "      tower et_l0  type=" << branchType(gep,"GEPCellsTowerEtaSK_et_l0") << "\n";
        f->Close(); delete f; return;
    }

    double fw=0,mw=0,fl=0,ml=0,ft=0,mt=0;
    vectorOccupancy<double>(wta, "Pt", kMaxEvents, fw, mw);
    vectorOccupancy<float> (lrj, "Et", kMaxEvents, fl, ml);
    vectorOccupancy<double>(tow, "Et", kMaxEvents, ft, mt);

    Long64_t nEvt = tow ? tow->GetEntries() : (wta ? wta->GetEntries() : -1);
    std::cout << "  " << label << "  [" << nEvt << " evt]\n"
              << "      WTACone  type=" << branchType(wta,"Pt")
              << "  evtWithJet=" << fw << "  <nJet>=" << mw << "\n"
              << "      LRJ      type=" << branchType(lrj,"Et")
              << "  evtWithJet=" << fl << "  <nJet>=" << ml
              << (fl == 0.0 ? "   <== EMPTY COLLECTION" : "") << "\n"
              << "      towers   type=" << branchType(tow,"Et")
              << "  evtWithTower=" << ft << "  <nTower>=" << mt << "\n";
    f->Close(); delete f;
}

void checkJZSliceContents(std::string gepFile = "") {
    if (!gepFile.empty()) { auditFile(gepFile, "GEP: " + gepFile); return; }

    std::cout << "\n=== caloShowerShape ntuple audit (first " << kMaxEvents
              << " events per file) ===\n";
    for (int jz = 0; jz < 10; ++jz) {
        std::string p = std::string(kNtupleDir) + "/caloShowerShape_dijet_JZ"
                      + std::to_string(jz) + ".root";
        auditFile(p, "JZ" + std::to_string(jz));
    }
    for (const char* s : {"displaced_dark_photon", "emerging_jets"}) {
        std::string p = std::string(kNtupleDir) + "/caloShowerShape_" + s + ".root";
        auditFile(p, s);
    }
    std::cout << "\nIf a slice shows LRJ evtWithJet=0 while WTACone is populated, the\n"
                 "JetTaggerLRJ collection is empty upstream -- rerun this macro on that\n"
                 "slice's GEPOutputReader file to confirm the ntuple is faithful:\n"
                 "  root -b -q 'caloShowerShapeStudy/checkJZSliceContents.C(\"<gepFile>\")'\n";
}