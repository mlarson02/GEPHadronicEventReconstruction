// Sum of event weights per JZ slice for the QCD_Dijet GEPOutputReader ntuples.
//
// C++/ROOT port of getSumOfWeights.py (much faster: reads only the `weight`
// branch and sums it at compiled speed via TTree::Draw, rather than a
// per-entry PyROOT loop).
//
// Run with:
//   lsetup root            # if not already set up
//   root -l -b -q getSumOfWeights.C
//
// Output: sumOfWeights.txt  ->  "JZ<n> <sumOfWeights> <nFiles>" per line.

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

void getSumOfWeights() {
    const TString baseDir =
        "/data/larsonma/GEPHadronicEventReconstruction/GEPOutputReaderNTuples_PU140/QCD_Dijet";

    std::ofstream out("sumOfWeights_PU140.txt");

    for (int jz = 0; jz < 10; ++jz) {
        // Version tag varies per JZ slice (v17/v18/v22, JZ1 is PU200), so match
        // on the dataset dir suffix rather than a fixed version string. Let the
        // shell expand the glob and give us one path per line.
        TString pattern =
            TString::Format("%s/JZ%d/user.mlarson.GEPNtupleJETM42.QCD_Dijet_JZ%d.*_EXT0/*.root",
                            baseDir.Data(), jz, jz);

        TString lsOut = gSystem->GetFromPipe(TString::Format("ls %s 2>/dev/null", pattern.Data()));

        std::vector<TString> files;
        {
            std::stringstream ss(lsOut.Data());
            std::string line;
            while (std::getline(ss, line)) {
                if (!line.empty()) files.emplace_back(line.c_str());
            }
        }

        if (files.empty()) {
            std::cout << "JZ" << jz << ": no files matched " << pattern << std::endl;
            out << "JZ" << jz << " 0.0 0\n";
            continue;
        }

        double thisSum = 0.0;
        Long64_t nEntries = 0;

        // Loop file-by-file so the temporary weight array (TTree::Draw) only ever
        // holds one file's worth of entries -> fast and bounded memory.
        for (const TString& fname : files) {
            TFile* f = TFile::Open(fname);
            if (!f || f->IsZombie()) {
                std::cerr << "  WARNING: could not open " << fname << std::endl;
                if (f) f->Close();
                continue;
            }

            TTree* t = dynamic_cast<TTree*>(f->Get("ntuple"));
            if (!t) {
                std::cerr << "  WARNING: no 'ntuple' tree in " << fname << std::endl;
                f->Close();
                continue;
            }

            Long64_t ne = t->GetEntries();
            nEntries += ne;

            if (ne > 0) {
                t->SetEstimate(ne + 1);                       // hold all values
                Long64_t n = t->Draw("weight", "", "goff");   // eval only `weight`
                const double* v = t->GetV1();
                for (Long64_t i = 0; i < n; ++i) thisSum += v[i];
            }

            f->Close();
        }

        out << "JZ" << jz << " " << thisSum << " " << files.size() << "\n";
        out.flush();  // write incrementally so progress is visible mid-run
        std::cout << "JZ" << jz << ": " << files.size() << " files, " << nEntries
                  << " entries, Sum of Weights = " << thisSum << std::endl;
    }

    out.close();
}
