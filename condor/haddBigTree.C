// Merge ROOT files whose combined TTrees exceed ROOT's default 100 GB tree-size
// limit (TTree::fgMaxTreeSize = 1e11 bytes).
//
// Plain `hadd` cannot do this: once the output tree crosses the limit, ROOT calls
// TTree::ChangeFile(), silently swaps the output TFile for {name}_1.root mid-merge
// and the TFileMerger aborts (core dump). What is left behind is a truncated file
// containing only the trees written before the abort — later trees are simply
// absent, which downstream shows up as a null TTree* from TFile::Get().
//
// Usage:
//   root -b -l -q 'haddBigTree.C("/path/out.root", "/path/inputs.txt")'
//
// inputs.txt holds one input file path per line (blank lines and #-comments
// ignored). A list file is used rather than a variadic argument list so that the
// ten long JZ-slice paths do not have to survive shell quoting.
//
// Single-threaded by design: `hadd -j` forks subprocesses, and the raised
// fgMaxTreeSize would not propagate into them.

#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TFileMerger.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

void haddBigTree(const char* outFile,
                 const char* inputListFile,
                 Long64_t    maxTreeSizeGB = 2000) {
    // Must be set before the merger creates the output trees: TTree's constructor
    // copies fgMaxTreeSize into the per-tree fMaxTreeSize.
    TTree::SetMaxTreeSize(maxTreeSizeGB * 1000000000LL);
    std::cout << "haddBigTree: TTree::GetMaxTreeSize() = "
              << TTree::GetMaxTreeSize() << " bytes\n";

    std::ifstream in(inputListFile);
    if (!in) {
        std::cerr << "haddBigTree: cannot open input list " << inputListFile << "\n";
        gSystem->Exit(1);
    }

    std::vector<std::string> inputs;
    std::string line;
    while (std::getline(in, line)) {
        // Trim whitespace so trailing \r or stray indentation does not break paths.
        const size_t b = line.find_first_not_of(" \t\r\n");
        if (b == std::string::npos) continue;
        const size_t e = line.find_last_not_of(" \t\r\n");
        line = line.substr(b, e - b + 1);
        if (line.empty() || line[0] == '#') continue;
        inputs.push_back(line);
    }

    if (inputs.empty()) {
        std::cerr << "haddBigTree: no input files listed in " << inputListFile << "\n";
        gSystem->Exit(1);
    }

    // Match hadd's default: inherit the compression settings of the first input
    // rather than silently recompressing the whole sample.
    Int_t comp = 101;
    if (TFile* f0 = TFile::Open(inputs.front().c_str(), "READ")) {
        if (!f0->IsZombie()) comp = f0->GetCompressionSettings();
        f0->Close();
        delete f0;
    } else {
        std::cerr << "haddBigTree: cannot open first input " << inputs.front() << "\n";
        gSystem->Exit(1);
    }

    TFileMerger merger(/*isLocal=*/kFALSE, /*histoOneGo=*/kFALSE);
    merger.SetPrintLevel(1);
    if (!merger.OutputFile(outFile, "RECREATE", comp)) {
        std::cerr << "haddBigTree: cannot create output " << outFile << "\n";
        gSystem->Exit(1);
    }

    for (const std::string& f : inputs) {
        if (!merger.AddFile(f.c_str())) {
            std::cerr << "haddBigTree: cannot add " << f << "\n";
            gSystem->Exit(1);
        }
        std::cout << "  + " << f << "\n";
    }

    std::cout << "haddBigTree: merging " << inputs.size() << " file(s) -> " << outFile << "\n";
    if (!merger.Merge()) {
        std::cerr << "haddBigTree: merge FAILED\n";
        gSystem->Exit(1);
    }

    // A surviving {out}_1.root means the size limit was still hit and the output
    // is truncated — fail loudly rather than leave a silently incomplete ntuple.
    std::string stub(outFile);
    const size_t dot = stub.rfind(".root");
    if (dot != std::string::npos) {
        stub = stub.substr(0, dot) + "_1.root";
        if (gSystem->AccessPathName(stub.c_str()) == kFALSE) {
            std::cerr << "haddBigTree: output was split — " << stub
                      << " exists, so " << outFile << " is TRUNCATED.\n"
                      << "  Raise maxTreeSizeGB (currently " << maxTreeSizeGB << ") and re-run.\n";
            gSystem->Exit(1);
        }
    }

    std::cout << "haddBigTree: OK -> " << outFile << "\n";
}
