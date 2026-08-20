#ifndef CHAIN_SOURCE_H
#define CHAIN_SOURCE_H

// Drop-in stand-in for a TFile* that transparently spans several input ntuples.
//
// Motivation: the ten QCD JZ slices used to be hadd'd into one file, but the
// combined trees exceed ROOT's default 100 GB TTree::fgMaxTreeSize, at which
// point hadd calls TTree::ChangeFile() mid-merge and aborts, leaving a truncated
// ntuple whose later trees are silently absent. Reading the slices directly as
// TChains removes both the merge step and the 100 GB duplicate on disk.
//
// The analysis macros only ever use their TFile* handles for IsZombie(),
// Get("treeName") and Close(), so this exposes exactly those three and every
// existing `(TTree*)f->Get("...")` call site is left untouched.
//
//   ChainSource* f = ChainSource::Open("/path/QCD_Dijet_JZ*_v4/...JZ*...root");
//   TTree* t = (TTree*)f->Get("eventInfoTree");   // nullptr if the tree is absent
//   f->Close();
//
// The argument is a single path or a glob. A glob keeps the callers' existing
// std::vector<std::pair<std::string,std::string>> input config intact — only the
// path literal changes.

#include <glob.h>

#include <algorithm>
#include <cctype>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TTree.h"

// Progress/timing printouts. Get() below opens EVERY file in the chain to count entries, which
// on a busy /data can take minutes per tree with no output at all — set this to see which file
// it is on and how long each tree took.
static const bool kChainSourceVerbose = true;

// --- Read-time branch pruning -----------------------------------------------------------------
// GetEntry() decompresses every ENABLED branch of a tree, whether or not the caller ever looks at
// it. These ntuples carry far more per tree than any one analysis reads — the tower and cell
// collections above all — so a large part of the read time goes on data that is decompressed and
// immediately discarded.
//
// The wanted set is deliberately NOT hand-listed: a hand-list rots silently the first time
// somebody adds a SetBranchAddress. It is read back off the tree instead, because
// SetBranchAddress is exactly what leaves a non-null address on a branch, so a branch without one
// is by construction never read by the caller.
//
// MUST be called after the caller's last SetBranchAddress for the tree and before its first
// GetEntry. A SetBranchAddress issued afterwards would point at a branch that is no longer being
// read, and the variable would silently keep whatever it held before — which is why the listing
// below prints the disabled names: anything the macro needs turning up under "disabled" is the
// symptom to look for.
//
// Returns the number of branches disabled.
inline int DisableUnusedBranches(TTree* t, const std::string& label = "", bool verbose = true) {
    if (!t) return 0;
    TObjArray* branches = t->GetListOfBranches();
    if (!branches) return 0;

    std::vector<std::string> used, unused;
    for (int i = 0; i < branches->GetEntriesFast(); ++i) {
        TBranch* b = (TBranch*)branches->At(i);
        if (!b) continue;
        if (b->GetAddress()) used.push_back(b->GetName());
        else                 unused.push_back(b->GetName());
    }

    const std::string name = label.empty() ? t->GetName() : label;

    // A tree with no addressed branch at all is one the caller only opened, or only calls
    // GetEntry on to advance. Disabling everything there would make GetEntry read zero bytes,
    // which is probably what is wanted but is too surprising to do silently, so leave it be.
    if (used.empty()) {
        if (verbose && !unused.empty())
            std::cout << "  [branches] " << name << ": no branch has an address set — left alone ("
                      << unused.size() << " branches)\n";
        return 0;
    }

    for (size_t i = 0; i < unused.size(); ++i)
        t->SetBranchStatus(unused[i].c_str(), 0);

    if (verbose) {
        std::cout << "  [branches] " << name << ": " << used.size() << " used / "
                  << (used.size() + unused.size()) << " total, " << unused.size() << " disabled\n";
        auto printList = [](const char* tag, const std::vector<std::string>& v) {
            if (v.empty()) return;
            std::cout << "             " << tag;
            size_t col = 0;
            for (size_t i = 0; i < v.size(); ++i) {
                if (col > 5) { std::cout << "\n                       "; col = 0; }
                std::cout << v[i] << (i + 1 < v.size() ? ", " : "");
                ++col;
            }
            std::cout << "\n";
        };
        printList("used    : ", used);
        printList("disabled: ", unused);
        std::cout << std::flush;
    }
    return (int)unused.size();
}

class ChainSource {
public:
    // Returns nullptr only if the pattern matches nothing, mirroring TFile::Open.
    // A file that exists but cannot be read leaves the object non-null and zombie,
    // so callers' `!f || f->IsZombie()` checks keep working unchanged.
    static ChainSource* Open(const char* pathOrGlob) {
        std::vector<std::string> files = ExpandGlob(pathOrGlob);
        if (files.empty()) {
            std::cerr << "ChainSource: no files match " << pathOrGlob << "\n";
            return nullptr;
        }

        ChainSource* src = new ChainSource();
        src->fFiles = files;

        // Probe the first file once for its tree names, so Get() can reproduce
        // TFile::Get's nullptr-for-absent-tree behaviour without reopening. The
        // macros rely on that nullptr to skip whole plot blocks (e.g. the
        // hasGFexSimMET guard), so returning an empty chain instead would
        // silently fill those plots with nothing.
        TFile* probe = TFile::Open(files.front().c_str(), "READ");
        if (!probe || probe->IsZombie()) {
            std::cerr << "ChainSource: cannot open " << files.front() << "\n";
            src->fZombie = true;
            delete probe;
            return src;
        }
        TIter next(probe->GetListOfKeys());
        while (TKey* key = (TKey*)next()) {
            if (TClass::GetClass(key->GetClassName()) &&
                TClass::GetClass(key->GetClassName())->InheritsFrom("TTree"))
                src->fTreeNames.insert(key->GetName());
        }
        probe->Close();
        delete probe;

        if (files.size() > 1)
            std::cout << "ChainSource: chaining " << files.size() << " file(s) from "
                      << pathOrGlob << "\n";
        return src;
    }

    bool IsZombie() const { return fZombie; }

    TTree* Get(const char* treeName) {
        const std::string name(treeName);

        std::map<std::string, TChain*>::iterator it = fChains.find(name);
        if (it != fChains.end()) return it->second;

        if (fTreeNames.find(name) == fTreeNames.end()) return nullptr;  // as TFile::Get

        TChain* chain = new TChain(treeName);
        for (size_t i = 0; i < fFiles.size(); ++i) chain->Add(fFiles[i].c_str());

        TStopwatch swGet;
        if (kChainSourceVerbose) {
            swGet.Start();
            std::cout << "ChainSource: Get(" << name << ") — opening " << fFiles.size()
                      << " file(s) to count entries...\n" << std::flush;
        }

        // TChain::FindBranch returns null until a tree is loaded, and callers use
        // it to detect optional branches (e.g. metTruthNonIntX). Load the first
        // tree up front so those probes behave as they did on a plain TFile.
        chain->LoadTree(0);
        if (kChainSourceVerbose) {
            // RealTime() stops the watch, so restart it for the GetEntries leg.
            std::cout << "ChainSource:   LoadTree(0) done at " << swGet.RealTime() << " s"
                      << " (first file: " << fFiles.front() << ")\n" << std::flush;
            swGet.Continue();
        }

        // All trees here are per-event, so a mismatch means the slices were not
        // written in lockstep and the shared iEvt index across trees is invalid.
        const Long64_t n = chain->GetEntries();
        if (kChainSourceVerbose)
            std::cout << "ChainSource:   " << name << " -> " << n << " entries, "
                      << swGet.RealTime() << " s total\n" << std::flush;
        if (fEntries < 0) {
            fEntries = n;
            fEntriesFrom = name;
        } else if (n != fEntries) {
            std::cerr << "ChainSource: WARNING entry-count mismatch — " << name
                      << " has " << n << " but " << fEntriesFrom << " has " << fEntries
                      << "; per-event indexing across trees will be wrong.\n";
        }

        fChains[name] = chain;
        return chain;
    }

    // Prune every tree this source has handed out — one call per input file, placed after the
    // macro's SetBranchAddress block. Only trees actually fetched via Get() exist here, so a tree
    // the macro never asked for costs nothing either way.
    int DisableUnusedBranches(bool verbose = true) {
        int nDisabled = 0;
        if (verbose && !fChains.empty())
            std::cout << "ChainSource: pruning unused branches across " << fChains.size()
                      << " tree(s)\n" << std::flush;
        for (std::map<std::string, TChain*>::iterator it = fChains.begin();
             it != fChains.end(); ++it)
            nDisabled += ::DisableUnusedBranches(it->second, it->first, verbose);
        if (verbose)
            std::cout << "ChainSource: " << nDisabled << " branch(es) disabled in total\n"
                      << std::flush;
        return nDisabled;
    }

    void Close() {
        for (std::map<std::string, TChain*>::iterator it = fChains.begin();
             it != fChains.end(); ++it)
            delete it->second;
        fChains.clear();
    }

    const std::vector<std::string>& Files() const { return fFiles; }

private:
    ChainSource() : fZombie(false), fEntries(-1) {}

    // Sorts digit runs numerically so JZ10 would follow JZ9 rather than JZ1.
    // Chain order sets the global entry index, and every tree must agree on it.
    static bool NaturalLess(const std::string& a, const std::string& b) {
        size_t i = 0, j = 0;
        while (i < a.size() && j < b.size()) {
            if (std::isdigit((unsigned char)a[i]) && std::isdigit((unsigned char)b[j])) {
                size_t i0 = i, j0 = j;
                while (i < a.size() && std::isdigit((unsigned char)a[i])) ++i;
                while (j < b.size() && std::isdigit((unsigned char)b[j])) ++j;
                const std::string na = a.substr(i0, i - i0), nb = b.substr(j0, j - j0);
                const unsigned long va = std::stoul(na), vb = std::stoul(nb);
                if (va != vb) return va < vb;
            } else {
                if (a[i] != b[j]) return a[i] < b[j];
                ++i; ++j;
            }
        }
        return a.size() < b.size();
    }

    static std::vector<std::string> ExpandGlob(const char* pattern) {
        std::vector<std::string> out;
        glob_t g;
        const int rc = glob(pattern, GLOB_NOSORT, nullptr, &g);
        if (rc == 0) {
            for (size_t i = 0; i < g.gl_pathc; ++i) out.push_back(g.gl_pathv[i]);
        }
        globfree(&g);
        std::sort(out.begin(), out.end(), NaturalLess);
        return out;
    }

    std::vector<std::string>      fFiles;
    std::set<std::string>         fTreeNames;
    std::map<std::string, TChain*> fChains;
    bool                          fZombie;
    Long64_t                      fEntries;
    std::string                   fEntriesFrom;
};

#endif  // CHAIN_SOURCE_H