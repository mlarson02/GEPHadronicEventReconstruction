#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <array>
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include <cmath>
#include <TMath.h>
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TROOT.h"
#include <cstdio>
#include <filesystem>
#include "emulationHelperFunctions_MET.h"


// Function for processing provided number of events with MET algorithm,
// then writing output MET data to a new TTree
void eventLoop(std::string inputNTuplePath, std::string outputNTuplePath,
               bool useSKObjects) {

    // GEPCellsTowers
    std::vector<double>* gepCellsTowersEtValues  = nullptr;
    std::vector<double>* gepCellsTowersEtaValues = nullptr;
    std::vector<double>* gepCellsTowersPhiValues = nullptr;

    // GEPCellsTowers (PU-suppressed)
    std::vector<double>* gepCellsTowersSKEtValues  = nullptr;
    std::vector<double>* gepCellsTowersSKEtaValues = nullptr;
    std::vector<double>* gepCellsTowersSKPhiValues = nullptr;

    // WTA-cone jets built from towers
    std::vector<double>*       gepWTAConeCellsTowersJetspTValues           = nullptr;
    std::vector<double>*       gepWTAConeCellsTowersJetsEtaValues          = nullptr;
    std::vector<double>*       gepWTAConeCellsTowersJetsPhiValues          = nullptr;
    std::vector<unsigned int>* gepWTAConeCellsTowersJetsNConstituentsValues = nullptr;

    // WTA-cone jets built from towers (PU-suppressed)
    std::vector<double>*       gepWTAConeCellsTowersSKJetspTValues           = nullptr;
    std::vector<double>*       gepWTAConeCellsTowersSKJetsEtaValues          = nullptr;
    std::vector<double>*       gepWTAConeCellsTowersSKJetsPhiValues          = nullptr;
    std::vector<unsigned int>* gepWTAConeCellsTowersSKJetsNConstituentsValues = nullptr;

    // Open input ROOT file
    TFile* inputFile = TFile::Open(inputNTuplePath.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputNTuplePath << std::endl;
        return;
    }

    // Open output ROOT file - assumed already copied from input, write new tree into it
    TFile* outputFile = TFile::Open(outputNTuplePath.c_str(), "UPDATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << outputNTuplePath << std::endl;
        return;
    }

    // Input TTrees
    TTree* gepCellsTowersTree              = (TTree*)inputFile->Get("gepCellsTowersTree");
    TTree* gepCellsTowersSKTree            = (TTree*)inputFile->Get("gepCellsTowersSKTree");
    TTree* gepWTAConeCellsTowersJetsTree   = (TTree*)inputFile->Get("gepWTAConeCellsTowersJetsTree");
    TTree* gepWTAConeCellsTowersSKJetsTree = (TTree*)inputFile->Get("gepWTAConeCellsTowersSKJetsTree");

    // Output MET values (one scalar set per event)
    double out_jetMetX    = 0.0;
    double out_jetMetY    = 0.0;
    double out_jetMet     = 0.0;
    double out_SumJetET   = 0.0;
    double out_towerMetX  = 0.0;
    double out_towerMetY  = 0.0;
    double out_towerMet   = 0.0;
    double out_SumTowerET = 0.0;
    double out_SumET      = 0.0;
    double out_totalMET   = 0.0;
    double out_totalMETX  = 0.0;
    double out_totalMETY  = 0.0;

    outputFile->cd();
    TTree* metTree = new TTree("metTree", "Tree storing event-wise MET data");
    metTree->Branch("JetMetX",    &out_jetMetX);
    metTree->Branch("JetMetY",    &out_jetMetY);
    metTree->Branch("JetMet",     &out_jetMet);
    metTree->Branch("SumJetET",   &out_SumJetET);
    metTree->Branch("TowerMetX",  &out_towerMetX);
    metTree->Branch("TowerMetY",  &out_towerMetY);
    metTree->Branch("TowerMet",   &out_towerMet);
    metTree->Branch("SumTowerET", &out_SumTowerET);
    metTree->Branch("TotalMETX",  &out_totalMETX);
    metTree->Branch("TotalMETY",  &out_totalMETY);
    metTree->Branch("SumET",      &out_SumET);
    metTree->Branch("TotalMET",   &out_totalMET);

    // === gepCellsTowersTree ===
    gepCellsTowersTree->SetBranchAddress("Et",  &gepCellsTowersEtValues);
    gepCellsTowersTree->SetBranchAddress("Eta", &gepCellsTowersEtaValues);
    gepCellsTowersTree->SetBranchAddress("Phi", &gepCellsTowersPhiValues);

    // === gepCellsTowersSKTree (PU-suppressed) ===
    gepCellsTowersSKTree->SetBranchAddress("Et",  &gepCellsTowersSKEtValues);
    gepCellsTowersSKTree->SetBranchAddress("Eta", &gepCellsTowersSKEtaValues);
    gepCellsTowersSKTree->SetBranchAddress("Phi", &gepCellsTowersSKPhiValues);

    // === gepWTAConeCellsTowersJetsTree ===
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("pT",           &gepWTAConeCellsTowersJetspTValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Eta",          &gepWTAConeCellsTowersJetsEtaValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Phi",          &gepWTAConeCellsTowersJetsPhiValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("NConstituents",&gepWTAConeCellsTowersJetsNConstituentsValues);

    // === gepWTAConeCellsTowersSKJetsTree (PU-suppressed) ===
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("pT",           &gepWTAConeCellsTowersSKJetspTValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Eta",          &gepWTAConeCellsTowersSKJetsEtaValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Phi",          &gepWTAConeCellsTowersSKJetsPhiValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("NConstituents",&gepWTAConeCellsTowersSKJetsNConstituentsValues);

    // Verifying constants
    std::cout << "half_pi_digitized_in_phi_: " << half_pi_digitized_in_phi_ << "\n";
    std::cout << "pi_digitized_in_phi_: " << pi_digitized_in_phi_ << "\n";

    // Event loop
    unsigned int eventsToProcess = gepCellsTowersTree->GetEntries();
    for (unsigned int iEvt = 0; iEvt < eventsToProcess; iEvt++) {

        // Reset output values for this event
        out_jetMetX    = 0.0;
        out_jetMetY    = 0.0;
        out_jetMet     = 0.0;
        out_SumJetET   = 0.0;
        out_towerMetX  = 0.0;
        out_towerMetY  = 0.0;
        out_towerMet   = 0.0;
        out_SumTowerET = 0.0;
        out_totalMETX  = 0.0;
        out_totalMETY  = 0.0;
        out_SumET      = 0.0;
        out_totalMET   = 0.0;

        // Load relevant input trees for this event
        if (useSKObjects) {
            gepCellsTowersSKTree->GetEntry(iEvt);
            gepWTAConeCellsTowersSKJetsTree->GetEntry(iEvt);
        } else {
            gepCellsTowersTree->GetEntry(iEvt);
            gepWTAConeCellsTowersJetsTree->GetEntry(iEvt);
        }

        // Select tower collection based on useSKObjects
        const std::vector<double>* towerEtVec  = useSKObjects ? gepCellsTowersSKEtValues  : gepCellsTowersEtValues;
        const std::vector<double>* towerPhiVec = useSKObjects ? gepCellsTowersSKPhiValues : gepCellsTowersPhiValues;

        // Clamp to maxTowersConsidered_ (runtime parameter) then to available count
        unsigned int towersProcessed = maxTowersConsidered_;
        if (towersProcessed > towerEtVec->size()) towersProcessed = towerEtVec->size();

        // -------------------------------------------------------
        // ---                   MET algorithm:                ---
        // -------------------------------------------------------

        // Tower loop (for TowerMetX/TowerMetY and SumET)
        int towerETxSum = 0;
        int towerETySum = 0;
        unsigned int towerSumEt = 0; 
        for (unsigned int iTower = 0; iTower < towersProcessed; iTower++) {
            std::cout << "towerET (GeV): " << towerEtVec->at(iTower) << "\n";
            std::cout << "towerPhi (undigi): " << towerPhiVec->at(iTower) << "\n";
            unsigned int towerEt  = digitize(towerEtVec->at(iTower),  et_bit_length_,  static_cast<double>(et_min_),  static_cast<double>(et_max_));
            unsigned int towerPhi = digitize(towerPhiVec->at(iTower), phi_bit_length_, phi_min_, phi_max_);
            std::cout << "towerET: " << towerEt << "\n";
            std::cout << "towerPhi: " << towerPhi << "\n";
            towerSumEt += towerEt;

            int towerCosPhi = sinLUT_[towerPhi + half_pi_digitized_in_phi_];
            int towerSinPhi = sinLUT_[towerPhi];

            std::cout << "towerCosPhi: " << towerCosPhi << "\n";
            std::cout << "towerSinPhi: " << towerSinPhi << "\n";

            int towerETx = (static_cast<int>(towerEt) * towerCosPhi) / (1 << (sin_bit_length_ - 1));
            int towerETy = (static_cast<int>(towerEt) * towerSinPhi) / (1 << (sin_bit_length_ - 1));
            std::cout << "towerETx: " << towerETx << "\n";
            std::cout << "towerETy: " << towerETy << "\n";
            towerETxSum += towerETx;
            towerETySum += towerETy;
        }
        int towerMETx = -towerETxSum; // Take negative of ET sums 
        int towerMETy = -towerETySum; // Take negative of ET sums 

        // Select jet collection based on useSKObjects
        const std::vector<double>* jetPtVec  = useSKObjects ? gepWTAConeCellsTowersSKJetspTValues  : gepWTAConeCellsTowersJetspTValues;
        const std::vector<double>* jetPhiVec = useSKObjects ? gepWTAConeCellsTowersSKJetsPhiValues : gepWTAConeCellsTowersJetsPhiValues;

        // Clamp to available count
        unsigned int jetsProcessed = maxJetsConsidered_;
        if (jetsProcessed > jetPtVec->size()) jetsProcessed = jetPtVec->size();

        // Jet loop (for JetMetX/JetMetY)
        int jetETxSum = 0;
        int jetETySum = 0;
        unsigned int jetSumEt = 0; 
        for (unsigned int iJet = 0; iJet < jetsProcessed; iJet++) {
            unsigned int jetEt  = digitize(jetPtVec->at(iJet),  et_bit_length_,  static_cast<double>(et_min_),  static_cast<double>(et_max_));
            unsigned int jetPhi = digitize(jetPhiVec->at(iJet), phi_bit_length_, phi_min_, phi_max_);

            jetSumEt += jetEt;

            int jetCosPhi = sinLUT_[jetPhi + half_pi_digitized_in_phi_];
            int jetSinPhi = sinLUT_[jetPhi];

            int jetETx = (static_cast<int>(jetEt) * jetCosPhi) / (1 << (sin_bit_length_ - 1));
            int jetETy = (static_cast<int>(jetEt) * jetSinPhi) / (1 << (sin_bit_length_ - 1));
            jetETxSum += jetETx;
            jetETySum += jetETy;
        }
        int jetMETx = -jetETxSum; // Take negative of ET sums
        int jetMETy = -jetETySum; // Take negative of ET sums

        int          totalMETx  = towerMETx + jetMETx;
        int          totalMETy  = towerMETy + jetMETy;
        unsigned int towerMET   = static_cast<unsigned int>(std::sqrt(towerMETx * towerMETx + towerMETy * towerMETy));
        unsigned int jetMET     = static_cast<unsigned int>(std::sqrt(jetMETx   * jetMETx   + jetMETy   * jetMETy));
        unsigned int totalMET   = static_cast<unsigned int>(std::sqrt(totalMETx * totalMETx + totalMETy * totalMETy));
        unsigned int totalSumET = towerSumEt + jetSumEt;

        // --- Pack into bitsets (bitwise-accurate representation) ---
        // Signed components (sign-magnitude: MSB = sign bit, remaining bits = magnitude)
        auto pack_signed_et = [](int val) -> uint64_t {
            uint64_t sign = (val < 0) ? 1u : 0u;
            uint64_t mag  = static_cast<uint64_t>(std::abs(val)) & maskN(signed_et_bit_length_ - 1);
            return (sign << (signed_et_bit_length_ - 1)) | mag;
        };
        std::bitset<signed_et_bit_length_> towerMETx_bitset(pack_signed_et(towerMETx));
        std::bitset<signed_et_bit_length_> towerMETy_bitset(pack_signed_et(towerMETy));
        std::bitset<signed_et_bit_length_> jetMETx_bitset  (pack_signed_et(jetMETx));
        std::bitset<signed_et_bit_length_> jetMETy_bitset  (pack_signed_et(jetMETy));
        std::bitset<signed_et_bit_length_> totalMETx_bitset(pack_signed_et(totalMETx));
        std::bitset<signed_et_bit_length_> totalMETy_bitset(pack_signed_et(totalMETy));
        // Unsigned magnitudes and sums
        std::bitset<et_bit_length_> towerMET_bitset  (towerMET   & maskN(et_bit_length_));
        std::bitset<et_bit_length_> jetMET_bitset    (jetMET     & maskN(et_bit_length_));
        std::bitset<et_bit_length_> totalMET_bitset  (totalMET   & maskN(et_bit_length_));
        std::bitset<et_bit_length_> towerSumEt_bitset(towerSumEt & maskN(et_bit_length_));
        std::bitset<et_bit_length_> jetSumEt_bitset  (jetSumEt   & maskN(et_bit_length_));
        std::bitset<et_bit_length_> totalSumET_bitset(totalSumET & maskN(et_bit_length_));

        // --- Undigitize to doubles for ntuple output ---
        out_towerMetX  = undigitize_signed_et(towerMETx_bitset);
        out_towerMetY  = undigitize_signed_et(towerMETy_bitset);
        out_jetMetX    = undigitize_signed_et(jetMETx_bitset);
        out_jetMetY    = undigitize_signed_et(jetMETy_bitset);
        out_totalMETX  = undigitize_signed_et(totalMETx_bitset);
        out_totalMETY  = undigitize_signed_et(totalMETy_bitset);
        out_towerMet   = undigitize_et(towerMET_bitset);
        out_jetMet     = undigitize_et(jetMET_bitset);
        out_totalMET   = undigitize_et(totalMET_bitset);
        out_SumTowerET = undigitize_et(towerSumEt_bitset);
        out_SumJetET   = undigitize_et(jetSumEt_bitset);
        out_SumET      = undigitize_et(totalSumET_bitset);

        metTree->Fill();
    } // Event loop

    outputFile->cd();
    std::cout << "Writing output file\n";
    metTree->Write("", TObject::kOverwrite);
    outputFile->Close();
    inputFile->Close();

} // eventLoop

// Main MET emulation function
// Use: root -b -l -q -e '.L metEmulation.cc+; metEmulation(true, false, "ZvvHbb")'
void metEmulation(bool signalBool,                 // true = signal sample, false = background (dijet)
                  bool useSKObjects,               // true = use PU-suppressed (SoftKiller) objects
                  std::string signalString = "ZvvHbb"         // Which signal sample: "VBF_hh_bbbb_cvv0/1", "ggF_hh_bbbb",
                                                              //   "ZvvHbb", "ttbar_had", "Zprime_ttbar", "ttbar_semilep", "ttbar_lep"
                  ) {

    if (signalBool) std::cout << "Processing signal: " << signalString << "\n";

    auto infile  = makeInputFileName(signalBool, signalString);
    auto outfile = makeOutputMETFileName(maxTowersConsidered_, signalBool, signalString, useSKObjects);

    std::cout << "infile:  " << infile  << "\n";
    std::cout << "outfile: " << outfile << "\n";

    try {
        std::filesystem::copy_file(infile, outfile,
                                   std::filesystem::copy_options::overwrite_existing);
        std::cout << "File copied successfully\n";
    } catch (std::filesystem::filesystem_error& e) {
        std::cerr << "Copy failed: " << e.what() << '\n';
    }
    gSystem->RedirectOutput("debuglog_MET.log", "w");
    std::cout << "Calling event loop\n";
    eventLoop(infile, outfile, useSKObjects);
}
