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
#include "emulationHelperFunctions_MET.h"


// Function for processing provided number of events with MET algorithm,
// then writing output MET data to a new TTree
void eventLoop(std::string inputNTuplePath, std::string outputNTuplePath,
               bool useSKObjects, double jetEtThreshold, double towerEtThreshold, bool doJetTowerOverlapRemoval,
               bool useEtaSKObjects = false,
               double towerScaleFactor = 1.0, double jetScaleFactor = 1.0) {
    // NOTE: towerScaleFactor / jetScaleFactor are scalar weights applied when
    // combining tower and jet MET into the total MET. They are intended as a
    // placeholder for what will eventually be eta-dependent, calibrated scale
    // factors derived to match truth MET. Once that calibration exists, replace
    // the scalar multiplication below with a per-object eta-binned lookup,
    // ideally fixed-point to match the FPGA implementation.

    // GEPCellsTowers
    std::vector<double>* gepCellsTowersEtValues  = nullptr;
    std::vector<double>* gepCellsTowersEtaValues = nullptr;
    std::vector<double>* gepCellsTowersPhiValues = nullptr;

    // GEPCellsTowers (PU-suppressed)
    std::vector<double>* gepCellsTowersSKEtValues  = nullptr;
    std::vector<double>* gepCellsTowersSKEtaValues = nullptr;
    std::vector<double>* gepCellsTowersSKPhiValues = nullptr;

    // WTA-cone jets built from towers
    std::vector<double>*       gepWTAConeCellsTowersJetsEtValues           = nullptr;
    std::vector<double>*       gepWTAConeCellsTowersJetsEtaValues          = nullptr;
    std::vector<double>*       gepWTAConeCellsTowersJetsPhiValues          = nullptr;
    std::vector<unsigned int>* gepWTAConeCellsTowersJetsNConstituentsValues = nullptr;

    // WTA-cone jets built from towers (PU-suppressed)
    std::vector<double>*       gepWTAConeCellsTowersSKJetsEtValues           = nullptr;
    std::vector<double>*       gepWTAConeCellsTowersSKJetsEtaValues          = nullptr;
    std::vector<double>*       gepWTAConeCellsTowersSKJetsPhiValues          = nullptr;
    std::vector<unsigned int>* gepWTAConeCellsTowersSKJetsNConstituentsValues = nullptr;

    // GEPCellsTowers (EtaSK PU-suppressed)
    std::vector<double>*       gepCellsTowersEtaSKEtValues  = nullptr;
    std::vector<double>*       gepCellsTowersEtaSKEtaValues = nullptr;
    std::vector<double>*       gepCellsTowersEtaSKPhiValues = nullptr;

    // WTA-cone jets built from towers (EtaSK PU-suppressed)
    std::vector<double>*       gepWTAConeCellsTowersEtaSKJetsEtValues           = nullptr;
    std::vector<double>*       gepWTAConeCellsTowersEtaSKJetsEtaValues          = nullptr;
    std::vector<double>*       gepWTAConeCellsTowersEtaSKJetsPhiValues          = nullptr;
    std::vector<unsigned int>* gepWTAConeCellsTowersEtaSKJetsNConstituentsValues = nullptr;

    // Run/event number passthrough (read from input, written to emulEventInfoTree for ordering validation)
    int gepRunNumberIn = 0, gepEventNumberIn = 0;
    int gepRunNumberOut = 0, gepEventNumberOut = 0;

    // Open input ROOT file
    TFile* inputFile = TFile::Open(inputNTuplePath.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputNTuplePath << std::endl;
        return;
    }

    // Create new output ROOT file containing only MET emulator trees
    TFile* outputFile = TFile::Open(outputNTuplePath.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << outputNTuplePath << std::endl;
        return;
    }

    // Input TTrees
    TTree* gepCellsTowersTree                  = (TTree*)inputFile->Get("gepCellsTowersTree");
    TTree* gepCellsTowersSKTree                = (TTree*)inputFile->Get("gepCellsTowersSKTree");
    TTree* gepCellsTowersEtaSKTree             = (TTree*)inputFile->Get("gepCellsTowersEtaSKTree");
    TTree* gepWTAConeCellsTowersJetsTree       = (TTree*)inputFile->Get("gepWTAConeCellsTowersJetsTree");
    TTree* gepWTAConeCellsTowersSKJetsTree     = (TTree*)inputFile->Get("gepWTAConeCellsTowersSKJetsTree");
    TTree* gepWTAConeCellsTowersEtaSKJetsTree  = (TTree*)inputFile->Get("gepWTAConeCellsTowersEtaSKJetsTree");
    TTree* eventInfoTreeIn                     = (TTree*)inputFile->Get("eventInfoTree");
    eventInfoTreeIn->SetBranchAddress("gepRunNumberOut",   &gepRunNumberIn);
    eventInfoTreeIn->SetBranchAddress("gepEventNumberOut", &gepEventNumberIn);

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
    // Event/run number passthrough tree (one entry per event, for ordering-alignment validation)
    TTree* emulEventInfoTree = new TTree("emulEventInfoTree", "Run/event number passthrough from HERNTupler input for ordering validation");
    emulEventInfoTree->Branch("gepRunNumberOut",   &gepRunNumberOut);
    emulEventInfoTree->Branch("gepEventNumberOut", &gepEventNumberOut);

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
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Et",           &gepWTAConeCellsTowersJetsEtValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Eta",          &gepWTAConeCellsTowersJetsEtaValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Phi",          &gepWTAConeCellsTowersJetsPhiValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("NConstituents",&gepWTAConeCellsTowersJetsNConstituentsValues);

    // === gepWTAConeCellsTowersSKJetsTree (PU-suppressed) ===
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Et",           &gepWTAConeCellsTowersSKJetsEtValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Eta",          &gepWTAConeCellsTowersSKJetsEtaValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Phi",          &gepWTAConeCellsTowersSKJetsPhiValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("NConstituents",&gepWTAConeCellsTowersSKJetsNConstituentsValues);

    // === gepCellsTowersEtaSKTree (EtaSK PU-suppressed) ===
    gepCellsTowersEtaSKTree->SetBranchAddress("Et",  &gepCellsTowersEtaSKEtValues);
    gepCellsTowersEtaSKTree->SetBranchAddress("Eta", &gepCellsTowersEtaSKEtaValues);
    gepCellsTowersEtaSKTree->SetBranchAddress("Phi", &gepCellsTowersEtaSKPhiValues);

    // === gepWTAConeCellsTowersEtaSKJetsTree (EtaSK PU-suppressed) ===
    gepWTAConeCellsTowersEtaSKJetsTree->SetBranchAddress("Et",           &gepWTAConeCellsTowersEtaSKJetsEtValues);
    gepWTAConeCellsTowersEtaSKJetsTree->SetBranchAddress("Eta",          &gepWTAConeCellsTowersEtaSKJetsEtaValues);
    gepWTAConeCellsTowersEtaSKJetsTree->SetBranchAddress("Phi",          &gepWTAConeCellsTowersEtaSKJetsPhiValues);
    gepWTAConeCellsTowersEtaSKJetsTree->SetBranchAddress("NConstituents",&gepWTAConeCellsTowersEtaSKJetsNConstituentsValues);

    // Verifying constants
    //std::cout << "half_pi_digitized_in_phi_: " << half_pi_digitized_in_phi_ << "\n";
    //std::cout << "pi_digitized_in_phi_: " << pi_digitized_in_phi_ << "\n";

    // Event loop
    unsigned int eventsToProcess = gepCellsTowersTree->GetEntries();
    for (unsigned int iEvt = 0; iEvt < eventsToProcess; iEvt++) {
        //std::cout << "-------------------------------------" << "\n";
        //std::cout << "iEvt: " << iEvt << "\n";
        //std::cout << "-------------------------------------" << "\n";
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

        // Run/event number passthrough for ordering-alignment validation
        eventInfoTreeIn->GetEntry(iEvt);
        gepRunNumberOut   = gepRunNumberIn;
        gepEventNumberOut = gepEventNumberIn;

        // Load relevant input trees for this event
        if (useEtaSKObjects) {
            gepCellsTowersEtaSKTree->GetEntry(iEvt);
            gepWTAConeCellsTowersEtaSKJetsTree->GetEntry(iEvt);
        } else if (useSKObjects) {
            gepCellsTowersSKTree->GetEntry(iEvt);
            gepWTAConeCellsTowersSKJetsTree->GetEntry(iEvt);
        } else {
            gepCellsTowersTree->GetEntry(iEvt);
            gepWTAConeCellsTowersJetsTree->GetEntry(iEvt);
        }

        // Select tower collection based on PU suppression mode
        const std::vector<double>* towerEtVec  = useEtaSKObjects ? gepCellsTowersEtaSKEtValues  : (useSKObjects ? gepCellsTowersSKEtValues  : gepCellsTowersEtValues);
        const std::vector<double>* towerPhiVec = useEtaSKObjects ? gepCellsTowersEtaSKPhiValues : (useSKObjects ? gepCellsTowersSKPhiValues : gepCellsTowersPhiValues);
        const std::vector<double>* towerEtaVec = useEtaSKObjects ? gepCellsTowersEtaSKEtaValues : (useSKObjects ? gepCellsTowersSKEtaValues : gepCellsTowersEtaValues);

        // Clamp to maxTowersConsidered_ (runtime parameter) then to available count
        unsigned int towersProcessed = maxTowersConsidered_;
        if (towersProcessed > towerEtVec->size()) towersProcessed = towerEtVec->size();

        // -------------------------------------------------------
        // ---                   MET algorithm:                ---
        // -------------------------------------------------------

        // Select jet collection based on PU suppression mode
        const std::vector<double>* jetPtVec  = useEtaSKObjects ? gepWTAConeCellsTowersEtaSKJetsEtValues   : (useSKObjects ? gepWTAConeCellsTowersSKJetsEtValues  : gepWTAConeCellsTowersJetsEtValues);
        const std::vector<double>* jetPhiVec = useEtaSKObjects ? gepWTAConeCellsTowersEtaSKJetsPhiValues  : (useSKObjects ? gepWTAConeCellsTowersSKJetsPhiValues : gepWTAConeCellsTowersJetsPhiValues);
        const std::vector<double>* jetEtaVec = useEtaSKObjects ? gepWTAConeCellsTowersEtaSKJetsEtaValues  : (useSKObjects ? gepWTAConeCellsTowersSKJetsEtaValues : gepWTAConeCellsTowersJetsEtaValues);

        std::vector<double> jetPhiOverlapVec;
        std::vector<double> jetEtaOverlapVec;

        // Clamp to available count
        unsigned int jetsProcessed = maxJetsConsidered_;
        if (jetsProcessed > jetPtVec->size()) jetsProcessed = jetPtVec->size();

        // Jet loop (for JetMetX/JetMetY)
        int jetETxSum = 0;
        int jetETySum = 0;
        unsigned int jetSumEt = 0; 
        for (unsigned int iJet = 0; iJet < jetsProcessed; iJet++) {
            if(jetPtVec->at(iJet) <= jetEtThreshold) continue;
            unsigned int jetEt  = digitize(jetPtVec->at(iJet),  et_bit_length_,  static_cast<double>(et_min_),  static_cast<double>(et_max_));
            unsigned int jetPhi = digitize_phi(jetPhiVec->at(iJet));
            unsigned int jetEta = digitize(jetEtaVec->at(iJet), eta_bit_length_, eta_min_, eta_max_);

            jetPhiOverlapVec.push_back(jetPhi);
            jetEtaOverlapVec.push_back(jetEta);

            jetSumEt += jetEt;

            int jetCosPhi = sinLUT_[wrapPhiUnsigned(jetPhi + half_pi_digitized_in_phi_)];
            int jetSinPhi = sinLUT_[jetPhi];

            int jetETx = (static_cast<int>(jetEt) * jetCosPhi) / (1 << (sin_bit_length_ - 1));
            int jetETy = (static_cast<int>(jetEt) * jetSinPhi) / (1 << (sin_bit_length_ - 1));
            jetETxSum += jetETx;
            jetETySum += jetETy;
        }
        int jetMETx = -jetETxSum; // Take negative of ET sums
        int jetMETy = -jetETySum; // Take negative of ET sums

        // Tower loop (for TowerMetX/TowerMetY and SumET)
        int towerETxSum = 0;
        int towerETySum = 0;
        unsigned int towerSumEt = 0; 
        for (unsigned int iTower = 0; iTower < towersProcessed; iTower++) {
            //std::cout << "iTower: " << iTower << "\n";
            //std::cout << "towerE_T:"  << towerEtVec->at(iTower) << " , towerEtThreshold: " << towerEtThreshold << "\n";
            if(towerEtVec->at(iTower) <= towerEtThreshold) continue;
            //std::cout << "towerET (GeV): " << towerEtVec->at(iTower) << "\n";
            //std::cout << "towerPhi (undigi): " << towerPhiVec->at(iTower) << "\n";
            unsigned int towerEt  = digitize(towerEtVec->at(iTower),  et_bit_length_,  static_cast<double>(et_min_),  static_cast<double>(et_max_));
            unsigned int towerPhi = digitize_phi(towerPhiVec->at(iTower));
            unsigned int towerEta = digitize(towerEtaVec->at(iTower), eta_bit_length_, eta_min_, eta_max_);
            // Check for overlap between jets and towers, if towers overlap, remove from computation of tower MET
           
            if(doJetTowerOverlapRemoval){
                bool foundTowerJetOverlap = false;
                for(unsigned int iOverlapJet = 0; iOverlapJet < jetPhiOverlapVec.size(); iOverlapJet++){
                    //std::cout << "iOverlapJet: " << iOverlapJet << "\n";
                    //std::cout << "jet eta: " << jetEtaOverlapVec[iOverlapJet] << " , jet phi: " << jetPhiOverlapVec[iOverlapJet] << "\n";
                    //std::cout << "towerEta: " << towerEta << " , towerPhi: " << towerPhi << "\n";
                    unsigned int digitizedDeltaR2TowerJet = digitizedDeltaR2(jetEtaOverlapVec[iOverlapJet], jetPhiOverlapVec[iOverlapJet], towerEta, towerPhi); 
                    //std::cout << "digitizedDeltaR2TowerJet: " << digitizedDeltaR2TowerJet << "\n";
                    //std::cout << "digitized_delta_R2Cut_: " << digitized_delta_R2Cut_ << "\n";
                    if(digitizedDeltaR2TowerJet <= digitized_delta_R2Cut_){
                        foundTowerJetOverlap = true;
                        break;
                    } 
                }
                if(foundTowerJetOverlap){
                    //std::cout << "found tower-jet overlap" << "\n";
                    continue;
                }  
            }            
            
            //std::cout << "towerET: " << towerEt << "\n";
            //std::cout << "towerPhi: " << towerPhi << "\n";
            towerSumEt += towerEt;

            int towerCosPhi = sinLUT_[wrapPhiUnsigned(towerPhi + half_pi_digitized_in_phi_)];
            int towerSinPhi = sinLUT_[towerPhi];

            //std::cout << "towerCosPhi: " << towerCosPhi << "\n";
            //std::cout << "towerSinPhi: " << towerSinPhi << "\n";

            int towerETx = (static_cast<int>(towerEt) * towerCosPhi) / (1 << (sin_bit_length_ - 1));
            int towerETy = (static_cast<int>(towerEt) * towerSinPhi) / (1 << (sin_bit_length_ - 1));
            //std::cout << "towerETx: " << towerETx << "\n";
            //std::cout << "towerETy: " << towerETy << "\n";
            towerETxSum += towerETx;
            towerETySum += towerETy;
        }
        int towerMETx = -towerETxSum; // Take negative of ET sums 
        int towerMETy = -towerETySum; // Take negative of ET sums 
        //std::cout << "towerETxsum: " << towerETxSum << " , towerMETx: " << towerMETx << "\n";
        //std::cout << "towerETysum: " << towerETySum << " , towerMETy: " << towerMETy << "\n";

        // Scalar tower/jet weighting. TODO: replace with eta-dependent, calibrated factors.
        
        int totalMETx = static_cast<int>(std::lround(towerScaleFactor * towerMETx + jetScaleFactor * jetMETx));
        int totalMETy = static_cast<int>(std::lround(towerScaleFactor * towerMETy + jetScaleFactor * jetMETy));
        //std::cout << "towerScaleFactor: " << towerScaleFactor << " , jetScaleFactor: " << jetScaleFactor << "\n";
        //std::cout << "towerMETx: " << towerMETx << ", jetMETx: " << jetMETx << ", totalMETx: " << totalMETx << "\n";
        //std::cout << "towerMETy: " << towerMETy << ", jetMETy: " << jetMETy << ", totalMETy: " << totalMETy << "\n";
        //std::cout << "totalMETx: " << totalMETx << " , totalMETy: " << totalMETy << "\n";
        unsigned int towerMET   = static_cast<unsigned int>(std::sqrt(towerMETx * towerMETx + towerMETy * towerMETy));
        unsigned int jetMET     = static_cast<unsigned int>(std::sqrt(jetMETx   * jetMETx   + jetMETy   * jetMETy));
        unsigned int totalMET   = static_cast<unsigned int>(std::sqrt(totalMETx * totalMETx + totalMETy * totalMETy));
        //std::cout << "towerMET: " << towerMET << " , jetMET: " << jetMET << " , totalMET: " << totalMET << "\n";
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
        //std::cout << "out_towerMetX: " << out_towerMetX << "\n";
        out_towerMetY  = undigitize_signed_et(towerMETy_bitset);
        //std::cout << "out_towerMetY: " << out_towerMetY << "\n";
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
        emulEventInfoTree->Fill();
    } // Event loop

    outputFile->cd();
    std::cout << "Writing output file\n";
    metTree->Write("", TObject::kOverwrite);
    emulEventInfoTree->Write("", TObject::kOverwrite);
    outputFile->Close();
    inputFile->Close();

} // eventLoop

// Main MET emulation function
// Use: root -b -l -q 'metEmulation.cc+(true, false, "ZvvHbb")'
void metEmulation(bool signalBool,                 // true = signal sample, false = background (dijet)
                  bool useSKObjects,               // true = use PU-suppressed (SoftKiller) objects
                  std::string signalString = "ZvvHbb",        // Which signal sample: "VBF_hh_bbbb_cvv0/1", "ggF_hh_bbbb",
                                                              //   "ZvvHbb", "ttbar_had", "Zprime_ttbar", "ttbar_semilep", "ttbar_lep"
                  double jetEtThreshold = 20.0,              // Minimum jet E_T [GeV] included in jet MET sum
                  bool doJetTowerOverlapRemoval = false,      // Remove towers overlapping with jets
                  double towerEtThreshold = 0.0,             // Minimum tower E_T [GeV] included in tower MET sum
                  bool useEtaSKObjects = false,              // true = use EtaSK PU-suppressed objects (towers + WTAConeJets)
                  std::string explicitInputPath = "",        // When non-empty, overrides makeInputFileName (used for per-file Condor parallelism)
                  int fileIndex = -1,                        // When >= 0, appended as _fileN to output name to avoid collisions across parallel jobs
                  double towerScaleFactor = 1.0,             // Scalar weight applied to tower MET in the totalMET sum (future: eta-binned, calibrated)
                  double jetScaleFactor = 1.0,               // Scalar weight applied to jet MET in the totalMET sum (future: eta-binned, calibrated)
                  unsigned int pileup = 200                  // Pileup scenario of the input sample; tags the output name (r16130 = PU200, r16129 = PU140)
                  ) {

    if (signalBool) std::cout << "Processing signal: " << signalString << "\n";

    auto infile  = explicitInputPath.empty() ? makeInputFileName(signalBool, signalString, "/data/larsonma/GEPHadronicEventReconstruction/ntuples/", pileup) : explicitInputPath;
    auto outfile = makeOutputMETFileName(maxTowersConsidered_, signalBool, signalString, useSKObjects, jetEtThreshold, towerEtThreshold, doJetTowerOverlapRemoval, "/data/larsonma/GEPMET/outputNTuplesDev_METv2/", useEtaSKObjects, towerScaleFactor, jetScaleFactor, pileup);
    if (fileIndex >= 0) {
        size_t pos = outfile.rfind(".root");
        if (pos != std::string::npos)
            outfile = outfile.substr(0, pos) + "_file" + std::to_string(fileIndex) + ".root";
    }

    std::cout << "infile:  " << infile  << "\n";
    std::cout << "outfile: " << outfile << "\n";

    gSystem->RedirectOutput("debuglog_MET.log", "w");
    std::cout << "Calling event loop\n";
    eventLoop(infile, outfile, useSKObjects, jetEtThreshold, towerEtThreshold, doJetTowerOverlapRemoval, useEtaSKObjects, towerScaleFactor, jetScaleFactor);
}


// =====================================================================================
// ===                        MET baseline validation path                           ===
// =====================================================================================
// Text-file in, text-file out. Reads the HERNTupler memory prints for towers and jets,
// runs the same digitized MET arithmetic the event loop above runs, and writes the
// result back out in the same memory-print format so an HLS C-sim testbench can be
// diffed against it object-for-object.
//
// Nothing here reads a ROOT file. The memory prints are already digitized, so the
// digitize() / digitize_phi() calls of the event loop are deliberately absent: the
// codes on disk are used as-is. That is the point of the path -- it isolates the MET
// arithmetic from the digitization, so a disagreement with the firmware can only come
// from the arithmetic.
//
// The baseline configuration is fixed on purpose and is not exposed as a parameter:
// no tower E_T threshold, no jet E_T threshold, no jet/tower overlap removal, and unit
// scale factors on both tower and jet MET. Anything else would make the reference
// vectors depend on a tuning choice that the firmware does not know about.
// =====================================================================================

// Set to true to print the per-event digitized quantities while writing the vectors.
constexpr bool met_baseline_debug_ = false;

const std::string met_baseline_mem_prints_path_ = "/eos/home-m/mlarson/TransferMemPrintsLUTs/data/MemPrints_v3/";
const std::string met_baseline_output_subdir_   = "METBaselineValidation/";
const std::string met_baseline_tower_output_suffix_ = "_METBaselineTower.dat";
const std::string met_baseline_jet_output_suffix_   = "_METBaselineJet.dat";
const std::string met_baseline_tower_subdir_    = "GEPCellsTowersSK/";
const std::string met_baseline_tower_suffix_    = "_gepcellstowerssk.dat";
const std::string met_baseline_jet_subdir_      = "GEPConeJetsCellsTowersSK/";
const std::string met_baseline_jet_suffix_      = "_gepconejetscellstowerssk.dat";
const std::string met_baseline_log_prefix_      = "debuglog_METBaselineValidation_";
const std::string met_baseline_log_suffix_      = ".log";

// Output word, 64 bits, laid out to match the Jet MET TOB map in the firmware
// (MET_Engine.v), MSB -> LSB:
//
//   [63:58]  reserved       6 bits, zero
//   [57:45]  ey_total      13 bits, signed accumulated Y component
//   [44:32]  ex_total      13 bits, signed accumulated X component
//   [31:23]  phi_reserved   9 bits, zero (reserved for phi_missing)
//   [22:13]  eta_reserved  10 bits, zero
//   [12:0]   et_missing    13 bits, MET after the square root
//
// The three reserved fields are written as zero but are still emitted in the binary
// column so a line reads the same way as the firmware word. Their widths are the
// firmware's, not this emulator's: the 9-bit phi and 10-bit eta of the TOB are wider
// than the phi_bit_length_ / eta_bit_length_ codes used on the GEP tower grid, so they
// are spelled out here rather than derived, and must not be swapped for the grid widths.
//
// ex_total / ey_total are sign-magnitude over signed_et_bit_length_ bits, matching
// pack_signed_et in the event loop above and what undigitize_signed_et expects.
//
// NOTE: the field layout is the firmware's but the Ex/Ey *encoding* is not -- MET_Engine.v
// carries them as two's complement. Sign-magnitude is the deliberate choice here so the
// validation path and the ROOT path agree with each other; a consumer diffing these
// vectors against firmware or C-sim output has to convert negative components first
// (-59 is 0x103B here, 0x1FC5 in the firmware). Positive values are identical either way.
constexpr unsigned int met_baseline_et_missing_width_   = et_bit_length_;        // 13
constexpr unsigned int met_baseline_eta_reserved_width_ = 10;
constexpr unsigned int met_baseline_phi_reserved_width_ = 9;
constexpr unsigned int met_baseline_ex_total_width_     = signed_et_bit_length_; // 13
constexpr unsigned int met_baseline_ey_total_width_     = signed_et_bit_length_; // 13
constexpr unsigned int met_baseline_reserved_width_     = 6;

constexpr unsigned int met_baseline_et_missing_low_   = 0;
constexpr unsigned int met_baseline_eta_reserved_low_ = met_baseline_et_missing_low_   + met_baseline_et_missing_width_;
constexpr unsigned int met_baseline_phi_reserved_low_ = met_baseline_eta_reserved_low_ + met_baseline_eta_reserved_width_;
constexpr unsigned int met_baseline_ex_total_low_     = met_baseline_phi_reserved_low_ + met_baseline_phi_reserved_width_;
constexpr unsigned int met_baseline_ey_total_low_     = met_baseline_ex_total_low_     + met_baseline_ex_total_width_;
constexpr unsigned int met_baseline_reserved_low_     = met_baseline_ey_total_low_     + met_baseline_ey_total_width_;
constexpr unsigned int met_baseline_total_bits_       = met_baseline_reserved_low_     + met_baseline_reserved_width_;
static_assert(met_baseline_total_bits_ == 64, "MET baseline output word must be exactly 64 bits");

// One object as it appears in a memory-print line: already-digitized codes.
struct MemPrintObject {
    unsigned int phi = 0;
    unsigned int eta = 0;
    unsigned int et  = 0;
};

// Read a whole memory-print file into a per-event list of objects.
//
// The binary column is treated as authoritative rather than the trailing hex word, for
// the same reason extract_values_from_file in fileRead.h does: it is the column that
// carries the field boundaries. Events are taken in file order, so index i of the
// returned vector is the i-th "Event :" block.
inline std::vector<std::vector<MemPrintObject>> readMemPrintFile(const std::string& fileName) {
    std::vector<std::vector<MemPrintObject>> eventObjects;

    std::ifstream inFile(fileName);
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open memory-print file " << fileName << std::endl;
        return eventObjects;
    }

    std::string line;
    while (std::getline(inFile, line)) {
        if (line.find("Event") != std::string::npos) {
            eventObjects.emplace_back();
            continue;
        }
        if (eventObjects.empty()) continue; // stray line before the first event header

        std::stringstream ss(line);
        std::string index, bin, hexWord;
        if (!(ss >> index >> bin >> hexWord)) continue; // blank or short line

        size_t firstPipe  = bin.find('|');
        size_t secondPipe = bin.rfind('|');
        if (firstPipe == std::string::npos || secondPipe == std::string::npos || firstPipe == secondPipe) {
            std::cerr << "Error: Malformed memory-print line -> " << line << std::endl;
            continue;
        }

        // Binary column is written MSB -> LSB as phi | eta | et
        MemPrintObject object;
        object.phi = static_cast<unsigned int>(std::stoul(bin.substr(0, firstPipe), nullptr, 2));
        object.eta = static_cast<unsigned int>(std::stoul(bin.substr(firstPipe + 1, secondPipe - firstPipe - 1), nullptr, 2));
        object.et  = static_cast<unsigned int>(std::stoul(bin.substr(secondPipe + 1), nullptr, 2));
        eventObjects.back().push_back(object);
    }

    inFile.close();
    return eventObjects;
}

// Physics part of the memory-print file name, shared by the tower and jet inputs and by
// the validation output. Mirrors fileName_ in fileRead.h so the testbench and this
// emulator name the same sample the same way.
inline std::string makeMemPrintBaseName(bool signalBool, std::string signalString, unsigned int jzSlice) {
    if (signalBool) {
        if (signalString == "ggF_hh_bbbb") return "mc21_14TeV_HHbbbb_HLLHC";
        std::cerr << "Error: no memory-print base name known for signal " << signalString << std::endl;
        return "";
    }
    if (jzSlice == 2) return "mc21_14TeV_jj_JZ2";
    if (jzSlice == 3) return "mc21_14TeV_jj_JZ3";
    if (jzSlice == 4) return "mc21_14TeV_jj_JZ4";
    std::cerr << "Error: no memory-print base name known for JZ slice " << jzSlice << std::endl;
    return "";
}

// Sign-magnitude packing, same convention as pack_signed_et in the event loop above.
inline uint64_t packSignedEtMemPrint(int value) {
    uint64_t sign = (value < 0) ? 1u : 0u;
    uint64_t mag  = static_cast<uint64_t>(std::abs(value)) & maskN(signed_et_bit_length_ - 1);
    return (sign << (signed_et_bit_length_ - 1)) | mag;
}

// Accumulate one collection of already-digitized objects into MET components.
// Identical arithmetic to the tower and jet loops of eventLoop: the sine lookup is
// indexed by the phi code, cosine is the same table read half_pi_digitized_in_phi_
// codes along, and the product is brought back to E_T counts by an integer divide by
// 1 << (sin_bit_length_ - 1). MET is the negative of the vector E_T sum.
inline void accumulateMemPrintMet(const std::vector<MemPrintObject>& objects,
                                  unsigned int maxObjects,
                                  int& metX, int& metY, unsigned int& met) {
    int etXSum = 0;
    int etYSum = 0;

    unsigned int objectsProcessed = maxObjects;
    if (objectsProcessed > objects.size()) objectsProcessed = objects.size();

    // No E_T threshold and no overlap removal by design -- every object in the memory
    // print contributes.
    for (unsigned int iObject = 0; iObject < objectsProcessed; iObject++) {
        const MemPrintObject& object = objects.at(iObject);

        int cosPhi = sinLUT_[wrapPhiUnsigned(object.phi + half_pi_digitized_in_phi_)];
        int sinPhi = sinLUT_[object.phi];

        etXSum += (static_cast<int>(object.et) * cosPhi) / (1 << (sin_bit_length_ - 1));
        etYSum += (static_cast<int>(object.et) * sinPhi) / (1 << (sin_bit_length_ - 1));
    }

    metX = -etXSum;
    metY = -etYSum;
    met  = static_cast<unsigned int>(std::sqrt(static_cast<double>(metX) * metX + static_cast<double>(metY) * metY));
}

// Write one memory-print line: index, the six binary fields of the TOB word, and the
// packed hex word. The reserved fields are all zero and contribute nothing to the hex,
// but are printed so the binary column mirrors the firmware layout field for field.
inline void writeMetBaselineLine(std::ofstream& outFile, unsigned int index,
                                 int metX, int metY, unsigned int met) {
    uint64_t exPacked      = packSignedEtMemPrint(metX);
    uint64_t eyPacked      = packSignedEtMemPrint(metY);
    uint64_t etMissingBits = static_cast<uint64_t>(met & maskN(met_baseline_et_missing_width_));

    std::bitset<met_baseline_reserved_width_>     reserved_bitset(0);
    std::bitset<met_baseline_ey_total_width_>     ey_total_bitset(eyPacked);
    std::bitset<met_baseline_ex_total_width_>     ex_total_bitset(exPacked);
    std::bitset<met_baseline_phi_reserved_width_> phi_reserved_bitset(0);
    std::bitset<met_baseline_eta_reserved_width_> eta_reserved_bitset(0);
    std::bitset<met_baseline_et_missing_width_>   et_missing_bitset(etMissingBits);

    uint64_t combined_value =
        (eyPacked      << met_baseline_ey_total_low_)   |
        (exPacked      << met_baseline_ex_total_low_)   |
        (etMissingBits << met_baseline_et_missing_low_);

    std::stringstream hex_stream;
    hex_stream << std::hex << std::nouppercase << std::setfill('0') << std::setw(met_baseline_total_bits_ / 4) << combined_value;

    outFile << "0x" << std::hex << std::setw(2) << std::setfill('0') << index << std::dec << std::setfill(' ')
            << " " << reserved_bitset.to_string()
            << "|" << ey_total_bitset.to_string()
            << "|" << ex_total_bitset.to_string()
            << "|" << phi_reserved_bitset.to_string()
            << "|" << eta_reserved_bitset.to_string()
            << "|" << et_missing_bitset.to_string()
            << " 0x" << hex_stream.str() << std::endl;
}

// MET baseline validation entry point.
//
// Use:
//   root -b -l
//   root [0] .L metEmulation.cc+
//   root [1] metBaselineValidation(true, "ggF_hh_bbbb")
//
// Writes two files per sample to met_baseline_output_subdir_, one for tower MET and one
// for jet MET, each with a single 0x00 line per event. Kept apart rather than as two
// lines of one file so each converts on its own to whatever the consumer wants without
// having to demultiplex by object index.
//
// Total MET is not written: with unit scale factors it is just the component-wise sum of
// the tower and jet words, so writing it would bake a redundant value into the reference.
//
// Per-event printouts go to met_baseline_log_prefix_ + sample + met_baseline_log_suffix_
// in the working directory; only the paths and the closing summary reach the terminal.
void metBaselineValidation(bool signalBool = true,          // true = signal sample, false = dijet background
                           std::string signalString = "ggF_hh_bbbb", // Which signal sample (memory prints exist for ggF_hh_bbbb)
                           unsigned int jzSlice = 3,        // JZ slice, used only when signalBool is false
                           int maxEvents = -1) {            // Process at most this many events; < 0 processes all

    std::string baseName = makeMemPrintBaseName(signalBool, signalString, jzSlice);
    if (baseName.empty()) return;

    std::string towerFile       = met_baseline_mem_prints_path_ + met_baseline_tower_subdir_ + baseName + met_baseline_tower_suffix_;
    std::string jetFile         = met_baseline_mem_prints_path_ + met_baseline_jet_subdir_   + baseName + met_baseline_jet_suffix_;
    std::string outputPath      = met_baseline_mem_prints_path_ + met_baseline_output_subdir_;
    std::string towerOutputFile = outputPath + baseName + met_baseline_tower_output_suffix_;
    std::string jetOutputFile   = outputPath + baseName + met_baseline_jet_output_suffix_;
    // Log stays in the working directory rather than next to the vectors on eos: it is a
    // run artifact, not a reference artifact, and it is rewritten on every run.
    std::string logFile         = met_baseline_log_prefix_ + baseName + met_baseline_log_suffix_;

    std::cout << "tower memory print: " << towerFile       << "\n";
    std::cout << "jet memory print:   " << jetFile         << "\n";
    std::cout << "tower output:       " << towerOutputFile << "\n";
    std::cout << "jet output:         " << jetOutputFile   << "\n";

    std::vector<std::vector<MemPrintObject>> towerEvents = readMemPrintFile(towerFile);
    std::vector<std::vector<MemPrintObject>> jetEvents   = readMemPrintFile(jetFile);

    if (towerEvents.empty() || jetEvents.empty()) {
        std::cerr << "Error: no events read, nothing written" << std::endl;
        return;
    }
    if (towerEvents.size() != jetEvents.size()) {
        // Not fatal, but the two prints are supposed to come from the same HERNTupler
        // pass over the same events, so a mismatch means they are not the same sample.
        std::cerr << "Warning: tower print has " << towerEvents.size() << " events, jet print has "
                  << jetEvents.size() << "; processing the overlap only" << std::endl;
    }

    unsigned int eventsToProcess = static_cast<unsigned int>(std::min(towerEvents.size(), jetEvents.size()));
    if (maxEvents >= 0 && static_cast<unsigned int>(maxEvents) < eventsToProcess)
        eventsToProcess = static_cast<unsigned int>(maxEvents);

    gSystem->mkdir(outputPath.c_str(), true);
    std::ofstream towerOutFile(towerOutputFile);
    if (!towerOutFile.is_open()) {
        std::cerr << "Error: Could not open file " << towerOutputFile << std::endl;
        return;
    }
    std::ofstream jetOutFile(jetOutputFile);
    if (!jetOutFile.is_open()) {
        std::cerr << "Error: Could not open file " << jetOutputFile << std::endl;
        return;
    }

    // Send the per-event printouts to a log rather than the terminal.
    //
    // The handle form is what makes this restorable: RedirectOutput saves the current
    // stdout/stderr into redirectHandle, and the paired call below with a null name and
    // the same handle puts them back. Redirecting without a handle -- as the ROOT path in
    // metEmulation above does -- leaves the rest of the session writing into the log file,
    // which is why nothing printed after it ever reaches the terminal.
    //
    // The redirect brackets the event loop only. Every early return sits above it, so
    // there is no path that leaves stdout redirected.
    RedirectHandle_t redirectHandle;
    gSystem->RedirectOutput(logFile.c_str(), "w", &redirectHandle);

    std::cout << "MET baseline validation: " << baseName << "\n";
    std::cout << "tower memory print: " << towerFile       << "\n";
    std::cout << "jet memory print:   " << jetFile         << "\n";
    std::cout << "tower output:       " << towerOutputFile << "\n";
    std::cout << "jet output:         " << jetOutputFile   << "\n";
    std::cout << "events:             " << std::dec << eventsToProcess << "\n";

    for (unsigned int iEvt = 0; iEvt < eventsToProcess; iEvt++) {
        towerOutFile << "Event : " << std::dec << iEvt << std::endl;
        jetOutFile   << "Event : " << std::dec << iEvt << std::endl;

        int towerMETx = 0, towerMETy = 0;
        unsigned int towerMET = 0;
        accumulateMemPrintMet(towerEvents.at(iEvt), maxTowersConsidered_, towerMETx, towerMETy, towerMET);

        int jetMETx = 0, jetMETy = 0;
        unsigned int jetMET = 0;
        accumulateMemPrintMet(jetEvents.at(iEvt), maxJetsConsidered_, jetMETx, jetMETy, jetMET);

        // Echo exactly the quantities that go into the output words, in digitized counts
        // with the GeV equivalent alongside, so a line of the .dat can be read back
        // without decoding it by hand.
        std::cout << "iEvt: " << std::dec << iEvt
                  << "  towers: " << towerEvents.at(iEvt).size()
                  << "  jets: "   << jetEvents.at(iEvt).size() << "\n";
        std::cout << "  tower  Ex: " << towerMETx << " (" << towerMETx * et_granularity_ << " GeV)"
                  << "  Ey: "        << towerMETy << " (" << towerMETy * et_granularity_ << " GeV)"
                  << "  MET: "       << towerMET  << " (" << towerMET  * et_granularity_ << " GeV)\n";
        std::cout << "  jet    Ex: " << jetMETx   << " (" << jetMETx   * et_granularity_ << " GeV)"
                  << "  Ey: "        << jetMETy   << " (" << jetMETy   * et_granularity_ << " GeV)"
                  << "  MET: "       << jetMET    << " (" << jetMET    * et_granularity_ << " GeV)\n";

        if (met_baseline_debug_) {
            for (unsigned int iObject = 0; iObject < towerEvents.at(iEvt).size(); iObject++) {
                const MemPrintObject& object = towerEvents.at(iEvt).at(iObject);
                std::cout << "    tower " << iObject << " phi: " << object.phi
                          << " eta: " << object.eta << " et: " << object.et
                          << "  sin: " << sinLUT_[object.phi]
                          << "  cos: " << sinLUT_[wrapPhiUnsigned(object.phi + half_pi_digitized_in_phi_)] << "\n";
            }
            for (unsigned int iObject = 0; iObject < jetEvents.at(iEvt).size(); iObject++) {
                const MemPrintObject& object = jetEvents.at(iEvt).at(iObject);
                std::cout << "    jet   " << iObject << " phi: " << object.phi
                          << " eta: " << object.eta << " et: " << object.et
                          << "  sin: " << sinLUT_[object.phi]
                          << "  cos: " << sinLUT_[wrapPhiUnsigned(object.phi + half_pi_digitized_in_phi_)] << "\n";
            }
        }

        writeMetBaselineLine(towerOutFile, 0, towerMETx, towerMETy, towerMET);
        writeMetBaselineLine(jetOutFile,   0, jetMETx,   jetMETy,   jetMET);
    }

    // Restore stdout/stderr before the summary, so the summary lands on the terminal.
    gSystem->RedirectOutput(nullptr, "", &redirectHandle);

    towerOutFile.close();
    jetOutFile.close();
    std::cout << "Wrote " << std::dec << eventsToProcess << " events to " << towerOutputFile << "\n";
    std::cout << "Wrote " << std::dec << eventsToProcess << " events to " << jetOutputFile   << "\n";
    std::cout << "Per-event printouts in " << logFile << "\n";
}
