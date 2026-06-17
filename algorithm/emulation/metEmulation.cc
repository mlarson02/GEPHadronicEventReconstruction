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
            unsigned int jetPhi = digitize(jetPhiVec->at(iJet), phi_bit_length_, phi_min_, phi_max_);
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
            unsigned int towerPhi = digitize(towerPhiVec->at(iTower), phi_bit_length_, phi_min_, phi_max_);
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
                  double jetScaleFactor = 1.0                // Scalar weight applied to jet MET in the totalMET sum (future: eta-binned, calibrated)
                  ) {

    if (signalBool) std::cout << "Processing signal: " << signalString << "\n";

    auto infile  = explicitInputPath.empty() ? makeInputFileName(signalBool, signalString) : explicitInputPath;
    auto outfile = makeOutputMETFileName(maxTowersConsidered_, signalBool, signalString, useSKObjects, jetEtThreshold, towerEtThreshold, doJetTowerOverlapRemoval, "/data/larsonma/GEPMET/outputNTuplesDev_METv2/", useEtaSKObjects, towerScaleFactor, jetScaleFactor);
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
