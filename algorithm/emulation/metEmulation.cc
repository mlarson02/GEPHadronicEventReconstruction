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
#include "TRandom3.h"
#include "emulationHelperFunctions_MET.h"

// Read the JES-calibrated WTACone CellsTowers jets (HERNTupler *CalibJetsTree, needs ntuples made with
// apply_gep_jes_calibration_ = true) instead of the uncalibrated ones, and write to their own output
// directory. Switch together with jetTaggerEmulation.cc, metAnalysisAndRates.C, largeRJetAnalysisCore.h,
// largeRJetStateRegistry.h and the condor scripts (hadd_emulator_outputs.sh).
const bool use_calibrated_wtacone_jets_ = true;
const std::string met_output_dir_ = use_calibrated_wtacone_jets_ ? "/data/larsonma/GEPMET/outputNTuplesDev_METv4_CalibJets/"
                                                                 : "/data/larsonma/GEPMET/outputNTuplesDev_METv3/";

// Input tree of the WTACone CellsTowers jets of one PU-suppression variant ("", "SK" or "EtaSK")
std::string wtaConeCellsTowersJetsTreeName(const std::string& variant) {
    return "gepWTAConeCellsTowers" + variant + (use_calibrated_wtacone_jets_ ? "CalibJetsTree" : "JetsTree");
}


// Function for processing provided number of events with MET algorithm,
// then writing output MET data to a new TTree
void eventLoop(std::string inputNTuplePath, std::string outputNTuplePath,
               bool useSKObjects, double jetEtThreshold, double towerEtThreshold, bool doJetTowerOverlapRemoval,
               bool useEtaSKObjects = false,
               double towerScaleFactor = 1.0, double jetScaleFactor = 1.0,
               bool useGEPJwoJ = false,
               GEPJwoJParams jwojParams = makeUniformGEPJwoJParams(10.0, 1.0, 1.0, 1)) {
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
    TTree* gepWTAConeCellsTowersJetsTree       = (TTree*)inputFile->Get(wtaConeCellsTowersJetsTreeName("").c_str());
    TTree* gepWTAConeCellsTowersSKJetsTree     = (TTree*)inputFile->Get(wtaConeCellsTowersJetsTreeName("SK").c_str());
    TTree* gepWTAConeCellsTowersEtaSKJetsTree  = (TTree*)inputFile->Get(wtaConeCellsTowersJetsTreeName("EtaSK").c_str());
    if (!gepWTAConeCellsTowersJetsTree || !gepWTAConeCellsTowersSKJetsTree || !gepWTAConeCellsTowersEtaSKJetsTree) {
        std::cerr << "Error: " << wtaConeCellsTowersJetsTreeName("") << " (or its SK / EtaSK variant) not found in "
                  << inputNTuplePath << (use_calibrated_wtacone_jets_ ? " -- calibrated jets need ntuples made with apply_gep_jes_calibration_ = true" : "")
                  << std::endl;
        return;
    }
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

    // GEP JwoJ MET outputs. Written only when useGEPJwoJ; see the algorithm block in the
    // event loop below. The hard and soft terms are written out UNCOEFFICIENTED, so the
    // coefficients can be re-derived downstream without re-running the emulation.
    double out_jwojMET       = 0.0, out_jwojMETX     = 0.0, out_jwojMETY     = 0.0;
    double out_jwojHardMET   = 0.0, out_jwojHardMETX = 0.0, out_jwojHardMETY = 0.0;
    double out_jwojSoftMET   = 0.0, out_jwojSoftMETX = 0.0, out_jwojSoftMETY = 0.0;
    double out_jwojSumHardET = 0.0, out_jwojSumSoftET = 0.0, out_jwojSumET    = 0.0;

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

    // GEP JwoJ branches, created ONLY in JwoJ mode. Absent rather than present-and-zero on a
    // standard run, so a consumer can tell the two apart from the tree alone -- which is what
    // metAnalysisAndRates.C's hasGEPJwoJ check relies on, alongside the filename tag.
    if (useGEPJwoJ) {
        metTree->Branch("GEPJwoJMET",       &out_jwojMET);
        metTree->Branch("GEPJwoJMETX",      &out_jwojMETX);
        metTree->Branch("GEPJwoJMETY",      &out_jwojMETY);
        metTree->Branch("GEPJwoJHardMET",   &out_jwojHardMET);
        metTree->Branch("GEPJwoJHardMETX",  &out_jwojHardMETX);
        metTree->Branch("GEPJwoJHardMETY",  &out_jwojHardMETY);
        metTree->Branch("GEPJwoJSoftMET",   &out_jwojSoftMET);
        metTree->Branch("GEPJwoJSoftMETX",  &out_jwojSoftMETX);
        metTree->Branch("GEPJwoJSoftMETY",  &out_jwojSoftMETY);
        metTree->Branch("GEPJwoJSumHardET", &out_jwojSumHardET);
        metTree->Branch("GEPJwoJSumSoftET", &out_jwojSumSoftET);
        metTree->Branch("GEPJwoJSumET",     &out_jwojSumET);
    }

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
    // --- GEP JwoJ block scratch, hoisted out of the event loop --------------------------
    // The tower collections are flat lists, so the blocks have to be accumulated into a grid
    // indexed by (eta block, phi block). Allocated once: re-allocating a 6272-entry grid on
    // every event would cost more than the tower loop it serves. Only the blocks an event
    // actually touched are reset afterwards -- jwojTouchedBlocks records them -- so the
    // per-event cost is O(occupied blocks) and not O(grid), which matters most at blockSize 1
    // where the two would otherwise differ by the tower occupancy.
    const unsigned int jwojBlockSize   = jwojParams.blockSize;
    const unsigned int nJwoJEtaBlocks  = eta_range_ / jwojBlockSize;   // exact: blockSize tiles the grid
    const unsigned int nJwoJPhiBlocks  = phi_range_ / jwojBlockSize;
    const unsigned int nJwoJBlocks     = nJwoJEtaBlocks * nJwoJPhiBlocks;
    // One entry per block: the summed digitized E_T that carries the vector sum, and the
    // summed undigitized E_T the hard threshold is compared against (the same convention the
    // per-tower threshold above uses, so a tower means the same thing to both cuts).
    std::vector<unsigned int>  jwojBlockEtCode(useGEPJwoJ ? nJwoJBlocks : 0, 0u);
    std::vector<double>        jwojBlockEtGeV (useGEPJwoJ ? nJwoJBlocks : 0, 0.0);
    std::vector<unsigned char> jwojBlockSeen  (useGEPJwoJ ? nJwoJBlocks : 0, 0u);
    std::vector<unsigned char> jwojBlockIsHard(useGEPJwoJ ? nJwoJBlocks : 0, 0u);
    std::vector<unsigned int>  jwojTouchedBlocks;
    // Each surviving tower's digitized codes and the block it landed in, kept from the first
    // pass so the soft-term pass does not digitize the same tower a second time.
    struct JwoJTowerCodes { unsigned int etCode, phiCode, etaCode, blockIdx; };
    std::vector<JwoJTowerCodes> jwojTowers;
    if (useGEPJwoJ) {
        jwojTouchedBlocks.reserve(maxTowersConsidered_);
        jwojTowers.reserve(maxTowersConsidered_);
    }

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
            // Tower-grid indices, through the jet TOB codes if input jets carry them
            unsigned int jetPhi = input_jets_use_jet_tob_encoding_ ? jet_phi_to_tower_index(digitize_jet_phi(jetPhiVec->at(iJet)))
                                                                   : digitize_phi(jetPhiVec->at(iJet));
            unsigned int jetEta = input_jets_use_jet_tob_encoding_ ? jet_eta_to_tower_index(digitize_jet_eta(jetEtaVec->at(iJet)))
                                                                   : digitize(jetEtaVec->at(iJet), eta_bit_length_, eta_min_, eta_max_, eta_range_);

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
            unsigned int towerEta = digitize(towerEtaVec->at(iTower), eta_bit_length_, eta_min_, eta_max_, eta_range_);
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

        // -------------------------------------------------------
        // ---              GEP JwoJ MET algorithm:            ---
        // -------------------------------------------------------
        // The same towers, the same tower E_T threshold and the same pileup-suppression
        // choice as the tower loop above, but split into a hard and a soft term at the
        // per-eta-bin threshold in jwojParams instead of being summed as one. The hard term
        // stands in for the jet term of the algorithm above -- it is the high-E_T part of the
        // event, reached without ever building a jet, which is what "jets without jets" means.
        //
        // Its OWN loop rather than an extra branch inside the tower loop above, because that
        // loop applies the jet/tower overlap removal: on an OR config it would silently delete
        // exactly the high-E_T towers this term is built from, and this algorithm has no jets
        // to remove against in the first place.
        std::cout << " why is this not showing up" << "\n";
        if (useGEPJwoJ) {
            std::cout << "using GEP JwoJ" << "\n";
            out_jwojMET       = 0.0; out_jwojMETX      = 0.0; out_jwojMETY     = 0.0;
            out_jwojHardMET   = 0.0; out_jwojHardMETX  = 0.0; out_jwojHardMETY = 0.0;
            out_jwojSoftMET   = 0.0; out_jwojSoftMETX  = 0.0; out_jwojSoftMETY = 0.0;
            out_jwojSumHardET = 0.0; out_jwojSumSoftET = 0.0; out_jwojSumET    = 0.0;

            const unsigned int nJwoJBins = jwojParams.nBins();
            // Accumulated per eta bin rather than as two running totals, so that a
            // per-eta-bin coefficient table can be applied at the recombination below. With
            // the single-bin default this is two totals with an extra index on them.
            std::vector<int> hardETxByBin(nJwoJBins, 0), hardETyByBin(nJwoJBins, 0);
            std::vector<int> softETxByBin(nJwoJBins, 0), softETyByBin(nJwoJBins, 0);
            unsigned int jwojHardSumEt = 0, jwojSoftSumEt = 0;

            // ---- pass 1: bin the surviving towers into blocks ----------------------------
            // A tower's block is integer division of its codes. blockSize tiles both axes
            // exactly (enforced by isSupportedGEPJwoJBlockSize), so every block is full and
            // no tower falls outside the grid.
            jwojTowers.clear();
            jwojTouchedBlocks.clear();
            for (unsigned int iTower = 0; iTower < towersProcessed; iTower++) {
                if (towerEtVec->at(iTower) <= towerEtThreshold) continue;
                unsigned int towerEt  = digitize(towerEtVec->at(iTower),  et_bit_length_,  static_cast<double>(et_min_),  static_cast<double>(et_max_));
                unsigned int towerPhi = digitize_phi(towerPhiVec->at(iTower));
                // eta_range_ as the code count, not 1 << eta_bit_length_: the block grid is
                // eta_range_ codes tall, and a 128-code eta would index past its end.
                unsigned int towerEta = digitize(towerEtaVec->at(iTower), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                std::cout << "iTower: " << iTower << "\n";
                const unsigned int blockEta = towerEta / jwojBlockSize;
                const unsigned int blockPhi = towerPhi / jwojBlockSize;
                const unsigned int blockIdx = blockEta * nJwoJPhiBlocks + blockPhi;

                std::cout << "jwojBlockSize: " << jwojBlockSize << "\n";
                std::cout << "blockEta, blocKphi: " << blockEta << " , " << blockPhi << "\n";

                if (!jwojBlockSeen[blockIdx]) {
                    jwojBlockSeen[blockIdx] = 1;
                    jwojTouchedBlocks.push_back(blockIdx);
                }
                jwojBlockEtCode[blockIdx] += towerEt;
                jwojBlockEtGeV[blockIdx]  += towerEtVec->at(iTower);
                std::cout << "jwojBlockEtCode[blockIdx]: " << jwojBlockEtCode[blockIdx] << "\n";
                jwojTowers.push_back({towerEt, towerPhi, towerEta, blockIdx});
            }

            // ---- pass 2: hard term, one contribution per block above threshold -----------
            // The block sits at its geometric centre. In HALF-codes that centre is
            // firstCode + lastCode on each axis -- 2*p0 + blockSize - 1 for a full block --
            // which is an integer whatever the block size, and whose parity says whether the
            // centre lands on a tower centre (odd blockSize) or on a tower corner (even).
            // sinLUTHalf_ covers both; see the note on it in constants.h.
            //
            // One multiply per block rather than one per tower is the point of blocking, and
            // it is also a real numerical difference from summing per-tower products: the
            // rounding of the E_T x sin/cos product now happens once on the block sum.
            for (unsigned int iBlock = 0; iBlock < jwojTouchedBlocks.size(); iBlock++) {
                const unsigned int blockIdx = jwojTouchedBlocks[iBlock];
                const unsigned int blockEta = blockIdx / nJwoJPhiBlocks;
                const unsigned int blockPhi = blockIdx % nJwoJPhiBlocks;

                const unsigned int firstEta = blockEta * jwojBlockSize;
                const unsigned int lastEta  = firstEta + jwojBlockSize - 1;
                const unsigned int firstPhi = blockPhi * jwojBlockSize;
                const unsigned int lastPhi  = firstPhi + jwojBlockSize - 1;

                // The eta bin is picked from the block's central tower code. For an even
                // block size that centre is between two towers, so integer division takes the
                // lower of them -- a fixed, documented convention, and immaterial in practice
                // since the bins are far coarser than a tower.
                const unsigned int iBin = jwojParams.binForEtaCode((firstEta + lastEta) / 2);

                if (jwojBlockEtGeV[blockIdx] <= jwojParams.hardEtThreshold[iBin]) continue;
                jwojBlockIsHard[blockIdx] = 1;

                const unsigned int blockPhiHalf = firstPhi + lastPhi;
                int blockCosPhi = sinLUTHalf_[wrapPhiHalfUnsigned(blockPhiHalf + half_pi_digitized_in_phi_half_)];
                int blockSinPhi = sinLUTHalf_[blockPhiHalf];

                const int blockEt = static_cast<int>(jwojBlockEtCode[blockIdx]);
                std::cout << "Hard E_x: " << (blockEt * blockCosPhi) / (1 << (sin_bit_length_ - 1)) << " , Hard E_y: " << (blockEt * blockSinPhi) / (1 << (sin_bit_length_ - 1)) << "\n";
                hardETxByBin[iBin] += (blockEt * blockCosPhi) / (1 << (sin_bit_length_ - 1));
                hardETyByBin[iBin] += (blockEt * blockSinPhi) / (1 << (sin_bit_length_ - 1));
                jwojHardSumEt += jwojBlockEtCode[blockIdx];
            }

            // ---- pass 3: soft term, the towers of every block that did not pass ----------
            // Per tower and at the tower's own centre, not at its block's: the soft term is
            // diffuse and there is nothing to gain by coarsening it, and this is what gFEX
            // does with the gTowers its gBlocks did not claim. Each tower's own eta bin is
            // used here, for the same reason.
            for (unsigned int iTower = 0; iTower < jwojTowers.size(); iTower++) {
                const JwoJTowerCodes& tower = jwojTowers[iTower];
                if (jwojBlockIsHard[tower.blockIdx]) continue;

                const unsigned int iBin = jwojParams.binForEtaCode(tower.etaCode);

                int towerCosPhi = sinLUT_[wrapPhiUnsigned(tower.phiCode + half_pi_digitized_in_phi_)];
                int towerSinPhi = sinLUT_[tower.phiCode];

                const int towerEt = static_cast<int>(tower.etCode);
                std::cout << "Hard E_x: " << (towerEt * towerCosPhi) / (1 << (sin_bit_length_ - 1)) << " , Hard E_y: " << (towerEt * towerSinPhi) / (1 << (sin_bit_length_ - 1)) << "\n";
                softETxByBin[iBin] += (towerEt * towerCosPhi) / (1 << (sin_bit_length_ - 1));
                softETyByBin[iBin] += (towerEt * towerSinPhi) / (1 << (sin_bit_length_ - 1));
                jwojSoftSumEt += tower.etCode;
            }

            // Clear only what this event touched, ready for the next one.
            for (unsigned int iBlock = 0; iBlock < jwojTouchedBlocks.size(); iBlock++) {
                const unsigned int blockIdx = jwojTouchedBlocks[iBlock];
                jwojBlockEtCode[blockIdx] = 0u;
                jwojBlockEtGeV[blockIdx]  = 0.0;
                jwojBlockSeen[blockIdx]   = 0u;
                jwojBlockIsHard[blockIdx] = 0u;
            }

            // Term MET is the negative of the vector E_T sum, as everywhere else in this file.
            // The two terms are accumulated across eta bins uncoefficiented -- that is what
            // gets written out -- while the recombination applies the coefficients bin by bin.
            // With one bin this is exactly the scalar form the total MET above uses:
            //   lround(jetScaleFactor * hardMET + towerScaleFactor * softMET)
            int jwojHardMETx = 0, jwojHardMETy = 0;
            int jwojSoftMETx = 0, jwojSoftMETy = 0;
            double jwojMETxAcc = 0.0, jwojMETyAcc = 0.0;
            for (unsigned int iBin = 0; iBin < nJwoJBins; iBin++) {
                jwojHardMETx += -hardETxByBin[iBin];
                jwojHardMETy += -hardETyByBin[iBin];
                jwojSoftMETx += -softETxByBin[iBin];
                jwojSoftMETy += -softETyByBin[iBin];
                jwojMETxAcc  += jwojParams.hardCoeff[iBin] * static_cast<double>(-hardETxByBin[iBin])
                              + jwojParams.softCoeff[iBin] * static_cast<double>(-softETxByBin[iBin]);
                jwojMETyAcc  += jwojParams.hardCoeff[iBin] * static_cast<double>(-hardETyByBin[iBin])
                              + jwojParams.softCoeff[iBin] * static_cast<double>(-softETyByBin[iBin]);
            }
            int jwojMETx = static_cast<int>(std::lround(jwojMETxAcc));
            int jwojMETy = static_cast<int>(std::lround(jwojMETyAcc));

            unsigned int jwojHardMET = static_cast<unsigned int>(std::sqrt(static_cast<double>(jwojHardMETx) * jwojHardMETx + static_cast<double>(jwojHardMETy) * jwojHardMETy));
            unsigned int jwojSoftMET = static_cast<unsigned int>(std::sqrt(static_cast<double>(jwojSoftMETx) * jwojSoftMETx + static_cast<double>(jwojSoftMETy) * jwojSoftMETy));
            unsigned int jwojMET     = static_cast<unsigned int>(std::sqrt(static_cast<double>(jwojMETx)     * jwojMETx     + static_cast<double>(jwojMETy)     * jwojMETy));
            unsigned int jwojSumEt   = jwojHardSumEt + jwojSoftSumEt;

            // Packed through the same sign-magnitude / bitset path as the outputs above, so a
            // JwoJ value and a total MET value are the same kind of number on disk.
            std::bitset<signed_et_bit_length_> jwojMETx_bitset    (pack_signed_et(jwojMETx));
            std::bitset<signed_et_bit_length_> jwojMETy_bitset    (pack_signed_et(jwojMETy));
            std::bitset<signed_et_bit_length_> jwojHardMETx_bitset(pack_signed_et(jwojHardMETx));
            std::bitset<signed_et_bit_length_> jwojHardMETy_bitset(pack_signed_et(jwojHardMETy));
            std::bitset<signed_et_bit_length_> jwojSoftMETx_bitset(pack_signed_et(jwojSoftMETx));
            std::bitset<signed_et_bit_length_> jwojSoftMETy_bitset(pack_signed_et(jwojSoftMETy));
            std::bitset<et_bit_length_> jwojMET_bitset      (jwojMET       & maskN(et_bit_length_));
            std::bitset<et_bit_length_> jwojHardMET_bitset  (jwojHardMET   & maskN(et_bit_length_));
            std::bitset<et_bit_length_> jwojSoftMET_bitset  (jwojSoftMET   & maskN(et_bit_length_));
            std::bitset<et_bit_length_> jwojHardSumEt_bitset(jwojHardSumEt & maskN(et_bit_length_));
            std::bitset<et_bit_length_> jwojSoftSumEt_bitset(jwojSoftSumEt & maskN(et_bit_length_));
            std::bitset<et_bit_length_> jwojSumEt_bitset    (jwojSumEt     & maskN(et_bit_length_));

            out_jwojMETX      = undigitize_signed_et(jwojMETx_bitset);
            out_jwojMETY      = undigitize_signed_et(jwojMETy_bitset);
            out_jwojHardMETX  = undigitize_signed_et(jwojHardMETx_bitset);
            out_jwojHardMETY  = undigitize_signed_et(jwojHardMETy_bitset);
            out_jwojSoftMETX  = undigitize_signed_et(jwojSoftMETx_bitset);
            out_jwojSoftMETY  = undigitize_signed_et(jwojSoftMETy_bitset);
            out_jwojMET       = undigitize_et(jwojMET_bitset);
            out_jwojHardMET   = undigitize_et(jwojHardMET_bitset);
            out_jwojSoftMET   = undigitize_et(jwojSoftMET_bitset);
            out_jwojSumHardET = undigitize_et(jwojHardSumEt_bitset);
            out_jwojSumSoftET = undigitize_et(jwojSoftSumEt_bitset);
            out_jwojSumET     = undigitize_et(jwojSumEt_bitset);
        }

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
                  unsigned int pileup = 200,                 // Pileup scenario of the input sample; tags the output name (r16130 = PU200, r16129 = PU140)
                  bool useGEPJwoJ = false,                   // true = also write the GEP JwoJ MET branches (see the algorithm block in eventLoop)
                  double jwojHardEtThreshold = 10.0,         // Block E_T [GeV] above which a tower block joins the GEP JwoJ hard term
                  unsigned int jwojBlockSize = 1             // Side of the GEP JwoJ tower block, in towers: 1 (0.1x0.1) or 2 (0.2x0.2)
                  ) {

    if (signalBool) std::cout << "Processing signal: " << signalString << "\n";

    // GEP JwoJ coefficients are the EXISTING two: in JwoJ mode the towers below the hard
    // threshold are the soft term and the ones above it stand in for the jets, so
    // towerScaleFactor keeps its meaning as the soft weight and jetScaleFactor as the hard
    // weight. Nothing new goes in the filename for them -- the _twrSF / _jetSF tags already
    // carry them -- which is why only the hard-term threshold gets a tag of its own.
    const GEPJwoJParams jwojParams = makeUniformGEPJwoJParams(jwojHardEtThreshold, jetScaleFactor, towerScaleFactor, jwojBlockSize);
    // Read back from the struct rather than echoing the argument: an unsupported block size is
    // reset to 1 there, and the output filename is tagged from the same value, so the log and
    // the file name cannot disagree about what was actually run.
    if (useGEPJwoJ)
        std::cout << "GEP JwoJ MET enabled: block " << jwojParams.blockSize << "x" << jwojParams.blockSize
                  << " (" << jwojParams.blockSize * eta_granularity_ << " x "
                  << jwojParams.blockSize * eta_granularity_ << "), hard-term E_T threshold "
                  << jwojHardEtThreshold << " GeV, hard coefficient " << jetScaleFactor
                  << ", soft coefficient " << towerScaleFactor << "\n";

    auto infile  = explicitInputPath.empty() ? makeInputFileName(signalBool, signalString, "/data/larsonma/GEPHadronicEventReconstruction/ntuples/", pileup) : explicitInputPath;
    auto outfile = makeOutputMETFileName(maxTowersConsidered_, signalBool, signalString, useSKObjects, jetEtThreshold, towerEtThreshold, doJetTowerOverlapRemoval, met_output_dir_, useEtaSKObjects, towerScaleFactor, jetScaleFactor, pileup, useGEPJwoJ, jwojHardEtThreshold, jwojParams.blockSize);
    if (fileIndex >= 0) {
        size_t pos = outfile.rfind(".root");
        if (pos != std::string::npos)
            outfile = outfile.substr(0, pos) + "_file" + std::to_string(fileIndex) + ".root";
    }

    std::cout << "infile:  " << infile  << "\n";
    std::cout << "outfile: " << outfile << "\n";

    gSystem->RedirectOutput("debuglog_MET.log", "w");
    std::cout << "Calling event loop\n";
    eventLoop(infile, outfile, useSKObjects, jetEtThreshold, towerEtThreshold, doJetTowerOverlapRemoval, useEtaSKObjects, towerScaleFactor, jetScaleFactor, useGEPJwoJ, jwojParams);
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
const std::string met_baseline_total_output_suffix_ = "_METBaselineTotal.dat";
// EtaSK pileup-suppressed memory prints. The plain-SK prints these once pointed at no
// longer exist on eos -- only the unsuppressed (GEPCellsTowers, GEPConeJetsCellsTowers)
// and EtaSK collections are produced now, and EtaSK is what the rest of the MET chain
// runs on. The suffix is the subdirectory name lowercased, which is the convention
// HERNTupler writes them under.
const std::string met_baseline_tower_subdir_    = "GEPCellsTowersEtaSK/";
const std::string met_baseline_tower_suffix_    = "_gepcellstowersetask.dat";
const std::string met_baseline_jet_subdir_      = "GEPConeJetsCellsTowersEtaSK/";
const std::string met_baseline_jet_suffix_      = "_gepconejetscellstowersetask.dat";
const std::string met_baseline_log_prefix_      = "debuglog_METBaselineValidation_";
const std::string met_baseline_log_suffix_      = ".log";

// ROOT ntuple companion to the .dat vectors, written to its OWN directory on /data.
//
// Kept away from both of the other outputs on purpose. It does not belong beside the .dat
// vectors on eos, which are reference artifacts a testbench diffs against and which should
// not grow a ROOT file alongside them; and it does not belong in the emulation's ntuple
// directory either, since nothing there was produced at this fixed configuration -- no
// thresholds, no overlap removal, unit scale factors -- and mixing it in would make the
// filename tags there mean two different things.
const std::string met_baseline_ntuple_path_   = "/data/larsonma/GEPMET/METBaselineValidation/";
const std::string met_baseline_ntuple_suffix_ = "_METBaselineValidation.root";

// Output word, 64 bits, the MET TOB data format. MSB -> LSB:
//
//   [63]     SumET overflow      SumET exceeded the 13-bit output range
//   [62]     upstream overflow   an accepted input object had already saturated its E_T
//   [61]     MET overflow        |MET| exceeded the 13-bit output range
//   [60:48]  SumET               unsigned 13-bit scalar energy, saturated at 8191
//   [47:35]  Ey_miss             13-bit two's complement, saturated at +-4095
//   [34:22]  Ex_miss             13-bit two's complement, saturated at +-4095
//   [21:13]  phi_miss            9 bits, {phi6, 3'b000}: the 6-bit metPhiIndex in the top bits
//   [12:0]   ET_miss             unsigned 13-bit magnitude, saturated at 8191
//
// This REPLACES the Table 4.5 layout (reserved bit, Ex/Ey overflow flags at [58]/[59],
// 6-bit phi at [18:13]). Every field but ET_miss has moved -- regenerate
// rather than diff. Ex/Ey still saturate at +-4095; their overflow conditions are no longer
// in the word and survive only as diagnostics in the validation ntuple and log.
//
// Saturation, not masking. The magnitude fields clamp at their maxima, which is what
// metTobTwosComplementEt / metTobSumEt / metIntegerRoot do; pack_signed_et masks instead.
constexpr unsigned int met_tob_et_missing_width_        = et_bit_length_;        // 13
constexpr unsigned int met_tob_phi_missing_padding_     = 3;                     // zero LSBs below phi6
constexpr unsigned int met_tob_phi_missing_width_       = met_phi_bit_length_ + met_tob_phi_missing_padding_; // 9
constexpr unsigned int met_tob_ex_total_width_          = signed_et_bit_length_; // 13
constexpr unsigned int met_tob_ey_total_width_          = signed_et_bit_length_; // 13
constexpr unsigned int met_tob_sum_et_width_            = et_bit_length_;        // 13
constexpr unsigned int met_tob_flag_width_              = 1;

constexpr unsigned int met_tob_et_missing_low_          = 0;
constexpr unsigned int met_tob_phi_missing_low_         = met_tob_et_missing_low_        + met_tob_et_missing_width_;   // 13
constexpr unsigned int met_tob_ex_total_low_            = met_tob_phi_missing_low_       + met_tob_phi_missing_width_;  // 22
constexpr unsigned int met_tob_ey_total_low_            = met_tob_ex_total_low_          + met_tob_ex_total_width_;     // 35
constexpr unsigned int met_tob_sum_et_low_              = met_tob_ey_total_low_          + met_tob_ey_total_width_;     // 48
constexpr unsigned int met_tob_met_overflow_low_        = met_tob_sum_et_low_            + met_tob_sum_et_width_;       // 61
constexpr unsigned int met_tob_upstream_overflow_low_   = met_tob_met_overflow_low_      + met_tob_flag_width_;         // 62
constexpr unsigned int met_tob_sum_et_overflow_low_     = met_tob_upstream_overflow_low_ + met_tob_flag_width_;         // 63
constexpr unsigned int met_tob_total_bits_              = met_tob_sum_et_overflow_low_   + met_tob_flag_width_;
static_assert(met_tob_total_bits_ == 64, "MET TOB output word must be exactly 64 bits");
static_assert(met_tob_et_missing_low_ == 0 && met_tob_phi_missing_low_ == 13 &&
              met_tob_ex_total_low_   == 22 && met_tob_ey_total_low_   == 35 &&
              met_tob_sum_et_low_     == 48 && met_tob_met_overflow_low_ == 61 &&
              met_tob_upstream_overflow_low_ == 62 && met_tob_sum_et_overflow_low_ == 63,
              "MET TOB field offsets must match the TOB format");

// One object as it appears in a memory-print line: already-digitized codes.
struct MemPrintObject {
    unsigned int phi = 0;
    unsigned int eta = 0;
    unsigned int et  = 0;
};

// Jet memory prints carry jet TOB codes when input jets do; decode them to the tower grid sinLUT_ indexes
inline std::vector<std::vector<MemPrintObject>> jetMemPrintsToTowerGrid(std::vector<std::vector<MemPrintObject>> events) {
    if (!input_jets_use_jet_tob_encoding_) return events;
    for (auto& event : events)
        for (auto& jet : event) {
            jet.eta = jet_eta_to_tower_index(jet.eta);
            jet.phi = jet_phi_to_tower_index(jet.phi);
        }
    return events;
}

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

// One MET term, fully digitized: everything the TOB word needs and nothing that depends
// on where the inputs came from. This is the boundary the GlobalSim port will use -- the
// float front end (digitize()/digitize_phi() on ROOT ntuple values) sits on one side of
// it, already-digitized inputs like the memory prints on the other, and the arithmetic
// below is shared by both rather than written out twice.
struct MetTermDigitized {
    int                metX  = 0; // full-width accumulator counts, MET = -(sum E_T)
    int                metY  = 0;
    unsigned int       met   = 0; // 13-bit magnitude from the LUT root, saturated
    unsigned int       phi   = 0; // 6-bit azimuth, MET_PHI_COMPARATOR convention
    unsigned long long sumEt = 0; // full-width scalar sum, saturated only at packing
    bool metOverflow    = false;  // the magnitude did not fit 13 bits
    bool inputSaturated = false;  // an accepted input object had already clipped its E_T
};

// Accumulate one collection of already-digitized objects into a MET term.
// Identical arithmetic to the tower and jet loops of eventLoop: the sine lookup is
// indexed by the phi code, cosine is the same table read half_pi_digitized_in_phi_
// codes along, and the product is brought back to E_T counts by an integer divide by
// 1 << (sin_bit_length_ - 1). MET is the negative of the vector E_T sum.
inline MetTermDigitized accumulateMemPrintMet(const std::vector<MemPrintObject>& objects,
                                              unsigned int maxObjects) {
    MetTermDigitized term;

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

        term.sumEt += object.et;
        // TOB bit [62]: this object arrived already clipped to the top of its E_T field,
        // so every sum built from it is a lower bound on the true energy.
        if (metInputEtSaturated(object.et)) term.inputSaturated = true;
    }

    term.metX = -etXSum;
    term.metY = -etYSum;
    term.met  = metIntegerRoot(term.metX, term.metY, term.metOverflow);
    term.phi  = metPhiIndex(term.metX, term.metY);
    return term;
}

// Total MET: the tower and jet terms combined.
//
// Combined at FULL WIDTH, before either term has been squeezed into a TOB field. The
// magnitude and the azimuth are then recomputed from the summed components rather than
// derived from the two input magnitudes, because MET is a vector sum -- |a + b| is not a
// function of |a| and |b|, and the phi ladder needs the summed components anyway.
//
// No scale factors. The baseline validation path is fixed at unit tower and jet weights
// by design (see the header of this section), so the total is the plain component-wise
// sum. The ROOT event loop, which does apply twrSF / jetSF, is where a weighted total
// lives.
//
// Flags:
//  - inputSaturated (TOB [62]) is inherited as an OR, which is what MET_Engine.v
//    specifies for that bit: "OR of accepted JET1 and TowerSum E_T overflow
//    indications". If either collection contained an object already clipped at the top
//    of its E_T field, every sum built from it is a lower bound, so the total has to say
//    so even when the other collection was clean.
//  - metOverflow is NOT inherited. It is recomputed here by metIntegerRoot on the summed
//    components: two terms that each fit can sum to one that does not, and two that each
//    overflow can cancel into one that fits.
//  - The Ex / Ey / SumET overflow flags are not carried on the struct at all. They are
//    properties of the output field and are raised by the packer, on these summed values.
inline MetTermDigitized combineMetTerms(const MetTermDigitized& towerTerm,
                                        const MetTermDigitized& jetTerm) {
    MetTermDigitized total;
    total.metX  = towerTerm.metX  + jetTerm.metX;
    total.metY  = towerTerm.metY  + jetTerm.metY;
    total.sumEt = towerTerm.sumEt + jetTerm.sumEt;
    total.inputSaturated = towerTerm.inputSaturated || jetTerm.inputSaturated;
    total.met = metIntegerRoot(total.metX, total.metY, total.metOverflow);
    total.phi = metPhiIndex(total.metX, total.metY);
    return total;
}

// Write one memory-print line: index, the eight binary fields of the TOB word in
// MSB -> LSB order, and the packed hex word.
//
// The flags are computed HERE, at packing time, rather than carried in from the
// arithmetic: overflow is a property of the output field, so the term keeps full-width
// values and only learns they did not fit at the moment they are squeezed into the TOB.
// The one exception is metOverflow, which MET_SUM_SQUARE raises inside the magnitude
// path, and inputSaturated, which is a fact about the inputs the packer cannot see.
inline void writeMetBaselineLine(std::ofstream& outFile, unsigned int index,
                                 const MetTermDigitized& term) {
    bool exOverflow = false, eyOverflow = false, sumEtOverflow = false; // Ex/Ey flags not in the word
    const uint64_t exPacked      = metTobTwosComplementEt(term.metX, exOverflow);
    const uint64_t eyPacked      = metTobTwosComplementEt(term.metY, eyOverflow);
    const uint64_t sumEtPacked   = metTobSumEt(term.sumEt, sumEtOverflow);
    const uint64_t etMissingBits = static_cast<uint64_t>(term.met & maskN(met_tob_et_missing_width_));
    const uint64_t phiMissBits   = static_cast<uint64_t>(term.phi & maskN(met_phi_bit_length_)) << met_tob_phi_missing_padding_;

    const uint64_t sumEtOverflowBit     = sumEtOverflow       ? 1ull : 0ull;
    const uint64_t upstreamOverflowBit  = term.inputSaturated ? 1ull : 0ull;
    const uint64_t metOverflowBit       = term.metOverflow    ? 1ull : 0ull;

    std::bitset<met_tob_flag_width_>        sum_et_overflow_bitset(sumEtOverflowBit);
    std::bitset<met_tob_flag_width_>        upstream_overflow_bitset(upstreamOverflowBit);
    std::bitset<met_tob_flag_width_>        met_overflow_bitset(metOverflowBit);
    std::bitset<met_tob_sum_et_width_>      sum_et_bitset(sumEtPacked);
    std::bitset<met_tob_ey_total_width_>    ey_total_bitset(eyPacked);
    std::bitset<met_tob_ex_total_width_>    ex_total_bitset(exPacked);
    std::bitset<met_tob_phi_missing_width_> phi_missing_bitset(phiMissBits);
    std::bitset<met_tob_et_missing_width_>  et_missing_bitset(etMissingBits);

    const uint64_t combined_value =
        (sumEtOverflowBit    << met_tob_sum_et_overflow_low_)   |
        (upstreamOverflowBit << met_tob_upstream_overflow_low_) |
        (metOverflowBit      << met_tob_met_overflow_low_)      |
        (sumEtPacked         << met_tob_sum_et_low_)            |
        (eyPacked            << met_tob_ey_total_low_)          |
        (exPacked            << met_tob_ex_total_low_)          |
        (phiMissBits         << met_tob_phi_missing_low_)       |
        (etMissingBits       << met_tob_et_missing_low_);

    std::stringstream hex_stream;
    hex_stream << std::hex << std::nouppercase << std::setfill('0') << std::setw(met_tob_total_bits_ / 4) << combined_value;

    outFile << "0x" << std::hex << std::setw(2) << std::setfill('0') << index << std::dec << std::setfill(' ')
            << " " << sum_et_overflow_bitset.to_string()
            << "|" << upstream_overflow_bitset.to_string()
            << "|" << met_overflow_bitset.to_string()
            << "|" << sum_et_bitset.to_string()
            << "|" << ey_total_bitset.to_string()
            << "|" << ex_total_bitset.to_string()
            << "|" << phi_missing_bitset.to_string()
            << "|" << et_missing_bitset.to_string()
            << " 0x" << hex_stream.str() << std::endl;
}

// One MET term's worth of branches in the validation ntuple.
//
// The ntuple is the plottable companion to the .dat vectors: same events, same arithmetic,
// same TOB encoding. It carries BOTH representations of every quantity --
//
//   * the digitized value, which is literally the field the TOB word carries, so a tree
//     entry can be read straight against the corresponding .dat line, and
//   * the GeV / radian equivalent, so the same tree can be histogrammed without the reader
//     having to know the LSB
//
// -- plus the five overflow flags, which are otherwise only visible by decoding the .dat.
//
// There is no risk of the two outputs drifting apart: fill() and writeMetBaselineLine()
// below are two call sites of the SAME encoders (metTobTwosComplementEt / metTobSumEt), not two
// implementations of the encoding.
struct MetBaselineNtupleTerm {
    // Digitized, exactly as the TOB carries them. Ex/Ey are the two's complement fields
    // decoded back to a signed count, i.e. already saturated.
    int          ex       = 0;
    int          ey       = 0;
    unsigned int met      = 0;
    unsigned int phiIndex = 0;
    unsigned int sumEt    = 0;
    // The same values in physical units.
    double exGeV    = 0.0;
    double eyGeV    = 0.0;
    double metGeV   = 0.0;
    double sumEtGeV = 0.0;
    double phiRad   = 0.0;   // bin center, folded into (-pi, pi]
    // Overflow flags; Ex/Ey ones are diagnostics only, the TOB carries SumET [63], upstream [62], MET [61].
    bool exOverflow       = false;
    bool eyOverflow       = false;
    bool sumEtOverflow    = false;
    bool metOverflow      = false;
    bool upstreamOverflow = false;

    // Branch names are <prefix>_<quantity>, e.g. Tower_Ex, Total_METGeV.
    void book(TTree* tree, const std::string& prefix) {
        tree->Branch((prefix + "_Ex").c_str(),       &ex);
        tree->Branch((prefix + "_Ey").c_str(),       &ey);
        tree->Branch((prefix + "_MET").c_str(),      &met);
        tree->Branch((prefix + "_PhiIndex").c_str(), &phiIndex);
        tree->Branch((prefix + "_SumET").c_str(),    &sumEt);

        tree->Branch((prefix + "_ExGeV").c_str(),    &exGeV);
        tree->Branch((prefix + "_EyGeV").c_str(),    &eyGeV);
        tree->Branch((prefix + "_METGeV").c_str(),   &metGeV);
        tree->Branch((prefix + "_SumETGeV").c_str(), &sumEtGeV);
        tree->Branch((prefix + "_PhiRad").c_str(),   &phiRad);

        tree->Branch((prefix + "_ExOverflow").c_str(),       &exOverflow);
        tree->Branch((prefix + "_EyOverflow").c_str(),       &eyOverflow);
        tree->Branch((prefix + "_SumETOverflow").c_str(),    &sumEtOverflow);
        tree->Branch((prefix + "_METOverflow").c_str(),      &metOverflow);
        tree->Branch((prefix + "_UpstreamOverflow").c_str(), &upstreamOverflow);
    }

    void fill(const MetTermDigitized& term) {
        // Through the TOB encoders, not straight off the accumulator: the .dat carries the
        // saturated field, so the ntuple has to as well or the two would disagree on
        // exactly the events where saturation bites.
        const unsigned int exPacked    = metTobTwosComplementEt(term.metX, exOverflow);
        const unsigned int eyPacked    = metTobTwosComplementEt(term.metY, eyOverflow);
        const unsigned int sumEtPacked = metTobSumEt(term.sumEt, sumEtOverflow);

        // Decode the two's complement field back to a signed count for the integer branch.
        ex = metTobTwosComplementToInt(exPacked);
        ey = metTobTwosComplementToInt(eyPacked);

        met      = term.met;
        phiIndex = term.phi; // 6-bit comparator index; the TOB field is phiIndex << met_tob_phi_missing_padding_
        sumEt    = sumEtPacked;

        exGeV    = ex    * et_granularity_;
        eyGeV    = ey    * et_granularity_;
        metGeV   = met   * et_granularity_;
        sumEtGeV = sumEt * et_granularity_;
        phiRad   = undigitize_met_phi(term.phi);

        metOverflow      = term.metOverflow;
        upstreamOverflow = term.inputSaturated;
    }
};

// MET baseline validation entry point.
//
// Use:
//   root -b -l
//   root [0] .L metEmulation.cc+
//   root [1] metBaselineValidation(true, "ggF_hh_bbbb")
//
// Writes three files per sample to met_baseline_output_subdir_ -- tower MET, jet MET and
// total MET -- each with a single 0x00 line per event. Kept apart rather than as three
// lines of one file so each converts on its own to whatever the consumer wants without
// having to demultiplex by object index.
//
// Total MET is a genuine third vector, not a value a consumer could reconstruct from the
// other two words: the tower and jet words carry Ex/Ey already SATURATED into 12-bit
// magnitudes, whereas the total is summed at full width before packing. It also has its
// own magnitude (the LUT root of the summed components, not a function of the two input
// magnitudes), its own azimuth, and its own overflow flags.
//
// Per-event printouts go to met_baseline_log_prefix_ + sample + met_baseline_log_suffix_
// in the working directory; only the paths and the closing summary reach the terminal.
void metBaselineValidation(bool signalBool = true,          // true = signal sample, false = dijet background
                           std::string signalString = "ggF_hh_bbbb", // Which signal sample (memory prints exist for ggF_hh_bbbb)
                           unsigned int jzSlice = 3,        // JZ slice, used only when signalBool is false
                           int maxEvents = 1000) {          // Process at most this many events; < 0 processes all

    std::string baseName = makeMemPrintBaseName(signalBool, signalString, jzSlice);
    if (baseName.empty()) return;

    std::string towerFile       = met_baseline_mem_prints_path_ + met_baseline_tower_subdir_ + baseName + met_baseline_tower_suffix_;
    std::string jetFile         = met_baseline_mem_prints_path_ + met_baseline_jet_subdir_   + baseName + met_baseline_jet_suffix_;
    std::string outputPath      = met_baseline_mem_prints_path_ + met_baseline_output_subdir_;
    std::string towerOutputFile = outputPath + baseName + met_baseline_tower_output_suffix_;
    std::string jetOutputFile   = outputPath + baseName + met_baseline_jet_output_suffix_;
    std::string totalOutputFile = outputPath + baseName + met_baseline_total_output_suffix_;
    // Log stays in the working directory rather than next to the vectors on eos: it is a
    // run artifact, not a reference artifact, and it is rewritten on every run.
    std::string logFile         = met_baseline_log_prefix_ + baseName + met_baseline_log_suffix_;
    std::string ntupleFile      = met_baseline_ntuple_path_ + baseName + met_baseline_ntuple_suffix_;

    std::cout << "tower memory print: " << towerFile       << "\n";
    std::cout << "jet memory print:   " << jetFile         << "\n";
    std::cout << "tower output:       " << towerOutputFile << "\n";
    std::cout << "jet output:         " << jetOutputFile   << "\n";
    std::cout << "total output:       " << totalOutputFile << "\n";
    std::cout << "ntuple output:      " << ntupleFile      << "\n";

    std::vector<std::vector<MemPrintObject>> towerEvents = readMemPrintFile(towerFile);
    std::vector<std::vector<MemPrintObject>> jetEvents   = jetMemPrintsToTowerGrid(readMemPrintFile(jetFile));

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
    std::ofstream totalOutFile(totalOutputFile);
    if (!totalOutFile.is_open()) {
        std::cerr << "Error: Could not open file " << totalOutputFile << std::endl;
        return;
    }

    // ROOT companion to the three .dat files. Opened alongside them, before the redirect
    // below, so a bad path is reported on the terminal like the others.
    gSystem->mkdir(met_baseline_ntuple_path_.c_str(), true);
    TFile* ntupleOutFile = TFile::Open(ntupleFile.c_str(), "RECREATE");
    if (!ntupleOutFile || ntupleOutFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << ntupleFile << std::endl;
        return;
    }
    ntupleOutFile->cd();
    TTree* metBaselineTree = new TTree("metBaselineTree",
                                       "Per-event MET baseline validation terms, digitized and in GeV");
    unsigned int out_eventIndex = 0;
    // The memory prints carry no run/event number, so the index into the print is the only
    // handle on an event -- and it is the same index the .dat "Event : N" headers use, so
    // an entry here and a block there line up by construction.
    metBaselineTree->Branch("eventIndex", &out_eventIndex);

    MetBaselineNtupleTerm towerNtuple, jetNtuple, totalNtuple;
    towerNtuple.book(metBaselineTree, "Tower");
    jetNtuple.book(metBaselineTree,   "Jet");
    totalNtuple.book(metBaselineTree, "Total");

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
    std::cout << "total output:       " << totalOutputFile << "\n";
    std::cout << "events:             " << std::dec << eventsToProcess << "\n";

    for (unsigned int iEvt = 0; iEvt < eventsToProcess; iEvt++) {
        towerOutFile << "Event : " << std::dec << iEvt << std::endl;
        jetOutFile   << "Event : " << std::dec << iEvt << std::endl;
        totalOutFile << "Event : " << std::dec << iEvt << std::endl;

        const MetTermDigitized towerTerm = accumulateMemPrintMet(towerEvents.at(iEvt), maxTowersConsidered_);
        const MetTermDigitized jetTerm   = accumulateMemPrintMet(jetEvents.at(iEvt),   maxJetsConsidered_);
        const MetTermDigitized totalTerm = combineMetTerms(towerTerm, jetTerm);

        // Echo exactly the quantities that go into the output words, in digitized counts
        // with the GeV equivalent alongside, so a line of the .dat can be read back
        // without decoding it by hand. phi is printed as its bin and as the bin centre in
        // radians; note that bin is measured from +x, NOT on the tower phi grid.
        auto printTerm = [](const char* label, const MetTermDigitized& term) {
            std::cout << "  " << label
                      << " Ex: "  << term.metX  << " (" << term.metX * et_granularity_ << " GeV)"
                      << "  Ey: " << term.metY  << " (" << term.metY * et_granularity_ << " GeV)"
                      << "  MET: "<< term.met   << " (" << term.met  * et_granularity_ << " GeV)"
                      << "  phi: "<< term.phi   << " (" << undigitize_met_phi(term.phi) << " rad)"
                      << "  SumET: " << term.sumEt << " (" << term.sumEt * et_granularity_ << " GeV)"
                      << (term.metOverflow    ? "  [MET overflow]"      : "")
                      << (term.inputSaturated ? "  [upstream overflow]" : "")
                      << "\n";
        };
        std::cout << "iEvt: " << std::dec << iEvt
                  << "  towers: " << towerEvents.at(iEvt).size()
                  << "  jets: "   << jetEvents.at(iEvt).size() << "\n";
        printTerm("tower ", towerTerm);
        printTerm("jet   ", jetTerm);
        printTerm("total ", totalTerm);

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

        writeMetBaselineLine(towerOutFile, 0, towerTerm);
        writeMetBaselineLine(jetOutFile,   0, jetTerm);
        writeMetBaselineLine(totalOutFile, 0, totalTerm);

        // Same three terms into the ntuple, from the same structs the .dat lines were
        // written from, so the two outputs cannot disagree about an event.
        out_eventIndex = iEvt;
        towerNtuple.fill(towerTerm);
        jetNtuple.fill(jetTerm);
        totalNtuple.fill(totalTerm);
        metBaselineTree->Fill();
    }

    // Restore stdout/stderr before the summary, so the summary lands on the terminal.
    gSystem->RedirectOutput(nullptr, "", &redirectHandle);

    towerOutFile.close();
    jetOutFile.close();
    totalOutFile.close();

    // Written after the redirect is restored, so the ROOT bookkeeping messages and the
    // summary below both land on the terminal rather than in the run log.
    ntupleOutFile->cd();
    metBaselineTree->Write("", TObject::kOverwrite);
    ntupleOutFile->Close();

    std::cout << "Wrote " << std::dec << eventsToProcess << " events to " << towerOutputFile << "\n";
    std::cout << "Wrote " << std::dec << eventsToProcess << " events to " << jetOutputFile   << "\n";
    std::cout << "Wrote " << std::dec << eventsToProcess << " events to " << totalOutputFile << "\n";
    std::cout << "Wrote " << std::dec << eventsToProcess << " events to " << ntupleFile      << "\n";
    std::cout << "Per-event printouts in " << logFile << "\n";
}


// =====================================================================================
// ===                 MET edge-case / random validation path                        ===
// =====================================================================================
// The baseline validation path above, but with made-up inputs instead of HERNTupler
// memory prints. A fixed block of edge cases comes first, then a seeded block of
// pseudo-random events. BOTH the inputs and the outputs are written as memory prints to
// one directory of their own, so a firmware testbench can be fed the input files and
// diffed against the output files directly.
//
// The arithmetic is not re-implemented here. The generated inputs are written to disk,
// read back with readMemPrintFileWithOverflow, and run through the same accumulateMemPrintMet ->
// combineMetTerms -> writeMetBaselineLine chain as the baseline path. The outputs are
// therefore computed from exactly what is in the input files, rather than from in-memory
// values that might not have survived the round trip through the text format.
//
// Same fixed configuration as the baseline: no thresholds, no overlap removal, unit
// scale factors. The fixed block never exceeds the object limits (maxTowersConsidered_ /
// maxJetsConsidered_); random events can carry more towers than maxTowersConsidered_,
// and those are truncated to the first maxTowersConsidered_ (see the random-block notes).
//
// Fixed block, in file order:
//   - all zero:           no objects at all; objects present at E_T = 0
//   - upstream overflow:  an E_T = 8191 tower, jet, and one of each back-to-back
//   - phi sweep:          one tower and one jet per phi code, all 64 codes
//   - Ex/Ey boundaries:   each component landing exactly on +-4095 (fits) and +-4096
//                         (flags)
//   - SumET:              exactly 8191 (fits), exactly 8192 (flags), and a large balanced
//                         sum with MET = 0
//   - MET magnitude:      Ex/Ey overflowing with MET still fitting; MET overflowing via
//                         the sum of squares; MET overflowing via a component >= 8192
//   - cancellation:       tower and jet terms that each overflow but cancel in the total,
//                         to zero and to a nonzero in-range value
//
// Some flags cannot be raised on their own, and no case pretends otherwise. A single
// object at E_T = 8191 cannot fit both components in +-4095, so the upstream cases also
// raise Ex/Ey overflow. MET can only overflow once a component has passed 4095, so every
// MET overflow case also raises Ex or Ey. The flags each case ACTUALLY raises are listed
// per event in the _Cases.txt file.
// =====================================================================================

// Set to true to print every input object into the run log as well as the per-event terms.
constexpr bool met_random_debug_ = false;

// Input memory-print object word, 33 bits, MSB -> LSB: overflow | phi | eta | et.
//
// The low 32 bits are the HERNTupler tower / jet print word, at the field widths
// HERNTupler writes them with (not the MET-internal phi / eta widths of constants.h --
// the codes are small enough for either). Bit [32] is the per-object overflow tag the
// firmware ORs into TOB bit [62]; HERNTupler prints do not carry it, so these files have
// one field more than the baseline-path inputs and are read by readMemPrintFileWithOverflow,
// not readMemPrintFile.
//
// The tag is set exactly when the object's E_T is saturated (metInputEtSaturated), which
// is the condition the emulator raises bit [62] on. Keeping the two tied means the tag in
// the file and the flag in the output always agree.
constexpr unsigned int met_mem_print_overflow_width_ = 1;
constexpr unsigned int met_mem_print_phi_width_      = 9;
constexpr unsigned int met_mem_print_eta_width_      = 10;
constexpr unsigned int met_mem_print_et_width_       = et_bit_length_; // 13
constexpr unsigned int met_mem_print_et_low_         = 0;
constexpr unsigned int met_mem_print_eta_low_        = met_mem_print_et_low_  + met_mem_print_et_width_;  // 13
constexpr unsigned int met_mem_print_phi_low_        = met_mem_print_eta_low_ + met_mem_print_eta_width_; // 23
constexpr unsigned int met_mem_print_overflow_low_   = met_mem_print_phi_low_ + met_mem_print_phi_width_; // 32
constexpr unsigned int met_mem_print_total_bits_     = met_mem_print_overflow_low_ + met_mem_print_overflow_width_;
static_assert(met_mem_print_overflow_low_ == 32 && met_mem_print_total_bits_ == 33,
              "memory-print object word must be the 32-bit HERNTupler word plus an overflow tag at [32]");

const std::string met_random_output_subdir_        = "METRandomValidation/";
const std::string met_random_base_prefix_          = "METRandomValidation_seed";
const std::string met_random_tower_input_suffix_   = "_TowerInput.dat";
const std::string met_random_jet_input_suffix_     = "_JetInput.dat";
const std::string met_random_cases_suffix_         = "_Cases.txt";
const std::string met_random_ntuple_path_          = "/data/larsonma/GEPMET/METRandomValidation/";
const std::string met_random_ntuple_suffix_        = ".root";
const std::string met_random_log_prefix_           = "debuglog_";
const std::string met_random_log_suffix_           = ".log";

// Saturated and largest unsaturated input E_T codes.
const unsigned int met_random_et_saturated_  = maskN(et_bit_length_);      // 8191
const unsigned int met_random_et_max_clean_  = maskN(et_bit_length_) - 1u; // 8190

// eta plays no part in the MET arithmetic; fixed cases just put everything mid-range.
constexpr unsigned int met_random_fixed_eta_ = eta_range_ / 2;

// Random block. Up to 6400 towers per event -- the total tower count -- with E_T
// exponential at a per-event mean drawn log-uniformly, so the block spans quiet events,
// events near the SumET boundary, and events well past it. Jets are mostly a softer
// exponential, with a fraction drawn flat across the whole E_T range to produce
// Ex/Ey/MET overflows, and a fraction forced to saturation.
//
// 6400 is above maxTowersConsidered_ (4096, = MAX_TOWERS in MET_Engine.v), so events
// with more than 4096 towers are truncated to the first 4096 by accumulateMemPrintMet.
// The input files still carry every tower.
constexpr unsigned int met_random_max_towers_        = 6400;
constexpr double       met_random_tower_mean_et_min_ = 1.0;    // counts, 0.25 GeV
constexpr double       met_random_tower_mean_et_max_ = 16.0;   // counts, 4 GeV
constexpr double       met_random_tower_saturate_prob_ = 1e-4; // per tower
constexpr double       met_random_jet_mean_et_       = 200.0;  // counts, 50 GeV
constexpr double       met_random_jet_hard_prob_     = 0.10;   // per jet, flat in [0, 8190]
constexpr double       met_random_jet_saturate_prob_ = 0.03;   // per jet, E_T = 8191

// One generated event: its label, what it is testing, and its two input collections.
struct MetValidationEvent {
    std::string name;
    std::string description;
    std::vector<MemPrintObject> towers;
    std::vector<MemPrintObject> jets;
};

inline MemPrintObject makeMemPrintObject(unsigned int phi, unsigned int eta, unsigned int et) {
    MemPrintObject object;
    object.phi = phi;
    object.eta = eta;
    object.et  = et;
    return object;
}

// Write one input memory-print line: index, binary overflow|phi|eta|et, packed hex word.
// Same layout as the HERNTupler prints with the overflow tag as an extra leading field;
// the hex word is 9 digits, the 33 bits rounded up to whole nibbles.
inline void writeMemPrintObjectLine(std::ofstream& outFile, unsigned int index, const MemPrintObject& object) {
    const uint64_t overflowBit = metInputEtSaturated(object.et) ? 1ull : 0ull;

    std::bitset<met_mem_print_overflow_width_> overflow_bitset(overflowBit);
    std::bitset<met_mem_print_phi_width_>      phi_bitset(object.phi);
    std::bitset<met_mem_print_eta_width_>      eta_bitset(object.eta);
    std::bitset<met_mem_print_et_width_>       et_bitset(object.et);

    const uint64_t combined_value =
        (overflowBit                       << met_mem_print_overflow_low_) |
        (static_cast<uint64_t>(object.phi) << met_mem_print_phi_low_)      |
        (static_cast<uint64_t>(object.eta) << met_mem_print_eta_low_)      |
        (static_cast<uint64_t>(object.et)  << met_mem_print_et_low_);

    std::stringstream hex_stream;
    hex_stream << std::hex << std::nouppercase << std::setfill('0')
               << std::setw((met_mem_print_total_bits_ + 3) / 4) << combined_value;

    outFile << "0x" << std::hex << std::setw(2) << std::setfill('0') << index << std::dec << std::setfill(' ')
            << " " << overflow_bitset.to_string()
            << "|" << phi_bitset.to_string()
            << "|" << eta_bitset.to_string()
            << "|" << et_bitset.to_string()
            << " 0x" << hex_stream.str() << std::endl;
}

// Reader for the files writeMemPrintObjectLine produces: readMemPrintFile, but for the
// four-field overflow|phi|eta|et binary column. Kept separate rather than folded into
// readMemPrintFile so the baseline path's parsing of the HERNTupler prints is untouched.
//
// The overflow tag is not stored -- MemPrintObject has no field for it, and the emulator
// raises TOB bit [62] from the E_T itself. Instead every line is checked: a tag that
// disagrees with the E_T means the file was not written by this path, and is reported.
inline std::vector<std::vector<MemPrintObject>> readMemPrintFileWithOverflow(const std::string& fileName) {
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

        std::vector<std::string> fields;
        std::stringstream binStream(bin);
        std::string field;
        while (std::getline(binStream, field, '|')) fields.push_back(field);
        if (fields.size() != 4) {
            std::cerr << "Error: Malformed memory-print line (expected overflow|phi|eta|et) -> " << line << std::endl;
            continue;
        }

        const unsigned int overflow = static_cast<unsigned int>(std::stoul(fields.at(0), nullptr, 2));
        MemPrintObject object;
        object.phi = static_cast<unsigned int>(std::stoul(fields.at(1), nullptr, 2));
        object.eta = static_cast<unsigned int>(std::stoul(fields.at(2), nullptr, 2));
        object.et  = static_cast<unsigned int>(std::stoul(fields.at(3), nullptr, 2));
        if ((overflow != 0u) != metInputEtSaturated(object.et))
            std::cerr << "Warning: overflow tag " << overflow << " disagrees with E_T " << object.et
                      << " -> " << line << std::endl;
        eventObjects.back().push_back(object);
    }

    inFile.close();
    return eventObjects;
}

// One MET component of a lone object, through the real accumulator rather than a copy of
// its arithmetic, so the boundary search below cannot disagree with what gets written.
inline int singleObjectMetComponent(unsigned int phi, unsigned int et, bool xComponent) {
    const MetTermDigitized term = accumulateMemPrintMet({makeMemPrintObject(phi, met_random_fixed_eta_, et)}, 1);
    return xComponent ? term.metX : term.metY;
}

// Smallest unsaturated E_T at which a lone object at this phi code gives exactly this MET
// component, or -1 if there is none. The per-object map E_T -> component has slope
// |LUT| / 4096 < 1, so every value up to ~8180 in magnitude is reachable.
inline int findEtForMetComponent(unsigned int phi, bool xComponent, int target) {
    for (unsigned int et = 0; et <= met_random_et_max_clean_; et++)
        if (singleObjectMetComponent(phi, et, xComponent) == target) return static_cast<int>(et);
    return -1;
}

// Towers only, one MET component landing EXACTLY on target.
//
// No phi code sits on an axis -- code k is the centre of tower k, so code 0 is pi/64
// -- so a lone object would also carry ~5% of its E_T in the other component. Instead the
// event is a pair straddling the axis: mirror images in the other component, which then
// cancels to zero or a count or two. MET is minus the vector sum, so a positive MET
// component needs the objects on the negative side of the axis.
inline MetValidationEvent makeComponentBoundaryEvent(bool xComponent, int target) {
    unsigned int phiA = 0, phiB = 0;
    if (xComponent) {
        // +Ex: objects straddling phi = pi; -Ex: straddling phi = 0.
        phiA = (target > 0) ? pi_digitized_in_phi_ - 1 : phi_range_ - 1;
        phiB = (target > 0) ? pi_digitized_in_phi_     : 0;
    } else {
        // +Ey: objects straddling phi = -pi/2; -Ey: straddling phi = +pi/2.
        phiA = (target > 0) ? pi_digitized_in_phi_ + half_pi_digitized_in_phi_ - 1 : half_pi_digitized_in_phi_ - 1;
        phiB = (target > 0) ? pi_digitized_in_phi_ + half_pi_digitized_in_phi_     : half_pi_digitized_in_phi_;
    }

    const int half = target / 2;
    const int etA  = findEtForMetComponent(phiA, xComponent, half);
    const int etB  = findEtForMetComponent(phiB, xComponent, target - half);

    MetValidationEvent event;
    const std::string axis = xComponent ? "Ex" : "Ey";
    event.name = axis + "_" + (target > 0 ? "+" : "") + std::to_string(target);
    event.description = "towers only, " + axis + " = " + std::to_string(target) + " exactly; " +
                        (std::abs(target) > static_cast<int>(maskN(met_tob_signed_mag_bit_length_))
                             ? "just outside +-4095, " + axis + " overflow"
                             : "largest magnitude that fits, no " + axis + " overflow");

    if (etA < 0 || etB < 0) {
        std::cerr << "Error: no E_T reaches " << axis << " = " << target << " for case " << event.name << std::endl;
        return event;
    }
    event.towers = {makeMemPrintObject(phiA, met_random_fixed_eta_, etA),
                    makeMemPrintObject(phiB, met_random_fixed_eta_, etB)};

    const MetTermDigitized check = accumulateMemPrintMet(event.towers, maxTowersConsidered_);
    const int reached = xComponent ? check.metX : check.metY;
    if (reached != target)
        std::cerr << "Error: case " << event.name << " reached " << axis << " = " << reached
                  << " instead of " << target << std::endl;
    return event;
}

// The fixed edge-case block. See the header of this section for the list.
inline std::vector<MetValidationEvent> makeMetEdgeCaseEvents() {
    std::vector<MetValidationEvent> events;
    const unsigned int eta = met_random_fixed_eta_;

    // Phi codes used below. Code k and k + pi_digitized_in_phi_ are exactly back-to-back,
    // so equal-E_T pairs on them cancel in both components.
    const unsigned int phiDiag     = pi_digitized_in_phi_ + half_pi_digitized_in_phi_ / 2; // 40, ~ -3pi/4
    const unsigned int phiDiagOpp  = half_pi_digitized_in_phi_ / 2;                        //  8, ~ +pi/4
    const unsigned int phiZeroLow  = phi_range_ - 1;                                       // 63, -pi/64
    const unsigned int phiZeroHigh = 0;                                                    //  0, +pi/64
    const unsigned int phiPiLow    = pi_digitized_in_phi_ - 1;                             // 31, pi - pi/64
    const unsigned int phiPiHigh   = pi_digitized_in_phi_;                                 // 32, -pi + pi/64
    // Codes straddling +pi/4. Objects here put MET near -3pi/4 with |Ex| ~= |Ey|.
    const unsigned int phiQuarterLow  = half_pi_digitized_in_phi_ / 2 - 1; //  7
    const unsigned int phiQuarterHigh = half_pi_digitized_in_phi_ / 2;     //  8

    // --- all zero ---
    {
        MetValidationEvent event;
        event.name = "zero_empty";
        event.description = "no towers, no jets: every field and flag zero";
        events.push_back(event);
    }
    {
        MetValidationEvent event;
        event.name = "zero_objects";
        event.description = "one E_T = 0 tower per phi code and a full set of E_T = 0 jets: every field and flag zero";
        for (unsigned int iPhi = 0; iPhi < phi_range_; iPhi++)
            event.towers.push_back(makeMemPrintObject(iPhi, eta, 0));
        for (unsigned int iJet = 0; iJet < maxJetsConsidered_; iJet++)
            event.jets.push_back(makeMemPrintObject((iJet * phi_range_) / maxJetsConsidered_, eta, 0));
        events.push_back(event);
    }

    // --- upstream overflow ---
    {
        MetValidationEvent event;
        event.name = "upstream_tower";
        event.description = "one tower at E_T = 8191: upstream flag on tower and total; "
                            "a lone 8191 cannot fit both components in +-4095, so Ex/Ey flags fire too";
        event.towers = {makeMemPrintObject(phiDiag, eta, met_random_et_saturated_)};
        events.push_back(event);
    }
    {
        MetValidationEvent event;
        event.name = "upstream_jet";
        event.description = "one jet at E_T = 8191: upstream flag on jet and total; Ex/Ey flags fire too";
        event.jets = {makeMemPrintObject(phiDiag, eta, met_random_et_saturated_)};
        events.push_back(event);
    }
    {
        MetValidationEvent event;
        event.name = "upstream_both_cancel";
        event.description = "an E_T = 8191 tower and jet back-to-back: total MET cancels, but the total "
                            "still inherits the upstream flag (OR) and overflows SumET";
        event.towers = {makeMemPrintObject(phiDiag,    eta, met_random_et_saturated_)};
        event.jets   = {makeMemPrintObject(phiDiagOpp, eta, met_random_et_saturated_)};
        events.push_back(event);
    }

    // --- phi sweep ---
    // Tower at 2000 counts and jet at 1000 counts on the same code: each term exercises its
    // own sin/cos lookup and the phi comparator, and the total (3000 counts) stays inside
    // +-4095 on every component at every angle, so no flag fires anywhere in the sweep.
    for (unsigned int iPhi = 0; iPhi < phi_range_; iPhi++) {
        MetValidationEvent event;
        event.name = "phi_sweep_" + std::to_string(iPhi);
        event.description = "one tower (E_T 2000) and one jet (E_T 1000) at phi code " + std::to_string(iPhi) + ", no flags";
        event.towers = {makeMemPrintObject(iPhi, iPhi % eta_range_, 2000)};
        event.jets   = {makeMemPrintObject(iPhi, iPhi % eta_range_, 1000)};
        events.push_back(event);
    }

    // --- Ex / Ey boundaries ---
    const int fitMax = static_cast<int>(maskN(met_tob_signed_mag_bit_length_)); // 4095
    for (bool xComponent : {true, false}) {
        for (int target : {fitMax, fitMax + 1, -fitMax, -(fitMax + 1)})
            events.push_back(makeComponentBoundaryEvent(xComponent, target));
    }

    // --- SumET ---
    // Back-to-back pairs, so MET stays at or within a count of zero and only SumET moves.
    {
        MetValidationEvent event;
        event.name = "sumet_8191";
        event.description = "towers back-to-back, SumET = 8191 exactly: largest value that fits, no SumET overflow";
        event.towers = {makeMemPrintObject(phiDiag,    eta, 4095),
                        makeMemPrintObject(phiDiagOpp, eta, 4096)};
        events.push_back(event);
    }
    {
        MetValidationEvent event;
        event.name = "sumet_8192";
        event.description = "towers back-to-back, SumET = 8192 exactly: SumET overflow, MET = 0";
        event.towers = {makeMemPrintObject(phiDiag,    eta, 4096),
                        makeMemPrintObject(phiDiagOpp, eta, 4096)};
        events.push_back(event);
    }
    {
        MetValidationEvent event;
        event.name = "sumet_balanced_large";
        event.description = "three back-to-back tower pairs at E_T 4000: SumET = 24000, overflow, MET = 0";
        for (unsigned int iPhi : {0u, 11u, 21u}) {
            event.towers.push_back(makeMemPrintObject(iPhi,                        eta, 4000));
            event.towers.push_back(makeMemPrintObject(iPhi + pi_digitized_in_phi_, eta, 4000));
        }
        events.push_back(event);
    }

    // --- MET magnitude ---
    // Pairs straddling +pi/4 give |Ex| ~= |Ey| ~= 1.41 * E_T. MET overflows once
    // Ex^2 + Ey^2 reaches 2^26, i.e. |Ex| = |Ey| ~= 5793.
    {
        MetValidationEvent event;
        event.name = "met_fits_components_overflow";
        event.description = "towers near pi/4, |Ex| ~= |Ey| ~= 5693: Ex and Ey overflow, "
                            "MET ~= 8051 still fits (sum of squares below 2^26)";
        event.towers = {makeMemPrintObject(phiQuarterLow,  eta, 4031),
                        makeMemPrintObject(phiQuarterHigh, eta, 4031)};
        events.push_back(event);
    }
    {
        MetValidationEvent event;
        event.name = "met_overflow_square";
        event.description = "towers near pi/4, |Ex| ~= |Ey| ~= 5992: MET overflow from the sum of squares "
                            "with both components below 8192; Ex, Ey, SumET overflow too";
        event.towers = {makeMemPrintObject(phiQuarterLow,  eta, 4243),
                        makeMemPrintObject(phiQuarterHigh, eta, 4243)};
        events.push_back(event);
    }
    {
        MetValidationEvent event;
        event.name = "met_overflow_component";
        event.description = "towers straddling phi = 0 at E_T 8190, Ex ~= -16360: MET overflow from a "
                            "component >= 8192; Ex and SumET overflow too";
        event.towers = {makeMemPrintObject(phiZeroLow,  eta, met_random_et_max_clean_),
                        makeMemPrintObject(phiZeroHigh, eta, met_random_et_max_clean_)};
        events.push_back(event);
    }

    // --- cancellation in the total ---
    // Flags are recomputed on the full-width total, not inherited from the terms.
    {
        MetValidationEvent event;
        event.name = "cancel_to_zero";
        event.description = "tower and jet terms each Ex ~= -/+16360 (Ex and MET overflow), "
                            "total Ex = 0: total has SumET overflow only";
        event.towers = {makeMemPrintObject(phiZeroLow, eta, met_random_et_max_clean_),
                        makeMemPrintObject(phiZeroHigh, eta, met_random_et_max_clean_)};
        event.jets   = {makeMemPrintObject(phiPiLow,   eta, met_random_et_max_clean_),
                        makeMemPrintObject(phiPiHigh,  eta, met_random_et_max_clean_)};
        events.push_back(event);
    }
    {
        MetValidationEvent event;
        event.name = "cancel_in_range";
        event.description = "tower term Ex ~= -16360 and jet term Ex ~= +13980, both overflowing; "
                            "total Ex ~= -2380 fits: total has SumET overflow only";
        event.towers = {makeMemPrintObject(phiZeroLow, eta, met_random_et_max_clean_),
                        makeMemPrintObject(phiZeroHigh, eta, met_random_et_max_clean_)};
        event.jets   = {makeMemPrintObject(phiPiLow,   eta, 7000),
                        makeMemPrintObject(phiPiHigh,  eta, 7000)};
        events.push_back(event);
    }

    return events;
}

// Clip a drawn E_T to the input field: anything at or past the top saturates to 8191,
// the same way an upstream digitization would.
inline unsigned int clipRandomEt(double et) {
    const long rounded = std::lround(et);
    if (rounded <= 0) return 0;
    if (rounded >= static_cast<long>(met_random_et_saturated_)) return met_random_et_saturated_;
    return static_cast<unsigned int>(rounded);
}

// The pseudo-random block, appended after the fixed cases. Fully determined by seed.
inline void appendMetRandomEvents(std::vector<MetValidationEvent>& events,
                                  unsigned int seed, unsigned int nRandomEvents) {
    TRandom3 rng(seed);
    for (unsigned int iRandom = 0; iRandom < nRandomEvents; iRandom++) {
        MetValidationEvent event;
        event.name = "random_" + std::to_string(iRandom);

        const unsigned int nTowers = rng.Integer(met_random_max_towers_ + 1);
        const double towerMeanEt = std::exp(rng.Uniform(std::log(met_random_tower_mean_et_min_),
                                                        std::log(met_random_tower_mean_et_max_)));
        for (unsigned int iTower = 0; iTower < nTowers; iTower++) {
            const unsigned int phi = rng.Integer(phi_range_);
            const unsigned int eta = rng.Integer(eta_range_);
            const unsigned int et  = (rng.Rndm() < met_random_tower_saturate_prob_)
                                         ? met_random_et_saturated_
                                         : clipRandomEt(rng.Exp(towerMeanEt));
            event.towers.push_back(makeMemPrintObject(phi, eta, et));
        }

        const unsigned int nJets = rng.Integer(maxJetsConsidered_ + 1);
        for (unsigned int iJet = 0; iJet < nJets; iJet++) {
            const unsigned int phi = rng.Integer(phi_range_);
            const unsigned int eta = rng.Integer(eta_range_);
            const double draw = rng.Rndm();
            unsigned int et = 0;
            if (draw < met_random_jet_saturate_prob_)
                et = met_random_et_saturated_;
            else if (draw < met_random_jet_saturate_prob_ + met_random_jet_hard_prob_)
                et = rng.Integer(met_random_et_max_clean_ + 1);
            else
                et = clipRandomEt(rng.Exp(met_random_jet_mean_et_));
            event.jets.push_back(makeMemPrintObject(phi, eta, et));
        }

        std::stringstream description;
        description << "random: " << nTowers << " towers (mean E_T " << std::fixed << std::setprecision(2)
                    << towerMeanEt << " counts), " << nJets << " jets";
        event.description = description.str();
        events.push_back(event);
    }
}

// The five TOB flags of one term, in TOB order [62] -> [58], for the _Cases.txt file.
// Raised through the same encoders writeMetBaselineLine uses.
inline std::string metTermFlagString(const MetTermDigitized& term) {
    bool exOverflow = false, eyOverflow = false, sumEtOverflow = false;
    metTobTwosComplementEt(term.metX, exOverflow);
    metTobTwosComplementEt(term.metY, eyOverflow);
    metTobSumEt(term.sumEt, sumEtOverflow);
    std::stringstream flags;
    flags << "up=" << term.inputSaturated << " met=" << term.metOverflow
          << " sum=" << sumEtOverflow << " ey=" << eyOverflow << " ex=" << exOverflow;
    return flags.str();
}

// MET edge-case / random validation entry point.
//
// Use:
//   root -b -l
//   root [0] .L metEmulation.cc+
//   root [1] metRandomValidation()            // seed 12345, 100 random events
//   root [2] metRandomValidation(777, 500)
//
// Writes to met_baseline_mem_prints_path_ + met_random_output_subdir_, all named
// METRandomValidation_seed<seed>:
//   _TowerInput.dat, _JetInput.dat    inputs, HERNTupler memory-print format plus an overflow tag at [32]
//   _METBaseline{Tower,Jet,Total}.dat outputs, same format as the baseline path
//   _Cases.txt                        per event: case name, object counts, the flags each
//                                     term raised, and what the case is testing
// plus a ROOT ntuple of the same terms in met_random_ntuple_path_, and a per-event log in
// the working directory.
//
// Event indices are shared by all of these: "Event : N" in every .dat, line N of the
// cases file, eventIndex N in the ntuple. Empty events still get their "Event : N" header
// in the input files, with no object lines under it.
void metRandomValidation(unsigned int seed = 12345,          // TRandom3 seed; must be nonzero
                         unsigned int nRandomEvents = 100) { // random events after the fixed block

    if (seed == 0) {
        // TRandom3(0) seeds from the clock, so the files could never be regenerated.
        std::cerr << "Error: seed 0 asks TRandom3 for a time-based seed; pick a nonzero seed" << std::endl;
        return;
    }

    std::vector<MetValidationEvent> events = makeMetEdgeCaseEvents();
    const unsigned int nFixedEvents = events.size();
    appendMetRandomEvents(events, seed, nRandomEvents);
    const unsigned int nEvents = events.size();

    const std::string baseName        = met_random_base_prefix_ + std::to_string(seed);
    const std::string outputPath      = met_baseline_mem_prints_path_ + met_random_output_subdir_;
    const std::string towerInputFile  = outputPath + baseName + met_random_tower_input_suffix_;
    const std::string jetInputFile    = outputPath + baseName + met_random_jet_input_suffix_;
    const std::string towerOutputFile = outputPath + baseName + met_baseline_tower_output_suffix_;
    const std::string jetOutputFile   = outputPath + baseName + met_baseline_jet_output_suffix_;
    const std::string totalOutputFile = outputPath + baseName + met_baseline_total_output_suffix_;
    const std::string casesFile       = outputPath + baseName + met_random_cases_suffix_;
    const std::string ntupleFile      = met_random_ntuple_path_ + baseName + met_random_ntuple_suffix_;
    const std::string logFile         = met_random_log_prefix_ + baseName + met_random_log_suffix_;

    std::cout << "tower input:   " << towerInputFile  << "\n";
    std::cout << "jet input:     " << jetInputFile    << "\n";
    std::cout << "tower output:  " << towerOutputFile << "\n";
    std::cout << "jet output:    " << jetOutputFile   << "\n";
    std::cout << "total output:  " << totalOutputFile << "\n";
    std::cout << "cases:         " << casesFile       << "\n";
    std::cout << "ntuple output: " << ntupleFile      << "\n";

    // --- write the inputs ---
    gSystem->mkdir(outputPath.c_str(), true);
    {
        std::ofstream towerInFile(towerInputFile);
        std::ofstream jetInFile(jetInputFile);
        if (!towerInFile.is_open() || !jetInFile.is_open()) {
            std::cerr << "Error: Could not open " << towerInputFile << " or " << jetInputFile << std::endl;
            return;
        }
        for (unsigned int iEvt = 0; iEvt < nEvents; iEvt++) {
            towerInFile << "Event : " << std::dec << iEvt << std::endl;
            jetInFile   << "Event : " << std::dec << iEvt << std::endl;
            for (unsigned int iObject = 0; iObject < events.at(iEvt).towers.size(); iObject++)
                writeMemPrintObjectLine(towerInFile, iObject, events.at(iEvt).towers.at(iObject));
            // Jets are generated on the tower grid; written as jet TOB codes if input jets carry them
            for (unsigned int iObject = 0; iObject < events.at(iEvt).jets.size(); iObject++) {
                MemPrintObject jet = events.at(iEvt).jets.at(iObject);
                if (input_jets_use_jet_tob_encoding_) {
                    jet.eta = tower_index_to_jet_eta(jet.eta);
                    jet.phi = tower_index_to_jet_phi(jet.phi);
                }
                writeMemPrintObjectLine(jetInFile, iObject, jet);
            }
        }
    }

    // --- read them back: the outputs are computed from the files, not from memory ---
    const std::vector<std::vector<MemPrintObject>> towerEvents = readMemPrintFileWithOverflow(towerInputFile);
    const std::vector<std::vector<MemPrintObject>> jetEvents   = jetMemPrintsToTowerGrid(readMemPrintFileWithOverflow(jetInputFile));
    if (towerEvents.size() != nEvents || jetEvents.size() != nEvents) {
        std::cerr << "Error: wrote " << nEvents << " events but read back " << towerEvents.size()
                  << " tower / " << jetEvents.size() << " jet events; nothing else written" << std::endl;
        return;
    }
    for (unsigned int iEvt = 0; iEvt < nEvents; iEvt++) {
        if (towerEvents.at(iEvt).size() != events.at(iEvt).towers.size() ||
            jetEvents.at(iEvt).size()   != events.at(iEvt).jets.size())
            std::cerr << "Warning: event " << iEvt << " (" << events.at(iEvt).name
                      << ") did not survive the round trip through the input files" << std::endl;
    }

    // --- open the outputs ---
    std::ofstream towerOutFile(towerOutputFile);
    std::ofstream jetOutFile(jetOutputFile);
    std::ofstream totalOutFile(totalOutputFile);
    std::ofstream casesOutFile(casesFile);
    if (!towerOutFile.is_open() || !jetOutFile.is_open() || !totalOutFile.is_open() || !casesOutFile.is_open()) {
        std::cerr << "Error: Could not open one of the output files in " << outputPath << std::endl;
        return;
    }

    gSystem->mkdir(met_random_ntuple_path_.c_str(), true);
    TFile* ntupleOutFile = TFile::Open(ntupleFile.c_str(), "RECREATE");
    if (!ntupleOutFile || ntupleOutFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << ntupleFile << std::endl;
        return;
    }
    ntupleOutFile->cd();
    TTree* metRandomTree = new TTree("metRandomValidationTree",
                                     "Per-event MET edge-case / random validation terms, digitized and in GeV");
    unsigned int out_eventIndex = 0;
    bool         out_isRandom   = false;
    std::string  out_caseName;
    unsigned int out_nTowers    = 0;
    unsigned int out_nJets      = 0;
    metRandomTree->Branch("eventIndex", &out_eventIndex);
    metRandomTree->Branch("isRandom",   &out_isRandom);
    metRandomTree->Branch("caseName",   &out_caseName);
    metRandomTree->Branch("nTowers",    &out_nTowers);
    metRandomTree->Branch("nJets",      &out_nJets);

    MetBaselineNtupleTerm towerNtuple, jetNtuple, totalNtuple;
    towerNtuple.book(metRandomTree, "Tower");
    jetNtuple.book(metRandomTree,   "Jet");
    totalNtuple.book(metRandomTree, "Total");

    casesOutFile << "# MET edge-case / random validation, seed " << seed << ": "
                 << nFixedEvents << " fixed events, then " << nRandomEvents << " random\n"
                 << "# flags in TOB order: up=[62] met=[61] sum=[60] ey=[59] ex=[58]\n";

    // Per-event printouts to the log; same restorable redirect as the baseline path, and
    // every early return is above it.
    RedirectHandle_t redirectHandle;
    gSystem->RedirectOutput(logFile.c_str(), "w", &redirectHandle);

    std::cout << "MET edge-case / random validation: " << baseName << "\n";
    std::cout << "events: " << nEvents << " (" << nFixedEvents << " fixed, " << nRandomEvents << " random)\n";

    auto printTerm = [](const char* label, const MetTermDigitized& term) {
        std::cout << "  " << label
                  << " Ex: "  << term.metX  << " (" << term.metX * et_granularity_ << " GeV)"
                  << "  Ey: " << term.metY  << " (" << term.metY * et_granularity_ << " GeV)"
                  << "  MET: "<< term.met   << " (" << term.met  * et_granularity_ << " GeV)"
                  << "  phi: "<< term.phi   << " (" << undigitize_met_phi(term.phi) << " rad)"
                  << "  SumET: " << term.sumEt << " (" << term.sumEt * et_granularity_ << " GeV)"
                  << "  [" << metTermFlagString(term) << "]\n";
    };

    for (unsigned int iEvt = 0; iEvt < nEvents; iEvt++) {
        towerOutFile << "Event : " << std::dec << iEvt << std::endl;
        jetOutFile   << "Event : " << std::dec << iEvt << std::endl;
        totalOutFile << "Event : " << std::dec << iEvt << std::endl;

        const MetTermDigitized towerTerm = accumulateMemPrintMet(towerEvents.at(iEvt), maxTowersConsidered_);
        const MetTermDigitized jetTerm   = accumulateMemPrintMet(jetEvents.at(iEvt),   maxJetsConsidered_);
        const MetTermDigitized totalTerm = combineMetTerms(towerTerm, jetTerm);

        std::cout << "iEvt: " << std::dec << iEvt << "  " << events.at(iEvt).name
                  << "  towers: " << towerEvents.at(iEvt).size()
                  << "  jets: "   << jetEvents.at(iEvt).size() << "\n";
        printTerm("tower ", towerTerm);
        printTerm("jet   ", jetTerm);
        printTerm("total ", totalTerm);

        if (met_random_debug_) {
            for (unsigned int iObject = 0; iObject < towerEvents.at(iEvt).size(); iObject++) {
                const MemPrintObject& object = towerEvents.at(iEvt).at(iObject);
                std::cout << "    tower " << iObject << " phi: " << object.phi
                          << " eta: " << object.eta << " et: " << object.et << "\n";
            }
            for (unsigned int iObject = 0; iObject < jetEvents.at(iEvt).size(); iObject++) {
                const MemPrintObject& object = jetEvents.at(iEvt).at(iObject);
                std::cout << "    jet   " << iObject << " phi: " << object.phi
                          << " eta: " << object.eta << " et: " << object.et << "\n";
            }
        }

        writeMetBaselineLine(towerOutFile, 0, towerTerm);
        writeMetBaselineLine(jetOutFile,   0, jetTerm);
        writeMetBaselineLine(totalOutFile, 0, totalTerm);

        casesOutFile << "Event : " << std::dec << iEvt << "  " << events.at(iEvt).name
                     << "  towers=" << towerEvents.at(iEvt).size() << " jets=" << jetEvents.at(iEvt).size()
                     << "  tower[" << metTermFlagString(towerTerm) << "]"
                     << "  jet["   << metTermFlagString(jetTerm)   << "]"
                     << "  total[" << metTermFlagString(totalTerm) << "]"
                     << "  -- " << events.at(iEvt).description << "\n";

        out_eventIndex = iEvt;
        out_isRandom   = (iEvt >= nFixedEvents);
        out_caseName   = events.at(iEvt).name;
        out_nTowers    = towerEvents.at(iEvt).size();
        out_nJets      = jetEvents.at(iEvt).size();
        towerNtuple.fill(towerTerm);
        jetNtuple.fill(jetTerm);
        totalNtuple.fill(totalTerm);
        metRandomTree->Fill();
    }

    gSystem->RedirectOutput(nullptr, "", &redirectHandle);

    towerOutFile.close();
    jetOutFile.close();
    totalOutFile.close();
    casesOutFile.close();

    ntupleOutFile->cd();
    metRandomTree->Write("", TObject::kOverwrite);
    ntupleOutFile->Close();

    std::cout << "Wrote " << std::dec << nEvents << " events (" << nFixedEvents << " fixed, "
              << nRandomEvents << " random) to " << outputPath << "\n";
    std::cout << "Wrote " << std::dec << nEvents << " events to " << ntupleFile << "\n";
    std::cout << "Per-event printouts in " << logFile << "\n";
}
