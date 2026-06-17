
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
#include "emulationHelperFunctions.h"

// Function for processing provided number of loops with JetTagger algorithm,
// then writing output large radius jets to a new TTree & text files
void eventLoop(std::string inputNTuplePath, std::string outputNTuplePath,std::string outputTextFilePath, std::string inputObjectType, std::string seedObjectType,
        bool useSKObjects, bool enableEtWeightedMidpoint, bool minEtSeedPosOptimization,
        double minEtSeedPosOptimizationCut, double subjetEtThreshold, // Both in GeV --> converted to digitized 125 MeV units later
        bool useEtaSKObjects = false
    ){
        std::cout << "inputobjecttype: " << inputObjectType << "\n";
        std::cout << "useSKObjects: " << useSKObjects << "\n";
    struct outputJet {
        // Main values
        unsigned int et;
        unsigned int eta;
        unsigned int phi;

        // Leading LRJ
        unsigned int nSubjets;
        unsigned int psiR;
        unsigned int tau_2;
        unsigned int tau_1;
        unsigned int massApprox;
    };

    struct inputObject {
        // Main values
        unsigned int et;
        unsigned int eta;
        unsigned int phi;
    };

    // Relevant data for algorithm read from inputNTuple

    // GEPCellsTowers
    std::vector<double>* gepCellsTowersEtValues = nullptr;
    std::vector<double>* gepCellsTowersEtaValues = nullptr;
    std::vector<double>* gepCellsTowersPhiValues = nullptr;

    // GEPCellsTowers (PU-suppressed)
    std::vector<double>* gepCellsTowersSKEtValues = nullptr;
    std::vector<double>* gepCellsTowersSKEtaValues = nullptr;
    std::vector<double>* gepCellsTowersSKPhiValues = nullptr;

    // GEPBasicClusters
    std::vector<double>* gepBasicClustersEtValues = nullptr;
    std::vector<double>* gepBasicClustersEtaValues = nullptr;
    std::vector<double>* gepBasicClustersPhiValues = nullptr;

    // GEPBasicClustersSK (PU-suppressed)
    std::vector<double>* gepBasicClustersSKEtValues = nullptr;
    std::vector<double>* gepBasicClustersSKEtaValues = nullptr;
    std::vector<double>* gepBasicClustersSKPhiValues = nullptr;

    // Cone jets from TrigGepPerf
    std::vector<double>* gepWTAConeCellsTowersJetsEtValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersJetsEtaValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersJetsPhiValues = nullptr;
    std::vector<unsigned int >* gepWTAConeCellsTowersJetsNConstituentsValues = nullptr;

    // Cone jets from TrigGepPerf (PU-suppressed)
    std::vector<double>* gepWTAConeCellsTowersSKJetsEtValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersSKJetsEtaValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersSKJetsPhiValues = nullptr;
    std::vector<unsigned int >* gepWTAConeCellsTowersSKJetsNConstituentsValues = nullptr;

    // GEPCellsTowers (EtaSK PU-suppressed)
    std::vector<double>* gepCellsTowersEtaSKEtValues = nullptr;
    std::vector<double>* gepCellsTowersEtaSKEtaValues = nullptr;
    std::vector<double>* gepCellsTowersEtaSKPhiValues = nullptr;

    // Cone jets from TrigGepPerf (EtaSK PU-suppressed)
    std::vector<double>* gepWTAConeCellsTowersEtaSKJetsEtValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersEtaSKJetsEtaValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersEtaSKJetsPhiValues = nullptr;
    std::vector<unsigned int>* gepWTAConeCellsTowersEtaSKJetsNConstituentsValues = nullptr;

    // gFEX small radius jets
    std::vector<unsigned int>* gFexSRJEtIndexValues = nullptr;
    std::vector<double>* gFexSRJEtValues = nullptr;
    std::vector<double>* gFexSRJEtaValues = nullptr;
    std::vector<double>* gFexSRJPhiValues = nullptr;

    // jFEX small radius jets
    std::vector<unsigned int>* jFexSRJEtIndexValues = nullptr;
    std::vector<double>* jFexSRJEtValues = nullptr;
    std::vector<double>* jFexSRJEtaValues = nullptr;
    std::vector<double>* jFexSRJPhiValues = nullptr;

    // Output JetTagger large radius jets
    std::vector<double>* jetTaggerLRJPsi_RValues = nullptr;
    std::vector<double>* jetTaggerLRJMassApproxValues = nullptr;
    std::vector<double>* jetTaggerLRJTau_1Values = nullptr;
    std::vector<double>* jetTaggerLRJTau_2Values = nullptr;
    std::vector<double>* jetTaggerLRJTau_21Values = nullptr;
    std::vector<unsigned int >* jetTaggerLRJSubjetMultiplicityValues = nullptr;
    std::vector<unsigned int>* jetTaggerLRJMergedIndicesValues = nullptr;
    std::vector<double>* jetTaggerLRJEtValues = nullptr;
    std::vector<double>* jetTaggerLRJEtaValues = nullptr;
    std::vector<double>* jetTaggerLRJPhiValues = nullptr;

    // Output JetTagger leading large radius jets
    std::vector<double>* jetTaggerLeadingLRJPsi_RValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJMassApproxValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJTau_1Values = nullptr;
    std::vector<double>* jetTaggerLeadingLRJTau_2Values = nullptr;
    std::vector<double>* jetTaggerLeadingLRJTau_21Values = nullptr;
    std::vector<unsigned int >* jetTaggerLeadingLRJSubjetMultiplicityValues = nullptr;
    std::vector<unsigned int>* jetTaggerLeadingLRJMergedIndicesValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJEtValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJEtaValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJPhiValues = nullptr;

    // Output JetTagger subleading large radius jets
    std::vector<double>* jetTaggerSubleadingLRJPsi_RValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJMassApproxValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJTau_1Values = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJTau_2Values = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJTau_21Values = nullptr;
    std::vector<unsigned int >* jetTaggerSubleadingLRJSubjetMultiplicityValues = nullptr;
    std::vector<unsigned int>* jetTaggerSubleadingLRJMergedIndicesValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJEtValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJEtaValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJPhiValues = nullptr;

    // Subjet Et, Eta, Phi (undigitized) for leading and subleading LRJs
    std::vector<double>* jetTaggerLeadingLRJSubjetEtValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJSubjetEtaValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJSubjetPhiValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJSubjetEtValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJSubjetEtaValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJSubjetPhiValues = nullptr;

    // Event/run number passthrough (for ordering-alignment validation)
    int gepRunNumberIn = 0, gepEventNumberIn = 0;

    int gepRunNumberOut = 0, gepEventNumberOut = 0;

    //std::cout << "test 2 " << "\n";

    // Open input ROOT file
    TFile* inputFile = TFile::Open(inputNTuplePath.c_str(), "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << inputNTuplePath << std::endl;
        return;
    }

    // Create new output ROOT file containing only jet tagger trees
    TFile* outputFile = TFile::Open(outputNTuplePath.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not open file " << outputNTuplePath << std::endl;
        return;
    }

    //std::cout << "test 3 " << "\n";

    // TTrees and branches required as input to algorithm
    TTree* caloTopoTowerTree = (TTree*)inputFile->Get("caloTopoTowerTree");
    TTree* topo422Tree = (TTree*)inputFile->Get("topo422Tree");
    TTree* gepBasicClustersTree = (TTree*)inputFile->Get("gepBasicClustersTree");
    TTree* gepCellsTowersTree = (TTree*)inputFile->Get("gepCellsTowersTree");
    TTree* gepBasicClustersSKTree = (TTree*)inputFile->Get("gepBasicClustersSKTree");
    TTree* gepCellsTowersSKTree = (TTree*)inputFile->Get("gepCellsTowersSKTree");
    TTree* gepWTAConeCellsTowersJetsTree = (TTree*)inputFile->Get("gepWTAConeCellsTowersJetsTree");
    TTree* gepWTAConeCellsTowersSKJetsTree = (TTree*)inputFile->Get("gepWTAConeCellsTowersSKJetsTree");
    TTree* gepCellsTowersEtaSKTree = (TTree*)inputFile->Get("gepCellsTowersEtaSKTree");
    TTree* gepWTAConeCellsTowersEtaSKJetsTree = (TTree*)inputFile->Get("gepWTAConeCellsTowersEtaSKJetsTree");
    TTree* gFexSRJTree = (TTree*)inputFile->Get("gFexSRJTree");
    TTree* jFexSRJTree = (TTree*)inputFile->Get("jFexSRJTree"); // Note that jFEX, gFEX trees are pre-sorted by Et in LRJNTupler.cc
    TTree* eventInfoTreeIn = (TTree*)inputFile->Get("eventInfoTree");
    eventInfoTreeIn->SetBranchAddress("gepRunNumberOut",   &gepRunNumberIn);
    eventInfoTreeIn->SetBranchAddress("gepEventNumberOut", &gepEventNumberIn);

    outputFile->cd();
    // Event/run number passthrough tree (one entry per event, for ordering-alignment validation)
    TTree* emulEventInfoTree = new TTree("emulEventInfoTree", "Run/event number passthrough from HERNTupler input for ordering validation");
    emulEventInfoTree->Branch("gepRunNumberOut",   &gepRunNumberOut);
    emulEventInfoTree->Branch("gepEventNumberOut", &gepEventNumberOut);

    // Create new output TTrees to be written to the output file
    TTree* jetTaggerLRJs = new TTree("jetTaggerLRJsTree", "Tree storing event-wise Substructure variable, Et, Eta, Phi");
    jetTaggerLRJs->Branch("Psi_R", &jetTaggerLRJPsi_RValues);
    jetTaggerLRJs->Branch("MergedIndices", &jetTaggerLRJMergedIndicesValues);
    jetTaggerLRJs->Branch("MassApprox", &jetTaggerLRJMassApproxValues);
    jetTaggerLRJs->Branch("Tau_1", &jetTaggerLRJTau_1Values);
    jetTaggerLRJs->Branch("Tau_2", &jetTaggerLRJTau_2Values);
    jetTaggerLRJs->Branch("Tau_21", &jetTaggerLRJTau_21Values);
    jetTaggerLRJs->Branch("SubjetMultiplicity", &jetTaggerLRJSubjetMultiplicityValues);
    jetTaggerLRJs->Branch("Et", &jetTaggerLRJEtValues);
    jetTaggerLRJs->Branch("Eta", &jetTaggerLRJEtaValues);
    jetTaggerLRJs->Branch("Phi", &jetTaggerLRJPhiValues);

    TTree* jetTaggerLeadingLRJs = new TTree("jetTaggerLeadingLRJsTree", "Tree storing event-wise Substructure variable, Et, Eta, Phi");
    jetTaggerLeadingLRJs->Branch("Psi_R", &jetTaggerLeadingLRJPsi_RValues);
    jetTaggerLeadingLRJs->Branch("MergedIndices", &jetTaggerLeadingLRJMergedIndicesValues);
    jetTaggerLeadingLRJs->Branch("MassApprox", &jetTaggerLeadingLRJMassApproxValues);
    jetTaggerLeadingLRJs->Branch("Tau_1", &jetTaggerLeadingLRJTau_1Values);
    jetTaggerLeadingLRJs->Branch("Tau_2", &jetTaggerLeadingLRJTau_2Values);
    jetTaggerLeadingLRJs->Branch("Tau_21", &jetTaggerLeadingLRJTau_21Values);
    jetTaggerLeadingLRJs->Branch("SubjetMultiplicity", &jetTaggerLeadingLRJSubjetMultiplicityValues);
    jetTaggerLeadingLRJs->Branch("Et", &jetTaggerLeadingLRJEtValues);
    jetTaggerLeadingLRJs->Branch("Eta", &jetTaggerLeadingLRJEtaValues);
    jetTaggerLeadingLRJs->Branch("Phi", &jetTaggerLeadingLRJPhiValues);
    jetTaggerLeadingLRJs->Branch("SubjetEt",  &jetTaggerLeadingLRJSubjetEtValues);
    jetTaggerLeadingLRJs->Branch("SubjetEta", &jetTaggerLeadingLRJSubjetEtaValues);
    jetTaggerLeadingLRJs->Branch("SubjetPhi", &jetTaggerLeadingLRJSubjetPhiValues);

    TTree* jetTaggerSubleadingLRJs = new TTree("jetTaggerSubleadingLRJsTree", "Tree storing event-wise Substructure variable, Et, Eta, Phi");
    jetTaggerSubleadingLRJs->Branch("Psi_R", &jetTaggerSubleadingLRJPsi_RValues);
    jetTaggerSubleadingLRJs->Branch("MassApprox", &jetTaggerSubleadingLRJMassApproxValues);
    jetTaggerSubleadingLRJs->Branch("MergedIndices", &jetTaggerSubleadingLRJMergedIndicesValues);
    jetTaggerSubleadingLRJs->Branch("Tau_1", &jetTaggerSubleadingLRJTau_1Values);
    jetTaggerSubleadingLRJs->Branch("Tau_2", &jetTaggerSubleadingLRJTau_2Values);
    jetTaggerSubleadingLRJs->Branch("Tau_21", &jetTaggerSubleadingLRJTau_21Values);
    jetTaggerSubleadingLRJs->Branch("SubjetMultiplicity", &jetTaggerSubleadingLRJSubjetMultiplicityValues);
    jetTaggerSubleadingLRJs->Branch("Et", &jetTaggerSubleadingLRJEtValues);
    jetTaggerSubleadingLRJs->Branch("Eta", &jetTaggerSubleadingLRJEtaValues);
    jetTaggerSubleadingLRJs->Branch("Phi", &jetTaggerSubleadingLRJPhiValues);
    jetTaggerSubleadingLRJs->Branch("SubjetEt",  &jetTaggerSubleadingLRJSubjetEtValues);
    jetTaggerSubleadingLRJs->Branch("SubjetEta", &jetTaggerSubleadingLRJSubjetEtaValues);
    jetTaggerSubleadingLRJs->Branch("SubjetPhi", &jetTaggerSubleadingLRJSubjetPhiValues);

    // === gepBasicClustersTree ===
    gepBasicClustersTree->SetBranchAddress("Et", &gepBasicClustersEtValues);
    gepBasicClustersTree->SetBranchAddress("Eta", &gepBasicClustersEtaValues);
    gepBasicClustersTree->SetBranchAddress("Phi", &gepBasicClustersPhiValues);

    // === gepBasicClustersSKTree (PU-suppressed) ===
    gepBasicClustersSKTree->SetBranchAddress("Et", &gepBasicClustersSKEtValues);
    gepBasicClustersSKTree->SetBranchAddress("Eta", &gepBasicClustersSKEtaValues);
    gepBasicClustersSKTree->SetBranchAddress("Phi", &gepBasicClustersSKPhiValues);

    // === gepCellsTowersTree ===
    gepCellsTowersTree->SetBranchAddress("Et", &gepCellsTowersEtValues);
    gepCellsTowersTree->SetBranchAddress("Eta", &gepCellsTowersEtaValues);
    gepCellsTowersTree->SetBranchAddress("Phi", &gepCellsTowersPhiValues);

    // === gepCellsTowersSKTree (PU-suppressed) ===
    gepCellsTowersSKTree->SetBranchAddress("Et", &gepCellsTowersSKEtValues);
    gepCellsTowersSKTree->SetBranchAddress("Eta", &gepCellsTowersSKEtaValues);
    gepCellsTowersSKTree->SetBranchAddress("Phi", &gepCellsTowersSKPhiValues);

    // gep WTA cone cells towers jets
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Et", &gepWTAConeCellsTowersJetsEtValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Eta", &gepWTAConeCellsTowersJetsEtaValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Phi", &gepWTAConeCellsTowersJetsPhiValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("NConstituents", &gepWTAConeCellsTowersJetsNConstituentsValues);

    // gep WTA cone cells towers SK jets (PU-suppressed)
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Et", &gepWTAConeCellsTowersSKJetsEtValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Eta", &gepWTAConeCellsTowersSKJetsEtaValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Phi", &gepWTAConeCellsTowersSKJetsPhiValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("NConstituents", &gepWTAConeCellsTowersSKJetsNConstituentsValues);

    // === gepCellsTowersEtaSKTree (EtaSK PU-suppressed) ===
    gepCellsTowersEtaSKTree->SetBranchAddress("Et",  &gepCellsTowersEtaSKEtValues);
    gepCellsTowersEtaSKTree->SetBranchAddress("Eta", &gepCellsTowersEtaSKEtaValues);
    gepCellsTowersEtaSKTree->SetBranchAddress("Phi", &gepCellsTowersEtaSKPhiValues);

    // gep WTA cone cells towers EtaSK jets (EtaSK PU-suppressed)
    gepWTAConeCellsTowersEtaSKJetsTree->SetBranchAddress("Et",           &gepWTAConeCellsTowersEtaSKJetsEtValues);
    gepWTAConeCellsTowersEtaSKJetsTree->SetBranchAddress("Eta",          &gepWTAConeCellsTowersEtaSKJetsEtaValues);
    gepWTAConeCellsTowersEtaSKJetsTree->SetBranchAddress("Phi",          &gepWTAConeCellsTowersEtaSKJetsPhiValues);
    gepWTAConeCellsTowersEtaSKJetsTree->SetBranchAddress("NConstituents",&gepWTAConeCellsTowersEtaSKJetsNConstituentsValues);

    // === gFexSRJTree ===
    gFexSRJTree->SetBranchAddress("EtIndex", &gFexSRJEtIndexValues);
    gFexSRJTree->SetBranchAddress("Et", &gFexSRJEtValues);
    gFexSRJTree->SetBranchAddress("Eta", &gFexSRJEtaValues);
    gFexSRJTree->SetBranchAddress("Phi", &gFexSRJPhiValues);

    // === jFexSRJTree ===
    jFexSRJTree->SetBranchAddress("EtIndex", &jFexSRJEtIndexValues);
    jFexSRJTree->SetBranchAddress("Et", &jFexSRJEtValues);
    jFexSRJTree->SetBranchAddress("Eta", &jFexSRJEtaValues);
    jFexSRJTree->SetBranchAddress("Phi", &jFexSRJPhiValues);

    std::ofstream f_output(outputTextFilePath);

    if (!f_output.is_open()) {
        std::cerr << "Error: One or more output files could not be opened for writing!" << std::endl;
        // handle error or return
    }

    // Event loop 
    unsigned int overlapRemovalCounter = 0;
    unsigned int eventsToProcess = gepBasicClustersTree->GetEntries();
    for (unsigned int iEvt = 0; iEvt < eventsToProcess; iEvt++) {
        //std::cout << "-----------------------" << "\n";
        //std::cout << "iEvt: " << iEvt << "\n";
        // Write event header to output text files
        f_output << "Event : " << std::dec << iEvt << "\n";

        //std::cout << "procesing event: " << std::dec << iEvt << "\n";

        std::bitset<num_subjets_length_ > num_subjets_bitset_array[nSeedsOutput_];
        std::bitset<N_subjetiness_bit_length_ > tau_1_bitset_array[nSeedsOutput_];
        std::bitset<N_subjetiness_bit_length_ > tau_2_bitset_array[nSeedsOutput_];
        std::bitset<mass_approx_bit_length_ > massApprox_bitset_array[nSeedsOutput_];
        std::bitset<psi_R_bit_length_ > psi_R_bitset_array[nSeedsOutput_];
        std::bitset<et_bit_length_  > et_bitset_array[nSeedsOutput_];
        std::bitset<eta_bit_length_ > eta_bitset_array[nSeedsOutput_];
        std::bitset<phi_bit_length_ > phi_bitset_array[nSeedsOutput_];

        double tau_21_array[nSeedsOutput_]; // NOTE: store tau_21 as just a double, since this wouldn't be computed within Jet Tagger algorithm anyways, but is still useful to have in the ntuple

        // Clear vectors for branches
        jetTaggerLRJPsi_RValues->clear();
        
        jetTaggerLRJMassApproxValues->clear();
        jetTaggerLRJTau_1Values->clear();
        jetTaggerLRJTau_2Values->clear();
        jetTaggerLRJTau_21Values->clear();
        jetTaggerLRJSubjetMultiplicityValues->clear();
        jetTaggerLRJEtValues->clear();
        jetTaggerLRJEtaValues->clear();
        jetTaggerLRJPhiValues->clear();
        jetTaggerLRJMergedIndicesValues->clear();

        jetTaggerLeadingLRJPsi_RValues->clear();
        jetTaggerLeadingLRJMergedIndicesValues->clear();
        jetTaggerLeadingLRJMassApproxValues->clear();
        jetTaggerLeadingLRJTau_1Values->clear();
        jetTaggerLeadingLRJTau_2Values->clear();
        jetTaggerLeadingLRJTau_21Values->clear();
        jetTaggerLeadingLRJSubjetMultiplicityValues->clear();
        jetTaggerLeadingLRJEtValues->clear();
        jetTaggerLeadingLRJEtaValues->clear();
        jetTaggerLeadingLRJPhiValues->clear();

        jetTaggerSubleadingLRJPsi_RValues->clear();
        jetTaggerSubleadingLRJMergedIndicesValues->clear();
        jetTaggerSubleadingLRJMassApproxValues->clear();
        jetTaggerSubleadingLRJTau_1Values->clear();
        jetTaggerSubleadingLRJTau_2Values->clear();
        jetTaggerSubleadingLRJTau_21Values->clear();
        jetTaggerSubleadingLRJSubjetMultiplicityValues->clear();
        jetTaggerSubleadingLRJEtValues->clear();
        jetTaggerSubleadingLRJEtaValues->clear();
        jetTaggerSubleadingLRJPhiValues->clear();

        jetTaggerLeadingLRJSubjetEtValues->clear();
        jetTaggerLeadingLRJSubjetEtaValues->clear();
        jetTaggerLeadingLRJSubjetPhiValues->clear();
        jetTaggerSubleadingLRJSubjetEtValues->clear();
        jetTaggerSubleadingLRJSubjetEtaValues->clear();
        jetTaggerSubleadingLRJSubjetPhiValues->clear();

        eventInfoTreeIn->GetEntry(iEvt); // get event/run number for passthrough

        gepRunNumberOut = gepRunNumberIn;
        gepEventNumberOut = gepEventNumberIn;

        if(seedObjectType == "jFEXSRJ"){
            jFexSRJTree->GetEntry(iEvt); // get seed data
        }
        else if(seedObjectType == "gFEXSRJ"){
            gFexSRJTree->GetEntry(iEvt); // get seed data
        }
        else if(seedObjectType == "gepWTAConeCellsTowersJets"){
            if(useEtaSKObjects){
                gepWTAConeCellsTowersEtaSKJetsTree->GetEntry(iEvt); // get seed object data (EtaSK PU-suppressed)
            }
            else if(useSKObjects){
                gepWTAConeCellsTowersSKJetsTree->GetEntry(iEvt); // get seed object data (PU-suppressed)
            }
            else{
                gepWTAConeCellsTowersJetsTree->GetEntry(iEvt); // get seed object data
            }
        }
        if(inputObjectType == "gepBasicClusters"){
            if(useSKObjects){
                gepBasicClustersSKTree->GetEntry(iEvt); // get input object data (PU-suppressed)
            }
            else{
                gepBasicClustersTree->GetEntry(iEvt); // get input object data
            }   
            
        }
        else if(inputObjectType == "gepCellsTowers"){
            if(useEtaSKObjects){
                gepCellsTowersEtaSKTree->GetEntry(iEvt); // get input object data (EtaSK PU-suppressed)
            }
            else if(useSKObjects){
                gepCellsTowersSKTree->GetEntry(iEvt); // get input object data (PU-suppressed)
            }
            else{
                gepCellsTowersTree->GetEntry(iEvt); // get input object data
            }
        }
        else if(inputObjectType == "gepWTAConeCellsTowersJets"){
            if(useEtaSKObjects){
                gepWTAConeCellsTowersEtaSKJetsTree->GetEntry(iEvt); // get input object data (EtaSK PU-suppressed)
            }
            else if(useSKObjects){
                gepWTAConeCellsTowersSKJetsTree->GetEntry(iEvt); // get input object data (PU-suppressed)
            }
            else{
                gepWTAConeCellsTowersJetsTree->GetEntry(iEvt); // get input object data
            }
        }
        else {
            std::cerr << "INPUT OBJECT UNKNOWN!" << "\n";
        }

        inputObject seedValuesOriginal[nSeedsInput_] = {}; // 3 = number of data types read in (Et, Eta, Phi)
        inputObject seedValues[nSeedsInput_] = {}; // 3 = number of data types read in (Et, Eta, Phi)
        inputObject inputObjectValues[maxObjectsConsidered_] = {}; // using struct now - stores number of data types read in (Et, Eta, Phi)
        outputJet outputJetValues[nSeedsOutput_] = {};         // using struct now - stores number of possible output data types (Et, Eta, Phi, # subjets (both LRJs), psi_R (both LRJs), 
                                                               // N subjetiness? (both LRJs), massApprox (both LRJs))
                                                               // NOTE: not all of these will be written out in firmware algo, but can be used for physics studies
        // Digitize and format data in same way expected for running in HLS
        for(unsigned int iSeed = 0; iSeed < nSeedsInput_; iSeed++){
            //std::cout << "iSeed :" << iSeed << "\n";
            //std::cout << "jFexSRJEtValues->size(): " << jFexSRJEtValues->size() << "\n";
            if(seedObjectType == "jFEXSRJ"){
                if(jFexSRJEtValues->size() < nSeedsInput_){
                    if(iSeed >= jFexSRJEtValues->size()){
                        seedValuesOriginal[iSeed].et = 0;
                        seedValuesOriginal[iSeed].eta = 0;
                        seedValuesOriginal[iSeed].phi = 0;
                        continue;
                    }
                }
                seedValuesOriginal[iSeed].et = digitize(jFexSRJEtValues->at(iSeed), et_bit_length_,
                                static_cast<double>(et_min_), static_cast<double>(et_max_));
                seedValuesOriginal[iSeed].eta = digitize(jFexSRJEtaValues->at(iSeed), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                seedValuesOriginal[iSeed].phi = digitize(jFexSRJPhiValues->at(iSeed), phi_bit_length_, phi_min_, phi_max_);
            }
            else if(seedObjectType == "gFEXSRJ"){
                if(gFexSRJEtValues->size() < nSeedsInput_){
                    if(iSeed >= gFexSRJEtValues->size()){
                        seedValuesOriginal[iSeed].et = 0;
                        seedValuesOriginal[iSeed].eta = 0;
                        seedValuesOriginal[iSeed].phi = 0;
                        continue;
                    }
                }
                seedValuesOriginal[iSeed].et = digitize(gFexSRJEtValues->at(iSeed), et_bit_length_,
                                static_cast<double>(et_min_), static_cast<double>(et_max_));
                seedValuesOriginal[iSeed].eta = digitize(gFexSRJEtaValues->at(iSeed), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                seedValuesOriginal[iSeed].phi = digitize(gFexSRJPhiValues->at(iSeed), phi_bit_length_, phi_min_, phi_max_);
            }
            else if(seedObjectType == "gepWTAConeCellsTowersJets"){
                if(useEtaSKObjects){
                    if(gepWTAConeCellsTowersEtaSKJetsEtValues->size() < nSeedsInput_){
                        if(iSeed >= gepWTAConeCellsTowersEtaSKJetsEtValues->size()){
                            seedValuesOriginal[iSeed].et = 0;
                            seedValuesOriginal[iSeed].eta = 0;
                            seedValuesOriginal[iSeed].phi = 0;
                            continue;
                        }
                    }
                    seedValuesOriginal[iSeed].et = digitize(gepWTAConeCellsTowersEtaSKJetsEtValues->at(iSeed), et_bit_length_,
                                    static_cast<double>(et_min_), static_cast<double>(et_max_));
                    seedValuesOriginal[iSeed].eta = digitize(gepWTAConeCellsTowersEtaSKJetsEtaValues->at(iSeed), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    seedValuesOriginal[iSeed].phi = digitize(gepWTAConeCellsTowersEtaSKJetsPhiValues->at(iSeed), phi_bit_length_, phi_min_, phi_max_);
                }
                else if(useSKObjects){
                    //std::cout << " seedobject type is correct" << "\n";
                    //std::cout << "gepWTAConeCellsTowersSKJetsEtValues->size() :" << gepWTAConeCellsTowersSKJetsEtValues->size() << "\n";
                    if(gepWTAConeCellsTowersSKJetsEtValues->size() < nSeedsInput_){
                        if(iSeed >= gepWTAConeCellsTowersSKJetsEtValues->size()){
                            seedValuesOriginal[iSeed].et = 0;
                            seedValuesOriginal[iSeed].eta = 0;
                            seedValuesOriginal[iSeed].phi = 0;
                            continue;
                        }
                    }
                    seedValuesOriginal[iSeed].et = digitize(gepWTAConeCellsTowersSKJetsEtValues->at(iSeed), et_bit_length_,
                                    static_cast<double>(et_min_), static_cast<double>(et_max_));
                    seedValuesOriginal[iSeed].eta = digitize(gepWTAConeCellsTowersSKJetsEtaValues->at(iSeed), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    seedValuesOriginal[iSeed].phi = digitize(gepWTAConeCellsTowersSKJetsPhiValues->at(iSeed), phi_bit_length_, phi_min_, phi_max_);
                }
                else{
                    //std::cout << " seedobject type is correct" << "\n";
                    //std::cout << "gepConeGEPBasicClustersJetsEtValues->size() :" << gepWTAConeCellsTowersJetsEtValues->size() << "\n";
                    if(gepWTAConeCellsTowersJetsEtValues->size() < nSeedsInput_){
                        if(iSeed >= gepWTAConeCellsTowersJetsEtValues->size()){
                            seedValuesOriginal[iSeed].et = 0;
                            seedValuesOriginal[iSeed].eta = 0;
                            seedValuesOriginal[iSeed].phi = 0;
                            continue;
                        }
                    }
                    seedValuesOriginal[iSeed].et = digitize(gepWTAConeCellsTowersJetsEtValues->at(iSeed), et_bit_length_,
                                    static_cast<double>(et_min_), static_cast<double>(et_max_));
                    seedValuesOriginal[iSeed].eta = digitize(gepWTAConeCellsTowersJetsEtaValues->at(iSeed), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    seedValuesOriginal[iSeed].phi = digitize(gepWTAConeCellsTowersJetsPhiValues->at(iSeed), phi_bit_length_, phi_min_, phi_max_);
                }
            }
        }
        std::memcpy(seedValues, seedValuesOriginal,
            sizeof(seedValuesOriginal));
        // Copy seed values prior to doing pre-proceessing (overlap removal & seed optimization)
        if(algoVersion_ != 2){ 
            // FIXME add a boolean to enable/disable overlap removal
            // Perform overlap removal (OR) ensuring that leading, subleading seeds don't overlap within deltaR < 2.0 
            //std::cout << "seed 1 OR eta: " << seedValues[0].eta << " , phi: " << seedValues[0].phi << "\n";
            //std::cout << "seed 2 OR eta: " << seedValues[1].eta << " , phi: " << seedValues[1].phi << "\n";
            unsigned int deltaR2LeadingSubleading = digitizedDeltaR2(seedValues[0].eta, seedValues[0].phi, seedValues[1].eta, seedValues[1].phi);
            //std::cout << "deltaR2leadingsubleading: " << deltaR2LeadingSubleading << "\n";
            //std::cout << "deltaR^2 cut: " << 2 * 2 * digitized_delta_R2Cut_ << "\n";
            if(deltaR2LeadingSubleading <= 2 * 2 * digitized_delta_R2Cut_){
                //std::cout << "OR triggered" << "\n";
                for(unsigned int iSeedOR = nSeedsOutput_; iSeedOR < nSeedsInput_; iSeedOR++){
                    //std::cout << "iSeedOR: " << iSeedOR << "\n";
                    if(seedValues[iSeedOR].et == 0 && seedValues[iSeedOR].eta == 0 && seedValues[iSeedOR].phi == 0) continue; // don't consider if et, eta, phi all = 0 (no jet)
                    unsigned int deltaR2LeadingSubleadingNthLeading = digitizedDeltaR2(seedValues[0].eta, seedValues[0].phi, seedValues[iSeedOR].eta, seedValues[iSeedOR].phi);
                    //std::cout << " deltaR2LeadingSubleadingNthLeading: " << deltaR2LeadingSubleadingNthLeading << "\n";
                    if(deltaR2LeadingSubleadingNthLeading > 2 * 2 * digitized_delta_R2Cut_){
                        //std::cout << "triggering for iSeedOR : "<< iSeedOR << "\n";
                        // Swap the entire (Et, eta, phi) triplet for original subleading, new subleading seed
                        std::swap(seedValues[1].et, seedValues[iSeedOR].et); // Et
                        std::swap(seedValues[1].eta, seedValues[iSeedOR].eta); // eta
                        std::swap(seedValues[1].phi, seedValues[iSeedOR].phi); // phi
                        overlapRemovalCounter++;
                        break; // break out of loop as we've found something to swap, thus preventing overlapping large-R jets
                    } // If deltaR^2 between original leading, other proto-seed farther than 2 times jet radius 
                } // Loop through potential additional seeds to find 
            } // If deltaR^2 between original leading, subleading closer than 2 times jet radius 
        } // if not using basic algorithm
        
        unsigned int objectsProcessed = maxObjectsConsidered_;
        if(inputObjectType == "gepBasicClusters"){
            if(useSKObjects){
                if (maxObjectsConsidered_ > gepBasicClustersSKEtValues->size()) objectsProcessed = gepBasicClustersSKEtValues->size();
            }
            else{
                if (maxObjectsConsidered_ > gepBasicClustersEtValues->size()) objectsProcessed = gepBasicClustersEtValues->size();
            }
        }
        else if(inputObjectType == "gepCellsTowers"){
            if(useEtaSKObjects){
                if (maxObjectsConsidered_ > gepCellsTowersEtaSKEtValues->size()) objectsProcessed = gepCellsTowersEtaSKEtValues->size();
            }
            else if(useSKObjects){
                if (maxObjectsConsidered_ > gepCellsTowersSKEtValues->size()) objectsProcessed = gepCellsTowersSKEtValues->size();
            }
            else{
                if (maxObjectsConsidered_ > gepCellsTowersEtValues->size()) objectsProcessed = gepCellsTowersEtValues->size();
            }
        }
        else if(inputObjectType == "gepWTAConeCellsTowersJets"){
            if(useEtaSKObjects){
                if (maxObjectsConsidered_ + nSeedsOutput_ >= gepWTAConeCellsTowersEtaSKJetsEtValues->size()){
                    if(gepWTAConeCellsTowersEtaSKJetsEtValues->size() >= nSeedsOutput_ ) objectsProcessed = gepWTAConeCellsTowersEtaSKJetsEtValues->size() - nSeedsOutput_;
                    else objectsProcessed = 0;
                }
            }
            else if(useSKObjects){
                if (maxObjectsConsidered_ + nSeedsOutput_ >= gepWTAConeCellsTowersSKJetsEtValues->size()){
                    if(gepWTAConeCellsTowersSKJetsEtValues->size() >= nSeedsOutput_ ) objectsProcessed = gepWTAConeCellsTowersSKJetsEtValues->size() - nSeedsOutput_;
                    else objectsProcessed = 0;
                }
            }
            else{
                if (maxObjectsConsidered_ + nSeedsOutput_ >= gepWTAConeCellsTowersJetsEtValues->size()){
                    if(gepWTAConeCellsTowersJetsEtValues->size() >= nSeedsOutput_ ) objectsProcessed = gepWTAConeCellsTowersJetsEtValues->size() - nSeedsOutput_;
                    else objectsProcessed = 0;
                }
            }
        }
        //std::cout << "--------------------------------------------" << "\n";
        //std::cout << "OBJECTS PROCESSED: " << objectsProcessed << "\n";
        for(unsigned int iIO = 0; iIO < objectsProcessed; iIO++){
            if(inputObjectType == "gepBasicClusters"){
                if(useSKObjects){
                    inputObjectValues[iIO].et = digitize(gepBasicClustersSKEtValues->at(iIO), et_bit_length_,
                              static_cast<double>(et_min_), static_cast<double>(et_max_));
                    inputObjectValues[iIO].eta = digitize(gepBasicClustersSKEtaValues->at(iIO), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    inputObjectValues[iIO].phi = digitize(gepBasicClustersSKPhiValues->at(iIO), phi_bit_length_, phi_min_, phi_max_);
                }
                else{
                    inputObjectValues[iIO].et = digitize(gepBasicClustersEtValues->at(iIO), et_bit_length_,
                              static_cast<double>(et_min_), static_cast<double>(et_max_));
                    inputObjectValues[iIO].eta = digitize(gepBasicClustersEtaValues->at(iIO), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    inputObjectValues[iIO].phi = digitize(gepBasicClustersPhiValues->at(iIO), phi_bit_length_, phi_min_, phi_max_);
                }
            }
            else if(inputObjectType == "gepCellsTowers"){
                if(useEtaSKObjects){
                    inputObjectValues[iIO].et = digitize(gepCellsTowersEtaSKEtValues->at(iIO), et_bit_length_,
                                static_cast<double>(et_min_), static_cast<double>(et_max_));
                    inputObjectValues[iIO].eta = digitize(gepCellsTowersEtaSKEtaValues->at(iIO), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    inputObjectValues[iIO].phi = digitize(gepCellsTowersEtaSKPhiValues->at(iIO), phi_bit_length_, phi_min_, phi_max_);
                }
                else if(useSKObjects){
                    inputObjectValues[iIO].et = digitize(gepCellsTowersSKEtValues->at(iIO), et_bit_length_,
                                static_cast<double>(et_min_), static_cast<double>(et_max_));
                    inputObjectValues[iIO].eta = digitize(gepCellsTowersSKEtaValues->at(iIO), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    inputObjectValues[iIO].phi = digitize(gepCellsTowersSKPhiValues->at(iIO), phi_bit_length_, phi_min_, phi_max_);
                }
                else{
                    inputObjectValues[iIO].et = digitize(gepCellsTowersEtValues->at(iIO), et_bit_length_,
                                static_cast<double>(et_min_), static_cast<double>(et_max_));
                    inputObjectValues[iIO].eta = digitize(gepCellsTowersEtaValues->at(iIO), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    inputObjectValues[iIO].phi = digitize(gepCellsTowersPhiValues->at(iIO), phi_bit_length_, phi_min_, phi_max_);
                }
            }
            else if(inputObjectType == "gepWTAConeCellsTowersJets"){
                //std::cout << "digitizing gepWTAConeCellsTowersJets" << "\n";
                if(useEtaSKObjects){
                    inputObjectValues[iIO].et = digitize(gepWTAConeCellsTowersEtaSKJetsEtValues->at(iIO + nSeedsOutput_), et_bit_length_,
                                static_cast<double>(et_min_), static_cast<double>(et_max_));
                    inputObjectValues[iIO].eta = digitize(gepWTAConeCellsTowersEtaSKJetsEtaValues->at(iIO + nSeedsOutput_), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    inputObjectValues[iIO].phi = digitize(gepWTAConeCellsTowersEtaSKJetsPhiValues->at(iIO + nSeedsOutput_), phi_bit_length_, phi_min_, phi_max_);
                }
                else if(useSKObjects){
                    inputObjectValues[iIO].et = digitize(gepWTAConeCellsTowersSKJetsEtValues->at(iIO + nSeedsOutput_), et_bit_length_,
                                static_cast<double>(et_min_), static_cast<double>(et_max_));
                    inputObjectValues[iIO].eta = digitize(gepWTAConeCellsTowersSKJetsEtaValues->at(iIO + nSeedsOutput_), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    inputObjectValues[iIO].phi = digitize(gepWTAConeCellsTowersSKJetsPhiValues->at(iIO + nSeedsOutput_), phi_bit_length_, phi_min_, phi_max_);
                }
                else{
                    inputObjectValues[iIO].et = digitize(gepWTAConeCellsTowersJetsEtValues->at(iIO + nSeedsOutput_), et_bit_length_,
                                static_cast<double>(et_min_), static_cast<double>(et_max_));
                    inputObjectValues[iIO].eta = digitize(gepWTAConeCellsTowersJetsEtaValues->at(iIO + nSeedsOutput_), eta_bit_length_, eta_min_, eta_max_, eta_range_);
                    inputObjectValues[iIO].phi = digitize(gepWTAConeCellsTowersJetsPhiValues->at(iIO + nSeedsOutput_), phi_bit_length_, phi_min_, phi_max_);
                }
            }

        }

        // --------------------------------------------
        // --- Pasted in and modified HLS algorithm ---
        // --------------------------------------------
        // Seed position optimization 
        if(algoVersion_ != 2){ // Only execute if using basic (v2) algorithm
            //std::cout << "PERFORMING SEED POSITION OPTIMIZATION" << "\n";
            unsigned int deltaRValuesLessThanSearchDistanceCounter[nSeedsOutput_] = {0, 0};
            int indices[nSeedsOutput_] = {-1, -1};
            bool indicesofProtoSeeds[nSeedsOutput_][nProtoSeeds_ - nSeedsOutput_] = {false};
            //std::cout << "indices[0]: " << indices[0] << " , indices[1]: " << indices[1] << "\n";
            //std::cout << " ----------------- SEED POS OPT -----------------" << "\n";
            for (unsigned int iSeed = 0; iSeed < nSeedsOutput_; iSeed++){ // loop through and calculate deltaR between leading, subleading & 3rd - 6th highest Et JFEX SRJ
                for (unsigned int iPreSeed = 0; iPreSeed < (nProtoSeeds_ - nSeedsOutput_); iPreSeed++){ // Loop through the number of considered
                    //std::cout << "iPreSeed: " << iPreSeed << "\n";
                    //std::cout << "-----------------------------------" << "\n";
                    //std::cout << "iSeed (seed opt): " << iSeed << "\n";
                    //std::cout << "iPreSeed (seed opt): " << iPreSeed << "\n";
                    //std::cout << "seed eta: " << seedValues[iSeed].eta << " , seed phi: " << seedValues[iSeed].phi << "\n";
                    //std::cout << "iPreSeed eta: " << seedValues[iPreSeed + nSeedsOutput_].eta << " , iPreSeed phi: " << seedValues[iPreSeed + nSeedsOutput_].phi << "\n";
                    
                    // Don't consider zero Et pre-seeds
                    if(seedValues[iPreSeed + nSeedsOutput_].et == 0) continue; 

                    // Don't consider certain preSeeds for seed position recalculation if below some mininmum E_T threshold
                    if(minEtSeedPosOptimization){
                        if(seedValues[iPreSeed + nSeedsOutput_].et <=  minEtSeedPosOptimizationCut / et_granularity_) continue;
                    }

                    //std::cout << "uDeltaPhi phi: " << uDeltaPhi << "\n";
                    //std::cout << "uDeltaPhi post wrap: " << uDeltaPhi << "\n";
                    unsigned int deltaR2SeedPreSeed = digitizedDeltaR2(seedValues[iSeed].eta, seedValues[iSeed].phi, seedValues[iPreSeed + nSeedsOutput_].eta, seedValues[iPreSeed + nSeedsOutput_].phi);
                    //std::cout << "deltaR2SeedPreSeed: " << deltaR2SeedPreSeed * 0.01 << "\n";
                    //std::cout << "deltaR2SeedPreSeed: " << deltaR2SeedPreSeed << " , digitized_d_search_squared_: " << digitized_d_search_squared_ << "\n";
                    if (deltaR2SeedPreSeed <= digitized_d_search_squared_){
                        //std::cout << "seed pos opt triggered" << "\n";
                        deltaRValuesLessThanSearchDistanceCounter[iSeed]++;
                        indicesofProtoSeeds[iSeed][iPreSeed] = true;
                        indices[iSeed] = iPreSeed;
                        //std::cout << "INDICES[iSEED] BEING SET TO : " << iPreSeed << "\n";
                    }
                }
            }
            //std::cout << "indices[0]: " << indices[0] << " , indices[1]: " << indices[1] << "\n";
            //std::cout << "pi_digitized_in_phi_: " << pi_digitized_in_phi_ << " , (2 * pi_digitized_in_phi_):  " << (2 * pi_digitized_in_phi_) << "\n";
            // Account for the case when multiple seeds are within deltaR customalizable value - use highest Et seed to compute halfway point
            // For seed 0 (leading seed)
            if (deltaRValuesLessThanSearchDistanceCounter[0] > 1) {
                //std::cout << "seed 0 why is this triggered" << "\n";
                unsigned int maxEt = 0;
                for (unsigned int iPreSeed = 0; iPreSeed < (nProtoSeeds_ - nSeedsOutput_); iPreSeed++) {
                    if (indicesofProtoSeeds[0][iPreSeed]) {
                        unsigned int et = seedValues[iPreSeed + nSeedsOutput_].et;
                        if (et > maxEt) {
                            maxEt = et;
                            indices[0] = iPreSeed;
                            //std::cout << "final ipreseed for jet 0: " << indices[0] << "\n";
                        }
                    }
                }
            } // Note: for case where only one potential pre-seed is found, index should already be set.

            // For seed 1
            if (deltaRValuesLessThanSearchDistanceCounter[1] > 1) {
                //std::cout << "seed 1 why is this triggered" << "\n";
                unsigned int maxEt = 0;
                for (unsigned int iPreSeed = 0; iPreSeed < (nProtoSeeds_ - nSeedsOutput_); iPreSeed++) {
                    if (indicesofProtoSeeds[1][iPreSeed]) {
                        unsigned int et = seedValues[iPreSeed + nSeedsOutput_].et;
                        if (et > maxEt) {
                            maxEt = et;
                            indices[1] = iPreSeed;
                            //std::cout << "final ipreseed for jet 1: " << indices[1] << "\n";
                        }
                    }
                }
            }

            // FIXME replace with a templated version of this s.t. you don't need to rewrite code for different number sorted (of index_of_min)
            
            // Account for the case when both original seeds would be shifted toward the same pre seeds, 
            // this would cause unnecessary overlap - this should be better thought out, 
            // but for now just prevent lower E_T seed from being shifted at all
            //std::cout << "indices[0]: " << indices[0] << "\n";
            //std::cout << "indices[1]: " << indices[1] << "\n";
            bool skipSecondSeed = false;
            if(indices[0] == indices[1] && indices[0] != -1){ 
                skipSecondSeed = true;
            }
            //std::cout << "skipSecondSeed: " << skipSecondSeed << "\n";
            // next update leading, subleading seed positions (and energy?) to be in between closest other seed 
            for (unsigned int iSeed = 0; iSeed < nSeedsOutput_; iSeed++){
                //std::cout << "iSeed: " << iSeed << " , skipSecondSeed: " << skipSecondSeed << "\n";
                //std::cout << "indices[iSeed]: " << indices[iSeed] << "\n";
                if(skipSecondSeed == true && iSeed == 1) continue; // case where don't want to shift 2nd seed since this would overlap with 1st seed
                if(indices[iSeed] == -1) continue; // if didn't find proto seed within search distance
                //std::cout << "iSeed MIDPOINT CALC: " << iSeed << "\n";
                // Get Et values for the two seeds (assumed already extracted)
                int et1 = seedValues[iSeed].et;
                int et2 = seedValues[indices[iSeed] + nSeedsOutput_].et;
                int etSum = et1 + et2; 

                // Convert eta to signed integers centered at 0
                int eta1 = seedValues[iSeed].eta - (1 << (eta_bit_length_ - 1));
                int eta2 = seedValues[indices[iSeed] + nSeedsOutput_].eta - (1 << (eta_bit_length_ - 1));

                // Convert phi to signed integers centered at 0
                int phi1s = seedValues[iSeed].phi - (1 << (phi_bit_length_ - 1));
                int phi2s = seedValues[indices[iSeed] + nSeedsOutput_].phi - (1 << (phi_bit_length_ - 1));
                //std::cout << "phi1s: " << phi1s << " , phi2s: " << phi2s << "\n";
                // --- Shortest-arc midpoint (unweighted) ---
                int dphi = phi2s - phi1s;

                dphi = wrapSym(dphi);                          // now in [-PI, PI)

                //std::cout << "PI_D: " << PI_D << "\n";
                //std::cout << "TWO_PI_D: " << TWO_PI_D << "\n";

                int phi_mid;
            
                // --- Compute midpoint in signed domain ---
                int eta_mid;
                if(!enableEtWeightedMidpoint || etSum == 0){ // Must account for case when there aren't any cone jets (or seeds in general) --> etSum == 0
                    eta_mid = (eta1 + eta2) >> 1; // Integer division by 2
                    
                    phi_mid = phi1s + (dphi >> 1);  // arithmetic shift (divide by 2)
                    //std::cout << "phi_mid before wrap: " << phi_mid << "\n";
                    phi_mid = wrapSym(phi_mid);
                    //std::cout << "phi_mid after wrap: " << phi_mid << "\n";
                } 
                else{
                    eta_mid = (et2 * eta2 + et1 * eta1) / etSum;
                    phi_mid = phi1s + ( et2 * dphi ) / etSum;
                    phi_mid = wrapSym(phi_mid);
                }

                // --- Convert midpoints back to digitized unsigned format ---
                unsigned int eta_mid_digitized = eta_mid + (1 << (eta_bit_length_ - 1));
                unsigned int phi_mid_digitized = phi_mid + (1 << (phi_bit_length_ - 1));
                //std::cout << "phi_mid_digitized: " << phi_mid_digitized << "\n";

                seedValues[iSeed].eta = eta_mid_digitized; // FIXME replace 1, 2 with etaIndexInput, phiIndexInput - fine for now
                seedValues[iSeed].phi = phi_mid_digitized;
            }
        }
        
        inputObject subjets[(1 << num_subjets_length_) - 1]; // Store the Et, Eta, Phi of subjets for use in substructure calculations (massApprox, tau_2, tau_1)
        inputObject subjetsSeedArray[nSeedsOutput_][(1 << num_subjets_length_) - 1] = {};
        unsigned int numSubjetsSeedArray[nSeedsOutput_] = {};
        std::vector<unsigned int > jet1MergedIndices; // store indices of merged input objects for each seed
        std::vector<unsigned int > jet2MergedIndices;

        // Combinations of seeds and input objects for merging 
        for (unsigned int iSeed = 0; iSeed < nSeedsOutput_; ++iSeed){
            std::vector<unsigned int > mergedInputObjectIndices; 
            //std::cout << "-----------------------------" << "\n";
            //std::cout << "iSeed: " << iSeed << "\n";
            //std::cout << "original seed eta: " << seedValuesOriginal[iSeed].eta << " , phi: " << seedValuesOriginal[iSeed].phi << " , et: " << seedValuesOriginal[iSeed].et << "\n";
            //std::cout << "seed eta: " << seedValues[iSeed].eta << " , seedphi: " << seedValues[iSeed].phi << " , et: " << seedValues[iSeed].et << "\n";
            unsigned int outputJetEt = 0;
            if(inputObjectType == "gepWTAConeCellsTowersJets") outputJetEt = seedValues[iSeed].et; // If clustering using jets, output jet E_T should start with seed E_T
            unsigned int numSubjets = 0; 
            unsigned int jet_psi_R = 0;
            unsigned int massApprox = 0;
            unsigned int tau_2 = 0;
            unsigned int tau_1 = 0;

            // Only do the main processing for non-zero seeds - otherwise end up reclustering input objects around 0,0 in eta-phi, and get non-zero output Et, other properties
            if(seedValues[iSeed].et != 0){
                for (unsigned int iInput = 0; iInput < objectsProcessed; ++iInput){ // loop through input objects to consider merging
                    //std::cout << "inputObjectValues[iInput].et: " << undigitize_et(inputObjectValues[iInput].et) << "\n";
                    if(inputObjectValues[iInput].et  == 0) break;
                    //std::cout << "iInput: " << iInput << "\n";
                    //std::cout << "input eta: " << inputObjectValues[iInput].eta << " , phi: " << inputObjectValues[iInput].phi << " , et: " << inputObjectValues[iInput].et << "\n";

                    unsigned int uDeltaEta = static_cast<unsigned int>(std::abs(int(seedValues[iSeed].eta) - int(inputObjectValues[iInput].eta)));
                    unsigned int uDeltaPhi = static_cast<unsigned int>(std::abs(int(seedValues[iSeed].phi) - int(inputObjectValues[iInput].phi)));

                    if (uDeltaPhi >= pi_digitized_in_phi_) uDeltaPhi = (2 * pi_digitized_in_phi_) - uDeltaPhi;
                    //std::cout << "deltaR2: " << uDeltaEta*uDeltaEta+uDeltaPhi*uDeltaPhi << " , digitized_delta_R2Cut_: " << digitized_delta_R2Cut_ << "\n";
                    if (uDeltaEta * uDeltaEta + uDeltaPhi * uDeltaPhi <= digitized_delta_R2Cut_){

                        if(outputJetEt + inputObjectValues[iInput].et >= int(et_max_/et_granularity_ - 1)) outputJetEt = int(et_max_/et_granularity_ - 1);
                        else outputJetEt += inputObjectValues[iInput].et; // add input object Et to seed Et for resultant output jet Et
                        //std::cout << "outputJetEt after sum: " << outputJetEt << "\n";
                        mergedInputObjectIndices.push_back(iInput);
                        // Compute psi_R using 8b 0,1 deltaR LUT - need to use a LUT since need deltaR rather than deltaR^2, would take too much latency to do sqrt computation directly!
                        unsigned int lutIndex = calcLutIndex(uDeltaEta, uDeltaPhi);
                        //std::cout << "inputObjectValues[iInput].et (psi_R): " << inputObjectValues[iInput].et << " , lutR_8b_[lutIndex]: " << lutR_8b_[lutIndex] << "\n";
                        jet_psi_R += inputObjectValues[iInput].et * lutR_8b_[lutIndex];

                        jetTaggerLRJMergedIndicesValues->push_back(iInput);
                        if(iSeed == 0) jet1MergedIndices.push_back(iInput);
                        if(iSeed == 1) jet2MergedIndices.push_back(iInput);

                        //std::cout << "jet_psi_R: " << jet_psi_R << "\n";
                        if(inputObjectType == "gepWTAConeCellsTowersJets"){
                            if(inputObjectValues[iInput].et > subjetEtThreshold/et_granularity_){ // check whether currently considered small-R jet input object is subjet (this covers all originally read in jets)
                                if(numSubjets < ((1 << num_subjets_length_) - 1)){
                                    subjets[numSubjets] = inputObjectValues[iInput]; // Assign subjet values from original seedValues 
                                    numSubjets++; // Clamp number of subjets to maximum value (max = 3 = 2^2 - 1) 
                                } 
                                else numSubjets = ((1 << num_subjets_length_) - 1); // Don't store for beyond 3 subjets - note assume the subjet input were sorted in Et (correct assumption at least for cone jets)

                            } 
                        }
                    }
                }

                if(inputObjectType == "gepWTAConeCellsTowersJets"){
                    for(unsigned int iSubjet = 0; iSubjet < nSeedsOutput_; iSubjet++){ // check whether ORIGINAL seeds are subjets
                        if(seedValuesOriginal[iSubjet].et > subjetEtThreshold/et_granularity_){
                            if(digitizedDeltaR2(seedValues[iSeed].eta, seedValues[iSeed].phi, seedValuesOriginal[iSubjet].eta, seedValuesOriginal[iSubjet].phi) <= digitized_delta_R2Cut_){
                                if(numSubjets < ((1 << num_subjets_length_) - 1)){
                                    //std::cout << "found seed subjet" << "\n";
                                    subjets[numSubjets] = seedValuesOriginal[iSubjet]; // Assign subjet values from original seedValues
                                    numSubjets++; // Clamp number of subjets to maximum value (max = 3 = 2^2 - 1)
                                }
                                else numSubjets = ((1 << num_subjets_length_) - 1); // Don't store for beyond 3 subjets - note assume the subjet input were sorted in Et (correct assumption at least for cone jets)
                            }
                        }
                    }
                }

                if(inputObjectType == "gepBasicClusters" || inputObjectType == "gepCellsTowers"){
                    for (unsigned int iSubjet = 0; iSubjet < nSeedsInput_ ; iSubjet++){ 
                        //std::cout << "iSubjet: " << iSubjet << "\n";
                        //std::cout << "seedValuesOriginal[iSubjet].et: " << seedValuesOriginal[iSubjet].et << "\n";
                        //std::cout << "subjetEtThreshold/et_granularity_: " << subjetEtThreshold/et_granularity_ << "\n";
                        if(seedValuesOriginal[iSubjet].et > subjetEtThreshold/et_granularity_){
                            //std::cout << " passes et threshold" << "\n";
                            //std::cout << "seedValues[iSeed].eta: " << undigitize_eta(seedValues[iSeed].eta)  << " ,seedValues[iSeed].phi: " << undigitize_phi(seedValues[iSeed].phi) << "\n";
                            //std::cout << "seedValuesOriginal[iSubjet].eta: " << undigitize_eta(seedValuesOriginal[iSubjet].eta)  << " ,seedValuesOriginal[iSubjet].phi: " << undigitize_phi(seedValuesOriginal[iSubjet].phi) << "\n";
                            //std::cout << "deltaR2: " << digitizedDeltaR2(...) << " , delta_R2Cut_: " << digitized_delta_R2Cut_ << "\n";
                            if(digitizedDeltaR2(seedValues[iSeed].eta, seedValues[iSeed].phi, seedValuesOriginal[iSubjet].eta, seedValuesOriginal[iSubjet].phi) <= digitized_delta_R2Cut_){
                                //std::cout << "FOUND SUBJET" << "\n";
                                if(numSubjets < ((1 << num_subjets_length_) - 1)){
                                    subjets[numSubjets] = seedValuesOriginal[iSubjet]; // Assign subjet values from original seedValues
                                    numSubjets++; // Clamp number of subjets to maximum value (max = 3 = 2^2 - 1)
                                }
                                else numSubjets = ((1 << num_subjets_length_) - 1); // Don't store for beyond 3 subjets - note assume the subjet input were sorted in Et (correct assumption at least for cone jets)
                            }
                        } // If passes subjet Et threshold
                    } // Loop through potential (non-seed) subjets
                } // If input object type isn't a small-R jet (only WTA-cone jets considered)
                //std::cout << "psi_R raw: " << jet_psi_R << "\n";
                //std::cout << "outputJetEt: " << outputJetEt << "\n";
                if(outputJetEt > 0) jet_psi_R = jet_psi_R / outputJetEt;
                else jet_psi_R = 0; 
                //std::cout << "numSubjets: " << numSubjets << "\n";
                //if(iSeed == 0) std::cout << "num lead subjets: " << numSubjets << "\n";
                //if(iSeed == 1) std::cout << "num subl subjets: " << numSubjets << "\n";
                //std::cout << "psi_R normalized: " << jet_psi_R << "\n";
                //if(iSeed == 0) std::cout << "num leading subjets: " << numSubjets << "\n";
                //if(iSeed == 1) std::cout << "num subleading subjets: " << numSubjets << "\n";
                //std::cout << "-----------------------" << "\n";
                // Compute massApproximation using deltaR between subjets & their Et
                if(numSubjets == 1){
                    // Compute distance between input objects and leading, subleading subjets for N-subjetiness
                    for (auto iInput : mergedInputObjectIndices){
                        //std::cout << "iInput subjets == 1: " << iInput << "\n";
                        //std::cout << "subjets[0].eta : " << undigitize_eta(subjets[0].eta) << " ,subjets[0].phi: " << undigitize_phi(subjets[0].phi) << "\n";
                        //std::cout << "inputObjectValues[iInput].eta : " << undigitize_eta(inputObjectValues[iInput].eta) << " , inputObjectValues[iInput].phi: " << undigitize_phi(inputObjectValues[iInput].phi) << "\n";
                        unsigned int uDeltaEta = static_cast<unsigned int>(std::abs(int(subjets[0].eta) - int(inputObjectValues[iInput].eta)));
                        unsigned int uDeltaPhi = static_cast<unsigned int>(std::abs(int(subjets[0].phi) - int(inputObjectValues[iInput].phi)));
                        if (uDeltaPhi >= pi_digitized_in_phi_) uDeltaPhi = (2 * pi_digitized_in_phi_) - uDeltaPhi;
                        //std::cout << "deltaR: " << lutR_8b_[calcLutIndex(uDeltaEta,uDeltaPhi)] * deltaR_granularity_ << " , inputObjectValues[iInput].et: " << undigitize_et(inputObjectValues[iInput].et) << "\n";
                        tau_1 += lutR_8b_[calcLutIndex(uDeltaEta, uDeltaPhi)] * inputObjectValues[iInput].et;
                        //std::cout << "tau_1: " << tau_1 << "\n";
                    }
                }
                else if(numSubjets >= 2){ // Require there to be at least 2 subjets for this computation to work
                    unsigned int subjetTotalEt = subjets[0].et + subjets[1].et;
                    //std::cout << "subjets[0].et: " << subjets[0].et << " , subjets[1].et: " << subjets[1].et << "\n";
                    //std::cout << "COMPUTING MASS APPROX subjet total et: " << subjetTotalEt << "\n";
                    //std::cout << "subjets[0].eta: " << subjets[0].eta << " , subjets[1].eta: " << subjets[1].eta << "\n";
                    //std::cout << "subjets[0].phi: " << subjets[0].phi << " , subjets[1].phi: " << subjets[1].phi << "\n";
                    // Compute deltaR between subjets for mass approximation
                    unsigned int uDeltaEta_sub = static_cast<unsigned int>(std::abs(int(subjets[0].eta) - int(subjets[1].eta)));
                    unsigned int uDeltaPhi_sub = static_cast<unsigned int>(std::abs(int(subjets[0].phi) - int(subjets[1].phi)));
                    if (uDeltaPhi_sub >= pi_digitized_in_phi_) uDeltaPhi_sub = (2 * pi_digitized_in_phi_) - uDeltaPhi_sub;
                    //std::cout << "subjetTotalEt: " << subjetTotalEt << " , lutR_8b_[calcLutIndex(...)]: " << lutR_8b_[calcLutIndex(uDeltaEta_sub,uDeltaPhi_sub)] << "\n";
                    unsigned int massApproxRaw = subjetTotalEt * lutR_8b_[calcLutIndex(uDeltaEta_sub, uDeltaPhi_sub)];
                    //std::cout << "massApproxRaw: " << massApproxRaw << "\n";
                    unsigned int massApproxTmp = massApproxRaw / massApproxDivisor_;
                    //std::cout << "massApproxTmp: " << massApproxTmp << "\n";
                    massApprox = (massApproxTmp > massApprox_max_) ? massApprox_max_ : massApproxTmp;
                    //std::cout << "massApprox: " << massApprox << "\n";

                    // Compute N-subjetiness: min distance from each constituent to the two subjets
                    for (auto iInput : mergedInputObjectIndices){
                        //std::cout << "iInput subjets >= 2: " << iInput << "\n";
                        //std::cout << "input eta; " << undigitize_eta(inputObjectValues[iInput].eta) << " , input phi : " << undigitize_phi(inputObjectValues[iInput].phi) << "\n";
                        unsigned int uDeltaEta_Subjet0 = static_cast<unsigned int>(std::abs(int(subjets[0].eta) - int(inputObjectValues[iInput].eta)));
                        unsigned int uDeltaPhi_Subjet0 = static_cast<unsigned int>(std::abs(int(subjets[0].phi) - int(inputObjectValues[iInput].phi)));
                        if (uDeltaPhi_Subjet0 >= pi_digitized_in_phi_) uDeltaPhi_Subjet0 = (2 * pi_digitized_in_phi_) - uDeltaPhi_Subjet0;
                        unsigned int deltaR_Subjet0 = lutR_8b_[calcLutIndex(uDeltaEta_Subjet0, uDeltaPhi_Subjet0)];
                        //std::cout << "subjet 0 eta; " << undigitize_eta(subjets[0].eta) << " , subjet 0 phi : " << undigitize_phi(subjets[0].phi) << "\n";
                        //std::cout << "deltaR_Subjet0: " << deltaR_Subjet0 * deltaR_granularity_ << "\n";
                        //std::cout << "inputObjectValues[iInput].et: " << undigitize_et(inputObjectValues[iInput].et) << "\n";
                        tau_1 += deltaR_Subjet0 * inputObjectValues[iInput].et;
                        //std::cout << "tau1 2 subjets: " << tau_1 << "\n";
                        unsigned int uDeltaEta_Subjet1 = static_cast<unsigned int>(std::abs(int(subjets[1].eta) - int(inputObjectValues[iInput].eta)));
                        unsigned int uDeltaPhi_Subjet1 = static_cast<unsigned int>(std::abs(int(subjets[1].phi) - int(inputObjectValues[iInput].phi)));
                        if (uDeltaPhi_Subjet1 >= pi_digitized_in_phi_) uDeltaPhi_Subjet1 = (2 * pi_digitized_in_phi_) - uDeltaPhi_Subjet1;
                        //std::cout << "subjet 1 eta; " << undigitize_eta(subjets[1].eta) << " , subjet 1 phi : " << undigitize_phi(subjets[1].phi) << "\n";
                        unsigned int deltaR_Subjet1 = lutR_8b_[calcLutIndex(uDeltaEta_Subjet1, uDeltaPhi_Subjet1)];
                        //std::cout << "deltaR_Subjet1: " << deltaR_Subjet1 * deltaR_granularity_ << "\n";
                        unsigned int min_deltaR = std::min(deltaR_Subjet0, deltaR_Subjet1);
                        //std::cout << "tau 2 min deltaR: " << min_deltaR * deltaR_granularity_  << "\n";
                        //std::cout << "tau 2 inputObjectValues[iInput].et: " << undigitize_et(inputObjectValues[iInput].et) << "\n";
                        //std::cout << "tau 2 contribution: " << min_deltaR * deltaR_granularity_ * undigitize_et(inputObjectValues[iInput].et) << "\n";
                        tau_2 += min_deltaR * inputObjectValues[iInput].et;
                        //std::cout << "tau_2: " << tau_2 << "\n";
                    }
                }
                if(outputJetEt > 0){
                    //std::cout << "tau_1 raw: " << tau_1 << "\n";
                    //std::cout << "normalization (Et * R jet): " << outputJetEt * rCut_ << "\n";
                    tau_1 = (tau_1 * ((1 << N_subjetiness_bit_length_) - 1)) / (outputJetEt * static_cast<unsigned int>(rCut_/(deltaR_granularity_) + 0.5)); // convert to integer between 0, 255 rather than double between 0, 1
                    //std::cout << "tau_1_after : " << tau_1 << "\n";
                    //std::cout << "tau_1 undigitized final: " << undigitize_N_subjetiness(tau_1) << "\n";
                    //if(undigitize_N_subjetiness(tau_1) > 1.0) std::cout << "TAU1 > 1" << "\n";
                    //std::cout << "tau_2 raw: " << tau_2 << "\n";
                    tau_2 = (tau_2 * ((1 << N_subjetiness_bit_length_) - 1)) / (outputJetEt * static_cast<unsigned int>(rCut_/(deltaR_granularity_) + 0.5)); // FIXME make this normalization a constant - 
                    //if(undigitize_N_subjetiness(tau_2) > 1.0) std::cout << "TAU2 > 1" << "\n";
                    // NOTE we normalize by 2 * since that is the max distance between subjet & constituent
                    //std::cout << "tau_2 after: " << tau_2 << "\n";
                    //std::cout << "tau_2 undigitized final: " << undigitize_N_subjetiness(tau_2) << "\n";
                }
                else{
                    tau_1 = 0;
                    tau_2 = 0;
                }
            }
            

            double tau_21;
            if(tau_1 != 0) tau_21 = double(tau_2) / double(tau_1);
            else tau_21 = 0;

            // Write (digitized) values to output memory files & (undigitized) values to ntuple

            // Values already computed:
            uint64_t num_subjets_value = static_cast<uint32_t>(numSubjets)        & maskN(num_subjets_length_);
            uint64_t tau_1_value = static_cast<uint32_t>(tau_1)        & maskN(N_subjetiness_bit_length_);
            uint64_t tau_2_value = static_cast<uint32_t>(tau_2)        & maskN(N_subjetiness_bit_length_);
            uint64_t massApprox_value = static_cast<uint32_t>(massApprox)        & maskN(mass_approx_bit_length_);
            uint64_t psi_R_value = static_cast<uint32_t>(jet_psi_R)        & maskN(psi_R_bit_length_);
            //std::cout << "dim_value: " << diam_value << "\n";
            uint64_t et_value   = static_cast<uint32_t>(outputJetEt)      & maskN(et_bit_length_);
            uint64_t eta_value  = static_cast<uint32_t>(seedValues[iSeed].eta) & maskN(eta_bit_length_);
            uint64_t phi_value  = static_cast<uint32_t>(seedValues[iSeed].phi) & maskN(phi_bit_length_);

            // Set output values - FIXME these aren't actually used currently, rather the bitset arrays below are used for setting final values
            outputJetValues[iSeed].et = et_value;   // Et
            outputJetValues[iSeed].eta = eta_value; // eta
            outputJetValues[iSeed].phi = phi_value; // phi
            outputJetValues[iSeed].nSubjets = numSubjets; // subjets above Et threshold
            outputJetValues[iSeed].psiR = psi_R_value; // psi_R jet profile
            outputJetValues[iSeed].tau_1 = tau_1_value; // 1 subjetiness
            outputJetValues[iSeed].tau_2 = tau_2_value; // 2 subjetiness 
            outputJetValues[iSeed].massApprox = massApprox_value; // massApprox = deltaR subjets * Et subjets
            
            // Bitset strings (unchanged)
            std::bitset<et_bit_length_  > et_bitset(et_value);
            std::bitset<eta_bit_length_ > eta_bitset(eta_value);
            std::bitset<phi_bit_length_ > phi_bitset(phi_value);
            std::bitset<num_subjets_length_> num_subjets_bitset(num_subjets_value);
            std::bitset<N_subjetiness_bit_length_> tau_1_bitset(tau_1_value);
            std::bitset<N_subjetiness_bit_length_> tau_2_bitset(tau_2_value);
            std::bitset<mass_approx_bit_length_> massApprox_bitset(massApprox_value);
            std::bitset<psi_R_bit_length_> psi_R_bitset(psi_R_value); // Note: psi_R only written to ntuple
            
            et_bitset_array[iSeed] = et_bitset;
            eta_bitset_array[iSeed] = eta_bitset;
            phi_bitset_array[iSeed] = phi_bitset;
            num_subjets_bitset_array[iSeed] = num_subjets_bitset;
            psi_R_bitset_array[iSeed] = psi_R_bitset;
            massApprox_bitset_array[iSeed] = massApprox_bitset;
            tau_1_bitset_array[iSeed] = tau_1_bitset;
            tau_2_bitset_array[iSeed] = tau_2_bitset;
            tau_21_array[iSeed] = tau_21;

            // Pack into one word: [nsubjets | et | eta | phi]
            // NOTE: psi_R now saved in output ntuple, but not in output test vector - kept for ease of using analysis script
            uint64_t combined_value;
            if(algoVersion_ == 2){
                combined_value =
                    ((et_value   & maskN(et_bit_length_  )) << (eta_bit_length_ + phi_bit_length_)) |
                    ((eta_value  & maskN(eta_bit_length_ )) <<  phi_bit_length_) |
                    ( phi_value  & maskN(phi_bit_length_ ));
            }
            else if(algoVersion_ == 3){
                combined_value =
                    ((num_subjets_value & maskN(num_subjets_length_)) << (et_bit_length_ + eta_bit_length_ + phi_bit_length_)) |
                    ((et_value   & maskN(et_bit_length_  )) << (eta_bit_length_ + phi_bit_length_)) |
                    ((eta_value  & maskN(eta_bit_length_ )) <<  phi_bit_length_) |
                    ( phi_value  & maskN(phi_bit_length_ ));
            }

            // Hex string with width = total_bits_/4 (=8 for 32 bits)
            std::stringstream hex_stream;
            hex_stream << std::hex << std::nouppercase << std::setfill('0') << std::setw(total_bits_output_ / 4) << combined_value;
            std::string hexValue = hex_stream.str();

            // Final writing of output test vectors - remove num subjets for v2 algo
            if(algoVersion_ == 2){
                f_output << "0x" << std::hex << std::setw(2) << std::setfill('0') << iSeed
                    << " " << et_bitset.to_string()
                    << "|" << eta_bitset.to_string()  << "|" << phi_bitset.to_string()
                    << " 0x" << hexValue << std::endl;
            }
            else if(algoVersion_ == 3){
                f_output << "0x" << std::hex << std::setw(2) << std::setfill('0') << iSeed
                    << " "  << num_subjets_bitset.to_string() << "|" << et_bitset.to_string()
                    << "|"  << eta_bitset.to_string()  << "|" << phi_bitset.to_string()
                    << " 0x" << hexValue << std::endl;
            }
            

            // Write data to all LRJ branch
            jetTaggerLRJPsi_RValues->push_back(undigitize_psi_R(psi_R_bitset));
            jetTaggerLRJMassApproxValues->push_back(undigitize_mass_approx(massApprox_bitset));
            jetTaggerLRJTau_1Values->push_back(undigitize_N_subjetiness(tau_1_bitset));
            jetTaggerLRJTau_2Values->push_back(undigitize_N_subjetiness(tau_2_bitset));
            jetTaggerLRJTau_21Values->push_back(tau_21); // Already undigitized! 
            jetTaggerLRJSubjetMultiplicityValues->push_back((numSubjets));
            //std::cout << "undigitized psi_R: " << undigitize_psi_R(psi_R_bitset) << "\n";
            jetTaggerLRJEtValues->push_back(undigitize_et(et_bitset));
            jetTaggerLRJEtaValues->push_back(undigitize_eta(eta_bitset));
            jetTaggerLRJPhiValues->push_back(undigitize_phi(phi_bitset));

            numSubjetsSeedArray[iSeed] = numSubjets;
            for(unsigned int iSubj = 0; iSubj < numSubjets; iSubj++){
                subjetsSeedArray[iSeed][iSubj] = subjets[iSubj];
            }
        } // Loop through seeds 

        // Write data to leading, subleading jet branches sorted by Et
        // jet1MergedIndices/jet2MergedIndices always correspond to iSeed==0 and iSeed==1 respectively
        const unsigned int leadIdx = (outputJetValues[0].et >= outputJetValues[1].et) ? 0 : 1;
        const unsigned int sublIdx = 1 - leadIdx;
        const auto& leadMergedIndices = (leadIdx == 0) ? jet1MergedIndices : jet2MergedIndices;
        const auto& sublMergedIndices = (leadIdx == 0) ? jet2MergedIndices : jet1MergedIndices;

        //std::cout << "leading num subjets post sort: " << num_subjets_bitset_array[leadIdx].to_ulong() << "\n";
        //std::cout << "subleading num subjets post sort: " << num_subjets_bitset_array[sublIdx].to_ulong() << "\n";
        jetTaggerLeadingLRJPsi_RValues->push_back(undigitize_psi_R(psi_R_bitset_array[leadIdx]));
        jetTaggerLeadingLRJSubjetMultiplicityValues->push_back(num_subjets_bitset_array[leadIdx].to_ulong());
        jetTaggerLeadingLRJMassApproxValues->push_back(undigitize_mass_approx(massApprox_bitset_array[leadIdx]));
        jetTaggerLeadingLRJTau_1Values->push_back(undigitize_N_subjetiness(tau_1_bitset_array[leadIdx]));
        jetTaggerLeadingLRJTau_2Values->push_back(undigitize_N_subjetiness(tau_2_bitset_array[leadIdx]));
        jetTaggerLeadingLRJTau_21Values->push_back(tau_21_array[leadIdx]);
        //std::cout << "et_bitset_array[leadIdx]: " << et_bitset_array[leadIdx] << "\n";
        //std::cout << "undigitize_et(et_bitset_array[leadIdx]): " << undigitize_et(et_bitset_array[leadIdx]) << "\n";
        jetTaggerLeadingLRJEtValues->push_back(undigitize_et(et_bitset_array[leadIdx]));
        jetTaggerLeadingLRJEtaValues->push_back(undigitize_eta(eta_bitset_array[leadIdx]));
        jetTaggerLeadingLRJPhiValues->push_back(undigitize_phi(phi_bitset_array[leadIdx]));

        for(unsigned int iSubj = 0; iSubj < numSubjetsSeedArray[leadIdx]; iSubj++){
            jetTaggerLeadingLRJSubjetEtValues->push_back(undigitize_et(subjetsSeedArray[leadIdx][iSubj].et));
            jetTaggerLeadingLRJSubjetEtaValues->push_back(undigitize_eta(subjetsSeedArray[leadIdx][iSubj].eta));
            jetTaggerLeadingLRJSubjetPhiValues->push_back(undigitize_phi(subjetsSeedArray[leadIdx][iSubj].phi));
        }

        // necessary to ensure merged indices match to leading, subleading correctly
        for(auto idx : leadMergedIndices) jetTaggerLeadingLRJMergedIndicesValues->push_back(idx);
        for(auto idx : sublMergedIndices) jetTaggerSubleadingLRJMergedIndicesValues->push_back(idx);

        jetTaggerSubleadingLRJPsi_RValues->push_back(undigitize_psi_R(psi_R_bitset_array[sublIdx]));
        jetTaggerSubleadingLRJSubjetMultiplicityValues->push_back(num_subjets_bitset_array[sublIdx].to_ulong());
        jetTaggerSubleadingLRJMassApproxValues->push_back(undigitize_mass_approx(massApprox_bitset_array[sublIdx]));
        jetTaggerSubleadingLRJTau_1Values->push_back(undigitize_N_subjetiness(tau_1_bitset_array[sublIdx]));
        jetTaggerSubleadingLRJTau_2Values->push_back(undigitize_N_subjetiness(tau_2_bitset_array[sublIdx]));
        jetTaggerSubleadingLRJTau_21Values->push_back(tau_21_array[sublIdx]);
        jetTaggerSubleadingLRJEtValues->push_back(undigitize_et(et_bitset_array[sublIdx]));
        jetTaggerSubleadingLRJEtaValues->push_back(undigitize_eta(eta_bitset_array[sublIdx]));
        jetTaggerSubleadingLRJPhiValues->push_back(undigitize_phi(phi_bitset_array[sublIdx]));

        for(unsigned int iSubj = 0; iSubj < numSubjetsSeedArray[sublIdx]; iSubj++){
            jetTaggerSubleadingLRJSubjetEtValues->push_back(undigitize_et(subjetsSeedArray[sublIdx][iSubj].et));
            jetTaggerSubleadingLRJSubjetEtaValues->push_back(undigitize_eta(subjetsSeedArray[sublIdx][iSubj].eta));
            jetTaggerSubleadingLRJSubjetPhiValues->push_back(undigitize_phi(subjetsSeedArray[sublIdx][iSubj].phi));
        }
        // Fill branch to output ntuple with
        jetTaggerLRJs->Fill();
        jetTaggerLeadingLRJs->Fill();
        jetTaggerSubleadingLRJs->Fill();
        emulEventInfoTree->Fill();
    } // Loop through events
    std::cout << "number of overlaps removed: " << overlapRemovalCounter << " for " << eventsToProcess << " events." << "\n";
    // Close output files
    outputFile->cd();
    std::cout << "writing output file" << "\n";
    jetTaggerLRJs->Write("", TObject::kOverwrite);
    jetTaggerLeadingLRJs->Write("", TObject::kOverwrite);
    jetTaggerSubleadingLRJs->Write("", TObject::kOverwrite);
    emulEventInfoTree->Write("", TObject::kOverwrite);
    outputFile->Close();
    inputFile->Close();

} // Event loop function

// Function for processing a different algorithm configuration
// Use(compiling, running with ROOT): .L jetTaggerEmulation.cc, jetTaggerEmulation(0.001, 128, 2, 1.0, true, 3, false)
void jetTaggerEmulation(double rMergeCut, // Distance in r-phi plane to look for other jFEX SRJs (from leading, subleading jets)
                                        // to recalculate seed position by taking midpoint between the two jets
                                        // to disable recalculation of seed position, set to 0.001
                        unsigned int numberIOs, // number of input objects (clusters or towers) considered for merging to seeds
                        unsigned int nSeeds, // Number of seeds and consequently number of output jets
                        double RSquaredCut, // radius squared of output jets
                        bool signalBool,  // whether processing a signal or background process
                        bool condorBool, // whether running using condor batch job submission (requires change in filepaths)
                        bool useSKObjects, // Whether to use PU-suppressed (with SoftKiller) objects
                        std::string signalString, // Which signal sample being used (functionality for: VBF_hh_bbbb, ggF_hh_bbbb, ZvvHbb, ttbar_had, Zprime_ttbar)
                        std::string inputObjectType = "gepCellsTowers", // Possibilities: "gepCellsTowers", "gepWTAConeCellsTowersJets", "gepBasicClusters"
                        std::string seedObjectType = "gepWTAConeCellsTowersJets", // Possibilities: "gepWTAConeCellsTowersJets" or "jFEXSRJ"
                        double subjetEtThreshold = 25.0, // Nominally 25 GeV for gepWTAConeCellsTowersJets, 50 GeV for jFEX SRJ, 35 GeV for gFEX SRJ based on their respective energy scales, not fully optimized
                        bool enableEtWeightedMidpoint = false, // Whether use E_T weighted midpoint calculation for seed optimization
                        bool minEtSeedPosOptimization = false, // Whether require a minimum E_T to allow a proto-seed to be considered for seed optimization
                        double minEtSeedPosOptimizationCut = 25.0, // Minimum E_T (in GeV) if minEtSeedPosOptimization is enabled for a proto-seed to be considered for seed optimization
                        bool writeMemPrints = true, // Whether to write memory-print text files (disable for batch/Condor runs)
                        std::string explicitInputPath = "", // When non-empty, overrides makeInputFileName (used for per-file Condor parallelism)
                        int fileIndex = -1, // When >= 0, appended as _fileN to output name to avoid collisions across parallel jobs
                        bool useEtaSKObjects = false // Whether to use EtaSK PU-suppressed objects (gepCellsTowers and WTAConeJets only)
                        ){
    if(signalBool) std::cout << "Processing signal of: " << signalString  << "\n";
    // Construct input and output ntuple, LUT paths based on configuration type
    auto infile = explicitInputPath.empty() ? makeInputFileName(signalBool, signalString) : explicitInputPath;
    auto outntuplefile = makeOutputFileName(rMergeCut, numberIOs, nSeeds, RSquaredCut, signalBool, signalString, inputObjectType, seedObjectType, useSKObjects, algoVersion_, subjetEtThreshold, enableEtWeightedMidpoint, minEtSeedPosOptimization, minEtSeedPosOptimizationCut, "/data/larsonma/LargeRadiusJets/outputNTuplesDev_CondorSubmission_NewSamples/", useEtaSKObjects);
    if (fileIndex >= 0) {
        size_t pos = outntuplefile.rfind(".root");
        if (pos != std::string::npos)
            outntuplefile = outntuplefile.substr(0, pos) + "_file" + std::to_string(fileIndex) + ".root";
    }
    auto outtextfile = writeMemPrints
        ? makeOutputTextFileName(rMergeCut, numberIOs, nSeeds, RSquaredCut, signalBool, signalString, inputObjectType, seedObjectType, useSKObjects, algoVersion_, "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrintsEmulation/", useEtaSKObjects)
        : std::string("/dev/null");
    std::cout << "infile: " << infile << "\n";
    std::cout << "outntuplefile: " << outntuplefile << "\n";
    if (writeMemPrints) std::cout << "outtextfile: " << outtextfile << "\n";
    std::cout << "calling event loop: " << "\n";
    eventLoop(infile, outntuplefile, outtextfile, inputObjectType, seedObjectType, useSKObjects,
        enableEtWeightedMidpoint, minEtSeedPosOptimization,
        minEtSeedPosOptimizationCut, subjetEtThreshold, useEtaSKObjects);
    //gSystem->Exit(0);
}
