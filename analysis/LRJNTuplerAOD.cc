// Save as: nTupler.C
// To execute: e.g., root ; .L nTupler.C ; nTupler(true, true, true) 

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <filesystem> // C++17
#include <algorithm>
#include <stdexcept>  // for std::runtime_error

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH2F.h"
#include "analysisHelperFunctions.h"

#define PROCESS_DAOD true

// Settings
const bool afBool = true;
const bool vbfBool = true;


unsigned int digitize(double value, int bit_length, double min_val, double max_val) {
    // Check if value is in range
    if (value < min_val || value > max_val) {
        throw std::runtime_error("Value " + std::to_string(value) +
                                 " is out of range (" + std::to_string(min_val) +
                                 ", " + std::to_string(max_val) + ")");
    }

    double scale = (std::pow(2, bit_length) - 1) / (max_val - min_val);
    return static_cast<unsigned int>(std::round((value - min_val) * scale));
}

// Define phi wrapping
float delta_phi(float phi1, float phi2) {
    float dphi = phi1 - phi2;
    while (dphi > M_PI) dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    return dphi;
}

// Define delta R
float delta_R(float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1 - eta2;
    float dphi = delta_phi(phi1, phi2);
    return sqrt(deta*deta + dphi*dphi);
}

// Recursively find first non-Higgs daughters
void find_non_higgs_daughters(const xAOD::TruthParticle* particle,
                              std::vector<const xAOD::TruthParticle*>& result) {
    if (!particle) return;

    // If the particle is NOT a Higgs, store it directly
    if (particle->pdgId() != 25) {
        result.push_back(particle);
        return;
    }

    // If it is a Higgs, recurse through its children
    for (unsigned int i = 0; i < particle->nChildren(); ++i) {
        find_non_higgs_daughters(particle->child(i), result);
    }
}


// Main function
void nTupler(bool signalBool, bool daodBool, bool vbfBool) {
    // Setup file paths based on whether processing signal or background, and vbf production or ggF production
    string fileDir;
    if (vbfBool){
        if (daodBool){
            if(signalBool){
                fileDir = "/data/larsonma/LargeRadiusJets/datasets/Signal_HHbbbb_DAODAOD/DAOD/VBF/mc21_14TeV.537536.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv0cv1.deriv.DAOD_JETM42.e8557_s4422_r16130_p6658";
            }
            else{
                fileDir = "/data/larsonma/LargeRadiusJets/datasets/Background_jj_JZ3_100kEv/mc21_14TeV.801168.Py8EG_A14NNPDF23LO_jj_JZ3.deriv.DAOD_JETM42.e8557_s4422_r16130_p6658";
            }
            
        }
        else{
            if (signalBool) {
                fileDir = afBool ?
                    "/data/larsonma/LargeRadiusJets/datasets/Signal_HHbbb_AOD_10kEv/mc21_14TeV.537536.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv0cv1.recon.AOD.e8557_s4422_r16130" :
                    "/home/larsonma/GEPHadronicEventReconstruction/datasets/Signal_HHbbbb/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130";
            } else {
                fileDir = afBool ?
                    "/data/larsonma/LargeRadiusJets/datasets/Background_jj_JZ3_100kEv/mc21_14TeV.801168.Py8EG_A14NNPDF23LO_jj_JZ3.recon.AOD.e8557_s4422_r16130" :
                    "/home/larsonma/GEPHadronicEventReconstruction/datasets/Background_jj_JZ3/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130";
            }
        }
    }
    else{
        if (daodBool){
            if (signalBool){
                fileDir = "/data/larsonma/LargeRadiusJets/datasets/Signal_HHbbbb_DAODAOD/DAOD/ggF/mc21_14TeV.603277.PhPy8EG_PDF4LHC21_HHbbbb_HLLHC_chhh1p0.deriv.DAOD_JETM42.e8564_s4422_r16130_p6658";
            }
            else{
                fileDir = "/data/larsonma/LargeRadiusJets/datasets/Background_jj_JZ3/mc21_14TeV.801168.Py8EG_A14NNPDF23LO_jj_JZ3.deriv.DAOD_JETM42.e8557_s4422_r16130_p6658";
            }
        }
        else{
            if (signalBool) {
                fileDir = afBool ?
                    "/data/larsonma/LargeRadiusJets/datasets/Signal_HHbbbb/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130":
                    "/home/larsonma/GEPHadronicEventReconstruction/datasets/Signal_HHbbbb/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130";
            } else {
                fileDir = afBool ?
                    "/data/larsonma/LargeRadiusJets/datasets/Background_jj_JZ3/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130" :
                    "/home/larsonma/GEPHadronicEventReconstruction/datasets/Background_jj_JZ3/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130";
            }
        }
    }

    std::string output_file_topo422;
    std::string output_file_calotopotowers;
    std::string output_file_gfex;
    std::string output_file_jfex;

    if (signalBool) {
        if (vbfBool) {
            output_file_topo422 = "/data/larsonma/LargeRadiusJets/MemPrints/CaloTopo_422/mc21_14TeV_hh_bbbb_vbf_novhh_topo422.dat";
            output_file_calotopotowers = "/data/larsonma/LargeRadiusJets/MemPrints/CaloTopoTowers/mc21_14TeV_hh_bbbb_vbf_novhh_calotopotowers.dat";
            output_file_gfex = "/data/larsonma/LargeRadiusJets/MemPrints/gFex/mc21_14TeV_hh_bbbb_vbf_novhh_gfex_smallrj.dat";
            output_file_jfex = "/data/larsonma/LargeRadiusJets/MemPrints/jFex/mc21_14TeV_hh_bbbb_vbf_novhh_jfex_smallrj.dat";
        } else {
            output_file_topo422 = "/data/larsonma/LargeRadiusJets/MemPrints/CaloTopo_422/mc21_14TeV_HHbbbb_HLLHC_topo422.dat";
            output_file_calotopotowers = "/data/larsonma/LargeRadiusJets/MemPrints/CaloTopoTowers/mc21_14TeV_HHbbbb_HLLHC_calotopotowers.dat";
            output_file_gfex = "/data/larsonma/LargeRadiusJets/MemPrints/gFex/mc21_14TeV_HHbbbb_HLLHC_gfex_smallrj.dat";
            output_file_jfex = "/data/larsonma/LargeRadiusJets/MemPrints/jFex/mc21_14TeV_HHbbbb_HLLHC_jfex_smallrj.dat";
        }
    } else {
        output_file_topo422 = "/data/larsonma/LargeRadiusJets/MemPrints/CaloTopo_422/mc21_14TeV_jj_JZ3_topo422.dat";
        output_file_calotopotowers = "/data/larsonma/LargeRadiusJets/MemPrints/CaloTopoTowers/mc21_14TeV_jj_JZ3_calotopotowers.dat";
        output_file_gfex = "/data/larsonma/LargeRadiusJets/MemPrints/gFex/mc21_14TeV_jj_JZ3_gfex_smallrj.dat";
        output_file_jfex = "/data/larsonma/LargeRadiusJets/MemPrints/jFex/mc21_14TeV_jj_JZ3_jfex_smallrj.dat";
    }

    std::ofstream f_topotower(output_file_calotopotowers);
    std::ofstream f_topo(output_file_topo422);
    std::ofstream f_gfex(output_file_gfex);
    std::ofstream f_jfex(output_file_jfex);

    if (!f_topotower.is_open() || !f_topo.is_open() || !f_gfex.is_open() || !f_jfex.is_open()) {
        std::cerr << "Error: One or more output files could not be opened for writing!" << std::endl;
        // handle error or return
    }
    

    TString outputFileName; 
    if (vbfBool){
        if (daodBool){
            if (signalBool) outputFileName = "outputRootFiles/mc21_14TeV_hh_bbbb_vbf_novhh_DAOD_NTUPLE_HIGGSPTCUT.root";
            else outputFileName = "outputRootFiles/mc21_14TeV_jj_JZ3_DAOD_NTUPLE.root";
        }
        else{
            if (signalBool) outputFileName = "outputRootFiles/mc21_14TeV_hh_bbbb_vbf_novhh_AOD_NTUPLE.root";
            else outputFileName = "outputRootFiles/mc21_14TeV_jj_JZ3_AOD_NTUPLE_100kEv.root";
        }
    }
    else{
        if (daodBool){
            if (signalBool) outputFileName = "outputRootFiles/mc21_14TeV_HHbbbb_HLLHC_DAOD_NTUPLE.root"; // FIXME change these files 
            else outputFileName = "outputRootFiles/mc21_14TeV_jj_JZ3_DAOD_NTUPLE.root";
        }
        else{
            if (signalBool) outputFileName = "outputRootFiles/mc21_14TeV_HHbbbb_HLLHC_AOD_NTUPLE.root";
            else outputFileName = "outputRootFiles/mc21_14TeV_jj_JZ3_AOD_NTUPLE.root";
        }
    }
    
    

    // Create ROOT output file
    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    // Create TTrees - commented out trees that cannot be filled with DAOD information.
    //TTree* truthParticleTree = new TTree("truthParticleTree", "Tree storing event-wise information about truth particles");
    TTree* truthbTree = new TTree("truthbTree", "Tree storing event-wise information about truth particles");
    TTree* truthHiggsTree = new TTree("truthHiggsTree", "Tree storing event-wise information about truth particles");
    //TTree* truthVBFQuark = new TTree("truthVBFQuark", "Tree storing event-wise information about truth particles");
    TTree* caloTopoTowerTree = new TTree("caloTopoTowerTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* topo422Tree = new TTree("topo422Tree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gFexSRJTree = new TTree("gFexSRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gFexLeadingSRJTree = new TTree("gFexLeadingSRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gFexSubleadingSRJTree = new TTree("gFexSubleadingSRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gFexLRJTree = new TTree("gFexLRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gFexLeadingLRJTree = new TTree("gFexLeadingLRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gFexSubleadingLRJTree = new TTree("gFexSubleadingLRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* inTimeAntiKt4TruthJetsTree = new TTree("inTimeAntiKt4TruthJetsTree", "Tree storing event-wise Et, eta phi");
    TTree* leadingInTimeAntiKt4TruthJetsTree = new TTree("leadingInTimeAntiKt4TruthJetsTree", "Tree storing event-wise Et, eta phi");
    TTree* subleadingInTimeAntiKt4TruthJetsTree = new TTree("subleadingInTimeAntiKt4TruthJetsTree", "Tree storing event-wise Et, eta phi");
    TTree* jFexSRJTree = new TTree("jFexSRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* jFexLeadingSRJTree = new TTree("jFexLeadingSRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* jFexSubleadingSRJTree = new TTree("jFexSubleadingSRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* jFexLRJTree = new TTree("jFexLRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* jFexLeadingLRJTree = new TTree("jFexLeadingLRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* jFexSubleadingLRJTree = new TTree("jFexSubleadingLRJTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* hltAntiKt4EMTopoJetsTree = new TTree("hltAntiKt4EMTopoJetsTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* leadingHltAntiKt4EMTopoJetsTree = new TTree("leadingHltAntiKt4EMTopoJetsTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* subleadingHltAntiKt4EMTopoJetsTree = new TTree("subleadingHltAntiKt4EMTopoJetsTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* recoAntiKt10UFOCSSKJets = new TTree("recoAntiKt10UFOCSSKJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* leadingRecoAntiKt10UFOCSSKJets = new TTree("leadingRecoAntiKt10UFOCSSKJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* subleadingRecoAntiKt10UFOCSSKJets = new TTree("subleadingRecoAntiKt10UFOCSSKJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* truthAntiKt4TruthDressedWZJets = new TTree("truthAntiKt4TruthDressedWZJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* leadingTruthAntiKt4TruthDressedWZJets = new TTree("leadingTruthAntiKt4TruthDressedWZJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* subleadingTruthAntiKt4TruthDressedWZJets = new TTree("subleadingTruthAntiKt4TruthDressedWZJets", "Tree storing event-wise Et, Eta, Phi");
    /*TTree* recoAntiKt10LCTopoJets = new TTree("recoAntiKt10LCTopoJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* leadingRecoAntiKt10LCTopoJets = new TTree("leadingRecoAntiKt10LCTopoJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* subleadingRecoAntiKt10LCTopoJets = new TTree("recoAntiKt10LCTopoJets", "Tree storing event-wise Et, Eta, Phi");*/
    //TTree* gFexRhoRoI = new TTree("gFexRhoRoI", "Tree storing event-wise Et, Eta, Phi");

    // Variables to store data to be written to TTrees

    // Truth Particle vectors
    //std::vector<int> truthParticlePDGId, truthParticleStatus; 
    //std::vector<double> truthParticleEtValues, truthParticleEnergyValues, truthParticlepTValues, truthParticlepxValues, truthParticlepyValues, truthParticlepzValues, truthParticleEtaValues, truthParticlePhiValues;
    std::vector<unsigned int> higgsIndexValues, indexOfHiggsValues;
    std::vector<double> truthbquarksEtValues, truthbquarksEnergyValues, truthbquarkspTValues, truthbquarkspxValues, truthbquarkspyValues, truthbquarkspzValues, truthbquarksEtaValues, truthbquarksPhiValues;
    std::vector<double> truthHiggsEtValues, truthHiggsEnergyValues, truthHiggspTValues, truthHiggspxValues, truthHiggspyValues, truthHiggspzValues, truthHiggsEtaValues, truthHiggsPhiValues, truthHiggsInvMassValues;
    //std::vector<double> truthVBFQuarkValues, truthVBFQuarkEnergyValues, truthVBFQuarkpTValues, truthVBFQuarkpxValues, truthVBFQuarkpyValues, truthVBFQuarkpzValues, truthVBFQuarkEtaValues, truthVBFQuarkPhiValues;

    // Tower / cluster vectors
    std::vector<double> caloTopoTowerEtValues, caloTopoTowerEtaValues, caloTopoTowerPhiValues;
    std::vector<double> topo422EtValues, topo422EtaValues, topo422PhiValues;

    // L1Calo jets vectors
    std::vector<unsigned int> gFexSRJEtIndexValues; // stores index in sorted by Et list of jets
    std::vector<double> gFexSRJEtValues, gFexSRJEtaValues, gFexSRJPhiValues;
    std::vector<double> gFexSRJLeadingEtValues, gFexSRJLeadingEtaValues, gFexSRJLeadingPhiValues;
    std::vector<double> gFexSRJSubleadingEtValues, gFexSRJSubleadingEtaValues, gFexSRJSubleadingPhiValues;
    std::vector<unsigned int> gFexLRJEtIndexValues; // stores index in sorted by Et list of jets
    std::vector<double> gFexLRJEtValues, gFexLRJEtaValues, gFexLRJPhiValues;
    std::vector<double> gFexLRJLeadingEtValues, gFexLRJLeadingEtaValues, gFexLRJLeadingPhiValues;
    std::vector<double> gFexLRJSubleadingEtValues, gFexLRJSubleadingEtaValues, gFexLRJSubleadingPhiValues;
    std::vector<unsigned int> jFexSRJEtIndexValues;
    std::vector<double> jFexSRJEtValues, jFexSRJEtaValues, jFexSRJPhiValues;
    std::vector<double> jFexSRJLeadingEtValues, jFexSRJLeadingEtaValues, jFexSRJLeadingPhiValues;
    std::vector<double> jFexSRJSubleadingEtValues, jFexSRJSubleadingEtaValues, jFexSRJSubleadingPhiValues;
    std::vector<unsigned int> jFexLRJEtIndexValues;
    std::vector<double> jFexLRJEtValues, jFexLRJEtaValues, jFexLRJPhiValues;
    std::vector<double> jFexLRJLeadingEtValues, jFexLRJLeadingEtaValues, jFexLRJLeadingPhiValues;
    std::vector<double> jFexLRJSubleadingEtValues, jFexLRJSubleadingEtaValues, jFexLRJSubleadingPhiValues;
    //std::vector<double> gFexRhoRoIEtValues, gFexRhoRoIEtaValues, gFexRhoRoIPhiValues; // skip adding these for now.

    // HLT jets vectors
    std::vector<unsigned int> hltAntiKt4SRJEtIndexValues; // stores index in sorted by Et list of jets
    std::vector<double> hltAntiKt4SRJEtValues, hltAntiKt4SRJEtaValues, hltAntiKt4SRJPhiValues;
    std::vector<double> hltAntiKt4SRJLeadingEtValues, hltAntiKt4SRJLeadingEtaValues, hltAntiKt4SRJLeadingPhiValues;
    std::vector<double> hltAntiKt4SRJSubleadingEtValues, hltAntiKt4SRJSubleadingEtaValues, hltAntiKt4SRJSubleadingPhiValues;

    // Reco offline jets vectors
    std::vector<unsigned int> recoAntiKt10LRJEtIndexValues;
    std::vector<double> recoAntiKt10LRJEtValues, recoAntiKt10LRJEtaValues, recoAntiKt10LRJPhiValues;
    std::vector<double> recoAntiKt10LRJLeadingEtValues, recoAntiKt10LRJLeadingEtaValues, recoAntiKt10LRJLeadingPhiValues;
    std::vector<double> recoAntiKt10LRJSubleadingEtValues, recoAntiKt10LRJSubleadingEtaValues, recoAntiKt10LRJSubleadingPhiValues;

    // Truth WZ Antikt4 jets vectors
    std::vector<unsigned int> truthAntiKt4WZSRJEtIndexValues;
    std::vector<double> truthAntiKt4WZSRJEtValues, truthAntiKt4WZSRJEtaValues, truthAntiKt4WZSRJPhiValues;
    std::vector<double> truthAntiKt4WZSRJLeadingEtValues, truthAntiKt4WZSRJLeadingEtaValues, truthAntiKt4WZSRJLeadingPhiValues;
    std::vector<double> truthAntiKt4WZSRJSubleadingEtValues, truthAntiKt4WZSRJSubleadingEtaValues, truthAntiKt4WZSRJSubleadingPhiValues;

    // In time anti-kt 4 truth jets vectors
    std::vector<unsigned int> inTimeAntiKt4TruthSRJEtIndexValues;
    std::vector<double> inTimeAntiKt4TruthSRJEtValues, inTimeAntiKt4TruthSRJEtaValues, inTimeAntiKt4TruthSRJPhiValues;
    std::vector<double> inTimeAntiKt4TruthSRJLeadingEtValues, inTimeAntiKt4TruthSRJLeadingEtaValues, inTimeAntiKt4TruthSRJLeadingPhiValues;
    std::vector<double> inTimeAntiKt4TruthSRJSubleadingEtValues, inTimeAntiKt4TruthSRJSubleadingEtaValues, inTimeAntiKt4TruthSRJSubleadingPhiValues;

    // Create branches for each TTree FIXME don't have access to full truth record for DAOD samples
    /*
    truthParticleTree->Branch("Et", &truthParticleEtValues);
    truthParticleTree->Branch("Eta", &truthParticleEtaValues);
    truthParticleTree->Branch("Phi", &truthParticlePhiValues);
    truthParticleTree->Branch("pdgId", &truthParticlePDGId);
    truthParticleTree->Branch("pdgStatus", &truthParticleStatus);
    truthParticleTree->Branch("pT", &truthParticlepTValues);
    truthParticleTree->Branch("px", &truthParticlepxValues);
    truthParticleTree->Branch("py", &truthParticlepyValues);
    truthParticleTree->Branch("pz", &truthParticlepzValues);*/

    // truthbTree
    truthbTree->Branch("higgsIndex", &higgsIndexValues);
    truthbTree->Branch("Et", &truthbquarksEtValues);
    truthbTree->Branch("Eta", &truthbquarksEtaValues);
    truthbTree->Branch("Phi", &truthbquarksPhiValues);
    truthbTree->Branch("pT", &truthbquarkspTValues);
    truthbTree->Branch("px", &truthbquarkspxValues);
    truthbTree->Branch("py", &truthbquarkspyValues);
    truthbTree->Branch("pz", &truthbquarkspzValues);
    truthbTree->Branch("Energy", &truthbquarksEnergyValues);

    // truthHiggsTree
    truthHiggsTree->Branch("indexOfHiggs", &indexOfHiggsValues);
    truthHiggsTree->Branch("invMass", &truthHiggsInvMassValues);
    truthHiggsTree->Branch("Et", &truthHiggsEtValues);
    truthHiggsTree->Branch("Eta", &truthHiggsEtaValues);
    truthHiggsTree->Branch("Phi", &truthHiggsPhiValues);
    truthHiggsTree->Branch("pT", &truthHiggspTValues);
    truthHiggsTree->Branch("px", &truthHiggspxValues);
    truthHiggsTree->Branch("py", &truthHiggspyValues);
    truthHiggsTree->Branch("pz", &truthHiggspzValues);
    truthHiggsTree->Branch("Energy", &truthHiggsEnergyValues);

    // truthVBFQuark FIXME don't have for DAOD samples - unless can find corresponding aod events for which full truth record is stored
    /*
    truthVBFQuark->Branch("Et", &truthVBFQuarkValues);
    truthVBFQuark->Branch("Eta", &truthVBFQuarkEtaValues);
    truthVBFQuark->Branch("Phi", &truthVBFQuarkPhiValues);
    truthVBFQuark->Branch("pT", &truthVBFQuarkpTValues);
    truthVBFQuark->Branch("px", &truthVBFQuarkpxValues);
    truthVBFQuark->Branch("py", &truthVBFQuarkpyValues);
    truthVBFQuark->Branch("pz", &truthVBFQuarkpzValues);
    truthVBFQuark->Branch("Energy", &truthVBFQuarkEnergyValues);*/

    // caloTopoTowerTree
    caloTopoTowerTree->Branch("Et", &caloTopoTowerEtValues);
    caloTopoTowerTree->Branch("Eta", &caloTopoTowerEtaValues);
    caloTopoTowerTree->Branch("Phi", &caloTopoTowerPhiValues);

    // topo422Tree
    topo422Tree->Branch("Et", &topo422EtValues);
    topo422Tree->Branch("Eta", &topo422EtaValues);
    topo422Tree->Branch("Phi", &topo422PhiValues);

    // gFexSRJTree
    gFexSRJTree->Branch("EtIndex", &gFexSRJEtIndexValues);
    gFexSRJTree->Branch("Et", &gFexSRJEtValues);
    gFexSRJTree->Branch("Eta", &gFexSRJEtaValues);
    gFexSRJTree->Branch("Phi", &gFexSRJPhiValues);

    // gFexLeadingSRJTree
    gFexLeadingSRJTree->Branch("Et", &gFexSRJLeadingEtValues);
    gFexLeadingSRJTree->Branch("Eta", &gFexSRJLeadingEtaValues);
    gFexLeadingSRJTree->Branch("Phi", &gFexSRJLeadingPhiValues);

    // gFexSubleadingSRJTree
    gFexSubleadingSRJTree->Branch("Et", &gFexSRJSubleadingEtValues);
    gFexSubleadingSRJTree->Branch("Eta", &gFexSRJSubleadingEtaValues);
    gFexSubleadingSRJTree->Branch("Phi", &gFexSRJSubleadingPhiValues);

    // gFexLRJTree
    gFexLRJTree->Branch("EtIndex", &gFexLRJEtIndexValues);
    gFexLRJTree->Branch("Et", &gFexLRJEtValues);
    gFexLRJTree->Branch("Eta", &gFexLRJEtaValues);
    gFexLRJTree->Branch("Phi", &gFexLRJPhiValues);

    // gFexLeadingLRJTree
    gFexLeadingLRJTree->Branch("Et", &gFexLRJLeadingEtValues);
    gFexLeadingLRJTree->Branch("Eta", &gFexLRJLeadingEtaValues);
    gFexLeadingLRJTree->Branch("Phi", &gFexLRJLeadingPhiValues);

    // gFexSubleadingLRJTree
    gFexSubleadingLRJTree->Branch("Et", &gFexLRJSubleadingEtValues);
    gFexSubleadingLRJTree->Branch("Eta", &gFexLRJSubleadingEtaValues);
    gFexSubleadingLRJTree->Branch("Phi", &gFexLRJSubleadingPhiValues);

    // inTimeAntiKt4TruthJetsTree
    inTimeAntiKt4TruthJetsTree->Branch("EtIndex", &inTimeAntiKt4TruthSRJEtIndexValues);
    inTimeAntiKt4TruthJetsTree->Branch("Et", &inTimeAntiKt4TruthSRJEtValues);
    inTimeAntiKt4TruthJetsTree->Branch("Eta", &inTimeAntiKt4TruthSRJEtaValues);
    inTimeAntiKt4TruthJetsTree->Branch("Phi", &inTimeAntiKt4TruthSRJPhiValues);

    // leadingInTimeAntiKt4TruthJetsTree
    leadingInTimeAntiKt4TruthJetsTree->Branch("Et", &inTimeAntiKt4TruthSRJLeadingEtValues);
    leadingInTimeAntiKt4TruthJetsTree->Branch("Eta", &inTimeAntiKt4TruthSRJLeadingEtaValues);
    leadingInTimeAntiKt4TruthJetsTree->Branch("Phi", &inTimeAntiKt4TruthSRJLeadingPhiValues);

    // subleadingInTimeAntiKt4TruthJetsTree
    subleadingInTimeAntiKt4TruthJetsTree->Branch("Et", &inTimeAntiKt4TruthSRJSubleadingEtValues);
    subleadingInTimeAntiKt4TruthJetsTree->Branch("Eta", &inTimeAntiKt4TruthSRJSubleadingEtaValues);
    subleadingInTimeAntiKt4TruthJetsTree->Branch("Phi", &inTimeAntiKt4TruthSRJSubleadingPhiValues);

    // jFexSRJTree
    jFexSRJTree->Branch("EtIndex", &jFexSRJEtIndexValues);
    jFexSRJTree->Branch("Et", &jFexSRJEtValues);
    jFexSRJTree->Branch("Eta", &jFexSRJEtaValues);
    jFexSRJTree->Branch("Phi", &jFexSRJPhiValues);

    // jFexLeadingSRJTree
    jFexLeadingSRJTree->Branch("Et", &jFexSRJLeadingEtValues);
    jFexLeadingSRJTree->Branch("Eta", &jFexSRJLeadingEtaValues);
    jFexLeadingSRJTree->Branch("Phi", &jFexSRJLeadingPhiValues);

    // jFexSubleadingSRJTree
    jFexSubleadingSRJTree->Branch("Et", &jFexSRJSubleadingEtValues);
    jFexSubleadingSRJTree->Branch("Eta", &jFexSRJSubleadingEtaValues);
    jFexSubleadingSRJTree->Branch("Phi", &jFexSRJSubleadingPhiValues);

    // jFexLRJTree
    jFexLRJTree->Branch("EtIndex", &jFexLRJEtIndexValues);
    jFexLRJTree->Branch("Et", &jFexLRJEtValues);
    jFexLRJTree->Branch("Eta", &jFexLRJEtaValues);
    jFexLRJTree->Branch("Phi", &jFexLRJPhiValues);

    // jFexLeadingLRJTree
    jFexLeadingLRJTree->Branch("Et", &jFexLRJLeadingEtValues);
    jFexLeadingLRJTree->Branch("Eta", &jFexLRJLeadingEtaValues);
    jFexLeadingLRJTree->Branch("Phi", &jFexLRJLeadingPhiValues);

    // jFexSubleadingLRJTree
    jFexSubleadingLRJTree->Branch("Et", &jFexLRJSubleadingEtValues);
    jFexSubleadingLRJTree->Branch("Eta", &jFexLRJSubleadingEtaValues);
    jFexSubleadingLRJTree->Branch("Phi", &jFexLRJSubleadingPhiValues);

    // hltAntiKt4EMTopoJetsTreehltAntiKt4SRJEtIndexValues
    hltAntiKt4EMTopoJetsTree->Branch("EtIndex", &recoAntiKt10LRJEtIndexValues);
    hltAntiKt4EMTopoJetsTree->Branch("Et", &hltAntiKt4SRJEtValues);
    hltAntiKt4EMTopoJetsTree->Branch("Eta", &hltAntiKt4SRJEtaValues);
    hltAntiKt4EMTopoJetsTree->Branch("Phi", &hltAntiKt4SRJPhiValues);

    // leadingHltAntiKt4EMTopoJetsTree
    leadingHltAntiKt4EMTopoJetsTree->Branch("Et", &hltAntiKt4SRJLeadingEtValues);
    leadingHltAntiKt4EMTopoJetsTree->Branch("Eta", &hltAntiKt4SRJLeadingEtaValues);
    leadingHltAntiKt4EMTopoJetsTree->Branch("Phi", &hltAntiKt4SRJLeadingPhiValues);

    // subleadingHltAntiKt4EMTopoJetsTree
    subleadingHltAntiKt4EMTopoJetsTree->Branch("Et", &hltAntiKt4SRJSubleadingEtValues);
    subleadingHltAntiKt4EMTopoJetsTree->Branch("Eta", &hltAntiKt4SRJSubleadingEtaValues);
    subleadingHltAntiKt4EMTopoJetsTree->Branch("Phi", &hltAntiKt4SRJSubleadingPhiValues);

    // recoAntiKt10UFOCSSKJets
    recoAntiKt10UFOCSSKJets->Branch("EtIndex", &recoAntiKt10LRJEtIndexValues);
    recoAntiKt10UFOCSSKJets->Branch("Et", &recoAntiKt10LRJEtValues);
    recoAntiKt10UFOCSSKJets->Branch("Eta", &recoAntiKt10LRJEtaValues);
    recoAntiKt10UFOCSSKJets->Branch("Phi", &recoAntiKt10LRJPhiValues);

    // leadingRecoAntiKt10UFOCSSKJets
    leadingRecoAntiKt10UFOCSSKJets->Branch("Et", &recoAntiKt10LRJLeadingEtValues);
    leadingRecoAntiKt10UFOCSSKJets->Branch("Eta", &recoAntiKt10LRJLeadingEtaValues);
    leadingRecoAntiKt10UFOCSSKJets->Branch("Phi", &recoAntiKt10LRJLeadingPhiValues);

    // subleadingRecoAntiKt10UFOCSSKJets
    subleadingRecoAntiKt10UFOCSSKJets->Branch("Et", &recoAntiKt10LRJSubleadingEtValues);
    subleadingRecoAntiKt10UFOCSSKJets->Branch("Eta", &recoAntiKt10LRJSubleadingEtaValues);
    subleadingRecoAntiKt10UFOCSSKJets->Branch("Phi", &recoAntiKt10LRJSubleadingPhiValues);

    // truthAntiKt4TruthDressedWZJets
    truthAntiKt4TruthDressedWZJets->Branch("EtIndex", &truthAntiKt4WZSRJEtIndexValues);
    truthAntiKt4TruthDressedWZJets->Branch("Et", &truthAntiKt4WZSRJEtValues);
    truthAntiKt4TruthDressedWZJets->Branch("Eta", &truthAntiKt4WZSRJEtaValues);
    truthAntiKt4TruthDressedWZJets->Branch("Phi", &truthAntiKt4WZSRJPhiValues);

    // leadingTruthAntiKt4TruthDressedWZJets
    leadingTruthAntiKt4TruthDressedWZJets->Branch("Et", &truthAntiKt4WZSRJLeadingEtValues);
    leadingTruthAntiKt4TruthDressedWZJets->Branch("Eta", &truthAntiKt4WZSRJLeadingEtaValues);
    leadingTruthAntiKt4TruthDressedWZJets->Branch("Phi", &truthAntiKt4WZSRJLeadingPhiValues);

    // subleadingTruthAntiKt4TruthDressedWZJets
    subleadingTruthAntiKt4TruthDressedWZJets->Branch("Et", &truthAntiKt4WZSRJSubleadingEtValues);
    subleadingTruthAntiKt4TruthDressedWZJets->Branch("Eta", &truthAntiKt4WZSRJSubleadingEtaValues);
    subleadingTruthAntiKt4TruthDressedWZJets->Branch("Phi", &truthAntiKt4WZSRJSubleadingPhiValues);

    // recoAntiKt10LCTopoJets - unsure if want to use these at all?? 
    /*
    recoAntiKt10LCTopoJets->Branch("Et", &recoAntiKt10LRJEtValues);
    recoAntiKt10LCTopoJets->Branch("Eta", &recoAntiKt10LRJEtaValues);
    recoAntiKt10LCTopoJets->Branch("Phi", &recoAntiKt10LRJPhiValues);

    // leadingRecoAntiKt10LCTopoJets
    leadingRecoAntiKt10LCTopoJets->Branch("Et", &recoAntiKt10LRJLeadingEtValues);
    leadingRecoAntiKt10LCTopoJets->Branch("Eta", &recoAntiKt10LRJLeadingEtaValues);
    leadingRecoAntiKt10LCTopoJets->Branch("Phi", &recoAntiKt10LRJLeadingPhiValues);

    // subleadingRecoAntiKt10LCTopoJets
    subleadingRecoAntiKt10LCTopoJets->Branch("Et", &recoAntiKt10LRJSubleadingEtValues);
    subleadingRecoAntiKt10LCTopoJets->Branch("Eta", &recoAntiKt10LRJSubleadingEtaValues);
    subleadingRecoAntiKt10LCTopoJets->Branch("Phi", &recoAntiKt10LRJSubleadingPhiValues);*/

    // gFexRhoRoI
    /*gFexRhoRoI->Branch("Et", &gFexRhoRoIEtValues);
    gFexRhoRoI->Branch("Eta", &gFexRhoRoIEtaValues);
    gFexRhoRoI->Branch("Phi", &gFexRhoRoIPhiValues);*/
    int higgsPassEventCounter = 0;
    // Collect input file names
    vector<string> fileNames;
    for (const auto& entry : filesystem::directory_iterator(fileDir)) {
        if (entry.path().extension() == ".1") { // Ends with .1
            fileNames.push_back(entry.path().string());
        }
    }

    if (fileNames.empty()) {
        cout << "No ROOT files found in directory." << endl;
        return;
    }

    cout << "Found " << fileNames.size() << " files." << endl;

    // Main processing loop
    int fileIt = 0;
    for (const auto& fileName : fileNames) {
        fileIt++; 
        //if (fileIt > 9) break; 
        cout << "Processing file: " << fileName << endl;

        TFile* f = TFile::Open(fileName.c_str());
        if (!f || f->IsZombie()) {
            cerr << "Could not open " << fileName << endl;
            continue;
        }

        xAOD::TEvent event(xAOD::TEvent::kClassAccess);
        if (!event.readFrom(f).isSuccess()) {
            cerr << "Cannot read xAOD from file." << endl;
            continue;
        }

        std::cout << "  Number of events: " << event.getEntries() << endl;

        
        
        for (Long64_t iEvt = 0; iEvt < event.getEntries(); ++iEvt) {
            event.getEntry(iEvt);

            if ((iEvt % 100) == 0) std::cout << "iEvt: " << iEvt << "\n";
            

            // -- retrieve collections from DOAD ---
            const xAOD::TruthParticleContainer* TruthParticles = nullptr;
            if (!event.retrieve(TruthParticles, "TruthParticles").isSuccess()) {
                cerr << "Cannot access TruthParticles" << endl;
                continue;
            }

            const xAOD::JetContainer* HLT_AntiKt4EMTopoJets_subjesIS = nullptr;
            if (!event.retrieve(HLT_AntiKt4EMTopoJets_subjesIS, "HLT_AntiKt4EMTopoJets_subjesIS").isSuccess()) {
                cerr << "Failed to retrieve HLT jets" << endl;
                continue;
            }

            const xAOD::gFexJetRoIContainer* L1_gFexSRJetRoI = nullptr;
            if (!event.retrieve(L1_gFexSRJetRoI, "L1_gFexSRJetRoI").isSuccess()) {
                std::cerr << "Failed to retrieve gFex SR jets" << std::endl;
                continue;
            }

            const xAOD::gFexJetRoIContainer* L1_gFexLRJetRoI = nullptr;
            if (!event.retrieve(L1_gFexLRJetRoI, "L1_gFexLRJetRoI").isSuccess()) {
                cerr << "Failed to retrieve gFex LR jets" << endl;
                continue;
            }

            const DataVector<xAOD::jFexSRJetRoI_v1>* L1_jFexSRJetRoI = nullptr;
            if (!event.retrieve(L1_jFexSRJetRoI, "L1_jFexSRJetRoI").isSuccess()) {
                std::cerr << "Failed to retrieve jFex SR jets" << std::endl;
                continue;
            }

            const DataVector<xAOD::jFexLRJetRoI_v1>* L1_jFexLRJetRoI = nullptr;
            if (!event.retrieve(L1_jFexLRJetRoI, "L1_jFexLRJetRoI").isSuccess()) {
                std::cerr << "Failed to retrieve jFex SR jets" << std::endl;
                continue;
            }

            /*
            const xAOD::gFexJetRoIContainer* L1_gFexRhoRoI = nullptr;
            if (!event.retrieve(L1_gFexRhoRoI, "L1_gFexRhoRoI").isSuccess()) {
                std::cerr << "Failed to retrieve gFex energy density" << std::endl;
                continue;
            }*/

            const xAOD::JetContainer* InTimeAntiKt4TruthJets = nullptr;
            if (!event.retrieve(InTimeAntiKt4TruthJets, "InTimeAntiKt4TruthJets").isSuccess()) {
                cerr << "Failed to retrieve Truth jets" << endl;
                continue;
            }


            // Retrieve the CaloTopoClusters422 container
            const DataVector<xAOD::CaloCluster_v1>* CaloTopoClusters422 = nullptr;
            if (!event.retrieve(CaloTopoClusters422, "CaloTopoClusters422").isSuccess()) {
                std::cerr << "Failed to retrieve CaloTopoClusters422" << std::endl;
                continue;
            }

            // Retrieve the CaloCalAllTopoTowers container
            const DataVector<xAOD::CaloCluster_v1>* CaloCalAllTopoTowers = nullptr;
            if (!event.retrieve(CaloCalAllTopoTowers, "CaloCalAllTopoTowers").isSuccess()) {
                std::cerr << "Failed to retrieve CaloCalAllTopoTowers" << std::endl;
                continue;
            }


            // FIXME add all this, fill trees and clear vectors:
            hltAntiKt4SRJEtaValues.clear();
            hltAntiKt4SRJEtValues.clear();
            hltAntiKt4SRJPhiValues.clear();
            hltAntiKt4SRJEtIndexValues.clear();
            hltAntiKt4SRJLeadingEtValues.clear();
            hltAntiKt4SRJLeadingEtaValues.clear();
            hltAntiKt4SRJLeadingPhiValues.clear();
            hltAntiKt4SRJSubleadingEtValues.clear();
            hltAntiKt4SRJSubleadingEtaValues.clear();
            hltAntiKt4SRJSubleadingPhiValues.clear();
            inTimeAntiKt4TruthSRJEtValues.clear();
            inTimeAntiKt4TruthSRJEtIndexValues.clear();
            inTimeAntiKt4TruthSRJEtaValues.clear();
            inTimeAntiKt4TruthSRJPhiValues.clear();
            inTimeAntiKt4TruthSRJLeadingEtValues.clear();
            inTimeAntiKt4TruthSRJLeadingEtaValues.clear();
            inTimeAntiKt4TruthSRJLeadingPhiValues.clear();
            inTimeAntiKt4TruthSRJSubleadingEtValues.clear();
            inTimeAntiKt4TruthSRJSubleadingEtaValues.clear();
            inTimeAntiKt4TruthSRJSubleadingPhiValues.clear();
            recoAntiKt10LRJEtIndexValues.clear();
            recoAntiKt10LRJEtValues.clear();
            recoAntiKt10LRJEtaValues.clear();
            recoAntiKt10LRJPhiValues.clear();
            recoAntiKt10LRJLeadingEtValues.clear();
            recoAntiKt10LRJLeadingEtaValues.clear();
            recoAntiKt10LRJLeadingPhiValues.clear();
            recoAntiKt10LRJSubleadingEtValues.clear();
            recoAntiKt10LRJSubleadingEtaValues.clear();
            recoAntiKt10LRJSubleadingPhiValues.clear();
            truthAntiKt4WZSRJEtIndexValues.clear();
            truthAntiKt4WZSRJEtValues.clear();
            truthAntiKt4WZSRJEtaValues.clear();
            truthAntiKt4WZSRJPhiValues.clear();
            truthAntiKt4WZSRJLeadingEtValues.clear();
            truthAntiKt4WZSRJLeadingEtaValues.clear();
            truthAntiKt4WZSRJLeadingPhiValues.clear();
            truthAntiKt4WZSRJSubleadingEtValues.clear();
            truthAntiKt4WZSRJSubleadingEtaValues.clear();
            truthAntiKt4WZSRJSubleadingPhiValues.clear();
            
            gFexLRJEtValues.clear();
            gFexLRJEtaValues.clear();
            gFexLRJPhiValues.clear();
            gFexLRJEtIndexValues.clear();
            gFexLRJLeadingEtValues.clear();
            gFexLRJLeadingEtaValues.clear();
            gFexLRJLeadingPhiValues.clear();
            gFexLRJSubleadingEtValues.clear();
            gFexLRJSubleadingEtaValues.clear();
            gFexLRJSubleadingPhiValues.clear();
            gFexSRJEtValues.clear();
            gFexSRJEtaValues.clear();
            gFexSRJPhiValues.clear();
            gFexSRJEtIndexValues.clear();
            gFexSRJLeadingEtValues.clear();
            gFexSRJLeadingEtaValues.clear();
            gFexSRJLeadingPhiValues.clear();
            gFexSRJSubleadingEtValues.clear();
            gFexSRJSubleadingEtaValues.clear();
            gFexSRJSubleadingPhiValues.clear();
            jFexSRJEtValues.clear();
            jFexSRJEtaValues.clear();
            jFexSRJPhiValues.clear();
            jFexSRJEtIndexValues.clear();
            jFexSRJLeadingEtValues.clear();
            jFexSRJLeadingEtaValues.clear();
            jFexSRJLeadingPhiValues.clear();
            jFexSRJSubleadingEtValues.clear();
            jFexSRJSubleadingEtaValues.clear();
            jFexSRJSubleadingPhiValues.clear();
            jFexLRJEtValues.clear();
            jFexLRJEtaValues.clear();
            jFexLRJPhiValues.clear();
            jFexLRJEtIndexValues.clear();
            jFexLRJLeadingEtValues.clear();
            jFexLRJLeadingEtaValues.clear();
            jFexLRJLeadingPhiValues.clear();
            jFexLRJSubleadingEtValues.clear();
            jFexLRJSubleadingEtaValues.clear();
            jFexLRJSubleadingPhiValues.clear();
            topo422EtValues.clear();
            topo422EtaValues.clear();
            topo422PhiValues.clear();
            caloTopoTowerEtValues.clear();
            caloTopoTowerEtaValues.clear();
            caloTopoTowerPhiValues.clear();
            truthHiggsInvMassValues.clear();
            truthHiggsEtValues.clear();
            truthHiggsEtaValues.clear();
            truthHiggsPhiValues.clear();
            truthHiggspTValues.clear();
            truthHiggspxValues.clear();
            truthHiggspyValues.clear();
            truthHiggspzValues.clear();
            truthHiggsEnergyValues.clear();
            indexOfHiggsValues.clear();
            higgsIndexValues.clear();
            truthbquarksEtValues.clear();
            truthbquarksEtaValues.clear();
            truthbquarksPhiValues.clear();
            truthbquarkspTValues.clear();
            truthbquarkspxValues.clear();
            truthbquarkspyValues.clear();
            truthbquarkspzValues.clear();
            truthbquarksEnergyValues.clear();

            // Initialize per-event counters
            unsigned int higgs_counter = -1;
            bool higgsPtCutPassed[2] = {false, false}; 
            // --- Loop over TruthParticles (for Higgs and B's) ---
            std::vector<std::vector<float > > allb_list;
            for (const auto& el : *TruthParticles) {
                if (el->pdgId() == 25 && el->status() == 22) {
                    higgs_counter++;

                    float pt = el->pt() / 1000.0;
                    float eta = el->eta();
                    float phi = el->phi();
                    float px = el->px() / 1000.0;
                    float py = el->py() / 1000.0;
                    float pz = el->pz() / 1000.0;
                    float energy = el->e() / 1000.0;
                    float et = energy / cosh(eta);

                    // Fill truthHiggsTree variables
                    truthHiggsEtValues.push_back(et);
                    truthHiggsEtaValues.push_back(eta);
                    truthHiggsPhiValues.push_back(phi);
                    truthHiggspTValues.push_back(pt);
                    truthHiggspxValues.push_back(px);
                    truthHiggspyValues.push_back(py);
                    truthHiggspzValues.push_back(pz);
                    truthHiggsEnergyValues.push_back(energy);
                    indexOfHiggsValues.push_back(higgs_counter);

                    if (higgs_counter == 0) { // fill b's for both higgs, and store index of which higgs they correspond to.  
                        if (pt > 0) higgsPtCutPassed[0] = true;
                        std::vector<const xAOD::TruthParticle*> b1_list;
                        find_non_higgs_daughters(el, b1_list);

                        if (b1_list.size() == 2) {
                            TLorentzVector b1, b2;
                            b1.SetPxPyPzE(b1_list[0]->px() / 1000.0, b1_list[0]->py() / 1000.0,
                                        b1_list[0]->pz() / 1000.0, b1_list[0]->e() / 1000.0);
                            b2.SetPxPyPzE(b1_list[1]->px() / 1000.0, b1_list[1]->py() / 1000.0,
                                        b1_list[1]->pz() / 1000.0, b1_list[1]->e() / 1000.0);
                            double mH = (b1 + b2).M();  // invariant mass in GeV
                            truthHiggsInvMassValues.push_back(mH);
                        } // compute invariant mass

                        for (const auto* b : b1_list) {
                            if (!b) continue;

                            float pt   = b->pt() / 1000.0;
                            float eta  = b->eta();
                            float phi  = b->phi();
                            float px   = b->px() / 1000.0;
                            float py   = b->py() / 1000.0;
                            float pz   = b->pz() / 1000.0;
                            float E    = b->e()  / 1000.0;
                            float Et   = E / std::cosh(eta);

                            higgsIndexValues.push_back(0);  // from 1st Higgs
                            truthbquarksEtValues.push_back(Et);
                            truthbquarksEtaValues.push_back(eta);
                            truthbquarksPhiValues.push_back(phi);
                            truthbquarkspTValues.push_back(pt);
                            truthbquarkspxValues.push_back(px);
                            truthbquarkspyValues.push_back(py);
                            truthbquarkspzValues.push_back(pz);
                            truthbquarksEnergyValues.push_back(E);
                        }
                    }

                    if (higgs_counter == 1) {
                        if (pt > 0) higgsPtCutPassed[1] = true;
                        std::vector<const xAOD::TruthParticle*> b2_list;
                        find_non_higgs_daughters(el, b2_list);

                        if (b2_list.size() == 2) {
                            TLorentzVector b1, b2;
                            b1.SetPxPyPzE(b2_list[0]->px() / 1000.0, b2_list[0]->py() / 1000.0,
                                        b2_list[0]->pz() / 1000.0, b2_list[0]->e() / 1000.0);
                            b2.SetPxPyPzE(b2_list[1]->px() / 1000.0, b2_list[1]->py() / 1000.0,
                                        b2_list[1]->pz() / 1000.0, b2_list[1]->e() / 1000.0);
                            double mH = (b1 + b2).M();  // invariant mass in GeV
                            truthHiggsInvMassValues.push_back(mH);
                        }

                        for (const auto* b : b2_list) {
                            if (!b) continue;

                            float pt   = b->pt() / 1000.0;
                            float eta  = b->eta();
                            float phi  = b->phi();
                            float px   = b->px() / 1000.0;
                            float py   = b->py() / 1000.0;
                            float pz   = b->pz() / 1000.0;
                            float E    = b->e()  / 1000.0;
                            float Et   = E / std::cosh(eta);

                            higgsIndexValues.push_back(1);  // from 2nd Higgs
                            truthbquarksEtValues.push_back(Et);
                            truthbquarksEtaValues.push_back(eta);
                            truthbquarksPhiValues.push_back(phi);
                            truthbquarkspTValues.push_back(pt);
                            truthbquarkspxValues.push_back(px);
                            truthbquarkspyValues.push_back(py);
                            truthbquarkspzValues.push_back(pz);
                            truthbquarksEnergyValues.push_back(E);
                        } // loop through b2 list
                    } // if 2nd higgs in event
                } // if higgs truth particle of interest
            } // loop through truth bosons with decay particles
            bool higgsPtCutsPassed = (higgsPtCutPassed[0] || higgsPtCutPassed[1]) || !signalBool;
            if (higgsPtCutsPassed) higgsPassEventCounter++;

            // Loop over clusters and fill Et, Eta, Phi
            unsigned int topotower_it = 0;
            for (const auto* cluster : *CaloCalAllTopoTowers) {
                if (!cluster) continue;

                double et = cluster->e() / cosh(cluster->eta()) / 1000.0; // Et in GeV
                caloTopoTowerEtValues.push_back(et);
                caloTopoTowerEtaValues.push_back(cluster->eta());
                caloTopoTowerPhiValues.push_back(cluster->phi());

                if (et < 0) continue; // don't store to digitized memories

                // Digitize each variable
                int phi_bin = digitize(cluster->phi(), phi_bit_length_, phi_min_, phi_max_);
                int eta_bin = digitize(cluster->eta(), eta_bit_length_, eta_min_, eta_max_);
                int et_bin  = digitize(et, et_bit_length_,
                              static_cast<double>(et_min_), static_cast<double>(et_max_));
                                        
                // 2. Build binary string (for debug or text output)
                std::stringstream binary_ss;
                binary_ss << std::bitset<et_bit_length_>(et_bin) << "|"
                        << std::bitset<eta_bit_length_>(eta_bin) << "|"
                        << std::bitset<phi_bit_length_>(phi_bin);
                std::string binary_word = binary_ss.str();

                // 3. Pack into 27-bit word (stored as uint32_t)
                uint32_t packed_word = (et_bin  << (eta_bit_length_ + phi_bit_length_)) |
                                    (eta_bin << phi_bit_length_) |
                                    (phi_bin);
                                    
                
                if (higgsPtCutsPassed || !signalBool){
                    if (topotower_it == 0) {
                        f_topotower << "Event : " << std::dec << iEvt << "\n";
                    }
                    // 4. Write to output file
                    f_topotower << "0x" << std::hex << std::setw(2) << std::setfill('0') << topotower_it
                                << " " << binary_word
                                << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    topotower_it++;
                }
            }

            unsigned int topocluster422_it = 0;
            // Loop over the clusters and store Et, Eta, Phi
            for (const auto* cluster : *CaloTopoClusters422) {
                if (!cluster) continue;

                double et = cluster->e() / cosh(cluster->eta()) / 1000.0; // Et in GeV
                topo422EtValues.push_back(et);
                topo422EtaValues.push_back(cluster->eta());
                topo422PhiValues.push_back(cluster->phi());

                

                if (et < 0) continue; // FIXME for now don't store negative Et to digitized memories
                // Digitize each variable
                int phi_bin = digitize(cluster->phi(), phi_bit_length_, phi_min_, phi_max_);
                int eta_bin = digitize(cluster->eta(), eta_bit_length_, eta_min_, eta_max_);
                int et_bin  = digitize(et, et_bit_length_,
                              static_cast<double>(et_min_), static_cast<double>(et_max_));

                //std::cout << "eta: " << cluster->eta() << " and eta_bin : " << eta_bin << "\n";
                                        
                // 2. Build binary string (for debug or text output)
                std::stringstream binary_ss;
                binary_ss << std::bitset<et_bit_length_>(et_bin) << "|"
                        << std::bitset<eta_bit_length_>(eta_bin) << "|"
                        << std::bitset<phi_bit_length_>(phi_bin);
                std::string binary_word = binary_ss.str();

                // 3. Pack into 27-bit word (stored as uint32_t)
                uint32_t packed_word = (et_bin  << (eta_bit_length_ + phi_bit_length_)) |
                                    (eta_bin << phi_bit_length_) |
                                    (phi_bin);

                //std::cout << "binary: " << binary_word << "\n";

                // 4. Write to output file
                
                if(higgsPtCutsPassed || !signalBool){
                    if (topocluster422_it == 0) {
                        f_topo << "Event : " << std::dec << iEvt << "\n";
                    }
                    f_topo << "0x" << std::hex << std::setw(2) << std::setfill('0') << topocluster422_it
                            << " " << binary_word
                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    topocluster422_it++; 
                }
            }

            // Temporary vector for sorting by Et
            std::vector<std::pair<size_t, double>> jFexSRJetEtWithIndex;

            for (size_t i = 0; i < L1_jFexSRJetRoI->size(); ++i) {
                const auto& jet = (*L1_jFexSRJetRoI)[i];
                double et = jet->et() / 1000.0; // Already in GeV
                jFexSRJetEtWithIndex.emplace_back(i, et);
                //std::cout << "i : " << i << " and Et : " << et << "\n";
            }

            // Sort descending by Et
            std::sort(jFexSRJetEtWithIndex.begin(), jFexSRJetEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });

            // Fill vectors in sorted order
            unsigned int jfex_it = 0;
            for (const auto& [index, et] : jFexSRJetEtWithIndex) {
                //std::cout << "index: " << index << "\n";
                //std::cout << "et : " << et << "\n";
                const auto& jet = (*L1_jFexSRJetRoI)[index];

                jFexSRJEtIndexValues.push_back(static_cast<unsigned int>(index));
                jFexSRJEtValues.push_back(et);
                jFexSRJEtaValues.push_back(jet->eta());
                jFexSRJPhiValues.push_back(jet->phi());

                if (et < 0) continue; // FIXME

                // Digitize each variable
                int phi_bin = digitize(jet->phi(), phi_bit_length_, phi_min_, phi_max_);
                int eta_bin = digitize(jet->eta(), eta_bit_length_, eta_min_, eta_max_);
                int et_bin  = digitize(et, et_bit_length_,
                              static_cast<double>(et_min_), static_cast<double>(et_max_));
                                        
                // 2. Build binary string (for debug or text output)
                std::stringstream binary_ss;
                binary_ss << std::bitset<et_bit_length_>(et_bin) << "|"
                        << std::bitset<eta_bit_length_>(eta_bin) << "|"
                        << std::bitset<phi_bit_length_>(phi_bin);
                std::string binary_word = binary_ss.str();

                // 3. Pack into 27-bit word (stored as uint32_t)
                uint32_t packed_word = (et_bin  << (eta_bit_length_ + phi_bit_length_)) |
                                    (eta_bin << phi_bit_length_) |
                                    (phi_bin);

                // 4. Write to output file
                
                if(higgsPtCutsPassed || !signalBool){
                    if (jfex_it == 0) {
                        f_jfex << "Event : " << std::dec << iEvt << "\n";
                    }
                    f_jfex << "0x" << std::hex << std::setw(2) << std::setfill('0') << jfex_it
                            << " " << binary_word
                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    jfex_it++;
                }
            }

            // Leading jet
            if (!jFexSRJetEtWithIndex.empty()) {
                const auto& leading = (*L1_jFexSRJetRoI)[jFexSRJetEtWithIndex[0].first];
                jFexSRJLeadingEtValues.push_back(jFexSRJEtValues[0]);
                jFexSRJLeadingEtaValues.push_back(leading->eta());
                jFexSRJLeadingPhiValues.push_back(leading->phi());
            }

            // Subleading jet
            if (jFexSRJetEtWithIndex.size() > 1) {
                const auto& subleading = (*L1_jFexSRJetRoI)[jFexSRJetEtWithIndex[1].first];
                jFexSRJSubleadingEtValues.push_back(jFexSRJEtValues[1]);
                jFexSRJSubleadingEtaValues.push_back(subleading->eta());
                jFexSRJSubleadingPhiValues.push_back(subleading->phi());
            }

            // Temporary vector for sorting by Et
            std::vector<std::pair<size_t, double>> jFexLRJetEtWithIndex;

            for (size_t i = 0; i < L1_jFexLRJetRoI->size(); ++i) {
                const auto& jet = (*L1_jFexLRJetRoI)[i];
                double et = jet->et() / 1000.0; // Already in GeV
                jFexLRJetEtWithIndex.emplace_back(i, et);
            }

            // Sort descending by Et
            std::sort(jFexLRJetEtWithIndex.begin(), jFexLRJetEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });

            // Fill vectors in sorted order
            for (const auto& [index, et] : jFexLRJetEtWithIndex) {
                const auto& jet = (*L1_jFexLRJetRoI)[index];

                jFexLRJEtIndexValues.push_back(static_cast<unsigned int>(index));
                jFexLRJEtValues.push_back(et);
                jFexLRJEtaValues.push_back(jet->eta());
                jFexLRJPhiValues.push_back(jet->phi());

            }

            // Leading jet
            if (!jFexLRJetEtWithIndex.empty()) {
                const auto& leading = (*L1_jFexLRJetRoI)[jFexLRJetEtWithIndex[0].first];
                jFexLRJLeadingEtValues.push_back(jFexLRJEtValues[0]);
                jFexLRJLeadingEtaValues.push_back(leading->eta());
                jFexLRJLeadingPhiValues.push_back(leading->phi());
            }

            // Subleading jet
            if (jFexLRJetEtWithIndex.size() > 1) {
                const auto& subleading = (*L1_jFexLRJetRoI)[jFexLRJetEtWithIndex[1].first];
                jFexLRJSubleadingEtValues.push_back(jFexLRJEtValues[1]);
                jFexLRJSubleadingEtaValues.push_back(subleading->eta());
                jFexLRJSubleadingPhiValues.push_back(subleading->phi());
            }

            // Temporary vector for sorting by Et
            std::vector<std::pair<size_t, double>> gFexSRJetEtWithIndex;

            for (size_t i = 0; i < L1_gFexSRJetRoI->size(); ++i) {
                const auto& jet = (*L1_gFexSRJetRoI)[i];
                double et = jet->et() / 1000.0; // Already in GeV
                gFexSRJetEtWithIndex.emplace_back(i, et);
            }

            // Sort descending by Et
            std::sort(gFexSRJetEtWithIndex.begin(), gFexSRJetEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });

            // Fill vectors in sorted order
            unsigned int gfex_it = 0;
            for (const auto& [index, et] : gFexSRJetEtWithIndex) {
                const auto& jet = (*L1_gFexSRJetRoI)[index];

                gFexSRJEtIndexValues.push_back(static_cast<unsigned int>(index));
                gFexSRJEtValues.push_back(et);
                gFexSRJEtaValues.push_back(jet->eta());
                gFexSRJPhiValues.push_back(jet->phi());

                if (et < 0) continue;

                // Digitize each variable
                int phi_bin = digitize(jet->phi(), phi_bit_length_, phi_min_, phi_max_);
                int eta_bin = digitize(jet->eta(), eta_bit_length_, eta_min_, eta_max_);
                int et_bin  = digitize(et, et_bit_length_,
                              static_cast<double>(et_min_), static_cast<double>(et_max_));
                                        
                // 2. Build binary string (for debug or text output)
                std::stringstream binary_ss;
                binary_ss << std::bitset<et_bit_length_>(et_bin) << "|"
                        << std::bitset<eta_bit_length_>(eta_bin) << "|"
                        << std::bitset<phi_bit_length_>(phi_bin);
                std::string binary_word = binary_ss.str();

                // 3. Pack into 27-bit word (stored as uint32_t)
                uint32_t packed_word = (et_bin  << (eta_bit_length_ + phi_bit_length_)) |
                                    (eta_bin << phi_bit_length_) |
                                    (phi_bin);

                // 4. Write to output file
                
                if(higgsPtCutsPassed || !signalBool){
                    if (gfex_it == 0) {
                        f_gfex << "Event : " << std::dec << iEvt << "\n";
                    }
                    f_gfex << "0x" << std::hex << std::setw(2) << std::setfill('0') << gfex_it
                            << " " << binary_word
                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    gfex_it++;
                }   
            }

            // Leading jet
            if (!gFexSRJetEtWithIndex.empty()) {
                const auto& leading = (*L1_gFexSRJetRoI)[gFexSRJetEtWithIndex[0].first];
                gFexSRJLeadingEtValues.push_back(gFexSRJEtValues[0]);
                gFexSRJLeadingEtaValues.push_back(leading->eta());
                gFexSRJLeadingPhiValues.push_back(leading->phi());
            }

            // Subleading jet
            if (gFexSRJetEtWithIndex.size() > 1) {
                const auto& subleading = (*L1_gFexSRJetRoI)[gFexSRJetEtWithIndex[1].first];
                gFexSRJSubleadingEtValues.push_back(gFexSRJEtValues[1]);
                gFexSRJSubleadingEtaValues.push_back(subleading->eta());
                gFexSRJSubleadingPhiValues.push_back(subleading->phi());
            }


            // Temporary vector for sorting by Et
            std::vector<std::pair<size_t, double>> gFexLRJetEtWithIndex;

            for (size_t i = 0; i < L1_gFexLRJetRoI->size(); ++i) {
                const auto& jet = (*L1_gFexLRJetRoI)[i];
                double et = jet->et() / 1000.0; // already in GeV
                gFexLRJetEtWithIndex.emplace_back(i, et);
            }

            // Sort descending by Et
            std::sort(gFexLRJetEtWithIndex.begin(), gFexLRJetEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });

            // Fill vectors in sorted order
            for (const auto& [index, et] : gFexLRJetEtWithIndex) {
                const auto& jet = (*L1_gFexLRJetRoI)[index];

                gFexLRJEtIndexValues.push_back(static_cast<unsigned int>(index));
                gFexLRJEtValues.push_back(et);
                gFexLRJEtaValues.push_back(jet->eta());
                gFexLRJPhiValues.push_back(jet->phi());
            }

            // Leading jet
            if (!gFexLRJetEtWithIndex.empty()) {
                const auto& leading = (*L1_gFexLRJetRoI)[gFexLRJetEtWithIndex[0].first];
                gFexLRJLeadingEtValues.push_back(gFexLRJEtValues[0]);
                gFexLRJLeadingEtaValues.push_back(leading->eta());
                gFexLRJLeadingPhiValues.push_back(leading->phi());
            }

            // Subleading jet
            if (gFexLRJetEtWithIndex.size() > 1) {
                const auto& subleading = (*L1_gFexLRJetRoI)[gFexLRJetEtWithIndex[1].first];
                gFexLRJSubleadingEtValues.push_back(gFexLRJEtValues[1]);
                gFexLRJSubleadingEtaValues.push_back(subleading->eta());
                gFexLRJSubleadingPhiValues.push_back(subleading->phi());
            }


            std::vector<std::pair<size_t, double>> hltJetEtWithIndex;
            for (size_t i = 0; i < HLT_AntiKt4EMTopoJets_subjesIS->size(); ++i) {
                const auto& el = (*HLT_AntiKt4EMTopoJets_subjesIS)[i];
                double et = el->e() / (1000.0 * cosh(el->eta()));
                hltJetEtWithIndex.emplace_back(i, et);  // Store index and Et for sorting
            }
            // Sort by descending Et
            std::sort(hltJetEtWithIndex.begin(), hltJetEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });
            // Now push back into vectors in sorted order
            for (const auto& [index, et] : hltJetEtWithIndex) {
                const auto& el = (*HLT_AntiKt4EMTopoJets_subjesIS)[index];
                hltAntiKt4SRJEtIndexValues.push_back(static_cast<unsigned int>(index));
                hltAntiKt4SRJEtValues.push_back(et);
                hltAntiKt4SRJEtaValues.push_back(el->eta());
                hltAntiKt4SRJPhiValues.push_back(el->phi());
            }

            // Store leading and subleading jets if available
            if (!hltJetEtWithIndex.empty()) {
                const auto& leading = (*HLT_AntiKt4EMTopoJets_subjesIS)[hltJetEtWithIndex[0].first];
                hltAntiKt4SRJLeadingEtValues.push_back(hltAntiKt4SRJEtValues[0]);
                hltAntiKt4SRJLeadingEtaValues.push_back(leading->eta());
                hltAntiKt4SRJLeadingPhiValues.push_back(leading->phi());
            }

            if (hltJetEtWithIndex.size() > 1) {
                const auto& subleading = (*HLT_AntiKt4EMTopoJets_subjesIS)[hltJetEtWithIndex[1].first];
                hltAntiKt4SRJSubleadingEtValues.push_back(hltAntiKt4SRJEtValues[1]);
                hltAntiKt4SRJSubleadingEtaValues.push_back(subleading->eta());
                hltAntiKt4SRJSubleadingPhiValues.push_back(subleading->phi());
            }

            // --- Loop over L1_gFexLRJetRoI ---
            //for (const auto& el : *gFexLRJets) {
            //    gfex_larger_jet_pt_values.push_back(el->et() / 1000.0);
            //    float gfex_larger_jet_Et = el->et() / 1000.0;
            //}

            // --- Loop over InTimeAntiKt4TruthJets ---
            // Temporary vector to hold (index, Et) for sorting
            std::vector<std::pair<size_t, double>> truthJetEtWithIndex;

            for (size_t i = 0; i < InTimeAntiKt4TruthJets->size(); ++i) {
                const auto& jet = (*InTimeAntiKt4TruthJets)[i];
                double et = jet->e() / (1000.0 * cosh(jet->eta()));
                truthJetEtWithIndex.emplace_back(i, et);  // Store original index and Et
            }

            // Sort by descending Et
            std::sort(truthJetEtWithIndex.begin(), truthJetEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });

            // Fill vectors in sorted order
            for (const auto& [index, et] : truthJetEtWithIndex) {
                const auto& jet = (*InTimeAntiKt4TruthJets)[index];
                inTimeAntiKt4TruthSRJEtIndexValues.push_back(static_cast<unsigned int>(index));
                inTimeAntiKt4TruthSRJEtValues.push_back(et);
                inTimeAntiKt4TruthSRJEtaValues.push_back(jet->eta());
                inTimeAntiKt4TruthSRJPhiValues.push_back(jet->phi());
            }

            // Store leading jet info if available
            if (!truthJetEtWithIndex.empty()) {
                const auto& leading = (*InTimeAntiKt4TruthJets)[truthJetEtWithIndex[0].first];
                inTimeAntiKt4TruthSRJLeadingEtValues.push_back(inTimeAntiKt4TruthSRJEtValues[0]);
                inTimeAntiKt4TruthSRJLeadingEtaValues.push_back(leading->eta());
                inTimeAntiKt4TruthSRJLeadingPhiValues.push_back(leading->phi());
            }

            // Store subleading jet info if available
            if (truthJetEtWithIndex.size() > 1) {
                const auto& subleading = (*InTimeAntiKt4TruthJets)[truthJetEtWithIndex[1].first];
                inTimeAntiKt4TruthSRJSubleadingEtValues.push_back(inTimeAntiKt4TruthSRJEtValues[1]);
                inTimeAntiKt4TruthSRJSubleadingEtaValues.push_back(subleading->eta());
                inTimeAntiKt4TruthSRJSubleadingPhiValues.push_back(subleading->phi());
            }

            
            if(higgsPtCutsPassed){
                truthbTree->Fill();
                truthHiggsTree->Fill();
                // truthVBFQuark->Fill();  // commented out as in your declaration
                caloTopoTowerTree->Fill();
                topo422Tree->Fill();
                gFexSRJTree->Fill();
                gFexLeadingSRJTree->Fill();
                gFexSubleadingSRJTree->Fill();
                gFexLRJTree->Fill();
                gFexLeadingLRJTree->Fill();
                gFexSubleadingLRJTree->Fill();
                inTimeAntiKt4TruthJetsTree->Fill();
                leadingInTimeAntiKt4TruthJetsTree->Fill();
                subleadingInTimeAntiKt4TruthJetsTree->Fill();
                jFexSRJTree->Fill();
                jFexLeadingSRJTree->Fill();
                jFexSubleadingSRJTree->Fill();
                jFexLRJTree->Fill();
                jFexLeadingLRJTree->Fill();
                jFexSubleadingLRJTree->Fill();
                hltAntiKt4EMTopoJetsTree->Fill();
                leadingHltAntiKt4EMTopoJetsTree->Fill();
                subleadingHltAntiKt4EMTopoJetsTree->Fill();
                recoAntiKt10UFOCSSKJets->Fill();
                leadingRecoAntiKt10UFOCSSKJets->Fill();
                subleadingRecoAntiKt10UFOCSSKJets->Fill();
                truthAntiKt4TruthDressedWZJets->Fill();
                leadingTruthAntiKt4TruthDressedWZJets->Fill();
                subleadingTruthAntiKt4TruthDressedWZJets->Fill();
            }
        } // loop through events
        
        f->Close();
    } // loop through filenames
    outputFile->cd();
    truthbTree->Write("", TObject::kOverwrite);
    truthHiggsTree->Write("", TObject::kOverwrite);
    // truthVBFQuark->Write();  // Optional, if used
    caloTopoTowerTree->Write("", TObject::kOverwrite);
    topo422Tree->Write("", TObject::kOverwrite);
    gFexSRJTree->Write("", TObject::kOverwrite);
    gFexLeadingSRJTree->Write("", TObject::kOverwrite);
    gFexSubleadingSRJTree->Write("", TObject::kOverwrite);
    gFexLRJTree->Write("", TObject::kOverwrite);
    gFexLeadingLRJTree->Write("", TObject::kOverwrite);
    gFexSubleadingLRJTree->Write("", TObject::kOverwrite);
    inTimeAntiKt4TruthJetsTree->Write("", TObject::kOverwrite);
    leadingInTimeAntiKt4TruthJetsTree->Write("", TObject::kOverwrite);
    subleadingInTimeAntiKt4TruthJetsTree->Write("", TObject::kOverwrite);
    jFexSRJTree->Write("", TObject::kOverwrite);
    jFexLeadingSRJTree->Write("", TObject::kOverwrite);
    jFexSubleadingSRJTree->Write("", TObject::kOverwrite);
    jFexLRJTree->Write("", TObject::kOverwrite);
    jFexLeadingLRJTree->Write("", TObject::kOverwrite);
    jFexSubleadingLRJTree->Write("", TObject::kOverwrite);
    hltAntiKt4EMTopoJetsTree->Write("", TObject::kOverwrite);
    leadingHltAntiKt4EMTopoJetsTree->Write("", TObject::kOverwrite);
    subleadingHltAntiKt4EMTopoJetsTree->Write("", TObject::kOverwrite);
    recoAntiKt10UFOCSSKJets->Write("", TObject::kOverwrite);
    leadingRecoAntiKt10UFOCSSKJets->Write("", TObject::kOverwrite);
    subleadingRecoAntiKt10UFOCSSKJets->Write("", TObject::kOverwrite);
    truthAntiKt4TruthDressedWZJets->Write("", TObject::kOverwrite);
    leadingTruthAntiKt4TruthDressedWZJets->Write("", TObject::kOverwrite);
    subleadingTruthAntiKt4TruthDressedWZJets->Write("", TObject::kOverwrite);
    outputFile->Close();
    std::cout << "Processing complete." << endl;
    std::cout << "higgsPassEventCounter: " << higgsPassEventCounter << "\n";
} // ntupler function
