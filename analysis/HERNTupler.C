// To execute: e.g., root ; .L nTupler.C ; nTupler(true, true, true) 
#include <algorithm>
#include <numeric>   // std::iota

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


// Ntupler functions
// Used for digitized data writing
struct OutputFiles {
    std::string topo422;
    std::string caloTopoTowers;
    std::string gepBasicClusters;
    std::string gepBasicClustersSK;
    std::string gepCellTowersSK; 
    std::string gepCellTowers;
    std::string gepCellTowersSort;
    std::string gepWTAConeJetsCellsTowers;
    std::string gepWTAConeJetsBasicClusters;
    std::string gepWTAConeJetsCellsTowersSK;
    std::string gepWTAConeJetsBasicClustersSK;
    std::string gFex;
    std::string jFex;
};

const unsigned int nJZSlices = 10;

// Constants required for reweighting

// In barns^-1 - 7.5*10^34 cm^-2 s^-1 * 1 s (HL-LHC 200 PU inst. lumi * 1 second) - use 1 second to make rates plots easy
const double reweightLuminosity = 7.5e10;

// Filter efficiencies from AMI
const double filterEffienciesByJZSlice[nJZSlices] = {0.9716436,    // JZ0
                                                     0.03777559,   // JZ1
                                                     0.01136654,   // JZ2
                                                     0.01367042,   // JZ3
                                                     0.01628158,   // JZ4
                                                     0.01905588,   // JZ5
                                                     0.01352844,   // JZ6
                                                     0.01764909,   // JZ7
                                                     0.01887484,   // JZ8
                                                     0.02827565};  // JZ9

// Cross sections from AMI [in b]
const double crossSectionsByJZSlice[nJZSlices] = {0.07893,      // JZ0
                                                  0.09679,      // JZ1
                                                  0.0026805,    // JZ2
                                                  0.000029984,  // JZ3
                                                  2.972e-7,     // JZ4
                                                  5.5384e-09,   // JZ5
                                                  3.2616e-10,   // JZ6
                                                  2.1734e-11,   // JZ7
                                                  9.2995e-13,   // JZ8
                                                  3.4519e-14};  // JZ9

// From CutBookkeeper printed from python script
const double sumOfEventWeightsByJZSlice[nJZSlices] = {10000.0,                // JZ0
                                                      467.7515499125784,      // JZ1
                                                      3.9288560244923474,     // JZ2
                                                      0.08040490627634789,    // JZ3
                                                      0.0025011268072319126,  // JZ4
                                                      0.00019515390387392362, // JZ5
                                                      3.325120293798389e-05,  // JZ6
                                                      9.518299562145949e-06,  // JZ7 
                                                      2.8526936132339253e-06, // JZ8
                                                      3.274964361399163e-09}; // JZ9 

// Helper to construct memprint / test vector filenames
inline OutputFiles makeMemPrintFilenames(std::string signalString, bool signalBool, int jzSlice, unsigned int algoVersion) {
    std::string base = "/data/larsonma/LargeRadiusJets/MemPrints";
    if(algoVersion == 2) base = "/data/larsonma/LargeRadiusJets/MemPrints_v2/";
    else if(algoVersion == 3) base = "/data/larsonma/LargeRadiusJets/MemPrints_v3/";
     
    std::cout << "BASE BEING USED FOR MEMPRINTS: " << base << "\n";
    std::string tag;
    if (signalBool) {
        if(signalString == "VBF_hh_bbbb") tag = "mc21_14TeV_hh_bbbb_vbf_novhh"; 
        else if (signalString == "ggF_hh_bbbb") tag = "mc21_14TeV_HHbbbb_HLLHC";
        else if (signalString == "ZvvHbb") tag = "mc21_14TeV_ZvvH125_bb";
        else if (signalString == "ttbar_had") tag = "mc21_14TeV_ttbar_hdamp258p75_allhad";
        else if (signalString == "Zprime_ttbar") tag = "mc21_14TeV_flatpT_Zprime_tthad";
    } else {
        // generic: works for all slices 0..9 etc
        tag = "mc21_14TeV_jj_JZ" + std::to_string(jzSlice);
    }

    OutputFiles out;
    out.topo422         = base + "CaloTopo_422/"    + tag + "_topo422.dat";
    out.caloTopoTowers  = base + "CaloTopoTowers/"  + tag + "_calotopotowers.dat";
    out.gepBasicClusters = base + "GEPBasicClusters/"+ tag + "_gepbasicclusters.dat";
    out.gepBasicClustersSK = base + "GEPBasicClustersSK/"+ tag + "_gepbasicclusters.dat";
    out.gepCellTowers = base + "GEPCellsTowers/"+ tag + "_gepcellstowers.dat";
    out.gepCellTowersSK = base + "GEPCellsTowersSK/"+ tag + "_gepcellstowerssk.dat";
    out.gepCellTowersSort =base + "GEPCellsTowers_Sort/" + tag + "_gepcellstowers.dat";
    out.gepWTAConeJetsCellsTowers = base + "GEPConeJetsCellsTowers/"+ tag + "_gepconejetscellstowers.dat";
    out.gepWTAConeJetsBasicClusters = base + "GEPConeJetsBasicClusters/"+ tag + "_gepconejetsbasicclusters.dat";
    out.gepWTAConeJetsCellsTowersSK = base + "GEPConeJetsCellsTowersSK/"+ tag + "_gepconejetscellstowerssk.dat";
    out.gepWTAConeJetsBasicClustersSK = base + "GEPConeJetsBasicClustersSK/"+ tag + "_gepconejetsbasicclusterssk.dat";
    out.gFex            = base + "gFex/"            + tag + "_gfex_smallrj.dat";
    out.jFex            = base + "jFex/"            + tag + "_jfex_smallrj.dat";
    return out;
}

// Helper to get input DAOD + GEP ntuple file names
std::string makeFileDir(std::string signalString, bool signalBool, int jzSlice) {
    static const std::string kRoot = "/data/larsonma/LargeRadiusJets/datasets";

    if (signalBool) {
        if(signalString == "VBF_hh_bbbb") return kRoot + "/DAOD_TrigGepPerf/Signal_HHbbbb_VBF_SM"; 
        else if (signalString == "ggF_hh_bbbb") return kRoot + "/DAOD_TrigGepPerf/Signal_HHbbbb_ggF";
        else if (signalString == "ZvvHbb") return kRoot + "/DAOD_TrigGepPerf/Signal_ZvvHbb";
        else if (signalString == "ttbar_had") return kRoot + "/DAOD_TrigGepPerf/Signal_ttbar_fullhad";
        else if (signalString == "Zprime_ttbar") return kRoot + "/DAOD_TrigGepPerf/Signal_Zprime_ttbar_flatpT";
    }

    // background path
    if (jzSlice < 0 || jzSlice > 9) {
        throw std::out_of_range("jzSlice must be in [0, 9]");
    }
    return kRoot + "/DAOD_TrigGepPerf/Background_jj_JZ" + std::to_string(jzSlice);
}

// Helper function to get output ntuple file name
TString makeOutputFileName(std::string signalString, bool signalBool, int jzSlice, unsigned int algoVersion) {
    TString outDir;
    if(algoVersion == 2) outDir = "outputRootFiles/v2/";
    else if(algoVersion == 3) outDir = "outputRootFiles/v3/";
   
    if(signalBool){
        if(signalString == "VBF_hh_bbbb") return outDir + "mc21_14TeV_hh_bbbb_vbf_novhh_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ggF_hh_bbbb") return outDir + "mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ZvvHbb") return outDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ttbar_had") return outDir + "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "Zprime_ttbar") return outDir + "mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else return "";
    }
    else {
        if (jzSlice < 0 || jzSlice > 9)
            throw std::out_of_range("jzSlice must be between 0 and 9");
        return outDir + "mc21_14TeV_jj_JZ" + TString(std::to_string(jzSlice))
            + "_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    }
 
}

unsigned int digitize(double value, int bit_length, double min_val, double max_val, unsigned int altRange = 0) {
    unsigned int range = (altRange == 0) ? (1u << bit_length) : altRange;
    double scale = double(range) / (max_val - min_val);
    //std::cout << "max_val - scale; " << max_val - (1/scale) << "\n";
    // Check if value is in range
    if (value < min_val) {
        value = min_val;
        //std::cout << "Warning: Value " << value
        //  << " is out of range (" << min_val
        //  << ", " << max_val << ")\n";
    }
    if (value >= max_val){
        return range - 1;
    }

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
void nTupler(bool signalBool, std::string signalString, unsigned int etaAltRange, unsigned int algoVersion, unsigned int jzSlice = 3) {
    //ServiceHandle<StoreGateSvc> inputMetaStore("StoreGateSvc/InputMetaDataStore","");
    // Setup file paths based on whether processing signal or background, and vbf production or ggF production
    std::string fileDir =  makeFileDir(signalString, signalBool, jzSlice);
    std::cout << "fileDir: " << fileDir << "\n";
    //xAOD::Init().ignore();

    // pick filenames
    OutputFiles out = makeMemPrintFilenames(signalString, signalBool, jzSlice, algoVersion);

    // open streams
    std::ofstream f_topotower(out.caloTopoTowers);
    std::ofstream f_gepbasicclusters(out.gepBasicClusters);
    std::ofstream f_gepcellstowers(out.gepCellTowers);
    std::ofstream f_gepbasicclusterssk(out.gepBasicClustersSK);
    std::ofstream f_gepcellstowerssk(out.gepCellTowersSK);
    std::ofstream f_gepcellstowers_sort(out.gepCellTowersSort);
    std::ofstream f_wtaconejetscellstowers(out.gepWTAConeJetsCellsTowers);
    std::ofstream f_wtaconejetsbasicclusters(out.gepWTAConeJetsBasicClusters);
    std::ofstream f_wtaconejetscellstowerssk(out.gepWTAConeJetsCellsTowersSK);
    std::ofstream f_wtaconejetsbasicclusterssk(out.gepWTAConeJetsBasicClustersSK);
    std::ofstream f_topo(out.topo422);
    std::ofstream f_gfex(out.gFex);
    std::ofstream f_jfex(out.jFex);

    if (!f_topotower.is_open() || !f_topo.is_open() || !f_gfex.is_open() || !f_jfex.is_open() || !f_gepbasicclusters.is_open() || !f_gepcellstowers.is_open() || !f_gepbasicclusterssk.is_open() || !f_gepcellstowerssk.is_open() || !f_wtaconejetscellstowers.is_open() || !f_wtaconejetsbasicclusters.is_open() || !f_wtaconejetscellstowerssk.is_open() || !f_wtaconejetsbasicclusterssk.is_open()) {
        std::cerr << "Error: One or more output files could not be opened for writing!" << std::endl;
        // handle error or return
    }

    TString outputFileName = makeOutputFileName(signalString, signalBool, jzSlice, algoVersion);
    
    // Create ROOT output file
    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    // Create TTrees - commented out trees that cannot be filled with DAOD information.
    //TTree* truthParticleTree = new TTree("truthParticleTree", "Tree storing event-wise information about truth particles");
    TTree* eventInfoTree = new TTree("eventInfoTree", "Tree storing event information (e.g., mcEventWeight) required for rate computation");
    TTree* truthbTree = new TTree("truthbTree", "Tree storing event-wise information about truth particles");
    TTree* truthHiggsTree = new TTree("truthHiggsTree", "Tree storing event-wise information about truth particles");
    TTree* truthTopTree = new TTree("truthTopTree", "Tree storing event-wise information about truth particles");
    //TTree* truthVBFQuark = new TTree("truthVBFQuark", "Tree storing event-wise information about truth particles");
    TTree* caloTopoTowerTree = new TTree("caloTopoTowerTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gepBasicClustersTree = new TTree("gepBasicClustersTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gepBasicClustersSKTree = new TTree("gepBasicClustersSKTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gepCellsTowersTree = new TTree("gepCellsTowersTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gepCellsTowersSKTree = new TTree("gepCellsTowersSKTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gepWTAConeCellsTowersJetsTree = new TTree("gepWTAConeCellsTowersJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepLeadingWTAConeCellsTowersJetsTree = new TTree("gepLeadingWTAConeCellsTowersJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepSubleadingWTAConeCellsTowersJetsTree = new TTree("gepSubleadingWTAConeCellsTowersJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepWTAConeCellsTowersSKJetsTree = new TTree("gepWTAConeCellsTowersSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepLeadingWTAConeCellsTowersSKJetsTree = new TTree("gepLeadingWTAConeCellsTowersSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepSubleadingWTAConeCellsTowersSKJetsTree = new TTree("gepSubleadingWTAConeCellsTowersSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepWTAConeBasicClustersJetsTree = new TTree("gepWTAConeBasicClustersJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepLeadingWTAConeBasicClustersJetsTree = new TTree("gepLeadingWTAConeBasicClustersJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepSubleadingWTAConeBasicClustersJetsTree = new TTree("gepSubleadingWTAConeBasicClustersJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepWTAConeBasicClustersSKJetsTree = new TTree("gepWTAConeBasicClustersSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepLeadingWTAConeBasicClustersSKJetsTree = new TTree("gepLeadingWTAConeBasicClustersSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepSubleadingWTAConeBasicClustersSKJetsTree = new TTree("gepSubleadingWTAConeBasicClustersSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
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
    TTree* recoAntiKt10UFOCSSKSoftDropJets = new TTree("recoAntiKt10UFOCSSKSoftDropJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* leadingRecoAntiKt10UFOCSSKSoftDropJets = new TTree("leadingRecoAntiKt10UFOCSSKSoftDropJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* subleadingRecoAntiKt10UFOCSSKSoftDropJets = new TTree("subleadingRecoAntiKt10UFOCSSKSoftDropJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* antiKt10TruthJetsTree = new TTree("antiKt10TruthJetsTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* leadingAntiKt10TruthJetsTree = new TTree("leadingAntiKt10TruthJetsTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* subleadingAntiKt10TruthJetsTree = new TTree("subleadingAntiKt10TruthJetsTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* antiKt10TruthSoftDropJetsTree = new TTree("antiKt10TruthSoftDropJetsTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* leadingAntiKt10TruthSoftDropJetsTree = new TTree("leadingAntiKt10TruthSoftDropJetsTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* subleadingAntiKt10TruthSoftDropJetsTree = new TTree("subleadingAntiKt10TruthSoftDropJetsTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* truthAntiKt4TruthDressedWZJets = new TTree("truthAntiKt4TruthDressedWZJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* leadingTruthAntiKt4TruthDressedWZJets = new TTree("leadingTruthAntiKt4TruthDressedWZJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* subleadingTruthAntiKt4TruthDressedWZJets = new TTree("subleadingTruthAntiKt4TruthDressedWZJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* gFexMETTree             = new TTree("gFexMETTree",             "Tree storing event-wise gFEX MET quantities");
    TTree* gFexMETNoiseCutTree     = new TTree("gFexMETNoiseCutTree",     "Tree storing event-wise gFEX MET quantities (NoiseCut)");
    TTree* gFexMETRmsTree          = new TTree("gFexMETRmsTree",          "Tree storing event-wise gFEX MET quantities (Rms)");
    TTree* metTruthTree            = new TTree("metTruthTree",            "Tree storing event-wise truth MET quantities");
    TTree* metCoreAntiKt4EMTopoTree  = new TTree("metCoreAntiKt4EMTopoTree",  "Tree storing event-wise reco core MET quantities (AntiKt4EMTopo)");
    TTree* metCoreAntiKt4EMPFlowTree = new TTree("metCoreAntiKt4EMPFlowTree", "Tree storing event-wise reco core MET quantities (AntiKt4EMPFlow)");
    /*TTree* recoAntiKt10LCTopoJets = new TTree("recoAntiKt10LCTopoJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* leadingRecoAntiKt10LCTopoJets = new TTree("leadingRecoAntiKt10LCTopoJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* subleadingRecoAntiKt10LCTopoJets = new TTree("recoAntiKt10LCTopoJets", "Tree storing event-wise Et, Eta, Phi");*/
    //TTree* gFexRhoRoI = new TTree("gFexRhoRoI", "Tree storing event-wise Et, Eta, Phi");

    // Variables to store data to be written to TTrees

    // Event info vectors
    std::vector<double> mcEventWeights;
    double sumOfWeightsForSample;
    int sampleJZSlice;
    std::vector<double> eventWeights;

    // Truth Particle vectors
    //std::vector<int> truthParticlePDGId, truthParticleStatus; 
    //std::vector<double> truthParticleEtValues, truthParticleEnergyValues, truthParticlepTValues, truthParticlepxValues, truthParticlepyValues, truthParticlepzValues, truthParticleEtaValues, truthParticlePhiValues;
    std::vector<unsigned int> higgsIndexValues, indexOfHiggsValues;
    std::vector<double> truthbquarksEtValues, truthbquarksEnergyValues, truthbquarkspTValues, truthbquarkspxValues, truthbquarkspyValues, truthbquarkspzValues, truthbquarksEtaValues, truthbquarksPhiValues;
    std::vector<double> truthHiggsEtValues, truthHiggsEnergyValues, truthHiggspTValues, truthHiggspxValues, truthHiggspyValues, truthHiggspzValues, truthHiggsEtaValues, truthHiggsPhiValues, truthHiggsInvMassValues;
    //std::vector<double> truthVBFQuarkValues, truthVBFQuarkEnergyValues, truthVBFQuarkpTValues, truthVBFQuarkpxValues, truthVBFQuarkpyValues, truthVBFQuarkpzValues, truthVBFQuarkEtaValues, truthVBFQuarkPhiValues;

    std::vector<unsigned int> topIndexValues, indexOfTopValues;
    std::vector<double> truthTopEtValues, truthTopEnergyValues, truthToppTValues, truthTopEtaValues, truthTopPhiValues;
    // Tower / cluster vectors
    std::vector<double> caloTopoTowerEtValues, caloTopoTowerEtaValues, caloTopoTowerPhiValues;
    std::vector<double> gepBasicClustersEtValues, gepBasicClustersEtaValues, gepBasicClustersPhiValues;
    std::vector<double> gepBasicClustersSKEtValues, gepBasicClustersSKEtaValues, gepBasicClustersSKPhiValues;
    std::vector<double> gepCellsTowersEtValues, gepCellsTowersEtaValues, gepCellsTowersPhiValues;
    std::vector<double> gepCellsTowersSKEtValues, gepCellsTowersSKEtaValues, gepCellsTowersSKPhiValues;
    std::vector<double> topo422EtValues, topo422EtaValues, topo422PhiValues;

    // Cone jets from TrigGepPerf
    std::vector<double> gepWTAConeCellsTowersJetspTValues, gepWTAConeCellsTowersJetsEtaValues, gepWTAConeCellsTowersJetsPhiValues;
    std::vector<unsigned int > gepWTAConeCellsTowersJetsNConstituentsValues;
    std::vector<double> gepWTAConeGEPBasicClustersJetspTValues, gepWTAConeGEPBasicClustersJetsEtaValues, gepWTAConeGEPBasicClustersJetsPhiValues;
    std::vector<unsigned int > gepWTAConeGEPBasicClustersJetsNConstituentsValues;

    // SK PU-suppressed cone jets from TrigGepPerf
    std::vector<double> gepWTAConeCellsTowersSKJetspTValues, gepWTAConeCellsTowersSKJetsEtaValues, gepWTAConeCellsTowersSKJetsPhiValues;
    std::vector<unsigned int > gepWTAConeCellsTowersSKJetsNConstituentsValues;
    std::vector<double> gepWTAConeGEPBasicClustersSKJetspTValues, gepWTAConeGEPBasicClustersSKJetsEtaValues, gepWTAConeGEPBasicClustersSKJetsPhiValues;
    std::vector<unsigned int > gepWTAConeGEPBasicClustersSKJetsNConstituentsValues;

    // Leading cone jets from TrigGepPerf
    std::vector<double> gepLeadingWTAConeCellsTowersJetspTValues, gepLeadingWTAConeCellsTowersJetsEtaValues, gepLeadingWTAConeCellsTowersJetsPhiValues;
    std::vector<unsigned int > gepLeadingWTAConeCellsTowersJetsNConstituentsValues;
    std::vector<double> gepLeadingWTAConeGEPBasicClustersJetspTValues, gepLeadingWTAConeGEPBasicClustersJetsEtaValues, gepLeadingWTAConeGEPBasicClustersJetsPhiValues;
    std::vector<unsigned int > gepLeadingWTAConeGEPBasicClustersJetsNConstituentsValues;

    // Leading SK PU-suppressed jets from TrigGepPerf
    std::vector<double> gepLeadingWTAConeCellsTowersSKJetspTValues, gepLeadingWTAConeCellsTowersSKJetsEtaValues, gepLeadingWTAConeCellsTowersSKJetsPhiValues;
    std::vector<unsigned int > gepLeadingWTAConeCellsTowersSKJetsNConstituentsValues;
    std::vector<double> gepLeadingWTAConeGEPBasicClustersSKJetspTValues, gepLeadingWTAConeGEPBasicClustersSKJetsEtaValues, gepLeadingWTAConeGEPBasicClustersSKJetsPhiValues;
    std::vector<unsigned int > gepLeadingWTAConeGEPBasicClustersSKJetsNConstituentsValues;

    // Subleading cone jets from TrigGepPerf
    std::vector<double> gepSubleadingWTAConeCellsTowersJetspTValues, gepSubleadingWTAConeCellsTowersJetsEtaValues, gepSubleadingWTAConeCellsTowersJetsPhiValues;
    std::vector<unsigned int > gepSubleadingWTAConeCellsTowersJetsNConstituentsValues;
    std::vector<double> gepSubleadingWTAConeGEPBasicClustersJetspTValues, gepSubleadingWTAConeGEPBasicClustersJetsEtaValues, gepSubleadingWTAConeGEPBasicClustersJetsPhiValues;
    std::vector<unsigned int > gepSubleadingWTAConeGEPBasicClustersJetsNConstituentsValues;

    // Subleading SK PU-suppressed cone jets from TrigGepPerf
    std::vector<double> gepSubleadingWTAConeCellsTowersSKJetspTValues, gepSubleadingWTAConeCellsTowersSKJetsEtaValues, gepSubleadingWTAConeCellsTowersSKJetsPhiValues;
    std::vector<unsigned int > gepSubleadingWTAConeCellsTowersSKJetsNConstituentsValues;
    std::vector<double> gepSubleadingWTAConeGEPBasicClustersSKJetspTValues, gepSubleadingWTAConeGEPBasicClustersSKJetsEtaValues, gepSubleadingWTAConeGEPBasicClustersSKJetsPhiValues;
    std::vector<unsigned int > gepSubleadingWTAConeGEPBasicClustersSKJetsNConstituentsValues;

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
    std::vector<double> recoAntiKt10LRJEtValues, recoAntiKt10LRJEtaValues, recoAntiKt10LRJPhiValues, recoAntiKt10LRJMassValues;
    std::vector<double> recoAntiKt10LRJLeadingEtValues, recoAntiKt10LRJLeadingEtaValues, recoAntiKt10LRJLeadingPhiValues, recoAntiKt10LRJLeadingMassValues;
    std::vector<double> recoAntiKt10LRJSubleadingEtValues, recoAntiKt10LRJSubleadingEtaValues, recoAntiKt10LRJSubleadingPhiValues, recoAntiKt10LRJSubleadingMassValues;

    // Reco AntiKt10 UFOCSSK SoftDrop jets vectors
    std::vector<unsigned int> recoAntiKt10UFOCSSKSDEtIndexValues;
    std::vector<double> recoAntiKt10UFOCSSKSDEtValues, recoAntiKt10UFOCSSKSDEtaValues, recoAntiKt10UFOCSSKSDPhiValues, recoAntiKt10UFOCSSKSDMassValues;
    std::vector<double> recoAntiKt10UFOCSSKSDLeadingEtValues, recoAntiKt10UFOCSSKSDLeadingEtaValues, recoAntiKt10UFOCSSKSDLeadingPhiValues, recoAntiKt10UFOCSSKSDLeadingMassValues;
    std::vector<double> recoAntiKt10UFOCSSKSDSubleadingEtValues, recoAntiKt10UFOCSSKSDSubleadingEtaValues, recoAntiKt10UFOCSSKSDSubleadingPhiValues, recoAntiKt10UFOCSSKSDSubleadingMassValues;

    // AntiKt10 Truth jets vectors
    std::vector<unsigned int> antiKt10TruthEtIndexValues;
    std::vector<double> antiKt10TruthEtValues, antiKt10TruthEtaValues, antiKt10TruthPhiValues, antiKt10TruthMassValues;
    std::vector<double> antiKt10TruthLeadingEtValues, antiKt10TruthLeadingEtaValues, antiKt10TruthLeadingPhiValues, antiKt10TruthLeadingMassValues;
    std::vector<double> antiKt10TruthSubleadingEtValues, antiKt10TruthSubleadingEtaValues, antiKt10TruthSubleadingPhiValues, antiKt10TruthSubleadingMassValues;

    // AntiKt10 Truth SoftDrop jets vectors
    std::vector<unsigned int> antiKt10TruthSDEtIndexValues;
    std::vector<double> antiKt10TruthSDEtValues, antiKt10TruthSDEtaValues, antiKt10TruthSDPhiValues, antiKt10TruthSDMassValues;
    std::vector<double> antiKt10TruthSDLeadingEtValues, antiKt10TruthSDLeadingEtaValues, antiKt10TruthSDLeadingPhiValues, antiKt10TruthSDLeadingMassValues;
    std::vector<double> antiKt10TruthSDSubleadingEtValues, antiKt10TruthSDSubleadingEtaValues, antiKt10TruthSDSubleadingPhiValues, antiKt10TruthSDSubleadingMassValues;

    // Truth WZ Antikt4 jets vectors
    std::vector<unsigned int> truthAntiKt4WZSRJEtIndexValues;
    std::vector<double> truthAntiKt4WZSRJEtValues, truthAntiKt4WZSRJEtaValues, truthAntiKt4WZSRJPhiValues, truthAntiKt4WZSRJMassValues;
    std::vector<double> truthAntiKt4WZSRJLeadingEtValues, truthAntiKt4WZSRJLeadingEtaValues, truthAntiKt4WZSRJLeadingPhiValues, truthAntiKt4WZSRJLeadingMassValues;
    std::vector<double> truthAntiKt4WZSRJSubleadingEtValues, truthAntiKt4WZSRJSubleadingEtaValues, truthAntiKt4WZSRJSubleadingPhiValues, truthAntiKt4WZSRJSubleadingMassValues;

    // gFEX MET quantities (scalar per event)
    double gMHX = 0.0, gMHY = 0.0;                 // hard term (L1_gMHTComponentsJwoj)
    double gMSX = 0.0, gMSY = 0.0;                 // soft term (L1_gMSTComponentsJwoj)
    double gMETX = 0.0, gMETY = 0.0;               // total MET components (L1_gMETComponentsJwoj)
    double gMET = 0.0, gSumET = 0.0;               // scalar MET and SumET (L1_gScalarEJwoj)
    double gMETX_NC = 0.0, gMETY_NC = 0.0;         // total MET components (L1_gMETComponentsNoiseCut)
    double gMET_NC = 0.0, gSumET_NC = 0.0;         // scalar MET and SumET (L1_gScalarENoiseCut)
    double gMETX_Rms = 0.0, gMETY_Rms = 0.0;       // total MET components (L1_gMETComponentsRms)
    double gMET_Rms = 0.0, gSumET_Rms = 0.0;       // scalar MET and SumET (L1_gScalarERms)
    // Reference MET quantities (scalar per event)
    double metTruthNonIntX = 0.0, metTruthNonIntY = 0.0;   // MET_Truth["NonInt"]
    double metTruthNonIntSumET = 0.0, metTruthNonInt = 0.0;
    double metTruthIntX = 0.0, metTruthIntY = 0.0;         // MET_Truth["Int"]
    double metTruthIntSumET = 0.0, metTruthInt = 0.0;
    double metTruthIntOutX = 0.0, metTruthIntOutY = 0.0;   // MET_Truth["IntOut"]
    double metTruthIntOutSumET = 0.0, metTruthIntOut = 0.0;
    // MET_Core_AntiKt4EMTopo — one set of variables per term
    double metCoreEMTopo_SoftClus_X = 0.0,   metCoreEMTopo_SoftClus_Y = 0.0,   metCoreEMTopo_SoftClus_SumET = 0.0,   metCoreEMTopo_SoftClus_MET = 0.0;
    double metCoreEMTopo_PVSoftTrk_X = 0.0,  metCoreEMTopo_PVSoftTrk_Y = 0.0,  metCoreEMTopo_PVSoftTrk_SumET = 0.0,  metCoreEMTopo_PVSoftTrk_MET = 0.0;
    double metCoreEMTopo_SoftClusEM_X = 0.0, metCoreEMTopo_SoftClusEM_Y = 0.0, metCoreEMTopo_SoftClusEM_SumET = 0.0, metCoreEMTopo_SoftClusEM_MET = 0.0;
    // MET_Core_AntiKt4EMPFlow — one set of variables per term
    double metCoreEMPFlow_SoftClus_X = 0.0,  metCoreEMPFlow_SoftClus_Y = 0.0,  metCoreEMPFlow_SoftClus_SumET = 0.0,  metCoreEMPFlow_SoftClus_MET = 0.0;
    double metCoreEMPFlow_PVSoftTrk_X = 0.0, metCoreEMPFlow_PVSoftTrk_Y = 0.0, metCoreEMPFlow_PVSoftTrk_SumET = 0.0, metCoreEMPFlow_PVSoftTrk_MET = 0.0;

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

    // event info tree
    eventInfoTree->Branch("mcEventWeight", &mcEventWeights);
    eventInfoTree->Branch("sumOfWeightsForSample", &sumOfWeightsForSample);
    eventInfoTree->Branch("eventWeights", &eventWeights);
    eventInfoTree->Branch("sampleJZSlice", &sampleJZSlice);

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

    // truthTopTree
    truthTopTree->Branch("indexOfHiggs", &indexOfTopValues);
    truthTopTree->Branch("Et", &truthTopEtValues);
    truthTopTree->Branch("Eta", &truthTopEtaValues);
    truthTopTree->Branch("Phi", &truthTopPhiValues);
    truthTopTree->Branch("pT", &truthToppTValues);

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

    // gepBasicClustersTree
    gepBasicClustersTree->Branch("Et", &gepBasicClustersEtValues);
    gepBasicClustersTree->Branch("Eta", &gepBasicClustersEtaValues);
    gepBasicClustersTree->Branch("Phi", &gepBasicClustersPhiValues);

    // gepBasicClustersSKTree
    gepBasicClustersSKTree->Branch("Et", &gepBasicClustersSKEtValues);
    gepBasicClustersSKTree->Branch("Eta", &gepBasicClustersSKEtaValues);
    gepBasicClustersSKTree->Branch("Phi", &gepBasicClustersSKPhiValues);

    // gepCellsTowersTree
    gepCellsTowersTree->Branch("Et", &gepCellsTowersEtValues);
    gepCellsTowersTree->Branch("Eta", &gepCellsTowersEtaValues);
    gepCellsTowersTree->Branch("Phi", &gepCellsTowersPhiValues);

    // gepCellsTowersSKTree
    gepCellsTowersSKTree->Branch("Et", &gepCellsTowersSKEtValues);
    gepCellsTowersSKTree->Branch("Eta", &gepCellsTowersSKEtaValues);
    gepCellsTowersSKTree->Branch("Phi", &gepCellsTowersSKPhiValues);

    // gep wta cone cells towers jets
    gepWTAConeCellsTowersJetsTree->Branch("pT", &gepWTAConeCellsTowersJetspTValues);
    gepWTAConeCellsTowersJetsTree->Branch("Eta", &gepWTAConeCellsTowersJetsEtaValues);
    gepWTAConeCellsTowersJetsTree->Branch("Phi", &gepWTAConeCellsTowersJetsPhiValues);
    gepWTAConeCellsTowersJetsTree->Branch("NConstituents", &gepWTAConeCellsTowersJetsNConstituentsValues);

    // gep wta cone cells towers jets
    gepLeadingWTAConeCellsTowersJetsTree->Branch("pT", &gepLeadingWTAConeCellsTowersJetspTValues);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Eta", &gepLeadingWTAConeCellsTowersJetsEtaValues);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Phi", &gepLeadingWTAConeCellsTowersJetsPhiValues);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("NConstituents", &gepLeadingWTAConeCellsTowersJetsNConstituentsValues);

    // gep wta cone cells towers jets
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("pT", &gepSubleadingWTAConeCellsTowersJetspTValues);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Eta", &gepSubleadingWTAConeCellsTowersJetsEtaValues);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Phi", &gepSubleadingWTAConeCellsTowersJetsPhiValues);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("NConstituents", &gepSubleadingWTAConeCellsTowersJetsNConstituentsValues);

     // gep wta cone cells towers SK jets
    gepWTAConeCellsTowersSKJetsTree->Branch("pT", &gepWTAConeCellsTowersSKJetspTValues);
    gepWTAConeCellsTowersSKJetsTree->Branch("Eta", &gepWTAConeCellsTowersSKJetsEtaValues);
    gepWTAConeCellsTowersSKJetsTree->Branch("Phi", &gepWTAConeCellsTowersSKJetsPhiValues);
    gepWTAConeCellsTowersSKJetsTree->Branch("NConstituents", &gepWTAConeCellsTowersSKJetsNConstituentsValues);

    // gep wta cone cells towers SK jets
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("pT", &gepLeadingWTAConeCellsTowersSKJetspTValues);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Eta", &gepLeadingWTAConeCellsTowersSKJetsEtaValues);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Phi", &gepLeadingWTAConeCellsTowersSKJetsPhiValues);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("NConstituents", &gepLeadingWTAConeCellsTowersSKJetsNConstituentsValues);

    // gep wta cone cells towers SK jets
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("pT", &gepSubleadingWTAConeCellsTowersSKJetspTValues);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Eta", &gepSubleadingWTAConeCellsTowersSKJetsEtaValues);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Phi", &gepSubleadingWTAConeCellsTowersSKJetsPhiValues);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("NConstituents", &gepSubleadingWTAConeCellsTowersSKJetsNConstituentsValues);

    // gep wta cone basic clusters jets
    gepWTAConeBasicClustersJetsTree->Branch("pT", &gepWTAConeGEPBasicClustersJetspTValues);
    gepWTAConeBasicClustersJetsTree->Branch("Eta", &gepWTAConeGEPBasicClustersJetsEtaValues);
    gepWTAConeBasicClustersJetsTree->Branch("Phi", &gepWTAConeGEPBasicClustersJetsPhiValues);
    gepWTAConeBasicClustersJetsTree->Branch("NConstituents", &gepWTAConeGEPBasicClustersJetsNConstituentsValues);

    // gep wta cone basic clusters jets
    gepLeadingWTAConeBasicClustersJetsTree->Branch("pT", &gepLeadingWTAConeGEPBasicClustersJetspTValues);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Eta", &gepLeadingWTAConeGEPBasicClustersJetsEtaValues);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Phi", &gepLeadingWTAConeGEPBasicClustersJetsPhiValues);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("NConstituents", &gepLeadingWTAConeGEPBasicClustersJetsNConstituentsValues);

    // gep wta cone basic clusters jets
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("pT", &gepSubleadingWTAConeGEPBasicClustersJetspTValues);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Eta", &gepSubleadingWTAConeGEPBasicClustersJetsEtaValues);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Phi", &gepSubleadingWTAConeGEPBasicClustersJetsPhiValues);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("NConstituents", &gepSubleadingWTAConeGEPBasicClustersJetsNConstituentsValues);

    // gep wta cone basic clusters SK jets
    gepWTAConeBasicClustersSKJetsTree->Branch("pT", &gepWTAConeGEPBasicClustersSKJetspTValues);
    gepWTAConeBasicClustersSKJetsTree->Branch("Eta", &gepWTAConeGEPBasicClustersSKJetsEtaValues);
    gepWTAConeBasicClustersSKJetsTree->Branch("Phi", &gepWTAConeGEPBasicClustersSKJetsPhiValues);
    gepWTAConeBasicClustersSKJetsTree->Branch("NConstituents", &gepWTAConeGEPBasicClustersSKJetsNConstituentsValues);

    // gep wta cone basic clusters SK jets
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("pT", &gepLeadingWTAConeGEPBasicClustersSKJetspTValues);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Eta", &gepLeadingWTAConeGEPBasicClustersSKJetsEtaValues);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Phi", &gepLeadingWTAConeGEPBasicClustersSKJetsPhiValues);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("NConstituents", &gepLeadingWTAConeGEPBasicClustersSKJetsNConstituentsValues);

    // gep wta cone basic clusters SK jets
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("pT", &gepSubleadingWTAConeGEPBasicClustersSKJetspTValues);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Eta", &gepSubleadingWTAConeGEPBasicClustersSKJetsEtaValues);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Phi", &gepSubleadingWTAConeGEPBasicClustersSKJetsPhiValues);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("NConstituents", &gepSubleadingWTAConeGEPBasicClustersSKJetsNConstituentsValues);

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
    recoAntiKt10UFOCSSKJets->Branch("Mass", &recoAntiKt10LRJMassValues);

    // leadingRecoAntiKt10UFOCSSKJets
    leadingRecoAntiKt10UFOCSSKJets->Branch("Et", &recoAntiKt10LRJLeadingEtValues);
    leadingRecoAntiKt10UFOCSSKJets->Branch("Eta", &recoAntiKt10LRJLeadingEtaValues);
    leadingRecoAntiKt10UFOCSSKJets->Branch("Phi", &recoAntiKt10LRJLeadingPhiValues);
    leadingRecoAntiKt10UFOCSSKJets->Branch("Mass", &recoAntiKt10LRJLeadingMassValues);

    // subleadingRecoAntiKt10UFOCSSKJets
    subleadingRecoAntiKt10UFOCSSKJets->Branch("Et", &recoAntiKt10LRJSubleadingEtValues);
    subleadingRecoAntiKt10UFOCSSKJets->Branch("Eta", &recoAntiKt10LRJSubleadingEtaValues);
    subleadingRecoAntiKt10UFOCSSKJets->Branch("Phi", &recoAntiKt10LRJSubleadingPhiValues);
    subleadingRecoAntiKt10UFOCSSKJets->Branch("Mass", &recoAntiKt10LRJSubleadingMassValues);

    // recoAntiKt10UFOCSSKSoftDropJets
    recoAntiKt10UFOCSSKSoftDropJets->Branch("EtIndex", &recoAntiKt10UFOCSSKSDEtIndexValues);
    recoAntiKt10UFOCSSKSoftDropJets->Branch("Et", &recoAntiKt10UFOCSSKSDEtValues);
    recoAntiKt10UFOCSSKSoftDropJets->Branch("Eta", &recoAntiKt10UFOCSSKSDEtaValues);
    recoAntiKt10UFOCSSKSoftDropJets->Branch("Phi", &recoAntiKt10UFOCSSKSDPhiValues);
    recoAntiKt10UFOCSSKSoftDropJets->Branch("Mass", &recoAntiKt10UFOCSSKSDMassValues);

    // leadingRecoAntiKt10UFOCSSKSoftDropJets
    leadingRecoAntiKt10UFOCSSKSoftDropJets->Branch("Et", &recoAntiKt10UFOCSSKSDLeadingEtValues);
    leadingRecoAntiKt10UFOCSSKSoftDropJets->Branch("Eta", &recoAntiKt10UFOCSSKSDLeadingEtaValues);
    leadingRecoAntiKt10UFOCSSKSoftDropJets->Branch("Phi", &recoAntiKt10UFOCSSKSDLeadingPhiValues);
    leadingRecoAntiKt10UFOCSSKSoftDropJets->Branch("Mass", &recoAntiKt10UFOCSSKSDLeadingMassValues);

    // subleadingRecoAntiKt10UFOCSSKSoftDropJets
    subleadingRecoAntiKt10UFOCSSKSoftDropJets->Branch("Et", &recoAntiKt10UFOCSSKSDSubleadingEtValues);
    subleadingRecoAntiKt10UFOCSSKSoftDropJets->Branch("Eta", &recoAntiKt10UFOCSSKSDSubleadingEtaValues);
    subleadingRecoAntiKt10UFOCSSKSoftDropJets->Branch("Phi", &recoAntiKt10UFOCSSKSDSubleadingPhiValues);
    subleadingRecoAntiKt10UFOCSSKSoftDropJets->Branch("Mass", &recoAntiKt10UFOCSSKSDSubleadingMassValues);

    // antiKt10TruthJetsTree
    antiKt10TruthJetsTree->Branch("EtIndex", &antiKt10TruthEtIndexValues);
    antiKt10TruthJetsTree->Branch("Et", &antiKt10TruthEtValues);
    antiKt10TruthJetsTree->Branch("Eta", &antiKt10TruthEtaValues);
    antiKt10TruthJetsTree->Branch("Phi", &antiKt10TruthPhiValues);
    antiKt10TruthJetsTree->Branch("Mass", &antiKt10TruthMassValues);

    // leadingAntiKt10TruthJetsTree
    leadingAntiKt10TruthJetsTree->Branch("Et", &antiKt10TruthLeadingEtValues);
    leadingAntiKt10TruthJetsTree->Branch("Eta", &antiKt10TruthLeadingEtaValues);
    leadingAntiKt10TruthJetsTree->Branch("Phi", &antiKt10TruthLeadingPhiValues);
    leadingAntiKt10TruthJetsTree->Branch("Mass", &antiKt10TruthLeadingMassValues);

    // subleadingAntiKt10TruthJetsTree
    subleadingAntiKt10TruthJetsTree->Branch("Et", &antiKt10TruthSubleadingEtValues);
    subleadingAntiKt10TruthJetsTree->Branch("Eta", &antiKt10TruthSubleadingEtaValues);
    subleadingAntiKt10TruthJetsTree->Branch("Phi", &antiKt10TruthSubleadingPhiValues);
    subleadingAntiKt10TruthJetsTree->Branch("Mass", &antiKt10TruthSubleadingMassValues);

    // antiKt10TruthSoftDropJetsTree
    antiKt10TruthSoftDropJetsTree->Branch("EtIndex", &antiKt10TruthSDEtIndexValues);
    antiKt10TruthSoftDropJetsTree->Branch("Et", &antiKt10TruthSDEtValues);
    antiKt10TruthSoftDropJetsTree->Branch("Eta", &antiKt10TruthSDEtaValues);
    antiKt10TruthSoftDropJetsTree->Branch("Phi", &antiKt10TruthSDPhiValues);
    antiKt10TruthSoftDropJetsTree->Branch("Mass", &antiKt10TruthSDMassValues);

    // leadingAntiKt10TruthSoftDropJetsTree
    leadingAntiKt10TruthSoftDropJetsTree->Branch("Et", &antiKt10TruthSDLeadingEtValues);
    leadingAntiKt10TruthSoftDropJetsTree->Branch("Eta", &antiKt10TruthSDLeadingEtaValues);
    leadingAntiKt10TruthSoftDropJetsTree->Branch("Phi", &antiKt10TruthSDLeadingPhiValues);
    leadingAntiKt10TruthSoftDropJetsTree->Branch("Mass", &antiKt10TruthSDLeadingMassValues);

    // subleadingAntiKt10TruthSoftDropJetsTree
    subleadingAntiKt10TruthSoftDropJetsTree->Branch("Et", &antiKt10TruthSDSubleadingEtValues);
    subleadingAntiKt10TruthSoftDropJetsTree->Branch("Eta", &antiKt10TruthSDSubleadingEtaValues);
    subleadingAntiKt10TruthSoftDropJetsTree->Branch("Phi", &antiKt10TruthSDSubleadingPhiValues);
    subleadingAntiKt10TruthSoftDropJetsTree->Branch("Mass", &antiKt10TruthSDSubleadingMassValues);

    // truthAntiKt4TruthDressedWZJets
    truthAntiKt4TruthDressedWZJets->Branch("EtIndex", &truthAntiKt4WZSRJEtIndexValues);
    truthAntiKt4TruthDressedWZJets->Branch("Et", &truthAntiKt4WZSRJEtValues);
    truthAntiKt4TruthDressedWZJets->Branch("Eta", &truthAntiKt4WZSRJEtaValues);
    truthAntiKt4TruthDressedWZJets->Branch("Phi", &truthAntiKt4WZSRJPhiValues);
    truthAntiKt4TruthDressedWZJets->Branch("Mass", &truthAntiKt4WZSRJMassValues);

    // leadingTruthAntiKt4TruthDressedWZJets
    leadingTruthAntiKt4TruthDressedWZJets->Branch("Et", &truthAntiKt4WZSRJLeadingEtValues);
    leadingTruthAntiKt4TruthDressedWZJets->Branch("Eta", &truthAntiKt4WZSRJLeadingEtaValues);
    leadingTruthAntiKt4TruthDressedWZJets->Branch("Phi", &truthAntiKt4WZSRJLeadingPhiValues);
    leadingTruthAntiKt4TruthDressedWZJets->Branch("Mass", &truthAntiKt4WZSRJLeadingMassValues);

    // subleadingTruthAntiKt4TruthDressedWZJets
    subleadingTruthAntiKt4TruthDressedWZJets->Branch("Et", &truthAntiKt4WZSRJSubleadingEtValues);
    subleadingTruthAntiKt4TruthDressedWZJets->Branch("Eta", &truthAntiKt4WZSRJSubleadingEtaValues);
    subleadingTruthAntiKt4TruthDressedWZJets->Branch("Phi", &truthAntiKt4WZSRJSubleadingPhiValues);
    subleadingTruthAntiKt4TruthDressedWZJets->Branch("Mass", &truthAntiKt4WZSRJSubleadingMassValues);

    // gFexMETTree branches (one scalar value per event)
    gFexMETTree->Branch("gMHX",   &gMHX);
    gFexMETTree->Branch("gMHY",   &gMHY);
    gFexMETTree->Branch("gMSX",   &gMSX);
    gFexMETTree->Branch("gMSY",   &gMSY);
    gFexMETTree->Branch("gMETX",  &gMETX);
    gFexMETTree->Branch("gMETY",  &gMETY);
    gFexMETTree->Branch("gMET",   &gMET);
    gFexMETTree->Branch("gSumET", &gSumET);

    // gFexMETNoiseCutTree branches
    gFexMETNoiseCutTree->Branch("gMETX",  &gMETX_NC);
    gFexMETNoiseCutTree->Branch("gMETY",  &gMETY_NC);
    gFexMETNoiseCutTree->Branch("gMET",   &gMET_NC);
    gFexMETNoiseCutTree->Branch("gSumET", &gSumET_NC);

    // gFexMETRmsTree branches
    gFexMETRmsTree->Branch("gMETX",  &gMETX_Rms);
    gFexMETRmsTree->Branch("gMETY",  &gMETY_Rms);
    gFexMETRmsTree->Branch("gMET",   &gMET_Rms);
    gFexMETRmsTree->Branch("gSumET", &gSumET_Rms);

    metTruthTree->Branch("metTruthNonIntX",     &metTruthNonIntX);
    metTruthTree->Branch("metTruthNonIntY",     &metTruthNonIntY);
    metTruthTree->Branch("metTruthNonIntSumET", &metTruthNonIntSumET);
    metTruthTree->Branch("metTruthNonInt",      &metTruthNonInt);
    metTruthTree->Branch("metTruthIntX",        &metTruthIntX);
    metTruthTree->Branch("metTruthIntY",        &metTruthIntY);
    metTruthTree->Branch("metTruthIntSumET",    &metTruthIntSumET);
    metTruthTree->Branch("metTruthInt",         &metTruthInt);
    metTruthTree->Branch("metTruthIntOutX",     &metTruthIntOutX);
    metTruthTree->Branch("metTruthIntOutY",     &metTruthIntOutY);
    metTruthTree->Branch("metTruthIntOutSumET", &metTruthIntOutSumET);
    metTruthTree->Branch("metTruthIntOut",      &metTruthIntOut);

    metCoreAntiKt4EMTopoTree->Branch("SoftClus_X",     &metCoreEMTopo_SoftClus_X);
    metCoreAntiKt4EMTopoTree->Branch("SoftClus_Y",     &metCoreEMTopo_SoftClus_Y);
    metCoreAntiKt4EMTopoTree->Branch("SoftClus_SumET", &metCoreEMTopo_SoftClus_SumET);
    metCoreAntiKt4EMTopoTree->Branch("SoftClus_MET",   &metCoreEMTopo_SoftClus_MET);
    metCoreAntiKt4EMTopoTree->Branch("PVSoftTrk_X",     &metCoreEMTopo_PVSoftTrk_X);
    metCoreAntiKt4EMTopoTree->Branch("PVSoftTrk_Y",     &metCoreEMTopo_PVSoftTrk_Y);
    metCoreAntiKt4EMTopoTree->Branch("PVSoftTrk_SumET", &metCoreEMTopo_PVSoftTrk_SumET);
    metCoreAntiKt4EMTopoTree->Branch("PVSoftTrk_MET",   &metCoreEMTopo_PVSoftTrk_MET);
    metCoreAntiKt4EMTopoTree->Branch("SoftClusEM_X",     &metCoreEMTopo_SoftClusEM_X);
    metCoreAntiKt4EMTopoTree->Branch("SoftClusEM_Y",     &metCoreEMTopo_SoftClusEM_Y);
    metCoreAntiKt4EMTopoTree->Branch("SoftClusEM_SumET", &metCoreEMTopo_SoftClusEM_SumET);
    metCoreAntiKt4EMTopoTree->Branch("SoftClusEM_MET",   &metCoreEMTopo_SoftClusEM_MET);

    metCoreAntiKt4EMPFlowTree->Branch("SoftClus_X",     &metCoreEMPFlow_SoftClus_X);
    metCoreAntiKt4EMPFlowTree->Branch("SoftClus_Y",     &metCoreEMPFlow_SoftClus_Y);
    metCoreAntiKt4EMPFlowTree->Branch("SoftClus_SumET", &metCoreEMPFlow_SoftClus_SumET);
    metCoreAntiKt4EMPFlowTree->Branch("SoftClus_MET",   &metCoreEMPFlow_SoftClus_MET);
    metCoreAntiKt4EMPFlowTree->Branch("PVSoftTrk_X",     &metCoreEMPFlow_PVSoftTrk_X);
    metCoreAntiKt4EMPFlowTree->Branch("PVSoftTrk_Y",     &metCoreEMPFlow_PVSoftTrk_Y);
    metCoreAntiKt4EMPFlowTree->Branch("PVSoftTrk_SumET", &metCoreEMPFlow_PVSoftTrk_SumET);
    metCoreAntiKt4EMPFlowTree->Branch("PVSoftTrk_MET",   &metCoreEMPFlow_PVSoftTrk_MET);

    int higgsPassEventCounter = 0;
    namespace fs = std::filesystem;
    // Collect input file names
    std::vector<std::string> fileNames;     // DAOD_JETM42*.root
    std::vector<std::string> gepFileNames;  // GEP*.root

    auto has_prefix = [](const std::string& s, const std::string& p) {
        return s.rfind(p, 0) == 0; // starts with
    };
    auto has_suffix = [](const std::string& s, const std::string& suf) {
        return s.size() >= suf.size() &&
            s.compare(s.size() - suf.size(), suf.size(), suf) == 0; // ends with
    };

    // Loop 1: DAOD_JETM42*.root (from derivation of an AOD)
    for (const auto& entry : fs::directory_iterator(fileDir)) {
        if (!entry.is_regular_file()) continue;
        const std::string fn = entry.path().filename().string();
        if (has_prefix(fn, "DAOD_JETM42") && has_suffix(fn, ".root")) {
            std::cout << "found DAOD file: " << entry.path().string() << "\n";
            fileNames.push_back(entry.path().string());
        }
    }

    // Loop 2: GEP*.root - load TrigGepPerf output ROOT file
    for (const auto& entry : fs::directory_iterator(fileDir)) {
        if (!entry.is_regular_file()) continue;
        const std::string fn = entry.path().filename().string();
        if (has_prefix(fn, "GEP") && has_suffix(fn, ".root")) {
            std::cout << "found TrigGepPerf file name: " << entry.path().string() << "\n";
            gepFileNames.push_back(entry.path().string());
        }
    }

    if (fileNames.empty()) {
        std::cout << "No ROOT files found in directory." << endl;
        return;
    }

    std::cout << "Found " << fileNames.size() << " files." << endl;

    // Main processing loop
    int fileIt = 0;
    for (const auto& fileName : fileNames) {
        fileIt++; 
        //if (fileIt > 9) break; 
        std::cout << "Processing file: " << fileName << endl;
        std::cout << "Processing GEP file: " << gepFileNames[fileIt - 1] << "\n";

        TFile* f = TFile::Open(fileName.c_str());
        if (!f || f->IsZombie()) {
            cerr << "Could not open " << fileName << endl;
            continue;
        }

        // Open the DAOD file
        xAOD::TEvent event(xAOD::TEvent::kClassAccess);
        if (!event.readFrom(f).isSuccess()) {
            cerr << "Cannot read xAOD from file." << endl;
            continue;
        }

        // --- Open GEP ROOT (non-xAOD) ---
        TFile* gf = TFile::Open(gepFileNames[fileIt - 1].c_str());
        if (!gf || gf->IsZombie()) {
            std::cerr << "Could not open GEP file: " << gepFileNames[fileIt - 1] << "\n";
            if (gf) { gf->Close(); delete gf; }
            f->Close(); delete f;
            continue;
        }

        TTree* gt = nullptr;
        gf->GetObject("ntuple", gt);
        if (!gt) {
            std::cerr << "GEP file has no TTree named 'ntuple'\n";
            gf->Close(); delete gf;
            f->Close(); delete f;
            continue;
        }

        // Enable needed branches from TrigGepPerf 
        gt->SetBranchStatus("*", 0);

        gt->SetBranchStatus("GEPBasicClusters_et",  1);
        gt->SetBranchStatus("GEPBasicClusters_eta", 1);
        gt->SetBranchStatus("GEPBasicClusters_phi", 1);

        gt->SetBranchStatus("GEPBasicClustersSK_et",  1);
        gt->SetBranchStatus("GEPBasicClustersSK_eta", 1);
        gt->SetBranchStatus("GEPBasicClustersSK_phi", 1);

        gt->SetBranchStatus("GEPCellsTower_et",  1);
        gt->SetBranchStatus("GEPCellsTower_eta", 1);
        gt->SetBranchStatus("GEPCellsTower_phi", 1);

        gt->SetBranchStatus("GEPCellsTowerSK_et",  1);
        gt->SetBranchStatus("GEPCellsTowerSK_eta", 1);
        gt->SetBranchStatus("GEPCellsTowerSK_phi", 1);

        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_pt",  1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_eta", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_phi", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_nConstituents", 1);

        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_pt",  1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_eta", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_phi", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_nConstituents", 1);

        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_pt",  1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_eta", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_phi", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_nConstituents", 1);

        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_pt",  1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_eta", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_phi", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_nConstituents", 1);

        std::vector<float>* gepBasicClustersEt  = nullptr;
        std::vector<float>* gepBasicClustersEta = nullptr;
        std::vector<float>* gepBasicClustersPhi = nullptr;

        std::vector<float>* gepBasicClustersSKEt  = nullptr;
        std::vector<float>* gepBasicClustersSKEta = nullptr;
        std::vector<float>* gepBasicClustersSKPhi = nullptr;

        std::vector<float>* gepCellsTowersEt  = nullptr;
        std::vector<float>* gepCellsTowersEta = nullptr;
        std::vector<float>* gepCellsTowersPhi = nullptr;

        std::vector<float>* gepCellsTowersSKEt  = nullptr;
        std::vector<float>* gepCellsTowersSKEta = nullptr;
        std::vector<float>* gepCellsTowersSKPhi = nullptr;

        std::vector<float>* gepWTAConeCellsTowersJetsPt  = nullptr;
        std::vector<float>* gepWTAConeCellsTowersJetsEta = nullptr;
        std::vector<float>* gepWTAConeCellsTowersJetsPhi = nullptr;
        std::vector<float>* gepWTAConeCellsTowersJetsNConstituents = nullptr;

        std::vector<float>* gepWTAConeCellsTowersSKJetsPt  = nullptr;
        std::vector<float>* gepWTAConeCellsTowersSKJetsEta = nullptr;
        std::vector<float>* gepWTAConeCellsTowersSKJetsPhi = nullptr;
        std::vector<float>* gepWTAConeCellsTowersSKJetsNConstituents = nullptr;

        std::vector<float>* gepWTAConeBasicClustersJetsPt  = nullptr;
        std::vector<float>* gepWTAConeBasicClustersJetsEta = nullptr;
        std::vector<float>* gepWTAConeBasicClustersJetsPhi = nullptr;
        std::vector<float>* gepWTAConeBasicClustersJetsNConstituents = nullptr;

        std::vector<float>* gepWTAConeBasicClustersSKJetsPt  = nullptr;
        std::vector<float>* gepWTAConeBasicClustersSKJetsEta = nullptr;
        std::vector<float>* gepWTAConeBasicClustersSKJetsPhi = nullptr;
        std::vector<float>* gepWTAConeBasicClustersSKJetsNConstituents = nullptr;

        gt->SetBranchAddress("GEPBasicClusters_et",  &gepBasicClustersEt);
        gt->SetBranchAddress("GEPBasicClusters_eta", &gepBasicClustersEta);
        gt->SetBranchAddress("GEPBasicClusters_phi", &gepBasicClustersPhi);

        gt->SetBranchAddress("GEPBasicClustersSK_et",  &gepBasicClustersSKEt);
        gt->SetBranchAddress("GEPBasicClustersSK_eta", &gepBasicClustersSKEta);
        gt->SetBranchAddress("GEPBasicClustersSK_phi", &gepBasicClustersSKPhi);

        gt->SetBranchAddress("GEPCellsTower_et",  &gepCellsTowersEt);
        gt->SetBranchAddress("GEPCellsTower_eta", &gepCellsTowersEta);
        gt->SetBranchAddress("GEPCellsTower_phi", &gepCellsTowersPhi);

        gt->SetBranchAddress("GEPCellsTowerSK_et",  &gepCellsTowersSKEt);
        gt->SetBranchAddress("GEPCellsTowerSK_eta", &gepCellsTowersSKEta);
        gt->SetBranchAddress("GEPCellsTowerSK_phi", &gepCellsTowersSKPhi);

        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_pt",  &gepWTAConeCellsTowersJetsPt);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_eta", &gepWTAConeCellsTowersJetsEta);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_phi", &gepWTAConeCellsTowersJetsPhi); 
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_nConstituents", &gepWTAConeCellsTowersJetsNConstituents); 

        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_pt",  &gepWTAConeCellsTowersSKJetsPt);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_eta", &gepWTAConeCellsTowersSKJetsEta);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_phi", &gepWTAConeCellsTowersSKJetsPhi); 
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_nConstituents", &gepWTAConeCellsTowersSKJetsNConstituents); 

        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_pt",  &gepWTAConeBasicClustersJetsPt);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_eta", &gepWTAConeBasicClustersJetsEta);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_phi", &gepWTAConeBasicClustersJetsPhi); 
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_nConstituents", &gepWTAConeBasicClustersJetsNConstituents); 

        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_pt",  &gepWTAConeBasicClustersSKJetsPt);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_eta", &gepWTAConeBasicClustersSKJetsEta);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_phi", &gepWTAConeBasicClustersSKJetsPhi); 
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_nConstituents", &gepWTAConeBasicClustersSKJetsNConstituents); 

        std::cout << "  Number of events: " << event.getEntries() << endl;

        unsigned int passedEventsCounter = 0;
        unsigned int skippedEventsEmptyTruth = 0;
        unsigned int skippedEventsHSTP = 0;
        for (Long64_t iEvt = 0; iEvt < event.getEntries(); ++iEvt) { // NOTE assume that # of events for GEP and DAOD files is the same, else will get a seg fault.
            //std::cout << "iEvt: " << iEvt << "\n";
            // Retrieve truth and in time pileup jets first to apply hard-scatter-softer-than-PU (HSTP) filter (described here: https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtMissMCSamples#Dijet_normalization_procedure_HS)
            // Also require that truth jet collection is not empty
            event.getEntry(iEvt);               // DAOD event iEvt
            gt->GetEntry(iEvt);                 // GEP event iEvt

            const xAOD::JetContainer* AntiKt4TruthDressedWZJets = nullptr;
            if (!event.retrieve(AntiKt4TruthDressedWZJets, "AntiKt4TruthDressedWZJets").isSuccess()) {
                std::cout << "Failed to retrieve reco Antik4 Truth Dressed WZ jets" << endl;
                skippedEventsEmptyTruth++;
                continue;
            }
            if(!signalBool){
                if ( AntiKt4TruthDressedWZJets->size() == 0 ){
                    skippedEventsEmptyTruth++;
                    continue;
                } 
            }

            const xAOD::JetContainer* InTimeAntiKt4TruthJets = nullptr;
            if (!event.retrieve(InTimeAntiKt4TruthJets, "InTimeAntiKt4TruthJets").isSuccess()) {
                std::cout << "Failed to retrieve Truth jets" << endl;
                continue;
            }

            // IMPORTANT NOTE: must first sort both jet collections by pT, most (but NOT ALL) are already sorted with element 0 representing leading jet, 
            // but failing to sort for ALL events will result in rate spikes (observed for ~200 GeV pileup jets for ~25 GeV truth jets!)
            // --- Sort AntiKt4TruthDressedWZJets by pT (descending)
            ConstDataVector<xAOD::JetContainer> sortedWZJets(SG::VIEW_ELEMENTS);
            sortedWZJets.reserve(AntiKt4TruthDressedWZJets->size());

            for (const xAOD::Jet* jet : *AntiKt4TruthDressedWZJets) {
                sortedWZJets.push_back(jet);
            }

            std::sort(sortedWZJets.begin(), sortedWZJets.end(),
                    [](const xAOD::Jet* a, const xAOD::Jet* b) {
                        return a->pt() > b->pt();
                    });

            // --- Sort InTimeAntiKt4TruthJets by pT (descending)
            ConstDataVector<xAOD::JetContainer> sortedTruthJets(SG::VIEW_ELEMENTS);
            sortedTruthJets.reserve(InTimeAntiKt4TruthJets->size());

            for (const xAOD::Jet* jet : *InTimeAntiKt4TruthJets) {
                sortedTruthJets.push_back(jet);
            }

            std::sort(sortedTruthJets.begin(), sortedTruthJets.end(),
                    [](const xAOD::Jet* a, const xAOD::Jet* b) {
                        return a->pt() > b->pt();
                    });

            // Filter out events where hard scatter is softer than PU
            if(!signalBool){
                if(sortedWZJets.at(0)->pt() <= sortedTruthJets.at(0)->pt()){
                    skippedEventsHSTP++;
                    continue; 
                } 
            }
            

            passedEventsCounter++;
            
            

            if(iEvt % 100 == 0) std::cout << "iEvt: " << iEvt << "\n";
            

            // -- retrieve collections from DAOD ---
            const xAOD::EventInfo_v1* EventInfo = nullptr;
            if (!event.retrieve(EventInfo, "EventInfo").isSuccess()) {
                cerr << "Cannot access EventInfo" << endl;
                continue;
            }
            
            const xAOD::TruthParticleContainer* TruthBosonsWithDecayParticles = nullptr;
            if (!event.retrieve(TruthBosonsWithDecayParticles, "TruthBosonsWithDecayParticles").isSuccess()) {
                cerr << "Cannot access TruthBosonsWithDecayParticles" << endl;
                continue;
            }

            const xAOD::TruthParticleContainer* TruthTop = nullptr;
            if (!event.retrieve(TruthTop, "TruthTop").isSuccess()) {
                cerr << "Cannot access TruthTop" << endl;
                continue;
            }
            
            const xAOD::TruthParticleContainer* TruthHFWithDecayParticles = nullptr;
            if (!event.retrieve(TruthHFWithDecayParticles, "TruthHFWithDecayParticles").isSuccess()) {
                cerr << "Cannot access TruthHFWithDecayParticles" << endl;
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

            // gFEX MET containers
            const DataVector<xAOD::gFexGlobalRoI_v1>* L1_gMHTComponentsJwoj = nullptr;
            if (!event.retrieve(L1_gMHTComponentsJwoj, "L1_gMHTComponentsJwoj").isSuccess()) {
                std::cerr << "Failed to retrieve gFEX MET hard term" << std::endl;
            }
            const DataVector<xAOD::gFexGlobalRoI_v1>* L1_gMSTComponentsJwoj = nullptr;
            if (!event.retrieve(L1_gMSTComponentsJwoj, "L1_gMSTComponentsJwoj").isSuccess()) {
                std::cerr << "Failed to retrieve gFEX MET soft term" << std::endl;
            }
            const DataVector<xAOD::gFexGlobalRoI_v1>* L1_gMETComponentsJwoj = nullptr;
            if (!event.retrieve(L1_gMETComponentsJwoj, "L1_gMETComponentsJwoj").isSuccess()) {
                std::cerr << "Failed to retrieve gFEX total MET" << std::endl;
            }
            const DataVector<xAOD::gFexGlobalRoI_v1>* L1_gScalarEJwoj = nullptr;
            if (!event.retrieve(L1_gScalarEJwoj, "L1_gScalarEJwoj").isSuccess()) {
                std::cerr << "Failed to retrieve gFEX scalar MET/SumET" << std::endl;
            }

            const DataVector<xAOD::gFexGlobalRoI_v1>* L1_gMETComponentsNoiseCut = nullptr;
            if (!event.retrieve(L1_gMETComponentsNoiseCut, "L1_gMETComponentsNoiseCut").isSuccess()) {
                std::cerr << "Failed to retrieve gFEX scalar MET/SumET" << std::endl;
            }

            const DataVector<xAOD::gFexGlobalRoI_v1>* L1_gScalarENoiseCut = nullptr;
            if (!event.retrieve(L1_gScalarENoiseCut, "L1_gScalarENoiseCut").isSuccess()) {
                std::cerr << "Failed to retrieve gFEX scalar MET/SumET" << std::endl;
            }

            const DataVector<xAOD::gFexGlobalRoI_v1>* L1_gMETComponentsRms = nullptr;
            if (!event.retrieve(L1_gMETComponentsRms, "L1_gMETComponentsRms").isSuccess()) {
                std::cerr << "Failed to retrieve gFEX scalar MET/SumET" << std::endl;
            }

            const DataVector<xAOD::gFexGlobalRoI_v1>* L1_gScalarERms = nullptr;
            if (!event.retrieve(L1_gScalarERms, "L1_gScalarERms").isSuccess()) {
                std::cerr << "Failed to retrieve gFEX scalar MET/SumET" << std::endl;
            }

            // Reference MET containers
            const xAOD::MissingETContainer_v1* MET_Truth = nullptr;
            if (!event.retrieve(MET_Truth, "MET_Truth").isSuccess()) {
                std::cerr << "Failed to retrieve MET_Truth" << std::endl;
            }

            const xAOD::MissingETContainer_v1* MET_Core_AntiKt4EMTopo = nullptr;
            if (!event.retrieve(MET_Core_AntiKt4EMTopo, "MET_Core_AntiKt4EMTopo").isSuccess()) {
                std::cerr << "Failed to retrieve MET_Core_AntiKt4EMTopo" << std::endl;
            }

            const xAOD::MissingETContainer_v1* MET_Core_AntiKt4EMPFlow = nullptr;
            if (!event.retrieve(MET_Core_AntiKt4EMPFlow, "MET_Core_AntiKt4EMPFlow").isSuccess()) {
                std::cerr << "Failed to retrieve MET_Core_AntiKt4EMPFlow" << std::endl;
            }

            const xAOD::JetContainer* AntiKt10UFOCSSKJets = nullptr;
            if (!event.retrieve(AntiKt10UFOCSSKJets, "AntiKt10UFOCSSKJets").isSuccess()) {
                cerr << "Failed to retrieve reco Antik10 UFOCSSK jets" << endl;
                continue;
            }

            const xAOD::JetContainer* AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets = nullptr;
            if (!event.retrieve(AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets, "AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets").isSuccess()) {
                cerr << "Failed to retrieve reco AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets" << endl;
                continue;
            }

            const xAOD::JetContainer* AntiKt10TruthJets = nullptr;
            if (!event.retrieve(AntiKt10TruthJets, "AntiKt10TruthJets").isSuccess()) {
                cerr << "Failed to retrieve reco AntiKt10TruthJets" << endl;
                continue;
            }

            const xAOD::JetContainer* AntiKt10TruthSoftDropBeta100Zcut10Jets = nullptr;
            if (!event.retrieve(AntiKt10TruthSoftDropBeta100Zcut10Jets, "AntiKt10TruthSoftDropBeta100Zcut10Jets").isSuccess()) {
                cerr << "Failed to retrieve reco AntiKt10TruthSoftDropBeta100Zcut10Jets" << endl;
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
            recoAntiKt10LRJMassValues.clear();
            recoAntiKt10LRJLeadingEtValues.clear();
            recoAntiKt10LRJLeadingEtaValues.clear();
            recoAntiKt10LRJLeadingPhiValues.clear();
            recoAntiKt10LRJLeadingMassValues.clear();
            recoAntiKt10LRJSubleadingEtValues.clear();
            recoAntiKt10LRJSubleadingEtaValues.clear();
            recoAntiKt10LRJSubleadingPhiValues.clear();
            recoAntiKt10LRJSubleadingMassValues.clear();
            recoAntiKt10UFOCSSKSDEtIndexValues.clear();
            recoAntiKt10UFOCSSKSDEtValues.clear();
            recoAntiKt10UFOCSSKSDEtaValues.clear();
            recoAntiKt10UFOCSSKSDPhiValues.clear();
            recoAntiKt10UFOCSSKSDMassValues.clear();
            recoAntiKt10UFOCSSKSDLeadingEtValues.clear();
            recoAntiKt10UFOCSSKSDLeadingEtaValues.clear();
            recoAntiKt10UFOCSSKSDLeadingPhiValues.clear();
            recoAntiKt10UFOCSSKSDLeadingMassValues.clear();
            recoAntiKt10UFOCSSKSDSubleadingEtValues.clear();
            recoAntiKt10UFOCSSKSDSubleadingEtaValues.clear();
            recoAntiKt10UFOCSSKSDSubleadingPhiValues.clear();
            recoAntiKt10UFOCSSKSDSubleadingMassValues.clear();
            antiKt10TruthEtIndexValues.clear();
            antiKt10TruthEtValues.clear();
            antiKt10TruthEtaValues.clear();
            antiKt10TruthPhiValues.clear();
            antiKt10TruthMassValues.clear();
            antiKt10TruthLeadingEtValues.clear();
            antiKt10TruthLeadingEtaValues.clear();
            antiKt10TruthLeadingPhiValues.clear();
            antiKt10TruthLeadingMassValues.clear();
            antiKt10TruthSubleadingEtValues.clear();
            antiKt10TruthSubleadingEtaValues.clear();
            antiKt10TruthSubleadingPhiValues.clear();
            antiKt10TruthSubleadingMassValues.clear();
            antiKt10TruthSDEtIndexValues.clear();
            antiKt10TruthSDEtValues.clear();
            antiKt10TruthSDEtaValues.clear();
            antiKt10TruthSDPhiValues.clear();
            antiKt10TruthSDMassValues.clear();
            antiKt10TruthSDLeadingEtValues.clear();
            antiKt10TruthSDLeadingEtaValues.clear();
            antiKt10TruthSDLeadingPhiValues.clear();
            antiKt10TruthSDLeadingMassValues.clear();
            antiKt10TruthSDSubleadingEtValues.clear();
            antiKt10TruthSDSubleadingEtaValues.clear();
            antiKt10TruthSDSubleadingPhiValues.clear();
            antiKt10TruthSDSubleadingMassValues.clear();
            truthAntiKt4WZSRJEtIndexValues.clear();
            truthAntiKt4WZSRJEtValues.clear();
            truthAntiKt4WZSRJEtaValues.clear();
            truthAntiKt4WZSRJPhiValues.clear();
            truthAntiKt4WZSRJMassValues.clear();
            truthAntiKt4WZSRJLeadingEtValues.clear();
            truthAntiKt4WZSRJLeadingEtaValues.clear();
            truthAntiKt4WZSRJLeadingPhiValues.clear();
            truthAntiKt4WZSRJLeadingMassValues.clear();
            truthAntiKt4WZSRJSubleadingEtValues.clear();
            truthAntiKt4WZSRJSubleadingEtaValues.clear();
            truthAntiKt4WZSRJSubleadingPhiValues.clear();
            truthAntiKt4WZSRJSubleadingMassValues.clear();
            
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
            gepBasicClustersEtValues.clear();
            gepBasicClustersEtaValues.clear();
            gepBasicClustersPhiValues.clear();
            gepCellsTowersEtValues.clear();
            gepCellsTowersEtaValues.clear();
            gepCellsTowersPhiValues.clear();
            gepBasicClustersSKEtValues.clear();
            gepBasicClustersSKEtaValues.clear();
            gepBasicClustersSKPhiValues.clear();
            gepCellsTowersSKEtValues.clear();
            gepCellsTowersSKEtaValues.clear();
            gepCellsTowersSKPhiValues.clear();
            gepWTAConeCellsTowersJetspTValues.clear();
            gepWTAConeCellsTowersJetsEtaValues.clear();
            gepWTAConeCellsTowersJetsPhiValues.clear();
            gepWTAConeCellsTowersJetsNConstituentsValues.clear();
            gepWTAConeGEPBasicClustersJetspTValues.clear();
            gepWTAConeGEPBasicClustersJetsEtaValues.clear();
            gepWTAConeGEPBasicClustersJetsPhiValues.clear();
            gepWTAConeGEPBasicClustersJetsNConstituentsValues.clear();
            gepLeadingWTAConeCellsTowersJetspTValues.clear();
            gepLeadingWTAConeCellsTowersJetsEtaValues.clear();
            gepLeadingWTAConeCellsTowersJetsPhiValues.clear();
            gepLeadingWTAConeCellsTowersJetsNConstituentsValues.clear();
            gepLeadingWTAConeGEPBasicClustersJetspTValues.clear();
            gepLeadingWTAConeGEPBasicClustersJetsEtaValues.clear();
            gepLeadingWTAConeGEPBasicClustersJetsPhiValues.clear();
            gepLeadingWTAConeGEPBasicClustersJetsNConstituentsValues.clear();
            gepSubleadingWTAConeCellsTowersJetspTValues.clear();
            gepSubleadingWTAConeCellsTowersJetsEtaValues.clear();
            gepSubleadingWTAConeCellsTowersJetsPhiValues.clear();
            gepSubleadingWTAConeCellsTowersJetsNConstituentsValues.clear();
            gepSubleadingWTAConeGEPBasicClustersJetspTValues.clear();
            gepSubleadingWTAConeGEPBasicClustersJetsEtaValues.clear();
            gepSubleadingWTAConeGEPBasicClustersJetsPhiValues.clear();
            gepSubleadingWTAConeGEPBasicClustersJetsNConstituentsValues.clear();
            gepWTAConeCellsTowersSKJetspTValues.clear();
            gepWTAConeCellsTowersSKJetsEtaValues.clear();
            gepWTAConeCellsTowersSKJetsPhiValues.clear();
            gepWTAConeCellsTowersSKJetsNConstituentsValues.clear();
            gepWTAConeGEPBasicClustersSKJetspTValues.clear();
            gepWTAConeGEPBasicClustersSKJetsEtaValues.clear();
            gepWTAConeGEPBasicClustersSKJetsPhiValues.clear();
            gepWTAConeGEPBasicClustersSKJetsNConstituentsValues.clear();
            gepLeadingWTAConeCellsTowersSKJetspTValues.clear();
            gepLeadingWTAConeCellsTowersSKJetsEtaValues.clear();
            gepLeadingWTAConeCellsTowersSKJetsPhiValues.clear();
            gepLeadingWTAConeCellsTowersSKJetsNConstituentsValues.clear();
            gepLeadingWTAConeGEPBasicClustersSKJetspTValues.clear();
            gepLeadingWTAConeGEPBasicClustersSKJetsEtaValues.clear();
            gepLeadingWTAConeGEPBasicClustersSKJetsPhiValues.clear();
            gepLeadingWTAConeGEPBasicClustersSKJetsNConstituentsValues.clear();
            gepSubleadingWTAConeCellsTowersSKJetspTValues.clear();
            gepSubleadingWTAConeCellsTowersSKJetsEtaValues.clear();
            gepSubleadingWTAConeCellsTowersSKJetsPhiValues.clear();
            gepSubleadingWTAConeCellsTowersSKJetsNConstituentsValues.clear();
            gepSubleadingWTAConeGEPBasicClustersSKJetspTValues.clear();
            gepSubleadingWTAConeGEPBasicClustersSKJetsEtaValues.clear();
            gepSubleadingWTAConeGEPBasicClustersSKJetsPhiValues.clear();
            gepSubleadingWTAConeGEPBasicClustersSKJetsNConstituentsValues.clear();
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
            mcEventWeights.clear();
            eventWeights.clear();
            sumOfWeightsForSample = 0;
            sampleJZSlice = 0;
            
            // Get sum of weights for sample, individual Monte Carlo event weight, and compute event weight used later for reweighting histograms
            sumOfWeightsForSample = sumOfEventWeightsByJZSlice[jzSlice];
            
            if (signalBool) sampleJZSlice = -1; 
            float mcEventWeight = EventInfo->mcEventWeight();
            //std::cout << "iEvt: " << iEvt << " and event weight: " << eventWeight << "\n";
            mcEventWeights.push_back(mcEventWeight);

            // Compute weight for histograms 
            double eventWeight = mcEventWeight * crossSectionsByJZSlice[jzSlice] * filterEffienciesByJZSlice[jzSlice] * reweightLuminosity / (sumOfEventWeightsByJZSlice[jzSlice]);
            eventWeights.push_back(eventWeight);
            
            //std::cout << "----------------- processing event --------------------------" << "\n";
            //std::cout << "iEvt: " << iEvt << "\n";
            for(const auto& el : *TruthTop){
                //std::cout << "----------------- looping within TruthTop ------------------- " << "\n";
                //std::cout << "particle pdg id : " << el->pdgId() << " , and status: " << el->status() << "\n";
                //std::cout << "particle pT: " << el->pt() / 1000.0 << "\n";
                //std::cout << "particle number of children: " << el->nChildren() << "\n";
                for (unsigned int i = 0; i < el->nChildren(); ++i) {
                    //std::cout << " particle child iterator i: " << i << "\n";
                    //std::cout << "pdg of particle child: " << el->child(i)->pdgId() << "\n";
                    //if(abs(el->child(i)->pdgId()) == 6) std::cout << "b pT: " << el->child(i)->pt() << "\n";
                    //if(abs(el->child(i)->pdgId()) == 24) std::cout << "W pT: " << el->child(i)->pt() << "\n";
                    
                }
                //if(el->pdgId() == 6) std::cout << "top quark within TruthHFWithDecayParticles " << "\n";
                
            }
            /*for(const auto& el : *TruthHFWithDecayParticles){
                std::cout << "----------------- looping within TruthHFWithDecayParticles ------------------- " << "\n";
                std::cout << "particle pdg id : " << el->pdgId() << " , and status: " << el->status() << "\n";
                std::cout << "particle pt: " << el->pt() / 1000.0 << "\n";
                std::cout << "particle number of children: " << el->nChildren() << "\n";
                for (unsigned int i = 0; i < el->nChildren(); ++i) {
                    std::cout << " particle child iterator i: " << i << "\n";
                    std::cout << "pdg of particle child: " << el->child(i)->pdgId() << "\n";
                }
                if(el->pdgId() == 6) std::cout << "top quark within TruthHFWithDecayParticles " << "\n";
                
            }*/

            // Initialize per-event counters
            unsigned int higgs_counter = -1;
            // --- Loop over TruthParticles (for Higgs and B's) ---
            std::vector<std::vector<float > > allb_list;
            for (const auto& el : *TruthBosonsWithDecayParticles) {
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

            unsigned int gepbasicclusters_it = 0;
            for (unsigned int i = 0; i < gepBasicClustersEt->size(); ++i) {
                float gepBasicClusterEt  = (*gepBasicClustersEt)[i] / 1000.0;
                float gepBasicClusterEta = (*gepBasicClustersEta)[i];
                float gepBasicClusterPhi = (*gepBasicClustersPhi)[i];

                gepBasicClustersEtValues.push_back(gepBasicClusterEt);
                gepBasicClustersEtaValues.push_back(gepBasicClusterEta);
                gepBasicClustersPhiValues.push_back(gepBasicClusterPhi);

                if (gepBasicClusterEt <= 0) continue; // don't store to digitized memories

                // Digitize each variable
                int phi_bin = digitize(gepBasicClusterPhi, phi_bit_length_, phi_min_, phi_max_);
                int eta_bin = digitize(gepBasicClusterEta, eta_bit_length_, eta_min_, eta_max_, etaAltRange);
                int et_bin  = digitize(gepBasicClusterEt, et_bit_length_,
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
                                    
                if (gepbasicclusters_it == 0) {
                    f_gepbasicclusters << "Event : " << std::dec << iEvt << "\n";
                }
                // 4. Write to output file
                f_gepbasicclusters << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepbasicclusters_it
                            << " " << binary_word
                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                gepbasicclusters_it++;

            }

            // SK basic clusters
            unsigned int gepbasicclusterssk_it = 0;
            for (unsigned int i = 0; i < gepBasicClustersSKEt->size(); ++i) {
                float gepBasicClusterSKEt  = (*gepBasicClustersSKEt)[i] / 1000.0;
                float gepBasicClusterSKEta = (*gepBasicClustersSKEta)[i];
                float gepBasicClusterSKPhi = (*gepBasicClustersSKPhi)[i];

                gepBasicClustersSKEtValues.push_back(gepBasicClusterSKEt);
                gepBasicClustersSKEtaValues.push_back(gepBasicClusterSKEta);
                gepBasicClustersSKPhiValues.push_back(gepBasicClusterSKPhi);

                if (gepBasicClusterSKEt <= 0) continue; // don't store to digitized memories - and remove 0 Et clusters, which were already "removed" by the soft killer algorith

                // Digitize each variable
                int phi_bin = digitize(gepBasicClusterSKPhi, phi_bit_length_, phi_min_, phi_max_);
                int eta_bin = digitize(gepBasicClusterSKEta, eta_bit_length_, eta_min_, eta_max_, etaAltRange);
                int et_bin  = digitize(gepBasicClusterSKEt, et_bit_length_,
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
                                    
                if (gepbasicclusterssk_it == 0) {
                    f_gepbasicclusterssk << "Event : " << std::dec << iEvt << "\n";
                }
                // 4. Write to output file
                f_gepbasicclusterssk << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepbasicclusterssk_it
                            << " " << binary_word
                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                gepbasicclusterssk_it++;
            }

            // ---------- gepCellsTowers ----------
            {
                // Build index order sorted by Et (descending)
                std::vector<unsigned int> order(gepCellsTowersEt->size());
                std::iota(order.begin(), order.end(), 0u);
                std::sort(order.begin(), order.end(),
                        [&](unsigned int a, unsigned int b){
                            return (*gepCellsTowersEt)[a] > (*gepCellsTowersEt)[b];
                        });

                unsigned int gepcellstowers_it = 0;
                for (unsigned int i : order) {
                    float gepCellsTowersEtValue  = (*gepCellsTowersEt)[i] / 1000.0f;
                    float gepCellsTowersEtaValue = (*gepCellsTowersEta)[i];
                    float gepCellsTowersPhiValue = (*gepCellsTowersPhi)[i];

                    // push in Et-sorted order
                    gepCellsTowersEtValues.push_back(gepCellsTowersEtValue);
                    gepCellsTowersEtaValues.push_back(gepCellsTowersEtaValue);
                    gepCellsTowersPhiValues.push_back(gepCellsTowersPhiValue);

                    if (gepCellsTowersEtValue <= 0) continue; // keep your existing rule

                    // Digitize
                    int phi_bin = digitize(gepCellsTowersPhiValue, phi_bit_length_, phi_min_, phi_max_);
                    int eta_bin = digitize(gepCellsTowersEtaValue, eta_bit_length_, eta_min_, eta_max_, etaAltRange);
                    int et_bin  = digitize(gepCellsTowersEtValue,  et_bit_length_,
                                        static_cast<double>(et_min_), static_cast<double>(et_max_));

                    int et_bin_sort_purposes = digitize(gepCellsTowersEtValue, 12, static_cast<double>(et_min_), static_cast<double>(et_max_));
                    int phi_bin_sort_purposes = digitize(gepCellsTowersPhiValue, 6, phi_min_, phi_max_);
                    int eta_bin_sort_purposes = digitize(gepCellsTowersEtValue, 7, static_cast<double>(et_min_), static_cast<double>(et_max_), etaAltRange);
                    int error_bin_sort_purposes = digitize(0, 7, 0, 1 << 7);

                    // Build binary string
                    {   
                        std::stringstream binary_ss;
                        binary_ss << std::bitset<7>(error_bin_sort_purposes)  << "|"
                                << std::bitset<7>(eta_bin_sort_purposes)  << "|"
                                << std::bitset<6>(phi_bin_sort_purposes) << "|"
                                << std::bitset<12>(et_bin_sort_purposes);
                        std::string binary_word = binary_ss.str();

                        // Pack into 27-bit word
                        uint32_t packed_word = (error_bin_sort_purposes  << (12 + 6 + 7)) |
                                            (eta_bin_sort_purposes  << (12 + 6)) |
                                            (phi_bin_sort_purposes <<  12) |
                                            (et_bin_sort_purposes);

                        if (gepcellstowers_it == 0) {
                            f_gepcellstowers_sort << "Event : " << std::dec << iEvt << "\n";
                        }
                        f_gepcellstowers_sort << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepcellstowers_it
                                        << " "  << binary_word
                                        << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    } // Fixing scope since using the same variable names here
                    

                    // Build binary string
                    std::stringstream binary_ss;
                    binary_ss << std::bitset<et_bit_length_>(et_bin)  << "|"
                            << std::bitset<eta_bit_length_>(eta_bin) << "|"
                            << std::bitset<phi_bit_length_>(phi_bin);
                    std::string binary_word = binary_ss.str();

                    // Pack into 27-bit word
                    uint32_t packed_word = (et_bin  << (eta_bit_length_ + phi_bit_length_)) |
                                        (eta_bin <<  phi_bit_length_) |
                                        (phi_bin);

                    if (gepcellstowers_it == 0) {
                        f_gepcellstowers << "Event : " << std::dec << iEvt << "\n";
                    }
                    f_gepcellstowers << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepcellstowers_it
                                    << " "  << binary_word
                                    << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    ++gepcellstowers_it;

                }
            }

            // ---------- gepCellsTowers SK PU Suppression applied ----------
            {
                // Build index order sorted by Et (descending)
                std::vector<unsigned int> order(gepCellsTowersSKEt->size());
                std::iota(order.begin(), order.end(), 0u);
                std::sort(order.begin(), order.end(),
                        [&](unsigned int a, unsigned int b){
                            return (*gepCellsTowersSKEt)[a] > (*gepCellsTowersSKEt)[b];
                        });

                unsigned int gepcellstowerssk_it = 0;
                for (unsigned int i : order) {
                    float gepCellsTowersSKEtValue  = (*gepCellsTowersSKEt)[i] / 1000.0f;
                    float gepCellsTowersSKEtaValue = (*gepCellsTowersSKEta)[i];
                    float gepCellsTowersSKPhiValue = (*gepCellsTowersSKPhi)[i];

                    // push in Et-sorted order
                    gepCellsTowersSKEtValues.push_back(gepCellsTowersSKEtValue);
                    gepCellsTowersSKEtaValues.push_back(gepCellsTowersSKEtaValue);
                    gepCellsTowersSKPhiValues.push_back(gepCellsTowersSKPhiValue);

                    if (gepCellsTowersSKEtValue <= 0) continue; // keep your existing rule

                    // Digitize
                    int phi_bin = digitize(gepCellsTowersSKPhiValue, phi_bit_length_, phi_min_, phi_max_);
                    int eta_bin = digitize(gepCellsTowersSKEtaValue, eta_bit_length_, eta_min_, eta_max_, etaAltRange);
                    int et_bin  = digitize(gepCellsTowersSKEtValue,  et_bit_length_,
                                        static_cast<double>(et_min_), static_cast<double>(et_max_));
                    

                    // Build binary string
                    std::stringstream binary_ss;
                    binary_ss << std::bitset<et_bit_length_>(et_bin)  << "|"
                            << std::bitset<eta_bit_length_>(eta_bin) << "|"
                            << std::bitset<phi_bit_length_>(phi_bin);
                    std::string binary_word = binary_ss.str();

                    // Pack into 27-bit word
                    uint32_t packed_word = (et_bin  << (eta_bit_length_ + phi_bit_length_)) |
                                        (eta_bin <<  phi_bit_length_) |
                                        (phi_bin);

                    if (gepcellstowerssk_it == 0) {
                        f_gepcellstowerssk << "Event : " << std::dec << iEvt << "\n";
                    }
                    f_gepcellstowerssk << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepcellstowerssk_it
                                    << " "  << binary_word
                                    << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    ++gepcellstowerssk_it;
                }
            }


            // Loop over clusters and fill Et, Eta, Phi
            unsigned int topotower_it = 0;
            for (const auto* cluster : *CaloCalAllTopoTowers) {
                if (!cluster) continue;

                double et = cluster->e() / cosh(cluster->eta()) / 1000.0; // Et in GeV
                caloTopoTowerEtValues.push_back(et);
                caloTopoTowerEtaValues.push_back(cluster->eta());
                caloTopoTowerPhiValues.push_back(cluster->phi());

                if (et <= 0) continue; // don't store to digitized memories

                // Digitize each variable
                int phi_bin = digitize(cluster->phi(), phi_bit_length_, phi_min_, phi_max_);
                int eta_bin = digitize(cluster->eta(), eta_bit_length_, eta_min_, eta_max_, etaAltRange);
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
                                    
                if (topotower_it == 0) {
                    f_topotower << "Event : " << std::dec << iEvt << "\n";
                }
                // 4. Write to output file
                f_topotower << "0x" << std::hex << std::setw(2) << std::setfill('0') << topotower_it
                            << " " << binary_word
                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                topotower_it++;
            }

            unsigned int topocluster422_it = 0;
            // Loop over the clusters and store Et, Eta, Phi
            for (const auto* cluster : *CaloTopoClusters422) {
                if (!cluster) continue;

                double et = cluster->e() / cosh(cluster->eta()) / 1000.0; // Et in GeV
                topo422EtValues.push_back(et);
                topo422EtaValues.push_back(cluster->eta());
                topo422PhiValues.push_back(cluster->phi());

                if (et <= 0) continue; // FIXME for now don't store negative Et to digitized memories
                // Digitize each variable
                int phi_bin = digitize(cluster->phi(), phi_bit_length_, phi_min_, phi_max_);
                int eta_bin = digitize(cluster->eta(), eta_bit_length_, eta_min_, eta_max_, etaAltRange);
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
                /*if (topocluster422_it == 0) {
                    f_topo << "Event : " << std::dec << iEvt << "\n";
                }
                f_topo << "0x" << std::hex << std::setw(2) << std::setfill('0') << topocluster422_it
                        << " " << binary_word
                        << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                topocluster422_it++; */
            }

            // GEP cone jets from GEPCellsTowers
            {
                std::vector<unsigned int> indices(gepWTAConeCellsTowersJetsPt->size());
                std::iota(indices.begin(), indices.end(), 0);

                std::sort(indices.begin(), indices.end(),
                        [&](unsigned int a, unsigned int b) {
                            return (*gepWTAConeCellsTowersJetsPt)[a] > (*gepWTAConeCellsTowersJetsPt)[b];
                        });
                unsigned int gepwtaconecellstowers_it = 0;
                for (unsigned int i = 0; i < indices.size(); i++) {
                    unsigned int idx = indices[i];
                    float WTAConeCellsTowersJetspT  = (*gepWTAConeCellsTowersJetsPt)[idx] / 1000.0;
                    float WTAConeCellsTowersJetsEta = (*gepWTAConeCellsTowersJetsEta)[idx];
                    float WTAConeCellsTowersJetsPhi = (*gepWTAConeCellsTowersJetsPhi)[idx];
                    unsigned int WTAConeCellsTowersJetsNConstituents = (*gepWTAConeCellsTowersJetsNConstituents)[idx];
                    if(i == 0){
                        gepLeadingWTAConeCellsTowersJetspTValues.push_back(WTAConeCellsTowersJetspT);
                        gepLeadingWTAConeCellsTowersJetsEtaValues.push_back(WTAConeCellsTowersJetsEta);
                        gepLeadingWTAConeCellsTowersJetsPhiValues.push_back(WTAConeCellsTowersJetsPhi);
                        gepLeadingWTAConeCellsTowersJetsNConstituentsValues.push_back(WTAConeCellsTowersJetsNConstituents);
                    }
                    if (i == 1){
                        gepSubleadingWTAConeCellsTowersJetspTValues.push_back(WTAConeCellsTowersJetspT);
                        gepSubleadingWTAConeCellsTowersJetsEtaValues.push_back(WTAConeCellsTowersJetsEta);
                        gepSubleadingWTAConeCellsTowersJetsPhiValues.push_back(WTAConeCellsTowersJetsPhi);
                        gepSubleadingWTAConeCellsTowersJetsNConstituentsValues.push_back(WTAConeCellsTowersJetsNConstituents);
                    }
                
                    gepWTAConeCellsTowersJetspTValues.push_back(WTAConeCellsTowersJetspT);
                    gepWTAConeCellsTowersJetsEtaValues.push_back(WTAConeCellsTowersJetsEta);
                    gepWTAConeCellsTowersJetsPhiValues.push_back(WTAConeCellsTowersJetsPhi);
                    gepWTAConeCellsTowersJetsNConstituentsValues.push_back(WTAConeCellsTowersJetsNConstituents);

                    if (WTAConeCellsTowersJetspT <= 0) continue;

                    int phi_bin = digitize(WTAConeCellsTowersJetsPhi, phi_bit_length_, phi_min_, phi_max_);
                    int eta_bin = digitize(WTAConeCellsTowersJetsEta, eta_bit_length_, eta_min_, eta_max_, etaAltRange);
                    int pt_bin  = digitize(WTAConeCellsTowersJetspT,  et_bit_length_, // digitize the pT the same as E_T would be digitized
                                        static_cast<double>(et_min_), static_cast<double>(et_max_));
                    //std::cout << "WTA pT [GeV]: " << WTAConeCellsTowersJetspT << "\n";
                    //std::cout << "WTA pT bin: " << pt_bin << "\n";
                    //std::cout << "et bit length: " << et_bit_length_ << " , et min: " << et_min_ << " , et max: " << et_max_ << "\n"; 
                    std::stringstream binary_ss;
                    binary_ss << std::bitset<et_bit_length_>(pt_bin)  << "|"
                            << std::bitset<eta_bit_length_>(eta_bin) << "|"
                            << std::bitset<phi_bit_length_>(phi_bin);
                    std::string binary_word = binary_ss.str();

                    uint32_t packed_word = (pt_bin  << (eta_bit_length_ + phi_bit_length_)) |
                                        (eta_bin <<  phi_bit_length_) |
                                        (phi_bin);

                    if (gepwtaconecellstowers_it == 0) {
                        f_wtaconejetscellstowers << "Event : " << std::dec << iEvt << "\n";
                    }
                    f_wtaconejetscellstowers << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepwtaconecellstowers_it
                                        << " "  << binary_word
                                        << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    ++gepwtaconecellstowers_it;

                }
            }

            // GEP cone jets from GEPCellsTowers with SoftKiller PU Suppression applied
            //std::cout << "gepWTAConeCellsTowersSKJetsPt->size(): " << gepWTAConeCellsTowersSKJetsPt->size() << "\n";
            if(gepWTAConeCellsTowersSKJetsPt->size() == 0) std::cout << "NO WTA CONE JETS FOR THIS EVENT!" << "\n";
            {
                std::vector<unsigned int> indices(gepWTAConeCellsTowersSKJetsPt->size());
                std::iota(indices.begin(), indices.end(), 0);

                std::sort(indices.begin(), indices.end(),
                        [&](unsigned int a, unsigned int b) {
                            return (*gepWTAConeCellsTowersSKJetsPt)[a] > (*gepWTAConeCellsTowersSKJetsPt)[b];
                        });
                unsigned int gepwtaconecellstowerssk_it = 0;
                for (unsigned int i = 0; i < indices.size(); i++) {
                    unsigned int idx = indices[i];
                    float WTAConeCellsTowersSKJetspT  = (*gepWTAConeCellsTowersSKJetsPt)[idx] / 1000.0;
                    float WTAConeCellsTowersSKJetsEta = (*gepWTAConeCellsTowersSKJetsEta)[idx];
                    float WTAConeCellsTowersSKJetsPhi = (*gepWTAConeCellsTowersSKJetsPhi)[idx];
                    unsigned int WTAConeCellsTowersSKJetsNConstituents = (*gepWTAConeCellsTowersSKJetsNConstituents)[idx];
                    if(i == 0){
                        gepLeadingWTAConeCellsTowersSKJetspTValues.push_back(WTAConeCellsTowersSKJetspT);
                        gepLeadingWTAConeCellsTowersSKJetsEtaValues.push_back(WTAConeCellsTowersSKJetsEta);
                        gepLeadingWTAConeCellsTowersSKJetsPhiValues.push_back(WTAConeCellsTowersSKJetsPhi);
                        gepLeadingWTAConeCellsTowersSKJetsNConstituentsValues.push_back(WTAConeCellsTowersSKJetsNConstituents);
                    }
                    if (i == 1){
                        gepSubleadingWTAConeCellsTowersSKJetspTValues.push_back(WTAConeCellsTowersSKJetspT);
                        gepSubleadingWTAConeCellsTowersSKJetsEtaValues.push_back(WTAConeCellsTowersSKJetsEta);
                        gepSubleadingWTAConeCellsTowersSKJetsPhiValues.push_back(WTAConeCellsTowersSKJetsPhi);
                        gepSubleadingWTAConeCellsTowersSKJetsNConstituentsValues.push_back(WTAConeCellsTowersSKJetsNConstituents);
                    }
                
                    gepWTAConeCellsTowersSKJetspTValues.push_back(WTAConeCellsTowersSKJetspT);
                    gepWTAConeCellsTowersSKJetsEtaValues.push_back(WTAConeCellsTowersSKJetsEta);
                    gepWTAConeCellsTowersSKJetsPhiValues.push_back(WTAConeCellsTowersSKJetsPhi);
                    gepWTAConeCellsTowersSKJetsNConstituentsValues.push_back(WTAConeCellsTowersSKJetsNConstituents);

                    if (WTAConeCellsTowersSKJetspT <= 0) continue; // don't want to skip over  == 0 Et cone jets 

                    int phi_bin = digitize(WTAConeCellsTowersSKJetsPhi, phi_bit_length_, phi_min_, phi_max_);
                    int eta_bin = digitize(WTAConeCellsTowersSKJetsEta, eta_bit_length_, eta_min_, eta_max_, etaAltRange);
                    //std::cout << "WTAConeCellsTowersSKJetsEta: " << WTAConeCellsTowersSKJetsEta << " , eta_bit_length_: " << eta_bit_length_ << "\n";
                    //std::cout << "eta min: " << eta_min_ << " , eta_max_: " << eta_max_ << " , etaAltRange: " << etaAltRange << "\n";
                    //std::cout << "eta_bin: " << eta_bin << "\n";
                    int pt_bin  = digitize(WTAConeCellsTowersSKJetspT,  et_bit_length_, // digitize the pT the same as E_T would be digitized
                                        static_cast<double>(et_min_), static_cast<double>(et_max_));

                    std::stringstream binary_ss;
                    binary_ss << std::bitset<et_bit_length_>(pt_bin)  << "|"
                            << std::bitset<eta_bit_length_>(eta_bin) << "|"
                            << std::bitset<phi_bit_length_>(phi_bin);
                    std::string binary_word = binary_ss.str();

                    uint32_t packed_word = (pt_bin  << (eta_bit_length_ + phi_bit_length_)) |
                                        (eta_bin <<  phi_bit_length_) |
                                        (phi_bin);


                    if (gepwtaconecellstowerssk_it == 0) {
                        f_wtaconejetscellstowerssk << "Event : " << std::dec << iEvt << "\n";
                    }
                    f_wtaconejetscellstowerssk << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepwtaconecellstowerssk_it
                                        << " "  << binary_word
                                        << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    ++gepwtaconecellstowerssk_it;

                }
            }

            // GEP cone jets from basic clusters
            {
                std::vector<unsigned int> indices(gepWTAConeBasicClustersJetsPt->size());
                std::iota(indices.begin(), indices.end(), 0);

                std::sort(indices.begin(), indices.end(),
                        [&](unsigned int a, unsigned int b) {
                            return (*gepWTAConeBasicClustersJetsPt)[a] > (*gepWTAConeBasicClustersJetsPt)[b];
                        });
                unsigned int gepwtaconebasicclusters_it = 0;
                for (unsigned int i = 0; i < indices.size(); i++) {
                    unsigned int idx = indices[i];
                    float coneWTAGEPBasicClustersJetspT  = (*gepWTAConeBasicClustersJetsPt)[idx] / 1000.0;
                    float coneWTAGEPBasicClustersJetsEta = (*gepWTAConeBasicClustersJetsEta)[idx];
                    float coneWTAGEPBasicClustersJetsPhi = (*gepWTAConeBasicClustersJetsPhi)[idx];
                    unsigned int coneWTAGEPBasicClustersJetsNConstituents = (*gepWTAConeBasicClustersJetsNConstituents)[idx];
                    if(i == 0){
                        gepLeadingWTAConeGEPBasicClustersJetspTValues.push_back(coneWTAGEPBasicClustersJetspT);
                        gepLeadingWTAConeGEPBasicClustersJetsEtaValues.push_back(coneWTAGEPBasicClustersJetsEta);
                        gepLeadingWTAConeGEPBasicClustersJetsPhiValues.push_back(coneWTAGEPBasicClustersJetsPhi);
                        gepLeadingWTAConeGEPBasicClustersJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersJetsNConstituents);
                    }
                    if (i == 1){
                        gepSubleadingWTAConeGEPBasicClustersJetspTValues.push_back(coneWTAGEPBasicClustersJetspT);
                        gepSubleadingWTAConeGEPBasicClustersJetsEtaValues.push_back(coneWTAGEPBasicClustersJetsEta);
                        gepSubleadingWTAConeGEPBasicClustersJetsPhiValues.push_back(coneWTAGEPBasicClustersJetsPhi);
                        gepSubleadingWTAConeGEPBasicClustersJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersJetsNConstituents);
                    }   
                    
                    gepWTAConeGEPBasicClustersJetspTValues.push_back(coneWTAGEPBasicClustersJetspT);
                    gepWTAConeGEPBasicClustersJetsEtaValues.push_back(coneWTAGEPBasicClustersJetsEta);
                    gepWTAConeGEPBasicClustersJetsPhiValues.push_back(coneWTAGEPBasicClustersJetsPhi);
                    gepWTAConeGEPBasicClustersJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersJetsNConstituents);

                    if (coneWTAGEPBasicClustersJetspT <= 0) continue;

                    int phi_bin = digitize(coneWTAGEPBasicClustersJetsPhi, phi_bit_length_, phi_min_, phi_max_);
                    int eta_bin = digitize(coneWTAGEPBasicClustersJetsEta, eta_bit_length_, eta_min_, eta_max_, etaAltRange);
                    int pt_bin  = digitize(coneWTAGEPBasicClustersJetspT,  et_bit_length_, // digitize the pT the same as E_T would be digitized
                                        static_cast<double>(et_min_), static_cast<double>(et_max_));

                    std::stringstream binary_ss;
                    binary_ss << std::bitset<et_bit_length_>(pt_bin)  << "|"
                            << std::bitset<eta_bit_length_>(eta_bin) << "|"
                            << std::bitset<phi_bit_length_>(phi_bin);
                    std::string binary_word = binary_ss.str();

                    uint32_t packed_word = (pt_bin  << (eta_bit_length_ + phi_bit_length_)) |
                                        (eta_bin <<  phi_bit_length_) |
                                        (phi_bin);


                    if (gepwtaconebasicclusters_it == 0) {
                        f_wtaconejetsbasicclusters << "Event : " << std::dec << iEvt << "\n";
                    }
                    f_wtaconejetsbasicclusters << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepwtaconebasicclusters_it
                                        << " "  << binary_word
                                        << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    ++gepwtaconebasicclusters_it;

                }
            }

            // GEP cone jets from basic clusters with SoftKiller PU Suppression applied
            {
                std::vector<unsigned int> indices(gepWTAConeBasicClustersSKJetsPt->size());
                std::iota(indices.begin(), indices.end(), 0);

                std::sort(indices.begin(), indices.end(),
                        [&](unsigned int a, unsigned int b) {
                            return (*gepWTAConeBasicClustersSKJetsPt)[a] > (*gepWTAConeBasicClustersSKJetsPt)[b];
                        });
                unsigned int gepwtaconebasicclusterssk_it = 0;
                for (unsigned int i = 0; i < indices.size(); i++) {
                    unsigned int idx = indices[i];
                    float coneWTAGEPBasicClustersSKJetspT  = (*gepWTAConeBasicClustersSKJetsPt)[idx] / 1000.0;
                    float coneWTAGEPBasicClustersSKJetsEta = (*gepWTAConeBasicClustersSKJetsEta)[idx];
                    float coneWTAGEPBasicClustersSKJetsPhi = (*gepWTAConeBasicClustersSKJetsPhi)[idx];
                    unsigned int coneWTAGEPBasicClustersSKJetsNConstituents = (*gepWTAConeBasicClustersJetsNConstituents)[idx];
                    if(i == 0){
                        gepLeadingWTAConeGEPBasicClustersSKJetspTValues.push_back(coneWTAGEPBasicClustersSKJetspT);
                        gepLeadingWTAConeGEPBasicClustersSKJetsEtaValues.push_back(coneWTAGEPBasicClustersSKJetsEta);
                        gepLeadingWTAConeGEPBasicClustersSKJetsPhiValues.push_back(coneWTAGEPBasicClustersSKJetsPhi);
                        gepLeadingWTAConeGEPBasicClustersSKJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersSKJetsNConstituents);
                    }
                    if (i == 1){
                        gepSubleadingWTAConeGEPBasicClustersSKJetspTValues.push_back(coneWTAGEPBasicClustersSKJetspT);
                        gepSubleadingWTAConeGEPBasicClustersSKJetsEtaValues.push_back(coneWTAGEPBasicClustersSKJetsEta);
                        gepSubleadingWTAConeGEPBasicClustersSKJetsPhiValues.push_back(coneWTAGEPBasicClustersSKJetsPhi);
                        gepSubleadingWTAConeGEPBasicClustersSKJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersSKJetsNConstituents);
                    }   
                    
                    gepWTAConeGEPBasicClustersSKJetspTValues.push_back(coneWTAGEPBasicClustersSKJetspT);
                    gepWTAConeGEPBasicClustersSKJetsEtaValues.push_back(coneWTAGEPBasicClustersSKJetsEta);
                    gepWTAConeGEPBasicClustersSKJetsPhiValues.push_back(coneWTAGEPBasicClustersSKJetsPhi);
                    gepWTAConeGEPBasicClustersSKJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersSKJetsNConstituents);

                    if (coneWTAGEPBasicClustersSKJetspT <= 0) continue;

                    int phi_bin = digitize(coneWTAGEPBasicClustersSKJetsPhi, phi_bit_length_, phi_min_, phi_max_);
                    int eta_bin = digitize(coneWTAGEPBasicClustersSKJetsEta, eta_bit_length_, eta_min_, eta_max_, etaAltRange);
                    int pt_bin  = digitize(coneWTAGEPBasicClustersSKJetspT,  et_bit_length_, // digitize the pT the same as E_T would be digitized
                                        static_cast<double>(et_min_), static_cast<double>(et_max_));

                    std::stringstream binary_ss;
                    binary_ss << std::bitset<et_bit_length_>(pt_bin)  << "|"
                            << std::bitset<eta_bit_length_>(eta_bin) << "|"
                            << std::bitset<phi_bit_length_>(phi_bin);
                    std::string binary_word = binary_ss.str();

                    uint32_t packed_word = (pt_bin  << (eta_bit_length_ + phi_bit_length_)) |
                                        (eta_bin <<  phi_bit_length_) |
                                        (phi_bin);

                    if (gepwtaconebasicclusterssk_it == 0) {
                        f_wtaconejetsbasicclusterssk << "Event : " << std::dec << iEvt << "\n";
                    }
                    f_wtaconejetsbasicclusterssk << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepwtaconebasicclusterssk_it
                                        << " "  << binary_word
                                        << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    ++gepwtaconebasicclusterssk_it;
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

                if (et <= 0) continue; // FIXME

                // Digitize each variable
                int phi_bin = digitize(jet->phi(), phi_bit_length_, phi_min_, phi_max_);
                int eta_bin = digitize(jet->eta(), eta_bit_length_, eta_min_, eta_max_, etaAltRange);
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
                if (jfex_it == 0) {
                    f_jfex << "Event : " << std::dec << iEvt << "\n";
                }
                f_jfex << "0x" << std::hex << std::setw(2) << std::setfill('0') << jfex_it
                        << " " << binary_word
                        << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                jfex_it++;

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

                if (et <= 0) continue;

                // Digitize each variable
                int phi_bin = digitize(jet->phi(), phi_bit_length_, phi_min_, phi_max_);
                int eta_bin = digitize(jet->eta(), eta_bit_length_, eta_min_, eta_max_, etaAltRange);
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
                if (gfex_it == 0) {
                    f_gfex << "Event : " << std::dec << iEvt << "\n";
                }
                f_gfex << "0x" << std::hex << std::setw(2) << std::setfill('0') << gfex_it
                        << " " << binary_word
                        << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                gfex_it++;
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

            // --- gFEX MET ---
            gMHX = 0.0; gMHY = 0.0;
            gMSX = 0.0; gMSY = 0.0;
            gMETX = 0.0; gMETY = 0.0;
            gMET = 0.0; gSumET = 0.0;
            gMETX_NC = 0.0; gMETY_NC = 0.0; gMET_NC = 0.0; gSumET_NC = 0.0;
            gMETX_Rms = 0.0; gMETY_Rms = 0.0; gMET_Rms = 0.0; gSumET_Rms = 0.0;
            metTruthNonIntX = 0.0; metTruthNonIntY = 0.0; metTruthNonIntSumET = 0.0; metTruthNonInt = 0.0;
            metTruthIntX    = 0.0; metTruthIntY    = 0.0; metTruthIntSumET    = 0.0; metTruthInt    = 0.0;
            metTruthIntOutX = 0.0; metTruthIntOutY = 0.0; metTruthIntOutSumET = 0.0; metTruthIntOut = 0.0;
            metCoreEMTopo_SoftClus_X   = 0.0; metCoreEMTopo_SoftClus_Y   = 0.0; metCoreEMTopo_SoftClus_SumET   = 0.0; metCoreEMTopo_SoftClus_MET   = 0.0;
            metCoreEMTopo_PVSoftTrk_X  = 0.0; metCoreEMTopo_PVSoftTrk_Y  = 0.0; metCoreEMTopo_PVSoftTrk_SumET  = 0.0; metCoreEMTopo_PVSoftTrk_MET  = 0.0;
            metCoreEMTopo_SoftClusEM_X = 0.0; metCoreEMTopo_SoftClusEM_Y = 0.0; metCoreEMTopo_SoftClusEM_SumET = 0.0; metCoreEMTopo_SoftClusEM_MET = 0.0;
            metCoreEMPFlow_SoftClus_X  = 0.0; metCoreEMPFlow_SoftClus_Y  = 0.0; metCoreEMPFlow_SoftClus_SumET  = 0.0; metCoreEMPFlow_SoftClus_MET  = 0.0;
            metCoreEMPFlow_PVSoftTrk_X = 0.0; metCoreEMPFlow_PVSoftTrk_Y = 0.0; metCoreEMPFlow_PVSoftTrk_SumET = 0.0; metCoreEMPFlow_PVSoftTrk_MET = 0.0;

            if (L1_gMHTComponentsJwoj) {
                for (size_t i = 0; i < L1_gMHTComponentsJwoj->size(); ++i) {
                    const auto& mht = (*L1_gMHTComponentsJwoj)[i];
                    gMHX = mht->METquantityOne() / 1000.0;
                    gMHY = mht->METquantityTwo() / 1000.0;
                }
            }
            if (L1_gMSTComponentsJwoj) {
                for (size_t i = 0; i < L1_gMSTComponentsJwoj->size(); ++i) {
                    const auto& mst = (*L1_gMSTComponentsJwoj)[i];
                    gMSX = mst->METquantityOne() / 1000.0;
                    gMSY = mst->METquantityTwo() / 1000.0;
                }
            }
            if (L1_gMETComponentsJwoj) {
                for (size_t i = 0; i < L1_gMETComponentsJwoj->size(); ++i) {
                    const auto& metTotal = (*L1_gMETComponentsJwoj)[i];
                    gMETX = metTotal->METquantityOne() / 1000.0;
                    gMETY = metTotal->METquantityTwo() / 1000.0;
                }
            }
            if (L1_gScalarEJwoj) {
                for (size_t i = 0; i < L1_gScalarEJwoj->size(); ++i) {
                    const auto& metScalar = (*L1_gScalarEJwoj)[i];
                    gMET   = metScalar->METquantityOne() / 1000.0;
                    gSumET = metScalar->SumEt() / 1000.0;
                }
            }

            // --- gFEX MET NoiseCut ---
            if (L1_gMETComponentsNoiseCut) {
                for (size_t i = 0; i < L1_gMETComponentsNoiseCut->size(); ++i) {
                    const auto& metTotal = (*L1_gMETComponentsNoiseCut)[i];
                    gMETX_NC = metTotal->METquantityOne() / 1000.0;
                    gMETY_NC = metTotal->METquantityTwo() / 1000.0;
                }
            }
            gMET_NC = std::sqrt(gMETX_NC * gMETX_NC + gMETY_NC * gMETY_NC);
            if (L1_gScalarENoiseCut) {
                for (size_t i = 0; i < L1_gScalarENoiseCut->size(); ++i) {
                    gSumET_NC = (*L1_gScalarENoiseCut)[i]->SumEt() / 1000.0;
                }
            }

            // --- gFEX MET Rms ---
            if (L1_gMETComponentsRms) {
                for (size_t i = 0; i < L1_gMETComponentsRms->size(); ++i) {
                    const auto& metTotal = (*L1_gMETComponentsRms)[i];
                    gMETX_Rms = metTotal->METquantityOne() / 1000.0;
                    gMETY_Rms = metTotal->METquantityTwo() / 1000.0;
                }
            }
            gMET_Rms = std::sqrt(gMETX_Rms * gMETX_Rms + gMETY_Rms * gMETY_Rms);
            if (L1_gScalarERms) {
                for (size_t i = 0; i < L1_gScalarERms->size(); ++i) {
                    gSumET_Rms = (*L1_gScalarERms)[i]->SumEt() / 1000.0;
                }
            }

            // --- MET_Truth ---
            if (MET_Truth) {
                auto fillTruthTerm = [&](const std::string& key,
                                         double& mx, double& my, double& sumet, double& mag) {
                    auto it = MET_Truth->find(key);
                    if (it != MET_Truth->end()) {
                        const xAOD::MissingET* m = *it;
                        mx    = m->mpx()   / 1000.0;
                        my    = m->mpy()   / 1000.0;
                        sumet = m->sumet() / 1000.0;
                        mag   = std::sqrt(mx * mx + my * my);
                    }
                };
                fillTruthTerm("NonInt", metTruthNonIntX, metTruthNonIntY, metTruthNonIntSumET, metTruthNonInt);
                fillTruthTerm("Int",    metTruthIntX,    metTruthIntY,    metTruthIntSumET,    metTruthInt);
                fillTruthTerm("IntOut", metTruthIntOutX, metTruthIntOutY, metTruthIntOutSumET, metTruthIntOut);
            }

            // Helper: fill per-term variables from a named term in a MissingET container
            auto fillCoreTerm = [](const xAOD::MissingETContainer* cont, const std::string& key,
                                   double& mx, double& my, double& sumet, double& mag, Long64_t iEvt) {
                auto it = cont->find(key);
                if (it != cont->end()) {
                    const xAOD::MissingET* m = *it;
                    mx    = m->mpx()   / 1000.0;
                    my    = m->mpy()   / 1000.0;
                    sumet = m->sumet() / 1000.0;
                    mag   = std::sqrt(mx * mx + my * my);
                    if(iEvt <= 10) std::cout << key << ": X=" << mx << " Y=" << my << " SumET=" << sumet << " MET=" << mag << "\n";
                }
            };

            // --- MET_Core_AntiKt4EMTopo ---
            if (MET_Core_AntiKt4EMTopo) {
                fillCoreTerm(MET_Core_AntiKt4EMTopo, "SoftClusCore",   metCoreEMTopo_SoftClus_X,   metCoreEMTopo_SoftClus_Y,   metCoreEMTopo_SoftClus_SumET,   metCoreEMTopo_SoftClus_MET, iEvt);
                fillCoreTerm(MET_Core_AntiKt4EMTopo, "PVSoftTrkCore",  metCoreEMTopo_PVSoftTrk_X,  metCoreEMTopo_PVSoftTrk_Y,  metCoreEMTopo_PVSoftTrk_SumET,  metCoreEMTopo_PVSoftTrk_MET, iEvt);
                fillCoreTerm(MET_Core_AntiKt4EMTopo, "SoftClusEMCore", metCoreEMTopo_SoftClusEM_X, metCoreEMTopo_SoftClusEM_Y, metCoreEMTopo_SoftClusEM_SumET, metCoreEMTopo_SoftClusEM_MET, iEvt);
            }

            // --- MET_Core_AntiKt4EMPFlow ---
            if (MET_Core_AntiKt4EMPFlow) {
                fillCoreTerm(MET_Core_AntiKt4EMPFlow, "SoftClusCore",  metCoreEMPFlow_SoftClus_X,  metCoreEMPFlow_SoftClus_Y,  metCoreEMPFlow_SoftClus_SumET,  metCoreEMPFlow_SoftClus_MET, iEvt);
                fillCoreTerm(MET_Core_AntiKt4EMPFlow, "PVSoftTrkCore", metCoreEMPFlow_PVSoftTrk_X, metCoreEMPFlow_PVSoftTrk_Y, metCoreEMPFlow_PVSoftTrk_SumET, metCoreEMPFlow_PVSoftTrk_MET, iEvt);
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

            // Temporary vector to hold (index, Et) for sorting
            std::vector<std::pair<size_t, double>> recoJetEtWithIndex;

            for (size_t i = 0; i < AntiKt10UFOCSSKJets->size(); ++i) {
                const auto& jet = (*AntiKt10UFOCSSKJets)[i];
                double et = jet->e() / (1000.0 * cosh(jet->eta()));
                recoJetEtWithIndex.emplace_back(i, et);
            }

            // Sort by descending Et
            std::sort(recoJetEtWithIndex.begin(), recoJetEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });

            // Fill vectors in sorted order
            for (const auto& [index, et] : recoJetEtWithIndex) {
                const auto& jet = (*AntiKt10UFOCSSKJets)[index];
                recoAntiKt10LRJEtIndexValues.push_back(static_cast<unsigned int>(index));
                recoAntiKt10LRJEtValues.push_back(et);
                recoAntiKt10LRJEtaValues.push_back(jet->eta());
                recoAntiKt10LRJPhiValues.push_back(jet->phi());
                recoAntiKt10LRJMassValues.push_back(jet->m() / 1000.0);
            }

            // Leading jet
            if (!recoJetEtWithIndex.empty()) {
                const auto& leading = (*AntiKt10UFOCSSKJets)[recoJetEtWithIndex[0].first];
                recoAntiKt10LRJLeadingEtValues.push_back(recoAntiKt10LRJEtValues[0]);
                recoAntiKt10LRJLeadingEtaValues.push_back(leading->eta());
                recoAntiKt10LRJLeadingPhiValues.push_back(leading->phi());
                recoAntiKt10LRJLeadingMassValues.push_back(leading->m() / 1000.0);
            }

            // Subleading jet
            if (recoJetEtWithIndex.size() > 1) {
                const auto& subleading = (*AntiKt10UFOCSSKJets)[recoJetEtWithIndex[1].first];
                recoAntiKt10LRJSubleadingEtValues.push_back(recoAntiKt10LRJEtValues[1]);
                recoAntiKt10LRJSubleadingEtaValues.push_back(subleading->eta());
                recoAntiKt10LRJSubleadingPhiValues.push_back(subleading->phi());
                recoAntiKt10LRJSubleadingMassValues.push_back(subleading->m() / 1000.0);
            }

            // AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets
            std::vector<std::pair<size_t, double>> recoSDJetEtWithIndex;
            for (size_t i = 0; i < AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets->size(); ++i) {
                const auto& jet = (*AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets)[i];
                double et = jet->e() / (1000.0 * cosh(jet->eta()));
                recoSDJetEtWithIndex.emplace_back(i, et);
            }
            std::sort(recoSDJetEtWithIndex.begin(), recoSDJetEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });
            for (const auto& [index, et] : recoSDJetEtWithIndex) {
                const auto& jet = (*AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets)[index];
                recoAntiKt10UFOCSSKSDEtIndexValues.push_back(static_cast<unsigned int>(index));
                recoAntiKt10UFOCSSKSDEtValues.push_back(et);
                recoAntiKt10UFOCSSKSDEtaValues.push_back(jet->eta());
                recoAntiKt10UFOCSSKSDPhiValues.push_back(jet->phi());
                recoAntiKt10UFOCSSKSDMassValues.push_back(jet->m() / 1000.0);
            }
            if (!recoSDJetEtWithIndex.empty()) {
                const auto& leading = (*AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets)[recoSDJetEtWithIndex[0].first];
                recoAntiKt10UFOCSSKSDLeadingEtValues.push_back(recoAntiKt10UFOCSSKSDEtValues[0]);
                recoAntiKt10UFOCSSKSDLeadingEtaValues.push_back(leading->eta());
                recoAntiKt10UFOCSSKSDLeadingPhiValues.push_back(leading->phi());
                recoAntiKt10UFOCSSKSDLeadingMassValues.push_back(leading->m() / 1000.0);
            }
            if (recoSDJetEtWithIndex.size() > 1) {
                const auto& subleading = (*AntiKt10UFOCSSKSoftDropBeta100Zcut10Jets)[recoSDJetEtWithIndex[1].first];
                recoAntiKt10UFOCSSKSDSubleadingEtValues.push_back(recoAntiKt10UFOCSSKSDEtValues[1]);
                recoAntiKt10UFOCSSKSDSubleadingEtaValues.push_back(subleading->eta());
                recoAntiKt10UFOCSSKSDSubleadingPhiValues.push_back(subleading->phi());
                recoAntiKt10UFOCSSKSDSubleadingMassValues.push_back(subleading->m() / 1000.0);
            }

            // AntiKt10TruthJets
            std::vector<std::pair<size_t, double>> truthLRJEtWithIndex;
            for (size_t i = 0; i < AntiKt10TruthJets->size(); ++i) {
                const auto& jet = (*AntiKt10TruthJets)[i];
                double et = jet->e() / (1000.0 * cosh(jet->eta()));
                truthLRJEtWithIndex.emplace_back(i, et);
            }
            std::sort(truthLRJEtWithIndex.begin(), truthLRJEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });
            for (const auto& [index, et] : truthLRJEtWithIndex) {
                const auto& jet = (*AntiKt10TruthJets)[index];
                antiKt10TruthEtIndexValues.push_back(static_cast<unsigned int>(index));
                antiKt10TruthEtValues.push_back(et);
                antiKt10TruthEtaValues.push_back(jet->eta());
                antiKt10TruthPhiValues.push_back(jet->phi());
                antiKt10TruthMassValues.push_back(jet->m() / 1000.0);
            }
            if (!truthLRJEtWithIndex.empty()) {
                const auto& leading = (*AntiKt10TruthJets)[truthLRJEtWithIndex[0].first];
                antiKt10TruthLeadingEtValues.push_back(antiKt10TruthEtValues[0]);
                antiKt10TruthLeadingEtaValues.push_back(leading->eta());
                antiKt10TruthLeadingPhiValues.push_back(leading->phi());
                antiKt10TruthLeadingMassValues.push_back(leading->m() / 1000.0);
            }
            if (truthLRJEtWithIndex.size() > 1) {
                const auto& subleading = (*AntiKt10TruthJets)[truthLRJEtWithIndex[1].first];
                antiKt10TruthSubleadingEtValues.push_back(antiKt10TruthEtValues[1]);
                antiKt10TruthSubleadingEtaValues.push_back(subleading->eta());
                antiKt10TruthSubleadingPhiValues.push_back(subleading->phi());
                antiKt10TruthSubleadingMassValues.push_back(subleading->m() / 1000.0);
            }

            // AntiKt10TruthSoftDropBeta100Zcut10Jets
            std::vector<std::pair<size_t, double>> truthSDLRJEtWithIndex;
            for (size_t i = 0; i < AntiKt10TruthSoftDropBeta100Zcut10Jets->size(); ++i) {
                const auto& jet = (*AntiKt10TruthSoftDropBeta100Zcut10Jets)[i];
                double et = jet->e() / (1000.0 * cosh(jet->eta()));
                truthSDLRJEtWithIndex.emplace_back(i, et);
            }
            std::sort(truthSDLRJEtWithIndex.begin(), truthSDLRJEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });
            for (const auto& [index, et] : truthSDLRJEtWithIndex) {
                const auto& jet = (*AntiKt10TruthSoftDropBeta100Zcut10Jets)[index];
                antiKt10TruthSDEtIndexValues.push_back(static_cast<unsigned int>(index));
                antiKt10TruthSDEtValues.push_back(et);
                antiKt10TruthSDEtaValues.push_back(jet->eta());
                antiKt10TruthSDPhiValues.push_back(jet->phi());
                antiKt10TruthSDMassValues.push_back(jet->m() / 1000.0);
            }
            if (!truthSDLRJEtWithIndex.empty()) {
                const auto& leading = (*AntiKt10TruthSoftDropBeta100Zcut10Jets)[truthSDLRJEtWithIndex[0].first];
                antiKt10TruthSDLeadingEtValues.push_back(antiKt10TruthSDEtValues[0]);
                antiKt10TruthSDLeadingEtaValues.push_back(leading->eta());
                antiKt10TruthSDLeadingPhiValues.push_back(leading->phi());
                antiKt10TruthSDLeadingMassValues.push_back(leading->m() / 1000.0);
            }
            if (truthSDLRJEtWithIndex.size() > 1) {
                const auto& subleading = (*AntiKt10TruthSoftDropBeta100Zcut10Jets)[truthSDLRJEtWithIndex[1].first];
                antiKt10TruthSDSubleadingEtValues.push_back(antiKt10TruthSDEtValues[1]);
                antiKt10TruthSDSubleadingEtaValues.push_back(subleading->eta());
                antiKt10TruthSDSubleadingPhiValues.push_back(subleading->phi());
                antiKt10TruthSDSubleadingMassValues.push_back(subleading->m() / 1000.0);
            }

            // Temporary vector to hold (index, Et) for sorting
            std::vector<std::pair<size_t, double>> truthWZJetEtWithIndex;

            for (size_t i = 0; i < AntiKt4TruthDressedWZJets->size(); ++i) {
                const auto& jet = (*AntiKt4TruthDressedWZJets)[i];
                double et = jet->e() / (1000.0 * cosh(jet->eta()));
                truthWZJetEtWithIndex.emplace_back(i, et);
            }

            // Sort by descending Et
            std::sort(truthWZJetEtWithIndex.begin(), truthWZJetEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });

            // Fill vectors in sorted order
            for (const auto& [index, et] : truthWZJetEtWithIndex) {
                const auto& jet = (*AntiKt4TruthDressedWZJets)[index];
                truthAntiKt4WZSRJEtIndexValues.push_back(static_cast<unsigned int>(index));
                truthAntiKt4WZSRJEtValues.push_back(et);
                truthAntiKt4WZSRJEtaValues.push_back(jet->eta());
                truthAntiKt4WZSRJPhiValues.push_back(jet->phi());
                truthAntiKt4WZSRJMassValues.push_back(jet->m() / 1000.0);
            }

            // Leading jet
            if (!truthWZJetEtWithIndex.empty()) {
                const auto& leading = (*AntiKt4TruthDressedWZJets)[truthWZJetEtWithIndex[0].first];
                truthAntiKt4WZSRJLeadingEtValues.push_back(truthAntiKt4WZSRJEtValues[0]);
                truthAntiKt4WZSRJLeadingEtaValues.push_back(leading->eta());
                truthAntiKt4WZSRJLeadingPhiValues.push_back(leading->phi());
                truthAntiKt4WZSRJLeadingMassValues.push_back(leading->m() / 1000.0);

            }

            // Subleading jet
            if (truthWZJetEtWithIndex.size() > 1) {
                const auto& subleading = (*AntiKt4TruthDressedWZJets)[truthWZJetEtWithIndex[1].first];
                truthAntiKt4WZSRJSubleadingEtValues.push_back(truthAntiKt4WZSRJEtValues[1]);
                truthAntiKt4WZSRJSubleadingEtaValues.push_back(subleading->eta());
                truthAntiKt4WZSRJSubleadingPhiValues.push_back(subleading->phi());
                truthAntiKt4WZSRJSubleadingMassValues.push_back(subleading->m() / 1000.0);
            }
            sampleJZSlice = jzSlice;
            eventInfoTree->Fill();
            truthbTree->Fill();
            truthHiggsTree->Fill();
            // truthVBFQuark->Fill();  // commented out as in your declaration
            caloTopoTowerTree->Fill();
            gepBasicClustersTree->Fill();
            gepCellsTowersTree->Fill();
            gepBasicClustersSKTree->Fill();
            gepCellsTowersSKTree->Fill();
            gepWTAConeCellsTowersJetsTree->Fill();
            gepWTAConeBasicClustersJetsTree->Fill();
            gepLeadingWTAConeCellsTowersJetsTree->Fill();
            gepLeadingWTAConeBasicClustersJetsTree->Fill();
            gepSubleadingWTAConeCellsTowersJetsTree->Fill();
            gepSubleadingWTAConeBasicClustersJetsTree->Fill();
            gepWTAConeCellsTowersSKJetsTree->Fill();
            gepWTAConeBasicClustersSKJetsTree->Fill();
            gepLeadingWTAConeCellsTowersSKJetsTree->Fill();
            gepLeadingWTAConeBasicClustersSKJetsTree->Fill();
            gepSubleadingWTAConeCellsTowersSKJetsTree->Fill();
            gepSubleadingWTAConeBasicClustersSKJetsTree->Fill();
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
            recoAntiKt10UFOCSSKSoftDropJets->Fill();
            leadingRecoAntiKt10UFOCSSKSoftDropJets->Fill();
            subleadingRecoAntiKt10UFOCSSKSoftDropJets->Fill();
            antiKt10TruthJetsTree->Fill();
            leadingAntiKt10TruthJetsTree->Fill();
            subleadingAntiKt10TruthJetsTree->Fill();
            antiKt10TruthSoftDropJetsTree->Fill();
            leadingAntiKt10TruthSoftDropJetsTree->Fill();
            subleadingAntiKt10TruthSoftDropJetsTree->Fill();
            truthAntiKt4TruthDressedWZJets->Fill();
            leadingTruthAntiKt4TruthDressedWZJets->Fill();
            subleadingTruthAntiKt4TruthDressedWZJets->Fill();
            gFexMETTree->Fill();
            gFexMETNoiseCutTree->Fill();
            gFexMETRmsTree->Fill();
            metTruthTree->Fill();
            metCoreAntiKt4EMTopoTree->Fill();
            metCoreAntiKt4EMPFlowTree->Fill();
        } // loop through events
        std::cout << "for jz: " << jzSlice << " these many events passed: " << passedEventsCounter << " out of: " << event.getEntries() << "\n";
        std::cout << " these many events skipped due to empty truth container: " << skippedEventsEmptyTruth << " these many events skipped due to pt hard < pt pileup: " << skippedEventsHSTP << "\n";
        //std::cout << "leading not 0 truth: " << truthLeadingNot0 << " , PU leading not 0: " << puLeadingNot0 << "\n";
        std::cout << "closing files" << "\n";
        gf->Close();
        f->Close();
    } // loop through filenames
    outputFile->cd();
    std::cout << "writing output file" << "\n";
    eventInfoTree->Write("", TObject::kOverwrite);
    truthbTree->Write("", TObject::kOverwrite);
    truthHiggsTree->Write("", TObject::kOverwrite);
    // truthVBFQuark->Write();  // Optional, if used
    caloTopoTowerTree->Write("", TObject::kOverwrite);
    gepBasicClustersTree->Write("", TObject::kOverwrite);
    gepBasicClustersSKTree->Write("", TObject::kOverwrite);
    gepCellsTowersTree->Write("", TObject::kOverwrite);
    gepCellsTowersSKTree->Write("", TObject::kOverwrite);
    gepWTAConeCellsTowersJetsTree->Write("", TObject::kOverwrite);
    gepWTAConeBasicClustersJetsTree->Write("", TObject::kOverwrite);
    gepWTAConeCellsTowersSKJetsTree->Write("", TObject::kOverwrite);
    gepWTAConeBasicClustersSKJetsTree->Write("", TObject::kOverwrite);
    gepLeadingWTAConeCellsTowersJetsTree->Write("", TObject::kOverwrite);
    gepLeadingWTAConeBasicClustersJetsTree->Write("", TObject::kOverwrite);
    gepLeadingWTAConeCellsTowersSKJetsTree->Write("", TObject::kOverwrite);
    gepLeadingWTAConeBasicClustersSKJetsTree->Write("", TObject::kOverwrite);
    gepSubleadingWTAConeCellsTowersJetsTree->Write("", TObject::kOverwrite);
    gepSubleadingWTAConeBasicClustersJetsTree->Write("", TObject::kOverwrite);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Write("", TObject::kOverwrite);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Write("", TObject::kOverwrite);
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
    recoAntiKt10UFOCSSKSoftDropJets->Write("", TObject::kOverwrite);
    leadingRecoAntiKt10UFOCSSKSoftDropJets->Write("", TObject::kOverwrite);
    subleadingRecoAntiKt10UFOCSSKSoftDropJets->Write("", TObject::kOverwrite);
    antiKt10TruthJetsTree->Write("", TObject::kOverwrite);
    leadingAntiKt10TruthJetsTree->Write("", TObject::kOverwrite);
    subleadingAntiKt10TruthJetsTree->Write("", TObject::kOverwrite);
    antiKt10TruthSoftDropJetsTree->Write("", TObject::kOverwrite);
    leadingAntiKt10TruthSoftDropJetsTree->Write("", TObject::kOverwrite);
    subleadingAntiKt10TruthSoftDropJetsTree->Write("", TObject::kOverwrite);
    truthAntiKt4TruthDressedWZJets->Write("", TObject::kOverwrite);
    leadingTruthAntiKt4TruthDressedWZJets->Write("", TObject::kOverwrite);
    subleadingTruthAntiKt4TruthDressedWZJets->Write("", TObject::kOverwrite);
    gFexMETTree->Write("", TObject::kOverwrite);
    gFexMETNoiseCutTree->Write("", TObject::kOverwrite);
    gFexMETRmsTree->Write("", TObject::kOverwrite);
    metTruthTree->Write("", TObject::kOverwrite);
    metCoreAntiKt4EMTopoTree->Write("", TObject::kOverwrite);
    metCoreAntiKt4EMPFlowTree->Write("", TObject::kOverwrite);
    outputFile->Close();
    std::cout << "Processing complete." << endl;
    //std::cout << "higgsPassEventCounter: " << higgsPassEventCounter << "\n";
} // ntupler function

// processes all jzSlices + signal
void HERNTupler(){
    //gSystem->Load("libxAODRootAccess");
    //xAOD::Init().ignore();
    

    gSystem->RedirectOutput("HERNTupler.log", "w");
    //std::vector<unsigned int > jzSlices = {3};
    std::vector<unsigned int > jzSlices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<unsigned int > algoVersions = {3}; // FIXME algo version needs to control eta, phi bit widths too! otherwise can't loop through them and do this all in one run
    //std::vector<unsigned int > jzSlices = {3};
    std::vector<std::string > jzOutputFilenames;

    for(auto algoVersion : algoVersions){
        unsigned int etaAltRange = 0;
        if(algoVersion == 2){
            etaAltRange = 784;
        }
        if(algoVersion == 3){
            etaAltRange = 98;
        }
        //nTupler(true, "VBF_hh_bbbb", etaAltRange, algoVersion); // call for signal (hh->4b VBF)
        //nTupler(true, "ggF_hh_bbbb", etaAltRange, algoVersion); // call for signal (hh->4b ggF)
        nTupler(true, "ZvvHbb", etaAltRange, algoVersion); // call for signal Z -> nu nu H -> bb
        //nTupler(true, "ttbar_had", etaAltRange, algoVersion); // call for sigal ttbar fully hadronic decays
        //nTupler(true, "Zprime_ttbar", etaAltRange, algoVersion); // call for signal Zprime -> ttbar (hadronic), flat in pT
    }

    for(auto jzSlice : jzSlices){
        for(auto algoVersion : algoVersions){
            unsigned int etaAltRange = 0;
            if(algoVersion == 2){
                etaAltRange = 784; // = 9.8 / 0.0125
            }
            if(algoVersion == 3){
                etaAltRange = 98; // = 9.8 / 0.1
            }
            nTupler(false, "", etaAltRange, algoVersion, jzSlice); // call for each background jzSlice
            TString out = makeOutputFileName("", false, jzSlice, algoVersion);
            jzOutputFilenames.push_back(out.Data());
        }
        
    }

    // FIXME allow hadd to work automatically 
    //std::ostringstream cmd;
    //cmd << "hadd -f -k outputRootFiles/background_merged_nojz7or9.root";
    //for (const auto& f : jzOutputFilenames) cmd << " '" << f << "'";
    //int rc = gSystem->Exec(cmd.str().c_str());
    //if (rc != 0) {
    //    Error("callNTupler", "hadd failed with code %d", rc);
    //}
    
    for(auto jzOutputFileName : jzOutputFilenames){ // for manual hadd'ing into a combined JZ slice ntuple for now
        std::cout << jzOutputFileName << "\n";
    }
    gSystem->Exit(0);
}
