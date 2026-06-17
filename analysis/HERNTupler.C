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
const double sumOfEventWeightsByJZSlice[nJZSlices] = {100000.0,                // JZ0
                                                      4692.711304682304,      // JZ1
                                                      40.668641448125456,     // JZ2
                                                      0.812919741273701,    // JZ3
                                                      0.012162432307614823,  // JZ4
                                                      0.00094489973084666, // JZ5
                                                      0.0001659370012750544,  // JZ6
                                                      4.6604455941866296e-05,  // JZ7 
                                                      1.3833386023877348e-05, // JZ8
                                                      3.292773754935112e-08}; // JZ9 

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
        else if (signalString == "ggF_hh_bbbb_resim") return "/data/larsonma/GEPHadronicEventReconstruction/DAOD_TrigGepPerf/ggF_HHbbbb";
        else if (signalString == "ZvvHbb") return kRoot + "/DAOD_TrigGepPerf/Signal_ZvvHbb";
        else if (signalString == "ZvvHbb_resim") return "/data/larsonma/GEPHadronicEventReconstruction/DAOD_TrigGepPerf/ZvvHbb";
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
TString makeOutputFileName(std::string signalString, bool signalBool, int jzSlice, unsigned int algoVersion, bool specialJZ0 = false, std::string outputDirOverride = "", std::string fileSuffix = "") {
    TString outDir;
    if (!outputDirOverride.empty()) {
        outDir = TString(outputDirOverride);
        if (outDir.Length() > 0 && outDir[outDir.Length()-1] != '/') outDir += "/";
    } else if(algoVersion == 2) outDir = "outputRootFiles/v2/";
    else if(algoVersion == 3) outDir = "outputRootFiles/v3/";

    TString suffix = TString(fileSuffix);
    if(signalBool){
        if(signalString == "VBF_hh_bbbb") return outDir + "mc21_14TeV_hh_bbbb_vbf_novhh_e8557_s4422_r16130_DAOD_NTUPLE_GEP" + suffix + ".root";
        else if (signalString == "ggF_hh_bbbb") return outDir + "mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_DAOD_NTUPLE_GEP" + suffix + ".root";
        else if (signalString == "ggF_hh_bbbb_resim") return outDir + "mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_resim_DAOD_NTUPLE_GEP" + suffix + ".root";
        else if (signalString == "ZvvHbb") return outDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_DAOD_NTUPLE_GEP" + suffix + ".root";
        else if (signalString == "ZvvHbb_resim") return outDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_resim_DAOD_NTUPLE_GEP" + suffix + ".root";
        else if (signalString == "ttbar_had") return outDir + "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422_r16130_DAOD_NTUPLE_GEP" + suffix + ".root";
        else if (signalString == "Zprime_ttbar") return outDir + "mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422_r16130_DAOD_NTUPLE_GEP" + suffix + ".root";
        else return "";
    }
    else {
        if (jzSlice < 0 || jzSlice > 9)
            throw std::out_of_range("jzSlice must be between 0 and 9");
        if(specialJZ0){
            return outDir + "mc21_14TeV_jj_noHSTP_JZ" + TString(std::to_string(jzSlice)) + "_e8557_s4422_r16130_DAOD_NTUPLE_GEP" + suffix + ".root";
        }
        else{
            return outDir + "mc21_14TeV_jj_JZ" + TString(std::to_string(jzSlice)) + "_e8557_s4422_r16130_DAOD_NTUPLE_GEP" + suffix + ".root";
        }
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
void nTupler(bool signalBool, std::string signalString, unsigned int etaAltRange, unsigned int algoVersion, unsigned int jzSlice = 3, bool specialJZ0Bool = false, std::string outputDirOverride = "", bool writeDatFiles = true, std::string specificDaodFile = "", std::string specificGepFile = "", std::string fileSuffix = "") {
    //ServiceHandle<StoreGateSvc> inputMetaStore("StoreGateSvc/InputMetaDataStore","");
    // Setup file paths based on whether processing signal or background, and vbf production or ggF production
    
    std::string fileDir =  makeFileDir(signalString, signalBool, jzSlice);
    std::cout << "fileDir: " << fileDir << "\n";
    //xAOD::Init().ignore();

    // Max number of events to write test vectors for
    const unsigned int maxDatEvents = 1000;

    // Declare .dat output streams (only opened when writeDatFiles == true)
    std::ofstream f_topotower, f_gepbasicclusters, f_gepcellstowers;
    std::ofstream f_gepbasicclusterssk, f_gepcellstowerssk, f_gepcellstowers_sort;
    std::ofstream f_wtaconejetscellstowers, f_wtaconejetsbasicclusters;
    std::ofstream f_wtaconejetscellstowerssk, f_wtaconejetsbasicclusterssk;
    std::ofstream f_topo, f_gfex, f_jfex;

    if (writeDatFiles) {
        OutputFiles out = makeMemPrintFilenames(signalString, signalBool, jzSlice, algoVersion);
        f_topotower.open(out.caloTopoTowers);
        f_gepbasicclusters.open(out.gepBasicClusters);
        f_gepcellstowers.open(out.gepCellTowers);
        f_gepbasicclusterssk.open(out.gepBasicClustersSK);
        f_gepcellstowerssk.open(out.gepCellTowersSK);
        f_gepcellstowers_sort.open(out.gepCellTowersSort);
        f_wtaconejetscellstowers.open(out.gepWTAConeJetsCellsTowers);
        f_wtaconejetsbasicclusters.open(out.gepWTAConeJetsBasicClusters);
        f_wtaconejetscellstowerssk.open(out.gepWTAConeJetsCellsTowersSK);
        f_wtaconejetsbasicclusterssk.open(out.gepWTAConeJetsBasicClustersSK);
        f_topo.open(out.topo422);
        f_gfex.open(out.gFex);
        f_jfex.open(out.jFex);
        if (!f_topotower.is_open() || !f_topo.is_open() || !f_gfex.is_open() || !f_jfex.is_open() || !f_gepbasicclusters.is_open() || !f_gepcellstowers.is_open() || !f_gepbasicclusterssk.is_open() || !f_gepcellstowerssk.is_open() || !f_wtaconejetscellstowers.is_open() || !f_wtaconejetsbasicclusters.is_open() || !f_wtaconejetscellstowerssk.is_open() || !f_wtaconejetsbasicclusterssk.is_open()) {
            std::cerr << "Error: One or more .dat output files could not be opened for writing!" << std::endl;
        }
    }

    TString outputFileName = makeOutputFileName(signalString, signalBool, jzSlice, algoVersion, specialJZ0Bool, outputDirOverride, fileSuffix);
    
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
    TTree* outOfTimeAntiKt4TruthJetsTree = new TTree("outOfTimeAntiKt4TruthJetsTree", "Tree storing event-wise Et, eta phi");
    TTree* leadingOutOfTimeAntiKt4TruthJetsTree = new TTree("leadingOutOfTimeAntiKt4TruthJetsTree", "Tree storing event-wise Et, eta phi");
    TTree* subleadingOutOfTimeAntiKt4TruthJetsTree = new TTree("subleadingOutOfTimeAntiKt4TruthJetsTree", "Tree storing event-wise Et, eta phi");
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
    // Sim gFEX/jFEX trees (from TrigGepPerf ntuple)
    TTree* jFexSRJSimTree = new TTree("jFexSRJSimTree", "Tree storing event-wise Et, Eta, Phi (sim jFEX SR jets)");
    TTree* jFexLeadingSRJSimTree = new TTree("jFexLeadingSRJSimTree", "Tree storing event-wise Et, Eta, Phi (sim jFEX SR jets)");
    TTree* jFexSubleadingSRJSimTree = new TTree("jFexSubleadingSRJSimTree", "Tree storing event-wise Et, Eta, Phi (sim jFEX SR jets)");
    TTree* gFexLRJSimTree = new TTree("gFexLRJSimTree", "Tree storing event-wise Et, Eta, Phi (sim gFEX LR jets)");
    TTree* gFexLeadingLRJSimTree = new TTree("gFexLeadingLRJSimTree", "Tree storing event-wise Et, Eta, Phi (sim gFEX LR jets)");
    TTree* gFexSubleadingLRJSimTree = new TTree("gFexSubleadingLRJSimTree", "Tree storing event-wise Et, Eta, Phi (sim gFEX LR jets)");
    // Sim gFEX/jFEX MET trees (from TrigGepPerf ntuple)
    TTree* gFexMETJwoJSimTree      = new TTree("gFexMETJwoJSimTree",      "Tree storing event-wise sim gFEX MET quantities (JwoJ)");
    TTree* gFexMETNoiseCutSimTree  = new TTree("gFexMETNoiseCutSimTree",  "Tree storing event-wise sim gFEX MET quantities (NoiseCut)");
    TTree* gFexMETRmsSimTree       = new TTree("gFexMETRmsSimTree",       "Tree storing event-wise sim gFEX MET quantities (Rms)");
    TTree* jFexMETTree             = new TTree("jFexMETTree",             "Tree storing event-wise jFEX MET quantities");
    // TODO: uncomment once L1_gFexSRJetRoISim branch exists in GEP ntuples
    // TTree* gFexSRJSimTree = new TTree("gFexSRJSimTree", "Tree storing event-wise Et, Eta, Phi (sim gFEX SR jets)");
    // TTree* gFexLeadingSRJSimTree = new TTree("gFexLeadingSRJSimTree", "Tree storing event-wise Et, Eta, Phi (sim gFEX SR jets)");
    // TTree* gFexSubleadingSRJSimTree = new TTree("gFexSubleadingSRJSimTree", "Tree storing event-wise Et, Eta, Phi (sim gFEX SR jets)");
    // EtaSK tower/cluster trees
    TTree* gepCellsTowersEtaSKTree = new TTree("gepCellsTowersEtaSKTree", "Tree storing event-wise Et, Eta, Phi");
    TTree* gepBasicClustersEtaSKTree = new TTree("gepBasicClustersEtaSKTree", "Tree storing event-wise Et, Eta, Phi");
    // EtaSK WTACone jet trees
    TTree* gepWTAConeCellsTowersEtaSKJetsTree = new TTree("gepWTAConeCellsTowersEtaSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepLeadingWTAConeCellsTowersEtaSKJetsTree = new TTree("gepLeadingWTAConeCellsTowersEtaSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepSubleadingWTAConeCellsTowersEtaSKJetsTree = new TTree("gepSubleadingWTAConeCellsTowersEtaSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepWTAConeBasicClustersEtaSKJetsTree = new TTree("gepWTAConeBasicClustersEtaSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepLeadingWTAConeBasicClustersEtaSKJetsTree = new TTree("gepLeadingWTAConeBasicClustersEtaSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    TTree* gepSubleadingWTAConeBasicClustersEtaSKJetsTree = new TTree("gepSubleadingWTAConeBasicClustersEtaSKJetsTree", "Tree storing event-wise Pt, Eta, Phi, Mass, NConstituents");
    /*TTree* recoAntiKt10LCTopoJets = new TTree("recoAntiKt10LCTopoJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* leadingRecoAntiKt10LCTopoJets = new TTree("leadingRecoAntiKt10LCTopoJets", "Tree storing event-wise Et, Eta, Phi");
    TTree* subleadingRecoAntiKt10LCTopoJets = new TTree("recoAntiKt10LCTopoJets", "Tree storing event-wise Et, Eta, Phi");*/
    //TTree* gFexRhoRoI = new TTree("gFexRhoRoI", "Tree storing event-wise Et, Eta, Phi");

    // Variables to store data to be written to TTrees

    // Event info vectors
    std::vector<double> mcEventWeights;
    double sumOfWeightsForSample;
    int sampleJZSlice;
    int gepRunNumberIn = 0, gepEventNumberIn = 0;
    int gepRunNumberOut = 0, gepEventNumberOut = 0;
    std::vector<double> eventWeights;
    bool passHSTP = true;

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
    std::vector<double> gepCellsTowersEtaSKEtValues, gepCellsTowersEtaSKEtaValues, gepCellsTowersEtaSKPhiValues;
    std::vector<double> gepBasicClustersEtaSKEtValues, gepBasicClustersEtaSKEtaValues, gepBasicClustersEtaSKPhiValues;
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

    // EtaSK WTACone jets from TrigGepPerf
    std::vector<double> gepWTAConeCellsTowersEtaSKJetspTValues, gepWTAConeCellsTowersEtaSKJetsEtaValues, gepWTAConeCellsTowersEtaSKJetsPhiValues;
    std::vector<unsigned int> gepWTAConeCellsTowersEtaSKJetsNConstituentsValues;
    std::vector<double> gepWTAConeBasicClustersEtaSKJetspTValues, gepWTAConeBasicClustersEtaSKJetsEtaValues, gepWTAConeBasicClustersEtaSKJetsPhiValues;
    std::vector<unsigned int> gepWTAConeBasicClustersEtaSKJetsNConstituentsValues;
    std::vector<double> gepLeadingWTAConeCellsTowersEtaSKJetspTValues, gepLeadingWTAConeCellsTowersEtaSKJetsEtaValues, gepLeadingWTAConeCellsTowersEtaSKJetsPhiValues;
    std::vector<unsigned int> gepLeadingWTAConeCellsTowersEtaSKJetsNConstituentsValues;
    std::vector<double> gepLeadingWTAConeBasicClustersEtaSKJetspTValues, gepLeadingWTAConeBasicClustersEtaSKJetsEtaValues, gepLeadingWTAConeBasicClustersEtaSKJetsPhiValues;
    std::vector<unsigned int> gepLeadingWTAConeBasicClustersEtaSKJetsNConstituentsValues;
    std::vector<double> gepSubleadingWTAConeCellsTowersEtaSKJetspTValues, gepSubleadingWTAConeCellsTowersEtaSKJetsEtaValues, gepSubleadingWTAConeCellsTowersEtaSKJetsPhiValues;
    std::vector<unsigned int> gepSubleadingWTAConeCellsTowersEtaSKJetsNConstituentsValues;
    std::vector<double> gepSubleadingWTAConeBasicClustersEtaSKJetspTValues, gepSubleadingWTAConeBasicClustersEtaSKJetsEtaValues, gepSubleadingWTAConeBasicClustersEtaSKJetsPhiValues;
    std::vector<unsigned int> gepSubleadingWTAConeBasicClustersEtaSKJetsNConstituentsValues;

    // Ring Et and TobN for WTACone CellsTower jets (all/leading/subleading)
    std::vector<float> gepWTAConeCellsTowersJetsRing0Et, gepWTAConeCellsTowersJetsRing1Et, gepWTAConeCellsTowersJetsRing2Et, gepWTAConeCellsTowersJetsRing3Et, gepWTAConeCellsTowersJetsRing4Et;
    std::vector<int>   gepWTAConeCellsTowersJetsTotalTobN, gepWTAConeCellsTowersJetsRing0TobN, gepWTAConeCellsTowersJetsRing1TobN, gepWTAConeCellsTowersJetsRing2TobN, gepWTAConeCellsTowersJetsRing3TobN, gepWTAConeCellsTowersJetsRing4TobN;
    std::vector<float> gepLeadingWTAConeCellsTowersJetsRing0Et, gepLeadingWTAConeCellsTowersJetsRing1Et, gepLeadingWTAConeCellsTowersJetsRing2Et, gepLeadingWTAConeCellsTowersJetsRing3Et, gepLeadingWTAConeCellsTowersJetsRing4Et;
    std::vector<int>   gepLeadingWTAConeCellsTowersJetsTotalTobN, gepLeadingWTAConeCellsTowersJetsRing0TobN, gepLeadingWTAConeCellsTowersJetsRing1TobN, gepLeadingWTAConeCellsTowersJetsRing2TobN, gepLeadingWTAConeCellsTowersJetsRing3TobN, gepLeadingWTAConeCellsTowersJetsRing4TobN;
    std::vector<float> gepSubleadingWTAConeCellsTowersJetsRing0Et, gepSubleadingWTAConeCellsTowersJetsRing1Et, gepSubleadingWTAConeCellsTowersJetsRing2Et, gepSubleadingWTAConeCellsTowersJetsRing3Et, gepSubleadingWTAConeCellsTowersJetsRing4Et;
    std::vector<int>   gepSubleadingWTAConeCellsTowersJetsTotalTobN, gepSubleadingWTAConeCellsTowersJetsRing0TobN, gepSubleadingWTAConeCellsTowersJetsRing1TobN, gepSubleadingWTAConeCellsTowersJetsRing2TobN, gepSubleadingWTAConeCellsTowersJetsRing3TobN, gepSubleadingWTAConeCellsTowersJetsRing4TobN;
    // Ring Et and TobN for WTACone CellsTower SK jets
    std::vector<float> gepWTAConeCellsTowersSKJetsRing0Et, gepWTAConeCellsTowersSKJetsRing1Et, gepWTAConeCellsTowersSKJetsRing2Et, gepWTAConeCellsTowersSKJetsRing3Et, gepWTAConeCellsTowersSKJetsRing4Et;
    std::vector<int>   gepWTAConeCellsTowersSKJetsTotalTobN, gepWTAConeCellsTowersSKJetsRing0TobN, gepWTAConeCellsTowersSKJetsRing1TobN, gepWTAConeCellsTowersSKJetsRing2TobN, gepWTAConeCellsTowersSKJetsRing3TobN, gepWTAConeCellsTowersSKJetsRing4TobN;
    std::vector<float> gepLeadingWTAConeCellsTowersSKJetsRing0Et, gepLeadingWTAConeCellsTowersSKJetsRing1Et, gepLeadingWTAConeCellsTowersSKJetsRing2Et, gepLeadingWTAConeCellsTowersSKJetsRing3Et, gepLeadingWTAConeCellsTowersSKJetsRing4Et;
    std::vector<int>   gepLeadingWTAConeCellsTowersSKJetsTotalTobN, gepLeadingWTAConeCellsTowersSKJetsRing0TobN, gepLeadingWTAConeCellsTowersSKJetsRing1TobN, gepLeadingWTAConeCellsTowersSKJetsRing2TobN, gepLeadingWTAConeCellsTowersSKJetsRing3TobN, gepLeadingWTAConeCellsTowersSKJetsRing4TobN;
    std::vector<float> gepSubleadingWTAConeCellsTowersSKJetsRing0Et, gepSubleadingWTAConeCellsTowersSKJetsRing1Et, gepSubleadingWTAConeCellsTowersSKJetsRing2Et, gepSubleadingWTAConeCellsTowersSKJetsRing3Et, gepSubleadingWTAConeCellsTowersSKJetsRing4Et;
    std::vector<int>   gepSubleadingWTAConeCellsTowersSKJetsTotalTobN, gepSubleadingWTAConeCellsTowersSKJetsRing0TobN, gepSubleadingWTAConeCellsTowersSKJetsRing1TobN, gepSubleadingWTAConeCellsTowersSKJetsRing2TobN, gepSubleadingWTAConeCellsTowersSKJetsRing3TobN, gepSubleadingWTAConeCellsTowersSKJetsRing4TobN;
    // Ring Et and TobN for WTACone CellsTower EtaSK jets
    std::vector<float> gepWTAConeCellsTowersEtaSKJetsRing0Et, gepWTAConeCellsTowersEtaSKJetsRing1Et, gepWTAConeCellsTowersEtaSKJetsRing2Et, gepWTAConeCellsTowersEtaSKJetsRing3Et, gepWTAConeCellsTowersEtaSKJetsRing4Et;
    std::vector<int>   gepWTAConeCellsTowersEtaSKJetsTotalTobN, gepWTAConeCellsTowersEtaSKJetsRing0TobN, gepWTAConeCellsTowersEtaSKJetsRing1TobN, gepWTAConeCellsTowersEtaSKJetsRing2TobN, gepWTAConeCellsTowersEtaSKJetsRing3TobN, gepWTAConeCellsTowersEtaSKJetsRing4TobN;
    std::vector<float> gepLeadingWTAConeCellsTowersEtaSKJetsRing0Et, gepLeadingWTAConeCellsTowersEtaSKJetsRing1Et, gepLeadingWTAConeCellsTowersEtaSKJetsRing2Et, gepLeadingWTAConeCellsTowersEtaSKJetsRing3Et, gepLeadingWTAConeCellsTowersEtaSKJetsRing4Et;
    std::vector<int>   gepLeadingWTAConeCellsTowersEtaSKJetsTotalTobN, gepLeadingWTAConeCellsTowersEtaSKJetsRing0TobN, gepLeadingWTAConeCellsTowersEtaSKJetsRing1TobN, gepLeadingWTAConeCellsTowersEtaSKJetsRing2TobN, gepLeadingWTAConeCellsTowersEtaSKJetsRing3TobN, gepLeadingWTAConeCellsTowersEtaSKJetsRing4TobN;
    std::vector<float> gepSubleadingWTAConeCellsTowersEtaSKJetsRing0Et, gepSubleadingWTAConeCellsTowersEtaSKJetsRing1Et, gepSubleadingWTAConeCellsTowersEtaSKJetsRing2Et, gepSubleadingWTAConeCellsTowersEtaSKJetsRing3Et, gepSubleadingWTAConeCellsTowersEtaSKJetsRing4Et;
    std::vector<int>   gepSubleadingWTAConeCellsTowersEtaSKJetsTotalTobN, gepSubleadingWTAConeCellsTowersEtaSKJetsRing0TobN, gepSubleadingWTAConeCellsTowersEtaSKJetsRing1TobN, gepSubleadingWTAConeCellsTowersEtaSKJetsRing2TobN, gepSubleadingWTAConeCellsTowersEtaSKJetsRing3TobN, gepSubleadingWTAConeCellsTowersEtaSKJetsRing4TobN;
    // Ring Et and TobN for WTACone BasicClusters jets
    std::vector<float> gepWTAConeBasicClustersJetsRing0Et, gepWTAConeBasicClustersJetsRing1Et, gepWTAConeBasicClustersJetsRing2Et, gepWTAConeBasicClustersJetsRing3Et, gepWTAConeBasicClustersJetsRing4Et;
    std::vector<int>   gepWTAConeBasicClustersJetsTotalTobN, gepWTAConeBasicClustersJetsRing0TobN, gepWTAConeBasicClustersJetsRing1TobN, gepWTAConeBasicClustersJetsRing2TobN, gepWTAConeBasicClustersJetsRing3TobN, gepWTAConeBasicClustersJetsRing4TobN;
    std::vector<float> gepLeadingWTAConeBasicClustersJetsRing0Et, gepLeadingWTAConeBasicClustersJetsRing1Et, gepLeadingWTAConeBasicClustersJetsRing2Et, gepLeadingWTAConeBasicClustersJetsRing3Et, gepLeadingWTAConeBasicClustersJetsRing4Et;
    std::vector<int>   gepLeadingWTAConeBasicClustersJetsTotalTobN, gepLeadingWTAConeBasicClustersJetsRing0TobN, gepLeadingWTAConeBasicClustersJetsRing1TobN, gepLeadingWTAConeBasicClustersJetsRing2TobN, gepLeadingWTAConeBasicClustersJetsRing3TobN, gepLeadingWTAConeBasicClustersJetsRing4TobN;
    std::vector<float> gepSubleadingWTAConeBasicClustersJetsRing0Et, gepSubleadingWTAConeBasicClustersJetsRing1Et, gepSubleadingWTAConeBasicClustersJetsRing2Et, gepSubleadingWTAConeBasicClustersJetsRing3Et, gepSubleadingWTAConeBasicClustersJetsRing4Et;
    std::vector<int>   gepSubleadingWTAConeBasicClustersJetsTotalTobN, gepSubleadingWTAConeBasicClustersJetsRing0TobN, gepSubleadingWTAConeBasicClustersJetsRing1TobN, gepSubleadingWTAConeBasicClustersJetsRing2TobN, gepSubleadingWTAConeBasicClustersJetsRing3TobN, gepSubleadingWTAConeBasicClustersJetsRing4TobN;
    // Ring Et and TobN for WTACone BasicClusters SK jets
    std::vector<float> gepWTAConeBasicClustersSKJetsRing0Et, gepWTAConeBasicClustersSKJetsRing1Et, gepWTAConeBasicClustersSKJetsRing2Et, gepWTAConeBasicClustersSKJetsRing3Et, gepWTAConeBasicClustersSKJetsRing4Et;
    std::vector<int>   gepWTAConeBasicClustersSKJetsTotalTobN, gepWTAConeBasicClustersSKJetsRing0TobN, gepWTAConeBasicClustersSKJetsRing1TobN, gepWTAConeBasicClustersSKJetsRing2TobN, gepWTAConeBasicClustersSKJetsRing3TobN, gepWTAConeBasicClustersSKJetsRing4TobN;
    std::vector<float> gepLeadingWTAConeBasicClustersSKJetsRing0Et, gepLeadingWTAConeBasicClustersSKJetsRing1Et, gepLeadingWTAConeBasicClustersSKJetsRing2Et, gepLeadingWTAConeBasicClustersSKJetsRing3Et, gepLeadingWTAConeBasicClustersSKJetsRing4Et;
    std::vector<int>   gepLeadingWTAConeBasicClustersSKJetsTotalTobN, gepLeadingWTAConeBasicClustersSKJetsRing0TobN, gepLeadingWTAConeBasicClustersSKJetsRing1TobN, gepLeadingWTAConeBasicClustersSKJetsRing2TobN, gepLeadingWTAConeBasicClustersSKJetsRing3TobN, gepLeadingWTAConeBasicClustersSKJetsRing4TobN;
    std::vector<float> gepSubleadingWTAConeBasicClustersSKJetsRing0Et, gepSubleadingWTAConeBasicClustersSKJetsRing1Et, gepSubleadingWTAConeBasicClustersSKJetsRing2Et, gepSubleadingWTAConeBasicClustersSKJetsRing3Et, gepSubleadingWTAConeBasicClustersSKJetsRing4Et;
    std::vector<int>   gepSubleadingWTAConeBasicClustersSKJetsTotalTobN, gepSubleadingWTAConeBasicClustersSKJetsRing0TobN, gepSubleadingWTAConeBasicClustersSKJetsRing1TobN, gepSubleadingWTAConeBasicClustersSKJetsRing2TobN, gepSubleadingWTAConeBasicClustersSKJetsRing3TobN, gepSubleadingWTAConeBasicClustersSKJetsRing4TobN;
    // Ring Et and TobN for WTACone BasicClusters EtaSK jets
    std::vector<float> gepWTAConeBasicClustersEtaSKJetsRing0Et, gepWTAConeBasicClustersEtaSKJetsRing1Et, gepWTAConeBasicClustersEtaSKJetsRing2Et, gepWTAConeBasicClustersEtaSKJetsRing3Et, gepWTAConeBasicClustersEtaSKJetsRing4Et;
    std::vector<int>   gepWTAConeBasicClustersEtaSKJetsTotalTobN, gepWTAConeBasicClustersEtaSKJetsRing0TobN, gepWTAConeBasicClustersEtaSKJetsRing1TobN, gepWTAConeBasicClustersEtaSKJetsRing2TobN, gepWTAConeBasicClustersEtaSKJetsRing3TobN, gepWTAConeBasicClustersEtaSKJetsRing4TobN;
    std::vector<float> gepLeadingWTAConeBasicClustersEtaSKJetsRing0Et, gepLeadingWTAConeBasicClustersEtaSKJetsRing1Et, gepLeadingWTAConeBasicClustersEtaSKJetsRing2Et, gepLeadingWTAConeBasicClustersEtaSKJetsRing3Et, gepLeadingWTAConeBasicClustersEtaSKJetsRing4Et;
    std::vector<int>   gepLeadingWTAConeBasicClustersEtaSKJetsTotalTobN, gepLeadingWTAConeBasicClustersEtaSKJetsRing0TobN, gepLeadingWTAConeBasicClustersEtaSKJetsRing1TobN, gepLeadingWTAConeBasicClustersEtaSKJetsRing2TobN, gepLeadingWTAConeBasicClustersEtaSKJetsRing3TobN, gepLeadingWTAConeBasicClustersEtaSKJetsRing4TobN;
    std::vector<float> gepSubleadingWTAConeBasicClustersEtaSKJetsRing0Et, gepSubleadingWTAConeBasicClustersEtaSKJetsRing1Et, gepSubleadingWTAConeBasicClustersEtaSKJetsRing2Et, gepSubleadingWTAConeBasicClustersEtaSKJetsRing3Et, gepSubleadingWTAConeBasicClustersEtaSKJetsRing4Et;
    std::vector<int>   gepSubleadingWTAConeBasicClustersEtaSKJetsTotalTobN, gepSubleadingWTAConeBasicClustersEtaSKJetsRing0TobN, gepSubleadingWTAConeBasicClustersEtaSKJetsRing1TobN, gepSubleadingWTAConeBasicClustersEtaSKJetsRing2TobN, gepSubleadingWTAConeBasicClustersEtaSKJetsRing3TobN, gepSubleadingWTAConeBasicClustersEtaSKJetsRing4TobN;

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
    // Sim jFEX SR jets (from TrigGepPerf ntuple)
    std::vector<unsigned int> jFexSRJSimEtIndexValues;
    std::vector<double> jFexSRJSimEtValues, jFexSRJSimEtaValues, jFexSRJSimPhiValues;
    std::vector<double> jFexSRJSimLeadingEtValues, jFexSRJSimLeadingEtaValues, jFexSRJSimLeadingPhiValues;
    std::vector<double> jFexSRJSimSubleadingEtValues, jFexSRJSimSubleadingEtaValues, jFexSRJSimSubleadingPhiValues;
    // Sim gFEX LR jets (from TrigGepPerf ntuple)
    std::vector<unsigned int> gFexLRJSimEtIndexValues;
    std::vector<double> gFexLRJSimEtValues, gFexLRJSimEtaValues, gFexLRJSimPhiValues;
    std::vector<double> gFexLRJSimLeadingEtValues, gFexLRJSimLeadingEtaValues, gFexLRJSimLeadingPhiValues;
    std::vector<double> gFexLRJSimSubleadingEtValues, gFexLRJSimSubleadingEtaValues, gFexLRJSimSubleadingPhiValues;
    // TODO: uncomment once L1_gFexSRJetRoISim branch exists in GEP ntuples
    // // Sim gFEX SR jets (from TrigGepPerf ntuple)
    // std::vector<unsigned int> gFexSRJSimEtIndexValues;
    // std::vector<double> gFexSRJSimEtValues, gFexSRJSimEtaValues, gFexSRJSimPhiValues;
    // std::vector<double> gFexSRJSimLeadingEtValues, gFexSRJSimLeadingEtaValues, gFexSRJSimLeadingPhiValues;
    // std::vector<double> gFexSRJSimSubleadingEtValues, gFexSRJSimSubleadingEtaValues, gFexSRJSimSubleadingPhiValues;
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
    double metTruthIntMuonsX = 0.0, metTruthIntMuonsY = 0.0;   // MET_Truth["IntMuons"]
    double metTruthIntMuonsSumET = 0.0, metTruthIntMuons = 0.0;
    // MET_Core_AntiKt4EMTopo — one set of variables per term
    double metCoreEMTopo_SoftClus_X = 0.0,   metCoreEMTopo_SoftClus_Y = 0.0,   metCoreEMTopo_SoftClus_SumET = 0.0,   metCoreEMTopo_SoftClus_MET = 0.0;
    double metCoreEMTopo_PVSoftTrk_X = 0.0,  metCoreEMTopo_PVSoftTrk_Y = 0.0,  metCoreEMTopo_PVSoftTrk_SumET = 0.0,  metCoreEMTopo_PVSoftTrk_MET = 0.0;
    double metCoreEMTopo_SoftClusEM_X = 0.0, metCoreEMTopo_SoftClusEM_Y = 0.0, metCoreEMTopo_SoftClusEM_SumET = 0.0, metCoreEMTopo_SoftClusEM_MET = 0.0;
    // MET_Core_AntiKt4EMPFlow — one set of variables per term
    double metCoreEMPFlow_SoftClus_X = 0.0,  metCoreEMPFlow_SoftClus_Y = 0.0,  metCoreEMPFlow_SoftClus_SumET = 0.0,  metCoreEMPFlow_SoftClus_MET = 0.0;
    double metCoreEMPFlow_PVSoftTrk_X = 0.0, metCoreEMPFlow_PVSoftTrk_Y = 0.0, metCoreEMPFlow_PVSoftTrk_SumET = 0.0, metCoreEMPFlow_PVSoftTrk_MET = 0.0;
    // Sim gFEX/jFEX MET quantities (from TrigGepPerf ntuple)
    double gFexMETJwoJSim = 0.0, gFexMETPhiJwoJSim = 0.0, gFexSumETJwoJSim = 0.0;
    double gFexMETNoiseCutSim = 0.0, gFexMETPhiNoiseCutSim = 0.0, gFexSumETNoiseCutSim = 0.0;
    double gFexMETRmsSim = 0.0, gFexMETPhiRmsSim = 0.0, gFexSumETRmsSim = 0.0;
    double jFexMET = 0.0;

    // In time anti-kt 4 truth jets vectors
    std::vector<unsigned int> inTimeAntiKt4TruthSRJEtIndexValues;
    std::vector<double> inTimeAntiKt4TruthSRJEtValues, inTimeAntiKt4TruthSRJEtaValues, inTimeAntiKt4TruthSRJPhiValues;
    std::vector<double> inTimeAntiKt4TruthSRJLeadingEtValues, inTimeAntiKt4TruthSRJLeadingEtaValues, inTimeAntiKt4TruthSRJLeadingPhiValues;
    std::vector<double> inTimeAntiKt4TruthSRJSubleadingEtValues, inTimeAntiKt4TruthSRJSubleadingEtaValues, inTimeAntiKt4TruthSRJSubleadingPhiValues;

    // Out of time anti-kt 4 truth jets vectors
    std::vector<unsigned int> outOfTimeAntiKt4TruthSRJEtIndexValues;
    std::vector<double> outOfTimeAntiKt4TruthSRJEtValues, outOfTimeAntiKt4TruthSRJEtaValues, outOfTimeAntiKt4TruthSRJPhiValues;
    std::vector<double> outOfTimeAntiKt4TruthSRJLeadingEtValues, outOfTimeAntiKt4TruthSRJLeadingEtaValues, outOfTimeAntiKt4TruthSRJLeadingPhiValues;
    std::vector<double> outOfTimeAntiKt4TruthSRJSubleadingEtValues, outOfTimeAntiKt4TruthSRJSubleadingEtaValues, outOfTimeAntiKt4TruthSRJSubleadingPhiValues;

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
    eventInfoTree->Branch("gepRunNumberOut", &gepRunNumberOut);
    eventInfoTree->Branch("gepEventNumberOut", &gepEventNumberOut);
    eventInfoTree->Branch("passHSTP", &passHSTP);

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

    // gepCellsTowersEtaSKTree
    gepCellsTowersEtaSKTree->Branch("Et", &gepCellsTowersEtaSKEtValues);
    gepCellsTowersEtaSKTree->Branch("Eta", &gepCellsTowersEtaSKEtaValues);
    gepCellsTowersEtaSKTree->Branch("Phi", &gepCellsTowersEtaSKPhiValues);

    // gepBasicClustersEtaSKTree
    gepBasicClustersEtaSKTree->Branch("Et", &gepBasicClustersEtaSKEtValues);
    gepBasicClustersEtaSKTree->Branch("Eta", &gepBasicClustersEtaSKEtaValues);
    gepBasicClustersEtaSKTree->Branch("Phi", &gepBasicClustersEtaSKPhiValues);

    // gep wta cone cells towers jets
    gepWTAConeCellsTowersJetsTree->Branch("Et", &gepWTAConeCellsTowersJetspTValues);
    gepWTAConeCellsTowersJetsTree->Branch("Eta", &gepWTAConeCellsTowersJetsEtaValues);
    gepWTAConeCellsTowersJetsTree->Branch("Phi", &gepWTAConeCellsTowersJetsPhiValues);
    gepWTAConeCellsTowersJetsTree->Branch("NConstituents", &gepWTAConeCellsTowersJetsNConstituentsValues);
    gepWTAConeCellsTowersJetsTree->Branch("Ring0Et", &gepWTAConeCellsTowersJetsRing0Et);
    gepWTAConeCellsTowersJetsTree->Branch("Ring1Et", &gepWTAConeCellsTowersJetsRing1Et);
    gepWTAConeCellsTowersJetsTree->Branch("Ring2Et", &gepWTAConeCellsTowersJetsRing2Et);
    gepWTAConeCellsTowersJetsTree->Branch("Ring3Et", &gepWTAConeCellsTowersJetsRing3Et);
    gepWTAConeCellsTowersJetsTree->Branch("Ring4Et", &gepWTAConeCellsTowersJetsRing4Et);
    gepWTAConeCellsTowersJetsTree->Branch("TotalTobN", &gepWTAConeCellsTowersJetsTotalTobN);
    gepWTAConeCellsTowersJetsTree->Branch("Ring0TobN", &gepWTAConeCellsTowersJetsRing0TobN);
    gepWTAConeCellsTowersJetsTree->Branch("Ring1TobN", &gepWTAConeCellsTowersJetsRing1TobN);
    gepWTAConeCellsTowersJetsTree->Branch("Ring2TobN", &gepWTAConeCellsTowersJetsRing2TobN);
    gepWTAConeCellsTowersJetsTree->Branch("Ring3TobN", &gepWTAConeCellsTowersJetsRing3TobN);
    gepWTAConeCellsTowersJetsTree->Branch("Ring4TobN", &gepWTAConeCellsTowersJetsRing4TobN);

    // gep wta cone cells towers leading jets
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Et", &gepLeadingWTAConeCellsTowersJetspTValues);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Eta", &gepLeadingWTAConeCellsTowersJetsEtaValues);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Phi", &gepLeadingWTAConeCellsTowersJetsPhiValues);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("NConstituents", &gepLeadingWTAConeCellsTowersJetsNConstituentsValues);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Ring0Et", &gepLeadingWTAConeCellsTowersJetsRing0Et);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Ring1Et", &gepLeadingWTAConeCellsTowersJetsRing1Et);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Ring2Et", &gepLeadingWTAConeCellsTowersJetsRing2Et);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Ring3Et", &gepLeadingWTAConeCellsTowersJetsRing3Et);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Ring4Et", &gepLeadingWTAConeCellsTowersJetsRing4Et);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("TotalTobN", &gepLeadingWTAConeCellsTowersJetsTotalTobN);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Ring0TobN", &gepLeadingWTAConeCellsTowersJetsRing0TobN);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Ring1TobN", &gepLeadingWTAConeCellsTowersJetsRing1TobN);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Ring2TobN", &gepLeadingWTAConeCellsTowersJetsRing2TobN);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Ring3TobN", &gepLeadingWTAConeCellsTowersJetsRing3TobN);
    gepLeadingWTAConeCellsTowersJetsTree->Branch("Ring4TobN", &gepLeadingWTAConeCellsTowersJetsRing4TobN);

    // gep wta cone cells towers subleading jets
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Et", &gepSubleadingWTAConeCellsTowersJetspTValues);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Eta", &gepSubleadingWTAConeCellsTowersJetsEtaValues);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Phi", &gepSubleadingWTAConeCellsTowersJetsPhiValues);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("NConstituents", &gepSubleadingWTAConeCellsTowersJetsNConstituentsValues);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Ring0Et", &gepSubleadingWTAConeCellsTowersJetsRing0Et);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Ring1Et", &gepSubleadingWTAConeCellsTowersJetsRing1Et);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Ring2Et", &gepSubleadingWTAConeCellsTowersJetsRing2Et);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Ring3Et", &gepSubleadingWTAConeCellsTowersJetsRing3Et);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Ring4Et", &gepSubleadingWTAConeCellsTowersJetsRing4Et);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("TotalTobN", &gepSubleadingWTAConeCellsTowersJetsTotalTobN);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Ring0TobN", &gepSubleadingWTAConeCellsTowersJetsRing0TobN);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Ring1TobN", &gepSubleadingWTAConeCellsTowersJetsRing1TobN);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Ring2TobN", &gepSubleadingWTAConeCellsTowersJetsRing2TobN);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Ring3TobN", &gepSubleadingWTAConeCellsTowersJetsRing3TobN);
    gepSubleadingWTAConeCellsTowersJetsTree->Branch("Ring4TobN", &gepSubleadingWTAConeCellsTowersJetsRing4TobN);

     // gep wta cone cells towers SK jets
    gepWTAConeCellsTowersSKJetsTree->Branch("Et", &gepWTAConeCellsTowersSKJetspTValues);
    gepWTAConeCellsTowersSKJetsTree->Branch("Eta", &gepWTAConeCellsTowersSKJetsEtaValues);
    gepWTAConeCellsTowersSKJetsTree->Branch("Phi", &gepWTAConeCellsTowersSKJetsPhiValues);
    gepWTAConeCellsTowersSKJetsTree->Branch("NConstituents", &gepWTAConeCellsTowersSKJetsNConstituentsValues);
    gepWTAConeCellsTowersSKJetsTree->Branch("Ring0Et", &gepWTAConeCellsTowersSKJetsRing0Et);
    gepWTAConeCellsTowersSKJetsTree->Branch("Ring1Et", &gepWTAConeCellsTowersSKJetsRing1Et);
    gepWTAConeCellsTowersSKJetsTree->Branch("Ring2Et", &gepWTAConeCellsTowersSKJetsRing2Et);
    gepWTAConeCellsTowersSKJetsTree->Branch("Ring3Et", &gepWTAConeCellsTowersSKJetsRing3Et);
    gepWTAConeCellsTowersSKJetsTree->Branch("Ring4Et", &gepWTAConeCellsTowersSKJetsRing4Et);
    gepWTAConeCellsTowersSKJetsTree->Branch("TotalTobN", &gepWTAConeCellsTowersSKJetsTotalTobN);
    gepWTAConeCellsTowersSKJetsTree->Branch("Ring0TobN", &gepWTAConeCellsTowersSKJetsRing0TobN);
    gepWTAConeCellsTowersSKJetsTree->Branch("Ring1TobN", &gepWTAConeCellsTowersSKJetsRing1TobN);
    gepWTAConeCellsTowersSKJetsTree->Branch("Ring2TobN", &gepWTAConeCellsTowersSKJetsRing2TobN);
    gepWTAConeCellsTowersSKJetsTree->Branch("Ring3TobN", &gepWTAConeCellsTowersSKJetsRing3TobN);
    gepWTAConeCellsTowersSKJetsTree->Branch("Ring4TobN", &gepWTAConeCellsTowersSKJetsRing4TobN);

    // gep wta cone cells towers SK leading jets
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Et", &gepLeadingWTAConeCellsTowersSKJetspTValues);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Eta", &gepLeadingWTAConeCellsTowersSKJetsEtaValues);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Phi", &gepLeadingWTAConeCellsTowersSKJetsPhiValues);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("NConstituents", &gepLeadingWTAConeCellsTowersSKJetsNConstituentsValues);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Ring0Et", &gepLeadingWTAConeCellsTowersSKJetsRing0Et);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Ring1Et", &gepLeadingWTAConeCellsTowersSKJetsRing1Et);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Ring2Et", &gepLeadingWTAConeCellsTowersSKJetsRing2Et);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Ring3Et", &gepLeadingWTAConeCellsTowersSKJetsRing3Et);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Ring4Et", &gepLeadingWTAConeCellsTowersSKJetsRing4Et);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("TotalTobN", &gepLeadingWTAConeCellsTowersSKJetsTotalTobN);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Ring0TobN", &gepLeadingWTAConeCellsTowersSKJetsRing0TobN);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Ring1TobN", &gepLeadingWTAConeCellsTowersSKJetsRing1TobN);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Ring2TobN", &gepLeadingWTAConeCellsTowersSKJetsRing2TobN);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Ring3TobN", &gepLeadingWTAConeCellsTowersSKJetsRing3TobN);
    gepLeadingWTAConeCellsTowersSKJetsTree->Branch("Ring4TobN", &gepLeadingWTAConeCellsTowersSKJetsRing4TobN);

    // gep wta cone cells towers SK subleading jets
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Et", &gepSubleadingWTAConeCellsTowersSKJetspTValues);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Eta", &gepSubleadingWTAConeCellsTowersSKJetsEtaValues);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Phi", &gepSubleadingWTAConeCellsTowersSKJetsPhiValues);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("NConstituents", &gepSubleadingWTAConeCellsTowersSKJetsNConstituentsValues);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Ring0Et", &gepSubleadingWTAConeCellsTowersSKJetsRing0Et);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Ring1Et", &gepSubleadingWTAConeCellsTowersSKJetsRing1Et);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Ring2Et", &gepSubleadingWTAConeCellsTowersSKJetsRing2Et);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Ring3Et", &gepSubleadingWTAConeCellsTowersSKJetsRing3Et);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Ring4Et", &gepSubleadingWTAConeCellsTowersSKJetsRing4Et);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("TotalTobN", &gepSubleadingWTAConeCellsTowersSKJetsTotalTobN);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Ring0TobN", &gepSubleadingWTAConeCellsTowersSKJetsRing0TobN);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Ring1TobN", &gepSubleadingWTAConeCellsTowersSKJetsRing1TobN);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Ring2TobN", &gepSubleadingWTAConeCellsTowersSKJetsRing2TobN);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Ring3TobN", &gepSubleadingWTAConeCellsTowersSKJetsRing3TobN);
    gepSubleadingWTAConeCellsTowersSKJetsTree->Branch("Ring4TobN", &gepSubleadingWTAConeCellsTowersSKJetsRing4TobN);

    // gep wta cone basic clusters jets
    gepWTAConeBasicClustersJetsTree->Branch("Et", &gepWTAConeGEPBasicClustersJetspTValues);
    gepWTAConeBasicClustersJetsTree->Branch("Eta", &gepWTAConeGEPBasicClustersJetsEtaValues);
    gepWTAConeBasicClustersJetsTree->Branch("Phi", &gepWTAConeGEPBasicClustersJetsPhiValues);
    gepWTAConeBasicClustersJetsTree->Branch("NConstituents", &gepWTAConeGEPBasicClustersJetsNConstituentsValues);
    gepWTAConeBasicClustersJetsTree->Branch("Ring0Et", &gepWTAConeBasicClustersJetsRing0Et);
    gepWTAConeBasicClustersJetsTree->Branch("Ring1Et", &gepWTAConeBasicClustersJetsRing1Et);
    gepWTAConeBasicClustersJetsTree->Branch("Ring2Et", &gepWTAConeBasicClustersJetsRing2Et);
    gepWTAConeBasicClustersJetsTree->Branch("Ring3Et", &gepWTAConeBasicClustersJetsRing3Et);
    gepWTAConeBasicClustersJetsTree->Branch("Ring4Et", &gepWTAConeBasicClustersJetsRing4Et);
    gepWTAConeBasicClustersJetsTree->Branch("TotalTobN", &gepWTAConeBasicClustersJetsTotalTobN);
    gepWTAConeBasicClustersJetsTree->Branch("Ring0TobN", &gepWTAConeBasicClustersJetsRing0TobN);
    gepWTAConeBasicClustersJetsTree->Branch("Ring1TobN", &gepWTAConeBasicClustersJetsRing1TobN);
    gepWTAConeBasicClustersJetsTree->Branch("Ring2TobN", &gepWTAConeBasicClustersJetsRing2TobN);
    gepWTAConeBasicClustersJetsTree->Branch("Ring3TobN", &gepWTAConeBasicClustersJetsRing3TobN);
    gepWTAConeBasicClustersJetsTree->Branch("Ring4TobN", &gepWTAConeBasicClustersJetsRing4TobN);

    // gep wta cone basic clusters leading jets
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Et", &gepLeadingWTAConeGEPBasicClustersJetspTValues);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Eta", &gepLeadingWTAConeGEPBasicClustersJetsEtaValues);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Phi", &gepLeadingWTAConeGEPBasicClustersJetsPhiValues);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("NConstituents", &gepLeadingWTAConeGEPBasicClustersJetsNConstituentsValues);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Ring0Et", &gepLeadingWTAConeBasicClustersJetsRing0Et);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Ring1Et", &gepLeadingWTAConeBasicClustersJetsRing1Et);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Ring2Et", &gepLeadingWTAConeBasicClustersJetsRing2Et);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Ring3Et", &gepLeadingWTAConeBasicClustersJetsRing3Et);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Ring4Et", &gepLeadingWTAConeBasicClustersJetsRing4Et);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("TotalTobN", &gepLeadingWTAConeBasicClustersJetsTotalTobN);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Ring0TobN", &gepLeadingWTAConeBasicClustersJetsRing0TobN);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Ring1TobN", &gepLeadingWTAConeBasicClustersJetsRing1TobN);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Ring2TobN", &gepLeadingWTAConeBasicClustersJetsRing2TobN);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Ring3TobN", &gepLeadingWTAConeBasicClustersJetsRing3TobN);
    gepLeadingWTAConeBasicClustersJetsTree->Branch("Ring4TobN", &gepLeadingWTAConeBasicClustersJetsRing4TobN);

    // gep wta cone basic clusters subleading jets
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Et", &gepSubleadingWTAConeGEPBasicClustersJetspTValues);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Eta", &gepSubleadingWTAConeGEPBasicClustersJetsEtaValues);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Phi", &gepSubleadingWTAConeGEPBasicClustersJetsPhiValues);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("NConstituents", &gepSubleadingWTAConeGEPBasicClustersJetsNConstituentsValues);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Ring0Et", &gepSubleadingWTAConeBasicClustersJetsRing0Et);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Ring1Et", &gepSubleadingWTAConeBasicClustersJetsRing1Et);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Ring2Et", &gepSubleadingWTAConeBasicClustersJetsRing2Et);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Ring3Et", &gepSubleadingWTAConeBasicClustersJetsRing3Et);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Ring4Et", &gepSubleadingWTAConeBasicClustersJetsRing4Et);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("TotalTobN", &gepSubleadingWTAConeBasicClustersJetsTotalTobN);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Ring0TobN", &gepSubleadingWTAConeBasicClustersJetsRing0TobN);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Ring1TobN", &gepSubleadingWTAConeBasicClustersJetsRing1TobN);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Ring2TobN", &gepSubleadingWTAConeBasicClustersJetsRing2TobN);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Ring3TobN", &gepSubleadingWTAConeBasicClustersJetsRing3TobN);
    gepSubleadingWTAConeBasicClustersJetsTree->Branch("Ring4TobN", &gepSubleadingWTAConeBasicClustersJetsRing4TobN);

    // gep wta cone basic clusters SK jets
    gepWTAConeBasicClustersSKJetsTree->Branch("Et", &gepWTAConeGEPBasicClustersSKJetspTValues);
    gepWTAConeBasicClustersSKJetsTree->Branch("Eta", &gepWTAConeGEPBasicClustersSKJetsEtaValues);
    gepWTAConeBasicClustersSKJetsTree->Branch("Phi", &gepWTAConeGEPBasicClustersSKJetsPhiValues);
    gepWTAConeBasicClustersSKJetsTree->Branch("NConstituents", &gepWTAConeGEPBasicClustersSKJetsNConstituentsValues);
    gepWTAConeBasicClustersSKJetsTree->Branch("Ring0Et", &gepWTAConeBasicClustersSKJetsRing0Et);
    gepWTAConeBasicClustersSKJetsTree->Branch("Ring1Et", &gepWTAConeBasicClustersSKJetsRing1Et);
    gepWTAConeBasicClustersSKJetsTree->Branch("Ring2Et", &gepWTAConeBasicClustersSKJetsRing2Et);
    gepWTAConeBasicClustersSKJetsTree->Branch("Ring3Et", &gepWTAConeBasicClustersSKJetsRing3Et);
    gepWTAConeBasicClustersSKJetsTree->Branch("Ring4Et", &gepWTAConeBasicClustersSKJetsRing4Et);
    gepWTAConeBasicClustersSKJetsTree->Branch("TotalTobN", &gepWTAConeBasicClustersSKJetsTotalTobN);
    gepWTAConeBasicClustersSKJetsTree->Branch("Ring0TobN", &gepWTAConeBasicClustersSKJetsRing0TobN);
    gepWTAConeBasicClustersSKJetsTree->Branch("Ring1TobN", &gepWTAConeBasicClustersSKJetsRing1TobN);
    gepWTAConeBasicClustersSKJetsTree->Branch("Ring2TobN", &gepWTAConeBasicClustersSKJetsRing2TobN);
    gepWTAConeBasicClustersSKJetsTree->Branch("Ring3TobN", &gepWTAConeBasicClustersSKJetsRing3TobN);
    gepWTAConeBasicClustersSKJetsTree->Branch("Ring4TobN", &gepWTAConeBasicClustersSKJetsRing4TobN);

    // gep wta cone basic clusters SK leading jets
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Et", &gepLeadingWTAConeGEPBasicClustersSKJetspTValues);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Eta", &gepLeadingWTAConeGEPBasicClustersSKJetsEtaValues);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Phi", &gepLeadingWTAConeGEPBasicClustersSKJetsPhiValues);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("NConstituents", &gepLeadingWTAConeGEPBasicClustersSKJetsNConstituentsValues);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Ring0Et", &gepLeadingWTAConeBasicClustersSKJetsRing0Et);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Ring1Et", &gepLeadingWTAConeBasicClustersSKJetsRing1Et);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Ring2Et", &gepLeadingWTAConeBasicClustersSKJetsRing2Et);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Ring3Et", &gepLeadingWTAConeBasicClustersSKJetsRing3Et);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Ring4Et", &gepLeadingWTAConeBasicClustersSKJetsRing4Et);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("TotalTobN", &gepLeadingWTAConeBasicClustersSKJetsTotalTobN);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Ring0TobN", &gepLeadingWTAConeBasicClustersSKJetsRing0TobN);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Ring1TobN", &gepLeadingWTAConeBasicClustersSKJetsRing1TobN);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Ring2TobN", &gepLeadingWTAConeBasicClustersSKJetsRing2TobN);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Ring3TobN", &gepLeadingWTAConeBasicClustersSKJetsRing3TobN);
    gepLeadingWTAConeBasicClustersSKJetsTree->Branch("Ring4TobN", &gepLeadingWTAConeBasicClustersSKJetsRing4TobN);

    // gep wta cone basic clusters SK subleading jets
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Et", &gepSubleadingWTAConeGEPBasicClustersSKJetspTValues);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Eta", &gepSubleadingWTAConeGEPBasicClustersSKJetsEtaValues);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Phi", &gepSubleadingWTAConeGEPBasicClustersSKJetsPhiValues);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("NConstituents", &gepSubleadingWTAConeGEPBasicClustersSKJetsNConstituentsValues);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Ring0Et", &gepSubleadingWTAConeBasicClustersSKJetsRing0Et);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Ring1Et", &gepSubleadingWTAConeBasicClustersSKJetsRing1Et);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Ring2Et", &gepSubleadingWTAConeBasicClustersSKJetsRing2Et);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Ring3Et", &gepSubleadingWTAConeBasicClustersSKJetsRing3Et);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Ring4Et", &gepSubleadingWTAConeBasicClustersSKJetsRing4Et);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("TotalTobN", &gepSubleadingWTAConeBasicClustersSKJetsTotalTobN);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Ring0TobN", &gepSubleadingWTAConeBasicClustersSKJetsRing0TobN);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Ring1TobN", &gepSubleadingWTAConeBasicClustersSKJetsRing1TobN);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Ring2TobN", &gepSubleadingWTAConeBasicClustersSKJetsRing2TobN);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Ring3TobN", &gepSubleadingWTAConeBasicClustersSKJetsRing3TobN);
    gepSubleadingWTAConeBasicClustersSKJetsTree->Branch("Ring4TobN", &gepSubleadingWTAConeBasicClustersSKJetsRing4TobN);

    // gep wta cone cells towers EtaSK jets
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Et", &gepWTAConeCellsTowersEtaSKJetspTValues);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Eta", &gepWTAConeCellsTowersEtaSKJetsEtaValues);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Phi", &gepWTAConeCellsTowersEtaSKJetsPhiValues);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("NConstituents", &gepWTAConeCellsTowersEtaSKJetsNConstituentsValues);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Ring0Et", &gepWTAConeCellsTowersEtaSKJetsRing0Et);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Ring1Et", &gepWTAConeCellsTowersEtaSKJetsRing1Et);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Ring2Et", &gepWTAConeCellsTowersEtaSKJetsRing2Et);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Ring3Et", &gepWTAConeCellsTowersEtaSKJetsRing3Et);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Ring4Et", &gepWTAConeCellsTowersEtaSKJetsRing4Et);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("TotalTobN", &gepWTAConeCellsTowersEtaSKJetsTotalTobN);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Ring0TobN", &gepWTAConeCellsTowersEtaSKJetsRing0TobN);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Ring1TobN", &gepWTAConeCellsTowersEtaSKJetsRing1TobN);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Ring2TobN", &gepWTAConeCellsTowersEtaSKJetsRing2TobN);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Ring3TobN", &gepWTAConeCellsTowersEtaSKJetsRing3TobN);
    gepWTAConeCellsTowersEtaSKJetsTree->Branch("Ring4TobN", &gepWTAConeCellsTowersEtaSKJetsRing4TobN);

    // gep wta cone cells towers EtaSK leading jets
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Et", &gepLeadingWTAConeCellsTowersEtaSKJetspTValues);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Eta", &gepLeadingWTAConeCellsTowersEtaSKJetsEtaValues);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Phi", &gepLeadingWTAConeCellsTowersEtaSKJetsPhiValues);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("NConstituents", &gepLeadingWTAConeCellsTowersEtaSKJetsNConstituentsValues);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring0Et", &gepLeadingWTAConeCellsTowersEtaSKJetsRing0Et);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring1Et", &gepLeadingWTAConeCellsTowersEtaSKJetsRing1Et);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring2Et", &gepLeadingWTAConeCellsTowersEtaSKJetsRing2Et);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring3Et", &gepLeadingWTAConeCellsTowersEtaSKJetsRing3Et);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring4Et", &gepLeadingWTAConeCellsTowersEtaSKJetsRing4Et);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("TotalTobN", &gepLeadingWTAConeCellsTowersEtaSKJetsTotalTobN);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring0TobN", &gepLeadingWTAConeCellsTowersEtaSKJetsRing0TobN);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring1TobN", &gepLeadingWTAConeCellsTowersEtaSKJetsRing1TobN);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring2TobN", &gepLeadingWTAConeCellsTowersEtaSKJetsRing2TobN);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring3TobN", &gepLeadingWTAConeCellsTowersEtaSKJetsRing3TobN);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring4TobN", &gepLeadingWTAConeCellsTowersEtaSKJetsRing4TobN);

    // gep wta cone cells towers EtaSK subleading jets
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Et", &gepSubleadingWTAConeCellsTowersEtaSKJetspTValues);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Eta", &gepSubleadingWTAConeCellsTowersEtaSKJetsEtaValues);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Phi", &gepSubleadingWTAConeCellsTowersEtaSKJetsPhiValues);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("NConstituents", &gepSubleadingWTAConeCellsTowersEtaSKJetsNConstituentsValues);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring0Et", &gepSubleadingWTAConeCellsTowersEtaSKJetsRing0Et);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring1Et", &gepSubleadingWTAConeCellsTowersEtaSKJetsRing1Et);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring2Et", &gepSubleadingWTAConeCellsTowersEtaSKJetsRing2Et);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring3Et", &gepSubleadingWTAConeCellsTowersEtaSKJetsRing3Et);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring4Et", &gepSubleadingWTAConeCellsTowersEtaSKJetsRing4Et);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("TotalTobN", &gepSubleadingWTAConeCellsTowersEtaSKJetsTotalTobN);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring0TobN", &gepSubleadingWTAConeCellsTowersEtaSKJetsRing0TobN);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring1TobN", &gepSubleadingWTAConeCellsTowersEtaSKJetsRing1TobN);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring2TobN", &gepSubleadingWTAConeCellsTowersEtaSKJetsRing2TobN);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring3TobN", &gepSubleadingWTAConeCellsTowersEtaSKJetsRing3TobN);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Branch("Ring4TobN", &gepSubleadingWTAConeCellsTowersEtaSKJetsRing4TobN);

    // gep wta cone basic clusters EtaSK jets
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Et", &gepWTAConeBasicClustersEtaSKJetspTValues);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Eta", &gepWTAConeBasicClustersEtaSKJetsEtaValues);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Phi", &gepWTAConeBasicClustersEtaSKJetsPhiValues);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("NConstituents", &gepWTAConeBasicClustersEtaSKJetsNConstituentsValues);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Ring0Et", &gepWTAConeBasicClustersEtaSKJetsRing0Et);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Ring1Et", &gepWTAConeBasicClustersEtaSKJetsRing1Et);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Ring2Et", &gepWTAConeBasicClustersEtaSKJetsRing2Et);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Ring3Et", &gepWTAConeBasicClustersEtaSKJetsRing3Et);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Ring4Et", &gepWTAConeBasicClustersEtaSKJetsRing4Et);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("TotalTobN", &gepWTAConeBasicClustersEtaSKJetsTotalTobN);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Ring0TobN", &gepWTAConeBasicClustersEtaSKJetsRing0TobN);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Ring1TobN", &gepWTAConeBasicClustersEtaSKJetsRing1TobN);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Ring2TobN", &gepWTAConeBasicClustersEtaSKJetsRing2TobN);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Ring3TobN", &gepWTAConeBasicClustersEtaSKJetsRing3TobN);
    gepWTAConeBasicClustersEtaSKJetsTree->Branch("Ring4TobN", &gepWTAConeBasicClustersEtaSKJetsRing4TobN);

    // gep wta cone basic clusters EtaSK leading jets
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Et", &gepLeadingWTAConeBasicClustersEtaSKJetspTValues);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Eta", &gepLeadingWTAConeBasicClustersEtaSKJetsEtaValues);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Phi", &gepLeadingWTAConeBasicClustersEtaSKJetsPhiValues);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("NConstituents", &gepLeadingWTAConeBasicClustersEtaSKJetsNConstituentsValues);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring0Et", &gepLeadingWTAConeBasicClustersEtaSKJetsRing0Et);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring1Et", &gepLeadingWTAConeBasicClustersEtaSKJetsRing1Et);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring2Et", &gepLeadingWTAConeBasicClustersEtaSKJetsRing2Et);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring3Et", &gepLeadingWTAConeBasicClustersEtaSKJetsRing3Et);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring4Et", &gepLeadingWTAConeBasicClustersEtaSKJetsRing4Et);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("TotalTobN", &gepLeadingWTAConeBasicClustersEtaSKJetsTotalTobN);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring0TobN", &gepLeadingWTAConeBasicClustersEtaSKJetsRing0TobN);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring1TobN", &gepLeadingWTAConeBasicClustersEtaSKJetsRing1TobN);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring2TobN", &gepLeadingWTAConeBasicClustersEtaSKJetsRing2TobN);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring3TobN", &gepLeadingWTAConeBasicClustersEtaSKJetsRing3TobN);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring4TobN", &gepLeadingWTAConeBasicClustersEtaSKJetsRing4TobN);

    // gep wta cone basic clusters EtaSK subleading jets
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Et", &gepSubleadingWTAConeBasicClustersEtaSKJetspTValues);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Eta", &gepSubleadingWTAConeBasicClustersEtaSKJetsEtaValues);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Phi", &gepSubleadingWTAConeBasicClustersEtaSKJetsPhiValues);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("NConstituents", &gepSubleadingWTAConeBasicClustersEtaSKJetsNConstituentsValues);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring0Et", &gepSubleadingWTAConeBasicClustersEtaSKJetsRing0Et);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring1Et", &gepSubleadingWTAConeBasicClustersEtaSKJetsRing1Et);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring2Et", &gepSubleadingWTAConeBasicClustersEtaSKJetsRing2Et);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring3Et", &gepSubleadingWTAConeBasicClustersEtaSKJetsRing3Et);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring4Et", &gepSubleadingWTAConeBasicClustersEtaSKJetsRing4Et);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("TotalTobN", &gepSubleadingWTAConeBasicClustersEtaSKJetsTotalTobN);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring0TobN", &gepSubleadingWTAConeBasicClustersEtaSKJetsRing0TobN);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring1TobN", &gepSubleadingWTAConeBasicClustersEtaSKJetsRing1TobN);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring2TobN", &gepSubleadingWTAConeBasicClustersEtaSKJetsRing2TobN);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring3TobN", &gepSubleadingWTAConeBasicClustersEtaSKJetsRing3TobN);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Branch("Ring4TobN", &gepSubleadingWTAConeBasicClustersEtaSKJetsRing4TobN);

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

    // outOfTimeAntiKt4TruthJetsTree
    outOfTimeAntiKt4TruthJetsTree->Branch("EtIndex", &outOfTimeAntiKt4TruthSRJEtIndexValues);
    outOfTimeAntiKt4TruthJetsTree->Branch("Et", &outOfTimeAntiKt4TruthSRJEtValues);
    outOfTimeAntiKt4TruthJetsTree->Branch("Eta", &outOfTimeAntiKt4TruthSRJEtaValues);
    outOfTimeAntiKt4TruthJetsTree->Branch("Phi", &outOfTimeAntiKt4TruthSRJPhiValues);

    // leadingOutOfTimeAntiKt4TruthJetsTree
    leadingOutOfTimeAntiKt4TruthJetsTree->Branch("Et", &outOfTimeAntiKt4TruthSRJLeadingEtValues);
    leadingOutOfTimeAntiKt4TruthJetsTree->Branch("Eta", &outOfTimeAntiKt4TruthSRJLeadingEtaValues);
    leadingOutOfTimeAntiKt4TruthJetsTree->Branch("Phi", &outOfTimeAntiKt4TruthSRJLeadingPhiValues);

    // subleadingOutOfTimeAntiKt4TruthJetsTree
    subleadingOutOfTimeAntiKt4TruthJetsTree->Branch("Et", &outOfTimeAntiKt4TruthSRJSubleadingEtValues);
    subleadingOutOfTimeAntiKt4TruthJetsTree->Branch("Eta", &outOfTimeAntiKt4TruthSRJSubleadingEtaValues);
    subleadingOutOfTimeAntiKt4TruthJetsTree->Branch("Phi", &outOfTimeAntiKt4TruthSRJSubleadingPhiValues);

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

    // jFexSRJSimTree
    jFexSRJSimTree->Branch("EtIndex", &jFexSRJSimEtIndexValues);
    jFexSRJSimTree->Branch("Et", &jFexSRJSimEtValues);
    jFexSRJSimTree->Branch("Eta", &jFexSRJSimEtaValues);
    jFexSRJSimTree->Branch("Phi", &jFexSRJSimPhiValues);

    // jFexLeadingSRJSimTree
    jFexLeadingSRJSimTree->Branch("Et", &jFexSRJSimLeadingEtValues);
    jFexLeadingSRJSimTree->Branch("Eta", &jFexSRJSimLeadingEtaValues);
    jFexLeadingSRJSimTree->Branch("Phi", &jFexSRJSimLeadingPhiValues);

    // jFexSubleadingSRJSimTree
    jFexSubleadingSRJSimTree->Branch("Et", &jFexSRJSimSubleadingEtValues);
    jFexSubleadingSRJSimTree->Branch("Eta", &jFexSRJSimSubleadingEtaValues);
    jFexSubleadingSRJSimTree->Branch("Phi", &jFexSRJSimSubleadingPhiValues);

    // gFexLRJSimTree
    gFexLRJSimTree->Branch("EtIndex", &gFexLRJSimEtIndexValues);
    gFexLRJSimTree->Branch("Et", &gFexLRJSimEtValues);
    gFexLRJSimTree->Branch("Eta", &gFexLRJSimEtaValues);
    gFexLRJSimTree->Branch("Phi", &gFexLRJSimPhiValues);

    // gFexLeadingLRJSimTree
    gFexLeadingLRJSimTree->Branch("Et", &gFexLRJSimLeadingEtValues);
    gFexLeadingLRJSimTree->Branch("Eta", &gFexLRJSimLeadingEtaValues);
    gFexLeadingLRJSimTree->Branch("Phi", &gFexLRJSimLeadingPhiValues);

    // gFexSubleadingLRJSimTree
    gFexSubleadingLRJSimTree->Branch("Et", &gFexLRJSimSubleadingEtValues);
    gFexSubleadingLRJSimTree->Branch("Eta", &gFexLRJSimSubleadingEtaValues);
    gFexSubleadingLRJSimTree->Branch("Phi", &gFexLRJSimSubleadingPhiValues);

    // TODO: uncomment once L1_gFexSRJetRoISim branch exists in GEP ntuples
    // // gFexSRJSimTree
    // gFexSRJSimTree->Branch("EtIndex", &gFexSRJSimEtIndexValues);
    // gFexSRJSimTree->Branch("Et", &gFexSRJSimEtValues);
    // gFexSRJSimTree->Branch("Eta", &gFexSRJSimEtaValues);
    // gFexSRJSimTree->Branch("Phi", &gFexSRJSimPhiValues);

    // // gFexLeadingSRJSimTree
    // gFexLeadingSRJSimTree->Branch("Et", &gFexSRJSimLeadingEtValues);
    // gFexLeadingSRJSimTree->Branch("Eta", &gFexSRJSimLeadingEtaValues);
    // gFexLeadingSRJSimTree->Branch("Phi", &gFexSRJSimLeadingPhiValues);

    // // gFexSubleadingSRJSimTree
    // gFexSubleadingSRJSimTree->Branch("Et", &gFexSRJSimSubleadingEtValues);
    // gFexSubleadingSRJSimTree->Branch("Eta", &gFexSRJSimSubleadingEtaValues);
    // gFexSubleadingSRJSimTree->Branch("Phi", &gFexSRJSimSubleadingPhiValues);

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
    metTruthTree->Branch("metTruthIntMuonsX",     &metTruthIntMuonsX);
    metTruthTree->Branch("metTruthIntMuonsY",     &metTruthIntMuonsY);
    metTruthTree->Branch("metTruthIntMuonsSumET", &metTruthIntMuonsSumET);
    metTruthTree->Branch("metTruthIntMuons",      &metTruthIntMuons);

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

    // gFexMETJwoJSimTree branches
    gFexMETJwoJSimTree->Branch("gMET",    &gFexMETJwoJSim);
    gFexMETJwoJSimTree->Branch("gMETPhi", &gFexMETPhiJwoJSim);
    gFexMETJwoJSimTree->Branch("gSumET",  &gFexSumETJwoJSim);

    // gFexMETNoiseCutSimTree branches
    gFexMETNoiseCutSimTree->Branch("gMET",    &gFexMETNoiseCutSim);
    gFexMETNoiseCutSimTree->Branch("gMETPhi", &gFexMETPhiNoiseCutSim);
    gFexMETNoiseCutSimTree->Branch("gSumET",  &gFexSumETNoiseCutSim);

    // gFexMETRmsSimTree branches
    gFexMETRmsSimTree->Branch("gMET",    &gFexMETRmsSim);
    gFexMETRmsSimTree->Branch("gMETPhi", &gFexMETPhiRmsSim);
    gFexMETRmsSimTree->Branch("gSumET",  &gFexSumETRmsSim);

    // jFexMETTree branches
    jFexMETTree->Branch("jMET", &jFexMET);

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

    if (!specificDaodFile.empty()) {
        // Condor single-file mode: use the provided DAOD file directly
        fileNames.push_back(specificDaodFile);
        std::cout << "Single-file mode: DAOD = " << specificDaodFile << "\n";
    } else {
        // Directory scan mode: collect DAOD_JETM42*.root files
        for (const auto& entry : fs::directory_iterator(fileDir)) {
            if (!entry.is_regular_file()) continue;
            const std::string fn = entry.path().filename().string();
            if (has_prefix(fn, "DAOD_JETM42") && has_suffix(fn, ".root")) {
                std::cout << "found DAOD file: " << entry.path().string() << "\n";
                fileNames.push_back(entry.path().string());
            }
        }
    }

    if (!specificGepFile.empty()) {
        // Condor single-file mode: use the provided GEP file directly
        gepFileNames.push_back(specificGepFile);
        std::cout << "Single-file mode: GEP  = " << specificGepFile << "\n";
    } else if (specificDaodFile.empty()) {
        // Directory scan mode: collect GEP*.root files (co-located with DAOD files)
        for (const auto& entry : fs::directory_iterator(fileDir)) {
            if (!entry.is_regular_file()) continue;
            const std::string fn = entry.path().filename().string();
            if (has_prefix(fn, "GEP") && has_suffix(fn, ".root")) {
                std::cout << "found TrigGepPerf file name: " << entry.path().string() << "\n";
                gepFileNames.push_back(entry.path().string());
            }
        }
    }

    if (fileNames.empty()) {
        std::cout << "No ROOT files found in directory." << endl;
        return;
    }

    std::cout << "Found " << fileNames.size() << " DAOD file(s), " << gepFileNames.size() << " GEP file(s)." << endl;

    // Main processing loop
    int fileIt = 0;
    for (const auto& fileName : fileNames) {
        fileIt++; 
        //if (fileIt > 9) break; 
        std::cout << "Processing file: " << fileName << endl;

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

        // --- Open GEP ROOT (non-xAOD); required to have a matching GEP file ---
        if (fileIt - 1 >= (int)gepFileNames.size()) {
            std::cerr << "Error: no matching GEP file for DAOD file: " << fileName << "\n";
            std::cerr << "  " << gepFileNames.size() << " GEP file(s) provided for " << fileNames.size() << " DAOD file(s).\n";
            f->Close(); delete f;
            continue;
        }
        TFile* gf = TFile::Open(gepFileNames[fileIt - 1].c_str());
        if (!gf || gf->IsZombie()) {
            std::cerr << "Could not open GEP file: " << gepFileNames[fileIt - 1] << "\n";
            if (gf) { gf->Close(); delete gf; }
            f->Close(); delete f;
            continue;
        }
        std::cout << "Processing GEP file: " << gepFileNames[fileIt - 1] << "\n";

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

        gt->SetBranchStatus("eventNumber", 1);
        gt->SetBranchStatus("runNumber",   1);
        if(!signalBool) gt->SetBranchStatus("filterHSTP", 1);

        gt->SetBranchStatus("GEPCellsTowerEtaSK_et",  1);
        gt->SetBranchStatus("GEPCellsTowerEtaSK_eta", 1);
        gt->SetBranchStatus("GEPCellsTowerEtaSK_phi", 1);

        gt->SetBranchStatus("GEPBasicClustersEtaSK_et",  1);
        gt->SetBranchStatus("GEPBasicClustersEtaSK_eta", 1);
        gt->SetBranchStatus("GEPBasicClustersEtaSK_phi", 1);

        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_ring0_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_ring1_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_ring2_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_ring3_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_ring4_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_total_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_ring0_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_ring1_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_ring2_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_ring3_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerJets_ring4_TobN", 1);

        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_ring0_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_ring1_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_ring2_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_ring3_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_ring4_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_total_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_ring0_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_ring1_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_ring2_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_ring3_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerSKJets_ring4_TobN", 1);

        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_pt",  1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_eta", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_phi", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_nConstituents", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_ring0_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_ring1_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_ring2_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_ring3_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_ring4_Et", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_total_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_ring0_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_ring1_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_ring2_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_ring3_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPCellsTowerEtaSKJets_ring4_TobN", 1);

        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_ring0_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_ring1_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_ring2_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_ring3_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_ring4_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_total_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_ring0_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_ring1_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_ring2_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_ring3_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersJets_ring4_TobN", 1);

        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_ring0_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_ring1_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_ring2_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_ring3_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_ring4_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_total_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_ring0_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_ring1_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_ring2_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_ring3_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersSKJets_ring4_TobN", 1);

        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_pt",  1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_eta", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_phi", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_nConstituents", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_ring0_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_ring1_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_ring2_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_ring3_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_ring4_Et", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_total_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_ring0_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_ring1_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_ring2_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_ring3_TobN", 1);
        gt->SetBranchStatus("WTAConeGEPBasicClustersEtaSKJets_ring4_TobN", 1);

        gt->SetBranchStatus("L1_jFexSRJetRoISim_et",  1);
        gt->SetBranchStatus("L1_jFexSRJetRoISim_eta", 1);
        gt->SetBranchStatus("L1_jFexSRJetRoISim_phi", 1);

        gt->SetBranchStatus("L1_gFexLRJetRoISim_et",  1);
        gt->SetBranchStatus("L1_gFexLRJetRoISim_eta", 1);
        gt->SetBranchStatus("L1_gFexLRJetRoISim_phi", 1);

        gt->SetBranchStatus("gFex_JwoJ_met",    1);
        gt->SetBranchStatus("gFex_JwoJ_metPhi", 1);
        gt->SetBranchStatus("gFex_JwoJ_sumEt",  1);
        gt->SetBranchStatus("gFex_NC_met",      1);
        gt->SetBranchStatus("gFex_NC_metPhi",   1);
        gt->SetBranchStatus("gFex_NC_sumEt",    1);
        gt->SetBranchStatus("gFex_RhoRMS_met",    1);
        gt->SetBranchStatus("gFex_RhoRMS_metPhi", 1);
        gt->SetBranchStatus("gFex_RhoRMS_sumEt",  1);
        gt->SetBranchStatus("jFex_met", 1);
        // TODO: uncomment once L1_gFexSRJetRoISim branch exists in GEP ntuples
        // gt->SetBranchStatus("L1_gFexSRJetRoISim_et",  1);
        // gt->SetBranchStatus("L1_gFexSRJetRoISim_eta", 1);
        // gt->SetBranchStatus("L1_gFexSRJetRoISim_phi", 1);

        // Local pointers for new GEP branches
        std::vector<float>* gepIn_BasicClustersEt  = nullptr;
        std::vector<float>* gepIn_BasicClustersEta = nullptr;
        std::vector<float>* gepIn_BasicClustersPhi = nullptr;

        std::vector<float>* gepIn_BasicClustersSKEt  = nullptr;
        std::vector<float>* gepIn_BasicClustersSKEta = nullptr;
        std::vector<float>* gepIn_BasicClustersSKPhi = nullptr;

        std::vector<float>* gepIn_CellsTowersEt  = nullptr;
        std::vector<float>* gepIn_CellsTowersEta = nullptr;
        std::vector<float>* gepIn_CellsTowersPhi = nullptr;

        std::vector<float>* gepIn_CellsTowersSKEt  = nullptr;
        std::vector<float>* gepIn_CellsTowersSKEta = nullptr;
        std::vector<float>* gepIn_CellsTowersSKPhi = nullptr;

        std::vector<float>* gepIn_WTAConeCellsTowersJetsPt  = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersJetsEta = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersJetsPhi = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersJetsNConstituents = nullptr;

        std::vector<float>* gepIn_WTAConeCellsTowersSKJetsPt  = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersSKJetsEta = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersSKJetsPhi = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersSKJetsNConstituents = nullptr;

        std::vector<float>* gepIn_WTAConeBasicClustersJetsPt  = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersJetsEta = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersJetsPhi = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersJetsNConstituents = nullptr;

        std::vector<float>* gepIn_WTAConeBasicClustersSKJetsPt  = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersSKJetsEta = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersSKJetsPhi = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersSKJetsNConstituents = nullptr;

        gt->SetBranchAddress("GEPBasicClusters_et",  &gepIn_BasicClustersEt);
        gt->SetBranchAddress("GEPBasicClusters_eta", &gepIn_BasicClustersEta);
        gt->SetBranchAddress("GEPBasicClusters_phi", &gepIn_BasicClustersPhi);

        gt->SetBranchAddress("GEPBasicClustersSK_et",  &gepIn_BasicClustersSKEt);
        gt->SetBranchAddress("GEPBasicClustersSK_eta", &gepIn_BasicClustersSKEta);
        gt->SetBranchAddress("GEPBasicClustersSK_phi", &gepIn_BasicClustersSKPhi);

        gt->SetBranchAddress("GEPCellsTower_et",  &gepIn_CellsTowersEt);
        gt->SetBranchAddress("GEPCellsTower_eta", &gepIn_CellsTowersEta);
        gt->SetBranchAddress("GEPCellsTower_phi", &gepIn_CellsTowersPhi);

        gt->SetBranchAddress("GEPCellsTowerSK_et",  &gepIn_CellsTowersSKEt);
        gt->SetBranchAddress("GEPCellsTowerSK_eta", &gepIn_CellsTowersSKEta);
        gt->SetBranchAddress("GEPCellsTowerSK_phi", &gepIn_CellsTowersSKPhi);

        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_pt",  &gepIn_WTAConeCellsTowersJetsPt);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_eta", &gepIn_WTAConeCellsTowersJetsEta);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_phi", &gepIn_WTAConeCellsTowersJetsPhi); 
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_nConstituents", &gepIn_WTAConeCellsTowersJetsNConstituents); 

        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_pt",  &gepIn_WTAConeCellsTowersSKJetsPt);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_eta", &gepIn_WTAConeCellsTowersSKJetsEta);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_phi", &gepIn_WTAConeCellsTowersSKJetsPhi); 
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_nConstituents", &gepIn_WTAConeCellsTowersSKJetsNConstituents); 

        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_pt",  &gepIn_WTAConeBasicClustersJetsPt);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_eta", &gepIn_WTAConeBasicClustersJetsEta);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_phi", &gepIn_WTAConeBasicClustersJetsPhi); 
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_nConstituents", &gepIn_WTAConeBasicClustersJetsNConstituents); 

        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_pt",  &gepIn_WTAConeBasicClustersSKJetsPt);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_eta", &gepIn_WTAConeBasicClustersSKJetsEta);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_phi", &gepIn_WTAConeBasicClustersSKJetsPhi);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_nConstituents", &gepIn_WTAConeBasicClustersSKJetsNConstituents);
        // Local pointers for new GEP branches
        std::vector<float>* gepIn_CellsTowersEtaSKEt  = nullptr;
        std::vector<float>* gepIn_CellsTowersEtaSKEta = nullptr;
        std::vector<float>* gepIn_CellsTowersEtaSKPhi = nullptr;

        std::vector<float>* gepIn_BasicClustersEtaSKEt  = nullptr;
        std::vector<float>* gepIn_BasicClustersEtaSKEta = nullptr;
        std::vector<float>* gepIn_BasicClustersEtaSKPhi = nullptr;

        std::vector<float>* gepIn_WTAConeCellsTowersJetsRing0Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersJetsRing1Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersJetsRing2Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersJetsRing3Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersJetsRing4Et = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersJetsTotalTobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersJetsRing0TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersJetsRing1TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersJetsRing2TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersJetsRing3TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersJetsRing4TobN = nullptr;

        std::vector<float>* gepIn_WTAConeCellsTowersSKJetsRing0Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersSKJetsRing1Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersSKJetsRing2Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersSKJetsRing3Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersSKJetsRing4Et = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersSKJetsTotalTobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersSKJetsRing0TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersSKJetsRing1TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersSKJetsRing2TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersSKJetsRing3TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersSKJetsRing4TobN = nullptr;

        std::vector<float>* gepIn_WTAConeCellsTowersEtaSKJetsPt  = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersEtaSKJetsEta = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersEtaSKJetsPhi = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersEtaSKJetsNConstituents = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersEtaSKJetsRing0Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersEtaSKJetsRing1Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersEtaSKJetsRing2Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersEtaSKJetsRing3Et = nullptr;
        std::vector<float>* gepIn_WTAConeCellsTowersEtaSKJetsRing4Et = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersEtaSKJetsTotalTobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersEtaSKJetsRing0TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersEtaSKJetsRing1TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersEtaSKJetsRing2TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersEtaSKJetsRing3TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeCellsTowersEtaSKJetsRing4TobN = nullptr;

        std::vector<float>* gepIn_WTAConeBasicClustersJetsRing0Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersJetsRing1Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersJetsRing2Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersJetsRing3Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersJetsRing4Et = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersJetsTotalTobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersJetsRing0TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersJetsRing1TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersJetsRing2TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersJetsRing3TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersJetsRing4TobN = nullptr;

        std::vector<float>* gepIn_WTAConeBasicClustersSKJetsRing0Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersSKJetsRing1Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersSKJetsRing2Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersSKJetsRing3Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersSKJetsRing4Et = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersSKJetsTotalTobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersSKJetsRing0TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersSKJetsRing1TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersSKJetsRing2TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersSKJetsRing3TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersSKJetsRing4TobN = nullptr;

        std::vector<float>* gepIn_WTAConeBasicClustersEtaSKJetsPt  = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersEtaSKJetsEta = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersEtaSKJetsPhi = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersEtaSKJetsNConstituents = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersEtaSKJetsRing0Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersEtaSKJetsRing1Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersEtaSKJetsRing2Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersEtaSKJetsRing3Et = nullptr;
        std::vector<float>* gepIn_WTAConeBasicClustersEtaSKJetsRing4Et = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersEtaSKJetsTotalTobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersEtaSKJetsRing0TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersEtaSKJetsRing1TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersEtaSKJetsRing2TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersEtaSKJetsRing3TobN = nullptr;
        std::vector<int>*   gepIn_WTAConeBasicClustersEtaSKJetsRing4TobN = nullptr;

        std::vector<float>* jFexSRJSimEt  = nullptr;
        std::vector<float>* jFexSRJSimEta = nullptr;
        std::vector<float>* jFexSRJSimPhi = nullptr;

        std::vector<float>* gFexLRJSimEt  = nullptr;
        std::vector<float>* gFexLRJSimEta = nullptr;
        std::vector<float>* gFexLRJSimPhi = nullptr;

        float gepIn_gFexJwoJMet = 0.0f,    gepIn_gFexJwoJMetPhi = 0.0f,    gepIn_gFexJwoJSumEt = 0.0f;
        float gepIn_gFexNCMet = 0.0f,      gepIn_gFexNCMetPhi = 0.0f,      gepIn_gFexNCSumEt = 0.0f;
        float gepIn_gFexRhoRMSMet = 0.0f,  gepIn_gFexRhoRMSMetPhi = 0.0f,  gepIn_gFexRhoRMSSumEt = 0.0f;
        float gepIn_jFexMet = 0.0f;
        // TODO: uncomment once L1_gFexSRJetRoISim branch exists in GEP ntuples
        // std::vector<float>* gFexSRJSimEt  = nullptr;
        // std::vector<float>* gFexSRJSimEta = nullptr;
        // std::vector<float>* gFexSRJSimPhi = nullptr;

        bool filterHSTP;

        gt->SetBranchAddress("eventNumber", &gepEventNumberIn);
        if(!signalBool) gt->SetBranchAddress("filterHSTP", &filterHSTP);
        gt->SetBranchAddress("runNumber",   &gepRunNumberIn);

        gt->SetBranchAddress("GEPCellsTowerEtaSK_et",  &gepIn_CellsTowersEtaSKEt);
        gt->SetBranchAddress("GEPCellsTowerEtaSK_eta", &gepIn_CellsTowersEtaSKEta);
        gt->SetBranchAddress("GEPCellsTowerEtaSK_phi", &gepIn_CellsTowersEtaSKPhi);

        gt->SetBranchAddress("GEPBasicClustersEtaSK_et",  &gepIn_BasicClustersEtaSKEt);
        gt->SetBranchAddress("GEPBasicClustersEtaSK_eta", &gepIn_BasicClustersEtaSKEta);
        gt->SetBranchAddress("GEPBasicClustersEtaSK_phi", &gepIn_BasicClustersEtaSKPhi);

        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_ring0_Et", &gepIn_WTAConeCellsTowersJetsRing0Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_ring1_Et", &gepIn_WTAConeCellsTowersJetsRing1Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_ring2_Et", &gepIn_WTAConeCellsTowersJetsRing2Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_ring3_Et", &gepIn_WTAConeCellsTowersJetsRing3Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_ring4_Et", &gepIn_WTAConeCellsTowersJetsRing4Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_total_TobN", &gepIn_WTAConeCellsTowersJetsTotalTobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_ring0_TobN", &gepIn_WTAConeCellsTowersJetsRing0TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_ring1_TobN", &gepIn_WTAConeCellsTowersJetsRing1TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_ring2_TobN", &gepIn_WTAConeCellsTowersJetsRing2TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_ring3_TobN", &gepIn_WTAConeCellsTowersJetsRing3TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerJets_ring4_TobN", &gepIn_WTAConeCellsTowersJetsRing4TobN);

        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_ring0_Et", &gepIn_WTAConeCellsTowersSKJetsRing0Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_ring1_Et", &gepIn_WTAConeCellsTowersSKJetsRing1Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_ring2_Et", &gepIn_WTAConeCellsTowersSKJetsRing2Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_ring3_Et", &gepIn_WTAConeCellsTowersSKJetsRing3Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_ring4_Et", &gepIn_WTAConeCellsTowersSKJetsRing4Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_total_TobN", &gepIn_WTAConeCellsTowersSKJetsTotalTobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_ring0_TobN", &gepIn_WTAConeCellsTowersSKJetsRing0TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_ring1_TobN", &gepIn_WTAConeCellsTowersSKJetsRing1TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_ring2_TobN", &gepIn_WTAConeCellsTowersSKJetsRing2TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_ring3_TobN", &gepIn_WTAConeCellsTowersSKJetsRing3TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerSKJets_ring4_TobN", &gepIn_WTAConeCellsTowersSKJetsRing4TobN);

        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_pt",  &gepIn_WTAConeCellsTowersEtaSKJetsPt);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_eta", &gepIn_WTAConeCellsTowersEtaSKJetsEta);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_phi", &gepIn_WTAConeCellsTowersEtaSKJetsPhi);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_nConstituents", &gepIn_WTAConeCellsTowersEtaSKJetsNConstituents);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_ring0_Et", &gepIn_WTAConeCellsTowersEtaSKJetsRing0Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_ring1_Et", &gepIn_WTAConeCellsTowersEtaSKJetsRing1Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_ring2_Et", &gepIn_WTAConeCellsTowersEtaSKJetsRing2Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_ring3_Et", &gepIn_WTAConeCellsTowersEtaSKJetsRing3Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_ring4_Et", &gepIn_WTAConeCellsTowersEtaSKJetsRing4Et);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_total_TobN", &gepIn_WTAConeCellsTowersEtaSKJetsTotalTobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_ring0_TobN", &gepIn_WTAConeCellsTowersEtaSKJetsRing0TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_ring1_TobN", &gepIn_WTAConeCellsTowersEtaSKJetsRing1TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_ring2_TobN", &gepIn_WTAConeCellsTowersEtaSKJetsRing2TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_ring3_TobN", &gepIn_WTAConeCellsTowersEtaSKJetsRing3TobN);
        gt->SetBranchAddress("WTAConeGEPCellsTowerEtaSKJets_ring4_TobN", &gepIn_WTAConeCellsTowersEtaSKJetsRing4TobN);

        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_ring0_Et", &gepIn_WTAConeBasicClustersJetsRing0Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_ring1_Et", &gepIn_WTAConeBasicClustersJetsRing1Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_ring2_Et", &gepIn_WTAConeBasicClustersJetsRing2Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_ring3_Et", &gepIn_WTAConeBasicClustersJetsRing3Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_ring4_Et", &gepIn_WTAConeBasicClustersJetsRing4Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_total_TobN", &gepIn_WTAConeBasicClustersJetsTotalTobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_ring0_TobN", &gepIn_WTAConeBasicClustersJetsRing0TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_ring1_TobN", &gepIn_WTAConeBasicClustersJetsRing1TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_ring2_TobN", &gepIn_WTAConeBasicClustersJetsRing2TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_ring3_TobN", &gepIn_WTAConeBasicClustersJetsRing3TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersJets_ring4_TobN", &gepIn_WTAConeBasicClustersJetsRing4TobN);

        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_ring0_Et", &gepIn_WTAConeBasicClustersSKJetsRing0Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_ring1_Et", &gepIn_WTAConeBasicClustersSKJetsRing1Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_ring2_Et", &gepIn_WTAConeBasicClustersSKJetsRing2Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_ring3_Et", &gepIn_WTAConeBasicClustersSKJetsRing3Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_ring4_Et", &gepIn_WTAConeBasicClustersSKJetsRing4Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_total_TobN", &gepIn_WTAConeBasicClustersSKJetsTotalTobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_ring0_TobN", &gepIn_WTAConeBasicClustersSKJetsRing0TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_ring1_TobN", &gepIn_WTAConeBasicClustersSKJetsRing1TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_ring2_TobN", &gepIn_WTAConeBasicClustersSKJetsRing2TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_ring3_TobN", &gepIn_WTAConeBasicClustersSKJetsRing3TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersSKJets_ring4_TobN", &gepIn_WTAConeBasicClustersSKJetsRing4TobN);

        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_pt",  &gepIn_WTAConeBasicClustersEtaSKJetsPt);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_eta", &gepIn_WTAConeBasicClustersEtaSKJetsEta);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_phi", &gepIn_WTAConeBasicClustersEtaSKJetsPhi);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_nConstituents", &gepIn_WTAConeBasicClustersEtaSKJetsNConstituents);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_ring0_Et", &gepIn_WTAConeBasicClustersEtaSKJetsRing0Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_ring1_Et", &gepIn_WTAConeBasicClustersEtaSKJetsRing1Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_ring2_Et", &gepIn_WTAConeBasicClustersEtaSKJetsRing2Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_ring3_Et", &gepIn_WTAConeBasicClustersEtaSKJetsRing3Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_ring4_Et", &gepIn_WTAConeBasicClustersEtaSKJetsRing4Et);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_total_TobN", &gepIn_WTAConeBasicClustersEtaSKJetsTotalTobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_ring0_TobN", &gepIn_WTAConeBasicClustersEtaSKJetsRing0TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_ring1_TobN", &gepIn_WTAConeBasicClustersEtaSKJetsRing1TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_ring2_TobN", &gepIn_WTAConeBasicClustersEtaSKJetsRing2TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_ring3_TobN", &gepIn_WTAConeBasicClustersEtaSKJetsRing3TobN);
        gt->SetBranchAddress("WTAConeGEPBasicClustersEtaSKJets_ring4_TobN", &gepIn_WTAConeBasicClustersEtaSKJetsRing4TobN);

        gt->SetBranchAddress("L1_jFexSRJetRoISim_et",  &jFexSRJSimEt);
        gt->SetBranchAddress("L1_jFexSRJetRoISim_eta", &jFexSRJSimEta);
        gt->SetBranchAddress("L1_jFexSRJetRoISim_phi", &jFexSRJSimPhi);

        gt->SetBranchAddress("L1_gFexLRJetRoISim_et",  &gFexLRJSimEt);
        gt->SetBranchAddress("L1_gFexLRJetRoISim_eta", &gFexLRJSimEta);
        gt->SetBranchAddress("L1_gFexLRJetRoISim_phi", &gFexLRJSimPhi);

        gt->SetBranchAddress("gFex_JwoJ_met",    &gepIn_gFexJwoJMet);
        gt->SetBranchAddress("gFex_JwoJ_metPhi", &gepIn_gFexJwoJMetPhi);
        gt->SetBranchAddress("gFex_JwoJ_sumEt",  &gepIn_gFexJwoJSumEt);
        gt->SetBranchAddress("gFex_NC_met",      &gepIn_gFexNCMet);
        gt->SetBranchAddress("gFex_NC_metPhi",   &gepIn_gFexNCMetPhi);
        gt->SetBranchAddress("gFex_NC_sumEt",    &gepIn_gFexNCSumEt);
        gt->SetBranchAddress("gFex_RhoRMS_met",    &gepIn_gFexRhoRMSMet);
        gt->SetBranchAddress("gFex_RhoRMS_metPhi", &gepIn_gFexRhoRMSMetPhi);
        gt->SetBranchAddress("gFex_RhoRMS_sumEt",  &gepIn_gFexRhoRMSSumEt);
        gt->SetBranchAddress("jFex_met", &gepIn_jFexMet);
        // TODO: uncomment once L1_gFexSRJetRoISim branch exists in GEP ntuples
        // gt->SetBranchAddress("L1_gFexSRJetRoISim_et",  &gFexSRJSimEt);
        // gt->SetBranchAddress("L1_gFexSRJetRoISim_eta", &gFexSRJSimEta);
        // gt->SetBranchAddress("L1_gFexSRJetRoISim_phi", &gFexSRJSimPhi);

        std::cout << "  Number of events: " << event.getEntries() << endl;

        unsigned int passedEventsCounter = 0;
        for (Long64_t iEvt = 0; iEvt < event.getEntries(); ++iEvt) { // NOTE assume that # of events for GEP and DAOD files is the same, else will get a seg fault.
            //std::cout << "iEvt: " << iEvt << "\n";
            // Retrieve truth and in time pileup jets first to apply hard-scatter-softer-than-PU (HSTP) filter (described here: https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEtMissMCSamples#Dijet_normalization_procedure_HS)
            // Also require that truth jet collection is not empty
            event.getEntry(iEvt);               // DAOD event iEvt
            gt->GetEntry(iEvt);                 // GEP event iEvt

            // Set HSTP pass flag from GEPOutputReader value (signal always passes)
            passHSTP = signalBool ? true : filterHSTP;

            gepEventNumberOut = gepEventNumberIn;
            gepRunNumberOut = gepRunNumberIn;

            const xAOD::JetContainer* AntiKt4TruthDressedWZJets = nullptr;
            if (!event.retrieve(AntiKt4TruthDressedWZJets, "AntiKt4TruthDressedWZJets").isSuccess()) {
                std::cout << "Failed to retrieve reco Antik4 Truth Dressed WZ jets" << endl;
                continue;
            }

            const xAOD::JetContainer* InTimeAntiKt4TruthJets = nullptr;
            if (!event.retrieve(InTimeAntiKt4TruthJets, "InTimeAntiKt4TruthJets").isSuccess()) {
                std::cout << "Failed to retrieve in time PU Truth jets" << endl;
                continue;
            }

            const xAOD::JetContainer* OutOfTimeAntiKt4TruthJets = nullptr;
            if (!event.retrieve(OutOfTimeAntiKt4TruthJets, "OutOfTimeAntiKt4TruthJets").isSuccess()) {
                std::cout << "Failed to retrieve out of time PU Truth jets" << endl;
                continue;
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
            outOfTimeAntiKt4TruthSRJEtValues.clear();
            outOfTimeAntiKt4TruthSRJEtIndexValues.clear();
            outOfTimeAntiKt4TruthSRJEtaValues.clear();
            outOfTimeAntiKt4TruthSRJPhiValues.clear();
            outOfTimeAntiKt4TruthSRJLeadingEtValues.clear();
            outOfTimeAntiKt4TruthSRJLeadingEtaValues.clear();
            outOfTimeAntiKt4TruthSRJLeadingPhiValues.clear();
            outOfTimeAntiKt4TruthSRJSubleadingEtValues.clear();
            outOfTimeAntiKt4TruthSRJSubleadingEtaValues.clear();
            outOfTimeAntiKt4TruthSRJSubleadingPhiValues.clear();
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
            gepCellsTowersEtaSKEtValues.clear();
            gepCellsTowersEtaSKEtaValues.clear();
            gepCellsTowersEtaSKPhiValues.clear();
            gepBasicClustersEtaSKEtValues.clear();
            gepBasicClustersEtaSKEtaValues.clear();
            gepBasicClustersEtaSKPhiValues.clear();
            gepWTAConeCellsTowersEtaSKJetspTValues.clear();
            gepWTAConeCellsTowersEtaSKJetsEtaValues.clear();
            gepWTAConeCellsTowersEtaSKJetsPhiValues.clear();
            gepWTAConeCellsTowersEtaSKJetsNConstituentsValues.clear();
            gepWTAConeBasicClustersEtaSKJetspTValues.clear();
            gepWTAConeBasicClustersEtaSKJetsEtaValues.clear();
            gepWTAConeBasicClustersEtaSKJetsPhiValues.clear();
            gepWTAConeBasicClustersEtaSKJetsNConstituentsValues.clear();
            gepLeadingWTAConeCellsTowersEtaSKJetspTValues.clear();
            gepLeadingWTAConeCellsTowersEtaSKJetsEtaValues.clear();
            gepLeadingWTAConeCellsTowersEtaSKJetsPhiValues.clear();
            gepLeadingWTAConeCellsTowersEtaSKJetsNConstituentsValues.clear();
            gepLeadingWTAConeBasicClustersEtaSKJetspTValues.clear();
            gepLeadingWTAConeBasicClustersEtaSKJetsEtaValues.clear();
            gepLeadingWTAConeBasicClustersEtaSKJetsPhiValues.clear();
            gepLeadingWTAConeBasicClustersEtaSKJetsNConstituentsValues.clear();
            gepSubleadingWTAConeCellsTowersEtaSKJetspTValues.clear();
            gepSubleadingWTAConeCellsTowersEtaSKJetsEtaValues.clear();
            gepSubleadingWTAConeCellsTowersEtaSKJetsPhiValues.clear();
            gepSubleadingWTAConeCellsTowersEtaSKJetsNConstituentsValues.clear();
            gepSubleadingWTAConeBasicClustersEtaSKJetspTValues.clear();
            gepSubleadingWTAConeBasicClustersEtaSKJetsEtaValues.clear();
            gepSubleadingWTAConeBasicClustersEtaSKJetsPhiValues.clear();
            gepSubleadingWTAConeBasicClustersEtaSKJetsNConstituentsValues.clear();
            gepWTAConeCellsTowersJetsRing0Et.clear(); gepWTAConeCellsTowersJetsRing1Et.clear(); gepWTAConeCellsTowersJetsRing2Et.clear(); gepWTAConeCellsTowersJetsRing3Et.clear(); gepWTAConeCellsTowersJetsRing4Et.clear();
            gepWTAConeCellsTowersJetsTotalTobN.clear(); gepWTAConeCellsTowersJetsRing0TobN.clear(); gepWTAConeCellsTowersJetsRing1TobN.clear(); gepWTAConeCellsTowersJetsRing2TobN.clear(); gepWTAConeCellsTowersJetsRing3TobN.clear(); gepWTAConeCellsTowersJetsRing4TobN.clear();
            gepLeadingWTAConeCellsTowersJetsRing0Et.clear(); gepLeadingWTAConeCellsTowersJetsRing1Et.clear(); gepLeadingWTAConeCellsTowersJetsRing2Et.clear(); gepLeadingWTAConeCellsTowersJetsRing3Et.clear(); gepLeadingWTAConeCellsTowersJetsRing4Et.clear();
            gepLeadingWTAConeCellsTowersJetsTotalTobN.clear(); gepLeadingWTAConeCellsTowersJetsRing0TobN.clear(); gepLeadingWTAConeCellsTowersJetsRing1TobN.clear(); gepLeadingWTAConeCellsTowersJetsRing2TobN.clear(); gepLeadingWTAConeCellsTowersJetsRing3TobN.clear(); gepLeadingWTAConeCellsTowersJetsRing4TobN.clear();
            gepSubleadingWTAConeCellsTowersJetsRing0Et.clear(); gepSubleadingWTAConeCellsTowersJetsRing1Et.clear(); gepSubleadingWTAConeCellsTowersJetsRing2Et.clear(); gepSubleadingWTAConeCellsTowersJetsRing3Et.clear(); gepSubleadingWTAConeCellsTowersJetsRing4Et.clear();
            gepSubleadingWTAConeCellsTowersJetsTotalTobN.clear(); gepSubleadingWTAConeCellsTowersJetsRing0TobN.clear(); gepSubleadingWTAConeCellsTowersJetsRing1TobN.clear(); gepSubleadingWTAConeCellsTowersJetsRing2TobN.clear(); gepSubleadingWTAConeCellsTowersJetsRing3TobN.clear(); gepSubleadingWTAConeCellsTowersJetsRing4TobN.clear();
            gepWTAConeCellsTowersSKJetsRing0Et.clear(); gepWTAConeCellsTowersSKJetsRing1Et.clear(); gepWTAConeCellsTowersSKJetsRing2Et.clear(); gepWTAConeCellsTowersSKJetsRing3Et.clear(); gepWTAConeCellsTowersSKJetsRing4Et.clear();
            gepWTAConeCellsTowersSKJetsTotalTobN.clear(); gepWTAConeCellsTowersSKJetsRing0TobN.clear(); gepWTAConeCellsTowersSKJetsRing1TobN.clear(); gepWTAConeCellsTowersSKJetsRing2TobN.clear(); gepWTAConeCellsTowersSKJetsRing3TobN.clear(); gepWTAConeCellsTowersSKJetsRing4TobN.clear();
            gepLeadingWTAConeCellsTowersSKJetsRing0Et.clear(); gepLeadingWTAConeCellsTowersSKJetsRing1Et.clear(); gepLeadingWTAConeCellsTowersSKJetsRing2Et.clear(); gepLeadingWTAConeCellsTowersSKJetsRing3Et.clear(); gepLeadingWTAConeCellsTowersSKJetsRing4Et.clear();
            gepLeadingWTAConeCellsTowersSKJetsTotalTobN.clear(); gepLeadingWTAConeCellsTowersSKJetsRing0TobN.clear(); gepLeadingWTAConeCellsTowersSKJetsRing1TobN.clear(); gepLeadingWTAConeCellsTowersSKJetsRing2TobN.clear(); gepLeadingWTAConeCellsTowersSKJetsRing3TobN.clear(); gepLeadingWTAConeCellsTowersSKJetsRing4TobN.clear();
            gepSubleadingWTAConeCellsTowersSKJetsRing0Et.clear(); gepSubleadingWTAConeCellsTowersSKJetsRing1Et.clear(); gepSubleadingWTAConeCellsTowersSKJetsRing2Et.clear(); gepSubleadingWTAConeCellsTowersSKJetsRing3Et.clear(); gepSubleadingWTAConeCellsTowersSKJetsRing4Et.clear();
            gepSubleadingWTAConeCellsTowersSKJetsTotalTobN.clear(); gepSubleadingWTAConeCellsTowersSKJetsRing0TobN.clear(); gepSubleadingWTAConeCellsTowersSKJetsRing1TobN.clear(); gepSubleadingWTAConeCellsTowersSKJetsRing2TobN.clear(); gepSubleadingWTAConeCellsTowersSKJetsRing3TobN.clear(); gepSubleadingWTAConeCellsTowersSKJetsRing4TobN.clear();
            gepWTAConeCellsTowersEtaSKJetsRing0Et.clear(); gepWTAConeCellsTowersEtaSKJetsRing1Et.clear(); gepWTAConeCellsTowersEtaSKJetsRing2Et.clear(); gepWTAConeCellsTowersEtaSKJetsRing3Et.clear(); gepWTAConeCellsTowersEtaSKJetsRing4Et.clear();
            gepWTAConeCellsTowersEtaSKJetsTotalTobN.clear(); gepWTAConeCellsTowersEtaSKJetsRing0TobN.clear(); gepWTAConeCellsTowersEtaSKJetsRing1TobN.clear(); gepWTAConeCellsTowersEtaSKJetsRing2TobN.clear(); gepWTAConeCellsTowersEtaSKJetsRing3TobN.clear(); gepWTAConeCellsTowersEtaSKJetsRing4TobN.clear();
            gepLeadingWTAConeCellsTowersEtaSKJetsRing0Et.clear(); gepLeadingWTAConeCellsTowersEtaSKJetsRing1Et.clear(); gepLeadingWTAConeCellsTowersEtaSKJetsRing2Et.clear(); gepLeadingWTAConeCellsTowersEtaSKJetsRing3Et.clear(); gepLeadingWTAConeCellsTowersEtaSKJetsRing4Et.clear();
            gepLeadingWTAConeCellsTowersEtaSKJetsTotalTobN.clear(); gepLeadingWTAConeCellsTowersEtaSKJetsRing0TobN.clear(); gepLeadingWTAConeCellsTowersEtaSKJetsRing1TobN.clear(); gepLeadingWTAConeCellsTowersEtaSKJetsRing2TobN.clear(); gepLeadingWTAConeCellsTowersEtaSKJetsRing3TobN.clear(); gepLeadingWTAConeCellsTowersEtaSKJetsRing4TobN.clear();
            gepSubleadingWTAConeCellsTowersEtaSKJetsRing0Et.clear(); gepSubleadingWTAConeCellsTowersEtaSKJetsRing1Et.clear(); gepSubleadingWTAConeCellsTowersEtaSKJetsRing2Et.clear(); gepSubleadingWTAConeCellsTowersEtaSKJetsRing3Et.clear(); gepSubleadingWTAConeCellsTowersEtaSKJetsRing4Et.clear();
            gepSubleadingWTAConeCellsTowersEtaSKJetsTotalTobN.clear(); gepSubleadingWTAConeCellsTowersEtaSKJetsRing0TobN.clear(); gepSubleadingWTAConeCellsTowersEtaSKJetsRing1TobN.clear(); gepSubleadingWTAConeCellsTowersEtaSKJetsRing2TobN.clear(); gepSubleadingWTAConeCellsTowersEtaSKJetsRing3TobN.clear(); gepSubleadingWTAConeCellsTowersEtaSKJetsRing4TobN.clear();
            gepWTAConeBasicClustersJetsRing0Et.clear(); gepWTAConeBasicClustersJetsRing1Et.clear(); gepWTAConeBasicClustersJetsRing2Et.clear(); gepWTAConeBasicClustersJetsRing3Et.clear(); gepWTAConeBasicClustersJetsRing4Et.clear();
            gepWTAConeBasicClustersJetsTotalTobN.clear(); gepWTAConeBasicClustersJetsRing0TobN.clear(); gepWTAConeBasicClustersJetsRing1TobN.clear(); gepWTAConeBasicClustersJetsRing2TobN.clear(); gepWTAConeBasicClustersJetsRing3TobN.clear(); gepWTAConeBasicClustersJetsRing4TobN.clear();
            gepLeadingWTAConeBasicClustersJetsRing0Et.clear(); gepLeadingWTAConeBasicClustersJetsRing1Et.clear(); gepLeadingWTAConeBasicClustersJetsRing2Et.clear(); gepLeadingWTAConeBasicClustersJetsRing3Et.clear(); gepLeadingWTAConeBasicClustersJetsRing4Et.clear();
            gepLeadingWTAConeBasicClustersJetsTotalTobN.clear(); gepLeadingWTAConeBasicClustersJetsRing0TobN.clear(); gepLeadingWTAConeBasicClustersJetsRing1TobN.clear(); gepLeadingWTAConeBasicClustersJetsRing2TobN.clear(); gepLeadingWTAConeBasicClustersJetsRing3TobN.clear(); gepLeadingWTAConeBasicClustersJetsRing4TobN.clear();
            gepSubleadingWTAConeBasicClustersJetsRing0Et.clear(); gepSubleadingWTAConeBasicClustersJetsRing1Et.clear(); gepSubleadingWTAConeBasicClustersJetsRing2Et.clear(); gepSubleadingWTAConeBasicClustersJetsRing3Et.clear(); gepSubleadingWTAConeBasicClustersJetsRing4Et.clear();
            gepSubleadingWTAConeBasicClustersJetsTotalTobN.clear(); gepSubleadingWTAConeBasicClustersJetsRing0TobN.clear(); gepSubleadingWTAConeBasicClustersJetsRing1TobN.clear(); gepSubleadingWTAConeBasicClustersJetsRing2TobN.clear(); gepSubleadingWTAConeBasicClustersJetsRing3TobN.clear(); gepSubleadingWTAConeBasicClustersJetsRing4TobN.clear();
            gepWTAConeBasicClustersSKJetsRing0Et.clear(); gepWTAConeBasicClustersSKJetsRing1Et.clear(); gepWTAConeBasicClustersSKJetsRing2Et.clear(); gepWTAConeBasicClustersSKJetsRing3Et.clear(); gepWTAConeBasicClustersSKJetsRing4Et.clear();
            gepWTAConeBasicClustersSKJetsTotalTobN.clear(); gepWTAConeBasicClustersSKJetsRing0TobN.clear(); gepWTAConeBasicClustersSKJetsRing1TobN.clear(); gepWTAConeBasicClustersSKJetsRing2TobN.clear(); gepWTAConeBasicClustersSKJetsRing3TobN.clear(); gepWTAConeBasicClustersSKJetsRing4TobN.clear();
            gepLeadingWTAConeBasicClustersSKJetsRing0Et.clear(); gepLeadingWTAConeBasicClustersSKJetsRing1Et.clear(); gepLeadingWTAConeBasicClustersSKJetsRing2Et.clear(); gepLeadingWTAConeBasicClustersSKJetsRing3Et.clear(); gepLeadingWTAConeBasicClustersSKJetsRing4Et.clear();
            gepLeadingWTAConeBasicClustersSKJetsTotalTobN.clear(); gepLeadingWTAConeBasicClustersSKJetsRing0TobN.clear(); gepLeadingWTAConeBasicClustersSKJetsRing1TobN.clear(); gepLeadingWTAConeBasicClustersSKJetsRing2TobN.clear(); gepLeadingWTAConeBasicClustersSKJetsRing3TobN.clear(); gepLeadingWTAConeBasicClustersSKJetsRing4TobN.clear();
            gepSubleadingWTAConeBasicClustersSKJetsRing0Et.clear(); gepSubleadingWTAConeBasicClustersSKJetsRing1Et.clear(); gepSubleadingWTAConeBasicClustersSKJetsRing2Et.clear(); gepSubleadingWTAConeBasicClustersSKJetsRing3Et.clear(); gepSubleadingWTAConeBasicClustersSKJetsRing4Et.clear();
            gepSubleadingWTAConeBasicClustersSKJetsTotalTobN.clear(); gepSubleadingWTAConeBasicClustersSKJetsRing0TobN.clear(); gepSubleadingWTAConeBasicClustersSKJetsRing1TobN.clear(); gepSubleadingWTAConeBasicClustersSKJetsRing2TobN.clear(); gepSubleadingWTAConeBasicClustersSKJetsRing3TobN.clear(); gepSubleadingWTAConeBasicClustersSKJetsRing4TobN.clear();
            gepWTAConeBasicClustersEtaSKJetsRing0Et.clear(); gepWTAConeBasicClustersEtaSKJetsRing1Et.clear(); gepWTAConeBasicClustersEtaSKJetsRing2Et.clear(); gepWTAConeBasicClustersEtaSKJetsRing3Et.clear(); gepWTAConeBasicClustersEtaSKJetsRing4Et.clear();
            gepWTAConeBasicClustersEtaSKJetsTotalTobN.clear(); gepWTAConeBasicClustersEtaSKJetsRing0TobN.clear(); gepWTAConeBasicClustersEtaSKJetsRing1TobN.clear(); gepWTAConeBasicClustersEtaSKJetsRing2TobN.clear(); gepWTAConeBasicClustersEtaSKJetsRing3TobN.clear(); gepWTAConeBasicClustersEtaSKJetsRing4TobN.clear();
            gepLeadingWTAConeBasicClustersEtaSKJetsRing0Et.clear(); gepLeadingWTAConeBasicClustersEtaSKJetsRing1Et.clear(); gepLeadingWTAConeBasicClustersEtaSKJetsRing2Et.clear(); gepLeadingWTAConeBasicClustersEtaSKJetsRing3Et.clear(); gepLeadingWTAConeBasicClustersEtaSKJetsRing4Et.clear();
            gepLeadingWTAConeBasicClustersEtaSKJetsTotalTobN.clear(); gepLeadingWTAConeBasicClustersEtaSKJetsRing0TobN.clear(); gepLeadingWTAConeBasicClustersEtaSKJetsRing1TobN.clear(); gepLeadingWTAConeBasicClustersEtaSKJetsRing2TobN.clear(); gepLeadingWTAConeBasicClustersEtaSKJetsRing3TobN.clear(); gepLeadingWTAConeBasicClustersEtaSKJetsRing4TobN.clear();
            gepSubleadingWTAConeBasicClustersEtaSKJetsRing0Et.clear(); gepSubleadingWTAConeBasicClustersEtaSKJetsRing1Et.clear(); gepSubleadingWTAConeBasicClustersEtaSKJetsRing2Et.clear(); gepSubleadingWTAConeBasicClustersEtaSKJetsRing3Et.clear(); gepSubleadingWTAConeBasicClustersEtaSKJetsRing4Et.clear();
            gepSubleadingWTAConeBasicClustersEtaSKJetsTotalTobN.clear(); gepSubleadingWTAConeBasicClustersEtaSKJetsRing0TobN.clear(); gepSubleadingWTAConeBasicClustersEtaSKJetsRing1TobN.clear(); gepSubleadingWTAConeBasicClustersEtaSKJetsRing2TobN.clear(); gepSubleadingWTAConeBasicClustersEtaSKJetsRing3TobN.clear(); gepSubleadingWTAConeBasicClustersEtaSKJetsRing4TobN.clear();
            jFexSRJSimEtValues.clear(); jFexSRJSimEtaValues.clear(); jFexSRJSimPhiValues.clear();
            jFexSRJSimLeadingEtValues.clear(); jFexSRJSimLeadingEtaValues.clear(); jFexSRJSimLeadingPhiValues.clear();
            jFexSRJSimSubleadingEtValues.clear(); jFexSRJSimSubleadingEtaValues.clear(); jFexSRJSimSubleadingPhiValues.clear();
            gFexLRJSimEtValues.clear(); gFexLRJSimEtaValues.clear(); gFexLRJSimPhiValues.clear();
            gFexLRJSimLeadingEtValues.clear(); gFexLRJSimLeadingEtaValues.clear(); gFexLRJSimLeadingPhiValues.clear();
            gFexLRJSimSubleadingEtValues.clear(); gFexLRJSimSubleadingEtaValues.clear(); gFexLRJSimSubleadingPhiValues.clear();
            // TODO: uncomment once L1_gFexSRJetRoISim branch exists in GEP ntuples
            // gFexSRJSimEtValues.clear(); gFexSRJSimEtaValues.clear(); gFexSRJSimPhiValues.clear();
            // gFexSRJSimLeadingEtValues.clear(); gFexSRJSimLeadingEtaValues.clear(); gFexSRJSimLeadingPhiValues.clear();
            // gFexSRJSimSubleadingEtValues.clear(); gFexSRJSimSubleadingEtaValues.clear(); gFexSRJSimSubleadingPhiValues.clear();
            // gepRunNumber / gepEventNumber are NOT cleared here — they are set each event by gt->GetEntry(iEvt)
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
            for (unsigned int i = 0; i < gepIn_BasicClustersEt->size(); ++i) {
                float gepBasicClusterEt  = (*gepIn_BasicClustersEt)[i] / 1000.0;
                float gepBasicClusterEta = (*gepIn_BasicClustersEta)[i];
                float gepBasicClusterPhi = (*gepIn_BasicClustersPhi)[i];

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
                                    
                if (iEvt < maxDatEvents) {
                    if (gepbasicclusters_it == 0) {
                        f_gepbasicclusters << "Event : " << std::dec << iEvt << "\n";
                    }
                    // 4. Write to output file
                    f_gepbasicclusters << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepbasicclusters_it
                                << " " << binary_word
                                << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    gepbasicclusters_it++;
                }

            }

            // SK basic clusters
            unsigned int gepbasicclusterssk_it = 0;
            for (unsigned int i = 0; i < gepIn_BasicClustersSKEt->size(); ++i) {
                float gepBasicClusterSKEt  = (*gepIn_BasicClustersSKEt)[i] / 1000.0;
                float gepBasicClusterSKEta = (*gepIn_BasicClustersSKEta)[i];
                float gepBasicClusterSKPhi = (*gepIn_BasicClustersSKPhi)[i];

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
                                    
                if (iEvt < maxDatEvents) {
                    if (gepbasicclusterssk_it == 0) {
                        f_gepbasicclusterssk << "Event : " << std::dec << iEvt << "\n";
                    }
                    // 4. Write to output file
                    f_gepbasicclusterssk << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepbasicclusterssk_it
                                << " " << binary_word
                                << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    gepbasicclusterssk_it++;
                }
            }

            // ---------- gepCellsTowers ----------
            {
                // Build index order sorted by Et (descending)
                std::vector<unsigned int> order(gepIn_CellsTowersEt->size());
                std::iota(order.begin(), order.end(), 0u);
                std::sort(order.begin(), order.end(),
                        [&](unsigned int a, unsigned int b){
                            return (*gepIn_CellsTowersEt)[a] > (*gepIn_CellsTowersEt)[b];
                        });

                unsigned int gepcellstowers_it = 0;
                for (unsigned int i : order) {
                    float gepCellsTowersEtValue  = (*gepIn_CellsTowersEt)[i] / 1000.0f;
                    float gepCellsTowersEtaValue = (*gepIn_CellsTowersEta)[i];
                    float gepCellsTowersPhiValue = (*gepIn_CellsTowersPhi)[i];

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

                        if (iEvt < maxDatEvents) {
                            if (gepcellstowers_it == 0) {
                                f_gepcellstowers_sort << "Event : " << std::dec << iEvt << "\n";
                            }
                            f_gepcellstowers_sort << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepcellstowers_it
                                            << " "  << binary_word
                                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                        }
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

                    if (iEvt < maxDatEvents) {
                        if (gepcellstowers_it == 0) {
                            f_gepcellstowers << "Event : " << std::dec << iEvt << "\n";
                        }
                        f_gepcellstowers << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepcellstowers_it
                                        << " "  << binary_word
                                        << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    }
                    ++gepcellstowers_it;

                }
            }

            // ---------- gepCellsTowers SK PU Suppression applied ----------
            {
                // Build index order sorted by Et (descending)
                std::vector<unsigned int> order(gepIn_CellsTowersSKEt->size());
                std::iota(order.begin(), order.end(), 0u);
                std::sort(order.begin(), order.end(),
                        [&](unsigned int a, unsigned int b){
                            return (*gepIn_CellsTowersSKEt)[a] > (*gepIn_CellsTowersSKEt)[b];
                        });

                unsigned int gepcellstowerssk_it = 0;
                for (unsigned int i : order) {
                    float gepCellsTowersSKEtValue  = (*gepIn_CellsTowersSKEt)[i] / 1000.0f;
                    float gepCellsTowersSKEtaValue = (*gepIn_CellsTowersSKEta)[i];
                    float gepCellsTowersSKPhiValue = (*gepIn_CellsTowersSKPhi)[i];

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

                    if (iEvt < maxDatEvents) {
                        if (gepcellstowerssk_it == 0) {
                            f_gepcellstowerssk << "Event : " << std::dec << iEvt << "\n";
                        }
                        f_gepcellstowerssk << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepcellstowerssk_it
                                        << " "  << binary_word
                                        << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    }
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
                                    
                if (iEvt < maxDatEvents) {
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
                std::vector<unsigned int> indices(gepIn_WTAConeCellsTowersJetsPt->size());
                std::iota(indices.begin(), indices.end(), 0);

                std::sort(indices.begin(), indices.end(),
                        [&](unsigned int a, unsigned int b) {
                            return (*gepIn_WTAConeCellsTowersJetsPt)[a] > (*gepIn_WTAConeCellsTowersJetsPt)[b];
                        });
                unsigned int gepwtaconecellstowers_it = 0;
                for (unsigned int i = 0; i < indices.size(); i++) {
                    unsigned int idx = indices[i];
                    float WTAConeCellsTowersJetspT  = (*gepIn_WTAConeCellsTowersJetsPt)[idx] / 1000.0;
                    float WTAConeCellsTowersJetsEta = (*gepIn_WTAConeCellsTowersJetsEta)[idx];
                    float WTAConeCellsTowersJetsPhi = (*gepIn_WTAConeCellsTowersJetsPhi)[idx];
                    unsigned int WTAConeCellsTowersJetsNConstituents = (*gepIn_WTAConeCellsTowersJetsNConstituents)[idx];
                    float ct_r0Et  = (*gepIn_WTAConeCellsTowersJetsRing0Et)[idx];
                    float ct_r1Et  = (*gepIn_WTAConeCellsTowersJetsRing1Et)[idx];
                    float ct_r2Et  = (*gepIn_WTAConeCellsTowersJetsRing2Et)[idx];
                    float ct_r3Et  = (*gepIn_WTAConeCellsTowersJetsRing3Et)[idx];
                    float ct_r4Et  = (*gepIn_WTAConeCellsTowersJetsRing4Et)[idx];
                    int   ct_totN  = (*gepIn_WTAConeCellsTowersJetsTotalTobN)[idx];
                    int   ct_r0N   = (*gepIn_WTAConeCellsTowersJetsRing0TobN)[idx];
                    int   ct_r1N   = (*gepIn_WTAConeCellsTowersJetsRing1TobN)[idx];
                    int   ct_r2N   = (*gepIn_WTAConeCellsTowersJetsRing2TobN)[idx];
                    int   ct_r3N   = (*gepIn_WTAConeCellsTowersJetsRing3TobN)[idx];
                    int   ct_r4N   = (*gepIn_WTAConeCellsTowersJetsRing4TobN)[idx];

                    if(i == 0){
                        gepLeadingWTAConeCellsTowersJetspTValues.push_back(WTAConeCellsTowersJetspT);
                        gepLeadingWTAConeCellsTowersJetsEtaValues.push_back(WTAConeCellsTowersJetsEta);
                        gepLeadingWTAConeCellsTowersJetsPhiValues.push_back(WTAConeCellsTowersJetsPhi);
                        gepLeadingWTAConeCellsTowersJetsNConstituentsValues.push_back(WTAConeCellsTowersJetsNConstituents);
                        gepLeadingWTAConeCellsTowersJetsRing0Et.push_back(ct_r0Et); gepLeadingWTAConeCellsTowersJetsRing1Et.push_back(ct_r1Et); gepLeadingWTAConeCellsTowersJetsRing2Et.push_back(ct_r2Et); gepLeadingWTAConeCellsTowersJetsRing3Et.push_back(ct_r3Et); gepLeadingWTAConeCellsTowersJetsRing4Et.push_back(ct_r4Et);
                        gepLeadingWTAConeCellsTowersJetsTotalTobN.push_back(ct_totN); gepLeadingWTAConeCellsTowersJetsRing0TobN.push_back(ct_r0N); gepLeadingWTAConeCellsTowersJetsRing1TobN.push_back(ct_r1N); gepLeadingWTAConeCellsTowersJetsRing2TobN.push_back(ct_r2N); gepLeadingWTAConeCellsTowersJetsRing3TobN.push_back(ct_r3N); gepLeadingWTAConeCellsTowersJetsRing4TobN.push_back(ct_r4N);
                    }
                    if (i == 1){
                        gepSubleadingWTAConeCellsTowersJetspTValues.push_back(WTAConeCellsTowersJetspT);
                        gepSubleadingWTAConeCellsTowersJetsEtaValues.push_back(WTAConeCellsTowersJetsEta);
                        gepSubleadingWTAConeCellsTowersJetsPhiValues.push_back(WTAConeCellsTowersJetsPhi);
                        gepSubleadingWTAConeCellsTowersJetsNConstituentsValues.push_back(WTAConeCellsTowersJetsNConstituents);
                        gepSubleadingWTAConeCellsTowersJetsRing0Et.push_back(ct_r0Et); gepSubleadingWTAConeCellsTowersJetsRing1Et.push_back(ct_r1Et); gepSubleadingWTAConeCellsTowersJetsRing2Et.push_back(ct_r2Et); gepSubleadingWTAConeCellsTowersJetsRing3Et.push_back(ct_r3Et); gepSubleadingWTAConeCellsTowersJetsRing4Et.push_back(ct_r4Et);
                        gepSubleadingWTAConeCellsTowersJetsTotalTobN.push_back(ct_totN); gepSubleadingWTAConeCellsTowersJetsRing0TobN.push_back(ct_r0N); gepSubleadingWTAConeCellsTowersJetsRing1TobN.push_back(ct_r1N); gepSubleadingWTAConeCellsTowersJetsRing2TobN.push_back(ct_r2N); gepSubleadingWTAConeCellsTowersJetsRing3TobN.push_back(ct_r3N); gepSubleadingWTAConeCellsTowersJetsRing4TobN.push_back(ct_r4N);
                    }

                    gepWTAConeCellsTowersJetspTValues.push_back(WTAConeCellsTowersJetspT);
                    gepWTAConeCellsTowersJetsEtaValues.push_back(WTAConeCellsTowersJetsEta);
                    gepWTAConeCellsTowersJetsPhiValues.push_back(WTAConeCellsTowersJetsPhi);
                    gepWTAConeCellsTowersJetsNConstituentsValues.push_back(WTAConeCellsTowersJetsNConstituents);
                    gepWTAConeCellsTowersJetsRing0Et.push_back(ct_r0Et); gepWTAConeCellsTowersJetsRing1Et.push_back(ct_r1Et); gepWTAConeCellsTowersJetsRing2Et.push_back(ct_r2Et); gepWTAConeCellsTowersJetsRing3Et.push_back(ct_r3Et); gepWTAConeCellsTowersJetsRing4Et.push_back(ct_r4Et);
                    gepWTAConeCellsTowersJetsTotalTobN.push_back(ct_totN); gepWTAConeCellsTowersJetsRing0TobN.push_back(ct_r0N); gepWTAConeCellsTowersJetsRing1TobN.push_back(ct_r1N); gepWTAConeCellsTowersJetsRing2TobN.push_back(ct_r2N); gepWTAConeCellsTowersJetsRing3TobN.push_back(ct_r3N); gepWTAConeCellsTowersJetsRing4TobN.push_back(ct_r4N);

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

                    if (iEvt < maxDatEvents) {
                        if (gepwtaconecellstowers_it == 0) {
                            f_wtaconejetscellstowers << "Event : " << std::dec << iEvt << "\n";
                        }
                        f_wtaconejetscellstowers << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepwtaconecellstowers_it
                                            << " "  << binary_word
                                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    }
                    ++gepwtaconecellstowers_it;

                }
            }

            // GEP cone jets from GEPCellsTowers with SoftKiller PU Suppression applied
            //std::cout << "gepIn_WTAConeCellsTowersSKJetsPt->size(): " << gepIn_WTAConeCellsTowersSKJetsPt->size() << "\n";
            if(gepIn_WTAConeCellsTowersSKJetsPt->size() == 0) std::cout << "NO WTA CONE JETS FOR THIS EVENT!" << "\n";
            {
                std::vector<unsigned int> indices(gepIn_WTAConeCellsTowersSKJetsPt->size());
                std::iota(indices.begin(), indices.end(), 0);

                std::sort(indices.begin(), indices.end(),
                        [&](unsigned int a, unsigned int b) {
                            return (*gepIn_WTAConeCellsTowersSKJetsPt)[a] > (*gepIn_WTAConeCellsTowersSKJetsPt)[b];
                        });
                unsigned int gepwtaconecellstowerssk_it = 0;
                for (unsigned int i = 0; i < indices.size(); i++) {
                    unsigned int idx = indices[i];
                    float WTAConeCellsTowersSKJetspT  = (*gepIn_WTAConeCellsTowersSKJetsPt)[idx] / 1000.0;
                    float WTAConeCellsTowersSKJetsEta = (*gepIn_WTAConeCellsTowersSKJetsEta)[idx];
                    float WTAConeCellsTowersSKJetsPhi = (*gepIn_WTAConeCellsTowersSKJetsPhi)[idx];
                    unsigned int WTAConeCellsTowersSKJetsNConstituents = (*gepIn_WTAConeCellsTowersSKJetsNConstituents)[idx];
                    float ctsk_r0Et = (*gepIn_WTAConeCellsTowersSKJetsRing0Et)[idx];
                    float ctsk_r1Et = (*gepIn_WTAConeCellsTowersSKJetsRing1Et)[idx];
                    float ctsk_r2Et = (*gepIn_WTAConeCellsTowersSKJetsRing2Et)[idx];
                    float ctsk_r3Et = (*gepIn_WTAConeCellsTowersSKJetsRing3Et)[idx];
                    float ctsk_r4Et = (*gepIn_WTAConeCellsTowersSKJetsRing4Et)[idx];
                    int   ctsk_totN = (*gepIn_WTAConeCellsTowersSKJetsTotalTobN)[idx];
                    int   ctsk_r0N  = (*gepIn_WTAConeCellsTowersSKJetsRing0TobN)[idx];
                    int   ctsk_r1N  = (*gepIn_WTAConeCellsTowersSKJetsRing1TobN)[idx];
                    int   ctsk_r2N  = (*gepIn_WTAConeCellsTowersSKJetsRing2TobN)[idx];
                    int   ctsk_r3N  = (*gepIn_WTAConeCellsTowersSKJetsRing3TobN)[idx];
                    int   ctsk_r4N  = (*gepIn_WTAConeCellsTowersSKJetsRing4TobN)[idx];

                    if(i == 0){
                        gepLeadingWTAConeCellsTowersSKJetspTValues.push_back(WTAConeCellsTowersSKJetspT);
                        gepLeadingWTAConeCellsTowersSKJetsEtaValues.push_back(WTAConeCellsTowersSKJetsEta);
                        gepLeadingWTAConeCellsTowersSKJetsPhiValues.push_back(WTAConeCellsTowersSKJetsPhi);
                        gepLeadingWTAConeCellsTowersSKJetsNConstituentsValues.push_back(WTAConeCellsTowersSKJetsNConstituents);
                        gepLeadingWTAConeCellsTowersSKJetsRing0Et.push_back(ctsk_r0Et); gepLeadingWTAConeCellsTowersSKJetsRing1Et.push_back(ctsk_r1Et); gepLeadingWTAConeCellsTowersSKJetsRing2Et.push_back(ctsk_r2Et); gepLeadingWTAConeCellsTowersSKJetsRing3Et.push_back(ctsk_r3Et); gepLeadingWTAConeCellsTowersSKJetsRing4Et.push_back(ctsk_r4Et);
                        gepLeadingWTAConeCellsTowersSKJetsTotalTobN.push_back(ctsk_totN); gepLeadingWTAConeCellsTowersSKJetsRing0TobN.push_back(ctsk_r0N); gepLeadingWTAConeCellsTowersSKJetsRing1TobN.push_back(ctsk_r1N); gepLeadingWTAConeCellsTowersSKJetsRing2TobN.push_back(ctsk_r2N); gepLeadingWTAConeCellsTowersSKJetsRing3TobN.push_back(ctsk_r3N); gepLeadingWTAConeCellsTowersSKJetsRing4TobN.push_back(ctsk_r4N);
                    }
                    if (i == 1){
                        gepSubleadingWTAConeCellsTowersSKJetspTValues.push_back(WTAConeCellsTowersSKJetspT);
                        gepSubleadingWTAConeCellsTowersSKJetsEtaValues.push_back(WTAConeCellsTowersSKJetsEta);
                        gepSubleadingWTAConeCellsTowersSKJetsPhiValues.push_back(WTAConeCellsTowersSKJetsPhi);
                        gepSubleadingWTAConeCellsTowersSKJetsNConstituentsValues.push_back(WTAConeCellsTowersSKJetsNConstituents);
                        gepSubleadingWTAConeCellsTowersSKJetsRing0Et.push_back(ctsk_r0Et); gepSubleadingWTAConeCellsTowersSKJetsRing1Et.push_back(ctsk_r1Et); gepSubleadingWTAConeCellsTowersSKJetsRing2Et.push_back(ctsk_r2Et); gepSubleadingWTAConeCellsTowersSKJetsRing3Et.push_back(ctsk_r3Et); gepSubleadingWTAConeCellsTowersSKJetsRing4Et.push_back(ctsk_r4Et);
                        gepSubleadingWTAConeCellsTowersSKJetsTotalTobN.push_back(ctsk_totN); gepSubleadingWTAConeCellsTowersSKJetsRing0TobN.push_back(ctsk_r0N); gepSubleadingWTAConeCellsTowersSKJetsRing1TobN.push_back(ctsk_r1N); gepSubleadingWTAConeCellsTowersSKJetsRing2TobN.push_back(ctsk_r2N); gepSubleadingWTAConeCellsTowersSKJetsRing3TobN.push_back(ctsk_r3N); gepSubleadingWTAConeCellsTowersSKJetsRing4TobN.push_back(ctsk_r4N);
                    }

                    gepWTAConeCellsTowersSKJetspTValues.push_back(WTAConeCellsTowersSKJetspT);
                    gepWTAConeCellsTowersSKJetsEtaValues.push_back(WTAConeCellsTowersSKJetsEta);
                    gepWTAConeCellsTowersSKJetsPhiValues.push_back(WTAConeCellsTowersSKJetsPhi);
                    gepWTAConeCellsTowersSKJetsNConstituentsValues.push_back(WTAConeCellsTowersSKJetsNConstituents);
                    gepWTAConeCellsTowersSKJetsRing0Et.push_back(ctsk_r0Et); gepWTAConeCellsTowersSKJetsRing1Et.push_back(ctsk_r1Et); gepWTAConeCellsTowersSKJetsRing2Et.push_back(ctsk_r2Et); gepWTAConeCellsTowersSKJetsRing3Et.push_back(ctsk_r3Et); gepWTAConeCellsTowersSKJetsRing4Et.push_back(ctsk_r4Et);
                    gepWTAConeCellsTowersSKJetsTotalTobN.push_back(ctsk_totN); gepWTAConeCellsTowersSKJetsRing0TobN.push_back(ctsk_r0N); gepWTAConeCellsTowersSKJetsRing1TobN.push_back(ctsk_r1N); gepWTAConeCellsTowersSKJetsRing2TobN.push_back(ctsk_r2N); gepWTAConeCellsTowersSKJetsRing3TobN.push_back(ctsk_r3N); gepWTAConeCellsTowersSKJetsRing4TobN.push_back(ctsk_r4N);

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


                    if (iEvt < maxDatEvents) {
                        if (gepwtaconecellstowerssk_it == 0) {
                            f_wtaconejetscellstowerssk << "Event : " << std::dec << iEvt << "\n";
                        }
                        f_wtaconejetscellstowerssk << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepwtaconecellstowerssk_it
                                            << " "  << binary_word
                                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    }
                    ++gepwtaconecellstowerssk_it;

                }
            }

            // GEP cone jets from basic clusters
            {
                std::vector<unsigned int> indices(gepIn_WTAConeBasicClustersJetsPt->size());
                std::iota(indices.begin(), indices.end(), 0);

                std::sort(indices.begin(), indices.end(),
                        [&](unsigned int a, unsigned int b) {
                            return (*gepIn_WTAConeBasicClustersJetsPt)[a] > (*gepIn_WTAConeBasicClustersJetsPt)[b];
                        });
                unsigned int gepwtaconebasicclusters_it = 0;
                for (unsigned int i = 0; i < indices.size(); i++) {
                    unsigned int idx = indices[i];
                    float coneWTAGEPBasicClustersJetspT  = (*gepIn_WTAConeBasicClustersJetsPt)[idx] / 1000.0;
                    float coneWTAGEPBasicClustersJetsEta = (*gepIn_WTAConeBasicClustersJetsEta)[idx];
                    float coneWTAGEPBasicClustersJetsPhi = (*gepIn_WTAConeBasicClustersJetsPhi)[idx];
                    unsigned int coneWTAGEPBasicClustersJetsNConstituents = (*gepIn_WTAConeBasicClustersJetsNConstituents)[idx];
                    float bc_r0Et  = (*gepIn_WTAConeBasicClustersJetsRing0Et)[idx];
                    float bc_r1Et  = (*gepIn_WTAConeBasicClustersJetsRing1Et)[idx];
                    float bc_r2Et  = (*gepIn_WTAConeBasicClustersJetsRing2Et)[idx];
                    float bc_r3Et  = (*gepIn_WTAConeBasicClustersJetsRing3Et)[idx];
                    float bc_r4Et  = (*gepIn_WTAConeBasicClustersJetsRing4Et)[idx];
                    int   bc_totN  = (*gepIn_WTAConeBasicClustersJetsTotalTobN)[idx];
                    int   bc_r0N   = (*gepIn_WTAConeBasicClustersJetsRing0TobN)[idx];
                    int   bc_r1N   = (*gepIn_WTAConeBasicClustersJetsRing1TobN)[idx];
                    int   bc_r2N   = (*gepIn_WTAConeBasicClustersJetsRing2TobN)[idx];
                    int   bc_r3N   = (*gepIn_WTAConeBasicClustersJetsRing3TobN)[idx];
                    int   bc_r4N   = (*gepIn_WTAConeBasicClustersJetsRing4TobN)[idx];
                    if(i == 0){
                        gepLeadingWTAConeGEPBasicClustersJetspTValues.push_back(coneWTAGEPBasicClustersJetspT);
                        gepLeadingWTAConeGEPBasicClustersJetsEtaValues.push_back(coneWTAGEPBasicClustersJetsEta);
                        gepLeadingWTAConeGEPBasicClustersJetsPhiValues.push_back(coneWTAGEPBasicClustersJetsPhi);
                        gepLeadingWTAConeGEPBasicClustersJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersJetsNConstituents);
                        gepLeadingWTAConeBasicClustersJetsRing0Et.push_back(bc_r0Et); gepLeadingWTAConeBasicClustersJetsRing1Et.push_back(bc_r1Et); gepLeadingWTAConeBasicClustersJetsRing2Et.push_back(bc_r2Et); gepLeadingWTAConeBasicClustersJetsRing3Et.push_back(bc_r3Et); gepLeadingWTAConeBasicClustersJetsRing4Et.push_back(bc_r4Et);
                        gepLeadingWTAConeBasicClustersJetsTotalTobN.push_back(bc_totN); gepLeadingWTAConeBasicClustersJetsRing0TobN.push_back(bc_r0N); gepLeadingWTAConeBasicClustersJetsRing1TobN.push_back(bc_r1N); gepLeadingWTAConeBasicClustersJetsRing2TobN.push_back(bc_r2N); gepLeadingWTAConeBasicClustersJetsRing3TobN.push_back(bc_r3N); gepLeadingWTAConeBasicClustersJetsRing4TobN.push_back(bc_r4N);
                    }
                    if (i == 1){
                        gepSubleadingWTAConeGEPBasicClustersJetspTValues.push_back(coneWTAGEPBasicClustersJetspT);
                        gepSubleadingWTAConeGEPBasicClustersJetsEtaValues.push_back(coneWTAGEPBasicClustersJetsEta);
                        gepSubleadingWTAConeGEPBasicClustersJetsPhiValues.push_back(coneWTAGEPBasicClustersJetsPhi);
                        gepSubleadingWTAConeGEPBasicClustersJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersJetsNConstituents);
                        gepSubleadingWTAConeBasicClustersJetsRing0Et.push_back(bc_r0Et); gepSubleadingWTAConeBasicClustersJetsRing1Et.push_back(bc_r1Et); gepSubleadingWTAConeBasicClustersJetsRing2Et.push_back(bc_r2Et); gepSubleadingWTAConeBasicClustersJetsRing3Et.push_back(bc_r3Et); gepSubleadingWTAConeBasicClustersJetsRing4Et.push_back(bc_r4Et);
                        gepSubleadingWTAConeBasicClustersJetsTotalTobN.push_back(bc_totN); gepSubleadingWTAConeBasicClustersJetsRing0TobN.push_back(bc_r0N); gepSubleadingWTAConeBasicClustersJetsRing1TobN.push_back(bc_r1N); gepSubleadingWTAConeBasicClustersJetsRing2TobN.push_back(bc_r2N); gepSubleadingWTAConeBasicClustersJetsRing3TobN.push_back(bc_r3N); gepSubleadingWTAConeBasicClustersJetsRing4TobN.push_back(bc_r4N);
                    }

                    gepWTAConeGEPBasicClustersJetspTValues.push_back(coneWTAGEPBasicClustersJetspT);
                    gepWTAConeGEPBasicClustersJetsEtaValues.push_back(coneWTAGEPBasicClustersJetsEta);
                    gepWTAConeGEPBasicClustersJetsPhiValues.push_back(coneWTAGEPBasicClustersJetsPhi);
                    gepWTAConeGEPBasicClustersJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersJetsNConstituents);
                    gepWTAConeBasicClustersJetsRing0Et.push_back(bc_r0Et); gepWTAConeBasicClustersJetsRing1Et.push_back(bc_r1Et); gepWTAConeBasicClustersJetsRing2Et.push_back(bc_r2Et); gepWTAConeBasicClustersJetsRing3Et.push_back(bc_r3Et); gepWTAConeBasicClustersJetsRing4Et.push_back(bc_r4Et);
                    gepWTAConeBasicClustersJetsTotalTobN.push_back(bc_totN); gepWTAConeBasicClustersJetsRing0TobN.push_back(bc_r0N); gepWTAConeBasicClustersJetsRing1TobN.push_back(bc_r1N); gepWTAConeBasicClustersJetsRing2TobN.push_back(bc_r2N); gepWTAConeBasicClustersJetsRing3TobN.push_back(bc_r3N); gepWTAConeBasicClustersJetsRing4TobN.push_back(bc_r4N);

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


                    if (iEvt < maxDatEvents) {
                        if (gepwtaconebasicclusters_it == 0) {
                            f_wtaconejetsbasicclusters << "Event : " << std::dec << iEvt << "\n";
                        }
                        f_wtaconejetsbasicclusters << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepwtaconebasicclusters_it
                                            << " "  << binary_word
                                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    }
                    ++gepwtaconebasicclusters_it;

                }
            }

            // GEP cone jets from basic clusters with SoftKiller PU Suppression applied
            {
                std::vector<unsigned int> indices(gepIn_WTAConeBasicClustersSKJetsPt->size());
                std::iota(indices.begin(), indices.end(), 0);

                std::sort(indices.begin(), indices.end(),
                        [&](unsigned int a, unsigned int b) {
                            return (*gepIn_WTAConeBasicClustersSKJetsPt)[a] > (*gepIn_WTAConeBasicClustersSKJetsPt)[b];
                        });
                unsigned int gepwtaconebasicclusterssk_it = 0;
                for (unsigned int i = 0; i < indices.size(); i++) {
                    unsigned int idx = indices[i];
                    float coneWTAGEPBasicClustersSKJetspT  = (*gepIn_WTAConeBasicClustersSKJetsPt)[idx] / 1000.0;
                    float coneWTAGEPBasicClustersSKJetsEta = (*gepIn_WTAConeBasicClustersSKJetsEta)[idx];
                    float coneWTAGEPBasicClustersSKJetsPhi = (*gepIn_WTAConeBasicClustersSKJetsPhi)[idx];
                    unsigned int coneWTAGEPBasicClustersSKJetsNConstituents = (*gepIn_WTAConeBasicClustersJetsNConstituents)[idx];
                    float bcsk_r0Et = (*gepIn_WTAConeBasicClustersSKJetsRing0Et)[idx];
                    float bcsk_r1Et = (*gepIn_WTAConeBasicClustersSKJetsRing1Et)[idx];
                    float bcsk_r2Et = (*gepIn_WTAConeBasicClustersSKJetsRing2Et)[idx];
                    float bcsk_r3Et = (*gepIn_WTAConeBasicClustersSKJetsRing3Et)[idx];
                    float bcsk_r4Et = (*gepIn_WTAConeBasicClustersSKJetsRing4Et)[idx];
                    int   bcsk_totN = (*gepIn_WTAConeBasicClustersSKJetsTotalTobN)[idx];
                    int   bcsk_r0N  = (*gepIn_WTAConeBasicClustersSKJetsRing0TobN)[idx];
                    int   bcsk_r1N  = (*gepIn_WTAConeBasicClustersSKJetsRing1TobN)[idx];
                    int   bcsk_r2N  = (*gepIn_WTAConeBasicClustersSKJetsRing2TobN)[idx];
                    int   bcsk_r3N  = (*gepIn_WTAConeBasicClustersSKJetsRing3TobN)[idx];
                    int   bcsk_r4N  = (*gepIn_WTAConeBasicClustersSKJetsRing4TobN)[idx];
                    if(i == 0){
                        gepLeadingWTAConeGEPBasicClustersSKJetspTValues.push_back(coneWTAGEPBasicClustersSKJetspT);
                        gepLeadingWTAConeGEPBasicClustersSKJetsEtaValues.push_back(coneWTAGEPBasicClustersSKJetsEta);
                        gepLeadingWTAConeGEPBasicClustersSKJetsPhiValues.push_back(coneWTAGEPBasicClustersSKJetsPhi);
                        gepLeadingWTAConeGEPBasicClustersSKJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersSKJetsNConstituents);
                        gepLeadingWTAConeBasicClustersSKJetsRing0Et.push_back(bcsk_r0Et); gepLeadingWTAConeBasicClustersSKJetsRing1Et.push_back(bcsk_r1Et); gepLeadingWTAConeBasicClustersSKJetsRing2Et.push_back(bcsk_r2Et); gepLeadingWTAConeBasicClustersSKJetsRing3Et.push_back(bcsk_r3Et); gepLeadingWTAConeBasicClustersSKJetsRing4Et.push_back(bcsk_r4Et);
                        gepLeadingWTAConeBasicClustersSKJetsTotalTobN.push_back(bcsk_totN); gepLeadingWTAConeBasicClustersSKJetsRing0TobN.push_back(bcsk_r0N); gepLeadingWTAConeBasicClustersSKJetsRing1TobN.push_back(bcsk_r1N); gepLeadingWTAConeBasicClustersSKJetsRing2TobN.push_back(bcsk_r2N); gepLeadingWTAConeBasicClustersSKJetsRing3TobN.push_back(bcsk_r3N); gepLeadingWTAConeBasicClustersSKJetsRing4TobN.push_back(bcsk_r4N);
                    }
                    if (i == 1){
                        gepSubleadingWTAConeGEPBasicClustersSKJetspTValues.push_back(coneWTAGEPBasicClustersSKJetspT);
                        gepSubleadingWTAConeGEPBasicClustersSKJetsEtaValues.push_back(coneWTAGEPBasicClustersSKJetsEta);
                        gepSubleadingWTAConeGEPBasicClustersSKJetsPhiValues.push_back(coneWTAGEPBasicClustersSKJetsPhi);
                        gepSubleadingWTAConeGEPBasicClustersSKJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersSKJetsNConstituents);
                        gepSubleadingWTAConeBasicClustersSKJetsRing0Et.push_back(bcsk_r0Et); gepSubleadingWTAConeBasicClustersSKJetsRing1Et.push_back(bcsk_r1Et); gepSubleadingWTAConeBasicClustersSKJetsRing2Et.push_back(bcsk_r2Et); gepSubleadingWTAConeBasicClustersSKJetsRing3Et.push_back(bcsk_r3Et); gepSubleadingWTAConeBasicClustersSKJetsRing4Et.push_back(bcsk_r4Et);
                        gepSubleadingWTAConeBasicClustersSKJetsTotalTobN.push_back(bcsk_totN); gepSubleadingWTAConeBasicClustersSKJetsRing0TobN.push_back(bcsk_r0N); gepSubleadingWTAConeBasicClustersSKJetsRing1TobN.push_back(bcsk_r1N); gepSubleadingWTAConeBasicClustersSKJetsRing2TobN.push_back(bcsk_r2N); gepSubleadingWTAConeBasicClustersSKJetsRing3TobN.push_back(bcsk_r3N); gepSubleadingWTAConeBasicClustersSKJetsRing4TobN.push_back(bcsk_r4N);
                    }

                    gepWTAConeGEPBasicClustersSKJetspTValues.push_back(coneWTAGEPBasicClustersSKJetspT);
                    gepWTAConeGEPBasicClustersSKJetsEtaValues.push_back(coneWTAGEPBasicClustersSKJetsEta);
                    gepWTAConeGEPBasicClustersSKJetsPhiValues.push_back(coneWTAGEPBasicClustersSKJetsPhi);
                    gepWTAConeGEPBasicClustersSKJetsNConstituentsValues.push_back(coneWTAGEPBasicClustersSKJetsNConstituents);
                    gepWTAConeBasicClustersSKJetsRing0Et.push_back(bcsk_r0Et); gepWTAConeBasicClustersSKJetsRing1Et.push_back(bcsk_r1Et); gepWTAConeBasicClustersSKJetsRing2Et.push_back(bcsk_r2Et); gepWTAConeBasicClustersSKJetsRing3Et.push_back(bcsk_r3Et); gepWTAConeBasicClustersSKJetsRing4Et.push_back(bcsk_r4Et);
                    gepWTAConeBasicClustersSKJetsTotalTobN.push_back(bcsk_totN); gepWTAConeBasicClustersSKJetsRing0TobN.push_back(bcsk_r0N); gepWTAConeBasicClustersSKJetsRing1TobN.push_back(bcsk_r1N); gepWTAConeBasicClustersSKJetsRing2TobN.push_back(bcsk_r2N); gepWTAConeBasicClustersSKJetsRing3TobN.push_back(bcsk_r3N); gepWTAConeBasicClustersSKJetsRing4TobN.push_back(bcsk_r4N);

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

                    if (iEvt < maxDatEvents) {
                        if (gepwtaconebasicclusterssk_it == 0) {
                            f_wtaconejetsbasicclusterssk << "Event : " << std::dec << iEvt << "\n";
                        }
                        f_wtaconejetsbasicclusterssk << "0x" << std::hex << std::setw(2) << std::setfill('0') << gepwtaconebasicclusterssk_it
                                            << " "  << binary_word
                                            << " 0x" << std::setw(8) << std::setfill('0') << packed_word << "\n";
                    }
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
                if (iEvt < maxDatEvents) {
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
                if (iEvt < maxDatEvents) {
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

            // GEP EtaSK tower collection
            for (size_t i = 0; i < gepIn_CellsTowersEtaSKEt->size(); ++i) {
                gepCellsTowersEtaSKEtValues.push_back((*gepIn_CellsTowersEtaSKEt)[i] / 1000.0);
                gepCellsTowersEtaSKEtaValues.push_back((*gepIn_CellsTowersEtaSKEta)[i]);
                gepCellsTowersEtaSKPhiValues.push_back((*gepIn_CellsTowersEtaSKPhi)[i]);
            }

            // GEP EtaSK cluster collection
            for (size_t i = 0; i < gepIn_BasicClustersEtaSKEt->size(); ++i) {
                gepBasicClustersEtaSKEtValues.push_back((*gepIn_BasicClustersEtaSKEt)[i] / 1000.0);
                gepBasicClustersEtaSKEtaValues.push_back((*gepIn_BasicClustersEtaSKEta)[i]);
                gepBasicClustersEtaSKPhiValues.push_back((*gepIn_BasicClustersEtaSKPhi)[i]);
            }

            // GEP EtaSK WTACone CellsTower jets
            {
                std::vector<unsigned int> indices(gepIn_WTAConeCellsTowersEtaSKJetsPt->size());
                std::iota(indices.begin(), indices.end(), 0);
                std::sort(indices.begin(), indices.end(),
                    [&](unsigned int a, unsigned int b) {
                        return (*gepIn_WTAConeCellsTowersEtaSKJetsPt)[a] > (*gepIn_WTAConeCellsTowersEtaSKJetsPt)[b];
                    });
                for (unsigned int i = 0; i < indices.size(); i++) {
                    unsigned int idx = indices[i];
                    float etaskct_pT  = (*gepIn_WTAConeCellsTowersEtaSKJetsPt)[idx] / 1000.0;
                    float etaskct_eta = (*gepIn_WTAConeCellsTowersEtaSKJetsEta)[idx];
                    float etaskct_phi = (*gepIn_WTAConeCellsTowersEtaSKJetsPhi)[idx];
                    unsigned int etaskct_nc = (*gepIn_WTAConeCellsTowersEtaSKJetsNConstituents)[idx];
                    float etaskct_r0Et = (*gepIn_WTAConeCellsTowersEtaSKJetsRing0Et)[idx];
                    float etaskct_r1Et = (*gepIn_WTAConeCellsTowersEtaSKJetsRing1Et)[idx];
                    float etaskct_r2Et = (*gepIn_WTAConeCellsTowersEtaSKJetsRing2Et)[idx];
                    float etaskct_r3Et = (*gepIn_WTAConeCellsTowersEtaSKJetsRing3Et)[idx];
                    float etaskct_r4Et = (*gepIn_WTAConeCellsTowersEtaSKJetsRing4Et)[idx];
                    int   etaskct_totN = (*gepIn_WTAConeCellsTowersEtaSKJetsTotalTobN)[idx];
                    int   etaskct_r0N  = (*gepIn_WTAConeCellsTowersEtaSKJetsRing0TobN)[idx];
                    int   etaskct_r1N  = (*gepIn_WTAConeCellsTowersEtaSKJetsRing1TobN)[idx];
                    int   etaskct_r2N  = (*gepIn_WTAConeCellsTowersEtaSKJetsRing2TobN)[idx];
                    int   etaskct_r3N  = (*gepIn_WTAConeCellsTowersEtaSKJetsRing3TobN)[idx];
                    int   etaskct_r4N  = (*gepIn_WTAConeCellsTowersEtaSKJetsRing4TobN)[idx];
                    if (i == 0) {
                        gepLeadingWTAConeCellsTowersEtaSKJetspTValues.push_back(etaskct_pT);
                        gepLeadingWTAConeCellsTowersEtaSKJetsEtaValues.push_back(etaskct_eta);
                        gepLeadingWTAConeCellsTowersEtaSKJetsPhiValues.push_back(etaskct_phi);
                        gepLeadingWTAConeCellsTowersEtaSKJetsNConstituentsValues.push_back(etaskct_nc);
                        gepLeadingWTAConeCellsTowersEtaSKJetsRing0Et.push_back(etaskct_r0Et); gepLeadingWTAConeCellsTowersEtaSKJetsRing1Et.push_back(etaskct_r1Et); gepLeadingWTAConeCellsTowersEtaSKJetsRing2Et.push_back(etaskct_r2Et); gepLeadingWTAConeCellsTowersEtaSKJetsRing3Et.push_back(etaskct_r3Et); gepLeadingWTAConeCellsTowersEtaSKJetsRing4Et.push_back(etaskct_r4Et);
                        gepLeadingWTAConeCellsTowersEtaSKJetsTotalTobN.push_back(etaskct_totN); gepLeadingWTAConeCellsTowersEtaSKJetsRing0TobN.push_back(etaskct_r0N); gepLeadingWTAConeCellsTowersEtaSKJetsRing1TobN.push_back(etaskct_r1N); gepLeadingWTAConeCellsTowersEtaSKJetsRing2TobN.push_back(etaskct_r2N); gepLeadingWTAConeCellsTowersEtaSKJetsRing3TobN.push_back(etaskct_r3N); gepLeadingWTAConeCellsTowersEtaSKJetsRing4TobN.push_back(etaskct_r4N);
                    }
                    if (i == 1) {
                        gepSubleadingWTAConeCellsTowersEtaSKJetspTValues.push_back(etaskct_pT);
                        gepSubleadingWTAConeCellsTowersEtaSKJetsEtaValues.push_back(etaskct_eta);
                        gepSubleadingWTAConeCellsTowersEtaSKJetsPhiValues.push_back(etaskct_phi);
                        gepSubleadingWTAConeCellsTowersEtaSKJetsNConstituentsValues.push_back(etaskct_nc);
                        gepSubleadingWTAConeCellsTowersEtaSKJetsRing0Et.push_back(etaskct_r0Et); gepSubleadingWTAConeCellsTowersEtaSKJetsRing1Et.push_back(etaskct_r1Et); gepSubleadingWTAConeCellsTowersEtaSKJetsRing2Et.push_back(etaskct_r2Et); gepSubleadingWTAConeCellsTowersEtaSKJetsRing3Et.push_back(etaskct_r3Et); gepSubleadingWTAConeCellsTowersEtaSKJetsRing4Et.push_back(etaskct_r4Et);
                        gepSubleadingWTAConeCellsTowersEtaSKJetsTotalTobN.push_back(etaskct_totN); gepSubleadingWTAConeCellsTowersEtaSKJetsRing0TobN.push_back(etaskct_r0N); gepSubleadingWTAConeCellsTowersEtaSKJetsRing1TobN.push_back(etaskct_r1N); gepSubleadingWTAConeCellsTowersEtaSKJetsRing2TobN.push_back(etaskct_r2N); gepSubleadingWTAConeCellsTowersEtaSKJetsRing3TobN.push_back(etaskct_r3N); gepSubleadingWTAConeCellsTowersEtaSKJetsRing4TobN.push_back(etaskct_r4N);
                    }
                    gepWTAConeCellsTowersEtaSKJetspTValues.push_back(etaskct_pT);
                    gepWTAConeCellsTowersEtaSKJetsEtaValues.push_back(etaskct_eta);
                    gepWTAConeCellsTowersEtaSKJetsPhiValues.push_back(etaskct_phi);
                    gepWTAConeCellsTowersEtaSKJetsNConstituentsValues.push_back(etaskct_nc);
                    gepWTAConeCellsTowersEtaSKJetsRing0Et.push_back(etaskct_r0Et); gepWTAConeCellsTowersEtaSKJetsRing1Et.push_back(etaskct_r1Et); gepWTAConeCellsTowersEtaSKJetsRing2Et.push_back(etaskct_r2Et); gepWTAConeCellsTowersEtaSKJetsRing3Et.push_back(etaskct_r3Et); gepWTAConeCellsTowersEtaSKJetsRing4Et.push_back(etaskct_r4Et);
                    gepWTAConeCellsTowersEtaSKJetsTotalTobN.push_back(etaskct_totN); gepWTAConeCellsTowersEtaSKJetsRing0TobN.push_back(etaskct_r0N); gepWTAConeCellsTowersEtaSKJetsRing1TobN.push_back(etaskct_r1N); gepWTAConeCellsTowersEtaSKJetsRing2TobN.push_back(etaskct_r2N); gepWTAConeCellsTowersEtaSKJetsRing3TobN.push_back(etaskct_r3N); gepWTAConeCellsTowersEtaSKJetsRing4TobN.push_back(etaskct_r4N);
                }
            }

            // GEP EtaSK WTACone BasicClusters jets
            {
                std::vector<unsigned int> indices(gepIn_WTAConeBasicClustersEtaSKJetsPt->size());
                std::iota(indices.begin(), indices.end(), 0);
                std::sort(indices.begin(), indices.end(),
                    [&](unsigned int a, unsigned int b) {
                        return (*gepIn_WTAConeBasicClustersEtaSKJetsPt)[a] > (*gepIn_WTAConeBasicClustersEtaSKJetsPt)[b];
                    });
                for (unsigned int i = 0; i < indices.size(); i++) {
                    unsigned int idx = indices[i];
                    float etaskbc_pT  = (*gepIn_WTAConeBasicClustersEtaSKJetsPt)[idx] / 1000.0;
                    float etaskbc_eta = (*gepIn_WTAConeBasicClustersEtaSKJetsEta)[idx];
                    float etaskbc_phi = (*gepIn_WTAConeBasicClustersEtaSKJetsPhi)[idx];
                    unsigned int etaskbc_nc = (*gepIn_WTAConeBasicClustersEtaSKJetsNConstituents)[idx];
                    float etaskbc_r0Et = (*gepIn_WTAConeBasicClustersEtaSKJetsRing0Et)[idx];
                    float etaskbc_r1Et = (*gepIn_WTAConeBasicClustersEtaSKJetsRing1Et)[idx];
                    float etaskbc_r2Et = (*gepIn_WTAConeBasicClustersEtaSKJetsRing2Et)[idx];
                    float etaskbc_r3Et = (*gepIn_WTAConeBasicClustersEtaSKJetsRing3Et)[idx];
                    float etaskbc_r4Et = (*gepIn_WTAConeBasicClustersEtaSKJetsRing4Et)[idx];
                    int   etaskbc_totN = (*gepIn_WTAConeBasicClustersEtaSKJetsTotalTobN)[idx];
                    int   etaskbc_r0N  = (*gepIn_WTAConeBasicClustersEtaSKJetsRing0TobN)[idx];
                    int   etaskbc_r1N  = (*gepIn_WTAConeBasicClustersEtaSKJetsRing1TobN)[idx];
                    int   etaskbc_r2N  = (*gepIn_WTAConeBasicClustersEtaSKJetsRing2TobN)[idx];
                    int   etaskbc_r3N  = (*gepIn_WTAConeBasicClustersEtaSKJetsRing3TobN)[idx];
                    int   etaskbc_r4N  = (*gepIn_WTAConeBasicClustersEtaSKJetsRing4TobN)[idx];
                    if (i == 0) {
                        gepLeadingWTAConeBasicClustersEtaSKJetspTValues.push_back(etaskbc_pT);
                        gepLeadingWTAConeBasicClustersEtaSKJetsEtaValues.push_back(etaskbc_eta);
                        gepLeadingWTAConeBasicClustersEtaSKJetsPhiValues.push_back(etaskbc_phi);
                        gepLeadingWTAConeBasicClustersEtaSKJetsNConstituentsValues.push_back(etaskbc_nc);
                        gepLeadingWTAConeBasicClustersEtaSKJetsRing0Et.push_back(etaskbc_r0Et); gepLeadingWTAConeBasicClustersEtaSKJetsRing1Et.push_back(etaskbc_r1Et); gepLeadingWTAConeBasicClustersEtaSKJetsRing2Et.push_back(etaskbc_r2Et); gepLeadingWTAConeBasicClustersEtaSKJetsRing3Et.push_back(etaskbc_r3Et); gepLeadingWTAConeBasicClustersEtaSKJetsRing4Et.push_back(etaskbc_r4Et);
                        gepLeadingWTAConeBasicClustersEtaSKJetsTotalTobN.push_back(etaskbc_totN); gepLeadingWTAConeBasicClustersEtaSKJetsRing0TobN.push_back(etaskbc_r0N); gepLeadingWTAConeBasicClustersEtaSKJetsRing1TobN.push_back(etaskbc_r1N); gepLeadingWTAConeBasicClustersEtaSKJetsRing2TobN.push_back(etaskbc_r2N); gepLeadingWTAConeBasicClustersEtaSKJetsRing3TobN.push_back(etaskbc_r3N); gepLeadingWTAConeBasicClustersEtaSKJetsRing4TobN.push_back(etaskbc_r4N);
                    }
                    if (i == 1) {
                        gepSubleadingWTAConeBasicClustersEtaSKJetspTValues.push_back(etaskbc_pT);
                        gepSubleadingWTAConeBasicClustersEtaSKJetsEtaValues.push_back(etaskbc_eta);
                        gepSubleadingWTAConeBasicClustersEtaSKJetsPhiValues.push_back(etaskbc_phi);
                        gepSubleadingWTAConeBasicClustersEtaSKJetsNConstituentsValues.push_back(etaskbc_nc);
                        gepSubleadingWTAConeBasicClustersEtaSKJetsRing0Et.push_back(etaskbc_r0Et); gepSubleadingWTAConeBasicClustersEtaSKJetsRing1Et.push_back(etaskbc_r1Et); gepSubleadingWTAConeBasicClustersEtaSKJetsRing2Et.push_back(etaskbc_r2Et); gepSubleadingWTAConeBasicClustersEtaSKJetsRing3Et.push_back(etaskbc_r3Et); gepSubleadingWTAConeBasicClustersEtaSKJetsRing4Et.push_back(etaskbc_r4Et);
                        gepSubleadingWTAConeBasicClustersEtaSKJetsTotalTobN.push_back(etaskbc_totN); gepSubleadingWTAConeBasicClustersEtaSKJetsRing0TobN.push_back(etaskbc_r0N); gepSubleadingWTAConeBasicClustersEtaSKJetsRing1TobN.push_back(etaskbc_r1N); gepSubleadingWTAConeBasicClustersEtaSKJetsRing2TobN.push_back(etaskbc_r2N); gepSubleadingWTAConeBasicClustersEtaSKJetsRing3TobN.push_back(etaskbc_r3N); gepSubleadingWTAConeBasicClustersEtaSKJetsRing4TobN.push_back(etaskbc_r4N);
                    }
                    gepWTAConeBasicClustersEtaSKJetspTValues.push_back(etaskbc_pT);
                    gepWTAConeBasicClustersEtaSKJetsEtaValues.push_back(etaskbc_eta);
                    gepWTAConeBasicClustersEtaSKJetsPhiValues.push_back(etaskbc_phi);
                    gepWTAConeBasicClustersEtaSKJetsNConstituentsValues.push_back(etaskbc_nc);
                    gepWTAConeBasicClustersEtaSKJetsRing0Et.push_back(etaskbc_r0Et); gepWTAConeBasicClustersEtaSKJetsRing1Et.push_back(etaskbc_r1Et); gepWTAConeBasicClustersEtaSKJetsRing2Et.push_back(etaskbc_r2Et); gepWTAConeBasicClustersEtaSKJetsRing3Et.push_back(etaskbc_r3Et); gepWTAConeBasicClustersEtaSKJetsRing4Et.push_back(etaskbc_r4Et);
                    gepWTAConeBasicClustersEtaSKJetsTotalTobN.push_back(etaskbc_totN); gepWTAConeBasicClustersEtaSKJetsRing0TobN.push_back(etaskbc_r0N); gepWTAConeBasicClustersEtaSKJetsRing1TobN.push_back(etaskbc_r1N); gepWTAConeBasicClustersEtaSKJetsRing2TobN.push_back(etaskbc_r2N); gepWTAConeBasicClustersEtaSKJetsRing3TobN.push_back(etaskbc_r3N); gepWTAConeBasicClustersEtaSKJetsRing4TobN.push_back(etaskbc_r4N);
                }
            }

            // Sim jFEX SR jets (from TrigGepPerf ntuple)
            {
                std::vector<std::pair<size_t, double>> jFexSRJSimEtWithIndex;
                for (size_t i = 0; i < jFexSRJSimEt->size(); ++i)
                    jFexSRJSimEtWithIndex.emplace_back(i, (*jFexSRJSimEt)[i] / 1000.0);
                std::sort(jFexSRJSimEtWithIndex.begin(), jFexSRJSimEtWithIndex.end(),
                    [](const std::pair<size_t,double>& a, const std::pair<size_t,double>& b){ return a.second > b.second; });
                for (const auto& [index, et] : jFexSRJSimEtWithIndex) {
                    jFexSRJSimEtIndexValues.push_back(static_cast<unsigned int>(index));
                    jFexSRJSimEtValues.push_back(et);
                    jFexSRJSimEtaValues.push_back((*jFexSRJSimEta)[index]);
                    jFexSRJSimPhiValues.push_back((*jFexSRJSimPhi)[index]);
                }
                if (!jFexSRJSimEtWithIndex.empty()) {
                    jFexSRJSimLeadingEtValues.push_back(jFexSRJSimEtWithIndex[0].second);
                    jFexSRJSimLeadingEtaValues.push_back((*jFexSRJSimEta)[jFexSRJSimEtWithIndex[0].first]);
                    jFexSRJSimLeadingPhiValues.push_back((*jFexSRJSimPhi)[jFexSRJSimEtWithIndex[0].first]);
                }
                if (jFexSRJSimEtWithIndex.size() > 1) {
                    jFexSRJSimSubleadingEtValues.push_back(jFexSRJSimEtWithIndex[1].second);
                    jFexSRJSimSubleadingEtaValues.push_back((*jFexSRJSimEta)[jFexSRJSimEtWithIndex[1].first]);
                    jFexSRJSimSubleadingPhiValues.push_back((*jFexSRJSimPhi)[jFexSRJSimEtWithIndex[1].first]);
                }
            }

            // Sim gFEX LR jets (from TrigGepPerf ntuple)
            {
                std::vector<std::pair<size_t, double>> gFexLRJSimEtWithIndex;
                for (size_t i = 0; i < gFexLRJSimEt->size(); ++i)
                    gFexLRJSimEtWithIndex.emplace_back(i, (*gFexLRJSimEt)[i] / 1000.0);
                std::sort(gFexLRJSimEtWithIndex.begin(), gFexLRJSimEtWithIndex.end(),
                    [](const std::pair<size_t,double>& a, const std::pair<size_t,double>& b){ return a.second > b.second; });
                for (const auto& [index, et] : gFexLRJSimEtWithIndex) {
                    gFexLRJSimEtIndexValues.push_back(static_cast<unsigned int>(index));
                    gFexLRJSimEtValues.push_back(et);
                    gFexLRJSimEtaValues.push_back((*gFexLRJSimEta)[index]);
                    gFexLRJSimPhiValues.push_back((*gFexLRJSimPhi)[index]);
                }
                if (!gFexLRJSimEtWithIndex.empty()) {
                    gFexLRJSimLeadingEtValues.push_back(gFexLRJSimEtWithIndex[0].second);
                    gFexLRJSimLeadingEtaValues.push_back((*gFexLRJSimEta)[gFexLRJSimEtWithIndex[0].first]);
                    gFexLRJSimLeadingPhiValues.push_back((*gFexLRJSimPhi)[gFexLRJSimEtWithIndex[0].first]);
                }
                if (gFexLRJSimEtWithIndex.size() > 1) {
                    gFexLRJSimSubleadingEtValues.push_back(gFexLRJSimEtWithIndex[1].second);
                    gFexLRJSimSubleadingEtaValues.push_back((*gFexLRJSimEta)[gFexLRJSimEtWithIndex[1].first]);
                    gFexLRJSimSubleadingPhiValues.push_back((*gFexLRJSimPhi)[gFexLRJSimEtWithIndex[1].first]);
                }
            }

            // TODO: uncomment once L1_gFexSRJetRoISim branch exists in GEP ntuples
            // // Sim gFEX SR jets (from TrigGepPerf ntuple)
            // {
            //     std::vector<std::pair<size_t, double>> gFexSRJSimEtWithIndex;
            //     for (size_t i = 0; i < gFexSRJSimEt->size(); ++i)
            //         gFexSRJSimEtWithIndex.emplace_back(i, (*gFexSRJSimEt)[i] / 1000.0);
            //     std::sort(gFexSRJSimEtWithIndex.begin(), gFexSRJSimEtWithIndex.end(),
            //         [](const std::pair<size_t,double>& a, const std::pair<size_t,double>& b){ return a.second > b.second; });
            //     for (const auto& [index, et] : gFexSRJSimEtWithIndex) {
            //         gFexSRJSimEtIndexValues.push_back(static_cast<unsigned int>(index));
            //         gFexSRJSimEtValues.push_back(et);
            //         gFexSRJSimEtaValues.push_back((*gFexSRJSimEta)[index]);
            //         gFexSRJSimPhiValues.push_back((*gFexSRJSimPhi)[index]);
            //     }
            //     if (!gFexSRJSimEtWithIndex.empty()) {
            //         gFexSRJSimLeadingEtValues.push_back(gFexSRJSimEtWithIndex[0].second);
            //         gFexSRJSimLeadingEtaValues.push_back((*gFexSRJSimEta)[gFexSRJSimEtWithIndex[0].first]);
            //         gFexSRJSimLeadingPhiValues.push_back((*gFexSRJSimPhi)[gFexSRJSimEtWithIndex[0].first]);
            //     }
            //     if (gFexSRJSimEtWithIndex.size() > 1) {
            //         gFexSRJSimSubleadingEtValues.push_back(gFexSRJSimEtWithIndex[1].second);
            //         gFexSRJSimSubleadingEtaValues.push_back((*gFexSRJSimEta)[gFexSRJSimEtWithIndex[1].first]);
            //         gFexSRJSimSubleadingPhiValues.push_back((*gFexSRJSimPhi)[gFexSRJSimEtWithIndex[1].first]);
            //     }
            // }

            // Sim gFEX/jFEX MET (from TrigGepPerf ntuple — scalar per event, already in GeV)
            gFexMETJwoJSim       = gepIn_gFexJwoJMet    / 1000.0;
            gFexMETPhiJwoJSim    = gepIn_gFexJwoJMetPhi;
            gFexSumETJwoJSim     = gepIn_gFexJwoJSumEt  / 1000.0;
            gFexMETNoiseCutSim   = gepIn_gFexNCMet      / 1000.0;
            gFexMETPhiNoiseCutSim = gepIn_gFexNCMetPhi;
            gFexSumETNoiseCutSim  = gepIn_gFexNCSumEt   / 1000.0;
            gFexMETRmsSim        = gepIn_gFexRhoRMSMet      / 1000.0;
            gFexMETPhiRmsSim     = gepIn_gFexRhoRMSMetPhi;
            gFexSumETRmsSim      = gepIn_gFexRhoRMSSumEt    / 1000.0;
            jFexMET              = gepIn_jFexMet        / 1000.0;

            // Event ID from GEP ntuple
            // gepRunNumber / gepEventNumber already set by gt->GetEntry(iEvt) above

            // --- gFEX MET ---
            gMHX = 0.0; gMHY = 0.0;
            gMSX = 0.0; gMSY = 0.0;
            gMETX = 0.0; gMETY = 0.0;
            gMET = 0.0; gSumET = 0.0;
            gMETX_NC = 0.0; gMETY_NC = 0.0; gMET_NC = 0.0; gSumET_NC = 0.0;
            gMETX_Rms = 0.0; gMETY_Rms = 0.0; gMET_Rms = 0.0; gSumET_Rms = 0.0;
            metTruthNonIntX    = 0.0; metTruthNonIntY    = 0.0; metTruthNonIntSumET    = 0.0; metTruthNonInt    = 0.0;
            metTruthIntX       = 0.0; metTruthIntY       = 0.0; metTruthIntSumET       = 0.0; metTruthInt       = 0.0;
            metTruthIntOutX    = 0.0; metTruthIntOutY    = 0.0; metTruthIntOutSumET    = 0.0; metTruthIntOut    = 0.0;
            metTruthIntMuonsX  = 0.0; metTruthIntMuonsY  = 0.0; metTruthIntMuonsSumET  = 0.0; metTruthIntMuons  = 0.0;
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
                std::cout << "MET_Truth terms available:\n";
                for (const xAOD::MissingET* m : *MET_Truth){
                    if(iEvt <= 5){
                        std::cout << "  name=" << m->name() << "  met=" << std::sqrt(m->mpx()*m->mpx()+m->mpy()*m->mpy())/1000.0 << " GeV\n";
                    }
                }
                    

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
                fillTruthTerm("NonInt",   metTruthNonIntX,   metTruthNonIntY,   metTruthNonIntSumET,   metTruthNonInt);
                fillTruthTerm("Int",      metTruthIntX,      metTruthIntY,      metTruthIntSumET,      metTruthInt);
                fillTruthTerm("IntOut",   metTruthIntOutX,   metTruthIntOutY,   metTruthIntOutSumET,   metTruthIntOut);
                fillTruthTerm("IntMuons", metTruthIntMuonsX, metTruthIntMuonsY, metTruthIntMuonsSumET, metTruthIntMuons);
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

            // --- Loop over OutOfTimeAntiKt4TruthJets ---
            std::vector<std::pair<size_t, double>> outOfTimeTruthJetEtWithIndex;
            for (size_t i = 0; i < OutOfTimeAntiKt4TruthJets->size(); ++i) {
                const auto& jet = (*OutOfTimeAntiKt4TruthJets)[i];
                double et = jet->e() / (1000.0 * cosh(jet->eta()));
                outOfTimeTruthJetEtWithIndex.emplace_back(i, et);
            }
            std::sort(outOfTimeTruthJetEtWithIndex.begin(), outOfTimeTruthJetEtWithIndex.end(),
                    [](const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
                        return a.second > b.second;
                    });
            for (const auto& [index, et] : outOfTimeTruthJetEtWithIndex) {
                const auto& jet = (*OutOfTimeAntiKt4TruthJets)[index];
                outOfTimeAntiKt4TruthSRJEtIndexValues.push_back(static_cast<unsigned int>(index));
                outOfTimeAntiKt4TruthSRJEtValues.push_back(et);
                outOfTimeAntiKt4TruthSRJEtaValues.push_back(jet->eta());
                outOfTimeAntiKt4TruthSRJPhiValues.push_back(jet->phi());
            }
            if (!outOfTimeTruthJetEtWithIndex.empty()) {
                const auto& leading = (*OutOfTimeAntiKt4TruthJets)[outOfTimeTruthJetEtWithIndex[0].first];
                outOfTimeAntiKt4TruthSRJLeadingEtValues.push_back(outOfTimeAntiKt4TruthSRJEtValues[0]);
                outOfTimeAntiKt4TruthSRJLeadingEtaValues.push_back(leading->eta());
                outOfTimeAntiKt4TruthSRJLeadingPhiValues.push_back(leading->phi());
            }
            if (outOfTimeTruthJetEtWithIndex.size() > 1) {
                const auto& subleading = (*OutOfTimeAntiKt4TruthJets)[outOfTimeTruthJetEtWithIndex[1].first];
                outOfTimeAntiKt4TruthSRJSubleadingEtValues.push_back(outOfTimeAntiKt4TruthSRJEtValues[1]);
                outOfTimeAntiKt4TruthSRJSubleadingEtaValues.push_back(subleading->eta());
                outOfTimeAntiKt4TruthSRJSubleadingPhiValues.push_back(subleading->phi());
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
            outOfTimeAntiKt4TruthJetsTree->Fill();
            leadingOutOfTimeAntiKt4TruthJetsTree->Fill();
            subleadingOutOfTimeAntiKt4TruthJetsTree->Fill();
            jFexSRJTree->Fill();
            jFexLeadingSRJTree->Fill();
            jFexSubleadingSRJTree->Fill();
            jFexLRJTree->Fill();
            jFexLeadingLRJTree->Fill();
            jFexSubleadingLRJTree->Fill();
            jFexSRJSimTree->Fill();
            jFexLeadingSRJSimTree->Fill();
            jFexSubleadingSRJSimTree->Fill();
            gFexLRJSimTree->Fill();
            gFexLeadingLRJSimTree->Fill();
            gFexSubleadingLRJSimTree->Fill();
            // TODO: uncomment once L1_gFexSRJetRoISim branch exists in GEP ntuples
            // gFexSRJSimTree->Fill();
            // gFexLeadingSRJSimTree->Fill();
            // gFexSubleadingSRJSimTree->Fill();
            gepCellsTowersEtaSKTree->Fill();
            gepBasicClustersEtaSKTree->Fill();
            gepWTAConeCellsTowersEtaSKJetsTree->Fill();
            gepLeadingWTAConeCellsTowersEtaSKJetsTree->Fill();
            gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Fill();
            gepWTAConeBasicClustersEtaSKJetsTree->Fill();
            gepLeadingWTAConeBasicClustersEtaSKJetsTree->Fill();
            gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Fill();
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
            gFexMETJwoJSimTree->Fill();
            gFexMETNoiseCutSimTree->Fill();
            gFexMETRmsSimTree->Fill();
            jFexMETTree->Fill();
        } // loop through events
        std::cout << "for jz: " << jzSlice << " these many events passed: " << passedEventsCounter << " out of: " << event.getEntries() << "\n";
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
    outOfTimeAntiKt4TruthJetsTree->Write("", TObject::kOverwrite);
    leadingOutOfTimeAntiKt4TruthJetsTree->Write("", TObject::kOverwrite);
    subleadingOutOfTimeAntiKt4TruthJetsTree->Write("", TObject::kOverwrite);
    jFexSRJTree->Write("", TObject::kOverwrite);
    jFexLeadingSRJTree->Write("", TObject::kOverwrite);
    jFexSubleadingSRJTree->Write("", TObject::kOverwrite);
    jFexLRJTree->Write("", TObject::kOverwrite);
    jFexLeadingLRJTree->Write("", TObject::kOverwrite);
    jFexSubleadingLRJTree->Write("", TObject::kOverwrite);
    jFexSRJSimTree->Write("", TObject::kOverwrite);
    jFexLeadingSRJSimTree->Write("", TObject::kOverwrite);
    jFexSubleadingSRJSimTree->Write("", TObject::kOverwrite);
    gFexLRJSimTree->Write("", TObject::kOverwrite);
    gFexLeadingLRJSimTree->Write("", TObject::kOverwrite);
    gFexSubleadingLRJSimTree->Write("", TObject::kOverwrite);
    gFexMETJwoJSimTree->Write("", TObject::kOverwrite);
    gFexMETNoiseCutSimTree->Write("", TObject::kOverwrite);
    gFexMETRmsSimTree->Write("", TObject::kOverwrite);
    jFexMETTree->Write("", TObject::kOverwrite);
    // TODO: uncomment once L1_gFexSRJetRoISim branch exists in GEP ntuples
    // gFexSRJSimTree->Write("", TObject::kOverwrite);
    // gFexLeadingSRJSimTree->Write("", TObject::kOverwrite);
    // gFexSubleadingSRJSimTree->Write("", TObject::kOverwrite);
    gepCellsTowersEtaSKTree->Write("", TObject::kOverwrite);
    gepBasicClustersEtaSKTree->Write("", TObject::kOverwrite);
    gepWTAConeCellsTowersEtaSKJetsTree->Write("", TObject::kOverwrite);
    gepLeadingWTAConeCellsTowersEtaSKJetsTree->Write("", TObject::kOverwrite);
    gepSubleadingWTAConeCellsTowersEtaSKJetsTree->Write("", TObject::kOverwrite);
    gepWTAConeBasicClustersEtaSKJetsTree->Write("", TObject::kOverwrite);
    gepLeadingWTAConeBasicClustersEtaSKJetsTree->Write("", TObject::kOverwrite);
    gepSubleadingWTAConeBasicClustersEtaSKJetsTree->Write("", TObject::kOverwrite);
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
        nTupler(true, "ggF_hh_bbbb_resim", etaAltRange, algoVersion); // call for signal (hh->4b ggF) (resimulated objects)
        nTupler(true, "ZvvHbb_resim", etaAltRange, algoVersion); // call for signal Z -> nu nu H -> bb (resimulated objects)
        //nTupler(true, "ZvvHbb", etaAltRange, algoVersion); // call for signal Z -> nu nu H -> bb
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
            if(jzSlice == 0){
                nTupler(false, "", etaAltRange, algoVersion, jzSlice, true); // call for each background jzSlice
            }
            TString out = makeOutputFileName("", false, jzSlice, algoVersion, false);
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

// Entry point for single-file Condor jobs.
// Called as: root -b -q 'HERNTupler.C+(signalBool, "signalString", algoVersion, jzSlice, specialJZ0, "outputDir", "daodFile", "gepFile", "fileSuffix")'
// etaAltRange is derived automatically from algoVersion.
void HERNTupler(bool signalBool, std::string signalString, unsigned int algoVersion,
                unsigned int jzSlice, bool specialJZ0,
                std::string outputDir, std::string daodFile, std::string gepFile,
                std::string fileSuffix) {
    unsigned int etaAltRange = 0;
    if (algoVersion == 2) etaAltRange = 784;  // 9.8 / 0.0125
    else if (algoVersion == 3) etaAltRange = 98;   // 9.8 / 0.1

    nTupler(signalBool, signalString, etaAltRange, algoVersion,
            jzSlice, specialJZ0,
            outputDir, /*writeDatFiles=*/false,
            daodFile, gepFile, fileSuffix);
    gSystem->Exit(0);
}
