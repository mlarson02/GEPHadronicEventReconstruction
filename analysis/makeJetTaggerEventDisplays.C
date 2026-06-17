#include "analysisHelperFunctions.h"

void callMakeEventDisplays(std::string jetTaggerInputFile, std::string herInputFile, unsigned int desiredJZSlice, bool signalBool, std::string signalString, std::string subjetType, bool skBool){
    SetPlotStyle();

    auto fileInfo = ParseFileName(jetTaggerInputFile);
    std::string inputObjectType = fileInfo.inputObjectType;
    std::string seedObjectType = fileInfo.seedObjectType;
    std::cout << "seedObjectType: " << seedObjectType << "\n";
    std::string rMergeValue = fileInfo.rMergeValue;
    double jetRadius = sqrt(ParseFileName(jetTaggerInputFile).jetRadiusSquared);

    unsigned int nInputObjectsAlgorithmConfiguration = 0;
    std::regex re("_IOs_(\\d+)_");
    std::smatch match;

    if (std::regex_search(jetTaggerInputFile, match, re)) {
        nInputObjectsAlgorithmConfiguration = std::stoi(match[1]);  // Convert to integer
        std::cout << "Extracted number of input objects for this algorithm configuration: " << nInputObjectsAlgorithmConfiguration << std::endl;
    } else {
        std::cout << "No match found." << std::endl;
    }
    std::vector<std::string > algorithmConfigurations;
    std::string key = "16130_";
    size_t pos = jetTaggerInputFile.find(key);

    if (pos != std::string::npos) {
        size_t start = pos + key.length();
        size_t end =  jetTaggerInputFile.rfind(".root");

        if (end != std::string::npos && end > start) {
            algorithmConfigurations.push_back( jetTaggerInputFile.substr(start, end - start));
        } else {
            std::cerr << ".root not found or too early in string.\n";
        }
    } else {
        std::cerr << "\"16130_\" not found in jet tagger file name.\n";
    }
    
   
    TString modifiedOutputFileDir = "eventDisplays/" + signalString + "_" + algorithmConfigurations[0] + "_seed" + seedObjectType + "_" + subjetType + "subjets/";
    gSystem->mkdir(modifiedOutputFileDir);

    std::vector<double>* mcEventWeightsValues = nullptr;
    double sumOfWeightsForSampleValues = 0.0;
    std::vector<double>* eventWeightsValues = nullptr;
    int sampleJZSliceValues = -1;
    std::vector<double>* jetTaggerLRJPsi_RValues = nullptr;
    std::vector<unsigned int>* jetTaggerLRJMergedIndicesValues = nullptr;
    std::vector<double>* jetTaggerLRJEtValues = nullptr;
    std::vector<double>* jetTaggerLRJEtaValues = nullptr;
    std::vector<double>* jetTaggerLRJPhiValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJPsi_RValues = nullptr;
    std::vector<unsigned int>* jetTaggerLeadingLRJMergedIndicesValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJEtValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJEtaValues = nullptr;
    std::vector<double>* jetTaggerLeadingLRJPhiValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJPsi_RValues = nullptr;
    std::vector<unsigned int>* jetTaggerSubleadingLRJMergedIndicesValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJEtValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJEtaValues = nullptr;
    std::vector<double>* jetTaggerSubleadingLRJPhiValues = nullptr;
    std::vector<unsigned int>* higgsIndexValues = nullptr;
    std::vector<unsigned int>* indexOfHiggsValues = nullptr;
    std::vector<double>* truthbquarksEtValues = nullptr;
    std::vector<double>* truthbquarksEnergyValues = nullptr;
    std::vector<double>* truthbquarkspTValues = nullptr;
    std::vector<double>* truthbquarkspxValues = nullptr;
    std::vector<double>* truthbquarkspyValues = nullptr;
    std::vector<double>* truthbquarkspzValues = nullptr;
    std::vector<double>* truthbquarksEtaValues = nullptr;
    std::vector<double>* truthbquarksPhiValues = nullptr;
    std::vector<double>* truthHiggsEtValues = nullptr;
    std::vector<double>* truthHiggsEnergyValues = nullptr;
    std::vector<double>* truthHiggspTValues = nullptr;
    std::vector<double>* truthHiggspxValues = nullptr;
    std::vector<double>* truthHiggspyValues = nullptr;
    std::vector<double>* truthHiggspzValues = nullptr;
    std::vector<double>* truthHiggsEtaValues = nullptr;
    std::vector<double>* truthHiggsPhiValues = nullptr;
    std::vector<double>* truthHiggsInvMassValues = nullptr;
    std::vector<double>* gepCellsTowersEtValues = nullptr;
    std::vector<double>* gepCellsTowersEtaValues = nullptr;
    std::vector<double>* gepCellsTowersPhiValues = nullptr;
    std::vector<double>* gepCellsTowersSKEtValues = nullptr;
    std::vector<double>* gepCellsTowersSKEtaValues = nullptr;
    std::vector<double>* gepCellsTowersSKPhiValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersJetspTValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersJetsEtaValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersJetsPhiValues = nullptr;
    std::vector<unsigned int >* gepWTAConeCellsTowersJetsNConstituentsValues = nullptr;
    std::vector<double>* gepLeadingWTAConeCellsTowersJetspTValues = nullptr;
    std::vector<double>* gepLeadingWTAConeCellsTowersJetsEtaValues = nullptr;
    std::vector<double>* gepLeadingWTAConeCellsTowersJetsPhiValues = nullptr;
    std::vector<unsigned int >* gepLeadingWTAConeCellsTowersJetsNConstituentsValues = nullptr;
    std::vector<double>* gepSubleadingWTAConeCellsTowersJetspTValues = nullptr;
    std::vector<double>* gepSubleadingWTAConeCellsTowersJetsEtaValues = nullptr;
    std::vector<double>* gepSubleadingWTAConeCellsTowersJetsPhiValues = nullptr;
    std::vector<unsigned int >* gepSubleadingWTAConeCellsTowersJetsNConstituentsValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersSKJetspTValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersSKJetsEtaValues = nullptr;
    std::vector<double>* gepWTAConeCellsTowersSKJetsPhiValues = nullptr;
    std::vector<unsigned int >* gepWTAConeCellsTowersSKJetsNConstituentsValues = nullptr;
    std::vector<double>* gepLeadingWTAConeCellsTowersSKJetspTValues = nullptr;
    std::vector<double>* gepLeadingWTAConeCellsTowersSKJetsEtaValues = nullptr;
    std::vector<double>* gepLeadingWTAConeCellsTowersSKJetsPhiValues = nullptr;
    std::vector<unsigned int >* gepLeadingWTAConeCellsTowersSKJetsNConstituentsValues = nullptr;
    std::vector<double>* gepSubleadingWTAConeCellsTowersSKJetspTValues = nullptr;
    std::vector<double>* gepSubleadingWTAConeCellsTowersSKJetsEtaValues = nullptr;
    std::vector<double>* gepSubleadingWTAConeCellsTowersSKJetsPhiValues = nullptr;
    std::vector<unsigned int >* gepSubleadingWTAConeCellsTowersSKJetsNConstituentsValues = nullptr;
    std::vector<unsigned int>* gFexSRJEtIndexValues = nullptr;
    std::vector<double>* gFexSRJEtValues = nullptr;
    std::vector<double>* gFexSRJEtaValues = nullptr;
    std::vector<double>* gFexSRJPhiValues = nullptr;
    std::vector<double>* gFexSRJLeadingEtValues = nullptr;
    std::vector<double>* gFexSRJLeadingEtaValues = nullptr;
    std::vector<double>* gFexSRJLeadingPhiValues = nullptr;
    std::vector<double>* gFexSRJSubleadingEtValues = nullptr;
    std::vector<double>* gFexSRJSubleadingEtaValues = nullptr;
    std::vector<double>* gFexSRJSubleadingPhiValues = nullptr;
    std::vector<unsigned int>* gFexLRJEtIndexValues = nullptr;
    std::vector<double>* gFexLRJEtValues = nullptr;
    std::vector<double>* gFexLRJEtaValues = nullptr;
    std::vector<double>* gFexLRJPhiValues = nullptr;
    std::vector<double>* gFexLRJLeadingEtValues = nullptr;
    std::vector<double>* gFexLRJLeadingEtaValues = nullptr;
    std::vector<double>* gFexLRJLeadingPhiValues = nullptr;
    std::vector<double>* gFexLRJSubleadingEtValues = nullptr;
    std::vector<double>* gFexLRJSubleadingEtaValues = nullptr;
    std::vector<double>* gFexLRJSubleadingPhiValues = nullptr;
    std::vector<unsigned int>* jFexSRJEtIndexValues = nullptr;
    std::vector<double>* jFexSRJEtValues = nullptr;
    std::vector<double>* jFexSRJEtaValues = nullptr;
    std::vector<double>* jFexSRJPhiValues = nullptr;
    std::vector<double>* jFexSRJLeadingEtValues = nullptr;
    std::vector<double>* jFexSRJLeadingEtaValues = nullptr;
    std::vector<double>* jFexSRJLeadingPhiValues = nullptr;
    std::vector<double>* jFexSRJSubleadingEtValues = nullptr;
    std::vector<double>* jFexSRJSubleadingEtaValues = nullptr;
    std::vector<double>* jFexSRJSubleadingPhiValues = nullptr;
    std::vector<unsigned int>* recoAntiKt10LRJEtIndexValues = nullptr;
    std::vector<double>* recoAntiKt10LRJEtValues = nullptr;
    std::vector<double>* recoAntiKt10LRJEtaValues = nullptr;
    std::vector<double>* recoAntiKt10LRJPhiValues = nullptr;
    std::vector<double>* recoAntiKt10LRJMassValues = nullptr;
    std::vector<double>* recoAntiKt10LRJLeadingEtValues = nullptr;
    std::vector<double>* recoAntiKt10LRJLeadingEtaValues = nullptr;
    std::vector<double>* recoAntiKt10LRJLeadingPhiValues = nullptr;
    std::vector<double>* recoAntiKt10LRJLeadingMassValues = nullptr;
    std::vector<double>* recoAntiKt10LRJSubleadingEtValues = nullptr;
    std::vector<double>* recoAntiKt10LRJSubleadingEtaValues = nullptr;
    std::vector<double>* recoAntiKt10LRJSubleadingPhiValues = nullptr;
    std::vector<double>* recoAntiKt10LRJSubleadingMassValues = nullptr;
    std::vector<unsigned int>* truthAntiKt4WZSRJEtIndexValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJEtValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJEtaValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJPhiValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJMassValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJLeadingEtValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJLeadingEtaValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJLeadingPhiValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJLeadingMassValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJSubleadingEtValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJSubleadingEtaValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJSubleadingPhiValues = nullptr;
    std::vector<double>* truthAntiKt4WZSRJSubleadingMassValues = nullptr;
    std::vector<unsigned int>* inTimeAntiKt4TruthSRJEtIndexValues = nullptr;
    std::vector<double>* inTimeAntiKt4TruthSRJEtValues = nullptr;
    std::vector<double>* inTimeAntiKt4TruthSRJEtaValues = nullptr;
    std::vector<double>* inTimeAntiKt4TruthSRJPhiValues = nullptr;
    std::vector<double>* inTimeAntiKt4TruthSRJLeadingEtValues = nullptr;
    std::vector<double>* inTimeAntiKt4TruthSRJLeadingEtaValues = nullptr;
    std::vector<double>* inTimeAntiKt4TruthSRJLeadingPhiValues = nullptr;
    std::vector<double>* inTimeAntiKt4TruthSRJSubleadingEtValues = nullptr;
    std::vector<double>* inTimeAntiKt4TruthSRJSubleadingEtaValues = nullptr;
    std::vector<double>* inTimeAntiKt4TruthSRJSubleadingPhiValues = nullptr;
    std::vector<double>* outOfTimeAntiKt4TruthSRJEtValues = nullptr;
    std::vector<double>* outOfTimeAntiKt4TruthSRJEtaValues = nullptr;
    std::vector<double>* outOfTimeAntiKt4TruthSRJPhiValues = nullptr;
    std::vector<double>* outOfTimeAntiKt4TruthSRJLeadingEtValues = nullptr;
    std::vector<double>* outOfTimeAntiKt4TruthSRJLeadingEtaValues = nullptr;
    std::vector<double>* outOfTimeAntiKt4TruthSRJLeadingPhiValues = nullptr;

    TFile* jetTaggerFile = TFile::Open(jetTaggerInputFile.c_str(), "READ");
    TFile* herFile = TFile::Open(herInputFile.c_str(), "READ");

    // jet tagger trees — from jet tagger NTuple
    TTree* jetTaggerLRJs = (TTree*)jetTaggerFile->Get("jetTaggerLRJsTree");
    TTree* jetTaggerLeadingLRJs = (TTree*)jetTaggerFile->Get("jetTaggerLeadingLRJsTree");
    TTree* jetTaggerSubleadingLRJs = (TTree*)jetTaggerFile->Get("jetTaggerSubleadingLRJsTree");
    // all other trees — from HER input NTuple
    TTree* eventInfoTree = (TTree*)herFile->Get("eventInfoTree");
    TTree* truthbTree = (TTree*)herFile->Get("truthbTree");
    TTree* truthHiggsTree = (TTree*)herFile->Get("truthHiggsTree");
    TTree* gepCellsTowersTree = (TTree*)herFile->Get("gepCellsTowersTree");
    TTree* gepCellsTowersSKTree = (TTree*)herFile->Get("gepCellsTowersSKTree");
    TTree* gepWTAConeCellsTowersJetsTree = (TTree*)herFile->Get("gepWTAConeCellsTowersJetsTree");
    TTree* gepLeadingWTAConeCellsTowersJetsTree = (TTree*)herFile->Get("gepLeadingWTAConeCellsTowersJetsTree");
    TTree* gepSubleadingWTAConeCellsTowersJetsTree = (TTree*)herFile->Get("gepSubleadingWTAConeCellsTowersJetsTree");
    TTree* gepWTAConeCellsTowersSKJetsTree = (TTree*)herFile->Get("gepWTAConeCellsTowersSKJetsTree");
    TTree* gepLeadingWTAConeCellsTowersSKJetsTree = (TTree*)herFile->Get("gepLeadingWTAConeCellsTowersSKJetsTree");
    TTree* gepSubleadingWTAConeCellsTowersSKJetsTree = (TTree*)herFile->Get("gepSubleadingWTAConeCellsTowersSKJetsTree");
    TTree* gFexSRJTree = (TTree*)herFile->Get("gFexSRJTree");
    TTree* gFexLeadingSRJTree = (TTree*)herFile->Get("gFexLeadingSRJTree");
    TTree* gFexSubleadingSRJTree = (TTree*)herFile->Get("gFexSubleadingSRJTree");
    TTree* gFexLRJTree = (TTree*)herFile->Get("gFexLRJTree");
    TTree* gFexLeadingLRJTree = (TTree*)herFile->Get("gFexLeadingLRJTree");
    TTree* gFexSubleadingLRJTree = (TTree*)herFile->Get("gFexSubleadingLRJTree");
    TTree* inTimeAntiKt4TruthJetsTree = (TTree*)herFile->Get("inTimeAntiKt4TruthJetsTree");
    TTree* leadingInTimeAntiKt4TruthJetsTree = (TTree*)herFile->Get("leadingInTimeAntiKt4TruthJetsTree");
    TTree* subleadingInTimeAntiKt4TruthJetsTree = (TTree*)herFile->Get("subleadingInTimeAntiKt4TruthJetsTree");
    TTree* outOfTimeAntiKt4TruthJetsTree = (TTree*)herFile->Get("outOfTimeAntiKt4TruthJetsTree");
    TTree* leadingOutOfTimeAntiKt4TruthJetsTree = (TTree*)herFile->Get("leadingOutOfTimeAntiKt4TruthJetsTree");
    TTree* jFexSRJTree = (TTree*)herFile->Get("jFexSRJTree");
    TTree* jFexLeadingSRJTree = (TTree*)herFile->Get("jFexLeadingSRJTree");
    TTree* jFexSubleadingSRJTree = (TTree*)herFile->Get("jFexSubleadingSRJTree");
    TTree* recoAntiKt10UFOCSSKJets = (TTree*)herFile->Get("recoAntiKt10UFOCSSKJets");
    TTree* leadingRecoAntiKt10UFOCSSKJets = (TTree*)herFile->Get("leadingRecoAntiKt10UFOCSSKJets");
    TTree* subleadingRecoAntiKt10UFOCSSKJets = (TTree*)herFile->Get("subleadingRecoAntiKt10UFOCSSKJets");
    TTree* truthAntiKt4TruthDressedWZJets = (TTree*)herFile->Get("truthAntiKt4TruthDressedWZJets");
    TTree* leadingTruthAntiKt4TruthDressedWZJets = (TTree*)herFile->Get("leadingTruthAntiKt4TruthDressedWZJets");
    TTree* subleadingTruthAntiKt4TruthDressedWZJets = (TTree*)herFile->Get("subleadingTruthAntiKt4TruthDressedWZJets");

    // === eventInfoTree ===
    eventInfoTree->SetBranchAddress("mcEventWeight", &mcEventWeightsValues);
    eventInfoTree->SetBranchAddress("sumOfWeightsForSample", &sumOfWeightsForSampleValues);
    eventInfoTree->SetBranchAddress("eventWeights", &eventWeightsValues);
    eventInfoTree->SetBranchAddress("sampleJZSlice", &sampleJZSliceValues);

    // === jetTaggerLRJs ===
    jetTaggerLRJs->SetBranchAddress("Psi_R", &jetTaggerLRJPsi_RValues);
    jetTaggerLRJs->SetBranchAddress("MergedIndices", &jetTaggerLRJMergedIndicesValues);
    jetTaggerLRJs->SetBranchAddress("Et", &jetTaggerLRJEtValues);
    jetTaggerLRJs->SetBranchAddress("Eta", &jetTaggerLRJEtaValues);
    jetTaggerLRJs->SetBranchAddress("Phi", &jetTaggerLRJPhiValues);

    // === jetTaggerLeadingLRJs ===
    jetTaggerLeadingLRJs->SetBranchAddress("Psi_R", &jetTaggerLeadingLRJPsi_RValues);
    jetTaggerLeadingLRJs->SetBranchAddress("MergedIndices", &jetTaggerLeadingLRJMergedIndicesValues);
    jetTaggerLeadingLRJs->SetBranchAddress("Et", &jetTaggerLeadingLRJEtValues);
    jetTaggerLeadingLRJs->SetBranchAddress("Eta", &jetTaggerLeadingLRJEtaValues);
    jetTaggerLeadingLRJs->SetBranchAddress("Phi", &jetTaggerLeadingLRJPhiValues);

    // === jetTaggerSubleadingLRJs ===
    jetTaggerSubleadingLRJs->SetBranchAddress("Psi_R", &jetTaggerSubleadingLRJPsi_RValues);
    jetTaggerSubleadingLRJs->SetBranchAddress("MergedIndices", &jetTaggerSubleadingLRJMergedIndicesValues);
    jetTaggerSubleadingLRJs->SetBranchAddress("Et", &jetTaggerSubleadingLRJEtValues);
    jetTaggerSubleadingLRJs->SetBranchAddress("Eta", &jetTaggerSubleadingLRJEtaValues);
    jetTaggerSubleadingLRJs->SetBranchAddress("Phi", &jetTaggerSubleadingLRJPhiValues);

    // === truthbTree ===
    truthbTree->SetBranchAddress("higgsIndex", &higgsIndexValues);
    truthbTree->SetBranchAddress("Et", &truthbquarksEtValues);
    truthbTree->SetBranchAddress("Energy", &truthbquarksEnergyValues);
    truthbTree->SetBranchAddress("Et", &truthbquarkspTValues);
    truthbTree->SetBranchAddress("px", &truthbquarkspxValues);
    truthbTree->SetBranchAddress("py", &truthbquarkspyValues);
    truthbTree->SetBranchAddress("pz", &truthbquarkspzValues);
    truthbTree->SetBranchAddress("Eta", &truthbquarksEtaValues);
    truthbTree->SetBranchAddress("Phi", &truthbquarksPhiValues);

    // === truthHiggsTree ===
    truthHiggsTree->SetBranchAddress("indexOfHiggs", &indexOfHiggsValues);
    truthHiggsTree->SetBranchAddress("invMass", &truthHiggsInvMassValues);
    truthHiggsTree->SetBranchAddress("Et", &truthHiggsEtValues);
    truthHiggsTree->SetBranchAddress("Energy", &truthHiggsEnergyValues);
    truthHiggsTree->SetBranchAddress("Et", &truthHiggspTValues);
    truthHiggsTree->SetBranchAddress("px", &truthHiggspxValues);
    truthHiggsTree->SetBranchAddress("py", &truthHiggspyValues);
    truthHiggsTree->SetBranchAddress("pz", &truthHiggspzValues);
    truthHiggsTree->SetBranchAddress("Eta", &truthHiggsEtaValues);
    truthHiggsTree->SetBranchAddress("Phi", &truthHiggsPhiValues);

    // === gepCellsTower tree ===
    gepCellsTowersTree->SetBranchAddress("Et", &gepCellsTowersEtValues);
    gepCellsTowersTree->SetBranchAddress("Eta", &gepCellsTowersEtaValues);
    gepCellsTowersTree->SetBranchAddress("Phi", &gepCellsTowersPhiValues);

    // === gepCellsTowerSK tree ===
    gepCellsTowersSKTree->SetBranchAddress("Et", &gepCellsTowersSKEtValues);
    gepCellsTowersSKTree->SetBranchAddress("Eta", &gepCellsTowersSKEtaValues);
    gepCellsTowersSKTree->SetBranchAddress("Phi", &gepCellsTowersSKPhiValues);

    // gep wta cone cells towers jets
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Et", &gepWTAConeCellsTowersJetspTValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Eta", &gepWTAConeCellsTowersJetsEtaValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("Phi", &gepWTAConeCellsTowersJetsPhiValues);
    gepWTAConeCellsTowersJetsTree->SetBranchAddress("NConstituents", &gepWTAConeCellsTowersJetsNConstituentsValues);

    // leading gep wta cone cells towers jets
    gepLeadingWTAConeCellsTowersJetsTree->SetBranchAddress("Et", &gepLeadingWTAConeCellsTowersJetspTValues);
    gepLeadingWTAConeCellsTowersJetsTree->SetBranchAddress("Eta", &gepLeadingWTAConeCellsTowersJetsEtaValues);
    gepLeadingWTAConeCellsTowersJetsTree->SetBranchAddress("Phi", &gepLeadingWTAConeCellsTowersJetsPhiValues);
    gepLeadingWTAConeCellsTowersJetsTree->SetBranchAddress("NConstituents", &gepLeadingWTAConeCellsTowersJetsNConstituentsValues);

    // subleading gep wta cone cells towers jets
    gepSubleadingWTAConeCellsTowersJetsTree->SetBranchAddress("Et", &gepSubleadingWTAConeCellsTowersJetspTValues);
    gepSubleadingWTAConeCellsTowersJetsTree->SetBranchAddress("Eta", &gepSubleadingWTAConeCellsTowersJetsEtaValues);
    gepSubleadingWTAConeCellsTowersJetsTree->SetBranchAddress("Phi", &gepSubleadingWTAConeCellsTowersJetsPhiValues);
    gepSubleadingWTAConeCellsTowersJetsTree->SetBranchAddress("NConstituents", &gepSubleadingWTAConeCellsTowersJetsNConstituentsValues);

    // gep wta cone cells towers SK jets
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Et", &gepWTAConeCellsTowersSKJetspTValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Eta", &gepWTAConeCellsTowersSKJetsEtaValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("Phi", &gepWTAConeCellsTowersSKJetsPhiValues);
    gepWTAConeCellsTowersSKJetsTree->SetBranchAddress("NConstituents", &gepWTAConeCellsTowersSKJetsNConstituentsValues);

    // leading gep wta cone cells towers SK jets
    gepLeadingWTAConeCellsTowersSKJetsTree->SetBranchAddress("Et", &gepLeadingWTAConeCellsTowersSKJetspTValues);
    gepLeadingWTAConeCellsTowersSKJetsTree->SetBranchAddress("Eta", &gepLeadingWTAConeCellsTowersSKJetsEtaValues);
    gepLeadingWTAConeCellsTowersSKJetsTree->SetBranchAddress("Phi", &gepLeadingWTAConeCellsTowersSKJetsPhiValues);
    gepLeadingWTAConeCellsTowersSKJetsTree->SetBranchAddress("NConstituents", &gepLeadingWTAConeCellsTowersSKJetsNConstituentsValues);

    // subleading gep wta cone cells towers SK jets
    gepSubleadingWTAConeCellsTowersSKJetsTree->SetBranchAddress("Et", &gepSubleadingWTAConeCellsTowersSKJetspTValues);
    gepSubleadingWTAConeCellsTowersSKJetsTree->SetBranchAddress("Eta", &gepSubleadingWTAConeCellsTowersSKJetsEtaValues);
    gepSubleadingWTAConeCellsTowersSKJetsTree->SetBranchAddress("Phi", &gepSubleadingWTAConeCellsTowersSKJetsPhiValues);
    gepSubleadingWTAConeCellsTowersSKJetsTree->SetBranchAddress("NConstituents", &gepSubleadingWTAConeCellsTowersSKJetsNConstituentsValues);

    // === gFexSRJTree ===
    gFexSRJTree->SetBranchAddress("EtIndex", &gFexSRJEtIndexValues);
    gFexSRJTree->SetBranchAddress("Et", &gFexSRJEtValues);
    gFexSRJTree->SetBranchAddress("Eta", &gFexSRJEtaValues);
    gFexSRJTree->SetBranchAddress("Phi", &gFexSRJPhiValues);

    // === gFexLeadingSRJTree ===
    gFexLeadingSRJTree->SetBranchAddress("Et", &gFexSRJLeadingEtValues);
    gFexLeadingSRJTree->SetBranchAddress("Eta", &gFexSRJLeadingEtaValues);
    gFexLeadingSRJTree->SetBranchAddress("Phi", &gFexSRJLeadingPhiValues);

    // === gFexSubleadingSRJTree ===
    gFexSubleadingSRJTree->SetBranchAddress("Et", &gFexSRJSubleadingEtValues);
    gFexSubleadingSRJTree->SetBranchAddress("Eta", &gFexSRJSubleadingEtaValues);
    gFexSubleadingSRJTree->SetBranchAddress("Phi", &gFexSRJSubleadingPhiValues);

    // === gFexLRJTree ===
    gFexLRJTree->SetBranchAddress("EtIndex", &gFexLRJEtIndexValues);
    gFexLRJTree->SetBranchAddress("Et", &gFexLRJEtValues);
    gFexLRJTree->SetBranchAddress("Eta", &gFexLRJEtaValues);
    gFexLRJTree->SetBranchAddress("Phi", &gFexLRJPhiValues);

    // === gFexLeadingLRJTree ===
    gFexLeadingLRJTree->SetBranchAddress("Et", &gFexLRJLeadingEtValues);
    gFexLeadingLRJTree->SetBranchAddress("Eta", &gFexLRJLeadingEtaValues);
    gFexLeadingLRJTree->SetBranchAddress("Phi", &gFexLRJLeadingPhiValues);

    // === gFexSubleadingLRJTree ===
    gFexSubleadingLRJTree->SetBranchAddress("Et", &gFexLRJSubleadingEtValues);
    gFexSubleadingLRJTree->SetBranchAddress("Eta", &gFexLRJSubleadingEtaValues);
    gFexSubleadingLRJTree->SetBranchAddress("Phi", &gFexLRJSubleadingPhiValues);

    // === jFexSRJTree ===
    jFexSRJTree->SetBranchAddress("EtIndex", &jFexSRJEtIndexValues);
    jFexSRJTree->SetBranchAddress("Et", &jFexSRJEtValues);
    jFexSRJTree->SetBranchAddress("Eta", &jFexSRJEtaValues);
    jFexSRJTree->SetBranchAddress("Phi", &jFexSRJPhiValues);

    // === jFexLeadingSRJTree ===
    jFexLeadingSRJTree->SetBranchAddress("Et", &jFexSRJLeadingEtValues);
    jFexLeadingSRJTree->SetBranchAddress("Eta", &jFexSRJLeadingEtaValues);
    jFexLeadingSRJTree->SetBranchAddress("Phi", &jFexSRJLeadingPhiValues);

    // === jFexSubleadingSRJTree ===
    jFexSubleadingSRJTree->SetBranchAddress("Et", &jFexSRJSubleadingEtValues);
    jFexSubleadingSRJTree->SetBranchAddress("Eta", &jFexSRJSubleadingEtaValues);
    jFexSubleadingSRJTree->SetBranchAddress("Phi", &jFexSRJSubleadingPhiValues);


    // === recoAntiKt10UFOCSSKJets ===
    recoAntiKt10UFOCSSKJets->SetBranchAddress("EtIndex", &recoAntiKt10LRJEtIndexValues);
    recoAntiKt10UFOCSSKJets->SetBranchAddress("Et", &recoAntiKt10LRJEtValues);
    recoAntiKt10UFOCSSKJets->SetBranchAddress("Eta", &recoAntiKt10LRJEtaValues);
    recoAntiKt10UFOCSSKJets->SetBranchAddress("Phi", &recoAntiKt10LRJPhiValues);
    recoAntiKt10UFOCSSKJets->SetBranchAddress("Phi", &recoAntiKt10LRJMassValues);

    // === leadingRecoAntiKt10UFOCSSKJets ===
    leadingRecoAntiKt10UFOCSSKJets->SetBranchAddress("Et", &recoAntiKt10LRJLeadingEtValues);
    leadingRecoAntiKt10UFOCSSKJets->SetBranchAddress("Eta", &recoAntiKt10LRJLeadingEtaValues);
    leadingRecoAntiKt10UFOCSSKJets->SetBranchAddress("Phi", &recoAntiKt10LRJLeadingPhiValues);
    leadingRecoAntiKt10UFOCSSKJets->SetBranchAddress("Mass", &recoAntiKt10LRJLeadingMassValues);

    // === subleadingRecoAntiKt10UFOCSSKJets ===
    subleadingRecoAntiKt10UFOCSSKJets->SetBranchAddress("Et", &recoAntiKt10LRJSubleadingEtValues);
    subleadingRecoAntiKt10UFOCSSKJets->SetBranchAddress("Eta", &recoAntiKt10LRJSubleadingEtaValues);
    subleadingRecoAntiKt10UFOCSSKJets->SetBranchAddress("Phi", &recoAntiKt10LRJSubleadingPhiValues);
    subleadingRecoAntiKt10UFOCSSKJets->SetBranchAddress("Mass", &recoAntiKt10LRJSubleadingMassValues);

    // === truthAntiKt4TruthDressedWZJets ===
    truthAntiKt4TruthDressedWZJets->SetBranchAddress("EtIndex", &truthAntiKt4WZSRJEtIndexValues);
    truthAntiKt4TruthDressedWZJets->SetBranchAddress("Et", &truthAntiKt4WZSRJEtValues);
    truthAntiKt4TruthDressedWZJets->SetBranchAddress("Eta", &truthAntiKt4WZSRJEtaValues);
    truthAntiKt4TruthDressedWZJets->SetBranchAddress("Phi", &truthAntiKt4WZSRJPhiValues);
    truthAntiKt4TruthDressedWZJets->SetBranchAddress("Mass", &truthAntiKt4WZSRJMassValues);

    // === leadingTruthAntiKt4TruthDressedWZJets ===
    leadingTruthAntiKt4TruthDressedWZJets->SetBranchAddress("Et", &truthAntiKt4WZSRJLeadingEtValues);
    leadingTruthAntiKt4TruthDressedWZJets->SetBranchAddress("Eta", &truthAntiKt4WZSRJLeadingEtaValues);
    leadingTruthAntiKt4TruthDressedWZJets->SetBranchAddress("Phi", &truthAntiKt4WZSRJLeadingPhiValues);
    leadingTruthAntiKt4TruthDressedWZJets->SetBranchAddress("Mass", &truthAntiKt4WZSRJLeadingMassValues);

    // === subleadingTruthAntiKt4TruthDressedWZJets ===
    subleadingTruthAntiKt4TruthDressedWZJets->SetBranchAddress("Et", &truthAntiKt4WZSRJSubleadingEtValues);
    subleadingTruthAntiKt4TruthDressedWZJets->SetBranchAddress("Eta", &truthAntiKt4WZSRJSubleadingEtaValues);
    subleadingTruthAntiKt4TruthDressedWZJets->SetBranchAddress("Phi", &truthAntiKt4WZSRJSubleadingPhiValues);
    subleadingTruthAntiKt4TruthDressedWZJets->SetBranchAddress("Mass", &truthAntiKt4WZSRJSubleadingMassValues);

    // === inTimeAntiKt4TruthJetsTree ===
    inTimeAntiKt4TruthJetsTree->SetBranchAddress("EtIndex", &inTimeAntiKt4TruthSRJEtIndexValues);
    inTimeAntiKt4TruthJetsTree->SetBranchAddress("Et", &inTimeAntiKt4TruthSRJEtValues);
    inTimeAntiKt4TruthJetsTree->SetBranchAddress("Eta", &inTimeAntiKt4TruthSRJEtaValues);
    inTimeAntiKt4TruthJetsTree->SetBranchAddress("Phi", &inTimeAntiKt4TruthSRJPhiValues);

    // === leadingInTimeAntiKt4TruthJetsTree ===
    leadingInTimeAntiKt4TruthJetsTree->SetBranchAddress("Et", &inTimeAntiKt4TruthSRJLeadingEtValues);
    leadingInTimeAntiKt4TruthJetsTree->SetBranchAddress("Eta", &inTimeAntiKt4TruthSRJLeadingEtaValues);
    leadingInTimeAntiKt4TruthJetsTree->SetBranchAddress("Phi", &inTimeAntiKt4TruthSRJLeadingPhiValues);

    // === subleadingInTimeAntiKt4TruthJetsTree ===
    subleadingInTimeAntiKt4TruthJetsTree->SetBranchAddress("Et", &inTimeAntiKt4TruthSRJSubleadingEtValues);
    subleadingInTimeAntiKt4TruthJetsTree->SetBranchAddress("Eta", &inTimeAntiKt4TruthSRJSubleadingEtaValues);
    subleadingInTimeAntiKt4TruthJetsTree->SetBranchAddress("Phi", &inTimeAntiKt4TruthSRJSubleadingPhiValues);

    // === outOfTimeAntiKt4TruthJetsTree ===
    outOfTimeAntiKt4TruthJetsTree->SetBranchAddress("Et", &outOfTimeAntiKt4TruthSRJEtValues);
    outOfTimeAntiKt4TruthJetsTree->SetBranchAddress("Eta", &outOfTimeAntiKt4TruthSRJEtaValues);
    outOfTimeAntiKt4TruthJetsTree->SetBranchAddress("Phi", &outOfTimeAntiKt4TruthSRJPhiValues);

    // === leadingOutOfTimeAntiKt4TruthJetsTree ===
    leadingOutOfTimeAntiKt4TruthJetsTree->SetBranchAddress("Et", &outOfTimeAntiKt4TruthSRJLeadingEtValues);
    leadingOutOfTimeAntiKt4TruthJetsTree->SetBranchAddress("Eta", &outOfTimeAntiKt4TruthSRJLeadingEtaValues);
    leadingOutOfTimeAntiKt4TruthJetsTree->SetBranchAddress("Phi", &outOfTimeAntiKt4TruthSRJLeadingPhiValues);

    //const int numProcessedEvents = eventInfoTree->GetEntries();
    TString pdf = modifiedOutputFileDir + "cellsTowersHighestEtSeedPositions.pdf";
    TCanvas c("c","c",800,700);

    // open multipage pdf
    c.Print(pdf + "(");

    TString pdf_ED = modifiedOutputFileDir + "cellsTowersHighestEtTruthAndPileupJets.pdf";
    TCanvas cEventDisplay("cEventDisplay","cEventDisplay",800,700);

    // open multipage pdf_ED
    cEventDisplay.Print(pdf_ED + "(");
    //gErrorIgnoreLevel = kWarning;
    const unsigned int maxDisplays = 200;
    unsigned int nAccepted = 0;

    const int numProcessedEvents = eventInfoTree->GetEntries();

    for (int iEvt = 154796; iEvt < numProcessedEvents && nAccepted < maxDisplays; ++iEvt) {
        eventInfoTree->GetEntry(iEvt);
        //std::cout << "iEvt: " << iEvt << "\n";
        //std::cout << "sampleJZSliceValues: " << desiredJZSlice << "\n";
        if (!signalBool && sampleJZSliceValues != desiredJZSlice)
            continue;

        // ---- event accepted ----
        ++nAccepted;
        //std::cout << "accepted!" << "\n";

        if(iEvt % 10 == 0) std::cout << "iEvt: " << iEvt << "\n";
        inTimeAntiKt4TruthJetsTree->GetEntry(iEvt); // NOTE all of these are pT / Et sorted!
        leadingInTimeAntiKt4TruthJetsTree->GetEntry(iEvt);
        outOfTimeAntiKt4TruthJetsTree->GetEntry(iEvt);
        leadingOutOfTimeAntiKt4TruthJetsTree->GetEntry(iEvt);
        gepCellsTowersTree->GetEntry(iEvt);
        gepCellsTowersSKTree->GetEntry(iEvt);
        jetTaggerLeadingLRJs->GetEntry(iEvt);
        eventInfoTree->GetEntry(iEvt);
        jetTaggerSubleadingLRJs->GetEntry(iEvt);

        truthbTree->GetEntry(iEvt);
        truthHiggsTree->GetEntry(iEvt);
        gepWTAConeCellsTowersJetsTree->GetEntry(iEvt);
        gepWTAConeCellsTowersSKJetsTree->GetEntry(iEvt);
        jFexSRJTree->GetEntry(iEvt);
        gFexSRJTree->GetEntry(iEvt);
        gFexLRJTree->GetEntry(iEvt);
        recoAntiKt10UFOCSSKJets->GetEntry(iEvt);
        leadingRecoAntiKt10UFOCSSKJets->GetEntry(iEvt);
        truthAntiKt4TruthDressedWZJets->GetEntry(iEvt);
        leadingTruthAntiKt4TruthDressedWZJets->GetEntry(iEvt);

        TH2F *cellsTowersHighestEtSeedPositions = new TH2F("cellsTowersHighestEtSeedPositions", "Sum of CellsTowers E_{T} in Each Bin; #eta;#phi", 
                            100, -5, 5,   
                            64, -3.2, 3.2);  

        TH2F *cellsTowersHighestEtTruthAndPileupJets = new TH2F("cellsTowersHighestEtTruthAndPileupJets", "Sum of CellsTowers E_{T} in Each Bin; #eta;#phi", 
                            100, -5, 5,   
                            64, -3.2, 3.2);  
        
        std::vector<std::pair<double, double > > seedJetPositions;
        std::vector<std::pair<double, double > > additionalSeedJetPositions;
        std::vector<std::pair<double, double > > newSeedPositions;
        std::vector<double > lrjEt;
        std::vector<double > lrjPsi_R;
        std::array<double, 2 > lrjDeltaRSubjets;
        std::array<double, 2 > lrjEtRatioSubjets;
        std::array<double, 2 > lrjMassApprox;
        std::array<double, 2 > lrjTau21;
        
        //std::cout << "seedObjectType: " << seedObjectType << "\n";
        
        double subjetEtThreshold;
        if(subjetType == "WTACone") subjetEtThreshold = 25.0;
        else if(subjetType == "jFEXSRJ") subjetEtThreshold = 50.0;
        else if(subjetType == "gFEXSRJ") subjetEtThreshold = 35.0;
        const double truthJetEtThreshold = 15.0;
        const double pileupJetEtThreshold = 15.0;

        std::vector<std::pair<double, double > > possibleSubjets;
        std::vector<std::pair<double, double > > truthJets;
        std::vector<std::pair<double, double > > pileupJets;
        std::vector<std::pair<double, double > > outOfTimePileupJets;
        std::array< std::vector<std::pair<double, double > >, 2 > bQuarks;
        std::array<double, 2 > bQuarksDeltaR;
        std::vector<std::pair<double, double > > higgses;
        std::vector<double > higgsPts;
        std::array<std::vector<double >, 2 > bEts; 

        for(unsigned int iCellsTower = 0; iCellsTower < nInputObjectsAlgorithmConfiguration; iCellsTower++){
            if(skBool){
                cellsTowersHighestEtSeedPositions->Fill(gepCellsTowersSKEtaValues->at(iCellsTower), gepCellsTowersSKPhiValues->at(iCellsTower), gepCellsTowersSKEtValues->at(iCellsTower));
                cellsTowersHighestEtTruthAndPileupJets->Fill(gepCellsTowersSKEtaValues->at(iCellsTower), gepCellsTowersSKPhiValues->at(iCellsTower), gepCellsTowersSKEtValues->at(iCellsTower));
            }
            else{
                cellsTowersHighestEtSeedPositions->Fill(gepCellsTowersEtaValues->at(iCellsTower), gepCellsTowersPhiValues->at(iCellsTower), gepCellsTowersEtValues->at(iCellsTower));
                cellsTowersHighestEtTruthAndPileupJets->Fill(gepCellsTowersEtaValues->at(iCellsTower), gepCellsTowersPhiValues->at(iCellsTower), gepCellsTowersEtValues->at(iCellsTower));
            }
            
        }
        std::cout << "seedObjectType: " << seedObjectType << "\n";
        if(seedObjectType == "jFEXSRJ"){
            seedJetPositions.push_back(std::make_pair(jFexSRJEtaValues->at(0), jFexSRJPhiValues->at(0)));
            seedJetPositions.push_back(std::make_pair(jFexSRJEtaValues->at(1), jFexSRJPhiValues->at(1)));
            additionalSeedJetPositions.push_back(std::make_pair(jFexSRJEtaValues->at(2), jFexSRJPhiValues->at(2)));
            additionalSeedJetPositions.push_back(std::make_pair(jFexSRJEtaValues->at(3), jFexSRJPhiValues->at(3)));
            additionalSeedJetPositions.push_back(std::make_pair(jFexSRJEtaValues->at(4), jFexSRJPhiValues->at(4)));
            additionalSeedJetPositions.push_back(std::make_pair(jFexSRJEtaValues->at(5), jFexSRJPhiValues->at(5)));
        }

        if(seedObjectType == "gFEXSRJ"){
            seedJetPositions.push_back(std::make_pair(gFexSRJEtaValues->at(0), gFexSRJPhiValues->at(0)));
            seedJetPositions.push_back(std::make_pair(gFexSRJEtaValues->at(1), gFexSRJPhiValues->at(1)));
            additionalSeedJetPositions.push_back(std::make_pair(gFexSRJEtaValues->at(2), gFexSRJPhiValues->at(2)));
            additionalSeedJetPositions.push_back(std::make_pair(gFexSRJEtaValues->at(3), gFexSRJPhiValues->at(3)));
            additionalSeedJetPositions.push_back(std::make_pair(gFexSRJEtaValues->at(4), gFexSRJPhiValues->at(4)));
            additionalSeedJetPositions.push_back(std::make_pair(gFexSRJEtaValues->at(5), gFexSRJPhiValues->at(5)));
        }
        
        else if(seedObjectType == "ConeJet"){
            if(skBool){
                for(unsigned int iWTACone = 0; iWTACone < gepWTAConeCellsTowersSKJetsEtaValues->size(); iWTACone++){
                    if(iWTACone <= 1){
                        seedJetPositions.push_back(std::make_pair(gepWTAConeCellsTowersSKJetsEtaValues->at(iWTACone), gepWTAConeCellsTowersSKJetsPhiValues->at(iWTACone)));
                    }
                    else if(iWTACone <= 5){
                        additionalSeedJetPositions.push_back(std::make_pair(gepWTAConeCellsTowersSKJetsEtaValues->at(iWTACone), gepWTAConeCellsTowersSKJetsPhiValues->at(iWTACone)));
                    }
                    else break;
                }
            }
            else{
                for(unsigned int iWTACone = 0; iWTACone < gepWTAConeCellsTowersJetsEtaValues->size(); iWTACone++){
                    if(iWTACone <= 1){
                        seedJetPositions.push_back(std::make_pair(gepWTAConeCellsTowersJetsEtaValues->at(iWTACone), gepWTAConeCellsTowersJetsPhiValues->at(iWTACone)));
                    }
                    else if(iWTACone <= 5){
                        additionalSeedJetPositions.push_back(std::make_pair(gepWTAConeCellsTowersJetsEtaValues->at(iWTACone), gepWTAConeCellsTowersJetsPhiValues->at(iWTACone)));
                    }
                    else break;
                }
            }
        }

        newSeedPositions.push_back(std::make_pair(jetTaggerLeadingLRJEtaValues->at(0), jetTaggerLeadingLRJPhiValues->at(0)));
        newSeedPositions.push_back(std::make_pair(jetTaggerSubleadingLRJEtaValues->at(0), jetTaggerSubleadingLRJPhiValues->at(0)));
        lrjEt.push_back(jetTaggerLeadingLRJEtValues->at(0));
        lrjEt.push_back(jetTaggerSubleadingLRJEtValues->at(0));
        lrjPsi_R.push_back(jetTaggerLeadingLRJPsi_RValues->at(0));
        lrjPsi_R.push_back(jetTaggerSubleadingLRJPsi_RValues->at(0));
        std::array<std::vector<unsigned int >, 2 > subjetIndices;
        std::array<unsigned int, 2 > NSubjets = {0};
        if(subjetType == "WTACone"){
            if(skBool){
                for(unsigned int iConeJet = 0; iConeJet < gepWTAConeCellsTowersSKJetspTValues->size(); iConeJet++){
                    if(gepWTAConeCellsTowersSKJetspTValues->at(iConeJet) > subjetEtThreshold){
                        double deltaRPossibleSubjetsSeedLeading = sqrt(calcDeltaR2(gepWTAConeCellsTowersSKJetsEtaValues->at(iConeJet), gepWTAConeCellsTowersSKJetsPhiValues->at(iConeJet), jetTaggerLeadingLRJEtaValues->at(0), jetTaggerLeadingLRJPhiValues->at(0)));
                        double deltaRPossibleSubjetsSeedSubleading = sqrt(calcDeltaR2(gepWTAConeCellsTowersSKJetsEtaValues->at(iConeJet), gepWTAConeCellsTowersSKJetsPhiValues->at(iConeJet), jetTaggerSubleadingLRJEtaValues->at(0), jetTaggerSubleadingLRJPhiValues->at(0)));
                        if(deltaRPossibleSubjetsSeedLeading < jetRadius){
                            NSubjets[0]++;
                            subjetIndices[0].push_back(iConeJet);
                        }
                        if(deltaRPossibleSubjetsSeedSubleading < jetRadius){
                            NSubjets[1]++;
                            subjetIndices[1].push_back(iConeJet);
                        }
                        possibleSubjets.push_back(std::make_pair(gepWTAConeCellsTowersSKJetsEtaValues->at(iConeJet), gepWTAConeCellsTowersSKJetsPhiValues->at(iConeJet)));
                    }
                }
            }
            else{
                for(unsigned int iConeJet = 0; iConeJet < gepWTAConeCellsTowersJetspTValues->size(); iConeJet++){
                    if(gepWTAConeCellsTowersJetspTValues->at(iConeJet) > subjetEtThreshold){
                        double deltaRPossibleSubjetsSeedLeading = sqrt(calcDeltaR2(gepWTAConeCellsTowersJetsEtaValues->at(iConeJet), gepWTAConeCellsTowersJetsPhiValues->at(iConeJet), jetTaggerLeadingLRJEtaValues->at(0), jetTaggerLeadingLRJPhiValues->at(0)));
                        double deltaRPossibleSubjetsSeedSubleading = sqrt(calcDeltaR2(gepWTAConeCellsTowersJetsEtaValues->at(iConeJet), gepWTAConeCellsTowersJetsPhiValues->at(iConeJet), jetTaggerSubleadingLRJEtaValues->at(0), jetTaggerSubleadingLRJPhiValues->at(0)));
                        if(deltaRPossibleSubjetsSeedLeading < jetRadius){
                            NSubjets[0]++;
                            subjetIndices[0].push_back(iConeJet);
                        }
                        if(deltaRPossibleSubjetsSeedSubleading < jetRadius){
                            NSubjets[1]++;
                            subjetIndices[1].push_back(iConeJet);
                        }
                        possibleSubjets.push_back(std::make_pair(gepWTAConeCellsTowersJetsEtaValues->at(iConeJet), gepWTAConeCellsTowersJetsPhiValues->at(iConeJet)));
                    }
                }
            }
        }
        else if(subjetType == "jFEXSRJ"){
            for(unsigned int ijFEX = 0; ijFEX < jFexSRJEtValues->size(); ijFEX++){
                if(jFexSRJEtValues->at(ijFEX) > subjetEtThreshold){
                    double deltaRPossibleSubjetsSeedLeading = sqrt(calcDeltaR2(jFexSRJEtaValues->at(ijFEX), jFexSRJPhiValues->at(ijFEX), jetTaggerLeadingLRJEtaValues->at(0), jetTaggerLeadingLRJPhiValues->at(0)));
                    double deltaRPossibleSubjetsSeedSubleading = sqrt(calcDeltaR2(jFexSRJEtaValues->at(ijFEX), jFexSRJPhiValues->at(ijFEX), jetTaggerSubleadingLRJEtaValues->at(0), jetTaggerSubleadingLRJPhiValues->at(0)));
                    if(deltaRPossibleSubjetsSeedLeading < jetRadius){
                        NSubjets[0]++;
                        subjetIndices[0].push_back(ijFEX);
                    }
                    if(deltaRPossibleSubjetsSeedSubleading < jetRadius){
                        NSubjets[1]++;
                        subjetIndices[1].push_back(ijFEX);
                    }
                    possibleSubjets.push_back(std::make_pair(jFexSRJEtaValues->at(ijFEX), jFexSRJPhiValues->at(ijFEX)));
                }
            }
        }

        else if(subjetType == "gFEXSRJ"){
            for(unsigned int igFEX = 0; igFEX < gFexSRJEtValues->size(); igFEX++){
                if(gFexSRJEtValues->at(igFEX) > subjetEtThreshold){
                    double deltaRPossibleSubjetsSeedLeading = sqrt(calcDeltaR2(gFexSRJEtaValues->at(igFEX), gFexSRJPhiValues->at(igFEX), jetTaggerLeadingLRJEtaValues->at(0), jetTaggerLeadingLRJPhiValues->at(0)));
                    double deltaRPossibleSubjetsSeedSubleading = sqrt(calcDeltaR2(gFexSRJEtaValues->at(igFEX), gFexSRJPhiValues->at(igFEX), jetTaggerSubleadingLRJEtaValues->at(0), jetTaggerSubleadingLRJPhiValues->at(0)));
                    if(deltaRPossibleSubjetsSeedLeading < jetRadius){
                        NSubjets[0]++;
                        subjetIndices[0].push_back(igFEX);
                    }
                    if(deltaRPossibleSubjetsSeedSubleading < jetRadius){
                        NSubjets[1]++;
                        subjetIndices[1].push_back(igFEX);
                    }
                    possibleSubjets.push_back(std::make_pair(gFexSRJEtaValues->at(igFEX), gFexSRJPhiValues->at(igFEX)));
                }
            }
        }

        if(NSubjets[0] >= 2){
            if(subjetType == "WTACone"){
                if(skBool){
                    // Note: this assumes that subjets were pre-sorted within ntupler
                    double deltaRSubjets = sqrt(calcDeltaR2(gepWTAConeCellsTowersSKJetsEtaValues->at(subjetIndices[0][0]), gepWTAConeCellsTowersSKJetsPhiValues->at(subjetIndices[0][0]), gepWTAConeCellsTowersSKJetsEtaValues->at(subjetIndices[0][1]), gepWTAConeCellsTowersSKJetsPhiValues->at(subjetIndices[0][1])));
                    lrjDeltaRSubjets[0] = deltaRSubjets;
                    lrjEtRatioSubjets[0] = gepWTAConeCellsTowersSKJetspTValues->at(subjetIndices[0][0])/gepWTAConeCellsTowersSKJetspTValues->at(subjetIndices[0][1]);
                    lrjMassApprox[0] = deltaRSubjets * (gepWTAConeCellsTowersSKJetspTValues->at(subjetIndices[0][0]) + gepWTAConeCellsTowersSKJetspTValues->at(subjetIndices[0][1]));
                }
                else{
                    // Note: this assumes that subjets were pre-sorted within ntupler
                    double deltaRSubjets = sqrt(calcDeltaR2(gepWTAConeCellsTowersJetsEtaValues->at(subjetIndices[0][0]), gepWTAConeCellsTowersJetsPhiValues->at(subjetIndices[0][0]), gepWTAConeCellsTowersJetsEtaValues->at(subjetIndices[0][1]), gepWTAConeCellsTowersJetsPhiValues->at(subjetIndices[0][1])));
                    lrjDeltaRSubjets[0] = deltaRSubjets;
                    lrjEtRatioSubjets[0] = gepWTAConeCellsTowersJetspTValues->at(subjetIndices[0][0])/gepWTAConeCellsTowersJetspTValues->at(subjetIndices[0][1]);
                    lrjMassApprox[0] = deltaRSubjets * (gepWTAConeCellsTowersJetspTValues->at(subjetIndices[0][0]) + gepWTAConeCellsTowersJetspTValues->at(subjetIndices[0][1]));
                }
                
            }
            else if(subjetType == "jFEXSRJ"){
                double deltaRSubjets = sqrt(calcDeltaR2(jFexSRJEtaValues->at(subjetIndices[0][0]), jFexSRJPhiValues->at(subjetIndices[0][0]), jFexSRJEtaValues->at(subjetIndices[0][1]), jFexSRJPhiValues->at(subjetIndices[0][1])));
                lrjDeltaRSubjets[0] = deltaRSubjets;
                lrjEtRatioSubjets[0] = jFexSRJEtValues->at(subjetIndices[0][0])/jFexSRJEtValues->at(subjetIndices[0][1]);
                lrjMassApprox[0] = deltaRSubjets * (jFexSRJEtValues->at(subjetIndices[0][0]) + jFexSRJEtValues->at(subjetIndices[0][1]));
            }
            else if(subjetType == "gFEXSRJ"){
                double deltaRSubjets = sqrt(calcDeltaR2(gFexSRJEtaValues->at(subjetIndices[0][0]), gFexSRJPhiValues->at(subjetIndices[0][0]), gFexSRJEtaValues->at(subjetIndices[0][1]), gFexSRJPhiValues->at(subjetIndices[0][1])));
                lrjDeltaRSubjets[0] = deltaRSubjets;
                lrjEtRatioSubjets[0] = gFexSRJEtValues->at(subjetIndices[0][0])/gFexSRJEtValues->at(subjetIndices[0][1]);
                lrjMassApprox[0] = deltaRSubjets * (gFexSRJEtValues->at(subjetIndices[0][0]) + gFexSRJEtValues->at(subjetIndices[0][1]));
            }
        }
        if(NSubjets[1] >= 2){
            if(subjetType == "WTACone"){
                if(skBool){
                    double deltaRSubjets = sqrt(calcDeltaR2(gepWTAConeCellsTowersSKJetsEtaValues->at(subjetIndices[1][0]), gepWTAConeCellsTowersSKJetsPhiValues->at(subjetIndices[1][0]), gepWTAConeCellsTowersSKJetsEtaValues->at(subjetIndices[1][1]), gepWTAConeCellsTowersSKJetsPhiValues->at(subjetIndices[1][1])));
                    lrjDeltaRSubjets[1] = deltaRSubjets;
                    lrjEtRatioSubjets[1] = gepWTAConeCellsTowersSKJetspTValues->at(subjetIndices[1][0])/gepWTAConeCellsTowersSKJetspTValues->at(subjetIndices[1][1]);
                    lrjMassApprox[1] = deltaRSubjets * (gepWTAConeCellsTowersSKJetspTValues->at(subjetIndices[1][0]) + gepWTAConeCellsTowersSKJetspTValues->at(subjetIndices[1][1]));
                }
                else{
                    double deltaRSubjets = sqrt(calcDeltaR2(gepWTAConeCellsTowersJetsEtaValues->at(subjetIndices[1][0]), gepWTAConeCellsTowersJetsPhiValues->at(subjetIndices[1][0]), gepWTAConeCellsTowersJetsEtaValues->at(subjetIndices[1][1]), gepWTAConeCellsTowersJetsPhiValues->at(subjetIndices[1][1])));
                    lrjDeltaRSubjets[1] = deltaRSubjets;
                    lrjEtRatioSubjets[1] = gepWTAConeCellsTowersJetspTValues->at(subjetIndices[1][0])/gepWTAConeCellsTowersJetspTValues->at(subjetIndices[1][1]);
                    lrjMassApprox[1] = deltaRSubjets * (gepWTAConeCellsTowersJetspTValues->at(subjetIndices[1][0]) + gepWTAConeCellsTowersJetspTValues->at(subjetIndices[1][1]));
                }
            }
            else if(subjetType == "jFEXSRJ"){
                double deltaRSubjets = sqrt(calcDeltaR2(jFexSRJEtaValues->at(subjetIndices[1][0]), jFexSRJPhiValues->at(subjetIndices[1][0]), jFexSRJEtaValues->at(subjetIndices[1][1]), jFexSRJPhiValues->at(subjetIndices[1][1])));
                lrjDeltaRSubjets[1] = deltaRSubjets;
                lrjEtRatioSubjets[1] = jFexSRJEtValues->at(subjetIndices[1][0])/jFexSRJEtValues->at(subjetIndices[1][1]);
                lrjMassApprox[1] = deltaRSubjets * (jFexSRJEtValues->at(subjetIndices[1][0]) + jFexSRJEtValues->at(subjetIndices[1][1]));
            }
            else if(subjetType == "gFEXSRJ"){
                double deltaRSubjets = sqrt(calcDeltaR2(gFexSRJEtaValues->at(subjetIndices[1][0]), gFexSRJPhiValues->at(subjetIndices[1][0]), gFexSRJEtaValues->at(subjetIndices[1][1]), gFexSRJPhiValues->at(subjetIndices[1][1])));
                lrjDeltaRSubjets[1] = deltaRSubjets;
                lrjEtRatioSubjets[1] = gFexSRJEtValues->at(subjetIndices[1][0])/gFexSRJEtValues->at(subjetIndices[1][1]);
                lrjMassApprox[1] = deltaRSubjets * (gFexSRJEtValues->at(subjetIndices[1][0]) + gFexSRJEtValues->at(subjetIndices[1][1]));
            }
        }

        for(unsigned int iTruthJet = 0; iTruthJet < truthAntiKt4WZSRJEtaValues->size(); iTruthJet++){
            if(truthAntiKt4WZSRJEtValues->at(iTruthJet) > truthJetEtThreshold){
                truthJets.push_back(std::make_pair(truthAntiKt4WZSRJEtaValues->at(iTruthJet), truthAntiKt4WZSRJPhiValues->at(iTruthJet)));
            }
        }

        for(unsigned int iPUJet = 0; iPUJet < inTimeAntiKt4TruthSRJEtValues->size(); iPUJet++){
            if(inTimeAntiKt4TruthSRJEtValues->at(iPUJet) > truthJetEtThreshold){
                pileupJets.push_back(std::make_pair(inTimeAntiKt4TruthSRJEtaValues->at(iPUJet), inTimeAntiKt4TruthSRJPhiValues->at(iPUJet)));
            }
        }

        for(unsigned int iOOTPUJet = 0; iOOTPUJet < outOfTimeAntiKt4TruthSRJEtValues->size(); iOOTPUJet++){
            if(outOfTimeAntiKt4TruthSRJEtValues->at(iOOTPUJet) > 200.0){
                outOfTimePileupJets.push_back(std::make_pair(outOfTimeAntiKt4TruthSRJEtaValues->at(iOOTPUJet), outOfTimeAntiKt4TruthSRJPhiValues->at(iOOTPUJet)));
            }
        }

        if(signalString == "VBF_hh_4b_cvv0" || signalString == "VBF_hh_4b_cvv1" || signalString == "ggF_hh_4b"){
            // Map: Higgs index -> list of b-quark indices
            std::unordered_map<int, std::vector<int>> higgsToB;
            higgsToB.reserve(higgsIndexValues->size());
            for (size_t ib = 0; ib < higgsIndexValues->size(); ++ib) {
                higgsToB[ (*higgsIndexValues)[ib] ].push_back((int)ib);
            }
            for (size_t ih = 0; ih < indexOfHiggsValues->size(); ++ih) {

                const int hid = (*indexOfHiggsValues)[ih];
                const double higgspT = (*truthHiggspTValues)[ih];
                const double higgsEta = (*truthHiggsEtaValues)[ih];
                const double higgsPhi = (*truthHiggsPhiValues)[ih];


                higgses.push_back(std::make_pair(higgsEta, higgsPhi));
                higgsPts.push_back(higgspT);

                auto it = higgsToB.find(hid);
                if (it == higgsToB.end() || it->second.size() < 2) continue;

                const int ib1 = it->second[0];
                const int ib2 = it->second[1];
                const double et1 = (*truthbquarksEtValues)[ib1];
                const double eta1 = (*truthbquarksEtaValues)[ib1];
                const double phi1 = (*truthbquarksPhiValues)[ib1];
                const double et2 = (*truthbquarksEtValues)[ib2];
                const double eta2 = (*truthbquarksEtaValues)[ib2];
                const double phi2 = (*truthbquarksPhiValues)[ib2];
                if(et1 >= et2){
                    bQuarks[ih].push_back(std::make_pair(eta1, phi1));
                    bQuarks[ih].push_back(std::make_pair(eta2, phi2));
                    bEts[ih].push_back(et1);
                    bEts[ih].push_back(et2);
                }
                else{
                    bQuarks[ih].push_back(std::make_pair(eta2, phi2));
                    bQuarks[ih].push_back(std::make_pair(eta1, phi1));
                    bEts[ih].push_back(et2);
                    bEts[ih].push_back(et1);
                }
                

                const double dr2 = calcDeltaR2(eta1, phi1, eta2, phi2);
                const double dR  = std::sqrt(dr2);

                bQuarksDeltaR[ih] = dR;



                
            }
        }

        {
            c.cd();
            c.Clear();
            cellsTowersHighestEtSeedPositions->GetZaxis()->SetTitle("E_{T} [GeV]");
            cellsTowersHighestEtSeedPositions->Draw("COLZ");
            //std::cout << "seedJetPositions.size(): " << seedJetPositions.size() << "\n";
            //std::cout << "additionalSeedJetPositions.size(): " << additionalSeedJetPositions.size() << "\n";
            for(unsigned int iSeedOriginal = 0; iSeedOriginal < seedJetPositions.size(); iSeedOriginal++){
                TEllipse *circle0 = new TEllipse(seedJetPositions[iSeedOriginal].first, seedJetPositions[iSeedOriginal].second, jetRadius, jetRadius); // R in both x and y
                circle0->SetLineColor(kRed);
                circle0->SetLineWidth(2);
                circle0->SetFillStyle(0); // no fill
                circle0->SetLineStyle(2);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                circle0->Draw("same");    // overlay on the existing plot
            }
            for(unsigned int iSeedAdditional = 0; iSeedAdditional < additionalSeedJetPositions.size(); iSeedAdditional++){
                TEllipse *circle2 = new TEllipse(additionalSeedJetPositions[iSeedAdditional].first, additionalSeedJetPositions[iSeedAdditional].second, jetRadius, jetRadius); // R in both x and y
                circle2->SetLineColor(kAzure+2);
                circle2->SetLineWidth(2);
                circle2->SetFillStyle(0); // no fill
                circle2->SetLineStyle(2);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                circle2->Draw("same");    // overlay on the existing plot
            }
            
            // NOTE there will always be 2 jet tagger jets.
            TEllipse *circle6 = new TEllipse(newSeedPositions[0].first, newSeedPositions[0].second, jetRadius, jetRadius); // R in both x and y
            circle6->SetLineColor(kGreen+2);
            circle6->SetLineWidth(2);
            circle6->SetFillStyle(0); // no fill
            circle6->SetLineStyle(1);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
            circle6->Draw("same");    // overlay on the existing plot

            TEllipse *circle7 = new TEllipse(newSeedPositions[1].first, newSeedPositions[1].second, jetRadius, jetRadius); // R in both x and y
            circle7->SetLineColor(kGreen+2);
            circle7->SetLineWidth(2);
            circle7->SetFillStyle(0); // no fill
            circle7->SetLineStyle(1);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
            circle7->Draw("same");    // overlay on the existing plot

            // Build the label text (1st LRJ)
            
            TString inputFileLabel1 = Form("Lead. E_{T}=%.1f GeV", lrjEt[0]);

            // Option 1: Use TLatex directly
            TLatex lat1;
            lat1.SetTextSize(0.025);     // Adjust size
            lat1.SetTextColor(kBlack); // Match circle color
            lat1.DrawLatex(newSeedPositions[0].first, newSeedPositions[0].second + 0.6, inputFileLabel1); // Slightly above the circle

            // Build the label text (2nd LRJ)
            
            TString inputFileLabel2 = Form("Subl. E_{T}=%.1f GeV", lrjEt[1]);

            // Option 1: Use TLatex directly
            TLatex Lat2;
            Lat2.SetTextSize(0.025);     // Adjust size
            Lat2.SetTextColor(kBlack); // Match circle color
            Lat2.DrawLatex(newSeedPositions[1].first, newSeedPositions[1].second + 0.6, inputFileLabel2); // Slightly above the circle

            c.Print(pdf);   // add a page

            c.Clear();

            TLatex lat;
            lat.SetTextSize(0.035);
            lat.SetTextAlign(13);  // left, top

            double y = 0.85;

            double deltaRSubjetsLead = -1; 
            double lrjEtRatioLead = -1;
            double massApproxLead = -1;
            if(NSubjets[0] >= 2){
                deltaRSubjetsLead = lrjDeltaRSubjets[0];
                lrjEtRatioLead = lrjEtRatioSubjets[0];
                massApproxLead = lrjMassApprox[0];
            }
            lat.DrawLatexNDC(
                0.10, 0.95,
                Form("Event = %.0d", iEvt)
            );
            lat.DrawLatexNDC(0.10, y,
                Form("Leading LRJ: E_{T}=%.1f GeV, #psi_{R}=%.2f, #DeltaR=%.2f, E_{T} ratio=%.2f",
                    lrjEt[0], lrjPsi_R[0], deltaRSubjetsLead, lrjEtRatioLead));

            y -= 0.08;
            
            double deltaRSubjetsSublead = -1; 
            double lrjEtRatioSublead = -1;
            double massApproxSublead = -1;
            if(NSubjets[1] >= 2){
                deltaRSubjetsSublead = lrjDeltaRSubjets[1];
                lrjEtRatioSublead = lrjEtRatioSubjets[1];
                massApproxSublead = lrjMassApprox[1]; 
            }

            lat.DrawLatexNDC(0.10, y,
                Form("Subleading LRJ: E_{T}=%.1f GeV, #psi_{R}=%.2f, #DeltaR=%.2f, E_{T} ratio=%.2f",
                    lrjEt[1], lrjPsi_R[1], deltaRSubjetsSublead, lrjEtRatioSublead));

            c.Print(pdf);   // ===== PAGE 2 =====
        } // within scope to allow for renaming variables. 

        { // generate event displays with cone subjets, jettagger LRJ, truth & pileup jets + b, higgs information
            cEventDisplay.cd();
            cEventDisplay.Clear();
            cellsTowersHighestEtTruthAndPileupJets->GetZaxis()->SetTitle("E_{T} [GeV]");
            cellsTowersHighestEtTruthAndPileupJets->Draw("COLZ");

            TEllipse *circle0 = new TEllipse(newSeedPositions[0].first, newSeedPositions[0].second, jetRadius, jetRadius); // R in both x and y
            circle0->SetLineColor(kBlue);
            circle0->SetLineWidth(2);
            circle0->SetFillStyle(0); // no fill
            circle0->SetLineStyle(1);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
            circle0->Draw("same");    // overlay on the existing plot

            TEllipse *circle1 = new TEllipse(newSeedPositions[1].first, newSeedPositions[1].second, jetRadius, jetRadius); // R in both x and y
            circle1->SetLineColor(kBlue);
            circle1->SetLineWidth(2);
            circle1->SetFillStyle(0); // no fill
            circle1->SetLineStyle(1);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
            circle1->Draw("same");    // overlay on the existing plot

            for(unsigned int iSubjet = 0; iSubjet < possibleSubjets.size(); iSubjet++){
                //if(subjetType == "WTACone"){
                    TEllipse *circle = new TEllipse(possibleSubjets[iSubjet].first, possibleSubjets[iSubjet].second, 0.4, 0.4); // R in both x and y
                    circle->SetLineColor(kBlack);
                    circle->SetLineWidth(2);
                    circle->SetFillStyle(0); // no fill
                    circle->SetLineStyle(2);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                    circle->Draw("same");    // overlay on the existing plot
                //}
                /*else if(subjetType == "jFEXSRJ"){
                    if()
                    TBox *square0 = new TBox(possibleSubjets[iSubjet].first - R, possibleSubjets[iSubjet].second - R, x + R, y + R);
                    square0->SetLineColor(kRed);
                    square0->SetLineWidth(2);
                    square0->SetFillStyle(0);   // no fill
                    square0->SetLineStyle(2);   // dashed
                    square0->Draw("same");
                }*/ // FIXME use actual jFEX specifications to generate square (or whatever it should be)
                
            }

            for(unsigned int iPU = 0; iPU < pileupJets.size(); iPU++){
                TEllipse *circle = new TEllipse(pileupJets[iPU].first, pileupJets[iPU].second, 0.4, 0.4); // R in both x and y
                circle->SetLineColor(kGreen);
                circle->SetLineWidth(2);
                circle->SetFillStyle(0); // no fill
                circle->SetLineStyle(3);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                circle->Draw("same");    // overlay on the existing plot
            }

            for(unsigned int iOOTPU = 0; iOOTPU < outOfTimePileupJets.size(); iOOTPU++){
                TEllipse *circle = new TEllipse(outOfTimePileupJets[iOOTPU].first, outOfTimePileupJets[iOOTPU].second, 0.4, 0.4);
                circle->SetLineColor(kViolet);
                circle->SetLineWidth(2);
                circle->SetFillStyle(0);
                circle->SetLineStyle(3);
                circle->Draw("same");
            }

            for(unsigned int iTruthJet = 0; iTruthJet < truthJets.size(); iTruthJet++){
                TEllipse *circle = new TEllipse(truthJets[iTruthJet].first, truthJets[iTruthJet].second, 0.4, 0.4); // R in both x and y
                circle->SetLineColor(kCyan);
                circle->SetLineWidth(2);
                circle->SetFillStyle(0); // no fill
                circle->SetLineStyle(1);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                circle->Draw("same");    // overlay on the existing plot
            }
            if(signalString == "VBF_hh_4b_cvv0" || signalString == "VBF_hh_4b_cvv1" || signalString == "ggF_hh_4b"){
                TEllipse *circle_h1 = new TEllipse(higgses[0].first, higgses[0].second, 0.1, 0.1); // R in both x and y
                circle_h1->SetLineColor(kOrange);
                circle_h1->SetLineWidth(2);
                circle_h1->SetFillStyle(0); // no fill
                circle_h1->SetLineStyle(1);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                circle_h1->Draw("same");    // overlay on the existing plot

                TEllipse *circle_b11 = new TEllipse(bQuarks[0][0].first, bQuarks[0][0].second, 0.1, 0.1); // R in both x and y
                circle_b11->SetLineColor(kRed);
                circle_b11->SetLineWidth(2);
                circle_b11->SetFillStyle(0); // no fill
                circle_b11->SetLineStyle(1);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                circle_b11->Draw("same");    // overlay on the existing plot

                TEllipse *circle_b12 = new TEllipse(bQuarks[0][1].first, bQuarks[0][1].second, 0.1, 0.1); // R in both x and y
                circle_b12->SetLineColor(kMagenta);
                circle_b12->SetLineWidth(2);
                circle_b12->SetFillStyle(0); // no fill
                circle_b12->SetLineStyle(1);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                circle_b12->Draw("same");    // overlay on the existing plot

                TEllipse *circle_h2 = new TEllipse(higgses[1].first, higgses[1].second, 0.1, 0.1); // R in both x and y
                circle_h2->SetLineColor(kOrange);
                circle_h2->SetLineWidth(2);
                circle_h2->SetFillStyle(0); // no fill
                circle_h2->SetLineStyle(2);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                circle_h2->Draw("same");    // overlay on the existing plot

                TEllipse *circle_b21 = new TEllipse(bQuarks[1][0].first, bQuarks[1][0].second, 0.1, 0.1); // R in both x and y
                circle_b21->SetLineColor(kRed);
                circle_b21->SetLineWidth(2);
                circle_b21->SetFillStyle(0); // no fill
                circle_b21->SetLineStyle(2);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                circle_b21->Draw("same");    // overlay on the existing plot

                TEllipse *circle_b22 = new TEllipse(bQuarks[1][1].first, bQuarks[1][1].second, 0.1, 0.1); // R in both x and y
                circle_b22->SetLineColor(kMagenta);
                circle_b22->SetLineWidth(2);
                circle_b22->SetFillStyle(0); // no fill
                circle_b22->SetLineStyle(2);  // dashed (1=solid, 2=dashed, 3=dotted, etc.)
                circle_b22->Draw("same");    // overlay on the existing plot

                // Build the label text (1st LRJ)
                TString higgsLabel1 = Form("#Delta R (b's) [h1]=%.2f", bQuarksDeltaR[0]);

                // Option 1: Use TLatex directly
                TLatex lat_Higgs1;
                lat_Higgs1.SetTextSize(0.025);     // Adjust size
                lat_Higgs1.SetTextColor(kBlack); // Match circle color
                lat_Higgs1.DrawLatex(higgses[0].first, higgses[0].second, higgsLabel1); // Slightly above the circle

                // Build the label text (1st LRJ)
                TString higgsLabel2 = Form("#Delta R (b's) [h2]=%.2f", bQuarksDeltaR[1]);

                // Option 1: Use TLatex directly
                TLatex lat_Higgs2;
                lat_Higgs2.SetTextSize(0.025);     // Adjust size
                lat_Higgs2.SetTextColor(kBlack); // Match circle color
                lat_Higgs2.DrawLatex(higgses[1].first, higgses[1].second, higgsLabel2); // Slightly above the circle

                
            }

            // Build the label text (1st LRJ)
            TString label1 = Form("Lead. E_{T}=%.1f GeV", lrjEt[0]);

            // Option 1: Use TLatex directly
            TLatex lat1;
            lat1.SetTextSize(0.025);     // Adjust size
            lat1.SetTextColor(kBlack); // Match circle color
            lat1.DrawLatex(newSeedPositions[0].first, newSeedPositions[0].second + 0.6, label1); // Slightly above the circle

            // Build the label text (2nd LRJ)
            TString label2 = Form("Subl. E_{T}=%.1f GeV", lrjEt[1]);

            // Option 1: Use TLatex directly
            TLatex Lat2;
            Lat2.SetTextSize(0.025);     // Adjust size
            Lat2.SetTextColor(kBlack); // Match circle color
            Lat2.DrawLatex(newSeedPositions[1].first, newSeedPositions[1].second + 0.6, label2); // Slightly above the circle
            
            

            cEventDisplay.Print(pdf_ED);

            cEventDisplay.Clear();

            TLatex lat;
            lat.SetTextSize(0.035);
            lat.SetTextAlign(13);  // left, top

            double y = 0.85;

            double deltaRSubjetsLead = -1; 
            double lrjEtRatioLead = -1;
            double massApproxLead = -1;
            if(NSubjets[0] >= 2){
                deltaRSubjetsLead = lrjDeltaRSubjets[0];
                lrjEtRatioLead = lrjEtRatioSubjets[0];
                massApproxLead = lrjMassApprox[0];
            }

            lat.DrawLatexNDC(0.05, y,
                Form("Lead.: E_{T}=%.1f GeV, #psi_{R}=%.2f, #DeltaR=%.2f, r_{E_{T}}=%.2f, #DeltaR #times r_{E_{T}}=%.2f",
                    lrjEt[0], lrjPsi_R[0], deltaRSubjetsLead, lrjEtRatioLead, deltaRSubjetsLead*lrjEtRatioLead));

            y -= 0.08;
            
            double deltaRSubjetsSublead = -1; 
            double lrjEtRatioSublead = -1;
            double massApproxSublead = -1; 
            if(NSubjets[1] >= 2){
                deltaRSubjetsSublead = lrjDeltaRSubjets[1];
                lrjEtRatioSublead = lrjEtRatioSubjets[1];
                massApproxSublead = lrjMassApprox[1];
            }

            lat.DrawLatexNDC(0.05, y,
                Form("Subl. LRJ: E_{T}=%.1f GeV, #psi_{R}=%.2f, #DeltaR=%.2f, r_{E_{T}}=%.2f, #DeltaR #times r_{E_{T}}=%.2f",
                    lrjEt[1], lrjPsi_R[1], deltaRSubjetsSublead, lrjEtRatioSublead, deltaRSubjetsSublead*lrjEtRatioSublead));
            y -= 0.10;        
            lat.DrawLatexNDC(0.05, y,
                    Form(" Leading Mass Approx=%.1f [GeV], Subleading Mass Approx =%.1f [GeV]",
                    massApproxLead, massApproxSublead));        
            if(signalString == "VBF_hh_4b_cvv0" || signalString == "VBF_hh_4b_cvv1" || signalString == "ggF_hh_4b"){
                //std::cout << "why not triggered" << "\n";
                //std::cout << "y1: " << y << "\n";
                y -= 0.10;
                lat.DrawLatexNDC(0.05, y,
                    Form("h_{1} p_{T}=%.1f, Lead. b E_{T}=%.1f, Subl. b E_{T}=%.1f [GeV], #DeltaR_{b's}=%.2f",
                        higgsPts[0], bEts[0][0], bEts[0][1], bQuarksDeltaR[0]));
                //std::cout << "y2: " << y << "\n";
                y -= 0.10;        
                lat.DrawLatexNDC(0.05, y,
                    Form("h_{2} p_{T}=%.1f, Lead. b E_{T}=%.1f, Subl. b E_{T}=%.1f [GeV], #DeltaR_{b's}=%.2f",
                        higgsPts[1], bEts[1][0], bEts[1][1], bQuarksDeltaR[1]));

                
                //std::cout << "y3: " << y << "\n";
            }
            if(!signalBool){
                y -= 0.10;
                lat.DrawLatexNDC(0.05, y,
                    Form("Event Weight (Rate Contribution):%.1f",
                        eventWeightsValues->at(0)));
                y -= 0.10;
                lat.DrawLatexNDC(0.05, y,
                    Form("JZ Slice: %d, Event Number: %d",
                        sampleJZSliceValues, iEvt));

            }

            y -= 0.10;
            lat.DrawLatexNDC(0.05, y,
                Form("Lead. Offline LRJ: E_{T}=%.1f GeV, #eta=%.2f, #phi=%.2f",
                    recoAntiKt10LRJLeadingEtValues->at(0), recoAntiKt10LRJLeadingEtaValues->at(0), recoAntiKt10LRJLeadingPhiValues->at(0)));
            y -= 0.10;
            lat.DrawLatexNDC(0.05, y,
                Form("Lead. Truth Jet: E_{T}=%.1f GeV, #eta=%.2f, #phi=%.2f",
                    truthAntiKt4WZSRJLeadingEtValues->at(0), truthAntiKt4WZSRJLeadingEtaValues->at(0), truthAntiKt4WZSRJLeadingPhiValues->at(0)));
            y -= 0.10;
            lat.DrawLatexNDC(0.05, y,
                Form("Lead. In-Time PU Truth Jet: E_{T}=%.1f GeV, #eta=%.2f, #phi=%.2f",
                    inTimeAntiKt4TruthSRJLeadingEtValues->at(0), inTimeAntiKt4TruthSRJLeadingEtaValues->at(0), inTimeAntiKt4TruthSRJLeadingPhiValues->at(0)));
            y -= 0.10;
            lat.DrawLatexNDC(0.05, y,
                Form("Lead. Out-of-Time PU Truth Jet: E_{T}=%.1f GeV, #eta=%.2f, #phi=%.2f",
                    outOfTimeAntiKt4TruthSRJLeadingEtValues->at(0), outOfTimeAntiKt4TruthSRJLeadingEtaValues->at(0), outOfTimeAntiKt4TruthSRJLeadingPhiValues->at(0)));

            cEventDisplay.Print(pdf_ED);   // ===== PAGE 2 =====
        }
    }
    
    // close multipage pdf
    c.Print(pdf + ")");
    cEventDisplay.Print(pdf_ED + ")");

}


void makeJetTaggerEventDisplays(){
    gErrorIgnoreLevel = kWarning;
    std::string inputFile = "";
    //std::string inputFile = "";
    unsigned int desiredJZSlice = -1;
    //unsigned int desiredJZSlice = 3;
    bool signalBool = true;
    //bool signalBool = false;
    std::string subjetType = "WTACone";

    /*std::string signalString;
    if(signalBool){
        if(vbfBool){
            signalString = "VBF_hh_4b";
        }
        else{
            signalString = "ggF_hh_4b";
        }
    } 
    else signalString = "jj_" + std::to_string(desiredJZSlice);*/
    // Call for signal, both subjet types
    // VBF hh->4b (Cone subjets)
    /*makeEventDisplays("/data/larsonma/LargeRadiusJets/outputNTuplesDev_VBFSM/mc21_14TeV_hh_bbbb_vbf_novhh_e8557_s4422_r16130_rMerge_2_IOs_128_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_NoSK.root", 
        -1, true, "VBF_hh_4b", "WTACone", false);*/
    
    // VBF hh->4b (cvv = 0) PU Suppressed (Cone subjets)
    //callMakeEventDisplays("<jet_tagger_file>", "<her_input_file>",
    //    -1, true, "VBF_hh_4b_cvv0", "WTACone", true);

    // VBF hh->4b (cvv = 1, SM) PU Suppressed (Cone subjets)
    //callMakeEventDisplays("<jet_tagger_file>", "<her_input_file>",
    //    -1, true, "VBF_hh_4b_cvv1", "WTACone", true);

    // ggF hh->4b (Cone subjets)
    //makeEventDisplays("/data/larsonma/LargeRadiusJets/outputNTuplesDev_VBFSM/mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_rMerge_2_IOs_128_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_NoSK.root", 
    //    -1, true, false, "ggF_hh_4b", false);
    
    // ggF hh->4b PU Suppressed (Cone subjets) — old path, kept for reference
    //callMakeEventDisplays("<jet_tagger_file>", "<her_input_file>",
    //    -1, true, "ggF_hh_4b", "WTACone", true);

    // ggF hh->4b, rMerge=2, R=1.1 (WTACone subjets, SK)
    /*callMakeEventDisplays(
        "/data/larsonma/LargeRadiusJets/outputNTuplesDev_CondorSubmission_NewSamples/mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_rMerge_2_IOs_128_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_SK_subjetEt25GeV_ewm0_mep1_mec20GeV_v3_file0.root",
        "/data/larsonma/GEPHadronicEventReconstruction/ntuples/ggF_HHbbbb_v3/mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_resim_DAOD_NTUPLE_GEP_000001.root",
        -1, true, "ggF_hh_4b", "WTACone", true);*/

    // ggF hh->4b, rMerge=0.001, R=1.1 (WTACone subjets, SK)
    /*callMakeEventDisplays(
        "/data/larsonma/LargeRadiusJets/outputNTuplesDev_CondorSubmission_NewSamples/mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_rMerge_0.001_IOs_128_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_SK_subjetEt25GeV_ewm0_mep1_mec20GeV_v3_file0.root",
        "/data/larsonma/GEPHadronicEventReconstruction/ntuples/ggF_HHbbbb_v3/mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_resim_DAOD_NTUPLE_GEP_000001.root",
        -1, true, "ggF_hh_4b", "WTACone", true);*/

    callMakeEventDisplays(
        "/data/larsonma/LargeRadiusJets/outputNTuplesDev_CondorSubmission_NewSamples/mc21_14TeV_jj_JZ_e8557_s4422_r16130_rMerge_0.001_IOs_128_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_SK_subjetEt25GeV_ewm0_mep1_mec20GeV_v3.root",
        "/data/larsonma/GEPHadronicEventReconstruction/ntuples/mc21_14TeV_jj_JZ_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root",
        1, false, "jj_1", "WTACone", true);

    //callMakeEventDisplays("<jet_tagger_file>", "<her_input_file>",
    //    -1, true, "ggF_hh_4b", "gFEXSRJ", true);

    //callMakeEventDisplays("<jet_tagger_file>", "<her_input_file>",
    //    -1, true, "ggF_hh_4b", "jFEXSRJ", true);

    // ttbar fully hadronic decays (Cone subjets)
    //makeEventDisplays("/data/larsonma/LargeRadiusJets/outputNTuplesDev_VBFSM/mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_rMerge_2_IOs_128_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_NoSK.root", 
    //    -1, true, false, "ttbar_had", false);
    
    // ttbar fully hadronic decays PU Suppressed (Cone subjets)
    //callMakeEventDisplays("<jet_tagger_file>", "<her_input_file>",
    //    -1, true, "ttbar_had", "WTACone", true);

    // Zprime->ttbar PU Suppressed (Cone subjets)
    //callMakeEventDisplays("<jet_tagger_file>", "<her_input_file>",
    //    -1, true, "Zprime_ttbar", "WTACone", true);
    
    //makeEventDisplays("/data/larsonma/LargeRadiusJets/outputNTuplesDev_VBFSM/mc21_14TeV_hh_bbbb_vbf_novhh_e8557_s4422_r16130_rMerge_2_IOs_128_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_OR.root", 
    //    -1, true, "jFEXSRJ");

    // JZ2 background (Cone subjets)
    //makeEventDisplays("/data/larsonma/LargeRadiusJets/outputNTuplesDev_VBFSM/mc21_14TeV_jj_JZ_e8557_s4422_r16130_rMerge_2_IOs_128_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_NoSK.root", 
    //    2, false, "jj_2", "WTACone", false);

    // JZ2 background PU suppressed (Cone subjets)
    //callMakeEventDisplays("<jet_tagger_file>", "/data/larsonma/GEPHadronicEventReconstruction/ntuples/QCD_Dijet_JZ2_v3/<her_file>",
    //    2, false, "jj_2", "WTACone", true);

    //callMakeEventDisplays("<jet_tagger_file>", "/data/larsonma/GEPHadronicEventReconstruction/ntuples/QCD_Dijet_JZ2_v3/<her_file>",
    //    2, false, "jj_2", "gFEXSRJ", true);

    //callMakeEventDisplays("<jet_tagger_file>", "/data/larsonma/GEPHadronicEventReconstruction/ntuples/QCD_Dijet_JZ2_v3/<her_file>",
    //    2, false, "jj_2", "jFEXSRJ", true);

    // JZ3 background (Cone subjets)
    //makeEventDisplays("/data/larsonma/LargeRadiusJets/outputNTuplesDev_VBFSM/mc21_14TeV_jj_JZ_e8557_s4422_r16130_rMerge_2_IOs_128_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_NoSK.root", 
    //    3, false, "jj_3", "WTACone", false);

    // JZ3 background PU suppressed (Cone subjets)
    //callMakeEventDisplays("<jet_tagger_file>", "/data/larsonma/GEPHadronicEventReconstruction/ntuples/QCD_Dijet_JZ3_v3/<her_file>",
    //    3, false, "jj_3", "WTACone", true);

    // Call for background (jz2, jz3)
    /*makeEventDisplays("/data/larsonma/LargeRadiusJets/outputNTuplesDev_VBFSM/mc21_14TeV_jj_JZ_e8557_s4422_r16130_rMerge_2_IOs_32_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_OR.root", 
        2, false, "WTACone");

    makeEventDisplays("/data/larsonma/LargeRadiusJets/outputNTuplesDev_VBFSM/mc21_14TeV_jj_JZ_e8557_s4422_r16130_rMerge_2_IOs_32_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_OR.root", 
        2, false, "jFEXSRJ");

    makeEventDisplays("/data/larsonma/LargeRadiusJets/outputNTuplesDev_VBFSM/mc21_14TeV_jj_JZ_e8557_s4422_r16130_rMerge_2_IOs_32_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_OR.root", 
        3, false, "WTACone");

    makeEventDisplays("/data/larsonma/LargeRadiusJets/outputNTuplesDev_VBFSM/mc21_14TeV_jj_JZ_e8557_s4422_r16130_rMerge_2_IOs_32_Seeds_2_R2_1.21_IO_gepCellsTowers_Seed_gepWTAConeCellsTowersJets_OR.root", 
        3, false, "jFEXSRJ");*/
    
    gSystem->Exit(0);
}