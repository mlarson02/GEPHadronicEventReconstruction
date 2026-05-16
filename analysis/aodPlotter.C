// Save as: process_xaod.C
// To execute: root -l -b -q process_xaod.C

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <filesystem> // C++17
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH2F.h"
#include "analysisHelperFunctions.h"

using namespace std;

// Settings
const bool signalBool = true;
const bool lxplusBool = true;
const bool higgsPtCut = true;
const bool jFexPlotsBool = false;
const bool daodBool = true;

const float higgs_min_pt = 100.0; // GeV
const float higgs_min_average_pt = 200.0; // GeV

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
void find_non_higgs_daughters(const xAOD::TruthParticle* particle, std::vector<std::vector<float>>& result) {
    if (!particle) return;

    // Check if the particle is not a Higgs
    if (particle->pdgId() != 25) { // Not Higgs
        float E = particle->e();
        float pt = particle->pt();
        float eta = particle->eta();
        float phi = particle->phi();
        float p = (std::abs(eta) < 100) ? pt * std::cosh(eta) : INFINITY;
        float Et = (p > 0) ? (E * (pt / p)) : 0.0;

        // Push the result as a vector of floats (instead of a tuple)
        result.push_back({Et / 1000.0f, eta, phi}); // Store in GeV
        return;
    }

    // If the particle is a Higgs, recursively find its non-Higgs daughters
    for (unsigned int i = 0; i < particle->nChildren(); ++i) {
        find_non_higgs_daughters(particle->child(i), result);
    }
}

// Main function
void process_xaod() {
    SetPlotStyle();
    // Setup file paths
    string fileDir;
    if (daodBool){
        fileDir = "/data/larsonma/LargeRadiusJets/datasets/Signal_HHbbbb_DAOD_JETM42/mc21_14TeV.537536.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv0cv1.deriv.DAOD_JETM42.e8557_s4422_r16130_p6658";
    }
    else{
        if (signalBool) {
            fileDir = lxplusBool ?
                "/data/larsonma/LargeRadiusJets/datasets/Signal_HHbbbb/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130" :
                "/home/larsonma/GEPHadronicEventReconstruction/datasets/Signal_HHbbbb/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130";
        } else {
            fileDir = lxplusBool ?
                "/data/larsonma/LargeRadiusJets/datasets/Background_jj_JZ3/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130" :
                "/home/larsonma/GEPHadronicEventReconstruction/datasets/Background_jj_JZ3/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130";
        }
        }
    

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

    // Declare vectors if not already
    std::vector<float> hlt_jet_pt_values;
    std::vector<float> gfex_larger_jet_pt_values;
    std::vector<float> jet_pt_values;
    std::vector<float> jet_log_pt_values;
    std::vector<float> deltaR_values;
    std::vector<float> higgs_pt_values_for_dr;
    std::vector<float> higgs_pt_values;
    std::vector<float> higgs_pt_values_by_event;
    std::vector<float> higgs_pt_values_after_cut;
    std::vector<float> b_Et_values;
    std::vector<float> b_Et_values_after_higgs_cut;

    std::vector<float> highestEt_b_closest_gfex_deltar_list;
    std::vector<float> highestEt_b_closest_truthjet_deltar_list;
    std::vector<float> second_highestEt_b_closest_gfex_deltar_list;
    std::vector<float> second_highestEt_b_closest_truthjet_deltar_list;
    std::vector<float> third_highestEt_b_closest_gfex_deltar_list;
    std::vector<float> third_highestEt_b_closest_truthjet_deltar_list;
    std::vector<float> fourth_highestEt_b_closest_gfex_deltar_list;
    std::vector<float> fourth_highestEt_b_closest_truthjet_deltar_list;

    std::vector<float> b_gfex_deltar_list;

    std::vector<float> b11_list, b21_list, b12_list, b22_list;

    // Main processing loop
    int fileIt = 0;
    for (const auto& fileName : fileNames) {
        fileIt++; 
        if (fileIt > 9) break; 
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

        cout << "  Number of events: " << event.getEntries() << endl;

        

        for (Long64_t entry = 0; entry < event.getEntries(); ++entry) {
            event.getEntry(entry);

            // -- ADD collections you want here ---
            const xAOD::TruthParticleContainer* truthParticles = nullptr;
            if (!event.retrieve(truthParticles, "TruthParticles").isSuccess()) {
                cerr << "Cannot access TruthParticles" << endl;
                continue;
            }

            const xAOD::JetContainer* HLTJets = nullptr;
            if (!event.retrieve(HLTJets, "HLT_AntiKt4EMTopoJets_subjesIS").isSuccess()) {
                cerr << "Failed to retrieve HLT jets" << endl;
                continue;
            }

            const xAOD::JetContainer* RecoJets = nullptr;
            if (!event.retrieve(RecoJets, "AntiKt10UFOCSSKJets").isSuccess()) {
                cerr << "Failed to retrieve HLT jets" << endl;
                continue;
            }


            const xAOD::gFexJetRoIContainer* gFexSRJets = nullptr;
            if (!event.retrieve(gFexSRJets, "L1_gFexSRJetRoI").isSuccess()) {
                std::cerr << "Failed to retrieve gFex SR jets" << std::endl;
                continue;
            }

            const xAOD::gFexJetRoIContainer* gFexLRJets = nullptr;
            if (!event.retrieve(gFexLRJets, "L1_gFexLRJetRoI").isSuccess()) {
                cerr << "Failed to retrieve gFex LR jets" << endl;
                continue;
            }

            const xAOD::JetContainer* TruthJets = nullptr;
            if (!event.retrieve(TruthJets, "InTimeAntiKt4TruthJets").isSuccess()) {
                cerr << "Failed to retrieve Truth jets" << endl;
                continue;
            }

            // Initialize per-event counters
            int higgs_counter = 0;
            bool higgs_passes_cut = false;

            float higgs_1_pt = 0, higgs_2_pt = 0;
            float higgs_1_eta = 0, higgs_1_phi = 0;
            float higgs_2_eta = 0, higgs_2_phi = 0;

            // --- Loop over HLT_AntiKt4EMTopoJets_subjesIS ---
            for (const auto& el : *HLTJets) {
                hlt_jet_pt_values.push_back(el->pt() / 1000.0);
                float hlt_jet_total_p = el->pt() * cosh(el->eta());
                float hlt_jet_ptoverp = el->pt() / hlt_jet_total_p;
                float hlt_jet_Et = (el->e() * hlt_jet_ptoverp) / 1000.0;
            }

            // --- Loop over L1_gFexLRJetRoI ---
            //for (const auto& el : *gFexLRJets) {
            //    gfex_larger_jet_pt_values.push_back(el->et() / 1000.0);
            //    float gfex_larger_jet_Et = el->et() / 1000.0;
            //}

            // --- Loop over InTimeAntiKt4TruthJets ---
            for (const auto& el : *TruthJets) {
                jet_pt_values.push_back(el->pt() / 1000.0);
                jet_log_pt_values.push_back(log10(el->pt() / 1000.0));
                float truth_jet_total_p = el->pt() * cosh(el->eta());
                float truth_jet_ptoverp = el->pt() / truth_jet_total_p;
                float truth_jet_Et = (el->e() * truth_jet_ptoverp) / 1000.0;
            }

            // --- Loop over TruthParticles (for Higgs and B's) ---
            std::vector<std::vector<float > > allb_list;
            for (const auto& el : *truthParticles) {
                if (el->pdgId() == 25 && el->status() == 22) {
                    higgs_counter++;
                    float total_p = el->pt() * cosh(el->eta());
                    float ptoverp = el->pt() / total_p;

                    if (higgs_counter == 1) {
                        higgs_1_pt = el->pt()/1000.0;
                        higgs_1_eta = el->eta();
                        higgs_1_phi = el->phi();
                        vector<vector<float>> b1_list;
                        find_non_higgs_daughters(el, b1_list);
                        b11_list = b1_list[0];
                        b21_list = b1_list[1];
                        float dR_1 = delta_R(b11_list[1], b11_list[2], b21_list[1], b21_list[2]);
                        deltaR_values.push_back(dR_1);
                        higgs_pt_values_for_dr.push_back(el->pt()/1000.0);
                    }
                    if (higgs_counter == 2) {
                        higgs_2_pt = el->pt()/1000.0;
                        higgs_2_eta = el->eta();
                        higgs_2_phi = el->phi();
                        vector<vector<float>> b2_list;
                        find_non_higgs_daughters(el, b2_list);
                        b12_list = b2_list[0];
                        b22_list = b2_list[1];
                        float dR_2 = delta_R(b12_list[1], b12_list[2], b22_list[1], b22_list[2]);
                        deltaR_values.push_back(dR_2);
                        higgs_pt_values_for_dr.push_back(el->pt()/1000.0);
                    }

                    float higgs_average_pt = (higgs_1_pt + higgs_2_pt) / 2.0;
                    higgs_pt_values.push_back(el->pt()/1000.0);
                    higgs_pt_values_by_event.push_back(el->pt()/1000.0);
                    std::cout << "higgs pt value: " << el->pt()/1000.0 << "\n";
                }
            }

            // --- After TruthParticles loop, build full b-list ---
            allb_list.push_back(b11_list);
            allb_list.push_back(b21_list);
            allb_list.push_back(b12_list);
            allb_list.push_back(b22_list);

            /*
            for (const auto& p : b11_list) allb_list.push_back({p[0], p[1], p[2]});
            for (const auto& p : b21_list) allb_list.push_back({p[0], p[1], p[2]});
            for (const auto& p : b12_list) allb_list.push_back({p[0], p[1], p[2]});
            for (const auto& p : b22_list) allb_list.push_back({p[0], p[1], p[2]});
            */

            std::sort(allb_list.begin(), allb_list.end(),
                [](const std::vector<float>& a, const std::vector<float>& b) {
                    return a[0] > b[0]; // Sort by Et descending
                });



            // --- Delta R matching between b-jets and truth jets / gFEX jets ---
            int b_iterator = 0;
            for (const auto& b : allb_list) {
                float b_Et  = b[0];
                float b_eta = b[1];
                float b_phi = b[2];
                std::vector<float> deltaR_truth_jets, deltaR_gFex_smallR_jets;

                for (const auto& jet : *TruthJets) {
                    deltaR_truth_jets.push_back(delta_R(b_eta, b_phi, jet->eta(), jet->phi()));
                }
                /*
                for (int i = 0 ; i < gFexSRJets->size(); i++){
                    deltaR_gFex_smallR_jets.push_back(delta_R(b_eta, b_phi, gFexSRJets[i].eta(), gFexSRJets[i].phi()));
                    b_gfex_deltar_list.push_back(delta_R(b_eta, b_phi, gFexSRJets[i].eta(), gFexSRJets[i].phi()));
                }*/

                for (const auto& jet : *gFexLRJets) {
                    deltaR_gFex_smallR_jets.push_back(delta_R(b_eta, b_phi, jet->eta(), jet->phi()));
                    b_gfex_deltar_list.push_back(delta_R(b_eta, b_phi, jet->eta(), jet->phi()));
                }

                std::sort(deltaR_truth_jets.begin(), deltaR_truth_jets.end());
                std::sort(deltaR_gFex_smallR_jets.begin(), deltaR_gFex_smallR_jets.end());

                if (!deltaR_truth_jets.empty() && !deltaR_gFex_smallR_jets.empty()) {
                    if (b_iterator == 0) {
                        highestEt_b_closest_truthjet_deltar_list.push_back(deltaR_truth_jets[0]);
                        highestEt_b_closest_gfex_deltar_list.push_back(deltaR_gFex_smallR_jets[0]);
                    }
                    if (b_iterator == 1) {
                        second_highestEt_b_closest_truthjet_deltar_list.push_back(deltaR_truth_jets[0]);
                        second_highestEt_b_closest_gfex_deltar_list.push_back(deltaR_gFex_smallR_jets[0]);
                    }
                    if (b_iterator == 2) {
                        third_highestEt_b_closest_truthjet_deltar_list.push_back(deltaR_truth_jets[0]);
                        third_highestEt_b_closest_gfex_deltar_list.push_back(deltaR_gFex_smallR_jets[0]);
                    }
                    if (b_iterator == 3) {
                        fourth_highestEt_b_closest_truthjet_deltar_list.push_back(deltaR_truth_jets[0]);
                        fourth_highestEt_b_closest_gfex_deltar_list.push_back(deltaR_gFex_smallR_jets[0]);
                    }
                }
                b_iterator++;
            }

            // --- Higgs Pt Cut ---
            if (higgsPtCut && signalBool) {
                if (higgs_pt_values_by_event[0] > 100.0 || higgs_pt_values_by_event[1] > 100.0){
                    higgs_passes_cut = true;
                    higgs_pt_values_after_cut.insert(higgs_pt_values_after_cut.end(), higgs_pt_values_by_event.begin(), higgs_pt_values_by_event.end());
                    for (const auto& b : allb_list) {
                        b_Et_values_after_higgs_cut.push_back(b[0]);
                    }
                }
                /*
                if (std::any_of(higgs_pt_values_by_event.begin(), higgs_pt_values_by_event.end(),
                                [&](float pt) { return (pt > higgs_min_pt); })) {
                    higgs_pt_values_after_cut.insert(higgs_pt_values_after_cut.end(),
                                                    higgs_pt_values_by_event.begin(), higgs_pt_values_by_event.end());
                    higgs_passes_cut = true;
                    for (const auto& b : allb_list) {
                        b_Et_values_after_higgs_cut.push_back(b[0]);
                    }*/
                        
                }
            }
        f->Close();
    }


    // --- PLOTTING STARTS ---

    // Create output directory if needed
    std::string outdir = signalBool ? "signalTruthPlotsROOT/" : "backgroundTruthPlotsROOT/";
    gSystem->mkdir(outdir.c_str(), true);

    if (signalBool) {
        // --- Higgs pT Distribution ---
        TH1F* h_higgs_pt = new TH1F("h_higgs_pt", "Higgs Truth p_{T} Distribution;p_{T} [GeV];# Higgs", 50, 0, 500);
        for (auto pt : higgs_pt_values) h_higgs_pt->Fill(pt);

        TCanvas* c_higgs_pt = new TCanvas("c_higgs_pt", "Higgs pT", 800, 600);
        h_higgs_pt->Draw("hist");
        c_higgs_pt->SaveAs((outdir + "higgs_pt_distribution.pdf").c_str());

        // --- Higgs pT after Cut ---
        if (higgsPtCut) {
            TH1F* h_higgs_pt_after_cut = new TH1F("h_higgs_pt_after_cut", "Higgs Truth p_{T} After Cut;p_{T} [GeV];# Higgs", 50, 0, 500);
            for (auto pt : higgs_pt_values_after_cut) h_higgs_pt_after_cut->Fill(pt);

            TCanvas* c_higgs_pt_cut = new TCanvas("c_higgs_pt_cut", "Higgs pT Cut", 800, 600);
            h_higgs_pt_after_cut->Draw("hist");
            c_higgs_pt_cut->SaveAs((outdir + "higgs_pt_distribution_after_cut.pdf").c_str());
        }

        // --- b-quark E_T Distribution ---
        TH1F* h_b_Et = new TH1F("h_b_Et", "b Truth E_{T} Distribution;E_{T} [GeV];# b quarks", 50, 0, 500);
        for (auto et : b_Et_values) h_b_Et->Fill(et);

        TCanvas* c_b_Et = new TCanvas("c_b_Et", "b ET", 800, 600);
        h_b_Et->Draw("hist");
        c_b_Et->SaveAs((outdir + "b_Et_distribution.pdf").c_str());

        // --- b-quark E_T after Higgs Cut ---
        TH1F* h_b_Et_after_cut = new TH1F("h_b_Et_after_cut", "b Truth E_{T} After Higgs Cut;E_{T} [GeV];# b quarks", 50, 0, 500);
        for (auto et : b_Et_values_after_higgs_cut) h_b_Et_after_cut->Fill(et);

        TCanvas* c_b_Et_cut = new TCanvas("c_b_Et_cut", "b ET Cut", 800, 600);
        h_b_Et_after_cut->Draw("hist");
        c_b_Et_cut->SaveAs((outdir + "b_Et_distribution_after_higgs_cut.pdf").c_str());

        // --- DeltaR_bb vs Higgs pT ---
        const int nBins = 22;
        double ptBins[nBins+1];
        for (int i = 0; i <= nBins; ++i) ptBins[i] = 40.0 + i * 20.0;

        std::vector<double> avgDeltaR;
        std::vector<double> errors;

        for (int i = 0; i < nBins; ++i) {
            std::vector<double> deltaRsInBin;
            for (size_t j = 0; j < higgs_pt_values_for_dr.size(); ++j) {
                if (higgs_pt_values_for_dr[j] >= ptBins[i] && higgs_pt_values_for_dr[j] < ptBins[i+1]) {
                    deltaRsInBin.push_back(deltaR_values[j]);
                }
            }
            if (!deltaRsInBin.empty()) {
                double mean = std::accumulate(deltaRsInBin.begin(), deltaRsInBin.end(), 0.0) / deltaRsInBin.size();
                double sq_sum = std::inner_product(deltaRsInBin.begin(), deltaRsInBin.end(), deltaRsInBin.begin(), 0.0);
                double std_dev = sqrt(sq_sum / deltaRsInBin.size() - mean * mean);
                avgDeltaR.push_back(mean);
                errors.push_back(std_dev / sqrt(deltaRsInBin.size()));
            } else {
                avgDeltaR.push_back(0);
                errors.push_back(0);
            }
        }

        double binCenters[nBins];
        for (int i = 0; i < nBins; ++i) binCenters[i] = (ptBins[i] + ptBins[i+1]) / 2;

        TGraphErrors* gr_deltaR = new TGraphErrors(nBins, binCenters, avgDeltaR.data(), 0, errors.data());
        gr_deltaR->SetTitle("Average #DeltaR vs Higgs p_{T};Higgs p_{T} [GeV];Average #DeltaR");
        gr_deltaR->SetMarkerStyle(20);
        gr_deltaR->SetMarkerSize(1);
        gr_deltaR->SetLineWidth(2);

        TCanvas* c_deltaR = new TCanvas("c_deltaR", "DeltaR vs Higgs pT", 800, 600);
        gr_deltaR->Draw("AP");
        c_deltaR->SaveAs((outdir + "deltaRbb_vs_higgs_pt.pdf").c_str());

        // --- DeltaR distributions ---
        auto drawDeltaRHist = [&](const std::vector<float>& data, const std::string& histname, const std::string& title, const std::string& filename) {
            TH1F* h = new TH1F(histname.c_str(), title.c_str(), 25, 0, 5);
            for (auto v : data) h->Fill(v);
            h->Scale(1.0 / h->Integral());
            TCanvas* c = new TCanvas((histname + "_c").c_str(), histname.c_str(), 800, 600);
            h->Draw("hist");
            c->SaveAs((outdir + filename).c_str());
        };

        drawDeltaRHist(b_gfex_deltar_list, "h_b_gfex", "#DeltaR (b quarks and gFEX smallR jets);#DeltaR;Normalized Events", "deltaRb_gFex_smallR_distribution.pdf");
        drawDeltaRHist(highestEt_b_closest_gfex_deltar_list, "h_highestEt_b_gfex", "#DeltaR (highest E_{T} b and closest gFEX smallR jet);#DeltaR;Normalized Events", "highestEt_b_closest_gFEX_dR.pdf");
        drawDeltaRHist(second_highestEt_b_closest_gfex_deltar_list, "h_second_b_gfex", "#DeltaR (second highest E_{T} b and closest gFEX smallR jet);#DeltaR;Normalized Events", "second_highestEt_b_closest_gFEX_dR.pdf");

        // --- Truth Jet pT ---
        TH1F* h_truthJet_pt = new TH1F("h_truthJet_pt", "Truth InTimeAntiKt4TruthJets p_{T};p_{T} [GeV];# Jets", 50, 0, 100);
        for (auto pt : jet_pt_values) h_truthJet_pt->Fill(pt);

        TCanvas* c_truthJet_pt = new TCanvas("c_truthJet_pt", "Truth Jet pt", 800, 600);
        h_truthJet_pt->Draw("hist");
        c_truthJet_pt->SaveAs((outdir + "truth_jet_pt_distribution.pdf").c_str());

        // --- Truth Jet pT with Log Scale ---
        TH1F* h_truthJet_pt_log = new TH1F("h_truthJet_pt_log", "Truth InTimeAntiKt4TruthJets p_{T} (Log Y);p_{T} [GeV];# Jets", 50, 0, 400);
        for (auto pt : jet_pt_values) h_truthJet_pt_log->Fill(pt);

        TCanvas* c_truthJet_log = new TCanvas("c_truthJet_log", "Truth Jet pt Log", 800, 600);
        c_truthJet_log->SetLogy();
        h_truthJet_pt_log->Draw("hist");
        c_truthJet_log->SaveAs((outdir + "truth_jet_pt_log_distribution.pdf").c_str());
    }
    else {
        // --- Background Plotting ---

        // Create output directory if needed
        std::string outdir = "backgroundTruthPlots/";
        gSystem->mkdir(outdir.c_str(), true);

        // --- Truth Jets pt Distribution ---
        TH1F* h_truthJet_pt_bg = new TH1F("h_truthJet_pt_bg", "Truth InTimeAntiKt4TruthJets p_{T} Distribution;p_{T} [GeV];# Jets", 50, 0, 100);
        for (auto pt : jet_pt_values) h_truthJet_pt_bg->Fill(pt);

        TCanvas* c_truthJet_pt_bg = new TCanvas("c_truthJet_pt_bg", "Truth Jet pt Background", 800, 600);
        h_truthJet_pt_bg->Draw("hist");
        c_truthJet_pt_bg->SaveAs((outdir + "truth_jet_pt_distribution.pdf").c_str());

        // --- HLT Jets pt Distribution ---
        TH1F* h_hltJet_pt_bg = new TH1F("h_hltJet_pt_bg", "HLT_AntiKt4EMTopoJets_subjesIS p_{T} Distribution;p_{T} [GeV];# Jets", 50, 0, 100);
        for (auto pt : hlt_jet_pt_values) h_hltJet_pt_bg->Fill(pt);

        TCanvas* c_hltJet_pt_bg = new TCanvas("c_hltJet_pt_bg", "HLT Jet pt Background", 800, 600);
        h_hltJet_pt_bg->Draw("hist");
        c_hltJet_pt_bg->SaveAs((outdir + "hlt_jet_pt_distribution.pdf").c_str());

        // --- HLT Jets pt Distribution (Log Scale) ---
        TH1F* h_hltJet_pt_log_bg = new TH1F("h_hltJet_pt_log_bg", "HLT_AntiKt4EMTopoJets_subjesIS p_{T} Distribution (Log Y);p_{T} [GeV];# Jets", 50, 0, 400);
        for (auto pt : hlt_jet_pt_values) h_hltJet_pt_log_bg->Fill(pt);

        TCanvas* c_hltJet_log_bg = new TCanvas("c_hltJet_log_bg", "HLT Jet pt Log Background", 800, 600);
        c_hltJet_log_bg->SetLogy();
        h_hltJet_pt_log_bg->Draw("hist");
        c_hltJet_log_bg->SaveAs((outdir + "hlt_jet_pt_log_distribution.pdf").c_str());
    }



    cout << "Processing complete." << endl;
}
