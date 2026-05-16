#!/usr/bin/env python

# Set up ROOT and RootCore
import ROOT
import os
import numpy as np
import matplotlib.pyplot as plt
import math
import subprocess
import heapq
import matplotlib.patches as patches

# Configuration booleans
signalBool = True
lxplusBool = True
higgsPtCut = True
jFexPlotsBool = False

# Clear plots
subprocess.run("rm -rf signalEventPlots/*", shell=True, check=True)

def find_non_higgs_daughters(particle, result=None, level=0):
    """ Recursively find the first non-Higgs daughters of a given particle. """
    if result is None:
        result = []  # Initialize the result list only once

    if particle.pdgId() != 25:  # If not a Higgs, store properties
        E = particle.e()  # Total energy
        p_T = particle.pt()  # Transverse momentum
        eta = particle.eta()
        phi = particle.phi()
        
        # Compute total momentum p = p_T * cosh(eta)
        p = p_T * math.cosh(eta) if abs(eta) < 100 else float('inf')

        # Compute transverse energy E_T
        E_T = E * (p_T / p) if p > 0 else 0  # Avoid division by zero

        result.append((E_T / 1000.0, eta, phi))  # Store in GeV
        return result  # Return the list with appended values
    
    # If it's a Higgs, continue to its children
    for i in range(particle.nChildren()):
        find_non_higgs_daughters(particle.child(i), result, level + 1)  # Pass result to child
    
    return result

def delta_phi(phi1, phi2):
    """Compute Δφ, ensuring it stays within [-π, π]"""
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2 * math.pi
    while dphi < -math.pi:
        dphi += 2 * math.pi
    return dphi

def delta_R(eta1, phi1, eta2, phi2):
    """Compute ΔR = sqrt((Δη)^2 + (Δφ)^2)"""
    return math.sqrt((eta1 - eta2) ** 2 + delta_phi(phi1, phi2) ** 2)

def reset_event_counter(eventCounter): # reset event counter to keep track of kept event number after higgs cuts
    if eventCounter == 0:
        eventCounter = 0
    else: 
        eventCounter -= 1

# Plotting the heatmaps
def plot_heatmap(objectType, data, eta_bins, phi_bins, title, filename, log=False,
                 Higgs1Pt=0, Higgs2Pt=0, Higgs1Eta=0, Higgs2Eta=0,
                 Higgs1Phi=0, Higgs2Phi=0, signalBool=False,
                 b1_list=[], b2_list=[], dr1=0, dr2=0):

    plt.figure(figsize=(8, 6))
    plt.imshow(data.T, origin='lower', aspect='auto',
               extent=[eta_bins[0], eta_bins[-1], phi_bins[0], phi_bins[-1]],
               cmap='viridis')
    plt.subplots_adjust(top=0.77)

    if log:
        plt.colorbar(label='Log₁₀(Sum of $E_T$)')
    else:
        plt.colorbar(label='Sum of $E_T$ [GeV]')

    plt.xlabel('Eta')
    plt.ylabel('Phi')
    plt.title(title)

    bin_eta_centers = 0.5 * (eta_bins[:-1] + eta_bins[1:])
    bin_phi_centers = 0.5 * (phi_bins[:-1] + phi_bins[1:])
    
    for i, eta in enumerate(bin_eta_centers):
        for j, phi in enumerate(bin_phi_centers):
            if objectType == "gfex_smallr":
                if data[i][j] >= 20:
                    rect = patches.Rectangle((eta - 0.3, phi - 0.3), 0.6, 0.6,
                                            edgecolor='white', facecolor='none', lw=1.5)
                    plt.gca().add_patch(rect)
            elif objectType == "gfex_larger":
                if data[i][j] >= 20:
                    circle = patches.Circle((eta, phi), radius=0.9,
                                            edgecolor='white', facecolor='none', lw=1.5)
                    plt.gca().add_patch(circle)
            elif objectType == "truthjet" or objectType == "hltjet":
                if data[i][j] >= 20:
                    circle = patches.Circle((eta, phi), radius=0.4,
                                            edgecolor='white', facecolor='none', lw=1.5)
                    plt.gca().add_patch(circle)

    if signalBool:
        b_colors = ['red', 'blue', 'green', 'orange']
        higgs_colors = ['magenta', 'cyan']  # Color for Higgs circles and legend

        # Unpack b quark info
        (b1_Et_1, b1_eta_1, b1_phi_1) = b1_list[0]
        (b2_Et_1, b2_eta_1, b2_phi_1) = b1_list[1]
        (b1_Et_2, b1_eta_2, b1_phi_2) = b2_list[0]
        (b2_Et_2, b2_eta_2, b2_phi_2) = b2_list[1]

        # Draw b-quark circles
        b_positions = [
            ((b1_eta_1, b1_phi_1), b_colors[0]),
            ((b2_eta_1, b2_phi_1), b_colors[1]),
            ((b1_eta_2, b1_phi_2), b_colors[2]),
            ((b2_eta_2, b2_phi_2), b_colors[3])
        ]
        for (eta, phi), color in b_positions:
            circle = patches.Circle((eta, phi), radius=0.1, edgecolor=color,
                                    facecolor='none', lw=1.5, linestyle='--')
            plt.gca().add_patch(circle)

        # Draw Higgs centers
        higgs_positions = [
            (Higgs1Eta, Higgs1Phi, higgs_colors[0]),
            (Higgs2Eta, Higgs2Phi, higgs_colors[1])
        ]
        for eta, phi, color in higgs_positions:
            circle = patches.Circle((eta, phi), radius=0.1, edgecolor=color,
                                    facecolor='none', lw=1.5, linestyle='-')
            plt.gca().add_patch(circle)

        # Create legend text with color coding
        higgs1_lines = [
            f"  Higgs 1: pT={Higgs1Pt:.2f} GeV, Eta={Higgs1Eta:.2f}, Phi={Higgs1Phi:.2f}",
            f"  1st b: Et={b1_Et_1:.2f}, Eta={b1_eta_1:.2f}, Phi={b1_phi_1:.2f}",
            f"  2nd b: Et={b2_Et_1:.2f}, Eta={b2_eta_1:.2f}, Phi={b2_phi_1:.2f}",
            f"  dr (b's)={dr1:.2f}"
        ]
        higgs1_colors = [higgs_colors[0], b_colors[0], b_colors[1], None]

        higgs2_lines = [
            f"  Higgs 2: pT={Higgs2Pt:.2f} GeV, Eta={Higgs2Eta:.2f}, Phi={Higgs2Phi:.2f}",
            f"  1st b: Et={b1_Et_2:.2f}, Eta={b1_eta_2:.2f}, Phi={b1_phi_2:.2f}",
            f"  2nd b: Et={b2_Et_2:.2f}, Eta={b2_eta_2:.2f}, Phi={b2_phi_2:.2f}",
            f"  dr (b's)={dr2:.2f}"
        ]
        higgs2_colors = [higgs_colors[1], b_colors[2], b_colors[3], None]

        # Place text annotations on figure
        y_start = 0.96
        y_step = 0.045

        for i, (line, color) in enumerate(zip(higgs1_lines, higgs1_colors)):
            plt.figtext(0.01, y_start - i * y_step, line, ha="left", fontsize=7.5,
                        color=color if color else "black")

        for i, (line, color) in enumerate(zip(higgs2_lines, higgs2_colors)):
            plt.figtext(0.51, y_start - i * y_step, line, ha="left", fontsize=7.5,
                        color=color if color else "black")

    plt.savefig(filename)
    plt.close()

# Initialize the xAOD infrastructure
ROOT.gSystem.Load("libxAODRootAccess")
ROOT.xAOD.Init()

# Set up the input file directory
if signalBool:
    if lxplusBool:
        fileDir = "/data/larsonma/LargeRadiusJets/datasets/Signal_HHbbbb/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130/"
    else:
        fileDir = "/data/larsonma/LargeRadiusJets/datasets/Signal_HHbbbb/mc21_14TeV.537540.MGPy8EG_hh_bbbb_vbf_novhh_5fs_l1cvv1cv1.recon.AOD.e8557_s4422_r16130/"
    
else: 
    if lxplusBool:
        fileDir = "/data/larsonma/LargeRadiusJets/datasets/Background_jj_JZ3/mc21_14TeV.801168.Py8EG_A14NNPDF23LO_jj_JZ3.recon.AOD.e8557_s4422_r16130"
    else:
        fileDir = "/data/larsonma/LargeRadiusJets/datasets/Background_jj_JZ3/mc21_14TeV.801168.Py8EG_A14NNPDF23LO_jj_JZ3.recon.AOD.e8557_s4422_r16130"
    
# Get a list of all ROOT files in the directory
fileNames = [os.path.join(fileDir, f) for f in os.listdir(fileDir) if f.endswith(".pool.root.1")]

if not fileNames:
    print("No ROOT files found in the specified directory.")
    exit(1)

# Initialize arrays to store phi, eta, and Et values

### ADD TO EVENT DISPLAY and kinematic distributions: HLT_AntiKt4EMTopoJets_subjesIS
#gfex_smallr_et_values = []
"""
topo_phi_values = []
topo_eta_values = []
topo_et_values = []

topo_422_phi_values = []
topo_422_eta_values = []
topo_422_et_values = []

gfex_smallr_phi_values = []
gfex_smallr_eta_values = []


jfex_smallr_phi_values = []
jfex_smallr_eta_values = []
jfex_smallr_et_values = []
"""
higgs_pt_values = []  # List to store plotted values
deltaR_values = []
higgs_pt_values_for_dr = []
jet_pt_values = []
jet_log_pt_values = []
hlt_jet_pt_values = []
gfex_larger_jet_pt_values = []
higgs_pt_values_after_cut = []
"""
highestEt_b_closest_gfex_deltar_list = []
second_highestEt_b_closest_gfex_deltar_list = []
third_highestEt_b_closest_gfex_deltar_list = []
fourth_highestEt_b_closest_gfex_deltar_list = []

highestEt_b_closest_truthjet_deltar_list = []
second_highestEt_b_closest_truthjet_deltar_list = []
third_highestEt_b_closest_truthjet_deltar_list = []
fourth_highestEt_b_closest_truthjet_deltar_list = []
"""
#b11_list = []
#b21_list = []
#b12_list = []
#b22_list = []

b_Et_values = []
b_Et_values_after_higgs_cut = []

# Heatmap bin definitions
phi_bins = np.linspace(-3.2, 3.2, 65)  # 64 bins
eta_bins = np.linspace(-5.0, 5.0, 101)  # 100 bins
if signalBool:
    if lxplusBool:
        output_file_topo422 = "/data/larsonma/LargeRadiusJets/MemPrints/CaloTopo_422/mc21_14TeV_hh_bbbb_vbf_novhh_topo422.dat"
        output_file_gfex = "/data/larsonma/LargeRadiusJets/MemPrints/gFex/mc21_14TeV_hh_bbbb_vbf_novhh_gfex_smallrj.dat"
        #output_file_jfex = "/data/larsonma/LargeRadiusJets/MemPrints/jFex/mc21_14TeV_hh_bbbb_vbf_novhh_jfex_smallrj.dat"
    else:
        output_file_topo422 = "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/CaloTopo_422/mc21_14TeV_hh_bbbb_vbf_novhh_topo422.dat"
        output_file_gfex = "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/gFex/mc21_14TeV_hh_bbbb_vbf_novhh_gfex_smallrj.dat"
        #output_file_jfex = "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/jFex/mc21_14TeV_hh_bbbb_vbf_novhh_jfex_smallrj.dat"
    
else:
    if lxplusBool:
        output_file_topo422 = "/data/larsonma/LargeRadiusJets/MemPrints/CaloTopo_422/mc21_14TeV_jj_JZ3_topo422.dat"
        output_file_gfex = "/data/larsonma/LargeRadiusJets/MemPrints/gFex/mc21_14TeV_jj_JZ3_gfex_smallrj.dat"
        #output_file_jfex = "/data/larsonma/LargeRadiusJets/MemPrints/jFex/mc21_14TeV_jj_JZ3_jfex_smallrj.dat"
    else:
        output_file_topo422 = "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/CaloTopo_422/mc21_14TeV_jj_JZ3_topo422.dat"
        output_file_gfex = "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/gFex/mc21_14TeV_jj_JZ3_gfex_smallrj.dat"
        #output_file_jfex = "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/jFex/mc21_14TeV_jj_JZ3_jfex_smallrj.dat"
    

# Function to scale and digitize a value
def digitize(value, bit_length, min_val, max_val):
    # Ensure the value is within range
    if value < min_val or value > max_val:
        raise ValueError(f"Value {value} is out of range ({min_val}, {max_val})")
    scale = (2**bit_length - 1) / (max_val - min_val)
    return int(round((value - min_val) * scale))

# Parameters
phi_bit_length = 6
phi_min = -3.2
phi_max = 3.2

eta_bit_length = 8
eta_min = -5.0
eta_max = 5.0

et_bit_length = 13
et_min = 0
et_max = 2048
et_granularity = 0.25  # GeV per bit

events_per_file = 100

higgs_min_pt = 100
higgs_max_pt = 5000000 # remove max pt

higgs_min_average_pt = 200

fileIt = 0
eventCounter = -1
higgsPlotCounter = 0
higgs_counter_average = 0
higgs_counter = 0

# Process each file
with open(output_file_topo422, "w") as f_topo:
    with open(output_file_gfex, "w") as f_gfex:
            
        for fileName in fileNames:
            print("fileit,", fileIt)
            if fileIt > 60:
                break
            print(f"Processing file: {fileName}")

            # Open the file and make the "transient tree"
            f = ROOT.TFile.Open(fileName)
            t = ROOT.xAOD.MakeTransientTree(f, 'CollectionTree')

            # Print the number of events in the file
            print('  Number of input events: %s' % t.GetEntries())
            # Loop through all events in the file
            for entry in range(t.GetEntries()):
                # Initialize variables used later in plotting
                higgs_1_pt = 0
                higgs_2_pt = 0
                higgs_1_eta = 0
                higgs_2_eta = 0
                higgs_1_phi = 0
                higgs_2_phi = 0
                higgs_counter_average += higgs_counter
                higgs_counter = 0
                higgs_average_pt = 0
                b1_Et_1 = b1_eta_1 = b1_phi_1 = b1_Et_2 = b1_eta_2 = b1_phi_2 =None
                b2_Et_1 = b2_eta_1 = b2_phi_1 = b2_Et_2 = b2_eta_2 = b2_phi_2 = None
                dR_1 = 0
                dR_2 = 0
                b1_list = [] # list for bbbar from 1st higgs
                b2_list = [] # list for bbbar from 2nd higgs
                allb_list = []
                sorted_b_list = []
                
                eventCounter += 1
                higgs_pt_values_by_event = []  # List to store pt values
                higgs_passes_cut = False
                t.GetEntry(entry)
                iEvt = eventCounter # 100 events per file
                #print("entry, fileIt, iEvt:", eventCounter, fileIt, iEvt)
                gfex_smallr_et_values_by_event = []
                """
                topo_phi_values_by_event = []
                topo_eta_values_by_event = []
                topo_et_values_by_event = []

                topo_422_phi_values_by_event = []
                topo_422_eta_values_by_event = []
                topo_422_et_values_by_event = []

                gfex_smallr_phi_values_by_event = []
                gfex_smallr_eta_values_by_event = []
                

                jfex_smallr_phi_values_by_event = []
                jfex_smallr_eta_values_by_event = []
                jfex_smallr_et_values_by_event = []
                """
                # Heatmap storage
                topo_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                topo_422_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                gfex_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                #jfex_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                truthInTimeJets_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                hltJets_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                gfex_largerJets_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))

                topo_log_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                topo_422_log_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                gfex_log_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                #jfex_log_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                truthInTimeJets_log_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                hltJets_log_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                gfex_largerJets_log_heatmap = np.zeros((len(eta_bins) - 1, len(phi_bins) - 1))
                
                if iEvt == 0: 
                    # Loop over all collections in the transient tree (only for the first event of the first file)
                    for key in t.GetListOfBranches():
                        collection_name = key.GetName()
                        collection = getattr(t, collection_name, None)

                        if collection:
                            print(f'  Collection: {collection_name}')

                sum_et = 0.0

                topo_it = 0

                #if iEvt <= 9:
                truthIt = 0

                for el in t.HLT_AntiKt4EMTopoJets_subjesIS:
                    hlt_jet_pt_values.append(el.pt() / 1000)
                    hlt_jet_total_p = el.pt() * math.cosh(el.eta())
                    hlt_jet_ptoverp = el.pt() / hlt_jet_total_p
                    hlt_jet_Et = (el.e() * hlt_jet_ptoverp )/1000.0
                    hltJets_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[hlt_jet_Et])[0]
                    hltJets_log_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[np.log10(hlt_jet_Et)])[0]

                for el in t.L1_gFexLRJetRoI:
                    gfex_larger_jet_pt_values.append(el.et() / 1000)
                    #gfex_larger_jet_total_p = el.pt() * math.cosh(el.eta())
                    #gfex_larger_jet_ptoverp = el.pt() / hlt_jet_total_p
                    #gfex_larger_jet_Et = (el.e() * hlt_jet_ptoverp )/1000.0
                    gfex_larger_jet_Et = el.et() / 1000
                    gfex_largerJets_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[gfex_larger_jet_Et])[0]
                    gfex_largerJets_log_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[np.log10(gfex_larger_jet_Et)])[0]

                for el in t.InTimeAntiKt4TruthJets:
                    #print("truth jet pt:", el.pt() / 1000)
                    jet_pt_values.append(el.pt() / 1000)
                    jet_log_pt_values.append(np.log10(el.pt() / 1000))
                    truth_jet_total_p = el.pt() * math.cosh(el.eta())
                    truth_jet_ptoverp = el.pt() / truth_jet_total_p
                    truth_jet_Et = (el.e() * truth_jet_ptoverp )/1000.0
                    truthInTimeJets_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[truth_jet_Et])[0]
                    truthInTimeJets_log_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[np.log10(truth_jet_Et)])[0]
                for el in t.TruthParticles: ### FIXME add check to ensure the b's are actually b's
                    if el.pdgId() == 25 and el.status() == 22:
                        higgs_counter += 1
                        if higgs_counter == 1:
                            total_p = el.pt() * math.cosh(el.eta())
                            ptoverp = el.pt() / total_p
                            #print(' Higgs 1 PDG ID = %d, Eta = %g, Phi = %g, Et [GeV] = %g, Pt [GeV] = %g, iEvt = %g, status = %g' %  
                            #(el.pdgId(), el.eta(), el.phi(), (el.e() * ptoverp )/1000.0, el.pt()/1000, iEvt, el.status()))
                            higgs_1_pt = el.pt()/1000
                            higgs_1_eta = el.eta()
                            higgs_1_phi = el.phi()
                            b1_list = find_non_higgs_daughters(el)
                            #print("len(b1_list):", len(b1_list))
                            (b1_Et_1, b1_eta_1, b1_phi_1) = b1_list[0]
                            #b11_list.append([b1_eta_1, b1_phi_1])
                            (b2_Et_1, b2_eta_1, b2_phi_1) = b1_list[1]
                            #b21_list.append([b2_eta_1, b2_phi_1])
                            dR_1 = delta_R(b1_eta_1, b1_phi_1, b2_eta_1, b2_phi_1)
                            #if (dR_1 < 1.0):
                            #    print("dr1 is < 1.0!")
                            #    print("dr:", dR_1)
                            deltaR_values.append(dR_1)
                            higgs_pt_values_for_dr.append(el.pt() / 1000)

                        if higgs_counter == 2:
                            total_p = el.pt() * math.cosh(el.eta())
                            ptoverp = el.pt() / total_p
                            #print(' Higgs 2 PDG ID = %d, Eta = %g, Phi = %g, Et [GeV] = %g, Pt [GeV] = %g, iEvt = %g, status = %g' %  
                            #(el.pdgId(), el.eta(), el.phi(), (el.e() * ptoverp )/1000.0, el.pt()/1000, iEvt, el.status()))
                            higgs_2_pt = el.pt()/1000
                            higgs_2_eta = el.eta()
                            higgs_2_phi = el.phi()
                            b2_list = find_non_higgs_daughters(el)
                            #print("len(b2_list):", len(b2_list))
                            (b1_Et_2, b1_eta_2, b1_phi_2) = b2_list[0]
                            #b12_list.append([b1_eta_2, b1_phi_2])
                            (b2_Et_2, b2_eta_2, b2_phi_2) = b2_list[1]
                            #b22_list.append([b2_eta_2, b2_phi_2])
                            dR_2 = delta_R(b1_eta_2, b1_phi_2, b2_eta_2, b2_phi_2)
                            #if (dR_2 < 1.0):
                            #    print("dr2 is < 1.0!")
                            #    print("dr:", dR_2)
                            deltaR_values.append(dR_2)
                            higgs_pt_values_for_dr.append(el.pt() / 1000)

                        higgs_average_pt = (higgs_1_pt + higgs_2_pt) / 2.0
                        higgs_pt_values.append(el.pt()/1000)
                        higgs_pt_values_by_event.append(el.pt()/1000)
                    truthIt += 1
                    
                if len(b1_list) > 1:
                    allb_list.append(b1_list[0])
                    allb_list.append(b1_list[1])
                if len(b2_list) > 1:
                    allb_list.append(b2_list[0])
                    allb_list.append(b2_list[1])

                sorted_b_list = sorted(allb_list, key=lambda x: x[0], reverse=True)

                b_iterator = 0
                for (b_Et, b_eta, b_phi) in sorted_b_list:
                    deltaR_gFex_smallR_jets = []
                    deltaR_truth_jets = []
                    for el in t.InTimeAntiKt4TruthJets:
                        deltaR_truth_jets.append(delta_R(b_eta, b_phi, el.eta(), el.phi()))
                    for el in t.L1_gFexSRJetRoI:
                        deltaR_gFex_smallR_jets.append(delta_R(b_eta, b_phi, el.eta(), el.phi()))
                    sorted_deltaR_truth_jets = sorted(deltaR_truth_jets)#, reverse=True)
                    sorted_deltaR_gFex_smallR_jets = sorted(deltaR_gFex_smallR_jets)#, reverse=True)
                    """
                    if len(sorted_deltaR_truth_jets) > 0: 
                        if b_iterator == 0:
                            highestEt_b_closest_gfex_deltar_list.append(sorted_deltaR_gFex_smallR_jets[0])
                            highestEt_b_closest_truthjet_deltar_list.append(sorted_deltaR_truth_jets[0])
                        if b_iterator == 1:
                            second_highestEt_b_closest_gfex_deltar_list.append(sorted_deltaR_gFex_smallR_jets[0])
                            second_highestEt_b_closest_truthjet_deltar_list.append(sorted_deltaR_truth_jets[0])
                        if b_iterator == 2:
                            third_highestEt_b_closest_gfex_deltar_list.append(sorted_deltaR_gFex_smallR_jets[0])
                            third_highestEt_b_closest_truthjet_deltar_list.append(sorted_deltaR_truth_jets[0])
                        if b_iterator == 3:
                            fourth_highestEt_b_closest_gfex_deltar_list.append(sorted_deltaR_gFex_smallR_jets[0])
                            fourth_highestEt_b_closest_truthjet_deltar_list.append(sorted_deltaR_truth_jets[0])
                        b_iterator += 1
                    """
                for (b_Et, b_eta, b_phi) in allb_list:
                    b_Et_values.append(b_Et)     
                if higgsPtCut and signalBool:
                        
                    if any(higgs_pt > higgs_min_pt and higgs_pt < higgs_max_pt for higgs_pt in higgs_pt_values_by_event):
                        higgs_pt_values_after_cut += higgs_pt_values_by_event
                        higgs_passes_cut = True
                        for (b_Et, b_eta, b_phi) in allb_list:
                            b_Et_values_after_higgs_cut.append(b_Et)    
                
                    if not higgs_passes_cut:
                        eventCounter -= 1
                        continue

                # Loop over the CaloCalTopoClusters collection
                for el in t.CaloCalTopoClusters:
                    #print('  topo eta = %g, phi = %g, Et = %g' %  (el.eta(), el.phi(), el.et()))
                    #topo_phi_values_by_event.append(el.phi())
                    #topo_eta_values_by_event.append(el.eta())
                    #topo_et_values_by_event.append(el.et())
                    if iEvt < 25:
                        topo_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[el.et()/1000])[0]
                        topo_log_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[np.log10(el.et() / 1000)])[0]
                    # Initialize sum of transverse energy for this event
                    topo_it += 1

                #topo_phi_values.append(topo_phi_values_by_event)
                #topo_eta_values.append(topo_eta_values_by_event)
                #topo_et_values.append(topo_et_values_by_event)
                
                # Loop over the CaloCalTopoClusters collection
                topo422_it = 0
                for el in t.CaloTopoClusters422:
                    #print('  topo422 eta = %g, phi = %g, Et = %g' %  (el.eta(), el.phi(), el.et()))
                    #topo_422_phi_values_by_event.append(el.phi())
                    #topo_422_eta_values_by_event.append(el.eta())
                    #topo_422_et_values_by_event.append(el.et() / 1000.0)

                    # Digitize values
                    if el.et() >= 0:
                        phi_bin = digitize(el.phi(), phi_bit_length, phi_min, phi_max)
                        eta_bin = digitize(el.eta(), eta_bit_length, eta_min, eta_max)
                        #print("eta_bit_length:", eta_bit_length)
                        #print("phi_bit_length:", phi_bit_length)
                        et_bin = digitize(el.et() / (et_granularity * 1000), et_bit_length, et_min / et_granularity, et_max / et_granularity)
                        #print("topo422 phi_bin, eta_bin, et_bin:", phi_bin, eta_bin, et_bin)

                        # Create binary word
                        binary_word = f"{et_bin:0{et_bit_length}b}|{eta_bin:0{eta_bit_length}b}|{phi_bin:0{phi_bit_length}b}"

                        # Create hex word
                        hex_word = f"{(et_bin << (eta_bit_length + phi_bit_length)) | (eta_bin << phi_bit_length) | phi_bin:08x}"

                        # Write to file
                        if topo422_it == 0:
                            f_topo.write(f"Event : {iEvt}\n")
                        f_topo.write(f"0x{topo422_it:02x} {binary_word} 0x{hex_word}\n")
                    if iEvt < 25:
                        topo_422_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[el.et()/1000])[0]
                        topo_422_log_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[np.log10(el.et() / 1000)])[0]
                    sum_et += el.et()  # Add transverse energy of each cluster
                    topo422_it += 1
                    # Initialize sum of transverse energy for this event
                
                #topo_422_phi_values.append(topo_422_phi_values_by_event)
                #topo_422_eta_values.append(topo_422_eta_values_by_event)
                #topo_422_et_values.append(topo_422_et_values_by_event)

                gfex_it = 0
                for el in t.L1_gFexSRJetRoI:
                    #print('  GFEX: L1_gFexSRJetRoI eta = %g, phi = %g, Et = %g' %  (el.eta(), el.phi(), el.et()))
                    #gfex_smallr_phi_values_by_event.append(el.phi())
                    #gfex_smallr_eta_values_by_event.append(el.eta())
                    #gfex_smallr_et_values_by_event.append(el.et())

                    if el.et() >= 0:
                        # Digitize values
                        phi_bin = digitize(el.phi(), phi_bit_length, phi_min, phi_max)
                        eta_bin = digitize(el.eta(), eta_bit_length, eta_min, eta_max)
                        et_bin = digitize(el.et() / (et_granularity * 1000), et_bit_length, et_min / et_granularity, et_max / et_granularity)
                        #print("gfex phi_bin, eta_bin, et_bin:", phi_bin, eta_bin, et_bin)

                        # Create binary word
                        binary_word = f"{et_bin:0{et_bit_length}b}|{eta_bin:0{eta_bit_length}b}|{phi_bin:0{phi_bit_length}b}"

                        # Create hex word
                        hex_word = f"{(et_bin << (eta_bit_length + phi_bit_length)) | (eta_bin << phi_bit_length) | phi_bin:08x}"

                        # Write to file
                        if gfex_it == 0: 
                            f_gfex.write(f"Event : {iEvt}\n")
                            #print("writing event")
                        #print("binary word:", binary_word)
                        f_gfex.write(f"0x{gfex_it:02x} {binary_word} 0x{hex_word}\n")
                    if iEvt < 25:
                        gfex_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[el.et()/1000])[0]
                        gfex_log_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[np.log10(el.et() / 1000)])[0]
                    gfex_it += 1

                #gfex_smallr_phi_values.append(gfex_smallr_phi_values_by_event)
                #gfex_smallr_eta_values.append(gfex_smallr_eta_values_by_event)
                #gfex_smallr_et_values.append(gfex_smallr_et_values_by_event)
                
                """
                #jfex_it = 0
                for el in t.L1_jFexSRJetRoI:
                    #print('  JFEX: L1_jFexSRJetRoI eta = %g, phi = %g, Et = %g' %  (el.eta(), el.phi(), el.et()))
                    #jfex_smallr_phi_values_by_event.append(el.phi())
                    #jfex_smallr_eta_values_by_event.append(el.eta())
                    #jfex_smallr_et_values_by_event.append(el.et())

                    if el.et() >= 0:
                        # Digitize values
                        phi_bin = digitize(el.phi(), phi_bit_length, phi_min, phi_max)
                        eta_bin = digitize(el.eta(), eta_bit_length, eta_min, eta_max)
                        et_bin = digitize(el.et() / (et_granularity * 1000), et_bit_length, et_min / et_granularity, et_max / et_granularity)
                        #print("jfex phi_bin, eta_bin, et_bin:", phi_bin, eta_bin, et_bin)

                        # Create binary word
                        binary_word = f"{et_bin:0{et_bit_length}b}|{eta_bin:0{eta_bit_length}b}|{phi_bin:0{phi_bit_length}b}"

                        # Create hex word
                        hex_word = f"{(et_bin << (eta_bit_length + phi_bit_length)) | (eta_bin << phi_bit_length) | phi_bin:08x}"

                        # Write to file
                        #if jfex_it == 0:
                        #    f_jfex.write(f"Event : {iEvt}\n")
                        #f_jfex.write(f"0x{jfex_it:02x} {binary_word} 0x{hex_word}\n")
                    #if iEvt < 25:
                    #    jfex_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[el.et()/1000])[0]
                    #    jfex_log_heatmap += np.histogram2d([el.eta()], [el.phi()], bins=[eta_bins, phi_bins], weights=[np.log10(el.et() / 1000)])[0]
                    jfex_it += 1
                
                #jfex_smallr_phi_values.append(jfex_smallr_phi_values_by_event)
                #jfex_smallr_eta_values_by_event.append(jfex_smallr_eta_values_by_event)
                #jfex_smallr_et_values_by_event.append(jfex_smallr_et_values_by_event)
                """
                # Print the sum of transverse energy for the event
                #print(f"  Event {iEvt}: Sum of transverse energy (ΣEt) = {sum_et/1000:.2f} GeV")
                # Save the heatmaps
                if signalBool:
                    plotsDir = "signalEventPlots"
                else: 
                    plotsDir = "backgroundEventPlots"
                if signalBool:
                    if higgs_average_pt > higgs_min_average_pt and higgsPlotCounter < 20:
                        higgsPlotCounter += 1
                        #plot_heatmap(topo_heatmap, eta_bins, phi_bins, 'TopoClusters: $E_T$ in Eta-Phi Plane', plotsDir + f'/topo_heatmap_event{iEvt}.pdf', False, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        plot_heatmap("topo422", topo_422_heatmap, eta_bins, phi_bins, 'Topo422Clusters: $E_T$ in Eta-Phi Plane', plotsDir + f'/topo_422_heatmap_event{iEvt}.pdf', False, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        plot_heatmap("gfex_smallr", gfex_heatmap, eta_bins, phi_bins, 'gFex SmallR Jets: $E_T$ in Eta-Phi Plane', plotsDir + f'/gfex_heatmap_event{iEvt}.pdf', False, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        #if jFexPlotsBool:
                        #    plot_heatmap("jfex", jfex_heatmap, eta_bins, phi_bins, 'jFex: $E_T$ in Eta-Phi Plane', plotsDir + f'/jfex_heatmap_event{iEvt}.pdf', False, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        plot_heatmap("truthjet", truthInTimeJets_heatmap, eta_bins, phi_bins, 'InTimeAntiKt4TruthJets: $E_T$ in Eta-Phi Plane', plotsDir + f'/truth_jets_heatmap_event{iEvt}.pdf', False, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        plot_heatmap("hltjet", hltJets_heatmap, eta_bins, phi_bins, 'HLT_AntiKt4EMTopoJets_subjesIS: $E_T$ in Eta-Phi Plane', plotsDir + f'/hlt_jets_heatmap_event{iEvt}.pdf', False, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        plot_heatmap("gfex_larger", gfex_largerJets_heatmap, eta_bins, phi_bins, 'L1_gFexLRJetRoI: $E_T$ in Eta-Phi Plane', plotsDir + f'/gfex_larger_jets_heatmap_event{iEvt}.pdf', False, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)

                        #plot_heatmap(topo_log_heatmap, eta_bins, phi_bins, 'TopoClusters: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/topo_log_heatmap_event{iEvt}.pdf', True, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        plot_heatmap("topo422", topo_422_log_heatmap, eta_bins, phi_bins, 'Topo422Clusters: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/topo_422_log_heatmap_event{iEvt}.pdf', True, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        plot_heatmap("gfex_smallr", gfex_log_heatmap, eta_bins, phi_bins, 'gFex SmallR Jets: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/gfex_log_heatmap_event{iEvt}.pdf', True, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        #if jFexPlotsBool:
                        #    plot_heatmap("jfex", jfex_log_heatmap, eta_bins, phi_bins, 'jFex: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/jfex_log_heatmap_event{iEvt}.pdf', True, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        plot_heatmap("truthjet", truthInTimeJets_log_heatmap, eta_bins, phi_bins, 'InTimeAntiKt4TruthJets: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/truth_jets_log_heatmap_event{iEvt}.pdf', True, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        plot_heatmap("hltjet", hltJets_log_heatmap, eta_bins, phi_bins, 'HLT_AntiKt4EMTopoJets_subjesIS: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/hlt_jets_log_heatmap_event{iEvt}.pdf', True, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                        plot_heatmap("gfex_larger", gfex_largerJets_log_heatmap, eta_bins, phi_bins, 'L1_gFexLRJetRoI: $E_T$ in Eta-Phi Plane', plotsDir + f'/gfex_larger_log_jets_heatmap_event{iEvt}.pdf', False, higgs_1_pt, higgs_2_pt, higgs_1_eta, higgs_2_eta, higgs_1_phi, higgs_2_phi, True, b1_list, b2_list, dR_1, dR_2)
                else: 
                    if iEvt < 25:
                        #plot_heatmap(topo_heatmap, eta_bins, phi_bins, 'TopoClusters: $E_T$ in Eta-Phi Plane', plotsDir + f'/topo_heatmap_event{iEvt}.pdf')
                        plot_heatmap("topo422", topo_422_heatmap, eta_bins, phi_bins, 'Topo422Clusters: $E_T$ in Eta-Phi Plane', plotsDir + f'/topo_422_heatmap_event{iEvt}.pdf')
                        plot_heatmap("gfex_smallr", gfex_heatmap, eta_bins, phi_bins, 'gFex SmallR Jets: $E_T$ in Eta-Phi Plane', plotsDir + f'/gfex_heatmap_event{iEvt}.pdf')
                        #if jFexPlotsBool:
                        #    plot_heatmap("jfex", jfex_heatmap, eta_bins, phi_bins, 'jFex: $E_T$ in Eta-Phi Plane', plotsDir + f'/jfex_heatmap_event{iEvt}.pdf')
                        #plot_heatmap(topo_log_heatmap, eta_bins, phi_bins, 'TopoClusters: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/topo_log_heatmap_event{iEvt}.pdf', True)
                        plot_heatmap("topo422", topo_422_log_heatmap, eta_bins, phi_bins, 'Topo422Clusters: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/topo_422_log_heatmap_event{iEvt}.pdf', True)
                        plot_heatmap("gfex_smallr", gfex_log_heatmap, eta_bins, phi_bins, 'gFex SmallR Jets: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/gfex_log_heatmap_event{iEvt}.pdf', True) # FIXME add two leading jets pt here
                        #if jFexPlotsBool:
                        #    plot_heatmap("jfex", jfex_log_heatmap, eta_bins, phi_bins, 'jFex: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/jfex_log_heatmap_event{iEvt}.pdf', True)
                        plot_heatmap("truthjet", truthInTimeJets_heatmap, eta_bins, phi_bins, 'InTimeAntiKt4TruthJets: $E_T$ in Eta-Phi Plane', plotsDir + f'/truth_jets_heatmap_event{iEvt}.pdf')
                        plot_heatmap("hltjet", hltJets_heatmap, eta_bins, phi_bins, 'HLT_AntiKt4EMTopoJets_subjesIS: $E_T$ in Eta-Phi Plane', plotsDir + f'/hlt_jets_heatmap_event{iEvt}.pdf', False)
                        plot_heatmap("truthjet", truthInTimeJets_log_heatmap, eta_bins, phi_bins, 'InTimeAntiKt4TruthJets: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/truth_jets_log_heatmap_event{iEvt}.pdf', True)
                        plot_heatmap("hltjet", hltJets_log_heatmap, eta_bins, phi_bins, 'HLT_AntiKt4EMTopoJets_subjesIS: $Log(E_T)$ in Eta-Phi Plane', plotsDir + f'/hlt_jets_log_heatmap_event{iEvt}.pdf', True)
                        
                    

            # Clean up the transient tree for the current file
            ROOT.xAOD.ClearTransientTrees()
            fileIt += 1
            del f
            del t
if signalBool:

    # After the loop (once you've collected all values):
    plt.hist(higgs_pt_values, bins=50, range=(0, 500), histtype='step', label=r"Higgs $p_T$ [GeV]")
    plt.xlabel(r"$p_T$ [GeV]")
    plt.ylabel('# Higgs')
    plt.title(r"Higgs Truth $p_T$ Distribution")
    plt.legend()

    # Save the histogram as an image
    plt.savefig('signalTruthPlots/higgs_pt_distribution.pdf')

    if higgsPtCut:
        plt.hist(higgs_pt_values_after_cut, bins=50, range=(0, 500), histtype='step', label=r"Higgs $p_T$ after cut [GeV]")
        plt.xlabel(r"$p_T$ [GeV]")
        plt.ylabel('# Higgs')
        plt.title(r"Higgs $p_T$ Distribution")
        plt.legend()

        # Save the histogram as an image
        plt.savefig('signalTruthPlots/higgs_pt_distribution_after_cut.pdf')
        plt.clf()
    
    # Step 3: Create the histogram of b_Et values
    plt.hist(b_Et_values, bins=50, range=(0, 500), histtype='step', label='b_Et')
    plt.xlabel(r"$E_T$ [GeV]")
    plt.ylabel('# b quarks')
    plt.title(r"$b$ Truth $E_T$ Distribution")
    plt.legend()

    # Step 4: Save the histogram as a PDF
    plt.savefig('signalTruthPlots/b_Et_distribution.pdf')

    # Step 3: Create the histogram of b_Et values
    plt.hist(b_Et_values_after_higgs_cut, bins=50, range=(0, 500), histtype='step', label=r"b $E_T$ after cut [GeV]")
    plt.xlabel(r"$E_T$ [GeV]")
    plt.ylabel('# b quarks')
    plt.title(r"$b$ Truth $E_T$ Distribution")
    plt.legend()

    # Step 4: Save the histogram as a PDF
    plt.savefig('signalTruthPlots/b_Et_distribution_after_higgs_cut.pdf')
    plt.clf()
    
    pt_bins = np.linspace(40, 500, 23).astype(float)
    average_deltaR = []
    errors = []
    bin_indices = np.digitize(higgs_pt_values_for_dr, pt_bins)
    for i in range(len(pt_bins) - 1):
        #print("i:", i)
        # Find the particles that fall into this pt bin
        mask = bin_indices == i # (higgs_pt_values_for_dr >= pt_bins[i]) & (higgs_pt_values_for_dr < pt_bins[i+1])
        mask = mask.astype(bool)
        #print("mask dtype:", type(mask))
        #print(mask)
        #print(mask.shape, deltaR_values.shape)
        deltaR_values = np.array(deltaR_values)
        # Get delta R values for these particles
        deltaR_in_bin = deltaR_values[mask]
        if len(deltaR_in_bin) > 0:
            # Compute the average delta R for the bin
            avg_deltaR = np.mean(deltaR_in_bin)

            # Compute the standard error (SE) for the bin
            std_error = np.std(deltaR_in_bin) / np.sqrt(len(deltaR_in_bin))
        else:
            avg_deltaR = 0  # If no particles in this bin, set average to 0
            std_error = 0  # Set the error to 0 if no particles are in the bin
        average_deltaR.append(avg_deltaR)
        errors.append(std_error)
    # Now, we have the average delta R for each pt bin, let's plot
    plt.figure(figsize=(8, 6))
    # Compute bin centers
    bin_centers = (pt_bins[:-1] + pt_bins[1:]) / 2
    plt.errorbar(
        bin_centers, average_deltaR, yerr=errors,
        fmt='o-', capsize=5, color='black', label='Average ΔR'
    )

    # Define and plot the function f(pT) = 250.22 / pT
    pt_fit = np.linspace(min(bin_centers), max(bin_centers), 100)  # Smooth x-axis range for the fit
    deltaR_fit = 250.22 / pt_fit  # Compute function values

    # Add solid blue line at ΔR = 1.0
    plt.axhline(y=1.0, color='blue', linestyle='--', label=r'$\Delta R = 1.0$')

    plt.plot(pt_fit, deltaR_fit, 'r--', label=r'Fit: $m_H$ / Higgs $p_T$')  # Red dashed line for fit
    plt.legend(loc='best', fontsize=12)  # Add legend with automatic placement

    # Customize plot labels
    plt.xlabel(r'Higgs $p_T$ [GeV]')
    plt.ylabel('Average ΔR (bb Higgs decay products) / 20 GeV')
    plt.title(r'Average ΔR (bb Higgs decay products) vs. Higgs $p_T$')
    plt.savefig('signalTruthPlots/deltaRbb_vs_higgs_pt.pdf')
    plt.clf()
    """
    b_gfex_deltar_list = []

    for i in range(len(gfex_smallr_et_values)): # loop through events 
        if len(gfex_smallr_et_values[i]) < 2:
            continue  # Skip if there are less than 2 values

        # Find the top two highest values and their indices
        top_two = heapq.nlargest(2, enumerate(gfex_smallr_et_values[i]), key=lambda x: x[1]) # get indices of two highest energy gfex smallr jets (currently treated as seeds)
        (j1, value1), (j2, value2) = top_two
        gfex_1_eta = gfex_smallr_eta_values[i][j1]
        gfex_1_phi = gfex_smallr_phi_values[i][j1]
        gfex_2_eta = gfex_smallr_eta_values[i][j2]
        gfex_2_phi = gfex_smallr_phi_values[i][j2]
        b_gfex_deltar_list.append(delta_R(b11_list[i][0], b11_list[i][1], gfex_1_eta, gfex_1_phi))
        b_gfex_deltar_list.append(delta_R(b21_list[i][0], b21_list[i][1], gfex_1_eta, gfex_1_phi))
        b_gfex_deltar_list.append(delta_R(b12_list[i][0], b12_list[i][1], gfex_1_eta, gfex_1_phi))
        b_gfex_deltar_list.append(delta_R(b22_list[i][0], b22_list[i][1], gfex_1_eta, gfex_1_phi))
        b_gfex_deltar_list.append(delta_R(b11_list[i][0], b11_list[i][1], gfex_2_eta, gfex_2_phi))
        b_gfex_deltar_list.append(delta_R(b21_list[i][0], b21_list[i][1], gfex_2_eta, gfex_2_phi))
        b_gfex_deltar_list.append(delta_R(b12_list[i][0], b12_list[i][1], gfex_2_eta, gfex_2_phi))
        b_gfex_deltar_list.append(delta_R(b22_list[i][0], b22_list[i][1], gfex_2_eta, gfex_2_phi))


    # After the loop (once you've collected all values):
    counts, bins = np.histogram(b_gfex_deltar_list, bins=25, range=(0, 5))
    counts = counts / np.sum(counts)  # normalize to sum = 1
    plt.plot(bins[:-1], counts, drawstyle='steps-post')
    plt.xlabel(r"ΔR (b's and 2 highest Et gFex small R Jets)")
    plt.ylabel('# b, gFex small R combinations / 0.2')
    plt.title(r"ΔR Distribution (b's and 2 highest Et gFex small R Jets)")
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/deltaRb_gFex_smallR_distribution.pdf')
    plt.clf()
    """
    """
    # After the loop (once you've collected all values):
    counts, bins = np.histogram(highestEt_b_closest_gfex_deltar_list, bins=25, range=(0, 5))
    counts = counts / np.sum(counts)  # normalize to sum = 1
    plt.plot(bins[:-1], counts, drawstyle='steps-post')
    plt.xlabel(r"ΔR")
    plt.ylabel('normalized # events / 0.2')
    plt.title(r"ΔR Distribution (highest E_T b quark, closest gFex small R Jet)")
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/deltaR_highestEt_b_closest_gFex_smallR_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    counts, bins = np.histogram(second_highestEt_b_closest_gfex_deltar_list, bins=25, range=(0, 5))
    counts = counts / np.sum(counts)  # normalize to sum = 1
    plt.plot(bins[:-1], counts, drawstyle='steps-post')
    plt.xlabel(r"ΔR (2nd highest Et b and closest gFex small R Jet)")
    plt.ylabel('normalized # events / 0.2')
    plt.title(r"ΔR Distribution (2nd highest Et b and closest gFex small R Jet)")
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/deltaR_2nd_highestEt_b_closest_gFex_smallR_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    counts, bins = np.histogram(third_highestEt_b_closest_gfex_deltar_list, bins=25, range=(0, 5))
    counts = counts / np.sum(counts)  # normalize to sum = 1
    plt.plot(bins[:-1], counts, drawstyle='steps-post')
    plt.xlabel(r"ΔR (3rd highest Et b and closest gFex small R Jet)")
    plt.ylabel('normalized # events / 0.2')
    plt.title(r"ΔR Distribution (3rd highest Et b and closest gFex small R Jet)")
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/deltaR_3rd_highestEt_b_closest_gFex_smallR_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    counts, bins = np.histogram(fourth_highestEt_b_closest_gfex_deltar_list, bins=25, range=(0, 5))
    counts = counts / np.sum(counts)  # normalize to sum = 1
    plt.plot(bins[:-1], counts, drawstyle='steps-post')
    plt.xlabel(r"ΔR (4th highest Et b and closest gFex small R Jet)")
    plt.ylabel('normalized # events / 0.2')
    plt.title(r"ΔR Distribution (4th highest E_T b quark, closest gFex small R Jet)")
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/deltaR_4th_highestEt_b_closest_gFex_smallR_distribution.pdf')
    plt.clf()
    
    # After the loop (once you've collected all values):
    counts, bins = np.histogram(highestEt_b_closest_truthjet_deltar_list, bins=25, range=(0, 5))
    counts = counts / np.sum(counts)  # normalize to sum = 1
    plt.plot(bins[:-1], counts, drawstyle='steps-post')
    plt.xlabel(r"ΔR (highest Et b and closest in time Truth Jet)")
    plt.ylabel('normalized # events / 0.2')
    plt.title(r"ΔR Distribution (highest Et b and closest in time truth Jet)")
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/deltaR_highestEt_b_closest_truth_jet_distribution.pdf')
    plt.clf()

    counts, bins = np.histogram(second_highestEt_b_closest_truthjet_deltar_list, bins=25, range=(0, 5))
    counts = counts / np.sum(counts)  # normalize to sum = 1
    plt.plot(bins[:-1], counts, drawstyle='steps-post')
    plt.xlabel(r"ΔR (2nd highest Et b and closest in time truth Jet)")
    plt.ylabel('normalized # events / 0.2')
    plt.title(r"ΔR Distribution (2nd highest Et b and closest in time truth Jet)")
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/deltaR_2nd_highestEt_b_closest_truth_jet_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    counts, bins = np.histogram(third_highestEt_b_closest_truthjet_deltar_list, bins=25, range=(0, 5))
    counts = counts / np.sum(counts)  # normalize to sum = 1
    plt.plot(bins[:-1], counts, drawstyle='steps-post')
    plt.xlabel(r"ΔR (3rd highest Et b and closest in time truth Jet)")
    plt.ylabel('normalized # events / 0.2')
    plt.title(r"ΔR Distribution (3rd highest Et b and closest in time truth Jet)")
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/deltaR_3rd_highestEt_b_closest_truth_jet_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    counts, bins = np.histogram(fourth_highestEt_b_closest_truthjet_deltar_list, bins=25, range=(0, 5))
    counts = counts / np.sum(counts)  # normalize to sum = 1
    plt.plot(bins[:-1], counts, drawstyle='steps-post')
    plt.xlabel(r"ΔR (4th highest Et b and closest in time truth Jet)")
    plt.ylabel('normalized # events / 0.2')
    plt.title(r"ΔR Distribution (4th highest Et b and closest in time truth Jet)")
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/deltaR_4th_highestEt_b_closest_truth_jet_distribution.pdf')
    plt.clf()
    """
    # After the loop (once you've collected all values):
    plt.hist(jet_pt_values, bins=50, range=(0, 100), histtype='step', label='Jet pt')
    plt.xlabel(r"$p_T$ [GeV]")
    plt.ylabel('# Jets')
    plt.title(r"Truth InTimeAntiKt4TruthJets $p_T$ Distribution")
    plt.legend()
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/truth_jet_pt_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    plt.hist(jet_pt_values, bins=50, range=(0, 400), histtype='step', label='Jet pt')
    plt.xlabel(r"$p_T$ [GeV]")
    plt.ylabel('# Jets')
    plt.yscale('log')
    plt.title(r"Truth InTimeAntiKt4TruthJets $(p_T)$ Distribution")
    plt.legend()
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/truth_jet_pt_log_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    plt.hist(hlt_jet_pt_values, bins=50, range=(0, 100), histtype='step', label='Jet pt')
    plt.xlabel(r"$p_T$ [GeV]")
    plt.ylabel('# Jets')
    plt.title(r"HLT_AntiKt4EMTopoJets_subjesIS $p_T$ Distribution")
    plt.legend()
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/hlt_jet_pt_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    plt.hist(hlt_jet_pt_values, bins=50, range=(0, 400), histtype='step', label='Jet pt')
    plt.xlabel(r"$p_T$ [GeV]")
    plt.ylabel('# Jets')
    plt.yscale('log')
    plt.title(r"HLT_AntiKt4EMTopoJets_subjesIS $(p_T)$ Distribution")
    plt.legend()
    # Save the histogram as an image
    plt.savefig('signalTruthPlots/hlt_jet_pt_log_distribution.pdf')
    plt.clf()
    plt.clf()
else:
    # After the loop (once you've collected all values):
    plt.hist(jet_pt_values, bins=50, range=(0, 100), histtype='step', label='Jet pt')
    plt.xlabel(r"$p_T$ [GeV]")
    plt.ylabel('# Jets')
    plt.title(r"Truth InTimeAntiKt4TruthJets $p_T$ Distribution")
    plt.legend()
    # Save the histogram as an image
    plt.savefig('backgroundTruthPlots/truth_jet_pt_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    plt.hist(jet_pt_values, bins=50, range=(0, 400), histtype='step', label='Jet pt')
    plt.xlabel(r"$p_T$ [GeV]")
    plt.ylabel('# Jets')
    plt.yscale('log')
    plt.title(r"Truth InTimeAntiKt4TruthJets $(p_T)$ Distribution")
    plt.legend()
    # Save the histogram as an image
    plt.savefig('backgroundTruthPlots/truth_jet_pt_log_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    plt.hist(hlt_jet_pt_values, bins=50, range=(0, 100), histtype='step', label='Jet pt')
    plt.xlabel(r"$p_T$ [GeV]")
    plt.ylabel('# Jets')
    plt.title(r"HLT_AntiKt4EMTopoJets_subjesIS $p_T$ Distribution")
    plt.legend()
    # Save the histogram as an image
    plt.savefig('backgroundTruthPlots/hlt_jet_pt_distribution.pdf')
    plt.clf()

    # After the loop (once you've collected all values):
    plt.hist(hlt_jet_pt_values, bins=50, range=(0, 400), histtype='step', label='Jet pt')
    plt.xlabel(r"$p_T$ [GeV]")
    plt.ylabel('# Jets')
    plt.yscale('log')
    plt.title(r"HLT_AntiKt4EMTopoJets_subjesIS $Log_(10)(p_T)$ Distribution")
    plt.legend()
    # Save the histogram as an image
    plt.savefig('backgroundTruthPlots/hlt_jet_pt_log_distribution.pdf')
    plt.clf()

"""
# Calculate and print min and max for phi, eta, and Et
if topo_phi_values and any(topo_phi_values):
    print("Topo Phi: min = %g, max = %g" % (
        min(min(event) for event in topo_phi_values if event), 
        max(max(event) for event in topo_phi_values if event)
    ))
else:
    print("Topo Phi: No data available.")

if topo_eta_values and any(topo_eta_values):
    print("Topo Eta: min = %g, max = %g" % (
        min(min(event) for event in topo_eta_values if event), 
        max(max(event) for event in topo_eta_values if event)
    ))
else:
    print("Topo Eta: No data available.")

if topo_et_values and any(topo_et_values):
    print("Topo Et: min = %g, max = %g" % (
        min(min(event) for event in topo_et_values if event), 
        max(max(event) for event in topo_et_values if event)
    ))
else:
    print("Topo Et: No data available.")

if topo_422_phi_values and any(topo_422_phi_values):
    print("Topo422 Phi: min = %g, max = %g" % (
        min(min(event) for event in topo_422_phi_values if event), 
        max(max(event) for event in topo_422_phi_values if event)
    ))
else:
    print("Topo422 Phi: No data available.")

if topo_422_eta_values and any(topo_422_eta_values):
    print("Topo422 Eta: min = %g, max = %g" % (
        min(min(event) for event in topo_422_eta_values if event), 
        max(max(event) for event in topo_422_eta_values if event)
    ))
else:
    print("Topo422 Eta: No data available.")

if topo_422_et_values and any(topo_422_et_values):
    print("Topo422 Et: min = %g, max = %g" % (
        min(min(event) for event in topo_422_et_values if event), 
        max(max(event) for event in topo_422_et_values if event)
    ))
else:
    print("Topo422 Et: No data available.")

if gfex_smallr_phi_values and any(gfex_smallr_phi_values):
    print("gFex Phi: min = %g, max = %g" % (
        min(min(event) for event in gfex_smallr_phi_values if event), 
        max(max(event) for event in gfex_smallr_phi_values if event)
    ))
else:
    print("gFex Phi: No data available.")

if gfex_smallr_eta_values and any(gfex_smallr_eta_values):
    print("gFex Eta: min = %g, max = %g" % (
        min(min(event) for event in gfex_smallr_eta_values if event), 
        max(max(event) for event in gfex_smallr_eta_values if event)
    ))
else:
    print("gFex Eta: No data available.")

if gfex_smallr_et_values and any(gfex_smallr_et_values):
    print("gFex Et: min = %g, max = %g" % (
        min(min(event) for event in gfex_smallr_et_values if event), 
        max(max(event) for event in gfex_smallr_et_values if event)
    ))
else:
    print("gFex Et: No data available.")

if jfex_smallr_phi_values and any(jfex_smallr_phi_values):
    print("jFex Phi: min = %g, max = %g" % (
        min(min(event) for event in jfex_smallr_phi_values if event), 
        max(max(event) for event in jfex_smallr_phi_values if event)
    ))
else:
    print("jFex Phi: No data available.")

if jfex_smallr_eta_values and any(jfex_smallr_eta_values):
    print("jFex Eta: min = %g, max = %g" % (
        min(min(event) for event in jfex_smallr_eta_values if event), 
        max(max(event) for event in jfex_smallr_eta_values if event)
    ))
else:
    print("jFex Eta: No data available.")

if jfex_smallr_et_values and any(jfex_smallr_et_values):
    print("jFex Et: min = %g, max = %g" % (
        min(min(event) for event in jfex_smallr_et_values if event), 
        max(max(event) for event in jfex_smallr_et_values if event)
    ))
else:
    print("jFex Et: No data available.")
"""
if signalBool and higgsPtCut:
    print("percent of events passing pt cut:", eventCounter/ (fileIt * events_per_file))
    print("average # intermediate higgs:", higgs_counter_average / (fileIt * events_per_file))



#L1_jFexSRJetRoI
#L1_gFexSRJetRoI