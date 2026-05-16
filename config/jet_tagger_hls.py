# To run: vitis-run --mode hls --py jet_tagger_hls.py
# vitis -s jet_tagger_hls.py 
import vitis
import os
import shutil
import sys

project_name = sys.argv[1]
hls_exec = int(sys.argv[2]) if len(sys.argv) > 1 else 1  # default to 1 if not specified
print("hls exec:", hls_exec)
cwd = os.getcwd() + '/'

# Component name
comp_name = project_name

# Define the directory path
comp_dir = f'./w/{comp_name}'

# Delete the directory if it exists
if os.path.exists(comp_dir):
    shutil.rmtree(comp_dir)

# Initialize Vitis client and set workspace
client = vitis.create_client()
client.set_workspace(path='./w')

# Create new HLS component with an empty template
comp = client.create_hls_component(name=comp_name, cfg_file=['hls_config.cfg'], template='empty')

# Modify config file
cfg_path = f'./w/{comp_name}/hls_config.cfg'
cfg = client.get_config_file(path=cfg_path)

cfg.set_value(key='part', value='xcvp1802-lsvc4072-2MP-e-S') #xcvc1902-vsva2197-2MP-e-S
cfg.set_value(section='hls', key='clock', value='4.16666666666666')  # 250 MHz
cfg.set_value(section='hls', key='flow_target', value='vivado')

cfg.set_value(section='hls',key='package.output.format', value='ip_catalog')
#cfg.set_value(key='package.output.syn', value='false')

# Set top-level function
cfg.set_value(section='hls', key='syn.top', value='jet_tag')

# Add design files
cfg.set_values(section='hls', key='syn.file', values=[
    
    #'/home/larsonma/GEPHadronicEventReconstruction/algorithm/jet_tagger_top.h',
    #'/home/larsonma/GEPHadronicEventReconstruction/algorithm/jet_tagger_top.cc',
    '/home/larsonma/GEPHadronicEventReconstruction/algorithm/jet_tag.h',
    '/home/larsonma/GEPHadronicEventReconstruction/algorithm/jet_tag.cc',
    '/home/larsonma/GEPHadronicEventReconstruction/algorithm/helperFunctions.h',
    '/home/larsonma/GEPHadronicEventReconstruction/data/LUTs/deltaR2Cut.h',
    '/home/larsonma/GEPHadronicEventReconstruction/algorithm/constants.h'
])

# Add testbench and data files
cfg.set_values(section='hls', key='tb.file', values=[
    '/home/larsonma/GEPHadronicEventReconstruction/algorithm/testbench/jet_tagger_testbench.cpp',
    '/home/larsonma/GEPHadronicEventReconstruction/algorithm/fileRead.h',
    '/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/gFex/mc21_14TeV_hh_bbbb_vbf_novhh_gfex_smallrj.dat',
    '/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/gFex/mc21_14TeV_jj_JZ3_gfex_smallrj.dat',
    '/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/CaloTopo_422/mc21_14TeV_hh_bbbb_vbf_novhh_topo422.dat',
    '/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/CaloTopo_422/mc21_14TeV_jj_JZ3_topo422.dat'
])

# Manually strip unwanted settings and add cosim waveform config
cfg_path = f'./w/{comp_name}/hls_config.cfg'
"""
with open(cfg_path, 'r') as f:
    lines = f.readlines()

# Filter out specific unwanted keys
filtered_lines = [
    line for line in lines
    if not (
        line.strip().startswith('package.output.format') or
        line.strip().startswith('package.output.syn')
    )
]

# Check if [cosim] section already exists
cosim_exists = any(line.strip().lower() == '[cosim]' for line in filtered_lines)

# If no [cosim] present, append it and the settings
if not cosim_exists:
    filtered_lines.append('\n[cosim]\n')
    #filtered_lines.append('cosim.tool=xsim\n')
    #filtered_lines.append('cosim.trace_level=port\n')
    filtered_lines.append('xsim.wave_debug = all\n')
    filtered_lines.append('xsim.dump = all\n')

with open(cfg_path, 'w') as f:
    f.writelines(filtered_lines)
"""
# Run flow
print("hls exec later:", hls_exec)
comp = client.get_component(name=comp_name)
if hls_exec == 0: # hls_exec determined by function input
    #print("this should be triggered")
    comp.run(operation='C_SIMULATION')
if hls_exec == 1: # hls_exec determined by function input
    #print("why is this triggered")
    #comp.run(operation='C_SIMULATION')
    comp.run(operation='SYNTHESIS')
elif hls_exec == 2:
    comp.run(operation='C_SIMULATION')
    comp.run(operation='SYNTHESIS')
    comp.run(operation='CO_SIMULATION')
elif hls_exec == 3:
    comp.run(operation='C_SIMULATION')
    comp.run(operation='SYNTHESIS')
    comp.run(operation='CO_SIMULATION')
    comp.run(operation='IMPLEMENTATION')