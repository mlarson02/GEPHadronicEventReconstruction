import ROOT
import glob

baseDir = "/data/larsonma/GEPHadronicEventReconstruction/GEPOutputReaderNTuples/QCD_Dijet"

with open("sumOfWeights.txt", "w") as out:
    for jz in range(10):
        pattern = f"{baseDir}/JZ{jz}/user.mlarson.GEPNtupleJETM42.QCD_Dijet_JZ{jz}.v2_TGP_JETM42_dijet_ZvvHbb_EXT0/*.root"
        files = sorted(glob.glob(pattern))

        if not files:
            print(f"JZ{jz}: no files matched {pattern}")
            out.write(f"JZ{jz} 0.0 0\n")
            continue

        tree = ROOT.TChain("ntuple")
        for file in files:
            tree.Add(file)

        this_sum = 0.0
        nEntries = tree.GetEntries()
        for i in range(nEntries):
            tree.GetEntry(i)
            this_sum += tree.weight

        out.write(f"JZ{jz} {this_sum} {len(files)}\n")
        print(f"JZ{jz}: {len(files)} files, {nEntries} entries, Sum of Weights = {this_sum}")
