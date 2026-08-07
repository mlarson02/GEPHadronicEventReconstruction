#!/usr/bin/env python3
#
# Compare the Athena TrigGepPerf JetTaggerLRJ output against the standalone
# emulation, for the trigGepPerfValidation round-trip.
#
#   Athena LRJ : athenaJetTaggerLRJTree  in  HERNTupler_trigGepPerfValidation_v<N>.root
#                (per event: Et/Eta/Phi/M/NSubjets/Psi_R/Tau_1/Tau_2/Tau_21/MassApprox
#                 as vectors, index 0 = leading, 1 = subleading; kinematics in MeV)
#   Emulation  : jetTaggerLeadingLRJsTree / jetTaggerSubleadingLRJsTree in
#                emulation_trigGepPerfValidation_v<N>.root  (per event: one value each;
#                 kinematics in GeV)
#
# Both trees are event-aligned (HERNTupler passed events == emulation events).
#
# Usage:  python3 compareValidation.py [version=3] [tol=1e-2]
#
import sys
import ROOT

RUN_DIR = "/home/larsonma/DevelopingTriggerSimulation/run/"
MEV2GEV = 1.0e-3


def as_list(entry, branch):
    try:
        return list(getattr(entry, branch))
    except Exception:
        return []


def main():
    version = sys.argv[1] if len(sys.argv) > 1 else "2"
    tol     = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0e-2

    ath_path = f"{RUN_DIR}HERNTupler_trigGepPerfValidation_v{version}.root"
    emu_path = f"{RUN_DIR}emulation_trigGepPerfValidation_v{version}.root"

    print("comparing:", ath_path)
    print("vs. :", emu_path)

    fa = ROOT.TFile.Open(ath_path)
    fe = ROOT.TFile.Open(emu_path)
    if not fa or fa.IsZombie():
        print(f"ERROR: cannot open {ath_path}"); return 1
    if not fe or fe.IsZombie():
        print(f"ERROR: cannot open {emu_path}"); return 1

    ta      = fa.Get("athenaJetTaggerLRJTree")
    te_lead = fe.Get("jetTaggerLeadingLRJsTree")
    te_subl = fe.Get("jetTaggerSubleadingLRJsTree")
    if not ta or not te_lead or not te_subl:
        print("ERROR: missing a tree "
              f"(athena={bool(ta)} lead={bool(te_lead)} subl={bool(te_subl)})"); return 1

    na, nl, ns = ta.GetEntries(), te_lead.GetEntries(), te_subl.GetEntries()
    print(f"Comparing v{version})")
    print(f"  entries: athena={na}  emu_leading={nl}  emu_subleading={ns}")
    n = min(na, nl, ns)
    if na != nl or na != ns:
        print("  WARNING: entry counts differ -> events may be misaligned; comparing first", n)

    # (athena_branch, emu_branch, athena_scale, is_int)
    quants = [
        ("Et",         "Et",                 MEV2GEV, False),
        ("Eta",        "Eta",                1.0,     False),
        ("Phi",        "Phi",                1.0,     False),
        ("NSubjets",   "SubjetMultiplicity", 1.0,     True),
        ("Psi_R",      "Psi_R",              1.0,     False),
        ("Tau_1",      "Tau_1",              1.0,     False),
        ("Tau_2",      "Tau_2",              1.0,     False),
        ("Tau_21",     "Tau_21",             1.0,     False),
        ("MassApprox", "MassApprox",         1.0,     False),
    ]

    # A jet with Et == 0 is an empty placeholder (no valid seed / no merged
    # constituents). Its eta/phi are meaningless and reported differently by the
    # two sides (Athena: 0 from a degenerate 4-vector; emulation: the seed pos),
    # so exclude such jets from the comparison. Threshold is below one Et LSB
    # (0.125 GeV), so it cleanly separates Et==0 from any real jet.
    EMPTY_ET = 0.06

    st = {q[0]: {"lm": 0, "lt": 0, "sm": 0, "st": 0, "maxd": 0.0, "ex": []} for q in quants}
    skipped = [0, 0]  # empty leading / subleading jets excluded

    for i in range(n):
        ta.GetEntry(i); te_lead.GetEntry(i); te_subl.GetEntry(i)

        # Decide which jets are empty (Et == 0 on both sides).
        ath_Et = as_list(ta, "Et")
        emu_Et = (as_list(te_lead, "Et"), as_list(te_subl, "Et"))
        empty = [False, False]
        for k in (0, 1):
            a_et = ath_Et[k] * MEV2GEV if len(ath_Et) > k else 0.0
            e_et = emu_Et[k][0] if len(emu_Et[k]) >= 1 else 0.0
            empty[k] = (a_et < EMPTY_ET) and (e_et < EMPTY_ET)
            if empty[k]:
                skipped[k] += 1

        for ab, eb, sc, isint in quants:
            av = as_list(ta, ab)        # [leading, subleading]
            el = as_list(te_lead, eb)   # [leading]
            es = as_list(te_subl, eb)   # [subleading]
            for k, (emu, key_m, key_t) in enumerate(((el, "lm", "lt"), (es, "sm", "st"))):
                if empty[k]:
                    continue            # empty jet: skip (no meaningful quantities)
                if len(av) > k and len(emu) >= 1:
                    a = av[k] * sc
                    e = emu[0]
                    d = abs(a - e)
                    st[ab][key_t] += 1
                    if d <= tol:
                        st[ab][key_m] += 1
                    elif len(st[ab]["ex"]) < 5:
                        st[ab]["ex"].append((i, "lead" if k == 0 else "subl", a, e, d))
                    if d > st[ab]["maxd"]:
                        st[ab]["maxd"] = d

    print(f"  excluded empty (Et=0) jets: leading={skipped[0]}  subleading={skipped[1]}")
    print()
    print(f"  {'quantity':<12} {'leading':>14} {'subleading':>14} {'max|diff|':>12}")
    all_ok = True
    for ab, _, _, _ in quants:
        s = st[ab]
        lead = f"{s['lm']}/{s['lt']}"
        subl = f"{s['sm']}/{s['st']}"
        ok = (s["lm"] == s["lt"]) and (s["sm"] == s["st"])
        all_ok = all_ok and ok
        flag = "" if ok else "  <-- MISMATCH"
        print(f"  {ab:<12} {lead:>14} {subl:>14} {s['maxd']:>12.4g}{flag}")

    # show a few examples for mismatching quantities
    for ab, _, _, _ in quants:
        if st[ab]["ex"]:
            print(f"\n  examples for {ab} (event, jet, athena, emulation, |diff|):")
            for ev, jet, a, e, d in st[ab]["ex"]:
                print(f"    evt {ev:>4} {jet}: athena={a:.5g}  emu={e:.5g}  |diff|={d:.5g}")

    print()
    print("  RESULT:", "ALL MATCH" if all_ok else "MISMATCHES FOUND")
    fa.Close(); fe.Close()
    return 0 if all_ok else 2


if __name__ == "__main__":
    sys.exit(main())
