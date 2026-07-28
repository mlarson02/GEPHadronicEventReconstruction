# Calorimeter shower-shape study (displaced-jet trigger)

Tooling for exploring a **displaced-jet trigger based on calorimeter shower
shape**, using GEP (Global Event Processor) EtaSK objects plus truth BSM
particles. Signals of interest: **displaced dark photons** (`HAHM ZdZd4e`) and
**emerging jets** (`Zprime2EJs`); **QCD dijet** is used as the background/prompt
comparison.

The physics idea: a long-lived BSM particle decays at a displaced vertex, so its
calorimeter energy has a different **longitudinal (per-layer) shower profile**
than a prompt jet. We distinguish the 7 calorimeter layer groups `l0..l6`
(PreSampler, EM1-3, Tile0-2) throughout.

## Files

| File | What it does |
|------|--------------|
| [`caloShowerShapeNTupler.C`](caloShowerShapeNTupler.C) | Reads a DAOD_JETM42 (xAOD) + its matching GEPOutputReader ntuple and writes a small focused ntuple (see trees below). |
| [`caloShowerEventDisplays.C`](caloShowerEventDisplays.C) | Per-event **event displays** (r-z, x-y, 3D) of the per-layer energy deposits of matched jets, with a truth-displacement overlay. Qualitative / for talks. |
| [`caloShowerShapePlots.C`](caloShowerShapePlots.C) | **Aggregate** per-layer shower-shape distributions (signal vs dijet). Quantitative — this is where the discrimination lives. |

Condor submission (under `../../condor/`):
`submit_caloShowerShape.py`, `run_caloShowerShape_job.sh`,
`submit_all_caloShowerShape.sh`.

## Data flow

```
DAOD_JETM42 (xAOD) ─┐
                    ├─► caloShowerShapeNTupler.C ─► caloShowerShape_<tag>.root ─┬─► caloShowerEventDisplays.C ─► *.pdf
GEPOutputReader ────┘                                                           └─► caloShowerShapePlots.C     ─► *.pdf
```

The ntupler is adapted from `../HERNTupler.C` and runs the same way (interpreted
ROOT macro, xAOD types via autoloading in an AnalysisBase environment).

### Ntuple contents (`caloShowerShape_<tag>.root`)

- `jetTaggerLRJEtaSKTree` — JetTaggerLRJ **EtaSK** jets, taken straight from
  GEPOutputReader (the objects used for TrigGepPerf validation).
- `wtaConeCellsTowersEtaSKTree` — WTACone **EtaSK** cells-tower jets (Et-sorted;
  pt, eta, phi, ring Et, TobN counts).
- `gepCellsTowersEtaSKTree` — GEPCellsTowerEtaSK towers with **per-layer Et**
  (`Et`, `Eta`, `Phi`, `Et_l0..Et_l6`).
- `truthBSMTree` — full `TruthBSM` collection: `pdgId`, `status`, 4-vector, and
  **production + decay vertices** (x/y/z/t, mm) so displacement is available.
  Filter by `pdgId` offline for each signal.
- `eventInfoTree` — DAOD event/run/mcChannel + weight, plus the GEP event/run
  numbers for DAOD↔GEP alignment checks.

## How to run

Ntuple (single file, e.g. under Condor):
```bash
root -b -q 'caloShowerShapeStudy/caloShowerShapeNTupler.C(1,"displaced_dark_photon","out/","daod.root","gep.root","_000001")'
```

Event displays — two PDFs per call, one per jet collection:
```bash
root -b -q 'caloShowerShapeStudy/caloShowerEventDisplays.C("signal.root","displaced_dark_photon",false,20,"plots/")'
root -b -q 'caloShowerShapeStudy/caloShowerEventDisplays.C("dijet.root","dijet_JZ4",true,20,"plots/")'   # isDijet=true
```

Aggregate shower-shape comparison (signal vs dijet), one PDF per jet collection:
```bash
root -b -q 'caloShowerShapeStudy/caloShowerShapePlots.C("signal.root","dijet.root","plots/")'
```

Grid/Condor ntuple production for all samples (signals + dijet JZ slices):
```bash
../../condor/submit_all_caloShowerShape.sh
```
`submit_caloShowerShape.py` supports `--signal {displaced_dark_photon,emerging_jets}`
or `--background --label dijet_JZ<N>` (the `--label` keeps JZ slices from
colliding and becomes the output-file tag).

## Important caveats

- **Projective towers.** GEP towers point at the IP: every layer of a tower
  shares one `(eta,phi)`, so the calo deposits in an event display always appear
  to point back to the origin. Displacement is shown only via the **truth
  overlay** (production/decay vertex + flight path drawn off-origin) and via the
  **per-layer profile** — not by the shower "bending." The quantitative
  discriminant is `caloShowerShapePlots.C`, not the displays.
- **Nominal geometry.** Event-display layer positions (`kRbarrel`/`kZendcap` in
  `caloShowerEventDisplays.C`) are approximate ATLAS radii/z for illustration,
  not measured positions.
- **No jet→tower links** in the ntuple: towers are associated to a jet by ΔR
  within `Rassoc` (top of the display/plot macros — WTACone `0.4`, LRJ `1.0`).
  **Tune these**; the WTACone cone radius in particular is a guess.
- **Per-layer Et depends on the upstream `GepClusteringAlg` sampling-energy
  fix.** Regenerate GEP ntuples with that build, or `Et_l*` reads zero for the
  EtaSK collections.
- **BSM match selection** (signal): BSM particles with a decay vertex and
  `pt > ptMinBSM` (default 5 GeV). Refine per sample once the `TruthBSM` pdgId
  content is known.
