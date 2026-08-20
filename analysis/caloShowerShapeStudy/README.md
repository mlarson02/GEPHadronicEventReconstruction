# Calorimeter shower-shape study (displaced-jet trigger)

Tooling for exploring a **displaced-jet trigger based on calorimeter shower
shape**, using GEP (Global Event Processor) EtaSK objects plus truth BSM
particles. Signals of interest: **displaced dark photons** (`HAHM ZdZd4e`),
**emerging jets** (`Zprime2EJs`) and **long-lived staus**
(`StauStauLLP_500_0_10ns`); **QCD dijet** is used as the background/prompt
comparison.

The physics idea: a long-lived BSM particle decays at a displaced vertex, so its
calorimeter energy has a different **longitudinal (per-layer) shower profile**
than a prompt jet. We distinguish the 7 calorimeter layer groups `l0..l6`
(PreSampler, EM1-3, Tile0-2) throughout.

**The goal this all serves:** lower the E_T threshold *for displaced objects
specifically* without paying for it in background rate. Whatever provides the
displacement handle — the regression network, a classifier, or the simple
per-layer-centroid line fit below — is only interesting if the QCD rate at the
lower threshold stays put. Page 7 of `caloShowerShapePlots.C` is that plot: QCD
rate (Hz) above a leading-jet E_T threshold, with and without the fitted-`dca3D`
requirement, beside the signal efficiency it costs.

## Files

| File | What it does |
|------|--------------|
| [`caloShowerShapeNTupler.C`](caloShowerShapeNTupler.C) | Reads a DAOD_JETM42 (xAOD) + its matching GEPOutputReader ntuple and writes a small focused ntuple (see trees below). |
| [`jzSliceWeights.h`](jzSliceWeights.h) | Per-JZ-slice cross sections / filter efficiencies / sums of weights, mirrored from [`../HERNTupler.C`](../HERNTupler.C) — keep the two in sync. |
| [`caloShowerPointing.h`](caloShowerPointing.h) | Nominal layer geometry + the Et-weighted per-layer-centroid line fit and its `dca3D`, shared by the display and plot macros. |
| [`caloShowerEventDisplays.C`](caloShowerEventDisplays.C) | Per-event **event displays** (r-z, x-y, 3D) of the per-layer energy deposits of matched jets, with a truth-displacement overlay. Qualitative / for talks. |
| [`caloShowerShapePlots.C`](caloShowerShapePlots.C) | **Aggregate** per-layer shower-shape distributions (signal vs dijet), plus the shower-pointing `dca3D` of the leading and subleading jet. Quantitative — this is where the discrimination lives. |

Condor submission (under `../../condor/`):
`submit_caloShowerShape.py`, `run_caloShowerShape_job.sh`,
`submit_all_caloShowerShape.sh`.

## Data flow

```
DAOD_JETM42 (xAOD) ─┐
                    ├─► caloShowerShapeNTupler.C ─► caloShowerShape_<tag>.root ─┬─► caloShowerEventDisplays.C ─► *.pdf
GEPOutputReader ────┘   (one job per file,          <tag> = signal key or       ├─► caloShowerShapePlots.C     ─► *.pdf
                         one <tag> per JZ slice)     dijet_JZ0..dijet_JZ9       └─► make_training_parquet.py   ─► *.parquet
```
Signals are one merged ntuple each; QCD dijet is ten, one per JZ slice, chained by
the consumers (see below).

The ntupler is adapted from `../HERNTupler.C` and runs the same way (interpreted
ROOT macro, xAOD types via autoloading in an AnalysisBase environment).

### Ntuple contents (`caloShowerShape_<tag>.root`)

- `jetTaggerLRJEtaSKTree` — JetTaggerLRJ **EtaSK** jets, taken straight from
  GEPOutputReader (the objects used for TrigGepPerf validation).
- `wtaConeCellsTowersEtaSKTree` — WTACone **EtaSK** cells-tower jets (Et-sorted;
  pt, eta, phi, ring Et, TobN counts).
- `gepCellsTowersEtaSKTree` — GEPCellsTowerEtaSK towers with **per-layer Et**
  (`Et`, `Eta`, `Phi`, `Et_l0..Et_l6`).
- `gepCellsTowersTree` / `gepCellsTowersSKTree` — the same towers with **no soft
  killer** and with plain SK. Both consumer macros now default to the unsuppressed
  tree (`kTowerTree`, configurable, with an automatic warned fallback to the EtaSK
  tree for older ntuples). EtaSK's dynamic O(1-2) GeV per-tower threshold is
  correct for jet finding but starves the shower measurement: a surviving tower
  spreads that E_T over up to 7 layers, so a QCD jet is left with 2-3 towers and
  1-2 lit layers — not enough for the per-layer-centroid fit, which is what
  produces the large-DCA₃D tail in the background. Measure the shower on these
  collections; keep the trigger jets EtaSK. **Requires re-running the ntuple
  production to appear.**
- `jetTaggerLRJEtaSKTree` — off by default (`kWriteJetTaggerLRJ`): the collection
  is empty in the JZ1 (v22 PU200) production and the study is WTACone-only.
- `truthBSMTree` — full `TruthBSM` collection: `pdgId`, `status`, 4-vector, and
  **production + decay vertices** (x/y/z/t, mm) so displacement is available.
  Filter by `pdgId` offline for each signal.
- `eventInfoTree` — DAOD event/run/mcChannel + weight, plus the GEP event/run
  numbers for DAOD↔GEP alignment checks, and the JZ bookkeeping (branch names
  match `HERNTupler.C`): `sampleJZSlice` (−1 for signal), `eventWeights` (the
  slice cross-section weight, see below) and `passHSTP`.

### QCD dijet background: ten JZ slices, chained

Background is produced one JZ slice at a time — `caloShowerShapeNTupler.C` takes
a `jzSlice` argument and writes `caloShowerShape_dijet_JZ<N>.root` — and the
consumers read all ten at once through a glob:

```
caloShowerShape_dijet_JZ[0-9].root  ──►  ChainSource (../chainSource.h)  ──►  one TChain per tree
```
The `[0-9]` matters: the per-job outputs (`caloShowerShape_dijet_JZ9_000510.root`)
remain in the same flat directory after the hadd, so a bare `JZ*.root` chains the
merged files **and** their ~1700 inputs — every event counted twice.

There is no combined hadd: this mirrors `metAnalysisAndRates.C` /
`largeRJetAnalysisAndRates.C`, which read the ten v4 QCD ntuples the same way
(the merged file would exceed ROOT's 100 GB `TTree::fgMaxTreeSize`).

Two consequences for anything that touches the background:

- **Weight by the slice weight.** The slices' cross sections span thirteen orders
  of magnitude, so raw event counts across a chain are meaningless. Every event
  carries `eventWeights[0] = mcEventWeight · σ_JZ · filterEff_JZ · L /
  sumOfWeights_JZ` (constants in `jzSliceWeights.h`).
  `caloShowerShapePlots.C` fills with it; `make_training_parquet.py` exposes it
  as the `jz_weight` column.
- **HSTP filter.** Background events with `passHSTP == 0` are dropped by default,
  as in the rate macros (`kApplyHSTPFilter` in the plots macro,
  `--no-hstp-filter` in the parquet script).

## How to run

Ntuple (single file, e.g. under Condor) — the last two arguments are the JZ slice
(−1 = signal) and the pileup scenario:
```bash
root -b -q 'caloShowerShapeStudy/caloShowerShapeNTupler.C(1,"displaced_dark_photon","out/","daod.root","gep.root","_000001",-1,200)'
root -b -q 'caloShowerShapeStudy/caloShowerShapeNTupler.C(0,"","out/","daod.root","gep.root","_000001",4,200)'   # -> caloShowerShape_dijet_JZ4_000001.root
```

Event displays — one PDF per jet collection (WTACone only unless `kRunLRJ`). With no
arguments all three signals **and** the QCD dijet chain are drawn; for the dijet chain the
pages are a fixed quota per slice — `kEventsPerSlice = 10` from each of
`kDisplaySlices = {1,2,3,4}`, i.e. 40 pages, with **JZ0 excluded** (the HSTP filter
removes essentially all of it, so its events carry no rate). Any input may be a glob; drawn events
are then spread over the whole chain and labeled with their JZ slice. Dijet pages
have no truth overlay (prompt by construction) and instead report the slice and the
event's rate contribution in Hz:
```bash
root -b -l -q 'caloShowerShapeStudy/caloShowerEventDisplays.C'   # signals + dijet
root -b -q 'caloShowerShapeStudy/caloShowerEventDisplays.C("signal.root","displaced_dark_photon",false,20,"plots/")'
root -b -q 'caloShowerShapeStudy/caloShowerEventDisplays.C("/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_dijet_JZ[0-9].root","dijet_JZ0to9",true,20,"plots/")'
```

Aggregate shower-shape comparison (signal vs the weighted JZ0-9 chain), one PDF
per jet collection, seven pages:

| Page | Content |
|------|---------|
| 1 | E_T fraction in each layer `l0..l6` |
| 2 | shower depth, EM fraction, n towers |
| 3 | **fitted** shower-pointing `dca3D`, leading and subleading jet |
| 4 | **truth** `dca3D` (leading/subleading), decay `Lxy`, decay `|r|` |
| 5 | truth shower-parent kinematics: `pT`, `eta`, `phi`, mass |
| 6 | fitted − truth `dca3D` residual and the fitted-vs-truth 2D correlation |
| 7 | QCD rate (Hz) vs leading-jet E_T threshold ± the `dca3D` cut, and the signal efficiency |
| 8 | `<DCA₃D>` profiled vs `|η|`, vs layers used in the fit, vs jet E_T — the diagnostics that separate geometry artefacts from shower physics |

With no arguments all three signals run against the JZ chain:
```bash
root -b -l -q 'caloShowerShapeStudy/caloShowerShapePlots.C'
root -b -q 'caloShowerShapeStudy/caloShowerShapePlots.C("signal.root","caloShowerShape_dijet_JZ[0-9].root","plots/")'
```

Grid/Condor ntuple production for all samples (signals + dijet JZ0-9):
```bash
../../condor/submit_all_caloShowerShape.sh
# then merge the per-job outputs into one file per tag (per JZ slice):
../../condor/hadd_emulator_outputs.sh --caloshowershape
```
`submit_caloShowerShape.py` supports
`--signal {displaced_dark_photon,emerging_jets,stau_stau}` or
`--background --jz <N>`. `--jz` gives the ntupler the slice weight, resolves the
`QCD_Dijet/JZ<N>` DAOD/GEP containers under
`/data/larsonma/GEPHadronicEventReconstruction/` (dijet is not part of the
CaloShowerShapeTriggers download), and tags the output `dijet_JZ<N>`. `stau_stau`
is in the same situation as dijet — its inputs are auto-resolved from
`StauStau/` under `/data/larsonma/GEPHadronicEventReconstruction/` (see
`SIGNAL_BASES`) — while the other two signals come from the
CaloShowerShapeTriggers layout.

Training parquet — signals + all ten JZ slices by default:
```bash
python make_training_parquet.py                  # add --no-background for signals only
```

## Truth DCA₃D is ~0 by construction — do not validate the fit against it

The LLP is produced at the IP and travels along its own direction **u**, so its decay
vertex is **v** = t·**u**. The line through **v** along **u** is therefore the same
line as the IP→vertex flight path, and its closest approach to the IP is zero. The
few tens of mm actually seen (e.g. truth DCA₃D = 50 mm for a decay at |r| = 3332 mm)
is just the production-vertex offset.

Consequences:

- The `truthShowerDca3D()` number, page 4's truth DCA₃D and page 6's fitted-vs-truth
  comparison are **not** a validation of the fit. What is physically displaced is the
  decay **radius** (`Lxy`, `|r|`) — which is already the regression target in
  `make_training_parquet.py` (`truth_decay_r3d_mm`).
- More fundamentally: a decay that is collinear with the parent produces energy along
  the *same* ray from the IP, so it leaves **no pointing signature at all**. Pointing
  only sees the decay opening angle, and the E_T-weighted daughter direction averages
  back toward the parent direction. This is a strong argument that the longitudinal
  profile (which layer the shower starts in — a decay at `Lxy` = 2.3 m is inside Tile,
  so the EM layers are empty) is the real handle for this topology, not DCA₃D.
- To validate a pointing fit properly the ntupler would need the **visible decay
  products**, so the truth target could be the impact parameter of their summed
  direction rather than the parent's. `TruthBSMWithDecayParticles` has them;
  `showerParentTree` currently keeps only the parents.

### Why DCA₃D has a large background tail — and the knobs for it

The miss distance is amplified by `r₁r₆/|P₆−P₁| ≈ 2.4 m` per **radian** of per-layer
angular drift, so **one 0.1 tower cell of drift ≈ 240 mm of DCA₃D**. Hundreds of mm
therefore needs only a fraction of a cell of layer-to-layer centroid movement, which
pileup and E_T fluctuation supply for free. Two geometry effects make it worse:

- the nominal geometry switches from fixed-r to fixed-z at `|η| = 1.5`, so a jet
  straddling it gets per-layer centroids from two different placement families;
- in the endcap `r = z/sinh η`, so `dr/dη ≈ −r·coth η` (≈ 1.6 m per unit η at η = 2)
  turns a small η-centroid shift into a large radial shift.

Removing the soft killer fixed the "fit had <2 layers" failures but can *increase*
the drift, because the cone then holds ~50 soft towers and in the outer layers the
E_T-weighted centroid gets dragged from the core toward the pileup average.

| Constant | Where | Effect |
|---|---|---|
| `kMaxJetAbsEta` | plots + displays | barrel-only. In the plots it is a **soft** cut: forward jets stay out of the physics distributions but still fill the page-8 `|η|` diagnostics, so one run shows both |
| `kFitEtWeightPower` | `caloShowerPointing.h` | 1 = E_T, 2 = E_T² — biases the centroid onto the core |
| `kFitRassoc` | `caloShowerPointing.h` | tighter cone for the **fit only** (try 0.15–0.2); ≤0 = same as the association cone |

### Configurations are self-labelling

Every knob above that changes a distribution is folded into the output directory
name, so runs sit side by side instead of overwriting each other:

```
plots/<sample>/<config>/caloShowerShapePlots_WTACone.pdf
plots/<sample>/<config>/caloShowerEventDisplays_WTACone_<label>.pdf

  noSK_twrEt0p00_eta1p2_etw2p00_rfitoff      <- tower collection, per-layer E_T cut,
  EtaSK_twrEt0p50_etaoff_etw1p00_rfitoff        |eta| cut, fit weight power, fit cone
  noSK_twrEt0p00_eta1p2_etw2p00_rfit0p15
```

Non-default background handling appends `_noJZwgt` / `_noHSTP`, `kRunLRJ` appends
`_withLRJ`, and an event cap appends `_capN`. The tag is also stamped on the plots
(page 1) and on every event-display page, so a stray PDF still says what produced it.
`kOutputTagOverride` (both macros) replaces the whole tag with a name of your own.

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
