#ifndef CALO_SHOWER_JZ_SLICE_WEIGHTS_H
#define CALO_SHOWER_JZ_SLICE_WEIGHTS_H

// Per-JZ-slice reweighting constants for the calo shower-shape study.
//
// The QCD dijet background is produced in ten JZ slices whose cross sections span
// thirteen orders of magnitude, so a chain of the ten per-slice ntuples is only a
// physical mixture once each event carries its slice weight:
//
//   w = mcEventWeight * sigma_JZ * filterEff_JZ * L / sumOfWeights_JZ
//
// caloShowerShapeNTupler.C stamps that weight into eventInfoTree ("eventWeights"),
// exactly as HERNTupler.C does, so the downstream macros just read it back.
//
// IMPORTANT: these tables mirror the ones at the top of ../HERNTupler.C (filter
// efficiencies + cross sections from AMI, sums of weights from getSumOfWeights.py).
// If you update one, update the other. They live in a namespace here so both can be
// loaded in the same ROOT session.

namespace JZSliceWeights {

const unsigned int nJZSlices = 10;

// In barns^-1 - 7.5*10^34 cm^-2 s^-1 * 1 s (HL-LHC 200 PU inst. lumi * 1 second)
const double reweightLuminosity = 7.5e10;

// PU140 counterpart: nominal HL-LHC 140 PU instantaneous luminosity.
const double reweightLuminosity_PU140 = 5.0e10;

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

// PU200 sum of weights, from getSumOfWeights.C over the PU200 ntuples.
const double sumOfEventWeightsByJZSlice[nJZSlices] = {100000.0,      // JZ0
                                                      9493.89,       // JZ1
                                                      40.6686,       // JZ2
                                                      0.81292,       // JZ3
                                                      0.0126565,     // JZ4
                                                      0.000982199,   // JZ5
                                                      0.000164118,   // JZ6
                                                      4.75189e-05,   // JZ7
                                                      1.42645e-05,   // JZ8
                                                      3.27642e-08};  // JZ9

// PU140 counterpart (carries HERNTupler.C's caveat: still placeholder-derived).
const double sumOfEventWeightsByJZSlice_PU140[nJZSlices] = {100000.0,                // JZ0
                                                            4692.711304682304,      // JZ1
                                                            40.668641448125456,     // JZ2
                                                            0.812919741273701,      // JZ3
                                                            0.012162432307614823,   // JZ4
                                                            0.00094489973084666,    // JZ5
                                                            0.0001659370012750544,  // JZ6
                                                            4.6604455941866296e-05, // JZ7
                                                            1.3833386023877348e-05, // JZ8
                                                            3.292773754935112e-08}; // JZ9

inline double sumOfEventWeightsForPU(int jzSlice, unsigned int pileup) {
    return (pileup == 140) ? sumOfEventWeightsByJZSlice_PU140[jzSlice]
                           : sumOfEventWeightsByJZSlice[jzSlice];
}

inline double reweightLuminosityForPU(unsigned int pileup) {
    return (pileup == 140) ? reweightLuminosity_PU140 : reweightLuminosity;
}

// Per-event histogram weight for a background event of slice jzSlice. Returns
// mcEventWeight unchanged for a slice outside [0, nJZSlices) (i.e. for signal).
inline double eventWeight(double mcEventWeight, int jzSlice, unsigned int pileup = 200) {
    if (jzSlice < 0 || jzSlice >= (int)nJZSlices) return mcEventWeight;
    return mcEventWeight * crossSectionsByJZSlice[jzSlice] * filterEffienciesByJZSlice[jzSlice]
           * reweightLuminosityForPU(pileup) / sumOfEventWeightsForPU(jzSlice, pileup);
}

}  // namespace JZSliceWeights

#endif  // CALO_SHOWER_JZ_SLICE_WEIGHTS_H
