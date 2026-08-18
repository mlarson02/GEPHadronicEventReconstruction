#ifndef CALO_SHOWER_POINTING_H
#define CALO_SHOWER_POINTING_H

// Nominal calorimeter layer geometry + the shower-pointing line fit, shared by
// caloShowerEventDisplays.C (which draws the fitted line) and
// caloShowerShapePlots.C (which histograms its dca3D). Both must use the same
// geometry and the same fit, or the distributions stop describing the displays.

#include <array>
#include <cmath>
#include <vector>

#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"

// ---------------------------------------------------------------------------
// Nominal calorimeter layer geometry (meters). Illustrative only.
// Layer grouping matches the ntupler (l0..l6):
//   l0 PreSampler | l1 EM1 | l2 EM2 | l3 EM3 | l4 Tile0 | l5 Tile1 | l6 Tile2
// ---------------------------------------------------------------------------
static const double kRbarrel[7] = { 1.45, 1.53, 1.63, 1.74, 2.45, 3.02, 3.63 };
static const double kZendcap[7] = { 3.60, 3.75, 3.87, 3.97, 4.40, 5.00, 5.60 };
static const double kEtaBarrelLimit = 1.5;

// --- fit tuning -------------------------------------------------------------
// The miss distance is amplified by r1*r6/|P6-P1| ~ 2.4 m per radian of per-layer
// angular drift, i.e. ONE 0.1 tower cell of drift is already ~240 mm of DCA3D. So
// the fit lives or dies on how stable the per-layer centroid is, and the enemy is
// the outer layers, where the jet's own Et is small and the Et-weighted centroid
// gets dragged from the shower core toward the pileup average over the cone.
//
// Two knobs to bias the fit toward the core:
//   kFitEtWeightPower  1 = Et weights (as before), 2 = Et^2 weights (core-biased).
//   kFitRassoc         > 0 uses this cone for the FIT only, independent of the
//                      jet-tower association cone used for the layer sums; <= 0
//                      keeps them the same. A tight fit cone (0.15-0.2) keeps the
//                      centroid on the core while the layer fractions still use
//                      the full cone.
static const double kFitEtWeightPower = 1.0;
static const double kFitRassoc        = -1.0;

static inline double caloShowerDR(double e1, double p1, double e2, double p2) {
    double dp = std::fabs(p1 - p2);
    while (dp > M_PI) dp = std::fabs(dp - 2*M_PI);
    double de = e1 - e2;
    return std::sqrt(de*de + dp*dp);
}

// (layer, eta, phi) -> nominal (x,y,z) in meters
static void layerXYZ(int l, double eta, double phi, double& x, double& y, double& z) {
    double r, zz;
    if (std::fabs(eta) < kEtaBarrelLimit) {
        r  = kRbarrel[l];
        zz = r * std::sinh(eta);
    } else {
        zz = (eta >= 0 ? 1.0 : -1.0) * kZendcap[l];
        double s = std::sinh(eta);
        r = (std::fabs(s) > 1e-6) ? std::fabs(kZendcap[l] / s) : kRbarrel[l];
    }
    x = r * std::cos(phi);
    y = r * std::sin(phi);
    z = zz;
}

// Et-weighted per-layer-centroid straight-line fit of a jet's calorimeter shower.
// For each layer l0..l6 we take the Et-weighted mean tower position (within Rassoc
// of the jet), then fit a 3D line (weighted principal axis) through those up-to-7
// layer centroids. Unlike the jet axis this line is NOT forced through the IP, so
// for a displaced shower it points back to an offset origin. Returns false if
// fewer than two layers have energy. Fills line anchor c[3], unit direction d[3],
// and dca3D = 3D closest approach of the fitted line to the origin (m).
//
// The direction is the eigenvector of the largest eigenvalue of the Et-weighted
// scatter matrix, i.e. total least squares: it minimizes the perpendicular
// distances to the line rather than residuals along one chosen axis, so the fit
// behaves the same for forward and central jets.
static bool jetShowerPointing(double je, double jp, double Rassoc, double etMinTower,
                              const std::vector<double>* towEta,
                              const std::vector<double>* towPhi,
                              std::vector<double>* const towEtL[7],
                              double c[3], double d[3], double& dca3D,
                              int* nLayersUsed = nullptr)
{
    if (nLayersUsed) *nLayersUsed = 0;
    if (!towEta || !towPhi) return false;
    // Fit cone: tighter than the association cone when kFitRassoc > 0.
    const double Rfit = (kFitRassoc > 0.0) ? kFitRassoc : Rassoc;
    double Lx[7]={0,0,0,0,0,0,0}, Ly[7]={0,0,0,0,0,0,0}, Lz[7]={0,0,0,0,0,0,0}, Lw[7]={0,0,0,0,0,0,0};
    for (size_t t=0; t<towEta->size(); ++t) {
        double te=(*towEta)[t], tp=(*towPhi)[t];
        if (caloShowerDR(te,tp,je,jp) >= Rfit) continue;
        for (int l=0;l<7;++l) {
            double et = towEtL[l] ? (*towEtL[l])[t] : 0.0;
            if (et < etMinTower) continue;
            double x,y,z; layerXYZ(l,te,tp,x,y,z);
            // Et^p weighting: p=2 pulls the centroid onto the core, which is what
            // keeps the per-layer drift (and hence the fake DCA3D) small.
            const double w = (kFitEtWeightPower == 1.0) ? et : std::pow(et, kFitEtWeightPower);
            Lx[l]+=w*x; Ly[l]+=w*y; Lz[l]+=w*z; Lw[l]+=w;
        }
    }
    std::vector<std::array<double,3>> P; std::vector<double> W;
    for (int l=0;l<7;++l) if (Lw[l]>0) { P.push_back({Lx[l]/Lw[l],Ly[l]/Lw[l],Lz[l]/Lw[l]}); W.push_back(Lw[l]); }
    if (nLayersUsed) *nLayersUsed = (int)P.size();
    if (P.size()<2) return false;

    double sw=0; c[0]=c[1]=c[2]=0;
    for (size_t i=0;i<P.size();++i){ sw+=W[i]; for(int k=0;k<3;++k) c[k]+=W[i]*P[i][k]; }
    for(int k=0;k<3;++k) c[k]/=sw;

    TMatrixDSym M(3);
    for(int a=0;a<3;++a) for(int b=0;b<3;++b) M(a,b)=0.0;
    for (size_t i=0;i<P.size();++i){
        double dp[3]={P[i][0]-c[0],P[i][1]-c[1],P[i][2]-c[2]};
        for(int a=0;a<3;++a) for(int b=0;b<3;++b) M(a,b)+=W[i]*dp[a]*dp[b];
    }
    TMatrixDSymEigen eig(M);
    TVectorD val=eig.GetEigenValues();
    TMatrixD vec=eig.GetEigenVectors();
    int imax=0; for(int i=1;i<3;++i) if(val[i]>val[imax]) imax=i;
    for(int k=0;k<3;++k) d[k]=vec(k,imax);
    double n=std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]); if(n>0) for(int k=0;k<3;++k) d[k]/=n;

    double cd=c[0]*d[0]+c[1]*d[1]+c[2]*d[2];
    double px=c[0]-cd*d[0], py=c[1]-cd*d[1], pz=c[2]-cd*d[2];
    dca3D=std::sqrt(px*px+py*py+pz*pz);
    return true;
}

// Truth counterpart of the fitted DCA3D, for the same jet: the LLP decays at vertex
// (vx,vy,vz) and its visible products travel on from there along the parent's
// direction (eta,phi), so the truth "shower line" is the ray v + t*u. Its closest
// approach to the IP is the displacement the fit is trying to recover:
//
//   DCA3D_truth = | v - (v.u) u |
//
// Same projection as the fit, so the two are directly comparable. Units follow the
// vertex, i.e. mm for the ntupler's decayVtx_* branches. A prompt parent (v = 0)
// gives exactly 0.
static inline double truthShowerDca3D(double vx, double vy, double vz,
                                      double eta, double phi) {
    const double ch = std::cosh(eta);
    const double u[3] = { std::cos(phi)/ch, std::sin(phi)/ch, std::tanh(eta) };
    const double vu = vx*u[0] + vy*u[1] + vz*u[2];
    const double px = vx - vu*u[0], py = vy - vu*u[1], pz = vz - vu*u[2];
    return std::sqrt(px*px + py*py + pz*pz);
}

#endif  // CALO_SHOWER_POINTING_H
