// digitized_dr2_map.C
// Usage:
//   root -l -q 'digitized_dr2_map.C(1.21,"digitized_dr2_map.pdf")'

#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TLine.h"
#include "TMath.h"
#include "TEllipse.h"
#include <cmath>
#include <string>
#include <numbers> // Required for std::numbers::pi

void digitized_dr2_map(double R2cut = 1.21, const char* outPdf = "digitized_dr2_map.pdf")
{
    // -----------------------------
    // Digitization formats (from your table)
    // -----------------------------
    const double etaMin = -4.9;
    const double etaMax =  4.9;
    const int    etaBits = 7;                 // 1<<8 = 256 codes
    const int    nEta    = 98;//(1 << etaBits);    // 256
    const double etaStep = 0.1;//(etaMax - etaMin) / nEta; // 10/256 = 0.0390625

    const double phiMin = -3.2;
    const double phiMax =  3.2;
    const int    phiBits = 6;                 // 1<<6 = 64 codes
    const int    nPhi    = (1 << phiBits);    // 64
    const double phiStep = (phiMax - phiMin) / nPhi; // 6.4/64 = 0.1

    // -----------------------------
    // Histogram: x = Δη, y = Δφ
    // -----------------------------
    TH2F* h = new TH2F("h_digitized_pass",
                      Form("Digitized pass map: (#Delta#eta,#Delta#phi) with #DeltaR^{2}<%.4g;#Delta#eta (digitized);#Delta#phi (digitized)",
                           R2cut),
                      nEta, etaMin, etaMax,
                      nPhi, phiMin, phiMax);

    // Fill bin content based on digitized bin centers
    for (int ix = 1; ix <= nEta; ++ix) {
        const double deta = h->GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= nPhi; ++iy) {
            const double dphi = h->GetYaxis()->GetBinCenter(iy);

            const double dr2 = deta*deta + dphi*dphi;
            const double pass = (dr2 < R2cut) ? 1.0 : 0.0;  // 1=green, 0=red
            h->SetBinContent(ix, iy, pass);
        }
    }

    // -----------------------------
    // Styling
    // -----------------------------
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(2);

    // Discrete 2-color palette (0=red, 1=green)
    Int_t pal[2] = { kRed+1, kGreen+2 };
    gStyle->SetPalette(2, pal);

    h->SetMinimum(-0.5);
    h->SetMaximum( 1.5);
    h->SetContour(2);

    // -----------------------------
    // Canvas/margins so ALL text fits
    // -----------------------------
    TCanvas* c = new TCanvas("c", "digitized dr2 map", 1200, 850);
    c->SetLeftMargin(0.11);
    c->SetRightMargin(0.16);
    c->SetBottomMargin(0.11);
    c->SetTopMargin(0.20);

    h->SetTitle(""); // draw our own title with TLatex

    h->GetXaxis()->SetTitle("#Delta#eta (digitized)");
    h->GetYaxis()->SetTitle("#Delta#phi (digitized)");
    h->GetZaxis()->SetTitle("Pass (1) / Fail (0)");

    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetZaxis()->SetTitleSize(0.040);

    h->GetXaxis()->SetLabelSize(0.038);
    h->GetYaxis()->SetLabelSize(0.038);
    h->GetZaxis()->SetLabelSize(0.035);

    // Base map
    h->Draw("COLZ");

    // -----------------------------
    // Draw outlines around EACH PASS bin (i.e. around the green bins themselves)
    // -----------------------------
    TH2F* hPassMask = (TH2F*)h->Clone("hPassMask");
    hPassMask->Reset("ICESM");
    hPassMask->SetFillStyle(0);          // no fill, just outlines
    hPassMask->SetMarkerSize(0);

    int nPassBins = 0;
    // populate only pass bins
    for (int ix = 1; ix <= nEta; ++ix) {
        for (int iy = 1; iy <= nPhi; ++iy) {
            if (h->GetBinContent(ix, iy) > 0.5) {
                nPassBins++;
                hPassMask->SetBinContent(ix, iy, 1.0);
            }
        }
    }
    std::cout << "Area of digitized jet cone: " << nPassBins * etaStep * phiStep << "\n";
    std::cout << "Area of full granularity jet circle : " << R2cut * M_PI << "\n";

    // style for bin outlines
    hPassMask->SetLineColorAlpha(kBlack, 0.75);
    hPassMask->SetLineWidth(1);          // minimum possible width

    // Draw a rectangle for every non-empty bin in hPassMask
    hPassMask->Draw("BOX SAME");

    // -----------------------------
    // Draw an R = 1.1 circle centered at origin
    // Note: axis scaling differs in x/y, so use TEllipse with different radii in x and y.
    // -----------------------------
    const double Rcircle = std::sqrt(R2cut);

    TEllipse* circle = new TEllipse(0.0, 0.0, Rcircle, Rcircle);
    circle->SetFillStyle(0);
    circle->SetLineColorAlpha(kBlack, 0.90);
    circle->SetLineWidth(2);
    circle->SetLineStyle(2); // dashed
    circle->Draw("SAME");

    // keep axes crisp
    gPad->RedrawAxis();

    // -----------------------------
    // Legend (all entries on ONE row)
    // -----------------------------
    TH1F* hPass = new TH1F("hPass","",1,0,1);
    TH1F* hFail = new TH1F("hFail","",1,0,1);

    hPass->SetFillColor(kGreen+2);
    hFail->SetFillColor(kRed+1);
    hPass->SetLineColor(kGreen+2);
    hFail->SetLineColor(kRed+1);

    // Legend slightly wider to fit 3 entries
    TLegend* leg = new TLegend(0.12, 0.80, 0.88, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.036);

    // 3 columns → everything on one row
    leg->SetNColumns(3);

    leg->AddEntry(hPass, Form("#DeltaR^{2} < %.2f (pass)", R2cut), "f");
    leg->AddEntry(hFail, Form("#DeltaR^{2} #geq %.2f (fail)", R2cut), "f");
    leg->AddEntry(circle, Form("R = %.1f circle", Rcircle), "l");

    leg->Draw();


    // -----------------------------
    // Custom title + step annotation
    // -----------------------------
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);

    lat.SetTextSize(0.050);
    lat.SetTextAlign(22);
    lat.DrawLatex(0.50, 0.965,
                  Form("Digitized pass map: (#Delta#eta,#Delta#phi) with  #DeltaR^{2} < %.2f", R2cut));

    lat.SetTextSize(0.036);
    lat.SetTextAlign(22);
    lat.DrawLatex(0.50, 0.91,
                  Form("#eta step = %.7f (256 bins),   #phi step = %.1f (64 bins)", etaStep, phiStep));

    c->SaveAs(outPdf);

    delete c;
}
