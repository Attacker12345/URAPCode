#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "vector"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "Rtypes.h"
#include "TLegend.h"
using namespace RooFit ;


void higgsFit(){

TFile *file = new TFile("invariantMassDistribution.root","READ");

//Getting Histograms from ROOT file
TH1D *h_mSidebands = (TH1D*)file->Get("h_Sidebands");
TH1D *h_mSignal = (TH1D*)file->Get("h_Signal");

//Obtaining Complete Invariant Mass Histogram
TH1D *h_m = new TH1D("h_m","h_m",100,105,160);
h_m->Add(h_mSidebands);
h_m->Add(h_mSignal);

//Performing a fit using RooBernstein to the sideband histogram, getting background fit

RooRealVar mgg("mgg","mgg",105,160);

RooRealVar coef1("coef1", "Coefficient 1", 0.5, 0, 1);
RooRealVar coef2("coef2", "Coefficient 2", 0.3, 0, 0.5);
RooRealVar coef3("coef3", "Coefficient 3", 0.1, 0, 0.5);
RooRealVar coef4("coef4", "Coefficient 4", 0.05, 0, 0.2);
RooArgList coefList(coef1, coef2, coef3, coef4);
RooBernstein background_pdf("bernstein", "Bernstein polynomial", mgg, coefList);

RooDataHist data("data", "data", mgg, Import(*h_m));

//Only fitting data to the sideband regions
mgg.setRange("low",105,120);
mgg.setRange("high",130,160);
mgg.setRange("total",105,160);
background_pdf.fitTo(data,RooFit::Range("low,high"),RooFit::Save());

//setting the polynomial coefficients so as they're not adjusted later
coef1.setConstant(true);
coef2.setConstant(true);
coef3.setConstant(true);
coef4.setConstant(true);
background_pdf.fitTo(data,RooFit::Range("total"),RooFit::Save());

//cout << coef1.getVal() << endl;
//cout << coef2.getVal() << endl;
//cout << coef3.getVal() << endl;
//cout << coef4.getVal() << endl;


//Defining the Guassian Signal PDF

RooRealVar mean("mean", "mean of Gaussian", 125, 120, 130);
RooRealVar sigma("sigma", "width of Gaussian", 1.5, 0.3, 3);
RooGaussian signal_pdf("gaussian", "Gaussian PDF", mgg, mean, sigma);

mean.setConstant(true);
sigma.setConstant(true);

//Fitting to an Extended PDF with Background + Signal Gaussian
float upperlim =  h_m->Integral();
RooRealVar Nb("Nb", "Number of background events", upperlim/100, 0, upperlim);
RooRealVar Ns("Ns", "Number of signal events", upperlim/100, 0, upperlim);
RooAddPdf model("model", "Extended PDF", RooArgList(signal_pdf, background_pdf), RooArgList(Ns, Nb));
model.fitTo(data,RooFit::Range("total"),RooFit::Save());

RooExtendPdf extended_pdf("extendedPdf", "Extended PDF with normalization", model, Ns);

//Computing Residuals

double dataValue,fittedValue,residual;
TH1D *h_residual = new TH1D("h_residuals","h_residuals",100,105,160);
double scaleFactor = h_m->Integral() / background_pdf.createIntegral(mgg,RooFit::Range(""))->getVal();
float width = h_m->GetBinWidth(1);
float center;
for (int i = 0; i < h_m->GetNbinsX(); ++i) {
        center = h_m->GetBinCenter(i+1);
        mgg.setRange("bin",center-width/2,center+width/2);
        dataValue = h_m->GetBinContent(i+1);
        fittedValue = background_pdf.createIntegral(mgg,RooFit::Range("bin"))->getVal()*scaleFactor;
        residual = dataValue - fittedValue;
       // cout << dataValue << " " << fittedValue << endl;
        //FIX ERRORS ON THIS
        h_residual->SetBinContent(i+1, residual);
        h_residual->SetBinError(i+1,h_m->GetBinError(i+1));
    }

RooDataHist data_residual("data", "data", mgg, Import(*h_residual));

//Plotting the invariant mass histogram and curve fit for the background data

double signalScale = Ns.getVal();

TCanvas *c = new TCanvas("c","Invariant Mass Distribution",800,600);
c->Divide(1, 2,0,0);

c->cd(1);
RooPlot* frame = mgg.frame();
signal_pdf.plotOn(frame,RooFit::Name("Signal PDF"),RooFit::Normalization(signalScale,RooAbsReal::NumEvent),RooFit::LineColor(kRed));
data.plotOn(frame,RooFit::Name("data"));
extended_pdf.plotOn(frame,RooFit::Name("Total PDF"));
background_pdf.plotOn(frame,RooFit::Name("Background PDF"),RooFit::LineStyle(kDotted), RooFit::LineColor(kBlack));
frame->GetXaxis()->SetTitle("");
frame->GetYaxis()->SetTitle("Events / GeV");
frame->SetTitle("Diphoton Invariant Mass Distribution");
frame->SetTitleSize(0.1, "t");
frame->Draw();
TLegend* legend = new TLegend(0.65, 0.75, 0.85, 0.85);
legend->AddEntry("data", "Data", "l");
legend->AddEntry("Total PDF", "Total PDF", "l");
legend->AddEntry("Signal PDF", "Signal PDF", "l");
legend->AddEntry("Background PDF", "Background PDF", "l");
legend->SetBorderSize(0);
legend->SetFillStyle(0);
legend->SetX1(0.65); // Adjust as needed
legend->SetY1(0.6); // Adjust as needed
legend->SetX2(0.95); // Adjust as needed
legend->SetY2(0.9);
//legend->SetTextSize(0.1);
legend->Draw();

c->cd(2);

RooPlot* frame1 = mgg.frame();
//data_residual.plotOn(frame1,RooFit::DataError(RooAbsData::None));
data_residual.plotOn(frame1);
mgg.setRange("signal1",120,130);
signal_pdf.fitTo(data_residual,RooFit::Range("signal1"),RooFit::Save());
signal_pdf.plotOn(frame1,RooFit::Range("total"));
frame1->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
frame1->GetYaxis()->SetTitle("Data - Background");
frame1->SetTitle("");
frame1->Draw();

c->Update();
c->SaveAs("InvariantMassFit.png");

}
