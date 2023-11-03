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
background_pdf.fitTo(data,RooFit::Range("low,high"),RooFit::Save());

//setting the polynomial coefficients so as they're not adjusted later

coef1.setConstant(true);
coef2.setConstant(true);
coef3.setConstant(true);
coef4.setConstant(true);

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

RooRealVar Nb("Nb", "Number of background events", 10000000, 0, 100000000000);
RooRealVar Ns("Ns", "Number of signal events", 1000000, 0, 1000000000);
RooAddPdf model("model", "Extended PDF", RooArgList(signal_pdf, background_pdf), RooArgList(Ns, Nb));
mgg.setRange("signal",120,130);
model.fitTo(data,RooFit::Range("signal"),RooFit::Save());

cout << Ns.getVal() << endl;
cout << Nb.getVal() << endl;
cout << sigma.getVal() << endl;

RooExtendPdf extended_pdf("extendedPdf", "Extended PDF with normalization", model, Ns);
//Computing Residuals

//RooRealVar x("x","Bin Centers",130,105,160);
double dataValue,fittedValue,residual;
int totalEntries = h_m->GetEntries();
TH1D *h_residual = new TH1D("h_residuals","h_residuals",100,105,160);
RooArgSet obs(mgg);
double scaleFactor = h_m->Integral() / extended_pdf.createIntegral(mgg)->getVal();

for (int i = 0; i < h_m->GetNbinsX(); ++i) {
        mgg.setVal(h_m->GetBinCenter(i+1));
        dataValue = h_m->GetBinContent(i+1);
        fittedValue = extended_pdf.getVal(RooArgSet(mgg))*scaleFactor;
// fittedValue = extended_pdf.getVal();
        residual = dataValue - fittedValue;
        cout << dataValue << " " << fittedValue << endl;
        h_residual->SetBinContent(i+1, residual);
    }

RooDataHist data_residual("data", "data", mgg, Import(*h_residual));

//Drawing Complete Invariant Mass Distribution

TCanvas *c2 = new TCanvas("c2","Residuals",800,600);
h_residual->Draw();
c2->SaveAs("ResidualsTEST.png");


//Plotting the invariant mass histogram and curve fit for the background data

TCanvas *c = new TCanvas("c","Invariant Mass Distribution",800,600);
//c->Divide(1, 2);


//c->cd(1);
RooPlot* frame = mgg.frame();
data.plotOn(frame);
extended_pdf.plotOn(frame);
frame->GetXaxis()->SetTitle("Diphoton Invariant Mass (GeV)");
frame->GetYaxis()->SetTitle("Number of Events");
frame->Draw();

//c->cd(2);

//RooPlot* frame1 = mgg.frame();
//data_residual.plotOn(frame1);
//signal_pdf.plotOn(frame1);
//frame1->Draw();

//c->Update();
c->SaveAs("InvariantMassFit.png");

}

