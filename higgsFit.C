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
#include "RooWorkspace.h"
#include "RooCrystalBall.h"
#include "RooFitResult.h"
using namespace RooFit ;


void higgsFit(){
RooRealVar mgg("mgg","mgg",105,160);

//Obtaining Binned Data
TFile *file = new TFile("Background.root","READ");
TFile *file1 = new TFile("Signal.root","READ");
//Getting Histograms from ROOT file
TH1D *h_mBackground = (TH1D*)file->Get("h_Background");
TH1D *h_mSignal = (TH1D*)file1->Get("h_Signal");

//Obtaining Complete Invariant Mass Histogram
TH1D *h_m = new TH1D("h_m","h_m",1100,105,160);
h_m->Add(h_mBackground);
h_m->Add(h_mSignal);
RooDataHist data("data", "data", mgg, Import(*h_m));
TH1D *h_mRebin = (TH1D*)h_m->Rebin(11,"Histogram Rebinned");
RooDataHist dataRebin("data", "data", mgg, Import(*h_mRebin));
//Obtaining Unbinned Data
/*
TFile *file = new TFile("BackgroundData.root","READ");
TFile *file1 = new TFile("SignalData.root","READ");
RooWorkspace *workspace = (RooWorkspace*)file->Get("workspace");
RooWorkspace *workspace1 = (RooWorkspace*)file1->Get("workspace");

RooDataSet *data_background = (RooDataSet*)workspace->data("Background Data");
RooDataSet *data_signal = (RooDataSet*)workspace1->data("Signal Data");

RooDataSet *data = new RooDataSet("Data", "Data", mgg);
data->append(*data_background);
data->append(*data_signal);
*/

//Defining Background PDF as Bernstein Polynomial
RooRealVar coef0("coef0", "Coefficient 0", 0.5, -100, 100);
RooRealVar coef1("coef1", "Coefficient 1", 0.5, -100, 100);
RooRealVar coef2("coef2", "Coefficient 2", 0.2, -100, 100);
RooRealVar coef3("coef3", "Coefficient 3", 0.1, -100, 100);
RooRealVar coef4("coef4", "Coefficient 4", 0.05, -100, 100);
//RooRealVar coef5("coef5", "Coefficient 5", 0.01, -100, 100);
RooArgList coefList(coef0,coef1, coef2, coef3, coef4);
RooBernstein background_pdf("bernstein", "Bernstein polynomial", mgg, coefList);

//Only fitting data to the sideband regions
mgg.setRange("low",105,120);
mgg.setRange("high",130,160);
mgg.setRange("total",105,160);
//background_pdf.fitTo(*data,RooFit::Range("low,high"),RooFit::Save());

//coef0.setConstant(true);
//coef1.setConstant(true);
//coef2.setConstant(true);
//coef3.setConstant(true);
//coef4.setConstant(true);

//Getting Parameters from the saved CrystalBallFit
TFile *params_file = new TFile("CrystalBallParameters.root","READ");
RooWorkspace *params_workspace = (RooWorkspace*)params_file->Get("ggH");
RooArgSet parameters = params_workspace->allVars();

RooRealVar* avg_mgg = dynamic_cast<RooRealVar*>(parameters.find("avg_mgg"));
RooRealVar* sigmaL = dynamic_cast<RooRealVar*>(parameters.find("sigmaL"));
RooRealVar* sigmaR = dynamic_cast<RooRealVar*>(parameters.find("sigmaR"));
RooRealVar* alphaL = dynamic_cast<RooRealVar*>(parameters.find("alphaL"));
RooRealVar* alphaR = dynamic_cast<RooRealVar*>(parameters.find("alphaR"));
RooRealVar* nL = dynamic_cast<RooRealVar*>(parameters.find("nL")); 
RooRealVar* nR = dynamic_cast<RooRealVar*>(parameters.find("nR"));	
RooCrystalBall signal_pdf("CrystalBall", "Crystal Ball Fit", mgg,*avg_mgg, *sigmaL, *sigmaR, *alphaL, *nL, *alphaR, *nR);

avg_mgg->setConstant(true);
sigmaL->setConstant(true);
sigmaR->setConstant(true);
alphaL->setConstant(true);
alphaR->setConstant(true);
nL->setConstant(true);
nR->setConstant(true);

//Fitting to an Extended PDF with Background + Signal Gaussian
float upperlim =  h_m->Integral();
RooRealVar Nb("Nb", "Number of background events", upperlim/100, 0, upperlim);
RooRealVar Ns("Ns", "Number of signal events", upperlim/100, 0, upperlim);
RooAddPdf model("model", "Extended PDF", RooArgList(signal_pdf, background_pdf), RooArgList(Ns, Nb));
RooFitResult *result = model.fitTo(data,RooFit::Range("total"),RooFit::Save());
RooExtendPdf extended_pdf("extendedPdf", "Extended PDF with normalization", model, Ns);

cout << "Ns Value: " << Ns.getVal() << " +- " << Ns.getError() << endl;
cout << "Nb Value: " << Nb.getVal() << " +- " << Nb.getError() << endl;
//Computing Residuals
double dataValue,fittedValue,residual,fittedValue1,residual1;
TH1D *h_residual = new TH1D("h_residuals","h_residuals",1100,105,160);
double scaleFactor = Nb.getVal() / background_pdf.createIntegral(mgg,RooFit::Range(""))->getVal();
float width = h_m->GetBinWidth(1);
float center;
for (int i = 0; i < h_m->GetNbinsX(); ++i) {
        center = h_m->GetBinCenter(i+1);
        mgg.setRange("bin",center-width/2,center+width/2);
        dataValue = h_m->GetBinContent(i+1);
        fittedValue = background_pdf.createIntegral(mgg,RooFit::Range("bin"))->getVal()*scaleFactor;
        residual = dataValue - fittedValue;

       // cout << dataValue << " " << fittedValue << endl;`
        h_residual->SetBinContent(i+1, residual);
        h_residual->SetBinError(i+1,h_m->GetBinError(i+1));
    }
TH1D *h_residualRebin =  (TH1D*)h_residual->Rebin(11,"residualshistogram");

RooDataHist data_residual("data", "data", mgg, Import(*h_residualRebin));

//Plotting the invariant mass histogram and curve fit for the background data

double signalScale = Ns.getVal();

TCanvas *c = new TCanvas("c","Invariant Mass Distribution",800,600);
c->Divide(1, 2,0,0);

c->cd(1);
RooPlot* frame = mgg.frame();
signal_pdf.plotOn(frame,RooFit::Name("Signal PDF"),RooFit::Normalization(signalScale,RooAbsReal::NumEvent),RooFit::LineColor(kAzure+6));
dataRebin.plotOn(frame,RooFit::Name("data"));
extended_pdf.plotOn(frame,RooFit::Name("Total PDF"));
background_pdf.plotOn(frame,RooFit::Name("Background PDF"),RooFit::LineStyle(kDotted), RooFit::LineColor(kBlack));
frame->GetXaxis()->SetTitle("");
frame->GetYaxis()->SetTitle("Events / GeV");
frame->GetYaxis()->SetTitleSize(0.05);
frame->SetTitle("Diphoton Invariant Mass Distribution");
frame->SetTitleSize(0.1, "t");
frame->Draw();
TLegend* legend = new TLegend(0.65, 0.75, 0.85, 0.85);
legend->AddEntry("data", "Data", "P");
legend->AddEntry("Total PDF", "Total PDF", "l");
legend->AddEntry("Signal PDF", "Signal PDF", "l");
legend->AddEntry("Background PDF", "Background PDF", "l");
legend->SetBorderSize(0);
legend->SetFillStyle(0);
legend->SetX1(0.65);
legend->SetY1(0.6); 
legend->SetX2(0.95); 
legend->SetY2(0.9);
legend->Draw();

c->cd(2);
RooPlot* frame1 = mgg.frame();
//data_residual.plotOn(frame1,RooFit::DataError(RooAbsData::None));
data_residual.plotOn(frame1,RooFit::Name("Residuals"));
mgg.setRange("signal1",120,130);
signal_pdf.plotOn(frame1,RooFit::Range("total"),RooFit::Normalization(signalScale,RooAbsReal::NumEvent),RooFit::Name("Signal PDF"),RooFit::LineColor(kAzure+6));
frame1->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
frame1->GetXaxis()->SetTitleSize(0.05);
frame1->GetYaxis()->SetTitle("Data - Background");
frame1->GetYaxis()->SetTitleSize(0.05);
frame1->SetTitle("");
frame1->Draw();
TLegend* legend1 = new TLegend(0.65, 0.75, 0.85, 0.85);
legend1->AddEntry("Residuals", "Residuals", "P");
legend1->AddEntry("Signal PDF", "Signal PDF", "l");
legend1->SetBorderSize(0);
legend1->SetFillStyle(0);
//legend1->SetX1(0.8); 
//legend1->SetY1(0.75); 
//legend1->SetX2(0.95); 
//legend1->SetY2(0.9);
legend1->SetX1(0.65);
legend1->SetY1(0.6);
legend1->SetX2(0.95);
legend1->SetY2(0.9);
legend1->Draw();

c->Update();
c->SaveAs("InvariantMassFit.png");

//Performing Fit for Null Hypothesis

Ns.setVal(0);
Ns.setConstant(true);
RooFitResult *resultNull = model.fitTo(data,RooFit::Range("total"),RooFit::Save());

double sig = sqrt(2*(resultNull->minNll()-result->minNll()));

cout << "Significance: " << sig << endl;



file->Close();
file1->Close();

}
