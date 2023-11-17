#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "vector"
#include "RooGlobalFunc.h"
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
#ifndef RooFit_RooFit_RooCrystalBall_h
#define RooFit_RooFit_RooCrystalBall_h
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCrystalBall.h"
#endif
using namespace RooFit ;

void CrystalBall(){
	
	//Obtaining Saved Histograms for each Production Channel
        TFile *file = new TFile("ProductionChannels.root","READ");
        TH1D *h_ggH = (TH1D*)file->Get("h_ggH");
        TH1D *h_VBF = (TH1D*)file->Get("h_VBF");
	TH1D *h_VH = (TH1D*)file->Get("h_VH");
	TH1D *h_ttH = (TH1D*)file->Get("h_ttH");

	TFile *file2 = new TFile("Signal.root","READ");
	TH1D *h_signal = (TH1D*)file2->Get("h_Signal");
	
	//Defining the parameters and PDF for the double sided crystal ball fit
	RooRealVar mgg("mgg","mgg", 105, 160);
        RooRealVar avg_mgg("avg_mgg","avg_mgg", 125, 120, 130);
        RooRealVar sigmaL("sigmaL", "sigmaL",1.5, 0.01, 5);
        RooRealVar sigmaR("sigmaR", "sigmaR",1.5, 0.01, 5);
        RooRealVar alphaL("alphaL", "alphaL",10, 0.01, 20);
        RooRealVar alphaR("alphaR", "alphaR",10, 0.01, 20);
        RooRealVar nL("nL","nL",5, 0.01, 15);
        RooRealVar nR("nR","nR",5, 0.01, 25);
	RooCrystalBall model("CrystalBall", "Crystal Ball Fit", mgg,avg_mgg, sigmaL, sigmaR, alphaL, nL, alphaR, nR);
	//avg_mgg.setConstant(true);
	
	
	TCanvas *c = new TCanvas("c", "CrystallBallFit",800,600);
	TCanvas *c1 = nullptr;
	c->Divide(2, 2,0.05,0);
	
	TFile *file1 = new TFile("CrystalBallParameters.root","RECREATE");
	for(int i = 0; i<5; ++i){
		//h_temp points to different production channel for each iteration of i
                TH1D* h_temp = nullptr;
		if(i==0) h_temp = h_ggH;
		if(i==1) h_temp = h_VBF;
		if(i==2) h_temp = h_VH;
		if(i==3) h_temp = h_ttH;
		if(i==4) h_temp = h_signal;
		if (h_temp == nullptr) cout << "RuhROh" << endl;
		
		std::string name = h_temp->GetName();
		std::string title = h_temp->GetTitle();
		
		//Create DataHist from histogram, perform fit
		RooDataHist data("data","data",mgg, Import(*h_temp));
		model.fitTo(data,RooFit::Save());
		
		//Saving all of the fit parameters to a RooWorspace and into the TFile
		RooWorkspace *workspace = new RooWorkspace(name.substr(2).c_str());
                workspace->import(avg_mgg);
                workspace->import(sigmaL);
                workspace->import(sigmaR);
                workspace->import(alphaL);
                workspace->import(alphaR);
                workspace->import(nL);
                workspace->import(nR);
                workspace->Write();
		workspace->Print();
                cout << avg_mgg.getVal() << " " << sigmaL.getVal() << " " << sigmaR.getVal() << endl;
                cout << alphaL.getVal() << " " << alphaR.getVal() << " " << nL.getVal() << " " << nR.getVal() << endl;
		delete workspace;
		
		//Saving the plot of four quadrants of each production channel fit
		if(i==4){
		   c->Update();
        	   c->SaveAs("CrystalBallFit.png");
		   c1 = new TCanvas("c1","SignalFit",800,600);
		} 
		//Making the plot in quadrant i+1 for each production channel
		else c->cd(i+1);
                RooPlot* frame = mgg.frame();
                data.plotOn(frame, RooFit::Name(title.c_str()));
                model.plotOn(frame);
                frame->SetTitle(title.c_str());
		if(i==1 or i==3) frame->GetYaxis()->SetLabelOffset(-0.05); 
		frame->Draw();
		//Saving the Signal Plot to a seperate png
		if(i==4) c1->SaveAs("CrystalBallFitSignal.png");
		
	}
	
	file->Close();
        file1->Close();
        file2->Close();
}
