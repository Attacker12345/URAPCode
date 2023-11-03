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
#include <TChain.h>
using namespace RooFit ;

string append;

//Declaring Variables that will be used as Ttree address holders
float y1_pt, y2_pt = 0;
float y1_eta,y2_eta = 0;
float y1_phi,y2_phi = 0;
float weight = 0;
Char_t ispassed = 0;

double myy;
//Defining two histograms, One for the Signal and one for the Sidebands, can be added together later
TH1D *h_mSignal = new TH1D("h_Signal","h_Signal",100,105,160);
TH1D *h_mSidebands = new TH1D("h_Sidebands","h_Sidebandsl",100,105,160);

void higgs(){

TChain *tree = new TChain("output");
tree->Add("h029_230114/background/mc21a.aMCPy8EG_aa_FxFx_2j_myy_90_175.MxAODDetailed.e8481_s3873_r13829_p5512.h029.*.root");
tree->Add("h029_230114/signal/*.root");


     //Setting Brnach Adresses 
tree->SetBranchAddress("y1_pt",&y1_pt);
tree->SetBranchAddress("y2_pt",&y2_pt);
tree->SetBranchAddress("y1_eta",&y1_eta);
tree->SetBranchAddress("y2_eta",&y2_eta);
tree->SetBranchAddress("y1_phi",&y1_phi);
tree->SetBranchAddress("y2_phi",&y2_phi);
tree->SetBranchAddress("HGamEventInfoAuxDyn.isPassed",&ispassed);
tree->SetBranchAddress("HGamEventInfoAuxDyn.weight",&weight);
TLorentzVector photon1,photon2;
     //Iterating over all entries in the Ttree
     for(int j=0;j<tree->GetEntries();++j){
             tree->GetEntry(j);
	     //verifying if event passes the cut requirments
              if (ispassed == 1){
	           //Calculating invariant mass using TLorentzVector Class      
                   photon1.SetPtEtaPhiM(y1_pt,y1_eta,y1_phi,0);
                   photon2.SetPtEtaPhiM(y2_pt,y2_eta,y2_phi,0);
                   myy = (photon1+photon2).M()/1000;
               //Filling Signal or Sideband histogram depending on if invariant mass is within signal region
               if(myy>120 and myy<130) h_mSignal->Fill(myy,weight);
               else h_mSidebands->Fill(myy,weight);
        }
     }
TCanvas *c1 = new TCanvas("c1","Invariant Mass Distribution (Sidebands)",800,600);
h_mSidebands->Draw();
c1->SaveAs("InvariantMassDistribution.png");

TCanvas *c2 = new TCanvas("c2","Invariant Mass Distribution (Sig)",800,600);
h_mSignal->Draw();
c2->SaveAs("InvariantMassDistributionTEST.png");


TFile *file = new TFile("invariantMassDistribution.root", "RECREATE");
h_mSidebands->Write();
h_mSignal->Write();
file->Close();
}

