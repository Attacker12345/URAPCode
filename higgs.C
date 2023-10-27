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
using namespace RooFit ;

string append;
float y1_pt, y2_pt = 0;
float y1_eta,y2_eta = 0;
float y1_phi,y2_phi = 0;
Char_t ispassed = 0;
double myy;
TH1D *h_mSignal = new TH1D("h_Signal","h_Signal",100,105,160);
TH1D *h_mSidebands = new TH1D("h_Sidebands","h_Sidebandsl",100,105,160);

void higgs(){
for (int i=1;i<497; ++i){
    if (i<10) append = "00";
    else  if (i<100) append = "0";
    else if (i==100) append ="";
    string filename = "h029_230114/background/mc21a.aMCPy8EG_aa_FxFx_2j_myy_90_175.MxAODDetailed.e8481_s3873_r13829_p5512.h029."+append+std::to_string(i)+".root"; 
 



    TFile* file = new TFile(filename.c_str());;
    TString rootFileName = filename;
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open ROOT file " << rootFileName << endl;
        continue;
    }
 

   // Get the TTree from the ROOT file
    TTree* tree = (TTree*)file->Get("output");
       if (!tree) {
           cerr << "Error: TTree 'output' not found in " << rootFileName << endl;
            return;
       }

     
     tree->SetBranchAddress("y1_pt",&y1_pt);
     tree->SetBranchAddress("y2_pt",&y2_pt);
     tree->SetBranchAddress("y1_eta",&y1_eta);
     tree->SetBranchAddress("y2_eta",&y2_eta);
     tree->SetBranchAddress("y1_phi",&y1_phi);
     tree->SetBranchAddress("y2_phi",&y2_phi);
     tree->SetBranchAddress("HGamEventInfoAuxDyn.isPassed",&ispassed);
     TLorentzVector photon1,photon2;

     for(int j=0;j<tree->GetEntries();++j){
          
             tree->GetEntry(j);
             

              if (ispassed == 1){
     
               

               photon1.SetPtEtaPhiM(y1_pt,y1_eta,y1_phi,0);
               photon2.SetPtEtaPhiM(y2_pt,y2_eta,y2_phi,0);

               myy = (photon1+photon2).M()/1000;
              // cout << myy << endl;
               if(myy>120 and myy<130) h_mSignal->Fill(myy);
               else h_mSidebands->Fill(myy);
        }
     }
    file->Close();
    }

TCanvas *c = new TCanvas("c","Invariant Mass Distribution (Sidebands)",800,600);
h_mSidebands->Draw();
c->SaveAs("InvariantMassDistribution.png");

  
}

