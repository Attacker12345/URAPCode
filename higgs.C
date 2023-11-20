#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "vector"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include <TChain.h>
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "Rtypes.h"
#include "TLegend.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
using namespace RooFit ;

void higgs(){

//Declaring Variables that will be used as Ttree address holders
float y1_pt, y2_pt = 0;
float y1_eta,y2_eta = 0;
float y1_phi,y2_phi = 0;
float weight,cross = 0;
UInt_t chnum = 0;
Char_t ispassed = 0;
vector <int> channels  = {364352,508784,601521,601482,601483,601484,601523,601522,601481};
vector <double> sum_weights = {357737100,95504710000,123287100,7593420,115394,184347,419728,3590.61,279324};
float myy;

//Defining histograms to store Binned Data
TH1D *h_mSignal = new TH1D("h_Signal","Signal Data",1100,105,160);
TH1D *h_mBackground = new TH1D("h_Background","Background Data",1100,105,160);
TH1D *h_ggH = new TH1D("h_ggH","Gluon Fusion Production",100,105,1160);
TH1D *h_VBF = new TH1D("h_VBF","Vector Boson Fusion Production",1100,105,160);
TH1D *h_VH = new TH1D("h_VH","Associated Production with Vector Boson",1100,105,160);
TH1D *h_ttH = new TH1D("h_ttH","Associated Production with Top Quarks",1100,105,160);

//Defining Data UNBINNED
/*
RooAbsData::setDefaultStorageType(RooAbsData::Tree); //For Storing Unbinned Data without Memory Overflow
RooRealVar mgg("mgg","mgg",105,160);
RooDataSet data_background("Background Data", "Background Data",mgg, RooFit::WeightVar("weight_background"));
RooDataSet data_signal("Signal Data", "Signal Data",mgg, RooFit::WeightVar("weight_signal"));
*/
//Each Iteration of a covers either background or a specific signal channel
for(int a=0;a<5;++a){
TChain *tree = new TChain("output");	
if(a==0) tree->Add("h029_230114/background/mc21a.aMCPy8EG_aa_FxFx_2j_myy_90_175.MxAODDetailed.e8481_s3873_r13829_p5512.h029.*.root");
//if(a==1) tree->Add("h029_230114/signal/*.root");
if(a==1) tree->Add("h029_230114/signal/mc21a.PhPy8EG_PDF4LHC21_gg*");
if(a==2) tree->Add("h029_230114/signal/mc21a.PhPy8EG_PDF4LHC21_VBF*");
if(a==3){ 
    tree->Add("h029_230114/signal/mc21a.PhPy8EG_PDF4LHC21_W*");
    tree->Add("h029_230114/signal/mc21a.PhPy8EG_PDF4LHC21_ZH125J_Zincl_MINLO.MxAODDetailedNoSkim.e8472_s3873_r13829_p5441_h029.root");
}
if(a==4) tree->Add("h029_230114/signal/mc21a.PhPy8EG_PDF4LHC21_ttH*");

//Setting Branch Adresses 
tree->SetBranchAddress("y1_pt",&y1_pt);
tree->SetBranchAddress("y2_pt",&y2_pt);
tree->SetBranchAddress("y1_eta",&y1_eta);
tree->SetBranchAddress("y2_eta",&y2_eta);
tree->SetBranchAddress("y1_phi",&y1_phi);
tree->SetBranchAddress("y2_phi",&y2_phi);
tree->SetBranchAddress("HGamEventInfoAuxDyn.isPassed",&ispassed);
tree->SetBranchAddress("HGamEventInfoAuxDyn.weight",&weight);
tree->SetBranchAddress("HGamEventInfoAuxDyn.crossSectionBRfilterEff",&cross);
tree->SetBranchAddress("EventInfoAuxDyn.mcChannelNumber",&chnum);
TLorentzVector photon1,photon2;
float weight_normalized;
int index;
     //Iterating over all entries in the Ttree
     for(int j=0;j<tree->GetEntries();++j){
             tree->GetEntry(j);
	     //verifying if event passes the cut requirments
              if (ispassed == 1){
	           //Calculating invariant mass using TLorentzVector Class      
                   photon1.SetPtEtaPhiM(y1_pt,y1_eta,y1_phi,0);
                   photon2.SetPtEtaPhiM(y2_pt,y2_eta,y2_phi,0);
                   myy = (photon1+photon2).M()/1000;
                   index = -1;
		   //Getting proper sum_weight value given the mcChannelNumber
                   for(int k=0;k<channels.size();++k){
                         if (chnum == channels[k]){
                            index = k;
		            break;
                        }
                     }
                   if (index<0){
                         cout << "ERROR CHANNEL INDEX INVALID" << endl;
                         continue;
                     }
               weight_normalized = 31400*cross*weight/sum_weights[index];
               //Filling Signal or Sideband histogram depending on if invariant mass is within signal region
               
	       if(a==0) h_mBackground->Fill(myy,weight_normalized);
               if(a>0)h_mSignal->Fill(myy,weight_normalized);
               if(a==1)h_ggH->Fill(myy,weight_normalized);  
               if(a==2)h_VBF->Fill(myy,weight_normalized);
	       if(a==3)h_VH->Fill(myy,weight_normalized);
	       if(a==4)h_ttH->Fill(myy,weight_normalized);   
		
		
	       //Filling Data UNBINNED
	       /*
		mgg.setVal(myy);
	       if(a==0) data_background.add(mgg,weight_normalized);
	       if(a>0) data_signal.add(mgg,weight_normalized);
	      */
		}
        }
  delete tree;
  /*
  RooWorkspace *workspace = new RooWorkspace("workspace");
  TFile *roofile = nullptr;
  if (a==0){ 
        workspace->import(data_background);
        roofile = new TFile("BackgroundData.root","RECREATE");
   }   
  
  if (a==4){  
      workspace->import(data_signal);
      roofile = new TFile("SignalData.root","RECREATE");
     }

  if (a!=1 and a!=2 and a!=3){ 
   workspace->Write();
   roofile->Close();
}
   delete workspace;
   delete roofile;
*/
}
/*
TCanvas *c2 = new TCanvas("c2","Unbinned Distribution",800,600);
RooDataSet data("Data", "Data", mgg);
data.append(data_background);
data.append(data_signal);
RooPlot* frame = mgg.frame();
data.plotOn(frame,RooFit::Name("data"));
frame->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");
frame->GetYaxis()->SetTitle("Events / GeV");
frame->SetTitle("Unbinned Data Distribution");
frame->SetTitleSize(0.1, "t");
frame->Draw();
c2->SaveAs("UnbinnedDistribution.png");
*/

//Writing all of the histograms to .root files
TFile *file = new TFile("Background.root", "RECREATE");
h_mBackground->Write();
file->Close();
TFile *file1 = new TFile("Signal.root", "RECREATE");
h_mSignal->Write();
file1->Close();
TFile *file2 = new TFile("ProductionChannels.root","RECREATE");
h_ggH->Write();
h_VBF->Write();
h_VH->Write();
h_ttH->Write();
file2->Close();

}

