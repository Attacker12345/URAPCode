//CUT FLOW EXERCISE
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <vector>


void cutflow(){
TFile* file = new TFile("photon.root");
    TString rootFileName = "photon.root";
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open ROOT file " << rootFileName << endl;
        return;
    }

    // Get the TTree from the ROOT file
TTree* tree = (TTree*)file->Get("output");
	if (!tree) {
             cerr << "Error: TTree 'output' not found in " << rootFileName << endl;
            return;
            }

vector<double> *ph_pt = nullptr;
vector<double> *ph_eta = nullptr;
vector<double> *ph_phi = nullptr;

//vector<double> *ph_pt,*ph_eta, *ph_phi;

tree->SetBranchAddress("ph_pt", &ph_pt);
tree->SetBranchAddress("ph_eta", &ph_eta);
tree->SetBranchAddress("ph_phi", &ph_phi);

TH1D *h_c = new TH1D("h_cutflow","h_cutflow",6,0,6);
int pt1,pt2;
int index1=0;
int index2=0;
double eta;
int count;
TLorentzVector photon1,photon2;
for(Long64_t i=0;i<tree->GetEntries(); ++i){
	
        tree->GetEntry(i);
        //Initial Sample Bin 
	h_c->Fill(0.5);
      
       // if (i>500) break;
	
        // Photon Selection 
	
	 pt1=0;
         pt2=0;
         count=0;
         

	if (ph_pt->size()<2) continue; 
	
        //Loop meant to find the two leading photons 
        for(size_t j=0;j<ph_pt->size();++j){
              eta = abs((*ph_eta)[j]);
               //Testing eta requirment for each photon
              if ((eta<2.37 and eta>1.52) or eta<1.37){
                   count = count + 1;
                    //Finding th eleading photon, second leading photon
	           if( (*ph_pt)[j] > pt1){
                      pt2 = pt1;
                      index2 = index1;
	              pt1 = (*ph_pt)[j];
                      index1 =j;	      
	            }  
	           else if( (*ph_pt)[j] > pt2){
	 	     pt2 = (*ph_pt)[j];
	              index2 = j;
	      }
        }
	 }
        
         //Must be at least two photons that meet eta requirments
         if (count<2) continue;

	 photon1.SetPtEtaPhiM((*ph_pt)[index1],(*ph_eta)[index1],(*ph_phi)[index1],0);
	 photon2.SetPtEtaPhiM((*ph_pt)[index2],(*ph_eta)[index2],(*ph_phi)[index2],0);
	 
	 //Must also both have pt above 25
        if (photon1.Pt()<25 or photon2.Pt()<25 ){
          continue;
	}
	
        h_c->Fill(1.5);
        
	//CUT2
	
	if (photon1.Pt()<40){ continue; }

	h_c->Fill(2.5);

        //CUT3

        if (photon2.Pt()<30){ continue; }

        h_c->Fill(3.5);

         TLorentzVector diphoton = photon1 + photon2;
	 double myy  = diphoton.M();
	
	//CUT4
	
	if(photon1.Pt()/myy<0.4 or photon2.Pt()/myy<0.3){
		continue;
	}
	h_c->Fill(4.5);


       //CUT5

       if(myy>160 or myy<105){ continue; }

	h_c->Fill(5.5);
}	
//Getting/Prinitng Cut FLow chart from histogram
double binContent,binError;
 for (int bin = 1; bin <= h_c->GetNbinsX(); ++bin) {
        binContent = h_c->GetBinContent(bin);
        binError = h_c->GetBinError(bin);
        if (bin == 1) std::cout<< "Initial Sample";
        else std::cout << "Cut " << bin ;
        std::cout<< ": Events = " << binContent << ", Error = " << binError << std::endl;
       } 


}        
