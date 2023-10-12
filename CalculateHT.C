// Load the necessary libraries
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"

/*
In particle physics:

HT = Scalar sum of the transverse momenta (p_T) of selected particles.

Definition:
HT = Σ |p_T|

For an event with three jets having transverse momenta p_{T1}, p_{T2}, and p_{T3}:
HT = p_{T1} + p_{T2} + p_{T3}

Use cases:
- Used in searches for new physics.
- Indicates significant undetected energy or missing transverse energy.
- Essential for studies involving potential new particles at colliders like the LHC.
*/

/*
Follow the central to make a plot for HT_central. It is the scalar sum of all the jets that has |eta| < 2.5

1. define HT_central outside the event loop
2. define its TH1 object
3. Set the jet eta branch properly
4. in the event loop, when you loop over jet, check jet eta, only jets that have |eta| < 2.5 enter the HT_central calculation
5. make the plot and save it to png file
*/

void CalculateHT() {
    // Open the ROOT file with the data
    TFile *file = new TFile("tWH.root");

    // Get the TTree named "nominal_Loose" from the file
    TTree *tree = (TTree*)file->Get("nominal_Loose");

    // Define a variable to store the calculated HT value for each event
    double HT = 0;
    double HTCentral;

    // Define a histogram to store HT values, with 100 bins ranging from 0 to 1000 (in MeV units)
    TH1D *h_HT = new TH1D("h_HT", "HT distribution", 100, 0, 1000);
    TH1D *h_HTCentral = new TH1D("h_HTCentral", "HT Central Distribution", 100, 0 , 1000);
        // Get jet momenta for this event
    std::vector<float> *jet_pt = 0, *jet_eta = 0;
    tree->SetBranchAddress("jet_pt", &jet_pt);
    tree->SetBranchAddress("jet_eta", &jet_eta);
    // Loop over all events in the TTree
    int nentries = tree->GetEntries();
    for(int i = 0; i < nentries; i++) {
        tree->GetEntry(i);  // Load data for the ith event




        // Reset HT value for this event
        HT = 0;
	HTCentral = 0;

        // Loop over all jets and sum their momenta
        for(unsigned int j = 0; j < jet_pt->size(); j++) {
	  //	  cout << (*jet_pt)[j]/1000. << endl;
            HT += (*jet_pt)[j]/1000.;
	    if (abs((*jet_eta)[j])<2.5){
	       HTCentral += (*jet_pt)[j]/1000.;

	    }
        }

        // Fill the histogram with the calculated HT value
        h_HT->Fill(HT);
	h_HTCentral->Fill(HTCentral);
    }

    // Draw the histogram
    TCanvas *c1 = new TCanvas("c1", "HT Plot", 800, 600);
    h_HT->Draw();

    // Save the plot to a file
    c1->SaveAs("HT_plot.png");

    TCanvas *c2 = new TCanvas("c2", "HT Central Plot",800,600);
    h_HTCentral->Draw();
    c2->SaveAs("HTCentral.png");
}
