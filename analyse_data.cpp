#include <iostream>
#include "TROOT.h"
#include <fstream>
#include <string>
#include <vector> 
#include "TGraph.h"


int analyse_data()
{

    // read the data in the stored tree
    std::unique_ptr<TFile> data( TFile::Open("stored_data.root") );
    auto tree = data->Get<TTree>("analysis");

    // Disable everything...
    tree->SetBranchStatus("*", false);
    // ...but the branches we need
    for (const auto& name : {"voltage", "channel", "cs"}) {
        tree->SetBranchStatus(name, true);
    }
    

    // Associate the branches to variables
    float voltage;
    tree->SetBranchAddress("voltage", &voltage);

    std::vector<int>* channel = nullptr;
    tree->SetBranchAddress("channel", &channel);

    std::vector<float>* cs = nullptr;
    tree->SetBranchAddress("cs", &cs);
    
    //define variables for the CV graph
    std::vector<float> x;
    std::vector<float> y;
    int dim; //number of different voltages tested

    for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
        // Load the data for the given tree entry
        tree->GetEntry(iEntry);

        //fill vaiables for the CV graph
        x.push_back(voltage);
        y.push_back(cs->at(0));

        // Now, the associated variables are set to the values of the corresponding branches in tree entry `iEntry`

        dim=iEntry+1;

    }
    // which channel to plot
    int ch = channel->at(0); 
    printf("Channel: %d\n", channel->at(0));
    
    // Create a TGraph object and fill it with the data
    TGraph *g = new TGraph(dim, &x[0], &y[0]);
    TString title = Form("Channel %d;X-Axis;Y-Axis", ch);
    g->SetTitle(title);
    g->Draw();
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();


    return 0;
}