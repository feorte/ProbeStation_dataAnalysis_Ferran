#include <iostream>
#include "TROOT.h"
#include <fstream>
#include <string>
#include <vector> 
#include "TGraph.h"
#include "TGraphErrors.h"



int analyse_data()
{

    // read the data in the stored tree
    std::unique_ptr<TFile> data( TFile::Open("stored_data.root") );
    auto tree = data->Get<TTree>("analysis");

    // Disable everything...
    tree->SetBranchStatus("*", false);
    // ...but the branches we need
    for (const auto& name : {"voltage", "channel", "cs", "cs_err"}) {
        tree->SetBranchStatus(name, true);
    }
    

    // Associate the branches to variables
    float voltage;
    tree->SetBranchAddress("voltage", &voltage);

    std::vector<int>* channel = nullptr;
    tree->SetBranchAddress("channel", &channel);

    std::vector<float>* cs = nullptr;
    tree->SetBranchAddress("cs", &cs);

    std::vector<float>* cs_err = nullptr;
    tree->SetBranchAddress("cs_err", &cs_err);
    
    //define variables for the CV graph
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> yerr;

    int dim; //number of different voltages tested

    for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
        // Load the data for the given tree entry
        tree->GetEntry(iEntry);

        //fill vaiables for the CV graph
        x.push_back(voltage);
        y.push_back(cs->at(0));
        yerr.push_back(cs_err->at(0));

        // Now, the associated variables are set to the values of the corresponding branches in tree entry `iEntry`

        dim=iEntry+1;

    }
    // which channel to plot
    int ch = channel->at(0); 
    printf("Channel: %d\n", channel->at(0));
    
    // Create a TGraphErrors object and fill it with the data
    std::vector<float> xerr(dim, 0); // no error on x-axis
    TGraph *g = new TGraphErrors(dim, &x[0], &y[0], &xerr[0], &yerr[0]);
    g->SetName(Form("Channel %d", ch));
    TString title = Form("Channel %d;Voltage [V]; Capacitance [pF]", ch);
    g->SetTitle(title);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1);
    g->SetMarkerColor(kBlue);
    g->SetLineColor(kBlue);
    g->SetLineWidth(2);
    g->SetLineStyle(1);
    g->SetFillColor(kBlue);
    g->SetFillStyle(3001);
    g->SetFillColorAlpha(kBlue, 0.1);
 
    g->Draw();
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();


    // Save the graphs
    g->SaveAs("CV_graph.png");
    std::unique_ptr<TFile> myFile( TFile::Open("CV_graphs.root", "RECREATE") );
    //TString graph_name = Form("Channel %d", ch);
    g->Write();



    return 0;
}