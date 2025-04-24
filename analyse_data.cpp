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

    // disable everything...
    tree->SetBranchStatus("*", false);
    // ...but the branches we need
    for (const auto& name : {"voltage", "channel", "cs", "cs_err"}) {
        tree->SetBranchStatus(name, true);
    }
    

    // associate the branches to variables
    float voltage;
    tree->SetBranchAddress("voltage", &voltage);

    std::vector<int>* channel = nullptr;
    tree->SetBranchAddress("channel", &channel);

    std::vector<float>* cs = nullptr;
    tree->SetBranchAddress("cs", &cs);

    std::vector<float>* cs_err = nullptr;
    tree->SetBranchAddress("cs_err", &cs_err);
    
    // define variables for the CV graph
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> yerr;

    int dim; //number of different voltages tested

    for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
        // load the data for the given tree entry
        tree->GetEntry(iEntry);

        // fill vaiables for the CV graph
        x.push_back(voltage);
        y.push_back(cs->at(0));
        yerr.push_back(cs_err->at(0));

        dim=iEntry+1;
    }

    // which channel to plot
    int ch = channel->at(0); 
    printf("Channel: %d\n", channel->at(0));
    
    // create a TGraphErrors object and fill it with the data
    std::vector<float> xerr(dim, 0); // no error on x-axis
    TGraph *g = new TGraphErrors(dim, &x[0], &y[0], &xerr[0], &yerr[0]);
    g->SetName(Form("Channel %d", ch));
    g->SetTitle(Form("Channel %d;Voltage [V]; Capacitance [pF]", ch));
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.75);
    g->SetMarkerColor(kBlue);
    g->Draw();
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();

    // create same graph but in log scale (necessary to get depletion voltage)
    // first create new vectors to hold the log-transformed data and the propagated errors:
    std::vector<float> x_log(dim), y_log(dim);
    std::vector<float> xerr_log(dim), yerr_log(dim);

    for (int i = 0; i < dim; ++i) {
        x_log[i] = std::log(x[i]);
        y_log[i] = std::log(y[i]);

        // Error propagation formula for log(x): simga_log(x) = sigma_x / x
        //xerr_log[i] = xerr[i] / x[i]; // still, al 0s
        yerr_log[i] = yerr[i] / y[i];
    }

    TGraphErrors* glog = new TGraphErrors(dim, x_log.data(), y_log.data(), &xerr[0], yerr_log.data());   
    glog->SetName(Form("Channel %d log scale", ch));
    glog->SetTitle(Form("Channel %d log scale;ln V; ln C", ch));
    glog->SetMarkerStyle(20);
    glog->SetMarkerSize(0.75);
    glog->SetMarkerColor(kBlue);
    glog->Draw();
    glog->GetXaxis()->CenterTitle();
    glog->GetYaxis()->CenterTitle();

    // Save the graphs
    //g->SaveAs("CV_graph.png");
    std::unique_ptr<TFile> myFile( TFile::Open("CV_graphs.root", "RECREATE") );
    g->Write();
    glog->Write();




    return 0;
}