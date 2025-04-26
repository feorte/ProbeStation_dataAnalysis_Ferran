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
    // create a canvas to write everithing together
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    glog->Draw();
    glog->GetXaxis()->CenterTitle();
    glog->GetYaxis()->CenterTitle();

    // fit two lines on the graph to get the depletion voltage
    // First fit: left region (rising region)
    glog->Fit("pol1", "0", "", x_log[0], x_log[4]); 
    TF1* lfit = glog->GetFunction("pol1");
    lfit->Draw("SAME"); 

    // Second fit: rigt region (plateau region)
    TF1* rfit = new TF1("rfit", "pol1", x_log[dim-6], x_log[dim-1]); // Define second fit
    glog->Fit(rfit, "R0"); // R = restrict to function range, 0 = no auto draw
    rfit->SetRange(3.5,5.5);
    rfit->Draw("SAME");

    //----------------------------------------
    // Extract parameters
    double p0_1 = lfit->GetParameter(0); // Intercept of first fit
    double p1_1 = lfit->GetParameter(1); // Slope of first fit

    double p0_2 = rfit->GetParameter(0); // Intercept of second fit
    double p1_2 = rfit->GetParameter(1); // Slope of second fit

    //----------------------------------------
    // Find intersection (depletion voltage V_dep)
    // The two lines cross at: p0_1 + p1_1 * x = p0_2 + p1_2 * x
    double V_dep = (p0_2 - p0_1) / (p1_1 - p1_2);

    // Print it
    std::cout << "Depletion voltage V_dep = " << V_dep << std::endl;

    //----------------------------------------
    // Draw a vertical line at V_dep
    TLine* line = new TLine(V_dep, glog->GetYaxis()->GetXmin(), V_dep, glog->GetYaxis()->GetXmax());
    line->SetLineColor(kGreen+2);
    line->SetLineStyle(2); // dashed
    line->SetLineWidth(2);
    line->Draw("SAME");

    // Save the graphs
    //g->SaveAs("CV_graph.png");
    std::unique_ptr<TFile> myFile( TFile::Open("CV_graphs.root", "RECREATE") );
    g->Write();
    c1->Write(); // log scale graph and fit in the same canvas
    glog->Write(); // justs log scale graph




    return 0;
}