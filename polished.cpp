#include <iostream>
#include "TROOT.h"
#include <fstream>
#include <string>
#include <vector> 
#include "TGraph.h"
#include "TGraphErrors.h"



int polished()
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
    printf("Channel: %d\n", ch);
    
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

    // Create canvas
    TCanvas* c1 = new TCanvas("c1", "CV Analysis", 800, 600);
    c1->SetGrid();
    c1->SetTicks();
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.15);

    // Draw graph
    glog->SetMarkerStyle(20);
    glog->SetMarkerSize(0.8);
    glog->SetMarkerColor(kBlue+2);
    glog->GetXaxis()->SetTitleFont(42);
    glog->GetYaxis()->SetTitleFont(42);
    glog->GetXaxis()->SetLabelFont(42);
    glog->GetYaxis()->SetLabelFont(42);
    glog->GetXaxis()->SetTitleSize(0.05);
    glog->GetYaxis()->SetTitleSize(0.05);
    glog->GetXaxis()->SetTitleOffset(1.2);
    glog->GetYaxis()->SetTitleOffset(1.4);
    glog->Draw("AP"); // Important: "AP" to redraw axis properly

    glog->GetXaxis()->CenterTitle();
    glog->GetYaxis()->CenterTitle();

    // First fit: left region
    glog->Fit("pol1", "0", "", x_log[0], x_log[4]);
    TF1* lfit = (TF1*)glog->GetFunction("pol1")->Clone("lfit");
    lfit->SetLineColor(kRed);
    lfit->SetLineWidth(2);
    lfit->SetLineStyle(2); // dashed
    lfit->SetRange(0,5);
    lfit->Draw("SAME");

    // Second fit: right region
    TF1* rfit = new TF1("rfit", "pol1", x_log[dim-6], x_log[dim-1]);
    glog->Fit(rfit, "R0");
    rfit->SetLineColor(kGreen+3);
    rfit->SetLineWidth(2);
    rfit->SetLineStyle(7); // dotted
    rfit->SetRange(3.5,5.5);
    rfit->Draw("SAME");

    // Calculate intersection
    double p0_1 = lfit->GetParameter(0);
    double p1_1 = lfit->GetParameter(1);
    double p0_2 = rfit->GetParameter(0);
    double p1_2 = rfit->GetParameter(1);

    double V_dep = (p0_2 - p0_1) / (p1_1 - p1_2);
    std::cout << "Depletion voltage V_dep = " << V_dep << std::endl;

    // Draw depletion voltage line
    TLine* line = new TLine(V_dep, glog->GetYaxis()->GetXmin(), V_dep, glog->GetYaxis()->GetXmax());
    line->SetLineColor(kMagenta+2);
    line->SetLineStyle(9); // long dashed
    line->SetLineWidth(2);
    line->Draw("SAME");

    // Add legend
    TLegend* legend = new TLegend(0.18,0.75,0.5,0.88);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(glog, "Data (ln C vs ln V)", "p");
    legend->AddEntry(lfit, "Rising fit", "l");
    legend->AddEntry(rfit, "Plateau fit", "l");
    legend->AddEntry(line, Form("V_{dep} = %.2f", V_dep), "l");
    legend->Draw("SAME");

    // Optional: Add text
    TLatex latex;
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatexNDC(0.18, 0.92, Form("Channel %d: Depletion voltage analysis", ch));

    // Save everything
    //c1->SaveAs("CV_graphs.pdf"); // good quality
    //c1->SaveAs("CV_graphs.png"); // for quick view
    c1->Write();
    glog->Write();

    return 0;
}
