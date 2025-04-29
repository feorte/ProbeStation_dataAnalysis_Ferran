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

    //----------------------------------------------------------------------------------

    // which channel to analyse
    int indx = 1; // index of the channel to study (0 = first channel, 8 channels in total)

    //----------------------------------------------------------------------------------

    int dim; //number of different voltages tested

    for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
        // load the data for the given tree entry
        tree->GetEntry(iEntry);

        // fill vaiables for the CV graph
        x.push_back(voltage);
        y.push_back(cs->at(indx));
        yerr.push_back(cs_err->at(indx));

        dim=iEntry+1;
    }

    // which channel to analyse
    int ch = channel->at(indx); 
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

    // // create a canvas to put everything together
    // TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    // glog->Draw();
    // glog->GetXaxis()->CenterTitle();
    // glog->GetYaxis()->CenterTitle();

    // // fit two lines on the graph to get the depletion voltage
    // // First fit: left region (rising region)
    // glog->Fit("pol1", "0", "", x_log[0], x_log[4]); 
    // TF1* lfit = (TF1*)glog->GetFunction("pol1")->Clone("lfit"); // <--- CLONE HERE
    // lfit->SetRange(0,5);
    // lfit->Draw("SAME"); 

    // // Second fit: rigt region (plateau region)
    // TF1* rfit = new TF1("rfit", "pol1", x_log[dim-6], x_log[dim-1]); // Define second fit
    // glog->Fit(rfit, "R0"); // R = restrict to function range, 0 = no auto draw
    // rfit->SetRange(3.5,5.5);
    // rfit->Draw("SAME");

    // //----------------------------------------
    // // Extract parameters
    // double p0_1 = lfit->GetParameter(0); // Intercept of first fit
    // double p1_1 = lfit->GetParameter(1); // Slope of first fit

    // double p0_2 = rfit->GetParameter(0); // Intercept of second fit
    // double p1_2 = rfit->GetParameter(1); // Slope of second fit


    // //----------------------------------------
    // // Find intersection (depletion voltage V_dep)
    // // The two lines cross at: p0_1 + p1_1 * x = p0_2 + p1_2 * x
    // double V_dep = (p0_2 - p0_1) / (p1_1 - p1_2);

    // // Print it
    // std::cout << "Depletion voltage V_dep = " << V_dep << std::endl;

    // //----------------------------------------
    // // Draw a vertical line at V_dep
    // TLine* line = new TLine(V_dep, glog->GetYaxis()->GetXmin(), V_dep, glog->GetYaxis()->GetXmax());
    // line->SetLineColor(kGreen+2);
    // line->SetLineStyle(2); // dashed
    // line->SetLineWidth(2);
    // line->Draw("SAME");

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
    glog->Fit("pol1", "0", "", x_log[1], x_log[6]);
    TF1* lfit = (TF1*)glog->GetFunction("pol1")->Clone("lfit");
    lfit->SetLineColor(kRed);
    lfit->SetLineWidth(2);
    lfit->SetLineStyle(2); // dashed
    lfit->SetRange(0,5);
    lfit->Draw("SAME");

    // Second fit: right region
    TF1* rfit = new TF1("rfit", "pol0", x_log[dim-6], x_log[dim-1]);
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
    // double p1_2 = rfit->GetParameter(1);
    double p1_2 = 0; // fit to a constaqnt, so slope is 0

    // intersection point (depletion voltage V_dep)
    // The two lines cross at: p0_1 + p1_1 * x = p0_2 + p1_2 * x
    double V_dep = std::exp((p0_2 - p0_1) / (p1_1 - p1_2));
    std::cout << "Depletion voltage V_dep = " << V_dep << " V" << std::endl;

    // Draw depletion voltage line
    TLine* line = new TLine(V_dep, glog->GetYaxis()->GetXmin(), V_dep, glog->GetYaxis()->GetXmax());
    line->SetLineColor(kMagenta+2);
    line->SetLineStyle(9); // long dashed
    line->SetLineWidth(2);
    line->Draw("SAME");

    // Add legend
    TLegend* legend = new TLegend(0.18, 0.18, 0.5, 0.31);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(glog, "Data (ln C vs ln V)", "p");
    legend->AddEntry(lfit, "Rising fit", "l");
    legend->AddEntry(rfit, "Plateau fit", "l");
    legend->AddEntry(line, Form("V_{dep} = %.2f V", V_dep), "l");
    legend->Draw("SAME");

    // Optional: Add text
    // TLatex latex;
    // latex.SetTextFont(42);
    // latex.SetTextSize(0.035);
    // latex.DrawLatexNDC(0.18, 0.92, Form("Channel %d: Depletion voltage analysis", ch));

    // //--------Donnor density----------

    // std::vector<float> y_new(dim);
    // std::vector<float> y_new_err(dim);

    // // calculate 1/cs^2
    // for (int i = 0; i < dim; ++i) {
    //     y_new[i] = std::pow(y[i], -2);

    //     // Error propagation formula for x^{-2}: simga_x^{-2} = 2*sigma_x / (x^3)
    //     y_new_err[i] = 2*yerr[i] / (pow(y[i],3));
    // }

    // TGraphErrors* gnew = new TGraphErrors(dim, &x[0], y_new.data(), &xerr[0], y_new_err.data());   
    // gnew->SetName(Form("Channel %d log scale", ch));
    // gnew->SetTitle(Form("Channel %d log scale;ln V; ln C", ch));
    // gnew->SetMarkerStyle(20);
    // gnew->SetMarkerSize(0.75);
    // gnew->SetMarkerColor(kBlue);

    // // // Create canvas
    // // TCanvas* c1 = new TCanvas("c1", "CV Analysis", 800, 600);
    // // c1->SetGrid();
    // // c1->SetTicks();
    // // c1->SetLeftMargin(0.15);
    // // c1->SetBottomMargin(0.15);

    // // Draw graph
    // gnew->SetMarkerStyle(20);
    // gnew->SetMarkerSize(0.8);
    // gnew->SetMarkerColor(kBlue+2);
    // gnew->GetXaxis()->SetTitleFont(42);
    // gnew->GetYaxis()->SetTitleFont(42);
    // gnew->GetXaxis()->SetLabelFont(42);
    // gnew->GetYaxis()->SetLabelFont(42);
    // gnew->GetXaxis()->SetTitleSize(0.05);
    // gnew->GetYaxis()->SetTitleSize(0.05);
    // gnew->GetXaxis()->SetTitleOffset(1.2);
    // gnew->GetYaxis()->SetTitleOffset(1.4);
    // gnew->Draw("AP"); // Important: "AP" to redraw axis properly

    // gnew->GetXaxis()->CenterTitle();
    // gnew->GetYaxis()->CenterTitle();

    // // First fit: left region
    // gnew->Fit("pol1", "", "", x[0], x[4]);
    // // TF1* lfit = (TF1*)gnew->GetFunction("pol1")->Clone("lfit");
    // // lfit->SetLineColor(kRed);
    // // lfit->SetLineWidth(2);
    // // lfit->SetLineStyle(2); // dashed
    // // lfit->SetRange(0,5);
    // // lfit->Draw("SAME");
 

    // Save the graphs
    //g->SaveAs("CV_graph.png");
    std::unique_ptr<TFile> myFile( TFile::Open("CV_graphs.root", "RECREATE") );
    g->Write();
    c1->Write(); // log scale graph and fit in the same canvas
    glog->Write(); // justs log scale graph
    //gnew->Write();




    return 0;
}