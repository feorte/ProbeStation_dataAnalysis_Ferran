#include <iostream>
#include <fstream>
#include <string>
#include <vector> 

#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
gROOT->SetBatch(kTRUE); // Disable graphics
gStyle->SetOptStat(0); // Disable statistics box

// used values/constants
const float e = 1.602176634e-19; // charge of an electron in C
const float eps0 = 8.8541878128e-12; // vacuum permittivity in F/m
const float epsSi = 11.7; // relative permittivity of silicon
const float eps = epsSi * eps0; // permittivity of silicon in F/m
const float A = 30.25e-6; // pad area in m^2 (1 mm^2 = 1e-6 m^2), pads are 5.5x5.5 mm^2


int analyse_data()
{
    // Load ROOT file and tree safely
    auto file = std::unique_ptr<TFile>(TFile::Open("stored_data.root"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file.\n";
        return -1;
    }

    auto tree = file->Get<TTree>("analysis");
    if (!tree) {
        std::cerr << "Error: Tree 'analysis' not found in the file.\n";
        return -1;
    }

    // Enable only required branches
    const std::vector<std::string> branches = {"voltage", "channel", "cs", "cs_err"};
    tree->SetBranchStatus("*", false);
    for (const auto& name : branches) {
        if (tree->GetBranch(name.c_str())) {
            tree->SetBranchStatus(name.c_str(), true);
        } else {
            std::cerr << "Warning: Branch '" << name << "' not found in the tree." << std::endl;
        }
    }

    // Bind branches to variables
    float voltage = 0.0f;
    std::vector<int>* channel = nullptr;
    std::vector<float>* cs = nullptr;
    std::vector<float>* cs_err = nullptr;

    tree->SetBranchAddress("voltage", &voltage);
    tree->SetBranchAddress("channel", &channel);
    tree->SetBranchAddress("cs", &cs);
    tree->SetBranchAddress("cs_err", &cs_err);
    
    // define variables for the CV graph
    std::vector<float> x ;
    std::vector<float> y ;
    std::vector<float> yerr ;

    // open file to store graphs and fits
    std::unique_ptr<TFile> myFile( TFile::Open("CV_graphs.root", "RECREATE") );

    TH1F* hVdepxch = new TH1F("hVdepxch", "Depletion Voltage per Channel;Channel;V_{dep} [V]", 7, 0, 7); // histogram to store depletion voltages per channel
    TH1F* hVdep = new TH1F("hVdep", "Depletion Voltage;V_{dep} [V];Entries", 10, 40, 60); // histogram to store depletion voltages distribution
    TH1F* hdonden = new TH1F("hdonden", "Donnor density;Donnor density [ne/cm^{3}];Entries", 10, 1e+10, 2e+11); // histogram to store depletion voltages distribution

    // canvas to store 2D map of sensor with results
    auto c0 = new TCanvas("c0", "Canvas", 800, 600);
    auto channel_name = new TH2F("channel_name","Sensor pixels;X;Y",
        16,0,16,  // X axis
        16,0,16); // Y axis

    // draw text with channel number
    for (int i = 0; i < 16; ++i) { 
        for (int j = 0; j < 16; ++j) {
            channel_name->Fill(i, j, i+16*j+1);
        }
    }

    channel_name->Draw("text");


    //------------------------------------MAIN LOOP----------------------------------------

    for (int indx = 0; indx < 7; ++indx) { // loop over all channels (0-7)
        int dim=0; //number of different voltages tested
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
        //g->Draw();
        g->GetXaxis()->CenterTitle();
        g->GetYaxis()->CenterTitle();

        // create same graph but in log scale (necessary to get depletion voltage)
        // first create new vectors to hold the log-transformed data and the propagated errors:
        std::vector<float> x_log(dim), y_log(dim);
        std::vector<float> xerr_log(dim), yerr_log(dim);

        for (int i = 0; i < dim; ++i) {
            x_log[i] = std::log(x[i]);
            y_log[i] = std::log(y[i]);

            // error propagation formula for log(x): simga_log(x) = sigma_x / x
            //xerr_log[i] = xerr[i] / x[i]; // still, al 0s
            yerr_log[i] = yerr[i] / y[i];
        }

        TGraphErrors* glog = new TGraphErrors(dim, x_log.data(), y_log.data(), &xerr[0], yerr_log.data());   
        glog->SetName(Form("Channel %d log scale", ch));
        glog->SetTitle(Form("Channel %d log scale;ln V; ln C", ch));
        glog->SetMarkerStyle(20);
        glog->SetMarkerSize(0.75);
        glog->SetMarkerColor(kBlue);

        // create canvas
        TCanvas* c1 = new TCanvas(Form("CV_dep_volt_channel_%d", ch), Form("CV_dep_volt_channel_%d", ch), 800, 600);
        c1->SetGrid();
        c1->SetTicks();
        c1->SetLeftMargin(0.15);
        c1->SetBottomMargin(0.15);

        // draw graph
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
        glog->Draw("AP"); // important: "AP" to redraw axis properly

        glog->GetXaxis()->CenterTitle();
        glog->GetYaxis()->CenterTitle();

        // first fit: left region
        glog->Fit("pol1", "Q0", "", x_log[1], x_log[6]);
        TF1* lfit = (TF1*)glog->GetFunction("pol1")->Clone("lfit");
        lfit->SetLineColor(kRed);
        lfit->SetLineWidth(2);
        lfit->SetLineStyle(2); // dashed
        lfit->SetRange(0,5);
        lfit->Draw("SAME");

        // second fit: right region
        TF1* rfit = new TF1("rfit", "pol0", x_log[dim-6], x_log[dim-1]);
        glog->Fit(rfit, "QR0");
        rfit->SetLineColor(kGreen+3);
        rfit->SetLineWidth(2);
        rfit->SetLineStyle(7); // dotted
        rfit->SetRange(3.5,5.5);
        rfit->Draw("SAME");

        // calculate intersection
        double p0_1 = lfit->GetParameter(0);
        double p1_1 = lfit->GetParameter(1);
        double p0_2 = rfit->GetParameter(0);
        // double p1_2 = rfit->GetParameter(1);
        double p1_2 = 0; // fit to a constaqnt, so slope is 0

        // intersection point (depletion voltage V_dep)
        // the two lines cross at: p0_1 + p1_1 * x = p0_2 + p1_2 * x
        double log_Vdep = (p0_2 - p0_1) / (p1_1 - p1_2);
        double V_dep = std::exp(log_Vdep); // convert back to linear scale
        std::cout << "Depletion voltage V_dep = " << V_dep << " V" << std::endl;

        // store depletion voltage in histograms
        hVdepxch->GetXaxis()->SetBinLabel(indx + 1, Form("Ch%d", ch)); // set bin label
        hVdepxch->SetBinContent(indx+1, V_dep); // channel index starts at 0
        hVdep->Fill(V_dep); // fill histogram with depletion voltage

        // draw depletion voltage line
        TLine* line = new TLine(log_Vdep, glog->GetYaxis()->GetXmin(), log_Vdep, glog->GetYaxis()->GetXmax());
        line->SetLineColor(kMagenta+2);
        line->SetLineStyle(9); // long dashed
        line->SetLineWidth(2);
        line->Draw("SAME");

        // add legend
        TLegend* legend = new TLegend(0.18, 0.18, 0.5, 0.31);
        legend->SetTextFont(42);
        legend->SetTextSize(0.03);
        legend->AddEntry(glog, "Data (ln C vs ln V)", "p");
        legend->AddEntry(lfit, "Rising fit", "l");
        legend->AddEntry(rfit, "Plateau fit", "l");
        legend->AddEntry(line, Form("V_{dep} = %.2f V", V_dep), "l");
        legend->Draw("SAME");


        //--------Donnor density----------
        std::vector<float> y_new(dim);
        std::vector<float> y_new_err(dim);

        // calculate 1/cs^2
        for (int i = 0; i < dim; ++i) {
            y_new[i] = std::pow(y[i], -2);

            // error propagation formula for x^{-2}: simga_x^{-2} = 2*sigma_x / (x^3)
            y_new_err[i] = 2*yerr[i] / (pow(y[i],3));
        }

        TGraphErrors* gnew = new TGraphErrors(dim, &x[0], y_new.data(), &xerr[0], y_new_err.data());   
        gnew->SetName(Form("Channel %d", ch));
        gnew->SetTitle(Form("Channel %d;Voltage [V]; 1/C^{2} [1/pF^{2}]", ch));
        gnew->SetMarkerStyle(20);
        gnew->SetMarkerSize(0.75);
        gnew->SetMarkerColor(kBlue);

        // create canvas
        TCanvas* c2 = new TCanvas(Form("donnor_density_channel_%d", ch), Form("donnor_density_channel_%d", ch), 800, 600);
        c2->SetGrid();
        c2->SetTicks();
        c2->SetLeftMargin(0.15);
        c2->SetBottomMargin(0.15);

        // draw graph
        gnew->SetMarkerStyle(20);
        gnew->SetMarkerSize(0.8);
        gnew->SetMarkerColor(kBlue+2);
        gnew->GetXaxis()->SetTitleFont(42);
        gnew->GetYaxis()->SetTitleFont(42);
        gnew->GetXaxis()->SetLabelFont(42);
        gnew->GetYaxis()->SetLabelFont(42);
        gnew->GetXaxis()->SetTitleSize(0.05);
        gnew->GetYaxis()->SetTitleSize(0.05);
        gnew->GetXaxis()->SetTitleOffset(1.2);
        gnew->GetYaxis()->SetTitleOffset(1.4);

        gnew->Draw("AP"); // important: "AP" to redraw axis properly
        gnew->GetXaxis()->CenterTitle();
        gnew->GetYaxis()->CenterTitle();

        // first fit: left region
        gnew->Fit("pol1", "Q0", "", x[1], x[6]);
        TF1* don_fit = (TF1*)gnew->GetFunction("pol1")->Clone("don_fit");
        don_fit->SetLineColor(kRed);
        don_fit->SetLineWidth(2);
        don_fit->SetLineStyle(2); // dashed
        don_fit->SetRange(x[0], x[8]);
        don_fit->Draw("SAME");

        //get donor density from slope of the fit
        double p1_don = don_fit->GetParameter(1); //slope of the fit in [V^{-1}pF^{-2}]
        p1_don = p1_don * std::pow(10,24); // convert from pF^{-2} to F^{-2}
        double donor_density = (2)/(e*eps*std::pow(A,2)*p1_don*std::pow(10,6)); //donor density in [number elctrons*cm^{-3}] 
        std::cout << "Donnor density = " << donor_density << " ne*cm^{-3}" << std::endl;
        hdonden->Fill(donor_density); // fill histogram with donor density

        //--------2D map------------
        

        // save the graphs
        //g->SaveAs("CV_graph.png");
        g->Write();
        c1->Write(); // log scale graph and fit in the same canvas
        glog->Write(); // justs log scale graph
        //gnew->Write();
        c2->Write();

        //clean data vectors
        x.clear();
        y.clear();
        yerr.clear();

    } // end of channel loop

    hVdepxch->Write();
    hVdep->Write(); // write histogram with depletion voltages distribution
    hdonden->Write(); // write histogram with donor density distribution
    c0->Write(); // write 2D histogram canvas



    //c0->cd();



    return 0;
}