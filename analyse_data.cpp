#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
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


int CVorIV(std::string filename) {
    // Check the file extension to determine if it's CV or IV
    if (filename.size() >= 6) {
        std::string ending = filename.substr(filename.size() - 7); // get last 7 chars

        if (ending == "CV.root") {
            std::cout << "File is CV type." << std::endl;
            return 1; // CV type
        } else if (ending == "IV.root") {
            std::cout << "File is IV type." << std::endl;
            return 0; // IV type
        } else {
            std::cerr << "Error: Unknown file type." << std::endl;
            return -1;
        }
    } else {
        std::cerr << "Error: Filename is too short to determine type." << std::endl;
        return -1;
    }
}

void choose_branches(TTree* tree, const std::vector<std::string>& branches) {
    // Enable only required branches
    tree->SetBranchStatus("*", false);
    for (const auto& name : branches) {
        if (tree->GetBranch(name.c_str())) {
            tree->SetBranchStatus(name.c_str(), true);
        } else {
            std::cerr << "Warning: Branch '" << name << "' not found in the tree." << std::endl;
        }
    }
}

void log_scale(int n_volt, std::vector<float> x, std::vector<float> y, std::vector<float> yerr,
               std::vector<float>& x_log, std::vector<float>& y_log, std::vector<float>& yerr_log) {
    for (int i = 0; i < n_volt; ++i) {
        x_log[i] = std::log(x[i]);
        y_log[i] = std::log(y[i]);


        // error propagation formula for log(x): simga_log(x) = sigma_x / x
        yerr_log[i] = yerr[i] / y[i];
        //xerr_log[i] = xerr[i] / x[i]; // still, al 0s
    }
}


int analyse_data()
{
    // ------------------------------------Load data from tree---------------------------------------------------

    std::string storingfile = "stored_data/stored_data_CV.root"; // replace with your file name

    // Load ROOT file
    auto file = std::unique_ptr<TFile>(TFile::Open(storingfile.c_str()));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file.\n";
        return -1;
    }

    // Load tree from file
    auto tree = file->Get<TTree>("analysis");
    if (!tree) {
        std::cerr << "Error: Tree 'analysis' not found in the file.\n";
        return -1;
    }

    // Check if the file is CV or IV type
    int isCV = CVorIV(storingfile); // 1 for CV, 0 for IV, -1 for error
    if (isCV == -1) {
        return -1;
    }

    // ------------------------------------Get necessary branches---------------------------------------------------

    if (isCV==true) {
        cout << "Processing CV data...\n";
        // Enable only required branches
        const std::vector<std::string> branches = {"voltage", "channel", "cs", "cs_err"};
        choose_branches(tree, branches);

        // Bind branches to variables
        float voltage = 0.0f;
        std::vector<int>* channel = 0;
        std::vector<float>* cs = 0;
        std::vector<float>* cs_err = 0;
        tree->SetBranchAddress("voltage", &voltage);
        tree->SetBranchAddress("channel", &channel);
        tree->SetBranchAddress("cs", &cs);
        tree->SetBranchAddress("cs_err", &cs_err);

        // ------------------------------------Create histograms---------------------------------------------------

        // 1D histograms to store depletion voltages and donor density distributions
        TH1F* hVdepxch = new TH1F("hVdepxch", "Depletion Voltage per Channel;Channel;V_{dep} [V]", 256, 0.5, 256.5); // histogram to store depletion voltages per channel
        TH1F* hVdep = new TH1F("hVdep", "Depletion Voltage;V_{dep} [V];Entries", 50, 17, 80); // histogram to store depletion voltages distribution
        TH1F* hndon = new TH1F("hndon", "Donnor density;Donnor density [ne/cm^{3}];Entries", 50, 3e+10, 4e+11); // histogram to store depletion voltages distribution

        // 2D histogram to store capacitance, depletion voltage and donor density of every channel
        auto hVdep_map = new TH2F("hVdep_map","Depletion Voltage;X;Y",
            16,0.5,16.5,  // X axis
            16,0.5,16.5); // Y axis
        auto hndon_map = new TH2F("hndon_map","Donnor density;X;Y", 16,0.5,16.5, 16,0.5,16.5);
        auto hcs_plateau_map = new TH2F("hcs_plateau_map","High voltage capacitance;X;Y", 16,0.5,16.5, 16,0.5,16.5);


        //2D histogram that maps channels to positions on the sensor
        auto channel_name = new TH2F("channel_name","Sensor pixels;X;Y", 16,0.5,16.5, 16,0.5,16.5);
        // fill it with channel numbers
        for (int i = 1; i < 17; ++i) {
            for (int j = 1; j < 17; ++j) {
                channel_name->Fill(i, j, i+16*(j-1));
            }
        }

        // --------------------Open file to store graphs, fits and histograms-------------------

        // folder to store results files
        std::filesystem::create_directories("results");
        // results root files
        std::unique_ptr<TFile> myFile( TFile::Open("results/CV_CSIS.root", "RECREATE") );
        // Create directories inside root file 
        TDirectory* dirCV = myFile->mkdir("CV_graphs");
        TDirectory* dirDepletion = myFile->mkdir("Depletion_voltage");
        TDirectory* dirDonnor = myFile->mkdir("Donnor_density");

        // ---------------------Get how many different voltages an channels were tested---------------------

        int n_volt=0; //number of different voltages tested
        int n_ch=0; //number of different voltages tested

        for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
            // load the data for the given tree entry
            tree->GetEntry(iEntry);
            n_volt=iEntry+1;
        }
        n_ch = channel->size();
        cout<<"Number of channels: " << n_ch << endl;

        // define variables for the CV graph
        std::vector<float> x ;
        std::vector<float> y ;
        std::vector<float> yerr ;

       //------------------------------------LOOP over all channels: CV curves and fill histograms----------------------------------------

        for (int indx = 0; indx < n_ch; ++indx) { // loop over all channels

            for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
                // load the data for the given tree entry
                tree->GetEntry(iEntry);
        
                // fill variables for the CV graph
                x.push_back(voltage);
                y.push_back(cs->at(indx));
                yerr.push_back(cs_err->at(indx));
            }

            // which channel analysing
            int ch = channel->at(indx);
            printf("Channel: %d\n", ch);

            // -------------------CV GRAPH--------------------
            std::vector<float> xerr(n_volt, 0); // no error on x-axis
            TGraph *g = new TGraphErrors(n_volt, &x[0], &y[0], &xerr[0], &yerr[0]);
            g->SetName(Form("Channel %d", ch));
            g->SetTitle(Form("Channel %d;Voltage [V]; Capacitance [pF]", ch));
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.75);
            g->SetMarkerColor(kBlue);
            //g->Draw();
            g->GetXaxis()->CenterTitle();
            g->GetYaxis()->CenterTitle();

            // -------------------CV GRAPH: log scale for depletion voltage--------------------
            // first create new vectors to hold the log-transformed data and the propagated errors:
            std::vector<float> x_log(n_volt), y_log(n_volt);
            std::vector<float> yerr_log(n_volt);

            log_scale(n_volt, x, y, yerr, x_log, y_log, yerr_log);// modifies log variables

            TGraphErrors* glog = new TGraphErrors(n_volt, x_log.data(), y_log.data(), &xerr[0], yerr_log.data());
            glog->SetName(Form("Channel %d log scale", ch));
            glog->SetTitle(Form("Channel %d log scale;ln V; ln C", ch));
            glog->SetMarkerStyle(20);
            glog->SetMarkerSize(0.75);
            glog->SetMarkerColor(kBlue);

            // create canvas to draw graph and fit
            TCanvas* c1 = new TCanvas("CV_dep_volt_channel_", Form("CV_dep_volt_channel_%d", ch), 800, 600);
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

            // -----------------Fit and calculate depletion voltage------------------

            // // first fit: left region (line with slope)
            // glog->Fit("pol1", "0E", "", x_log[1], x_log[6]);
            // TF1* lfit = (TF1*)glog->GetFunction("pol1")->Clone("lfit");
            // lfit->SetLineColor(kRed);
            // lfit->SetLineWidth(2);
            // lfit->SetLineStyle(2); // dashed
            // lfit->SetRange(x_log[1], x_log[6]);
            // lfit->Draw("SAME");

            // second fit: right region (horizontal line)
            TF1* lfit = new TF1("lfit", "pol1", x_log[1], x_log[6]);
            glog->Fit(lfit, "QR0");
            lfit->SetLineColor(kGreen+3);
            lfit->SetLineWidth(2);
            lfit->SetLineStyle(7); // dotted
            lfit->SetRange(x_log[1], x_log[(n_volt-1)/3]);
            lfit->Draw("SAME");

            // second fit: right region (horizontal line)
            TF1* rfit = new TF1("rfit", "pol0", x_log[n_volt-(n_volt-1)/3], x_log[n_volt-1]);
            glog->Fit(rfit, "QR0");
            rfit->SetLineColor(kGreen+3);
            rfit->SetLineWidth(2);
            rfit->SetLineStyle(7); // dotted
            rfit->SetRange(x_log[n_volt-6],x_log[n_volt-1]);
            rfit->Draw("SAME");

            // calculate intersection (intersection point is depletion voltage V_dep)
            double p0_1 = lfit->GetParameter(0);
            double p1_1 = lfit->GetParameter(1);
            double p0_2 = rfit->GetParameter(0);
            double p1_2 = 0; // fit to a constaqnt, so slope is 0
            double log_Vdep = (p0_2 - p0_1) / (p1_1 - p1_2); // the lines cross at: p0_1 + p1_1 * x = p0_2 + p1_2 * x
            double V_dep = std::exp(log_Vdep); // convert back to linear scale
            // std::cout << "Depletion voltage V_dep = " << V_dep << " V" << std::endl;

            double chi2_lfit = lfit->GetChisquare();         // chi-squared
            int ndf_lfit = lfit->GetNDF();                   // number of degrees of freedom
            double pval_lfit = TMath::Prob(chi2_lfit, ndf_lfit);

            double chi2_rfit = rfit->GetChisquare();         // chi-squared
            int ndf_rfit = rfit->GetNDF();                   // number of degrees of freedom
            double pval_rfit = TMath::Prob(chi2_rfit, ndf_rfit);

            // store capacitance of the right plateau of CV in 2D histogram
            hcs_plateau_map->Fill(((ch-1)%16)+1 // X position (from 1 to 16)
                        , ((ch-1)/16)+1 // Y position
                        , std::exp(p0_2)); // value of constant line in V

            // // lets try to do the fit directly on CV (no log scale)
            // TF1* rfit2 = new TF1("rfit2", "pol0", x[n_volt-6], x[n_volt-1]);
            // g->Fit(rfit2, "QR0");

            // double p0_3 = rfit2->GetParameter(0);
            // cout << "fit to log scale: " << std::exp(p0_2) << std::endl << "fit to linear scale: " << p0_3 << std::endl;

            // store depletion voltage in histograms
            // hVdepxch->GetXaxis()->SetBinLabel(indx + 1, Form("Ch%d", ch)); // set bin label
            hVdepxch->SetBinContent(indx+1, V_dep); // channel index starts at 0
            hVdep->Fill(V_dep); // fill histogram with depletion voltage
            hVdep_map->Fill(((ch-1)%16)+1 // X position (from 1 to 16)
                            , ((ch-1)/16)+1 // Y position
                            , V_dep); // fill 2D histogram with depletion voltage

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


            //----------------------Donnor density-------------------------

            std::vector<float> y_new(n_volt);
            std::vector<float> y_new_err(n_volt);

            // calculate 1/cs^2
            for (int i = 0; i < n_volt; ++i) {
                y_new[i] = std::pow(y[i], -2);

                // error propagation formula for x^{-2}: simga_x^{-2} = 2*sigma_x / (x^3)
                y_new_err[i] = 2*yerr[i] / (pow(y[i],3));
            }

            TGraphErrors* gnew = new TGraphErrors(n_volt, &x[0], y_new.data(), &xerr[0], y_new_err.data());
            gnew->SetName(Form("Channel %d", ch));
            gnew->SetTitle(Form("Channel %d;Voltage [V]; 1/C^{2} [1/pF^{2}]", ch));
            gnew->SetMarkerStyle(20);
            gnew->SetMarkerSize(0.75);
            gnew->SetMarkerColor(kBlue);

            // create canvas
            TCanvas* c2 = new TCanvas("donnor_density_channel", Form("donnor_density_channel_%d", ch), 800, 600);
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

            // get donor density from slope of the fit
            double p1_don = don_fit->GetParameter(1); //slope of the fit in [V^{-1}pF^{-2}]
            p1_don = p1_don * std::pow(10,24); // convert from pF^{-2} to F^{-2}
            double donor_density = (2)/(e*eps*std::pow(A,2)*p1_don*std::pow(10,6)); //donor density in [number elctrons*cm^{-3}]
            // std::cout << "Donnor density = " << donor_density << " ne*cm^{-3}" << std::endl;

            // fill histogram with donor density
            hndon->Fill(donor_density);
            hndon_map->Fill(((ch-1)%16)+1 // X position (from 1 to 16)
                            , ((ch-1)/16)+1 // Y position
                            , donor_density); // fill 2D histogram with depletion voltage

            dirCV->cd();
            g->Write();

            dirDepletion->cd();
            c1->Write(); // log scale graph and fit in the same canvas

            dirDonnor->cd();
            c2->Write();

            //clean data vectors
            x.clear();
            y.clear();
            yerr.clear();
        } // end of channel loop

        // ------------------------------------Create 2D map of sensor---------------------------------------------------

        // depletion voltage

        gStyle->SetPalette(kBlueRedYellow); // Set default color palette
       
        auto cVdep = new TCanvas("cVdep", "Canvas", 600, 600);
        cVdep->cd();
        hVdep_map->SetMinimum(43);
        hVdep_map->SetMaximum(58);
        hVdep_map->Draw("COLZ");
        hVdep_map->GetZaxis()->SetTitle("Vdep [V]");
        channel_name->Draw("text same");

        // donor density
        auto cndon = new TCanvas("cndon_map", "Canvas", 600, 600);
        cndon->cd();
        hndon_map->SetMinimum(1.9e+11);
        hndon_map->SetMaximum(2.5e+11);
        hndon_map->Draw("COLZ");
        hndon_map->GetZaxis()->SetTitle("n_{don} [ne/cm^{3}]");
        channel_name->Draw("text same");

        // capacitance
        auto ccsplat = new TCanvas("ccsplat", "Canvas", 600, 600);
        ccsplat->cd();
        // hcs_plateau_map->SetMinimum(5.5);
        // hcs_plateau_map->SetMaximum(6.5);
        hcs_plateau_map->SetMinimum(5.4);
        hcs_plateau_map->SetMaximum(6);
        hcs_plateau_map->Draw("COLZ");
        hcs_plateau_map->GetZaxis()->SetTitle("Capacitance [pF]");
        channel_name->Draw("text same");

        myFile->cd();
        hVdepxch->Write(); // write histogram with depletion voltages per channel
        hVdep->Write(); // write histogram with depletiozn voltages distribution
        hndon->Write(); // write histogram with donor density distribution
        cVdep->Write();
        cndon->Write();
        ccsplat->Write();
    } // end isCV
    


    
    else {
        cout << "Processing IV data...\n";
        // Enable only required branches
        const std::vector<std::string> branches = {"voltage", "channel", "current", "current_err"};
        choose_branches(tree, branches);

        // Bind branches to variables
        float voltage = 0.0f;
        std::vector<int>* channel = 0;
        std::vector<float>* current = 0;
        std::vector<float>* current_err = 0;
        tree->SetBranchAddress("voltage", &voltage);
        tree->SetBranchAddress("channel", &channel);
        tree->SetBranchAddress("current", &current);
        tree->SetBranchAddress("current_err", &current_err);

        // --------------------Open file to store graphs, fits and histograms-------------------

        // folder to store results files
        std::filesystem::create_directories("results");
        // results root files
        std::unique_ptr<TFile> myFile( TFile::Open("results/IV_CSIS.root", "RECREATE") );

        //------------------------------------LOOP over all channels: CV curves and fill histograms----------------------------------------

        // define variables for the CV graph
        std::vector<float> x ;
        std::vector<float> y ;
        std::vector<float> yerr ;

        // ---------------------Get how many different voltages an channels were tested---------------------

        int n_volt=0; //number of different voltages tested
        int n_ch; //number of different voltages tested
        for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
            // load the data for the given tree entry
            tree->GetEntry(iEntry);
            n_volt=iEntry+1;
        }
        n_ch = channel->size();
        cout<<"Number of channels: " << n_ch << endl;
        
        for (int indx = 0; indx < n_ch; ++indx) { // loop over all channels

            //------------------------------------Load all data from given channel----------------------------------------
            for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
                // load the data for the given tree entry
                tree->GetEntry(iEntry);
                // fill variables for the CV graph
                x.push_back(voltage);
                y.push_back(current->at(indx));
                yerr.push_back(current_err->at(indx));
            }

            // which channel analysing
            int ch = channel->at(indx);
            printf("Channel: %d\n", ch);

            // -------------------CV GRAPH--------------------
            std::vector<float> xerr(n_volt, 0); // no error on x-axis
            TGraph *g = new TGraphErrors(n_volt, &x[0], &y[0], &xerr[0], &yerr[0]);
            g->SetName(Form("Channel %d", ch));
            g->SetTitle(Form("Channel %d;Voltage [V]; Current [nA]", ch));
            g->SetMarkerStyle(20);
            g->SetMarkerSize(0.75);
            g->SetMarkerColor(kBlue);
            //g->Draw();
            g->GetXaxis()->CenterTitle();
            g->GetYaxis()->CenterTitle();    

            myFile->cd();
            g->Write();

            //clean data vectors
            x.clear();
            y.clear();
            yerr.clear();
        }    
    }

    return 0;
}