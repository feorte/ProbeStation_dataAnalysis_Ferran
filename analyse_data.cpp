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
const float A = 30.25e-6; // pad area in m^2; 5.5 × 5.5 mm^2 = 30.25 mm^2 = 30.25e-6 m^2

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


int analyse_data(std::string number_sensor = "20", std::string type = "IV")
{
    std::string storingfile = "stored_data/stored_data_" + number_sensor + "_" + type + ".root";
    
    // ------------------------------------Load data from tree---------------------------------------------------

    // std::string storingfile = "stored_data/stored_data_CV.root"; // replace with your file name

    bool CSIS = false; // set to true if the data is from CSIS, false if from CSIS2
    int CSIS_ch_map[264] = {48,208,192,240,224,144,128,176,160,80,64,112,96,256,16,32,47,207,191,239,223,143,127,175,159,79,63,111,95,255,15,31,46,206,190,238,222,142,126,174,158,78,62,110,94,254,14,30,45,205,189,237,221,141,125,173,157,77,61,109,93,253,13,29,44,204,188,236,220,140,124,172,156,76,60,108,92,252,12,28,43,203,187,235,219,139,123,171,155,75,59,107,91,251,11,27,42,202,186,234,218,138,122,170,154,74,58,106,90,250,10,26,41,201,185,233,217,137,121,169,153,73,57,105,89,249,9,25,40,200,184,232,216,136,120,168,152,72,56,104,88,248,8,24,39,199,183,231,215,135,119,167,151,71,55,103,87,247,7,23,38,198,182,230,214,134,118,166,150,70,54,102,86,246,6,22,37,197,181,229,213,133,117,165,149,69,53,101,85,245,5,21,36,196,180,228,212,132,116,164,148,68,52,100,84,244,4,20,35,195,179,227,211,131,115,163,147,67,51,99,83,243,3,19,34,194,178,226,210,130,114,162,146,66,50,98,82,242,2,18,33,193,177,225,209,129,113,161,145,65,49,97,81,241,1,17};
    
    std::vector<int> border_channels = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
                                        17,32,33,48,64,65,80,81,96,97,112,113,128,129,144,145,160,161,176,177,192,193,208,209,224,225,240,241,
                                        241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256}; // borders of channels in the sensor
    
                                        // Load ROOT file with data
    auto file = std::unique_ptr<TFile>(TFile::Open(storingfile.c_str()));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file with data.\n";
        return -1;
    }

    // Load ROOT file with system capacitance
    auto file_systCap = std::unique_ptr<TFile>(TFile::Open("stored_data/systemCapacitanceCV.root"));
    if (!file || file_systCap->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file with system capacitance.\n";
        return -1;
    }

    // Load tree with data from file
    auto tree_data = file->Get<TTree>("analysis");
    if (!tree_data) {
        std::cerr << "Error: Tree 'analysis' not found in the file with data.\n";
        return -1;
    }

    // Load tree with system capacitance
    auto tree_systCap = file_systCap->Get<TTree>("analysis");
    if (!tree_data) {
        std::cerr << "Error: Tree 'analysis' not found in the file with system capacitance.\n";
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
        const std::vector<std::string> branches = {"voltage", "channel", "cs_uncorr", "cs_err"};
        choose_branches(tree_data, branches);
        choose_branches(tree_systCap, branches);

        // Bind branches to variables
        float voltage = 0.0f;
        std::vector<int>* channel = 0;
        std::vector<float>* cs = 0;
        std::vector<float>* cs_err = 0;
        std::vector<float>* cs_syst = 0;
        std::vector<float>* cs_syst_err = 0;
        tree_data->SetBranchAddress("voltage", &voltage);
        tree_data->SetBranchAddress("channel", &channel);
        tree_data->SetBranchAddress("cs_uncorr", &cs);
        tree_data->SetBranchAddress("cs_err", &cs_err);
        tree_systCap->SetBranchAddress("cs_uncorr", &cs_syst);
        tree_systCap->SetBranchAddress("cs_err", &cs_syst_err);

        // ------------------------------------Create histograms---------------------------------------------------

        // 1D histograms to store depletion voltages and donor density distributions
        TH1F* hVdepxch = new TH1F("hVdepxch", "Depletion Voltage per Channel;Channel;V_{dep} [V]", 256, 0.5, 256.5); // histogram to store depletion voltages per channel
        // TH1F* hVdep = new TH1F("hVdep", "Depletion Voltage;V_{dep} [V];Entries", 200, 17, 80); // histogram to store depletion voltages distribution
        TH1F* hndonxch = new TH1F("hndonxch", "Donnor density per Channel;Channel;Donnor density [ne/cm^{3}]", 256, 0.5, 256.5); // histogram to store depletion voltages per channel
        // TH1F* hndon = new TH1F("hndon", "Donnor density;Donnor density [ne/cm^{3}];Entries", 200, 3e+10, 1e+12); // histogram to store depletion voltages distribution
        // TH1F* hndon_border = new TH1F("hndon_border", "Donnor density at borders;Donnor density [ne/cm^{3}];Entries", 200, 3e+10, 1e+12); // histogram to store depletion voltages distribution at borders
        // TH1F* hndon_inner = new TH1F("hndon_inner", "Donnor density in inner channels;Donnor density [ne/cm^{3}];Entries", 200, 3e+10, 1e+12); // histogram to store depletion voltages distribution in inner channels
        TH1F* hcsxch = new TH1F("hcsxch", "Full depletion capacitance per Channel;Channel;Capacitance [pF]", 256, 0.5, 256.5); // histogram to store depletion voltages per channel
        // TH1F* hcs = new TH1F("hcs", "Full depletion capacitance;Capacitance [pF];Entries", 50, 5.2, 6.2); // histogram to store depletion voltages distribution
        // TH1F* hcs_border = new TH1F("hcs_border", "Full depletion capacitance at borders;Capacitance [pF];Entries", 50, 5.2, 6.2); // histogram to store depletion voltages distribution at borders
        // TH1F* hcs_inner = new TH1F("hcs_inner", "Full depletion capacitance in inner channels;Capacitance [pF];Entries", 50, 5.2, 6.2); // histogram to store depletion voltages distribution in inner channels
        TH1F* hchi2xch = new TH1F("hchi2xch", "Global chi2 per Channel;Channel;Chi2", 256, 0.5, 256.5); // histogram to store depletion voltages per channel
        TH1F* hchi2 = new TH1F("hchi2","Global chi2;Chi2;Entries", 150, 0, 2e5); //histogram to store global chi2 values
        std::vector<float> chi2_glob; // vector to store global chi2 values

        // 2D histogram to store capacitance, depletion voltage and donor density of every channel
        auto hVdep_map = new TH2F("hVdep_map",Form("Depletion Voltage for sensor %s;X;Y",number_sensor.c_str()),
            16,0.5,16.5,  // X axis
            16,0.5,16.5); // Y axis
        auto hndon_map = new TH2F("hndon_map",Form("Donnor density for sensor %s;X;Y",number_sensor.c_str()),
            16,0.5,16.5,  // X axis
            16,0.5,16.5); // Y axis
        auto hcs_plateau_map = new TH2F("hcs_plateau_map",Form("Full depletion capacitance for sensor %s;X;Y",number_sensor.c_str()),
            16,0.5,16.5,  // X axis
            16,0.5,16.5); // Y axis
        auto hchi2sus = new TH2F("hchi2sus","Weird chi2;X;Y", 16,0.5,16.5, 16,0.5,16.5); //chi2 suspicius


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
        std:: string results_file = "results/CV_" + number_sensor + ".root";
        std::unique_ptr<TFile> myFile( TFile::Open(results_file.c_str(), "RECREATE") );
        // Create directories inside root file 
        TDirectory* dirCV = myFile->mkdir("CV_graphs");
        TDirectory* dirDepletion = myFile->mkdir("Depletion_voltage");
        TDirectory* dirDonnor = myFile->mkdir("Donnor_density");

        // ---------------------Get how many different voltages an channels were tested---------------------

        int n_volt=0; //number of different voltages tested
        int n_ch=0; //number of different voltages tested

        for (int iEntry = 0; tree_data->LoadTree(iEntry) >= 0; ++iEntry) {
            // load the data for the given tree entry
            tree_data->GetEntry(iEntry);
            n_volt=iEntry+1;
        }
        n_ch = channel->size();
        // cout<<"Number of channels: " << n_ch << endl;

        // define variables for the CV graph
        std::vector<float> x ;
        std::vector<float> y ;
        std::vector<float> yerr ;
        float y_syst; // no error on x-axis
        float yerr_syst; // no error on x-axis

        // vectors to store main values
        std::vector<float> vect_vdep; // vector to store depletion voltages
        std::vector<float> vect_donor_density; // vector to store donor densities
        std::vector<float> vect_donor_density_inner; // vector to store donor densities
        std::vector<float> vect_donor_density_border; // vector to store donor densities at borders
        std::vector<float> vect_capacitance; // vector to store capacitances
        std::vector<float> vect_capacitance_inner; // vector to store capacitances
        std::vector<float> vect_capacitance_border; // vector to store capacitances at borders


       //------------------------------------LOOP over all channels: CV curves and fill histograms----------------------------------------

        for (int indx = 0; indx < n_ch; ++indx) { // loop over all channels

            for (int iEntry = 0; tree_systCap->LoadTree(iEntry) >= 0; ++iEntry) {
                // load the data for the given tree entry
                tree_systCap->GetEntry(iEntry);
                y_syst = cs_syst->at(indx);
                yerr_syst = cs_syst_err->at(indx);
            }

            for (int iEntry = 0; tree_data->LoadTree(iEntry) >= 0; ++iEntry) {
                // load the data for the given tree entry
                tree_data->GetEntry(iEntry);
                // fill variables for the CV graph
                x.push_back(voltage);
                y.push_back(cs->at(indx));
                yerr.push_back(cs_err->at(indx));
            }

            int ch; // current channel analysing
            // which channel analysing
            if (CSIS) {
                // map CSIS channels to sensor channels
                ch = CSIS_ch_map[channel->at(indx)-1];
                // printf("Channel: %d\n", ch);
            } else {
                // use channel number as is
                ch = channel->at(indx);
                // printf("Channel: %d\n", ch);
            }

            if (ch < 1 || ch > 256) {
                std::cerr << "Error: Channel number out of range (1-256): " << ch << std::endl;
                continue; // skip this channel
            }

            // Substract system capacitance from the measured capacitance
            for (int i = 0; i < n_volt; ++i) {
                y[i] -= y_syst;
                yerr[i] = std::sqrt(yerr[i]*yerr[i] + yerr_syst*yerr_syst); // propagate error
            }

            // -------------------CV GRAPH--------------------
            std::vector<float> xerr(n_volt, 0); // no error on x-axis
            TGraph *g = new TGraphErrors(n_volt, &x[0], &y[0], &xerr[0], &yerr[0]);
            g->SetName(Form("Channel_%d", ch));
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
            TCanvas* c1 = new TCanvas(Form("Channel_%d", ch), Form("Channel_%d", ch), 800, 600);
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

            // first fit: left region (line with slope)
            // TF1* lfit = new TF1("lfit", "pol1", x_log[1], x_log[(n_volt-1)/3]);
            TF1* lfit = new TF1("lfit", "pol1", x_log[1], std::log(35)); // fit to a line with slope, range from first point to 35 V
            glog->Fit(lfit, "EQR0");
            lfit->SetLineColor(kGreen+3);
            lfit->SetLineWidth(2);
            lfit->SetLineStyle(7); // dotted

            // second fit: right region (horizontal line)
            // TF1* rfit = new TF1("rfit", "pol0", x_log[n_volt-(n_volt-1)/3], x_log[n_volt-1]);
            TF1* rfit = new TF1("rfit", "pol0", std::log(60), x_log[n_volt-1]);// fit to a constant, range from 60 V to last point
            glog->Fit(rfit, "EQR0");
            rfit->SetLineColor(kRed);
            rfit->SetLineWidth(2);
            rfit->SetLineStyle(7); // dotted

            // calculate intersection (intersection point is depletion voltage V_dep)
            double p0_1 = lfit->GetParameter(0);
            double p1_1 = lfit->GetParameter(1);
            double p0_2 = rfit->GetParameter(0);
            double p1_2 = 0; // fit to a constaqnt, so slope is 0
            double log_Vdep = (p0_2 - p0_1) / (p1_1 - p1_2); // the lines cross at: p0_1 + p1_1 * x = p0_2 + p1_2 * x
            double V_dep = std::exp(log_Vdep); // convert back to linear scale
            // std::cout << "Depletion voltage V_dep = " << V_dep << " V" << std::endl;

            // draw fit lines
            lfit->SetRange(x_log[1], log_Vdep+1);
            rfit->SetRange(log_Vdep-1, x_log[n_volt-1]);
            lfit->Draw("SAME");
            rfit->Draw("SAME");

            double chi2_lfit = lfit->GetChisquare();         // chi-squared
            // int ndf_lfit = lfit->GetNDF();                   // number of degrees of freedom
            // double pval_lfit = TMath::Prob(chi2_lfit, ndf_lfit);
            double chi2_rfit = rfit->GetChisquare();         // chi-squared
            // int ndf_rfit = rfit->GetNDF();                   // number of degrees of freedom
            // double pval_rfit = TMath::Prob(chi2_rfit, ndf_rfit);
            chi2_glob.push_back(chi2_lfit + chi2_rfit); // sum of chi-squared values
            if (!std::isnan(V_dep) && std::isfinite(V_dep)) {
                hchi2xch->SetBinContent(indx+1, chi2_lfit + chi2_rfit); // store chi2 value in histogram
            } else {
                std::cerr << "Warning: Skipping invalid chi2_glob for channel index " << indx+1 << std::endl;
                //hchi2sus->Fill(((ch-1)%16)+1, ((ch-1)/16)+1, 1); 
            }
            hchi2->Fill(chi2_lfit + chi2_rfit); // store chi2 value in histogram
            
            //store capacitance in histograms
            hcsxch->SetBinContent(indx+1, std::exp(p0_2)); // channel index starts at 0
            // hcs->Fill(std::exp(p0_2));
            vect_capacitance.push_back(std::exp(p0_2));
            // store capacitance of the right plateau of CV in 2D histogram
            hcs_plateau_map->Fill(((ch-1)%16)+1 // X position (from 1 to 16)
                        , ((ch-1)/16)+1 // Y position
                        , std::exp(p0_2)); // value of constant line in V

            // check if channel is at the border or in the inner part of the sensor
            bool is_border = std::find(border_channels.begin(), border_channels.end(), ch) != border_channels.end();
            if (is_border) {
                // hcs_border->Fill(std::exp(p0_2)); // fill histogram with capacitance at borders
                vect_capacitance_border.push_back(std::exp(p0_2));
            }
            else {
                // hcs_inner->Fill(std::exp(p0_2)); // fill histogram with capacitance in inner channels
                vect_capacitance_inner.push_back(std::exp(p0_2));
            }

                    // // lets try to do the fit directly on CV (no log scale)
                    // TF1* rfit2 = new TF1("rfit2", "pol0", x[n_volt-6], x[n_volt-1]);
                    // g->Fit(rfit2, "QR0");

                    // double p0_3 = rfit2->GetParameter(0);
                    // cout << "fit to log scale: " << std::exp(p0_2) << std::endl << "fit to linear scale: " << p0_3 << std::endl;

            // store depletion voltage in histograms
            // hVdepxch->GetXaxis()->SetBinLabel(indx + 1, Form("Ch%d", ch)); // set bin label
            if (!std::isnan(V_dep) && std::isfinite(V_dep)) {
                hVdepxch->SetBinContent(indx+1, V_dep);
            } else {
                std::cerr << "Warning: Skipping invalid V_dep for channel index " << indx+1 << std::endl; 
            }
            // hVdep->Fill(V_dep); // fill histogram with depletion voltage
            vect_vdep.push_back(V_dep); // store depletion voltage in vector
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
            legend->AddEntry(rfit, Form("Plateau fit = %.2f pF", std::exp(p0_2)), "l");
            legend->AddEntry(line, Form("V_{dep} = %.2f V", V_dep), "l");
            legend->Draw("SAME");

            // add


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
            dirDonnor->cd();
            // TCanvas* c2 = new TCanvas("donnor_density_channel", Form("donnor_density_channel_%d", ch), 800, 600);
            TCanvas* c2 = new TCanvas(Form("donn_Channel_%d", ch), Form("Channel_%d", ch), 800, 600);


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

            // fit: left region
            gnew->Fit("pol1", "Q0", "", x[1], 35);
            TF1* don_fit = (TF1*)gnew->GetFunction("pol1")->Clone("don_fit");
            don_fit->SetLineColor(kRed);
            don_fit->SetLineWidth(2);
            don_fit->SetLineStyle(2); // dashed
            don_fit->SetRange(x[0], x[8]);
            don_fit->Draw("SAME");

            // get donor density from slope of the fit
            double p1_don = don_fit->GetParameter(1); //slope of the fit in [V^{-1}pF^{-2}]
            p1_don = p1_don * std::pow(10,24); // convert from pF^{-2} to F^{-2}
            double donor_density = (2)/(e*eps*std::pow(A,2)*p1_don)*1e-6; //donor density in [number elctrons*cm^{-3}]
            // std::cout << "Donnor density = " << donor_density << " ne*cm^{-3}" << std::endl;

            // fill histogram with donor density
            hndonxch->SetBinContent(indx+1, donor_density); // channel index starts at 0
            // hndon->Fill(donor_density);
            vect_donor_density.push_back(donor_density); // store donor density in vector
            hndon_map->Fill(((ch-1)%16)+1 // X position (from 1 to 16)
                            , ((ch-1)/16)+1 // Y position
                            , donor_density); // fill 2D histogram with depletion voltage

            // already checked if channel is at the border or in the inner part of the sensor
            if (is_border) {
                // hndon_border->Fill(donor_density); // fill histogram with donor density at borders
                vect_donor_density_border.push_back(donor_density); // store donor density at borders in vector
            } else {
                // hndon_inner->Fill(donor_density); // fill histogram with donor density in inner channels
                vect_donor_density_inner.push_back(donor_density); // store donor density in inner channels in vector
            }

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

        //------------------------Set ranges for histograms-------------------------

        //range for 1D histograms
        // find min and max values in the vectors
        double min_vdep = *std::min_element(vect_vdep.begin(), vect_vdep.end());
        double max_vdep = *std::max_element(vect_vdep.begin(), vect_vdep.end());
        double min_ndon = *std::min_element(vect_donor_density.begin(), vect_donor_density.end());
        double max_ndon = *std::max_element(vect_donor_density.begin(), vect_donor_density.end());
        double min_ndon_border = *std::min_element(vect_donor_density_border.begin(), vect_donor_density_border.end());
        double max_ndon_border = *std::max_element(vect_donor_density_border.begin(), vect_donor_density_border.end());
        double min_ndon_inner = *std::min_element(vect_donor_density_inner.begin(), vect_donor_density_inner.end());
        double max_ndon_inner = *std::max_element(vect_donor_density_inner.begin(), vect_donor_density_inner.end());
        double min_cs = *std::min_element(vect_capacitance.begin(), vect_capacitance.end());
        double max_cs = *std::max_element(vect_capacitance.begin(), vect_capacitance.end());
        double min_cs_border = *std::min_element(vect_capacitance_border.begin(), vect_capacitance_border.end());
        double max_cs_border = *std::max_element(vect_capacitance_border.begin(), vect_capacitance_border.end());
        double min_cs_inner = *std::min_element(vect_capacitance_inner.begin(), vect_capacitance_inner.end());
        double max_cs_inner = *std::max_element(vect_capacitance_inner.begin(), vect_capacitance_inner.end());

        TH1F* hVdep = new TH1F("hVdep", "Depletion Voltage;V_{dep} [V];Entries", 200, min_vdep, max_vdep); // histogram to store depletion voltages distribution
        TH1F* hndon = new TH1F("hndon", "Donnor density;Donnor density [ne/cm^{3}];Entries", 200, min_ndon, max_ndon); // histogram to store depletion voltages distribution
        TH1F* hndon_border = new TH1F("hndon_border", "Donnor density at borders;Donnor density [ne/cm^{3}];Entries", 200, min_ndon_border, max_ndon_border); // histogram to store depletion voltages distribution at borders
        TH1F* hndon_inner = new TH1F("hndon_inner", "Donnor density in inner channels;Donnor density [ne/cm^{3}];Entries", 200, min_ndon_inner, max_ndon_inner); // histogram to store depletion voltages distribution in inner channels
        TH1F* hcs = new TH1F("hcs", "Full depletion capacitance;Capacitance [pF];Entries", 200, min_cs, max_cs); // histogram to store depletion voltages distribution
        TH1F* hcs_border = new TH1F("hcs_border", "Full depletion capacitance at borders;Capacitance [pF];Entries", 200, min_cs_border, max_cs_border); // histogram to store depletion voltages distribution at borders
        TH1F* hcs_inner = new TH1F("hcs_inner", "Full depletion capacitance in inner channels;Capacitance [pF];Entries", 200, min_cs_inner, max_cs_inner); // histogram to store depletion voltages distribution in inner channels

        for (const auto& val : vect_vdep) {
            hVdep->Fill(val);
        }
        for (const auto& val : vect_donor_density) {
            hndon->Fill(val);
        }
        for (const auto& val : vect_donor_density_border) {
            hndon_border->Fill(val);
        }
        for (const auto& val : vect_donor_density_inner) {
            hndon_inner->Fill(val);
        }
        for (const auto& val : vect_capacitance) {
            hcs->Fill(val);
        }
        for (const auto& val : vect_capacitance_border) {
            hcs_border->Fill(val);
        }
        for (const auto& val : vect_capacitance_inner) {
            hcs_inner->Fill(val);
        }
        // set ranges for histogram
        float Vdep_mean = hVdep->GetMean();
        float Vdep_std = hVdep->GetStdDev();
        // cout<< "Vdep mean: " << Vdep_mean << ", Vdep std: " << hVdep->GetStdDev() << std::endl;
        // hVdep->GetXaxis()->SetRangeUser(Vdep_mean - 1*Vdep_std, Vdep_mean + 1*Vdep_std);
        // cout<< "Vdep mean after range: " << Vdep_mean << std::endl; // this does not change the mean, just the visualizationn range
        float ndon_border_mean = hndon_border->GetMean();
        float ndon_border_std = hndon_border->GetStdDev();
        // cout<< "ndon border mean: " << ndon_border_mean << ", ndon border std: " << hndon_border->GetStdDev() << std::endl;
        // hndon_border->GetXaxis()->SetRangeUser(ndon_border_mean - 2*ndon_border_std, ndon_border_mean + 2*ndon_border_std);
        float ndon_inner_mean = hndon_inner->GetMean();
        float ndon_inner_std = hndon_inner->GetStdDev();
        // cout<< "ndon inner mean: " << ndon_inner_mean << ", ndon inner std: " << hndon_inner->GetStdDev() << std::endl;
        // hndon_inner->GetXaxis()->SetRangeUser(ndon_inner_mean - 2*ndon_inner_std, ndon_inner_mean + 2*ndon_inner_std);
        float cs_border_mean = hcs_border->GetMean();
        float cs_border_std = hcs_border->GetStdDev();
        // cout<< "cs border mean: " << cs_border_mean << ", cs border std: " << hcs_border->GetStdDev() << std::endl;
        // hcs_border->GetXaxis()->SetRangeUser(cs_border_mean - 2*cs_border_std, cs_border_mean + 2*cs_border_std);
        float cs_inner_mean = hcs_inner->GetMean();
        float cs_inner_std = hcs_inner->GetStdDev();
        // cout<< "cs inner mean: " << cs_inner_mean << ", cs inner std: " << hcs_inner->GetStdDev() << std::endl;
        // hcs_inner->GetXaxis()->SetRangeUser(cs_inner_mean - 2*cs_inner_std, cs_inner_mean + 2*cs_inner_std);

        TH1F* hVdep_dist = new TH1F("hVdep_dist", "Depletion Voltage;V_{dep} [V];Entries", 25, Vdep_mean - 1*Vdep_std, Vdep_mean + 1*Vdep_std); // histogram to store depletion voltages distribution
        TH1F* hndon_dist = new TH1F("hndon_dist", "Donnor density;Donnor density [ne/cm^{3}];Entries", 50, ndon_inner_mean - 1*ndon_inner_std, ndon_border_mean + 1*ndon_border_std); // histogram to store depletion voltages distribution
        TH1F* hndon_border_dist = new TH1F("hndon_border_dist", "Donnor density at borders;Donnor density [ne/cm^{3}];Entries", 50, ndon_inner_mean - 1*ndon_inner_std, ndon_border_mean + 1*ndon_border_std); // histogram to store depletion voltages distribution at borders
        TH1F* hndon_inner_dist = new TH1F("hndon_inner_dist", "Donnor density in inner channels;Donnor density [ne/cm^{3}];Entries", 50, ndon_inner_mean - 1*ndon_inner_std, ndon_inner_mean + 1*ndon_inner_std); // histogram to store depletion voltages distribution in inner channels
        TH1F* hcs_dist = new TH1F("hcs_dist", "Full depletion capacitance;Capacitance [pF];Entries", 50, cs_inner_mean - 1*cs_inner_std, cs_border_mean + 1*cs_border_std); // histogram to store depletion voltages distribution
        TH1F* hcs_border_dist = new TH1F("hcs_border_dist", "Full depletion capacitance at borders;Capacitance [pF];Entries", 50, cs_border_mean - 1*cs_border_std, cs_border_mean + 1*cs_border_std); // histogram to store depletion voltages distribution at borders
        TH1F* hcs_inner_dist = new TH1F("hcs_inner_dist", "Full depletion capacitance in inner channels;Capacitance [pF];Entries", 50, cs_inner_mean - 1*cs_inner_std, cs_inner_mean + 1*cs_inner_std); // histogram to store depletion voltages distribution in inner channels

        // Fill the histogram with each value in the vector
        for (const auto& val : vect_vdep) {
            hVdep_dist->Fill(val);
        }
        for (const auto& val : vect_donor_density) {
            hndon_dist->Fill(val);
        }
        for (const auto& val : vect_donor_density_border) {
            hndon_border_dist->Fill(val);
        }
        for (const auto& val : vect_donor_density_inner) {
            hndon_inner_dist->Fill(val);
        }
        for (const auto& val : vect_capacitance) {
            hcs_dist->Fill(val);
        }
        for (const auto& val : vect_capacitance_border) {
            hcs_border_dist->Fill(val);
        }
        for (const auto& val : vect_capacitance_inner) {
            hcs_inner_dist->Fill(val);
        }

        float Vdep_dist_mean = hVdep_dist->GetMean();
        float Vdep_dist_std = hVdep_dist->GetStdDev();
        cout << "Vdep distribution mean: " << Vdep_dist_mean << ", Vdep distribution std: " << Vdep_dist_std << std::endl;
        float ndon_border_dist_mean = hndon_border_dist->GetMean();
        float ndon_border_dist_std = hndon_border_dist->GetStdDev();
        cout << "ndon border distribution mean: " << ndon_border_dist_mean << ", ndon border distribution std: " << ndon_border_dist_std << std::endl;
        float ndon_inner_dist_mean = hndon_inner_dist->GetMean();
        float ndon_inner_dist_std = hndon_inner_dist->GetStdDev();
        cout << "ndon inner distribution mean: " << ndon_inner_dist_mean << ", ndon inner distribution std: " << ndon_inner_dist_std << std::endl;
        float cs_border_dist_mean = hcs_border_dist->GetMean();
        float cs_border_dist_std = hcs_border_dist->GetStdDev();
        // cout << "cs border distribution mean: " << cs_border_dist_mean << ", cs border distribution std: " << cs_border_dist_std << std::endl;
        float cs_inner_dist_mean = hcs_inner_dist->GetMean();
        float cs_inner_dist_std = hcs_inner_dist->GetStdDev();
        // cout << "cs inner distribution mean: " << cs_inner_dist_mean << ", cs inner distribution std: " << cs_inner_dist_std << std::endl;

        // TH1F* hndon_border_dist2 = new TH1F("hndon_border_dist2", "Donnor density at borders;Donnor density [ne/cm^{3}];Entries", 20, ndon_border_dist_mean - 1*ndon_border_dist_std, ndon_border_dist_mean + 1*ndon_border_dist_std); // histogram to store depletion voltages distribution at borders
        // TH1F* hndon_inner_dist2 = new TH1F("hndon_inner_dist2", "Donnor density in inner channels;Donnor density [ne/cm^{3}];Entries", 50, ndon_inner_dist_mean - 1*ndon_inner_dist_std, ndon_inner_dist_mean + 1*ndon_inner_dist_std); // histogram to store depletion voltages distribution in inner channels
        TH1F* hcs_dist2 = new TH1F("hcs_dist2", "Full depletion capacitance;Capacitance [pF];Entries", 50, cs_inner_dist_mean - 1*cs_inner_dist_std, cs_border_dist_mean + 1*cs_border_dist_std); // histogram to store depletion voltages distribution
        TH1F* hcs_border_dist2 = new TH1F("hcs_border_dist2", "Full depletion capacitance at borders;Capacitance [pF];Entries", 50, cs_border_dist_mean - 1*cs_border_dist_std, cs_border_dist_mean + 1*cs_border_dist_std); // histogram to store depletion voltages distribution at borders
        TH1F* hcs_inner_dist2 = new TH1F("hcs_inner_dist2", "Full depletion capacitance in inner channels;Capacitance [pF];Entries", 50, cs_inner_dist_mean - 1*cs_inner_dist_std, cs_inner_dist_mean + 1*cs_inner_dist_std); // histogram to store depletion voltages distribution in inner channels

        // for (const auto& val : vect_donor_density_border) {
        //     hndon_border_dist2->Fill(val);
        // }
        // for (const auto& val : vect_donor_density_inner) {
        //     hndon_inner_dist2->Fill(val);
        // }
        for (const auto& val : vect_capacitance) {
            hcs_dist2->Fill(val);
        }
        for (const auto& val : vect_capacitance_border) {
            hcs_border_dist2->Fill(val);
        }
        for (const auto& val : vect_capacitance_inner) {
            hcs_inner_dist2->Fill(val);
        }
        // float ndon_border_dist_mean2 = hndon_border_dist2->GetMean();
        // float ndon_border_dist_std2 = hndon_border_dist2->GetStdDev();
        // cout << "ndon border distribution mean2: " << ndon_border_dist_mean2
        // << ", ndon border distribution std2: " << ndon_border_dist_std2 << std::endl;
        // float ndon_inner_dist_mean2 = hndon_inner_dist2->GetMean();
        // float ndon_inner_dist_std2 = hndon_inner_dist2->GetStdDev();
        // cout << "ndon inner distribution mean2: " << ndon_inner_dist_mean2
        // << ", ndon inner distribution std2: " << ndon_inner_dist_std2 << std::endl;
        float cs_border_dist_mean2 = hcs_border_dist2->GetMean();
        float cs_border_dist_std2 = hcs_border_dist2->GetStdDev();
        cout << "cs border distribution mean2: " << cs_border_dist_mean2
        << ", cs border distribution std2: " << cs_border_dist_std2 << std::endl;
        float cs_inner_dist_mean2 = hcs_inner_dist2->GetMean();
        float cs_inner_dist_std2 = hcs_inner_dist2->GetStdDev();
        cout << "cs inner distribution mean2: " << cs_inner_dist_mean2
        << ", cs inner distribution std2: " << cs_inner_dist_std2 << std::endl;

        
        hVdep_map->SetMinimum(Vdep_mean - 1*Vdep_std);
        hVdep_map->SetMaximum(Vdep_mean + 1*Vdep_std);
        hndon_map->SetMinimum(ndon_inner_dist_mean - 2*ndon_inner_dist_std);
        hndon_map->SetMaximum(ndon_border_dist_mean + 1*ndon_border_dist_std);
        hcs_plateau_map->SetMinimum(cs_inner_dist_mean - 1*cs_inner_dist_std);
        hcs_plateau_map->SetMaximum(cs_border_dist_mean + 1*cs_border_dist_std);





        //Apart from ranges, also find weird chi2 values and fill histogram with them
        // float chi2_min = hchi2->GetXaxis()->GetXmin();
        // float chi2_max = hchi2->GetXaxis()->GetXmax();
        // hchi2->GetXaxis()->SetRangeUser(-9, 9); // set range for chi2 histogram
        float chi2_mean = hchi2->GetMean();
        float chi2_std = hchi2->GetStdDev();
        cout << "Chi2 mean: " << chi2_mean << ", Chi2 std: " << chi2_std << std::endl;

        // ---------------------------Check for weird chi2 values and fill histogram with them--------------------------
        for (int i = 0; i<chi2_glob.size(); ++i) {
            int ch; // current channel analysing
            // which channel analysing
            if (CSIS) {
                // map CSIS channels to sensor channels
                ch = CSIS_ch_map[channel->at(i)-1];
            } else {
                // use channel number as is
                ch = channel->at(i);
            }
            if (ch < 1 || ch > 256) {
                continue; // skip this channel
            }
            // if (chi2_glob[i] < chi2_mean - 2*chi2_std || chi2_glob[i] > chi2_mean + 2*chi2_std) {
            if (chi2_glob[i] < 10) {
                // fill histogram with weird chi2 values
                hchi2sus->Fill(((ch-1)%16)+1, ((ch-1)/16)+1, 1); 
            }
        }

        // ------------------------------------Create 2D map of sensor---------------------------------------------------

        gStyle->SetPalette(kBlueRedYellow); // Set default color palette

        // depletion voltage
        auto cVdep = new TCanvas("cVdep", "Canvas", 600, 600);
        cVdep->cd();
        hVdep_map->Draw("COLZ");
        hVdep_map->GetZaxis()->SetTitle("Vdep [V]");
        channel_name->Draw("text same");

        // donor density
        auto cndon = new TCanvas("cndon_map", "Canvas", 600, 600);
        cndon->cd();
        hndon_map->Draw("COLZ");
        hndon_map->GetZaxis()->SetTitle("n_{don} [ne/cm^{3}]");
        channel_name->Draw("text same");

        // capacitance
        auto ccsplat = new TCanvas("ccsplat", "Canvas", 600, 600);
        ccsplat->cd();
        hcs_plateau_map->Draw("COLZ");
        hcs_plateau_map->GetZaxis()->SetTitle("Capacitance [pF]");
        channel_name->Draw("text same");

        // chi2 suspicius
        auto cchi2sus = new TCanvas("cchi2sus", "Canvas", 600, 600);
        cchi2sus->cd();
        hchi2sus->SetMinimum(0);
        hchi2sus->SetMaximum(1);
        hchi2sus->Draw("COLZ");
        hchi2sus->GetZaxis()->SetTitle("Weird chi2");
        channel_name->Draw("text same");

        myFile->cd();
        hVdepxch->Write(); // write histogram with depletion voltages per channel
        hVdep->Write(); // write histogram with depletiozn voltages distribution
        hndonxch->Write(); // write histogram with donor density per channel
        hndon->Write(); // write histogram with donor density distribution
        hndon_border->Write(); // write histogram with donor density at borders
        hndon_inner->Write(); // write histogram with donor density in inner channels
        hcsxch->Write(); // write histogram with capacitance per channel
        hcs->Write(); // write histogram with capacitance distribution
        hcs_border->Write(); // write histogram with capacitance at borders
        hcs_inner->Write(); // write histogram with capacitance in inner channels
        hchi2xch->Write(); // write histogram with chi2 values per channel
        hchi2->Write(); // write histogram with chi2 values
        cVdep->Write();
        cndon->Write();
        ccsplat->Write();
        cchi2sus->Write();

        hVdep_dist->Write(); // write histogram with depletion voltages distribution
        hndon_dist->Write(); // write histogram with donor density distribution
        hndon_border_dist->Write(); // write histogram with donor density at borders distribution
        hndon_inner_dist->Write(); // write histogram with donor density in inner channels distribution
        hcs_dist->Write(); // write histogram with capacitance distribution
        hcs_border_dist->Write(); // write histogram with capacitance at borders distribution
        hcs_inner_dist->Write(); // write histogram with capacitance in inner channels distribution
        // hndon_border_dist2->Write(); // write histogram with donor density at borders distribution
        // hndon_inner_dist2->Write(); // write histogram with donor density in inner channels
        hcs_dist2->Write(); // write histogram with capacitance distribution
        hcs_border_dist2->Write(); // write histogram with capacitance at borders distribution
        hcs_inner_dist2->Write(); // write histogram with capacitance in inner channels distribution
    } // end isCV
    



    
    
    else {
        cout << "Processing IV data...\n";
        // Enable only required branches
        const std::vector<std::string> branches = {"voltage", "channel", "current", "current_err"};
        choose_branches(tree_data, branches);

        // Bind branches to variables
        float voltage = 0.0f;
        std::vector<int>* channel = 0;
        std::vector<float>* current = 0;
        std::vector<float>* current_err = 0;
        tree_data->SetBranchAddress("voltage", &voltage);
        tree_data->SetBranchAddress("channel", &channel);
        tree_data->SetBranchAddress("current", &current);
        tree_data->SetBranchAddress("current_err", &current_err);

        // ------------------------------------Create histograms---------------------------------------------------
        // Histograms for current at certain voltage
        int voltage_check = 120; // voltage for which current is measured
        auto hcurrxch = new TH1F("hcurrxch", Form("Current per Channel at %d V;Channel;Current [nA]", voltage_check), 256, 0.5, 256.5); // histogram to store current per channel
        auto hcurr_map = new TH2F("hcurr_map",Form("Current at %d V for sensor %s;X;Y", voltage_check, number_sensor.c_str()), 16,0.5,16.5, 16,0.5,16.5);


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
        std:: string results_file = "results/IV_" + number_sensor + ".root";
        std::unique_ptr<TFile> myFile( TFile::Open(results_file.c_str(), "RECREATE") );
        // Create directories inside root file 
        TDirectory* dirIV = myFile->mkdir("IV_graphs");

        //------------------------------------LOOP over all channels: CV curves and fill histograms----------------------------------------

        // define variables for the CV graph
        std::vector<float> x ;
        std::vector<float> y ;
        std::vector<float> yerr ;

        // ---------------------Get how many different voltages an channels were tested. Also find position of voltage used to map current---------------------

        int n_volt=0; //number of different voltages tested
        int n_ch; //number of different channels tested
        int map_indx = -1; // index of the voltage used to map current
        for (int iEntry = 0; tree_data->LoadTree(iEntry) >= 0; ++iEntry) {
            // load the data for the given tree entry
            tree_data->GetEntry(iEntry);
            if (voltage==voltage_check){
                map_indx = iEntry;
            }
            n_volt=iEntry+1;
        }
        if (map_indx == -1) {
            std::cerr << "Error: Voltage " << voltage_check << " not found in the data. Check data or choose another voltage to do the current 2D map.\n";
            return -1;
        }
        n_ch = channel->size();
        // cout<<"Number of channels: " << n_ch << endl;
        
        std::vector<float> vect_current; // vector to store current values at given voltage
        std::vector<float> vect_current_inner; // vector to store current values at given voltage in inner channels
        std::vector<float> vect_current_border; // vector to store current values at given voltage at borders

        for (int indx = 0; indx < n_ch; ++indx) { // loop over all channels

            //------------------------------------Load all data from given channel----------------------------------------
            for (int iEntry = 0; tree_data->LoadTree(iEntry) >= 0; ++iEntry) {
                // load the data for the given tree entry
                tree_data->GetEntry(iEntry);
                // fill variables for the CV graph
                x.push_back(voltage);
                y.push_back(current->at(indx));
                yerr.push_back(current_err->at(indx));
            }

            int ch; // current channel analysing
            // which channel analysing
            if (CSIS) {
                // map CSIS channels to sensor channels
                ch = CSIS_ch_map[channel->at(indx)-1];
                // printf("Channel: %d\n", ch);
            } else {
                // use channel number as is
                ch = channel->at(indx);
                // printf("Channel: %d\n", ch);
            }
            if (ch < 1 || ch > 256) {
                std::cerr << "Error: Channel number out of range (1-256): " << ch << std::endl;
                continue; // skip this channel
            }

            // check if channel is at the border or in the inner part of the sensor
            bool is_border = std::find(border_channels.begin(), border_channels.end(), ch) != border_channels.end();

            // -------------------IV GRAPH--------------------
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

            // Fill histograms with current at given voltage
            hcurrxch->SetBinContent(indx+1, y[map_indx]); // store current at given voltage in histogram
            // hcurr->Fill(y[map_indx]); // store current at given voltage in histogram    
            // store capacitance of the right plateau of CV in 2D histogram

            vect_current.push_back(y[map_indx]); // store current at given voltage in vector
            if (is_border) {
                vect_current_border.push_back(y[map_indx]); // store current at given voltage in vector for borders
            } else {
                vect_current_inner.push_back(y[map_indx]); // store current at given voltage in vector for inner channels
            }

            hcurr_map->Fill(((ch-1)%16)+1 // X position (from 1 to 16)
                        , ((ch-1)/16)+1 // Y position
                        , y[map_indx]); // value of current at given voltage

            //             // Fill histograms with current at given voltage
            // hcurrxch->SetBinContent(indx+1, TMath::Mean(y.size(), &y[0])); // store current at given voltage in histogram
            // hcurr->Fill(TMath::Mean(y.size(), &y[0])); // store current at given voltage in histogram    
            // // store capacitance of the right plateau of CV in 2D histogram
            // hcurr_map->Fill(((ch-1)%16)+1 // X position (from 1 to 16)
            //             , ((ch-1)/16)+1 // Y position
            //             , TMath::Mean(y.size(), &y[0])); // value of current at given voltage


            dirIV->cd();
            g->Write();

            //clean data vectors
            x.clear();
            y.clear();
            yerr.clear();
        }   

        auto min_curr = *std::min_element(vect_current.begin(), vect_current.end());
        auto max_curr = *std::max_element(vect_current.begin(), vect_current.end());
        auto min_curr_border = *std::min_element(vect_current_border.begin(), vect_current_border.end());
        auto max_curr_border = *std::max_element(vect_current_border.begin(), vect_current_border.end());
        auto min_curr_inner = *std::min_element(vect_current_inner.begin(), vect_current_inner.end());
        auto max_curr_inner = *std::max_element(vect_current_inner.begin(), vect_current_inner.end());

        auto hcurr = new TH1F("hcurr", "Current distribution;Current [nA];Entries", 200, min_curr, max_curr); // histogram to store current distribution
        auto hcurr_border = new TH1F("hcurr_border", "Current distribution at borders;Current [nA];Entries", 100, min_curr_border, max_curr_border); // histogram to store current distribution at borders
        auto hcurr_inner = new TH1F("hcurr_inner", "Current distribution in inner channels;Current [nA];Entries", 200, min_curr_inner, max_curr_inner); // histogram to store current distribution in inner channels
        for (const auto& val : vect_current) {
            hcurr->Fill(val);
        }
        for (const auto& val : vect_current_border) {
            hcurr_border->Fill(val);
        }
        for (const auto& val : vect_current_inner) {
            hcurr_inner->Fill(val);
        }
        
        float curr_mean_border = hcurr_border->GetMean();
        float curr_std_border = hcurr_border->GetStdDev();
        // cout << "Current at borders mean: " << curr_mean_border << ", Current at borders std: " << curr_std_border << std::endl;
        float curr_mean_inner = hcurr_inner->GetMean();
        float curr_std_inner = hcurr_inner->GetStdDev();
        // cout << "Current in inner channels mean: " << curr_mean_inner << ", Current in inner channels std: " << curr_std_inner << std::endl;
  
        auto hcurr_zoom = new TH1F("hcurr_zoom", "Current distribution zoomed;Current [nA];Entries", 50, curr_mean_inner - 1*curr_std_inner, curr_mean_border + 1*curr_std_border); // histogram to store current distribution zoomed
        auto hcurr_border_zoom = new TH1F("hcurr_border_zoom", "Current distribution at borders zoomed;Current [nA];Entries", 25, curr_mean_border - 1*curr_std_border, curr_mean_border + 1*curr_std_border); // histogram to store current distribution at borders zoomed
        auto hcurr_inner_zoom = new TH1F("hcurr_inner_zoom", "Current distribution in inner channels zoomed;Current [nA];Entries", 25, curr_mean_inner - 1*curr_std_inner, curr_mean_inner + 1*curr_std_inner); // histogram to store current distribution in inner channels zoomed
        for (const auto& val : vect_current) {
            hcurr_zoom->Fill(val);
        }
        for (const auto& val : vect_current_border) {
            hcurr_border_zoom->Fill(val);
        }
        for (const auto& val : vect_current_inner) {
            hcurr_inner_zoom->Fill(val);
        }

        float curr_mean_border_zoom = hcurr_border_zoom->GetMean();
        float curr_std_border_zoom = hcurr_border_zoom->GetStdDev();
        cout << "Current at borders zoomed mean: " << curr_mean_border_zoom << ", Current at borders zoomed std: " << curr_std_border_zoom << std::endl;
        float curr_mean_inner_zoom = hcurr_inner_zoom->GetMean();
        float curr_std_inner_zoom = hcurr_inner_zoom->GetStdDev();
        cout << "Current in inner channels zoomed mean: " << curr_mean_inner_zoom << ", Current in inner channels zoomed std: " << curr_std_inner_zoom << std::endl;


        hcurr_map->SetMinimum(curr_mean_inner - 1*curr_std_inner);
        hcurr_map->SetMaximum(curr_mean_border + 1*curr_std_border);

        TH1F* hcurr_zoom2 = new TH1F("hcurr_zoom2", "Current distribution zoomed2;Current [nA];Entries", 50, curr_mean_inner_zoom - 1*curr_std_inner_zoom, curr_mean_border_zoom + 1*curr_std_border_zoom); // histogram to store current distribution zoomed
        TH1F* hcurr_border_zoom2 = new TH1F("hcurr_border_zoom2", "Current distribution at borders zoomed2;Current [nA];Entries", 25, curr_mean_border_zoom - 1*curr_std_border_zoom, curr_mean_border_zoom + 1*curr_std_border_zoom); // histogram to store current distribution at borders zoomed
        TH1F* hcurr_inner_zoom2 = new TH1F("hcurr_inner_zoom2", "Current distribution in inner channels zoomed2;Current [nA];Entries", 25, curr_mean_inner_zoom - 1*curr_std_inner_zoom, curr_mean_inner_zoom + 1*curr_std_inner_zoom); // histogram to store current distribution in inner channels zoomed
        for (const auto& val : vect_current) {
            hcurr_zoom2->Fill(val);
        }
        for (const auto& val : vect_current_border) {
            hcurr_border_zoom2->Fill(val);
        }
        for (const auto& val : vect_current_inner) {
            hcurr_inner_zoom2->Fill(val);
        }

        float curr_mean_border_zoom2 = hcurr_border_zoom2->GetMean();
        float curr_std_border_zoom2 = hcurr_border_zoom2->GetStdDev();
        cout << "Current at borders zoomed2 mean: " << curr_mean_border_zoom2 << ", Current at borders zoomed2 std: " << curr_std_border_zoom2 << std::endl;
        float curr_mean_inner_zoom2 = hcurr_inner_zoom2->GetMean();
        float curr_std_inner_zoom2 = hcurr_inner_zoom2->GetStdDev();
        cout << "Current in inner channels zoomed2 mean: " << curr_mean_inner_zoom2 << ", Current in inner channels zoomed2 std: " << curr_std_inner_zoom2 << std::endl;
        
        gStyle->SetPalette(kBlueRedYellow); // Set default color palette
        // Current at given voltage
        auto ccurr = new TCanvas("ccurr", "Canvas", 600, 600);
        ccurr->cd();
        // hcurr_map->SetMinimum(0.025);
        // hcurr_map->SetMaximum(0.18);
        hcurr_map->Draw("COLZ");
        hcurr_map->GetZaxis()->SetTitle("Current [nA]");
        channel_name->Draw("text same");

        myFile->cd();
        hcurrxch->Write(); // write histogram with current per channel
        hcurr->Write(); // write histogram with current distribution
        ccurr->Write(); 
        hcurr_border->Write(); // write histogram with current at borders
        hcurr_inner->Write(); // write histogram with current in inner channels
        hcurr_zoom->Write(); // write histogram with current distribution zoomed
        hcurr_border_zoom->Write(); // write histogram with current at borders zoomed
        hcurr_inner_zoom->Write(); // write histogram with current in inner channels zoomed
        hcurr_zoom2->Write(); // write histogram with current distribution zoomed2
        hcurr_border_zoom2->Write(); // write histogram with current at borders zoomed2
        hcurr_inner_zoom2->Write(); // write histogram with current in inner channels zoomed

    }

    return 0;
}