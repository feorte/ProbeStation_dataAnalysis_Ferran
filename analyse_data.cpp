#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

namespace {

constexpr double kElectronCharge = 1.602176634e-19;
constexpr double kVacuumPermittivity = 8.8541878128e-12;
constexpr double kSiliconRelativePermittivity = 11.7;
constexpr double kSiliconPermittivity = kSiliconRelativePermittivity * kVacuumPermittivity;
constexpr double kPadArea = 30.25e-6;
constexpr int kSensorChannels = 256;
constexpr int kTotalChannels = 264;
constexpr int kColumns = 16;
constexpr int kRows = 16;
constexpr int kIvMapVoltage = 120;

const int kCsisChannelMap[kTotalChannels] = {
    240, 256, 16,  176, 160, 208, 192, 112, 96,  144, 128, 48,  32,  80,  64,  224, 239, 255, 15,
    175, 159, 207, 191, 111, 95,  143, 127, 47,  31,  79,  63,  223, 238, 254, 14,  174, 158, 206,
    190, 110, 94,  142, 126, 46,  30,  78,  62,  222, 237, 253, 13,  173, 157, 205, 189, 109, 93,
    141, 125, 45,  29,  77,  61,  221, 236, 252, 12,  172, 156, 204, 188, 108, 92,  140, 124, 44,
    28,  76,  60,  220, 235, 251, 11,  171, 155, 203, 187, 107, 91,  139, 123, 43,  27,  75,  59,
    219, 234, 250, 10,  170, 154, 202, 186, 106, 90,  138, 122, 42,  26,  74,  58,  218, 233, 249,
    9,   169, 153, 201, 185, 105, 89,  137, 121, 41,  25,  73,  57,  217, 232, 248, 8,   168, 152,
    200, 184, 104, 88,  136, 120, 40,  24,  72,  56,  216, 231, 247, 7,   167, 151, 199, 183, 103,
    87,  135, 119, 39,  23,  71,  55,  215, 230, 246, 6,   166, 150, 198, 182, 102, 86,  134, 118,
    38,  22,  70,  54,  214, 229, 245, 5,   165, 149, 197, 181, 101, 85,  133, 117, 37,  21,  69,
    53,  213, 228, 244, 4,   164, 148, 196, 180, 100, 84,  132, 116, 36,  20,  68,  52,  212, 227,
    243, 3,   163, 147, 195, 179, 99,  83,  131, 115, 35,  19,  67,  51,  211, 226, 242, 2,   162,
    146, 194, 178, 98,  82,  130, 114, 34,  18,  66,  50,  210, 225, 241, 1,   161, 145, 193, 177,
    97,  81,  129, 113, 33,  17,  65,  49,  209};

const std::vector<int> kBorderChannels = {
    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,  15,
    16,  17,  32,  33,  48,  49,  64,  65,  80,  81,  96,  97,  112, 113, 128,
    129, 144, 145, 160, 161, 176, 177, 192, 193, 208, 209, 224, 225, 240, 241,
    242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256};

const std::vector<int> kLeftBorderChannels = {17,  33,  49,  65,  81,  97,  113, 129,
                                              145, 161, 177, 193, 209, 225, 241};

const std::vector<int> kRightBorderChannels = {16,  32,  48,  64,  80,  96,  112, 128,
                                               144, 160, 176, 192, 208, 224, 240};

int file_type(const std::string &filename) {
    if (filename.size() < 7) {
        std::cerr << "Error: Filename is too short to determine type." << std::endl;
        return -1;
    }

    const std::string ending = filename.substr(filename.size() - 7);
    if (ending == "CV.root") {
        std::cout << "File is CV type." << std::endl;
        return 1;
    }
    if (ending == "IV.root") {
        std::cout << "File is IV type." << std::endl;
        return 0;
    }

    std::cerr << "Error: Unknown file type." << std::endl;
    return -1;
}

void choose_branches(TTree *tree, const std::vector<std::string> &branches) {
    tree->SetBranchStatus("*", false);
    for (const std::string &name : branches) {
        if (tree->GetBranch(name.c_str())) {
            tree->SetBranchStatus(name.c_str(), true);
        } else {
            std::cerr << "Warning: Branch '" << name << "' not found in the tree." << std::endl;
        }
    }
}

void log_scale(std::size_t number_of_points, const std::vector<float> &x,
               const std::vector<float> &y, const std::vector<float> &yerr,
               std::vector<float> &x_log, std::vector<float> &y_log, std::vector<float> &yerr_log) {
    for (std::size_t i = 0; i < number_of_points; ++i) {
        x_log[i] = std::log(x[i]);
        y_log[i] = std::log(y[i]);
        yerr_log[i] = yerr[i] / y[i];
    }
}

bool contains(const std::vector<int> &values, int value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

bool in_range(int value, int minimum, int maximum) {
    return value >= minimum && value <= maximum;
}

int mapped_channel(int index, const std::vector<int> *channel, bool csis) {
    const int channel_number = channel->at(index);
    return csis ? kCsisChannelMap[channel_number - 1] : channel_number;
}

void pixel_position(int channel, int &x, int &y) {
    if (in_range(channel, kSensorChannels + 1, kTotalChannels)) {
        x = kColumns + 1 + (channel - (kSensorChannels + 1));
        y = kRows;
        return;
    }

    x = ((channel - 1) % kColumns) + 1;
    y = ((channel - 1) / kColumns) + 1;
}

std::pair<double, double> value_range(const std::vector<float> &values) {
    if (values.empty()) {
        return {0.0, 1.0};
    }
    const auto range = std::minmax_element(values.begin(), values.end());
    return {*range.first, *range.second};
}

void style_graph(TGraph &graph) {
    graph.SetMarkerStyle(20);
    graph.SetMarkerSize(0.8);
    graph.SetMarkerColor(kBlue + 2);
    graph.GetXaxis()->SetTitleFont(42);
    graph.GetYaxis()->SetTitleFont(42);
    graph.GetXaxis()->SetLabelFont(42);
    graph.GetYaxis()->SetLabelFont(42);
    graph.GetXaxis()->SetTitleSize(0.05);
    graph.GetYaxis()->SetTitleSize(0.05);
    graph.GetXaxis()->SetTitleOffset(1.2);
    graph.GetYaxis()->SetTitleOffset(1.4);
    graph.GetXaxis()->CenterTitle();
    graph.GetYaxis()->CenterTitle();
}

void style_canvas(TCanvas &canvas) {
    canvas.SetGrid();
    canvas.SetTicks();
    canvas.SetLeftMargin(0.15);
    canvas.SetBottomMargin(0.15);
}

void write_channel_map(TH2F &channel_map) {
    for (int row = 1; row <= kRows; ++row) {
        for (int column = 1; column <= kColumns; ++column) {
            channel_map.Fill(column, row, column + kColumns * (row - 1));
        }
    }
}

int analyse_cv(const std::string &sensor_number, TTree *tree_data, TTree *tree_systCap) {
    std::cout << "Processing CV data...\n";

    const bool csis = true;
    const std::vector<std::string> branches = {"voltage", "channel", "cs_uncorr", "cs_err"};
    choose_branches(tree_data, branches);
    choose_branches(tree_systCap, branches);

    float voltage = 0.0f;
    std::vector<int> *channel = nullptr;
    std::vector<float> *cs = nullptr;
    std::vector<float> *cs_err = nullptr;
    std::vector<float> *cs_syst = nullptr;
    std::vector<float> *cs_syst_err = nullptr;
    tree_data->SetBranchAddress("voltage", &voltage);
    tree_data->SetBranchAddress("channel", &channel);
    tree_data->SetBranchAddress("cs_uncorr", &cs);
    tree_data->SetBranchAddress("cs_err", &cs_err);
    tree_systCap->SetBranchAddress("cs_uncorr", &cs_syst);
    tree_systCap->SetBranchAddress("cs_err", &cs_syst_err);

    auto *hVdepxch = new TH1F("hVdepxch", "Depletion Voltage per Channel;Channel;V_{dep} [V]",
                              kTotalChannels, 0.5, kTotalChannels + 0.5);
    auto *hndonxch =
        new TH1F("hndonxch", "Donnor density per Channel;Channel;Donnor density [ne/cm^{3}]",
                 kTotalChannels, 0.5, kTotalChannels + 0.5);
    auto *hcsxch =
        new TH1F("hcsxch", "Full depletion capacitance per Channel;Channel;Capacitance [pF]",
                 kTotalChannels, 0.5, kTotalChannels + 0.5);
    auto *hchi2xch = new TH1F("hchi2xch", "Global chi2 per Channel;Channel;Chi2", kTotalChannels,
                              0.5, kTotalChannels + 0.5);
    auto *hchi2 = new TH1F("hchi2", "Global chi2;Chi2;Entries", 150, 0, 2e5);
    std::vector<float> chi2_glob;

    auto *hVdep_map =
        new TH2F("hVdep_map", Form("Depletion Voltage for sensor %s;X;Y", sensor_number.c_str()),
                 kColumns, 0.5, kColumns + 0.5, kRows, 0.5, kRows + 0.5);
    auto *hndon_map =
        new TH2F("hndon_map", Form("Donnor density for sensor %s;X;Y", sensor_number.c_str()),
                 kColumns, 0.5, kColumns + 0.5, kRows, 0.5, kRows + 0.5);
    auto *hcs_plateau_map =
        new TH2F("hcs_plateau_map",
                 Form("Full depletion capacitance for sensor %s;X;Y", sensor_number.c_str()),
                 kColumns, 0.5, kColumns + 0.5, kRows, 0.5, kRows + 0.5);
    auto *hchi2sus = new TH2F("hchi2sus", "Weird chi2;X;Y", kColumns, 0.5, kColumns + 0.5, kRows,
                              0.5, kRows + 0.5);

    auto *channel_name = new TH2F("channel_name", "Sensor pixels;X;Y", kColumns, 0.5,
                                  kColumns + 0.5, kRows, 0.5, kRows + 0.5);
    write_channel_map(*channel_name);

    std::filesystem::create_directories("results");
    const std::string results_file = "results/CV_" + sensor_number + ".root";
    auto myFile = std::make_unique<TFile>(results_file.c_str(), "RECREATE");
    if (!myFile || myFile->IsZombie()) {
        std::cerr << "Error: Cannot create results file: " << results_file << std::endl;
        return -1;
    }

    TDirectory *dirCV = myFile->mkdir("CV_graphs");
    TDirectory *dirDepletion = myFile->mkdir("Depletion_voltage");
    TDirectory *dirDonnor = myFile->mkdir("Donnor_density");

    int n_volt = 0;
    for (int iEntry = 0; tree_data->LoadTree(iEntry) >= 0; ++iEntry) {
        tree_data->GetEntry(iEntry);
        n_volt = iEntry + 1;
    }
    const int n_ch = static_cast<int>(channel->size());

    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> yerr;
    float y_syst = 0.0f;
    float yerr_syst = 0.0f;

    std::vector<float> vect_vdep;
    std::vector<float> vect_donor_density;
    std::vector<float> vect_donor_density_inner;
    std::vector<float> vect_donor_density_border;
    std::vector<float> vect_capacitance;
    std::vector<float> vect_capacitance_inner;
    std::vector<float> vect_capacitance_border;
    std::vector<float> vect_capacitance_uncorrected;
    std::vector<float> vect_donor_density_uncorrected;
    std::vector<int> invalid_ch;

    for (int indx = 0; indx < n_ch; ++indx) {
        x.clear();
        y.clear();
        yerr.clear();

        for (int iEntry = 0; tree_systCap->LoadTree(iEntry) >= 0; ++iEntry) {
            tree_systCap->GetEntry(iEntry);
            y_syst = cs_syst->at(indx);
            yerr_syst = cs_syst_err->at(indx);
        }

        for (int iEntry = 0; tree_data->LoadTree(iEntry) >= 0; ++iEntry) {
            tree_data->GetEntry(iEntry);
            x.push_back(voltage);
            y.push_back(cs->at(indx));
            yerr.push_back(cs_err->at(indx));
        }

        const int ch = mapped_channel(indx, channel, csis);
        if (!in_range(ch, 1, kTotalChannels)) {
            std::cerr << "Error: Channel number out of range (1-264): " << ch << std::endl;
            continue;
        }

        if (sensor_number == "15") {
            y_syst += 0.6046573f;
        }
        if (sensor_number == "16") {
            y_syst += 0.6010227f;
        }
        if (sensor_number == "17") {
            y_syst += 0.5761993f;
        }
        if (sensor_number == "18") {
            y_syst += 0.5857236f;
        }

        for (int i = 0; i < n_volt; ++i) {
            y[i] -= y_syst;
            yerr[i] = std::sqrt(yerr[i] * yerr[i] + yerr_syst * yerr_syst);
        }

        std::vector<float> xerr(n_volt, 0.0f);
        auto *g = new TGraphErrors(n_volt, x.data(), y.data(), xerr.data(), yerr.data());
        g->SetName(Form("Channel_%d", ch));
        g->SetTitle(Form("Channel %d;Voltage [V]; Capacitance [pF]", ch));
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.75);
        g->SetMarkerColor(kBlue);
        g->GetXaxis()->CenterTitle();
        g->GetYaxis()->CenterTitle();

        std::vector<float> x_log(n_volt);
        std::vector<float> y_log(n_volt);
        std::vector<float> yerr_log(n_volt);
        log_scale(n_volt, x, y, yerr, x_log, y_log, yerr_log);

        auto *glog =
            new TGraphErrors(n_volt, x_log.data(), y_log.data(), xerr.data(), yerr_log.data());
        glog->SetName(Form("Channel %d log scale", ch));
        glog->SetTitle(Form("Channel %d log scale;ln V; ln C", ch));
        glog->SetMarkerStyle(20);
        glog->SetMarkerSize(0.75);
        glog->SetMarkerColor(kBlue);

        auto *c1 = new TCanvas(Form("Channel_%d", ch), Form("Channel_%d", ch), 800, 600);
        style_canvas(*c1);
        style_graph(*glog);
        glog->Draw("AP");

        auto *lfit = new TF1("lfit", "pol1", x_log[1], std::log(35));
        glog->Fit(lfit, "EQR0");
        lfit->SetLineColor(kGreen + 3);
        lfit->SetLineWidth(2);
        lfit->SetLineStyle(7);

        auto *rfit = new TF1("rfit", "pol0", std::log(60), x_log[n_volt - 1]);
        glog->Fit(rfit, "EQR0");
        rfit->SetLineColor(kRed);
        rfit->SetLineWidth(2);
        rfit->SetLineStyle(7);

        const double p0_1 = lfit->GetParameter(0);
        const double p1_1 = lfit->GetParameter(1);
        const double p0_2 = rfit->GetParameter(0);
        const double p1_2 = 0.0;
        const double log_Vdep = (p0_2 - p0_1) / (p1_1 - p1_2);
        const double V_dep = std::exp(log_Vdep);

        lfit->SetRange(x_log[1], log_Vdep + 1);
        rfit->SetRange(log_Vdep - 1, x_log[n_volt - 1]);
        lfit->Draw("SAME");
        rfit->Draw("SAME");

        const double chi2_lfit = lfit->GetChisquare();
        const double chi2_rfit = rfit->GetChisquare();
        chi2_glob.push_back(chi2_lfit + chi2_rfit);

        if (std::isfinite(V_dep)) {
            hchi2xch->SetBinContent(indx + 1, chi2_lfit + chi2_rfit);
        } else {
            std::cerr << "Warning: Skipping invalid chi2_glob for channel index " << indx + 1
                      << std::endl;
        }
        hchi2->Fill(chi2_lfit + chi2_rfit);

        int x_pos = 0;
        int y_pos = 0;
        pixel_position(ch, x_pos, y_pos);

        bool is_invalid_channel =
            (ch == 3 || ch == 77 || ch == 78 || ch == 164) || (y.back() < 0 || y[0] < 4);
        if (csis) {
            is_invalid_channel = (y.back() < 0 || y[0] < 4);
        }
        if (sensor_number == "20") {
            is_invalid_channel = is_invalid_channel || (ch == 37 || ch == 80 || ch == 115 ||
                                                        ch == 162 || ch == 170 || ch == 184);
        }

        if (is_invalid_channel && in_range(ch, 1, kSensorChannels)) {
            invalid_ch.push_back(ch);
            if (ch - kColumns >= 1) {
                invalid_ch.push_back(ch - kColumns);
            }
            if (ch - 1 >= 1 && !contains(kLeftBorderChannels, ch)) {
                invalid_ch.push_back(ch - 1);
            }
            if (ch + 1 <= kSensorChannels && !contains(kRightBorderChannels, ch)) {
                invalid_ch.push_back(ch + 1);
            }
            if (ch + kColumns <= kSensorChannels) {
                invalid_ch.push_back(ch + kColumns);
            }
        }

        vect_capacitance_uncorrected.push_back(std::exp(p0_2));

        if ((sensor_number == "17" && (ch == 210 || ch == 226)) ||
            (sensor_number == "18" && (ch == 120 || ch == 136))) {
            is_invalid_channel = true;
        }

        float Vdep_value = -1.0f;
        if (in_range(ch, 1, kSensorChannels)) {
            if (!is_invalid_channel) {
                Vdep_value = static_cast<float>(V_dep);
                vect_vdep.push_back(Vdep_value);
            }
            hVdep_map->Fill(x_pos, y_pos, Vdep_value);
        }

        if (std::isfinite(V_dep)) {
            hVdepxch->SetBinContent(indx + 1, V_dep);
        } else {
            std::cerr << "Warning: Skipping invalid V_dep for channel index " << indx + 1
                      << std::endl;
        }

        auto *line =
            new TLine(log_Vdep, glog->GetYaxis()->GetXmin(), log_Vdep, glog->GetYaxis()->GetXmax());
        line->SetLineColor(kMagenta + 2);
        line->SetLineStyle(9);
        line->SetLineWidth(2);
        line->Draw("SAME");

        auto *legend = new TLegend(0.18, 0.18, 0.5, 0.31);
        legend->SetTextFont(42);
        legend->SetTextSize(0.03);
        legend->AddEntry(glog, "Data (ln C vs ln V)", "p");
        legend->AddEntry(lfit, "Rising fit", "l");
        legend->AddEntry(rfit, Form("Plateau fit = %.2f pF", std::exp(p0_2)), "l");
        legend->AddEntry(line, Form("V_{dep} = %.2f V", V_dep), "l");
        legend->Draw("SAME");

        std::vector<float> y_new(n_volt);
        std::vector<float> y_new_err(n_volt);
        for (int i = 0; i < n_volt; ++i) {
            y_new[i] = std::pow(y[i], -2);
            y_new_err[i] = 2 * yerr[i] / std::pow(y[i], 3);
        }

        auto *gnew =
            new TGraphErrors(n_volt, x.data(), y_new.data(), xerr.data(), y_new_err.data());
        gnew->SetName(Form("Channel %d", ch));
        gnew->SetTitle(Form("Channel %d;Voltage [V]; 1/C^{2} [1/pF^{2}]", ch));

        dirDonnor->cd();
        auto *c2 = new TCanvas(Form("donn_Channel_%d", ch), Form("Channel_%d", ch), 800, 600);
        style_canvas(*c2);
        style_graph(*gnew);
        gnew->Draw("AP");

        gnew->Fit("pol1", "Q0", "", x[1], 35);
        auto *don_fit = static_cast<TF1 *>(gnew->GetFunction("pol1")->Clone("don_fit"));
        don_fit->SetLineColor(kRed);
        don_fit->SetLineWidth(2);
        don_fit->SetLineStyle(2);
        don_fit->SetRange(x[0], x[8]);
        don_fit->Draw("SAME");

        double p1_don = don_fit->GetParameter(1);
        p1_don = p1_don * std::pow(10, 24);
        const double donor_density =
            (2) / (kElectronCharge * kSiliconPermittivity * std::pow(kPadArea, 2) * p1_don) * 1e-6;
        vect_donor_density_uncorrected.push_back(static_cast<float>(donor_density));

        dirCV->cd();
        g->Write();

        dirDepletion->cd();
        c1->Write();

        dirDonnor->cd();
        c2->Write();
    }

    if (sensor_number == "17") {
        invalid_ch.push_back(210);
        invalid_ch.push_back(226);
    }
    if (sensor_number == "18") {
        invalid_ch.push_back(120);
        invalid_ch.push_back(136);
    }

    float cs_value = -1.0f;
    float ndon_value = -1.0f;
    for (int i = 0; i < n_ch; ++i) {
        const int ch = mapped_channel(i, channel, csis);
        if (contains(invalid_ch, ch)) {
            cs_value = -1.0f;
            ndon_value = -1.0f;
        } else {
            const bool is_border = contains(kBorderChannels, ch);
            cs_value = vect_capacitance_uncorrected[i];
            ndon_value = vect_donor_density_uncorrected[i];
            if (in_range(ch, 1, kSensorChannels)) {
                vect_capacitance.push_back(cs_value);
                vect_donor_density.push_back(ndon_value);
                if (is_border) {
                    vect_capacitance_border.push_back(cs_value);
                    vect_donor_density_border.push_back(ndon_value);
                } else {
                    vect_capacitance_inner.push_back(cs_value);
                    vect_donor_density_inner.push_back(ndon_value);
                }
            }
        }

        if (in_range(ch, 1, kSensorChannels)) {
            int x_pos = 0;
            int y_pos = 0;
            pixel_position(ch, x_pos, y_pos);
            hcs_plateau_map->Fill(x_pos, y_pos, cs_value);
            hndon_map->Fill(x_pos, y_pos, ndon_value);
        }

        hcsxch->SetBinContent(ch, cs_value);
        hndonxch->SetBinContent(ch, ndon_value);
    }

    const auto vdep_range = value_range(vect_vdep);
    const auto ndon_range = value_range(vect_donor_density);
    const auto ndon_border_range = value_range(vect_donor_density_border);
    const auto ndon_inner_range = value_range(vect_donor_density_inner);
    const auto cs_range = value_range(vect_capacitance);
    const auto cs_border_range = value_range(vect_capacitance_border);
    const auto cs_inner_range = value_range(vect_capacitance_inner);

    auto *hVdep = new TH1F("hVdep", "Depletion Voltage;V_{dep} [V];Entries", 30, vdep_range.first,
                           vdep_range.second);
    auto *hndon = new TH1F("hndon", "Donnor density;Donnor density [ne/cm^{3}];Entries", 30,
                           ndon_range.first - 1e9, 1.5e11);
    auto *hndon_border =
        new TH1F("hndon_border", "Donnor density at borders;Donnor density [ne/cm^{3}];Entries",
                 200, ndon_border_range.first, ndon_border_range.second);
    auto *hndon_inner = new TH1F(
        "hndon_inner", "Donnor density in inner channels;Donnor density [ne/cm^{3}];Entries", 200,
        ndon_inner_range.first, ndon_inner_range.second);
    auto *hcs = new TH1F("hcs", "Full depletion capacitance;Capacitance [pF];Entries", 20,
                         cs_range.first - 0.05, 5.2);
    auto *hcs_border =
        new TH1F("hcs_border", "Full depletion capacitance at borders;Capacitance [pF];Entries",
                 200, cs_border_range.first, cs_border_range.second);
    auto *hcs_inner = new TH1F(
        "hcs_inner", "Full depletion capacitance in inner channels;Capacitance [pF];Entries", 200,
        cs_inner_range.first, cs_inner_range.second);

    for (const float value : vect_vdep) {
        hVdep->Fill(value);
    }
    for (const float value : vect_donor_density) {
        hndon->Fill(value);
    }
    for (const float value : vect_donor_density_border) {
        hndon_border->Fill(value);
    }
    for (const float value : vect_donor_density_inner) {
        hndon_inner->Fill(value);
    }
    for (const float value : vect_capacitance) {
        hcs->Fill(value);
    }
    for (const float value : vect_capacitance_border) {
        hcs_border->Fill(value);
    }
    for (const float value : vect_capacitance_inner) {
        hcs_inner->Fill(value);
    }

    auto *gaus_vdep = new TF1("gaus_vdep", "gaus", vdep_range.first, 40.5);
    hVdep->Fit(gaus_vdep, "REM");

    auto *gaus_hcs_inner = new TF1("gaus_hcs_inner", "gaus", 4.66, 4.76);
    auto *gaus_hcs_border = new TF1("gaus_hcs_border", "gaus", 4.85, 5.02);
    hcs->Fit(gaus_hcs_inner, "REM");
    hcs->Fit(gaus_hcs_border, "REM+");

    auto *gaus_hndon_inner = new TF1("gaus_hndon_inner", "gaus", 1.1e11, 1.2e11);
    auto *gaus_hndon_border = new TF1("gaus_hndon_border", "gaus", 1.25e11, 1.35e11);
    hndon->Fit(gaus_hndon_inner, "REM");
    hndon->Fit(gaus_hndon_border, "REM+");

    const float Vdep_mean = hVdep->GetMean();
    const float Vdep_std = hVdep->GetStdDev();
    const float Vdep_sum_of_weights = hVdep->GetSumOfWeights();
    std::cout << "Vdep mean: " << Vdep_mean << ", Vdep std: " << Vdep_std
              << ", Vdep sum of weights: " << Vdep_sum_of_weights << std::endl;

    if (sensor_number == "16" || sensor_number == "17" || sensor_number == "18") {
        std::cout << "last sensors" << std::endl;
        hVdep_map->SetMinimum(vdep_range.first);
        hVdep_map->SetMaximum(vdep_range.second);
    } else if (sensor_number == "20") {
        std::cout << "sensor 20" << std::endl;
        hVdep_map->SetMinimum(vdep_range.first - 0.5);
        hVdep_map->SetMaximum(vdep_range.second + 0.5);
    } else {
        hVdep_map->SetMinimum(vdep_range.first - 0.5);
        hVdep_map->SetMaximum(vdep_range.second + 0.5);
    }

    hndon_map->SetMinimum(ndon_range.first);
    hndon_map->SetMaximum(ndon_range.second);
    hcs_plateau_map->SetMinimum(cs_range.first);
    hcs_plateau_map->SetMaximum(cs_range.second);

    const float chi2_mean = hchi2->GetMean();
    const float chi2_std = hchi2->GetStdDev();
    std::cout << "Chi2 mean: " << chi2_mean << ", Chi2 std: " << chi2_std << std::endl;

    for (std::size_t i = 0; i < chi2_glob.size(); ++i) {
        const int ch = mapped_channel(static_cast<int>(i), channel, csis);
        if (!in_range(ch, 1, kSensorChannels)) {
            continue;
        }
        if (chi2_glob[i] < 10) {
            int x_pos = 0;
            int y_pos = 0;
            pixel_position(ch, x_pos, y_pos);
            hchi2sus->Fill(x_pos, y_pos, 1);
        }
    }

    gStyle->SetPalette(kBlueRedYellow);

    auto *cVdep = new TCanvas("cVdep", "Canvas", 600, 600);
    cVdep->cd();
    hVdep_map->Draw("COLZ");
    hVdep_map->GetZaxis()->SetTitle("Vdep [V]");
    channel_name->Draw("text same");

    auto *cndon = new TCanvas("cndon_map", "Canvas", 600, 600);
    cndon->cd();
    hndon_map->Draw("COLZ");
    hndon_map->GetZaxis()->SetTitle("n_{don} [ne/cm^{3}]");
    channel_name->Draw("text same");

    auto *ccsplat = new TCanvas("ccsplat", "Canvas", 600, 600);
    ccsplat->cd();
    hcs_plateau_map->Draw("COLZ");
    hcs_plateau_map->GetZaxis()->SetTitle("Capacitance [pF]");
    channel_name->Draw("text same");

    myFile->cd();
    hVdepxch->Write();
    hVdep->Write();
    hndonxch->Write();
    hndon->Write();
    hcsxch->Write();
    hcs->Write();
    cVdep->Write();
    cndon->Write();
    ccsplat->Write();

    double sum_ndon = 0.0;
    int count = 0;
    for (int ch = kSensorChannels + 1; ch <= kTotalChannels; ++ch) {
        const double ndon = hndonxch->GetBinContent(ch);
        if (ndon > 0) {
            sum_ndon += ndon;
            ++count;
        }
    }
    if (count > 0) {
        std::cout << "Mean donor density (channels 257-264): " << sum_ndon / count << " ne/cm^3"
                  << std::endl;
    } else {
        std::cerr << "No valid data found for channels 257-264." << std::endl;
    }

    double sum_cs_high_voltage = 0.0;
    int count_cs_high_voltage = 0;
    for (int ch_index = 0; ch_index < n_ch; ++ch_index) {
        const int ch = mapped_channel(ch_index, channel, csis);
        if (!in_range(ch, kSensorChannels + 1, kTotalChannels)) {
            continue;
        }

        for (int iEntry = 0; tree_data->LoadTree(iEntry) >= 0; ++iEntry) {
            tree_data->GetEntry(iEntry);
            if (voltage >= 80 && voltage <= 200) {
                const float capacitance = cs->at(ch_index) - cs_syst->at(ch_index);
                sum_cs_high_voltage += capacitance;
                ++count_cs_high_voltage;
            }
        }
    }
    if (count_cs_high_voltage > 0) {
        std::cout << "Mean capacitance (80-200 V) for channels 257-264: "
                  << sum_cs_high_voltage / count_cs_high_voltage << " pF" << std::endl;
    } else {
        std::cerr << "No valid capacitance data found for channels 257-264 in 80-200 V range."
                  << std::endl;
    }

    return 0;
}

int analyse_iv(const std::string &sensor_number, TTree *tree_data) {
    std::cout << "Processing IV data...\n";

    const bool csis = true;
    const std::vector<std::string> branches = {"voltage", "channel", "current", "current_err"};
    choose_branches(tree_data, branches);

    float voltage = 0.0f;
    std::vector<int> *channel = nullptr;
    std::vector<float> *current = nullptr;
    std::vector<float> *current_err = nullptr;
    tree_data->SetBranchAddress("voltage", &voltage);
    tree_data->SetBranchAddress("channel", &channel);
    tree_data->SetBranchAddress("current", &current);
    tree_data->SetBranchAddress("current_err", &current_err);

    auto *hcurrxch = new TH1F(
        "hcurrxch", Form("Current per Channel at %d V;Channel;Current [nA]", kIvMapVoltage),
        kTotalChannels, 0.5, kTotalChannels + 0.5);
    auto *hcurr_map =
        new TH2F("hcurr_map",
                 Form("Current at %d V for sensor %s;X;Y", kIvMapVoltage, sensor_number.c_str()),
                 kColumns, 0.5, kColumns + 0.5, kRows, 0.5, kRows + 0.5);

    auto *channel_name = new TH2F("channel_name", "Sensor pixels;X;Y", kColumns, 0.5,
                                  kColumns + 0.5, kRows, 0.5, kRows + 0.5);
    write_channel_map(*channel_name);

    std::filesystem::create_directories("results");
    const std::string results_file = "results/IV_" + sensor_number + ".root";
    auto myFile = std::make_unique<TFile>(results_file.c_str(), "RECREATE");
    if (!myFile || myFile->IsZombie()) {
        std::cerr << "Error: Cannot create results file: " << results_file << std::endl;
        return -1;
    }

    TDirectory *dirIV = myFile->mkdir("IV_graphs");

    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> yerr;

    int n_volt = 0;
    int map_indx = -1;
    for (int iEntry = 0; tree_data->LoadTree(iEntry) >= 0; ++iEntry) {
        tree_data->GetEntry(iEntry);
        if (voltage == kIvMapVoltage) {
            map_indx = iEntry;
        }
        n_volt = iEntry + 1;
    }
    if (map_indx == -1) {
        std::cerr << "Error: Voltage " << kIvMapVoltage
                  << " not found in the data. Check data or choose another voltage to do the "
                     "current 2D map.\n";
        return -1;
    }
    const int n_ch = static_cast<int>(channel->size());

    std::vector<float> vect_current;
    std::vector<float> vect_current_inner;
    std::vector<float> vect_current_border;

    for (int indx = 0; indx < n_ch; ++indx) {
        x.clear();
        y.clear();
        yerr.clear();

        for (int iEntry = 0; tree_data->LoadTree(iEntry) >= 0; ++iEntry) {
            tree_data->GetEntry(iEntry);
            x.push_back(voltage);
            y.push_back(current->at(indx));
            yerr.push_back(current_err->at(indx));
        }

        const int ch = mapped_channel(indx, channel, csis);
        if (!in_range(ch, 1, kTotalChannels)) {
            std::cerr << "Error: Channel number out of range (1-264): " << ch << std::endl;
            continue;
        }

        const bool is_border = contains(kBorderChannels, ch);

        std::vector<float> xerr(n_volt, 0.0f);
        auto *g = new TGraphErrors(n_volt, x.data(), y.data(), xerr.data(), yerr.data());
        g->SetName(Form("Channel %d", ch));
        g->SetTitle(Form("Channel %d;Voltage [V]; Current [nA]", ch));
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.75);
        g->SetMarkerColor(kBlue);
        g->GetXaxis()->CenterTitle();
        g->GetYaxis()->CenterTitle();

        int x_pos = 0;
        int y_pos = 0;
        pixel_position(ch, x_pos, y_pos);

        float current_value = -1.0f;

        bool is_invalid_channel = (ch == 3 || ch == 77 || ch == 78 || ch == 164) ||
                                  (y[map_indx] > 0.15 && is_border) ||
                                  (y[map_indx] > 0.06 && !is_border) || (y[map_indx] < 0);

        if (csis) {
            is_invalid_channel = (y[map_indx] > 0.15 && is_border) ||
                                 (y[map_indx] > 0.06 && !is_border) || (y[map_indx] < 0);
        }

        if (ch == 1 || ch == 16 || ch == 241 || ch == 256) {
            is_invalid_channel = (y[map_indx] > 0.3) || (y[map_indx] < 0);
        }

        if (in_range(ch, 1, kSensorChannels)) {
            if (!is_invalid_channel) {
                current_value = y[map_indx];
                vect_current.push_back(current_value);
                if (is_border) {
                    vect_current_border.push_back(current_value);
                } else {
                    vect_current_inner.push_back(current_value);
                }
            }
            hcurr_map->Fill(x_pos, y_pos, current_value);
        }

        hcurrxch->SetBinContent(indx + 1, current_value);

        dirIV->cd();
        g->Write();
    }

    const auto current_range = value_range(vect_current);
    const auto current_border_range = value_range(vect_current_border);
    const auto current_inner_range = value_range(vect_current_inner);

    auto *hcurr = new TH1F("hcurr", "Current distribution;Current [nA];Entries", 20,
                           current_range.first - 0.02, current_range.second + 0.01);
    auto *hcurr_border =
        new TH1F("hcurr_border", "Current distribution at borders;Current [nA];Entries", 100,
                 current_border_range.first, current_border_range.second);
    auto *hcurr_inner =
        new TH1F("hcurr_inner", "Current distribution in inner channels;Current [nA];Entries", 200,
                 current_inner_range.first, current_inner_range.second);
    for (const float value : vect_current) {
        hcurr->Fill(value);
    }
    for (const float value : vect_current_border) {
        hcurr_border->Fill(value);
    }
    for (const float value : vect_current_inner) {
        hcurr_inner->Fill(value);
    }

    auto *gaus_hcurr_inner = new TF1("gaus_hcurr_inner", "gaus", 0.01, 0.043);
    auto *gaus_hcurr_border = new TF1("gaus_hcurr_border", "gaus", 0.06, 0.084);
    hcurr->Fit(gaus_hcurr_inner, "REM");
    hcurr->Fit(gaus_hcurr_border, "REM+");

    hcurr_map->SetMinimum(current_range.first);
    hcurr_map->SetMaximum(current_range.second);

    gStyle->SetPalette(kBlueRedYellow);

    auto *ccurr = new TCanvas("ccurr", "Canvas", 600, 600);
    ccurr->cd();
    hcurr_map->Draw("COLZ");
    hcurr_map->GetZaxis()->SetTitle("Current [nA]");
    channel_name->Draw("text same");

    myFile->cd();
    hcurrxch->Write();
    hcurr->Write();
    ccurr->Write();

    return 0;
}

} // namespace

int analyse_data(std::string number_sensor = "CSIS_001", std::string type = "CV") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    std::cout << "Analysing data for sensor: " << number_sensor << ", type: " << type << std::endl;
    const std::string storing_file =
        "stored_data/stored_data_" + number_sensor + "_" + type + ".root";

    auto file = std::unique_ptr<TFile>(TFile::Open(storing_file.c_str()));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file with data: " << storing_file << std::endl;
        return -1;
    }

    auto *tree_data = file->Get<TTree>("analysis");
    if (!tree_data) {
        std::cerr << "Error: Tree 'analysis' not found in the data file." << std::endl;
        return -1;
    }

    const int is_cv = file_type(storing_file);
    if (is_cv == -1) {
        return -1;
    }

    if (is_cv == 1) {
        auto system_capacitance_file =
            std::unique_ptr<TFile>(TFile::Open("stored_data/systemCapacitanceCV.root"));
        if (!system_capacitance_file || system_capacitance_file->IsZombie()) {
            std::cerr << "Error: Cannot open ROOT file with system capacitance." << std::endl;
            return -1;
        }

        auto *tree_system_capacitance = system_capacitance_file->Get<TTree>("analysis");
        if (!tree_system_capacitance) {
            std::cerr << "Error: Tree 'analysis' not found in the system capacitance file."
                      << std::endl;
            return -1;
        }

        return analyse_cv(number_sensor, tree_data, tree_system_capacitance);
    }

    return analyse_iv(number_sensor, tree_data);
}
