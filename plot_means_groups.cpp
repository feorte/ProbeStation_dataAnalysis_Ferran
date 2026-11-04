#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"

namespace {

struct Series {
    std::string label;
    std::vector<double> values;
    std::vector<double> errors;
    Color_t color;
    Style_t marker;
};

double mean(const std::vector<double> &values) {
    double sum = 0.0;
    for (const double value : values) {
        sum += value;
    }
    return sum / static_cast<double>(values.size());
}

double standard_error_of_mean(const std::vector<double> &values) {
    const double average = mean(values);
    double sum_of_squares = 0.0;
    for (const double value : values) {
        const double difference = value - average;
        sum_of_squares += difference * difference;
    }
    const double standard_deviation =
        std::sqrt(sum_of_squares / static_cast<double>(values.size() - 1));
    return standard_deviation / std::sqrt(static_cast<double>(values.size()));
}

double minimum(const std::vector<double> &values) {
    return *std::min_element(values.begin(), values.end());
}

double maximum(const std::vector<double> &values) {
    return *std::max_element(values.begin(), values.end());
}

double minimum(const std::vector<double> &first, const std::vector<double> &second) {
    return std::min(minimum(first), minimum(second));
}

double maximum(const std::vector<double> &first, const std::vector<double> &second) {
    return std::max(maximum(first), maximum(second));
}

std::string value_label(const std::string &text, const std::string &series_label,
                        const std::string &number_format, double value, bool include_series) {
    const std::string prefix = include_series ? text + " (" + series_label + ") = " : text + " = ";
    return Form((prefix + number_format).c_str(), value);
}

void add_mean_band(TLegend &legend, const std::vector<double> &x, double average, double error,
                   Color_t color, int line_style, double fill_alpha, const std::string &mean_label,
                   const std::string &error_label) {
    auto *mean_line = new TLine(x.front(), average, x.back(), average);
    mean_line->SetLineColor(color);
    mean_line->SetLineStyle(line_style);
    mean_line->Draw("SAME");

    auto *error_box = new TBox(x.front(), average - error, x.back(), average + error);
    error_box->SetFillColorAlpha(color, fill_alpha);
    error_box->SetLineColor(color);
    error_box->Draw("SAME");

    legend.AddEntry(mean_line, mean_label.c_str(), "l");
    legend.AddEntry(error_box, error_label.c_str(), "f");
}

TCanvas *draw_series_plot(const std::string &canvas_name, const std::string &canvas_title,
                          const std::string &graph_title, const std::vector<double> &x,
                          const std::vector<Series> &series, double y_min, double y_max,
                          double legend_x1, double legend_y1, double legend_x2, double legend_y2,
                          const std::string &number_format) {
    auto *canvas = new TCanvas(canvas_name.c_str(), canvas_title.c_str(), 700, 500);
    canvas->cd();

    std::vector<TGraphErrors *> graphs;
    graphs.reserve(series.size());

    for (const Series &item : series) {
        auto *graph = new TGraphErrors(static_cast<int>(x.size()), x.data(), item.values.data(),
                                       nullptr, item.errors.data());
        graph->SetTitle(graph_title.c_str());
        graph->SetMarkerColor(item.color);
        graph->SetLineColor(item.color);
        graph->SetMarkerStyle(item.marker);
        graphs.push_back(graph);
    }

    graphs.front()->GetYaxis()->SetRangeUser(y_min, y_max);
    graphs.front()->Draw("AP");
    for (std::size_t i = 1; i < graphs.size(); ++i) {
        graphs[i]->Draw("P SAME");
    }

    auto *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2);
    for (std::size_t i = 0; i < series.size(); ++i) {
        legend->AddEntry(graphs[i], series[i].label.c_str(), "lep");
    }

    const bool include_series = series.size() > 1;
    for (const Series &item : series) {
        const double average = mean(item.values);
        const double error = standard_error_of_mean(item.values);
        add_mean_band(
            *legend, x, average, error, item.color, 2, 0.2,
            value_label("all sensors mean", item.label, number_format, average, include_series),
            value_label("error of mean", item.label, number_format, error, include_series));
    }
    legend->Draw("SAME");
    return canvas;
}

} // namespace

void plot_means_groups() {
    const std::vector<double> x = {9, 10, 13, 14, 15, 16, 17, 18, 20};

    const Series I_border{
        "border",
        {0.102106, 0.124096, 0.141025, 0.0996983, 0.133513, 0.130092, 0.236907, 0.116035, 0.114175},
        {0.00136566, 0.00207136, 0.003234042, 0.001537522, 0.002476146, 0.003149085, 0.083538352,
         0.001503999, 0.003237171},
        kBlue,
        20};

    const Series I_inner{"inner",
                         {0.0277687, 0.0323819, 0.0422819, 0.0261969, 0.02996, 0.0325467, 0.0487106,
                          0.0498156, 0.0385838},
                         {0.00066764, 0.00116804, 0.00121416, 0.0005778, 0.0006884, 0.00111217,
                          0.00298892, 0.00326956, 0.00165396},
                         kRed,
                         21};

    const Series Vdep{
        "mean values",
        {39.2443, 42.1407, 41.5955, 39.5250, 39.2711, 40.1523, 39.2806, 39.9765, 39.9105},
        {0.518085, 0.593510, 0.350522, 0.398827, 0.590837, 0.613164, 0.511247, 0.363510, 0.557616},
        kMagenta + 2,
        21};

    const Series cs_border{
        "border",
        {4.92100, 4.93655, 4.93350, 4.94820, 4.96133, 4.98191, 4.98339, 4.96732, 4.96264},
        {0.00595503, 0.00768074, 0.00693998, 0.0141305, 0.00911386, 0.0102376, 0.0171502, 0.0132169,
         0.0169416},
        kGreen + 2,
        22};

    const Series cs_inner{
        "inner",
        {4.73487, 4.74862, 4.74626, 4.75559, 4.77339, 4.79067, 4.79506, 4.77681, 4.76937},
        {0.00579868, 0.00832303, 0.00672275, 0.00895443, 0.00985736, 0.00844902, 0.0112740,
         0.0105223, 0.0125169},
        kRed + 2,
        23};

    const Series ndon_border{"border",
                             {1.27814e11, 1.35721e11, 1.34992e11, 1.29874e11, 1.29836e11,
                              1.33602e11, 1.30404e11, 1.31473e11, 1.31337e11},
                             {1.11698e9, 1.83940e9, 2.01939e9, 1.58712e9, 1.68886e9, 1.32846e9,
                              2.23421e9, 1.97902e9, 2.21336e9},
                             kBlack,
                             29};

    const Series ndon_inner{"inner",
                            {1.16177e11, 1.24833e11, 1.23430e11, 1.15856e11, 1.17994e11, 1.20852e11,
                             1.19276e11, 1.17619e11, 1.19376e11},
                            {1.77822e9, 1.51702e9, 1.82049e9, 1.17523e9, 1.53691e9, 1.39873e9,
                             1.72993e9, 1.25877e9, 1.41285e9},
                            kPink + 6,
                            30};

    const Series I_border_zoom{"border",
                               {1.00802e-1, 1.19649e-1, 1.34433e-1, 9.61273e-2, 1.29500e-1,
                                1.12517e-1, 1.32401e-1, 1.14193e-1, 1.04858e-1},
                               {4.10529e-3, 8.23437e-3, 9.52896e-3, 4.44298e-3, 5.24834e-3,
                                4.62978e-3, 6.22953e-3, 7.28483e-3, 6.05784e-3},
                               kBlue + 2,
                               22};

    const Series I_inner_zoom{"inner",
                              {2.27959e-2, 2.69181e-2, 3.61673e-2, 2.33366e-2, 2.42778e-2,
                               2.59393e-2, 3.36645e-2, 2.96695e-2, 2.82294e-2},
                              {6.33633e-3, 1.95250e-3, 3.66714e-3, 2.34243e-3, 4.53124e-3,
                               2.06813e-3, 2.99997e-3, 3.25653e-3, 2.46318e-3},
                              kOrange + 7,
                              23};

    TCanvas *c1 = draw_series_plot(
        "c1", "I (nA) one zoom", "Mean current at 120 V (one zoom);Sensor ID;Mean current [nA]", x,
        {I_border, I_inner}, minimum(I_border.values, I_inner.values) * 0.9,
        maximum(I_border.values, I_inner.values) * 1.1, 0.15, 0.72, 0.45, 0.9, "%.2f");

    TCanvas *c2 = draw_series_plot(
        "c2", "I (nA) two zoom", "Mean current at 120 V;Sensor ID;Mean current [nA]", x,
        {I_border_zoom, I_inner_zoom}, minimum(I_border_zoom.values, I_inner_zoom.values) * 0.8,
        maximum(I_border_zoom.values, I_inner_zoom.values) * 1.2, 0.15, 0.72, 0.45, 0.9, "%.3f");

    TCanvas *c3 = draw_series_plot(
        "c3", "Vdep", "Mean depletion voltage;Sensor ID;Mean V_{dep} [V]", x, {Vdep},
        minimum(Vdep.values) * 0.9, maximum(Vdep.values) * 1.1, 0.15, 0.77, 0.45, 0.9, "%.2f");

    TCanvas *c4 = draw_series_plot(
        "c4", "C_fd (pF)", "Mean full depletion capacitance;Sensor ID;C_{fd} [pF]", x,
        {cs_border, cs_inner}, minimum(cs_border.values, cs_inner.values) * 0.9,
        maximum(cs_border.values, cs_inner.values) * 1.1, 0.15, 0.72, 0.48, 0.9, "%.3f");

    TCanvas *c5 = draw_series_plot(
        "c5", "ndon (ne/cm^3)", "Mean donor density;Sensor ID;Mean n_{don} [ne/cm^{3}]", x,
        {ndon_border, ndon_inner}, minimum(ndon_border.values, ndon_inner.values) * 0.9,
        maximum(ndon_border.values, ndon_inner.values) * 1.1, 0.55, 0.72, 0.9, 0.9, "%.2e");

    auto *results_file = new TFile("all_sensors_results.root", "RECREATE");
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    results_file->Close();
}
