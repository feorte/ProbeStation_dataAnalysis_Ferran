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

TCanvas *draw_vdep_plot(const std::vector<double> &x, const Series &series) {
    auto *canvas = new TCanvas("c3", "Vdep", 700, 500);
    canvas->cd();

    auto *graph = new TGraphErrors(static_cast<int>(x.size()), x.data(), series.values.data(),
                                   nullptr, series.errors.data());
    graph->SetTitle("Mean depletion voltage;Sensor ID;Mean V_{dep} [V]");
    graph->SetMarkerColor(series.color);
    graph->SetLineColor(series.color);
    graph->SetMarkerStyle(series.marker);

    graph->GetYaxis()->SetRangeUser(minimum(series.values) * 0.9, maximum(series.values) * 1.1);
    graph->Draw("AP");

    auto *legend = new TLegend(0.15, 0.77, 0.45, 0.9);
    legend->AddEntry(graph, "mean values", "lep");

    const std::vector<double> high(series.values.begin(), series.values.begin() + 3);
    const std::vector<double> low(series.values.begin() + 3, series.values.end());

    const double high_mean = mean(high);
    const double high_error = standard_error_of_mean(high);
    const double low_mean = mean(low);
    const double low_error = standard_error_of_mean(low);
    const std::string high_mean_label = Form("high V_{dep} group mean = %.2f", high_mean);
    const std::string high_error_label = Form("error of mean (high V_{dep}) = %.2f", high_error);
    const std::string low_mean_label = Form("low V_{dep} group mean = %.2f", low_mean);
    const std::string low_error_label = Form("error of mean (low V_{dep}) = %.2f", low_error);

    add_mean_band(*legend, x, high_mean, high_error, series.color, 3, 0.1, high_mean_label,
                  high_error_label);
    add_mean_band(*legend, x, low_mean, low_error, series.color, 3, 0.1, low_mean_label,
                  low_error_label);

    legend->Draw("SAME");
    return canvas;
}

} // namespace

void plot_means_groups_CSIS() {
    const std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8};

    const Series I_border{
        "border",
        {0.102106, 0.124096, 0.141025, 0.0996983, 0.133513, 0.130092, 0.236907, 0.116035},
        {0.00136566, 0.00207136, 0.003234042, 0.001537522, 0.002476146, 0.003149085, 0.083538352,
         0.001503999},
        kBlue,
        20};

    const Series I_inner{
        "inner",
        {0.0277687, 0.0323819, 0.0422819, 0.0261969, 0.02996, 0.0325467, 0.0487106, 0.0498156},
        {0.00066764, 0.00116804, 0.00121416, 0.0005778, 0.0006884, 0.00111217, 0.00298892,
         0.00326956},
        kRed,
        21};

    const Series Vdep{
        "mean values",
        {50.6603, 46.1005, 49.2717, 39.4758, 39.7025, 39.5328, 39.4188, 39.7485},
        {0.881010, 1.18942, 1.47291, 0.271412, 0.464647, 0.430186, 0.499530, 0.453411},
        kMagenta + 2,
        21};

    const Series cs_border{
        "border",
        {4.99895, 4.99010, 4.96726, 4.94736, 4.94510, 4.97209, 4.89675, 4.90434},
        {0.0649411, 0.0609032, 0.0659735, 0.0754100, 0.00911386, 0.0701133, 0.0981383, 0.0557984},
        kGreen + 2,
        22};

    const Series cs_inner{
        "inner",
        {4.81578, 4.80319, 4.80479, 4.79513, 4.81693, 4.82159, 4.72190, 4.72744},
        {0.0333381, 0.0241327, 0.0317854, 0.0293461, 0.0235649, 0.0175699, 0.0228334, 0.0152061},
        kRed + 2,
        23};

    const Series ndon_border{
        "border",
        {1.64833e11, 1.50792e11, -1.0, 1.27893e11, 1.32475e11, 1.25328e11, 1.27645e11, 1.27761e11},
        {5.22334e9, 5.19116e9, 0.0, 9.34799e9, 4.97980e9, 2.97154e9, 3.72298e9, 4.40407e9},
        kBlack,
        29};

    const Series ndon_inner{
        "inner",
        {1.55153e11, 1.42815e11, -1.0, 1.18213e11, 1.19795e11, 1.19034e11, 1.16531e11, 1.16940e11},
        {2.09555e9, 2.39920e9, 0.0, 2.19398e9, 1.95682e9, 2.20134e9, 2.13760e9, 1.94701e9},
        kPink + 6,
        30};

    const Series I_border_zoom{"border",
                               {1.06852e-1, 8.28184e-2, 7.37621e-2, 1.03035e-1, 7.29869e-2,
                                1.01432e-1, 9.08203e-2, 1.04323e-1},
                               {1.37942e-2, 7.41617e-3, 5.92693e-3, 1.31972e-2, 6.08683e-3,
                                8.88252e-3, 6.13501e-3, 7.16592e-3},
                               kBlue + 2,
                               22};

    const Series I_inner_zoom{"inner",
                              {3.35711e-2, 2.19322e-2, 2.27658e-2, 2.82497e-2, 2.11307e-2,
                               2.77180e-2, 2.29835e-2, 2.55917e-2},
                              {4.56371e-3, 7.14088e-3, 2.84803e-3, 4.41030e-3, 4.63249e-3,
                               4.83072e-3, 4.72603e-3, 2.72624e-3},
                              kOrange + 7,
                              23};

    TCanvas *c1 = draw_series_plot(
        "c1", "I (nA) one zoom", "Mean current at 120 V (one zoom);Sensor ID;Mean current [nA]", x,
        {I_border, I_inner}, minimum(I_border.values, I_inner.values) * 0.9,
        maximum(I_border.values, I_inner.values) * 1.1, 0.15, 0.72, 0.45, 0.9, "%.2f");

    TCanvas *c2 = draw_series_plot(
        "c2", "I (nA) two zoom", "Mean current at 120 V;Sensor ID;Mean current [nA]", x,
        {I_border_zoom, I_inner_zoom}, 0.01,
        maximum(I_border_zoom.values, I_inner_zoom.values) * 1.5, 0.15, 0.72, 0.45, 0.9, "%.3f");

    TCanvas *c3 = draw_vdep_plot(x, Vdep);

    TCanvas *c4 = draw_series_plot(
        "c4", "C_fd (pF)", "Mean full depletion capacitance;Sensor ID;C_{fd} [pF]", x,
        {cs_border, cs_inner}, minimum(cs_border.values, cs_inner.values) * 0.9,
        maximum(cs_border.values, cs_inner.values) * 1.1, 0.15, 0.72, 0.48, 0.9, "%.3f");

    TCanvas *c5 = draw_series_plot(
        "c5", "ndon (ne/cm^3)", "Mean donor density;Sensor ID;Mean n_{don} [ne/cm^{3}]", x,
        {ndon_border, ndon_inner}, minimum(ndon_border.values, ndon_inner.values) * 0.9,
        maximum(ndon_border.values, ndon_inner.values) * 1.1, 0.55, 0.72, 0.9, 0.9, "%.2e");

    auto *results_file = new TFile("all_sensors_results_CSIS.root", "RECREATE");
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    results_file->Close();
}
