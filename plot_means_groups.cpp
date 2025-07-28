#include <TCanvas.h>
#include <TGraphErrors.h>

void plot_means_groups() {
    // Channel/ID for X (from your first column, e.g., 9, 10, 13, ...)
    const int n = 9;
    double x[n]      = {9, 10, 13, 14, 15, 16, 17, 18, 20};

    // "I border" group: mean and error
    double y1[n]     = {0.102106, 0.124096, 0.141025, 0.0996983, 0.133513, 0.130092, 0.236907, 0.116035, 0.114175};
    double yerr1[n]  = {0.00946157, 0.0143508, 0.0224061, 0.0109801, 0.0169756, 0.0218175, 0.619537, 0.0107407, 0.0226602};

    // "I inner" group: mean and error
    double y2[n]     = {0.0277687, 0.0323819, 0.0422819, 0.0261969, 0.02996, 0.0325467, 0.0487106, 0.0498156, 0.0385838};
    double yerr2[n]  = {0.00883199, 0.0155836, 0.0163799, 0.00759972, 0.00915857, 0.0146705, 0.0410909, 0.0454222, 0.0217544};

    TCanvas *c1 = new TCanvas("c1","Means of I border and I inner",800,600);

    // Create TGraphErrors for each group
    TGraphErrors *g1 = new TGraphErrors(n, x, y1, 0, yerr1);
    TGraphErrors *g2 = new TGraphErrors(n, x, y2, 0, yerr2);

    g1->SetMarkerStyle(20); g1->SetMarkerColor(kBlue); g1->SetLineColor(kBlue);
    g2->SetMarkerStyle(21); g2->SetMarkerColor(kRed);  g2->SetLineColor(kRed);

    g1->SetTitle("I border and I inner Means;Channel;Current (nA)");
    g1->GetYaxis()->SetRangeUser(0, 0.3);


    g1->Draw("AP");
    g2->Draw("P SAME");
    
    // Add a legend
    auto leg = new TLegend(0.65,0.75,0.88,0.88);
    leg->AddEntry(g1,"I border","lep");
    leg->AddEntry(g2,"I inner","lep");
    leg->Draw();

    c1->Update();

////////////// More zoom
// X-labels for each group/measurement
    const int num_points = 9;

    double indices[num_points] = {9, 10, 13, 14, 15, 16, 17, 18, 20};

    // Means and errors for "I border"
    double border_avg[num_points]     = {0.100806, 0.12062, 0.135656, 0.0963891, 0.129376, 0.120814, 0.152946, 0.115018, 0.107442};
    double border_stddev[num_points]  = {0.00366852, 0.00661943, 0.00933218, 0.00410977, 0.00547003, 0.0101554, 0.0566367, 0.00532697, 0.00733312};

    // Means and errors for "I inner"
    double inner_avg[num_points]     = {0.0250042, 0.0279725, 0.0382536, 0.0238612, 0.0276031, 0.0285825, 0.0373671, 0.0359152, 0.0311854};
    double inner_stddev[num_points]  = {0.00403847, 0.00380338, 0.00593112, 0.0024266, 0.00292108, 0.00487144, 0.00969866, 0.0142474, 0.00688912};

    TCanvas *canvas_zoom = new TCanvas("canvas_zoom", "Zoomed Group Means", 800, 600);

    TGraphErrors *graph_border = new TGraphErrors(num_points, indices, border_avg, 0, border_stddev);
    TGraphErrors *graph_inner  = new TGraphErrors(num_points, indices, inner_avg, 0, inner_stddev);

    graph_border->SetMarkerStyle(22); graph_border->SetMarkerColor(kGreen+2); graph_border->SetLineColor(kGreen+3);
    graph_inner->SetMarkerStyle(26);  graph_inner->SetMarkerColor(kViolet+2);  graph_inner->SetLineColor(kMagenta+2);

    graph_border->SetTitle("Group Means (Zoomed);Point Index;Current [nA]");
    graph_border->GetYaxis()->SetRangeUser(0, 0.2);

    graph_border->Draw("AP");
    graph_inner->Draw("P SAME");

    auto legend_zoom = new TLegend(0.18, 0.75, 0.42, 0.88);
    legend_zoom->AddEntry(graph_border, "I border", "lep");
    legend_zoom->AddEntry(graph_inner, "I inner", "lep");
    legend_zoom->Draw();

    canvas_zoom->Update();
}


