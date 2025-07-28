#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLine.h>
#include <iostream>
#include <vector>

void plot_means_groups() {
    const int n = 9; // 9 points (channels)
    double x[n] = {9, 10, 13, 14, 15, 16, 17, 18, 20};

    // 1) I (nA) FIRST "ONE ZOOM" block
    double I_border[n] = {0.102106, 0.124096, 0.141025, 0.0996983, 0.133513, 0.130092, 0.236907, 0.116035, 0.114175};
    double I_border_std[n] = {
        0.00136566, 0.00207136, 0.003234042, 0.001537522,
        0.002476146, 0.003149085, 0.083538352, 0.001503999, 0.003237171
    };

    double I_inner[n] = {0.0277687, 0.0323819, 0.0422819, 0.0261969, 0.02996, 0.0325467, 0.0487106, 0.0498156, 0.0385838};
    double I_inner_std[n] = {
        0.00066764, 0.00116804, 0.00121416, 0.0005778,
        0.0006884, 0.00111217, 0.00298892, 0.00326956, 0.00165396
    };

    // 2) Vdep
    double Vdep[n] = {39.2555, 43.04, 42.603, 39.5754, 38.8523, 39.66, 38.82, 40.06, 40.06};
    double Vdep_std[n] = {
        0.03780501, 0.060241389, 0.086240722, 0.030549367,
        0.059404365, 0.040913548, 0.034136491, 0.035638625, 0.054092008
    };

    // 3) cs (pF)
    double cs_border[n] = {4.94972, 4.97321, 4.96237, 4.97211, 5.61206, 5.62857, 5.5846, 5.58511, 5.00635};
    double cs_border_std[n] = {
        0.00903608, 0.01114395, 0.00900834, 0.00900115,
        0.01250869, 0.01380738, 0.00893851, 0.00988191, 0.01301919
    };
    double cs_inner[n] = {4.76043, 4.77507, 4.7658, 4.77434, 5.3996, 5.41376, 5.38917, 5.38202, 4.81161};
    double cs_inner_std[n] = {
        0.00537691, 0.00549771, 0.00492559, 0.00488354,
        0.00518196, 0.00521576, 0.00493167, 0.00503737, 0.0070233
    };

    // 4) ndon (ne/cm^3)
    double ndon_border[n] = {
        128593e6, 137073e6, 136001e6, 131361e6, 168669e6,
        171926e6, 164653e6, 173880e6, 133589e6
    };
    double ndon_border_std[n] = {
        544838973.3, 614812157.8, 551941338.3, 689705421.5,
        973215776.7, 985635406.7, 612865095.2, 912180268.2, 756498980.9
    };
    double ndon_inner[n] = {
        117585e6, 132229e6, 128448e6, 116860e6, 152679e6,
        155203e6, 151460e6, 156954e6, 122504e6
    };
    double ndon_inner_std[n] = {
        308947799, 483042437, 398802782, 308230863,
        329007832, 324315428, 310813613, 387044398, 392019059
    };

    // 5) "TWO ZOOMS": I (nA) second measurement (zoomed)
    double I_border_zoom[n] = {0.100806, 0.12062, 0.135656, 0.0963891, 0.129376, 0.120814, 0.152946, 0.115018, 0.107442};
    double I_border_zoom_std[n] = {
        0.00059511, 0.00099792, 0.00145744, 0.00060595,
        0.00082464, 0.00162617, 0.00770728, 0.00082197, 0.00110551
    };
    double I_inner_zoom[n] = {0.0250042, 0.0279725, 0.0382536, 0.0238612, 0.0276031, 0.0285825, 0.0373671, 0.0359152, 0.0311854};
    double I_inner_zoom_std[n] = {
        0.00032649, 0.00030258, 0.00047037, 0.00019491,
        0.00023539, 0.00039255, 0.00074827, 0.00109921, 0.00057211
    };

    // Lambda for mean
    auto mean = [](const double* arr, int length) {
        double sum = 0;
        for (int i = 0; i < length; i++) sum += arr[i];
        return sum / length;
    };

    // Helper: calculate mean over a subset of indices
    auto mean_subset = [&](const double* arr, const std::vector<int>& indices) {
        double sum = 0;
        for (auto idx : indices) sum += arr[idx];
        return sum / indices.size();
    };

    // Create canvases
    TCanvas *c1 = new TCanvas("c1", "I (nA) one zoom", 700, 500);
    TCanvas *c2 = new TCanvas("c2", "I (nA) two zoom", 700, 500);
    TCanvas *c3 = new TCanvas("c3", "Vdep", 700, 500);
    TCanvas *c4 = new TCanvas("c4", "C_fd (pF)", 700, 500);
    TCanvas *c5 = new TCanvas("c5", "ndon (ne/cm^3)", 700, 500);

    // 1. I one zoom
    c1->cd();
    TGraphErrors *gI_border = new TGraphErrors(n, x, I_border, 0, I_border_std);
    TGraphErrors *gI_inner  = new TGraphErrors(n, x, I_inner, 0, I_inner_std);
    gI_border->SetMarkerColor(kBlue);     gI_border->SetLineColor(kBlue); gI_border->SetMarkerStyle(20);
    gI_inner->SetMarkerColor(kRed);       gI_inner->SetLineColor(kRed);   gI_inner->SetMarkerStyle(21);
    gI_border->SetTitle("Mean current at 120 V (one zoom);Channel;I [nA]");

    // Adjust y-axis range to include all points + margin
    double I_min = std::min(*std::min_element(I_border, I_border + n), *std::min_element(I_inner, I_inner + n));
    double I_max = std::max(*std::max_element(I_border, I_border + n), *std::max_element(I_inner, I_inner + n));
    gI_border->GetYaxis()->SetRangeUser(I_min * 0.9, I_max * 1.1);

    gI_border->Draw("AP");
    gI_inner->Draw("P SAME");
    double m1b = mean(I_border, n), m1i = mean(I_inner, n);
    TLine *l1b = new TLine(x[0], m1b, x[n - 1], m1b);
    l1b->SetLineColor(kBlue);
    l1b->SetLineStyle(2);
    l1b->Draw();
    TLine *l1i = new TLine(x[0], m1i, x[n - 1], m1i);
    l1i->SetLineColor(kRed);
    l1i->SetLineStyle(2);
    l1i->Draw();
    auto leg1 = new TLegend(0.15, 0.77, 0.4, 0.9);
    leg1->AddEntry(gI_border, "border", "lep");
    leg1->AddEntry(gI_inner, "inner", "lep");
    leg1->AddEntry(l1b, "mean border", "l");
    leg1->AddEntry(l1i, "mean inner", "l");
    leg1->Draw();

    // 2. I two zooms
    c2->cd();
    TGraphErrors *gI_borderz = new TGraphErrors(n, x, I_border_zoom, 0, I_border_zoom_std);
    TGraphErrors *gI_innerz  = new TGraphErrors(n, x, I_inner_zoom, 0, I_inner_zoom_std);
    gI_borderz->SetMarkerColor(kBlue + 2);
    gI_borderz->SetLineColor(kBlue + 2);
    gI_borderz->SetMarkerStyle(22);
    gI_innerz->SetMarkerColor(kOrange + 7);
    gI_innerz->SetLineColor(kOrange + 7);
    gI_innerz->SetMarkerStyle(23);
    gI_borderz->SetTitle("Mean current at 120 V (two zooms);Channel;I [nA]");

    double I2_min = std::min(*std::min_element(I_border_zoom, I_border_zoom + n), *std::min_element(I_inner_zoom, I_inner_zoom + n));
    double I2_max = std::max(*std::max_element(I_border_zoom, I_border_zoom + n), *std::max_element(I_inner_zoom, I_inner_zoom + n));
    gI_borderz->GetYaxis()->SetRangeUser(I2_min * 0.9, I2_max * 1.1);

    gI_borderz->Draw("AP");
    gI_innerz->Draw("P SAME");
    double m2b = mean(I_border_zoom, n), m2i = mean(I_inner_zoom, n);
    TLine *l2b = new TLine(x[0], m2b, x[n - 1], m2b);
    l2b->SetLineColor(kBlue + 2);
    l2b->SetLineStyle(2);
    l2b->Draw();
    TLine *l2i = new TLine(x[0], m2i, x[n - 1], m2i);
    l2i->SetLineColor(kOrange + 7);
    l2i->SetLineStyle(2);
    l2i->Draw();
    auto leg2 = new TLegend(0.15, 0.77, 0.4, 0.9);
    leg2->AddEntry(gI_borderz, "border zoom", "lep");
    leg2->AddEntry(gI_innerz, "inner zoom", "lep");
    leg2->AddEntry(l2b, "mean border", "l");
    leg2->AddEntry(l2i, "mean inner", "l");
    leg2->Draw();

    // 3. Vdep
    c3->cd();
    TGraphErrors *gVdep = new TGraphErrors(n, x, Vdep, 0, Vdep_std);
    gVdep->SetMarkerColor(kMagenta + 2);
    gVdep->SetLineColor(kMagenta + 2);
    gVdep->SetMarkerStyle(21);
    gVdep->SetTitle("Mean depletion voltage;Channel;V_{dep} [V]");
    // Adjust y-axis range for Vdep
    gVdep->GetYaxis()->SetRangeUser(*std::min_element(Vdep, Vdep + n) * 0.9, *std::max_element(Vdep, Vdep + n) * 1.1);

    gVdep->Draw("AP");
    double mVdep = mean(Vdep, n);
    TLine *l3v = new TLine(x[0], mVdep, x[n - 1], mVdep);
    l3v->SetLineColor(kMagenta + 2);
    l3v->SetLineStyle(2);
    l3v->Draw();
    auto leg3 = new TLegend(0.15, 0.83, 0.35, 0.93);
    leg3->AddEntry(gVdep, "mean values", "lep");
    leg3->AddEntry(l3v, "mean", "l");
    leg3->Draw();

    // 4. cs (pF)
    c4->cd();
    TGraphErrors *gcsb = new TGraphErrors(n, x, cs_border, 0, cs_border_std);
    TGraphErrors *gcsi = new TGraphErrors(n, x, cs_inner, 0, cs_inner_std);
    gcsb->SetMarkerColor(kGreen + 2);
    gcsb->SetLineColor(kGreen + 2);
    gcsb->SetMarkerStyle(22);
    gcsi->SetMarkerColor(kRed + 2);
    gcsi->SetLineColor(kRed + 2);
    gcsi->SetMarkerStyle(23);
    gcsb->SetTitle("Mean full depletion capacitance;Channel;C_{fd} [pF]");
    // Set y-axis range to include both groups fully
    double cs_min = std::min(*std::min_element(cs_border, cs_border + n), *std::min_element(cs_inner, cs_inner + n));
    double cs_max = std::max(*std::max_element(cs_border, cs_border + n), *std::max_element(cs_inner, cs_inner + n));
    gcsb->GetYaxis()->SetRangeUser(cs_min * 0.9, cs_max * 1.1);

    gcsb->Draw("AP");
    gcsi->Draw("P SAME");

    double mcsb = mean(cs_border, n), mcsi = mean(cs_inner, n);

    // For cs capacitance draw two means: channels (0,1,2,3,8) and (4,5,6,7) zero-based indices
    std::vector<int> groupA = {0, 1, 2, 3, 8}; // channels 9,10,13,14,20
    std::vector<int> groupB = {4, 5, 6, 7};    // channels 15,16,17,18

    double mcsbA = mean_subset(cs_border, groupA);
    double mcsbB = mean_subset(cs_border, groupB);
    double mcsiA = mean_subset(cs_inner, groupA);
    double mcsiB = mean_subset(cs_inner, groupB);

    TLine *l4bA = new TLine(x[0], mcsbA, x[n - 1], mcsbA);
    l4bA->SetLineColor(kGreen + 2);
    l4bA->SetLineStyle(2);
    l4bA->Draw();
    TLine *l4bB = new TLine(x[0], mcsbB, x[n - 1], mcsbB);
    l4bB->SetLineColor(kGreen + 2);
    l4bB->SetLineStyle(3);
    l4bB->Draw();

    TLine *l4iA = new TLine(x[0], mcsiA, x[n - 1], mcsiA);
    l4iA->SetLineColor(kRed + 2);
    l4iA->SetLineStyle(2);
    l4iA->Draw();
    TLine *l4iB = new TLine(x[0], mcsiB, x[n - 1], mcsiB);
    l4iB->SetLineColor(kRed + 2);
    l4iB->SetLineStyle(3);
    l4iB->Draw();

    auto leg4 = new TLegend(0.15, 0.77, 0.48, 0.9);
    leg4->AddEntry(gcsb, "border", "lep");
    leg4->AddEntry(gcsi, "inner", "lep");
    leg4->AddEntry(l4bA, "mean border group 9,10,13,14,20", "l");
    leg4->AddEntry(l4bB, "mean border group 15,16,17,18", "l");
    leg4->AddEntry(l4iA, "mean inner group 9,10,13,14,20", "l");
    leg4->AddEntry(l4iB, "mean inner group 15,16,17,18", "l");
    leg4->Draw();

    // 5. ndon (ne/cm^3)
    c5->cd();
    TGraphErrors *gndb = new TGraphErrors(n, x, ndon_border, 0, ndon_border_std);
    TGraphErrors *gndi = new TGraphErrors(n, x, ndon_inner, 0, ndon_inner_std);
    gndb->SetMarkerColor(kBlack);
    gndb->SetLineColor(kBlack);
    gndb->SetMarkerStyle(29);
    gndi->SetMarkerColor(kPink + 6);
    gndi->SetLineColor(kPink + 6);
    gndi->SetMarkerStyle(30);
    gndb->SetTitle("Mean donor density;Channel;n_{don} [ne/cm^{3}]");

    // Adjust y-axis range for donor density (include all points)
    double nd_min = std::min(*std::min_element(ndon_border, ndon_border + n), *std::min_element(ndon_inner, ndon_inner + n));
    double nd_max = std::max(*std::max_element(ndon_border, ndon_border + n), *std::max_element(ndon_inner, ndon_inner + n));
    gndb->GetYaxis()->SetRangeUser(nd_min * 0.9, nd_max * 1.1);

    gndb->Draw("AP");
    gndi->Draw("P SAME");

    // Two groups mean lines: channels {9,10,13,14,20} and {15,16,17,18}
    double mndbA = mean_subset(ndon_border, groupA);
    double mndbB = mean_subset(ndon_border, groupB);
    double mndiA = mean_subset(ndon_inner, groupA);
    double mndiB = mean_subset(ndon_inner, groupB);

    TLine *l5bA = new TLine(x[0], mndbA, x[n - 1], mndbA);
    l5bA->SetLineColor(kBlack);
    l5bA->SetLineStyle(2);
    l5bA->Draw();
    TLine *l5bB = new TLine(x[0], mndbB, x[n - 1], mndbB);
    l5bB->SetLineColor(kBlack);
    l5bB->SetLineStyle(3);
    l5bB->Draw();

    TLine *l5iA = new TLine(x[0], mndiA, x[n - 1], mndiA);
    l5iA->SetLineColor(kPink + 6);
    l5iA->SetLineStyle(2);
    l5iA->Draw();
    TLine *l5iB = new TLine(x[0], mndiB, x[n - 1], mndiB);
    l5iB->SetLineColor(kPink + 6);
    l5iB->SetLineStyle(3);
    l5iB->Draw();

    auto leg5 = new TLegend(0.55, 0.77, 0.9, 0.9);
    leg5->AddEntry(gndb, "border", "lep");
    leg5->AddEntry(gndi, "inner", "lep");
    leg5->AddEntry(l5bA, "mean border group 9,10,13,14,20", "l");
    leg5->AddEntry(l5bB, "mean border group 15,16,17,18", "l");
    leg5->AddEntry(l5iA, "mean inner group 9,10,13,14,20", "l");
    leg5->AddEntry(l5iB, "mean inner group 15,16,17,18", "l");
    leg5->Draw();

    TFile* f = new TFile("all_sensors_results.root", "RECREATE");
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    f->Close();



}
