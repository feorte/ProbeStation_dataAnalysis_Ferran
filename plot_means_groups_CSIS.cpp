#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLine.h>
#include <iostream>
#include <vector>

void plot_means_groups_CSIS() {
    const int n = 8; // 9 points (sensor ids)
    double x[n] = {1, 2, 3, 4, 5, 6, 7, 8}; // sensor IDs

    // 1) I (nA) FIRST "ONE ZOOM" block
    double I_border[n] = {0.102106, 0.124096, 0.141025, 0.0996983, 0.133513, 0.130092, 0.236907, 0.116035};
    double I_border_std[n] = {
        0.00136566, 0.00207136, 0.003234042, 0.001537522,
        0.002476146, 0.003149085, 0.083538352, 0.001503999
    };

    double I_inner[n] = {0.0277687, 0.0323819, 0.0422819, 0.0261969, 0.02996, 0.0325467, 0.0487106, 0.0498156};
    double I_inner_std[n] = {
        0.00066764, 0.00116804, 0.00121416, 0.0005778,
        0.0006884, 0.00111217, 0.00298892, 0.00326956
    };

    // 2) Vdep
    double Vdep[n] = {5.06603e+01, 4.61005e+01, 4.92717e+01, 3.94758e+01, 
                      3.97025e+01, 3.95328e+01, 3.94188e+01, 3.97485e+01};

    double Vdep_std[n] = {
        3.67608e-02, 9.37175e-02, 3.18743e-02, 2.90011e-02,
        5.11238e-02, 4.08274e-02, 3.28996e-02, 3.83613e-02
    };
    double Vdep_sigma[n] = {
        8.81010e-01, 1.18942e+00, 1.47291e+00, 2.71412e-01,
        4.64647e-01, 4.30186e-01, 4.99530e-01, 4.53411e-01
    };

    // 3) cs (pF)

    double cs_inner[n] = {4.81578e+00, 4.80319e+00, 4.80479e+00, 4.79513e+00, 
                          4.81693e+00, 4.82159e+00, 4.72190e+00, 4.72744e+00};

    double cs_inner_std[n] = {
        4.79163e-04, 7.37504e-04,  8.88164e-04, 1.07228e-03,
        8.62581e-04, 7.27618e-04, 9.70663e-04, 8.66491e-04};

    double cs_inner_sigma[n] = {
        3.33381e-02, 2.41327e-02, 3.17854e-02, 2.93461e-02,
        2.35649e-02, 1.75699e-02, 2.28334e-02, 1.52061e-02
    };

    double cs_border[n] = {4.99895e+00, 4.99010e+00, 4.96726e+00, 4.94736e+00, 
                           4.94510e+00, 4.97209e+00, 4.89675e+00, 4.90434e+00};

    double cs_border_std[n] = {
        1.06944e-03, 1.80462e-03, 1.18650e-03, 2.37933e-03,
        8.54184e-02, 1.72863e-03, 2.85368e-03, 2.25959e-03
    };

    double cs_border_sigma[n] = {
        6.49411e-02, /*1.37339e-01*/6.09032e-02,  6.59735e-02, 7.54100e-02,
        9.11386e-03, 7.01133e-02, 9.81383e-02, 5.57984e-02
    };



    // 4) ndon (ne/cm^3)

    double ndon_inner[n] = {
        1.55153e+11, 1.42815e+11, -1, 1.18213e+11,
        1.19795e+11, 1.19034e+11, 1.16531e+11, 1.16940e+11
    };
    double ndon_inner_std[n] = {
        1.65130e+08, 6.00411e+08, 3.49825e+08, 8.98481e+07,
        1.36326e+08, 1.18543e+08, 1.34866e+08, 1.22382e+08
    };
    double ndon_inner_sigma[n] = {
        2.09555e+09, 2.39920e+09,  0, 2.19398e+09,
        1.95682e+09, 2.20134e+09, 2.13760e+09, 1.94701e+09
    };


    double ndon_border[n] = {
        1.64833e+11, 1.50792e+11, -1, 1.27893e+11,
       1.32475e+11, 1.25328e+11, 1.27645e+111, 1.27761e+11
    };
    double ndon_border_std[n] = {
        2.24065e+08, 4.81926e+08, 4.31075e+08, 2.92711e+08,
        3.29726e+08, 2.49179e+08, 3.64356e+08, 3.28163e+08
    };
    double ndon_border_sigma[n] = {
        5.22334e+09, 5.19116e+09, 0, 9.34799e+09,
        4.97980e+09, 2.97154e+09, 3.72298e+09, 4.40407e+09
    };




    // 5) "TWO ZOOMS": I (nA) second measurement (zoomed)

    double I_inner_zoom[n] = {3.35711e-02, 2.19322e-02, 2.27658e-02, 2.82497e-02, 
                              2.11307e-02, 2.77180e-02, 2.29835e-02, 2.55917e-02};

    double I_inner_zoom_sigma[n] = {
        4.56371e-03, 7.14088e-03, 2.84803e-03, 4.41030e-03,
        4.63249e-03, 4.83072e-03, 4.72603e-03, 2.72624e-03
    };


    double I_border_zoom[n] = {1.06852e-01, 8.28184e-02, 7.37621e-02, 1.03035e-01, 
        7.29869e-02, 1.01432e-01, 9.08203e-02, 1.04323e-01};

    double I_border_zoom_sigma[n] = {
        1.37942e-02, 7.41617e-03, 5.92693e-03, 1.31972e-02,
        6.08683e-03, 8.88252e-03, 6.13501e-03, 7.16592e-03
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
    gI_border->SetTitle("Mean current at 120 V (one zoom);Sensor ID;Mean current [nA]");

    // Adjust y-axis range to include all points + margin
    double I_min = std::min(*std::min_element(I_border, I_border + n), *std::min_element(I_inner, I_inner + n));
    double I_max = std::max(*std::max_element(I_border, I_border + n), *std::max_element(I_inner, I_inner + n));
    gI_border->GetYaxis()->SetRangeUser(I_min * 0.9, I_max * 1.1);

    gI_border->Draw("AP");
    gI_inner->Draw("P SAME");
    double m1b = mean(I_border, n), m1i = mean(I_inner, n);

    // --- compute error of the mean for border
    double sumsq_b = 0;
    for (int i = 0; i < n; i++) sumsq_b += (I_border[i] - m1b) * (I_border[i] - m1b);
    double sigma_b = sqrt(sumsq_b / (n - 1));   // standard deviation
    double err_b = sigma_b / sqrt(n);           // error of the mean

    // --- compute error of the mean for inner
    double sumsq_i = 0;
    for (int i = 0; i < n; i++) sumsq_i += (I_inner[i] - m1i) * (I_inner[i] - m1i);
    double sigma_i = sqrt(sumsq_i / (n - 1));
    double err_i = sigma_i / sqrt(n);

    TLine *l1b = new TLine(x[0], m1b, x[n - 1], m1b);
    l1b->SetLineColor(kBlue);
    l1b->SetLineStyle(2);
    l1b->Draw();

    TLine *l1i = new TLine(x[0], m1i, x[n - 1], m1i);
    l1i->SetLineColor(kRed);
    l1i->SetLineStyle(2);
    l1i->Draw();

    // --- draw error band for border mean
    TBox *errBoxB = new TBox(x[0], m1b - err_b, x[n - 1], m1b + err_b);
    errBoxB->SetFillColorAlpha(kBlue, 0.2); // semi-transparent blue
    errBoxB->SetLineColor(kBlue);
    errBoxB->Draw();

    // --- draw error band for inner mean
    TBox *errBoxI = new TBox(x[0], m1i - err_i, x[n - 1], m1i + err_i);
    errBoxI->SetFillColorAlpha(kRed, 0.2); // semi-transparent red
    errBoxI->SetLineColor(kRed);
    errBoxI->Draw();

    // --- legend
    auto leg1 = new TLegend(0.15, 0.72, 0.45, 0.9);
    leg1->AddEntry(gI_border, "border", "lep");
    leg1->AddEntry(gI_inner, "inner", "lep");
    leg1->AddEntry(l1b, Form("all sensors mean (border) = %.2f", m1b), "l");
    leg1->AddEntry(errBoxB, Form("error of mean (border) = %.2f", err_b), "f");
    leg1->AddEntry(l1i, Form("all sensors mean (inner) = %.2f", m1i), "l");
    leg1->AddEntry(errBoxI, Form("error of mean (inner) = %.2f", err_i), "f");
    leg1->Draw();

    // // 2. I two zooms
    // c2->cd();
    // TGraphErrors *gI_borderz = new TGraphErrors(n, x, I_border_zoom, 0, I_border_zoom_sigma);
    // TGraphErrors *gI_innerz  = new TGraphErrors(n, x, I_inner_zoom, 0, I_inner_zoom_sigma);
    // gI_borderz->SetMarkerColor(kBlue + 2);
    // gI_borderz->SetLineColor(kBlue + 2);
    // gI_borderz->SetMarkerStyle(22);
    // gI_innerz->SetMarkerColor(kOrange + 7);
    // gI_innerz->SetLineColor(kOrange + 7);
    // gI_innerz->SetMarkerStyle(23);
    // gI_borderz->SetTitle("Mean current at 120 V;Sensor ID;Mean current [nA]");

    // double I2_min = std::min(*std::min_element(I_border_zoom, I_border_zoom + n), *std::min_element(I_inner_zoom, I_inner_zoom + n));
    // double I2_max = std::max(*std::max_element(I_border_zoom, I_border_zoom + n), *std::max_element(I_inner_zoom, I_inner_zoom + n));
    // gI_borderz->GetYaxis()->SetRangeUser(I2_min * 0.8, I2_max * 1.2);

    // gI_borderz->Draw("AP");
    // gI_innerz->Draw("P SAME");
    // double m2b = mean(I_border_zoom, n), m2i = mean(I_inner_zoom, n);
    // TLine *l2b = new TLine(x[0], m2b, x[n - 1], m2b);
    // l2b->SetLineColor(kBlue + 2);
    // l2b->SetLineStyle(2);
    // l2b->Draw();
    // TLine *l2i = new TLine(x[0], m2i, x[n - 1], m2i);
    // l2i->SetLineColor(kOrange + 7);
    // l2i->SetLineStyle(2);
    // l2i->Draw();
    // auto leg2 = new TLegend(0.15, 0.77, 0.4, 0.9);
    // leg2->AddEntry(gI_borderz, "border", "lep");
    // leg2->AddEntry(gI_innerz, "inner", "lep");
    // leg2->AddEntry(l2b, Form("all sensors mean (border) = %.2f", m2b), "l");
    // leg2->AddEntry(l2i, Form("all sensors mean (inner) = %.2f", m2i), "l");
    // leg2->Draw();

    // // 3. Vdep
    // c3->cd();
    // TGraphErrors *gVdep = new TGraphErrors(n, x, Vdep, 0, Vdep_sigma);
    // gVdep->SetMarkerColor(kMagenta + 2);
    // gVdep->SetLineColor(kMagenta + 2);
    // gVdep->SetMarkerStyle(21);
    // gVdep->SetTitle("Mean depletion voltage;Sensor ID;Mean V_{dep} [V]");
    // // Adjust y-axis range for Vdep
    // gVdep->GetYaxis()->SetRangeUser(*std::min_element(Vdep, Vdep + n) * 0.9, *std::max_element(Vdep, Vdep + n) * 1.1);

    // gVdep->Draw("AP");
    // double mVdep = mean(Vdep, n);
    // TLine *l3v = new TLine(x[0], mVdep, x[n - 1], mVdep);
    // l3v->SetLineColor(kMagenta + 2);
    // l3v->SetLineStyle(2);
    // l3v->Draw();
    // auto leg3 = new TLegend(0.15, 0.83, 0.35, 0.93);
    // leg3->AddEntry(gVdep, "mean values", "lep");
    // leg3->AddEntry(l3v, Form("all sensors mean = %.2f",mVdep), "l");
    // leg3->Draw();

    // // 4. cs (pF)
    // c4->cd();
    // TGraphErrors *gcsb = new TGraphErrors(n, x, cs_border, 0, cs_border_sigma);
    // TGraphErrors *gcsi = new TGraphErrors(n, x, cs_inner, 0, cs_inner_sigma);
    // gcsb->SetMarkerColor(kGreen + 2);
    // gcsb->SetLineColor(kGreen + 2);
    // gcsb->SetMarkerStyle(22);
    // gcsi->SetMarkerColor(kRed + 2);
    // gcsi->SetLineColor(kRed + 2);
    // gcsi->SetMarkerStyle(23);
    // gcsb->SetTitle("Mean full depletion capacitance;Sensor ID;C_{fd} [pF]");
    // // Set y-axis range to include both groups fully
    // double cs_min = std::min(*std::min_element(cs_border, cs_border + n), *std::min_element(cs_inner, cs_inner + n));
    // double cs_max = std::max(*std::max_element(cs_border, cs_border + n), *std::max_element(cs_inner, cs_inner + n));
    // gcsb->GetYaxis()->SetRangeUser(cs_min * 0.9, cs_max * 1.1);

    // gcsb->Draw("AP");
    // gcsi->Draw("P SAME");

    // double mcsb = mean(cs_border, n), mcsi = mean(cs_inner, n);

    // // For cs capacitance draw two means: channels (0,1,2,3,8) and (4,5,6,7) zero-based indices
    // std::vector<int> groupA = {0, 1, 2, 3, 8}; // channels 9,10,13,14,20
    // std::vector<int> groupB = {4, 5, 6, 7};    // channels 15,16,17,18

    // double mcsbA = mean_subset(cs_border, groupA);
    // double mcsbB = mean_subset(cs_border, groupB);
    // double mcsiA = mean_subset(cs_inner, groupA);
    // double mcsiB = mean_subset(cs_inner, groupB);

    // TLine *l4bA = new TLine(x[0], mcsb, x[n - 1], mcsb);
    // l4bA->SetLineColor(kGreen + 2);
    // l4bA->SetLineStyle(2);
    // l4bA->Draw();
    // // TLine *l4bB = new TLine(x[0], mcsbB, x[n - 1], mcsbB);
    // // l4bB->SetLineColor(kGreen + 2);
    // // l4bB->SetLineStyle(3);
    // // l4bB->Draw();

    // TLine *l4iA = new TLine(x[0], mcsi, x[n - 1], mcsi);
    // l4iA->SetLineColor(kRed + 2);
    // l4iA->SetLineStyle(2);
    // l4iA->Draw();
    // // TLine *l4iB = new TLine(x[0], mcsiB, x[n - 1], mcsiB);
    // // l4iB->SetLineColor(kRed + 2);
    // // l4iB->SetLineStyle(3);
    // // l4iB->Draw();

    // auto leg4 = new TLegend(0.15, 0.77, 0.48, 0.9);
    // leg4->AddEntry(gcsb, "border", "lep");
    // leg4->AddEntry(gcsi, "inner", "lep");
    // leg4->AddEntry(l4bA, Form("all sensors mean (border) = %.2f",mcsb), "l");
    // // leg4->AddEntry(l4bB, Form("sensors 15,16,17,18 mean (border) = %.2f",mcsbB), "l");
    // // leg4->AddEntry(l4iB, Form("sensors 9,10,13,14,20 mean (inner) = %.2f",mcsiA), "l");
    // leg4->AddEntry(l4iA, Form("all sensors mean (inner) = %.2f",mcsi), "l");
    // leg4->Draw();

    // // 5. ndon (ne/cm^3)
    // c5->cd();
    // TGraphErrors *gndb = new TGraphErrors(n, x, ndon_border, 0, ndon_border_sigma);
    // TGraphErrors *gndi = new TGraphErrors(n, x, ndon_inner, 0, ndon_inner_sigma);
    // gndb->SetMarkerColor(kBlack);
    // gndb->SetLineColor(kBlack);
    // gndb->SetMarkerStyle(29);
    // gndi->SetMarkerColor(kPink + 6);
    // gndi->SetLineColor(kPink + 6);
    // gndi->SetMarkerStyle(30);
    // gndb->SetTitle("Mean donor density;Sensor ID;Mean n_{don} [ne/cm^{3}]");

    // // Adjust y-axis range for donor density (include all points)
    // double nd_min = std::min(*std::min_element(ndon_border, ndon_border + n), *std::min_element(ndon_inner, ndon_inner + n));
    // double nd_max = std::max(*std::max_element(ndon_border, ndon_border + n), *std::max_element(ndon_inner, ndon_inner + n));
    // gndb->GetYaxis()->SetRangeUser(nd_min * 0.9, nd_max * 1.1);

    // gndb->Draw("AP");
    // gndi->Draw("P SAME");

    // double mndb = mean(ndon_border, n), mndi = mean(ndon_inner, n);

    // // Two groups mean lines: channels {9,10,13,14,20} and {15,16,17,18}
    // double mndbA = mean_subset(ndon_border, groupA);
    // double mndbB = mean_subset(ndon_border, groupB);
    // double mndiA = mean_subset(ndon_inner, groupA);
    // double mndiB = mean_subset(ndon_inner, groupB);

    // TLine *l5bA = new TLine(x[0], mndb, x[n - 1], mndb);
    // l5bA->SetLineColor(kBlack);
    // l5bA->SetLineStyle(2);
    // l5bA->Draw();
    // // TLine *l5bB = new TLine(x[0], mndbB, x[n - 1], mndbB);
    // // l5bB->SetLineColor(kBlack);
    // // l5bB->SetLineStyle(3);
    // // l5bB->Draw();

    // // TLine *l5iA = new TLine(x[0], mndiA, x[n - 1], mndiA);
    // // l5iA->SetLineColor(kPink + 6);
    // // l5iA->SetLineStyle(2);
    // // l5iA->Draw();
    // TLine *l5iB = new TLine(x[0], mndi, x[n - 1], mndi);
    // l5iB->SetLineColor(kPink + 6);
    // l5iB->SetLineStyle(3);
    // l5iB->Draw();

    // auto leg5 = new TLegend(0.55, 0.77, 0.9, 0.9);
    // leg5->AddEntry(gndb, "border", "lep");
    // leg5->AddEntry(gndi, "inner", "lep");

    // leg5->AddEntry(l5bA, Form("all sensors mean (border) = %.2e",mndb), "l");
    // // leg5->AddEntry(l5bB, Form("sensors 15,16,17,18 mean (border) = %.2e",mndbB), "l");
    // // leg5->AddEntry(l5iA, Form("sensors 9,10,13,14,20 mean (inner) = %.2e",mndiA), "l");
    // leg5->AddEntry(l5iB, Form("all sensors mean (inner) = %.2e",mndi), "l");
    // leg5->Draw();

    // TFile* f = new TFile("all_sensors_results.root", "RECREATE");
    // c1->Write();
    // c2->Write();
    // c3->Write();
    // c4->Write();
    // c5->Write();
    // f->Close();

        // 2. I two zooms
    c2->cd();
    TGraphErrors *gI_borderz = new TGraphErrors(n, x, I_border_zoom, 0, I_border_zoom_sigma);
    TGraphErrors *gI_innerz  = new TGraphErrors(n, x, I_inner_zoom, 0, I_inner_zoom_sigma);
    gI_borderz->SetMarkerColor(kBlue + 2);
    gI_borderz->SetLineColor(kBlue + 2);
    gI_borderz->SetMarkerStyle(22);
    gI_innerz->SetMarkerColor(kOrange + 7);
    gI_innerz->SetLineColor(kOrange + 7);
    gI_innerz->SetMarkerStyle(23);
    gI_borderz->SetTitle("Mean current at 120 V;Sensor ID;Mean current [nA]");

    double I2_min = std::min(*std::min_element(I_border_zoom, I_border_zoom + n), *std::min_element(I_inner_zoom, I_inner_zoom + n));
    double I2_max = std::max(*std::max_element(I_border_zoom, I_border_zoom + n), *std::max_element(I_inner_zoom, I_inner_zoom + n));
    gI_borderz->GetYaxis()->SetRangeUser(0.01, I2_max * 1.5);

    gI_borderz->Draw("AP");
    gI_innerz->Draw("P SAME");

    // mean + error for border zoom
    double m2b = mean(I_border_zoom, n), m2i = mean(I_inner_zoom, n);
    double sumsq2b=0, sumsq2i=0;
    for (int i=0;i<n;i++) {
        sumsq2b += (I_border_zoom[i]-m2b)*(I_border_zoom[i]-m2b);
        sumsq2i += (I_inner_zoom[i]-m2i)*(I_inner_zoom[i]-m2i);
    }
    double err2b = sqrt(sumsq2b/(n-1))/sqrt(n);
    double err2i = sqrt(sumsq2i/(n-1))/sqrt(n);

    TLine *l2b = new TLine(x[0], m2b, x[n - 1], m2b);
    l2b->SetLineColor(kBlue + 2);
    l2b->SetLineStyle(2);
    l2b->Draw();
    TBox *errBox2b = new TBox(x[0], m2b-err2b, x[n-1], m2b+err2b);
    errBox2b->SetFillColorAlpha(kBlue+2,0.2);
    errBox2b->SetLineColor(kBlue+2);
    errBox2b->Draw();

    TLine *l2i = new TLine(x[0], m2i, x[n - 1], m2i);
    l2i->SetLineColor(kOrange + 7);
    l2i->SetLineStyle(2);
    l2i->Draw();
    TBox *errBox2i = new TBox(x[0], m2i-err2i, x[n-1], m2i+err2i);
    errBox2i->SetFillColorAlpha(kOrange+7,0.2);
    errBox2i->SetLineColor(kOrange+7);
    errBox2i->Draw();

    auto leg2 = new TLegend(0.15, 0.72, 0.45, 0.9);
    leg2->AddEntry(gI_borderz, "border", "lep");
    leg2->AddEntry(gI_innerz, "inner", "lep");
    leg2->AddEntry(l2b, Form("all sensors mean (border) = %.3f", m2b), "l");
    leg2->AddEntry(errBox2b, Form("error of mean (border) = %.3f", err2b), "f");
    leg2->AddEntry(l2i, Form("all sensors mean (inner) = %.3f", m2i), "l");
    leg2->AddEntry(errBox2i, Form("error of mean (inner) = %.3f", err2i), "f");
    leg2->Draw();

    // 3. Vdep
    c3->cd();
    TGraphErrors *gVdep = new TGraphErrors(n, x, Vdep, 0, Vdep_sigma);
    gVdep->SetMarkerColor(kMagenta + 2);
    gVdep->SetLineColor(kMagenta + 2);
    gVdep->SetMarkerStyle(21);
    gVdep->SetTitle("Mean depletion voltage;Sensor ID;Mean V_{dep} [V]");
    gVdep->GetYaxis()->SetRangeUser(*std::min_element(Vdep, Vdep + n) * 0.9, *std::max_element(Vdep, Vdep + n) * 1.1);
    gVdep->Draw("AP");

    // 2) Vdep
    double Vdephigh[3] = {5.06603e+01, 4.61005e+01, 4.92717e+01};
    double Vdeplow[5] = {3.94758e+01, 3.97025e+01, 3.95328e+01, 3.94188e+01, 3.97485e+01};

    double Vdephigh_sigma[3] = {
        8.81010e-01, 1.18942e+00, 1.47291e+00};
    double Vdeplow_sigma[5] = {2.71412e-01,
    4.64647e-01, 4.30186e-01, 4.99530e-01, 4.53411e-01};

    double mVdep = mean(Vdep, n);
    double mVdeph = mean(Vdephigh, 3);
    double mVdepl = mean(Vdeplow, 5);
    double sumsqV=0; for(int i=0;i<n;i++) sumsqV+=(Vdep[i]-mVdep)*(Vdep[i]-mVdep);
    double sumsqVh=0; for(int i=0;i<3;i++) sumsqVh+=(Vdephigh[i]-mVdeph)*(Vdephigh[i]-mVdeph);
    double sumsqVl=0; for(int i=0;i<5;i++) sumsqVl+=(Vdeplow[i]-mVdepl)*(Vdeplow[i]-mVdepl);
    double errVdep = sqrt(sumsqV/(n-1))/sqrt(n);
    double errVdeph = sqrt(sumsqVh/(3-1))/sqrt(3);
    double errVdepl = sqrt(sumsqVl/(5-1))/sqrt(5);

    // TLine *l3v = new TLine(x[0], mVdep, x[n - 1], mVdep);
    // l3v->SetLineColor(kMagenta + 2);
    // l3v->SetLineStyle(2);
    // l3v->Draw();
    // TBox *errBox3v = new TBox(x[0], mVdep-errVdep, x[n-1], mVdep+errVdep);
    // errBox3v->SetFillColorAlpha(kMagenta+2,0.2);
    // errBox3v->SetLineColor(kMagenta+2);
    // errBox3v->Draw();

    TLine *l3vh = new TLine(x[0], mVdeph, x[n - 1], mVdeph);
    l3vh->SetLineColor(kMagenta + 2);
    l3vh->SetLineStyle(3);
    l3vh->Draw();
    TBox *errBox3vh = new TBox(x[0], mVdeph-errVdeph, x[n-1], mVdeph+errVdeph);
    errBox3vh->SetFillColorAlpha(kMagenta+2,0.1);
    errBox3vh->SetLineColor(kMagenta+2);
    errBox3vh->Draw();

    TLine *l3vl = new TLine(x[0], mVdepl, x[n - 1], mVdepl);
    l3vl->SetLineColor(kMagenta + 2);
    l3vl->SetLineStyle(3);
    l3vl->Draw();
    TBox *errBox3vl = new TBox(x[0], mVdepl-errVdepl, x[n-1], mVdepl+errVdepl);
    errBox3vl->SetFillColorAlpha(kMagenta+2,0.1);
    errBox3vl->SetLineColor(kMagenta+2);
    errBox3vl->Draw();


    auto leg3 = new TLegend(0.15, 0.77, 0.45, 0.9);
    leg3->AddEntry(gVdep, "mean values", "lep");
    // leg3->AddEntry(l3v, Form("all sensors mean = %.2f",mVdep), "l");
    // leg3->AddEntry(errBox3v, Form("error of mean = %.2f",errVdep), "f");
    leg3->AddEntry(l3vh, Form("high V_{dep} group mean = %.2f",mVdeph), "l");
    leg3->AddEntry(errBox3vh, Form("error of mean (high V_{dep}) = %.2f",errVdeph), "f");
    leg3->AddEntry(l3vl, Form("low V_{dep} group mean = %.2f",mVdepl), "l");
    leg3->AddEntry(errBox3vl, Form("error of mean (low V_{dep}) = %.2f",errVdepl), "f");
    leg3->Draw();

    // 4. cs (pF)
    c4->cd();
    TGraphErrors *gcsb = new TGraphErrors(n, x, cs_border, 0, cs_border_sigma);
    TGraphErrors *gcsi = new TGraphErrors(n, x, cs_inner, 0, cs_inner_sigma);
    gcsb->SetMarkerColor(kGreen + 2);
    gcsb->SetLineColor(kGreen + 2);
    gcsb->SetMarkerStyle(22);
    gcsi->SetMarkerColor(kRed + 2);
    gcsi->SetLineColor(kRed + 2);
    gcsi->SetMarkerStyle(23);
    gcsb->SetTitle("Mean full depletion capacitance;Sensor ID;C_{fd} [pF]");
    double cs_min = std::min(*std::min_element(cs_border, cs_border + n), *std::min_element(cs_inner, cs_inner + n));
    double cs_max = std::max(*std::max_element(cs_border, cs_border + n), *std::max_element(cs_inner, cs_inner + n));
    gcsb->GetYaxis()->SetRangeUser(cs_min * 0.9, cs_max * 1.1);
    gcsb->Draw("AP");
    gcsi->Draw("P SAME");

    double mcsb = mean(cs_border, n), mcsi = mean(cs_inner, n);
    double sumsqcsb=0,sumsqcsi=0;
    for(int i=0;i<n;i++){sumsqcsb+=(cs_border[i]-mcsb)*(cs_border[i]-mcsb); sumsqcsi+=(cs_inner[i]-mcsi)*(cs_inner[i]-mcsi);}
    double errcsb = sqrt(sumsqcsb/(n-1))/sqrt(n);
    double errcsi = sqrt(sumsqcsi/(n-1))/sqrt(n);

    TLine *l4bA = new TLine(x[0], mcsb, x[n - 1], mcsb);
    l4bA->SetLineColor(kGreen + 2);
    l4bA->SetLineStyle(2);
    l4bA->Draw();
    TBox *errBox4b = new TBox(x[0], mcsb-errcsb, x[n-1], mcsb+errcsb);
    errBox4b->SetFillColorAlpha(kGreen+2,0.2);
    errBox4b->SetLineColor(kGreen+2);
    errBox4b->Draw("same");

    TLine *l4iA = new TLine(x[0], mcsi, x[n - 1], mcsi);
    l4iA->SetLineColor(kRed + 2);
    l4iA->SetLineStyle(2);
    l4iA->Draw();
    TBox *errBox4i = new TBox(x[0], mcsi-errcsi, x[n-1], mcsi+errcsi);
    errBox4i->SetFillColorAlpha(kRed+2,0.2);
    errBox4i->SetLineColor(kRed+2);
    errBox4i->Draw("same");

    auto leg4 = new TLegend(0.15, 0.72, 0.48, 0.9);
    leg4->AddEntry(gcsb, "border", "lep");
    leg4->AddEntry(gcsi, "inner", "lep");
    leg4->AddEntry(l4bA, Form("all sensors mean (border) = %.3f",mcsb), "l");
    leg4->AddEntry(errBox4b, Form("error of mean (border) = %.3f",errcsb), "f");
    leg4->AddEntry(l4iA, Form("all sensors mean (inner) = %.3f",mcsi), "l");
    leg4->AddEntry(errBox4i, Form("error of mean (inner) = %.3f",errcsi), "f");
    leg4->Draw();

    // 5. ndon (ne/cm^3)
    c5->cd();
    TGraphErrors *gndb = new TGraphErrors(n, x, ndon_border, 0, ndon_border_sigma);
    TGraphErrors *gndi = new TGraphErrors(n, x, ndon_inner, 0, ndon_inner_sigma);
    gndb->SetMarkerColor(kBlack);
    gndb->SetLineColor(kBlack);
    gndb->SetMarkerStyle(29);
    gndi->SetMarkerColor(kPink + 6);
    gndi->SetLineColor(kPink + 6);
    gndi->SetMarkerStyle(30);
    gndb->SetTitle("Mean donor density;Sensor ID;Mean n_{don} [ne/cm^{3}]");

    double nd_min = std::min(*std::min_element(ndon_border, ndon_border + n), *std::min_element(ndon_inner, ndon_inner + n));
    double nd_max = std::max(*std::max_element(ndon_border, ndon_border + n), *std::max_element(ndon_inner, ndon_inner + n));
    gndb->GetYaxis()->SetRangeUser(nd_min * 0.9, nd_max * 1.1);

    gndb->Draw("AP");
    gndi->Draw("P SAME");

    double mndb = mean(ndon_border, n), mndi = mean(ndon_inner, n);
    double sumsqndb=0,sumsqndi=0;
    for(int i=0;i<n;i++){sumsqndb+=(ndon_border[i]-mndb)*(ndon_border[i]-mndb); sumsqndi+=(ndon_inner[i]-mndi)*(ndon_inner[i]-mndi);}
    double errndb = sqrt(sumsqndb/(n-1))/sqrt(n);
    double errndi = sqrt(sumsqndi/(n-1))/sqrt(n);

    TLine *l5bA = new TLine(x[0], mndb, x[n - 1], mndb);
    l5bA->SetLineColor(kBlack);
    l5bA->SetLineStyle(2);
    l5bA->Draw();
    TBox *errBox5b = new TBox(x[0], mndb-errndb, x[n-1], mndb+errndb);
    errBox5b->SetFillColorAlpha(kBlack,0.2);
    errBox5b->SetLineColor(kBlack);
    errBox5b->Draw("same");

    TLine *l5iB = new TLine(x[0], mndi, x[n - 1], mndi);
    l5iB->SetLineColor(kPink + 6);
    l5iB->SetLineStyle(2);
    l5iB->Draw();
    TBox *errBox5i = new TBox(x[0], mndi-errndi, x[n-1], mndi+errndi);
    errBox5i->SetFillColorAlpha(kPink+6,0.2);
    errBox5i->SetLineColor(kPink+6);
    errBox5i->Draw("same");

    auto leg5 = new TLegend(0.55, 0.72, 0.9, 0.9);
    leg5->AddEntry(gndb, "border", "lep");
    leg5->AddEntry(gndi, "inner", "lep");
    leg5->AddEntry(l5bA, Form("all sensors mean (border) = %.2e",mndb), "l");
    leg5->AddEntry(errBox5b, Form("error of mean (border) = %.2e",errndb), "f");
    leg5->AddEntry(l5iB, Form("all sensors mean (inner) = %.2e",mndi), "l");
    leg5->AddEntry(errBox5i, Form("error of mean (inner) = %.2e",errndi), "f");
    leg5->Draw();

    TFile* f = new TFile("all_sensors_results_CSIS.root", "RECREATE");
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    f->Close();


}
