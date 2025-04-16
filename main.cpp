#include <iostream>
#include <TROOT.h>
#include <fstream>
#include <string>
#include <vector> 


int main()
{
    // read the data
    std::ifstream data ("CALICE_6in_256ch_77_20250221_1_CV_customscan.txt"); // open the file directly when initializing the stream object
    if (data.is_open()) {

        std::string line;
        
        // read first the header lines which are not tabular data 
        for (int i=0; i<30; i++) {
            std::getline(data,line);
        }
         
        //// read and show the tabular data 
        // while (std::getline(data,line)) {
        //     std::cout << line << endl;
        // }

        // create the storing file and a tree
        std::unique_ptr<TFile> myFile( TFile::Open("file.root", "RECREATE") );
        auto tree = std::make_unique<TTree>("tree", "The Tree Title");

        // add branches 

        //// use the following if you want to take the title of the branch from the file as the variable x (char *)
        // char *x="VOLTAGE";
        // tree->Branch(Form("%s",x), &volt);
        
        float voltage;
        tree->Branch("voltage", &voltage);

        int channel;
        tree->Branch("channel", &channel);

        float cs;
        tree->Branch("cs", &cs);

        float cs_err;
        tree->Branch("cs_err", &cs_err);
        
        float tot_curr;
        tree->Branch("tot_curr", &tot_curr);

        float act_vlt;
        tree->Branch("act_vlt", &act_vlt);

        float time;
        tree->Branch("time", &time);

        float temp;
        tree->Branch("temp", &temp);

        float hum;
        tree->Branch("hum", &hum);

        float cp;
        tree->Branch("cp", &cp);

        float cp_err;
        tree->Branch("cp_err", &cp_err);

        float impedance;
        tree->Branch("impedance", &impedance);

        float impedance_err;
        tree->Branch("impedance_err", &impedance_err);

        float phase;
        tree->Branch("phase", &phase);

        float phase_err;
        tree->Branch("phase_err", &phase_err);

        float cs_uncorr;
        tree->Branch("cs_uncorr", &cs_uncorr); //uncorrected

        float cp_uncorr;
        tree->Branch("cp_uncorr", &cp_uncorr);


        // create the tree for analysis
        auto tree_anl = std::make_unique<TTree>("tree_analysis", "Tree analysis");

        int n_ch = 8; // number of channels

        // create a vector of channels
        std::vector <int> used_ch = {1,48,57,65,120,193,241,256};

        // define the branches;
        tree_anl->Branch("voltage", &voltage);

        std::vector<float> cs_anl(n_ch);
        tree_anl->Branch("cs", &cs_anl);

        std::vector<float> cs_err_anl(n_ch);
        tree_anl->Branch("cs_err", &cs_err_anl);

        std::vector<float> tot_curr_anl(n_ch);
        tree_anl->Branch("tot_curr", &tot_curr_anl);

        std::vector<float> act_vlt_anl(n_ch);
        tree_anl->Branch("act_vlt", &act_vlt_anl);

        std::vector<float> time_anl(n_ch);
        tree_anl->Branch("time", &time_anl);

        std::vector<float> temp_anl(n_ch);
        tree_anl->Branch("temp", &temp_anl);

        std::vector<float> hum_anl(n_ch);
        tree_anl->Branch("hum", &hum_anl);

        std::vector<float> cp_anl(n_ch);
        tree_anl->Branch("cp", &cp_anl);

        std::vector<float> cp_err_anl(n_ch);
        tree_anl->Branch("cp_err", &cp_err_anl);

        std::vector<float> impedance_anl(n_ch);
        tree_anl->Branch("impedance", &impedance_anl);

        std::vector<float> impedance_err_anl(n_ch);
        tree_anl->Branch("impedance_err", &impedance_err_anl);

        std::vector<float> phase_anl(n_ch);
        tree_anl->Branch("phase", &phase_anl);

        std::vector<float> phase_err_anl(n_ch);
        tree_anl->Branch("phase_err", &phase_err_anl);

        std::vector<float> cs_uncorr_anl(n_ch);
        tree_anl->Branch("cs_uncorr", &cs_uncorr_anl); // uncorrected

        std::vector<float> cp_uncorr_anl(n_ch);
        tree_anl->Branch("cp_uncorr", &cp_uncorr_anl);

        // mean values for the temperature and humidity in one voltage
        float mean_temp;
        tree_anl->Branch("mean_temp", &mean_temp);

        float mean_hum;
        tree_anl->Branch("mean_hum", &mean_hum);

        //standard deviation for the temperature and humidity in one voltage
        float std_temp;
        tree_anl->Branch("std_temp", &std_temp);

        float std_hum;
        tree_anl->Branch("std_hum", &std_hum);


        float dummy; // to jump over columns we don't want to store 

        while (true){
            bool success = true;

            for (int i=0; i<n_ch; i++) { 
                // read the tabular data
                // the data is separated by tabs, so we can use >> to read it
                if (!(data >> voltage >> channel >> cs >> cs_err >> tot_curr >> act_vlt 
                    >> time >> temp >> hum >> cp >> cp_err >> impedance >> impedance_err 
                    >> phase >> phase_err >> cs_uncorr >> cp_uncorr)) {
                  success = false;
                  break;
              }

                tree->Fill();

                cs_anl[i] = cs;
                cs_err_anl[i] = cs_err;
                tot_curr_anl[i] = tot_curr;
                act_vlt_anl[i] = act_vlt; 
                time_anl[i] = time;
                temp_anl[i] = temp;
                hum_anl[i] = hum;
                cp_anl[i] = cp;
                cp_err_anl[i] = cp_err;
                impedance_anl[i] = impedance;
                impedance_err_anl[i] = impedance_err;
                phase_anl[i] = phase;
                phase_err_anl[i] = phase_err;
                cs_uncorr_anl[i] = cs_uncorr;
                cp_uncorr_anl[i] = cp_uncorr;
            }

            if (!success) break;

            // calculate the mean and standard deviation for temperature and humidity
            mean_temp = 0;
            mean_hum = 0;
            std_temp = 0;
            std_hum = 0;
            for (int i=0; i<n_ch; i++) {
                mean_temp += temp_anl[i];
                mean_hum += hum_anl[i];
            }
            mean_temp /= n_ch;
            mean_hum /= n_ch;
            for (int i=0; i<n_ch; i++) {
                std_temp += (temp_anl[i] - mean_temp) * (temp_anl[i] - mean_temp);
                std_hum += (hum_anl[i] - mean_hum) * (hum_anl[i] - mean_hum);
            }
            std_temp = sqrt(std_temp / n_ch);
            std_hum = sqrt(std_hum / n_ch);
            
            tree_anl->Fill();

            // stop reading if we reach the end of the file
            if (data.eof()) break;
            if (!data) {
                // Optional: handle error (e.g. I/O error)
                break;
            }

        }
        
        tree_anl->Scan();

        tree->Write();
        tree_anl->Write();
        // myFile->Close();
        // data.close();
        std::cout << "Storing file created successfully." << std::endl;
    }
 

    else {
        cout << "File could not be opened" << endl;
    }    


    return 0;
} 