#include <iostream>
#include "TROOT.h"
#include <fstream>
#include <string>
#include <vector> 


int alt_read_store()
{
    // read the data
    string filename = "CALICE_6in_256ch_77_20250509_3_CV.txt"; // replace with your file name
    std::ifstream data (filename); // open the file directly when initializing the stream object
    cout<<"Reading file: " << filename << endl;

    std::vector<int> customChannels; 
    std::string measurementType;

    if (data.is_open()) {
        std::string line;

        // Read the first 30 lines (these include headers and metadata)
        for (int i = 0; i < 30; ++i) {
            std::getline(data, line);

             // Check first line for measurement type (IV or CV)
             if (i == 0) {
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    measurementType = line.substr(pos + 1);
                    // Trim whitespace
                    measurementType.erase(0, measurementType.find_first_not_of(" \t\r\n"));
                    measurementType.erase(measurementType.find_last_not_of(" \t\r\n") + 1);
                    std::cout << "Measurement type: " << measurementType << std::endl;
                }
            }

            // Line 16 (index 15) contains the custom channels info
            if (i == 15) {
                size_t pos = line.find(":");
                if (pos != std::string::npos) {
                    // Extract everything after the colon
                    std::string channelsStr = line.substr(pos + 1);

                    // Trim leading and trailing whitespace
                    channelsStr.erase(0, channelsStr.find_first_not_of(" \t\r\n"));
                    channelsStr.erase(channelsStr.find_last_not_of(" \t\r\n") + 1); 

                    if (channelsStr == "None") {
                        // If "None", that means all 256 channels are used
                        for (int ch = 1; ch <= 256; ++ch) {
                            customChannels.push_back(ch);
                        }
                    } else {
                        // Otherwise, parse the comma-separated list of custom channels
                        std::stringstream ss(channelsStr);
                        std::string used_channels;

                        while (std::getline(ss, used_channels, ',')) {
                            int ch = std::stoi(used_channels); // convert string to int
                            customChannels.push_back(ch);
                        }
                    }
                }
            }
        }

        // Optional: print the list of channels for verification
        std::cout << "Custom channels (" << customChannels.size() << "): ";
        for (int ch : customChannels) {
            std::cout << ch << " ";
        }
        std::cout << std::endl;

        int n_ch = customChannels.size(); // number of channels

        // create the storing file and a tree to store the data
        std::unique_ptr<TFile> myFile( TFile::Open("stored_data.root", "RECREATE") );
        auto raw_measurements = std::make_unique<TTree>("raw_measurements", "Raw measurements");

        // create the tree for the analysis
        auto analysis = std::make_unique<TTree>("analysis", "Analysis");

        // add branches present in both CV and IV for raw measurements...
        float voltage;
        raw_measurements->Branch("voltage", &voltage);

        int channel;
        raw_measurements->Branch("channel", &channel);

        float tot_curr;
        raw_measurements->Branch("tot_curr", &tot_curr);

        float act_vlt;
        raw_measurements->Branch("act_vlt", &act_vlt);

        float time;
        raw_measurements->Branch("time", &time);

        float temp;
        raw_measurements->Branch("temp", &temp);

        float hum;
        raw_measurements->Branch("hum", &hum);

        //... and for analysis data
        analysis->Branch("voltage", &voltage);

        analysis->Branch("channel", &customChannels);

        std::vector<float> tot_curr_anl(n_ch);
        analysis->Branch("tot_curr", &tot_curr_anl);

        std::vector<float> act_vlt_anl(n_ch);
        analysis->Branch("act_vlt", &act_vlt_anl);

        std::vector<float> time_anl(n_ch);
        analysis->Branch("time", &time_anl);

        std::vector<float> temp_anl(n_ch);
        analysis->Branch("temp", &temp_anl);

        std::vector<float> hum_anl(n_ch);
        analysis->Branch("hum", &hum_anl);

        // mean values for the temperature and humidity in one voltage
        float mean_temp;
        analysis->Branch("mean_temp", &mean_temp);

        float mean_hum;
        analysis->Branch("mean_hum", &mean_hum);

        //standard deviation for the temperature and humidity in one voltage
        float std_temp;
        analysis->Branch("std_temp", &std_temp);

        float std_hum;
        analysis->Branch("std_hum", &std_hum);
 

        // Specific processing based on measurement type
        if (measurementType == "IV") {
            std::cout << "Processing IV data...\n";

            float current;
            raw_measurements->Branch("current", &current);
    
            float current_err;
            raw_measurements->Branch("current_err", &current_err);

            std::vector<float> current_anl(n_ch);
            analysis->Branch("current", &current_anl);

            std::vector<float> current_err_anl(n_ch);
            analysis->Branch("current_err", &current_err_anl);

            float dummy; // to jump over columns we don't want to store 

            while (true){
                bool success = true;
    
                for (int i=0; i<n_ch; i++) { 
                    // read the tabular data
                    // the data is separated by tabs, so we can use >> to read it
                    if (!(data >> voltage >> channel >> current >> current_err >> tot_curr >> act_vlt 
                        >> time >> temp >> hum >> dummy >> dummy >> dummy >> dummy 
                        >> dummy >> dummy >> dummy >> dummy)) {
                      success = false;
                      break;
                  }
    
                    // Fill measurements into the tree
                    raw_measurements->Fill();
                    analysis->Fill();

                    current_anl[i] = current;
                    current_err_anl[i] = current_err;
                    tot_curr_anl[i] = tot_curr;
                    act_vlt_anl[i] = act_vlt; 
                    time_anl[i] = time;
                    temp_anl[i] = temp;
                    hum_anl[i] = hum;
                }  
                if (!success) break;
            }      

        } else if (measurementType == "CV") {
            std::cout << "Processing CV data...\n";

            float cs;
            raw_measurements->Branch("cs", &cs);
    
            float cs_err;
            raw_measurements->Branch("cs_err", &cs_err);

            float cp;
            raw_measurements->Branch("cp", &cp);

            float cp_err;
            raw_measurements->Branch("cp_err", &cp_err);

            float impedance;
            raw_measurements->Branch("impedance", &impedance);

            float impedance_err;
            raw_measurements->Branch("impedance_err", &impedance_err);

            float phase;
            raw_measurements->Branch("phase", &phase);

            float phase_err;
            raw_measurements->Branch("phase_err", &phase_err);

            float cs_uncorr;
            raw_measurements->Branch("cs_uncorr", &cs_uncorr); //uncorrected

            float cp_uncorr;
            raw_measurements->Branch("cp_uncorr", &cp_uncorr);

            std::vector<float> cs_anl(n_ch);
            analysis->Branch("cs", &cs_anl);

            std::vector<float> cs_err_anl(n_ch);
            analysis->Branch("cs_err", &cs_err_anl);

            std::vector<float> cp_anl(n_ch);
            analysis->Branch("cp", &cp_anl);

            std::vector<float> cp_err_anl(n_ch);
            analysis->Branch("cp_err", &cp_err_anl);

            std::vector<float> impedance_anl(n_ch);
            analysis->Branch("impedance", &impedance_anl);

            std::vector<float> impedance_err_anl(n_ch);
            analysis->Branch("impedance_err", &impedance_err_anl);

            std::vector<float> phase_anl(n_ch);
            analysis->Branch("phase", &phase_anl);

            std::vector<float> phase_err_anl(n_ch);
            analysis->Branch("phase_err", &phase_err_anl);

            std::vector<float> cs_uncorr_anl(n_ch);
            analysis->Branch("cs_uncorr", &cs_uncorr_anl); // uncorrected

            std::vector<float> cp_uncorr_anl(n_ch);
            analysis->Branch("cp_uncorr", &cp_uncorr_anl);

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
    
                    // Fill measurements into the tree
                    raw_measurements->Fill();
                    analysis->Fill();

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

            } 

        } else {
            std::cerr << "Unknown measurement type: " << measurementType << std::endl;
        }
        

        // add branches 

        //// use the following if you want to take the title of the branch from the file as the variable x (char *)
        // char *x="VOLTAGE";
        // raw_measurements->Branch(Form("%s",x), &volt);
        
        // float voltage;
        // raw_measurements->Branch("voltage", &voltage);

        // int channel;
        // raw_measurements->Branch("channel", &channel);

        // float cs;
        // raw_measurements->Branch("cs", &cs);

        // float cs_err;
        // raw_measurements->Branch("cs_err", &cs_err);
        
        // float tot_curr;
        // raw_measurements->Branch("tot_curr", &tot_curr);

        // float act_vlt;
        // raw_measurements->Branch("act_vlt", &act_vlt);

        // float time;
        // raw_measurements->Branch("time", &time);

        // float temp;
        // raw_measurements->Branch("temp", &temp);

        // float hum;
        // raw_measurements->Branch("hum", &hum);

        // float cp;
        // raw_measurements->Branch("cp", &cp);

        // float cp_err;
        // raw_measurements->Branch("cp_err", &cp_err);

        // float impedance;
        // raw_measurements->Branch("impedance", &impedance);

        // float impedance_err;
        // raw_measurements->Branch("impedance_err", &impedance_err);

        // float phase;
        // raw_measurements->Branch("phase", &phase);

        // float phase_err;
        // raw_measurements->Branch("phase_err", &phase_err);

        // float cs_uncorr;
        // raw_measurements->Branch("cs_uncorr", &cs_uncorr); //uncorrected

        // float cp_uncorr;
        // raw_measurements->Branch("cp_uncorr", &cp_uncorr);


        // // create the tree for the analysis
        // auto analysis = std::make_unique<TTree>("analysis", "Analysis");

        // // define the branches;
        // analysis->Branch("voltage", &voltage);

        // analysis->Branch("channel", &customChannels);

        // std::vector<float> cs_anl(n_ch);
        // analysis->Branch("cs", &cs_anl);

        // std::vector<float> cs_err_anl(n_ch);
        // analysis->Branch("cs_err", &cs_err_anl);

        // std::vector<float> tot_curr_anl(n_ch);
        // analysis->Branch("tot_curr", &tot_curr_anl);

        // std::vector<float> act_vlt_anl(n_ch);
        // analysis->Branch("act_vlt", &act_vlt_anl);

        // std::vector<float> time_anl(n_ch);
        // analysis->Branch("time", &time_anl);

        // std::vector<float> temp_anl(n_ch);
        // analysis->Branch("temp", &temp_anl);

        // std::vector<float> hum_anl(n_ch);
        // analysis->Branch("hum", &hum_anl);

        // std::vector<float> cp_anl(n_ch);
        // analysis->Branch("cp", &cp_anl);

        // std::vector<float> cp_err_anl(n_ch);
        // analysis->Branch("cp_err", &cp_err_anl);

        // std::vector<float> impedance_anl(n_ch);
        // analysis->Branch("impedance", &impedance_anl);

        // std::vector<float> impedance_err_anl(n_ch);
        // analysis->Branch("impedance_err", &impedance_err_anl);

        // std::vector<float> phase_anl(n_ch);
        // analysis->Branch("phase", &phase_anl);

        // std::vector<float> phase_err_anl(n_ch);
        // analysis->Branch("phase_err", &phase_err_anl);

        // std::vector<float> cs_uncorr_anl(n_ch);
        // analysis->Branch("cs_uncorr", &cs_uncorr_anl); // uncorrected

        // std::vector<float> cp_uncorr_anl(n_ch);
        // analysis->Branch("cp_uncorr", &cp_uncorr_anl);

        // // mean values for the temperature and humidity in one voltage
        // float mean_temp;
        // analysis->Branch("mean_temp", &mean_temp);

        // float mean_hum;
        // analysis->Branch("mean_hum", &mean_hum);

        // //standard deviation for the temperature and humidity in one voltage
        // float std_temp;
        // analysis->Branch("std_temp", &std_temp);

        // float std_hum;
        // analysis->Branch("std_hum", &std_hum);


        // float dummy; // to jump over columns we don't want to store 

        // while (true){
        //     bool success = true;

        //     for (int i=0; i<n_ch; i++) { 
        //         // read the tabular data
        //         // the data is separated by tabs, so we can use >> to read it
        //         if (!(data >> voltage >> channel >> cs >> cs_err >> tot_curr >> act_vlt 
        //             >> time >> temp >> hum >> cp >> cp_err >> impedance >> impedance_err 
        //             >> phase >> phase_err >> cs_uncorr >> cp_uncorr)) {
        //           success = false;
        //           break;
        //       }

        //         // Fill measurements into the tree
        //         raw_measurements->Fill();

        //         cs_anl[i] = cs;
        //         cs_err_anl[i] = cs_err;
        //         tot_curr_anl[i] = tot_curr;
        //         act_vlt_anl[i] = act_vlt; 
        //         time_anl[i] = time;
        //         temp_anl[i] = temp;
        //         hum_anl[i] = hum;
        //         cp_anl[i] = cp;
        //         cp_err_anl[i] = cp_err;
        //         impedance_anl[i] = impedance;
        //         impedance_err_anl[i] = impedance_err;
        //         phase_anl[i] = phase;
        //         phase_err_anl[i] = phase_err;
        //         cs_uncorr_anl[i] = cs_uncorr;
        //         cp_uncorr_anl[i] = cp_uncorr;
        //     }

        //     if (!success) break;

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

            // fill the analysis tree
            analysis->Fill();

            // // stop reading if we reach the end of the file
            // if (data.eof()) {break;}
            // if (!data) {
            //     // Optional: handle error (e.g. I/O error)
            //     break;
            // }

        
        
        //analysis->Scan();

        raw_measurements->Write();
        analysis->Write();
        // myFile->Close();
        // data.close();
        std::cout << "Storing file created successfully." << std::endl;
    }
 

    else {
        cout << "File could not be opened" << endl;
    }    


    return 0;
} 