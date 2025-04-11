#include <iostream>
#include <TROOT.h>
#include <fstream>
#include <string>

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

        float capacitance_series;
        tree->Branch("capacitance_series", &capacitance_series);

        float error_capacitance_series;
        tree->Branch("error_capacitance_series", &error_capacitance_series);
        
        float total_current;
        tree->Branch("total_current", &total_current);

        float active_voltage;
        tree->Branch("active_voltage", &active_voltage);

        float time;
        tree->Branch("time", &time);

        float temperature;
        tree->Branch("temperature", &temperature);

        float humidity;
        tree->Branch("humidity", &humidity);

        float capacitance_parallel;
        tree->Branch("capacitance_parallel", &capacitance_parallel);

        float error_capacitance_parallel;
        tree->Branch("error_capacitance_parallel", &error_capacitance_parallel);

        float impedance;
        tree->Branch("impedance", &impedance);

        float impedance_error;
        tree->Branch("impedance_error", &impedance_error);

        float phase;
        tree->Branch("phase", &phase);

        float phase_error;
        tree->Branch("phase_error", &phase_error);

        float capacitance_series_uncorrelated;
        tree->Branch("capacitance_series_uncorrelated", &capacitance_series_uncorrelated);

        float capacitance_parellel_uncorrelated;
        tree->Branch("capacitance_parellel_uncorrelated", &capacitance_parellel_uncorrelated);

        
     
        // read the tabular data
        // the data is separated by tabs, so we can use >> to read it
        while (data >> voltage >> channel >> capacitance_series >> error_capacitance_series >> total_current >> active_voltage >> time >> temperature >> humidity >> capacitance_parallel >> error_capacitance_parallel >> impedance >> impedance_error >> phase >> phase_error >> capacitance_series_uncorrelated >> capacitance_parellel_uncorrelated) {
            
            // Now write the header
            tree->Fill();
        }
    }
 

    else {
        cout << "No se pudo abrir el archivo." << endl;
    }    




    return 0;
} 