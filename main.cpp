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
        for (int i=0; i<28; i++) {
            std::getline(data,line);
        }
         
        // read the tabular data 
        while (std::getline(data,line)) {
            std::cout << line << endl;
        }
    }
 
    else {
        cout << "No se pudo abrir el archivo." << endl;
    }    
        
    ////////////////////////////////////////////////////////////////////////////////

    // create a tree and store it 

    // create the storing file and the tree
    std::unique_ptr<TFile> myFile( TFile::Open("file.root", "RECREATE") );
    auto tree = std::make_unique<TTree>("tree", "The Tree Title");

    // add branches 
    // use the following if you want to take the title of the branch from the file as the variable x (char *)
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

    
    // store data
    for (int iEntry = 0; iEntry < 1000; ++iEntry) {
        voltage = 0.3 * iEntry;
        // Fill the current value of `voltage` into branch `voltage`
        tree->Fill();
     }
     
     // Now write the header
     tree->Write();




    return 0;
} 