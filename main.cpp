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
        while (std::getline(data, line)) {
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
    float volt;
    tree->Branch("Voltage", &volt);


    // store data
    for (int iEntry = 0; iEntry < 1000; ++iEntry) {
        volt = 0.3 * iEntry;
        // Fill the current value of `volt` into `Voltage`
        tree->Fill();
     }
     
     // Now write the header
     tree->Write();




    return 0;
} 