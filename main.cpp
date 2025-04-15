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

        float cs_error;
        tree->Branch("cs_error", &cs_error);
        
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

        float cp_error;
        tree->Branch("cp_error", &cp_error);

        float impedance;
        tree->Branch("impedance", &impedance);

        float impedance_error;
        tree->Branch("impedance_error", &impedance_error);

        float phase;
        tree->Branch("phase", &phase);

        float phase_error;
        tree->Branch("phase_error", &phase_error);

        float cs_uncorr;
        tree->Branch("cs_uncorr", &cs_uncorr); //uncorrected

        float cp_uncorr;
        tree->Branch("cp_uncorr", &cp_uncorr);
        
     
        // read the tabular data
        // the data is separated by tabs, so we can use >> to read it
        while (data >> voltage >> channel >> cs >> cs_error >> tot_curr >> act_vlt 
            >> time >> temp >> hum >> cp >> cp_error >> impedance >> impedance_error 
            >> phase >> phase_error >> cs_uncorr >> cp_uncorr) {
            
            // Now write the header
            tree->Fill();
        }
        tree->Scan("voltage:cs_error");
        tree->Show(0);
        tree->Write();


        // create the tree for analysis
        auto tree_analy = std::make_unique<TTree>("tree_analysis", "Tree analysis");

        int n_ch = 8; // number of channels

        // create a vector of channels
        std::vector <int> used_ch = {1,48,57,65,120,193,241,256};

        // define the branches;
        tree->Branch("voltage", &voltage);

        std::vector <float> cs_anl
        tree->Branch("cs", &cs);

        float cs_error;
        tree->Branch("cs_error", &cs_error);
        
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

        float cp_error;
        tree->Branch("cp_error", &cp_error);

        float impedance;
        tree->Branch("impedance", &impedance);

        float impedance_error;
        tree->Branch("impedance_error", &impedance_error);

        float phase;
        tree->Branch("phase", &phase);

        float phase_error;
        tree->Branch("phase_error", &phase_error);

        float cs_uncorr;
        tree->Branch("cs_uncorr", &cs_uncorr); //uncorrected

        float cp_uncorr;
        tree->Branch("cp_uncorr", &cp_uncorr);

        float dummy; // to jump over columns we don't want to store 

        while (true){
            data >> voltage >> dummy;

            for (int i=0; i<n_ch; i++) {
                cs >> cs_error >> tot_curr >> act_vlt 
                >> time >> temp >> hum >> cp >> cp_error >> impedance >> impedance_error 
                >> phase >> phase_error >> cs_uncorr >> cp_uncorr
            }

        }



        while (data >> voltage >> channel >> cs >> cs_error >> tot_curr >> act_vlt 
            >> time >> temp >> hum >> cp >> cp_error >> impedance >> impedance_error 
            >> phase >> phase_error >> cs_uncorr >> cp_uncorr) {
            
            // Now write the header
            tree_analy->Fill();
        }
        tree_analy->Scan("voltage:cs_error");
        tree_analy->Show(0);
        tree_analy->Write();





    }
 

    else {
        cout << "No se pudo abrir el archivo." << endl;
    }    

    return 0;
} 