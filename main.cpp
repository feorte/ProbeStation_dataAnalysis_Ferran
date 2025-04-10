#include <iostream>
#include <TROOT.h>
#include <fstream>
#include <string>

int main()
{
    std::ifstream CV_data ("CALICE_6in_256ch_77_20250221_1_CV_customscan.txt"); // open the file directly when initializing the stream object
    if (CV_data.is_open()) {
        std::string line;
        
        // read first the header lines which are not tabular data 
        for (int i=0; i<28; i++) {
            std::getline(CV_data,line);
        }

        // read the tabular data 
        while (std::getline(CV_data, line)) {
            std::cout << line << endl;
        }
    }

    else {
        cout << "No se pudo abrir el archivo." << endl;
    }    

    return 0;
} 