#include <iostream>
#include <TROOT.h>
#include <fstream>
#include <string>

int main()
{
    std::ifstream CV_data ("CALICE_6in_256ch_77_20250221_1_CV_customscan.txt");
    if (!CV_data.is_open()) {
        cout << "No se pudo abrir el archivo." << endl;
    }
    
    string line;
    while (getline(CV_data, line)) {
        cout << line << endl;
}
    



}