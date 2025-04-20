#include <iostream>
#include "TROOT.h"
#include <fstream>
#include <string>
#include <vector> 


int analyse_data()
{

    // read the data in the stored tree
    std::unique_ptr<TFile> data( TFile::Open("stored_data.root") );
    auto tree = data->Get<TTree>("analysis");

    // Disable everything...
    tree->SetBranchStatus("*", false);
    // ...but the branches we need
    for (const auto& name : {"voltage", "channel"}) {
        tree->SetBranchStatus(name, true);
    }
    
    // Associate the branches to variables
    float voltage;
    tree->SetBranchAddress("voltage", &voltage);

    std::vector<int>* channel = nullptr;
    tree->SetBranchAddress("channel", &channel);
    


    for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
        // Load the data for the given tree entry
        tree->GetEntry(iEntry);

        // Now, `variable` is set to the value of the branch
        // "branchName" in tree entry `iEntry`
        printf("%f\n", voltage);
        printf("%d\n", channel->at(1));  // or (*channel)[0]

    }





    return 0;
}