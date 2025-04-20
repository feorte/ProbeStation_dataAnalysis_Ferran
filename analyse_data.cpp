#include <iostream>
#include "TROOT.h"
#include <fstream>
#include <string>
#include <vector> 


int analyse_data()
{
    // read the data in the stored tree
    std::unique_ptr<TFile> data( TFile::Open("stored_data.root") );
    auto tree = myFile->Get<TTree>("TreeName");

    int variable;
    tree->SetBranchAddress("branchName", &variable);

    for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
    // Load the data for the given tree entry
    tree->GetEntry(iEntry);

    // Now, `variable` is set to the value of the branch
    // "branchName" in tree entry `iEntry`
    printf("%d\n", variable);
    }





    return 0;
]