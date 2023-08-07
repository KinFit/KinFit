#include "KFitRootAnalyzer.h"

Int_t analysis_user(TString infileList="~/KinFit/input.root", TString outfile = "fitted_test.root", Int_t nEvents=1000){

    // Initialize Analyzer 
    // Arguments are input file(list), output file name and number of events to analyze
    KFitRootAnalyzer RootAnalyzer(infileList, outfile, nEvents);

    // Create vector that contains the PIDs of the particles that should be analyzed
    std::vector<int> pids;
    //pids.push_back(14); pids.push_back(11); pids.push_back(14); pids.push_back(9);
    pids.push_back(14); pids.push_back(14); pids.push_back(9);
    //pids.push_back(14); pids.push_back(9);

    // Other input that might be needed, depending on constraint
    TLorentzVector ppSystem(0,0,5.3567,2*0.938272+4.500);
    Double_t mLambda = 1.11568;
    Double_t mK = 0.493677;

    //RootAnalyzer.doFitterTask("4C", pids, -1, ppSystem);
    //RootAnalyzer.doFitterTask("Mass", pids, mLambda);  
    //RootAnalyzer.doFitterTask("MassVtx", pids, mLambda);  //no candidate found
    //RootAnalyzer.doFitterTask("MM", pids, mK, ppSystem);  // no candidate found
    RootAnalyzer.doFitterTask("Mom", pids, mK, ppSystem); 
    //RootAnalyzer.doFitterTask("Vertex", pids, -1);    // symbol lookup error in no mLambda is given??
    
    return 0;
}
