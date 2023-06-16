//#include "TLorentzVector.h"
/*
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <iomanip>
#include <math.h>
*/

#include "KFitRootAnalyzer.h"

//using namespace std;

Int_t analysis_user(TString infileList="~/KinFit/input.root", TString outfile = "fitted_test.root", Int_t nEvents=1000){

    KFitRootAnalyzer RootAnalyzer(infileList, outfile, nEvents);
    std::vector<int> pids;
    //std::vector<int> pidsPrimVertex; // Jenny: User must set this for 3C fit
    //stc::vector<int> pidsDecayVertex; // Jenny: User must set this for 3C fit
    //pids.push_back(14); pids.push_back(11); pids.push_back(14); pids.push_back(9);
    pids.push_back(14); pids.push_back(14); pids.push_back(9);
    //pids.push_back(14); pids.push_back(9);
    TLorentzVector ppSystem(0,0,5.3567,2*0.938272+4.500);
    Double_t mass = 1.11568;
    Double_t mK = 0.493677;

    //RootAnalyzer.doFitterTask("4C", pids, -1, ppSystem);
    //RootAnalyzer.doFitterTask("Mass", pids, mass);  
    //RootAnalyzer.doFitterTask("MassVtx", pids, mass);  //no candidate found
    //RootAnalyzer.doFitterTask("MM", pids, mK, ppSystem);  // no candidate found
    RootAnalyzer.doFitterTask("Mom", pids, mK, ppSystem); 
    //RootAnalyzer.doFitterTask("Vertex", pids, -1);    // symbol lookup error in no mass is given??
    //cout<<"fitter task performed"<<endl;
    //RootFitter.finish();
    //cout<<"finished"<<endl;
    
    return 0;
}
