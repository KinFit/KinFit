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

Int_t analysis_user(TString infileList="~/KinFit/input.root", TString outfile = "fitted_mass.root", Int_t nEvents=10){

    KFitRootAnalyzer RootAnalyzer(infileList, outfile, nEvents);
    std::vector<int> pids;
    //std::vector<int> pidsPrimVertex; // Jenny: User must set this for 3C fit
    //stc::vector<int> pidsDecayVertex; // Jenny: User must set this for 3C fit
    //pids.push_back(14); pids.push_back(11); pids.push_back(14); pids.push_back(9);
    pids.push_back(14); pids.push_back(9);
    //TLorentzVector ppSystem(0,0,4337.96,2*938.272+3500);
    Double_t mass = 1.11568;

    //RootFitter.addFitterTask("4c", pids, ppSystem);
    RootAnalyzer.doFitterTask("Mass", pids, mass);
    //cout<<"fitter task performed"<<endl;
    //RootFitter.finish();
    //cout<<"finished"<<endl;
    
    return 0;
}
