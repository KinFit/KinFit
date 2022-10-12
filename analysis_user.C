#include "TLorentzVector.h"

#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <iomanip>
#include <math.h>


#include "hrootfitter.h"

using namespace std;

Int_t analysis_user(TString infileList="input.root", TString outfile = "fitted.root", Int_t nEvents=10){

    HRootFitter RootFitter(infileList, outfile, nEvents);
    std::vector<int> pids;
    //std::vector<int> pidsPrimVertex; // Jenny: User must set this for 3C fit
    //stc::vector<int> pidsDecayVertex; // Jenny: User must set this for 3C fit
    //pids.push_back(14); pids.push_back(11); pids.push_back(14); pids.push_back(9);
    pids.push_back(14); pids.push_back(9);
    //TLorentzVector ppSystem(0,0,4337.96,2*938.272+3500);
    Double_t mass = 0.938272;

    //RootFitter.addFitterTask("4c", pids, ppSystem);
    RootFitter.doFitterTask("Mass", pids, mass);
    RootFitter.finish();
    
    return 0;
}
