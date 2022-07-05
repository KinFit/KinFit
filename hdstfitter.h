/**
 * HDSTFitter.h
 *
 *
 */

#ifndef HDSTFITTER_H
#define HDSTFITTER_H

// system includes
#include <iostream>
#include <vector>
#include <cmath>

// framework includes
//#include "TLorentzVector.h"

#include "hades.h"
#include "hcategorymanager.h"
#include "hloop.h"
#include "htool.h"

#include "hcategory.h"
#include "hparticlecandsim.h"

#include "hgeantkine.h"
#include "hparticlegeant.h"
#include "hparticlegeantdecay.h"
#include "hparticlegeantevent.h"

#include "hdecaybuilder.h"
#include "hcovariancekinfit.h"

using std::cout;
using std::endl;

class HDSTFitter
{
private:
    TString fInfileList;
    bool fIncludeFw;
    std::vector< std::vector<HRefitCand> > fCandsFit;
    std::vector<Int_t> fPids;
    // Variables used for setting the covariance matrix
    bool fMomDepErrors;
    Int_t fEvents;

    Int_t fVerbose;

    //Read in data
    HCategory *fcatParticle;
    HCategory *fcatFwParticle;

    // Method to fill the data from a HRefitCand for simulations
    void FillData(HParticleCandSim *cand, HRefitCand &outcand, double arr[5], double mass);
    //void FillDataFW(HFwDetCandSim *cand, HRefitCand *outcand, double arr[], double mass); //adjust to HForwardCand for newer Hydra
    
public:
    HDSTFitter(TString infileList, bool includeFw = false, bool momDepErrors=false, Int_t nEvents=-1);
    ~HDSTFitter(){};

    //User functions
    void addFitterTask(TString task, std::vector<Int_t> pids, TLorentzVector lv = TLorentzVector(), HRefitCand mother = HRefitCand(), Double_t mm=0.);
    void addFitterTask(TString task, std::vector<Int_t> primPids, std::vector<Int_t> decayPids); // Jenny, for 3C fit

    void addBuilderTask(TString task, std::vector<Int_t> pids, TLorentzVector lv);

    void setIncludeFw(bool val){ fIncludeFw = val; }
    void setErrors();
    void setPids(std::vector<Int_t> val){ fPids = val; }
    void setVerbosity(Int_t val){ fVerbose = val; }

    std::vector<Int_t> getPids(){ return fPids; }


    void selectCandidates();
};

#endif /* HDSTFITTER_H */

// How to handle several tasks???
