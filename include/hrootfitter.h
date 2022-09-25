/**
 * HRootFitter.h
 *
 *
 */

#ifndef HROOTFITTER_H
#define HROOTFITTER_H

// system includes
#include <iostream>
#include <vector>
#include <cmath>

// framework includes
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "hrefitcand.h"
#include "hdecaybuilder.h"

using std::cout;
using std::endl;

class HRootFitter
{
private:
    //TString fInfileList, fOutprefix;
    std::vector< std::vector<HRefitCand> > fCandsFit;
    std::vector<Int_t> fPids;
    Int_t fEvents;

    Int_t fVerbose;

    //Read in data
    TFile *finFile;
    TTree *fTree;
    TClonesArray *fCands_in = new TClonesArray("HRefitCand");

    //Output data
    TFile *foutFile;
    TTree *fTree_out;

    // Fill HRefitCand with necessary information
    //void FillData(TLorentzVector *cand, Double_t R, Double_t Z, HRefitCand &outcand, TMatrix_D cov, double mass);
    // put to HRefitCand, here cands should be fully filled already
    
public:
    HRootFitter(TString inFileName, TString outFileName, Int_t nEvents=-1);
    ~HRootFitter(){};

    //User functions
    void addFitterTask(TString task, std::vector<Int_t> pids, TLorentzVector lv = TLorentzVector(), HRefitCand mother = HRefitCand(), Double_t mm=0.);
    //void addFitterTask(TString task, std::vector<Int_t> primPids, std::vector<Int_t> decayPids); // Jenny, for 3C fit
    //void addBuilderTask(TString task, std::vector<Int_t> pids, TLorentzVector lv);

    void setPids(std::vector<Int_t> val){ fPids = val; }
    void setVerbosity(Int_t val){ fVerbose = val; }

    std::vector<Int_t> getPids(){ return fPids; }
    TTree* getFittedTree(){ return fTree_out; }

    // select and sort candidates according to their PID
    void selectCandidates();

    void finish(); //write output
};

#endif /* HROOTFITTER_H */