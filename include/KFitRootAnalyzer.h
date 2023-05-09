/**
 * KFitRootAnalyzer.h
 *
 *
 */

#ifndef KFITROOTANALYZER_H
#define KFITROOTANALYZER_H

// system includes
#include <iostream>
#include <vector>
#include <cmath>

// framework includes
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "KFitParticle.h"
#include "KFitDecayBuilder.h"

using std::cout;
using std::endl;

class KFitRootAnalyzer
{
private:
    //TString fInfileList, fOutprefix;
    std::vector< std::vector<KFitParticle> > fCandsFit;
    std::vector<int> fPids;
    int fEvents;

    int fVerbose;

    //Read in data
    TFile *finFile;
    TTree *fTree;
    TClonesArray *fCands_in = new TClonesArray("KFitParticle");

    //Output data
    TFile *foutFile;
    TTree *fTree_out;

    // Fill KFitParticle with necessary information
    //void FillData(TLorentzVector *cand, double R, double Z, KFitParticle &outcand, TMatrix_D cov, double mass);
    // put to KFitParticle, here cands should be fully filled already
    
public:
    KFitRootAnalyzer(TString inFileName, TString outFileName, int nEvents=-1);
    ~KFitRootAnalyzer(){};

    //User functions
    //void doFitterTask(TString task, std::vector<int> pids, TLorentzVector lv = TLorentzVector(), KFitParticle mother = KFitParticle(), double mm=0.);
    void doFitterTask(TString task, std::vector<int> pids, double mm=0., TLorentzVector lv = TLorentzVector(), KFitParticle mother = KFitParticle());
    //void addFitterTask(TString task, std::vector<int> primPids, std::vector<int> decayPids); // Jenny, for 3C fit
    //void addBuilderTask(TString task, std::vector<int> pids, TLorentzVector lv);

    void setPids(std::vector<int> val){ fPids = val; }
    void setVerbosity(int val){ fVerbose = val; }

    std::vector<int> getPids(){ return fPids; }
    TTree* getFittedTree(){ return fTree_out; }

    // select and sort candidates according to their PID
    void selectCandidates();

    void finish(); //write output
};

#endif /* KFITROOTANALYZER_H */