/**
 * HDecayBuilder.h
 *
 *
 */

#ifndef HDECAYBUILDER_H
#define HDECAYBUILDER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// framework includes
#include "KFitParticle.h"
#include "hkinfitter.h"
//#include "hvertexfinder.h"
//#include "hneutralcandfinder.h"

//#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"

#include <algorithm>
#include <iterator>

using std::cout;
using std::endl;

class HDecayBuilder
{
private:
    // Working Particles
    std::vector<KFitParticle> fFitCands;
    std::vector<std::vector<KFitParticle>> fCands;

    // Output particles after fitting
    std::vector<KFitParticle> fOutputCands;

    // Fitter input variables
    TString fTask;
    std::vector<Int_t> fPids;
    //std::vector<int> fPidsPrim;
    //std::vector<int> fPidsDecay;
    TLorentzVector fIniSys;
    KFitParticle fMother;
    Double_t fMass;

    // For combinatorics
    Int_t fTotalCombos;
    Int_t fCombiCounter;
    std::vector<Int_t> particleCounter;
    Bool_t doubleParticle;

    // Probability
    Double_t fProb;
    Double_t fBestProb;
    Double_t fBestChi2;

    Int_t fVerbose;

public:
    HDecayBuilder(std::vector<std::vector<KFitParticle>> &cands, TString &task, std::vector<Int_t> &pids, TLorentzVector lv = TLorentzVector(), KFitParticle mother = KFitParticle(), Double_t mass = 0.);
    HDecayBuilder(TString &task, std::vector<Int_t> &pids, TLorentzVector lv = TLorentzVector(), KFitParticle mother = KFitParticle(), Double_t mass = 0.);
    ~HDecayBuilder(){};

    void setVerbosity(int val) { fVerbose = val; }

    // setters
    void setInputCands(std::vector<std::vector<KFitParticle>> cands) {fCands = cands; }
    void setIniSys(TLorentzVector val) { fIniSys = val; }
    void setMother(KFitParticle val) { fMother = val; }
    void setMass(Double_t val) { fMass = val; }
/*
    void setPidsInVertices(std::vector<int> val1, std::vector<int> val2)
    {
        fPidsPrim = val1;
        fPidsDecay = val2;
    }*/

    void buildDecay();

    //void createNeutralCandidate();

    Bool_t doFit();

    void countCombis();
    void fillFitCands();
    void checkDoubleParticle(size_t i);

    void createOutputParticle(KFitParticle FittedCand);
    void getFitCands(std::vector<KFitParticle> &cands) { cands = fOutputCands; }
    Double_t getChi2() { return fBestChi2; }
    Double_t getProbability() { return fBestProb; }
};

#endif /* HDECAYBUILDER_H */
