/**
 * KFitDecayBuilder.h
 *
 *
 */

#ifndef KFITDECAYBUILDER_H
#define KFITDECAYBUILDER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// framework includes
#include "KFitParticle.h"
#include "KinFitter.h"
//#include "KFitVertexFinder.h"
//#include "hneutralcandfinder.h"

//#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"

#include <algorithm>
#include <iterator>

using std::cout;
using std::endl;

class KFitDecayBuilder
{
private:
    // Working Particles
    std::vector<KFitParticle> fFitCands;
    std::vector<std::vector<KFitParticle>> fCands;

    // Output particles after fitting
    std::vector<KFitParticle> fOutputCands;

    // Fitter input variables
    TString fTask;
    std::vector<int> fPids;
    //std::vector<int> fPidsPrim;
    //std::vector<int> fPidsDecay;
    TLorentzVector fIniSys;
    KFitParticle fMother;
    double fMass;

    // For combinatorics
    int fTotalCombos;
    int fCombiCounter;
    std::vector<int> particleCounter;
    bool doubleParticle;

    // Probability
    double fProb;
    double fBestProb;
    double fBestChi2;

    int fVerbose;

public:
    KFitDecayBuilder(std::vector<std::vector<KFitParticle>> &cands, TString &task, std::vector<int> &pids, TLorentzVector lv = TLorentzVector(), KFitParticle mother = KFitParticle(), double mass = 0.);
    KFitDecayBuilder(TString &task, std::vector<int> &pids, TLorentzVector lv = TLorentzVector(), KFitParticle mother = KFitParticle(), double mass = 0.);
    ~KFitDecayBuilder(){};

    void setVerbosity(int val) { fVerbose = val; }

    // setters
    void setInputCands(std::vector<std::vector<KFitParticle>> cands) {fCands = cands; }
    void setIniSys(TLorentzVector val) { fIniSys = val; }
    void setMother(KFitParticle val) { fMother = val; }
    void setMass(double val) { fMass = val; }
/*
    void setPidsInVertices(std::vector<int> val1, std::vector<int> val2)
    {
        fPidsPrim = val1;
        fPidsDecay = val2;
    }*/

    void buildDecay();

    //void createNeutralCandidate();

    bool doFit();

    void countCombis();
    void fillFitCands();
    void checkDoubleParticle(size_t i);

    void createOutputParticle(KFitParticle FittedCand);
    void getFitCands(std::vector<KFitParticle> &cands) { cands = fOutputCands; }
    double getChi2() { return fBestChi2; }
    double getProbability() { return fBestProb; }
};

#endif /* KFITDECAYBUILDER_H */
