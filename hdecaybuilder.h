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
#include "hrefitcand.h"

#include "TH1F.h"

#include "hkinfitter.h"
#include "hvertexfinder.h"
#include "hneutralcandfinder.h"

#include <algorithm>
#include <iterator>

using std::cout;
using std::endl;

class HDecayBuilder
{
private:
    // Working Particles
    std::vector<HRefitCand> fFitCands;
    std::vector< std::vector<HRefitCand> > fCands;

    // Output particles after fitting
    std::vector<HRefitCand> fOutputCands;


    //Fitter input variables
    TString fTask;
    std::vector<Int_t> fPids;
    std::vector<int> fPidsPrim;
    std::vector<int> fPidsDecay;
    TLorentzVector fIniSys;
    HRefitCand fMother;
    Double_t fMass;
    
    //For combinatorics
    Int_t fTotalCombos;
    Int_t fCombiCounter;
    std::vector<Int_t> particleCounter;
    bool doubleParticle;
    
    //Probability
    Double_t fProb;
    Double_t fBestProb;

    int fVerbose;

public:
    HDecayBuilder(std::vector< std::vector<HRefitCand> > &cands, TString &task, std::vector<Int_t> &pids, TLorentzVector lv = TLorentzVector(), HRefitCand mother = HRefitCand(), Double_t mass=0.);
    ~HDecayBuilder(){};

    void setVerbosity(int val) { fVerbose = val; }
    
    //setters
    void setIniSys(TLorentzVector val) {fIniSys = val;}
    void setMother(HRefitCand val) {fMother = val;}
    void setMass(Double_t val) {fMass = val;}

    void setPidsInVertices(std::vector<int> val1, std::vector<int> val2){
        fPidsPrim=val1;
        fPidsDecay=val2;
    }

    void buildDecay();

    void createNeutralCandidate();

    bool do4cFit();
    bool do3cFit();
    bool doMissMomFit();

    void fillFitCands();
    void checkDoubleParticle(size_t i);

    // Functions for getting the pulls

    void createOutputParticle(HRefitCand FittedCand);
    void getFitCands(std::vector<HRefitCand> &cands) { cands = fOutputCands; }
};

#endif /* HDECAYBUILDER_H */
