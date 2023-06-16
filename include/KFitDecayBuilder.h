/**
 * KFitRootAnalyzer.h
 *
 * @updated 02.06.2023
 * @version v1.0.0
 * 
 * Class responsible for doing the combinatorics for each event in the automated
 * fitting procedure, calls KinFitter, chooses best combination
 *
 */

#ifndef KFITDECAYBUILDER_H
#define KFITDECAYBUILDER_H

// framework includes
#include "KFitParticle.h"
#include "KinFitter.h"

//#include "TH1F.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>

using std::cout;
using std::endl;

class KFitDecayBuilder
{
private:
    // Working Particles
    std::vector<std::vector<KFitParticle>> fCands;  // Vector of vector of particle candidates for each PID
    std::vector<KFitParticle> fFitCands;    // Vector of particles the fit should be performed on

    std::vector<KFitParticle> fOutputCands; // Output particles after fitting

    // Fitter input variables
    TString fTask;  // Task to be performed
    std::vector<int> fPids;  // Vector of PIDs of all particles included in fit
    TLorentzVector fIniSys; // Four-vector as input to fit
    KFitParticle fMother;   // Mother particle as input to fit
    double fMass;   // Mass as input to fit

    // For combinatorics
    int fTotalCombos;   // Total number of combinations for event
    int fCombiCounter;  // Number of combination that is evaluated 
    std::vector<int> particleCounter;   // Vector, each entry is a counter for respective PID of PID vector
    bool doubleParticle;    // True if same particle was used twice in same combination

    double fBestProb;   // Probability of best combination
    double fBestChi2;   // Chi2 of best combination

    int fVerbose;   // Verbosity 0-1-high
    
    /** @brief Fills vector of particles to be fitted with current combination */
    void fillFitCands();

    /** @brief Calls KFitter */
    bool doFit();

    /** @brief Checks if same particle was used twice in same combination */
    void checkDoubleParticle(size_t i);

public:
    /** @brief Constructor
     * @param task String defining the task to be performed. Possibilities: 4C, Vertex, Mass
     * @param pids Vector of PIDs of all particles included in fit
     * @param lv optinal, if 4-vector input is needed for fit
     * @param mother optional, if mother particle is needed for fit
     * @param mm optional, mass or missing mass of particle if needed for fit
    */
    KFitDecayBuilder(TString task, std::vector<int> pids, TLorentzVector lv = TLorentzVector(), KFitParticle mother = KFitParticle(), double mass = -1.);
    
    /** @brief Default deconstructor */
    ~KFitDecayBuilder(){};

    void setVerbosity(int val) { fVerbose = val; }

    // setters
    /** @brief Set input candidates
     * @param cands Vector of vector of particle candidates for each PID
    */
    void setInputCands(std::vector<std::vector<KFitParticle>> cands) {fCands = cands; }

    void setIniSys(TLorentzVector val) { fIniSys = val; }
    void setMother(KFitParticle val) { fMother = val; }
    void setMass(double val) { fMass = val; }

    /** @brief Count number of combination of particles in event */
    void countCombis();

    /** @brief Do combinatorics, call Fitter, choose best combination */
    void buildDecay();

    /** @brief Access fitted particles chosen */
    void getFitCands(std::vector<KFitParticle> &cands) { cands = fOutputCands; }

    /** @brief Returns chi2 of best combination */
    double getChi2() { return fBestChi2; }

    /** @brief Returns probability of best combination */
    double getProbability() { return fBestProb; }
};

#endif /* KFITDECAYBUILDER_H */
