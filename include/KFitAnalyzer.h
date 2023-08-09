/**
 * KFitAnalyzer.h
 *
 * @updated 01.08.2023
 * @version v1.0.0
 * 
 * "User interface" class to automatically analyze and fit root file 
 * Particle candidates must be stored in a TClones array of KFitParticles
 *
 */

#ifndef KFITANALYZER_H
#define KFITANALYZER_H

#include "KFitParticle.h"
#include "KFitDecayBuilder.h"

// framework includes
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

// system includes
#include <iostream>
#include <vector>
#include <cmath>

using std::cout;
using std::endl;

class KFitAnalyzer
{
private:
    std::vector< std::vector<KFitParticle> > fCandsFit; // Vector of vector of particle candidates for each PID
    std::vector<int> fPids; // Vector of allPIDs of particles included in fit
    int fEvents;    // Number of events to be analyzed

    int fVerbose;   // Verbosity

    //Read in data
    TTree *fTree;   // Input data tree
    TClonesArray *fCands_in = new TClonesArray("KFitParticle"); // Input TClones array of KFitParticles

    //Output data
    TFile *foutFile;    // Output file
    TTree *fTree_out;   // Output data tree

    /** @brief Function that returns the output tree
    */
    TTree* getFittedTree(){ return fTree_out; }

    /** @brief select and sort candidates according to their PID
    */
    void selectCandidates();

    /** @brief write output tree to output file, close
    */
    void finish(); 
    
public:
    /** @brief Constructor
     * @param inFileName Full file name of input root-file
     * @param outFileName Hull file name of ouput root-file
     * @param nEvents Number of events to be analyzed
    */
    KFitAnalyzer(TString inFileName, TString outFileName, int nEvents=-1);

    /** @brief Default deconstructor */
    ~KFitAnalyzer(){};
    
    void setVerbosity(int val){ fVerbose = val; }

    //---------------User functions---------------------------------------

    /** @brief Call this function to perform the automated fitting procedure
     *          reads input tree,
     *          initialized DecayBuilder,
     *          performs event loop,
     *          runs decay builder,
     *          fills output tree with fit result
     * @param task String defining the task to be performed. Possibilities: 4C, Vertex, Mass
     * @param pids Vector of PIDs of all particles included in fit
     * @param mass optional, mass or missing mass of particle if needed for fit
     * @param lv optinal, if 4-vector input is needed for fit
    */
    void doFitterTask(TString task, std::vector<int> pids, double mass=-1., TLorentzVector lv = TLorentzVector(), KFitParticle mother = KFitParticle());
    //void addFitterTask(TString task, std::vector<int> primPids, std::vector<int> decayPids); // Jenny, for 3C fit
    //void addBuilderTask(TString task, std::vector<int> pids, TLorentzVector lv);

    /** @brief Set PIDs for fit
     * @param val Vector of PIDs of particles included in fit
    */
    void setPids(std::vector<int> val){ fPids = val; }

    /** @brief Function that returns the PID vector
    */
    std::vector<int> getPids(){ return fPids; }

};

#endif /* KFITANALYZER_H */