/**
 * KFitNeutralCandFinder.h
 *
 *@updated 26.04.2023
 *@version 1.0 
 *
 * Class to calculate the neutral candidate
 * from two decay products 
 * 
 */

#ifndef KFITNEUTRALCANDFINDER_H
#define KFITNEUTRALCANDFINDER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
// framework includes
#include "KFitParticle.h"

using std::cout;
using std::endl;

class KFitNeutralCandFinder
{
private:
    std::vector<KFitParticle> fCands; // Vector of decay products

    TVector3 fPrimaryVertex; // Vector pointing to the primary vertex
    TVector3 fDecayVertex; // Vector pointing to the decay vertex
    Int_t fVerbose = 0;

    Double_t fMomentumBeforeDecay; // Estimated momentum of the neutral particle
    Double_t fNeutralCandMass; // Assuption for the mass of the neutral particle
    
    KFitParticle fNeutralMotherCandidate;

    //Double_t fDistParticle1Vertex;
    //Double_t fDistParticle2Vertex;
    //Double_t fDistParticle1Origin;
    //Double_t fDistParticle2Origin; 

    double fPrimVtxResX; // Primary vertex resolution in x-direction
    double fPrimVtxResY; // Primary vertex resolution in y-direction
    double fPrimVtxResZ; // Primary vertex resolution in z-direction

    double fDecVtxResX; // Decay vertex resolution in x-direction
    double fDecVtxResY; // Decay vertex resolution in y-direction
    double fDecVtxResZ; // Decay vertex resolution in z-direction

    /** Covariance matrix of the neutral mother candidate
    /* Diagonal entries correspond to the covariances 
    /* in the parameters in the following order
    /* 
    /* -----------------------------
    /* | 1/p                       |
    /* |     theta                 |
    /* |            phi            |
    /* |                   R       |
    /* |                        Z  |
    /* -----------------------------
    /* 
    /* Off diagonal elements corresponds to the
    /* correlations between the parameters
    */
    TMatrixD fCovarianceNeutralMother; 

public:
    /** @brief Constructor
    * @param cands - vector of decay particles
    * @param neutralCandMass - mass asumption for the neutral decaying particle
    * @param primaryVertex - vector with primary vertex positions
    * @param decayVertex - vector with decay vertex positions
    * @param primVtxResX - primary vertex resolution in x-direction
    * @param primVtxResY - primary vertex resolution in y-direction
    * @param primVtxResZ - primary vertex resolution in z-direction   
    * @param decVtxResX - decay vertex resolution in x-direction
    * @param decVtxResY - decay vertex resolution in y-direction
    * @param decVtxResZ - decay vertex resolution in z-direction
    */
    KFitNeutralCandFinder(const std::vector<KFitParticle> &cands, double neutralCandMass, TVector3 primaryVertex, TVector3 decayVertex, double primVtxResX, double primVtxResY, double primVtxResZ, double decVtxResX, double decVtxResY, double decVtxResZ);
    
    /** @brief Constructor
    * Decault vertex resolutions are used and lambda mass hypothesis
    * @param cands - vector of decay particles
    * @param primaryVertex - vector with primary vertex positions
    * @param decayVertex - vector with decay vertex positions
    */
    KFitNeutralCandFinder(const std::vector<KFitParticle> &cands, TVector3 primaryVertex, TVector3 decayVertex);
    
    /** Default Constructor **/
    ~KFitNeutralCandFinder(){};

    void setVerbosity(Int_t val) { fVerbose = val; }

    /** @brief Function that returns the neutral mother candidate
    */
    KFitParticle getNeutralMotherCandidate() { return fNeutralMotherCandidate; }

    /** @brief Function that returns the covariance matrix of the neutral mother candidate
    */
    TMatrixD getCovarianceMatrixNeutralMother() { return fCovarianceNeutralMother; }
    
private:

   /** @brief Function that calculates the neutral mother candidate
   */
    void calculateNeutralMotherCand();

};

#endif /* KFITNEUTRALCANDFINDER_H */
