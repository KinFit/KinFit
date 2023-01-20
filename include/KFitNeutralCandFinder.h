/**
 * KFitNeutralCandFinder.h
 *
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
//#include "KFitVertexFinder.h"
//#include "hparticletool.h"

using std::cout;
using std::endl;

class KFitNeutralCandFinder
{
private:
    std::vector<KFitParticle> fCands;

    TVector3 fDecayVertex;
    TVector3 fPrimaryVertex;
    Int_t fVerbose = 0;

    Double_t fMomentumAfterDecay;
    Double_t fNeutralCandMass;

    KFitParticle fNeutralMotherCandidate;

    Double_t fDistParticle1Vertex;
    Double_t fDistParticle2Vertex;
    Double_t fDistParticle1Origin;
    Double_t fDistParticle2Origin;

    double fPrimVtxResX;
    double fPrimVtxResY;
    double fPrimVtxResZ;

    double fDecVtxResX;
    double fDecVtxResY;
    double fDecVtxResZ;

    TMatrixD fCovarianceNeutralMother;

public:
    KFitNeutralCandFinder(const std::vector<KFitParticle> &cands, double neutralCandMass, TVector3 primaryVertex, TVector3 decayVertex, double primVtxResX, double primVtxResY, double primVtxResZ, double decVtxResX, double decVtxResY, double decVtxResZ);
    KFitNeutralCandFinder(const std::vector<KFitParticle> &cands, TVector3 primaryVertex, TVector3 decayVertex);
    ~KFitNeutralCandFinder(){};

    void setVerbosity(Int_t val) { fVerbose = val; }

    // The first function is for creating a neutral mother candidate if only information of the decay vertex is available
    // The second function is for creating the neutral candidate if information about the primary vertex is also available
    //void setNeutralMotherCand(Double_t valMomentum, Double_t valTheta, Double_t valPhi, Double_t valR, Double_t ValZ, TVector3 decayVertex);
    void setNeutralMotherCand();
    //void setMassNutralCand(Double_t val) { fNeutralCandMass = val; }

    KFitParticle getNeutralMotherCandidate() { return fNeutralMotherCandidate; }

    TMatrixD getCovarianceMatrixNeutralMother() { return fCovarianceNeutralMother; }
};

#endif /* KFITNEUTRALCANDFINDER_H */
