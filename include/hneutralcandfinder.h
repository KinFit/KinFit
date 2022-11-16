/**
 * HNeutralCandFinder.h
 *
 *
 */

#ifndef HNEUTRALCANDFINDER_H
#define HNEUTRALCANDFINDER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
// framework includes
#include "hrefitcand.h"
//#include "hvertexfinder.h"
//#include "hparticletool.h"

using std::cout;
using std::endl;

class HNeutralCandFinder
{
private:
    std::vector<HRefitCand> fCands;

    TVector3 fDecayVertex;
    TVector3 fPrimaryVertex;
    Int_t fVerbose = 0;

    Double_t fMomentumAfterDecay;
    Double_t fNeutralCandMass;

    HRefitCand fNeutralMotherCandidate;

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
    HNeutralCandFinder(const std::vector<HRefitCand> &cands, double neutralCandMass, TVector3 decayVertex, TVector3 primaryVertex, double primVtxResX, double primVtxResY, double primVtxResZ, double decVtxResX, double decVtxResY, double decVtxResZ);
    HNeutralCandFinder(const std::vector<HRefitCand> &cands, TVector3 decayVertex, TVector3 primaryVertex);
    ~HNeutralCandFinder(){};

    void setVerbosity(Int_t val) { fVerbose = val; }

    // The first function is for creating a neutral mother candidate if only information of the decay vertex is available
    // The second function is for creating the neutral candidate if information about the primary vertex is also available
    //void setNeutralMotherCand(Double_t valMomentum, Double_t valTheta, Double_t valPhi, Double_t valR, Double_t ValZ, TVector3 decayVertex);
    void setNeutralMotherCand();
    //void setMassNutralCand(Double_t val) { fNeutralCandMass = val; }

    HRefitCand getNeutralMotherCandidate() { return fNeutralMotherCandidate; }

    TMatrixD getCovarianceMatrixNeutralMother() { return fCovarianceNeutralMother; }
};

#endif /* HNEUTRALCANDFINDER_H */
