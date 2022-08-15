/**
 * HVertexFinder.h
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
#include "hvertexfinder.h"
#include "hparticletool.h"

using std::cout;
using std::endl;

class HNeutralCandFinder
{
private:
    std::vector<HRefitCand> fCands;

    TVector3 fVertex;
    TVector3 fPrimaryVertex;
    Int_t fVerbose = 0;

    Double_t fMomentumAfterDecay;
    Double_t fNeutralCandMass;

    HRefitCand fNeutralMotherCandidate;

    Double_t fDistParticle1Vertex;
    Double_t fDistParticle2Vertex;
    Double_t fDistParticle1Origin;
    Double_t fDistParticle2Origin;

    TMatrixD fCovarianceNeutralMother;
    Bool_t fPrimaryVertexFound;

public:
    HNeutralCandFinder(const std::vector<HRefitCand> &cands);
    ~HNeutralCandFinder(){};

    void setVerbosity(Int_t val) { fVerbose = val; }

    // The first function is for creating a neutral mother candidate if only information of the decay vertex is available
    // The second function is for creating the neutral candidate if information about the primary vertex is also available
    void setNeutralMotherCand(Double_t valMomentum, Double_t valTheta, Double_t valPhi, Double_t valR, Double_t ValZ, TVector3 decayVertex);
    void setNeutralMotherCand(TVector3 primVtx, TVector3 decayVtx);
    void setMassNutralCand(Double_t val) { fNeutralCandMass = val; }

    HRefitCand getNeutralMotherCandidate() { return fNeutralMotherCandidate; }

    TMatrixD getCovarianceMatrixNeutralMother() { return fCovarianceNeutralMother; }
};

#endif /* HNEUTRALCANDFINDER_H */
