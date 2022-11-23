/**
 * HVertexTools.h
 *
 *
 */

#ifndef HVERTEXTOOLS_H
#define HVERTEXTOOLS_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// framework includes
#include "KFitParticle.h"
#include "hparticletool.h"
#include "hgeomvertexfit.h"

using std::cout;
using std::endl;

class HVertexTools

{
private:
    std::vector<KFitParticle> fCands;

    TVector3 fVertex;
    TVector3 fPrimaryVertex;
    Int_t fVerbose;

    Double_t fDistanceParticleToParticle;
    Double_t fDistanceParticleToVertex;

    Double_t fDistParticle1Vertex;
    Double_t fDistParticle2Vertex;
    Double_t fDistParticle1Origin;
    Double_t fDistParticle2Origin;

    Bool_t fPrimaryVertexFound;

    // Properties related to difference between primary and decay vertex
    TVector3 fVecPrimToDecayVertex;
    Double_t fDistPrimToDecayVertex;

    // Checks if the decay vertex is located further downstream in z-position
    // and at larger values of R compared to the primary vertex
    Bool_t fPrimaryVertexIsBetforeDecayVertex;
    Bool_t fPrimaryVertexIsInsideDecayVertex;

    Bool_t fUsePrimaryVertexInNeutralCandidateCalculation;

public:
    HVertexTools();
    ~HVertexTools(){};

    void setVerbosity(Int_t val) { fVerbose = val; }

    TVector3 findVertex(const std::vector<KFitParticle> &cands);
    TVector3 findPrimaryVertex(const std::vector<KFitParticle> &cands);
    std::vector<KFitParticle> UpdateTrackParameters(std::vector<KFitParticle> &cands, TVector3 &VertexPos);
    void calculateVertexProperties(TVector3 primaryVertex, TVector3 decayVertex);

    TVector3 getVertex() const { return fVertex; }               // Function that the user should use in the analysis macro
    TVector3 gePrimarytVertex() const { return fPrimaryVertex; } // Function that the user should use in the analysis macro

    Double_t getDistanceBetweenFittedParticles() const { return fDistanceParticleToParticle; }
    Double_t getDistanceFirstParticleVertex() const { return fDistParticle1Vertex; }
    Double_t getDistanceSecondParticleVertex() const { return fDistParticle2Vertex; }
    Double_t getDistanceFirstParticleOrigin() const { return fDistParticle1Origin; }
    Double_t getDistanceSecondParticleOrigin() const { return fDistParticle2Origin; }

    Double_t getDistBetweenVertices() { return fDistPrimToDecayVertex; }

    Bool_t isPrimVertexBeforeDecayVertex() { return fPrimaryVertexIsBetforeDecayVertex; }
    Bool_t isPrimVertexInsideDecayVertex() { return fPrimaryVertexIsInsideDecayVertex; }

    TVector3 getPrimaryVertex() { return fPrimaryVertex; }
};

#endif /* HVERTEXTOOLS_H */
