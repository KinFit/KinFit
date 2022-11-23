/**
 * HVertexFinder.h
 *
 *
 */

#ifndef HVERTEXFINDER_H
#define HVERTEXFINDER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// ROOT includes
#include "TMatrixD.h"

// framework includes
#include "KFitParticle.h"

using std::cout;
using std::endl;

class HVertexFinder
{

private:
    std::vector<KFitParticle> fCands;
    // std::vector<int> fPids; // Index of a Pid in this vector correspond to that of the cand on the corresponding index in fCands.
    // std::vector<double> fWeight;

    TVector3 fVertex;
    int fVerbose;

    TMatrixD fM; // Temporal matrix for calculations

    TVector3 fDir;

    TVector3 fBase;

protected:
    TMatrixD fSys; // LSM system inverse matrix
    TVector3 fB;         // LSM independent term

public:
    HVertexFinder(std::vector<KFitParticle> &);
    ~HVertexFinder(){};

    void setVerbosity(int val) { fVerbose = val; }

    void addLinesToVertex(const TVector3 &r, const TVector3 &alpha, const Double_t w = 1.0);

    void reset();

    void findVertex();

    TVector3 getVertex() const { return fVertex; }
};

#endif /* HVERTEXFINDER_H */
