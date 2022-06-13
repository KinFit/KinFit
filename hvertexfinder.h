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

// framework includes
#include "hrefitcand.h"
#include "hgeomvector.h"
#include "hparticletool.h"
#include "hgeomvertexfit.h"

using std::cout;
using std::endl;

class HVertexFinder
{

private:
    std::vector<HRefitCand> fCands;
    //std::vector<int> fPids; // Index of a Pid in this vector correspond to that of the cand on the corresponding index in fCands.
    //std::vector<double> fWeight;

    TVector3 fVertex;
    int fVerbose;

public:
    HVertexFinder();
    ~HVertexFinder(){};

    void setVerbosity(int val) { fVerbose = val; }

    TVector3 findVertex(const std::vector<HRefitCand> &cands);

    TVector3 getVertex() const { return fVertex; }

};

#endif /* HVERTEXFINDER_H */
