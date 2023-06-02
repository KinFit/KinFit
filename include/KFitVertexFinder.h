/**
 * KFitVertexFinder.h
 *
 *@updated 27.04.2023
 *@version v1.0.0 
 *
 * Class 
 */

#ifndef KFITVERTEXFINDER_H
#define KFITVERTEXFINDER_H

// framework includes
#include "KFitParticle.h"

// ROOT includes
#include "TMatrixD.h"

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using std::cout;
using std::endl;

class KFitVertexFinder
{

private:
    std::vector<KFitParticle> fCands;

    TVector3 fVertex;
    int fVerbose;

    TMatrixD fM; // Temporal matrix for calculations

    TVector3 fDir; // Direction vector used for each track in the fitting

    TVector3 fBase; // Base vector used for each track in the fitting
 
    void addLinesToVertex(const TVector3 &r, const TVector3 &alpha, const double w = 1.0);

    /** @brief Function that finds the vertex via matrix multiplications
    */
    void findVertex();    
    
    void reset();

protected:
    TMatrixD fSys; // LSM system inverse matrix
    TVector3 fB;   // LSM independent term

public:
    
    /**  Constructor 
    */
    KFitVertexFinder(std::vector<KFitParticle> &);

    /** Default Destructor **/
    ~KFitVertexFinder(){};

    void setVerbosity(int val) { fVerbose = val; }

    /** @brief Function that returns the vertex
    * @return A Tvector3 with the vertex X, Y and Z positions
    */
    TVector3 getVertex() const { return fVertex; }
};

#endif /* KFITVERTEXFINDER_H */
