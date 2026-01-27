//****************************************************************************
//*                    This file is part of KinFit.                          *
//*                                                                          *
//*            	KinFit is distributed under the terms of the                 *
//*              GNU General Public License (GPL) version 3,                 *
//*                 copied verbatim in the file "LICENSE".                   *
//*                                                                          *
//*  				Copyright 2024                               *
//*		GSI Helmholtzzentrum für Schwerionenforschung                *
//* 	     This software is distributed under the terms of the             *
//*	     GNU General Public Licence version 3 (GPL Version 3)            *
//*		      			     				     *
//*     The copyright holders are listed in the file "COPYRIGHTHOLDERS".     *
//*               The authors are listed in the file "AUTHORS".              *
//****************************************************************************

/**
 * KFitVertexFinder.h
 *
 *@updated 03.08.2023
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

    int fVerbose; // Verbosity level

    std::vector<KFitParticle> fCands;

    TVector3 fVertex; // Vertex after finding

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
    KFitVertexFinder(std::vector<KFitParticle> &cands);

    /** Default Destructor **/
    ~KFitVertexFinder(){};

    void setVerbosity(int val) { fVerbose = val; }

    /** @brief Function that returns the vertex
    * @return A Tvector3 with the vertex X, Y and Z positions
    */
    TVector3 getVertex() const { return fVertex; }
};

#endif /* KFITVERTEXFINDER_H */
