//****************************************************************************
//*                    This file is part of KinFit.                          *
//*                                                                          *
//*             KinFit is distributed under the terms of the                 *
//*              GNU General Public License (GPL) version 3,                 *
//*                 copied verbatim in the file "LICENSE".                   *
//*                                                                          *
//*                             Copyright 2024                               *
//*             GSI Helmholtzzentrum für Schwerionenforschung                *
//*          This software is distributed under the terms of the             *
//*          GNU General Public Licence version 3 (GPL Version 3)            *
//*                                                                          *
//*     The copyright holders are listed in the file "COPYRIGHTHOLDERS".     *
//*               The authors are listed in the file "AUTHORS".              *
//****************************************************************************

/**
 * KFitVertexFinder.h
 *
 * @updated 03.08.2023
 * @version v1.0.0
 *
 * Vertex finder using a linear least-squares solution of track-line
 * intersections.
 */

#pragma once

// framework includes
#include "KinFitParticle.h"

// ROOT includes

// system includes
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

class KFitVertexFinder {

private:
  int fVerbose; // Verbosity level

  TVector3 fVertex; // Vertex after finding

  TMatrixD fM; // Temporal matrix for calculations

  TVector3 fDir; // Direction vector used for each track in the fitting

  TVector3 fBase; // Base vector used for each track in the fitting

  void addLinesToVertex(const TVector3 &r, const TVector3 &alpha,
                        const double w = 1.0);

  /** @brief Find the vertex via matrix operations.
   */
  void findVertex();

  /** @brief Reset internal matrices and vectors used in the fit.
   */
  void reset();

protected:
  TMatrixD fSys; // LSM system inverse matrix
  TVector3 fB;   // LSM independent term

public:
  /** @brief Constructor.
   * @param cands Input track candidates.
   */
  KFitVertexFinder(const std::vector<KinFitParticle> &cands);

  /** @brief Default destructor. */
  ~KFitVertexFinder() {};

  void setVerbosity(int val) { fVerbose = val; }

  /** @brief Return the fitted vertex.
   * @return A TVector3 with vertex x, y, and z positions.
   */
  TVector3 getVertex() const { return fVertex; }
};

