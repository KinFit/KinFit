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
 * KFitDecayCandFinder.h
 *
 * @updated 03.08.2023
 * @version v1.0.0
 *
 * Class to calculate a decay candidate
 * from its decay products.
 *
 */

#pragma once

// framework includes
#include "KinFitParticle.h"

// system includes
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

class KFitDecayCandFinder {
private:
  std::vector<KinFitParticle> fCands; // Vector of decay products

  TVector3 fPrimaryVertex; // Vector pointing to the primary vertex
  TVector3 fDecayVertex;   // Vector pointing to the decay vertex
  int fVerbose = 0;

  double fMomentumBeforeDecay; // Estimated momentum of the decay particle
  double fDecayCandMass;       // Assumption for the mass of the decay particle

  KinFitParticle fDecayCand;

  double fPrimVtxResX; // Primary vertex resolution in x-direction
  double fPrimVtxResY; // Primary vertex resolution in y-direction
  double fPrimVtxResZ; // Primary vertex resolution in z-direction

  double fDecVtxResX; // Decay vertex resolution in x-direction
  double fDecVtxResY; // Decay vertex resolution in y-direction
  double fDecVtxResZ; // Decay vertex resolution in z-direction

  double fCorrPrimXY = 0; // Correlations between primary vertex uncertainties
                          // in X and Y direction
  double fCorrPrimXZ = 0; // Correlations between primary vertex uncertainties
                          // in X and Z direction
  double fCorrPrimYZ = 0; // Correlations between primary vertex uncertainties
                          // in Y and Z direction

  double fCorrDecXY =
      0; // Correlations between decay vertex uncertainties in X and Y direction
  double fCorrDecXZ =
      0; // Correlations between decay vertex uncertainties in X and Z direction
  double fCorrDecYZ =
      0; // Correlations between decay vertex uncertainties in Y and Z direction

  /** Covariance matrix of the decay candidate
   * Diagonal entries correspond to the covariances
   * in the parameters in the following order
   *
   * -----------------------------
   * | 1/p                       |
   * |     theta                 |
   * |            phi            |
   * |                   R       |
   * |                        Z  |
   * -----------------------------
   *
   * Off diagonal elements correspond to the
   * covariances between the parameters.
   */
  TMatrixD fCovarianceDecayCand;

public:
  /** @brief Constructor
   * @param cands Vector of decay particles.
   * @param decayCandMass Mass assumption for the decay particle.
   * @param primaryVertex Primary vertex position.
   * @param decayVertex Decay vertex position.
   * @param primVtxResX Primary vertex resolution in x-direction.
   * @param primVtxResY Primary vertex resolution in y-direction.
   * @param primVtxResZ Primary vertex resolution in z-direction.
   * @param decVtxResX Decay vertex resolution in x-direction.
   * @param decVtxResY Decay vertex resolution in y-direction.
   * @param decVtxResZ Decay vertex resolution in z-direction.
   */
  KFitDecayCandFinder(const std::vector<KinFitParticle> &cands,
                      double decayCandMass, TVector3 primaryVertex,
                      TVector3 decayVertex, double primVtxResX,
                      double primVtxResY, double primVtxResZ, double decVtxResX,
                      double decVtxResY, double decVtxResZ);

  /** @brief Constructor
   * Default vertex resolutions are used with a lambda-mass hypothesis.
   * @param cands Vector of decay particles.
   * @param primaryVertex Primary vertex position.
   * @param decayVertex Decay vertex position.
   */
  KFitDecayCandFinder(const std::vector<KinFitParticle> &cands,
                      TVector3 primaryVertex, TVector3 decayVertex);

  /** @brief Default destructor. */
  ~KFitDecayCandFinder() {};

  void setVerbosity(int val) { fVerbose = val; }

  /** @brief Set correlations between x, y, and z uncertainties of the primary
   * and decay vertices.
   * @param valPrimXY x-y correlation of primary-vertex uncertainties.
   * @param valPrimXZ x-z correlation of primary-vertex uncertainties.
   * @param valPrimYZ y-z correlation of primary-vertex uncertainties.
   * @param valDecXY x-y correlation of decay-vertex uncertainties.
   * @param valDecXZ x-z correlation of decay-vertex uncertainties.
   * @param valDecYZ y-z correlation of decay-vertex uncertainties.
   * Correlation coefficients in [-1, 1]; the candidate is recomputed.
   */
  void setVertexCorrelations(double valPrimXY, double valPrimXZ,
                             double valPrimYZ, double valDecXY, double valDecXZ,
                             double valDecYZ) {

    fCorrPrimXY = valPrimXY;
    fCorrPrimXZ = valPrimXZ;
    fCorrPrimYZ = valPrimYZ;

    fCorrDecXY = valDecXY;
    fCorrDecXZ = valDecXZ;
    fCorrDecYZ = valDecYZ;

    calculateDecayCand(); // used to have no effect after construction
  }

  /** @brief Return the calculated decay candidate.
   * @return Decay candidate as KinFitParticle.
   */
  KinFitParticle getDecayCand() { return fDecayCand; }

  /** @brief Return covariance matrix of the decay candidate.
   * @return Covariance matrix in (1/p, theta, phi, R, Z) ordering.
   */
  TMatrixD getCovarianceMatrixDecayCand() { return fCovarianceDecayCand; }

private:
  /** @brief Calculate decay-candidate kinematics and covariance.
   */
  void calculateDecayCand();
};

