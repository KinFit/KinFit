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

#include "KFitDecayCandFinder.h"

#include <cmath>

/**
 * @file KFitDecayCandFinder.cxx
 * @brief Implementation of decay-candidate kinematics and covariance estimate.
 */

KFitDecayCandFinder::KFitDecayCandFinder(
    const std::vector<KinFitParticle> &cands, TVector3 primaryVertex,
    TVector3 decayVertex)
    : fCands(cands), fVerbose(0), fMomentumBeforeDecay(-1.),
      fPrimaryVertex(primaryVertex), fDecayVertex(decayVertex),
      fDecayCandMass(1115.683), fPrimVtxResX(1.78590), fPrimVtxResY(1.75516),
      fPrimVtxResZ(3.00431), fDecVtxResX(5.75369), fDecVtxResY(5.57198),
      fDecVtxResZ(10.2602) {
  if (fVerbose > 0) {
    std::cout << "--------------- KFitDecayCandFinder -----------------"
              << std::endl;
  }

  calculateDecayCand();
}

KFitDecayCandFinder::KFitDecayCandFinder(
    const std::vector<KinFitParticle> &cands, double DecayCandMass,
    TVector3 primaryVertex, TVector3 decayVertex, double primVtxResX,
    double primVtxResY, double primVtxResZ, double decVtxResX,
    double decVtxResY, double decVtxResZ)
    : fCands(cands), fVerbose(0), fMomentumBeforeDecay(-1.),
      fPrimaryVertex(primaryVertex), fDecayVertex(decayVertex),
      fDecayCandMass(DecayCandMass), fPrimVtxResX(primVtxResX),
      fPrimVtxResY(primVtxResY), fPrimVtxResZ(primVtxResZ),
      fDecVtxResX(decVtxResX), fDecVtxResY(decVtxResY),
      fDecVtxResZ(decVtxResZ) {
  if (fVerbose > 0) {
    std::cout << "--------------- KFitDecayCandFinder -----------------"
              << std::endl;
  }

  calculateDecayCand();
}

void KFitDecayCandFinder::calculateDecayCand() {
  if (fVerbose > 0) {
    std::cout
        << " ----------- KFitDecayCandFinder::calculateDecayCand() -----------"
        << std::endl;
  }

  // --- Momentum: |p| = sqrt((sum E_i)^2 - m^2) with the mass hypothesis. ---
  // (A local variable used to shadow fMomentumBeforeDecay, which stayed -1.)
  double totalEnergy = 0;
  TVector3 sumP;
  for (const KinFitParticle &cand : fCands) {
    totalEnergy += cand.E();
    sumP += cand.Vect();
  }
  const double p2 = totalEnergy * totalEnergy - fDecayCandMass * fDecayCandMass;
  // Below threshold (sum E < m) fall back to |sum p| instead of NaN.
  fMomentumBeforeDecay = p2 > 0 ? std::sqrt(p2) : sumP.Mag();

  // --- Direction: primary -> decay vertex ---
  const TVector3 &D = fDecayVertex;
  const TVector3 v = fDecayVertex - fPrimaryVertex;
  const double theta = v.Theta();
  const double phi = v.Phi();
  const double st = std::sin(theta), ct = std::cos(theta);
  const double sp = std::sin(phi), cp = std::cos(phi);

  // --- R, Z of the flight line through the decay vertex ---
  // Same convention as KinFitParticle: the point of closest approach to the
  // z axis is (-R sin(phi), R cos(phi), Z). With L = D.(cos phi, sin phi):
  //   R = -sin(phi) Dx + cos(phi) Dy,   Z = Dz - cot(theta) L.
  // (The line used to be anchored at D - P, and the sign of R was taken from
  // the primary vertex y.)
  const double cotTheta = ct / st;
  const double L = cp * D.X() + sp * D.Y();
  const double valR = -sp * D.X() + cp * D.Y();
  const double valZ = D.Z() - cotTheta * L;

  fDecayCand.SetXYZM(fMomentumBeforeDecay * st * cp,
                     fMomentumBeforeDecay * st * sp, fMomentumBeforeDecay * ct,
                     fDecayCandMass);
  fDecayCand.setMomentum(fMomentumBeforeDecay);
  fDecayCand.setThetaRad(theta);
  fDecayCand.setPhiRad(phi);
  fDecayCand.setR(valR);
  fDecayCand.setZ(valZ);

  // --- Covariance of (theta, phi, R, Z) ---
  // Linear propagation from the primary (P) and decay (D) vertex covariances:
  //   J_P = -J_v,  J_D = J_v + J_direct,  C = J_P C_P J_P^T + J_D C_D J_D^T
  // where J_v is the derivative w.r.t. v = D - P through (theta, phi) and
  // J_direct the explicit dependence of R, Z on D. The old code mixed sqrt
  // of correlations, asymmetric signs and wrong derivatives (dtheta/dx used
  // z/r^2, dphi had |x|).
  const double rho2 = v.Perp2();
  const double rho = std::sqrt(rho2);
  const double r2 = v.Mag2();
  TMatrixD dAngles(2, 3); // d(theta, phi) / d(vx, vy, vz)
  dAngles(0, 0) = v.X() * v.Z() / (r2 * rho);
  dAngles(0, 1) = v.Y() * v.Z() / (r2 * rho);
  dAngles(0, 2) = -rho / r2;
  dAngles(1, 0) = -v.Y() / rho2;
  dAngles(1, 1) = v.X() / rho2;
  dAngles(1, 2) = 0;

  // d(R, Z)/d(theta, phi) at fixed D
  const double dR_dphi = -L;
  const double dZ_dtheta = L / (st * st);
  const double dZ_dphi = -cotTheta * valR;

  TMatrixD Jv(4, 3); // d(theta, phi, R, Z) / dv
  for (int k = 0; k < 3; ++k) {
    Jv(0, k) = dAngles(0, k);
    Jv(1, k) = dAngles(1, k);
    Jv(2, k) = dR_dphi * dAngles(1, k);
    Jv(3, k) = dZ_dtheta * dAngles(0, k) + dZ_dphi * dAngles(1, k);
  }
  TMatrixD JD(Jv); // + explicit dependence on D
  JD(2, 0) += -sp;
  JD(2, 1) += cp;
  JD(3, 0) += -cotTheta * cp;
  JD(3, 1) += -cotTheta * sp;
  JD(3, 2) += 1.;
  const TMatrixD JP = -1. * Jv;

  auto vertexCov = [](double sx, double sy, double sz, double cxy, double cxz,
                      double cyz) {
    TMatrixD c(3, 3);
    c(0, 0) = sx * sx;
    c(1, 1) = sy * sy;
    c(2, 2) = sz * sz;
    c(0, 1) = c(1, 0) = cxy * sx * sy;
    c(0, 2) = c(2, 0) = cxz * sx * sz;
    c(1, 2) = c(2, 1) = cyz * sy * sz;
    return c;
  };
  const TMatrixD CP = vertexCov(fPrimVtxResX, fPrimVtxResY, fPrimVtxResZ,
                                fCorrPrimXY, fCorrPrimXZ, fCorrPrimYZ);
  const TMatrixD CD = vertexCov(fDecVtxResX, fDecVtxResY, fDecVtxResZ,
                                fCorrDecXY, fCorrDecXZ, fCorrDecYZ);
  const TMatrixD C = JP * CP * TMatrixD(TMatrixD::kTransposed, JP) +
                     JD * CD * TMatrixD(TMatrixD::kTransposed, JD);

  fCovarianceDecayCand.ResizeTo(5, 5);
  fCovarianceDecayCand.Zero();
  fCovarianceDecayCand(0, 0) = 9999999; // 1/p: free parameter in the 3C fit
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      fCovarianceDecayCand(i + 1, j + 1) = 0.5 * (C(i, j) + C(j, i));

  fDecayCand.setCovariance(fCovarianceDecayCand);

  if (fVerbose > 0) {
    std::cout << "calculateDecayCandidate, fDecayCand: theta= "
              << fDecayCand.Theta() << " phi= " << fDecayCand.Phi()
              << " R= " << valR << " Z= " << valZ << std::endl;
  }
}
