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

// Tests for KFitDecayCandFinder: building a Lambda candidate from its
// daughters and the primary and decay vertices.

#include "KinFitTest.h"
#include "ToyTracks.h"

#include "KFitDecayCandFinder.h"
#include "KinFitter.h"

#include <TMatrixDSymEigen.h>
#include <TVector2.h>

namespace {

KFitDecayCandFinder makeFinder(const toy::LambdaEvent &ev) {
  std::vector<KinFitParticle> daughters = {
      toy::exactTrack(ev.proton, ev.decayVertex, toy::Resolution()),
      toy::exactTrack(ev.pion, ev.decayVertex, toy::Resolution())};
  return KFitDecayCandFinder(daughters, toy::kMassLambda, ev.primaryVertex,
                             ev.decayVertex, 0.5, 0.5, 1.0, 2.0, 2.0, 4.0);
}

} // namespace

// |p| = sqrt((sum E)^2 - m^2) reproduces the true Lambda momentum.
KF_TEST(momentum_from_daughter_energies) {
  TRandom3 rnd(1);
  for (int i = 0; i < 100; ++i) {
    const auto ev = toy::lambdaEvent(rnd);
    const KinFitParticle cand = makeFinder(ev).getDecayCand();
    KF_CHECK_REL(cand.P(), ev.lambda.P(), 1e-9, 1);
    KF_CHECK_NEAR(cand.M(), toy::kMassLambda, 1e-6);
  }
}

// Direction = primary -> decay vertex.
KF_TEST(direction_from_vertices) {
  TRandom3 rnd(2);
  for (int i = 0; i < 100; ++i) {
    const auto ev = toy::lambdaEvent(rnd);
    const KinFitParticle cand = makeFinder(ev).getDecayCand();
    KF_CHECK_NEAR(cand.Theta(), ev.lambda.Theta(), 1e-9);
    KF_CHECK_NEAR(cand.Phi(), ev.lambda.Phi(), 1e-9);
  }
}

// R, Z must be those of the flight line (through the decay vertex along the
// flight direction). The line used to be anchored at (decay - primary) and
// the sign of R was taken from the primary vertex y.
KF_TEST(R_Z_of_flight_line) {
  TRandom3 rnd(3);
  for (int i = 0; i < 100; ++i) {
    const auto ev = toy::lambdaEvent(rnd);
    double R, Z;
    toy::straightRZ(ev.decayVertex, ev.lambda.Vect(), R, Z);
    const KinFitParticle cand = makeFinder(ev).getDecayCand();
    KF_CHECK_NEAR(cand.getR(), R, 1e-6);
    KF_CHECK_NEAR(cand.getZ(), Z, 1e-6);
  }
}

KF_TEST(covariance_accessors_agree) {
  TRandom3 rnd(4);
  const auto ev = toy::lambdaEvent(rnd);
  auto finder = makeFinder(ev);
  const TMatrixD a = finder.getCovarianceMatrixDecayCand();
  const TMatrixD b = finder.getDecayCand().getCovariance();
  KF_CHECK(a.GetNrows() == 5 && b.GetNrows() == 5);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      KF_CHECK_NEAR(a(i, j), b(i, j), 0);
}

// A covariance matrix must be symmetric and positive semi-definite, also
// with vertex correlations.
KF_TEST(covariance_symmetric_positive) {
  TRandom3 rnd(5);
  for (int i = 0; i < 50; ++i) {
    const auto ev = toy::lambdaEvent(rnd);
    std::vector<KinFitParticle> daughters = {
        toy::exactTrack(ev.proton, ev.decayVertex, toy::Resolution()),
        toy::exactTrack(ev.pion, ev.decayVertex, toy::Resolution())};
    KFitDecayCandFinder finder(daughters, toy::kMassLambda, ev.primaryVertex,
                               ev.decayVertex, 0.5, 0.5, 1.0, 2.0, 2.0, 4.0);
    finder.setVertexCorrelations(0.3, 0.1, -0.2, 0.3, 0.1, -0.2);
    const TMatrixD cov = finder.getCovarianceMatrixDecayCand();
    TMatrixDSym sym(5);
    for (int r = 0; r < 5; ++r)
      for (int c = 0; c < 5; ++c) {
        KF_CHECK_NEAR(cov(r, c), cov(c, r), 1e-12 * std::abs(cov(r, r)) + 1e-15);
        sym(r, c) = 0.5 * (cov(r, c) + cov(c, r));
      }
    TMatrixDSymEigen eig(sym);
    const TVectorD ev_ = eig.GetEigenValues();
    for (int k = 0; k < 5; ++k)
      KF_CHECK(ev_(k) >= -1e-12 * std::abs(ev_(0)));
  }
}

// The (theta, phi, R, Z) covariance agrees with the spread obtained by
// smearing both vertices (with correlations) and recomputing the candidate.
KF_TEST(covariance_matches_vertex_smearing) {
  TRandom3 rnd(7);
  const auto ev = toy::lambdaEvent(rnd);
  std::vector<KinFitParticle> daughters = {
      toy::exactTrack(ev.proton, ev.decayVertex, toy::Resolution()),
      toy::exactTrack(ev.pion, ev.decayVertex, toy::Resolution())};
  // Small resolutions keep the propagation linear.
  const double sP[3] = {0.05, 0.05, 0.1}, sD[3] = {0.2, 0.2, 0.4};
  const double corr[3] = {0.3, 0.1, -0.2}; // xy, xz, yz
  KFitDecayCandFinder finder(daughters, toy::kMassLambda, ev.primaryVertex,
                             ev.decayVertex, sP[0], sP[1], sP[2], sD[0], sD[1],
                             sD[2]);
  finder.setVertexCorrelations(corr[0], corr[1], corr[2], corr[0], corr[1],
                               corr[2]);
  const TMatrixD cov = finder.getCovarianceMatrixDecayCand();

  // Correlated Gaussian via Cholesky of the 3x3 correlation matrix.
  const double c0 = corr[0], c1 = corr[1], c2 = corr[2];
  const double l11 = std::sqrt(1 - c0 * c0);
  const double l21 = (c2 - c1 * c0) / l11;
  const double l22 = std::sqrt(1 - c1 * c1 - l21 * l21);
  auto smear = [&](const TVector3 &x, const double *s) {
    const double g0 = rnd.Gaus(), g1 = rnd.Gaus(), g2 = rnd.Gaus();
    return TVector3(x.X() + s[0] * g0, x.Y() + s[1] * (c0 * g0 + l11 * g1),
                    x.Z() + s[2] * (c1 * g0 + l21 * g1 + l22 * g2));
  };

  const int n = 20000;
  double sum[4] = {0}, sum2[4][4] = {{0}};
  const KinFitParticle ref = finder.getDecayCand();
  const double par0[4] = {ref.Theta(), ref.Phi(), ref.getR(), ref.getZ()};
  for (int i = 0; i < n; ++i) {
    KFitDecayCandFinder f(daughters, toy::kMassLambda,
                          smear(ev.primaryVertex, sP),
                          smear(ev.decayVertex, sD), sP[0], sP[1], sP[2],
                          sD[0], sD[1], sD[2]);
    const KinFitParticle c = f.getDecayCand();
    const double par[4] = {c.Theta() - par0[0],
                           TVector2::Phi_mpi_pi(c.Phi() - par0[1]),
                           c.getR() - par0[2], c.getZ() - par0[3]};
    for (int a = 0; a < 4; ++a) {
      sum[a] += par[a];
      for (int b = 0; b < 4; ++b)
        sum2[a][b] += par[a] * par[b];
    }
  }
  const char *name[4] = {"theta", "phi", "R", "Z"};
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b <= a; ++b) {
      const double mc = sum2[a][b] / n - sum[a] * sum[b] / (double(n) * n);
      const double an = cov(a + 1, b + 1);
      const double scale = std::sqrt(cov(a + 1, a + 1) * cov(b + 1, b + 1));
      kftest::info("cov(%s,%s): analytic %.4g, MC %.4g", name[a], name[b], an,
                   mc);
      KF_CHECK_NEAR(an, mc, 0.05 * scale);
    }
}

// The candidate can be used as mother in a 3C fit of exact daughters.
// Depends on correct R, Z (see R_Z_of_flight_line).
KF_TEST(usable_as_3C_mother) {
  TRandom3 rnd(6);
  int nConverged = 0;
  toy::Stats chi2;
  for (int i = 0; i < 20; ++i) {
    const auto ev = toy::lambdaEvent(rnd);
    std::vector<KinFitParticle> daughters = {
        toy::exactTrack(ev.proton, ev.decayVertex, toy::Resolution()),
        toy::exactTrack(ev.pion, ev.decayVertex, toy::Resolution())};
    KFitDecayCandFinder finder(daughters, toy::kMassLambda, ev.primaryVertex,
                               ev.decayVertex, 0.5, 0.5, 1.0, 2.0, 2.0, 4.0);
    KinFitter fitter(daughters);
    fitter.add3Constraint(finder.getDecayCand());
    if (fitter.fit()) {
      ++nConverged;
      chi2.add(fitter.getChi2());
    }
  }
  kftest::info("converged %d/20, <chi2> = %g", nConverged, chi2.mean());
  KF_CHECK(nConverged == 20);
  KF_CHECK(chi2.mean() < 1e-3);
}

KF_TEST_MAIN()
