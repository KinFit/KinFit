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

// Tests for CoordinateConversion: (p, theta, phi, R, Z) <-> Cartesian and
// error propagation.

#include "KinFitTest.h"
#include "ToyTracks.h"

#include "CoordinateConversion.h"

namespace {

struct Track {
  double p, theta, phi, R, Z;
};

Track randomTrack(TRandom3 &rnd, double thetaMin, double thetaMax) {
  return {rnd.Uniform(100, 1000), rnd.Uniform(thetaMin, thetaMax),
          rnd.Uniform(-TMath::Pi(), TMath::Pi()), rnd.Uniform(-10, 10),
          rnd.Uniform(-50, 50)};
}

void checkCartToSph(double thetaMin, double thetaMax, unsigned seed) {
  TRandom3 rnd(seed);
  for (int i = 0; i < 200; ++i) {
    const Track t = randomTrack(rnd, thetaMin, thetaMax);
    // Cartesian input: momentum and the PCA as point on the track.
    TVector3 mom;
    mom.SetMagThetaPhi(t.p, t.theta, t.phi);
    double cart[6] = {-t.R * std::sin(t.phi), t.R * std::cos(t.phi), t.Z,
                      mom.X(), mom.Y(), mom.Z()};
    CoordinateConversion cc;
    cc.setparametersCart(cart);
    cc.setvarsSpher();
    KF_CHECK_REL(cc.getmom(), t.p, 1e-9, 1);
    KF_CHECK_NEAR(cc.gettheta(), t.theta, 1e-9);
    KF_CHECK_NEAR(cc.getphi(), t.phi, 1e-9);
    KF_CHECK_NEAR(cc.getR(), t.R, 1e-6);
    KF_CHECK_NEAR(cc.getZ(), t.Z, 1e-6);
  }
}

} // namespace

KF_TEST(sph_to_cart) {
  TRandom3 rnd(1);
  for (int i = 0; i < 200; ++i) {
    const Track t = randomTrack(rnd, 0.1, 3.0);
    double sph[5] = {t.p, t.theta, t.phi, t.R, t.Z};
    CoordinateConversion cc;
    cc.setparametersSpher(sph);
    cc.setvarsCart();
    KF_CHECK_NEAR(cc.getPx(), t.p * std::sin(t.theta) * std::cos(t.phi), 1e-6);
    KF_CHECK_NEAR(cc.getPy(), t.p * std::sin(t.theta) * std::sin(t.phi), 1e-6);
    KF_CHECK_NEAR(cc.getPz(), t.p * std::cos(t.theta), 1e-6);
    KF_CHECK_NEAR(cc.getX(), -t.R * std::sin(t.phi), 1e-6);
    KF_CHECK_NEAR(cc.getY(), t.R * std::cos(t.phi), 1e-6);
  }
}

KF_TEST(cart_to_sph_forward) { checkCartToSph(0.1, 1.5, 2); }

// theta used to be atan(pT/pz), negative for backward tracks.
KF_TEST(cart_to_sph_backward) { checkCartToSph(1.65, 3.0, 3); }

// Error propagation spherical -> Cartesian for *uncorrelated* input errors
// against a numerical Jacobian. The formulas used to add cross terms
// 2 (df/da)(df/db) sigma_a sigma_b as if the inputs were fully correlated.
KF_TEST(errors_sph_to_cart) {
  TRandom3 rnd(4);
  int nBad = 0;
  for (int i = 0; i < 100; ++i) {
    const Track t = randomTrack(rnd, 0.2, 2.9);
    const double err[5] = {0.02 * t.p, 0.002, 0.003, 0.5, 1.0};
    double sph[5] = {t.p, t.theta, t.phi, t.R, t.Z};
    CoordinateConversion cc;
    cc.setparametersSpher(sph);
    double e[5] = {err[0], err[1], err[2], err[3], err[4]};
    cc.setparametererrorsSpher(e);
    cc.setvarsCart();
    cc.seterrorsCart();

    // Numerical propagation: sigma_f^2 = sum_k (df/dx_k)^2 sigma_k^2.
    auto cart = [](const double *x, int which) {
      switch (which) {
      case 0: return x[0] * std::sin(x[1]) * std::cos(x[2]);
      case 1: return x[0] * std::sin(x[1]) * std::sin(x[2]);
      case 2: return x[0] * std::cos(x[1]);
      case 3: return -x[3] * std::sin(x[2]);
      default: return x[3] * std::cos(x[2]);
      }
    };
    const double got[5] = {cc.errPx, cc.errPy, cc.errPz, cc.errX, cc.errY};
    for (int w = 0; w < 5; ++w) {
      double var = 0;
      for (int k = 0; k < 5; ++k) {
        double xp[5] = {t.p, t.theta, t.phi, t.R, t.Z};
        double xm[5] = {t.p, t.theta, t.phi, t.R, t.Z};
        const double h = 1e-6 * std::max(std::abs(xp[k]), 1e-3);
        xp[k] += h;
        xm[k] -= h;
        const double d = (cart(xp, w) - cart(xm, w)) / (2 * h);
        var += d * d * err[k] * err[k];
      }
      const double expected = std::sqrt(var);
      if (!(std::abs(got[w] - expected) <= 1e-4 * expected + 1e-9)) {
        if (++nBad <= 5)
          kftest::info("component %d: got %g, expected %g", w, got[w],
                       expected);
      }
    }
  }
  KF_CHECK(nBad == 0);
}

KF_TEST_MAIN()
