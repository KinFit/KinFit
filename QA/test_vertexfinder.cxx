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

// Tests for KFitVertexFinder (straight lines) and KFitVertexFinderConstB
// (helices in constant Bz).

#include "KinFitTest.h"
#include "ToyTracks.h"

#include "KFitVertexFinder.h"
#include "KFitVertexFinderConstB.h"

namespace {

bool isFailureMarker(const TVector3 &v) {
  return v.X() < -999 && v.Y() < -999 && v.Z() < -999;
}

std::vector<KinFitParticle> straightTracksFrom(const TVector3 &v, int n,
                                               TRandom3 &rnd) {
  std::vector<KinFitParticle> cands;
  for (int t = 0; t < n; ++t) {
    TVector3 p;
    p.SetMagThetaPhi(rnd.Uniform(200, 1000), rnd.Uniform(0.3, 1.3),
                     rnd.Uniform(-TMath::Pi(), TMath::Pi()));
    TLorentzVector lv;
    lv.SetVectM(p, toy::kMassPion);
    cands.push_back(toy::exactTrack(lv, v, toy::Resolution()));
  }
  return cands;
}

std::vector<KinFitParticle> helixTracksFrom(const TVector3 &v, int n,
                                            double bz, double ptMin,
                                            double ptMax, TRandom3 &rnd) {
  std::vector<KinFitParticle> cands;
  for (int t = 0; t < n; ++t) {
    TVector3 p;
    p.SetPtThetaPhi(rnd.Uniform(ptMin, ptMax), rnd.Uniform(0.2, 2.9),
                    rnd.Uniform(-TMath::Pi(), TMath::Pi()));
    cands.push_back(
        toy::helixTrackAtPCA(v, p, (t % 2) ? -1 : 1, bz, toy::kMassPion));
  }
  return cands;
}

} // namespace

// ---------------------------------------------------------------------------
// KFitVertexFinder
// ---------------------------------------------------------------------------

KF_TEST(lines_exact_intersection) {
  TRandom3 rnd(10);
  double worst = 0;
  for (int ev = 0; ev < 200; ++ev) {
    const TVector3 v(rnd.Uniform(-30, 30), rnd.Uniform(-30, 30),
                     rnd.Uniform(-60, 60));
    auto cands = straightTracksFrom(v, 2 + ev % 4, rnd);
    KFitVertexFinder finder(cands);
    worst = std::max(worst, (finder.getVertex() - v).Mag());
  }
  kftest::info("worst |vertex - truth| = %g mm", worst);
  KF_CHECK(worst < 1e-4);
}

KF_TEST(lines_parallel_tracks_fail_cleanly) {
  std::vector<KinFitParticle> cands;
  for (int k = 0; k < 2; ++k) {
    TLorentzVector lv;
    lv.SetXYZM(0, 0, 500, toy::kMassPion);
    cands.emplace_back(lv, 1.0 + k, 0.0);
  }
  KFitVertexFinder finder(cands);
  KF_CHECK(isFailureMarker(finder.getVertex()));
}

// Nearly parallel tracks (not exactly along z): the determinant is tiny but
// not exactly 0, which used to give a vertex far outside the detector.
KF_TEST(lines_nearly_parallel_fail_cleanly) {
  std::vector<KinFitParticle> cands;
  for (int k = 0; k < 2; ++k) {
    TLorentzVector lv;
    lv.SetXYZM(200, 100 + 1e-7 * k, 500, toy::kMassPion);
    cands.emplace_back(lv, 1.0 + k, 3.0);
  }
  KFitVertexFinder finder(cands);
  KF_CHECK(isFailureMarker(finder.getVertex()));
}

// ---------------------------------------------------------------------------
// KFitVertexFinderConstB
// ---------------------------------------------------------------------------

KF_TEST(helix_exact_vertex) {
  TRandom3 rnd(42);
  double worstHelix = 0, bestLine = 1e9;
  for (int ev = 0; ev < 200; ++ev) {
    const TVector3 v(rnd.Uniform(-30, 30), rnd.Uniform(-30, 30),
                     rnd.Uniform(-50, 50));
    auto cands = helixTracksFrom(v, 3, 1.0, 150, 800, rnd);
    KFitVertexFinderConstB helix(cands, 1.0);
    KFitVertexFinder line(cands);
    worstHelix = std::max(worstHelix, (helix.getVertex() - v).Mag());
    bestLine = std::min(bestLine, (line.getVertex() - v).Mag());
  }
  kftest::info("worst helix |d| = %g mm, best straight-line |d| = %g mm",
               worstHelix, bestLine);
  KF_CHECK(worstHelix < 1e-4);
}

KF_TEST(helix_low_pt_loopers) {
  TRandom3 rnd(7);
  double worst = 0;
  for (int ev = 0; ev < 300; ++ev) {
    const TVector3 v(rnd.Uniform(-5, 5), rnd.Uniform(-5, 5),
                     rnd.Uniform(-20, 20));
    auto cands = helixTracksFrom(v, 3, 1.0, 20, 120, rnd);
    KFitVertexFinderConstB finder(cands, 1.0);
    worst = std::max(worst, (finder.getVertex() - v).Mag());
  }
  kftest::info("worst |d| = %g mm", worst);
  KF_CHECK(worst < 1e-4);
}

KF_TEST(helix_negative_field) {
  TRandom3 rnd(43);
  double worst = 0;
  for (int ev = 0; ev < 100; ++ev) {
    const TVector3 v(rnd.Uniform(-30, 30), rnd.Uniform(-30, 30),
                     rnd.Uniform(-50, 50));
    auto cands = helixTracksFrom(v, 3, -0.7, 150, 800, rnd);
    KFitVertexFinderConstB finder(cands, -0.7);
    worst = std::max(worst, (finder.getVertex() - v).Mag());
  }
  KF_CHECK(worst < 1e-4);
}

KF_TEST(helix_zero_field_equals_lines) {
  TRandom3 rnd(11);
  for (int ev = 0; ev < 50; ++ev) {
    const TVector3 v(rnd.Uniform(-30, 30), rnd.Uniform(-30, 30),
                     rnd.Uniform(-50, 50));
    auto cands = straightTracksFrom(v, 3, rnd);
    for (auto &c : cands)
      c.setCharge(1);
    KFitVertexFinderConstB helix(cands, 0.0);
    KFitVertexFinder line(cands);
    KF_CHECK_NEAR((helix.getVertex() - line.getVertex()).Mag(), 0, 1e-6);
  }
}

KF_TEST(helix_neutral_tracks_are_lines) {
  TRandom3 rnd(12);
  const TVector3 v(10, -5, 20);
  auto cands = straightTracksFrom(v, 3, rnd); // charge 0 by default
  KFitVertexFinderConstB finder(cands, 1.0);
  KF_CHECK_NEAR((finder.getVertex() - v).Mag(), 0, 1e-4);
}

KF_TEST(helix_refit_after_setters) {
  TRandom3 rnd(13);
  const TVector3 v(12, -8, 30);
  auto cands = helixTracksFrom(v, 3, 1.0, 150, 800, rnd);
  KFitVertexFinderConstB finder(cands, -1.0); // wrong field sign
  const double before = (finder.getVertex() - v).Mag();
  finder.setMagneticFieldBz(1.0);
  finder.setMaxIterations(20);
  finder.setConvergenceTolerance(1e-5);
  finder.fit();
  const double after = (finder.getVertex() - v).Mag();
  kftest::info("wrong field: %g mm, after refit: %g mm", before, after);
  KF_CHECK(before > 1e-2);
  KF_CHECK(after < 1e-4);
}

KF_TEST(helix_parallel_tracks_fail_cleanly) {
  std::vector<KinFitParticle> cands;
  for (int k = 0; k < 2; ++k) {
    TLorentzVector lv;
    lv.SetXYZM(0, 0, 500, toy::kMassPion);
    cands.emplace_back(lv, 1.0 + k, 0.0);
    cands.back().setCharge(1);
  }
  KFitVertexFinderConstB finder(cands, 0.0);
  KF_CHECK(isFailureMarker(finder.getVertex()));
}

KF_TEST_MAIN()
