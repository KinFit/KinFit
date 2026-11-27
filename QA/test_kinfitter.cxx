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

// Tests for KinFitter, per constraint type:
//   exact_<C>    : unsmeared input -> converges, chi2 ~ 0, parameters kept
//   pulls_<C>    : smeared input   -> constraint fulfilled after the fit,
//                                     pulls ~ N(0,1), <prob> ~ 0.5
//   jacobian_<C> : analytic Feta_eval / Fxi_eval == finite differences of
//                  f_eval
// plus a few behavioural edge cases.
//
// Physics: Lambda -> p pi- from a displaced vertex (see ToyTracks.h).

#include "KinFitTest.h"
#include "ToyTracks.h"

#include "KinFitter.h"

#include <string>

namespace {

enum class C { FourC, Mass, Vertex, MassVtx, MissingMass, MissingParticle, ThreeC };

const char *name(C c) {
  switch (c) {
  case C::FourC: return "4C";
  case C::Mass: return "Mass";
  case C::Vertex: return "Vertex";
  case C::MassVtx: return "MassVtx";
  case C::MissingMass: return "MissingMass";
  case C::MissingParticle: return "MissingParticle";
  case C::ThreeC: return "3C";
  }
  return "?";
}

/** Everything a fit needs plus the truth. */
struct Scenario {
  std::vector<KinFitParticle> cands; // measured
  KinFitParticle mother;             // 3C only (measured direction/position)
  TLorentzVector init;               // initial system (4C, Missing*)
  double mass = 0;                   // constraint mass
  TLorentzVector missingTruth;       // MissingParticle only
};

const toy::Resolution kTrackRes;

// Resolution of the Lambda direction/position for the 3C fit.
constexpr double kMotherSigmaTheta = 0.01, kMotherSigmaPhi = 0.01;
constexpr double kMotherSigmaR = 1.0, kMotherSigmaZ = 2.0;

KinFitParticle motherCandidate(const toy::LambdaEvent &ev, bool smear,
                               TRandom3 &rnd) {
  double R, Z;
  toy::straightRZ(ev.primaryVertex, ev.lambda.Vect(), R, Z);
  double theta = ev.lambda.Theta(), phi = ev.lambda.Phi();
  if (smear) {
    theta = rnd.Gaus(theta, kMotherSigmaTheta);
    phi = rnd.Gaus(phi, kMotherSigmaPhi);
    R = rnd.Gaus(R, kMotherSigmaR);
    Z = rnd.Gaus(Z, kMotherSigmaZ);
  }
  // |p| is a free parameter in the 3C fit; start 5 % off.
  TLorentzVector lv;
  TVector3 p;
  p.SetMagThetaPhi(1.05 * ev.lambda.P(), theta, phi);
  lv.SetVectM(p, toy::kMassLambda);
  KinFitParticle m(lv, R, Z);
  TMatrixD cov(5, 5);
  cov(0, 0) = 1e6; // not used: 1/p of the mother is unmeasured
  cov(1, 1) = kMotherSigmaTheta * kMotherSigmaTheta;
  cov(2, 2) = kMotherSigmaPhi * kMotherSigmaPhi;
  cov(3, 3) = kMotherSigmaR * kMotherSigmaR;
  cov(4, 4) = kMotherSigmaZ * kMotherSigmaZ;
  m.setCovariance(cov);
  return m;
}

Scenario makeScenario(C c, bool smear, TRandom3 &rnd) {
  const toy::LambdaEvent ev = toy::lambdaEvent(rnd);
  auto track = [&](const TLorentzVector &lv) {
    return smear ? toy::smearedTrack(lv, ev.decayVertex, kTrackRes, rnd)
                 : toy::exactTrack(lv, ev.decayVertex, kTrackRes);
  };

  Scenario s;
  switch (c) {
  case C::FourC:
    s.cands = {track(ev.proton), track(ev.pion)};
    s.init = ev.lambda;
    break;
  case C::Mass:
  case C::Vertex:
  case C::MassVtx:
    s.cands = {track(ev.proton), track(ev.pion)};
    s.mass = toy::kMassLambda;
    break;
  case C::MissingMass:
    // Lambda known, only the proton measured, pion missing.
    s.cands = {track(ev.proton)};
    s.init = ev.lambda;
    s.mass = toy::kMassPion;
    break;
  case C::MissingParticle: {
    // Initial system = Lambda + K+; p and pi- measured, K+ missing.
    TVector3 pK;
    pK.SetMagThetaPhi(rnd.Uniform(300, 800), rnd.Uniform(0.3, 1.2),
                      rnd.Uniform(-TMath::Pi(), TMath::Pi()));
    s.missingTruth.SetVectM(pK, toy::kMassKaon);
    s.cands = {track(ev.proton), track(ev.pion)};
    s.init = ev.lambda + s.missingTruth;
    s.mass = toy::kMassKaon;
    break;
  }
  case C::ThreeC:
    s.cands = {track(ev.proton), track(ev.pion)};
    s.mother = motherCandidate(ev, smear, rnd);
    break;
  }
  return s;
}

void addConstraint(KinFitter &f, C c, const Scenario &s) {
  switch (c) {
  case C::FourC: f.add4Constraint(s.init); break;
  case C::Mass: f.addMassConstraint(s.mass); break;
  case C::Vertex: f.addVertexConstraint(); break;
  case C::MassVtx: f.addMassVtxConstraint(s.mass); break;
  case C::MissingMass: f.addMissingMassConstraint(s.init, s.mass); break;
  case C::MissingParticle: f.addMissingParticleConstraint(s.init, s.mass); break;
  case C::ThreeC: f.add3Constraint(s.mother); break;
  }
}

/** Indices (5*i + j) of the measured parameters that enter the constraint. */
std::vector<int> constrainedParameters(C c, int nCands) {
  std::vector<int> idx;
  for (int i = 0; i < nCands; ++i) {
    for (int j = 0; j < 5; ++j) {
      const bool momentum = j < 3;
      const bool vertexPar = j >= 1;
      bool used = false;
      switch (c) {
      case C::Vertex: used = vertexPar && i < 2; break;
      case C::MassVtx: used = momentum || (vertexPar && i < 2); break;
      default: used = momentum; break;
      }
      if (used)
        idx.push_back(5 * i + j);
    }
  }
  if (c == C::ThreeC) // mother theta, phi (R, Z do not enter the 3C equations)
    for (int j = 0; j < 2; ++j)
      idx.push_back(5 * nCands + j);
  return idx;
}

TVector3 basePoint(const KinFitParticle &p) {
  return TVector3(-p.getR() * std::sin(p.Phi()), p.getR() * std::cos(p.Phi()),
                  p.getZ());
}

double lineDistance(const KinFitParticle &a, const KinFitParticle &b) {
  const TVector3 n = a.Vect().Unit().Cross(b.Vect().Unit());
  return std::abs(n.Dot(basePoint(a) - basePoint(b))) / n.Mag();
}

TLorentzVector sum(const std::vector<KinFitParticle> &v) {
  TLorentzVector s;
  for (const auto &p : v)
    s += p;
  return s;
}

/** Check that the fitted particles fulfil the constraint. */
void checkConstraintFulfilled(C c, const Scenario &s, KinFitter &f) {
  std::vector<KinFitParticle> d;
  f.getDaughters(d);
  const double tolP = 1e-2;  // MeV
  const double tolX = 1e-3;  // mm
  switch (c) {
  case C::FourC: {
    const TLorentzVector diff = sum(d) - s.init;
    KF_CHECK_NEAR(diff.Px(), 0, tolP);
    KF_CHECK_NEAR(diff.Py(), 0, tolP);
    KF_CHECK_NEAR(diff.Pz(), 0, tolP);
    KF_CHECK_NEAR(diff.E(), 0, tolP);
    break;
  }
  case C::Mass:
    KF_CHECK_NEAR(sum(d).M(), s.mass, tolP);
    break;
  case C::Vertex:
    KF_CHECK_NEAR(lineDistance(d[0], d[1]), 0, tolX);
    break;
  case C::MassVtx:
    KF_CHECK_NEAR(sum(d).M(), s.mass, tolP);
    KF_CHECK_NEAR(lineDistance(d[0], d[1]), 0, tolX);
    break;
  case C::MissingMass:
    KF_CHECK_NEAR((s.init - sum(d)).M(), s.mass, tolP);
    break;
  case C::MissingParticle: {
    const TLorentzVector diff = sum(d) + f.getMissingDaughter() - s.init;
    KF_CHECK_NEAR(diff.Px(), 0, tolP);
    KF_CHECK_NEAR(diff.Py(), 0, tolP);
    KF_CHECK_NEAR(diff.Pz(), 0, tolP);
    KF_CHECK_NEAR(diff.E(), 0, tolP);
    KF_CHECK_NEAR(f.getMissingDaughter().M(), s.mass, tolP);
    break;
  }
  case C::ThreeC: {
    const TLorentzVector diff = sum(d) - f.getMother();
    KF_CHECK_NEAR(diff.Px(), 0, tolP);
    KF_CHECK_NEAR(diff.Py(), 0, tolP);
    KF_CHECK_NEAR(diff.Pz(), 0, tolP);
    KF_CHECK_NEAR(diff.E(), 0, tolP);
    KF_CHECK_NEAR(f.getMother().M(), toy::kMassLambda, tolP);
    break;
  }
  }
}

// ---------------------------------------------------------------------------
// Generic test bodies
// ---------------------------------------------------------------------------

void runExact(C c) {
  TRandom3 rnd(100);
  int nConverged = 0;
  const int nEvents = 50;
  double worstChi2 = 0, worstShift = 0;
  for (int ev = 0; ev < nEvents; ++ev) {
    Scenario s = makeScenario(c, false, rnd);
    KinFitter f(s.cands);
    addConstraint(f, c, s);
    if (!f.fit())
      continue;
    ++nConverged;
    worstChi2 = std::max(worstChi2, std::abs(f.getChi2()));
    for (size_t i = 0; i < s.cands.size(); ++i)
      worstShift = std::max(
          worstShift, (f.getDaughter(i).Vect() - s.cands[i].Vect()).Mag());
    if (c == C::MissingParticle)
      KF_CHECK_NEAR((f.getMissingDaughter().Vect() - s.missingTruth.Vect()).Mag(),
                    0, 1e-3);
  }
  kftest::info("%s exact: converged %d/%d, worst chi2 %g, worst |dp| %g MeV/c",
               name(c), nConverged, nEvents, worstChi2, worstShift);
  KF_CHECK(nConverged == nEvents);
  KF_CHECK(worstChi2 < 1e-6);
  KF_CHECK(worstShift < 1e-3);
}

void runPulls(C c) {
  TRandom3 rnd(200);
  const int nEvents = 1000;
  int nConverged = 0;
  toy::Stats prob;
  int nProbBelow5 = 0;
  std::vector<toy::Stats> pulls;
  std::vector<int> idx;

  for (int ev = 0; ev < nEvents; ++ev) {
    Scenario s = makeScenario(c, true, rnd);
    KinFitter f(s.cands);
    addConstraint(f, c, s);
    if (!f.fit())
      continue;
    ++nConverged;
    if (nConverged <= 20) // spot-check the constraint on a subset
      checkConstraintFulfilled(c, s, f);
    prob.add(f.getProb());
    if (f.getProb() < 0.05)
      ++nProbBelow5;
    if (idx.empty()) {
      idx = constrainedParameters(c, s.cands.size());
      pulls.resize(idx.size());
    }
    for (size_t k = 0; k < idx.size(); ++k)
      pulls[k].add(f.getPull(idx[k]));
  }

  kftest::info("%s: converged %d/%d, <prob> = %.3f +- %.3f, P(prob<0.05) = %.3f",
               name(c), nConverged, nEvents, prob.mean(), prob.meanError(),
               nConverged ? double(nProbBelow5) / nConverged : 0.);
  KF_CHECK(nConverged >= 0.99 * nEvents);
  if (nConverged < 100)
    return;

  KF_CHECK_NEAR(prob.mean(), 0.5, 0.04);
  KF_CHECK_NEAR(double(nProbBelow5) / nConverged, 0.05, 0.025);
  const char *par[5] = {"1/p", "theta", "phi", "R", "Z"};
  for (size_t k = 0; k < idx.size(); ++k) {
    const int i = idx[k] / 5, j = idx[k] % 5;
    kftest::info("  pull particle %d %-5s: mean %+.3f width %.3f", i, par[j],
                 pulls[k].mean(), pulls[k].rms());
    KF_CHECK_NEAR(pulls[k].mean(), 0, 0.12);
    KF_CHECK_NEAR(pulls[k].rms(), 1, 0.12);
  }
}

/** Measured-parameter vector y and unmeasured xi as KinFitter builds them. */
void fitterVectors(C c, const Scenario &s, TMatrixD &y, TMatrixD &xi) {
  const TMatrixD yDaughters = toy::parameterVector(s.cands);
  y.ResizeTo(yDaughters);
  y = yDaughters;
  xi.ResizeTo(1, 1);
  xi.Zero();
  if (c == C::ThreeC) {
    const int n = y.GetNrows();
    y.ResizeTo(n + 4, 1);
    y(n + 0, 0) = s.mother.Theta();
    y(n + 1, 0) = s.mother.Phi();
    y(n + 2, 0) = s.mother.getR();
    y(n + 3, 0) = s.mother.getZ();
    xi(0, 0) = 1.0 / s.mother.P();
  }
  if (c == C::MissingParticle) {
    xi.ResizeTo(3, 1);
    const TLorentzVector miss = s.init - sum(s.cands);
    xi(0, 0) = miss.Px();
    xi(1, 0) = miss.Py();
    xi(2, 0) = miss.Pz();
  }
}

/** Compare an analytic Jacobian with central finite differences of f. */
void compareJacobian(const char *what, const TMatrixD &analytic,
                     const std::function<TMatrixD(const TMatrixD &)> &f,
                     const TMatrixD &x0) {
  const TMatrixD f0 = f(x0);
  KF_CHECK(analytic.GetNrows() == f0.GetNrows());
  KF_CHECK(analytic.GetNcols() == x0.GetNrows());
  if (analytic.GetNrows() != f0.GetNrows() ||
      analytic.GetNcols() != x0.GetNrows())
    return;

  int nBad = 0;
  for (int k = 0; k < x0.GetNrows(); ++k) {
    const double h = 1e-6 * std::max(std::abs(x0(k, 0)), 1e-3);
    TMatrixD xp(x0), xm(x0);
    xp(k, 0) += h;
    xm(k, 0) -= h;
    const TMatrixD fp = f(xp), fm = f(xm);
    for (int r = 0; r < f0.GetNrows(); ++r) {
      // Scale: largest analytic entry of this row (derivatives of one
      // constraint share its units).
      double scale = 0;
      for (int cc = 0; cc < analytic.GetNcols(); ++cc)
        scale = std::max(scale, std::abs(analytic(r, cc)));
      const double numeric = (fp(r, 0) - fm(r, 0)) / (2 * h);
      const double tol = 1e-5 * std::max(scale, 1e-12) +
                         1e-6 * std::abs(numeric);
      if (std::abs(numeric - analytic(r, k)) > tol) {
        if (++nBad <= 5)
          kftest::info("%s: d f[%d] / d x[%d]: analytic %.6g, numeric %.6g",
                       what, r, k, analytic(r, k), numeric);
      }
    }
  }
  if (nBad)
    kftest::reportFailure(__FILE__, __LINE__,
                          std::string(what) + ": " + std::to_string(nBad) +
                              " Jacobian entries differ");
}

void runJacobian(C c) {
  TRandom3 rnd(300);
  for (int ev = 0; ev < 5; ++ev) {
    Scenario s = makeScenario(c, true, rnd);
    KinFitter fitter(s.cands);
    addConstraint(fitter, c, s);
    TMatrixD y, xi;
    fitterVectors(c, s, y, xi);

    compareJacobian("Feta", fitter.Feta_eval(y, xi),
                    [&](const TMatrixD &yy) { return fitter.f_eval(yy, xi); },
                    y);
    if (c == C::ThreeC || c == C::MissingParticle)
      compareJacobian(
          "Fxi", fitter.Fxi_eval(y, xi),
          [&](const TMatrixD &xx) { return fitter.f_eval(y, xx); }, xi);
  }
}

} // namespace

// ---------------------------------------------------------------------------
// Cases. Known failures are registered as DISABLED in QA/CMakeLists.txt.
// ---------------------------------------------------------------------------

KF_TEST(exact_4C) { runExact(C::FourC); }
KF_TEST(exact_Mass) { runExact(C::Mass); }
KF_TEST(exact_Vertex) { runExact(C::Vertex); }
KF_TEST(exact_MassVtx) { runExact(C::MassVtx); }
KF_TEST(exact_MissingMass) { runExact(C::MissingMass); }
KF_TEST(exact_MissingParticle) { runExact(C::MissingParticle); }
KF_TEST(exact_3C) { runExact(C::ThreeC); }

KF_TEST(pulls_4C) { runPulls(C::FourC); }
KF_TEST(pulls_Mass) { runPulls(C::Mass); }
KF_TEST(pulls_Vertex) { runPulls(C::Vertex); }
KF_TEST(pulls_MassVtx) { runPulls(C::MassVtx); }
KF_TEST(pulls_MissingMass) { runPulls(C::MissingMass); }
KF_TEST(pulls_MissingParticle) { runPulls(C::MissingParticle); }
KF_TEST(pulls_3C) { runPulls(C::ThreeC); }

KF_TEST(jacobian_4C) { runJacobian(C::FourC); }
KF_TEST(jacobian_Mass) { runJacobian(C::Mass); }
KF_TEST(jacobian_Vertex) { runJacobian(C::Vertex); }
KF_TEST(jacobian_MassVtx) { runJacobian(C::MassVtx); }
KF_TEST(jacobian_MissingMass) { runJacobian(C::MissingMass); }
KF_TEST(jacobian_MissingParticle) { runJacobian(C::MissingParticle); }
KF_TEST(jacobian_3C) { runJacobian(C::ThreeC); }

// Parameters that do not enter the constraint have an undefined pull,
// reported as -10000 (documented in docs/kinfitter.md).
KF_TEST(pull_of_unconstrained_parameter) {
  TRandom3 rnd(400);
  Scenario s = makeScenario(C::Mass, true, rnd);
  KinFitter f(s.cands);
  f.addMassConstraint(s.mass);
  KF_CHECK(f.fit());
  KF_CHECK_NEAR(f.getPull(3), -10000, 0); // R of particle 0
  KF_CHECK_NEAR(f.getPull(4), -10000, 0); // Z of particle 0
}

// Mass constraint on the first n particles only.
KF_TEST(mass_exclusive_candidates) {
  TRandom3 rnd(500);
  Scenario s = makeScenario(C::Mass, true, rnd);
  TLorentzVector extra;
  extra.SetXYZM(100, 50, 400, toy::kMassPion);
  s.cands.push_back(toy::smearedTrack(extra, TVector3(1, 1, 1), kTrackRes, rnd));
  KinFitter f(s.cands);
  f.setNumberOfExclusiveCandidates(2);
  f.addMassConstraint(s.mass);
  KF_CHECK(f.fit());
  KF_CHECK_NEAR((f.getDaughter(0) + f.getDaughter(1)).M(), s.mass, 1e-2);
  // The third particle is not constrained and must stay where it was.
  KF_CHECK_NEAR((f.getDaughter(2).Vect() - s.cands[2].Vect()).Mag(), 0, 1e-6);
}

// With a single measured track the 4x4 system D V D^T has rank <= 3:
// the fit must fail cleanly instead of returning garbage.
KF_TEST(missing_particle_single_track_fails_cleanly) {
  TRandom3 rnd(600);
  Scenario s = makeScenario(C::MissingParticle, true, rnd);
  s.cands.pop_back();
  KinFitter f(s.cands);
  f.addMissingParticleConstraint(s.init, s.mass);
  KF_CHECK(!f.fit());
}

// A measured parameter that is exactly zero (R = 0 for a track from the
// origin) must not prevent convergence (the relative convergence criterion
// used to divide by it).
KF_TEST(converges_with_zero_parameters) {
  TRandom3 rnd(700);
  int nConverged = 0;
  for (int ev = 0; ev < 20; ++ev) {
    Scenario s = makeScenario(C::Mass, true, rnd);
    for (auto &c : s.cands) {
      c.setR(0);
      c.setZ(0);
    }
    KinFitter f(s.cands);
    f.addMassConstraint(s.mass);
    nConverged += f.fit();
  }
  KF_CHECK(nConverged == 20);
}

// One iteration only: fit() reports success, the constraint is linearised
// once.
KF_TEST(single_iteration_returns_true) {
  TRandom3 rnd(800);
  Scenario s = makeScenario(C::Mass, true, rnd);
  KinFitter f(s.cands);
  f.addMassConstraint(s.mass);
  f.setNumberOfIterations(1);
  KF_CHECK(f.fit());
}

// Mass constraint with three particles: the Jacobian used to zero
// H(row, 4 + 4q), i.e. the phi derivative of particle 2.
KF_TEST(jacobian_Mass_three_particles) {
  TRandom3 rnd(900);
  for (int ev = 0; ev < 5; ++ev) {
    Scenario s = makeScenario(C::Mass, true, rnd);
    TLorentzVector extra;
    extra.SetXYZM(120, -60, 350, toy::kMassPion);
    s.cands.push_back(
        toy::smearedTrack(extra, TVector3(1, 1, 1), kTrackRes, rnd));
    KinFitter fitter(s.cands);
    fitter.addMassConstraint(1500.);
    TMatrixD y, xi;
    fitterVectors(C::Mass, s, y, xi);
    compareJacobian("Feta", fitter.Feta_eval(y, xi),
                    [&](const TMatrixD &yy) { return fitter.f_eval(yy, xi); },
                    y);
  }
}

// Two mass constraints on particle subsets via setUseCandInFit(): the path
// used to be unreachable (selection cleared before it was tested).
KF_TEST(two_mass_constraints_on_subsets) {
  TRandom3 rnd(901);
  int nConverged = 0;
  for (int ev = 0; ev < 20; ++ev) {
    const toy::LambdaEvent l1 = toy::lambdaEvent(rnd);
    const toy::LambdaEvent l2 = toy::lambdaEvent(rnd);
    std::vector<KinFitParticle> cands = {
        toy::smearedTrack(l1.proton, l1.decayVertex, kTrackRes, rnd),
        toy::smearedTrack(l2.proton, l2.decayVertex, kTrackRes, rnd),
        toy::smearedTrack(l1.pion, l1.decayVertex, kTrackRes, rnd),
        toy::smearedTrack(l2.pion, l2.decayVertex, kTrackRes, rnd)};
    KinFitter f(cands);
    f.setUseCandInFit(0);
    f.setUseCandInFit(2);
    f.addMassConstraint(toy::kMassLambda);
    f.setUseCandInFit(1);
    f.setUseCandInFit(3);
    f.addMassConstraint(toy::kMassLambda);
    if (!f.fit())
      continue;
    ++nConverged;
    KF_CHECK(f.getNDF() == 2);
    KF_CHECK_NEAR((f.getDaughter(0) + f.getDaughter(2)).M(), toy::kMassLambda,
                  1e-2);
    KF_CHECK_NEAR((f.getDaughter(1) + f.getDaughter(3)).M(), toy::kMassLambda,
                  1e-2);
  }
  KF_CHECK(nConverged == 20);
}

// A mass constraint is appended to another constraint (used to replace it).
KF_TEST(vertex_plus_mass) {
  TRandom3 rnd(902);
  int nConverged = 0;
  for (int ev = 0; ev < 20; ++ev) {
    Scenario s = makeScenario(C::Vertex, true, rnd);
    KinFitter f(s.cands);
    f.addVertexConstraint();
    f.addMassConstraint(toy::kMassLambda);
    if (!f.fit())
      continue;
    ++nConverged;
    KF_CHECK(f.getNDF() == 2);
    KF_CHECK_NEAR((f.getDaughter(0) + f.getDaughter(1)).M(), toy::kMassLambda,
                  1e-2);
    KF_CHECK_NEAR(lineDistance(f.getDaughter(0), f.getDaughter(1)), 0, 1e-3);
  }
  KF_CHECK(nConverged == 20);
}

// No constraint / two exclusive constraint types: fit() returns false
// (used to abort() resp. silently use only the last one).
KF_TEST(invalid_configuration_returns_false) {
  TRandom3 rnd(903);
  Scenario s = makeScenario(C::FourC, true, rnd);
  KinFitter none(s.cands);
  KF_CHECK(!none.fit());
  KinFitter both(s.cands);
  both.add4Constraint(s.init);
  both.addVertexConstraint();
  KF_CHECK(!both.fit());
}

// NDF = number of equations - number of unmeasured parameters.
KF_TEST(ndf_per_constraint) {
  TRandom3 rnd(904);
  const std::pair<C, int> expected[] = {
      {C::FourC, 4}, {C::Mass, 1},   {C::Vertex, 1},         {C::MassVtx, 2},
      {C::ThreeC, 3}, {C::MissingMass, 1}, {C::MissingParticle, 1}};
  for (const auto &e : expected) {
    Scenario s = makeScenario(e.first, true, rnd);
    KinFitter f(s.cands);
    addConstraint(f, e.first, s);
    f.fit();
    KF_CHECK(f.getNDF() == e.second);
  }
}

// The fitted covariance of each daughter is a symmetric 5x5 block with
// correlations, and fitting twice gives the same result (fit() used to
// overwrite the measurements with the fitted values).
KF_TEST(covariance_symmetric_and_refit_stable) {
  TRandom3 rnd(905);
  Scenario s = makeScenario(C::FourC, true, rnd);
  KinFitter f(s.cands);
  f.add4Constraint(s.init);
  KF_CHECK(f.fit());
  const double chi2 = f.getChi2();
  const TMatrixD cov = f.getDaughter(0).getCovariance();
  double offDiag = 0;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) {
      KF_CHECK_NEAR(cov(i, j), cov(j, i), 1e-9 * std::abs(cov(i, i)));
      if (i != j)
        offDiag = std::max(offDiag, std::abs(cov(i, j)));
    }
  KF_CHECK(offDiag > 0);
  KF_CHECK(f.fit());
  KF_CHECK_NEAR(f.getChi2(), chi2, 1e-9 * std::max(1., chi2));
}

KF_TEST_MAIN()
