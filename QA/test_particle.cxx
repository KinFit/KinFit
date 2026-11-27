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

// Tests for KinFitParticle: construction, (R, Z) convention, covariance and
// ROOT I/O round trip.

#include "KinFitTest.h"
#include "ToyTracks.h"

#include "KinFitParticle.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

#include <memory>

KF_TEST(construct_from_R_Z) {
  TLorentzVector lv;
  lv.SetXYZM(100, -200, 500, toy::kMassPion);
  KinFitParticle p(lv, 1.5, -3.0);

  KF_CHECK_NEAR(p.P(), lv.P(), 1e-9);
  KF_CHECK_NEAR(p.M(), toy::kMassPion, 1e-6);
  KF_CHECK_NEAR(p.getR(), 1.5, 0);
  KF_CHECK_NEAR(p.getZ(), -3.0, 0);
  // Separately stored copies must agree with the four-vector.
  KF_CHECK_NEAR(p.getMomentum(), lv.P(), 1e-9);
  KF_CHECK_NEAR(p.getThetaRad(), lv.Theta(), 1e-12);
  KF_CHECK_NEAR(p.getPhiRad(), lv.Phi(), 1e-12);
  KF_CHECK(p.getPid() == -1);
  KF_CHECK(p.getTrackId() == -1);
}

// The (lv, X, Y, Z) constructor must give R and Z of the point of closest
// approach of the straight track through (X, Y, Z) to the z axis, with the
// documented sign convention PCA = (-R sin phi, R cos phi, Z).
KF_TEST(construct_from_point_gives_PCA) {
  TRandom3 rnd(1);
  for (int i = 0; i < 1000; ++i) {
    const TVector3 v(rnd.Uniform(-50, 50), rnd.Uniform(-50, 50),
                     rnd.Uniform(-100, 100));
    TVector3 mom;
    mom.SetMagThetaPhi(500, rnd.Uniform(0.2, 2.9),
                       rnd.Uniform(-TMath::Pi(), TMath::Pi()));
    TLorentzVector lv;
    lv.SetVectM(mom, toy::kMassPion);

    double R, Z;
    toy::straightRZ(v, mom, R, Z);
    KinFitParticle p(lv, v.X(), v.Y(), v.Z());
    KF_CHECK_NEAR(p.getR(), R, 1e-9);
    KF_CHECK_NEAR(p.getZ(), Z, 1e-9);
  }
}

KF_TEST(sign_of_R) {
  // Track along +x passing the beam axis at y = +2: PCA = (0, 2, 0).
  // phi = 0 -> PCA = (-R sin 0, R cos 0) = (0, R) -> R = +2.
  TLorentzVector lv;
  lv.SetXYZM(300, 0, 0, toy::kMassPion);
  KinFitParticle above(lv, 10.0, 2.0, 0.0);
  KF_CHECK_NEAR(above.getR(), 2.0, 1e-9);
  KinFitParticle below(lv, 10.0, -2.0, 0.0);
  KF_CHECK_NEAR(below.getR(), -2.0, 1e-9);
}

KF_TEST(covariance_5x5) {
  KinFitParticle p;
  TMatrixD cov(5, 5);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      cov(i, j) = (i == j) ? 1.0 + i : 0.1 * (i + j);
  p.setCovariance(cov);
  const TMatrixD got = p.getCovariance();
  KF_CHECK(got.GetNrows() == 5 && got.GetNcols() == 5);
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      KF_CHECK_NEAR(got(i, j), cov(i, j), 0);
}

KF_TEST(covariance_wrong_size_ignored) {
  KinFitParticle p;
  TMatrixD cov(5, 5);
  cov(0, 0) = 7;
  p.setCovariance(cov);
  TMatrixD wrong(3, 3);
  wrong(0, 0) = 99;
  p.setCovariance(wrong);
  KF_CHECK_NEAR(p.getCovariance()(0, 0), 7, 0);
}

// Write a TClonesArray of KinFitParticles to a tree and read it back.
// The (x, y, z) constructor used to leave the covariance 0x0.
KF_TEST(construct_from_point_has_5x5_covariance) {
  TLorentzVector lv;
  lv.SetXYZM(100, 200, 300, 938.272);
  const KinFitParticle p(lv, 1., 2., 3.);
  KF_CHECK(p.getCovariance().GetNrows() == 5);
  KF_CHECK(p.getCovariance().GetNcols() == 5);
}

// setPhiRad/Deg used to overwrite theta of the internal vector.
KF_TEST(setters_store_values) {
  KinFitParticle p;
  p.setThetaRad(0.3);
  p.setPhiRad(1.2);
  KF_CHECK_NEAR(p.getThetaRad(), 0.3, 0);
  KF_CHECK_NEAR(p.getPhiRad(), 1.2, 0);
  p.setPhiDeg(45.);
  KF_CHECK_NEAR(p.getPhiDeg(), 45., 1e-12);
  p.setTOF(12.5); // used to be truncated to int
  p.setRapidity(0.75);
  KF_CHECK_NEAR(p.getTOF(), 12.5, 0);
  KF_CHECK_NEAR(p.getRapidity(), 0.75, 0);
}

KF_TEST(io_roundtrip_tclonesarray) {
  const char *fname = "test_particle_io.root";
  TRandom3 rnd(2);
  std::vector<KinFitParticle> written;
  {
    TFile f(fname, "RECREATE");
    TTree t("data", "io test");
    auto *arr = new TClonesArray("KinFitParticle");
    t.Branch("KinFitParticle", &arr);
    for (int ev = 0; ev < 3; ++ev) {
      arr->Clear();
      for (int k = 0; k < 2; ++k) {
        TLorentzVector lv;
        lv.SetXYZM(rnd.Uniform(-300, 300), rnd.Uniform(-300, 300),
                   rnd.Uniform(100, 900), toy::kMassProton);
        auto *p = new ((*arr)[k]) KinFitParticle(lv, rnd.Gaus(0, 1),
                                                 rnd.Gaus(0, 5));
        p->setCovariance(toy::covariance(lv.P(), toy::Resolution()));
        p->setPid(14 + k);
        p->setTrackId(100 * ev + k);
        p->setCharge(k ? -1 : 1);
        written.push_back(*p);
      }
      t.Fill();
    }
    t.Write();
    delete arr;
  }

  TFile f(fname, "READ");
  auto *t = f.Get<TTree>("data");
  KF_CHECK(t != nullptr);
  if (!t)
    return;
  KF_CHECK(t->GetEntries() == 3);
  auto *arr = new TClonesArray("KinFitParticle");
  t->SetBranchAddress("KinFitParticle", &arr);
  size_t idx = 0;
  for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
    t->GetEntry(ev);
    KF_CHECK(arr->GetEntries() == 2);
    for (int k = 0; k < arr->GetEntries(); ++k, ++idx) {
      const auto *p = static_cast<KinFitParticle *>(arr->At(k));
      const KinFitParticle &w = written[idx];
      KF_CHECK_NEAR(p->Px(), w.Px(), 1e-9);
      KF_CHECK_NEAR(p->Py(), w.Py(), 1e-9);
      KF_CHECK_NEAR(p->Pz(), w.Pz(), 1e-9);
      KF_CHECK_NEAR(p->M(), w.M(), 1e-6);
      KF_CHECK_NEAR(p->getR(), w.getR(), 0);
      KF_CHECK_NEAR(p->getZ(), w.getZ(), 0);
      KF_CHECK(p->getPid() == w.getPid());
      KF_CHECK(p->getTrackId() == w.getTrackId());
      KF_CHECK_NEAR(p->getCharge(), w.getCharge(), 0);
      for (int i = 0; i < 5; ++i)
        KF_CHECK_NEAR(p->getCovariance()(i, i), w.getCovariance()(i, i), 0);
    }
  }
  t->ResetBranchAddresses();
  delete arr;
}

KF_TEST_MAIN()
