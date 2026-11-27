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

// Tests for the analysis layer: KFitDecayBuilder combinatorics and the
// KFitAnalyzer file-to-file workflow.

#include "KinFitTest.h"
#include "ToyTracks.h"

#include "KFitAnalyzer.h"
#include "KFitDecayBuilder.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

namespace {

constexpr int kPidProton = 14, kPidPion = 9;
const toy::Resolution kRes;

/** Write `nEvents` Lambda -> p pi- events in the KFitAnalyzer input format. */
void writeLambdaFile(const char *fname, int nEvents, unsigned seed) {
  TRandom3 rnd(seed);
  TFile f(fname, "RECREATE");
  TTree t("data", "KinFit QA input");
  auto *arr = new TClonesArray("KinFitParticle");
  t.Branch("KinFitParticle", &arr);
  for (int ev = 0; ev < nEvents; ++ev) {
    arr->Clear();
    const auto le = toy::lambdaEvent(rnd);
    auto *p = new ((*arr)[0])
        KinFitParticle(toy::smearedTrack(le.proton, le.decayVertex, kRes, rnd));
    p->setPid(kPidProton);
    p->setTrackId(1000 * ev + 0);
    auto *pi = new ((*arr)[1])
        KinFitParticle(toy::smearedTrack(le.pion, le.decayVertex, kRes, rnd));
    pi->setPid(kPidPion);
    pi->setTrackId(1000 * ev + 1);
    t.Fill();
  }
  t.Write();
  delete arr;
}

struct OutputSummary {
  Long64_t entries = 0;
  int firstEvent = -1;
  int nBadMass = 0, nBadPid = 0, nBadProb = 0;
};

OutputSummary readOutput(const char *fname) {
  OutputSummary s;
  TFile f(fname, "READ");
  auto *t = f.Get<TTree>("data_fitted");
  if (!t)
    return s;
  int event;
  double prob, chi2;
  auto *arr = new TClonesArray("KinFitParticle");
  t->SetBranchAddress("Event", &event);
  t->SetBranchAddress("Prob", &prob);
  t->SetBranchAddress("Chi2", &chi2);
  t->SetBranchAddress("KinFitParticle", &arr);
  s.entries = t->GetEntries();
  for (Long64_t i = 0; i < s.entries; ++i) {
    t->GetEntry(i);
    if (i == 0)
      s.firstEvent = event;
    if (arr->GetEntries() != 2) {
      ++s.nBadPid;
      continue;
    }
    auto *p = static_cast<KinFitParticle *>(arr->At(0));
    auto *pi = static_cast<KinFitParticle *>(arr->At(1));
    if (p->getPid() != kPidProton || pi->getPid() != kPidPion)
      ++s.nBadPid;
    if (std::abs((*p + *pi).M() - toy::kMassLambda) > 1e-2)
      ++s.nBadMass;
    if (!(prob > 0 && prob <= 1))
      ++s.nBadProb;
  }
  t->ResetBranchAddresses();
  delete arr;
  return s;
}

} // namespace

// ---------------------------------------------------------------------------
// KFitDecayBuilder
// ---------------------------------------------------------------------------

// One proton, the true pion and a fake pion: the best combination
// (highest probability) should be the true one in most events.
KF_TEST(builder_selects_true_combination) {
  TRandom3 rnd(1);
  const int nEvents = 200;
  int nCorrect = 0, nResult = 0;
  for (int ev = 0; ev < nEvents; ++ev) {
    const auto le = toy::lambdaEvent(rnd);
    KinFitParticle p = toy::smearedTrack(le.proton, le.decayVertex, kRes, rnd);
    KinFitParticle piTrue = toy::smearedTrack(le.pion, le.decayVertex, kRes, rnd);
    TVector3 fakeMom;
    fakeMom.SetMagThetaPhi(rnd.Uniform(100, 400), rnd.Uniform(0.3, 1.2),
                           rnd.Uniform(-TMath::Pi(), TMath::Pi()));
    TLorentzVector fakeLv;
    fakeLv.SetVectM(fakeMom, toy::kMassPion);
    KinFitParticle piFake = toy::smearedTrack(fakeLv, TVector3(0, 0, 0), kRes, rnd);
    p.setTrackId(0);
    piTrue.setTrackId(1);
    piFake.setTrackId(2);

    KFitDecayBuilder builder("Mass", {kPidProton, kPidPion}, TLorentzVector(),
                             toy::kMassLambda);
    builder.setInputCands({{p}, {piFake, piTrue}});
    builder.countCombis();
    builder.buildDecay();
    std::vector<KinFitParticle> result;
    builder.getFitCands(result);
    if (result.size() == 2) {
      ++nResult;
      nCorrect += (result[1].getTrackId() == 1);
    }
  }
  kftest::info("result in %d/%d events, true pion chosen in %d", nResult,
               nEvents, nCorrect);
  KF_CHECK(nResult == nEvents);
  KF_CHECK(nCorrect >= 0.8 * nEvents);
}

// Two identical PIDs: a track must never be combined with itself.
KF_TEST(builder_no_double_use) {
  TRandom3 rnd(2);
  for (int ev = 0; ev < 50; ++ev) {
    const auto le = toy::lambdaEvent(rnd);
    // Use the "proton" PID for both tracks to exercise the identical-PID path.
    KinFitParticle a = toy::smearedTrack(le.proton, le.decayVertex, kRes, rnd);
    KinFitParticle b = toy::smearedTrack(le.pion, le.decayVertex, kRes, rnd);
    a.setTrackId(0);
    b.setTrackId(1);
    KFitDecayBuilder builder("Mass", {kPidProton, kPidProton}, TLorentzVector(),
                             toy::kMassLambda);
    builder.setInputCands({{a, b}, {a, b}});
    builder.countCombis();
    builder.buildDecay();
    std::vector<KinFitParticle> result;
    builder.getFitCands(result);
    KF_CHECK(result.size() == 2);
    if (result.size() == 2)
      KF_CHECK(result[0].getTrackId() != result[1].getTrackId());
  }
}

// ---------------------------------------------------------------------------
// KFitAnalyzer
// ---------------------------------------------------------------------------

KF_TEST(analyzer_mass_fit_output) {
  const int nEvents = 30;
  writeLambdaFile("test_analyzer_in.root", nEvents, 10);
  {
    KFitAnalyzer analyzer("test_analyzer_in.root", "test_analyzer_out.root", -1);
    analyzer.doFitterTask("Mass", {kPidProton, kPidPion}, toy::kMassLambda);
  }
  const OutputSummary s = readOutput("test_analyzer_out.root");
  kftest::info("output entries %lld (input %d), first Event %d", s.entries,
               nEvents, s.firstEvent);
  KF_CHECK(s.entries >= nEvents - 1);
  KF_CHECK(s.nBadPid == 0);
  KF_CHECK(s.nBadMass == 0);
  KF_CHECK(s.nBadProb == 0);
}

// Every input event must be analysed (the event loop used to start at 1).
KF_TEST(analyzer_processes_first_event) {
  const int nEvents = 10;
  writeLambdaFile("test_analyzer_first_in.root", nEvents, 11);
  {
    KFitAnalyzer analyzer("test_analyzer_first_in.root",
                          "test_analyzer_first_out.root", -1);
    analyzer.doFitterTask("Mass", {kPidProton, kPidPion}, toy::kMassLambda);
  }
  const OutputSummary s = readOutput("test_analyzer_first_out.root");
  KF_CHECK(s.entries == nEvents);
  KF_CHECK(s.firstEvent == 0);
}

// nEvents limits the number of analysed events.
KF_TEST(analyzer_event_limit) {
  writeLambdaFile("test_analyzer_limit_in.root", 20, 12);
  {
    KFitAnalyzer analyzer("test_analyzer_limit_in.root",
                          "test_analyzer_limit_out.root", 5);
    analyzer.doFitterTask("Mass", {kPidProton, kPidPion}, toy::kMassLambda);
  }
  const OutputSummary s = readOutput("test_analyzer_limit_out.root");
  KF_CHECK(s.entries == 5);
}

// An unknown task is reported and skipped; it used to reach
// KinFitter::fit() without a constraint, which called abort().
KF_TEST(builder_unknown_task_no_abort) {
  TRandom3 rnd(13);
  const auto ev = toy::lambdaEvent(rnd);
  std::vector<std::vector<KinFitParticle>> cands = {
      {toy::exactTrack(ev.proton, ev.decayVertex, toy::Resolution())},
      {toy::exactTrack(ev.pion, ev.decayVertex, toy::Resolution())}};
  KFitDecayBuilder builder("NoSuchTask", {kPidProton, kPidPion});
  builder.setInputCands(cands);
  builder.countCombis();
  builder.buildDecay();
  KF_CHECK_NEAR(builder.getProbability(), 0, 0);
}

KF_TEST_MAIN()
