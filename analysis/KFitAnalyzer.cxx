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

#include "KFitAnalyzer.h"
#include "KFitDecayBuilder.h"

#include "TFile.h"

/**
 * @file KFitAnalyzer.cxx
 * @brief Implementation of the high-level analysis workflow wrapper.
 */

KFitAnalyzer::KFitAnalyzer(TString inFileName, TString outFileName, int nEvents)
    : fEvents(nEvents), fVerbose(0),
      fCands_in(new TClonesArray("KinFitParticle")) {
  if (fVerbose > 0) {
    std::cout << "--------------- KFitAnalyzer() -----------------"
              << std::endl;
  }

  fInFile.reset(TFile::Open(inFileName, "READ"));
  if (!fInFile || fInFile->IsZombie()) {
    std::cerr << "KFitAnalyzer: cannot open input file " << inFileName
              << std::endl;
    fInFile.reset();
  } else {
    fInFile->GetObject("data", fTree);
    if (!fTree)
      std::cerr << "KFitAnalyzer: no TTree 'data' in " << inFileName
                << std::endl;
  }

  foutFile.reset(TFile::Open(outFileName, "RECREATE"));
  if (!foutFile || foutFile->IsZombie()) {
    std::cerr << "KFitAnalyzer: cannot open output file " << outFileName
              << std::endl;
    foutFile.reset();
  } else {
    // Created in the output file's directory, which owns and deletes it.
    fTree_out = new TTree("data_fitted", "output_tree");
  }
}

KFitAnalyzer::~KFitAnalyzer() {
  // Detach the input tree from fCands_in before either is destroyed.
  if (fTree)
    fTree->ResetBranchAddresses();
  if (foutFile && foutFile->IsOpen())
    foutFile->Close(); // deletes fTree_out
  if (fInFile)
    fInFile->Close(); // deletes fTree
}

// Select and sort particles according to their PID.
void KFitAnalyzer::selectCandidates() {
  if (fVerbose > 0) {
    std::cout << " ----------- KFitAnalyzer::selectCandidates() -----------"
              << std::endl;
    std::cout << "" << std::endl;
  }

  int ntracks = fCands_in->GetEntries();

  fCandsFit.clear();
  std::vector<KinFitParticle> tempVec;

  for (size_t it = 0; it < fPids.size(); it++) {
    tempVec.clear();
    for (int j = 0; j < ntracks; j++) {
      KinFitParticle *cand = (KinFitParticle *)fCands_in->At(j);

      if (cand->getPid() == fPids[it]) {
        tempVec.push_back(*cand);
      }
    }
    fCandsFit.push_back(tempVec);

  } // end of PIDs loop
}

void KFitAnalyzer::doFitterTask(TString task, std::vector<int> pids,
                                double mass, TLorentzVector lv,
                                KinFitParticle mother) {
  if (fVerbose > 0) {
    std::cout << " ----------- KFitAnalyzer::selectCandidates() -----------"
              << std::endl;
    std::cout << "" << std::endl;
    std::cout << "Task added: " << task << std::endl;
    std::cout << "" << std::endl;
  }

  if (!fTree || !fTree_out) {
    std::cerr << "KFitAnalyzer::doFitterTask: no input tree or output already "
                 "written (call doFitterTask once per KFitAnalyzer)"
              << std::endl;
    return;
  }

  // Read input tree
  int Event;
  TClonesArray *candsIn = fCands_in.get();
  fTree->SetBranchAddress("KinFitParticle", &candsIn);

  // Create output tree
  fTree_out->SetName("data_fitted");
  TString out_branchname = "cands_fitted";
  // Owned here; the output tree only stores its address and is written and
  // deleted (finish()) before this goes out of scope.
  std::unique_ptr<TClonesArray> fittedCandsOwner(
      new TClonesArray("KinFitParticle"));
  TClonesArray *fitted_cands = fittedCandsOwner.get();
  TClonesArray &fit_arrayRef = *fitted_cands;
  double Chi2, Prob;

  fTree_out->Branch("Event", &Event, "Event/I");
  fTree_out->Branch("KinFitParticle", "TClonesArray", &fitted_cands);
  fTree_out->Branch("Chi2", &Chi2, "Chi2/D");
  fTree_out->Branch("Prob", &Prob, "Prob/D");

  setPids(pids);

  int entries = fTree->GetEntries();
  if (fEvents < entries && fEvents > 0)
    entries = fEvents;
  if (fVerbose > 1) {
    std::cout << "EntEvents to be analyzed: " << entries << std::endl;
  }

  KFitDecayBuilder builder(task, fPids, lv, mass);
  builder.setVerbosity(fVerbose);

  // start of the event loop
  if (fVerbose > 1) {
    std::cout << "Start of event loop" << std::endl;
  }
  for (int i = 1; i < entries; i++) {
    fTree->GetEntry(i);
    Event = i;

    selectCandidates();

    // if not all particles are found, skip event
    bool isIncomplete = false;
    for (size_t it = 0; it < fPids.size(); it++) {
      if (fCandsFit[it].size() == 0) {
        isIncomplete = true;
        break;
      }
    }
    if (isIncomplete)
      continue;

    // initialize DecayBuilder
    builder.setInputCands(fCandsFit);
    builder.countCombis();

    builder.buildDecay();

    // get result from DecayBuilder
    std::vector<KinFitParticle> result;
    if (result.size() > 0)
      result.clear();
    builder.getFitCands(result);
    Chi2 = builder.getChi2();
    Prob = builder.getProbability();

    // Fill output TClonesArray with fit result
    int ii = 0;
    if (fitted_cands->GetEntries() > 0)
      fitted_cands->Clear();
    if (Prob > 0) {
      for (int k = 0; k < result.size(); k++) {
        KinFitParticle *fitted_cand = new (fit_arrayRef[ii])
            KinFitParticle(result[k], result[k].getR(), result[k].getZ());
        fitted_cand->setPid(result[k].getPid());
        fitted_cand->setTrackId(result[k].getTrackId());
        fitted_cand->setCovariance(result[k].getCovariance());

        ii++;
      }

      fTree_out->Fill();
      fitted_cands->Clear();
    } else {
      cout << "no candidate found" << endl;
      continue;
    }

  } // end of event loop

  if (fVerbose > 0) {
    std::cout << "Event loop finished!" << std::endl;
  }
  finish();
  fTree->ResetBranchAddresses();
}

// Write output and close all resources.
void KFitAnalyzer::finish() {
  foutFile->cd();
  fTree_out->Write();
  foutFile->Close(); // deletes fTree_out
  fTree_out = nullptr;
}