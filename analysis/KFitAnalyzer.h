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
 * KFitAnalyzer.h
 *
 * @updated 01.08.2023
 * @version v1.0.0
 *
 * "User interface" class to automatically analyze and fit ROOT files.
 * Particle candidates must be stored in a TClonesArray of KinFitParticle
 * objects.
 *
 */

#pragma once

#include "KinFitParticle.h"

// framework includes
#include "TClonesArray.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"

// system includes
#include <iostream>
#include <memory>
#include <vector>

using std::cout;
using std::endl;

class KFitAnalyzer {
private:
  std::vector<std::vector<KinFitParticle>>
      fCandsFit; // Vector of vector of particle candidates for each PID
  std::vector<int> fPids; // Vector of allPIDs of particles included in fit
  int fEvents;            // Number of events to be analyzed

  int fVerbose; // Verbosity

  // Read in data
  std::unique_ptr<TFile> fInFile;         // Input file (owns fTree)
  TTree *fTree = nullptr;                 // Input data tree
  std::unique_ptr<TClonesArray> fCands_in; // Input KinFitParticles

  // Output data
  std::unique_ptr<TFile> foutFile; // Output file (owns fTree_out)
  TTree *fTree_out = nullptr;      // Output data tree, nullptr once written

  /** @brief Return the output tree with fit results.
   * @return Pointer to the output TTree.
   */
  TTree *getFittedTree() { return fTree_out; }

  /** @brief Select and sort candidates according to their PID.
   */
  void selectCandidates();

  /** @brief Write the output tree to file and close resources.
   */
  void finish();

public:
  /** @brief Constructor
   * @param inFileName Full file name of the input ROOT file.
   * @param outFileName Full file name of the output ROOT file.
   * @param nEvents Number of events to analyze. If negative, process all
   * available events.
   */
  KFitAnalyzer(TString inFileName, TString outFileName, int nEvents = -1);

  /** @brief Closes the input and (if still open) the output file. */
  ~KFitAnalyzer();

  KFitAnalyzer(const KFitAnalyzer &) = delete;
  KFitAnalyzer &operator=(const KFitAnalyzer &) = delete;

  void setVerbosity(int val) { fVerbose = val; }

  //---------------User functions---------------------------------------

  /** @brief Perform the automated fitting procedure.
   *
   * Reads the input tree, initializes the decay builder, runs the event loop,
   * executes the fit, and fills the output tree with the fit result.
   *
   * @param task Task type to run (for example: "4C", "Vertex", "Mass").
   * @param pids Vector of PIDs of all particles included in the fit.
   * @param mass Optional mass or missing mass used by mass-constrained tasks.
   * @param lv Optional four-vector input required by specific tasks.
   * @param mother Optional mother particle used by specific tasks.
   */
  void doFitterTask(TString task, std::vector<int> pids, double mass = -1.,
                    TLorentzVector lv = TLorentzVector(),
                    KinFitParticle mother = KinFitParticle());
  // void addFitterTask(TString task, std::vector<int> primPids,
  // std::vector<int> decayPids); // Jenny, for 3C fit void
  // addBuilderTask(TString task, std::vector<int> pids, TLorentzVector lv);

  /** @brief Set PIDs used for the fit.
   * @param val Vector of PIDs of particles included in the fit.
   */
  void setPids(std::vector<int> val) { fPids = val; }

  /** @brief Return the configured PID vector.
   * @return Vector of PIDs used in the fit.
   */
  std::vector<int> getPids() { return fPids; }
};
