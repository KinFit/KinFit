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
 * KFitDecayBuilder.h
 *
 * @updated 01.08.2023
 * @version v1.0.0
 *
 * Class responsible for combinatorics in each event of the automated
 * fitting procedure. Calls KinFitter and chooses the best combination.
 *
 */

#pragma once

// framework includes
#include "KinFitParticle.h"

#include "TString.h"

// system includes
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

class KFitDecayBuilder {
private:
  // Working Particles
  std::vector<std::vector<KinFitParticle>>
      fCands; // Vector of vector of particle candidates for each PID
  std::vector<KinFitParticle>
      fFitCands; // Vector of particles the fit should be performed on

  std::vector<KinFitParticle> fOutputCands; // Output particles after fitting

  // Fitter input variables
  TString fTask;          // Task to be performed
  std::vector<int> fPids; // Vector of PIDs of all particles included in fit
  TLorentzVector fIniSys; // Four-vector as input to fit
  KinFitParticle fMother; // Mother particle as input to fit
  double fMass;           // Mass as input to fit

  // For combinatorics
  int fTotalCombis;                 // Total number of combinations for event
  int fCombiCounter;                // Number of combination that is evaluated
  std::vector<int> particleCounter; // Vector, each entry is a counter for
                                    // respective PID of PID vector
  bool doubleParticle; // True if same particle was used twice in same
                       // combination

  double fBestProb = 0;   // Probability of best combination
  double fBestChi2 = 1e6; // Chi2 of best combination

  int fVerbose; // Verbosity 0-1-high

  /** @brief Fills vector of particles to be fitted with current combination */
  void fillFitCands();

  /** @brief Calls KFitter */
  bool doFit();

  /** @brief Checks if same particle was used twice in same combination */
  void checkDoubleParticle(size_t i);

public:
  /** @brief Constructor
   * @param task String defining the task to be performed. Examples: 4C,
   * Vertex, Mass.
   * @param pids Vector of PIDs of all particles included in the fit.
   * @param lv Optional input four-vector needed by specific tasks.
   * @param mass Optional mass or missing mass used by specific tasks.
   */
  KFitDecayBuilder(TString task, std::vector<int> pids,
                   TLorentzVector lv = TLorentzVector(), double mass = -1.);

  /** @brief Default destructor. */
  ~KFitDecayBuilder() {};

  void setVerbosity(int val) { fVerbose = val; }

  // setters
  /** @brief Set input candidates
   * @param cands Vector of vector of particle candidates for each PID
   */
  void setInputCands(const std::vector<std::vector<KinFitParticle>> &cands) {
    fCands = cands;
  }

  /** @brief Set beam + target four-momentum.
   * @param val Beam + target four-momentum.
   */
  void setIniSys(TLorentzVector val) { fIniSys = val; }

  /** @brief Set mother particle.
   * @param val Mother KinFitParticle.
   */
  void setMother(KinFitParticle val) { fMother = val; }

  /** @brief Set input mass.
   * @param val Mass value.
   */
  void setMass(double val) { fMass = val; }

  /** @brief Count the number of particle combinations in the event. */
  void countCombis();

  /** @brief Run combinatorics, call the fitter, and choose the best result. */
  void buildDecay();

  /** @brief Access the selected fitted particles.
   * @param cands Output vector filled with fitted particles.
   */
  void getFitCands(std::vector<KinFitParticle> &cands) { cands = fOutputCands; }

  /** @brief Return chi2 of the best combination.
   * @return Best-fit chi2.
   */
  double getChi2() { return fBestChi2; }

  /** @brief Return probability of the best combination.
   * @return Best-fit probability.
   */
  double getProbability() { return fBestProb; }
};

