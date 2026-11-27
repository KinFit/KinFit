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
 * KinFitter.h
 *
 * @updated 03.08.2023
 * @version v1.0.0
 *
 * Main class that performs the kinematic fit
 *
 */

#pragma once

// framework includes
#include "KinFitParticle.h"

// ROOT includes
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObject.h"

// system includes
#include <iomanip>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

const double pi2 = TMath::PiOver2();

template <typename T> void Print(T const &matrix) {
  int nrows = matrix.GetNrows();
  int ncols = matrix.GetNcols();

  cout << "shape(" << nrows << "," << ncols << ")" << endl;

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      double element = matrix(i, j);
      if (TMath::Abs(element) < 1e-10)
        element = 0.;
      if (element >= 0.)
        cout << " " << std::fixed << std::setw(8) << std::scientific << element
             << " ";
      else
        cout << std::fixed << std::setw(8) << std::scientific << element << " ";
    }
    cout << endl;
  }

  cout << endl;
}

class KinFitter : public TObject {
private:
  TMatrixD y;              // Vector of measured variables (input, kept)
  TMatrixD x;              // Vector of unmeasured variables (fitted)
  TMatrixD V;              // Covariance matrix of measured variables (input)
  TMatrixD Vx;             // Covariance matrix of unmeasured variables
  TMatrixD fYFit;          //! Fitted measured variables
  TMatrixD fVFit;          //! Fitted covariance of the measured variables
  TMatrixD fVxy;           //! Fitted covariance unmeasured x measured
  std::vector<double> fPull; // Pull value for all measured variables
  double fChi2 = 1e6;      // Chi2 of fit
  double fProb = 1e6;      // Probability of fit
  bool fConverged = false; // True if fit has converged
  int fIteration = 0;      // Iterations needed until convergence
  int fN;                  // Number of input candidates
  int fyDim;               // Dimension of y
  std::vector<KinFitParticle> fCands; // Vector of input candidates
  KinFitParticle fMother;             // Decaying particle
  TLorentzVector fMissDaughter;       // Missing decay product
  std::vector<double> fMassVec; // One mass per addMassConstraint() call

  // Data members for constraints
  int fNdf = 0;                          // Number of degrees of freedom
  std::vector<double> fM;                // Vector of particle masses
  TLorentzVector fInit;                  // 4-vector used for constraint
  double fMass = 99999.9;                // Mass used for constraint
  double fMassMissingParticle = 99999.9; // Mass used for constraint

  // Constraints, true if constraint is set, only one at a time
  bool fMassConstraint = false;
  bool fMissingMassConstraint = false;
  bool fMassVtxConstraint = false;
  bool fVtxConstraint = false;
  bool f3Constraint = false;
  bool f4Constraint = false;
  bool fMissingParticleConstraint = false;

  int fNumIterations = 20; // Maximum number of iterations

  // Convergence criteria: difference in chi2, constraint equation d, difference
  // in track parameters between iterations If default parameters are used only
  // the chi2 will be used for convergence
  double fConvergenceCriterionChi2 = 1e-4;
  double fConvergenceCriterionD = 1e6;
  double fConvergenceCriterionAlpha = 1e6;

  std::vector<int> fFlexiParticlesInFit; // selection for the next mass constraint
  std::vector<std::vector<int>> fMassFitPair; // particles per mass constraint,
                                              // empty = default selection

  int fVerbose = 0;

  int fNExclusive = -1;

  TMatrixD calcMissingMom(const TMatrixD &m_iter);

  /** Particles entering mass constraint k (explicit selection, else the
   * first fNExclusive particles, else all). */
  std::vector<int> massParticles(size_t k) const;
  /** Number of unmeasured parameters (rows of x) used by the fit. */
  int nUnmeasured() const;

public:
  /** @brief Constructor
   * @param cands Vector of particles to be fitted.
   */
  KinFitter(const std::vector<KinFitParticle> &cands);

  /** @brief Default destructor. */
  ~KinFitter() {};

  /** @brief Evaluation of constraint equations
   * @param m_iter Measured track parameters.
   * @param xi_iter Unmeasured track parameters.
   * @return Constraint vector evaluated at the current iteration.
   */
  TMatrixD f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter);

  /** @brief Evaluation of Jacobian w.r.t. measured track parameters
   * @param miter Measured track parameters.
   * @param xi_iter Unmeasured track parameters.
   * @return Jacobian matrix with respect to measured parameters.
   */
  TMatrixD Feta_eval(const TMatrixD &miter, const TMatrixD &xi_iter);

  /** @brief Evaluation of Jacobian w.r.t. unmeasured track parameters
   * @param miter Measured track parameters.
   * @param xi_iter Unmeasured track parameters.
   * @return Jacobian matrix with respect to unmeasured parameters.
   */
  TMatrixD Fxi_eval(const TMatrixD &miter, const TMatrixD &xi_iter);

  /** @brief Functions for selecting the constraint equation.
   * @param mother Mother particle.
   * @param lv Initial beam-target four-momentum.
   * @param mass Mass of particle or missing mass of system (depending on
   * constraint)
   */
  void add3Constraint(
      KinFitParticle mother); // 4-momentum constraint in a decay vertex
  void
  add4Constraint(TLorentzVector lv); // 4-momentum constraint of final state
                                     // particles to the initial system, lv
  void addVertexConstraint();        // Geometrical vertex constraint
  void addMissingParticleConstraint(
      TLorentzVector lv,
      double
          mass); // Constraint of the final state particles + one undetected to
                 // the initial system, lv, and a missing particle with mass m
  void addMassConstraint(
      double mass); // Constraining decay products from decaying particle to the
                    // mass, m, of the decaying particle
  void addMissingMassConstraint(
      TLorentzVector lv,
      double mass); // Constrains all final state particles to the initial
                    // system, lv, and a missing mass , m
  void addMassVtxConstraint(
      double mass); // Constrains all particles to a common vertex and a mass,
                    // m, of a decaying particle

  /** @brief Mass constraints without explicit selection use only the first
   * `number` particles. */
  void setNumberOfExclusiveCandidates(int number) { fNExclusive = number; }
  /** @brief Select particle `number` for the *next* addMassConstraint() call.
   * Call once per particle, then addMassConstraint(); repeat for further
   * mass constraints on other particle subsets. */
  void setUseCandInFit(int number) { fFlexiParticlesInFit.push_back(number); };

  /** @brief Set maximum number of iterations, default = 20
   */
  void setNumberOfIterations(int val) { fNumIterations = val; }

  /** @brief Set convergence criteria
   * @param val1 Difference in chi2 of consecutive iterations.
   * @param val2 Norm of all constraint equations.
   * @param val3 Difference in norm of track-parameter vectors of consecutive
   * iterations
   */
  void setConvergenceCriteria(double val1, double val2, double val3) {

    fConvergenceCriterionChi2 = val1;
    fConvergenceCriterionD = val2;
    fConvergenceCriterionAlpha = val3;
  }

  void setConvergenceCriterionChi2(double val) {
    fConvergenceCriterionChi2 = val;
  }

  /** @brief Return chi2 of the fit.
   * @return Fit chi2.
   */
  double getChi2() const { return fChi2; }

  /** @brief Return fit probability.
   * @return Fit probability.
   */
  double getProb() const { return fProb; }

  /** @brief Number of degrees of freedom (constraints - unmeasured
   * parameters), set by fit(). */
  int getNDF() const { return fNdf; }

  /** @brief Return pull value for one fit parameter.
   * val = 0 + 5 * n : 1/p
   * val = 1 + 5 * n : theta
   * val = 2 + 5 * n : phi
   * val = 3 + 5 * n : R
   * val + 4 + 5 * n : Z
   * n = particle at entry 0, 1, 2... in the input vector
   * @param val Index of the pull entry.
   * @return Pull value (diagonal element of the pull matrix).
   */
  double getPull(int val = 0) const {
    return (val >= 0 && val < (int)fPull.size()) ? fPull[val] : -10000.;
  }
  int getIteration() const { return fIteration; }

  /** @brief Returns true if fit converged within max number of iterations
   */
  bool isConverged() const { return fConverged; }

  /** @brief Main fit function.
   *
   * Applies an iterative fit with Lagrange multipliers, updates covariance
   * and track parameters, and computes pull values.
   * @return True if the fit converged successfully.
   */
  bool fit();

  void setVerbosity(int val) { fVerbose = val; }

  /** @brief Return fitted daughter particle with updated parameters.
   * @param val Index of daughter particle.
   * @return Updated daughter particle.
   */
  KinFitParticle getDaughter(int val) { return fCands[val]; }
  /** @brief Return all fitted daughter particles.
   * @param daughters Output vector filled with all fitted daughters.
   */
  void getDaughters(std::vector<KinFitParticle> &daughters) {
    daughters = fCands;
  }
  /** @brief Return fitted mother particle.
   * @return Fitted mother particle.
   */
  KinFitParticle getMother() { return fMother; }
  /** @brief Return fitted missing particle four-vector if available.
   * @return Missing daughter four-vector.
   */
  TLorentzVector getMissingDaughter() { return fMissDaughter; }

private:
  void updateDaughters();
  void updateMother();

  ClassDef(KinFitter, 2)
};
