/**
 * KinFitter.h
 *
 * @updated 31.05.2023
 * @version v1.0.0
 * 
 * Main class that performs the kinematic fit
 *
 */

#ifndef KINFITTER_H
#define KINFITTER_H

// framework includes
#include "KFitParticle.h"

// ROOT includes
#include "TObject.h"

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using std::cout;
using std::endl;

const double pi2 = TMath::PiOver2();

template <typename T>
void Print(T const &matrix)
{
    int nrows = matrix.GetNrows();
    int ncols = matrix.GetNcols();

    cout << "shape(" << nrows << "," << ncols << ")" << endl;

    for (int i = 0; i < nrows; i++)
    {
        for (int j = 0; j < ncols; j++)
        {
            double element = matrix(i, j);
            if (TMath::Abs(element) < 1e-10)
                element = 0.;
            if (element >= 0.)
                cout << " " << std::fixed << std::setw(8) << std::scientific << element << " ";
            else
                cout << std::fixed << std::setw(8) << std::scientific << element << " ";
        }
        cout << endl;
    }

    cout << endl;
}

class KinFitter : public TObject
{
private:
    TMatrixD y; // Vector of measured variables
    TMatrixD x;  // Vector of unmeasured variables
    TMatrixD V;  // Covariance matrix of measured variables
    TMatrixD Vx;  // Covariance matrix of unmeasured variables
    TMatrixD fPull;  // Pull value for all measured variables
    double fChi2;    // Chi2 of fit
    double fProb;   // Probability of fit
    bool fConverged;    // True if fit has converged
    int fIteration; // Iterations needed until convergence
    int fN;     // Number of input candidates
    int fyDim;  // Dimension of y
    std::vector<KFitParticle> fCands;
    KFitParticle fMother;
    TLorentzVector fMissDaughter;

    // data members for constraints
    int fNdf;   // Number of degrees of freedom
    std::vector<double> fM; // Vector of particle masses
    TLorentzVector fInit;   // 4-vector used for constraint
    double fMass;   // Mass used for constraint

    // Constraints, true if constraint is set, only one at a time
    bool fMassConstraint, fMMConstraint, fMassVtxConstraint, fVtxConstraint, f3Constraint, f4Constraint, fMomConstraint;

    int fNumIterations; //Maximum number of iterations
    // Convergence criteria: difference in chi2, constraint equation d, difference in track parameters between iterations
    double fConvergenceCriterionChi2, fConvergenceCriterionD, fConvergenceCriterionAlpha;

    int fVerbose;

    TMatrixD calcMissingMom(const TMatrixD &m_iter);


public:
    /** @brief Constructor
    * @param cands - vector of particles to be fitted
    */
    KinFitter(const std::vector<KFitParticle> &cands);

    /** Default Deconstructor **/
    ~KinFitter(){};

    /** @brief Evaluation of constraint equations
    * @param m_iter - measured track parameters
    * @param xi_iter - unmeasured track parameters    
    */
    TMatrixD f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter);

    /** @brief Evaluation of Jacobian w.r.t. measured track parameters
    * @param m_iter - measured track parameters
    * @param xi_iter - unmeasured track parameters    
    */
    TMatrixD Feta_eval(const TMatrixD &miter, const TMatrixD &xi_iter);

    /** @brief Evaluation of Jacobian w.r.t. unmeasured track parameters
    * @param m_iter - measured track parameters
    * @param xi_iter - unmeasured track parameters    
    */
    TMatrixD Fxi_eval(const TMatrixD &miter, const TMatrixD &xi_iter);


    /** @brief Functions for choice of constraint equation
    * @param mother - mother particle
    * @param lv - initial beam-target 4-momentum  
    * @param mass - mass of particle / missing mass of system (depending on constraint)  
    */
    void add3Constraint(KFitParticle mother);
    void add4Constraint(TLorentzVector lv);
    void addVertexConstraint();
    void addMomConstraint(TLorentzVector lv, double mass);
    void addMassConstraint(double mass);
    void addMMConstraint(double mass, TLorentzVector lv);
    void addMassVtxConstraint(double mass);


    /** @brief Set maximum number of iterations, default = 20
    */
    void setNumberOfIterations(int val) { fNumIterations = val; }

    /** @brief Set convergence criteria
     * @param val1 - difference in chi2 of consecutive iterations
     * @param val2 - Norm of all constraint equations
     * @param val2 - Difference in norm of track parameter vector of consecutive iterations
    */
    void setConvergenceCriteria(double val1, double val2, double val3);
    //void setCovariance(TMatrixD &val) { V = val; }
    //void setMeasurement(TMatrixD &val) { y = val; }

    double getChi2() const { return fChi2; }
    double getProb() const { return fProb; }
    double getPull(int val = 0) { return fPull(val, val); }
    int getIteration() const { return fIteration; }

    /** @brief Returns true if fit converged within max number of iterations
    */
    bool isConverged() const { return fConverged; }

    /** @brief Main fit function, iterative fitting procedure using 
     * Lagrange multiplyers is applied, covariance and track parameters 
     * are updated, pulls are calculated
    */
    bool fit();

    void setVerbosity(int val) { fVerbose = val; }

    /** @brief Fitted particle number val is returned with updated track parameters
    */
    KFitParticle getDaughter(int val);
    /** @brief All fitted particles are returned 
    */
    void getDaughters(std::vector<KFitParticle> &daughters) { daughters = fCands; }
    /** @brief Returns fitted mother particle
    */
    KFitParticle getMother();
    /** @brief Returns TLorentzVector fitted missing particle if existing
    */
    TLorentzVector getMissingDaughter();

    //void update();

protected:
    void updateDaughters();
    void updateMother();
    ClassDef(KinFitter, 1)
};

#endif /* KINFITTER_H */
