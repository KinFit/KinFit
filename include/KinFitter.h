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
    KinFitter(const std::vector<KFitParticle> &cands);
    ~KinFitter(){};

    TMatrixD f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter);
    TMatrixD Feta_eval(const TMatrixD &miter, const TMatrixD &xi_iter);
    TMatrixD Fxi_eval(const TMatrixD &miter, const TMatrixD &xi_iter);

    void add3Constraint(KFitParticle mother);
    void add4Constraint(TLorentzVector lv);
    void addVertexConstraint();
    void addMomConstraint(TLorentzVector lv, double mass);
    void addMassConstraint(double mass);
    void addMMConstraint(double mm, TLorentzVector init);
    void addMassVtxConstraint(double mass);

    void setNumberOfIterations(int val) { fNumIterations = val; }
    void setConvergenceCriteria(double val1, double val2, double val3);
    void setCovariance(TMatrixD &val) { V = val; }
    void setMeasurement(TMatrixD &val) { y = val; }

    double getChi2() const { return fChi2; }
    double getProb() const { return fProb; }
    double getPull(int val = 0) { return fPull(val, val); }
    int getIteration() const { return fIteration; }

    bool isConverged() const { return fConverged; }

    bool fit();

    void setVerbosity(int val) { fVerbose = val; }

    KFitParticle getDaughter(int val);
    void getDaughters(std::vector<KFitParticle> &daughters) { daughters = fCands; }
    KFitParticle getMother();
    TLorentzVector getMissingDaughter();

    void update();

protected:
    void updateDaughters();
    void updateMother();
    ClassDef(KinFitter, 1)
};

#endif /* KINFITTER_H */
