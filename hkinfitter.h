/**
 * hkinfitter.h
 *
 *
 */

#ifndef HKINFITTER_H
#define HKINFITTER_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// framework includes
#include "hrefitcand.h"
#include "hgeomvector.h"
#include "hparticletool.h"

using std::cout;
using std::endl;

const double pi2 = TMath::PiOver2();

template <typename T>
void Print(T const &matrix)
{
    Int_t nrows = matrix.GetNrows();
    Int_t ncols = matrix.GetNcols();

    cout << "shape(" << nrows << "," << ncols << ")" << endl;

    for (Int_t i = 0; i < nrows; i++)
    {
        for (Int_t j = 0; j < ncols; j++)
        {
            Double_t element = matrix(i, j);
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

class HKinFitter
{
private:
    TMatrixD y, x, V, Vx, fPull;
    double fChi2, fProb;
    bool fConverged;
    int fIteration, fN, fyDim;   
    std::vector<HRefitCand> fCands;
    HRefitCand fMother;
    TLorentzVector fMissDaughter;

    // data members for constraints
    int fNdf;
    std::vector<double> fM;
    TLorentzVector fInit;
    Double_t fMass;

    bool fVtxConstraint, f3Constraint, f4Constraint, fMomConstraint;
    int fVerbose;

    double fLearningRate;
    int fNumIterations; 
    double fConvergenceCriterion;

public:
    HKinFitter(const std::vector<HRefitCand> &cands);
    HKinFitter(const std::vector<HRefitCand> &cands, HRefitCand &mother);
    HKinFitter(const std::vector<HRefitCand> &cands, TLorentzVector &lv);
    HKinFitter(const std::vector<HRefitCand> &cands, TLorentzVector &lv, Double_t mass);
    ~HKinFitter(){};

    TMatrixD calcMissingMom(const TMatrixD &m_iter);

    TMatrixD f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter);
    TMatrixD Feta_eval(const TMatrixD &miter, const TMatrixD &xi_iter);
    TMatrixD Fxi_eval(const TMatrixD &miter, const TMatrixD &xi_iter);

    void add3Constraint();
    void add4Constraint();
    void addVertexConstraint();
    void addMomConstraint();

    void setLearningRate(double val) { fLearningRate = val; }
    void setNumberOfIterations(int val) { fNumIterations = val; }
    void setConvergenceCriteria(double val) { fConvergenceCriterion = val; }

    double getChi2() const { return fChi2; }
    double getProb() const { return fProb; }
    double getPull(int val = 0) { return fPull(val, val); }

    bool isConverged() const { return fConverged; }
    int getIteration() const { return fIteration; }
    void setCovariance(TMatrixD &val) { V = val; }
    void setMeasurement(TMatrixD &val) { y = val; }

    bool fit();

    void setVerbosity(int val) { fVerbose = val; }

    HRefitCand getDaughter(int val);
    void getDaughters(std::vector<HRefitCand> &daughters) { daughters = fCands; }
    HRefitCand getMother();
    TLorentzVector getMissingDaughter();

    void update();

protected:
    void updateDaughters();
    void updateMother();
};

#endif /* HKINFITTER_H */
