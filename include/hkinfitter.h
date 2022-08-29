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

// ROOT includes
#include "TObject.h"

// framework includes
#include "hrefitcand.h"

using std::cout;
using std::endl;

const Double_t pi2 = TMath::PiOver2();

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

class HKinFitter : public TObject
{
private:
    TMatrixD y, x, V, Vx, fPull;
    Double_t fChi2, fProb;
    Bool_t fConverged;
    Int_t fIteration, fN, fyDim;
    std::vector<HRefitCand> fCands;
    HRefitCand fMother;
    TLorentzVector fMissDaughter;

    // data members for constraints
    Int_t fNdf;
    std::vector<Double_t> fM;
    TLorentzVector fInit;
    Double_t fMass;

    bool fMassConstraint, fMMConstraint, fMassVtxConstraint, fVtxConstraint, f3Constraint, f4Constraint, fMomConstraint;
    int fVerbose;

    Double_t fLearningRate;
    Int_t fNumIterations;
    Double_t fConvergenceCriterion;

public:
    HKinFitter(const std::vector<HRefitCand> &cands);
    ~HKinFitter(){};

    TMatrixD calcMissingMom(const TMatrixD &m_iter);

    TMatrixD f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter);
    TMatrixD Feta_eval(const TMatrixD &miter, const TMatrixD &xi_iter);
    TMatrixD Fxi_eval(const TMatrixD &miter, const TMatrixD &xi_iter);

    void add3Constraint(HRefitCand mother);
    void add4Constraint(TLorentzVector lv);
    void addVertexConstraint();
    void addMomConstraint(TLorentzVector lv, Double_t mass);
    void addMassConstraint(Double_t mass);
    void addMMConstraint(Double_t mm, TLorentzVector init);
    void addMassVtxConstraint(Double_t mass);

    void setLearningRate(Double_t val) { fLearningRate = val; }
    void setNumberOfIterations(Int_t val) { fNumIterations = val; }
    void setConvergenceCriterion(Double_t val) { fConvergenceCriterion = val; }
    void setCovariance(TMatrixD &val) { V = val; }
    void setMeasurement(TMatrixD &val) { y = val; }

    Double_t getChi2() const { return fChi2; }
    Double_t getProb() const { return fProb; }
    Double_t getPull(Int_t val = 0) { return fPull(val, val); }
    Int_t getIteration() const { return fIteration; }

    Bool_t isConverged() const { return fConverged; }

    Bool_t fit();

    void setVerbosity(Int_t val) { fVerbose = val; }

    HRefitCand getDaughter(Int_t val);
    void getDaughters(std::vector<HRefitCand> &daughters) { daughters = fCands; }
    HRefitCand getMother();
    TLorentzVector getMissingDaughter();

    void update();

protected:
    void updateDaughters();
    void updateMother();
    ClassDef(HKinFitter, 0)
};

#endif /* HKINFITTER_H */