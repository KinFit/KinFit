#include "hkinfitter.h"

const size_t cov_dim = 5;

//To be used for vertex fitter
HKinFitter::HKinFitter(const std::vector<HRefitCand> &cands) : fCands(cands),
                                                                     fVerbose(0),
                                                                     fLearningRate(1),
                                                                     fNumIterations(10),
                                                                     fConvergenceCriteria(1)
{
    // fN is the number of daughters e.g. (L->ppi-) n=2
    fN = cands.size();
    fyDim = fN * cov_dim; //Dimension of full covariance matrix (number of measured variables x cov_dim)

    y.ResizeTo(fyDim, 1);
    x.ResizeTo(1, 1);
    V.ResizeTo(fyDim, fyDim);
    Vx.ResizeTo(1, 1);

    y.Zero();
    x.Zero();
    V.Zero();
    Vx.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;
    f3Constraint = false;
    f4Constraint = false;
    fVtxConstraint = false;
    fMomConstraint = false;

    // set 'y=alpha' measurements
    // and the covariance
    for (int ix = 0; ix < fN; ix++)
    {
        HRefitCand cand = cands[ix];
        y(0 + ix * cov_dim, 0) = 1. / cand.P();
        y(1 + ix * cov_dim, 0) = cand.Theta();
        y(2 + ix * cov_dim, 0) = cand.Phi();
        y(3 + ix * cov_dim, 0) = cand.getR();
        y(4 + ix * cov_dim, 0) = cand.getZ();
        fM.push_back(cand.M());

        // FIX ME: only for diagonal elements
        TMatrixD covariance = cand.getCovariance();
        V(0 + ix * cov_dim, 0 + ix * cov_dim) = covariance(0, 0);
        V(1 + ix * cov_dim, 1 + ix * cov_dim) = covariance(1, 1);
        V(2 + ix * cov_dim, 2 + ix * cov_dim) = covariance(2, 2);
        V(3 + ix * cov_dim, 3 + ix * cov_dim) = covariance(3, 3);
        V(4 + ix * cov_dim, 4 + ix * cov_dim) = covariance(4, 4);

    }
}

//To be used for 3C fitter
HKinFitter::HKinFitter(const std::vector<HRefitCand> &cands, HRefitCand &mother) : fCands(cands),
                                                                                         fMother(mother),
                                                                                         fVerbose(0),
                                                                                         fLearningRate(1),
                                                                                         fNumIterations(10),
                                                                                         fConvergenceCriteria(1)
{
    // fN is the number of daughters e.g. (L->ppi-) n=2
    fN = cands.size();
    fyDim = (fN + 1) * cov_dim - 1; //Dimension of full covariance matrix (number of measured variables x cov_dim). Mother momentum is not measured

    y.ResizeTo(fyDim, 1);
    x.ResizeTo(1, 1);
    V.ResizeTo(fyDim, fyDim);
    Vx.ResizeTo(1, 1);

    y.Zero();
    V.Zero();
    x.Zero();
    Vx.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;
    f3Constraint = false;
    f4Constraint = false;
    fVtxConstraint = false;
    fMomConstraint = false;

    // set y to measurements and the covariance, set mass
    for (int ix = 0; ix < fN; ix++) //for daughters
    {
        HRefitCand cand = cands[ix];

        y(0 + ix * cov_dim, 0) = 1. / cand.P();
        y(1 + ix * cov_dim, 0) = cand.Theta();
        y(2 + ix * cov_dim, 0) = cand.Phi();
        y(3 + ix * cov_dim, 0) = cand.getR();
        y(4 + ix * cov_dim, 0) = cand.getZ();
        fM.push_back(cand.M());

        // FIX ME: only for diagonal elements
        TMatrixD covariance = cand.getCovariance();
        V(0 + ix * cov_dim, 0 + ix * cov_dim) = covariance(0, 0);
        V(1 + ix * cov_dim, 1 + ix * cov_dim) = covariance(1, 1);
        V(2 + ix * cov_dim, 2 + ix * cov_dim) = covariance(2, 2);
        V(3 + ix * cov_dim, 3 + ix * cov_dim) = covariance(3, 3);
        V(4 + ix * cov_dim, 4 + ix * cov_dim) = covariance(4, 4);
    }

    //for mother
    HRefitCand cand = fMother;

    y(fN * cov_dim, 0) = cand.Theta();
    y(1 + fN * cov_dim, 0) = cand.Phi();
    y(2 + fN * cov_dim, 0) = cand.getR();
    y(3 + fN * cov_dim, 0) = cand.getZ();
    fM.push_back(cand.M());

    TMatrixD covariance = cand.getCovariance();
    V(0 + fN * cov_dim, 0 + fN * cov_dim) = covariance(1, 1);
    V(1 + fN * cov_dim, 1 + fN * cov_dim) = covariance(2, 2);
    V(2 + fN * cov_dim, 2 + fN * cov_dim) = covariance(3, 3);
    V(3 + fN * cov_dim, 3 + fN * cov_dim) = covariance(4, 4);

    V(0 + fN * cov_dim, 1 + fN * cov_dim)=covariance(1, 2);
    V(1 + fN * cov_dim, 0 + fN * cov_dim)=covariance(2, 1);
    // Jenny: Comments below are for testing so that the covariance matrix is read in correctly
    //std::cout << "Fitter" << std::endl;
    //std::cout << "Diag elements: " << V(0 + fN * cov_dim, 0 + fN * cov_dim) << " " << V(1 + fN * cov_dim, 1 + fN * cov_dim) << std::endl;
    //std::cout << "Off diag elements: " << V(0 + fN * cov_dim, 1 + fN * cov_dim) << " " << V(1 + fN * cov_dim, 0 + fN * cov_dim) << std::endl;


}

//To be used for 4C fitter
HKinFitter::HKinFitter(const std::vector<HRefitCand> &cands, TLorentzVector& lv) : fCands(cands),
                                                                                         fInit(lv),
                                                                                         fVerbose(0),
                                                                                         fLearningRate(1),
                                                                                         fNumIterations(10),
                                                                                         fConvergenceCriteria(1)
{
    // fN is the number of daughters e.g. (L->ppi-) n=2
    fN = cands.size();
    fyDim = fN * cov_dim; //Dimension of full covariance matrix (number of measured variables x cov_dim).

    y.ResizeTo(fyDim, 1);
    x.ResizeTo(1, 1);
    V.ResizeTo(fyDim, fyDim);
    Vx.ResizeTo(1, 1);

    y.Zero();
    V.Zero();
    x.Zero();
    Vx.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;
    f3Constraint = false;
    f4Constraint = false;
    fVtxConstraint = false;
    fMomConstraint = false;

    // set y to measurements and the covariance, set mass
    for (int ix = 0; ix < fN; ix++) //for daughters
    {
        HRefitCand cand = cands[ix];

        y(0 + ix * cov_dim, 0) = 1. / cand.P();
        y(1 + ix * cov_dim, 0) = cand.Theta();
        y(2 + ix * cov_dim, 0) = cand.Phi();
        y(3 + ix * cov_dim, 0) = cand.getR();
        y(4 + ix * cov_dim, 0) = cand.getZ();
        fM.push_back(cand.M());

        // FIX ME: only for diagonal elements
        TMatrixD covariance = cand.getCovariance();
        V(0 + ix * cov_dim, 0 + ix * cov_dim) = covariance(0, 0);
        V(1 + ix * cov_dim, 1 + ix * cov_dim) = covariance(1, 1);
        V(2 + ix * cov_dim, 2 + ix * cov_dim) = covariance(2, 2);
        V(3 + ix * cov_dim, 3 + ix * cov_dim) = covariance(3, 3);
        V(4 + ix * cov_dim, 4 + ix * cov_dim) = covariance(4, 4);
    }
}

//To be used for missing particle fitter
HKinFitter::HKinFitter(const std::vector<HRefitCand> &cands, TLorentzVector& lv, Double_t mass) : fCands(cands),
                                                                                         fInit(lv),
                                                                                         fMass(mass),
                                                                                         fVerbose(0),
                                                                                         fLearningRate(1),
                                                                                         fNumIterations(10),
                                                                                         fConvergenceCriteria(1)
{
    // fN is the number of daughters e.g. (L->ppi-) n=2
    fN = cands.size();
    fyDim = fN * cov_dim; //Dimension of full covariance matrix (number of measured variables x cov_dim).

    y.ResizeTo(fyDim, 1);
    x.ResizeTo(3, 1);
    V.ResizeTo(fyDim, fyDim);
    Vx.ResizeTo(3, 3);

    y.Zero();
    V.Zero();
    x.Zero();
    Vx.Zero();

    fConverged = false;
    fIteration = 0;
    fNdf = 0;
    f3Constraint = false;
    f4Constraint = false;
    fVtxConstraint = false;
    fMomConstraint = false;

    // set y to measurements and the covariance, set mass
    for (int ix = 0; ix < fN; ix++) //for daughters
    {
        HRefitCand cand = cands[ix];

        y(0 + ix * cov_dim, 0) = 1. / cand.P();
        y(1 + ix * cov_dim, 0) = cand.Theta();
        y(2 + ix * cov_dim, 0) = cand.Phi();
        y(3 + ix * cov_dim, 0) = cand.getR();
        y(4 + ix * cov_dim, 0) = cand.getZ();
        fM.push_back(cand.M());

        // FIX ME: only for diagonal elements
        TMatrixD covariance = cand.getCovariance();
        V(0 + ix * cov_dim, 0 + ix * cov_dim) = covariance(0, 0);
        V(1 + ix * cov_dim, 1 + ix * cov_dim) = covariance(1, 1);
        V(2 + ix * cov_dim, 2 + ix * cov_dim) = covariance(2, 2);
        V(3 + ix * cov_dim, 3 + ix * cov_dim) = covariance(3, 3);
        V(4 + ix * cov_dim, 4 + ix * cov_dim) = covariance(4, 4);
    }
}

void HKinFitter::add3Constraint()
{
    if (!f3Constraint)
        fNdf += 3;
    f3Constraint = true;
}

void HKinFitter::add4Constraint()
{
    if (!f4Constraint)
        fNdf += 4;
    f4Constraint = true;
}

void HKinFitter::addVertexConstraint()
{
    if (!fVtxConstraint)
        fNdf += 1;
    fVtxConstraint = true;
}

void HKinFitter::addMomConstraint()
{
    if (!fMomConstraint){
        fM.push_back(fMass);
        fNdf += 1;
    }
    fMomConstraint = true;
}

TMatrixD HKinFitter::calcMissingMom(const TMatrixD& m_iter)
{
    TMatrix xi(3,1);

    xi(0, 0) = fInit.Px();
    xi(1, 0) = fInit.Py();
    xi(2, 0) = fInit.Pz();

    for(int q=0; q<fN; q++){
        xi(0, 0) -= 1. / m_iter(0 + q * cov_dim, 0)*sin(m_iter(1 + q * cov_dim, 0))*cos(m_iter(2 + q * cov_dim, 0));
        xi(1, 0) -= 1. / m_iter(0 + q * cov_dim, 0)*sin(m_iter(1 + q * cov_dim, 0))*sin(m_iter(2 + q * cov_dim, 0));
        xi(2, 0) -= 1. / m_iter(0 + q * cov_dim, 0)*cos(m_iter(1 + q * cov_dim, 0));
    }

    return xi;
}

//Constraint equations
TMatrixD HKinFitter::f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter)
{
    TMatrixD d;

    if (fVtxConstraint)
    {
        d.ResizeTo(1, 1);
        TVector3 base_1, base_2, dir_1, dir_2;

        base_1.SetXYZ(
            m_iter(3 + 0 * cov_dim, 0) *
                std::cos(m_iter(2 + 0 * cov_dim, 0) + TMath::PiOver2()),
            m_iter(3 + 0 * cov_dim, 0) *
                std::sin(m_iter(2 + 0 * cov_dim, 0) + TMath::PiOver2()),
            m_iter(4 + 0 * cov_dim, 0));
        base_2.SetXYZ(
            m_iter(3 + 1 * cov_dim, 0) *
                std::cos(m_iter(2 + 1 * cov_dim, 0) + TMath::PiOver2()),
            m_iter(3 + 1 * cov_dim, 0) *
                std::sin(m_iter(2 + 1 * cov_dim, 0) + TMath::PiOver2()),
            m_iter(4 + 1 * cov_dim, 0));

        dir_1.SetXYZ(std::sin(m_iter(1 + 0 * cov_dim, 0)) *
                         std::cos(m_iter(2 + 0 * cov_dim, 0)),
                     std::sin(m_iter(1 + 0 * cov_dim, 0)) *
                         std::sin(m_iter(2 + 0 * cov_dim, 0)),
                     std::cos(m_iter(1 + 0 * cov_dim, 0)));
        dir_2.SetXYZ(std::sin(m_iter(1 + 1 * cov_dim, 0)) *
                         std::cos(m_iter(2 + 1 * cov_dim, 0)),
                     std::sin(m_iter(1 + 1 * cov_dim, 0)) *
                         std::sin(m_iter(2 + 1 * cov_dim, 0)),
                     std::cos(m_iter(1 + 1 * cov_dim, 0)));

        d(0, 0) = std::fabs((dir_1.Cross(dir_2)).Dot((base_1 - base_2)));
    }

    if (f3Constraint)
    {
        //mother
        d.ResizeTo(4, 1);
        d(0, 0) = -1. / xi_iter(0, 0) * sin(m_iter(0 + fN * cov_dim, 0)) * cos(m_iter(1 + fN * cov_dim, 0));
        d(1, 0) = -1. / xi_iter(0, 0) * sin(m_iter(0 + fN * cov_dim, 0)) * sin(m_iter(1 + fN * cov_dim, 0));
        d(2, 0) = -1. / xi_iter(0, 0) * cos(m_iter(0 + fN * cov_dim, 0));
        d(3, 0) = -sqrt(pow((1. / xi_iter(0, 0)), 2) + pow(fM[fN], 2));

        //daughters
        for (int q = 0; q < fN; q++)
        {
            d(0, 0) += 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            d(1, 0) += 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            d(2, 0) += 1. / m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0));
            d(3, 0) += sqrt(pow((1. / m_iter(0 + q * cov_dim, 0)), 2) + pow(fM[q], 2));
        }
    }

    if (f4Constraint)
    {
        d.ResizeTo(4, 1);

        //initial system
        d(0, 0) = -fInit.Px();
        d(1, 0) = -fInit.Py();
        d(2, 0) = -fInit.Pz();
        d(3, 0) = -fInit.E();

        //daughters
        for (int q = 0; q < fN; q++)
        {
            d(0, 0) += 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            d(1, 0) += 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            d(2, 0) += 1. / m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0));
            d(3, 0) += sqrt(pow((1. / m_iter(0 + q * cov_dim, 0)), 2) + pow(fM[q], 2));
        }
    }

    if (fMomConstraint)
    {
	
	d.ResizeTo(4, 1);
        d(0, 0) = -fInit.Px() + xi_iter(0,0);
        d(1, 0) = -fInit.Py() + xi_iter(1,0);
        d(2, 0) = -fInit.Pz() + xi_iter(2,0);
        d(3, 0) = -fInit.E() + sqrt(pow(xi_iter(0,0),2)+pow(xi_iter(1,0),2)+pow(xi_iter(2,0),2)+pow(fM[fN],2));
	
        for(int q=0; q<fN; q++){
            d(0,0) += 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            d(1,0) += 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            d(2,0) += 1. / m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0));
            d(3,0) += sqrt(pow((1. / m_iter(0 + q * cov_dim, 0)), 2) + pow(fM[q], 2));
        }
    }

    return d;
}

//Jacobian (derivative of constraint equations with respect to measured variables)
TMatrixD HKinFitter::Feta_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter)
{
    TMatrixD H;

    if (fVtxConstraint)
    {
        H.ResizeTo(1, fyDim);
        H.Zero();

        H(0, 0) = 0.;
        H(0, 5) = 0.;

        H(0, 1) = m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                      std::cos(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) +
                  m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                      std::sin(m_iter(2, 0)) -
                  m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                      std::cos(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) -
                  m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                      std::sin(m_iter(1, 0)) -
                  m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                      std::cos(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) -
                  m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::sin(m_iter(1, 0)) +
                  m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                      std::cos(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) +
                  m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::sin(m_iter(1, 0)) +
                  (m_iter(4, 0) - m_iter(9, 0)) *
                      (std::cos(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) -
                       std::cos(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)));

        H(0, 6) = -1 * m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) -
                  m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                      std::cos(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) + //
                  m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) +
                  m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                      std::cos(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) +
                  m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) +
                  m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                      std::cos(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) -
                  m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) -
                  m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                      std::cos(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) +
                  (m_iter(4, 0) - m_iter(9, 0)) *
                      (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::cos(m_iter(6, 0)) * std::sin(m_iter(7, 0)) -
                       std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                           std::cos(m_iter(6, 0)) * std::cos(m_iter(7, 0)));

        H(0, 2) = m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) -
                  m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) +
                  m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) -
                  m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) +
                  m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) -
                  m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) +
                  m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) - //
                  m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) +
                  (m_iter(9, 0) - m_iter(4, 0)) *
                      (std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) -
                       std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)));

        H(0, 7) = -1 * m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) +
                  m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) +
                  m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) -
                  m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) -
                  m_iter(3, 0) * std::sin(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) +
                  m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) +
                  m_iter(8, 0) * std::sin(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) -
                  m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) +
                  (m_iter(4, 0) - m_iter(9, 0)) *
                      (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) -
                       std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)));

        H(0, 3) = std::cos(m_iter(2, 0) + pi2) *
                      (std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                           std::cos(m_iter(6, 0)) -
                       std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                           std::cos(m_iter(1, 0))) -
                  std::sin(m_iter(2, 0) + pi2) *
                      (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::cos(m_iter(6, 0)) -
                       std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                           std::cos(m_iter(1, 0)));

        H(0, 8) = -1 * std::cos(m_iter(7, 0) + pi2) *
                      (std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                           std::cos(m_iter(6, 0)) -
                       std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                           std::cos(m_iter(1, 0))) +
                  std::sin(m_iter(7, 0) + pi2) *
                      (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::cos(m_iter(6, 0)) -
                       std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                           std::cos(m_iter(1, 0)));

        H(0, 4) = std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) -
                  std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0));

        H(0, 9) = -1 * std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) +
                  std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0));
    }

    if (f3Constraint)
    {
        H.ResizeTo(4, fyDim);
        H.Zero();

        //Daughter variables
        for (int q = 0; q < fN; q++)
        {
            //d(1/p)
            H(0, q * cov_dim) = -1. / pow(m_iter(0 + q * cov_dim, 0), 2) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            H(1, q * cov_dim) = -1. / pow(m_iter(0 + q * cov_dim, 0), 2) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            H(2, q * cov_dim) = -1. / pow(m_iter(0 + q * cov_dim, 0), 2) * cos(m_iter(1 + q * cov_dim, 0));
            H(3, q * cov_dim) = -1. / pow(m_iter(0 + q * cov_dim, 0), 3) * 1. / sqrt(pow(1. / (m_iter(0 + q * cov_dim, 0)), 2) + pow(fM[q], 2));

            //dtht
            H(0, 1 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            H(1, 1 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            H(2, 1 + q * cov_dim) = -1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0));

            //dphi
            H(0, 2 + q * cov_dim) = -1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            H(1, 2 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
        }

        //Mother variables
        //dtht
        H(0, fN * cov_dim) = -1. / xi_iter(0, 0) * cos(m_iter(0 + fN * cov_dim, 0)) * cos(m_iter(1 + fN * cov_dim, 0));
        H(1, fN * cov_dim) = -1. / xi_iter(0, 0) * cos(m_iter(0 + fN * cov_dim, 0)) * sin(m_iter(1 + fN * cov_dim, 0));
        H(2, fN * cov_dim) = 1. / xi_iter(0, 0) * sin(m_iter(0 + fN * cov_dim, 0));
        //dphi
        H(0, 1 + fN * cov_dim) = 1. / xi_iter(0, 0) * sin(m_iter(0 + fN * cov_dim, 0)) * sin(m_iter(1 + fN * cov_dim, 0));
        H(1, 1 + fN * cov_dim) = -1. / xi_iter(0, 0) * sin(m_iter(0 + fN * cov_dim, 0)) * cos(m_iter(1 + fN * cov_dim, 0));
    }

    if (f4Constraint || fMomConstraint)
    {
        //H.ResizeTo(4, fN * cov_dim);
        H.ResizeTo(4, fyDim);
        H.Zero();

        for(int q=0; q<fN; q++){

            //d(1/p)
            H(0, q * cov_dim) = -1. / pow(m_iter(0 + q * cov_dim, 0),2) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            H(1, q * cov_dim) = -1. / pow(m_iter(0 + q * cov_dim, 0),2) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            H(2, q * cov_dim) = -1. / pow(m_iter(0 + q * cov_dim, 0),2) * cos(m_iter(1 + q * cov_dim, 0));
            H(3, q * cov_dim) = -1. / pow(m_iter(0 + q * cov_dim, 0),3) * 1. / sqrt(pow(1. / (m_iter(0 + q * cov_dim, 0)),2) + pow(fM[q], 2));
            
            //dtht
            H(0, 1 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            H(1, 1 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            H(2, 1 + q * cov_dim) = -1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0));

            //dphi
            H(0, 2 + q * cov_dim) = -1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            H(1, 2 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
        }
    }

    return H;
}

//Jacobian (derivative of constraint equations with respect to unmeasured variables)
TMatrixD HKinFitter::Fxi_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter)
{
    TMatrixD H;

    if (f3Constraint)
    {
        H.ResizeTo(4, 1);
        H.Zero();

        //d(1/p)
        H(0, 0) = 1. / pow(xi_iter(0, 0), 2) * sin(m_iter(0 + fN * cov_dim, 0)) * cos(m_iter(1 + fN * cov_dim, 0));
        H(1, 0) = 1. / pow(xi_iter(0, 0), 2) * sin(m_iter(0 + fN * cov_dim, 0)) * sin(m_iter(1 + fN * cov_dim, 0));
        H(2, 0) = 1. / pow(xi_iter(0, 0), 2) * cos(m_iter(0 + fN * cov_dim, 0));
        H(3, 0) = 1. / pow(xi_iter(0, 0), 3) * 1. / sqrt(pow(1. / xi_iter(0, 0), 2) + pow(fM[fN], 2));
    }

    if (fMomConstraint)
    {
        H.ResizeTo(4, 3);
        H.Zero();

        H(0, 0) = 1;
        H(1, 1) = 1;
        H(2, 2) = 1;

        H(3, 0) = xi_iter(0,0) / sqrt( pow(xi_iter(0,0),2) + pow(xi_iter(1,0),2) + pow(xi_iter(2,0),2) + pow(fM[fN],2) );
        H(3, 1) = xi_iter(1,0) / sqrt( pow(xi_iter(0,0),2) + pow(xi_iter(1,0),2) + pow(xi_iter(2,0),2) + pow(fM[fN],2) );
        H(3, 2) = xi_iter(2,0) / sqrt( pow(xi_iter(0,0),2) + pow(xi_iter(1,0),2) + pow(xi_iter(2,0),2) + pow(fM[fN],2) );
    }

    return H;
}

bool HKinFitter::fit()
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HKinFitter::fit() -----------" << std::endl;
        std::cout << "Vertex constraint set: " << fVtxConstraint << std::endl;
        std::cout << "3C set: " << f3Constraint << std::endl;
        std::cout << "4C set: " << f4Constraint << std::endl;
        std::cout << "Momentum constraint set: " << fMomConstraint << std::endl;
        std::cout << "" << std::endl;
    }

    double lr = fLearningRate;
    TMatrixD alpha0(fyDim, 1), alpha(fyDim, 1);
    TMatrixD xi0(x), xi(x), neu_xi(x);
    TMatrixD A0(y), V0(V);
    alpha0 = y;
    alpha = alpha0;

    // Calculating the inverse of the original covariance matrix that is not changed in the iterations
    TMatrixD V0_inv(V);
    V0_inv.Invert();
/*
    xi0.Zero();
    xi.Zero();
    neu_xi.Zero();
*/
    if (f4Constraint == false && f3Constraint == false && fVtxConstraint == false && fMomConstraint == false)
    {
        Error("fit()", Form("No constraint is chosen, please add constraint"));
        abort();
    }

    if (f3Constraint)
    {
        xi0(0, 0) = 1 / fMother.P();
        xi = xi0;
    }

    if (fMomConstraint)
    {
        xi0 = calcMissingMom(alpha0);
        xi = xi0;
    }

    double chi2 = 1e6;
    if (fVerbose > 1)
    {
        cout << " calc Feta" << endl;
    }
    TMatrixD D = Feta_eval(alpha, xi);
    TMatrixD DT(D.GetNcols(), D.GetNrows());
    if (fVerbose > 1)
    {
        cout << " calc f" << endl;
    }
    TMatrixD d = f_eval(alpha, xi);
    TMatrixD D_xi(d.GetNrows(), 1), DT_xi(1, d.GetNrows());
    if (d.GetNrows()-fNdf>0) D_xi.ResizeTo(d.GetNrows(), d.GetNrows()-fNdf); DT_xi.ResizeTo(d.GetNrows()-fNdf, d.GetNrows()); //check dimension if other fitters are added
    D_xi.Zero();
    DT_xi.Zero();
    if (f3Constraint || fMomConstraint){
        if (fVerbose > 1)
        {
            cout << " calc Fxi" << endl;
        }
        D_xi = Fxi_eval(alpha, xi);
    }
    TMatrixD VD(D.GetNrows(), D.GetNrows());
    VD.Zero();
    TMatrixD VDD(D_xi.GetNcols(), D_xi.GetNcols());
    VDD.Zero();

    for (int q = 0; q < fNumIterations; q++)
    {
        TMatrixD delta_alpha = alpha0 - alpha;
        //calc r
        if (fVerbose > 1)
        {
            cout << " calc r" << endl;
        }
        TMatrixD r = d + D * delta_alpha;
        DT.Transpose(D);
        //calc S
        if (fVerbose > 1)
        {
            cout << " calc S" << endl;
        }
        VD = D * V0 * DT;
        VD.Invert();
        if (f3Constraint || fMomConstraint)
        {
            DT_xi.Transpose(D_xi);
            if (fVerbose > 1)
            {
                cout << " calc Sxi" << endl;
            }
            VDD = DT_xi * VD * D_xi;
            VDD.Invert();
        }

        //calculate values for next iteration
        TMatrixD lambda(d); //Lagrange multiplier
        lambda.Zero();
        if (f3Constraint || fMomConstraint)
        {
            if (fVerbose > 1)
            {
                cout << " calc neuxi" << endl;
            }
            neu_xi = xi - lr * VDD * DT_xi * VD * r;
        }
        TMatrixD delta_xi = neu_xi - xi;
        if (fVerbose > 1)
        {
            cout << " calc lambda" << endl;
        }
        lambda = VD * (r + D_xi * delta_xi);
        TMatrixD lambdaT(lambda.GetNcols(), lambda.GetNrows());
        lambdaT.Transpose(lambda);
        TMatrixD neu_alpha(fyDim, 1);
        if (fVerbose > 1)
        {
            cout << " calc neueta" << endl;
        }
        neu_alpha = alpha0 - lr * V0 * DT * lambda;

        //Calculate new chi2
        TMatrixD chisqrd(1, 1);
        TMatrixD delta_alphaT(delta_alpha.GetNcols(), delta_alpha.GetNrows());
        delta_alphaT.Transpose(delta_alpha);
        TMatrixD two(1, 1);
        two(0, 0) = 2;
        if (fVerbose > 1)
        {
            cout << " calc chi2" << endl;
        }
        chisqrd = delta_alphaT * V0_inv * delta_alpha + two * lambdaT * d;
        //chisqrd = delta_alphaT * V0_inv * delta_alpha + two * lambdaT * f_eval(neu_alpha, neu_xi);
        
        //for checking convergence
        fIteration = q;
        if (fabs(chi2 - chisqrd(0, 0)) < fConvergenceCriteria)
        {
            fConverged = true;
            chi2 = chisqrd(0, 0);
            alpha = neu_alpha;
            if (f3Constraint || fMomConstraint)
                xi = neu_xi;
            break;
        }

        if (fVerbose > 0)
        {
            std::cout << "Iteration: " << q << std::endl;
            std::cout << "Printing d: " << std::endl;
            d.Print();
        }

        chi2 = chisqrd(0, 0);
        alpha = neu_alpha;
        if (f3Constraint || fMomConstraint)
            xi = neu_xi;
        D = Feta_eval(alpha, xi);
        if (f3Constraint || fMomConstraint)
            D_xi = Fxi_eval(alpha, xi);
        d = f_eval(alpha, xi);

    }

    if (fVerbose > 0)
    {
        cout << " exit iterations" << endl;
    }
    y = alpha;
    x = xi;
    if (fMomConstraint)
        fMissDaughter.SetXYZM(x(0,0), x(1,0), x(2,0), fM[fN]);
    fChi2 = chi2;
    fProb = TMath::Prob(chi2, fNdf);

    //Update covariance
    if (fVerbose > 1)
    {
        cout << " calc Vneu" << endl;
    }
    if (f3Constraint || fMomConstraint)
    {
        TMatrixD matrix = DT * VD * D_xi;
        TMatrixD matrixT(matrix.GetNcols(), matrix.GetNrows());
        matrixT.Transpose(matrix);
        TMatrixD invertedMatrix = DT_xi * VD * D_xi;
        invertedMatrix.Invert();
        V = V0 - lr * V0 * (DT * VD * D - (matrix * invertedMatrix * matrixT)) * V0;
        Vx = invertedMatrix;
    }
    if (fVtxConstraint || f4Constraint)
        V = V0 - lr * V0 * DT * VD * D * V0;

    // -----------------------------------------
    // Pull
    // -----------------------------------------
    fPull.ResizeTo(fyDim, fyDim);
    for (int b = 0; b < (fyDim); b++)
        fPull(b, b) = -10000;

    if (true)
    {
        for (int b = 0; b < (fyDim); b++)
        {
            double num = A0(b, 0) - alpha(b, 0);
            double dem = V0(b, b) - V(b, b);
            if (dem > 0)
            {
                fPull(b, b) = num / std::sqrt(dem);
            }
        }
    }

    updateDaughters();
    if (f3Constraint)
        updateMother();

    return fConverged; // for number of iterations greater than 1
    //return true; // for number of iterations equal to 1
}

HRefitCand HKinFitter::getDaughter(int val)
{
    return fCands[val];
}

HRefitCand HKinFitter::getMother()
{
    return fMother;
}

TLorentzVector HKinFitter::getMissingDaughter()
{
    return fMissDaughter;
}

void HKinFitter::updateDaughters()
{
    for (int val = 0; val < fN; ++val)
    {
        HRefitCand &cand = fCands[val];
        double Px = (1. / y(0 + val * cov_dim, 0)) *
                    std::sin(y(1 + val * cov_dim, 0)) *
                    std::cos(y(2 + val * cov_dim, 0));
        double Py = (1. / y(0 + val * cov_dim, 0)) *
                    std::sin(y(1 + val * cov_dim, 0)) *
                    std::sin(y(2 + val * cov_dim, 0));
        double Pz =
            (1. / y(0 + val * cov_dim, 0)) * std::cos(y(1 + val * cov_dim, 0));
        double M = fM[val];
        cand.SetXYZM(Px, Py, Pz, M);
        cand.setR(y(3 + val * cov_dim, 0));
        cand.setZ(y(4 + val * cov_dim, 0));

        // ---------------------------------------------------------------------------
        // set covariance
        // ---------------------------------------------------------------------------
        TMatrixD cov(5, 5);
        cov(0, 0) = V(0 + val * cov_dim, 0 + val * cov_dim);
        cov(1, 1) = V(1 + val * cov_dim, 1 + val * cov_dim);
        cov(2, 2) = V(2 + val * cov_dim, 2 + val * cov_dim);
        cov(3, 3) = V(3 + val * cov_dim, 3 + val * cov_dim);
        cov(4, 4) = V(4 + val * cov_dim, 4 + val * cov_dim);
        cand.setCovariance(cov);
        // ---------------------------------------------------------------------------
    }
}

void HKinFitter::updateMother()
{
    HRefitCand &mother = fMother;
    double Px = (1. / x(0, 0)) *
                std::sin(y(0 + fN * cov_dim, 0)) *
                std::cos(y(1 + fN * cov_dim, 0));
    double Py = (1. / x(0, 0)) *
                std::sin(y(0 + fN * cov_dim, 0)) *
                std::sin(y(1 + fN * cov_dim, 0));
    double Pz =
        (1. / x(0, 0)) * std::cos(y(0 + fN * cov_dim, 0));
    double M = fM[fN];
    mother.SetXYZM(Px, Py, Pz, M);
    mother.setR(y(2 + fN * cov_dim, 0));
    mother.setZ(y(3 + fN * cov_dim, 0));

    // ---------------------------------------------------------------------------
    // set covariance
    // ---------------------------------------------------------------------------
    TMatrixD cov(5, 5);
    cov(0, 0) = Vx(0, 0);
    cov(1, 1) = V(0 + fN * cov_dim, 0 + fN * cov_dim);
    cov(2, 2) = V(1 + fN * cov_dim, 1 + fN * cov_dim);
    cov(3, 3) = V(2 + fN * cov_dim, 2 + fN * cov_dim);
    cov(4, 4) = V(3 + fN * cov_dim, 3 + fN * cov_dim);
    mother.setCovariance(cov);
    // ---------------------------------------------------------------------------
}

void HKinFitter::update()
{
    for (int val = 0; val < fN; ++val)
    {
        HRefitCand &cand = fCands[val];
        cand.update();
    }
}
