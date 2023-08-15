#include "KinFitter.h"

const size_t cov_dim = 5;

// general KinFitter constructor
KinFitter::KinFitter(const std::vector<KFitParticle> &cands) : fCands(cands)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::KinFitter() -----------------" << std::endl;
    }

    // fN is the number of daughters e.g. (L->ppi-) n=2
    fN = fCands.size();
    fyDim = fN * cov_dim; // Dimension of full covariance matrix (number of measured variables x cov_dim)

    y.ResizeTo(fyDim, 1);
    x.ResizeTo(1, 1);
    V.ResizeTo(fyDim, fyDim);
    Vx.ResizeTo(1, 1);

    y.Zero();
    x.Zero();
    V.Zero();
    Vx.Zero();

    // set 'y=alpha' measurements
    // and the covariance
    for (int ix = 0; ix < fN; ix++)
    {
        KFitParticle cand = fCands[ix];
        y(0 + ix * cov_dim, 0) = 1. / cand.P();
        y(1 + ix * cov_dim, 0) = cand.Theta();
        y(2 + ix * cov_dim, 0) = cand.Phi();
        y(3 + ix * cov_dim, 0) = cand.getR();
        y(4 + ix * cov_dim, 0) = cand.getZ();
        fM.push_back(cand.M());

        TMatrixD covariance = cand.getCovariance();
        V(0 + ix * cov_dim, 0 + ix * cov_dim) = covariance(0, 0);
        V(1 + ix * cov_dim, 1 + ix * cov_dim) = covariance(1, 1);
        V(2 + ix * cov_dim, 2 + ix * cov_dim) = covariance(2, 2);
        V(3 + ix * cov_dim, 3 + ix * cov_dim) = covariance(3, 3);
        V(4 + ix * cov_dim, 4 + ix * cov_dim) = covariance(4, 4);
    }
}

void KinFitter::addMassConstraint(double mass)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::addMassConstraint() -----------------" << std::endl;
    }

    fMass = mass;
    if (!fMassConstraint)
        fNdf += 1;
    fMassConstraint = true;
}

void KinFitter::addMissingMassConstraint(TLorentzVector lv, double mass)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::addMissingMassConstraint() -----------------" << std::endl;
    }

    fMass = mass;
    fInit = lv;
    if (!fMissingMassConstraint)
        fNdf += 1;
    fMissingMassConstraint = true;
}

void KinFitter::addMassVtxConstraint(double mass)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::addMassVtxConstraint() -----------------" << std::endl;
    }

    if (!fMassVtxConstraint)
        fNdf += 2;
    fMass = mass;
    fMassVtxConstraint = true;
}

void KinFitter::add4Constraint(TLorentzVector lv)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::add4Constraint() -----------------" << std::endl;
    }

    fInit = lv;

    if (!f4Constraint)
        fNdf += 4;
    f4Constraint = true;
}

void KinFitter::add3Constraint(KFitParticle mother)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::add3Constraint() -----------------" << std::endl;
    }

    fyDim = (fN + 1) * cov_dim - 1; // Dimension of full covariance matrix (number of measured variables x cov_dim). Mother momentum is not measured

    y.ResizeTo(fyDim, 1);
    V.ResizeTo(fyDim, fyDim);

    y.Zero();
    V.Zero();

    fM.clear();

    // Set y to measurements and the covariance, set mass
    for (int ix = 0; ix < fN; ix++) // for daughters
    {
        KFitParticle cand = fCands[ix];

        y(0 + ix * cov_dim, 0) = 1. / cand.P();
        y(1 + ix * cov_dim, 0) = cand.Theta();
        y(2 + ix * cov_dim, 0) = cand.Phi();
        y(3 + ix * cov_dim, 0) = cand.getR();
        y(4 + ix * cov_dim, 0) = cand.getZ();
        fM.push_back(cand.M());

        TMatrixD covariance = cand.getCovariance();
        V(0 + ix * cov_dim, 0 + ix * cov_dim) = covariance(0, 0);
        V(1 + ix * cov_dim, 1 + ix * cov_dim) = covariance(1, 1);
        V(2 + ix * cov_dim, 2 + ix * cov_dim) = covariance(2, 2);
        V(3 + ix * cov_dim, 3 + ix * cov_dim) = covariance(3, 3);
        V(4 + ix * cov_dim, 4 + ix * cov_dim) = covariance(4, 4);
    }

    // Set the covariance matrix for the mother particle
    TMatrixD test = fMother.getCovariance();
    test.ResizeTo(5, 5);
    fMother.setCovariance(test);
    fMother = mother;

    y(fN * cov_dim, 0) = fMother.Theta();
    y(1 + fN * cov_dim, 0) = fMother.Phi();
    y(2 + fN * cov_dim, 0) = fMother.getR();
    y(3 + fN * cov_dim, 0) = fMother.getZ();
    fM.push_back(fMother.M());

    TMatrixD covariance = fMother.getCovariance();

    if (fVerbose > 0)
    {
        std::cout << "Decaying Particle Covariance: " << covariance(0, 0) << " " << covariance(1, 1) << " " << covariance(2, 2) << " " << covariance(3, 3) << " " << covariance(4, 4) << std::endl;
    }

    V(0 + fN * cov_dim, 0 + fN * cov_dim) = covariance(1, 1);
    V(1 + fN * cov_dim, 1 + fN * cov_dim) = covariance(2, 2);
    V(2 + fN * cov_dim, 2 + fN * cov_dim) = covariance(3, 3);
    V(3 + fN * cov_dim, 3 + fN * cov_dim) = covariance(4, 4);

    V(0 + fN * cov_dim, 1 + fN * cov_dim) = covariance(1, 2);
    V(1 + fN * cov_dim, 0 + fN * cov_dim) = covariance(2, 1);

    if (fVerbose > 0)
    {
        std::cout << "Fitter" << std::endl;
        std::cout << "Diag elements: " << V(0 + fN * cov_dim, 0 + fN * cov_dim) << " " << V(1 + fN * cov_dim, 1 + fN * cov_dim) << std::endl;
        std::cout << "Off diag elements: " << V(0 + fN * cov_dim, 1 + fN * cov_dim) << " " << V(1 + fN * cov_dim, 0 + fN * cov_dim) << std::endl;
    }

    if (!f3Constraint)
        fNdf += 3;
    f3Constraint = true;
}

void KinFitter::addPvaConstraint(KFitParticle mother)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::addPvaConstraint() -----------------" << std::endl;
    }

    fyDim = (fN + 1) * cov_dim - 1; // Dimension of full covariance matrix (number of measured variables x cov_dim). Mother momentum is not measured

    y.ResizeTo(fyDim, 1);
    V.ResizeTo(fyDim, fyDim);

    y.Zero();
    V.Zero();

    fM.clear();

    // Set y to measurements and the covariance, set mass
    for (int ix = 0; ix < fN; ix++) // for daughters
    {
        KFitParticle cand = fCands[ix];

        y(0 + ix * cov_dim, 0) = 1. / cand.P();
        y(1 + ix * cov_dim, 0) = cand.Theta();
        y(2 + ix * cov_dim, 0) = cand.Phi();
        y(3 + ix * cov_dim, 0) = cand.getR();
        y(4 + ix * cov_dim, 0) = cand.getZ();
        fM.push_back(cand.M());

        TMatrixD covariance = cand.getCovariance();
        V(0 + ix * cov_dim, 0 + ix * cov_dim) = covariance(0, 0);
        V(1 + ix * cov_dim, 1 + ix * cov_dim) = covariance(1, 1);
        V(2 + ix * cov_dim, 2 + ix * cov_dim) = covariance(2, 2);
        V(3 + ix * cov_dim, 3 + ix * cov_dim) = covariance(3, 3);
        V(4 + ix * cov_dim, 4 + ix * cov_dim) = covariance(4, 4);
    }

    // Set the covariance matrix for the mother particle
    TMatrixD test = fMother.getCovariance();
    test.ResizeTo(5, 5);
    fMother.setCovariance(test);
    fMother = mother;

    y(fN * cov_dim, 0) = fMother.Theta();
    y(1 + fN * cov_dim, 0) = fMother.Phi();
    y(2 + fN * cov_dim, 0) = fMother.getR();
    y(3 + fN * cov_dim, 0) = fMother.getZ();
    fM.push_back(fMother.M());

    TMatrixD covariance = fMother.getCovariance();

    if (fVerbose > 0)
    {
        std::cout << "Decaying Particle Covariance: " << covariance(0, 0) << " " << covariance(1, 1) << " " << covariance(2, 2) << " " << covariance(3, 3) << " " << covariance(4, 4) << std::endl;
    }

    V(0 + fN * cov_dim, 0 + fN * cov_dim) = covariance(1, 1);
    V(1 + fN * cov_dim, 1 + fN * cov_dim) = covariance(2, 2);
    V(2 + fN * cov_dim, 2 + fN * cov_dim) = covariance(3, 3);
    V(3 + fN * cov_dim, 3 + fN * cov_dim) = covariance(4, 4);

    V(0 + fN * cov_dim, 1 + fN * cov_dim) = covariance(1, 2);
    V(1 + fN * cov_dim, 0 + fN * cov_dim) = covariance(2, 1);

    if (fVerbose > 0)
    {
        std::cout << "Fitter" << std::endl;
        std::cout << "Diag elements: " << V(0 + fN * cov_dim, 0 + fN * cov_dim) << " " << V(1 + fN * cov_dim, 1 + fN * cov_dim) << std::endl;
        std::cout << "Off diag elements: " << V(0 + fN * cov_dim, 1 + fN * cov_dim) << " " << V(1 + fN * cov_dim, 0 + fN * cov_dim) << std::endl;
    }

    if (!fPvaConstraint)
        fNdf += 2;
    fPvaConstraint = true;
}

void KinFitter::addVertexConstraint()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::addVertexConstraint() -----------------" << std::endl;
    }

    if (!fVtxConstraint)
        fNdf += 1;
    fVtxConstraint = true;
}

// For fit with missing particle
void KinFitter::addMissingParticleConstraint(TLorentzVector lv, double mass)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::addMissingParticleConstraint() -----------------" << std::endl;
    }

    fInit = lv;
    fMass = mass;

    x.ResizeTo(3, 1);
    Vx.ResizeTo(3, 3);
    x.Zero();
    Vx.Zero();

    if (!fMissingParticleConstraint)
    {
        fM.push_back(fMass);
        fNdf += 1;
    }

    fMissingParticleConstraint = true;
}

TMatrixD KinFitter::calcMissingMom(const TMatrixD &m_iter)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::calcMissingMom() -----------------" << std::endl;
    }

    TMatrix xi(3, 1);

    xi(0, 0) = fInit.Px();
    xi(1, 0) = fInit.Py();
    xi(2, 0) = fInit.Pz();

    for (int q = 0; q < fN; q++)
    {
        xi(0, 0) -= 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
        xi(1, 0) -= 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
        xi(2, 0) -= 1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0));
    }

    return xi;
}

// Constraint equations
TMatrixD KinFitter::f_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::f_eval() -----------------" << std::endl;
    }

    TMatrixD d;

    if (fMassConstraint)
    {
        d.ResizeTo(1, 1);
        double Px = 0., Py = 0., Pz = 0., E = 0.;

        for (int q = 0; q < fN; q++)
        {
            E += std::sqrt((1. / m_iter(0 + q * cov_dim, 0)) *
                               (1. / m_iter(0 + q * cov_dim, 0)) +
                           fM[q] * fM[q]);
            Px += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::sin(m_iter(1 + q * cov_dim, 0)) *
                  std::cos(m_iter(2 + q * cov_dim, 0));
            Py += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::sin(m_iter(1 + q * cov_dim, 0)) *
                  std::sin(m_iter(2 + q * cov_dim, 0));
            Pz += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::cos(m_iter(1 + q * cov_dim, 0));
        }

        d(0, 0) = std::pow(E, 2) - std::pow(Px, 2) - std::pow(Py, 2) -
                  std::pow(Pz, 2) - fMass * fMass;
    }

    // Invariant mass + vertex constraint
    if (fMassVtxConstraint)
    {
        d.ResizeTo(2, 1);
        double Px = 0., Py = 0., Pz = 0., E = 0.;

        for (int q = 0; q < fN; q++)
        {
            E += std::sqrt((1. / m_iter(0 + q * cov_dim, 0)) *
                               (1. / m_iter(0 + q * cov_dim, 0)) +
                           fM[q] * fM[q]);
            Px += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::sin(m_iter(1 + q * cov_dim, 0)) *
                  std::cos(m_iter(2 + q * cov_dim, 0));
            Py += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::sin(m_iter(1 + q * cov_dim, 0)) *
                  std::sin(m_iter(2 + q * cov_dim, 0));
            Pz += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::cos(m_iter(1 + q * cov_dim, 0));
        }

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

        d(0, 0) = std::pow(E, 2) - std::pow(Px, 2) - std::pow(Py, 2) -
                  std::pow(Pz, 2) - fMass * fMass;
        d(1, 0) = (dir_1.Cross(dir_2)).Dot(base_2 - base_1);
    }

    // Missing mass constraint
    if (fMissingMassConstraint)
    {
        d.ResizeTo(1, 1);

        // initial system
        double Px = fInit.Px();
        double Py = fInit.Py();
        double Pz = fInit.Pz();
        double E = fInit.E();

        for (int q = 0; q < fN; q++)
        {
            Px += 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            Py += 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            Pz += 1. / m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0));
            E += sqrt(pow((1. / m_iter(0 + q * cov_dim, 0)), 2) + pow(fM[q], 2));
        }
        d(0, 0) = std::pow(E, 2) - std::pow(Px, 2) - std::pow(Py, 2) - std::pow(Pz, 2) - fMass * fMass;
    }

    // Vertex constraint
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

        d(0, 0) = (dir_1.Cross(dir_2)).Dot(base_1 - base_2);
    }

    // For 4-momentum fit in vertex with fixed mass, momentum unmeasured
    if (f3Constraint)
    {
        // Mother
        d.ResizeTo(4, 1);
        d(0, 0) = -1. / xi_iter(0, 0) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::cos(m_iter(1 + fN * cov_dim, 0));
        d(1, 0) = -1. / xi_iter(0, 0) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::sin(m_iter(1 + fN * cov_dim, 0));
        d(2, 0) = -1. / xi_iter(0, 0) * std::cos(m_iter(0 + fN * cov_dim, 0));
        d(3, 0) = -sqrt(pow((1. / xi_iter(0, 0)), 2) + std::pow(fM[fN], 2));

        // Daughters
        for (int q = 0; q < fN; q++)
        {
            d(0, 0) += 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
            d(1, 0) += 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            d(2, 0) += 1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0));
            d(3, 0) += sqrt(pow((1. / m_iter(0 + q * cov_dim, 0)), 2) + std::pow(fM[q], 2));
        }
    }

    // For 3-momentum fit in vertex, momentum magnitude unmeasured
    if (fPvaConstraint)
    {
        // Mother
        d.ResizeTo(3, 1);
        d(0, 0) = -1. / xi_iter(0, 0) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::cos(m_iter(1 + fN * cov_dim, 0));
        d(1, 0) = -1. / xi_iter(0, 0) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::sin(m_iter(1 + fN * cov_dim, 0));
        d(2, 0) = -1. / xi_iter(0, 0) * std::cos(m_iter(0 + fN * cov_dim, 0));

        // Daughters
        for (int q = 0; q < fN; q++)
        {
            d(0, 0) += 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
            d(1, 0) += 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            d(2, 0) += 1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0));
        }
    }

    // 4C fit
    if (f4Constraint)
    {
        d.ResizeTo(4, 1);

        // Initial system
        d(0, 0) = -fInit.Px();
        d(1, 0) = -fInit.Py();
        d(2, 0) = -fInit.Pz();
        d(3, 0) = -fInit.E();

        // Daughters
        for (int q = 0; q < fN; q++)
        {
            d(0, 0) += 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
            d(1, 0) += 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            d(2, 0) += 1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0));
            d(3, 0) += sqrt(pow((1. / m_iter(0 + q * cov_dim, 0)), 2) + std::pow(fM[q], 2));
        }
    }

    // For missing particle fit, all 4-momenta constrained to that of initial system
    if (fMissingParticleConstraint)
    {
        d.ResizeTo(4, 1);

        d(0, 0) = -fInit.Px() + xi_iter(0, 0);
        d(1, 0) = -fInit.Py() + xi_iter(1, 0);
        d(2, 0) = -fInit.Pz() + xi_iter(2, 0);
        d(3, 0) = -fInit.E() + sqrt(pow(xi_iter(0, 0), 2) + pow(xi_iter(1, 0), 2) + pow(xi_iter(2, 0), 2) + pow(fM[fN], 2));

        for (int q = 0; q < fN; q++)
        {
            d(0, 0) += 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * cos(m_iter(2 + q * cov_dim, 0));
            d(1, 0) += 1. / m_iter(0 + q * cov_dim, 0) * sin(m_iter(1 + q * cov_dim, 0)) * sin(m_iter(2 + q * cov_dim, 0));
            d(2, 0) += 1. / m_iter(0 + q * cov_dim, 0) * cos(m_iter(1 + q * cov_dim, 0));
            d(3, 0) += sqrt(pow((1. / m_iter(0 + q * cov_dim, 0)), 2) + pow(fM[q], 2));
        }
    }

    return d;
}

// Jacobian (derivative of constraint equations with respect to measured variables)
TMatrixD KinFitter::Feta_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::Feta_eval() -----------------" << std::endl;
    }

    TMatrixD H;

    if (fMassVtxConstraint)
    {
        H.ResizeTo(2, fyDim);
        H.Zero();

        double Px = 0., Py = 0., Pz = 0., E = 0.;

        for (int q = 0; q < fN; q++)
        {
            E += std::sqrt((1. / m_iter(0 + q * cov_dim, 0)) *
                               (1. / m_iter(0 + q * cov_dim, 0)) +
                           fM[q] * fM[q]);
            Px += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::sin(m_iter(1 + q * cov_dim, 0)) *
                  std::cos(m_iter(2 + q * cov_dim, 0));
            Py += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::sin(m_iter(1 + q * cov_dim, 0)) *
                  std::sin(m_iter(2 + q * cov_dim, 0));
            Pz += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::cos(m_iter(1 + q * cov_dim, 0));
        }

        for (int q = 0; q < fN; q++)
        {
            double Pi = 1. / m_iter(0 + q * cov_dim, 0);
            double Ei = std::sqrt(Pi * Pi + fM[q] * fM[q]);
            H(0, 0 + q * cov_dim) =
                -2 * E * (std::pow(Pi, 3) / Ei) +
                2 * std::pow(Pi, 2) * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::cos(m_iter(2 + q * cov_dim, 0)) * Px +
                2 * std::pow(Pi, 2) * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::sin(m_iter(2 + q * cov_dim, 0)) * Py +
                2 * std::pow(Pi, 2) * std::cos(m_iter(1 + q * cov_dim, 0)) * Pz;

            H(0, 1 + q * cov_dim) =
                -2 * Pi * std::cos(m_iter(1 + q * cov_dim, 0)) *
                    std::cos(m_iter(2 + q * cov_dim, 0)) * Px -
                2 * Pi * std::cos(m_iter(1 + q * cov_dim, 0)) *
                    std::sin(m_iter(2 + q * cov_dim, 0)) * Py +
                2 * Pi * std::sin(m_iter(1 + q * cov_dim, 0)) * Pz;

            H(0, 2 + q * cov_dim) =
                2 * Pi * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::sin(m_iter(2 + q * cov_dim, 0)) * Px -
                2 * Pi * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::cos(m_iter(2 + q * cov_dim, 0)) * Py;
            H(0, 3 + q * cov_dim) = 0.;
            H(0, 4 + q * cov_dim) = 0.;
            H(1, 0 + q * cov_dim) = 0.;
        }

        H(1, 1) = m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                      std::cos(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::cos(m_iter(6, 0)) +
                  m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                      std::sin(m_iter(1, 0)) -
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

        H(1, 6) = -1 * m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
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

        H(1, 2) = m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
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
                  (m_iter(4, 0) - m_iter(9, 0)) *
                      (-std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) -
                       std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)));

        H(1, 7) = -1 * m_iter(3, 0) * std::cos(m_iter(2, 0) + pi2) *
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
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) -
                  m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) +
                  (m_iter(4, 0) - m_iter(9, 0)) *
                      (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) +
                       std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)));

        H(1, 3) = std::cos(m_iter(2, 0) + pi2) *
                      (std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                           std::cos(m_iter(6, 0)) -
                       std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                           std::cos(m_iter(1, 0))) -
                  std::sin(m_iter(2, 0) + pi2) *
                      (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::cos(m_iter(6, 0)) -
                       std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                           std::cos(m_iter(1, 0)));

        H(1, 8) = -1 * std::cos(m_iter(7, 0) + pi2) *
                      (std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                           std::cos(m_iter(6, 0)) -
                       std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                           std::cos(m_iter(1, 0))) +
                  std::sin(m_iter(7, 0) + pi2) *
                      (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::cos(m_iter(6, 0)) -
                       std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                           std::cos(m_iter(1, 0)));

        H(1, 4) = std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) -
                  std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0));

        H(1, 9) = -1 * std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) +
                  std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0));
    }

    // Missing mass constraint
    if (fMissingMassConstraint)
    {
        H.ResizeTo(1, fyDim);
        H.Zero();

        double Px = 0., Py = 0., Pz = 0., E = 0.;
        Px += fInit.Px();
        Py += fInit.Py();
        Pz += fInit.Pz();
        E += fInit.E();
        for (int q = 0; q < fN; q++)
        {
            E -= std::sqrt((1. / m_iter(0 + q * cov_dim, 0)) *
                               (1. / m_iter(0 + q * cov_dim, 0)) +
                           fM[q] * fM[q]);
            Px -= (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::sin(m_iter(1 + q * cov_dim, 0)) *
                  std::cos(m_iter(2 + q * cov_dim, 0));
            Py -= (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::sin(m_iter(1 + q * cov_dim, 0)) *
                  std::sin(m_iter(2 + q * cov_dim, 0));
            Pz -= (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::cos(m_iter(1 + q * cov_dim, 0));
        }

        for (int q = 0; q < fN; q++)
        {
            double Pi = 1. / m_iter(0 + q * cov_dim, 0);
            double Ei = std::sqrt(Pi * Pi + fM[q] * fM[q]);
            H(0, 0 + q * cov_dim) =
                2 * E * (std::pow(Pi, 3) / Ei) -
                2 * std::pow(Pi, 2) * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::cos(m_iter(2 + q * cov_dim, 0)) * Px -
                2 * std::pow(Pi, 2) * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::sin(m_iter(2 + q * cov_dim, 0)) * Py -
                2 * std::pow(Pi, 2) * std::cos(m_iter(1 + q * cov_dim, 0)) * Pz;

            H(0, 1 + q * cov_dim) =
                2 * Pi * std::cos(m_iter(1 + q * cov_dim, 0)) *
                    std::cos(m_iter(2 + q * cov_dim, 0)) * Px +
                2 * Pi * std::cos(m_iter(1 + q * cov_dim, 0)) *
                    std::sin(m_iter(2 + q * cov_dim, 0)) * Py -
                2 * Pi * std::sin(m_iter(1 + q * cov_dim, 0)) * Pz;

            H(0, 2 + q * cov_dim) =
                -2 * Pi * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::sin(m_iter(2 + q * cov_dim, 0)) * Px +
                2 * Pi * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::cos(m_iter(2 + q * cov_dim, 0)) * Py;
        }
    }

    // Vertex constraint
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
                      std::sin(m_iter(1, 0)) -
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
                  (m_iter(4, 0) - m_iter(9, 0)) *
                      (-std::sin(m_iter(1, 0)) * std::sin(m_iter(2, 0)) *
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
                      std::sin(m_iter(6, 0)) * std::sin(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) -
                  m_iter(8, 0) * std::cos(m_iter(7, 0) + pi2) *
                      std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) *
                      std::cos(m_iter(1, 0)) +
                  (m_iter(4, 0) - m_iter(9, 0)) *
                      (std::sin(m_iter(1, 0)) * std::cos(m_iter(2, 0)) *
                           std::sin(m_iter(6, 0)) * std::cos(m_iter(7, 0)) +
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

    // 3C constraint
    if (f3Constraint)
    {
        H.ResizeTo(4, fyDim);
        H.Zero();

        // Daughter variables
        for (int q = 0; q < fN; q++)
        {
            // d(1/p)
            H(0, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 2) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
            H(1, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 2) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            H(2, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 2) * std::cos(m_iter(1 + q * cov_dim, 0));
            H(3, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 3) * 1. / sqrt(pow(1. / (m_iter(0 + q * cov_dim, 0)), 2) + std::pow(fM[q], 2));

            // dtht
            H(0, 1 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
            H(1, 1 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            H(2, 1 + q * cov_dim) = -1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0));

            // dphi
            H(0, 2 + q * cov_dim) = -1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            H(1, 2 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
        }

        // Mother variables
        // dtht
        H(0, fN * cov_dim) = -1. / xi_iter(0, 0) * std::cos(m_iter(0 + fN * cov_dim, 0)) * std::cos(m_iter(1 + fN * cov_dim, 0));
        H(1, fN * cov_dim) = -1. / xi_iter(0, 0) * std::cos(m_iter(0 + fN * cov_dim, 0)) * std::sin(m_iter(1 + fN * cov_dim, 0));
        H(2, fN * cov_dim) = 1. / xi_iter(0, 0) * std::sin(m_iter(0 + fN * cov_dim, 0));
        // dphi
        H(0, 1 + fN * cov_dim) = 1. / xi_iter(0, 0) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::sin(m_iter(1 + fN * cov_dim, 0));
        H(1, 1 + fN * cov_dim) = -1. / xi_iter(0, 0) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::cos(m_iter(1 + fN * cov_dim, 0));
    }

    // PVA constraint
    if (fPvaConstraint)
    {
        H.ResizeTo(3, fyDim);
        H.Zero();

        // Daughter variables
        for (int q = 0; q < fN; q++)
        {
            // d(1/p)
            H(0, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 2) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
            H(1, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 2) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            H(2, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 2) * std::cos(m_iter(1 + q * cov_dim, 0));

            // dtht
            H(0, 1 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
            H(1, 1 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            H(2, 1 + q * cov_dim) = -1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0));

            // dphi
            H(0, 2 + q * cov_dim) = -1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            H(1, 2 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
        }

        // Mother variables
        // dtht
        H(0, fN * cov_dim) = -1. / xi_iter(0, 0) * std::cos(m_iter(0 + fN * cov_dim, 0)) * std::cos(m_iter(1 + fN * cov_dim, 0));
        H(1, fN * cov_dim) = -1. / xi_iter(0, 0) * std::cos(m_iter(0 + fN * cov_dim, 0)) * std::sin(m_iter(1 + fN * cov_dim, 0));
        H(2, fN * cov_dim) = 1. / xi_iter(0, 0) * std::sin(m_iter(0 + fN * cov_dim, 0));
        // dphi
        H(0, 1 + fN * cov_dim) = 1. / xi_iter(0, 0) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::sin(m_iter(1 + fN * cov_dim, 0));
        H(1, 1 + fN * cov_dim) = -1. / xi_iter(0, 0) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::cos(m_iter(1 + fN * cov_dim, 0));
    }

    // 4C or missing particle constraint
    if (f4Constraint || fMissingParticleConstraint)
    {

        H.ResizeTo(4, fyDim);
        H.Zero();

        for (int q = 0; q < fN; q++)
        {

            // d(1/p)
            H(0, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 2) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
            H(1, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 2) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            H(2, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 2) * std::cos(m_iter(1 + q * cov_dim, 0));
            H(3, q * cov_dim) = -1. / std::pow(m_iter(0 + q * cov_dim, 0), 3) * 1. / sqrt(pow(1. / (m_iter(0 + q * cov_dim, 0)), 2) + std::pow(fM[q], 2));

            // dtht
            H(0, 1 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
            H(1, 1 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * std::cos(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            H(2, 1 + q * cov_dim) = -1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0));

            // dphi
            H(0, 2 + q * cov_dim) = -1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::sin(m_iter(2 + q * cov_dim, 0));
            H(1, 2 + q * cov_dim) = 1. / m_iter(0 + q * cov_dim, 0) * std::sin(m_iter(1 + q * cov_dim, 0)) * std::cos(m_iter(2 + q * cov_dim, 0));
        }
    }

    // Mass constraint
    if (fMassConstraint)
    {
        H.ResizeTo(1, fyDim);
        H.Zero();

        double Px = 0., Py = 0., Pz = 0., E = 0.;

        for (int q = 0; q < fN; q++)
        {
            E += std::sqrt((1. / m_iter(0 + q * cov_dim, 0)) *
                               (1. / m_iter(0 + q * cov_dim, 0)) +
                           fM[q] * fM[q]);
            Px += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::sin(m_iter(1 + q * cov_dim, 0)) *
                  std::cos(m_iter(2 + q * cov_dim, 0));
            Py += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::sin(m_iter(1 + q * cov_dim, 0)) *
                  std::sin(m_iter(2 + q * cov_dim, 0));
            Pz += (1. / m_iter(0 + q * cov_dim, 0)) *
                  std::cos(m_iter(1 + q * cov_dim, 0));
        }

        for (int q = 0; q < fN; q++)
        {
            double Pi = 1. / m_iter(0 + q * cov_dim, 0);
            double Ei = std::sqrt(Pi * Pi + fM[q] * fM[q]);

            H(0, 0 + q * cov_dim) =
                -2 * E * (std::pow(Pi, 3) / Ei) +
                2 * std::pow(Pi, 2) * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::cos(m_iter(2 + q * cov_dim, 0)) * Px +
                2 * std::pow(Pi, 2) * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::sin(m_iter(2 + q * cov_dim, 0)) * Py +
                2 * std::pow(Pi, 2) * std::cos(m_iter(1 + q * cov_dim, 0)) * Pz;

            H(0, 1 + q * cov_dim) =
                -2 * Pi * std::cos(m_iter(1 + q * cov_dim, 0)) *
                    std::cos(m_iter(2 + q * cov_dim, 0)) * Px -
                2 * Pi * std::cos(m_iter(1 + q * cov_dim, 0)) *
                    std::sin(m_iter(2 + q * cov_dim, 0)) * Py +
                2 * Pi * std::sin(m_iter(1 + q * cov_dim, 0)) * Pz;

            H(0, 2 + q * cov_dim) =
                2 * Pi * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::sin(m_iter(2 + q * cov_dim, 0)) * Px -
                2 * Pi * std::sin(m_iter(1 + q * cov_dim, 0)) *
                    std::cos(m_iter(2 + q * cov_dim, 0)) * Py;

            H(0, 3 + q * cov_dim) = 0.;
            H(0, 4 + q * 4) = 0.;
        }
    }

    return H;
}

// Jacobian (derivative of constraint equations with respect to unmeasured variables)
TMatrixD KinFitter::Fxi_eval(const TMatrixD &m_iter, const TMatrixD &xi_iter)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KinFitter::Fxi_eval() -----------------" << std::endl;
    }

    TMatrixD H;

    if (f3Constraint)
    {
        H.ResizeTo(4, 1);
        H.Zero();

        // d(1/p)
        H(0, 0) = 1. / std::pow(xi_iter(0, 0), 2) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::cos(m_iter(1 + fN * cov_dim, 0));
        H(1, 0) = 1. / std::pow(xi_iter(0, 0), 2) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::sin(m_iter(1 + fN * cov_dim, 0));
        H(2, 0) = 1. / std::pow(xi_iter(0, 0), 2) * std::cos(m_iter(0 + fN * cov_dim, 0));
        H(3, 0) = 1. / std::pow(xi_iter(0, 0), 3) * 1. / sqrt(pow(1. / xi_iter(0, 0), 2) + std::pow(fM[fN], 2));
    }

    if (fPvaConstraint)
    {
        H.ResizeTo(3, 1);
        H.Zero();

        // d(1/p)
        H(0, 0) = 1. / std::pow(xi_iter(0, 0), 2) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::cos(m_iter(1 + fN * cov_dim, 0));
        H(1, 0) = 1. / std::pow(xi_iter(0, 0), 2) * std::sin(m_iter(0 + fN * cov_dim, 0)) * std::sin(m_iter(1 + fN * cov_dim, 0));
        H(2, 0) = 1. / std::pow(xi_iter(0, 0), 2) * std::cos(m_iter(0 + fN * cov_dim, 0));
    }

    if (fMissingParticleConstraint)
    {
        H.ResizeTo(4, 3);
        H.Zero();

        H(0, 0) = 1;
        H(1, 1) = 1;
        H(2, 2) = 1;

        H(3, 0) = xi_iter(0, 0) / sqrt(pow(xi_iter(0, 0), 2) + std::pow(xi_iter(1, 0), 2) + std::pow(xi_iter(2, 0), 2) + std::pow(fM[fN], 2));
        H(3, 1) = xi_iter(1, 0) / sqrt(pow(xi_iter(0, 0), 2) + std::pow(xi_iter(1, 0), 2) + std::pow(xi_iter(2, 0), 2) + std::pow(fM[fN], 2));
        H(3, 2) = xi_iter(2, 0) / sqrt(pow(xi_iter(0, 0), 2) + std::pow(xi_iter(1, 0), 2) + std::pow(xi_iter(2, 0), 2) + std::pow(fM[fN], 2));
    }

    return H;
}

Bool_t KinFitter::fit()
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- KinFitter::fit() -----------" << std::endl;
        std::cout << "Vertex constraint set: " << fVtxConstraint << std::endl;
        std::cout << "3C set: " << f3Constraint << std::endl;
        std::cout << "PVA constraint set: " << fPvaConstraint << std::endl;
        std::cout << "4C set: " << f4Constraint << std::endl;
        std::cout << "Momentum constraint set: " << fMissingParticleConstraint << std::endl;
        std::cout << "Mass constraint set: " << fMassConstraint << std::endl;
        std::cout << "Missing Mass constraint set: " << fMissingMassConstraint << std::endl;
        std::cout << "Mass + Vertex constraint set: " << fMassVtxConstraint << std::endl;
        std::cout << "" << std::endl;
    }

    TMatrixD alpha0(fyDim, 1), alpha(fyDim, 1);
    TMatrixD xi0(x), xi(x), neu_xi(x);
    TMatrixD A0(y), V0(V);
    alpha0 = y;
    alpha = alpha0;
    double chi2 = 1e6;

    // Calculating the inverse of the original covariance matrix that is not changed in the iterations
    TMatrixD V0_inv(V);
    V0_inv.Invert();

    if (f4Constraint == false && f3Constraint == false && fVtxConstraint == false && fMissingParticleConstraint == false &&
        fMassConstraint == false && fMissingMassConstraint == false && fMassVtxConstraint == false && fPvaConstraint == false)
    {
        std::cout << "FATAL: No constraint is chosen, please add constraint!" << std::endl;
        abort();
    }

    if (f3Constraint || fPvaConstraint)
    {
        xi0(0, 0) = 1 / fMother.P();
        xi = xi0;
    }

    if (fMissingParticleConstraint)
    {
        xi0 = calcMissingMom(alpha0);
        xi = xi0;
    }

    if (fVerbose > 1)
    {
        cout << " Calculate Feta" << endl;
    }
    TMatrixD D = Feta_eval(alpha, xi);
    TMatrixD DT(D.GetNcols(), D.GetNrows());
    if (fVerbose > 1)
    {
        cout << " Calculate f" << endl;
    }
    TMatrixD d = f_eval(alpha, xi);
    TMatrixD D_xi(d.GetNrows(), 1), DT_xi(1, d.GetNrows());
    if (d.GetNrows() - fNdf > 0)
        D_xi.ResizeTo(d.GetNrows(), d.GetNrows() - fNdf);
    DT_xi.ResizeTo(d.GetNrows() - fNdf, d.GetNrows()); // check dimension if other fitters are added
    D_xi.Zero();
    DT_xi.Zero();
    if (f3Constraint || fMissingParticleConstraint || fPvaConstraint)
    {
        if (fVerbose > 1)
        {
            cout << " Calculate Fxi" << endl;
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
        // Calculate r
        if (fVerbose > 1)
        {
            cout << " Calculate r" << endl;
        }
        TMatrixD r = d + D * delta_alpha;
        DT.Transpose(D);
        // Calculate S
        if (fVerbose > 1)
        {
            cout << " Calculate S" << endl;
        }
        VD = D * V0 * DT;
        VD.Invert();
        if (f3Constraint || fMissingParticleConstraint || fPvaConstraint)
        {
            DT_xi.Transpose(D_xi);
            if (fVerbose > 1)
            {
                cout << " Calculate Sxi" << endl;
            }
            VDD = DT_xi * VD * D_xi;
            VDD.Invert();
        }

        // Calculate values for next iteration
        TMatrixD lambda(d); // Lagrange multiplier
        lambda.Zero();
        if (f3Constraint || fMissingParticleConstraint || fPvaConstraint)
        {
            if (fVerbose > 1)
            {
                cout << " Calculate neuxi" << endl;
            }
            neu_xi = xi - VDD * DT_xi * VD * r;
        }
        TMatrixD delta_xi = neu_xi - xi;
        if (fVerbose > 1)
        {
            cout << " Calculate lambda" << endl;
        }
        lambda = VD * (r + D_xi * delta_xi);
        TMatrixD lambdaT(lambda.GetNcols(), lambda.GetNrows());
        lambdaT.Transpose(lambda);
        TMatrixD neu_alpha(fyDim, 1);
        if (fVerbose > 1)
        {
            cout << " Calculate neueta" << endl;
        }
        neu_alpha = alpha0 - V0 * DT * lambda;
        delta_alpha = alpha0 - neu_alpha;

        // Calculate new chi2
        TMatrixD chisqrd(1, 1);
        TMatrixD delta_alphaT(delta_alpha.GetNcols(), delta_alpha.GetNrows());
        delta_alphaT.Transpose(delta_alpha);
        TMatrixD two(1, 1);
        two(0, 0) = 2;
        if (fVerbose > 1)
        {
            cout << " Calculate chi2" << endl;
        }
        chisqrd = delta_alphaT * V0_inv * delta_alpha + two * lambdaT * f_eval(neu_alpha, neu_xi);

        // For testing the convergence
        fIteration = q;
        double dNorm = 0;
        double alphaNorm = 0;
        TMatrixD delta_alpha_it = alpha - neu_alpha;
        for (int i = 0; i < d.GetNrows(); i++)
            dNorm += pow(d(i, 0), 2);
        dNorm = sqrt(dNorm);
        for (int i = 0; i < delta_alpha.GetNrows(); i++)
            alphaNorm += pow((delta_alpha_it(i, 0)/alpha0(i, 0)), 2);
        if (f3Constraint || fMissingParticleConstraint || fPvaConstraint)
        {
            for (int i = 0; i < delta_xi.GetNrows(); i++)
                alphaNorm += pow((delta_xi(i, 0)/xi0(i, 0)), 2);
        }
        alphaNorm = sqrt(alphaNorm);
        if (fVerbose > 2)
        {
            cout << "Chi2: " << fabs(chi2 - chisqrd(0, 0)) << endl;
            cout << "Delta apha norm: " << alphaNorm << endl;
            cout << "Constraints norm: " << dNorm << endl;
        }
        if (fabs(chi2 - chisqrd(0, 0)) < fConvergenceCriterionChi2 && dNorm < fConvergenceCriterionD && alphaNorm < fConvergenceCriterionAlpha)
        {
            fConverged = true;
            chi2 = chisqrd(0, 0);
            alpha = neu_alpha;
            if (f3Constraint || fMissingParticleConstraint || fPvaConstraint)
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
        if (f3Constraint || fMissingParticleConstraint || fPvaConstraint)
            xi = neu_xi;
        D = Feta_eval(alpha, xi);
        if (f3Constraint || fMissingParticleConstraint || fPvaConstraint)
            D_xi = Fxi_eval(alpha, xi);
        d = f_eval(alpha, xi);
    }

    if (fVerbose > 0)
    {
        cout << " Exit iterations" << endl;
    }
    y = alpha;
    x = xi;
    if (fMissingParticleConstraint)
        fMissDaughter.SetXYZM(x(0, 0), x(1, 0), x(2, 0), fM[fN]);
    fChi2 = chi2;
    fProb = TMath::Prob(chi2, fNdf);

    // Update covariance
    if (fVerbose > 1)
    {
        cout << " Calculate Vneu" << endl;
    }
    if (f3Constraint || fMissingParticleConstraint || fPvaConstraint)
    {
        TMatrixD matrix = DT * VD * D_xi;
        TMatrixD matrixT(matrix.GetNcols(), matrix.GetNrows());
        matrixT.Transpose(matrix);
        TMatrixD invertedMatrix = DT_xi * VD * D_xi;
        invertedMatrix.Invert();
        V = V0 - V0 * (DT * VD * D - (matrix * invertedMatrix * matrixT)) * V0;
        Vx = invertedMatrix;
    }
    if (fVtxConstraint || f4Constraint || fMassConstraint || fMissingMassConstraint || fMassVtxConstraint)
        V = V0 - V0 * DT * VD * D * V0;

    // -----------------------------------------
    // Pulls
    // -----------------------------------------
    if (fVerbose > 1)
    {
        cout << " Calculate pulls" << endl;
    }
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
    if (f3Constraint || fPvaConstraint)
        updateMother();

    if (fNumIterations == 1)
        return kTRUE; // For number of iterations greater than 1
    else
        return fConverged; // For number of iterations equal to 1
}

void KinFitter::updateDaughters()
{
    if (fVerbose > 1)
    {
        cout << " Update daughters" << endl;
    }
    for (int val = 0; val < fN; ++val)
    {
        KFitParticle &cand = fCands[val];
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
        cand.setMomentum(1. / y(0 + val * cov_dim, 0));
        cand.setThetaRad(y(1 + val * cov_dim, 0));
        cand.setPhiRad(y(2 + val * cov_dim, 0));
        cand.setR(y(3 + val * cov_dim, 0));
        cand.setZ(y(4 + val * cov_dim, 0));

        // ---------------------------------------------------
        // Set covariance
        // ---------------------------------------------------
        TMatrixD cov(5, 5);
        cov(0, 0) = V(0 + val * cov_dim, 0 + val * cov_dim);
        cov(1, 1) = V(1 + val * cov_dim, 1 + val * cov_dim);
        cov(2, 2) = V(2 + val * cov_dim, 2 + val * cov_dim);
        cov(3, 3) = V(3 + val * cov_dim, 3 + val * cov_dim);
        cov(4, 4) = V(4 + val * cov_dim, 4 + val * cov_dim);
        cand.setCovariance(cov);
        // ---------------------------------------------------
    }
}

void KinFitter::updateMother()
{
    if (fVerbose > 1)
    {
        cout << " Update mother" << endl;
    }

    KFitParticle &mother = fMother;
    double Px = (1. / x(0, 0)) *
                std::sin(y(0 + fN * cov_dim, 0)) *
                std::cos(y(1 + fN * cov_dim, 0));
    double Py = (1. / x(0, 0)) *
                std::sin(y(0 + fN * cov_dim, 0)) *
                std::sin(y(1 + fN * cov_dim, 0));
    double Pz =
        (1. / x(0, 0)) * std::cos(y(0 + fN * cov_dim, 0));
    double M = fabs(fM[fN]);
    mother.SetXYZM(Px, Py, Pz, M);
    mother.setR(y(2 + fN * cov_dim, 0));
    mother.setZ(y(3 + fN * cov_dim, 0));

    // ---------------------------------------------------------------------------
    // Set covariance
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

void KinFitter::setConvergenceCriteria(double val1, double val2, double val3)
{
    fConvergenceCriterionChi2 = val1;
    fConvergenceCriterionD = val2;
    fConvergenceCriterionAlpha = val3;
}

ClassImp(KinFitter)
