// ROOT includes
#include "hrefitcand.h"

HRefitCand::HRefitCand(TLorentzVector *cand, Double_t R, Double_t Z)
    : TLorentzVector(*cand), cand(cand), fMomentum(cand->P()), fTheta(cand->Theta()), fPhi(cand->Phi()), fR(R), fZ(Z)
{
    fPid = -1;
}

HRefitCand::HRefitCand(TLorentzVector *cand, Double_t X, Double_t Y, Double_t Z)
    : TLorentzVector(*cand), cand(cand), fMomentum(cand->P()), fTheta(cand->Theta()), fPhi(cand->Phi())
{
    Double_t deg2rad = TMath::DegToRad();

    TVector3 base(X, Y, Z);
    TVector3 beam_base, beam_dir;
    TVector3 dir(TMath::Sin(fTheta) * TMath::Cos(fPhi),
                 TMath::Sin(fTheta) * TMath::Sin(fPhi),
                 TMath::Cos(fTheta));

    HParticleTool::calcSegVector(1, 0, 0, 0, beam_base, beam_dir);
    TVector3 POCA = HParticleTool::calculatePointOfClosestApproach(base, dir, beam_base, beam_dir);
    fR = base.Y() * TMath::Cos(fPhi) - base.X() * TMath::Sin(fPhi);
    fZ = POCA.Z();

    fPid = -1;
}

HRefitCand::HRefitCand()
    : TLorentzVector()
{
    cand = new TLorentzVector();
    fR = 0;
    fZ = 0;
    fPid = -1;
}

void HRefitCand::setCovariance(const TMatrixD &cov)
{
    // HADES default track parametrization
    // 0 = 1/p
    // 1 = theta
    // 2 = phi
    // 3 = R
    // 4 = z

    // An alternative parametrization
    // 0 = px
    // 1 = py
    // 2 = pz
    // 3 = x (vertex x position)
    // 4 = y (vertex y position)
    // 5 = z (vertex z position)

    if (cov.GetNoElements() == 5 * 5)
    {
        fCov.ResizeTo(5, 5);
        fCov = cov;
    }

    else if (cov.GetNoElements() == 6 * 6)
    {
        // TO-DO: do the error propagtion calculations
    }
}

void HRefitCand::reset()
{
    *(TLorentzVector *)this = *cand;
}

void HRefitCand::update()
{
    *((TLorentzVector *)cand) = *this;
}
