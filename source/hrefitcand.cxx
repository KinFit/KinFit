// ROOT includes
#include "hrefitcand.h"

HRefitCand::HRefitCand(TLorentzVector* cand, Double_t R, Double_t Z)
    : TLorentzVector(*cand), cand(cand), fMomentum(cand->P()), fTheta(cand->Theta()), fPhi(cand->Phi()), fR(R), fZ(Z)
{
}

HRefitCand::HRefitCand()
    : TLorentzVector()
{
    cand = new TLorentzVector();
    fR = 0;
    fZ = 0;
}

void HRefitCand::setCovariance(const TMatrixD& cov)
{
    // 0 = 1/p
    // 1 = theta
    // 2 = phi
    // 3 = R
    // 4 = z
    fCov.ResizeTo(5, 5);
    fCov = cov;
}

void HRefitCand::reset() {
    *(TLorentzVector*)this = *cand;
}

void HRefitCand::update() {
    *((TLorentzVector*)cand) = *this;
}
