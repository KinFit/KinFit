// ROOT includes
#include "hrefitcand.h"

HRefitCand::HRefitCand()
    : TLorentzVector(), fMomentum(-9.0), fTheta(), fPhi(-9.0), fR(-9.0), fZ(-9.0), fIsForward(false)
{
    //HRefitCand cand = new HRefitCand();
    //fR = cand->getR();
    //fZ = cand->getZ();
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

void HRefitCand::reset(HRefitCand* cand) {
    *(TLorentzVector*)this = *cand;
}

void HRefitCand::update(HRefitCand* cand) {
    *((TLorentzVector*)cand) = *this;
}
