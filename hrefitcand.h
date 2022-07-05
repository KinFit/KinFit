#ifndef HREFITCAND_H
#define HREFITCAND_H

// ROOT includes
#include <TMath.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>

class HRefitCand : public TLorentzVector
{
private:
    Double_t fMomentum, fR, fZ, fTheta, fPhi;
    Bool_t fIsForward, fIsUsed;
    TMatrixD fCov;

  
public:
    HRefitCand(); //Does that work like this? See https://os.mbed.com/users/fpucher/code/HIM0Board/wiki/Vererbung-in-C%2B%2B
    ~HRefitCand(){};
    void setR(Double_t val) { fR = val; }
    void setZ(Double_t val) { fZ = val; }
    void setIsUsed(Double_t val) { fIsUsed = val; }
    void setIsForward(Bool_t val) { fIsForward = val; }
    void setTheta(Double_t val) { fTheta = val; }
    void setPhi(Double_t val) { fPhi = val; }
    void setMomentum(Double_t val){ fMomentum = val; }

    void setCovariance(const TMatrixD &cov);
    Double_t getR() const { return fR; }
    Double_t getZ() const { return fZ; }
    Double_t getTheta() const { return fTheta; }
    Double_t getPhi() const { return fPhi; }
    Double_t getMomentum() const { return fMomentum; }
    Bool_t getIsForward() const { return fIsForward; }
    Bool_t getIsUsed() const { return fIsUsed; };
    TMatrixD getCovariance() const { return fCov; }

    void reset(HRefitCand* cand);
    void update(HRefitCand* cand);
};

#endif /* HREFITCAND_H */
