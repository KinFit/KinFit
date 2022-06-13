#ifndef HREFITCAND_H
#define HREFITCAND_H

// ROOT includes
#include <TMath.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>

class HRefitCand : public TLorentzVector
{
private:
    TLorentzVector *cand;
    Double_t fMomentum, fTheta, fPhi, fR, fZ;
    TMatrixD fCov;

public:
    HRefitCand(TLorentzVector *cand, Double_t R, Double_t Z);
    HRefitCand(); //Does that work like this? See https://os.mbed.com/users/fpucher/code/HIM0Board/wiki/Vererbung-in-C%2B%2B
    ~HRefitCand(){};
    void setMomentum(Double_t val){ fMomentum = val; }
    void setTheta(Double_t val) { fTheta = val; }
    void setPhi(Double_t val) { fPhi = val; }
    void setR(Double_t val) { fR = val; }
    void setZ(Double_t val) { fZ = val; }
    void setCovariance(const TMatrixD &cov);

    Double_t getMomentum() const { return fMomentum; }
    Double_t getTheta() const { return fTheta; }
    Double_t getPhi() const { return fPhi; }
    Double_t getR() const { return fR; }
    Double_t getZ() const { return fZ; }
    TMatrixD getCovariance() const { return fCov; }

    void reset();
    void update();
};

#endif /* HREFITCAND_H */
