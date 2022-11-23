#ifndef KFITPARTICLE_H
#define KFITPARTICLE_H

// ROOT includes
#include <TMath.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>

// framework includes
//#include "hparticletool.h"

class KFitParticle : public TLorentzVector
{

    ClassDef(KFitParticle, 1);

private:
    TLorentzVector *cand;
    Double_t fMomentum, fTheta, fPhi, fR, fZ;
    Int_t fPid;
    TMatrixD fCov;

public:
    KFitParticle(TLorentzVector *cand, Double_t R, Double_t Z);
    KFitParticle(TLorentzVector *cand, Double_t X, Double_t Y, Double_t Z);
    KFitParticle(); 
    ~KFitParticle(){};
    void setMomentum(Double_t val) { fMomentum = val; }
    void setTheta(Double_t val) { fTheta = val; }
    void setPhi(Double_t val) { fPhi = val; }
    void setR(Double_t val) { fR = val; }
    void setZ(Double_t val) { fZ = val; }
    void setCovariance(const TMatrixD &cov);
    void setPid(Int_t val) { fPid = val; }

    Double_t getMomentum() const { return fMomentum; }
    Double_t getTheta() const { return fTheta; }
    Double_t getPhi() const { return fPhi; }
    Double_t getR() const { return fR; }
    Double_t getZ() const { return fZ; }
    TMatrixD getCovariance() const { return fCov; }
    Int_t getPid() const { return fPid; }

    void reset();
    void update();
};

#endif /* KFITPARTICLE_H */
