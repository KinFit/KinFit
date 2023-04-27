/**
 * KFitParticle.h
 *
 *@updated 27.04.2023
 *@version v1.0.0 
 *
 * Object containing track parameters
 * as well as the covariance matrix
 */

#ifndef KFITPARTICLE_H
#define KFITPARTICLE_H

// ROOT includes
#include <TMath.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>

class KFitParticle : public TLorentzVector
{

    ClassDef(KFitParticle, 1);

private:
    TLorentzVector *cand; // TLorentzVector the KFitParticle
    // inherits from
    Double_t fMomentum; // Momentum [MeV]
    Double_t fTheta; // Polar angle [rad]
    Double_t fPhi; // Azimuthal angle [rad]
    Double_t fR; // Closest distance to beamline [cm]
    Double_t fZ; // Point along beamline where track is closest to it [cm]
    Int_t fPid; // PID code for the particle spices
    Int_t fTrackId; // Track id

    /** Covariance matrix of the KFitParticle
    /* Diagonal entries correspond to the covariances 
    /* in the parameters in the following order
    /* 
    /* -----------------------------
    /* | 1/p                       |
    /* |     theta                 |
    /* |            phi            |
    /* |                   R       |
    /* |                        Z  |
    /* -----------------------------
    /* 
    /* Off diagonal elements corresponds to the
    /* correlations between the parameters
    */
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
    void setTrackId(Int_t val) { fTrackId = val; }

    Double_t getMomentum() const { return fMomentum; }
    Double_t getTheta() const { return fTheta; }
    Double_t getPhi() const { return fPhi; }
    Double_t getR() const { return fR; }
    Double_t getZ() const { return fZ; }
    TMatrixD getCovariance() const { return fCov; }
    Int_t getPid() const { return fPid; }
    Int_t getTrackId() const { return fTrackId; }

    void reset();
    void update();
};

#endif /* KFITPARTICLE_H */
