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

    ClassDef(KFitParticle, 2);

private:
    TLorentzVector cand; // TLorentzVector the KFitParticle inherits from
    double fMomentum; // Momentum [MeV]
    double fTheta; // Polar angle [rad]
    double fPhi; // Azimuthal angle [rad]
    double fR; // Closest distance to beamline [cm]
    double fZ; // Point along beamline where track is closest to it [cm]
    int fPid; // PID code for the particle spices
    int fTrackId; // Track id

    /** Covariance matrix of the KFitParticle
    * Diagonal entries correspond to the covariances 
    * in the parameters in the following order
    * 
    * -----------------------------
    * | 1/p                       |
    * |     theta                 |
    * |            phi            |
    * |                   R       |
    * |                        Z  |
    * -----------------------------
    * 
    * Off diagonal elements corresponds to the
    * correlations between the parameters
    */
    TMatrixD fCov;

public:
    KFitParticle(TLorentzVector cand, double R, double Z);
    KFitParticle(TLorentzVector cand, double X, double Y, double Z);
    KFitParticle(); 
    ~KFitParticle(){};
    void setMomentum(double val) { fMomentum = val; }
    void setThetaRad(double val) { fTheta = val; cand.SetTheta(val);}
    void setPhiRad(double val) { fPhi = val; cand.SetTheta(val);}
    void setThetaDeg(double val) { fTheta = val*TMath::DegToRad(); cand.SetTheta(val*TMath::DegToRad()); }
    void setPhiDeg(double val) { fPhi = val*TMath::DegToRad(); cand.SetTheta(val*TMath::DegToRad()); }
    void setR(double val) { fR = val; }
    void setZ(double val) { fZ = val; }
    void setCovariance(const TMatrixD &cov);
    void setPid(int val) { fPid = val; }
    void setTrackId(int val) { fTrackId = val; }

    double getMomentum() const { return fMomentum; }
    double getThetaRad() const { return fTheta; }
    double getPhiRad() const { return fPhi; }
    double getThetaDeg() const { return fTheta*TMath::RadToDeg(); }
    double getPhiDeg() const { return fPhi*TMath::RadToDeg(); }
    double getR() const { return fR; }
    double getZ() const { return fZ; }
    TMatrixD getCovariance() const { return fCov; }
    int getPid() const { return fPid; }
    int getTrackId() const { return fTrackId; }
};

#endif /* KFITPARTICLE_H */
