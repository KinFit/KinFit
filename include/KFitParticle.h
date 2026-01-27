//****************************************************************************
//*                    This file is part of KinFit.                          *
//*                                                                          *
//*            	KinFit is distributed under the terms of the                 *
//*              GNU General Public License (GPL) version 3,                 *
//*                 copied verbatim in the file "LICENSE".                   *
//*                                                                          *
//*  				Copyright 2024                               *
//*		GSI Helmholtzzentrum für Schwerionenforschung                *
//* 	     This software is distributed under the terms of the             *
//*	     GNU General Public Licence version 3 (GPL Version 3)            *
//*		      			     				     *
//*     The copyright holders are listed in the file "COPYRIGHTHOLDERS".     *
//*               The authors are listed in the file "AUTHORS".              *
//****************************************************************************

/**
 * KFitParticle.h
 *
 *@updated 03.08.2023
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
    double fMomentum; // Momentum [MeV/c]
    double fTheta; // Polar angle [rad]
    double fPhi; // Azimuthal angle [rad]
    double fR = 9.0; // Closest distance to beamline [mm]
    double fZ = 9.0; // Point along beamline where track is closest to it [mm]

    double fGenMomentum;
    double fGenTheta;
    double fGenPhi;
    double fGenR;
    double fGenZ;

    double fVtxX = -999999999.9;
    double fVtxY = -999999999.9;
    double fVtxZ = -999999999.9;
    double fGenVtxX = -999999999.9;
    double fGenVtxY = -999999999.9;
    double fGenVtxZ = -999999999.9;

    double fSector; // HADES sector 0-5

    double fRapidity = -9999999999.9;
    double fPt = -9999999999.9;
    double fPx = -9999999999.9; 
    double fPy = -9999999999.9; 
    double fPz = -9999999999.9;  

    double fBeta = -9999999999.9;
    double fChi2 = -9999999999.9;
    double fTOF = -9999999999.9;

    int fPID = -99999;
    int fParentPID = -99999;
    int fGrandParentPID = -99999;   

    int fCharge = 0;

    int fPid; // PID code for the particle spices
    int fTrackId; // Track id, can be used to differentiate different track types, e.g. reconstructed in different parts of the detector

    /** Covariance matrix of the KFitParticle
    * Diagonal entries correspond to the variances 
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
    * covariances of the parameters
    */
    TMatrixD fCov;

public:
    /** Default Constructor **/
    KFitParticle(); 

    /** @brief Constructor
    * @param cand - TLorentzVector with theta, phi and momentum set
    * @param R - Closest distance between particle cand and beam line
    * @param Z - Point along the beamline where the particle passes closest to it
    */
    KFitParticle(TLorentzVector cand, double R, double Z);

    /** @brief Constructor
    * @param cand - TLorentzVector with theta, phi and momentum set
    * @param X - Starting X-position of track
    * @param Y - Starting Y-position of track
    * @param Z - Starting Z-position of track
    */
    KFitParticle(TLorentzVector cand, double X, double Y, double Z);
    
    /** Default Destructor **/
    ~KFitParticle(){};

    void setMomentum(double val) { fMomentum = val; } // Sets the momentum [MeV/c]
    void setThetaRad(double val) { fTheta = val; cand.SetTheta(val); } // Sets the polar angle [radians]
    void setPhiRad(double val) { fPhi = val; cand.SetTheta(val); } // Sets the azimuthal angle [radians]
    void setThetaDeg(double val) { fTheta = val*TMath::DegToRad(); cand.SetTheta(val*TMath::DegToRad()); } // Sets the polar angle [deg]
    void setPhiDeg(double val) { fPhi = val*TMath::DegToRad(); cand.SetTheta(val*TMath::DegToRad()); } // Sets the azimuthal angle [deg]
    void setR(double val) { fR = val; } // Sets R [mm]
    void setZ(double val) { fZ = val; } // Sets Z [mm]
    void setCharge(double val) { fCharge = val; } // Sets Charge
    void setBeta(double val) { fBeta = val; } // Sets Beta
    void setChi2(double val) { fChi2 = val; } // Sets Chi2
    void setCovariance(const TMatrixD &cov); // Sets the covariance matrix
    void setPid(int val) { fPid = val; } // Sets the PID of the particle (user defined)
    void setTrackId(int val) { fTrackId = val; } // Sets the track id
    void setTOF(int val) { fTOF = val; } // Sets the tof
    void setRapidity(int val) { fRapidity = val; } // Sets the rapidity
    void setPt(double val) { fPt = val; } // Sets the transverse momentum [MeV/c]
    void setPx(double val) { fPx = val; } // Sets the x momentum [MeV/c]
    void setPy(double val) { fPy = val; } // Sets the y momentum [MeV/c]
    void setPz(double val) { fPz = val; } // Sets the z momentum [MeV/c]
    void setVertexX(double val) { fVtxX = val; } // Sets the x-vertex [mm]
    void setVertexY(double val) { fVtxY = val; } // Sets the y-vertex [mm]
    void setVertexZ(double val) { fVtxZ = val; } // Sets the z-vertex [mm]
    void setGenVertexX(double val) { fGenVtxX = val; } // Sets the generated x-vertex [mm]
    void setGenVertexY(double val) { fGenVtxY = val; } // Sets the generated y-vertex [mm]
    void setGenVertexZ(double val) { fGenVtxZ = val; } // Sets the generated z-vertex [mm]

    // MC Truth Matching Information
    void setPID(int val) { fPID = val; }
    void setParentPID(int val) { fParentPID = val; }
    void setGrandParentPID(int val) {fGrandParentPID = val; }

    void setGenMomentum(double val) { fGenMomentum = val; } // Sets the momentum [MeV/c]
    void setGenThetaDeg(double val) { fGenTheta = val; } // Sets the polar angle [deg]
    void setGenPhiDeg(double val) { fGenPhi = val; } // Sets the azimuthal angle [deg]
    void setGenR(double val) { fGenR = val; } // Sets R [mm]
    void setGenZ(double val) { fGenZ = val; } // Sets Z [mm]

    void setSector(int val) { fSector = val; } // Sets the sector

    double getMomentum() const { return fMomentum; } // Returns the momentum [MeV/c]
    double getThetaRad() const { return fTheta; } // Returns the polar angle [radians]
    double getPhiRad() const { return fPhi; } // Returns the azimuthal angle [radians]
    double getThetaDeg() const { return fTheta*TMath::RadToDeg(); } // Returns the polar angle [deg]
    double getPhiDeg() const { return fPhi*TMath::RadToDeg(); } // Returns the azimuthal angle [deg]
    double getGenThetaDeg() const { return fGenTheta; } // Returns the generated polar angle [deg]
    double getGenPhiDeg() const { return fGenPhi; } // Returns the generated azimuthal angle [deg]
    double getR() const { return fR; } // Returns R [mm]
    double getZ() const { return fZ; } // Returns Z [mm]
    double getCharge() const { return fCharge; } // Returns  Charge
    double getBeta() const { return fBeta; } // Returns Beta
    double getChi2() const { return fChi2; } // Returns Chi2
    TMatrixD getCovariance() const { return fCov; } // Returns the covariance matrix
    int getPid() const { return fPid; } // Returns the PID of the particle (user defined)
    int getTrackId() const { return fTrackId; } // Returns the track id
    int getTOF() const { return fTOF; } // Returns the TOF
    double getPt() const { return fPt; } // Returns the transverse momentum [MeV/c]
    double getPx() const { return fPx; } // Returns the x momentum [MeV/c]
    double getPy() const { return fPy; } // Returns the y momentum [MeV/c]
    double getPz() const { return fPz; } // Returns the z momentum [MeV/c]
    double getVertexX() const { return fVtxX; } // Returns the x-vertex [mm]
    double getVertexY() const { return fVtxY; } // Returns the y-vertex [mm] 
    double getVertexZ() const { return fVtxZ; } // Returns the z-vertex [mm]
    double getGenVertexX() const { return fGenVtxX; } // Returns the generated x-vertex [mm]
    double getGenVertexY() const { return fGenVtxY; } // Returns the generated y-vertex [mm] 
    double getGenVertexZ() const { return fGenVtxZ; } // Returns the generated z-vertex [mm]
    double getRapidity() const { return fRapidity; } // Returns the rapidity

    // MC Truth Matching Information
    int getPID() const { return fPID; }
    int getParentPID() const { return fParentPID; }
    int getGrandParentPID() const { return fGrandParentPID; }

    double getGenMomentum() const { return fGenMomentum; }
    double getGenTheta() const { return fGenTheta; }
    double getGenPhi() const { return fGenPhi; }
    double getGenR() const { return fGenR; }
    double getGenZ() const { return fGenZ; }

    int getSector() const { return fSector; } // Returns the sector
};

#endif /* KFITPARTICLE_H */
