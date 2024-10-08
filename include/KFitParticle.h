//****************************************************************************
//*                    This file is part of KinFit.                          *
//*                                                                          *
//*            	KinFit is distributed under the terms of the                 *
//*              GNU General Public License (GPL) version 3,                 *
//*                 copied verbatim in the file "LICENSE".                   *
//*                                                                          *
//*  				           Copyright 2024                                *
//*		         GSI Helmholtzzentrum für Schwerionenforschung               *
//* 	      This software is distributed under the terms of the            *
//*	          GNU General Public Licence version 3 (GPL Version 3)           *
//*		      			     				                                 *
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
    double fR; // Closest distance to beamline [mm]
    double fZ; // Point along beamline where track is closest to it [mm]
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
    void setCovariance(const TMatrixD &cov); // Sets the covariance matrix
    void setPid(int val) { fPid = val; } // Sets the PID of the particle (user defined)
    void setTrackId(int val) { fTrackId = val; } // Sets the track id

    double getMomentum() const { return fMomentum; } // Returns the momentum [MeV/c]
    double getThetaRad() const { return fTheta; } // Returns the polar angle [radians]
    double getPhiRad() const { return fPhi; } // Returns the azimuthal angle [radians]
    double getThetaDeg() const { return fTheta*TMath::RadToDeg(); } // Returns the polar angle [deg]
    double getPhiDeg() const { return fPhi*TMath::RadToDeg(); } // Returns the azimuthal angle [deg]
    double getR() const { return fR; } // Returns R [mm]
    double getZ() const { return fZ; } // Returns Z [mm]
    TMatrixD getCovariance() const { return fCov; } // Returns the covariance matrix
    int getPid() const { return fPid; } // Returns the PID of the particle (user defined)
    int getTrackId() const { return fTrackId; } // Returns the track id
};

#endif /* KFITPARTICLE_H */
