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

#include "KFitParticle.h"

KFitParticle::KFitParticle(TLorentzVector cand, double R, double Z)
    : TLorentzVector(cand), cand(cand), fMomentum(cand.P()), fTheta(cand.Theta()), fPhi(cand.Phi()), fR(R), fZ(Z)
{
    fPid = -1;
    fTrackId = -1;
    fCov.ResizeTo(5, 5);
}

KFitParticle::KFitParticle(TLorentzVector cand, double X, double Y, double Z)
    : TLorentzVector(cand), cand(cand), fMomentum(cand.P()), fTheta(cand.Theta()), fPhi(cand.Phi())
{
    double deg2rad = TMath::DegToRad();

    // Base and direction vettor of particle candidate
    TVector3 base(X, Y, Z);
    TVector3 dir(TMath::Sin(fTheta) * TMath::Cos(fPhi),
                 TMath::Sin(fTheta) * TMath::Sin(fPhi),
                 TMath::Cos(fTheta));
    
    // Base and direction vector of beamline
    TVector3 beam_base(0,0,1); 
    TVector3 beam_dir(0,0,1);

    TVector3 cross = dir.Cross(beam_dir);
    TVector3 diff=base-beam_base;

    // fR = diff.Dot(cross)/(cross.Mag());
    fR = -((cross.Dot(diff)) / cross.Mag()); // R is signed

    double a = beam_base.Dot(beam_dir);
    double b = beam_dir.Dot(beam_dir);
    double c = base.Dot(beam_dir);
    double d = (beam_base.Dot(dir)) * (dir.Dot(beam_dir)) / dir.Dot(dir);
    double e = (beam_dir.Dot(dir)) * (dir.Dot(beam_dir)) / dir.Dot(dir);
    double f = (base.Dot(dir)) * (dir.Dot(beam_dir)) / dir.Dot(dir);
    double u1 = (-a + c + d - f) / (b - e);

    fZ = beam_base.Z() + beam_dir.Z() * u1;

    double y = beam_base.Y() + dir.Y() * u1;

    fPid = -1;
    fTrackId = -1;
}

KFitParticle::KFitParticle()
    : TLorentzVector()
{
    cand.SetXYZM(0,0,0,0);
    fR = 0;
    fZ = 0;
    fPid = -1;
    fTrackId = -1;
    fCov.ResizeTo(5, 5);
}

void KFitParticle::setCovariance(const TMatrixD &cov)
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
