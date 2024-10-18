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

//#include "KFitParticle.h"
#include "KFitCoordinateConversion.h"

KFitCoordinateConversion::KFitCoordinateConversion()
: TObject()
{

}

void KFitCoordinateConversion::setVarsCart()
{

    fPx = (fMom * (sin(fTheta)) * (cos(fPhi)));
    fPy = (fMom * sin(fTheta) * sin(fPhi));
    fPz = (fMom * cos(fTheta));
    fX = (fR * cos((fPhi) + fPi2));
    fY = (fR * sin((fPhi) + fPi2));
}

void KFitCoordinateConversion::setVarsSpher()
{

    fMom = sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2));
    fTheta = atan((sqrt((pow(fPx, 2)) + (pow(fPy, 2)))) / fPz);
    fPhi = atan2(fPy, fPx);

    TVector3 trackdist;
    TVector3 trackbase;
    TVector3 beamdist;
    TVector3 beambase;

    trackdist.SetXYZ((sin(fTheta)) * (cos(fPhi)), (sin(fTheta)) * (sin(fPhi)), (cos(fTheta)));
    trackbase.SetXYZ(fX, fY, fZ);
    beamdist.SetXYZ(0, 0, 1);
    beambase.SetXYZ(0, 0, 1);

    TVector3 cross = trackdist.Cross(beamdist);
    TVector3 diff = trackbase - beambase;

    fR = -(((cross.Dot(diff))) / (cross.Mag()));
    // Rtest = (((fX * sin(fTheta) * sin(fPhi)) - (fY * sin(fTheta) * cos(fPhi)))) / (sqrt(pow((sin(fTheta) * sin(fPhi)), 2) + pow((sin(fTheta) * cos(fPhi)), 2)));
}

void KFitCoordinateConversion::setErrorsCart()
{

    fErrPx = sqrt(pow((fMom * cos(fTheta) * cos(fPhi) * fErrTheta), 2) + pow((fMom * sin(fTheta) * sin(fPhi) * fErrPhi), 2) + pow((sin(fTheta) * cos(fPhi) * fErrMom), 2) + 2 * (fMom * cos(fTheta) * cos(fPhi)) * ((fMom * sin(fTheta) * sin(fPhi))) * (fErrTheta * fErrPhi) + 2 * (fMom * cos(fTheta) * cos(fPhi)) * (sin(fTheta) * cos(fPhi)) * (fErrTheta * fErrMom) + 2 * ((fMom * sin(fTheta) * sin(fPhi)) * ((sin(fTheta) * cos(fPhi)) * (fErrMom * fErrPhi))));
    fErrPy = sqrt(pow((fMom * cos(fTheta) * sin(fPhi) * fErrTheta), 2) + pow((fMom * sin(fTheta) * cos(fPhi) * fErrPhi), 2) + pow((sin(fTheta) * sin(fPhi) * fErrMom), 2) + 2 * (fMom * cos(fTheta) * sin(fPhi) * ((fMom * sin(fTheta) * cos(fPhi)))) * (fErrTheta * fErrPhi) + 2 * (fMom * cos(fTheta) * sin(fPhi)) * (sin(fTheta) * sin(fPhi)) * (fErrTheta * fErrMom) + 2 * (fMom * sin(fTheta) * cos(fPhi)) * (sin(fTheta) * sin(fPhi)) * (fErrPhi * fErrMom));
    fErrPz = sqrt(pow((fMom * sin(fTheta) * fErrTheta), 2) + pow((cos(fTheta) * fErrMom), 2) + 2 * (fMom * sin(fTheta)) * (cos(fTheta) * (fErrTheta * fErrMom)));
    fErrX = sqrt((pow((fR * sin(fPhi + fPi2) * fErrPhi), 2) + pow((cos(fPhi + fPi2) * fErrR), 2)) + 2 * (fR * sin(fPhi + fPi2)) * (cos(fPhi + fPi2)) * (fErrPhi * fErrR));
    fErrY = sqrt((pow((fR * cos(fPhi + fPi2) * fErrPhi), 2) + pow((sin(fPhi + fPi2) * fErrR), 2)) + 2 * (fR * cos(fPhi + fPi2)) * (sin(fPhi + fPi2)) * (fErrPhi * fErrR));
}

void KFitCoordinateConversion::setErrorsSpher()
{

    fErrMom = sqrt((pow((fPx / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2)))), 2) * (pow(fErrPx, 2))) + (pow((fPy / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2)))), 2) * (pow(fErrPy, 2))) + (pow((fPz / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2)))), 2) * (pow(fErrPz, 2))) + 2 * ((fPx / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * ((fPy / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * (fErrPx * fErrPy) + 2 * ((fPx / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * ((fPz / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * (fErrPx * fErrPz) + 2 * ((fPz / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * ((fPy / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * (fErrPy * fErrPz));
    fErrTheta = sqrt(pow(((fPz * fPx) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2)))), 2) * (pow(fErrPx, 2)) + pow(((fPz * fPy) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2)))), 2) * (pow(fErrPy, 2)) + pow(((sqrt(((pow(fPx, 2)) + (pow(fPy, 2))))) / (((pow(fPx, 2)) + (pow(fPy, 2))) + (pow(fPz, 2)))), 2) * (pow(fErrPz, 2)) + 2 * (((fPz * fPx) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2))))) * (((fPz * fPy) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2))))) * (fErrPx * fErrPy) + 2 * (((fPz * fPx) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2))))) * ((sqrt(((pow(fPx, 2)) + (pow(fPy, 2))))) / (((pow(fPx, 2)) + (pow(fPy, 2))) + (pow(fPz, 2)))) * (fErrPx * fErrPz) + 2 * (((fPz * fPy) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2))))) * (((sqrt(((pow(fPx, 2)) + (pow(fPy, 2))))) / (((pow(fPx, 2)) + (pow(fPy, 2))) + (pow(fPz, 2))))) * (fErrPz * fErrPy)); // with correlation
    fErrPhi = sqrt((pow((fPy / ((pow(fPx, 2)) + (pow(fPy, 2)))), 2)) * (pow(fErrPx, 2)) + (pow((fPx / ((pow(fPx, 2)) + (pow(fPy, 2)))), 2) * (pow(fErrPy, 2))) + 2 * (fPy / ((pow(fPx, 2)) + (pow(fPy, 2)))) * (fPx / ((pow(fPx, 2)) + (pow(fPy, 2)))) * (fErrPx * fErrPy));

    double dRdX = (((sin(fPhi)) * (sin(fTheta))) / (sqrt(pow((sin(fPhi) * sin(fTheta)), 2) + pow((cos(fPhi) * sin(fTheta)), 2))));
    double dRdY = (-(cos(fPhi) * sin(fTheta)) / (sqrt(pow((sin(fPhi) * sin(fTheta)), 2) + pow((cos(fPhi) * sin(fTheta)), 2))));
    double dRdphi = ((sin(fTheta)) * (fY * sin(fPhi) + fX * cos(fPhi))) / (sqrt(pow((sin(fPhi) * sin(fTheta)), 2) + pow((cos(fPhi) * sin(fTheta)), 2)));
    double fErrR = (sqrt(pow((dRdX * fErrX), 2) + pow((dRdY * fErrY), 2) + pow((dRdphi * fErrPhi), 2) + 2 * (dRdX) * (dRdY) * (fErrX * fErrY) + 2 * (dRdX) * (dRdphi) * (fErrX * fErrPhi) + 2 * (dRdY) * (dRdphi) * (fErrY * fErrPhi)));
} 

void KFitCoordinateConversion::convertToSpherical(double cartParameters[11], double spherParameters[9])
{
    fX = cartParameters[0];
    fY = cartParameters[1];
    fZ = cartParameters[2];
    fPx = cartParameters[3];
    fPy = cartParameters[4];
    fPz = cartParameters[5];

    fErrX = cartParameters[6];
    fErrY = cartParameters[7];
    fErrZ = cartParameters[8];
    fErrPx = cartParameters[9];
    fErrPy = cartParameters[10];
    fErrPz = cartParameters[11];

    setVarsSpher();
    
    setErrorsSpher();

    spherParameters[0] = fMom;
    spherParameters[1] = fTheta;
    spherParameters[2] = fPhi;
    spherParameters[3] = fR;
    spherParameters[4] = fZ;
    spherParameters[5] = fErrMom;
    spherParameters[6] = fErrTheta;
    spherParameters[7] = fErrPhi;
    spherParameters[8] = fErrR;
    spherParameters[9] = fErrZ;

}

void KFitCoordinateConversion::convertToCart(double spherParameters[9], double cartParameters[11])
{

    fMom = spherParameters[0];
    fTheta = spherParameters[1];
    fPhi = spherParameters[2];
    fR = spherParameters[3];
    fZ = spherParameters[4];

    fErrMom = spherParameters[5];
    fErrTheta = spherParameters[6];
    fErrPhi = spherParameters[7];
    fErrR = spherParameters[8];
    fErrZ = spherParameters[9];
    
    setVarsCart();

    setErrorsCart();

    cartParameters[0]=fX;
    cartParameters[1]=fY;
    cartParameters[2]=fZ;
    cartParameters[3]=fPx;
    cartParameters[4]=fPy;
    cartParameters[5]=fPz;
    cartParameters[6]=fErrX;
    cartParameters[7]=fErrY;
    cartParameters[8]=fErrZ;
    cartParameters[9]=fErrPx;
    cartParameters[10]=fErrPy;
    cartParameters[11]=fErrPz;
}