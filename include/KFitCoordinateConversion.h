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
 * KFitCoordinateConversion.h
 *
 *@updated 17.10.2024
 *@version v1.0.0 
 *
 * Class for performing coordinate conversion 
 * between cartesian and spherical coordinates
 * 
 */

#ifndef KFITCOORDINATECONVERSION_H
#define KFITCOORDINATECONVERSION_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// root includes
#include <TMath.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>

// using std::cout;
// using std::endl;

using namespace std;

class KFitCoordinateConversion
{

    ClassDef(KFitCoordinateConversion, 2);

private:

    double fPi2 = TMath::PiOver2();
    double fMom = 9.9e6;
    double fTheta = 9.9e6;
    double fPhi = 9.9e6;
    double fR = 9.9e6;
    double fZ = 9.9e6;
    double fPx = 9.9e6;
    double fPy = 9.9e6;
    double fPz = 9.9e6;
    double fX = 9.9e6;
    double fY = 9.9e6;
    double fErrMom = 9.9e6;
    double fErrTheta = 9.9e6;
    double fErrPhi = 9.9e6;
    double fErrR = 9.9e6;
    double fErrZ = 9.9e6;
    double fErrPx = 9.9e6;
    double fErrPy = 9.9e6;
    double fErrPz = 9.9e6;
    double fErrX = 9.9e6;
    double fErrY = 9.9e6;

public:

    /** Default Constructor **/
    KFitCoordinateConversion();

    /** Default Destructor **/
    ~KFitCoordinateConversion() {};
     

    void setMom(double val, double errval)
    {
        fMom = val;
        fErrMom = errval;
    }
    void setTheta(double val, double errval)
    {
        fTheta = val;
        fErrTheta = errval;
    }
    void setPhi(double val, double errval)
    {
        fPhi = val;
        fErrPhi = errval;
    }
    void setR(double val, double errval)
    {
        fR = val;
        fErrR = errval;
    }

    void setPx(double val, double errval)
    {
        fPx = val;
        fErrPx = errval;
    }
    void setPy(double val, double errval)
    {
        fPy = val;
        fErrPy = errval;
    }
    void setPz(double val, double errval)
    {
        fPz = val;
        fErrPz = errval;
    }
    void setX(double val, double errval)
    {
        fX = val;
        fErrX = errval;
    }
    void setY(double val, double errval)
    {
        fY = val;
        fErrY = errval;
    }
    void setZ(float val, double errval)
    {
        fZ = val;
        fErrZ = errval;
    }

    void setParametersCart(double cartParameters[])
    {
        fX = cartParameters[0];
        fY = cartParameters[1];
        fZ = cartParameters[2];
        fPx = cartParameters[3];
        fPy = cartParameters[4];
        fPz = cartParameters[5];
    }
    void setParameterErrorsCart(double cartErrorParameters[])
    {
        fErrX = cartErrorParameters[0];
        fErrY = cartErrorParameters[1];
        fErrZ = cartErrorParameters[2];
        fErrPx = cartErrorParameters[3];
        fErrPy = cartErrorParameters[4];
        fErrPz = cartErrorParameters[5];
    }

    void setParametersSpher(double spherParameters[])
    {
        fMom = spherParameters[0];
        fTheta = spherParameters[1];
        fPhi = spherParameters[2];
        fR = spherParameters[3];
        fZ = spherParameters[4];
    }
    void setParameterErrorsSpher(double spherErrorParameters[])
    {
        fErrMom = spherErrorParameters[0];
        fErrTheta = spherErrorParameters[1];
        fErrPhi = spherErrorParameters[2];
        fErrR = spherErrorParameters[3];
        fErrZ = spherErrorParameters[4];
    }

    double getMom() const { return fMom; }
    double getTheta() const { return fTheta; }
    double getPhi() const { return fPhi; }
    double getR() const { return fR; }

    double getPx() const { return fPx; }
    double getPy() const { return fPy; }
    double getPz() const { return fPz; }
    double getX() const { return fX; }
    double getY() const { return fY; }
    double getZ() const { return fZ; }

    void setVarsCart()
    {
        fPx = (fMom * (sin(fTheta)) * (cos(fPhi)));
        fPy = (fMom * sin(fTheta) * sin(fPhi));
        fPz = (fMom * cos(fTheta));
        fX = (fR * cos((fPhi) + fPi2));
        fY = (fR * sin((fPhi) + fPi2));
    }

    void setVarsSpher()
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

    void setErrorsCart()
    {

        fErrPx = sqrt(pow((fMom * cos(fTheta) * cos(fPhi) * fErrTheta), 2) + pow((fMom * sin(fTheta) * sin(fPhi) * fErrPhi), 2) + pow((sin(fTheta) * cos(fPhi) * fErrMom), 2) + 2 * (fMom * cos(fTheta) * cos(fPhi)) * ((fMom * sin(fTheta) * sin(fPhi))) * (fErrTheta * fErrPhi) + 2 * (fMom * cos(fTheta) * cos(fPhi)) * (sin(fTheta) * cos(fPhi)) * (fErrTheta * fErrMom) + 2 * ((fMom * sin(fTheta) * sin(fPhi)) * ((sin(fTheta) * cos(fPhi)) * (fErrMom * fErrPhi))));
        fErrPy = sqrt(pow((fMom * cos(fTheta) * sin(fPhi) * fErrTheta), 2) + pow((fMom * sin(fTheta) * cos(fPhi) * fErrPhi), 2) + pow((sin(fTheta) * sin(fPhi) * fErrMom), 2) + 2 * (fMom * cos(fTheta) * sin(fPhi) * ((fMom * sin(fTheta) * cos(fPhi)))) * (fErrTheta * fErrPhi) + 2 * (fMom * cos(fTheta) * sin(fPhi)) * (sin(fTheta) * sin(fPhi)) * (fErrTheta * fErrMom) + 2 * (fMom * sin(fTheta) * cos(fPhi)) * (sin(fTheta) * sin(fPhi)) * (fErrPhi * fErrMom));
        fErrPz = sqrt(pow((fMom * sin(fTheta) * fErrTheta), 2) + pow((cos(fTheta) * fErrMom), 2) + 2 * (fMom * sin(fTheta)) * (cos(fTheta) * (fErrTheta * fErrMom)));
        fErrX = sqrt((pow((fR * sin(fPhi + fPi2) * fErrPhi), 2) + pow((cos(fPhi + fPi2) * fErrR), 2)) + 2 * (fR * sin(fPhi + fPi2)) * (cos(fPhi + fPi2)) * (fErrPhi * fErrR));
        fErrY = sqrt((pow((fR * cos(fPhi + fPi2) * fErrPhi), 2) + pow((sin(fPhi + fPi2) * fErrR), 2)) + 2 * (fR * cos(fPhi + fPi2)) * (sin(fPhi + fPi2)) * (fErrPhi * fErrR));
    }

    void setErrorsSpher()
    {

        fErrMom = sqrt((pow((fPx / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2)))), 2) * (pow(fErrPx, 2))) + (pow((fPy / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2)))), 2) * (pow(fErrPy, 2))) + (pow((fPz / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2)))), 2) * (pow(fErrPz, 2))) + 2 * ((fPx / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * ((fPy / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * (fErrPx * fErrPy) + 2 * ((fPx / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * ((fPz / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * (fErrPx * fErrPz) + 2 * ((fPz / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * ((fPy / (sqrt(pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2))))) * (fErrPy * fErrPz));
        fErrTheta = sqrt(pow(((fPz * fPx) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2)))), 2) * (pow(fErrPx, 2)) + pow(((fPz * fPy) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2)))), 2) * (pow(fErrPy, 2)) + pow(((sqrt(((pow(fPx, 2)) + (pow(fPy, 2))))) / (((pow(fPx, 2)) + (pow(fPy, 2))) + (pow(fPz, 2)))), 2) * (pow(fErrPz, 2)) + 2 * (((fPz * fPx) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2))))) * (((fPz * fPy) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2))))) * (fErrPx * fErrPy) + 2 * (((fPz * fPx) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2))))) * ((sqrt(((pow(fPx, 2)) + (pow(fPy, 2))))) / (((pow(fPx, 2)) + (pow(fPy, 2))) + (pow(fPz, 2)))) * (fErrPx * fErrPz) + 2 * (((fPz * fPy) / (sqrt(pow(fPx, 2) + pow(fPy, 2)) * (pow(fPx, 2) + pow(fPz, 2) + pow(fPy, 2))))) * (((sqrt(((pow(fPx, 2)) + (pow(fPy, 2))))) / (((pow(fPx, 2)) + (pow(fPy, 2))) + (pow(fPz, 2))))) * (fErrPz * fErrPy)); // with correlation
        fErrPhi = sqrt((pow((fPy / ((pow(fPx, 2)) + (pow(fPy, 2)))), 2)) * (pow(fErrPx, 2)) + (pow((fPx / ((pow(fPx, 2)) + (pow(fPy, 2)))), 2) * (pow(fErrPy, 2))) + 2 * (fPy / ((pow(fPx, 2)) + (pow(fPy, 2)))) * (fPx / ((pow(fPx, 2)) + (pow(fPy, 2)))) * (fErrPx * fErrPy));

        double dRdX = (((sin(fPhi)) * (sin(fTheta))) / (sqrt(pow((sin(fPhi) * sin(fTheta)), 2) + pow((cos(fPhi) * sin(fTheta)), 2))));
        double dRdY = (-(cos(fPhi) * sin(fTheta)) / (sqrt(pow((sin(fPhi) * sin(fTheta)), 2) + pow((cos(fPhi) * sin(fTheta)), 2))));
        double dRdphi = ((sin(fTheta)) * (fY * sin(fPhi) + fX * cos(fPhi))) / (sqrt(pow((sin(fPhi) * sin(fTheta)), 2) + pow((cos(fPhi) * sin(fTheta)), 2)));
        double fErrR = (sqrt(pow((dRdX * fErrX), 2) + pow((dRdY * fErrY), 2) + pow((dRdphi * fErrPhi), 2) + 2 * (dRdX) * (dRdY) * (fErrX * fErrY) + 2 * (dRdX) * (dRdphi) * (fErrX * fErrPhi) + 2 * (dRdY) * (dRdphi) * (fErrY * fErrPhi)));
    }
};
#endif /* KFITCOORDINATECONVERSION_H */
