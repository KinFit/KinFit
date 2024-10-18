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

class KFitCoordinateConversion : public TObject
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

    void setVarsCart();

    void setVarsSpher();

    void setErrorsCart();

    void setErrorsSpher();

    void convertToSpherical(double arr[11], double arr2[9]);
    
    void convertToCart(double arr[9], double arr2[11]);
    
};
#endif /* KFITCOORDINATECONVERSION_H */
