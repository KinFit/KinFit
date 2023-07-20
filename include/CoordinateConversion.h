/**
 * CoordinateConversion.h
 *
 *
 */

#ifndef COORDINATECONVERSION_H
#define COORDINATECONVERSION_H

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

//using std::cout;
//using std::endl;

using namespace std;

class CoordinateConversion{

public:

double mom = (9*10^(-999)), theta = (9*10^(-999)), phi = (9*10^(-999)), R = (9*10^(-999)), Z = (9*10^(-999)), Px = (9*10^(-999)), Py = (9*10^(-999)), Pz = (9*10^(-999)), X = (9*10^(-999)), Y = (9*10^(-999)), errmom = (9*10^(-999)), errtheta = (9*10^(-999)), errphi = (9*10^(-999)), errR = (9*10^(-999)), errZ = (9*10^(-999)), errPx = (9*10^(-999)), errPy = (9*10^(-999)), errPz = (9*10^(-999)), errX = (9*10^(-999)), errY = (9*10^(-999)), Rtest, Rtemp;

    void setmom(double val, double errval){mom = val; errmom = errval;}
    void settheta(double val, double errval){theta = val; errtheta = errval;}
    void setphi(double val, double errval){phi = val; errphi = errval;}
    void setR(double val, double errval){R = val; errR = errval;}

    void setPx(double val, double errval){Px = val; errPx = errval;}
    void setPy(double val, double errval){Py = val; errPy = errval;}
    void setPz(double val, double errval){Pz = val; errPz = errval;}
    void setX(double val, double errval){X = val; errX = errval;}
    void setY(double val, double errval){Y = val; errY = errval;}
    void setZ(float val, double errval){Z = val; errZ = errval;}

    void setparametersCart(double CartParameters[]){X = CartParameters[0]; Y = CartParameters[1]; Z = CartParameters[2]; Px = CartParameters[3]; Py = CartParameters[4]; Pz  = CartParameters[5];}
    void setparametererrorsCart(double CartErrorParameters[]){errX = CartErrorParameters[0]; errY = CartErrorParameters[1]; errZ = CartErrorParameters[2]; errPx = CartErrorParameters[3]; errPy = CartErrorParameters[4]; errPz = CartErrorParameters[5];
    }

    void setparametersSpher(double SpherParameters[]){mom = SpherParameters[0]; theta = SpherParameters[1]; phi = SpherParameters[2]; R = SpherParameters[3]; Z = SpherParameters[4];}
    void setparametererrorsSpher(double SpherErrorParameters[]){errmom = SpherErrorParameters[0]; errtheta = SpherErrorParameters[1]; errphi = SpherErrorParameters[2]; errR = SpherErrorParameters[3]; errZ = SpherErrorParameters[4];}

    double getmom() const{return mom;}
    double gettheta() const{return theta;}
    double getphi() const{return phi;}
    double getR() const{return R;}

    double getPx() const{return Px;}
    double getPy() const{return Py;}
    double getPz() const{return Pz;}
    double getX() const{return X;}
    double getY() const{return Y;}
    float getZ() const{return Z;}

    void setvarsCart(){
        Px = (mom*(sin(theta))*(cos(phi)));
        Py = (mom*sin(theta)*sin(phi));
        Pz = (mom*cos(theta));
        X = (R*cos((phi)+(3.14159265/2)));
        Y = (R*sin((phi)+(3.14159265/2)));
    }

    void seterrorsCart(){
        errPx = sqrt(pow((mom*cos(theta)*cos(phi)*errtheta), 2)+pow((mom*sin(theta)*sin(phi)*errphi), 2)+pow((sin(theta)*cos(phi)*errmom), 2)+2*(mom*cos(theta)*cos(phi))*((mom*sin(theta)*sin(phi)))*(errtheta*errphi)+2*(mom*cos(theta)*cos(phi))*(sin(theta)*cos(phi))*(errtheta*errmom)+2*((mom*sin(theta)*sin(phi))*((sin(theta)*cos(phi))*(errmom*errphi))));
//        errPx = sqrt(pow((mom*cos(theta)*cos(phi)*errtheta), 2)+pow((mom*sin(theta)*sin(phi)*errphi), 2)+pow((sin(theta)*cos(phi)*errmom), 2));
        errPy = sqrt(pow((mom*cos(theta)*sin(phi)*errtheta), 2)+pow((mom*sin(theta)*cos(phi)*errphi), 2)+pow((sin(theta)*sin(phi)*errmom), 2)+2*(mom*cos(theta)*sin(phi)*((mom*sin(theta)*cos(phi))))*(errtheta*errphi)+2*(mom*cos(theta)*sin(phi))*(sin(theta)*sin(phi))*(errtheta*errmom)+2*(mom*sin(theta)*cos(phi))*(sin(theta)*sin(phi))*(errphi*errmom));
//        errPy = sqrt(pow((mom*cos(theta)*sin(phi)*errtheta), 2)+pow((mom*sin(theta)*cos(phi)*errphi), 2)+pow((sin(theta)*sin(phi)*errmom), 2));
        errPz = sqrt(pow((mom*sin(theta)*errtheta), 2)+pow((cos(theta)*errmom), 2)+2*(mom*sin(theta))*(cos(theta)*(errtheta*errmom)));
//        errPz = sqrt(pow((mom*sin(theta)*errtheta), 2)+pow((cos(theta)*errmom), 2));
        errX = sqrt((pow((R*sin(phi+(3.14159265/2))*errphi), 2)+pow((cos(phi+(3.14159265/2))*errR), 2))+2*(R*sin(phi+(3.14159265/2)))*(cos(phi+(3.14159265/2)))*(errphi*errR));
//        errX = sqrt((pow((R*sin(phi+(3.14159265/2))*errphi), 2)+pow((cos(phi+(3.14159265/2))*errR), 2)));
        errY = sqrt((pow((R*cos(phi+(3.14159265/2))*errphi), 2)+pow((sin(phi+(3.14159265/2))*errR), 2))+2*(R*cos(phi+(3.14159265/2)))*(sin(phi+(3.14159265/2)))*(errphi*errR));
//        errY = sqrt((pow((R*cos(phi+(3.14159265/2))*errphi), 2)+pow((sin(phi+(3.14159265/2))*errR), 2)));

    ;}

    void setvarsSpher(){

        mom = sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2));
        theta = atan((sqrt((pow(Px, 2))+(pow(Py, 2))))/Pz);
        phi = atan2(Py,Px);
//        phi = acos(Px/(sqrt((pow(Px, 2))+(pow(Py, 2)))));
        //phi = asin(Py/(sqrt((pow(Px, 2))+(pow(Py, 2)))));

        TVector3 trackdist;
        TVector3 trackbase;
        TVector3 beamdist;
        TVector3 beambase;

        trackdist.SetXYZ((sin(theta))*(cos(phi)), (sin(theta))*(sin(phi)), (cos(theta)));
        trackbase.SetXYZ(X, Y, Z);
        beamdist.SetXYZ(0, 0, 1);
        beambase.SetXYZ(0, 0, 1);

        TVector3 cross = trackdist.Cross(beamdist);
        TVector3 diff = trackbase-beambase;

        R = -(((cross.Dot(diff)))/(cross.Mag()));

                //R and Y - same sign
/*
if(Px > 0){
    R = ((abs((cross.Dot(diff))))/(cross.Mag()));
}
else if(Px < 0){
    R = -((abs((cross.Dot(diff))))/(cross.Mag()));
}
*/

//        R = ((abs((cross.Dot(diff))))/(cross.Mag()));
//        Rtest = (abs(((X*sin(theta)*sin(phi))-(Y*sin(theta)*cos(phi)))))/(sqrt(pow((sin(theta)*sin(phi)), 2)+pow((sin(theta)*cos(phi)), 2)));
        Rtest = (((X*sin(theta)*sin(phi))-(Y*sin(theta)*cos(phi))))/(sqrt(pow((sin(theta)*sin(phi)), 2)+pow((sin(theta)*cos(phi)), 2)));
/*
if(Y < 0){
    if(Rtemp > 0){
    R = -Rtest;}
    else if(Rtemp > 0){
        R = Rtest;
    }
}
else if(Y > 0){
    if (Rtemp < 0){
        R = -Rtest;}
    else if (Rtemp > 0){
        R = Rtest;
    }
}
*/
    ;}

    void seterrorsSpher(){

        errmom = sqrt((pow((Px/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))), 2)*(pow(errPx, 2)))+(pow((Py/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))), 2)*(pow(errPy, 2)))+(pow((Pz/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))), 2)*(pow(errPz, 2)))+2*((Px/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))))*((Py/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))))*(errPx*errPy)+2*((Px/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))))*((Pz/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))))*(errPx*errPz)+2*((Pz/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))))*((Py/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))))*(errPy*errPz));
//        errmom = sqrt((pow((Px/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))), 2)*(pow(errPx, 2)))+(pow((Py/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))), 2)*(pow(errPy, 2)))+(pow((Pz/(sqrt(pow(Px, 2)+pow(Py, 2)+pow(Pz, 2)))), 2)*(pow(errPz, 2))));
//        errtheta = sqrt(pow(((Pz*Px)/(sqrt(pow(Px, 2)+pow(Py, 2))*(pow(Px, 2)+pow(Pz, 2)+pow(Py, 2)))), 2)*(pow(errPx, 2))+pow(((Pz*Py)/(sqrt(pow(Px, 2)+pow(Py, 2))*(pow(Px, 2)+pow(Pz, 2)+pow(Py, 2)))), 2)*(pow(errPy, 2))+pow(((sqrt(((pow(Px, 2))+(pow(Py, 2)))))/(((pow(Px, 2))+(pow(Py, 2)))+(pow(Pz, 2)))), 2)*(pow(errPz, 2)));
        errtheta = sqrt(pow(((Pz*Px)/(sqrt(pow(Px, 2)+pow(Py, 2))*(pow(Px, 2)+pow(Pz, 2)+pow(Py, 2)))), 2)*(pow(errPx, 2))+pow(((Pz*Py)/(sqrt(pow(Px, 2)+pow(Py, 2))*(pow(Px, 2)+pow(Pz, 2)+pow(Py, 2)))), 2)*(pow(errPy, 2))+pow(((sqrt(((pow(Px, 2))+(pow(Py, 2)))))/(((pow(Px, 2))+(pow(Py, 2)))+(pow(Pz, 2)))), 2)*(pow(errPz, 2))+2*(((Pz*Px)/(sqrt(pow(Px, 2)+pow(Py, 2))*(pow(Px, 2)+pow(Pz, 2)+pow(Py, 2)))))*(((Pz*Py)/(sqrt(pow(Px, 2)+pow(Py, 2))*(pow(Px, 2)+pow(Pz, 2)+pow(Py, 2)))))*(errPx*errPy)+2*(((Pz*Px)/(sqrt(pow(Px, 2)+pow(Py, 2))*(pow(Px, 2)+pow(Pz, 2)+pow(Py, 2)))))*((sqrt(((pow(Px, 2))+(pow(Py, 2)))))/(((pow(Px, 2))+(pow(Py, 2)))+(pow(Pz, 2))))*(errPx*errPz)+2*(((Pz*Py)/(sqrt(pow(Px, 2)+pow(Py, 2))*(pow(Px, 2)+pow(Pz, 2)+pow(Py, 2)))))*(((sqrt(((pow(Px, 2))+(pow(Py, 2)))))/(((pow(Px, 2))+(pow(Py, 2)))+(pow(Pz, 2)))))*(errPz*errPy)); // with correlation
        errphi = sqrt((pow((Py/((pow(Px, 2))+(pow(Py, 2)))), 2))*(pow(errPx, 2))+(pow((Px/((pow(Px, 2))+(pow(Py, 2)))), 2)*(pow(errPy, 2)))+2*(Py/((pow(Px, 2))+(pow(Py, 2))))*(Px/((pow(Px, 2))+(pow(Py, 2))))*(errPx*errPy)); // with correlation
//        errphi = sqrt((pow((Py/((pow(Px, 2))+(pow(Py, 2)))), 2))*(pow(errPx, 2))+(pow((Px/((pow(Px, 2))+(pow(Py, 2)))), 2)*(pow(errPy, 2))));


//        errphi = sqrt((pow((Py/((pow(Px, 2))+(pow(Py, 2)))), 2))*(pow(errPx, 2))+(pow((Px/((pow(Px, 2))+(pow(Py, 2)))), 2)*(pow(errPy, 2)))+pow((errPx/Px), 2)+pow((errPy/Py), 2)+2*(errPx/Px)*(errPy/Py)*(-0.000752836/(errPx*errPy)));
// wrong        double dRdX = ((pow((sin(phi)*sin(theta)), 2)*(X-Y))/((sqrt((pow((sin(phi)*sin(theta)), 2)+pow((cos(phi)*sin(theta)), 2))))*abs((sin(phi)*sin(theta))*(X-Y))));
// 2nd it.        double dRdX = (((sin(phi)*sin(theta))*(sin(phi)*sin(theta)*X-cos(phi)*sin(theta)*Y))/((sqrt(pow((sin(phi)*sin(theta)), 2)+pow((cos(phi)*sin(theta)), 2)))*abs((sin(phi)*sin(theta)*X)-(cos(phi)*sin(theta)*Y))));
double dRdX = (((sin(phi))*(sin(theta)))/(sqrt(pow((sin(phi)*sin(theta)), 2)+pow((cos(phi)*sin(theta)), 2))));
// wrong       double dRdY = ((pow((sin(phi)*sin(theta)), 2)*(Y-X))/((sqrt((pow((sin(phi)*sin(theta)), 2)+pow((cos(phi)*sin(theta)), 2))))*abs((sin(phi)*sin(theta))*(Y-X))));
// 2nd it.        double dRdY = (((cos(phi)*sin(theta))*(cos(phi)*sin(theta)*Y-sin(phi)*sin(theta)*X))/((sqrt(pow((sin(phi)*sin(theta)), 2)+pow((cos(phi)*sin(theta)), 2)))*abs((cos(phi)*sin(theta)*Y)-(sin(phi)*sin(theta)*X))));
double dRdY = (-(cos(phi)*sin(theta))/(sqrt(pow((sin(phi)*sin(theta)), 2)+pow((cos(phi)*sin(theta)), 2))));
// 2nd it.        double dRdphi = ((pow((sin(theta)), 2))*(X*sin(phi)-Y*cos(phi))*(Y*sin(phi)+X*cos(phi)))/((sqrt(pow((sin(phi)*sin(theta)), 2)+pow((cos(phi)*sin(theta)), 2)))*(abs((sin(phi)*sin(theta)*X)-(cos(phi)*sin(theta)*Y))));
// wrong       double dRdphi = (pow((sin(theta)*(Y-X)), 2)*cos(phi)*sin(phi))/((sqrt(pow((sin(theta)*sin(phi)), 2)+pow((sin(theta)*cos(phi)), 2)))*abs((sin(theta)*sin(phi)*(Y-X))));
double dRdphi = ((sin(theta))*(Y*sin(phi)+X*cos(phi)))/(sqrt(pow((sin(phi)*sin(theta)), 2)+pow((cos(phi)*sin(theta)), 2)));

        errR = (sqrt(pow((dRdX*errX), 2)+pow((dRdY*errY), 2)+pow((dRdphi*errphi), 2)+2*(dRdX)*(dRdY)*(errX*errY)+2*(dRdX)*(dRdphi)*(errX*errphi)+2*(dRdY)*(dRdphi)*(errY*errphi)));
//        errR = (sqrt(pow((dRdX*errX), 2)+pow((dRdY*errY), 2)+pow((dRdphi*errphi), 2)));


    ;}

/*
    void setvectors(){

TVector3 trackdist;
TVector3 trackbase;
TVector3 beamdist;
TVector3 beambase;

trackdist.SetXYZ((sin(theta))*(cos(phi)), (sin(theta))*(sin(phi)), (cos(theta)));
trackbase.SetXYZ(X, Y, Z);
beamdist.SetXYZ(0, 0, 1);
beambase.SetXYZ(0, 0, 1);

TVector3 cross = trackdist.Cross(beamdist);
TVector3 diff = trackbase-beambase;

R = ((abs((cross.Dot(diff))))/(cross.Mag()));

    ;}
*/


};
#endif
