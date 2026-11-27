//****************************************************************************
//*                    This file is part of KinFit.                          *
//*                                                                          *
//*             KinFit is distributed under the terms of the                 *
//*              GNU General Public License (GPL) version 3,                 *
//*                 copied verbatim in the file "LICENSE".                   *
//*                                                                          *
//*                             Copyright 2024                               *
//*             GSI Helmholtzzentrum für Schwerionenforschung                *
//*          This software is distributed under the terms of the             *
//*          GNU General Public Licence version 3 (GPL Version 3)            *
//*                                                                          *
//*     The copyright holders are listed in the file "COPYRIGHTHOLDERS".     *
//*               The authors are listed in the file "AUTHORS".              *
//****************************************************************************

/**
 * CoordinateConversion.h
 *
 * Utility class to convert track parameters between spherical
 * (p, theta, phi, R, Z) and Cartesian (X, Y, Z, Px, Py, Pz) forms,
 * including first-order uncertainty propagation.
 */

#pragma once

// system includes
#include <cmath>

// root includes
#include <TMath.h>

/**
 * Track parameters: spherical (p, theta, phi, R, Z) with the point of
 * closest approach (PCA) to the z axis at (-R sin(phi), R cos(phi), Z), and
 * Cartesian (X, Y, Z, Px, Py, Pz) with (X, Y, Z) a point on the track
 * (the PCA after setvarsCart()).
 *
 * Errors are propagated to first order assuming *uncorrelated* inputs (only
 * the uncertainties are known to this class). The old formulas added
 * cross terms 2 (df/da)(df/db) sigma_a sigma_b, i.e. assumed full
 * correlation.
 */
class CoordinateConversion {

public:
  // (The old defaults were 9 * 10 ^ (-999), an integer XOR, not 9e-999.)
  double mom = 0, theta = 0, phi = 0, R = 0, Z = 0, Px = 0, Py = 0, Pz = 0,
         X = 0, Y = 0, errmom = 0, errtheta = 0, errphi = 0, errR = 0,
         errZ = 0, errPx = 0, errPy = 0, errPz = 0, errX = 0, errY = 0;

  void setmom(double val, double errval) {
    mom = val;
    errmom = errval;
  }
  void settheta(double val, double errval) {
    theta = val;
    errtheta = errval;
  }
  void setphi(double val, double errval) {
    phi = val;
    errphi = errval;
  }
  void setR(double val, double errval) {
    R = val;
    errR = errval;
  }

  void setPx(double val, double errval) {
    Px = val;
    errPx = errval;
  }
  void setPy(double val, double errval) {
    Py = val;
    errPy = errval;
  }
  void setPz(double val, double errval) {
    Pz = val;
    errPz = errval;
  }
  void setX(double val, double errval) {
    X = val;
    errX = errval;
  }
  void setY(double val, double errval) {
    Y = val;
    errY = errval;
  }
  void setZ(double val, double errval) { // was float
    Z = val;
    errZ = errval;
  }

  void setparametersCart(const double CartParameters[]) {
    X = CartParameters[0];
    Y = CartParameters[1];
    Z = CartParameters[2];
    Px = CartParameters[3];
    Py = CartParameters[4];
    Pz = CartParameters[5];
  }
  void setparametererrorsCart(const double CartErrorParameters[]) {
    errX = CartErrorParameters[0];
    errY = CartErrorParameters[1];
    errZ = CartErrorParameters[2];
    errPx = CartErrorParameters[3];
    errPy = CartErrorParameters[4];
    errPz = CartErrorParameters[5];
  }

  void setparametersSpher(const double SpherParameters[]) {
    mom = SpherParameters[0];
    theta = SpherParameters[1];
    phi = SpherParameters[2];
    R = SpherParameters[3];
    Z = SpherParameters[4];
  }
  void setparametererrorsSpher(const double SpherErrorParameters[]) {
    errmom = SpherErrorParameters[0];
    errtheta = SpherErrorParameters[1];
    errphi = SpherErrorParameters[2];
    errR = SpherErrorParameters[3];
    errZ = SpherErrorParameters[4];
  }

  double getmom() const { return mom; }
  double gettheta() const { return theta; }
  double getphi() const { return phi; }
  double getR() const { return R; }

  double getPx() const { return Px; }
  double getPy() const { return Py; }
  double getPz() const { return Pz; }
  double getX() const { return X; }
  double getY() const { return Y; }
  double getZ() const { return Z; } // was float

  /** Spherical -> Cartesian: momentum components and the PCA. */
  void setvarsCart() {
    const double st = std::sin(theta), ct = std::cos(theta);
    const double sp = std::sin(phi), cp = std::cos(phi);
    Px = mom * st * cp;
    Py = mom * st * sp;
    Pz = mom * ct;
    X = -R * sp;
    Y = R * cp;
  }

  /** Errors of setvarsCart() outputs from errmom, errtheta, errphi, errR. */
  void seterrorsCart() {
    const double st = std::sin(theta), ct = std::cos(theta);
    const double sp = std::sin(phi), cp = std::cos(phi);
    errPx = quad(st * cp * errmom, mom * ct * cp * errtheta,
                 mom * st * sp * errphi);
    errPy = quad(st * sp * errmom, mom * ct * sp * errtheta,
                 mom * st * cp * errphi);
    errPz = quad(ct * errmom, mom * st * errtheta, 0);
    errX = quad(sp * errR, R * cp * errphi, 0);
    errY = quad(cp * errR, R * sp * errphi, 0);
  }

  /** Cartesian -> spherical: (p, theta, phi) from the momentum, (R, Z) of
   * the PCA of the line through (X, Y, Z). */
  void setvarsSpher() {
    const double pT2 = Px * Px + Py * Py;
    mom = std::sqrt(pT2 + Pz * Pz);
    theta = std::atan2(std::sqrt(pT2), Pz); // was atan(pT/pz): < 0 backwards
    phi = std::atan2(Py, Px);
    const double sp = std::sin(phi), cp = std::cos(phi);
    R = -X * sp + Y * cp;
    if (pT2 > 0) // move (X, Y, Z) along the track to the PCA
      Z += -(X * Px + Y * Py) / pT2 * Pz;
  }

  /** Errors of setvarsSpher() outputs from errPx, errPy, errPz, errX, errY.
   * errZ is not propagated. */
  void seterrorsSpher() {
    const double pT2 = Px * Px + Py * Py;
    const double pT = std::sqrt(pT2);
    const double p2 = pT2 + Pz * Pz;
    const double p = std::sqrt(p2);
    errmom = quad(Px / p * errPx, Py / p * errPy, Pz / p * errPz);
    errtheta = quad(Px * Pz / (pT * p2) * errPx, Py * Pz / (pT * p2) * errPy,
                    pT / p2 * errPz);
    errphi = quad(Py / pT2 * errPx, Px / pT2 * errPy, 0);
    // R = -X sin(phi) + Y cos(phi); phi depends only on the momentum.
    const double sp = std::sin(phi), cp = std::cos(phi);
    const double dRdphi = -X * cp - Y * sp;
    errR = std::sqrt(sq(sp * errX) + sq(cp * errY) + sq(dRdphi * errphi));
  }

private:
  static double sq(double v) { return v * v; }
  static double quad(double a, double b, double c) {
    return std::sqrt(a * a + b * b + c * c);
  }
};
