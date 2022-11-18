// This code was borrowed completely
// from Hydra software (https://hades.gsi.de/?q=computing)
// HGeomVector was replace by TVector3

#ifndef HPARTICLETOOL_H
#define HPARTICLETOOL_H

// ROOT includes
#include "TVector3.h"
#include "TMath.h"

// system includes
#include <math.h>

// Create namespace HParticleTool for compatibility with Hydra
namespace HParticleTool
{
  void calcSegVector(Double_t z, Double_t rho, Double_t phi, Double_t theta,
                     TVector3 &base, TVector3 &dir);

  Double_t calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,
                                                   TVector3 &point);

  Double_t calculateMinimumDistance(TVector3 &base1, TVector3 &dir1,
                                    TVector3 &base2, TVector3 &dir2);

  TVector3 calculateCrossPoint(TVector3 &base1, TVector3 &dir1, TVector3 &base2,
                               TVector3 &dir2);

  Double_t calcDeterminant(TVector3 &v1, TVector3 &v2, TVector3 &v3);

  TVector3 calculatePointOfClosestApproach(TVector3 &base1, TVector3 &dir1,
                                           TVector3 &base2, TVector3 &dir2);

  TVector3 calcVertexAnalytical(TVector3 &base1, TVector3 &dir1, TVector3 &base2,
                                TVector3 &dir2);
} // end of the namespace

#endif /* HPARTICLETOOL_H */