//*-- Author : M. Sanchez
//*-- Modified : 13.03.2001 (M. Sanchez)
//*-- Modified : 04.08.2022 (W. Esmail)
// This code was borrowed completely
// from Hydra software (https://hades.gsi.de/?q=computing)
// HGeomVector is replace by TVector3
// HGeomMatrix is replace by TMatrixD
#ifndef HGEOMVERTEXFIT_H
#define HGEOMVERTEXFIT_H

#include "TMatrixD.h"
#include "TVector3.h"

class HGeomVertexFit
{
private:
    TMatrixD fM; // Temporal matrix for calculations
protected:
    TMatrixD fSys; // LSM system inverse matrix
    TVector3 fB;   // LSM independent term
public:
    HGeomVertexFit(void);
    ~HGeomVertexFit(void);
    void addLine(const TVector3 &r, const TVector3 &alpha,
                 const Double_t w = 1.0);
    void getVertex(TVector3 &out);
    void reset(void);
};
#endif