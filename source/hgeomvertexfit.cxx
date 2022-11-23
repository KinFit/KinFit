#include "hgeomvertexfit.h"

//_KinFit_CLASS_DESCRIPTION
////////////////////////////////////////////////
// HGeomVertexFit
//
//   Calculates the point of maximun approach to any given
// set of n tracks with user specified weights.
//
//   To start the operation the function reset() is called, then
// addLine() should be used on all tracks on the sample and
// finally GetVertex() will provide the desired point.
//
//   Note that GetVertex() is not destructive, you can continue
// adding lines after GetVertex() was called
////////////////////////////////////////////////

HGeomVertexFit::HGeomVertexFit(void)
{
    // The default constructor
    fM.ResizeTo(3, 3);
    fSys.ResizeTo(3, 3);
}

HGeomVertexFit::~HGeomVertexFit(void)
{
    // Everything that is constructed has to be destructed
}

void HGeomVertexFit::addLine(const TVector3 &r, const TVector3 &alpha,
                             const Double_t w)
{
    // Function to add lines to the fit.
    // Input
    //   r --> A point in the straight line
    //   alpha --> normalized direction vector
    //   w --> weight of this line in the fit

    fM(0, 0) = w * (alpha.Y() * alpha.Y() + alpha.Z() * alpha.Z());
    fM(0, 1) = -w * (alpha.X() * alpha.Y());
    fM(0, 2) = -w * (alpha.X() * alpha.Z());
    //  fM(1,0)=-(alpha.X() * alpha.Y());
    fM(1, 1) = w * (alpha.X() * alpha.X() + alpha.Z() * alpha.Z());
    fM(1, 2) = -w * (alpha.Y() * alpha.Z());
    // fM(2,0)=-(alpha.X() * alpha.Z());
    // fM(2,1)=-(alpha.Y() * alpha.Z());
    fM(2, 2) = w * (alpha.Y() * alpha.Y() + alpha.X() * alpha.X());

    fSys(0, 0) += fM(0, 0);
    fSys(0, 1) += fM(0, 1);
    fSys(0, 2) += fM(0, 2);
    fSys(1, 1) += fM(1, 1);
    fSys(1, 2) += fM(1, 2);
    fSys(2, 2) += fM(2, 2);

    fB.X() += fM(0, 0) * r.X() + fM(0, 1) * r.Y() + fM(0, 2) * r.Z();
    fB.Y() += fM(0, 1) * r.X() + fM(1, 1) * r.Y() + fM(1, 2) * r.Z();
    fB.Z() += fM(0, 2) * r.X() + fM(1, 2) * r.Y() + fM(2, 2) * r.Z();
}

void HGeomVertexFit::getVertex(TVector3 &out);
{
    // This method fills the vector "out" with the coordinates
    // of the point of maximun approach to the lines added with
    // addLine()
    Double_t det = 0;

    det = fSys.Determinant();

    fM(0, 0) = fSys(1, 1) * fSys(2, 2) - fSys(1, 2) * fSys(1, 2);
    fM(0, 1) = fSys(1, 2) * fSys(0, 2) - fSys(0, 1) * fSys(2, 2);
    fM(0, 2) = fSys(0, 1) * fSys(1, 2) - fSys(1, 1) * fSys(0, 2);
    fM(1, 1) = fSys(0, 0) * fSys(2, 2) - fSys(0, 2) * fSys(0, 2);
    fM(1, 2) = fSys(0, 1) * fSys(0, 2) - fSys(0, 0) * fSys(1, 2);
    fM(2, 2) = fSys(0, 0) * fSys(1, 1) - fSys(0, 1) * fSys(0, 1);
    fM(1, 0) = fM(0, 1);
    fM(2, 0) = fM(0, 2);
    fM(2, 1) = fM(1, 2);

    if (det == 0)
    {
        out.SetXYZ(-1000., -1000., -1000.);
        return;
    }

    fM *= (1. / det);
    out = fM * fB;
}

void HGeomVertexFit::reset(void)
{
    // Resets the fitting procedure. That is, the class goes to the
    // state where no line was added
    for (Int_t i = 0; i < 3; i++)
    {
        for (Int_t j = 0; j < 3; j++)
            fM(i, j) = fSys(i, j) = 0.0;
    }
    fB.SetXYZ(0.0, 0.0, 0.0);
}