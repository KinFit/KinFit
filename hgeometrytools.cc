#include "hgeometrytools.h"

// HGeometryTools::HGeometryTools(void)  {
//    // The default constructor
// }

// HGeometryTools::~HGeometryTools(void) {
//    // Everything that is constructed has to be destructed
// }

void HGeometryTools::addLine(const HGeomVector &r, const HGeomVector &alpha, const Double_t w)
{

    // Function to add lines to the fit.
    // Input
    //   r --> A point in the straight line
    //   alpha --> normalized direction vector
    //   w --> weight of this line in the fit

    fM(0, 0) = w * (alpha.getY() * alpha.getY() + alpha.getZ() * alpha.getZ());
    fM(0, 1) = -w * (alpha.getX() * alpha.getY());
    fM(0, 2) = -w * (alpha.getX() * alpha.getZ());
    //  fM(1,0)=-(alpha.X() * alpha.Y());
    fM(1, 1) = w * (alpha.getX() * alpha.getX() + alpha.getZ() * alpha.getZ());
    fM(1, 2) = -w * (alpha.getY() * alpha.getZ());
    // fM(2,0)=-(alpha.X() * alpha.Z());
    // fM(2,1)=-(alpha.Y() * alpha.Z());
    fM(2, 2) = w * (alpha.getY() * alpha.getY() + alpha.getX() * alpha.getX());

    fSys(0, 0) += fM(0, 0);
    fSys(0, 1) += fM(0, 1);
    fSys(0, 2) += fM(0, 2);
    fSys(1, 1) += fM(1, 1);
    fSys(1, 2) += fM(1, 2);
    fSys(2, 2) += fM(2, 2);

    fB.X() += fM(0, 0) * r.getX() + fM(0, 1) * r.getY() + fM(0, 2) * r.getZ();
    fB.Y() += fM(0, 1) * r.getX() + fM(1, 1) * r.getY() + fM(1, 2) * r.getZ();
    fB.Z() += fM(0, 2) * r.getX() + fM(1, 2) * r.getY() + fM(2, 2) * r.getZ();
}

void HGeometryTools::getVertex(HGeomVector &out)
{
    // This method fills the vector "out" with the coordinates
    // of the point of maximun approach to the lines added with
    // addLine()
    Double_t det = 0;

    det = fSys(0, 0) * fSys(1, 1) * fSys(2, 2) + fSys(0, 1) * fSys(1, 2) * fSys(0, 2) + fSys(0, 1) * fSys(1, 2) * fSys(0, 2) - fSys(0, 2) * fSys(1, 1) * fSys(0, 2) - fSys(0, 1) * fSys(0, 1) * fSys(2, 2) - fSys(1, 2) * fSys(1, 2) * fSys(0, 0);

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
        out.setXYZ(-1000., -1000., -1000.);
        return;
    }
}
