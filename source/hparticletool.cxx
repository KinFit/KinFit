// framework includes
#include "hparticletool.h"

void HParticleTool::calcSegVector(Double_t z, Double_t rho, Double_t phi, Double_t theta,
                                  TVector3 &base, TVector3 &dir)
{
    base.SetXYZ(rho * std::cos(phi + TMath::PiOver2()),
                rho * std::sin(phi + TMath::PiOver2()), z);
    dir.SetXYZ(sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
}
Double_t HParticleTool::calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,
                                                                TVector3 &point)
{
    // calculates the minimum distance of a point to a straight
    // given as parametric straight x = base + n * dir

    if (!(dir.Mag() > 0))
    {
        return -1000000.;
    }

    TVector3 diff = base - point;

    TVector3 cross = dir.Cross(diff);
    return cross.Mag() / dir.Mag();
}

Double_t HParticleTool::calculateMinimumDistance(TVector3 &base1, TVector3 &dir1,
                                                 TVector3 &base2, TVector3 &dir2)
{
    // calculates the minimum distance of two tracks given
    // as parametric straights x = base + n * dir

    TVector3 cross = dir1.Cross(dir2);

    TVector3 ab = base1 - base2;

    if (!(std::fabs(cross.Mag()) > 0.)) // dir1 || dir2
    {
        return calculateMinimumDistanceStraightToPoint(base1, dir1, base2);
    }
    return fabs(ab.Dot(cross) / cross.Mag());
}

TVector3 HParticleTool::calculateCrossPoint(TVector3 &base1, TVector3 &dir1, TVector3 &base2,
                                            TVector3 &dir2)
{
    Double_t d1d1 = dir1(0) * dir1(0) + dir1(1) * dir1(1) + dir1(2) * dir1(2);
    Double_t d2d2 = dir2(0) * dir2(0) + dir2(1) * dir2(1) + dir2(2) * dir2(2);
    Double_t d1d2 = dir1(0) * dir2(0) + dir1(1) * dir2(1) + dir1(2) * dir2(2);

    Double_t D = d1d1 * d2d2 - (d1d2 * d1d2);

    if (!(fabs(D) > 0.))
    {
        ::Warning("calculateCrossPoint",
                  "Error while calculating cross point ... eqns are lin. "
                  "dependent:returning default Vertex (-20000,-20000,-20000)");

        return TVector3(-20000., -20000., -20000.);
    }

    Double_t d1diff = dir1(0) * (base2(0) - base1(0)) +
                      dir1(1) * (base2(1) - base1(1)) +
                      dir1(2) * (base2(2) - base1(2));
    Double_t d2diff = dir2(0) * (base2(0) - base1(0)) +
                      dir2(1) * (base2(1) - base1(1)) +
                      dir2(2) * (base2(2) - base1(2));

    Double_t Dlambda = d1diff * d2d2 - d1d2 * d2diff;

    Double_t lambda = Dlambda / D;

    TVector3 vertex;
    vertex += dir1;
    vertex *= lambda;
    vertex += base1;

    return TVector3(vertex);
}
Double_t HParticleTool::calcDeterminant(TVector3 &v1, TVector3 &v2, TVector3 &v3)
{
    return (v1(0) * v2(1) * v3(2) + v2(0) * v3(1) * v1(2) +
            v3(0) * v1(1) * v2(2) - v3(0) * v2(1) * v1(2) -
            v1(0) * v3(1) * v2(2) - v2(0) * v1(1) * v3(2));
}
TVector3 HParticleTool::calculatePointOfClosestApproach(TVector3 &base1, TVector3 &dir1,
                                                        TVector3 &base2, TVector3 &dir2)
{

    TVector3 cross = dir1.Cross(dir2); // cross product: dir1 x dir2

    // straight lines are either skew or have a cross point

    TVector3 diff = base1;
    diff -= base2; // Difference of two base vectors base1 - base2

    Double_t D;
    D = calcDeterminant(dir2, dir1, cross);

    if (!(fabs(D) > 0.))
    {
        ::Warning(":calculatePointOfClosestApproach",
                  "Dirs and cross-product are lin. dependent: returning default "
                  "Vertex (-20000,-20000,-20000)");

        return TVector3(-20000., -20000., -20000.);
    }

    Double_t Dm = calcDeterminant(diff, dir1, cross);
    Double_t Dl = -calcDeterminant(diff, dir2, cross);

    TVector3 vertex;
    TVector3 dm;
    TVector3 dl;

    dm = dir2;
    dm *= Dm;

    dl = dir1;
    dl *= Dl;

    vertex = dm - dl;

    vertex *= ((1.) / D);

    vertex += base1;
    vertex += base2;
    vertex *= 0.5;

    return TVector3(vertex);
}
TVector3 HParticleTool::calcVertexAnalytical(TVector3 &base1, TVector3 &dir1, TVector3 &base2,
                                             TVector3 &dir2)
{
    // 1. exists a unique solution ?

    if ((dir1.Cross(dir2)).Mag() > 0.) // dir1 and dir2 linear independent
    {
        // straight lines are either skew or have a cross point

        TVector3 diff = base1;
        diff -= base2; // Difference of two base vectors base1 - base2

        // 2. skew or intersecting ?

        if (fabs(calcDeterminant(dir2, dir1, diff)) > 0.)
        {
            // 3. (b) skew
            return TVector3(
                calculatePointOfClosestApproach(base1, dir1, base2, dir2));
        }
        else
        {
            // 3. (a) intersection
            return TVector3(calculateCrossPoint(base1, dir1, base2, dir2));
        }
    }
    else
    {
        // dir1 and dir2 linear dependent -> g1 and g2 identical or parallel
        return TVector3(-10000000., -10000000., -10000000.);
    }
    return TVector3(-10000000., -10000000., -10000000.);
}
