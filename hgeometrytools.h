/**
 * HGeometryTools.h
 *
 *
 */

#ifndef HGEOMETRYTOOLS_H
#define HGEOMETRYTOOLS_H

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// framework includes
#include "hrefitcand.h"
#include "hgeomvector.h"
#include "hparticletool.h"
#include "hgeomvertexfit.h"
#include "hgeommatrix.h"

using std::cout;
using std::endl;

class HGeometryTools : public TObject

{
private:

    HGeomMatrix fM; // Temporal matrix for calculations
protected:
    HGeomMatrix fSys; // LSM system inverse matrix
    HGeomVector fB;   // LSM independent term
public:
    HGeometryTools();
    ~HGeometryTools(){};
    
    void getVertex(HGeomVector &out);
    void addLine(const HGeomVector &r, const HGeomVector &alpha, const Double_t w);
    void getVertex(HGeomVector &out);
    void reset(void);

   ClassDef(HGeometryTools,0)
};

#endif /* HGEOMETRYTOOLS_H */
