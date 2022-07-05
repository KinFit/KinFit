#ifndef HCOVARIANCEKINFIT_H
#define HCOVARIANCEKINFIT_H

#include "TMath.h"
#include "TString.h"
#include "hphysicsconstants.h"  //knows conversion from id to mass
//#include "hparamlist.h"
#include "TF1.h"
#include "TFile.h"

#include <stdlib.h>
#include <cmath>

using namespace TMath;

//#include "hparcond.h"    is this useful here?

//class HCovarianceKinFit : public HParCond {
class HCovarianceKinFit {
private:
	TString fSetup;
	bool fMomDepErrors;
	Int_t fVerbose;

  //static   HCovarianceKinFit* gCovariance;
  
public:
  /*HCovarianceKinFit(const Char_t* name    = "CovarianceKinFit",
                     const Char_t* title   = "Covariance matrix estimation for kin fit",
                     const Char_t* context = "CovarianceMatrixEstimate");
  */
  HCovarianceKinFit();
  ~HCovarianceKinFit(void) {}
  //static HCovarianceKinFit* getObject(void) {return gCovariance;}
  
  void estimateCov(Int_t pid, Double_t mom, double (&covariance)[5]);
  void setSetup(TString run){ fSetup=run; }
  void setMomDepErrors(bool val){ fMomDepErrors=val; }
/*
private:
  ClassDef(HCovarianceKinFit,1) // Parameter container for energy loss correction. Wozu?
*/
};

#endif  /*HCOVARIANCEKINFIT_H */
