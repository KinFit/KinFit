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
	bool fMomDepErrors = false;
	Int_t fVerbose = 0;

  static   HCovarianceKinFit* gCovariance;
  ClassDef(HCovarianceKinFit,1) // Parameter container for energy loss correction. Wozu?
  
public:
  /*HCovarianceKinFit(const Char_t* name    = "CovarianceKinFit",
                     const Char_t* title   = "Covariance matrix estimation for kin fit",
                     const Char_t* context = "CovarianceMatrixEstimate");
  */
  HCovarianceKinFit();
  ~HCovarianceKinFit(void) {}
  static HCovarianceKinFit* getObject(void) {return gCovariance;}
  
  void estimateCov(Int_t pid, Double_t mom, double (&covariance)[5]);
  void setSetup(TString run){ fSetup=run; }
  void setMomDepErrors(bool val){ fMomDepErrors=val; }
};

#endif  /*HCOVARIANCEKINFIT_H */
