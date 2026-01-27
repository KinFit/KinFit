// A million includes for looping DSTs that should get out of the way
#ifndef __INCLUDES_H__
#define __INCLUDES_H__

#include "TObject.h"
#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "TList.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TCutG.h"
#include "TF1.h"


#include "hloop.h"
#include "hcategorymanager.h"
#include "hparticletool.h" // Jenny added
#include "hgeomvector.h" // Jenny added
#include "hgeomvertexfit.h" // Jenny added
#include "htool.h"
#include "hphysicsconstants.h"
#include "hparticletracksorter.h"
#include "hfrpchit.h"
#include "hforwardtools.h"

#include "heventheader.h"
#include "hparticleevtinfo.h"
#include "hgeantheader.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hgeantkine.h"
#include "hforwardcand.h"
#include "forwarddef.h" // Jenny added

#include "/lustre/hades/user/jregina/ExternalLibsHydra/lib/KinFitter.h"
//#include "/lustre/hades/user/jregina/ExternalLibsHydra/lib/KFitVertexFinder.h"
//#include "/lustre/hades/user/jregina/ExternalLibsHydra/lib/KFitDecayCandFinder.h"
//#include "KinFitter.h"
//#include "KFitVertexFinder.h"
//#include "KFitDecayCandFinder.h"

#include "hparticlebeamtilt.h"
#include "henergylosscorrpar.h"

#include "hmdcsizescells.h"
// For Start detector information: Added from gen02 onwards
#include "hstart2cal.h"
#include "hstart2hit.h"
#include "hstartdef.h"

#include "hitofdef.h"
#include "hitofcal.h"
#include "hitofdetector.h"
// #include "hstart2clusterfinder.h"

#include "hstsgeompar.h"
#include "hgeomvolume.h"
#include "hgeomcompositevolume.h"
#include "hspectrometer.h"
#include "hfrpcdetector.h"
#include "hstsdetector.h"
#include "forwarddef.h"
#include "stsdef.h"
#include "frpcdef.h"

#include "hdetector.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#endif // __INCLUDES_H__
