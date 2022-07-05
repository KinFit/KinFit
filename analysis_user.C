/*
#include "hades.h"
#include "hloop.h"
#include "htool.h"
#include "hcategorymanager.h"
#include "hparticleanglecor.h"
#include "hparticlepairmaker.h"
#include "hparticletool.h"
#include "hphysicsconstants.h"
#include "hhistmap.h"
#include "hparticletracksorter.h"
#include "henergylosscorrpar.h"


#include "hcategory.h"
#include "hlinearcategory.h"
#include "hrichhit.h"
#include "hrichhitsim.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticlepair.h"
#include "hparticlegeantpair.h"

#include "hgeantkine.h"
#include "hparticledef.h"
#include "hstartdef.h"
#include "richdef.h"

#include "hparticlegeant.h"
#include "hparticlegeantdecay.h"
#include "hparticlegeantevent.h"
#include "hparticlecutrange.h"

#include "TTree.h"
*/

//#include "TLorentzVector.h"

#include <vector>
#include <algorithm>
#include <map>
#include <iostream>
#include <iomanip>
#include <math.h>

/*
#include "hkinfitter.h"
#include "hvertexfinder.h"
#include "hneutralcandfinder.h"
*/
#include "hdstfitter.h"

using namespace std;
using namespace Particle;

Int_t analysis_user(TString infileList="/lustre/hades/user/jrieger/pp_pKLambda/sim/pp_pKlambda_100000evts1_dst_apr12.root", Int_t nEvents=1000){

    HDSTFitter DSTFitter(infileList, false, false, nEvents);
    //HDSTFitter DSTFitter(infileList);
    //std::vector<Int_t> pids{ 14, 11, 14, 9 };
    std::vector<int> pids;
    std::vector<int> pidsPrimVertex; // Jenny: User must set this for 3C fit
    stc::vector<int> pidsDecayVertex; // Jenny: User must set this for 3C fit
    pids.push_back(14); pids.push_back(11); pids.push_back(14); pids.push_back(9);
    TLorentzVector ppSystem(0,0,4337.96,2*938.272+3500);

    DSTFitter.addFitterTask("4c", pids, ppSystem);
    
    return 0;
}
