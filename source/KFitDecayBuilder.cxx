#include "KFitDecayBuilder.h"

KFitDecayBuilder::KFitDecayBuilder(TString task, std::vector<int> pids, TLorentzVector lv, KFitParticle mother, double mass) : fTask(task),
                                                                                                                            fPids(pids),
                                                                                                                            fVerbose(0)

{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder() -----------------" << std::endl;
    }

    setIniSys(lv);
    setMother(mother);
    setMass(mass);

}

void KFitDecayBuilder::countCombis()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::countCombis() -----------------" << std::endl;
    }
    fCombiCounter = 0;
    // Determine the number of combinations of the input particles
    fTotalCombos = 1;
    if(particleCounter.size()>0) particleCounter.clear();
    for (size_t i = 0; i < fPids.size(); i++)
    {
        fTotalCombos *= fCands[i].size();
        particleCounter.push_back(0);
    }

}

void KFitDecayBuilder::buildDecay()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::buildDecay() -----------------" << std::endl;
    }
    if (fVerbose > 1)
    {
        std::cout << "Combi counter: " << fCombiCounter << " total Combos: " << fTotalCombos << std::endl;
        std::cout<< "Task is " << fTask <<std::endl;
    }

    // For selection of best event according to probability
    fBestProb = 0;
    fBestChi2 = 1e6;

    while (fCombiCounter < (fTotalCombos ))
    {
        cout << "fill fit cands" << endl;
        // Take one combination of particles
        fillFitCands();
        
        // Check if correct number of particles
        if (fFitCands.size() != fPids.size())
        {
            cout << "wrong number of particles" << endl;
            fCombiCounter++;
            continue;
        }
        // Check for double usage of same particle
        if (doubleParticle)
        {
            cout << "particle used twice" << endl;
            continue;
        }
        // Do task that was chosen by user
        doFit();

        if (fVerbose > 1)
        {
            std::cout << "Combi counter: " << fCombiCounter << " total Combos: " << fTotalCombos << std::endl;
            std::cout << "fit finished" <<std::endl;
        }
        
        
        
    }
    
}

/*
void KFitDecayBuilder::createNeutralCandidate()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::createNeutralCandidate() -----------------" << std::endl;
    }

    std::vector<KFitParticle> candsPrim;
    std::vector<KFitParticle> candsDecay;
    int tempPid;

    for (int i_pids = 0; i_pids = (int)fPids.size(); i_pids++)
    {

        tempPid = fPids[i_pids];

        for (int i_cand = 0; i_cand < (int)fCands[i_pids].size(); i_cand++)
        {

            std::vector<int>::iterator it_prim = std::find(fPidsPrim.begin(), fPidsPrim.end(), tempPid);

            if (it_prim != fPidsPrim.end())
            {
                candsPrim.push_back(fCands[i_pids][i_cand]);
            }

            std::vector<int>::iterator it_decay = std::find(fPidsDecay.begin(), fPidsDecay.end(), tempPid);

            if (it_decay != fPidsDecay.end())
            {
                candsDecay.push_back(fCands[i_pids][i_cand]);
            }
        }
    }

    TVector3 primVertex;
    TVector3 decayVertex;

    // Find first vertex from all particle combinations if user set PIDs

    KFitVertexFinder *vtxFinderPrim = new KFitVertexFinder();
    primVertex = vtxFinderPrim->findVertex(candsPrim);

    // Find second vertex from all particle combinations of user set PIDs

    KFitVertexFinder *vtxFinderSec = new KFitVertexFinder();
    decayVertex = vtxFinderSec->findVertex(candsDecay);

    KFitNeutralCandFinder candFinder(candsDecay);
    candFinder.setNeutralMotherCand(primVertex, decayVertex);

    fFitCands = candsDecay;
    fMother = candFinder.getNeutralMotherCandidate();

    do3cFit();
}
*/

Bool_t KFitDecayBuilder::doFit()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::doFit() -----------------" << std::endl;
    }
    KinFitter Fitter(fFitCands);
    if (fTask == "4C") Fitter.add4Constraint(fIniSys);
    else if (fTask == "Vertex") Fitter.addVertexConstraint();
    else if (fTask == "MassVtx") Fitter.addMassVtxConstraint(fMass);
    else if (fTask == "MM" || fTask == "MissingMass" ) Fitter.addMissingMassConstraint(fIniSys, fMass);
    else if (fTask == "Mom" || fTask == "MissingParticle") Fitter.addMissingParticleConstraint(fIniSys, fMass);
    else if (fTask == "Mass"){
        Fitter.addMassConstraint(fMass);
        /*
        cout << "constraint added, Mass = " << fMass << endl;
        cout << "proton momentum = " << fFitCands[0].getMomentum() << endl;
        cout << "pion momentum = " << fFitCands[1].getMomentum() << endl;
        TMatrixD cov_p = fFitCands[0].getCovariance();
        TMatrixD cov_pi = fFitCands[1].getCovariance();
        /*
        cout << "proton momentum error = " << cov_p(0,0) << endl;
        cout << "pion momentum error = " << cov_pi(0,0) << endl;
*/
    } 
    else cout << "Task not available" << endl;

    if (Fitter.fit())
    {
        if (fVerbose > 1)
        {
            std::cout << "fit successful" << std::endl;
        }
        if (Fitter.getProb() > fBestProb)
        {
            fBestProb = Fitter.getProb();
            fBestChi2 = Fitter.getChi2();
            if(fOutputCands.size() != 0) fOutputCands.clear();
            Fitter.getDaughters(fOutputCands);
        }
        return kTRUE;
    }
    else
    {
        if (fVerbose > 1)
        {
            std::cout << "fit not successful" << std::endl;
        }
        return kFALSE;
    }
}

void KFitDecayBuilder::fillFitCands()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::fillFitCands() -----------------" << std::endl;
    }

    if(fFitCands.size()>0) fFitCands.clear();
    doubleParticle = false;
    
    for (size_t i = 0; i < fPids.size(); i++)
    {
        checkDoubleParticle(i);
        cout << "particle counter " << i << " is " << particleCounter[i] << " now" << endl;
        fFitCands.push_back(fCands[i][particleCounter[i]]);
        cout << "fill in pid " << fPids[i] << endl;
    }
    int a = fPids.size() - 1;
    cout << "particle counter " << a << " is " << particleCounter[a] << " now" << endl;
    while (particleCounter[a] == (fCands[a].size() - 1))
    {
        particleCounter[a] = 0;
        cout << "particle counter " << a << " is " << particleCounter[a] << " now" << endl;
        a--;
        if (a < 0)
        {
            if (!(fCombiCounter == fTotalCombos-1))
                cout << "counted wrong: " << fCombiCounter << " != " << fTotalCombos << endl;
            break;

            a = 1e6;
        }
    }
    particleCounter[a]++;
    cout << "particle counter " << a << " is " << particleCounter[a] << " now" << endl;
    fCombiCounter++;
    

    if (doubleParticle && (fCombiCounter < fTotalCombos))
        fillFitCands(); // If some particle has been filled more than once into fFitCands, repeat the procedure with the next combination
}

void KFitDecayBuilder::checkDoubleParticle(size_t i)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::checkDoubleParticle() -----------------" << std::endl;
    }
    for (size_t j = 0; j < i; j++)
    {
        if ((fPids[j] == fPids[i]) && (particleCounter[j] == particleCounter[i]))
        {
            cout << "double particle" << endl;
            doubleParticle = true;
        }
    }
    if (fVerbose > 1)
    {
        std::cout << "no double particle" << std::endl;
    }
}