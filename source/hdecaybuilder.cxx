#include "hdecaybuilder.h"

HDecayBuilder::HDecayBuilder(std::vector<std::vector<HRefitCand>> &cands, TString &task, std::vector<Int_t> &pids, TLorentzVector lv, HRefitCand mother, Double_t mass) : fCands(cands),
                                                                                                                                                                          fTask(task),
                                                                                                                                                                          fPids(pids),
                                                                                                                                                                          fCombiCounter(0),
                                                                                                                                                                          fProb(0.01),
                                                                                                                                                                          fVerbose(0)

{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder() -----------------" << std::endl;
    }

    setIniSys(lv);
    setMother(mother);
    setMass(mass);

    // Determine the number of combinations of the input particles
    fTotalCombos = 1;
    if(particleCounter.size()>0) particleCounter.clear();
    for (size_t i = 0; i < fPids.size(); i++)
    {
        fTotalCombos *= fCands[i].size();
        particleCounter.push_back(0);
    }
}

void HDecayBuilder::buildDecay()
{

    cout << "Combi counter: " << fCombiCounter << " total Combos: " << fTotalCombos << endl;
    cout << "Task is " << fTask << endl;

    // For selection of best event according to probability
    fBestProb = 0;
    fBestChi2 = 1e6;

    while (fCombiCounter < (fTotalCombos ))
    { // double check this
        cout << "fill fit cands" << endl;
        // Take one combination of particles
        fillFitCands();
        // Check if correct number of particles
        if (fFitCands.size() != fPids.size())
        {
            cout << "wrong number of particles" << endl;
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
        cout << "Combi counter: " << fCombiCounter << " total Combos: " << fTotalCombos << endl;
        cout << "fit finished" <<endl;
    }
}

/*
void HDecayBuilder::createNeutralCandidate()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::createNeutralCandidate() -----------------" << std::endl;
    }

    std::vector<HRefitCand> candsPrim;
    std::vector<HRefitCand> candsDecay;
    Int_t tempPid;

    for (Int_t i_pids = 0; i_pids = (Int_t)fPids.size(); i_pids++)
    {

        tempPid = fPids[i_pids];

        for (Int_t i_cand = 0; i_cand < (Int_t)fCands[i_pids].size(); i_cand++)
        {

            std::vector<Int_t>::iterator it_prim = std::find(fPidsPrim.begin(), fPidsPrim.end(), tempPid);

            if (it_prim != fPidsPrim.end())
            {
                candsPrim.push_back(fCands[i_pids][i_cand]);
            }

            std::vector<Int_t>::iterator it_decay = std::find(fPidsDecay.begin(), fPidsDecay.end(), tempPid);

            if (it_decay != fPidsDecay.end())
            {
                candsDecay.push_back(fCands[i_pids][i_cand]);
            }
        }
    }

    TVector3 primVertex;
    TVector3 decayVertex;

    // Find first vertex from all particle combinations if user set PIDs

    HVertexFinder *vtxFinderPrim = new HVertexFinder();
    primVertex = vtxFinderPrim->findVertex(candsPrim);

    // Find second vertex from all particle combinations of user set PIDs

    HVertexFinder *vtxFinderSec = new HVertexFinder();
    decayVertex = vtxFinderSec->findVertex(candsDecay);

    HNeutralCandFinder candFinder(candsDecay);
    candFinder.setNeutralMotherCand(primVertex, decayVertex);

    fFitCands = candsDecay;
    fMother = candFinder.getNeutralMotherCandidate();

    do3cFit();
}
*/

Bool_t HDecayBuilder::doFit()
{
    HKinFitter Fitter(fFitCands);
    if (fTask == "4C") Fitter.add4Constraint(fIniSys);
    else if (fTask == "Mass"){
        Fitter.addMassConstraint(fMass);
        cout << "constraint added, Mass = " << fMass << endl;
        cout << "proton momentum = " << fFitCands[0].getMomentum() << endl;
        cout << "pion momentum = " << fFitCands[1].getMomentum() << endl;
        TMatrixD cov_p = fFitCands[0].getCovariance();
        TMatrixD cov_pi = fFitCands[1].getCovariance();
        cout << "proton momentum error = " << cov_p(0,0) << endl;
        cout << "pion momentum error = " << cov_pi(0,0) << endl;

    } 
    else cout << "Task not available" << endl;

    Fitter.setNumberOfIterations(10);
    Fitter.setConvergenceCriterion(0.01);
    //Fitter.setVerbosity(2);

    if (Fitter.fit() && Fitter.getProb() > fProb)
    {
        cout << "fit successful" << endl;
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
        cout << "fit not successful" << endl;
        return kFALSE;
    }
}
/*
bool HDecayBuilder::do3cFit()
{
    createNeutralCandidate();
    HKinFitter Fitter(fFitCands, fMother);
    Fitter.add3Constraint();
    Fitter.fit();
}

bool HDecayBuilder::doMissMomFit()
{
    HKinFitter Fitter(fFitCands, fIniSys, fMass);
    Fitter.addMomConstraint();
    Fitter.fit();
}
*/
void HDecayBuilder::fillFitCands()
{

    if(fFitCands.size()>0) fFitCands.clear();
    doubleParticle = false;
    for (size_t i = 0; i < fPids.size(); i++)
    {
        checkDoubleParticle(i);
        cout << "particle counter " << i << " is " << particleCounter[i] << " now" << endl;
        fFitCands.push_back(fCands[i][particleCounter[i]]);
        cout << "fill in pid " << fPids[i] << endl;
    }
    Int_t a = fPids.size() - 1;
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

void HDecayBuilder::checkDoubleParticle(size_t i)
{
    for (size_t j = 0; j < i; j++)
    {
        if ((fPids[j] == fPids[i]) && (particleCounter[j] == particleCounter[i]))
        {
            cout << "double particle" << endl;
            doubleParticle = true;
        }
    }
    cout << "no double particle" << endl;
}

/*
void HDecayBuilder::createOutputParticle(HRefitCand refitCand)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::createOutputParticle() -----------------" << std::endl;
    }
    HRefitCand newParticle;

    newParticle.setMomentum(refitCand.P());
    newParticle.setTheta(refitCand.Theta());
    newParticle.setPhi(refitCand.Phi());
    newParticle.setR(refitCand.getR());
    newParticle.setZ(refitCand.getZ());
    newParticle.setCovariance(refitCand.getCovariance());
    newParticle.setPid(refitCand.getPid());

    fOutputCands.push_back(newParticle);
}
*/