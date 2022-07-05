#include "hdecaybuilder.h"

HDecayBuilder::HDecayBuilder(std::vector< std::vector<HRefitCand> > &cands, TString &task, std::vector<Int_t> &pids, TLorentzVector lv, HRefitCand mother, Double_t mass) : fCands(cands), 
                                                                                                                                                            fTask(task),
                                                                                                                                                            fPids(pids),
                                                                                                                                                            fCombiCounter(0),
                                                                                                                                                            fProb(0),
                                                                                                                                                            fVerbose(0)

{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder() -----------------" << std::endl;
    }
    
    setIniSys(lv);
    setMother(mother);
    setMass(mass);
    
    fTotalCombos = 1;
    particleCounter.clear();
    for (size_t i=0; i<fPids.size(); i++){
		fTotalCombos *= fCands[i].size();
		particleCounter.push_back(0);
	}
}

void HDecayBuilder::buildDecay(){
    
    cout<<"Combi counter: " << fCombiCounter << " total Combos: " << fTotalCombos <<endl;
    cout<<"Task is "<< fTask << endl;
	for (size_t i=0; i<particleCounter.size(); i++){
		cout<< "particle counter "<<i<<" = "<<particleCounter[i]<<endl;
	}
    while(fCombiCounter<(fTotalCombos-1)){	//double check this
        cout<<"fill fit cands"<<endl;
        fillFitCands();
        if (fFitCands.size()!=fPids.size()){
			cout<<"wrong number of particles"<<endl;
			continue;
		}
		if (doubleParticle){
			cout<<"particle used twice"<<endl;
			continue;
		}
        if(fTask == "createNeutral"){
            createNeutralCandidate();
        }else if(fTask == "3C"){
            if(do3cFit()) cout<<"3C task successful"<<endl;
        }else if(fTask = "4C"){
			cout<<"4C task received"<<endl;
            if(do4cFit()) cout<<"4C task successful"<<endl;
        }else if(fTask == "missMom"){
            if(doMissMomFit()) cout<<"miss mom task successful"<<endl;
        }else{
            cout<<"Task not available"<<endl;
        }
        cout<<"Combi counter: " << fCombiCounter << " total Combos: " << fTotalCombos <<endl;
    }
    

}

void HDecayBuilder::createNeutralCandidate()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::createNeutralCandidate() -----------------" << std::endl;
    }

    std::vector<HRefitCand> candsPrim;
    std::vector<HRefitCand> candsDecay;
    int tempPid;

    for (int i_pids = 0; i_pids = (int)fPids.size(); i_pids++)
    {

        tempPid=fPids[i_pids];

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

    HVertexFinder *vtxFinderPrim = new HVertexFinder();
    primVertex = vtxFinderPrim->findVertex(candsPrim);

    // Find second vertex from all particle combinations of user set PIDs

    HVertexFinder *vtxFinderSec = new HVertexFinder();
    decayVertex = vtxFinderSec->findVertex(candsDecay);

    HNeutralCandFinder candFinder(candsDecay);
    candFinder.setNeutralMotherCand(primVertex, decayVertex);

    fFitCands=candsDecay;
    fMother=candFinder.getNeutralMotherCandidate();

    do3cFit();

}

/*  
    


    bool bestPrimVertexFound=false;
    bool bestDecayVertexFound=false;

    TVector3 decayVertex;
    TVector3 primVertex;
    
    int indexPrimaryProton=-1, indexSecondBestPrimaryProton=-1;
    int indexDecayProton=-1, indexSecondBestDecayProton=-1;

    std::vector<HRefitCand> cands3c;
    cands3c.clear();
    fOutputCands.clear();

    double probPrim = -99999, probSec = -99999;
    double probSecondBestPrim = -99999, probSecondBestDecay = -99999;
    double probDecayVertex_Temp = -1;
    double probPrimVertex_Temp = -1;
    
    // Perform the analysis
    for (size_t n = 0; n < fProtons.size(); n++)
    {
        HRefitCand* cand1Ptr = fProtons[n];
        HRefitCand cand1= *cand1Ptr;
        
        for (size_t m = 0; m < fPions.size(); m++)
        {
            HRefitCand* cand2Ptr = fPions[m]; 
            HRefitCand cand2= *cand2Ptr;               
            
            std::vector<HRefitCand> candsSec;
            candsSec.clear();
            candsSec.push_back(cand1);
            candsSec.push_back(cand2);
            HVertexFinder *vtxFinderSec = new HVertexFinder();
            decayVertex = vtxFinderSec->findVertex(candsSec);
            
            // Kinematic fitting with vertex constraint
            HKinFitter vtxFitterSecCands(candsSec); 

            vtxFitterSecCands.addVertexConstraint();
            vtxFitterSecCands.fit();
            probSecondBestDecay=probSec;
            probSec = vtxFitterSecCands.getProb();

            if (probSec > probDecayVertex_Temp)
            {
                bestDecayVertexFound = true;
                probDecayVertex_Temp = probSec;
                indexSecondBestDecayProton = indexDecayProton;
                indexDecayProton = n;

                cands3c.clear();
                // Get the proton daughter and pass it to the 3C fit later
                cands3c.push_back(vtxFitterSecCands.getDaughter(0));
                // Get the pion daughter and pass it to the 3C fit later
                cands3c.push_back(vtxFitterSecCands.getDaughter(1));
                //std::cout << "vertex position: " << decayVertex.X() << " " << decayVertex.Y() << " " << decayVertex.Z() << std::endl;
            }
            
        }

        for (size_t p = 0; p < fKaons.size(); p++)
        {
            HRefitCand* cand3Ptr = fKaons[p];
            HRefitCand cand3= *cand3Ptr;

            std::vector<HRefitCand> candsPrim;
            candsPrim.clear();
            candsPrim.push_back(cand1);
            candsPrim.push_back(cand3);
            HVertexFinder *vtxFinderPrim = new HVertexFinder();
            primVertex = vtxFinderPrim->findVertex(candsPrim);            
            
            HKinFitter vtxFitterPrimCands(candsPrim);
            vtxFitterPrimCands.addVertexConstraint();
            vtxFitterPrimCands.fit();
            probSecondBestPrim=probPrim;
            probPrim = vtxFitterPrimCands.getProb();

            if (probPrim > probPrimVertex_Temp)
            {
                bestPrimVertexFound = true;
                probPrimVertex_Temp = probPrim;
                indexSecondBestPrimaryProton = indexPrimaryProton;
                indexPrimaryProton = n;                    
                
                fCandsMissPos.push_back(cand1);
                fCandsMissPos.push_back(cand3);
            }
        }
    }

    if (bestPrimVertexFound == true && bestDecayVertexFound == true)
    {

        if (indexDecayProton != indexPrimaryProton && probPrimVertex_Temp>fPrimVertexProbabilityCut && probDecayVertex_Temp>fDecayVertexProbabilityCut)
        {
            HNeutralCandFinder lambdaCandFinder(cands3c);

            lambdaCandFinder.setUsePrimaryVertexInNeutralMotherCalculation(true);

            lambdaCandFinder.setNeutralMotherCandFromPrimaryVtxInfo(primVertex, decayVertex);

            HVirtualCand lambdaCand = lambdaCandFinder.getNeutralMotherCandidate();

            HRefitCand lambdaCandRefit(&lambdaCand);
            lambdaCandRefit.SetXYZM(lambdaCand.getMomentum() * std::sin(lambdaCand.getTheta() * deg2rad) *
                                        std::cos(lambdaCand.getPhi() * deg2rad),
                                    lambdaCand.getMomentum() * std::sin(lambdaCand.getTheta() * deg2rad) *
                                        std::sin(lambdaCand.getPhi() * deg2rad),
                                    lambdaCand.getMomentum() * std::cos(lambdaCand.getTheta() * deg2rad),
                                    1115.683);

            TMatrixD lambdaCov(5, 5);
            lambdaCov = lambdaCandFinder.getCovarianceMatrixNeutralMother();
            lambdaCandRefit.setCovariance(lambdaCov);
            lambdaCandRefit.setR(lambdaCand.getR());
            lambdaCandRefit.setZ(lambdaCand.getZ());
            lambdaCandRefit.SetTheta(lambdaCand.getTheta() * deg2rad);
            lambdaCandRefit.SetPhi(lambdaCand.getPhi() * deg2rad);

            HKinFitter Fitter3c(cands3c, lambdaCandRefit);
            Fitter3c.add3Constraint();
            Fitter3c.setNumberOfIterations(20);

            Fitter3c.fit();

            HRefitCand cand13C = Fitter3c.getDaughter(0); // proton
            createOutputParticle(cand13C);
            HRefitCand cand23C = Fitter3c.getDaughter(1); // pion
            createOutputParticle(cand23C);
            
            HRefitCand lambdaCand3C = Fitter3c.getMother();
            fCandsMissPos.push_back(lambdaCand3C);

        }
    }
}*/

bool HDecayBuilder::do4cFit()
{
    HKinFitter Fitter(fFitCands, fIniSys);
    Fitter.add4Constraint();
    cout<<"constraint added"<<endl;
    if(Fitter.fit() && Fitter.getProb()>fProb){
		cout<<"fit successful"<<endl;
        fOutputCands.clear();
        Fitter.getDaughters(fOutputCands);
        return true;
        /*
        for(iterator it = fPids.begin(); it != fPids.end(); ++it){
            fOutputCands.push_back(Fitter.getDaughter(it));
        }*/
        //fillHistograms();
    } else {
		return false;
	}
}

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

void HDecayBuilder::fillFitCands()
{
    fFitCands.clear();
    doubleParticle = false;
    for (size_t i=0; i<fPids.size(); i++){
		checkDoubleParticle(i);
		cout<<"particle counter "<<i<< " is "<< particleCounter[i] << " now"<<endl;
        fFitCands.push_back(fCands[i][particleCounter[i]]);
        cout<<"fill in pid "<<fPids[i]<<endl;
	}
	Int_t a = fPids.size()-1;
	cout<<"particle counter "<<a<< " is "<< particleCounter[a] << " now"<<endl;
	while( particleCounter[a]==(fCands[a].size()-1) ){ 
		particleCounter[a] = 0;
		cout<<"particle counter "<<a<< " is "<< particleCounter[a] << " now"<<endl;
		a--;
		if(a<0){
			if(!(fCombiCounter==fTotalCombos)) cout<<"counted wrong: "<<fCombiCounter<<" != "<<fTotalCombos<<endl;
			break;
		}
	}
	particleCounter[a]++;
	cout<<"particle counter "<<a<< " is "<< particleCounter[a] << " now"<<endl;
	fCombiCounter++;
	
    if(doubleParticle && (fCombiCounter<fTotalCombos)) fillFitCands();	//If some particle has been filled more than once into fFitCands, repeat the procedure with the next combination
}

void HDecayBuilder::checkDoubleParticle(size_t i)
{
	for (size_t j=0; j<i; j++){
		if((fPids[j]==fPids[i]) && (particleCounter[j]==particleCounter[i])) 
		{
			cout<< "double particle" << endl;
			doubleParticle = true;
		}
	}
	cout<< "no double particle" << endl;
}

void HDecayBuilder::createOutputParticle(HRefitCand refitCand)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::createOutputParticle() -----------------" << std::endl;
    }
    HRefitCand newParticle;

    newParticle.setPhi(refitCand.Theta());//!!!
    newParticle.setR(refitCand.getR());
    newParticle.setZ(refitCand.getZ());
    newParticle.setMomentum(refitCand.P());

    fOutputCands.push_back(newParticle);
}

void HDecayBuilder::createOutputCategory()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::createOutputCategory() -----------------" << std::endl;
    }

    HCategoryManager catManager;
    HCategory *cat = catManager.addCategory(1, "HParticleCandSimAfterFit");
}
/*
void HDecayBuilder::fillHistograms()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HDecayBuilder::writeHistograms() -----------------" << std::endl;
    }

    hmLam_prefit->Fill((fFitCands[2]+fFitCands[3]).M());
    hmLam_post4C->Fill((fOutputCands[2]+fOutputCands[3]).M());
}*/
