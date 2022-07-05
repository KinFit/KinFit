#include "hdstfitter.h"
#include "hphysicsconstants.h"  //knows conversion from id to mass

HDSTFitter::HDSTFitter(TString infileList, bool includeFw, bool momDepErrors, Int_t nEvents) : fInfileList(infileList),
																								fIncludeFw(includeFw),
																								fMomDepErrors(momDepErrors),
																								fEvents(nEvents),
																								fVerbose(0)
{
}

void HDSTFitter::selectCandidates()
{
    // for each event there are a number of tracks
    Int_t ntracks = fcatParticle->getEntries();
    if(fIncludeFw) Int_t nFwTracks = fcatFwParticle->getEntries();

    //Object for covariance matrix estimation
    HCovarianceKinFit cov;
    cov.setSetup("pp35");

    fCandsFit.clear();

    for (size_t it = 0; it < fPids.size(); it++)
    {
        std::vector<HRefitCand> tempVec;

        for (Int_t j = 0; j < ntracks; j++)
        {
            //HParticleCand* cand = HCategoryManager::getObject(cand, catParticle, j);
            HParticleCandSim *cand = HCategoryManager::getObject(cand, fcatParticle, j);
            // skip ghost tracks (only avalible for MC events)
            if (cand->isGhostTrack())
                continue;
            // select "good" tracks
            if (!cand->isFlagBit(Particle::kIsUsed))
                continue;

            HRefitCand candidate(cand);

            //for(auto it = std::begin(fPids); it != std::end(fPids); ++it) {
            if (cand->getGeantPID() == fPids[it])
            {
                Double_t mom = cand->P();
                //vector<double> errors;
                //getErrors(fPids[it], mom, errors);
                Double_t errors[5];
                cov.estimateCov(fPids[it], mom, errors);
                Double_t mCand = HPhysicsConstants::mass(fPids[it]);
                //std::cout <<  << std::endl;
                FillData(cand, candidate, errors, mCand);
                tempVec.push_back(candidate);
            }
        }
        //cout << "Iter: " << it << endl;
        //cout << "tempvec size " << tempVec.size() << endl; 
        fCandsFit.push_back(tempVec);
        //cout << "fCandsFit size " << fCandsFit[it].size() << endl; 

    } // end of HADES track loop
    /*
    if(fIncludeFw){ //find correct index of proton HRefitCand vector
        for(Int_t j=0; j<nFwTracks; j++){
            HFwDetCandSim *cand = HCategoryManager::getObject(cand,catFwParticle,j);
            cand->calc4vectorProperties(938.272);
                
            HRefitCand candidate(cand);
            
            // select particles based on MC info
            if (cand->getGeantPID()==14){
                Double_t mom = cand->P();
                Double_t errors[];
                cov.setSetup("FwDet");
                cov.setMomDepErrors(true);
                cov.etstimateCov(fPid[it], mom, errors);
                FillDataFw(cand, candidate, errors, 938.272);
                fCandsFit[it that is 14].push_back(candidate);
            }
            else continue;
        }
    } // end fwTrack loop*/
}

void HDSTFitter::addBuilderTask(TString val, std::vector<Int_t> pids, TLorentzVector lv = TLorentzVector()){

    //fCandsFit.Clear();
    selectCandidates();

    //initialize DecayBuilder

    //Write output category

     //end of the event loop
}

//return bool?
void HDSTFitter::addFitterTask(TString task, std::vector<Int_t> pids, TLorentzVector lv, HRefitCand mother, Double_t mm)
{

    cout << "Task added" << endl;
    TFile *outfile = new TFile("test_userfit.root", "recreate");

    TH1F *hmLam_prefit = new TH1F("hLambdaMassPreFit", "", 100, 1070, 1170);
    hmLam_prefit->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    hmLam_prefit->SetYTitle(" events ");
    hmLam_prefit->GetXaxis()->SetTitleSize(0.05);
    hmLam_prefit->GetXaxis()->SetLabelSize(0.05);
    hmLam_prefit->GetYaxis()->SetTitleSize(0.05);
    hmLam_prefit->GetYaxis()->SetLabelSize(0.05);
    hmLam_prefit->SetLineColor(kBlack);
    TH1F *hmLam_post4C = (TH1F *)hmLam_prefit->Clone("hmLam_post4C");
    hmLam_post4C->SetLineColor(kBlue);

    setPids(pids);

    TStopwatch timer;
    timer.Start();

    HLoop loop(kTRUE);
    Bool_t ret = loop.addFiles(fInfileList);
    if (ret == 0)
    {
        cout << "READBACK: ERROR : cannot find inputfiles : "
             << fInfileList.Data() << endl;
        //return 1;
    }

    // select categories here
    if (fIncludeFw)
    {
        if (!loop.setInput("-*,+HParticleCandSim,+HFwDetCandSim"))
        {
            cout << "READBACK: ERROR : cannot read input !" << endl;
            exit(1);
        } // read all categories

        fcatFwParticle = loop.getCategory("HFwDetCandSim");
        if (!fcatFwParticle)
        {
            std::cout << "No FWparticleCat in input!" << std::endl;
            exit(1);
        }
    }
    else
    {
        if (!loop.setInput("-*,+HParticleCandSim"))
        {
            cout << "READBACK: ERROR : cannot read input !" << endl;
            exit(1);
        } // read all categories
    }

    loop.printCategories();
    loop.printChain();

    fcatParticle = loop.getCategory("HParticleCandSim");
    if (!fcatParticle)
    {
        std::cout << "No particleCat in input!" << std::endl;
        exit(1);
    }

    Int_t entries = loop.getEntries();
    if (fEvents > entries || fEvents <= 0)
        fEvents = entries;

    cout << "events: " << fEvents << endl;

    // start of the event loop
    for (Int_t i = 1; i < fEvents; i++)
    {
        //----------break if last event is reached-------------
        if (loop.nextEvent(i) <= 0)
        {
            cout << " end recieved " << endl;
            break;
        } // last event reached
        HTool::printProgress(i, fEvents, 1, "Analysing evt# :");

        //fCandsFit.Clear();
        //fCandsFit = {};
        selectCandidates();

        //if not all particles are found, skip event
        bool isIncomplete = false;
        for (size_t it = 0; it < fPids.size(); it++)
        {
            cout << "fCandsFit size " << fCandsFit[it].size() << endl;
            if(fCandsFit[it].size()==0){ 
				isIncomplete = true;
				break;
			}
        }
        if (isIncomplete) continue;
        
        //initialize DecayBuilder
        cout << "ini Decay Builder" << endl;
        HDecayBuilder builder(fCandsFit, task, fPids, lv, mother, mm);
        cout << "build decay" << endl;
        builder.buildDecay();
        std::vector<HRefitCand> result;
        cout << "get result" << endl;
        builder.getFitCands(result);
        cout << result.size() << endl;
        cout << "fill histos" << endl;
        if (result.size() > 2)
        {
            hmLam_prefit->Fill((result[2] + result[3]).M());
            hmLam_post4C->Fill((result[2] + result[3]).M());
        }

        //Get output particles
    } // end of event loop

    cout << "write output file" << endl;
    outfile->cd();
    hmLam_prefit->Write();
    hmLam_post4C->Write();
    outfile->Close();

    //Write output category
}

void HDSTFitter::FillData(HParticleCandSim* cand, HRefitCand &outcand, double arr[5], double mass)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- HDSTFitter::FillData() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    double deg2rad = TMath::DegToRad();

    TMatrixD cov(5, 5);
    cov(0, 0) = std::pow(arr[0], 2);
    cov(1, 1) = std::pow(arr[1], 2);
    cov(2, 2) = std::pow(arr[2], 2);
    cov(3, 3) = std::pow(arr[3], 2);
    cov(4, 4) = std::pow(arr[4], 2);

    outcand.SetXYZM(cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::cos(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::sin(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::cos(cand->getTheta() * deg2rad),
                    mass);
    outcand.setR(cand->getR());
    outcand.setZ(cand->getZ());
    outcand.setCovariance(cov);   
}
/*
void HDSTFitter::FillDataFw(HFwDetCandSim* cand, HRefitCand& outcand, double arr[], double mass)
{
    //double deg2rad = TMath::DegToRad();

    TMatrixD cov(5, 5);
    cov(0, 0) = std::pow(arr[0], 2);
    cov(1, 1) = std::pow(arr[1], 2);
    cov(2, 2) = std::pow(arr[2], 2);
    cov(3, 3) = std::pow(arr[3], 2);
    cov(4, 4) = std::pow(arr[4], 2);
    
    cand -> calc4vectorProperties(mass);
    
    outcand.SetXYZM(cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::cos(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                        std::sin(cand->getPhi() * deg2rad),
                    cand->getMomentum() * std::cos(cand->getTheta() * deg2rad),
                    mass);
    outcand.setR(cand->getR());
    outcand.setZ(cand->getZ());
    outcand.setIsForward(true);
    outcand.setCovariance( cov );
}*/
