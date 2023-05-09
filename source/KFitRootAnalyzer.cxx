#include "KFitRootAnalyzer.h"

KFitRootAnalyzer::KFitRootAnalyzer(TString inFileName, TString outFileName, int nEvents) : fEvents(nEvents),
                                                                                    fVerbose(0)
{
    finFile = new TFile(inFileName, "READ");
    fTree = (TTree*)finFile->Get("data");

    //foutFile = new TFile(outFileName.Data(),"RECREATE");
    foutFile = new TFile(outFileName,"RECREATE");
    fTree_out = new TTree("data_fitted", "output_tree");

}

// Select and sort particles according to their PID
void KFitRootAnalyzer::selectCandidates()
{
    int ntracks = fCands_in->GetEntries();

    fCandsFit.clear();
    std::vector<KFitParticle > tempVec;
    //KFitParticle *cand = new KFitParticle();

    for (size_t it = 0; it < fPids.size(); it++)
    {
        tempVec.clear();
        for (int j = 0; j < ntracks; j++)
        {
            KFitParticle *cand = (KFitParticle *)fCands_in->At(j);

            if (cand->getPid() == fPids[it])
            {
                tempVec.push_back(*cand);
            }
        }
        fCandsFit.push_back(tempVec);

    } // end of PIDs loop
}

void KFitRootAnalyzer::doFitterTask(TString task, std::vector<int> pids, double mm, TLorentzVector lv, KFitParticle mother)
{
    cout << "Task added: " << task << endl;

    // Read input tree
    int Event;
    fTree->SetBranchAddress("KFitParticle", &fCands_in);
    
    // Create output tree
    fTree_out->SetName("data_fitted");
    TString out_branchname = "cands_fitted";
    TClonesArray *fitted_cands = new TClonesArray("KFitParticle");
    TClonesArray &fit_arrayRef = *fitted_cands;
    double Chi2, Prob;

    fTree_out->Branch("Event", &Event, "Event/I");
    fTree_out->Branch("KFitParticle", "TClonesArray", &fitted_cands);
    fTree_out->Branch("Chi2", &Chi2, "Chi2/D");
    fTree_out->Branch("Prob", &Prob, "Prob/D");

    setPids(pids);

    int entries = fTree->GetEntries();
    if (fEvents < entries && fEvents > 0)
        entries = fEvents;

    cout<<"entries: "<<entries<<endl;

    KFitDecayBuilder builder(task, fPids, lv, mother, mm);

    // start of the event loop
    for (int i = 1; i < entries; i++)
    {
        //cout<<"Event: "<<i<<endl;
        fTree->GetEntry(i);
        Event = i;

        selectCandidates();

        // if not all particles are found, skip event
        bool isIncomplete = false;
        for (size_t it = 0; it < fPids.size(); it++)
        {
            if (fCandsFit[it].size() == 0)
            {
                isIncomplete = true;
                break;
            }
        }
        if (isIncomplete)
            continue;

        // initialize DecayBuilder
        builder.setInputCands(fCandsFit);
        builder.countCombis();
        
        builder.buildDecay();

        std::vector<KFitParticle> result;
        if(result.size()>0) result.clear();
        builder.getFitCands(result);
        Chi2 = builder.getChi2();
        Prob = builder.getProbability();

        // Fill output TClonesArray with fit result
        int ii = 0;
        if(fitted_cands->GetEntries()>0) fitted_cands->Clear();
        if(Prob>0){
        for (int k = 0; k < result.size(); k++)
        {
            KFitParticle *fitted_cand = new (fit_arrayRef[ii]) KFitParticle(result[k], result[k].getR(), result[k].getZ());
            fitted_cand->setPid(result[k].getPid());
            fitted_cand->setTrackId(result[k].getTrackId());
            fitted_cand->setCovariance(result[k].getCovariance());

            ii++;
        }

        fTree_out->Fill();
        fitted_cands->Clear();
        }
        else 
        {
            cout << "no candidate found" << endl;
            continue;
        }
        
    } // end of event loop
    //cout<<"Event loop ended"<<endl;
    finish();
    //cout<<"finished"<<endl;
}

// Close everything
void KFitRootAnalyzer::finish()
{
    foutFile->cd();
    fTree_out->Write();
    foutFile->Save();
    foutFile->Close();

    finFile->Close();
}

/*
// Function to fill KFitParticle -- not needed here if input is already KFitParticle
void KFitRootAnalyzer::FillData(HParticleCandSim *cand, KFitParticle &outcand, double arr[5], double mass)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- KFitRootAnalyzer::FillData() -----------" << std::endl;
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
*/