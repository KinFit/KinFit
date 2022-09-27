#include "hrootfitter.h"

HRootFitter::HRootFitter(TString inFileName, TString outFileName, Int_t nEvents) : fEvents(nEvents),
                                                                                    fVerbose(0)
{
    finFile = new TFile(inFileName, "open");
    fTree = (TTree*)finFile->Get("data");

    foutFile = new TFile(outFileName.Data(),"RECREATE");
    fTree_out = new TTree();

}

// Select and sort particles according to their PID
void HRootFitter::selectCandidates()
{
    Int_t ntracks = fCands_in->GetEntries();

    fCandsFit.clear();
    std::vector<HRefitCand > tempVec;
    HRefitCand *cand = new HRefitCand();

    for (size_t it = 0; it < fPids.size(); it++)
    {
        tempVec.clear();
        for (Int_t j = 0; j < ntracks; j++)
        {
            HRefitCand *cand = (HRefitCand *)fCands_in->At(j);

            if (cand->getPid() == fPids[it])
            {
                tempVec.push_back(*cand);
            }
        }
        fCandsFit.push_back(tempVec);

    } // end of PIDs loop
}

/*
void HRootFitter::addBuilderTask(TString val, std::vector<Int_t> pids, TLorentzVector lv = TLorentzVector())
{

    // fCandsFit.Clear();
    selectCandidates();

    // initialize DecayBuilder

    // Write output category

    // end of the event loop
}
*/

void HRootFitter::doFitterTask(TString task, std::vector<Int_t> pids, TLorentzVector lv, HRefitCand mother, Double_t mm)
{
    cout << "Task added: " << task << endl;

    // Read input tree
    //TClonesArray *input_cands = new TClonesArray("HRefitCand");
    Int_t Event;
    fTree->SetBranchAddress("cands", &fCands_in);
    
    // Create output tree
    //fTree_out = fTree->CopyTree();
    fTree_out->SetName("data_fitted");
    //fTree_out->BuildIndex("Event");
    TString out_branchname = "cands_fitted";
    TClonesArray *fitted_cands = new TClonesArray("HRefitCand");
    TClonesArray &fit_arrayRef = *fitted_cands;
    Double_t Chi2, Prob;
    //TClonesArray &fit_arrayRef = *fit_array;

    fTree_out->Branch("Event", &Event, "Event/I");
    fTree_out->Branch(out_branchname, "TClonesArray", &fitted_cands);
    fTree_out->Branch("Chi2", &Chi2, "Chi2/D");
    fTree_out->Branch("Prob", &Prob, "Prob/D");

    setPids(pids);
/*
    TFile tree_file(fInfileList, "read");
    fTree = (TTree *)tree_file.Get("data");
*/

    Int_t entries = fTree->GetEntries();
    if (fEvents < entries && fEvents > 0)
        entries = fEvents;

    // start of the event loop
    for (Int_t i = 0; i < fEvents; i++)
    {
        fTree->GetEntry(i);
        Event = i;

        selectCandidates();

        // if not all particles are found, skip event
        bool isIncomplete = false;
        for (size_t it = 0; it < fPids.size(); it++)
        {
            cout << "fCandsFit size " << fCandsFit[it].size() << endl;
            if (fCandsFit[it].size() == 0)
            {
                isIncomplete = true;
                break;
            }
        }
        if (isIncomplete)
            continue;

        // initialize DecayBuilder
        cout << "ini Decay Builder" << endl;
        HDecayBuilder builder(fCandsFit, task, fPids, lv, mother, mm);
        cout << "build decay" << endl;
        builder.buildDecay();
        std::vector<HRefitCand> result;
        cout << "get result" << endl;
        cout << result.size() << endl;
        builder.getFitCands(result);
        Chi2 = builder.getChi2();
        Prob = builder.getProbability();

        // Fill output TClonesArray with fit result
        Int_t ii = 0;
        fitted_cands->Clear();
        for (Int_t k = 0; k < result.size(); k++)
        {
            HRefitCand *fitted_cand = new (fit_arrayRef[ii]) HRefitCand();
            fitted_cand = &result[k];

            ii++;
        }
        fTree_out->Fill();
/*
        cout << "fill histos" << endl;
        if (result.size() > 2)
        {
            hmLam_prefit->Fill((result[2] + result[3]).M());
            hmLam_post4C->Fill((result[2] + result[3]).M());
        }
*/
    } // end of event loop
/*
    cout << "write output file" << endl;
    outfile->cd();
    hmLam_prefit->Write();
    hmLam_post4C->Write();
    outfile->Close();
*/
    // Write tree to outfile and close
}

// Close everything
void HRootFitter::finish()
{
    foutFile->cd();
    fTree_out->Write();
    foutFile->Save();
    foutFile->Close();

    finFile->Close();
}

/*
// Function to fill HRefitCand -- not needed here if input is already HRefitCand
void HRootFitter::FillData(HParticleCandSim *cand, HRefitCand &outcand, double arr[5], double mass)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- HRootFitter::FillData() -----------" << std::endl;
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