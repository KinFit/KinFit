//****************************************************************************
//*                    This file is part of KinFit.                          *
//*                                                                          *
//*            	KinFit is distributed under the terms of the                 *
//*              GNU General Public License (GPL) version 3,                 *
//*                 copied verbatim in the file "LICENSE".                   *
//*                                                                          *
//*  				           Copyright 2024                                *
//*		         GSI Helmholtzzentrum für Schwerionenforschung               *
//* 	      This software is distributed under the terms of the            *
//*	          GNU General Public Licence version 3 (GPL Version 3)           *
//*		      			     				                                 *
//*     The copyright holders are listed in the file "COPYRIGHTHOLDERS".     *
//*               The authors are listed in the file "AUTHORS".              *
//****************************************************************************

#include "KFitAnalyzer.h"

KFitAnalyzer::KFitAnalyzer(TString inFileName, TString outFileName, int nEvents) : fEvents(nEvents),
                                                                                    fVerbose(0)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitAnalyzer() -----------------" << std::endl;
    }

    TFile *inFile = new TFile(inFileName, "READ");
    fTree = (TTree*)inFile->Get("data");

    foutFile = new TFile(outFileName,"RECREATE");
    fTree_out = new TTree("data_fitted", "output_tree");

}

// Select and sort particles according to their PID
void KFitAnalyzer::selectCandidates()
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- KFitAnalyzer::selectCandidates() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    int ntracks = fCands_in->GetEntries();

    fCandsFit.clear();
    std::vector<KFitParticle > tempVec;

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

void KFitAnalyzer::doFitterTask(TString task, std::vector<int> pids, double mass, TLorentzVector lv, KFitParticle mother)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- KFitAnalyzer::selectCandidates() -----------" << std::endl;
        std::cout << "" << std::endl;
        std::cout << "Task added: " << task << std::endl;
        std::cout << "" << std::endl;
    }

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
    if (fVerbose > 1 )
    {
        std::cout<<"EntEvents to be analyzed: "<<entries<<std::endl;
    }
    
    KFitDecayBuilder builder(task, fPids, lv, mass);
    builder.setVerbosity(fVerbose);

    // start of the event loop
    if (fVerbose > 1 )
    {
        std::cout<<"Start of event loop"<<std::endl;
    }
    for (int i = 1; i < entries; i++)
    {
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

        // get result from DecayBuilder
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
    
    if (fVerbose > 0 )
    {
        std::cout<<"Event loop finished!"<<std::endl;
    }
    finish();
}

// Close everything
void KFitAnalyzer::finish()
{
    foutFile->cd();
    fTree_out->Write();
    foutFile->Save();
    foutFile->Close();

}