//****************************************************************************
//*                    This file is part of KinFit.                          *
//*                                                                          *
//*            	KinFit is distributed under the terms of the                 *
//*              GNU General Public License (GPL) version 3,                 *
//*                 copied verbatim in the file "LICENSE".                   *
//*                                                                          *
//*  				Copyright 2024                               *
//*		GSI Helmholtzzentrum für Schwerionenforschung                *
//* 	     This software is distributed under the terms of the             *
//*	     GNU General Public Licence version 3 (GPL Version 3)            *
//*		      			     				     *
//*     The copyright holders are listed in the file "COPYRIGHTHOLDERS".     *
//*               The authors are listed in the file "AUTHORS".              *
//****************************************************************************

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "KinFitter.h"

using namespace std;

void FillData(KFitParticle& outcand, double arr[])
{
    double deg2rad = TMath::DegToRad();

    TMatrixD cov(5, 5);
    cov(0, 0) = std::pow(arr[0], 2);
    cov(1, 1) = std::pow(arr[1], 2);
    cov(2, 2) = std::pow(arr[2], 2);
    cov(3, 3) = std::pow(arr[3], 2);
    cov(4, 4) = std::pow(arr[4], 2);

    outcand.setCovariance(cov);
}

Int_t fit_toyMC_fromPluto(TString infile, Int_t nEvents)
{
    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile* outfile = new TFile("testFit_toyMC_fromPluto_vtxfit.root", "recreate");

    // create histograms
    TH1F* h01 = new TH1F("hLambdaMassPreFit", "", 100, 1.070, 1.250);
    h01->SetXTitle(" M_{p#pi^{-}} [GeV/c^{2}]");
    h01->SetYTitle(" events ");
    TH1F* h012 = new TH1F("hLambdaMomentumPreFit", "", 100, 2, 4);
    h01->SetXTitle(" p_{#Lambda} [MGV/c]");
    h01->SetYTitle(" events ");

    TH1F* h02 = new TH1F("hChi2", "", 100, 0, 100);
    h02->SetXTitle("#chi^{2}");
    h02->SetYTitle(" counts ");
    TH1F *h022 = (TH1F*)h02->Clone("hChi2_probCut");
    h022 -> SetLineColor(kGreen);
    TH1F *h023 = (TH1F*)h02->Clone("hChi2_converged");
    h023 -> SetLineColor(kBlue);

    TH1F* h03 = new TH1F("hPChi2", "", 100, 0, 1);
    h03->SetXTitle("P(#chi^{2})");
    h03->SetYTitle(" counts ");

    TH1F* h04 = new TH1F("hLambdaMassPostFit", "", 100, 1.070, 1.250);
    h04->SetXTitle(" M_{p#pi^{-}} [GeV/c^{2}]");
    h04->SetYTitle(" events ");
    h04 -> SetLineColor(kRed);
    TH1F *h042 = (TH1F*)h04->Clone("hLambdaMassPostFit_probCut");
    h042 -> SetLineColor(kGreen);
    TH1F *h043 = (TH1F*)h04->Clone("hLambdaMassPostFit_converged");
    h043 -> SetLineColor(kBlue);

    TH1F* h09 = new TH1F("hLambdaMomPostFit", "", 100, 2, 4);
    h04->SetXTitle(" M_{p#pi^{-}} [GeV/c]");
    h04->SetYTitle(" events ");
    h04 -> SetLineColor(kRed);
    TH1F *h092 = (TH1F*)h09->Clone("hLambdaMomPostFit_probCut");
    h092 -> SetLineColor(kGreen);

    TH1F* h05 = new TH1F("hPull", "", 100, -5, 5);
    h05->SetXTitle("Pull(1/P_{p})");
    h05->SetYTitle(" counts ");
    TH1F *h052 = (TH1F*)h05->Clone("hPull_probCut");
    h052 -> SetLineColor(kGreen);
    TH1F *h053 = (TH1F*)h05->Clone("hPull_converged");
    h053 -> SetLineColor(kBlue);
    TH1F *h054 = (TH1F*)h05->Clone("hPull_tht");
    TH1F *h055 = (TH1F*)h05->Clone("hPull_phi");
    TH1F *h056 = (TH1F*)h05->Clone("hPull_p_pi");
    TH1F *h057 = (TH1F*)h05->Clone("hPull_tht_pi");
    TH1F *h058 = (TH1F*)h05->Clone("hPull_phi_pi");

    TH1F* h08 = new TH1F("hNIterations", "", 20, 0, 20);
    h08->SetXTitle(" Iteration");
    h08->SetYTitle(" events ");
    
    
    // -----------------------------------------------------------------------
    
    
   Float_t p1CandTrueP, p1CandTrueTheta, p1CandTruePhi, p1CandTrueR, p1CandTrueZ, p1CandRecoP, p1CandRecoTheta, p1CandRecoPhi, p1CandRecoR, p1CandRecoZ,
            KCandTrueP, KCandTrueTheta, KCandTruePhi, KCandTrueR, KCandTrueZ, KCandRecoP, KCandRecoTheta, KCandRecoPhi,  KCandRecoR, KCandRecoZ,
            p2CandTrueP, p2CandTrueTheta, p2CandTruePhi, p2CandTrueR, p2CandTrueZ, p2CandRecoP, p2CandRecoTheta, p2CandRecoPhi,  p2CandRecoR, p2CandRecoZ,
            piCandTrueP, piCandTrueTheta, piCandTruePhi, piCandTrueR, piCandTrueZ, piCandRecoP, piCandRecoTheta, piCandRecoPhi,  piCandRecoR, piCandRecoZ;
    
    TFile tree_file(infile, "READ");
    TTree *t = (TTree*)tree_file.Get("data");

    t->SetBranchAddress("p1CandTrueP", &p1CandTrueP);
    t->SetBranchAddress("p1CandRecoP", &p1CandRecoP);
    t->SetBranchAddress("p1CandTrueTheta", &p1CandTrueTheta);
    t->SetBranchAddress("p1CandRecoTheta", &p1CandRecoTheta);
    t->SetBranchAddress("p1CandTruePhi", &p1CandTruePhi);
    t->SetBranchAddress("p1CandRecoPhi", &p1CandRecoPhi);
    t->SetBranchAddress("p1CandTrueR", &p1CandTrueR);
    t->SetBranchAddress("p1CandRecoR", &p1CandRecoR);
    t->SetBranchAddress("p1CandTrueZ", &p1CandTrueZ);
    t->SetBranchAddress("p1CandRecoZ", &p1CandRecoZ);
    t->SetBranchAddress("KCandTrueP", &KCandTrueP);
    t->SetBranchAddress("KCandRecoP", &KCandRecoP);
    t->SetBranchAddress("KCandTrueTheta", &KCandTrueTheta);
    t->SetBranchAddress("KCandRecoTheta", &KCandRecoTheta);
    t->SetBranchAddress("KCandTruePhi", &KCandTruePhi);
    t->SetBranchAddress("KCandRecoPhi", &KCandRecoPhi);
    t->SetBranchAddress("KCandTrueR", &KCandTrueR);
    t->SetBranchAddress("KCandRecoR", &KCandRecoR);
    t->SetBranchAddress("KCandTrueZ", &KCandTrueZ);
    t->SetBranchAddress("KCandRecoZ", &KCandRecoZ);
    t->SetBranchAddress("p2CandTrueP", &p2CandTrueP);
    t->SetBranchAddress("p2CandRecoP", &p2CandRecoP);
    t->SetBranchAddress("p2CandTrueTheta", &p2CandTrueTheta);
    t->SetBranchAddress("p2CandRecoTheta", &p2CandRecoTheta);
    t->SetBranchAddress("p2CandTruePhi", &p2CandTruePhi);
    t->SetBranchAddress("p2CandRecoPhi", &p2CandRecoPhi);
    t->SetBranchAddress("p2CandTrueR", &p2CandTrueR);
    t->SetBranchAddress("p2CandRecoR", &p2CandRecoR);
    t->SetBranchAddress("p2CandTrueZ", &p2CandTrueZ);
    t->SetBranchAddress("p2CandRecoZ", &p2CandRecoZ);
    t->SetBranchAddress("piCandTrueP", &piCandTrueP);
    t->SetBranchAddress("piCandRecoP", &piCandRecoP);
    t->SetBranchAddress("piCandTrueTheta", &piCandTrueTheta);
    t->SetBranchAddress("piCandRecoTheta", &piCandRecoTheta);
    t->SetBranchAddress("piCandTruePhi", &piCandTruePhi);
    t->SetBranchAddress("piCandRecoPhi", &piCandRecoPhi);
    t->SetBranchAddress("piCandTrueR", &piCandTrueR);
    t->SetBranchAddress("piCandRecoR", &piCandRecoR);
    t->SetBranchAddress("piCandTrueZ", &piCandTrueZ);
    t->SetBranchAddress("piCandRecoZ", &piCandRecoZ);

    Double_t mp = 0.9382720813;
    Double_t p01 = sqrt(pow((4.500+mp),2)-pow(mp,2));
    TLorentzVector ini(0, 0, p01, 4.500+2*mp);  //M=3459

    Long64_t nevts = t->GetEntries();
    if (nEvents <= 0 || nEvents > nevts)
        nEvents = nevts;
    // event loop
    for (Long64_t ev = 0; ev < nEvents; ev++)
    {
        t->GetEntry(ev);

        TLorentzVector *proton1 = new TLorentzVector();
        proton1->SetXYZM(p1CandRecoP * std::sin(p1CandRecoTheta) * std::cos(p1CandRecoPhi),
                    p1CandRecoP * std::sin(p1CandRecoTheta) * std::sin(p1CandRecoPhi),
                    p1CandRecoP * std::cos(p1CandRecoTheta), 0.938272);
        double proton1_errors[] = {0.025*(1/p1CandRecoP), 0.0009, 0.0009, 0.5, 1.};

        TLorentzVector *kaon = new TLorentzVector();
        kaon->SetXYZM(KCandRecoP * std::sin(KCandRecoTheta) * std::cos(KCandRecoPhi),
                    KCandRecoP * std::sin(KCandRecoTheta) * std::sin(KCandRecoPhi),
                    KCandRecoP * std::cos(KCandRecoTheta), 0.493677);
        double kaon_errors[] = {0.025*(1/KCandRecoP), 0.0009, 0.0009, 0.5, 1.};

        TLorentzVector *proton2 = new TLorentzVector();
        proton2->SetXYZM(p2CandRecoP * std::sin(p2CandRecoTheta) * std::cos(p2CandRecoPhi),
                    p2CandRecoP * std::sin(p2CandRecoTheta) * std::sin(p2CandRecoPhi),
                    p2CandRecoP * std::cos(p2CandRecoTheta), 0.938272);
        double proton2_errors[] = {0.025*(1/p2CandRecoP), 0.0009, 0.0009, 0.5, 1.};

        TLorentzVector *pion = new TLorentzVector();
        pion->SetXYZM(piCandRecoP * std::sin(piCandRecoTheta) * std::cos(piCandRecoPhi),
                    piCandRecoP * std::sin(piCandRecoTheta) * std::sin(piCandRecoPhi),
                    piCandRecoP * std::cos(piCandRecoTheta), 0.13957);
        double pion_errors[] = {0.025*(1/piCandRecoP), 0.0009, 0.0009, 0.5, 1.};

        TLorentzVector lambda = *proton2 + *pion;
        h01->Fill(lambda.M());
        h012->Fill(lambda.P());

        KFitParticle proton1_fit(*proton1,p1CandRecoR,p1CandRecoZ);
        FillData(proton1_fit, proton1_errors);
        KFitParticle kaon_fit(*kaon,KCandRecoR,KCandRecoZ);
        FillData(kaon_fit, kaon_errors);
        KFitParticle proton2_fit(*proton2,p2CandRecoR,p2CandRecoZ);
        FillData(proton2_fit, proton2_errors);
        KFitParticle pion_fit(*pion,piCandRecoR, piCandRecoZ);
        FillData(pion_fit, pion_errors);

        // ---------------------------------------------------------------------------------
        // begin kinfit here
        // ---------------------------------------------------------------------------------
        std::vector<KFitParticle> cands;
        cands.clear();
        cands.push_back(proton2_fit);
        cands.push_back(pion_fit);

        KinFitter fitter(cands);
        fitter.setVerbosity(0);
        fitter.setNumberOfIterations(20);
        fitter.setConvergenceCriteria(0.01, 1e5, 1e6);
        fitter.addVertexConstraint();
        if(fitter.fit()){

            KFitParticle fcand1 = fitter.getDaughter(0); // proton
            KFitParticle fcand2 = fitter.getDaughter(1); // pion

            TLorentzVector daughter = fitter.getMissingDaughter();
            cout<<daughter.Px()<<endl;

            h02->Fill(fitter.getChi2());
            h03->Fill(fitter.getProb());
            TLorentzVector lambda_fit = fcand1 + fcand2;
            h04->Fill(lambda_fit.M());
            // get Pull example (1/P for the fitted proton)
            h05->Fill(fitter.getPull(0));
            h054->Fill(fitter.getPull(1));
            h055->Fill(fitter.getPull(2));
            h056->Fill(fitter.getPull(5));
            h057->Fill(fitter.getPull(6));
            h058->Fill(fitter.getPull(7));
            h09->Fill(lambda_fit.P());
            
            if(fitter.getProb()>0.01){
                h022->Fill(fitter.getChi2());
                h042->Fill(lambda_fit.M());

                // get Pull example (1/P for the fitted proton)
                h052->Fill(fitter.getPull(0));
                h092->Fill(lambda_fit.P());
            }

            h08->Fill(fitter.getIteration());

        }
    }

    // write histograms to the output file
    outfile->cd();
    h01->Write();
    h012->Write();
    h02->Write();
    h022->Write();
    h023->Write();
    h03->Write();
    h04->Write();
    h042->Write();
    h043->Write();
    h05->Write();
    h052->Write();
    h053->Write();
    h054->Write();
    h055->Write();
    h056->Write();
    h057->Write();
    h058->Write();
    h08->Write();
    h09->Write();
    h092->Write();
    outfile->Close();

    return 0;
}