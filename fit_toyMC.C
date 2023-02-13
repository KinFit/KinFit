#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "/home/jana/KinFit/include/KinFitter.h"

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

Int_t fit_toyMC(TString infile, Int_t nEvents)
{
    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile* outfile = new TFile("testFit_toyMC_M.root", "recreate");
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
    /*
    TH1F* h06 = new TH1F("hTotMomPreFit", "", 100, 3800, 4800);
    h06->SetXTitle(" p [MeV/c]");
    h06->SetYTitle(" events ");
    
    TH1F* h07 = new TH1F("hTotMomPostFit", "", 100, 3800, 4800);
    h07->SetXTitle(" p [MeV/c]");
    h07->SetYTitle(" events ");
    h07 -> SetLineColor(kRed);
    TH1F *h072 = (TH1F*)h07->Clone("hTotMomPostFit_probCut");
    h072 -> SetLineColor(kGreen);
    TH1F *h073 = (TH1F*)h07->Clone("hTotMomPostFit_converged");
    h073 -> SetLineColor(kBlue);
    */
    TH1F* h08 = new TH1F("hNIterations", "", 10, 0, 10);
    h08->SetXTitle(" Iteration");
    h08->SetYTitle(" events ");
    
    
    // -----------------------------------------------------------------------
    
    Float_t pCandTrueP, pCandRecoP, pCandTrueTheta, pCandRecoTheta, pCandTruePhi, pCandRecoPhi, 
                piCandTrueP, piCandRecoP, piCandTrueTheta, piCandRecoTheta, piCandTruePhi, piCandRecoPhi;
    
    TFile tree_file(infile, "READ");
    TTree *t = (TTree*)tree_file.Get("data");

    t->SetBranchAddress("pCandTrueP", &pCandTrueP);
    t->SetBranchAddress("pCandRecoP", &pCandRecoP);
    t->SetBranchAddress("pCandTrueTheta", &pCandTrueTheta);
    t->SetBranchAddress("pCandRecoTheta", &pCandRecoTheta);
    t->SetBranchAddress("pCandTruePhi", &pCandTruePhi);
    t->SetBranchAddress("pCandRecoPhi", &pCandRecoPhi);
    t->SetBranchAddress("piCandTrueP", &piCandTrueP);
    t->SetBranchAddress("piCandRecoP", &piCandRecoP);
    t->SetBranchAddress("piCandTrueTheta", &piCandTrueTheta);
    t->SetBranchAddress("piCandRecoTheta", &piCandRecoTheta);
    t->SetBranchAddress("piCandTruePhi", &piCandTruePhi);
    t->SetBranchAddress("piCandRecoPhi", &piCandRecoPhi);

    TLorentzVector ini;
    ini.SetXYZM(0.1, .1, 3., 1.11568);

    Long64_t nevts = t->GetEntries();
    if (nEvents <= 0 || nEvents > nevts)
        nEvents = nevts;
    // event loop
    for (Long64_t ev = 0; ev < nEvents; ev++)
    {
        t->GetEntry(ev);

        TLorentzVector *proton = new TLorentzVector();
        proton->SetXYZM(pCandRecoP * std::sin(pCandRecoTheta) * std::cos(pCandRecoPhi),
                    pCandRecoP * std::sin(pCandRecoTheta) * std::sin(pCandRecoPhi),
                    pCandRecoP * std::cos(pCandRecoTheta), 0.938272);
        double proton_errors[] = {0.025*(1/pCandRecoP), 0.0009, 0.0009, 0.0001, 0.0001};

        TLorentzVector *pion = new TLorentzVector();
        pion->SetXYZM(piCandRecoP * std::sin(piCandRecoTheta) * std::cos(piCandRecoPhi),
                    piCandRecoP * std::sin(piCandRecoTheta) * std::sin(piCandRecoPhi),
                    piCandRecoP * std::cos(piCandRecoTheta), 0.13957);
        double pion_errors[] = {0.025*(1/piCandRecoP), 0.0009, 0.0009, 0.0001, 0.0001};

        TLorentzVector lambda = *proton + *pion;
        h01->Fill(lambda.M());
        h012->Fill(lambda.P());

        KFitParticle proton_fit(proton,0,0);
        FillData(proton_fit, proton_errors);
        KFitParticle pion_fit(pion,0,0);
        FillData(pion_fit, pion_errors);

        // ---------------------------------------------------------------------------------
        // begin kinfit here
        // ---------------------------------------------------------------------------------
        std::vector<KFitParticle> cands;
        cands.clear();
        cands.push_back(proton_fit);
        cands.push_back(pion_fit);

        KinFitter fitter(cands);
        fitter.setVerbosity(0);
        fitter.setNumberOfIterations(10);
        //fitter.setLearningRate(0.5);
        //fitter.setConvergenceCriterion(0.01);
        fitter.addMassConstraint(1.11568);
        //fitter.add4Constraint(ini);
        if(fitter.fit()){

            KFitParticle fcand1 = fitter.getDaughter(0); // proton
            KFitParticle fcand2 = fitter.getDaughter(1); // pion

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
               // h072->Fill(all_fit.P());

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
    /*
    h06->Write();
    h07->Write();
    h072->Write();
    h073->Write();
    */
    h08->Write();
    h09->Write();
    h092->Write();
    outfile->Close();

    return 0;
}