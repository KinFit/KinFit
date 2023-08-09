#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCutG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "/home/jana/KinFit_official/KinFit/include/KinFitter.h"

#include "plotsStyleMacro.C"

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

Int_t fit4C_toyMC_fromPluto(TString infile, Int_t nEvents)
{
    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile* outfile = new TFile("fit4C_toyMC_myRZ.root", "recreate");

    TH1F* hLambdaMassPreFit = new TH1F("hLambdaMassPreFit", "", 100, 1.070, 1.250);
    hLambdaMassPreFit->SetXTitle(" M_{p#pi^{-}} [GeV/c^{2}]");
    hLambdaMassPreFit->SetYTitle(" events ");
    TH1F* hKMomentumPreFit = new TH1F("hKMomentumPreFit", "", 100, 2, 4);
    hKMomentumPreFit->SetXTitle(" p_{K} [GeV/c]");
    hKMomentumPreFit->SetYTitle(" events ");
    TH1F* hKThetaPreFit = new TH1F("hKThetaPreFit", "", 100, 0, 3.2);
    hKThetaPreFit->SetXTitle(" #vartheta_{K} [rad]");
    hKThetaPreFit->SetYTitle(" events ");
    TH1F* hKPhiPreFit = new TH1F("hKPhiPreFit", "", 100, -3.2, 3.2);
    hKPhiPreFit->SetXTitle(" #varphi_{K} [rad]");
    hKPhiPreFit->SetYTitle(" events ");
    TH1F* hKMomentumPreFitRes = new TH1F("hKMomentumPreFitRes", "", 200, -0.1, 0.1);
    hKMomentumPreFitRes->SetXTitle(" (1/p_{gen}-1/p_{reco})/(1/p_{gen})(K)");
    hKMomentumPreFitRes->SetYTitle(" events ");
    TH1F* hKThetaPreFitRes = new TH1F("hKThetaPreFitRes", "", 200, -0.01, 0.01);
    hKThetaPreFitRes->SetXTitle(" #vartheta_{gen}-#vartheta_{reco}(K) [rad]");
    hKThetaPreFitRes->SetYTitle(" events ");
    TH1F* hKPhiPreFitRes = new TH1F("hKPhiPreFitRes", "", 200, -0.01, 0.01);
    hKPhiPreFitRes->SetXTitle(" #varphi_{gen}-#varphi_{reco}(K) [rad]");
    hKPhiPreFitRes->SetYTitle(" events ");
    
    TH1F* hPiMomentumPreFit = new TH1F("hPiMomentumPreFit", "", 100, 2, 4);
    hPiMomentumPreFit->SetXTitle(" p_{#pi^{-}} [GeV/c]");
    hPiMomentumPreFit->SetYTitle(" events ");
    TH1F* hPiThetaPreFit = new TH1F("hPiThetaPreFit", "", 100, 0, 3.2);
    hPiThetaPreFit->SetXTitle(" #vartheta_{#pi^{-}} [rad]");
    hPiThetaPreFit->SetYTitle(" events ");
    TH1F* hPiPhiPreFit = new TH1F("hPiPhiPreFit", "", 100, -3.2, 3.2);
    hPiPhiPreFit->SetXTitle(" #varphi_{#pi^{-}} [rad]");
    hPiPhiPreFit->SetYTitle(" events ");
    TH1F* hPiMomentumPreFitRes = new TH1F("hPiMomentumPreFitRes", "", 200, -0.1, 0.1);
    hPiMomentumPreFitRes->SetXTitle(" (1/p_{gen}-1/p_{reco})/(1/p_{gen})(#pi^{-})");
    hPiMomentumPreFitRes->SetYTitle(" events ");
    TH1F* hPiThetaPreFitRes = new TH1F("hPiThetaPreFitRes", "", 200, -0.01, 0.01);
    hPiThetaPreFitRes->SetXTitle(" #vartheta_{gen}-#vartheta_{reco}(#pi^{-}) [rad]");
    hPiThetaPreFitRes->SetYTitle(" events ");
    TH1F* hPiPhiPreFitRes = new TH1F("hPiPhiPreFitRes", "", 200, -0.01, 0.01);
    hPiPhiPreFitRes->SetXTitle(" #varphi_{gen}-#varphi_{reco}(#pi^{-}) [rad]");
    hPiPhiPreFitRes->SetYTitle(" events ");

    TH1F* hChi2 = new TH1F("hChi2", "", 100, 0, 100);
    hChi2->SetXTitle("#chi^{2}");
    hChi2->SetYTitle(" counts ");
    TH1F *hChi2_probCut = (TH1F*)hChi2->Clone("hChi2_probCut");
    hChi2_probCut -> SetLineColor(kGreen);
    TH1F *hChi2_converged = (TH1F*)hChi2->Clone("hChi2_converged");
    hChi2_converged -> SetLineColor(kBlue);

    TH1F* hPChi2 = new TH1F("hPChi2", "", 100, 0, 1);
    hPChi2->SetXTitle("P(#chi^{2})");
    hPChi2->SetYTitle(" counts ");
    hPChi2->SetMinimum(0);

    TH1F* hLambdaMassPostFit = new TH1F("hLambdaMassPostFit", "", 100, 1.070, 1.250);
    hLambdaMassPostFit->SetXTitle(" M_{p#pi^{-}} [GeV/c^{2}]");
    hLambdaMassPostFit->SetYTitle(" events ");
    hLambdaMassPostFit -> SetLineColor(kRed);
    TH1F *hLambdaMassPostFit_probCut = (TH1F*)hLambdaMassPostFit->Clone("hLambdaMassPostFit_probCut");
    hLambdaMassPostFit_probCut -> SetLineColor(kGreen);
    TH1F *hLambdaMassPostFit_converged = (TH1F*)hLambdaMassPostFit->Clone("hLambdaMassPostFit_converged");
    hLambdaMassPostFit_converged -> SetLineColor(kBlue);

    TH1F* hKMomPostFit = new TH1F("hKMomPostFit", "", 100, 2, 4);
    hKMomPostFit->SetXTitle(" P_{K} [GeV/c]");
    hKMomPostFit->SetYTitle(" events ");
    hKMomPostFit -> SetLineColor(kRed);
    TH1F *hKMomPostFit_probCut = (TH1F*)hKMomPostFit->Clone("hKMomPostFit_probCut");
    hKMomPostFit_probCut -> SetLineColor(kGreen);
    TH1F* hKThetaPostFit = new TH1F("hKThetaPostFit", "", 100, 0, 3.2);
    hKThetaPostFit->SetXTitle(" #vartheta_{K} [rad]");
    hKThetaPostFit->SetYTitle(" events ");
    TH1F* hKPhiPostFit = new TH1F("hKPhiPostFit", "", 100, -3.2, 3.2);
    hKPhiPostFit->SetXTitle(" #varphi_{K} [rad]");
    hKPhiPostFit->SetYTitle(" events ");
    TH1F* hKMomPostFitRes = new TH1F("hKMomPostFitRes", "", 200, -0.1, 0.1);
    hKMomPostFitRes->SetXTitle(" (1/p_{gen}-1/p_{reco})/(1/p_{gen})(K)");
    hKMomPostFitRes->SetYTitle(" events ");
    TH1F* hKThetaPostFitRes = new TH1F("hKThetaPostFitRes", "", 200, -0.01, 0.01);
    hKThetaPostFitRes->SetXTitle(" #vartheta_{gen}-#vartheta_{reco}(K) [rad]");
    hKThetaPostFitRes->SetYTitle(" events ");
    TH1F* hKPhiPostFitRes = new TH1F("hKPhiPostFitRes", "", 200, -0.01, 0.01);
    hKPhiPostFitRes->SetXTitle(" #varphi_{gen}-#varphi_{reco}(K) [rad]");
    hKPhiPostFitRes->SetYTitle(" events ");

    TH1F* hPiMomPostFit = new TH1F("hPiMomPostFit", "", 100, 2, 4);
    hPiMomPostFit->SetXTitle(" P_{#pi^{-}} [GeV/c]");
    hPiMomPostFit->SetYTitle(" events ");
    hPiMomPostFit -> SetLineColor(kRed);
    TH1F *hPiMomPostFit_probCut = (TH1F*)hPiMomPostFit->Clone("hPiMomPostFit_probCut");
    hPiMomPostFit_probCut -> SetLineColor(kGreen);
    TH1F* hPiThetaPostFit = new TH1F("hPiThetaPostFit", "", 100, 0, 3.2);
    hPiThetaPostFit->SetXTitle(" #vartheta_{#pi^{-}} [rad]");
    hPiThetaPostFit->SetYTitle(" events ");
    TH1F* hPiPhiPostFit = new TH1F("hPiPhiPostFit", "", 100, -3.2, 3.2);
    hPiPhiPostFit->SetXTitle(" #varphi_{#pi^{-}} [rad]");
    hPiPhiPostFit->SetYTitle(" events ");
    TH1F* hPiMomPostFitRes = new TH1F("hPiMomPostFitRes", "", 200, -0.1, 0.1);
    hPiMomPostFitRes->SetXTitle(" (1/p_{gen}-1/p_{reco})/(1/p_{gen})(#pi^{-})");
    hPiMomPostFitRes->SetYTitle(" events ");
    TH1F* hPiThetaPostFitRes = new TH1F("hPiThetaPostFitRes", "", 200, -0.01, 0.01);
    hPiThetaPostFitRes->SetXTitle(" #vartheta_{gen}-#vartheta_{reco}(#pi^{-}) [rad]");
    hPiThetaPostFitRes->SetYTitle(" events ");
    TH1F* hPiPhiPostFitRes = new TH1F("hPiPhiPostFitRes", "", 200, -0.01, 0.01);
    hPiPhiPostFitRes->SetXTitle(" #varphi_{gen}-#varphi_{reco}(#pi^{-}) [rad]");
    hPiPhiPostFitRes->SetYTitle(" events ");

    TH1F* hPull = new TH1F("hPull", "", 100, -5, 5);
    hPull->SetXTitle("Pull(1/P_{p2})");
    hPull->SetYTitle(" counts ");
    TH1F *hPull_probCut = (TH1F*)hPull->Clone("hPull_probCut");
    hPull_probCut -> SetLineColor(kGreen);
    TH1F *hPull_converged = (TH1F*)hPull->Clone("hPull_converged");
    hPull_converged -> SetLineColor(kBlue);
    TH1F *hPull_tht = (TH1F*)hPull->Clone("hPull_tht");
    hPull_tht->SetXTitle("Pull(#theta_{p2})");
    TH1F *hPull_phi = (TH1F*)hPull->Clone("hPull_phi");
    hPull_phi->SetXTitle("Pull(#phi_{p2})");
    TH1F *hPull_p_pi = (TH1F*)hPull->Clone("hPull_p_pi");
    hPull_p_pi->SetXTitle("Pull(1/P_{#pi^{-}})");
    TH1F *hPull_tht_pi = (TH1F*)hPull->Clone("hPull_tht_pi");
    hPull_tht_pi->SetXTitle("Pull(#theta_{#pi^{-}})");
    TH1F *hPull_phi_pi = (TH1F*)hPull->Clone("hPull_phi_pi");
    hPull_phi_pi->SetXTitle("Pull(#phi_{#pi^{-}})");
    TH1F *hPull_p_K = (TH1F*)hPull->Clone("hPull_p_K");
    hPull_p_K->SetXTitle("Pull(1/P_{K^{+}})");
    TH1F *hPull_tht_K = (TH1F*)hPull->Clone("hPull_tht_K");
    hPull_tht_K->SetXTitle("Pull(#theta_{K^{+}})");
    TH1F *hPull_phi_K = (TH1F*)hPull->Clone("hPull_phi_K");
    hPull_phi_K->SetXTitle("Pull(#phi_{K^{+}})");
    TH1F *hPull_p_p1 = (TH1F*)hPull->Clone("hPull_p_p1");
    hPull_p_p1->SetXTitle("Pull(1/P_{p1})");
    TH1F *hPull_tht_p1 = (TH1F*)hPull->Clone("hPull_tht_p1");
    hPull_tht_p1->SetXTitle("Pull(#theta_{p1})");
    TH1F *hPull_phi_p1 = (TH1F*)hPull->Clone("hPull_phi_p1");
    hPull_phi_p1->SetXTitle("Pull(#phi_{p1})");


    TH1F* hNIterations = new TH1F("hNIterations", "", 20, 0, 20);
    hNIterations->SetXTitle(" Iteration");
    hNIterations->SetYTitle(" events ");
    
    
    
    // -----------------------------------------------------------------------
    
    
   Float_t p1CandTrueP, p1CandTrueTheta, p1CandTruePhi, p1CandTrueR, p1CandTrueZ, p1CandRecoP, p1CandRecoTheta, p1CandRecoPhi, p1CandRecoR, p1CandRecoZ,
            KCandTrueP, KCandTrueTheta, KCandTruePhi, KCandTrueR, KCandTrueZ, KCandRecoP, KCandRecoTheta, KCandRecoPhi,  KCandRecoR, KCandRecoZ,
            p2CandTrueP, p2CandTrueTheta, p2CandTruePhi, p2CandTrueR, p2CandTrueZ, p2CandRecoP, p2CandRecoTheta, p2CandRecoPhi,  p2CandRecoR, p2CandRecoZ,
            piCandTrueP, piCandTrueTheta, piCandTruePhi, piCandTrueR, piCandTrueZ, piCandRecoP, piCandRecoTheta, piCandRecoPhi,  piCandRecoR, piCandRecoZ,
            trueDecVtxX, trueDecVtxY, trueDecVtxZ,
            recoPrimVtxX, recoPrimVtxY, recoPrimVtxZ, recoDecVtxX, recoDecVtxY, recoDecVtxZ;
    
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
    t->SetBranchAddress("trueDecVtxX", &trueDecVtxX);
    t->SetBranchAddress("trueDecVtxY", &trueDecVtxY);
    t->SetBranchAddress("trueDecVtxZ", &trueDecVtxZ);
    t->SetBranchAddress("recoPrimVtxX", &recoPrimVtxX);
    t->SetBranchAddress("recoPrimVtxY", &recoPrimVtxY);
    t->SetBranchAddress("recoPrimVtxZ", &recoPrimVtxZ);
    t->SetBranchAddress("recoDecVtxX", &recoDecVtxX);
    t->SetBranchAddress("recoDecVtxY", &recoDecVtxY);
    t->SetBranchAddress("recoDecVtxZ", &recoDecVtxZ);

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

        TLorentzVector *proton2_gen = new TLorentzVector();
        proton2_gen->SetXYZM(p2CandTrueP * std::sin(p2CandTrueTheta) * std::cos(p2CandTruePhi),
                    p2CandTrueP * std::sin(p2CandTrueTheta) * std::sin(p2CandTruePhi),
                    p2CandTrueP * std::cos(p2CandTrueTheta), 0.938272);

        TLorentzVector *pion_gen = new TLorentzVector();
        pion_gen->SetXYZM(piCandTrueP * std::sin(piCandTrueTheta) * std::cos(piCandTruePhi),
                    piCandTrueP * std::sin(piCandTrueTheta) * std::sin(piCandTruePhi),
                    piCandTrueP * std::cos(piCandTrueTheta), 0.13957);

        TLorentzVector *kaon_gen = new TLorentzVector();
        kaon_gen->SetXYZM(KCandTrueP * std::sin(KCandTrueTheta) * std::cos(KCandTruePhi),
                    KCandTrueP * std::sin(KCandTrueTheta) * std::sin(KCandTruePhi),
                    KCandTrueP * std::cos(KCandTrueTheta), 0.493677);

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

        TLorentzVector lambda_gen = *proton2_gen + *pion_gen;

        TLorentzVector lambda_cand = *proton2 + *pion;
        hLambdaMassPreFit->Fill((lambda_cand).M());

        KFitParticle proton1_fit(*proton1,p1CandRecoR,p1CandRecoZ);
        FillData(proton1_fit, proton1_errors);
        KFitParticle kaon_fit(*kaon,KCandRecoR,KCandRecoZ);
        FillData(kaon_fit, kaon_errors);
        KFitParticle proton2_fit(*proton2,p2CandRecoR,p2CandRecoZ);
        //KFitParticle proton2_fit(proton2,trueDecVtxX,trueDecVtxY,trueDecVtxZ);      // Use this for calculating R and Z with KFitParticle functions instead
        FillData(proton2_fit, proton2_errors);
        KFitParticle pion_fit(*pion,piCandRecoR, piCandRecoZ);
        //KFitParticle pion_fit(pion,trueDecVtxX, trueDecVtxY, trueDecVtxZ);          // Use this for calculating R and Z with KFitParticle functions instead 
        FillData(pion_fit, pion_errors);
        //cout<<"pion momentum filled: "<<pion_fit.P()<<endl;


        hKMomentumPreFitRes->Fill((1/kaon->P()-1/kaon_gen->P())/(1/kaon_gen->P()));
        hKThetaPreFitRes->Fill(kaon->Theta()-kaon_gen->Theta());
        hKPhiPreFitRes->Fill(kaon->Phi()-kaon_gen->Phi());
        hPiMomentumPreFitRes->Fill((1/pion->P()-1/pion_gen->P())/(1/pion_gen->P()));
        hPiThetaPreFitRes->Fill(pion->Theta()-pion_gen->Theta());
        hPiPhiPreFitRes->Fill(pion->Phi()-pion_gen->Phi());

        std::vector<KFitParticle> cands;
        cands.clear();
        cands.push_back(proton1_fit);
        cands.push_back(proton2_fit);
        cands.push_back(pion_fit);    
        cands.push_back(kaon_fit);    
        //cout<<"pion momentum filled: "<<cands[2].P()<<endl;   
        
        // ---------------------------------------------------------------------------------
        // do Missing Particle Fit
        // ---------------------------------------------------------------------------------
        KinFitter fitter(cands);
        fitter.setVerbosity(0);
        //fitter.setLearningRate(0.5);
        //fitter.setConvergenceCriteria(0.01, 1e6, 1e6);
        fitter.setNumberOfIterations(20);
        fitter.add4Constraint(ini);
        if(fitter.fit()){

            KFitParticle fcand1 = fitter.getDaughter(0); // proton1
            KFitParticle fcand2 = fitter.getDaughter(1); // proton2
            KFitParticle fcand3 = fitter.getDaughter(2); // pion
            KFitParticle fcand4 = fitter.getDaughter(3); // kaon

            TLorentzVector lambda_fit = fcand2+fcand3;

            hChi2->Fill(fitter.getChi2());
            hPChi2->Fill(fitter.getProb());
            //TLorentzVector lambda_fit = fcand1 + fcand2;
            hLambdaMassPostFit->Fill(lambda_fit.M());
            // get Pull example (1/P for the fitted proton)
            hPull->Fill(fitter.getPull(5));
            hPull_tht->Fill(fitter.getPull(6));
            hPull_phi->Fill(fitter.getPull(7));
            hPull_p_pi->Fill(fitter.getPull(10));
            hPull_tht_pi->Fill(fitter.getPull(11));
            hPull_phi_pi->Fill(fitter.getPull(12));
            hPull_p_p1->Fill(fitter.getPull(0));
            hPull_tht_p1->Fill(fitter.getPull(1));
            hPull_phi_p1->Fill(fitter.getPull(2));
            hPull_p_K->Fill(fitter.getPull(15));
            hPull_tht_K->Fill(fitter.getPull(16));
            hPull_phi_K->Fill(fitter.getPull(17));
            hKMomPostFit->Fill(fcand4.P());
            hKThetaPostFit->Fill(fcand4.Theta());
            hKPhiPostFit->Fill(fcand4.Phi());
            hKMomPostFitRes->Fill((1/fcand4.P()-1/kaon_gen->P())/(1/kaon_gen->P()));
            hKThetaPostFitRes->Fill(fcand4.Theta()-kaon_gen->Theta());
            hKPhiPostFitRes->Fill(fcand4.Phi()-kaon_gen->Phi());
            hPiMomPostFit->Fill(fcand3.getMomentum());
            hPiThetaPostFit->Fill(fcand3.Theta());
            hPiPhiPostFit->Fill(fcand3.Phi());
            hPiMomPostFitRes->Fill((1/fcand3.P()-1/pion_gen->P())/(1/pion_gen->P()));
            hPiThetaPostFitRes->Fill(fcand3.Theta()-pion_gen->Theta());
            hPiPhiPostFitRes->Fill(fcand3.Phi()-pion_gen->Phi());
            
            if(fitter.getProb()>0.01){
                hChi2_probCut->Fill(fitter.getChi2());
                hLambdaMassPostFit_probCut->Fill(lambda_fit.M());

                // get Pull example (1/P for the fitted proton)
                hPull_probCut->Fill(fitter.getPull(10));
            }

            hNIterations->Fill(fitter.getIteration());

        }
        
    }
    
    //Outut plots
    outfile->cd();
    plotsStyleMacro(hPChi2, "prob_4C");
    plotsStyleMacro(hKMomentumPreFitRes, hKMomPostFitRes, "mom_res_K_4C", "before KinFit", "after KinFit");
    plotsStyleMacro(hKThetaPreFitRes, hKThetaPostFitRes, "tht_res_K_4C", "before KinFit", "after KinFit");
    plotsStyleMacro(hKPhiPreFitRes, hKPhiPostFitRes, "phi_res_K_4C", "before KinFit", "after KinFit");
    plotsStyleMacro(hPiMomentumPreFitRes, hPiMomPostFitRes, "mom_res_Pi_4C", "before KinFit", "after KinFit");
    plotsStyleMacro(hPiThetaPreFitRes, hPiThetaPostFitRes, "tht_res_Pi_4C", "before KinFit", "after KinFit");
    plotsStyleMacro(hPiPhiPreFitRes, hPiPhiPostFitRes, "phi_res_Pi_4C", "before KinFit", "after KinFit");
    plotsStyleMacro_andFit(hPull, "pull_p2_p_4C");
    plotsStyleMacro_andFit(hPull_tht, "pull_p2_tht_4C");
    plotsStyleMacro_andFit(hPull_phi, "pull_p2_phi_4C");
    plotsStyleMacro_andFit(hPull_p_p1, "pull_p1_p_4C");
    plotsStyleMacro_andFit(hPull_tht_p1, "pull_p1_tht_4C");
    plotsStyleMacro_andFit(hPull_phi_p1, "pull_p1_phi_4C");
    plotsStyleMacro_andFit(hPull_p_pi, "pull_pi_p_4C");
    plotsStyleMacro_andFit(hPull_tht_pi, "pull_pi_tht_4C");
    plotsStyleMacro_andFit(hPull_phi_pi, "pull_pi_phi_4C");
    plotsStyleMacro_andFit(hPull_p_K, "pull_K_p_4C");
    plotsStyleMacro_andFit(hPull_tht_K, "pull_K_tht_4C");
    plotsStyleMacro_andFit(hPull_phi_K, "pull_K_phi_4C");
    

    // write histograms to the output file

    hLambdaMassPreFit->Write();
    hKMomentumPreFit->Write();
    hKThetaPreFit->Write();
    hKPhiPreFit->Write();
    hKMomentumPreFitRes->Write();
    hKThetaPreFitRes->Write();
    hKPhiPreFitRes->Write();
    hPiMomentumPreFit->Write();
    hPiThetaPreFit->Write();
    hPiPhiPreFit->Write();
    hPiMomentumPreFitRes->Write();
    hPiThetaPreFitRes->Write();
    hPiPhiPreFitRes->Write();


    hChi2->Write();
    hChi2_probCut->Write();
    hPChi2->Write();
    hLambdaMassPostFit->Write();
    hLambdaMassPostFit_probCut->Write();
    hKMomPostFit->Write();
    hKMomPostFit_probCut->Write();
    hKThetaPostFit->Write();
    hKPhiPostFit->Write();
    hKMomPostFitRes->Write();
    hKThetaPostFitRes->Write();
    hKPhiPostFitRes->Write();
    hPiMomPostFit->Write();
    hPiMomPostFit_probCut->Write();
    hPiThetaPostFit->Write();
    hPiPhiPostFit->Write();
    hPiMomPostFitRes->Write();
    hPiThetaPostFitRes->Write();
    hPiPhiPostFitRes->Write();
    hPull->Write();
    hPull_probCut->Write();
    hPull_tht->Write();
    hPull_phi->Write();
    hPull_p_pi->Write();
    hPull_tht_pi->Write();
    hPull_phi_pi->Write();
    hNIterations->Write();

    outfile->Close();

    return 0;
}