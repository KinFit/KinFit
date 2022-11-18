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

#include "/home/jana/KinFit/include/hkinfitter.h"
#include "/home/jana/KinFit/include/hvertexfinder.h"
#include "/home/jana/KinFit/include/hneutralcandfinder.h"

using namespace std;

void FillData(HRefitCand& outcand, double arr[])
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

Int_t fitLambda3C_toyMC_fromPluto(TString infile, Int_t nEvents)
{
    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile* outfile = new TFile("fitLambda3C_toyMC_fromPluto_myRZ.root", "recreate");

    TH3F* hvertex1xyz = new TH3F("hvertex1xyz", "", 100, -2, 2, 100, -2, 2, 100, -5, 5);
    hvertex1xyz->SetXTitle(" Vx ");
    hvertex1xyz->SetYTitle(" Vy ");
    hvertex1xyz->SetYTitle(" Vz ");
    TH3F* hvertex2xyz = new TH3F("hvertex2xyz", "", 100, -2, 2, 100, -2, 2, 100, -5, 5);
    hvertex2xyz->SetXTitle(" Vx ");
    hvertex2xyz->SetYTitle(" Vy ");
    hvertex2xyz->SetYTitle(" Vz ");
    TH3F* hvertex2xyzres = new TH3F("hvertex2xyzres", "", 100, -2, 2, 100, -2, 2, 100, -5, 5);
    hvertex2xyzres->SetXTitle(" Vx_{gen}-Vx_{reco} ");
    hvertex2xyzres->SetYTitle(" Vy_{gen}-Vy_{reco} ");
    hvertex2xyzres->SetYTitle(" Vz_{gen}-Vz_{reco} ");

    
    TH1F* hLambdaMassPreFit = new TH1F("hLambdaMassPreFit", "", 100, 1.070, 1.250);
    hLambdaMassPreFit->SetXTitle(" M_{p#pi^{-}} [GeV/c^{2}]");
    hLambdaMassPreFit->SetYTitle(" events ");
    TH1F* hLambdaMomentumPreFit = new TH1F("hLambdaMomentumPreFit", "", 100, 2, 4);
    hLambdaMomentumPreFit->SetXTitle(" p_{#Lambda} [GeV/c]");
    hLambdaMomentumPreFit->SetYTitle(" events ");
    TH1F* hLambdaThetaPreFit = new TH1F("hLambdaThetaPreFit", "", 100, 0, 3.2);
    hLambdaThetaPreFit->SetXTitle(" #vartheta_{#Lambda} [rad]");
    hLambdaThetaPreFit->SetYTitle(" events ");
    TH1F* hLambdaPhiPreFit = new TH1F("hLambdaPhiPreFit", "", 100, -3.2, 3.2);
    hLambdaPhiPreFit->SetXTitle(" #varphi_{#Lambda} [rad]");
    hLambdaPhiPreFit->SetYTitle(" events ");
    TH1F* hLambdaMomentumPreFitRes = new TH1F("hLambdaMomentumPreFitRes", "", 100, -0.5, 0.5);
    hLambdaMomentumPreFitRes->SetXTitle(" p_{gen}-p_{reco}(#Lambda) [GeV/c]");
    hLambdaMomentumPreFitRes->SetYTitle(" events ");
    TH1F* hLambdaThetaPreFitRes = new TH1F("hLambdaThetaPreFitRes", "", 100, -0.5, 0.5);
    hLambdaThetaPreFitRes->SetXTitle(" #vartheta_{gen}-#vartheta_{reco}(#Lambda) [rad]");
    hLambdaThetaPreFitRes->SetYTitle(" events ");
    TH1F* hLambdaPhiPreFitRes = new TH1F("hLambdaPhiPreFitRes", "", 100, -0.5, 0.5);
    hLambdaPhiPreFitRes->SetXTitle(" #varphi_{gen}-#varphi_{reco}(#Lambda) [rad]");
    hLambdaPhiPreFitRes->SetYTitle(" events ");

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

    TH1F* hLambdaMassPostFit = new TH1F("hLambdaMassPostFit", "", 100, 1.070, 1.250);
    hLambdaMassPostFit->SetXTitle(" M_{p#pi^{-}} [GeV/c^{2}]");
    hLambdaMassPostFit->SetYTitle(" events ");
    hLambdaMassPostFit -> SetLineColor(kRed);
    TH1F *hLambdaMassPostFit_probCut = (TH1F*)hLambdaMassPostFit->Clone("hLambdaMassPostFit_probCut");
    hLambdaMassPostFit_probCut -> SetLineColor(kGreen);
    TH1F *hLambdaMassPostFit_converged = (TH1F*)hLambdaMassPostFit->Clone("hLambdaMassPostFit_converged");
    hLambdaMassPostFit_converged -> SetLineColor(kBlue);

    TH1F* hLambdaMomPostFit = new TH1F("hLambdaMomPostFit", "", 100, 2, 4);
    hLambdaMomPostFit->SetXTitle(" P_{#Lambda} [GeV/c]");
    hLambdaMomPostFit->SetYTitle(" events ");
    hLambdaMomPostFit -> SetLineColor(kRed);
    TH1F *hLambdaMomPostFit_probCut = (TH1F*)hLambdaMomPostFit->Clone("hLambdaMomPostFit_probCut");
    hLambdaMomPostFit_probCut -> SetLineColor(kGreen);
    TH1F* hLambdaThetaPostFit = new TH1F("hLambdaThetaPostFit", "", 100, 0, 3.2);
    hLambdaThetaPostFit->SetXTitle(" #vartheta_{#Lambda} [rad]");
    hLambdaThetaPostFit->SetYTitle(" events ");
    TH1F* hLambdaPhiPostFit = new TH1F("hLambdaPhiPostFit", "", 100, -3.2, 3.2);
    hLambdaPhiPostFit->SetXTitle(" #varphi_{#Lambda} [rad]");
    hLambdaPhiPostFit->SetYTitle(" events ");
    TH1F* hLambdaMomPostFitRes = new TH1F("hLambdaMomPostFitRes", "", 100, -0.5, 0.5);
    hLambdaMomPostFitRes->SetXTitle(" p_{gen}-p_{reco}(#Lambda) [GeV/c]");
    hLambdaMomPostFitRes->SetYTitle(" events ");
    TH1F* hLambdaThetaPostFitRes = new TH1F("hLambdaThetaPostFitRes", "", 100, -0.05, 0.05);
    hLambdaThetaPostFitRes->SetXTitle(" #vartheta_{gen}-#vartheta_{reco}(#Lambda) [rad]");
    hLambdaThetaPostFitRes->SetYTitle(" events ");
    TH1F* hLambdaPhiPostFitRes = new TH1F("hLambdaPhiPostFitRes", "", 100, -0.05, 0.05);
    hLambdaPhiPostFitRes->SetXTitle(" #varphi_{gen}-#varphi_{reco}(#Lambda) [rad]");
    hLambdaPhiPostFitRes->SetYTitle(" events ");

    TH1F* hPull = new TH1F("hPull", "", 100, -5, 5);
    hPull->SetXTitle("Pull(1/P_{p})");
    hPull->SetYTitle(" counts ");
    TH1F *hPull_probCut = (TH1F*)hPull->Clone("hPull_probCut");
    hPull_probCut -> SetLineColor(kGreen);
    TH1F *hPull_converged = (TH1F*)hPull->Clone("hPull_converged");
    hPull_converged -> SetLineColor(kBlue);
    TH1F *hPull_tht = (TH1F*)hPull->Clone("hPull_tht");
    TH1F *hPull_phi = (TH1F*)hPull->Clone("hPull_phi");
    TH1F *hPull_p_pi = (TH1F*)hPull->Clone("hPull_p_pi");
    TH1F *hPull_tht_pi = (TH1F*)hPull->Clone("hPull_tht_pi");
    TH1F *hPull_phi_pi = (TH1F*)hPull->Clone("hPull_phi_pi");

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
    TH1F* hNIterations = new TH1F("hNIterations", "", 10, 0, 10);
    hNIterations->SetXTitle(" Iteration");
    hNIterations->SetYTitle(" events ");
    
    
    
    // -----------------------------------------------------------------------
    
    
   Float_t p1CandTrueP, p1CandTrueTheta, p1CandTruePhi, p1CandTrueR, p1CandTrueZ, p1CandRecoP, p1CandRecoTheta, p1CandRecoPhi, p1CandRecoR, p1CandRecoZ,
            KCandTrueP, KCandTrueTheta, KCandTruePhi, KCandTrueR, KCandTrueZ, KCandRecoP, KCandRecoTheta, KCandRecoPhi,  KCandRecoR, KCandRecoZ,
            p2CandTrueP, p2CandTrueTheta, p2CandTruePhi, p2CandTrueR, p2CandTrueZ, p2CandRecoP, p2CandRecoTheta, p2CandRecoPhi,  p2CandRecoR, p2CandRecoZ,
            piCandTrueP, piCandTrueTheta, piCandTruePhi, piCandTrueR, piCandTrueZ, piCandRecoP, piCandRecoTheta, piCandRecoPhi,  piCandRecoR, piCandRecoZ,
            trueDecVtxX, trueDecVtxY, trueDecVtxZ;
    
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
        //h01->Fill(lambda_gen.M());
        //h012->Fill(lambda_gen.P());

        HRefitCand proton1_fit(proton1,p1CandRecoR,p1CandRecoZ);
        FillData(proton1_fit, proton1_errors);
        HRefitCand kaon_fit(kaon,KCandRecoR,KCandRecoZ);
        FillData(kaon_fit, kaon_errors);
        HRefitCand proton2_fit(proton2,p2CandRecoR,p2CandRecoZ);
        //HRefitCand proton2_fit(proton2,trueDecVtxX,trueDecVtxY,trueDecVtxZ);
        FillData(proton2_fit, proton2_errors);
        HRefitCand pion_fit(pion,piCandRecoR, piCandRecoZ);
        //HRefitCand pion_fit(pion,trueDecVtxX, trueDecVtxY, trueDecVtxZ);
        FillData(pion_fit, pion_errors);

        // ---------------------------------------------------------------------------------
        // find the primary and decay vertex
        // ---------------------------------------------------------------------------------
        std::vector<HRefitCand> cands1, cands2;
        cands1.clear();
        cands1.push_back(proton1_fit);
        cands1.push_back(kaon_fit);
        cands2.clear();
        cands2.push_back(proton2_fit);
        cands2.push_back(pion_fit);

        HVertexFinder vtx1finder(cands1);
        HVertexFinder vtx2finder(cands2);

        TVector3 vtx1 = vtx1finder.getVertex();
        TVector3 vtx2 = vtx2finder.getVertex();

        hvertex1xyz->Fill(vtx1.X(), vtx1.Y(), vtx1.Z());
        hvertex2xyz->Fill(vtx2.X(), vtx2.Y(), vtx2.Z());
        hvertex2xyzres->Fill(trueDecVtxX-vtx2.X(), trueDecVtxY-vtx2.Y(), trueDecVtxZ-vtx2.Z());

        // ---------------------------------------------------------------------------------
        // find the neutral candidate
        // ---------------------------------------------------------------------------------

        HNeutralCandFinder lambdafinder(cands2, 1.115683, vtx2, vtx1, 0.26, 0.26, 0.69, 0.26, 0.26, 0.7);
        lambdafinder.setVerbosity(10);
        lambdafinder.setNeutralMotherCand();

        HRefitCand lambda_cand = lambdafinder.getNeutralMotherCandidate();
        cout<<"lambda momentum"<<lambda_cand.P()<<endl;
        cout<<"lambda Z"<<lambda_cand.getZ()<<endl;
        
        //hLambdaMassPreFit->Fill((proton2_fit + pion_fit).M());
        hLambdaMassPreFit->Fill((lambda_cand).M());
        hLambdaMomentumPreFit->Fill(lambda_cand.P());
        hLambdaThetaPreFit->Fill(lambda_cand.Theta());
        hLambdaPhiPreFit->Fill(lambda_cand.Phi());

        hLambdaMomentumPreFitRes->Fill(lambda_cand.P()-lambda_gen.P());
        hLambdaThetaPreFitRes->Fill(lambda_cand.Theta()-lambda_gen.Theta());
        hLambdaPhiPreFitRes->Fill(lambda_cand.Phi()-lambda_gen.Phi());
/*
     TMatrixD covariance = lambda_cand.getCovariance();
     std::cout << "Lambda covariances: " << covariance(0, 0) << " "<< covariance(1, 1) << " " << covariance(2, 2) << " " << covariance(3, 3) << " " << covariance(4, 4) <<  std::endl;
     covariance.Zero();
     covariance(0,0) = 1e7;
     covariance(1,1) = 0.075*0.075;
     covariance(2,2) = 0.1*0.1;
     covariance(3,3) = 0.1*0.1;
     covariance(4,4) = 0.1*0.1;
     lambda_cand.setCovariance(covariance);
     */
        
        // ---------------------------------------------------------------------------------
        // do 3C fit in decay vertex
        // ---------------------------------------------------------------------------------
        HKinFitter fitter(cands2);
        fitter.setVerbosity(0);
        //fitter.setLearningRate(0.5);
        fitter.setConvergenceCriterion(0.01);
        fitter.setNumberOfIterations(10);
        fitter.add3Constraint(lambda_cand);
        //fitter.addVertexConstraint();
        //fitter.add4Constraint(ini);
        fitter.fit();
        if(1){

            HRefitCand fcand1 = fitter.getDaughter(0); // proton
            HRefitCand fcand2 = fitter.getDaughter(1); // pion
            HRefitCand lambda_fit;
            TMatrixD test = lambda_fit.getCovariance();
            test.ResizeTo(5,5);
            lambda_fit.setCovariance(test);
            lambda_fit = fitter.getMother();

            hChi2->Fill(fitter.getChi2());
            hPChi2->Fill(fitter.getProb());
            //TLorentzVector lambda_fit = fcand1 + fcand2;
            hLambdaMassPostFit->Fill(lambda_fit.M());
            cout<<"LAmbda Mass: "<<lambda_fit.M()<<endl;
            // get Pull example (1/P for the fitted proton)
            hPull->Fill(fitter.getPull(0));
            hPull_tht->Fill(fitter.getPull(1));
            hPull_phi->Fill(fitter.getPull(2));
            hPull_p_pi->Fill(fitter.getPull(5));
            hPull_tht_pi->Fill(fitter.getPull(6));
            hPull_phi_pi->Fill(fitter.getPull(7));
            hLambdaMomPostFit->Fill(lambda_fit.P());
            hLambdaThetaPostFit->Fill(lambda_fit.Theta());
            hLambdaPhiPostFit->Fill(lambda_fit.Phi());
            hLambdaMomPostFitRes->Fill(lambda_fit.P()-lambda_gen.P());
            hLambdaThetaPostFitRes->Fill(lambda_fit.Theta()-lambda_gen.Theta());
            hLambdaPhiPostFitRes->Fill(lambda_fit.Phi()-lambda_gen.Phi());
            
            if(fitter.getProb()>0.01){
                hChi2_probCut->Fill(fitter.getChi2());
                hLambdaMassPostFit_probCut->Fill(lambda_fit.M());
               // h072->Fill(all_fit.P());

                // get Pull example (1/P for the fitted proton)
                hPull_probCut->Fill(fitter.getPull(0));
                hLambdaMomPostFit_probCut->Fill(lambda_fit.P());
            }

            hNIterations->Fill(fitter.getIteration());

        }
        
    }

    // write histograms to the output file
    outfile->cd();

    hvertex1xyz->Write();
    hvertex2xyz->Write();
    hvertex2xyzres->Write();

    hLambdaMassPreFit->Write();
    hLambdaMomentumPreFit->Write();
    hLambdaThetaPreFit->Write();
    hLambdaPhiPreFit->Write();
    hLambdaMomentumPreFitRes->Write();
    hLambdaThetaPreFitRes->Write();
    hLambdaPhiPreFitRes->Write();


    hChi2->Write();
    hChi2_probCut->Write();
    hPChi2->Write();
    hLambdaMassPostFit->Write();
    hLambdaMassPostFit_probCut->Write();
    hLambdaMomPostFit->Write();
    hLambdaMomPostFit_probCut->Write();
    hLambdaThetaPostFit->Write();
    hLambdaPhiPostFit->Write();
    hLambdaMomPostFitRes->Write();
    hLambdaThetaPostFitRes->Write();
    hLambdaPhiPostFitRes->Write();
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