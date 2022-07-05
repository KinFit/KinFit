#include "hades.h"
#include "hcategorymanager.h"
#include "henergylosscorrpar.h"
#include "hhistmap.h"
#include "hloop.h"
#include "hparticleanglecor.h"
#include "hparticlepairmaker.h"
#include "hparticletool.h"
#include "hparticletracksorter.h"
#include "hphysicsconstants.h"
#include "htool.h"

#include "hcategory.h"
#include "hlinearcategory.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticlegeantpair.h"
#include "hparticlepair.h"
#include "hrichhit.h"
#include "hrichhitsim.h"

#include "hgeantkine.h"
#include "hparticledef.h"
#include "hstartdef.h"
#include "richdef.h"

#include "hparticlecutrange.h"
#include "hparticlegeant.h"
#include "hparticlegeantdecay.h"
#include "hparticlegeantevent.h"

#include "TTree.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "hdecaybuilder.h"
#include "hkinfitter.h"
#include "hvertexfinder.h"
#include "hneutralcandfinder.h"

using namespace std;
using namespace Particle;

void FillData(HParticleCand *cand, HRefitCand &outcand, double arr[],
              double mass)
{
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

Bool_t selectHadrons(HParticleCand *pcand)
{
    // build in selection function for hadron candidates.
    // Requires besides an RK + META and fitted
    // inner+outer segment.

    Bool_t test = kFALSE;
    if (pcand->isFlagAND(4, Particle::kIsAcceptedHitInnerMDC,
                         Particle::kIsAcceptedHitOuterMDC,
                         Particle::kIsAcceptedHitMETA,
                         Particle::kIsAcceptedRK) &&
        pcand->getInnerSegmentChi2() > 0 && pcand->getChi2() < 10000 // RK
    )
        test = kTRUE;

    if (!test)
        return kFALSE;

    if (test)
        test = pcand->getMetaMatchQuality() < 3 ? kTRUE : kFALSE;

    return test;
}

Int_t decayBuilderAnalysis(TString infileList = "/lustre/hades/user/jregina/Pluto/DST/OutputFolder/pplpk_million_Hgeant_*1_dst_apr12.root", Int_t nEvents = 1500000)
{

    // Original input: pp_pKlambda_100000evts1_dst_apr12.root
    // /lustre/hades/user/jregina/Pluto/DST/OutputFolder/pplpk_million_Hgeant_[1-9]1_dst_apr12.root,/lustre/hades/user/jregina/Pluto/DST/OutputFolder/pplpk_million_Hgeant_1[0-9]1_dst_apr12.root
    TStopwatch timer;
    timer.Start();

    //Momentum dependent uncertainty estimation input
    TFile *momErr = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p3500p_momDepMomErrors_allParticles_errFunc.root", "read");
    TF1 *momErrP = (TF1 *)momErr->Get("f_pP");
    TF1 *momErrPi = (TF1 *)momErr->Get("f_pPi");
    TF1 *momErrK = (TF1 *)momErr->Get("f_pK");

    TFile *thtErr = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p3500p_momDepThtErrors_allParticles_errFunc.root", "read");
    TFile *thtErr_PK = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p3500p_momDepThtErrors_PK_errFunc.root", "read");
    TF1 *thtErrP = (TF1 *)thtErr->Get("f_thtP");
    TF1 *thtErrPi = (TF1 *)thtErr->Get("f_thtPi");
    TF1 *thtErrK = (TF1 *)thtErr_PK->Get("f_thtK");

    TFile *phiErr = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepPhiErrors_PPi_errFunc.root", "read"); //if 4500 in name this is just a name issue
    TFile *phiErr_PK = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepPhiErrors_PK_errFunc.root", "read");
    TF1 *phiErrP = (TF1 *)phiErr->Get("f_phiP");
    TF1 *phiErrPi = (TF1 *)phiErr->Get("f_phiPi");
    TF1 *phiErrK = (TF1 *)phiErr_PK->Get("f_phiK");

    TFile *RZErr_PPi = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepRZErrors_PPi_errFunc.root", "read");
    TFile *RZErr_PK = new TFile("/lustre/hades/user/jrieger/pp_pKLambda/sim/ana/pp_pKLambda_p4500p_momDepRZErrors_PK_errFunc.root", "read");
    TF1 *RErrP = (TF1 *)RZErr_PPi->Get("f_RP");
    TF1 *RErrPi = (TF1 *)RZErr_PPi->Get("f_RPi");
    TF1 *RErrK = (TF1 *)RZErr_PK->Get("f_RK");
    TF1 *ZErrP = (TF1 *)RZErr_PPi->Get("f_ZP");
    TF1 *ZErrPi = (TF1 *)RZErr_PPi->Get("f_ZPi");
    TF1 *ZErrK = (TF1 *)RZErr_PK->Get("f_ZK");

    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile *outfile = new TFile("vertexfit_Realistic.root", "recreate");

    TH1F *h01 = new TH1F("hLambdaMassPreFit", "", 100, 1070, 1170);
    h01->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h01->SetYTitle(" events ");

    TH1F *h02 = new TH1F("hChi2", "", 1000, 0, 10);
    h02->SetXTitle("#chi^{2}");
    h02->SetYTitle(" Counts ");
    TH1F *h021 = (TH1F *)h02->Clone("hChi2_3c");
    TH1F *h023 = (TH1F *)h02->Clone("hChi2_3cConv");

    TH1F *h03Zoomed = new TH1F("hPChi2Zoomed", "", 10000, 0, pow(10, -200));
    h03Zoomed->SetXTitle("P(#chi^{2})");
    h03Zoomed->SetYTitle(" Counts ");

    TH1F *h03 = new TH1F("hPChi2", "", 1000, 0, 1);
    h03->SetXTitle("P(#chi^{2})");
    h03->SetYTitle(" Counts ");
    TH1F *h031 = (TH1F *)h03->Clone("hPChi2_3c");
    TH1F *h033 = (TH1F *)h03->Clone("hPChi2_3cConv");

    TH1F *h02SecondaryVtx = new TH1F("hChi2SecondaryVtx", "", 1000, 0, 10);
    h02SecondaryVtx->SetXTitle("#chi^{2}");
    h02SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h03SecondaryVtx = new TH1F("hPChi2SecondaryVtx", "", 1000, 0, 1);
    h03SecondaryVtx->SetXTitle("P(#chi^{2})");
    h03SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h03SecondaryVtxZoomed = new TH1F("hPChi2SecondaryVtxZoomed", "", 10000, 0, 0.0001);
    h03SecondaryVtxZoomed->SetXTitle("P(#chi^{2})");
    h03SecondaryVtxZoomed->SetYTitle(" Counts ");

    TH1F *h04 = new TH1F("hLambdaMassPostFit", "", 100, 1070, 1170);
    h04->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    h04->SetYTitle(" Events ");
    TH1F *hLambdaMassBeforeFit = (TH1F *)h04->Clone("hLambdaMassBeforeFit");
    TH1F *h040 = (TH1F *)h04->Clone("hLambdaMassPostFitCut");
    TH1F *h041 = (TH1F *)h04->Clone("hLambdaMassPostFit_3c");
    TH1F *h044 = (TH1F *)h04->Clone("hLambdaMassPostFit_3cConvCut");
    TH1F *h045 = (TH1F *)h04->Clone("hLambdaMassPreFit");
    //TH1F *h17 = new TH1F("hLambdaMassPostFitCut", "", 100, 1070, 1170);
    //h17->SetXTitle(" M_{p#pi^{-}} [MeV/c^{2}]");
    //h17->SetYTitle(" Events ");

    TH1F *h05 = new TH1F("hPullPInv", "", 100, -5, 5);
    h05->SetXTitle("Pull(1/P_{p})");
    h05->SetYTitle(" Counts ");
    TH1F *h050 = (TH1F *)h05->Clone("hPullPInvCut");
    TH1F *h051 = (TH1F *)h05->Clone("hPullPInv_3c");
    TH1F *h054 = (TH1F *)h05->Clone("hPullPInv_3cConvCut");

    TH1F *h06 = new TH1F("hPullTheta", "", 100, -5, 5);
    h06->SetXTitle("Pull(#theta)");
    h06->SetYTitle(" Counts ");

    TH1F *h07 = new TH1F("hPullPhi", "", 100, -5, 5);
    h07->SetXTitle("Pull(#phi)");
    h07->SetYTitle(" Counts ");

    TH1F *h08 = new TH1F("hPullR", "", 100, -5, 5);
    h08->SetXTitle("Pull(R)");
    h08->SetYTitle(" Counts ");

    TH1F *h09 = new TH1F("hPullZ", "", 100, -5, 5);
    h09->SetXTitle("Pull(Z)");
    h09->SetYTitle(" Counts ");

    TH1F *h06SecondaryVtx = new TH1F("hPullThetaSecondaryVtx", "", 100, -5, 5);
    h06SecondaryVtx->SetXTitle("Pull(#theta)");
    h06SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h07SecondaryVtx = new TH1F("hPullPhiSecondaryVtx", "", 100, -5, 5);
    h07SecondaryVtx->SetXTitle("Pull(#phi)");
    h07SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h08SecondaryVtx = new TH1F("hPullRSecondaryVtx", "", 100, -5, 5);
    h08SecondaryVtx->SetXTitle("Pull(R)");
    h08SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h09SecondaryVtx = new TH1F("hPullZSecondaryVtx", "", 100, -5, 5);
    h09SecondaryVtx->SetXTitle("Pull(Z)");
    h09SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h10 = new TH1F("hSuccessfulConvergence", "", 2, 0, 2);
    h10->SetXTitle("Fit has converged");
    h10->SetYTitle(" Counts ");

    TH1F *h11 = new TH1F("hIterations", "", 100, 0, 20);
    h11->SetXTitle("Number of Iterations");
    h11->SetYTitle(" Counts ");
    TH1F *h111 = (TH1F *)h11->Clone("hNIterations_3c");
    TH1F *h114 = (TH1F *)h11->Clone("hNIterations_3cCut");

    TH1F *h11SecVtx = new TH1F("hIterationsSecVtx", "", 100, 0, 20);
    h11SecVtx->SetXTitle("Number of Iterations");
    h11SecVtx->SetYTitle(" Counts ");

    // ------------------------- Histos after probability cut    --------------------

    TH1F *h12 = new TH1F("hPullPInvCut", "", 100, -5, 5);
    h12->SetXTitle("Pull(1/P_{p})");
    h12->SetYTitle(" Counts ");

    TH1F *h13 = new TH1F("hPullThetaCut", "", 100, -5, 5);
    h13->SetXTitle("Pull(#theta)");
    h13->SetYTitle(" Counts ");

    TH1F *h14 = new TH1F("hPullPhiCut", "", 100, -5, 5);
    h14->SetXTitle("Pull(#phi)");
    h14->SetYTitle(" Counts ");

    TH1F *h15 = new TH1F("hPullRCut", "", 100, -5, 5);
    h15->SetXTitle("Pull(R)");
    h15->SetYTitle(" Counts ");

    TH1F *h16 = new TH1F("hPullZCut", "", 100, -5, 5);
    h16->SetXTitle("Pull(Z)");
    h16->SetYTitle(" Counts ");

    TH1F *h12SecondaryVtx = new TH1F("hPullPInvCutSecondaryVtx", "", 100, -5, 5);
    h12SecondaryVtx->SetXTitle("Pull(1/P_{p})");
    h12SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h13SecondaryVtx = new TH1F("hPullThetaCutSecondaryVtx", "", 100, -5, 5);
    h13SecondaryVtx->SetXTitle("Pull(#theta)");
    h13SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h14SecondaryVtx = new TH1F("hPullPhiCutSecondaryVtx", "", 100, -5, 5);
    h14SecondaryVtx->SetXTitle("Pull(#phi)");
    h14SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h15SecondaryVtx = new TH1F("hPullRCutSecondaryVtx", "", 100, -5, 5);
    h15SecondaryVtx->SetXTitle("Pull(R)");
    h15SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h16SecondaryVtx = new TH1F("hPullZCutSecondaryVtx", "", 100, -5, 5);
    h16SecondaryVtx->SetXTitle("Pull(Z)");
    h16SecondaryVtx->SetYTitle(" Counts ");

    // -------------------- Pulls Pions -----------------------

    TH1F *h17 = new TH1F("hPullThetaPion", "", 100, -5, 5);
    h17->SetXTitle("Pull(#theta)");
    h17->SetYTitle(" Counts ");

    TH1F *h18 = new TH1F("hPullPhiPion", "", 100, -5, 5);
    h18->SetXTitle("Pull(#phi)");
    h18->SetYTitle(" Counts ");

    TH1F *h19 = new TH1F("hPullRPion", "", 100, -5, 5);
    h19->SetXTitle("Pull(R)");
    h19->SetYTitle(" Counts ");

    TH1F *h20 = new TH1F("hPullZPion", "", 100, -5, 5);
    h20->SetXTitle("Pull(Z)");
    h20->SetYTitle(" Counts ");

    // ------------------------- Histos after probability cut    --------------------

    TH1F *h21 = new TH1F("hPullThetaCutPion", "", 100, -5, 5);
    h21->SetXTitle("Pull(#theta)");
    h21->SetYTitle(" Counts ");

    TH1F *h22 = new TH1F("hPullPhiCutPion", "", 100, -5, 5);
    h22->SetXTitle("Pull(#phi)");
    h22->SetYTitle(" Counts ");

    TH1F *h23 = new TH1F("hPullRCutPion", "", 100, -5, 5);
    h23->SetXTitle("Pull(R)");
    h23->SetYTitle(" Counts ");

    TH1F *h24 = new TH1F("hPullZCutPion", "", 100, -5, 5);
    h24->SetXTitle("Pull(Z)");
    h24->SetYTitle(" Counts ");

    TH1F *h21SecondaryVtx = new TH1F("hPullThetaCutPionSecondaryVtx", "", 100, -5, 5);
    h21SecondaryVtx->SetXTitle("Pull(#theta)");
    h21SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h22SecondaryVtx = new TH1F("hPullPhiCutPionSecondaryVtx", "", 100, -5, 5);
    h22SecondaryVtx->SetXTitle("Pull(#phi)");
    h22SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h23SecondaryVtx = new TH1F("hPullRCutPionSecondaryVtx", "", 100, -5, 5);
    h23SecondaryVtx->SetXTitle("Pull(R)");
    h23SecondaryVtx->SetYTitle(" Counts ");

    TH1F *h24SecondaryVtx = new TH1F("hPullZCutPionSecondaryVtx", "", 100, -5, 5);
    h24SecondaryVtx->SetXTitle("Pull(Z)");
    h24SecondaryVtx->SetYTitle(" Counts ");
    // -----------------------------------------------------------------------
    // ------------------------- Histos after convergence    --------------------

    TH1F *h25 = new TH1F("hPullPInvConverged", "", 100, -5, 5);
    h25->SetXTitle("Pull(1/P_{p})");
    h25->SetYTitle(" Counts ");

    TH1F *h26 = new TH1F("hPullThetaConverged", "", 100, -5, 5);
    h26->SetXTitle("Pull(#theta)");
    h26->SetYTitle(" Counts ");

    TH1F *h27 = new TH1F("hPullPhiConverged", "", 100, -5, 5);
    h27->SetXTitle("Pull(#phi)");
    h27->SetYTitle(" Counts ");

    TH1F *h28 = new TH1F("hPullRConverged", "", 100, -5, 5);
    h28->SetXTitle("Pull(R)");
    h28->SetYTitle(" Counts ");

    TH1F *h29 = new TH1F("hPullZConverged", "", 100, -5, 5);
    h29->SetXTitle("Pull(Z)");
    h29->SetYTitle(" Counts ");

    TH1F *h30 = new TH1F("hPullThetaPionConverged", "", 100, -5, 5);
    h30->SetXTitle("Pull(#theta)");
    h30->SetYTitle(" Counts ");

    TH1F *h31 = new TH1F("hPullPhiPioConvergedn", "", 100, -5, 5);
    h31->SetXTitle("Pull(#phi)");
    h31->SetYTitle(" Counts ");

    TH1F *h32 = new TH1F("hPullRPionConverged", "", 100, -5, 5);
    h32->SetXTitle("Pull(R)");
    h32->SetYTitle(" Counts ");

    TH1F *h33 = new TH1F("hPullZPionConverged", "", 100, -5, 5);
    h33->SetXTitle("Pull(Z)");
    h33->SetYTitle(" Counts ");

    // -----------------------------------------------------------------------
    // ------- Vertex histograms pre fit -----------------------

    TH1F *hVertexXPreFit = new TH1F("hVertexXPreFit", "", 1000, -100, 100);
    hVertexXPreFit->SetXTitle("Vertex, X / mm");
    hVertexXPreFit->SetYTitle(" Counts ");

    TH1F *hVertexYPreFit = new TH1F("hVertexYPreFit", "", 1000, -100, 100);
    hVertexYPreFit->SetXTitle("Vertex, Y / mm");
    hVertexYPreFit->SetYTitle(" Counts ");

    TH1F *hVertexZPreFit = new TH1F("hVertexZPreFit", "", 1000, -60, 1000);
    hVertexZPreFit->SetXTitle("Vertex, Z / mm");
    hVertexZPreFit->SetYTitle(" Counts ");

    TH2F *hVertexPreFit = new TH2F("hVertexPreFit", "", 1000, -100, 100, 1000, -100, 100);
    hVertexPreFit->SetXTitle("Vertex, Z / mm");
    hVertexPreFit->SetYTitle("Vertex, R / mm");

    TH1F *hPrimaryVertexXPreFit = new TH1F("hPrimaryVertexXPreFit", "", 1000, -100, 100);
    hPrimaryVertexXPreFit->SetXTitle("Vertex, X / mm");
    hPrimaryVertexXPreFit->SetYTitle(" Counts ");

    TH1F *hPrimaryVertexYPreFit = new TH1F("hPrimaryVertexYPreFit", "", 1000, -100, 100);
    hPrimaryVertexYPreFit->SetXTitle("Vertex, Y / mm");
    hPrimaryVertexYPreFit->SetYTitle(" Counts ");

    TH1F *hPrimaryVertexZPreFit = new TH1F("hPrimaryVertexZPreFit", "", 1000, -60, 1000);
    hPrimaryVertexZPreFit->SetXTitle("Vertex, Z / mm");
    hPrimaryVertexZPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceBetweenProtonAndPionPreFit = new TH1F("hDistanceBetweenProtonAndPionPreFit", "", 500, 0, 50);
    hDistanceBetweenProtonAndPionPreFit->SetXTitle("Distance between particles / mm");
    hDistanceBetweenProtonAndPionPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceBetweenPrimaryProtonAndKaonPreFit = new TH1F("hDistanceBetweenPrimaryProtonAndKaonPreFit", "", 500, 0, 50);
    hDistanceBetweenPrimaryProtonAndKaonPreFit->SetXTitle("Distance between particles / mm");
    hDistanceBetweenPrimaryProtonAndKaonPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToVertexProtonPreFit = new TH1F("hDistanceToVertexProtonPreFit", "", 1000, 0, 100);
    hDistanceToVertexProtonPreFit->SetXTitle("Distance to vertex / mm");
    hDistanceToVertexProtonPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToVertexPionPreFit = new TH1F("hDistanceToVertexPionPreFit", "", 1000, 0, 100);
    hDistanceToVertexPionPreFit->SetXTitle("Distance to vertex / mm");
    hDistanceToVertexPionPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToVertexPrimaryProtonPreFit = new TH1F("hDistanceToVertexPrimaryProtonPreFit", "", 1000, 0, 100);
    hDistanceToVertexPrimaryProtonPreFit->SetXTitle("Distance to vertex / mm");
    hDistanceToVertexPrimaryProtonPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToVertexKaonPreFit = new TH1F("hDistanceToVertexKaonPreFit", "", 1000, 0, 100);
    hDistanceToVertexKaonPreFit->SetXTitle("Distance to vertex / mm");
    hDistanceToVertexKaonPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToOriginProtonPreFit = new TH1F("hDistanceToOriginProtonPreFit", "", 1000, 0, 100);
    hDistanceToOriginProtonPreFit->SetXTitle("Distance to origin / mm");
    hDistanceToOriginProtonPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToOriginPionPreFit = new TH1F("hDistanceToOriginPionPreFit", "", 1000, 0, 100);
    hDistanceToOriginPionPreFit->SetXTitle("Distance to Origin / mm");
    hDistanceToOriginPionPreFit->SetYTitle(" Counts ");

    TH1F *hDistanceToOriginSumPreFit = new TH1F("hDistanceToOriginSumPreFit", "", 1000, 0, 100);
    hDistanceToOriginSumPreFit->SetXTitle("Distance to Origin Sum / mm");
    hDistanceToOriginSumPreFit->SetYTitle(" Counts ");

    // --------- Vertex histograms post fit --------------------

    TH1F *hVertexXPostFit = new TH1F("hVertexXPostFit", "", 1000, -100, 100);
    hVertexXPostFit->SetXTitle("Vertex, X / mm");
    hVertexXPostFit->SetYTitle(" Counts ");

    TH1F *hVertexYPostFit = new TH1F("hVertexYPostFit", "", 1000, -100, 100);
    hVertexYPostFit->SetXTitle("Vertex, Y / mm");
    hVertexYPostFit->SetYTitle(" Counts ");

    TH1F *hVertexZPostFit = new TH1F("hVertexZPostFit", "", 1000, -100, 100);
    hVertexZPostFit->SetXTitle("Vertex, Z / mm");
    hVertexZPostFit->SetYTitle(" Counts ");

    TH2F *hVertexPostFit = new TH2F("hVertexPostFit", "", 1000, -100, 100, 1000, -100, 100);
    hVertexPostFit->SetXTitle("Vertex, Z / mm");
    hVertexPostFit->SetYTitle("Vertex, R /mm");

    TH1F *hPrimVertexXPostFit = new TH1F("hPrimVertexXPostFit", "", 1000, -100, 100);
    hPrimVertexXPostFit->SetXTitle("Vertex, X / mm");
    hPrimVertexXPostFit->SetYTitle(" Counts ");

    TH1F *hPrimVertexYPostFit = new TH1F("hPrimVertexYPostFit", "", 1000, -100, 100);
    hPrimVertexYPostFit->SetXTitle("Vertex, Y / mm");
    hPrimVertexYPostFit->SetYTitle(" Counts ");

    TH1F *hPrimVertexZPostFit = new TH1F("hPrimVertexZPostFit", "", 1000, -100, 100);
    hPrimVertexZPostFit->SetXTitle("Vertex, Z / mm");
    hPrimVertexZPostFit->SetYTitle(" Counts ");

    TH1F *hDistanceBetweenProtonAndPionPostFit = new TH1F("hDistanceBetweenProtonAndPionPostFit", "", 1000, 0, 100);
    hDistanceBetweenProtonAndPionPostFit->SetXTitle("Distance between particles / mm");
    hDistanceBetweenProtonAndPionPostFit->SetYTitle(" Counts ");

    TH1F *hDistanceToVertexProtonPostFit = new TH1F("hDistanceToVertexProtonPostFit", "", 1000, 0, 100);
    hDistanceToVertexProtonPostFit->SetXTitle("Distance to vertex / mm");
    hDistanceToVertexProtonPostFit->SetYTitle(" Counts ");

    TH1F *hDistanceToVertexPionPostFit = new TH1F("hDistanceToVertexPionPostFit", "", 1000, 0, 100);
    hDistanceToVertexPionPostFit->SetXTitle("Distance to vertex / mm");
    hDistanceToVertexPionPostFit->SetYTitle(" Counts ");

    //---------------------- MC (Geant) information ------------------------

    TH1F *hGeantVertexXAll = new TH1F("hGeantVertexXAll", "", 1000, -100, 100);
    hGeantVertexXAll->SetXTitle("Geant Vertex, X / mm");
    hGeantVertexXAll->SetYTitle(" Counts ");

    TH1F *hGeantVertexYAll = new TH1F("hGeantVertexYAll", "", 1000, -100, 100);
    hGeantVertexYAll->SetXTitle("Geant Vertex, Y / mm");
    hGeantVertexYAll->SetYTitle(" Counts ");

    TH1F *hGeantVertexZAll = new TH1F("hGeantVertexZAll", "", 1000, -100, 100);
    hGeantVertexZAll->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertexZAll->SetYTitle(" Counts ");

    TH1F *hGeantVertexXLambda = new TH1F("hGeantVertexXLambda", "", 1000, -100, 100);
    hGeantVertexXLambda->SetXTitle("Geant Vertex, X / mm");
    hGeantVertexXLambda->SetYTitle(" Counts ");

    TH1F *hGeantVertexYLambda = new TH1F("hGeantVertexYLambda", "", 1000, -100, 100);
    hGeantVertexYLambda->SetXTitle("Geant Vertex, Y / mm");
    hGeantVertexYLambda->SetYTitle(" Counts ");

    TH1F *hGeantVertexZLambda = new TH1F("hGeantVertexZLambda", "", 1000, -100, 100);
    hGeantVertexZLambda->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertexZLambda->SetYTitle(" Counts ");

    TH1F *hGeantVertexXCand1 = new TH1F("hGeantVertexXCand1", "", 1000, -100, 100);
    hGeantVertexXCand1->SetXTitle("Geant Vertex, X / mm");
    hGeantVertexXCand1->SetYTitle(" Counts ");

    TH1F *hGeantVertexYCand1 = new TH1F("hGeantVertexYCand1", "", 1000, -100, 100);
    hGeantVertexYCand1->SetXTitle("Geant Vertex, Y / mm");
    hGeantVertexYCand1->SetYTitle(" Counts ");

    TH1F *hGeantVertexZCand1 = new TH1F("hGeantVertexZCand1", "", 1000, -60, 1000);
    hGeantVertexZCand1->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertexZCand1->SetYTitle(" Counts ");

    TH1F *hGeantVertexXCand2 = new TH1F("hGeantVertexXCand2", "", 1000, -100, 100);
    hGeantVertexXCand2->SetXTitle("Geant Vertex, X / mm");
    hGeantVertexXCand2->SetYTitle(" Counts ");

    TH1F *hGeantVertexYCand2 = new TH1F("hGeantVertexYCand2", "", 1000, -100, 100);
    hGeantVertexYCand2->SetXTitle("Geant Vertex, Y / mm");
    hGeantVertexYCand2->SetYTitle(" Counts ");

    TH1F *hGeantVertexZCand2 = new TH1F("hGeantVertexZCand2", "", 1000, -60, 1000);
    hGeantVertexZCand2->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertexZCand2->SetYTitle(" Counts ");

    TH1F *hGeantVertexXCand1Prim = new TH1F("hGeantVertexXCand1Prim", "", 1000, -100, 100);
    hGeantVertexXCand1Prim->SetXTitle("Geant Vertex, X / mm");
    hGeantVertexXCand1Prim->SetYTitle(" Counts ");

    TH1F *hGeantVertexYCand1Prim = new TH1F("hGeantVertexYCand1Prim", "", 1000, -100, 100);
    hGeantVertexYCand1Prim->SetXTitle("Geant Vertex, Y / mm");
    hGeantVertexYCand1Prim->SetYTitle(" Counts ");

    TH1F *hGeantVertexZCand1Prim = new TH1F("hGeantVertexZCand1Prim", "", 1000, -60, 1000);
    hGeantVertexZCand1Prim->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertexZCand1Prim->SetYTitle(" Counts ");

    TH1F *hGeantVertexXCand3 = new TH1F("hGeantVertexXCand3", "", 1000, -100, 100);
    hGeantVertexXCand3->SetXTitle("Geant Vertex, X / mm");
    hGeantVertexXCand3->SetYTitle(" Counts ");

    TH1F *hGeantVertexYCand3 = new TH1F("hGeantVertexYCand3", "", 1000, -100, 100);
    hGeantVertexYCand3->SetXTitle("Geant Vertex, Y / mm");
    hGeantVertexYCand3->SetYTitle(" Counts ");

    TH1F *hGeantVertexZCand3 = new TH1F("hGeantVertexZCand3", "", 1000, -60, 1000);
    hGeantVertexZCand3->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertexZCand3->SetYTitle(" Counts ");

    TH2F *hGeantVertex = new TH2F("hGeantVertex", "", 1000, -100, 100, 1000, -100, 100);
    hGeantVertex->SetXTitle("Geant Vertex, Z / mm");
    hGeantVertex->SetYTitle("Geant Vertex, R / mm");

    TH1F *hVertexXDiff = new TH1F("hVertexXDiff", "", 1000, -100, 100);
    hVertexXDiff->SetXTitle("Difference Vertex, X / mm");
    hVertexXDiff->SetYTitle(" Counts ");

    TH1F *hVertexYDiff = new TH1F("hVertexYDiff", "", 1000, -100, 100);
    hVertexYDiff->SetXTitle("Difference Vertex, Y / mm");
    hVertexYDiff->SetYTitle(" Counts ");

    TH1F *hVertexZDiff = new TH1F("hVertexZDiff", "", 1000, -100, 100);
    hVertexZDiff->SetXTitle("Difference Vertex, Z / mm");
    hVertexZDiff->SetYTitle(" Counts ");

    TH1F *hVertexXDiffPrim = new TH1F("hVertexXDiffPrim", "", 1000, -100, 100);
    hVertexXDiffPrim->SetXTitle("Difference Vertex, X / mm");
    hVertexXDiffPrim->SetYTitle(" Counts ");

    TH1F *hVertexYDiffPrim = new TH1F("hVertexYDiffPrim", "", 1000, -100, 100);
    hVertexYDiffPrim->SetXTitle("Difference Vertex, Y / mm");
    hVertexYDiffPrim->SetYTitle(" Counts ");

    TH1F *hVertexZDiffPrim = new TH1F("hVertexZDiffPrim", "", 1000, -100, 100);
    hVertexZDiffPrim->SetXTitle("Difference Vertex, Z / mm");
    hVertexZDiffPrim->SetYTitle(" Counts ");

    TH1F *hVertexXDiff_ProbCut = (TH1F *)hVertexXDiff->Clone("hVertexXDiff_ProbCut");
    TH1F *hVertexYDiff_ProbCut = (TH1F *)hVertexYDiff->Clone("hVertexYDiff_ProbCut");
    TH1F *hVertexZDiff_ProbCut = (TH1F *)hVertexZDiff->Clone("hVertexZDiff_ProbCut");
    TH1F *hVertexXDiffPrim_ProbCut = (TH1F *)hVertexXDiffPrim->Clone("hVertexXDiffPrim_ProbCut");
    TH1F *hVertexYDiffPrim_ProbCut = (TH1F *)hVertexYDiffPrim->Clone("hVertexYDiffPrim_ProbCut");
    TH1F *hVertexZDiffPrim_ProbCut = (TH1F *)hVertexZDiffPrim->Clone("hVertexZDiffPrim_ProbCut");

    TH1F *hVertexXDiff_BestComb = (TH1F *)hVertexXDiff->Clone("hVertexXDiff_BestComb");
    TH1F *hVertexYDiff_BestComb = (TH1F *)hVertexYDiff->Clone("hVertexYDiff_BestComb");
    TH1F *hVertexZDiff_BestComb = (TH1F *)hVertexZDiff->Clone("hVertexZDiff_BestComb");

    TH1F *hVertexXDiffPrim_BestComb = (TH1F *)hVertexXDiffPrim->Clone("hVertexXDiffPrim_BestComb");
    TH1F *hVertexYDiffPrim_BestComb = (TH1F *)hVertexYDiffPrim->Clone("hVertexYDiffPrim_BestComb");
    TH1F *hVertexZDiffPrim_BestComb = (TH1F *)hVertexZDiffPrim->Clone("hVertexZDiffPrim_BestComb");

    TH1F *hVertexXDiff_ProbCut_BestComb = (TH1F *)hVertexXDiff->Clone("hVertexXDiff_ProbCut_BestComb");
    TH1F *hVertexYDiff_ProbCut_BestComb = (TH1F *)hVertexYDiff->Clone("hVertexYDiff_ProbCut_BestComb");
    TH1F *hVertexZDiff_ProbCut_BestComb = (TH1F *)hVertexZDiff->Clone("hVertexZDiff_ProbCut_BestComb");

    TH1F *hVertexXDiffPrim_ProbCut_BestComb = (TH1F *)hVertexXDiffPrim->Clone("hVertexXDiffPrim_ProbCut_BestComb");
    TH1F *hVertexYDiffPrim_ProbCut_BestComb = (TH1F *)hVertexYDiffPrim->Clone("hVertexYDiffPrim_ProbCut_BestComb");
    TH1F *hVertexZDiffPrim_ProbCut_BestComb = (TH1F *)hVertexZDiffPrim->Clone("hVertexZDiffPrim_ProbCut_BestComb");

    TH1F *hVertexXDiff_BothVerticesFound = (TH1F *)hVertexXDiff->Clone("hVertexXDiff_BothVerticesFound");
    TH1F *hVertexYDiff_BothVerticesFound = (TH1F *)hVertexYDiff->Clone("hVertexYDiff_BothVerticesFound");
    TH1F *hVertexZDiff_BothVerticesFound = (TH1F *)hVertexZDiff->Clone("hVertexZDiff_BothVerticesFound");

    TH1F *hVertexXDiffPrim_BothVerticesFound = (TH1F *)hVertexXDiffPrim->Clone("hVertexXDiffPrim_BothVerticesFound");
    TH1F *hVertexYDiffPrim_BothVerticesFound = (TH1F *)hVertexYDiffPrim->Clone("hVertexYDiffPrim_BothVerticesFound");
    TH1F *hVertexZDiffPrim_BothVerticesFound = (TH1F *)hVertexZDiffPrim->Clone("hVertexZDiffPrim_BothVerticesFound");

    // ---- Momentum -----
    TH1F *hGeantTotMomentumProtons = new TH1F("hGeantTotMomentumProtons", "", 1000, 0, 3000);
    hGeantTotMomentumProtons->SetXTitle("Momentum, X / MeV/c");
    hGeantTotMomentumProtons->SetYTitle(" Counts ");
    TH1F *hGeantXMomentumProtons = new TH1F("hGeantXMomentumProtons", "", 1000, -1000, 1000);
    hGeantXMomentumProtons->SetXTitle("Momentum, X / MeV/c");
    hGeantXMomentumProtons->SetYTitle(" Counts ");
    TH1F *hGeantYMomentumProtons = new TH1F("hGeantYMomentumProtons", "", 1000, -1000, 1000);
    hGeantYMomentumProtons->SetXTitle("Momentum, Y / MeV/c");
    hGeantYMomentumProtons->SetYTitle(" Counts ");
    TH1F *hGeantZMomentumProtons = new TH1F("hGeantZMomentumProtons", "", 1000, -100, 3000);
    hGeantZMomentumProtons->SetXTitle("Momentum, Z / MeV/c");
    hGeantZMomentumProtons->SetYTitle(" Counts ");

    TH1F *hGeantTotMomentumPions = new TH1F("hGeantTotMomentumPions", "", 1000, 0, 3000);
    hGeantTotMomentumPions->SetXTitle("Momentum, X / MeV/c");
    hGeantTotMomentumPions->SetYTitle(" Counts ");
    TH1F *hGeantXMomentumPions = new TH1F("hGeantXMomentumPions", "", 1000, -1000, 1000);
    hGeantXMomentumPions->SetXTitle("Momentum, X / MeV/c");
    hGeantXMomentumPions->SetYTitle(" Counts ");
    TH1F *hGeantYMomentumPions = new TH1F("hGeantYMomentumPions", "", 1000, -1000, 1000);
    hGeantYMomentumPions->SetXTitle("Momentum, Y / MeV/c");
    hGeantYMomentumPions->SetYTitle(" Counts ");
    TH1F *hGeantZMomentumPions = new TH1F("hGeantZMomentumPions", "", 1000, -100, 3000);
    hGeantZMomentumPions->SetXTitle("Momentum, Z / MeV/c");
    hGeantZMomentumPions->SetYTitle(" Counts ");

    TH1F *hGeantTotMomentumPrimProtons = new TH1F("hGeantTotMomentumPrimProtons", "", 1000, 0, 3000);
    hGeantTotMomentumPrimProtons->SetXTitle("Momentum, X / MeV/c");
    hGeantTotMomentumPrimProtons->SetYTitle(" Counts ");
    TH1F *hGeantXMomentumPrimProtons = new TH1F("hGeantXMomentumPrimProtons", "", 1000, -1000, 1000);
    hGeantXMomentumPrimProtons->SetXTitle("Momentum, X / MeV/c");
    hGeantXMomentumPrimProtons->SetYTitle(" Counts ");
    TH1F *hGeantYMomentumPrimProtons = new TH1F("hGeantYMomentumPrimProtons", "", 1000, -1000, 1000);
    hGeantYMomentumPrimProtons->SetXTitle("Momentum, Y / MeV/c");
    hGeantYMomentumPrimProtons->SetYTitle(" Counts ");
    TH1F *hGeantZMomentumPrimProtons = new TH1F("hGeantZMomentumPrimProtons", "", 1000, -100, 3000);
    hGeantZMomentumPrimProtons->SetXTitle("Momentum, Z / MeV/c");
    hGeantZMomentumPrimProtons->SetYTitle(" Counts ");

    TH1F *hGeantTotMomentumKaons = new TH1F("hGeantTotMomentumKaons", "", 1000, 0, 3000);
    hGeantTotMomentumKaons->SetXTitle("Momentum, X / MeV/c");
    hGeantTotMomentumKaons->SetYTitle(" Counts ");
    TH1F *hGeantXMomentumKaons = new TH1F("hGeantXMomentumKaons", "", 1000, -1000, 1000);
    hGeantXMomentumKaons->SetXTitle("Momentum, X / MeV/c");
    hGeantXMomentumKaons->SetYTitle(" Counts ");
    TH1F *hGeantYMomentumKaons = new TH1F("hGeantYMomentumKaons", "", 1000, -1000, 1000);
    hGeantYMomentumKaons->SetXTitle("Momentum, Y / MeV/c");
    hGeantYMomentumKaons->SetYTitle(" Counts ");
    TH1F *hGeantZMomentumKaons = new TH1F("hGeantZMomentumKaons", "", 1000, -100, 3000);
    hGeantZMomentumKaons->SetXTitle("Momentum, Z / MeV/c");
    hGeantZMomentumKaons->SetYTitle(" Counts ");

    // ------- Reconstructed quantities -------
    TH1F *hRecoRProtons = new TH1F("hRecoRProtons", "", 1000, -100, 100);
    hRecoRProtons->SetXTitle("R / mm");
    hRecoRProtons->SetYTitle(" Counts ");
    TH1F *hRecoZProtons = new TH1F("hRecoZProtons", "", 1000, -100, 100);
    hRecoZProtons->SetXTitle("Z / mm");
    hRecoZProtons->SetYTitle(" Counts ");
    TH1F *hRecoThetaProtons = new TH1F("hRecoThetaProtons", "", 500, 0, 3);
    hRecoThetaProtons->SetXTitle("#theta / rad");
    hRecoThetaProtons->SetYTitle(" Counts ");
    TH1F *hRecoPhiProtons = new TH1F("hRecoPhiProtons", "", 500, -4, 4);
    hRecoPhiProtons->SetXTitle("#phi / rad");
    hRecoPhiProtons->SetYTitle(" Counts ");
    TH1F *hRecoMomentumProtons = new TH1F("hRecoMomentumProtons", "", 1000, 0, 3000);
    hRecoMomentumProtons->SetXTitle("Momentum / MeV/c");
    hRecoMomentumProtons->SetYTitle(" Counts ");
    TH1F *hRecoBetaProtons = new TH1F("hRecoBetaProtons", "", 100, -0.1, 1);
    hRecoBetaProtons->SetXTitle("#beta");
    hRecoBetaProtons->SetYTitle(" Counts ");

    TH1F *hRecoRPions = new TH1F("hRecoRPions", "", 1000, -100, 100);
    hRecoRPions->SetXTitle("R / mm");
    hRecoRPions->SetYTitle(" Counts ");
    TH1F *hRecoZPions = new TH1F("hRecoZPions", "", 1000, -100, 100);
    hRecoZPions->SetXTitle("Z / mm");
    hRecoZPions->SetYTitle(" Counts ");
    TH1F *hRecoThetaPions = new TH1F("hRecoThetaPions", "", 500, 0, 3);
    hRecoThetaPions->SetXTitle("#theta / rad");
    hRecoThetaPions->SetYTitle(" Counts ");
    TH1F *hRecoPhiPions = new TH1F("hRecoPhiPions", "", 500, -4, 4);
    hRecoPhiPions->SetXTitle("#phi / rad");
    hRecoPhiPions->SetYTitle(" Counts ");
    TH1F *hRecoMomentumPions = new TH1F("hRecoMomentumPions", "", 1000, 0, 3000);
    hRecoMomentumPions->SetXTitle("Momentum / MeV/c");
    hRecoMomentumPions->SetYTitle(" Counts ");
    TH1F *hRecoBetaPions = new TH1F("hRecoBetaPions", "", 100, -0.1, 1);
    hRecoBetaPions->SetXTitle("#beta");
    hRecoBetaPions->SetYTitle(" Counts ");

    // ------- Reconstructed quantities Primaries -----------
    TH1F *hRecoRPrimProtons = new TH1F("hRecoRPrimProtons", "", 1000, -100, 100);
    hRecoRPrimProtons->SetXTitle("R / mm");
    hRecoRPrimProtons->SetYTitle(" Counts ");
    TH1F *hRecoZPrimProtons = new TH1F("hRecoZPrimProtons", "", 1000, -100, 100);
    hRecoZPrimProtons->SetXTitle("Z / mm");
    hRecoZPrimProtons->SetYTitle(" Counts ");
    TH1F *hRecoThetaPrimProtons = new TH1F("hRecoThetaPrimProtons", "", 500, 0, 3);
    hRecoThetaPrimProtons->SetXTitle("#theta / rad");
    hRecoThetaPrimProtons->SetYTitle(" Counts ");
    TH1F *hRecoPhiPrimProtons = new TH1F("hRecoPhiPrimProtons", "", 500, -4, 4);
    hRecoPhiPrimProtons->SetXTitle("#phi / rad");
    hRecoPhiPrimProtons->SetYTitle(" Counts ");
    TH1F *hRecoMomentumPrimProtons = new TH1F("hRecoMomentumPrimProtons", "", 1000, 0, 3000);
    hRecoMomentumPrimProtons->SetXTitle("Momentum / MeV/c");
    hRecoMomentumPrimProtons->SetYTitle(" Counts ");

    TH1F *hRecoRKaons = new TH1F("hRecoRKaons", "", 1000, -100, 100);
    hRecoRKaons->SetXTitle("R / mm");
    hRecoRKaons->SetYTitle(" Counts ");
    TH1F *hRecoZKaons = new TH1F("hRecoZKaons", "", 1000, -100, 100);
    hRecoZKaons->SetXTitle("Z / mm");
    hRecoZKaons->SetYTitle(" Counts ");
    TH1F *hRecoThetaKaons = new TH1F("hRecoThetaKaons", "", 500, 0, 3);
    hRecoThetaKaons->SetXTitle("#theta / rad");
    hRecoThetaKaons->SetYTitle(" Counts ");
    TH1F *hRecoPhiKaons = new TH1F("hRecoPhiKaons", "", 500, -4, 4);
    hRecoPhiKaons->SetXTitle("#phi / rad");
    hRecoPhiKaons->SetYTitle(" Counts ");
    TH1F *hRecoMomentumKaons = new TH1F("hRecoMomentumKoans", "", 1000, 0, 3000);
    hRecoMomentumKaons->SetXTitle("Momentum / MeV/c");
    hRecoMomentumKaons->SetYTitle(" Counts ");

    // ---------------- LAMBDA PLOTS --------------------
    TH1F *hMomLambda = new TH1F("hMomLambda", "", 1000, 0, 4000);
    hMomLambda->SetXTitle("Momentum / MeV/c");
    hMomLambda->SetYTitle(" Counts ");
    TH1F *hMomLambda_AfterFit = (TH1F *)hMomLambda->Clone("hMomLambda_AfterFit");
    TH1F *hRecoThetaLambda = new TH1F("hRecoThetaLambda", "", 500, 0, 3);
    hRecoThetaLambda->SetXTitle("#theta / rad");
    hRecoThetaLambda->SetYTitle(" Counts ");
    TH1F *hRecoThetaLambda_AfterFit = (TH1F *)hRecoThetaLambda->Clone("hRecoThetaLambda_AfterFit");
    TH1F *hRecoThetaLambdaZCut = new TH1F("hRecoThetaLambdaZCut", "", 500, 0, 3);
    hRecoThetaLambdaZCut->SetXTitle("#theta / rad");
    hRecoThetaLambdaZCut->SetYTitle(" Counts ");

    TH1F *hRecoPhiLambda = new TH1F("hRecoPhiLambda", "", 500, -4, 4);
    hRecoPhiLambda->SetXTitle("#phi / rad");
    hRecoPhiLambda->SetYTitle(" Counts ");
    TH1F *hRecoPhiLambda_AfterFit = (TH1F *)hRecoPhiLambda->Clone("hRecoPhiLambda_AfterFit");
    TH1F *hRecoPhiLambdaZCut = new TH1F("hRecoPhiLambdaZCut", "", 500, -4, 4);
    hRecoPhiLambdaZCut->SetXTitle("#phi / rad");
    hRecoPhiLambdaZCut->SetYTitle(" Counts ");

    TH1F *hRecoRLambda = new TH1F("hRecoRLambda", "", 1000, -100, 100);
    hRecoRLambda->SetXTitle("R / mm");
    hRecoRLambda->SetYTitle(" Counts ");
    TH1F *hRecoZLambda = new TH1F("hRecoZLambda", "", 1000, -100, 100);
    hRecoZLambda->SetXTitle("Z / mm");
    hRecoZLambda->SetYTitle(" Counts ");
    TH1F *hRecoRLambdaZCut = new TH1F("hRecoRLambdaZCut", "", 1000, -100, 100);
    hRecoRLambdaZCut->SetXTitle("R / mm");
    hRecoRLambdaZCut->SetYTitle(" Counts ");
    TH1F *hRecoZLambdaZCut = new TH1F("hRecoZLambdaZCut", "", 1000, -100, 100);
    hRecoZLambdaZCut->SetXTitle("Z / mm");
    hRecoZLambdaZCut->SetYTitle(" Counts ");

    TH1F *hErrorRLambda = new TH1F("hErrorRLambda", "", 1000, -100, 100);
    hErrorRLambda->SetXTitle("R / mm");
    hErrorRLambda->SetYTitle(" Counts ");

    TH1F *hErrorZLambda = new TH1F("hErrorZLambda", "", 1000, -100, 100);
    hErrorZLambda->SetXTitle("Z / mm");
    hErrorZLambda->SetYTitle(" Counts ");
    TH1F *hErrorThetaLambda = new TH1F("hErrorThetaLambda", "", 500, 0, 3);
    hErrorThetaLambda->SetXTitle("Error #theta / rad");
    hErrorThetaLambda->SetYTitle(" Counts ");
    //TH1F *hErrorThetaLambda_AfterFit = (TH1F*)hErrorThetaLambda->Clone("hErrorThetaLambda_AfterFit");
    TH1F *hErrorThetaLambda_AfterFit = new TH1F("hErrorThetaLambda_AfterFit", "", 1000, 0, 0.01);
    hErrorThetaLambda_AfterFit->SetXTitle("Error #theta / rad");
    hErrorThetaLambda_AfterFit->SetYTitle(" Counts ");
    TH1F *hErrorPhiLambda = new TH1F("hErrorPhiLambda", "", 500, 0, 3);
    hErrorPhiLambda->SetXTitle("Error #phi / rad");
    hErrorPhiLambda->SetYTitle(" Counts ");
    //TH1F *hErrorPhiLambda_AfterFit = (TH1F*)hErrorPhiLambda->Clone("hErrorPhiLambda_AfterFit");
    TH1F *hErrorPhiLambda_AfterFit = new TH1F("hErrorPhiLambda_AfterFit", "", 1000, 0, 0.01);
    hErrorPhiLambda_AfterFit->SetXTitle("Error #phi / rad");
    hErrorPhiLambda_AfterFit->SetYTitle(" Counts ");

    // ---------------- LAMBDA PLOTS CUT --------------------
    TH1F *hMomLambdaCut = new TH1F("hMomLambdaCut", "", 1000, 0, 3000);
    hMomLambdaCut->SetXTitle("Momentum / MeV/c");
    hMomLambdaCut->SetYTitle(" Counts ");
    TH1F *hRecoThetaLambdaCut = new TH1F("hRecoThetaLambdaCut", "", 500, 0, 3);
    hRecoThetaLambdaCut->SetXTitle("#theta / rad");
    hRecoThetaLambdaCut->SetYTitle(" Counts ");
    TH1F *hRecoPhiLambdaCut = new TH1F("hRecoPhiLambdaCut", "", 500, -4, 4);
    hRecoPhiLambdaCut->SetXTitle("#phi / rad");
    hRecoPhiLambdaCut->SetYTitle(" Counts ");

    TH1F *hErrorRLambdaCut = new TH1F("hErrorRLambdaCut", "", 1000, -100, 100);
    hErrorRLambdaCut->SetXTitle("R / mm");
    hErrorRLambdaCut->SetYTitle(" Counts ");
    TH1F *hErrorZLambdaCut = new TH1F("hErrorZLambdaCut", "", 1000, -100, 100);
    hErrorZLambdaCut->SetXTitle("Z / mm");
    hErrorZLambdaCut->SetYTitle(" Counts ");
    TH1F *hErrorThetaLambdaCut = new TH1F("hErrorThetaLambdaCut", "", 500, 0, 3);
    hErrorThetaLambdaCut->SetXTitle("#theta / rad");
    hErrorThetaLambdaCut->SetYTitle(" Counts ");
    TH1F *hErrorPhiLambdaCut = new TH1F("hErrorPhiLambdaCut", "", 500, -4, 4);
    hErrorPhiLambdaCut->SetXTitle("#phi / rad");
    hErrorPhiLambdaCut->SetYTitle(" Counts ");

    // Vertex prim and decay info
    TH1F *hDistPrimToDecayVertex = new TH1F("hDistPrimToDecayVertex", "", 100, 0, 150);
    hDistPrimToDecayVertex->SetXTitle(" Distance Between Primary and Decay Vertex ");
    hDistPrimToDecayVertex->SetYTitle(" Counts ");

    TH1F *hvertex_z = new TH1F("hvertex_z", "", 200, -100, 100);
    TH1F *hvertex_x = new TH1F("hvertex_x", "", 200, -50, 50);
    TH1F *hvertex_y = new TH1F("hvertex_y", "", 200, -50, 50);

    TH1F *hErrorMomProton_AfterFit = new TH1F("hErrorMomProton_AfterFit", "", 1000, 0, 0.0001);
    hErrorMomProton_AfterFit->SetXTitle("Error Momentum / MeV/c");
    hErrorMomProton_AfterFit->SetYTitle(" Counts ");

    TH1F *hErrorMomPion_AfterFit = new TH1F("hErrorMomPion_AfterFit", "", 1000, 0, 0.0001);
    hErrorMomPion_AfterFit->SetXTitle("Error Momentum / MeV/c");
    hErrorMomPion_AfterFit->SetYTitle(" Counts ");

    TH1F *hErrorThetaProton_AfterFit = new TH1F("hErrorThetaProton_AfterFit", "", 10000, 0, 0.01);
    hErrorThetaProton_AfterFit->SetXTitle("Error #theta / rad");
    hErrorThetaProton_AfterFit->SetYTitle(" Counts ");

    TH1F *hErrorThetaPion_AfterFit = new TH1F("hErrorThetaPion_AfterFit", "", 10000, 0, 0.01);
    hErrorThetaPion_AfterFit->SetXTitle("Error #theta / rad");
    hErrorThetaPion_AfterFit->SetYTitle(" Counts ");

    TH1F *hErrorPhiProton_AfterFit = new TH1F("hErrorPhiProton_AfterFit", "", 10000, 0, 0.01);
    hErrorPhiProton_AfterFit->SetXTitle("Error #theta / rad");
    hErrorPhiProton_AfterFit->SetYTitle(" Counts ");

    TH1F *hErrorPhiPion_AfterFit = new TH1F("hErrorPhiPion_AfterFit", "", 10000, 0, 0.1);
    hErrorPhiPion_AfterFit->SetXTitle("Error #theta / rad");
    hErrorPhiPion_AfterFit->SetYTitle(" Counts ");

    HLoop loop(kTRUE);
    Bool_t ret = loop.addFiles(infileList);
    if (ret == 0)
    {
        cout << "READBACK: ERROR : cannot find inputfiles : "
             << infileList.Data() << endl;
        return 1;
    }

    // select categories here
    if (!loop.setInput("-*,+HParticleCandSim,+HGeantKine"))
    {
        cout << "READBACK: ERROR : cannot read input !" << endl;
        exit(1);
    } // read all categories

    loop.printCategories();
    loop.printChain();

    HCategory *catParticle = loop.getCategory("HParticleCandSim");
    if (!catParticle)
    {
        std::cout << "No particleCat in input!" << std::endl;
        exit(1);
    }
    HCategory *catGeant = loop.getCategory("HGeantKine");
    if (!catGeant)
    {
        std::cout << "No kineCat in input!" << std::endl;
        exit(1);
    }

    Int_t entries = loop.getEntries();
    if (nEvents < entries && nEvents >= 0)
        entries = nEvents;

    int primVertexBeforeDecayVertex = 0;
    int decayVertexBeforePrimVertex = 0;
    int primVertexInsideDecayVertex = 0;
    int decayVertexInsidePrimVertex = 0;
    bool primVtxFound, secVtxFound;

    double deg2rad = TMath::DegToRad();

    // start of the event loop
    for (Int_t i = 1; i < nEvents; i++)
    {
        //----------break if last event is reached-------------
        if (loop.nextEvent(i) <= 0)
        {
            cout << " end recieved " << endl;
            break;
        } // last event reached
        HTool::printProgress(i, nEvents, 1, "Analysing evt# :");

        // for each event there are a number of tracks
        Int_t ntracks = catParticle->getEntries();

        std::vector<HRefitCand> protons, pions, kaons;
        std::vector<HParticleCandSim *> virtualCandLambdas, virtualCandProtons, virtualCandPions, virtualCandKaons;

        HEventHeader *eventheader = gHades->getCurrentEvent()->getHeader();
        HVertex evtVertex = eventheader->getVertexReco();
        hvertex_x->Fill(evtVertex.getX());
        hvertex_y->Fill(evtVertex.getY());
        hvertex_z->Fill(evtVertex.getZ());

        std::vector<HParticleCandSim *> candVector;
        candVector.clear();

        for (Int_t k = 0; k < ntracks; k++)
        {
            HParticleCandSim *cand =
                HCategoryManager::getObject(cand, catParticle, k);

            if (cand->isGhostTrack())
                continue;
            // select "good" tracks
            if (!cand->isFlagBit(Particle::kIsUsed))
                continue;

            candVector.push_back(cand);
        }

        HDecayBuilder decayBuilder(candVector);

        for (Int_t j = 0; j < ntracks; j++)
        {
            HParticleCandSim *cand =
                HCategoryManager::getObject(cand, catParticle, j);
            // skip ghost tracks (only avalible for MC events)
            if (cand->isGhostTrack())
                continue;
            // select "good" tracks
            if (!cand->isFlagBit(Particle::kIsUsed))
                continue;

            hGeantVertexXAll->Fill(cand->getGeantxVertex());
            hGeantVertexYAll->Fill(cand->getGeantyVertex());
            hGeantVertexZAll->Fill(cand->getGeantzVertex());

            HRefitCand candidate(cand);

            // select particles based on MC info
            // proton pdg==14, pion pdg==9, k+ pdg==11, lambda pdg==18
            // error values obtained from resoultion plots

            // The code below makes sure that all protons and pions in the analysis come from the Lambda decay
            // No combinatorial background

            /*if (cand->getGeantPID() == 14 && cand->getGeantParentPID() == 18) //Proton found
            {
                virtualCandProtons.push_back(cand);
                double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                                   1.188, 2.652};
                FillData(cand, candidate, errors, 938.272);
                protons.push_back(candidate);
            }
            else if (cand->getGeantPID() == 9 && cand->getGeantParentPID() == 18) // Pion found
            {
                virtualCandPions.push_back(cand);
                double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                                   4.006, 7.629};
                FillData(cand, candidate, errors, 139.570);
                pions.push_back(candidate);
            }*/

            if (cand->getGeantPID() == 14) //Proton found
            {
                //std::cout << "Parent ID: " << cand->getGeantParentPID() << std::endl;
                if (selectHadrons(cand) == true)
                {
                    virtualCandProtons.push_back(cand);
                    double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                                       1.188, 2.652};
                    FillData(cand, candidate, errors, 938.272);
                    protons.push_back(candidate);
                }
            }
            else if (cand->getGeantPID() == 9) // Pion found
            {
                if (selectHadrons(cand) == true)
                {
                    virtualCandPions.push_back(cand);
                    double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                                       4.006, 7.629};
                    FillData(cand, candidate, errors, 139.570);
                    pions.push_back(candidate);
                }
            }
            else if (cand->getGeantPID() == 11) // Kaon found
            {
                if (selectHadrons(cand) == true)
                {
                    virtualCandKaons.push_back(cand);
                    double errors[] = {1.947 * 1e-5, 2.296 * 1e-3, 6.312 * 1e-3,
                                       1.404, 2.723};
                    FillData(cand, candidate, errors, 493.7);
                    kaons.push_back(candidate);
                }
            }
            if (cand->getGeantPID() == 18) //Lambda found
            {
                // There are probably no candidates for Lambda due to track requirement
                virtualCandLambdas.push_back(cand);
            }
            else
                continue;
        } // end track loop

        /*if (cand->getGeantPID() == 14) //Proton found
            {
                virtualCandProtons.push_back(cand);
                Double_t mom = cand->getGeantTotalMom();
                //double errors[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                 //                  1.188, 2.652};   //rough error estimates
                double errors[] = {momErrP->Eval(mom), thtErrP->Eval(mom), phiErrP->Eval(mom),
                                   RErrP->Eval(mom), ZErrP->Eval(mom)};  //momentum dependent error estimates                 
                FillData(cand, candidate, errors, 938.272);
                protons.push_back(candidate);
            }
            else if (cand->getGeantPID() == 9) // Pion found
            {
                virtualCandPions.push_back(cand);
                Double_t mom = cand->getGeantTotalMom();
                //double errors[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                 //                  4.006, 7.629};   //rough error estimates
                double errors[] = {momErrPi->Eval(mom), thtErrPi->Eval(mom), phiErrPi->Eval(mom),
                                   RErrPi->Eval(mom), ZErrPi->Eval(mom)};    //momentum dependent error estimates                
                FillData(cand, candidate, errors, 139.570);
                pions.push_back(candidate);
            }
            else if (cand->getGeantPID() == 11) // Kaon found
            {
                virtualCandKaons.push_back(cand);
                Double_t mom = cand->getGeantTotalMom();
                //double errors[] = {1.947 * 1e-5, 2.296 * 1e-3, 6.312 * 1e-3,
                  //                 1.404, 2.723};   //rough error estimates
                double errors[] = {momErrK->Eval(mom), thtErrK->Eval(mom), phiErrK->Eval(mom),
                                   RErrK->Eval(mom), ZErrK->Eval(mom)};  //momentum dependent error estimates
                FillData(cand, candidate, errors, 493.7);
                kaons.push_back(candidate);
            }
            if (cand->getGeantPID() == 18) //Lambda found
            {
                // There are probably no candidates for Lambda due to track requirement
                virtualCandLambdas.push_back(cand);
            }
            else
                continue; 
        } // end track loop  */

        // -----------------------------------------------------------------------
        // looking at Lambda invariant mass here
        // -----------------------------------------------------------------------

        //std::cout << "Event number: "  << i << std::endl;
        //std::cout << " " << std::endl;
        //std::cout << "Number of Kaons: " << kaons.size() << std::endl;
        //std::cout << "Number of protons: " << protons.size() << std::endl;
        //std::cout << "Number of Pions: " << pions.size() << std::endl;

        primVtxFound = false;
        secVtxFound = false;

        //if (protons.size() < 2)
        //  continue;

        TVector3 primaryVertex, decayVertex;
        std::vector<HRefitCand> eventCandidates;
        eventCandidates.clear();
        double probPrim = -99999, probSec = -99999;

        bool bestDecayVertexFound = false;
        bool bestPrimVertexFound = false;

        TVector3 decayVertex_Temp, primVertex_Temp;
        TVector3 decayVertex_TempUpdated, tempVertex_TempUpdated;
        TVector3 decayVertexBestFit, primVertexBestFit;

        double probDecayVertex_Temp = -1;
        double probPrimVertex_Temp = -1;
        double probDecayVertex_TempUpdated = -1;
        double probPrimVertex_TempUpdated = -1;

        double bestDiffXDecay = -1;
        double bestDiffYDecay = -1;
        double bestDiffZDecay = -1;
        double bestDiffXPrim = -1;
        double bestDiffYPrim = -1;
        double bestDiffZPrim = -1;

        double indexPrimaryProton;
        double indexDecayProton;

        std::vector<HRefitCand> cands3c;
        cands3c.clear();

        for (size_t n = 0; n < protons.size(); n++)
        {

            HRefitCand cand1 = protons[n];
            HParticleCandSim *virtualCand1 = virtualCandProtons[n];

            for (size_t m = 0; m < pions.size(); m++)
            {

                HRefitCand cand2 = pions[m];
                HParticleCandSim *virtualCand2 = virtualCandPions[m];

                secVtxFound = true;
                std::vector<HRefitCand> candsSec;
                candsSec.clear();
                candsSec.push_back(cand1);
                candsSec.push_back(cand2);
                TLorentzVector lambda = cand1 + cand2;
                h045->Fill(lambda.M());
                //eventCandidates.push_back(cand1);
                //eventCandidates.push_back(cand2);

                HVertexFinder *vtxFinderSec = new HVertexFinder();
                decayVertex = vtxFinderSec->findVertex(candsSec);
                std::cout << "vertex position, old way: " << decayVertex.X() << " " << decayVertex.Y() << " " << decayVertex.Z() << std::endl;
                
                //decayVertex_Temp = decayVertex;

                hVertexXPreFit->Fill(decayVertex.X());
                hVertexYPreFit->Fill(decayVertex.Y());
                hVertexZPreFit->Fill(decayVertex.Z());
                hDistanceToVertexProtonPreFit->Fill(vtxFinderSec->getDistanceFirstParticleVertex());
                hDistanceToVertexPionPreFit->Fill(vtxFinderSec->getDistanceSecondParticleVertex());
                hDistanceBetweenProtonAndPionPreFit->Fill(vtxFinderSec->getDistanceBetweenFittedParticles());
                // Perform fitting of secondary vertex
                HKinFitter vtxFitterSecCands(candsSec);
                vtxFitterSecCands.addVertexConstraint();
                vtxFitterSecCands.fit();
                h02SecondaryVtx->Fill(vtxFitterSecCands.getChi2());
                h03SecondaryVtx->Fill(vtxFitterSecCands.getProb());
                h03Zoomed->Fill(vtxFitterSecCands.getProb());
                
                //if(vtxFitterSecCands.isConverged()){
                h06SecondaryVtx->Fill(vtxFitterSecCands.getPull(1));
                h07SecondaryVtx->Fill(vtxFitterSecCands.getPull(2));
                h08SecondaryVtx->Fill(vtxFitterSecCands.getPull(3));
                h09SecondaryVtx->Fill(vtxFitterSecCands.getPull(4));
                //}

                probSec = vtxFitterSecCands.getProb();
                //std::cout << "Secondary vertex probability: " << probSec << std::endl;

                eventCandidates.push_back(vtxFitterSecCands.getDaughter(0));
                eventCandidates.push_back(vtxFitterSecCands.getDaughter(1));

                if (probSec > probDecayVertex_Temp)
                {

                    bestDecayVertexFound = true;
                    probDecayVertex_Temp = probSec;
                    decayVertexBestFit = decayVertex;
                    indexDecayProton = n;
                    bestDiffXDecay = decayVertex.X() - virtualCand2->getGeantxVertex();
                    bestDiffYDecay = decayVertex.Y() - virtualCand2->getGeantyVertex();
                    bestDiffZDecay = decayVertex.Z() - virtualCand2->getGeantzVertex();

                    cands3c.clear();

                    // Get the proton daughter and pass it to the 3C fit later
                    cands3c.push_back(vtxFitterSecCands.getDaughter(0));

                    // Get the pion daughter and pass it to the 3C fit later
                    cands3c.push_back(vtxFitterSecCands.getDaughter(1));
                }

                h11SecVtx->Fill(vtxFitterSecCands.getIteration());

                //if (probSec > 0.000001)
                if (probSec > 0)
                {

                    h13SecondaryVtx->Fill(vtxFitterSecCands.getPull(1));
                    h14SecondaryVtx->Fill(vtxFitterSecCands.getPull(2));
                    h15SecondaryVtx->Fill(vtxFitterSecCands.getPull(3));
                    h16SecondaryVtx->Fill(vtxFitterSecCands.getPull(4));

                    hVertexXDiff_ProbCut->Fill(decayVertex.X() - virtualCand2->getGeantxVertex());
                    hVertexYDiff_ProbCut->Fill(decayVertex.Y() - virtualCand2->getGeantyVertex());
                    hVertexZDiff_ProbCut->Fill(decayVertex.Z() - virtualCand2->getGeantzVertex());

                    if (vtxFitterSecCands.isConverged())
                    {
                        TLorentzVector lambdaVtxFit = vtxFitterSecCands.getDaughter(0) + vtxFitterSecCands.getDaughter(1);

                        h040->Fill(lambdaVtxFit.M());
                    }
                }

                hGeantVertexXCand1->Fill(virtualCand1->getGeantxVertex());
                hGeantVertexYCand1->Fill(virtualCand1->getGeantyVertex());
                hGeantVertexZCand1->Fill(virtualCand1->getGeantzVertex());

                hGeantVertexXCand2->Fill(virtualCand2->getGeantxVertex());
                hGeantVertexYCand2->Fill(virtualCand2->getGeantyVertex());
                hGeantVertexZCand2->Fill(virtualCand2->getGeantzVertex());

                hVertexXDiff->Fill(decayVertex.X() - virtualCand2->getGeantxVertex());
                hVertexYDiff->Fill(decayVertex.Y() - virtualCand2->getGeantyVertex());
                hVertexZDiff->Fill(decayVertex.Z() - virtualCand2->getGeantzVertex());

                hGeantTotMomentumProtons->Fill(virtualCand1->getGeantTotalMom());
                hGeantXMomentumProtons->Fill(virtualCand1->getGeantxMom());
                hGeantYMomentumProtons->Fill(virtualCand1->getGeantyMom());
                hGeantZMomentumProtons->Fill(virtualCand1->getGeantzMom());
                hGeantTotMomentumPions->Fill(virtualCand2->getGeantTotalMom());
                hGeantXMomentumPions->Fill(virtualCand2->getGeantxMom());
                hGeantYMomentumPions->Fill(virtualCand2->getGeantyMom());
                hGeantZMomentumPions->Fill(virtualCand2->getGeantzMom());

                hRecoRProtons->Fill(cand1.getR());
                hRecoZProtons->Fill(cand1.getZ());
                hRecoThetaProtons->Fill(cand1.Theta());
                hRecoPhiProtons->Fill(cand1.Phi());
                hRecoMomentumProtons->Fill(cand1.P());
                //hRecoBetaProtons->Fill(cand1.getBeta());

                hRecoRPions->Fill(cand2.getR());
                hRecoZPions->Fill(cand2.getZ());
                hRecoThetaPions->Fill(cand2.Theta());
                hRecoPhiPions->Fill(cand2.Phi());
                hRecoMomentumPions->Fill(cand2.P());
            }

            for (size_t p = 0; p < kaons.size(); p++)
            {

                HRefitCand cand3 = kaons[p];
                HParticleCandSim *virtualCand3 = virtualCandKaons[p];

                primVtxFound = true;
                std::vector<HRefitCand> candsPrim;
                candsPrim.clear();
                candsPrim.push_back(cand1);
                candsPrim.push_back(cand3);

                //eventCandidates.push_back(cand1);
                //eventCandidates.push_back(cand3);

                HVertexFinder *vtxFinderPrim = new HVertexFinder();
                primaryVertex = vtxFinderPrim->findVertex(candsPrim);

                hPrimaryVertexXPreFit->Fill(primaryVertex.X());
                hPrimaryVertexYPreFit->Fill(primaryVertex.Y());
                hPrimaryVertexZPreFit->Fill(primaryVertex.Z());
                hDistanceBetweenPrimaryProtonAndKaonPreFit->Fill(vtxFinderPrim->getDistanceBetweenFittedParticles());
                hDistanceToVertexPrimaryProtonPreFit->Fill(vtxFinderPrim->getDistanceFirstParticleVertex());
                hDistanceToVertexKaonPreFit->Fill(vtxFinderPrim->getDistanceSecondParticleVertex());

                hGeantTotMomentumPrimProtons->Fill(virtualCand1->getGeantTotalMom());
                hGeantXMomentumPrimProtons->Fill(virtualCand1->getGeantxMom());
                hGeantYMomentumPrimProtons->Fill(virtualCand1->getGeantyMom());
                hGeantZMomentumPrimProtons->Fill(virtualCand1->getGeantzMom());
                hGeantTotMomentumKaons->Fill(virtualCand3->getGeantTotalMom());
                hGeantXMomentumKaons->Fill(virtualCand3->getGeantxMom());
                hGeantYMomentumKaons->Fill(virtualCand3->getGeantyMom());
                hGeantZMomentumKaons->Fill(virtualCand3->getGeantzMom());

                hRecoRPrimProtons->Fill(cand1.getR());
                hRecoZPrimProtons->Fill(cand1.getZ());
                hRecoThetaPrimProtons->Fill(cand1.Theta());
                hRecoPhiPrimProtons->Fill(cand1.Phi());
                hRecoMomentumPrimProtons->Fill(cand1.P());

                hRecoRKaons->Fill(cand3.getR());
                hRecoZKaons->Fill(cand3.getZ());
                hRecoThetaKaons->Fill(cand3.Theta());
                hRecoPhiKaons->Fill(cand3.Phi());
                hRecoMomentumKaons->Fill(cand3.P());
                // Perform fitting of primary vertex
                HKinFitter vtxFitterPrimCands(candsPrim);
                //vtxFitterPrimCands.setNumberOfIterations(1);
                vtxFitterPrimCands.addVertexConstraint();
                vtxFitterPrimCands.fit();

                h02->Fill(vtxFitterPrimCands.getChi2());
                h03->Fill(vtxFitterPrimCands.getProb());

                probPrim = vtxFitterPrimCands.getProb();

                eventCandidates.push_back(vtxFitterPrimCands.getDaughter(0));
                eventCandidates.push_back(vtxFitterPrimCands.getDaughter(1));

                if (probPrim > probPrimVertex_Temp)
                {

                    bestPrimVertexFound = true;
                    probPrimVertex_Temp = probPrim;
                    primVertexBestFit = primaryVertex;
                    indexPrimaryProton = n;
                    bestDiffXPrim = primaryVertex.X() - virtualCand3->getGeantxVertex();
                    bestDiffYPrim = primaryVertex.Y() - virtualCand3->getGeantyVertex();
                    bestDiffZPrim = primaryVertex.Z() - virtualCand3->getGeantzVertex();
                }

                h11->Fill(vtxFitterPrimCands.getIteration());

                //if(vtxFitterPrimCands.isConverged()){

                h06->Fill(vtxFitterPrimCands.getPull(1));
                h07->Fill(vtxFitterPrimCands.getPull(2));
                h08->Fill(vtxFitterPrimCands.getPull(3));
                h09->Fill(vtxFitterPrimCands.getPull(4));

                //}

                hGeantVertexXCand1Prim->Fill(virtualCand1->getGeantxVertex());
                hGeantVertexYCand1Prim->Fill(virtualCand1->getGeantyVertex());
                hGeantVertexZCand1Prim->Fill(virtualCand1->getGeantzVertex());

                hGeantVertexXCand3->Fill(virtualCand3->getGeantxVertex());
                hGeantVertexYCand3->Fill(virtualCand3->getGeantyVertex());
                hGeantVertexZCand3->Fill(virtualCand3->getGeantzVertex());

                hVertexXDiffPrim->Fill(primaryVertex.X() - virtualCand3->getGeantxVertex());
                hVertexYDiffPrim->Fill(primaryVertex.Y() - virtualCand3->getGeantyVertex());
                hVertexZDiffPrim->Fill(primaryVertex.Z() - virtualCand3->getGeantzVertex());

                //if (probPrim > 0.000001)
                if (probPrim > 0)
                {

                    hVertexXDiffPrim_ProbCut->Fill(primaryVertex.X() - virtualCand3->getGeantxVertex());
                    hVertexYDiffPrim_ProbCut->Fill(primaryVertex.Y() - virtualCand3->getGeantyVertex());
                    hVertexZDiffPrim_ProbCut->Fill(primaryVertex.Z() - virtualCand3->getGeantzVertex());

                    h13->Fill(vtxFitterPrimCands.getPull(1));
                    h14->Fill(vtxFitterPrimCands.getPull(2));
                    h15->Fill(vtxFitterPrimCands.getPull(3));
                    h16->Fill(vtxFitterPrimCands.getPull(4));
                }
            }
        }

        if (bestDecayVertexFound == true)
        {

            hVertexXDiff_BestComb->Fill(bestDiffXDecay);
            hVertexYDiff_BestComb->Fill(bestDiffYDecay);
            hVertexZDiff_BestComb->Fill(bestDiffZDecay);

            if (probDecayVertex_Temp > 0.000001)
            {

                hVertexXDiff_ProbCut_BestComb->Fill(bestDiffXDecay);
                hVertexYDiff_ProbCut_BestComb->Fill(bestDiffYDecay);
                hVertexZDiff_ProbCut_BestComb->Fill(bestDiffZDecay);
            }
        }

        if (bestPrimVertexFound == true)
        {

            hVertexXDiffPrim_BestComb->Fill(bestDiffXPrim);
            hVertexYDiffPrim_BestComb->Fill(bestDiffYPrim);
            hVertexZDiffPrim_BestComb->Fill(bestDiffZPrim);

            if (probPrimVertex_Temp > 0.000001)
            {

                hVertexXDiffPrim_ProbCut_BestComb->Fill(bestDiffXPrim);
                hVertexYDiffPrim_ProbCut_BestComb->Fill(bestDiffYPrim);
                hVertexZDiffPrim_ProbCut_BestComb->Fill(bestDiffZPrim);
            }
        }

        if (bestPrimVertexFound == true && bestDecayVertexFound == true)
        {
            if (indexDecayProton != indexPrimaryProton)
            {

                double distPrimToDecayVertex = sqrt((decayVertexBestFit.X() - primVertexBestFit.X()) * (decayVertexBestFit.X() - primVertexBestFit.X()) + (decayVertexBestFit.Y() - primVertexBestFit.Y()) * (decayVertexBestFit.Y() - primVertexBestFit.Y()) + (decayVertexBestFit.Z() - primVertexBestFit.Z()) * (decayVertexBestFit.Z() - primVertexBestFit.Z()));

                hDistPrimToDecayVertex->Fill(distPrimToDecayVertex);

                hVertexXDiff_BothVerticesFound->Fill(bestDiffXDecay);
                hVertexYDiff_BothVerticesFound->Fill(bestDiffYDecay);
                hVertexZDiff_BothVerticesFound->Fill(bestDiffZDecay);

                hVertexXDiffPrim_BothVerticesFound->Fill(bestDiffXPrim);
                hVertexYDiffPrim_BothVerticesFound->Fill(bestDiffYPrim);
                hVertexZDiffPrim_BothVerticesFound->Fill(bestDiffZPrim);

                double R_primaryVertex, R_decayVertex;

                R_primaryVertex = sqrt(primVertexBestFit.X() * primVertexBestFit.X() + primVertexBestFit.Y() * primVertexBestFit.Y());
                R_decayVertex = sqrt(decayVertexBestFit.X() * decayVertexBestFit.X() + decayVertexBestFit.Y() * decayVertexBestFit.Y());

                hVertexXPostFit->Fill(decayVertexBestFit.X());
                hVertexYPostFit->Fill(decayVertexBestFit.Y());
                hVertexZPostFit->Fill(decayVertexBestFit.Z());

                hPrimVertexXPostFit->Fill(primVertexBestFit.X());
                hPrimVertexYPostFit->Fill(primVertexBestFit.Y());
                hPrimVertexZPostFit->Fill(primVertexBestFit.Z());

                HNeutralCandFinder lambdaCandFinder(cands3c);

                lambdaCandFinder.setUsePrimaryVertexInNeutralMotherCalculation(true);

                lambdaCandFinder.setNeutralMotherCandFromPrimaryVtxInfo(primaryVertex, decayVertex);

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
                hMomLambda->Fill(lambdaCandRefit.P());

                lambdaCandRefit.setR(lambdaCand.getR());
                lambdaCandRefit.setZ(lambdaCand.getZ());

                lambdaCandRefit.SetTheta(lambdaCand.getTheta() * deg2rad);
                lambdaCandRefit.SetPhi(lambdaCand.getPhi() * deg2rad);

                hRecoRLambda->Fill(lambdaCand.getR());
                hRecoZLambda->Fill(lambdaCand.getZ());
                hRecoThetaLambda->Fill(lambdaCand.getTheta() * deg2rad);
                hRecoPhiLambda->Fill(lambdaCand.getPhi() * deg2rad);

                // Take the lambda mass before the fit
                TLorentzVector lambdaCandBefore3C = cands3c[0] + cands3c[1];
                hLambdaMassBeforeFit->Fill(lambdaCandBefore3C.M());

                if (cands3c.size() == 2)
                {
                    HKinFitter Fitter3c(cands3c, lambdaCandRefit);
                    Fitter3c.add3Constraint();
                    Fitter3c.setNumberOfIterations(20);

                    Fitter3c.fit();
                    if (Fitter3c.isConverged())
                    {
                        //HRefitCand fcand1 = Fitter3c.getDaughter(0); // proton
                        //HRefitCand fcand2 = Fitter3c.getDaughter(1); // pion
                        HRefitCand flambda = Fitter3c.getMother();

                        HRefitCand cand13C = Fitter3c.getDaughter(0); // proton
                        HRefitCand cand23C = Fitter3c.getDaughter(1); // pion
                        TLorentzVector lambdaCand3C = cand13C + cand23C;
                        h04->Fill(lambdaCand3C.M());
                    }

                    h021->Fill(Fitter3c.getChi2());
                    h031->Fill(Fitter3c.getProb());
                }

                if (R_primaryVertex < R_decayVertex)
                {

                    primVertexInsideDecayVertex++;
                }
                else
                {
                    decayVertexInsidePrimVertex++;
                }

                if (primVertexBestFit.Z() < decayVertexBestFit.Z())
                {

                    hRecoRLambdaZCut->Fill(lambdaCand.getR());
                    hRecoZLambdaZCut->Fill(lambdaCand.getZ());
                    hRecoThetaLambdaZCut->Fill(lambdaCand.getTheta() * deg2rad);
                    hRecoPhiLambdaZCut->Fill(lambdaCand.getPhi() * deg2rad);

                    primVertexBeforeDecayVertex++;
                }
                else
                {
                    decayVertexBeforePrimVertex++;
                }
            }
        }

    } // end of the events loop

    std::cout << "Number of prim vertices before the decay vertex: " << primVertexBeforeDecayVertex << ", and after: " << decayVertexBeforePrimVertex << std::endl;
    std::cout << "Number of prim vertices inside the decay vertex: " << primVertexInsideDecayVertex << ", and outside: " << decayVertexInsidePrimVertex << std::endl;

    // write histograms to the output file
    outfile->cd();
    h01->Write();
    h040->Write();
    h02->Write();
    h03->Write();
    h03Zoomed->Write();
    h021->Write();
    h031->Write();
    h041->Write();
    h044->Write();
    h054->Write(); // get Pull example (1/P for the fitted proton)
    h114->Write();
    h023->Write();
    h033->Write();
    h111->Write();
    h02SecondaryVtx->Write();
    h03SecondaryVtx->Write();
    h045->Write();
    h04->Write();
    hLambdaMassBeforeFit->Write();
    h05->Write();
    h06SecondaryVtx->Write();
    h07SecondaryVtx->Write();
    h08SecondaryVtx->Write();
    h09SecondaryVtx->Write();
    h06->Write();
    h07->Write();
    h08->Write();
    h09->Write();

    h10->Write();
    h11->Write();
    h11SecVtx->Write();
    h12->Write();
    h13->Write();
    h14->Write();
    h15->Write();
    h16->Write();
    h13SecondaryVtx->Write();
    h14SecondaryVtx->Write();
    h15SecondaryVtx->Write();
    h16SecondaryVtx->Write();
    h17->Write();
    h18->Write();
    h19->Write();
    h20->Write();
    h21->Write();
    h22->Write();
    h23->Write();
    h24->Write();

    hVertexXPreFit->Write();
    hVertexYPreFit->Write();
    hVertexZPreFit->Write();
    hPrimaryVertexXPreFit->Write();
    hPrimaryVertexYPreFit->Write();
    hPrimaryVertexZPreFit->Write();

    hVertexPreFit->Write();
    hDistanceBetweenProtonAndPionPreFit->Write();
    hDistanceBetweenPrimaryProtonAndKaonPreFit->Write();
    hDistanceToVertexProtonPreFit->Write();
    hDistanceToVertexPionPreFit->Write();
    hDistanceToVertexPrimaryProtonPreFit->Write();
    hDistanceToVertexKaonPreFit->Write();
    hDistanceToOriginProtonPreFit->Write();
    hDistanceToOriginPionPreFit->Write();
    hDistanceToOriginSumPreFit->Write();
    hPrimVertexXPostFit->Write();
    hPrimVertexYPostFit->Write();
    hPrimVertexZPostFit->Write();

    hVertexXPostFit->Write();
    hVertexYPostFit->Write();
    hVertexZPostFit->Write();
    hVertexPostFit->Write();

    hDistanceBetweenProtonAndPionPostFit->Write();
    hDistanceToVertexProtonPostFit->Write();
    hDistanceToVertexPionPostFit->Write();

    hGeantVertexXAll->Write();
    hGeantVertexYAll->Write();
    hGeantVertexZAll->Write();

    hGeantVertexXLambda->Write();
    hGeantVertexYLambda->Write();
    hGeantVertexZLambda->Write();

    hGeantVertexXCand1->Write();
    hGeantVertexYCand1->Write();
    hGeantVertexZCand1->Write();

    hGeantVertexXCand2->Write();
    hGeantVertexYCand2->Write();
    hGeantVertexZCand2->Write();

    hVertexXDiff->Write();
    hVertexYDiff->Write();
    hVertexZDiff->Write();

    hVertexXDiff_ProbCut->Write();
    hVertexYDiff_ProbCut->Write();
    hVertexZDiff_ProbCut->Write();

    hVertexXDiff_BestComb->Write();
    hVertexYDiff_BestComb->Write();
    hVertexZDiff_BestComb->Write();

    hVertexXDiff_ProbCut_BestComb->Write();
    hVertexYDiff_ProbCut_BestComb->Write();
    hVertexZDiff_ProbCut_BestComb->Write();

    hVertexXDiff_BothVerticesFound->Write();
    hVertexYDiff_BothVerticesFound->Write();
    hVertexZDiff_BothVerticesFound->Write();

    hGeantVertexXCand1Prim->Write();
    hGeantVertexYCand1Prim->Write();
    hGeantVertexZCand1Prim->Write();

    hGeantVertexXCand3->Write();
    hGeantVertexYCand3->Write();
    hGeantVertexZCand3->Write();

    hVertexXDiffPrim->Write();
    hVertexYDiffPrim->Write();
    hVertexZDiffPrim->Write();

    hVertexXDiffPrim_ProbCut->Write();
    hVertexYDiffPrim_ProbCut->Write();
    hVertexZDiffPrim_ProbCut->Write();

    hVertexXDiffPrim_BestComb->Write();
    hVertexYDiffPrim_BestComb->Write();
    hVertexZDiffPrim_BestComb->Write();

    hVertexXDiffPrim_ProbCut_BestComb->Write();
    hVertexYDiffPrim_ProbCut_BestComb->Write();
    hVertexZDiffPrim_ProbCut_BestComb->Write();

    hVertexXDiffPrim_BothVerticesFound->Write();
    hVertexYDiffPrim_BothVerticesFound->Write();
    hVertexZDiffPrim_BothVerticesFound->Write();

    hGeantTotMomentumProtons->Write();
    hGeantXMomentumProtons->Write();
    hGeantYMomentumProtons->Write();
    hGeantZMomentumProtons->Write();
    hGeantTotMomentumPions->Write();
    hGeantXMomentumPions->Write();
    hGeantYMomentumPions->Write();
    hGeantZMomentumPions->Write();

    hGeantTotMomentumPrimProtons->Write();
    hGeantXMomentumPrimProtons->Write();
    hGeantYMomentumPrimProtons->Write();
    hGeantZMomentumPrimProtons->Write();
    hGeantTotMomentumKaons->Write();
    hGeantXMomentumKaons->Write();
    hGeantYMomentumKaons->Write();
    hGeantZMomentumKaons->Write();

    hRecoRProtons->Write();
    hRecoZProtons->Write();
    hRecoThetaProtons->Write();
    hRecoPhiProtons->Write();
    hRecoMomentumProtons->Write();
    hRecoBetaProtons->Write();

    hRecoRPions->Write();
    hRecoZPions->Write();
    hRecoThetaPions->Write();
    hRecoPhiPions->Write();
    hRecoMomentumPions->Write();
    hRecoBetaPions->Write();

    hRecoRPrimProtons->Write();
    hRecoZPrimProtons->Write();
    hRecoThetaPrimProtons->Write();
    hRecoPhiPrimProtons->Write();
    hRecoMomentumPrimProtons->Write();

    hRecoRKaons->Write();
    hRecoZKaons->Write();
    hRecoThetaKaons->Write();
    hRecoPhiKaons->Write();
    hRecoMomentumKaons->Write();

    h25->Write();
    h26->Write();
    h27->Write();
    h28->Write();
    h29->Write();
    h30->Write();
    h31->Write();
    h32->Write();
    h33->Write();

    hMomLambda->Write();
    hRecoThetaLambda->Write();
    hRecoPhiLambda->Write();
    hRecoRLambda->Write();
    hRecoZLambda->Write();
    hRecoThetaLambdaZCut->Write();
    hRecoPhiLambdaZCut->Write();
    hRecoRLambdaZCut->Write();
    hRecoZLambdaZCut->Write();

    hErrorRLambda->Write();
    hErrorZLambda->Write();
    hErrorThetaLambda->Write();
    hErrorThetaLambda_AfterFit->Write();
    hErrorPhiLambda->Write();
    hErrorPhiLambda_AfterFit->Write();

    hMomLambdaCut->Write();
    hRecoThetaLambdaCut->Write();
    hRecoPhiLambdaCut->Write();

    hErrorRLambdaCut->Write();
    hErrorZLambdaCut->Write();
    hErrorThetaLambdaCut->Write();
    hErrorPhiLambdaCut->Write();

    hDistPrimToDecayVertex->Write();

    hvertex_x->Write();
    hvertex_y->Write();
    hvertex_z->Write();

    hErrorMomProton_AfterFit->Write();
    hErrorThetaProton_AfterFit->Write();
    hErrorPhiProton_AfterFit->Write();

    hErrorMomPion_AfterFit->Write();
    hErrorThetaPion_AfterFit->Write();
    hErrorPhiPion_AfterFit->Write();

    hMomLambda_AfterFit->Write();
    hRecoThetaLambda_AfterFit->Write();
    hRecoPhiLambda_AfterFit->Write();

    // hMomProton_AFterFit->Write();
    // hRecoThetaProton_AFterFit->Write();
    // hRecoPhiProton_AFterFit->Write();
    // hMomPion_AFterFit->Write();
    // hRecoThetaPion_AFterFit->Write();
    // hRecoPhiPion_AFterFit->Write();

    outfile->Close();

    return 0;
} // end of the macro
