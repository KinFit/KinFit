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

using namespace std;

Int_t test_toyMC(TString infile, Int_t nEvents)
{
    // -----------------------------------------------------------------------
    // define output file and some histograms
    // -----------------------------------------------------------------------
    // set ouput file
    TFile* outfile = new TFile("test_toyMC.root", "recreate");

    TH1F* h_p_1op = new TH1F("h_p_1op", "", 100, -0.1, 0.1);
    h_p_1op->SetXTitle(" 1/p_{gen}-1/p_{reco} [c/GeV]");
    h_p_1op->SetYTitle(" events ");
    TH1F* h_p_tht = new TH1F("h_p_tht", "", 100, -0.01, 0.01);
    h_p_tht->SetXTitle(" #theta_{gen}-#theta_{rec}");
    h_p_tht->SetYTitle(" events ");
    TH1F* h_p_phi = new TH1F("h_p_phi", "", 100, -0.01, 0.01);
    h_p_phi->SetXTitle(" #phi_{gen}-#phi_{rec}");
    h_p_phi->SetYTitle(" events ");


    TH1F* h_pi_1op = new TH1F("h_pi_1op", "", 100, -0.6, 0.6);
    h_pi_1op->SetXTitle(" 1/p_{gen}-1/p_{reco} [c/GeV]");
    h_pi_1op->SetYTitle(" events ");
    TH1F* h_pi_tht = new TH1F("h_pi_tht", "", 100, -0.01, 0.01);
    h_pi_tht->SetXTitle(" #theta_{gen}-#theta_{rec}");
    h_pi_tht->SetYTitle(" events ");
    TH1F* h_pi_phi = new TH1F("h_pi_phi", "", 100, -0.01, 0.01);
    h_pi_phi->SetXTitle(" #phi_{gen}-#phi_{rec}");
    h_pi_phi->SetYTitle(" events ");

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

    Long64_t nevts = t->GetEntries();
    if (nEvents <= 0 || nEvents > nevts)
        nEvents = nevts;
    // event loop
    for (Long64_t ev = 0; ev < nEvents; ev++)
    {
        t->GetEntry(ev);

        h_p_1op->Fill(1/pCandTrueP - 1/pCandRecoP);
        h_p_tht->Fill(pCandTrueTheta - pCandRecoTheta);
        h_p_phi->Fill(pCandTruePhi - pCandRecoPhi);

        h_pi_1op->Fill(1/piCandTrueP - 1/piCandRecoP);
        h_pi_tht->Fill(piCandTrueTheta - piCandRecoTheta);
        h_pi_phi->Fill(piCandTruePhi - piCandRecoPhi);

    }

    // write histograms to the output file
    outfile->cd();

    h_p_1op->Write();
    h_p_tht->Write();
    h_p_phi->Write();

    h_pi_1op->Write();
    h_pi_tht->Write();
    h_pi_phi->Write();

    outfile->Close();

    return 0;
}