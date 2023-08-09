#include "plotsStyleMacro.C"

void plot_fitted(){

    TFile *infile = new TFile("~/KinFit_official/KinFit/fitted_vtx1.root", "read");
    TFile* outfile = new TFile("plotsoffitted_mass.root", "recreate");

    TH1F *prob = new TH1F("prob", "", 100, 0, 1);
    /*
    TH1F *combi = new TH1F("combi", "", 4, 1, 5);

    combi->GetXaxis()->SetBinLabel(1, "proton 1");
    combi->GetXaxis()->SetBinLabel(2, "kaon");
    combi->GetXaxis()->SetBinLabel(3, "pion");
    combi->GetXaxis()->SetBinLabel(4, "proton 2");
    combi->SetYTitle("counts");
    combi->SetBarWidth(0.9*combi->GetBinWidth(1));
*/
    TH1F *combi = new TH1F("combi", "", 2, 1, 3);

    combi->GetXaxis()->SetBinLabel(1, "proton 1");
    combi->GetXaxis()->SetBinLabel(2, "proton 2");
    combi->SetYTitle("counts");
    combi->SetBarWidth(0.9*combi->GetBinWidth(1));
    prob->SetXTitle("P(#chi^{2})");
    prob->SetYTitle("counts");
    prob->SetMinimum(0);

    TTree *tree = (TTree*)infile->Get("data_fitted");

    tree->Draw("KFitParticle.getTrackId()>>combi");
    tree->Draw("Prob>>prob");

    outfile->cd();
    plotsStyleMacro(prob, "prob_mass");
    plotsStyleMacro(combi, "selected_particles_mass", "bar");
    outfile->Close();
}