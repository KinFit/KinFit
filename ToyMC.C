
void ToyMC(Int_t nEvents = 500000)
{
   TFile *outFile = new TFile("toy_montecarlo.root", "recreate");

   TNtuple *ntuple = new TNtuple("data", "output flat tree", "pCandTrueP:pCandTrueTheta:pCandTruePhi:"
                                                             "piCandTrueP:piCandTrueTheta:piCandTruePhi:"
                                                             "pCandRecoP:pCandRecoTheta:pCandRecoPhi:"
                                                             "piCandRecoP:piCandRecoTheta:piCandRecoPhi");

   // create Lambda object and assign random values
   // for momentum components
   // (Momentum, Energy units are Gev/C, GeV)
   TLorentzVector Lambda;
   Lambda.SetXYZM(0.1, .1, 3., 1.11568);
   // masses of the decay products (p, pi-)
   const Double_t masses[] = {0.938272, 0.13957};

   // define the event object
   // and the decay chain
   TGenPhaseSpace event;
   event.SetDecay(Lambda, 2, masses);
   // add noise to momentum
   TRandom3 *noise = new TRandom3();
   // smearing values (these values are up to the user)
   // they correspong to p, theta, phi
   Double_t smear[3] = {0.025, 0.0009, 0.0009};

   // event loop
   for (Int_t j = 0; j < nEvents; j++)
   {
      if (j % 100000 == 0)
         std::cout << "Processing Event " << j << " ... " << std::endl;
      Double_t weight = event.Generate();
      TLorentzVector *pCand = event.GetDecay(0);
      TLorentzVector *piCand = event.GetDecay(1);

      auto pCand_p = (1. * (pCand->P())) + (noise->Gaus(0, (1/pCand->P()) * smear[0]));
      auto piCand_p = (1. * (piCand->P())) + (noise->Gaus(0, (1/piCand->P()) * smear[0]));

      auto pCand_theta = (1. * (pCand->Theta())) + noise->Gaus(0, smear[1]);
      auto piCand_theta = (1. * (piCand->Theta())) + noise->Gaus(0, smear[1]);

      auto pCand_phi = (1. * (pCand->Phi())) + noise->Gaus(0, smear[2]);
      auto piCand_phi = (1. * (piCand->Phi())) + noise->Gaus(0, smear[2]);

      // store in the tree
      // arranged as Proton MC values (p, theta, phi), Pion MC values (p, theta, phi)
      // then Proton "reco" values and Pion "reco values"
      ntuple->Fill(pCand->P(), pCand->Theta(), pCand->Phi(), piCand->P(), piCand->Theta(), piCand->Phi(),
                   pCand_p, pCand_theta, pCand_phi, piCand_p, piCand_theta, piCand_phi);

   } // end event loop
   outFile->cd();
   ntuple->Write();
   outFile->Close();
}