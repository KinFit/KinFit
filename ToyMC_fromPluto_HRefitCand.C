
void ToyMC_fromPluto_HRefitCand(TStrng inFile, Int_t nEvents = 500000)
{
   TFile *outFile = new TFile("toy_montecarlo.root", "recreate");

   TNtuple *ntuple = new TNtuple("data", "output flat tree", "pCandTrueP:pCandTrueTheta:pCandTruePhi:"
                                                             "piCandTrueP:piCandTrueTheta:piCandTruePhi:"
                                                             "pCandRecoP:pCandRecoTheta:pCandRecoPhi:"
                                                             "piCandRecoP:piCandRecoTheta:piCandRecoPhi");

   
   TClonesArray *cands = new TClonesArray("Particles");
   Float_t p1CandTrueP, p1CandTrueTheta, p1CandTrueR, p1CandTrueZ, p1CandRecoP, p1CandRecoTheta, p1CandRecoPhi, p1CandRecoR, p1CandRecoZ,
            KCandTrueP, KCandTrueTheta, KCandTrueR, KCandTrueZ, KCandRecoP, KCandRecoTheta, KCandRecoPhi,  KCandRecoR, KCandRecoZ,
            p2CandTrueP, p2CandTrueTheta, p2CandTrueR, p2CandTrueZ, p2CandRecoP, p2CandRecoTheta, p2CandRecoPhi,  p2CandRecoR, p2CandRecoZ,
            piCandTrueP, piCandTrueTheta, piCandTrueR, piCandTrueZ, piCandRecoP, piCandRecoTheta, piCandRecoPhi,  piCandRecoR, piCandRecoZ;
    
   TFile tree_file(inFile, "READ");
   TTree *t = (TTree*)tree_file.Get("data");

   t->SetBranchAddress("Particles", &cands);
   
   // masses of the measured particles (p, K, p, pi-)
   const Double_t masses[] = {0.938272, 0.493677, 0.938272, 0.13957};

   Int_t entries = t->GetEntries();
   cout<<"Entries: "<<entries<<endl;
   if(nEvents<entries && nEvents>0) entries = nEvents;

   //Tree loop
   for(Int_t i=0; i<entries; i++){

      t->GetEntry(i);

      Particles p1 = cands->At(0);
      Particles K = cands->At(1);
      Particles p2 = cands->At(2);
      Particles pi = cands->At(3);
   }

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

      auto pCand_p = (1. * (pCand->P())) / (1 + (pCand->P() * noise->Gaus(0, (1/pCand->P()) * smear[0])));
      auto piCand_p = (1. * (piCand->P())) / (1 + (piCand->P() * noise->Gaus(0, (1/piCand->P()) * smear[0])));
      //auto pCand_p = (1. * (pCand->P())) / (1 + (pCand->P() * noise->Gaus(0, smear[0])));
      //auto piCand_p = (1. * (piCand->P())) / (1 + (piCand->P() * noise->Gaus(0, smear[0])));

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