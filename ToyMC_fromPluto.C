
TVector3 calcPoca(TVector3 vertex, TVector3 mom){
  Float_t px, py, pz;
  Float_t u = 0., t=0.;
  
  px = mom.X()/mom.Mag();
  py = mom.Y()/mom.Mag();
  pz = mom.Z()/mom.Mag();
  u = - (px*vertex.X() + py*vertex.Y())/(px*px + py*py);
  t = vertex.Z() + u*pz;
  TVector3 result(vertex.X()+u*px, vertex.Y()+u*py, t);
  return result;
}


void ToyMC_fromPluto(TString inFile, Int_t nEvents = 500000)
{
   TFile *outFile = new TFile("toy_montecarlo_rz_vtx.root", "recreate");

   TNtuple *ntuple = new TNtuple("data", "output flat tree", "p1CandTrueP:p1CandTrueTheta:p1CandTruePhi:p1CandTrueR:p1CandTrueZ:"
                                                             "p1CandRecoP:p1CandRecoTheta:p1CandRecoPhi:p1CandRecoR:p1CandRecoZ:"
                                                             "KCandTrueP:KCandTrueTheta:KCandTruePhi:KCandTrueR:KCandTrueZ:"
                                                             "KCandRecoP:KCandRecoTheta:KCandRecoPhi:KCandRecoR:KCandRecoZ:"
                                                             "p2CandTrueP:p2CandTrueTheta:p2CandTruePhi:p2CandTrueR:p2CandTrueZ:"
                                                             "p2CandRecoP:p2CandRecoTheta:p2CandRecoPhi:p2CandRecoR:p2CandRecoZ:"
                                                             "piCandTrueP:piCandTrueTheta:piCandTruePhi:piCandTrueR:piCandTrueZ:"
                                                             "piCandRecoP:piCandRecoTheta:piCandRecoPhi:piCandRecoR:piCandRecoZ:"
                                                             "trueDecVtxX:trueDecVtxY:trueDecVtxZ");

   
   TClonesArray *cands = new TClonesArray("PParticle");
   Float_t p1CandTrueP, p1CandTrueTheta, p1CandTrueR, p1CandTrueZ, p1CandRecoP, p1CandRecoTheta, p1CandRecoPhi, p1CandRecoR, p1CandRecoZ,
            KCandTrueP, KCandTrueTheta, KCandTrueR, KCandTrueZ, KCandRecoP, KCandRecoTheta, KCandRecoPhi,  KCandRecoR, KCandRecoZ,
            p2CandTrueP, p2CandTrueTheta, p2CandTrueR, p2CandTrueZ, p2CandRecoP, p2CandRecoTheta, p2CandRecoPhi,  p2CandRecoR, p2CandRecoZ,
            piCandTrueP, piCandTrueTheta, piCandTrueR, piCandTrueZ, piCandRecoP, piCandRecoTheta, piCandRecoPhi,  piCandRecoR, piCandRecoZ,
            trueDecVtxX, trueDecVtxY, trueDecVtxZ;
    
   TFile tree_file(inFile, "READ");
   TTree *t = (TTree*)tree_file.Get("data");

   t->SetBranchAddress("Particles", &cands);

   PParticle *p1 = new PParticle();
   PParticle *K = new PParticle();
   PParticle *p2 = new PParticle();
   PParticle *pi = new PParticle();
   
   // masses of the measured particles (p, K, p, pi-)
   //const Double_t masses[] = {0.938272, 0.493677, 0.938272, 0.13957};
   
   // add noise
   TRandom3 *noise = new TRandom3();
   // smearing values (these values are up to the user)
   // they correspong to p, theta, phi
   Double_t smear[5] = {0.025, 0.0009, 0.0009, 0.5, 1}; //x 0.1 more realistic, 

   Int_t entries = t->GetEntries();
   cout<<"Entries: "<<entries<<endl;
   if(nEvents<entries && nEvents>0) entries = nEvents;

   //Tree loop
   for(Int_t i=0; i<entries; i++){

      t->GetEntry(i);

      p1 = (PParticle *)cands->At(0);
      K = (PParticle *)cands->At(1);
      pi = (PParticle *)cands->At(2);
      p2 = (PParticle *)cands->At(3);

      Float_t p1Cand_p = (1. * (p1->P()) )/ (1 + (p1->P() * noise->Gaus(0, (1/p1->P()) * smear[0])));
      Float_t KCand_p = (1. * (K->P())) / (1 + (K->P() * noise->Gaus(0, (1/K->P()) * smear[0])));
      Float_t p2Cand_p = (1. * (p2->P())) / (1 + (p2->P() * noise->Gaus(0, (1/p2->P()) * smear[0])));
      Float_t piCand_p = (1. * (pi->P())) / (1 + (pi->P() * noise->Gaus(0, (1/pi->P()) * smear[0])));

      Float_t p1Cand_theta = (1. * (p1->Theta())) + noise->Gaus(0, smear[1]);
      Float_t KCand_theta = (1. * (K->Theta())) + noise->Gaus(0, smear[1]);
      Float_t p2Cand_theta = (1. * (p2->Theta())) + noise->Gaus(0, smear[1]);
      Float_t piCand_theta = (1. * (pi->Theta())) + noise->Gaus(0, smear[1]);

      Float_t p1Cand_phi = (1. * (p1->Phi())) + noise->Gaus(0, smear[2]);
      Float_t KCand_phi = (1. * (K->Phi())) + noise->Gaus(0, smear[2]);
      Float_t p2Cand_phi = (1. * (p2->Phi())) + noise->Gaus(0, smear[2]);
      Float_t piCand_phi = (1. * (pi->Phi())) + noise->Gaus(0, smear[2]);
      
      TVector3 pi_mom(pi->Px(), pi->Py(), pi->Pz());
      TVector3 pi_vtx(pi->X(), pi->Y(), pi->Z());
      TVector3 p2_mom(p2->Px(), p2->Py(), p2->Pz());
      TVector3 p2_vtx(p2->X(), p2->Y(), p2->Z());
      TVector3 poca_p2 = calcPoca(p2_vtx, p2_mom);
      TVector3 poca_pi = calcPoca(pi_vtx, pi_mom);

      Float_t p1Cand_R = (1. * 0) + noise->Gaus(0, smear[3]);
      Float_t KCand_R = (1. * 0) + noise->Gaus(0, smear[3]);
      Float_t p2Cand_R = (1. * (poca_p2.X()/cos(p2Cand_phi+M_PI/2))) + noise->Gaus(0, smear[3]);
      Float_t piCand_R = (1. * (poca_pi.X()/cos(piCand_phi+M_PI/2))) + noise->Gaus(0, smear[3]);

      Float_t p1Cand_Z = (1. * 0) + noise->Gaus(0, smear[4]);
      Float_t KCand_Z = (1. * 0) + noise->Gaus(0, smear[4]);
      Float_t p2Cand_Z = (1. * (poca_p2.Z())) + noise->Gaus(0, smear[4]);
      Float_t piCand_Z = (1. * (poca_pi.Z())) + noise->Gaus(0, smear[4]);


      // store in the tree
      // arranged as first all MC values then all reco values for all particles
      Float_t array[43] = {Float_t(p1->P()), Float_t(p1->Theta()), Float_t(p1->Phi()), 0., 0.,  p1Cand_p, p1Cand_theta, p1Cand_phi, p1Cand_R, p1Cand_Z,
                     Float_t(K->P()), Float_t(K->Theta()), Float_t(K->Phi()), 0., 0.,  KCand_p, KCand_theta, KCand_phi, KCand_R, KCand_Z,
                     Float_t(p2->P()), Float_t(p2->Theta()), Float_t(p2->Phi()), Float_t(poca_p2.X()/cos(p2Cand_phi+M_PI/2)), Float_t(poca_p2.Z()),  p2Cand_p, p2Cand_theta, p2Cand_phi, p2Cand_R, p2Cand_Z,
                     Float_t(pi->P()), Float_t(pi->Theta()), Float_t(pi->Phi()), Float_t(poca_pi.X()/cos(piCand_phi+M_PI/2)), Float_t(poca_pi.Z()),  piCand_p, piCand_theta, piCand_phi, piCand_R, piCand_Z,
                     Float_t(pi->X()), Float_t(pi->Y()), Float_t(pi->Z())};
      /*ntuple->Fill(p1->P(), p1->Theta(), p1->Phi(), 0, 0,  p1Cand_p, p1Cand_theta, p1Cand_phi, p1Cand_R, p1Cand_Z,
                     K->P(), K->Theta(), K->Phi(), 0, 0,  KCand_p, KCand_theta, KCand_phi, KCand_R, KCand_Z,
                     p2->P(), p2->Theta(), p2->Phi(), p2->R(), p2->Z(),  p2Cand_p, p2Cand_theta, p2Cand_phi, p2Cand_R, p2Cand_Z,
                     pi->P(), pi->Theta(), pi->Phi(), pi->R(), pi->Z(),  piCand_p, piCand_theta, piCand_phi, piCand_R, piCand_Z);*/
      ntuple->Fill(array);
   }

   outFile->cd();
   ntuple->Write();
   outFile->Close();
}