//#include "KFitParticle.h"

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

void FillData(KFitParticle **outcand, Double_t mom, Double_t tht, Double_t phi, Double_t R, Double_t Z, Double_t mass, Int_t pid, Int_t id)
{
    double deg2rad = TMath::DegToRad();

    (*outcand)->SetXYZM(mom * std::sin(tht) * std::cos(phi),
                     mom * std::sin(tht) * std::sin(phi),
                     mom * std::cos(tht), mass);

    (*outcand)->setMomentum(mom);
    (*outcand)->setTheta(tht);
    (*outcand)->setPhi(phi);
    (*outcand)->setR(R);
    (*outcand)->setZ(Z);
    (*outcand)->setPid(pid);
    (*outcand)->setTrackId(id);

    TMatrixD cov(5, 5);
    if(mom>0) cov(0, 0) = std::pow(0.025*(1/mom), 2);
    cov(1, 1) = std::pow(0.0009, 2);
    cov(2, 2) = std::pow(0.0009, 2);
    cov(3, 3) = std::pow(0.5, 2);
    cov(4, 4) = std::pow(1, 2);

    (*outcand)->setCovariance(cov);
    
}

void ToyMC_fromPluto_KFitParticle_workaround(TString inFile, Int_t nEvents = 500000)
{
   
   Float_t p1CandTrueP, p1CandTrueTheta, p1CandTruePhi, p1CandTrueR, p1CandTrueZ, p1CandRecoP, p1CandRecoTheta, p1CandRecoPhi, p1CandRecoR, p1CandRecoZ,
            KCandTrueP, KCandTrueTheta, KCandTruePhi, KCandTrueR, KCandTrueZ, KCandRecoP, KCandRecoTheta, KCandRecoPhi,  KCandRecoR, KCandRecoZ,
            p2CandTrueP, p2CandTrueTheta, p2CandTruePhi, p2CandTrueR, p2CandTrueZ, p2CandRecoP, p2CandRecoTheta, p2CandRecoPhi,  p2CandRecoR, p2CandRecoZ,
            piCandTrueP, piCandTrueTheta, piCandTruePhi, piCandTrueR, piCandTrueZ, piCandRecoP, piCandRecoTheta, piCandRecoPhi,  piCandRecoR, piCandRecoZ;
    
    TFile tree_file(inFile, "READ");
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

   //output tree
   TFile *outFile = new TFile("/home/jana/KinFit/input.root", "recreate");
    TTree *tree = new TTree("data", "input data for fitter");

    TClonesArray *p_array = new TClonesArray("KFitParticle");
    TClonesArray &p_arrayRef = *p_array;
    Int_t Event;

    TNtuple *ntuple = new TNtuple("data_true", "output flat tree", "p1CandTrueP:p1CandTrueTheta:p1CandTruePhi:p1CandTrueR:p1CandTrueZ:"
                                                             "KCandTrueP:KCandTrueTheta:KCandTruePhi:KCandTrueR:KCandTrueZ:"
                                                             "p2CandTrueP:p2CandTrueTheta:p2CandTruePhi:p2CandTrueR:p2CandTrueZ:"
                                                             "piCandTrueP:piCandTrueTheta:piCandTruePhi:piCandTrueR:piCandTrueZ:");

    tree->Branch("Event", &Event, "Event/I");
    tree->Branch("KFitParticle", "TClonesArray", &p_array);

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
   for(Int_t i=0; i<entries; i++)
   {
      Event = i;

      if (i % 10000 == 0)
         std::cout << "Processing Event " << i << " ... " << std::endl;
      
      t->GetEntry(i);

      KFitParticle *p1_fit = new (p_arrayRef[0]) KFitParticle();
      KFitParticle *K_fit = new (p_arrayRef[1]) KFitParticle();
      KFitParticle *pi_fit = new (p_arrayRef[2]) KFitParticle();
      KFitParticle *p2_fit = new (p_arrayRef[3]) KFitParticle();

      FillData(&p1_fit, p1CandRecoP, p1CandRecoTheta, p1CandRecoPhi, p1CandRecoR, p1CandRecoZ, 0.938272, 14, 1);
      FillData(&K_fit, KCandRecoP, KCandRecoTheta, KCandRecoPhi, KCandRecoR, KCandRecoZ, 0.493677, 11, 2);
      FillData(&pi_fit, piCandRecoP, piCandRecoTheta, piCandRecoPhi, piCandRecoR, piCandRecoZ, 0.13957, 9, 3);
      FillData(&p2_fit, p2CandRecoP, p2CandRecoTheta, p2CandRecoPhi, p2CandRecoR, p2CandRecoZ, 0.938272, 14, 4);

      // store in the tree
      // arranged as Proton MC values (p, theta, phi), Pion MC values (p, theta, phi)
      Float_t array[20] = {p1CandTrueP, p1CandTrueTheta, p1CandTruePhi, 0., 0.,
                     KCandTrueP, KCandTrueTheta, KCandTruePhi, 0., 0.,
                     piCandTrueP, piCandTrueTheta, piCandTruePhi, piCandTrueR, piCandTrueZ,
                     p2CandTrueP, p2CandTrueTheta, p2CandTruePhi, p2CandTrueR, p2CandTrueZ};

      ntuple->Fill(array);

      tree->Fill();
      p_array->Clear();

   } // end event loop

   outFile->cd();
   ntuple->Write();
   tree->Write();
   outFile->Close();
}