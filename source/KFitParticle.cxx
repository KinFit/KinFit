// ROOT includes
#include "KFitParticle.h"

KFitParticle::KFitParticle(TLorentzVector *cand, Double_t R, Double_t Z)
    : TLorentzVector(*cand), cand(cand), fMomentum(cand->P()), fTheta(cand->Theta()), fPhi(cand->Phi()), fR(R), fZ(Z)
{
    fPid = -1;
    fTrackId = -1;
    fCov.ResizeTo(5, 5);
}

KFitParticle::KFitParticle(TLorentzVector *cand, Double_t X, Double_t Y, Double_t Z)
    : TLorentzVector(*cand), cand(cand), fMomentum(cand->P()), fTheta(cand->Theta()), fPhi(cand->Phi())
{
    Double_t deg2rad = TMath::DegToRad();

    // Base and direction vettor of particle cand
    TVector3 base(X, Y, Z);
    TVector3 dir(TMath::Sin(fTheta) * TMath::Cos(fPhi),
                 TMath::Sin(fTheta) * TMath::Sin(fPhi),
                 TMath::Cos(fTheta));
    
    // Base and direction vector of beamline
    TVector3 beam_base(0,0,1); 
    TVector3 beam_dir(0,0,1);

    TVector3 cross = dir.Cross(beam_dir);
    TVector3 diff=base-beam_base;

    fR = abs(diff.Dot(cross)/(cross.Mag()));
    
    Double_t a = beam_base.Dot(beam_dir);
    Double_t b = beam_dir.Dot(beam_dir);
    Double_t c = base.Dot(beam_dir);
    Double_t d = (beam_base.Dot(dir)) * (dir.Dot(beam_dir)) / dir.Dot(dir);
    Double_t e = (beam_dir.Dot(dir)) * (dir.Dot(beam_dir)) / dir.Dot(dir);
    Double_t f = (base.Dot(dir)) * (dir.Dot(beam_dir)) / dir.Dot(dir);
    Double_t u1 = (-a + c + d - f) / (b - e);

    fZ = beam_base.Z() + beam_dir.Z() * u1;

    double y = beam_base.Y() + dir.Y() * u1;

    if(y < 0){

        fR = -1 * fR;
    }

    //HParticleTool::calcSegVector(1, 0, 0, 0, beam_base, beam_dir);
    //TVector3 POCA = HParticleTool::calculatePointOfClosestApproach(base, dir, beam_base, beam_dir);
    //fR = base.Y() * TMath::Cos(fPhi) - base.X() * TMath::Sin(fPhi);
    //fZ = POCA.Z();

    fPid = -1;
    fTrackId = -1;
}

KFitParticle::KFitParticle()
    : TLorentzVector()
{
    cand = new TLorentzVector();
    fR = 0;
    fZ = 0;
    fPid = -1;
    fTrackId = -1;
    fCov.ResizeTo(5, 5);
}

void KFitParticle::setCovariance(const TMatrixD &cov)
{
    // HADES default track parametrization
    // 0 = 1/p
    // 1 = theta
    // 2 = phi
    // 3 = R
    // 4 = z

    // An alternative parametrization
    // 0 = px
    // 1 = py
    // 2 = pz
    // 3 = x (vertex x position)
    // 4 = y (vertex y position)
    // 5 = z (vertex z position)

    if (cov.GetNoElements() == 5 * 5)
    {
        fCov.ResizeTo(5, 5);
        fCov = cov;
    }

    else if (cov.GetNoElements() == 6 * 6)
    {
        // TO-DO: do the error propagtion calculations
    }
}

void KFitParticle::reset()
{
    *(TLorentzVector *)this = *cand;
}

void KFitParticle::update()
{
    *((TLorentzVector *)cand) = *this;
}
