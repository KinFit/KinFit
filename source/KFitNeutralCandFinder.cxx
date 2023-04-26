#include "KFitNeutralCandFinder.h"

KFitNeutralCandFinder::KFitNeutralCandFinder(const std::vector<KFitParticle> &cands, TVector3 primaryVertex, TVector3 decayVertex) : fCands(cands), fVerbose(0), fMomentumBeforeDecay(-1.), fPrimaryVertex(primaryVertex), fDecayVertex(decayVertex), fNeutralCandMass(1115.683), fPrimVtxResX(1.78590), fPrimVtxResY(1.75516), fPrimVtxResZ(3.00431), fDecVtxResX(5.75369), fDecVtxResY(5.57198), fDecVtxResZ(10.2602)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitNeutralCandFinder -----------------" << std::endl;
    }
}

KFitNeutralCandFinder::KFitNeutralCandFinder(const std::vector<KFitParticle> &cands, double neutralCandMass, TVector3 primaryVertex, TVector3 decayVertex, double primVtxResX, double primVtxResY, double primVtxResZ, double decVtxResX, double decVtxResY, double decVtxResZ) : fCands(cands), fVerbose(0), fMomentumBeforeDecay(-1.), fPrimaryVertex(primaryVertex), fDecayVertex(decayVertex), fNeutralCandMass(neutralCandMass), fPrimVtxResX(primVtxResX), fPrimVtxResY(primVtxResY), fPrimVtxResZ(primVtxResZ), fDecVtxResX(decVtxResX), fDecVtxResY(decVtxResY), fDecVtxResZ(decVtxResZ)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitNeutralCandFinder -----------------" << std::endl;
    }
}

void KFitNeutralCandFinder::calculateNeutralMotherCand()
{
    Double_t param_p1, param_p2;

    KFitParticle cand1 = fCands[0];

    param_p1 = cand1.P(); // Not the inverse, this momentum is used for estimating the momentum of the Lambda Candidate

    KFitParticle cand2 = fCands[1];

    param_p2 = cand2.P(); // Not the inverse, this momentum is used for estimating the momentum of the Lambda Candidate

    Double_t energy_cand1, energy_cand2;
    energy_cand1 = std::sqrt(param_p1 * param_p1 + cand1.M() * cand1.M());
    energy_cand2 = std::sqrt(param_p2 * param_p2 + cand2.M() * cand2.M());

    fMomentumBeforeDecay = std::sqrt(energy_cand1 * energy_cand1 + 2 * energy_cand1 * energy_cand2 + energy_cand2 * energy_cand2 - fNeutralCandMass * fNeutralCandMass);

    TVector3 primaryVertex = fPrimaryVertex;
    TVector3 decayVertex = fDecayVertex;

    if (fVerbose > 0)
    {
        std::cout << " ----------- KFitNeutralCandFinder::calculateNeutralMotherCand() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    TVector3 vecPrimToDecayVertex = decayVertex - primaryVertex;

    TVector3 geom_dir_Z, vtx_geom_dir, geom_base_Z, vtx_geom_base;

    vtx_geom_base.SetX(vecPrimToDecayVertex.X());
    vtx_geom_base.SetY(vecPrimToDecayVertex.Y());
    vtx_geom_base.SetZ(vecPrimToDecayVertex.Z());
/*
    vtx_geom_dir.SetX(std::sin(TMath::RadToDeg() * vecPrimToDecayVertex.Theta()) * std::cos(TMath::RadToDeg() * vecPrimToDecayVertex.Phi()));
    vtx_geom_dir.SetY(std::sin(TMath::RadToDeg() * vecPrimToDecayVertex.Theta()) * std::sin(TMath::RadToDeg() * vecPrimToDecayVertex.Phi()));
    vtx_geom_dir.SetZ(std::cos(TMath::RadToDeg() * vecPrimToDecayVertex.Theta()));
*/
    vtx_geom_dir.SetX(std::sin(vecPrimToDecayVertex.Theta()) * std::cos(vecPrimToDecayVertex.Phi()));
    vtx_geom_dir.SetY(std::sin(vecPrimToDecayVertex.Theta()) * std::sin(vecPrimToDecayVertex.Phi()));
    vtx_geom_dir.SetZ(std::cos(vecPrimToDecayVertex.Theta()));

    geom_dir_Z.SetX(0);
    geom_dir_Z.SetY(0);
    geom_dir_Z.SetZ(1);

    geom_base_Z.SetX(0);
    geom_base_Z.SetY(0);
    geom_base_Z.SetZ(1);

    Double_t a = geom_base_Z.Dot(geom_dir_Z);
    Double_t b = geom_dir_Z.Dot(geom_dir_Z);
    Double_t c = vtx_geom_base.Dot(geom_dir_Z);
    Double_t d = (geom_base_Z.Dot(vtx_geom_dir)) * (vtx_geom_dir.Dot(geom_dir_Z)) / vtx_geom_dir.Dot(vtx_geom_dir);
    Double_t e = (geom_dir_Z.Dot(vtx_geom_dir)) * (vtx_geom_dir.Dot(geom_dir_Z)) / vtx_geom_dir.Dot(vtx_geom_dir);
    Double_t f = (vtx_geom_base.Dot(vtx_geom_dir)) * (vtx_geom_dir.Dot(geom_dir_Z)) / vtx_geom_dir.Dot(vtx_geom_dir);

    Double_t u1 = (-a + c + d - f) / (b - e);

    Double_t valZ = geom_base_Z.Z() + geom_dir_Z.Z() * u1;

    TVector3 cross = vtx_geom_dir.Cross(geom_dir_Z);
    TVector3 diff=vtx_geom_base-geom_base_Z;

    double valR = abs(diff.Dot(cross)/(cross.Mag()));

    if (primaryVertex.Y() < 0)
    {
        valR = -1 * valR;
    }

    Double_t thetaPrimaryToSecondaryVertex, phiPrimaryToSecondaryVertex;

    thetaPrimaryToSecondaryVertex = vecPrimToDecayVertex.Theta();
    phiPrimaryToSecondaryVertex = vecPrimToDecayVertex.Phi();

    fNeutralMotherCandidate.setR(valR);
    fNeutralMotherCandidate.setZ(valZ);

    Double_t Px = (fMomentumBeforeDecay *
                   std::sin(thetaPrimaryToSecondaryVertex) *
                   std::cos(phiPrimaryToSecondaryVertex));
    Double_t Py = (fMomentumBeforeDecay *
                   std::sin(thetaPrimaryToSecondaryVertex) *
                   std::sin(phiPrimaryToSecondaryVertex));
    Double_t Pz =
        (fMomentumBeforeDecay * std::cos(thetaPrimaryToSecondaryVertex));
    Double_t M = fNeutralCandMass;

    fNeutralMotherCandidate.SetXYZM(Px, Py, Pz, M);

    if (fVerbose > 0)
    {
        std::cout << "calculateNeutralMotherCandidate, fNeutralMotherCandidate: theta= " << fNeutralMotherCandidate.Theta() << " and phi = " << fNeutralMotherCandidate.Phi() << std::endl;
    }

    Double_t x_vertex = vecPrimToDecayVertex.X();
    Double_t y_vertex = vecPrimToDecayVertex.Y();
    Double_t z_vertex = vecPrimToDecayVertex.Z();

    double sigma_x = sqrt(fPrimVtxResX * fPrimVtxResX + fDecVtxResX * fDecVtxResX);
    double sigma_y = sqrt(fPrimVtxResY * fPrimVtxResY + fDecVtxResY * fDecVtxResY);
    double sigma_z = sqrt(fPrimVtxResZ * fPrimVtxResZ + fDecVtxResZ * fDecVtxResZ);

    // Use coordinate transformation cartesian->polar to estimate error in theta and phi

    // Calculate the error in theta
    Double_t r = std::sqrt(x_vertex * x_vertex + y_vertex * y_vertex + z_vertex * z_vertex);

    Double_t dr_dx = x_vertex / r;
    Double_t dr_dy = y_vertex / r;
    Double_t dr_dz = z_vertex / r;

    Double_t dtheta_dx = x_vertex * z_vertex / (r * r * r * std::sqrt(abs(1 - z_vertex / (r * r))));
    Double_t dtheta_dy = y_vertex * z_vertex / (r * r * r * std::sqrt(abs(1 - z_vertex / (r * r))));
    Double_t dtheta_dz = (1 / r - z_vertex * z_vertex / (r * r * r)) / std::sqrt(abs(1 - z_vertex * z_vertex / (r * r)));

    // Double_t sigma_theta =std::sqrt(dtheta_dx * dtheta_dx * sigma_x * sigma_x + dtheta_dy * dtheta_dy * sigma_y * sigma_y + dtheta_dz * dtheta_dz * sigma_z * sigma_z);

    // Calculate the error in phi
    Double_t r_2D = std::sqrt(x_vertex * x_vertex + y_vertex * y_vertex);

    Double_t dphi_dx = -x_vertex * y_vertex / (sqrt(x_vertex * x_vertex / (r_2D * r_2D)) * r_2D * r_2D * r_2D);
    Double_t dphi_dy = std::sqrt(x_vertex * x_vertex / (r_2D * r_2D)) / r_2D;

    Double_t dphi_dz = 0;

    // Double_t sigma_phi =std::sqrt(dphi_dx * dphi_dx * sigma_x * sigma_x + dphi_dy * dphi_dy * sigma_y * sigma_y);

    // Calculate the error in R
    Double_t dR_dx = x_vertex / r_2D;
    Double_t dR_dy = y_vertex / r_2D;

    Double_t sigma_R = std::sqrt(dR_dx * dR_dx * sigma_x * sigma_x + dR_dy * dR_dy * sigma_y * sigma_y);

    TMatrixD covarianceDetectorSystem;
    covarianceDetectorSystem.ResizeTo(3, 3);
    covarianceDetectorSystem.Zero();

    covarianceDetectorSystem(0, 0) = sigma_x * sigma_x;
    covarianceDetectorSystem(1, 1) = sigma_y * sigma_y;
    covarianceDetectorSystem(2, 2) = sigma_z * sigma_z;

    TMatrixD derivativeMatrix;
    derivativeMatrix.ResizeTo(3, 3);
    derivativeMatrix.Zero();

    // First row derivatives of r
    derivativeMatrix(0, 0) = dr_dx;
    derivativeMatrix(0, 1) = dr_dy;
    derivativeMatrix(0, 2) = dr_dz;

    // Second row derivative of theta
    derivativeMatrix(1, 0) = dtheta_dx;
    derivativeMatrix(1, 1) = dtheta_dy;
    derivativeMatrix(1, 2) = dtheta_dz;

    // Last row, derviatives of phi
    derivativeMatrix(2, 0) = dphi_dx;
    derivativeMatrix(2, 1) = dphi_dy;
    derivativeMatrix(2, 2) = dphi_dz;

    TMatrixD derivativeMatrixTranspose;
    derivativeMatrixTranspose.ResizeTo(3, 3);
    derivativeMatrixTranspose.Zero();
    derivativeMatrixTranspose.Transpose(derivativeMatrix);

    TMatrixD covarianceNeutralMother;
    covarianceNeutralMother.ResizeTo(3, 3);
    covarianceNeutralMother = derivativeMatrix * covarianceDetectorSystem * derivativeMatrixTranspose;

    fCovarianceNeutralMother.ResizeTo(5, 5);
    fCovarianceNeutralMother(0, 0) = 9999999;
    fCovarianceNeutralMother(1, 1) = covarianceNeutralMother(1, 1);
    fCovarianceNeutralMother(2, 2) = covarianceNeutralMother(2, 2);
    fCovarianceNeutralMother(1, 2) = covarianceNeutralMother(1, 2);
    fCovarianceNeutralMother(2, 1) = covarianceNeutralMother(2, 1);
    fCovarianceNeutralMother(3, 3) = sigma_R * sigma_R;
    fCovarianceNeutralMother(4, 4) = sigma_z * sigma_z;

    fNeutralMotherCandidate.setCovariance(fCovarianceNeutralMother);

    fNeutralMotherCandidate.setMomentum(fMomentumBeforeDecay);

    // Set angles

    fNeutralMotherCandidate.setTheta(thetaPrimaryToSecondaryVertex);
    fNeutralMotherCandidate.setPhi(phiPrimaryToSecondaryVertex);
}
