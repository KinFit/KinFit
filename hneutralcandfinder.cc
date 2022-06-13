#include "hneutralcandfinder.h"

HNeutralCandFinder::HNeutralCandFinder(const std::vector<HRefitCand> &cands) : fCands(cands),  fVerbose(0), fMomentumAfterDecay(-1.), fNeutralCandMass(1115.683), fPrimaryVertexFound(false)
{
    if(fVerbose>0){
    std::cout << "--------------- HNeutralCandFinder -----------------" << std::endl;
    }
    double param_p1, param_p2;
    
    HRefitCand cand1 = cands[0];
    

    param_p1 = cand1.P(); // Not the inverse, this momentum is used for estimating the momentum of the Lambda Candidate

    HRefitCand cand2 = cands[1];

    param_p2 = cand2.P(); // Not the inverse, this momentum is used for estimating the momentum of the Lambda Candidate
    
    double energy_cand1, energy_cand2;
    energy_cand1 = sqrt(param_p1 * param_p1 + cand1.M() * cand1.M());
    energy_cand2 = sqrt(param_p2 * param_p2 + cand2.M() * cand2.M());

    fMomentumAfterDecay = sqrt(energy_cand1 * energy_cand1 + 2 * energy_cand1 * energy_cand2 + energy_cand2 * energy_cand2 - fNeutralCandMass*fNeutralCandMass);
    
}

void HNeutralCandFinder::setNeutralMotherCand(TVector3 primaryVertex, TVector3 decayVertex)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- HNeutralCandFinder::setNeutralMotherCand() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    TVector3 vecPrimToDecayVertex = decayVertex - primaryVertex;

    HGeomVector geom_dir_Z, vtx_geom_dir, geom_base_Z, vtx_geom_base;

    vtx_geom_base.setX(vecPrimToDecayVertex.X());
    vtx_geom_base.setY(vecPrimToDecayVertex.Y());
    vtx_geom_base.setZ(vecPrimToDecayVertex.Z());

    vtx_geom_dir.setX(std::sin(TMath::RadToDeg() *vecPrimToDecayVertex.Theta()) * std::cos(TMath::RadToDeg() *vecPrimToDecayVertex.Phi()));
    vtx_geom_dir.setY(std::sin(TMath::RadToDeg() *vecPrimToDecayVertex.Theta()) * std::sin(TMath::RadToDeg() *vecPrimToDecayVertex.Phi()));
    vtx_geom_dir.setZ(std::cos(TMath::RadToDeg() *vecPrimToDecayVertex.Theta()));

    geom_dir_Z.setX(0);
    geom_dir_Z.setY(0);
    geom_dir_Z.setZ(1);

    geom_base_Z.setX(0);
    geom_base_Z.setY(0);
    geom_base_Z.setZ(1);

    double a = geom_base_Z .scalarProduct(geom_dir_Z);
    double b = geom_dir_Z .scalarProduct(geom_dir_Z);
    double c = vtx_geom_base. scalarProduct(geom_dir_Z);
    double d = (geom_base_Z .scalarProduct(vtx_geom_dir)) * (vtx_geom_dir .scalarProduct(geom_dir_Z)) / vtx_geom_dir .scalarProduct(vtx_geom_dir);
    double e = (geom_dir_Z.scalarProduct(vtx_geom_dir)) * (vtx_geom_dir .scalarProduct(geom_dir_Z)) / vtx_geom_dir .scalarProduct(vtx_geom_dir);
    double f = (vtx_geom_base.scalarProduct(vtx_geom_dir)) * (vtx_geom_dir .scalarProduct(geom_dir_Z)) / vtx_geom_dir .scalarProduct(vtx_geom_dir);
    
    Double_t u1 = (-a + c + d - f) / (b - e); 

    double valZ = geom_base_Z.getZ() + geom_dir_Z.getZ()*u1;

    Double_t valR = HParticleTool::calculateMinimumDistance(vtx_geom_base,vtx_geom_dir,geom_base_Z,geom_dir_Z);
    
    if(primaryVertex.Y()<0){
        valR = -1 * valR;
    }

    double thetaPrimaryToSecondaryVertex, phiPrimaryToSecondaryVertex;

    thetaPrimaryToSecondaryVertex = vecPrimToDecayVertex.Theta();
    phiPrimaryToSecondaryVertex = vecPrimToDecayVertex.Phi();

    //fNeutralMotherCandidate.setTheta(TMath::RadToDeg() * thetaPrimaryToSecondaryVertex);
    //fNeutralMotherCandidate.setPhi(TMath::RadToDeg() * phiPrimaryToSecondaryVertex);
    //fNeutralMotherCandidate.SetTheta(thetaPrimaryToSecondaryVertex);
    //fNeutralMotherCandidate.SetPhi(phiPrimaryToSecondaryVertex);
    fNeutralMotherCandidate.setR(valR);
    fNeutralMotherCandidate.setZ(valZ); 

    //std::cout << "Neutral cand finder " << std::endl;
    //std::cout << "Lambda cand angles: " << fNeutralMotherCandidate.Theta() << " " << fNeutralMotherCandidate.Phi() << std::endl;
    
    //double Px = (fMomentumAfterDecay*
    //            std::sin(thetaPrimaryToSecondaryVertex*TMath::DegToRad()) *
    //            std::cos(phiPrimaryToSecondaryVertex*TMath::DegToRad()));
    //double Py = (fMomentumAfterDecay *
    //             std::sin(thetaPrimaryToSecondaryVertex*TMath::DegToRad()) *
    //            std::sin(phiPrimaryToSecondaryVertex*TMath::DegToRad()));
    //double Pz =
    //    (fMomentumAfterDecay * std::cos(thetaPrimaryToSecondaryVertex*TMath::DegToRad()));
    //double M = fNeutralCandMass;
    
    double Px = (fMomentumAfterDecay*
                std::sin(thetaPrimaryToSecondaryVertex) *
                std::cos(phiPrimaryToSecondaryVertex));
    double Py = (fMomentumAfterDecay *
                std::sin(thetaPrimaryToSecondaryVertex) *
                std::sin(phiPrimaryToSecondaryVertex));
    double Pz =
        (fMomentumAfterDecay * std::cos(thetaPrimaryToSecondaryVertex));
    double M = fNeutralCandMass;
    
    fNeutralMotherCandidate.SetXYZM(Px, Py, Pz, M);

    //std::cout << "Lambda cand angles 2: " << fNeutralMotherCandidate.Theta() << " " << fNeutralMotherCandidate.Phi() << std::endl;

    //fNeutralMotherCandidate.setMomentum(fMomentumAfterDecay);

    if (fVerbose > 0)
    {
        std::cout << "setNeutralMotherCandidate, fNeutralMotherCandidate: theta= " << fNeutralMotherCandidate.Theta() << " and phi = " << fNeutralMotherCandidate.Phi() << std::endl;
    }

    double x_vertex = vecPrimToDecayVertex.X();
    double y_vertex = vecPrimToDecayVertex.Y();
    double z_vertex = vecPrimToDecayVertex.Z();

    // the errors below are estimated from difference distributions between reconstructed - MC truth for the vertex
    // The errors are estimated from the histograms where both vertices were found in an event
    double sigma_x = sqrt(1.78590*1.78590+5.75369*5.75369); // when fwhm were used: sqrt(16.97*16.97+14.56*14.56);  // In mm
    double sigma_y = sqrt(1.75516*1.75516+5.57198*5.57198); // when fwhm were used: sqrt(16.80*16.80+14.59*14.59);  // In mm
    double sigma_z = sqrt(3.00431*3.00431+10.2602*10.2602); // when fwhm were used: sqrt(25.81*25.81+19.84*19.84);  // In mm

    // Use coordinate transformation cartesian->polar to estimate error in theta and phi

    // Calculate the error in theta
    double r = sqrt(x_vertex * x_vertex + y_vertex * y_vertex + z_vertex * z_vertex);

    double dr_dx=x_vertex/r;
    double dr_dy=y_vertex/r;
    double dr_dz=z_vertex/r;

    double dtheta_dx = x_vertex * z_vertex / (r * r * r * sqrt(1 - z_vertex / (r * r)));
    double dtheta_dy = y_vertex * z_vertex / (r * r * r * sqrt(1 - z_vertex / (r * r)));
    double dtheta_dz = (1 / r - z_vertex * z_vertex / (r * r * r)) / sqrt(1 - z_vertex * z_vertex / (r * r));

    //double sigma_theta = sqrt(dtheta_dx * dtheta_dx * sigma_x * sigma_x + dtheta_dy * dtheta_dy * sigma_y * sigma_y + dtheta_dz * dtheta_dz * sigma_z * sigma_z);

    // Calculate the error in phi
    double r_2D = sqrt(x_vertex * x_vertex + y_vertex * y_vertex);

    double dphi_dx = -x_vertex * y_vertex / (sqrt(x_vertex * x_vertex / (r_2D * r_2D)) * r_2D * r_2D * r_2D);
    double dphi_dy = sqrt(x_vertex * x_vertex / (r_2D * r_2D)) / r_2D;

    double dphi_dz=0;

    //double sigma_phi = sqrt(dphi_dx * dphi_dx * sigma_x * sigma_x + dphi_dy * dphi_dy * sigma_y * sigma_y);

    // Calculate the error in R
    double dR_dx = x_vertex / r_2D;
    double dR_dy = y_vertex / r_2D;

    double sigma_R = sqrt(dR_dx * dR_dx * sigma_x * sigma_x + dR_dy * dR_dy * sigma_y * sigma_y);
    
    TMatrixD covarianceDetectorSystem;
    covarianceDetectorSystem.ResizeTo(3, 3);
    covarianceDetectorSystem.Zero();

    covarianceDetectorSystem(0,0)= sigma_x * sigma_x;
    covarianceDetectorSystem(1,1)= sigma_y * sigma_y;
    covarianceDetectorSystem(2,2)= sigma_z * sigma_z;
    
    TMatrixD derivativeMatrix;
    derivativeMatrix.ResizeTo(3, 3);
    derivativeMatrix.Zero();
    
    // First row derivatives of r
    derivativeMatrix(0,0)=dr_dx;
    derivativeMatrix(0,1)=dr_dy;
    derivativeMatrix(0,2)=dr_dz;

    // Second row derivative of theta
    derivativeMatrix(1,0)=dtheta_dx;
    derivativeMatrix(1,1)=dtheta_dy;
    derivativeMatrix(1,2)=dtheta_dz;

    // Last row, derviatives of phi
    derivativeMatrix(2,0)=dphi_dx;
    derivativeMatrix(2,1)=dphi_dy;
    derivativeMatrix(2,2)=dphi_dz;

    TMatrixD derivativeMatrixTranspose;
    derivativeMatrixTranspose.ResizeTo(3, 3);
    derivativeMatrixTranspose.Zero();
    derivativeMatrixTranspose.Transpose(derivativeMatrix);

    TMatrixD covarianceNeutralMother;
    covarianceNeutralMother.ResizeTo(3, 3);
    covarianceNeutralMother=derivativeMatrix*covarianceDetectorSystem*derivativeMatrixTranspose;

    fCovarianceNeutralMother.ResizeTo(5, 5);
    fCovarianceNeutralMother(0, 0) = 9999999;
    //fCovarianceNeutralMother(1, 1) = sigma_theta * sigma_theta;
    //fCovarianceNeutralMother(2, 2) = sigma_phi * sigma_phi;
    fCovarianceNeutralMother(1, 1) = covarianceNeutralMother(1,1);
    fCovarianceNeutralMother(2, 2) = covarianceNeutralMother(2,2);
    fCovarianceNeutralMother(1, 2) = covarianceNeutralMother(1,2);
    fCovarianceNeutralMother(2, 1) = covarianceNeutralMother(2,1);
    fCovarianceNeutralMother(3, 3) = sigma_R * sigma_R;
    fCovarianceNeutralMother(4, 4) = sigma_z * sigma_z;

    fNeutralMotherCandidate.setCovariance(fCovarianceNeutralMother);
    
    //TMatrixD covariance = fNeutralMotherCandidate.getCovariance();
    //std::cout << "Neutral cand finder " << std::endl;
    //std::cout << "Lambda covariances: " << covariance(0, 0) << " "<< covariance(1, 1) << " " << covariance(2, 2) << " " << covariance(3, 3) << " " << covariance(4, 4) <<  std::endl;
    
    // Jenny: Comments below are for testing so that the covariance matrix is read in correctly
    // std::cout << "Nautral Cand Finder" << std::endl;
    // std::cout << "Diag elements: " << fCovarianceNeutralMother(1, 1) << " " << fCovarianceNeutralMother(2, 2) << std::endl;
    // std::cout << "Off diag elements: " << fCovarianceNeutralMother(1,2) << " " <<  fCovarianceNeutralMother(2,1) << std::endl;

}
