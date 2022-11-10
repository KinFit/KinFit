#include "hvertextools.h"

HVertexTools::HVertexTools() : fVerbose(0), fPrimaryVertexFound(false), fUsePrimaryVertexInNeutralCandidateCalculation(false)
{
}

TVector3 HVertexTools::findVertex(const std::vector<HRefitCand> &cands)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexTools::findVertex() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    Double_t param_theta1, param_phi1, param_R1, param_Z1;
    Double_t param_theta2, param_phi2, param_R2, param_Z2;

    if (cands.size() < 2)
    {

        std::cout << "WARNING: No vertex can be found, not enough charged particles in the event" << std::endl;
        fVertex.SetXYZ(-2000., -2000., -2000.);
        return fVertex;
    }

    HRefitCand cand1 = cands[0];

    param_theta1 = cand1.Theta();
    param_phi1 = cand1.Phi();
    param_R1 = cand1.getR();
    param_Z1 = cand1.getZ();

    HRefitCand cand2 = cands[1];

    param_theta2 = cand2.Theta();
    param_phi2 = cand2.Phi();
    param_R2 = cand2.getR();
    param_Z2 = cand2.getZ();

    // Calculate the base and direction vectors of the two candidates
    TVector3 vtx_base_1, vtx_base_2, vtx_dir_1, vtx_dir_2;

    // Base vectors
    vtx_base_1.SetXYZ(param_R1 * std::cos(param_phi1 + TMath::PiOver2()),
                      param_R1 * std::sin(param_phi1 + TMath::PiOver2()),
                      param_Z1);

    vtx_base_2.SetXYZ(param_R2 * std::cos(param_phi2 + TMath::PiOver2()),
                      param_R2 * std::sin(param_phi2 + TMath::PiOver2()),
                      param_Z2);

    // Direction vectors
    vtx_dir_1.SetXYZ(std::sin(param_theta1) * std::cos(param_phi1),
                     std::sin(param_theta1) * std::sin(param_phi1),
                     std::cos(param_theta1));

    vtx_dir_2.SetXYZ(std::sin(param_theta2) * std::cos(param_phi2),
                     std::sin(param_theta2) * std::sin(param_phi2),
                     std::cos(param_theta2));

    // Calculate the distance between the two tracks
    // Double_t distance = std::fabs((vtx_dir_1.Cross(vtx_dir_2)).Dot((vtx_base_1 - vtx_base_2)))/(vtx_dir_1.Cross(vtx_dir_2).Mag());

    TVector3 vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_base_1, vtx_geom_base_2;

    // Direction vectors
    vtx_geom_dir_1.SetX(std::sin(param_theta1) * std::cos(param_phi1));
    vtx_geom_dir_1.SetY(std::sin(param_theta1) * std::sin(param_phi1));
    vtx_geom_dir_1.SetZ(std::cos(param_theta1));
    vtx_geom_dir_2.SetX(std::sin(param_theta2) * std::cos(param_phi2));
    vtx_geom_dir_2.SetY(std::sin(param_theta2) * std::sin(param_phi2));
    vtx_geom_dir_2.SetZ(std::cos(param_theta2));

    // Base vectors
    vtx_geom_base_1.SetX(param_R1 * std::cos(param_phi1 + TMath::PiOver2()));
    vtx_geom_base_1.SetY(param_R1 * std::sin(param_phi1 + TMath::PiOver2()));
    vtx_geom_base_1.SetZ(param_Z1);
    vtx_geom_base_2.SetX(param_R2 * std::cos(param_phi2 + TMath::PiOver2()));
    vtx_geom_base_2.SetY(param_R2 * std::sin(param_phi2 + TMath::PiOver2()));
    vtx_geom_base_2.SetZ(param_Z2);

    if (fVerbose > 0)
    {
        std::cout << " " << std::endl;
        std::cout << " ---------------------------------------" << std::endl;
        std::cout << " " << std::endl;
        std::cout << "Printing base vector 1: " << std::endl;
        vtx_geom_base_1.print();

        std::cout << " " << std::endl;
        std::cout << "Printing direction vector 1: " << std::endl;
        vtx_geom_dir_1.print();
        std::cout << " " << std::endl;

        std::cout << " " << std::endl;
        std::cout << "Printing base vector 2: " << std::endl;
        vtx_geom_base_2.print();

        std::cout << " " << std::endl;
        std::cout << "Printing direction vector 2: " << std::endl;
        vtx_geom_dir_2.print();
        std::cout << " " << std::endl;
        std::cout << " ---------------------------------------" << std::endl;
        std::cout << " " << std::endl;
    }

    Double_t dist = HParticleTool::calculateMinimumDistance(vtx_geom_base_1, vtx_geom_dir_1, vtx_geom_base_2, vtx_geom_dir_2);

    fDistanceParticleToParticle = dist;

    // TVector3 vertex = HParticleTool::calculatePointOfClosestApproach(vtx_geom_base_1, vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_dir_2);

    TVector3 vertex;
    HGeomVertexFit *vtxFit = new HGeomVertexFit();
    vtxFit->addLine(vtx_geom_base_1, vtx_geom_dir_1, 1);
    vtxFit->addLine(vtx_geom_base_2, vtx_geom_dir_2, 1);
    vtxFit->getVertex(vertex);

    /* TVector3 cross = vtx_geom_dir_1.vectorProduct(vtx_geom_dir_2); // cross product: dir1 x dir2

    // straight lines are either skew or have a cross point

    TVector3 diff = vtx_geom_base_1;
    diff-=vtx_geom_base_2; // Difference of two base vectors base1 - base2

    Double_t D;
    D =  HParticleTool::calcDeterminant(vtx_geom_dir_2, vtx_geom_dir_1 ,cross);

    Double_t Dm =  HParticleTool::calcDeterminant(diff , vtx_geom_dir_1, cross);
    Double_t Dl = -HParticleTool::calcDeterminant(diff , vtx_geom_dir_2, cross);

    TVector3 vertex;
    TVector3 dm;
    TVector3 dl;

    dm = vtx_geom_dir_2;
    dm *= Dm;

    dl = vtx_geom_dir_1;
    dl *= Dl;

    vertex = dm - dl;

    vertex *= 2*((1.)/D); // Jenny added 2

    vertex+=vtx_geom_base_1;
    vertex+=vtx_geom_base_2;
    vertex*=0.5; */

        HGeomVector cross = dir1.vectorProduct(dir2); // cross product: dir1 x dir2

    // straight lines are either skew or have a cross point

    HGeomVector diff = base1;
    diff-=base2; // Difference of two base vectors base1 - base2

    Double_t D;
    D =  calcDeterminant(dir2, dir1 ,cross);

    if (!(fabs(D) > 0.))
    {
	::Warning(":calculatePointOfClosestApproach","Dirs and cross-product are lin. dependent: returning default Vertex (-20000,-20000,-20000)");

	return HGeomVector(-20000.,-20000.,-20000.);
    }

    Double_t Dm =  calcDeterminant(diff , dir1, cross);
    Double_t Dl = -calcDeterminant(diff , dir2, cross);

    HGeomVector vertex;
    HGeomVector dm;
    HGeomVector dl;

    dm = dir2;
    dm *= Dm;

    dl = dir1;
    dl *= Dl;

    vertex = dm - dl;

    vertex *= ((1.)/D);

    vertex+=base1;
    vertex+=base2;
    vertex*=0.5;

    fVertex.SetXYZ(vertex.X(), vertex.Y(), vertex.Z());

    Double_t distanceFromParticleToVertex_1 = HParticleTool::calculateMinimumDistanceStraightToPoint(vtx_geom_base_1, vtx_geom_dir_1, vertex);
    Double_t distanceFromParticleToVertex_2 = HParticleTool::calculateMinimumDistanceStraightToPoint(vtx_geom_base_2, vtx_geom_dir_2, vertex);

    TVector3 originVertex;
    originVertex.SetX(0.0);
    originVertex.SetY(0.0);
    originVertex.SetZ(0.0);

    Double_t distanceFromParticleToOrigin_1 = HParticleTool::calculateMinimumDistanceStraightToPoint(vtx_geom_base_1, vtx_geom_dir_1, originVertex);
    Double_t distanceFromParticleToOrigin_2 = HParticleTool::calculateMinimumDistanceStraightToPoint(vtx_geom_base_2, vtx_geom_dir_2, originVertex);

    fDistParticle1Vertex = distanceFromParticleToVertex_1;
    fDistParticle2Vertex = distanceFromParticleToVertex_2;
    fDistParticle1Origin = distanceFromParticleToOrigin_1;
    fDistParticle2Origin = distanceFromParticleToOrigin_2;

    // Find the primary vertex

    if (fVerbose > 0)
    {
        std::cout << "Vertex: theta: " << fVertex.Theta() << " and phi: " << fVertex.Phi() << std::endl;
    }

    fPrimaryVertexFound = false;

    if (cands.size() > 2)
    {
        findPrimaryVertex(cands);
    }

    if (fPrimaryVertexFound == true)
    {
        calculateVertexProperties(fPrimaryVertex, fVertex);
    }

    return fVertex;
}

std::vector<HRefitCand> HVertexTools::UpdateTrackParameters(std::vector<HRefitCand> &cands, TVector3 &vertexPos)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexTools::UpdateTrackParameters() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    Double_t param_theta1, param_phi1, param_R1, param_Z1;
    Double_t param_theta2, param_phi2, param_R2, param_Z2;

    HRefitCand cand1 = cands[0];

    // param_p_inv1 = 1. / cand1.P();
    param_theta1 = cand1.Theta();
    param_phi1 = cand1.Phi();
    param_R1 = cand1.getR();
    param_Z1 = cand1.getZ();

    HRefitCand cand2 = cands[1];

    // param_p_inv2 = 1. / cand2.P();
    param_theta2 = cand2.Theta();
    param_phi2 = cand2.Phi();
    param_R2 = cand2.getR();
    param_Z2 = cand2.getZ();

    // Calculate the base and direction vectors of the two candidates
    TVector3 vtx_base_1, vtx_base_2, vtx_dir_1, vtx_dir_2;

    // Base vectors
    vtx_base_1.SetXYZ(param_R1 * std::cos(param_phi1 + TMath::PiOver2()),
                      param_R1 * std::sin(param_phi1 + TMath::PiOver2()),
                      param_Z1);

    vtx_base_2.SetXYZ(param_R2 * std::cos(param_phi2 + TMath::PiOver2()),
                      param_R2 * std::sin(param_phi2 + TMath::PiOver2()),
                      param_Z2);

    // Direction vectors
    vtx_dir_1.SetXYZ(std::sin(param_theta1) * std::cos(param_theta1),
                     std::sin(param_theta1) * std::sin(param_phi1),
                     std::cos(param_theta1));

    vtx_dir_2.SetXYZ(std::sin(param_theta2) * std::cos(param_theta2),
                     std::sin(param_theta2) * std::sin(param_phi2),
                     std::cos(param_theta2));

    // Vectors pointing from vertex to POCA to Beam axis
    TVector3 vtx_base_1_updated, vtx_base_2_updated;
    vtx_base_1_updated = vtx_base_1 - vertexPos;
    vtx_base_2_updated = vtx_base_2 - vertexPos;

    if (fVerbose > 0)
    {
        std::cout << " " << std::endl;
        std::cout << "Position of base vector 1: x=" << vtx_base_1.X() << ", y=" << vtx_base_1.Y() << ", z=" << vtx_base_1.Z() << std::endl;
        std::cout << "Position of updated base vector 1: x=" << vtx_base_1_updated.X() << ", y=" << vtx_base_1_updated.Y() << ", z=" << vtx_base_1_updated.Z() << std::endl;
        std::cout << "Theta, original base vector 1: " << vtx_base_1.Theta() << ", Phi, original base vector: " << vtx_base_1.Phi() << std::endl;
        std::cout << "Theta, updated base vector 1: " << vtx_base_1_updated.Theta() << ", Phi, updated base vector: " << vtx_base_1_updated.Phi() << std::endl;
        std::cout << " " << std::endl;

        std::cout << " " << std::endl;
        std::cout << "Position of base vector 2: x=" << vtx_base_2.X() << ", y=" << vtx_base_2.Y() << ", z=" << vtx_base_2.Z() << std::endl;
        std::cout << "Position of updated base vector 2: x=" << vtx_base_2_updated.X() << ", y=" << vtx_base_2_updated.Y() << ", z=" << vtx_base_2_updated.Z() << std::endl;
        std::cout << "Theta, original base vector 2: " << vtx_base_2.Theta() << ", Phi, original base vector: " << vtx_base_2.Phi() << std::endl;
        std::cout << "Theta, updated base vector 2: " << vtx_base_2_updated.Theta() << ", Phi, updated base vector: " << vtx_base_2_updated.Phi() << std::endl;
        std::cout << " " << std::endl;
    }

    Double_t theta_secondary1 = vtx_base_1_updated.Theta();
    Double_t theta_secondary2 = vtx_base_2_updated.Theta();

    Double_t phi_secondary1 = vtx_base_1_updated.Phi();
    Double_t phi_secondary2 = vtx_base_2_updated.Phi();

    TVector3 vtx_dir_1_updated, vtx_dir_2_updated;

    vtx_dir_1_updated.SetXYZ(std::sin(theta_secondary1) * std::cos(phi_secondary1),
                             std::sin(theta_secondary1) * std::sin(phi_secondary1),
                             std::cos(theta_secondary1));
    vtx_dir_2_updated.SetXYZ(std::sin(theta_secondary2) * std::cos(phi_secondary2),
                             std::sin(theta_secondary2) * std::sin(phi_secondary2),
                             std::cos(theta_secondary2));

    // Calculate the distance between the two tracks
    Double_t dist_new = std::fabs((vtx_dir_1_updated.Cross(vtx_dir_2_updated)).Dot((vtx_base_1_updated - vtx_base_2_updated)));

    if (fVerbose > 0)
    {
        std::cout << "Minimum NEW distance between tracks: " << dist_new << std::endl;
        std::cout << " " << std::endl;
    }

    if (fVerbose > 0)
    {
        std::cout << "Before update " << std::endl;

        std::cout << "Cand1, theta: " << cand1.Theta() << std::endl;
        std::cout << "Cand1, phi: " << cand1.Phi() << std::endl;
        std::cout << "Cand2, theta: " << cand2.Theta() << std::endl;
        std::cout << "Cand2, phi: " << cand2.Phi() << std::endl;
    }

    cand1.SetTheta(theta_secondary1);
    cand1.SetPhi(phi_secondary1);

    cand2.SetTheta(theta_secondary2);
    cand2.SetPhi(phi_secondary2);

    if (fVerbose > 0)
    {
        std::cout << "After update " << std::endl;

        std::cout << "Cand1, theta: " << cand1.Theta() << std::endl;
        std::cout << "Cand1, phi: " << cand1.Phi() << std::endl;
        std::cout << "Cand2, theta: " << cand2.Theta() << std::endl;
        std::cout << "Cand2, phi: " << cand2.Phi() << std::endl;
    }

    std::vector<HRefitCand> newCands;
    newCands.clear();
    newCands.push_back(cand1);
    newCands.push_back(cand2);

    return newCands;
}

TVector3 HVertexTools::findPrimaryVertex(const std::vector<HRefitCand> &cands)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexTools::findPrimaryVertex() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    Double_t param_theta1, param_phi1, param_R1, param_Z1;
    Double_t param_theta2, param_phi2, param_R2, param_Z2;

    // All found protons in the event
    HRefitCand primaryCand1 = cands[0];

    param_theta1 = primaryCand1.Theta();
    param_phi1 = primaryCand1.Phi();
    param_R1 = primaryCand1.getR();
    param_Z1 = primaryCand1.getZ();

    HRefitCand primaryCand2 = cands[2];
    param_theta2 = primaryCand2.Theta();
    param_phi2 = primaryCand2.Phi();
    param_R2 = primaryCand2.getR();
    param_Z2 = primaryCand2.getZ();

    // Calculate the base and direction vectors of the two candidates
    TVector3 vtx_base_1, vtx_base_2, vtx_dir_1, vtx_dir_2;

    // Base vectors
    vtx_base_1.SetXYZ(param_R1 * std::cos(param_phi1 + TMath::PiOver2()),
                      param_R1 * std::sin(param_phi1 + TMath::PiOver2()),
                      param_Z1);

    vtx_base_2.SetXYZ(param_R2 * std::cos(param_phi2 + TMath::PiOver2()),
                      param_R2 * std::sin(param_phi2 + TMath::PiOver2()),
                      param_Z2);

    // Direction vectors
    vtx_dir_1.SetXYZ(std::sin(param_theta1) * std::cos(param_theta1),
                     std::sin(param_theta1) * std::sin(param_phi1),
                     std::cos(param_theta1));

    vtx_dir_2.SetXYZ(std::sin(param_theta2) * std::cos(param_theta2),
                     std::sin(param_theta2) * std::sin(param_phi2),
                     std::cos(param_theta2));

    TVector3 vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_base_1, vtx_geom_base_2;

    // Direction vectors
    vtx_geom_dir_1.SetX(std::sin(param_theta1) * std::cos(param_theta1));
    vtx_geom_dir_1.SetY(std::sin(param_theta1) * std::sin(param_phi1));
    vtx_geom_dir_1.SetZ(std::cos(param_theta1));
    vtx_geom_dir_2.SetX(std::sin(param_theta2) * std::cos(param_theta2));
    vtx_geom_dir_2.SetY(std::sin(param_theta2) * std::sin(param_phi2));
    vtx_geom_dir_2.SetZ(std::cos(param_theta2));

    // Base vectors
    vtx_geom_base_1.SetX(param_R1 * std::cos(param_phi1 + TMath::PiOver2()));
    vtx_geom_base_1.SetY(param_R1 * std::sin(param_phi1 + TMath::PiOver2()));
    vtx_geom_base_1.SetZ(param_Z1);
    vtx_geom_base_2.SetX(param_R2 * std::cos(param_phi2 + TMath::PiOver2()));
    vtx_geom_base_2.SetY(param_R2 * std::sin(param_phi2 + TMath::PiOver2()));
    vtx_geom_base_2.SetZ(param_Z2);

    TVector3 primaryVertex = HParticleTool::calculatePointOfClosestApproach(vtx_geom_base_1, vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_dir_2);

    fPrimaryVertex.SetXYZ(primaryVertex.X(), primaryVertex.Y(), primaryVertex.Z());

    // if the primary vertex was not found, each coordinate x,y,z is set to -20000
    if (primaryVertex.X() != -2000)
    {
        fPrimaryVertexFound = true;
    }

    return fPrimaryVertex;
}

void HVertexTools::calculateVertexProperties(TVector3 primaryVertex, TVector3 decayVertex)
{

    fVecPrimToDecayVertex = decayVertex - primaryVertex;

    fDistPrimToDecayVertex = sqrt((decayVertex.X() - primaryVertex.X()) * (decayVertex.X() - primaryVertex.X()) + (decayVertex.Y() - primaryVertex.Y()) * (decayVertex.Y() - primaryVertex.Y()) + (decayVertex.Z() - primaryVertex.Z()) * (decayVertex.Z() - primaryVertex.Z()));

    if (decayVertex.Z() > primaryVertex.Z())
    {
        fPrimaryVertexIsBetforeDecayVertex = true;
    }
    else
    {
        fPrimaryVertexIsBetforeDecayVertex = false;
    }

    Double_t R_primaryVertex, R_decayVertex;

    R_primaryVertex = sqrt(primaryVertex.X() * primaryVertex.X() + primaryVertex.Y() * primaryVertex.Y());
    R_decayVertex = sqrt(decayVertex.X() * decayVertex.X() + decayVertex.Y() * decayVertex.Y());

    if (R_primaryVertex < R_decayVertex)
    {
        fPrimaryVertexIsInsideDecayVertex = true;
    }
    else
    {
        fPrimaryVertexIsInsideDecayVertex = false;
    }
}
