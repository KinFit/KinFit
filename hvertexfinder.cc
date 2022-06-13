#include "hvertexfinder.h"

HVertexFinder::HVertexFinder() : fVerbose(0)
{
}

TVector3 HVertexFinder::findVertex(const std::vector<HRefitCand> &cands)
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HVertexFinder::findVertex() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    double param_theta, param_phi, param_R, param_Z;
    // Calculate the base and direction vectors of the two candidates
    TVector3 vtx_base, vtx_dir;
    HGeomVector vtx_geom_dir, vtx_geom_base;
    

    if (cands.size() < 2)
    {

        std::cout << "WARNING: No vertex can be found, not enough charged particles in the event" << std::endl;
        fVertex.SetXYZ(-2000., -2000., -2000.);
        return fVertex;
    }

    HGeomVector vertex;
    HGeomVertexFit *vtxFit = new HGeomVertexFit();

    for (int i_cands = 0; i_cands < (int)cands.size(); i_cands++)
    {

        HRefitCand cand = cands[i_cands];

        param_theta = cand.Theta();
        param_phi = cand.Phi();
        param_R = cand.getR();
        param_Z = cand.getZ();

        // Base vectors
        vtx_base.SetXYZ(param_R * std::cos(param_phi + TMath::PiOver2()),
                          param_R * std::sin(param_phi + TMath::PiOver2()),
                          param_Z);

        // Direction vectors
        vtx_dir.SetXYZ(std::sin(param_theta) * std::cos(param_phi),
                         std::sin(param_theta) * std::sin(param_phi),
                         std::cos(param_theta));

        // Direction vectors
        vtx_geom_dir.setX(std::sin(param_theta) * std::cos(param_phi));
        vtx_geom_dir.setY(std::sin(param_theta) * std::sin(param_phi));
        vtx_geom_dir.setZ(std::cos(param_theta));

        // Base vectors
        vtx_geom_base.setX(param_R * std::cos(param_phi + TMath::PiOver2()));
        vtx_geom_base.setY(param_R * std::sin(param_phi + TMath::PiOver2()));
        vtx_geom_base.setZ(param_Z);

        vtxFit->addLine(vtx_geom_base, vtx_geom_dir, 1);
    }

    vtxFit->getVertex(vertex);

    fVertex.SetXYZ(vertex.X(), vertex.Y(), vertex.Z());

    return fVertex;
}

// std::vector<HRefitCand> HVertexFinder::UpdateTrackParameters(std::vector<HRefitCand> &cands, TVector3 &vertexPos)
// {
//     if (fVerbose > 0)
//     {
//         std::cout << " ----------- HVertexFinder::UpdateTrackParameters() -----------" << std::endl;
//         std::cout << "" << std::endl;
//     }

//     double param_theta1, param_phi1, param_R1, param_Z1;
//     double param_theta2, param_phi2, param_R2, param_Z2;

//     HRefitCand cand1 = cands[0];

//     //param_p_inv1 = 1. / cand1.P();
//     param_theta1 = cand1.Theta();
//     param_phi1 = cand1.Phi();
//     param_R1 = cand1.getR();
//     param_Z1 = cand1.getZ();

//     HRefitCand cand2 = cands[1];

//     //param_p_inv2 = 1. / cand2.P();
//     param_theta2 = cand2.Theta();
//     param_phi2 = cand2.Phi();
//     param_R2 = cand2.getR();
//     param_Z2 = cand2.getZ();

//     // Calculate the base and direction vectors of the two candidates
//     TVector3 vtx_base_1, vtx_base_2, vtx_dir_1, vtx_dir_2;

//     // Base vectors
//     vtx_base_1.SetXYZ(param_R1 * std::cos(param_phi1 + TMath::PiOver2()),
//                       param_R1 * std::sin(param_phi1 + TMath::PiOver2()),
//                       param_Z1);

//     vtx_base_2.SetXYZ(param_R2 * std::cos(param_phi2 + TMath::PiOver2()),
//                       param_R2 * std::sin(param_phi2 + TMath::PiOver2()),
//                       param_Z2);

//     // Direction vectors
//     vtx_dir_1.SetXYZ(std::sin(param_theta1) * std::cos(param_theta1),
//                      std::sin(param_theta1) * std::sin(param_phi1),
//                      std::cos(param_theta1));

//     vtx_dir_2.SetXYZ(std::sin(param_theta2) * std::cos(param_theta2),
//                      std::sin(param_theta2) * std::sin(param_phi2),
//                      std::cos(param_theta2));

//     //Vectors pointing from vertex to POCA to Beam axis
//     TVector3 vtx_base_1_updated, vtx_base_2_updated;
//     vtx_base_1_updated = vtx_base_1 - vertexPos;
//     vtx_base_2_updated = vtx_base_2 - vertexPos;

//     if (fVerbose > 0)
//     {
//         std::cout << " " << std::endl;
//         std::cout << "Position of base vector 1: x=" << vtx_base_1.X() << ", y=" << vtx_base_1.Y() << ", z=" << vtx_base_1.Z() << std::endl;
//         std::cout << "Position of updated base vector 1: x=" << vtx_base_1_updated.X() << ", y=" << vtx_base_1_updated.Y() << ", z=" << vtx_base_1_updated.Z() << std::endl;
//         std::cout << "Theta, original base vector 1: " << vtx_base_1.Theta() << ", Phi, original base vector: " << vtx_base_1.Phi() << std::endl;
//         std::cout << "Theta, updated base vector 1: " << vtx_base_1_updated.Theta() << ", Phi, updated base vector: " << vtx_base_1_updated.Phi() << std::endl;
//         std::cout << " " << std::endl;

//         std::cout << " " << std::endl;
//         std::cout << "Position of base vector 2: x=" << vtx_base_2.X() << ", y=" << vtx_base_2.Y() << ", z=" << vtx_base_2.Z() << std::endl;
//         std::cout << "Position of updated base vector 2: x=" << vtx_base_2_updated.X() << ", y=" << vtx_base_2_updated.Y() << ", z=" << vtx_base_2_updated.Z() << std::endl;
//         std::cout << "Theta, original base vector 2: " << vtx_base_2.Theta() << ", Phi, original base vector: " << vtx_base_2.Phi() << std::endl;
//         std::cout << "Theta, updated base vector 2: " << vtx_base_2_updated.Theta() << ", Phi, updated base vector: " << vtx_base_2_updated.Phi() << std::endl;
//         std::cout << " " << std::endl;
//     }

//     double theta_secondary1 = vtx_base_1_updated.Theta();
//     double theta_secondary2 = vtx_base_2_updated.Theta();

//     double phi_secondary1 = vtx_base_1_updated.Phi();
//     double phi_secondary2 = vtx_base_2_updated.Phi();

//     TVector3 vtx_dir_1_updated, vtx_dir_2_updated;

//     vtx_dir_1_updated.SetXYZ(std::sin(theta_secondary1) * std::cos(phi_secondary1),
//                              std::sin(theta_secondary1) * std::sin(phi_secondary1),
//                              std::cos(theta_secondary1));
//     vtx_dir_2_updated.SetXYZ(std::sin(theta_secondary2) * std::cos(phi_secondary2),
//                              std::sin(theta_secondary2) * std::sin(phi_secondary2),
//                              std::cos(theta_secondary2));

//     // Calculate the distance between the two tracks
//     double dist_new = std::fabs((vtx_dir_1_updated.Cross(vtx_dir_2_updated)).Dot((vtx_base_1_updated - vtx_base_2_updated)));

//     if (fVerbose > 0)
//     {
//         std::cout << "Minimum NEW distance between tracks: " << dist_new << std::endl;
//         std::cout << " " << std::endl;
//     }

//     if (fVerbose > 0)
//     {
//         std::cout << "Before update " << std::endl;

//         std::cout << "Cand1, theta: " << cand1.Theta() << std::endl;
//         std::cout << "Cand1, phi: " << cand1.Phi() << std::endl;
//         std::cout << "Cand2, theta: " << cand2.Theta() << std::endl;
//         std::cout << "Cand2, phi: " << cand2.Phi() << std::endl;
//     }

//     cand1.SetTheta(theta_secondary1);
//     cand1.SetPhi(phi_secondary1);

//     cand2.SetTheta(theta_secondary2);
//     cand2.SetPhi(phi_secondary2);

//     if (fVerbose > 0)
//     {
//         std::cout << "After update " << std::endl;

//         std::cout << "Cand1, theta: " << cand1.Theta() << std::endl;
//         std::cout << "Cand1, phi: " << cand1.Phi() << std::endl;
//         std::cout << "Cand2, theta: " << cand2.Theta() << std::endl;
//         std::cout << "Cand2, phi: " << cand2.Phi() << std::endl;
//     }

//     std::vector<HRefitCand> newCands;
//     newCands.clear();
//     newCands.push_back(cand1);
//     newCands.push_back(cand2);

//     return newCands;
// }

// TVector3 HVertexFinder::findPrimaryVertex(const std::vector<HRefitCand> &cands)
// {
//     if (fVerbose > 0)
//     {
//         std::cout << " ----------- HVertexFinder::findPrimaryVertex() -----------" << std::endl;
//         std::cout << "" << std::endl;
//     }

//     double param_theta1, param_phi1, param_R1, param_Z1;
//     double param_theta2, param_phi2, param_R2, param_Z2;

//     // All found protons in the event
//     HRefitCand primaryCand1 = cands[0];

//     param_theta1 = primaryCand1.Theta();
//     param_phi1 = primaryCand1.Phi();
//     param_R1 = primaryCand1.getR();
//     param_Z1 = primaryCand1.getZ();

//     HRefitCand primaryCand2 = cands[2];
//     param_theta2 = primaryCand2.Theta();
//     param_phi2 = primaryCand2.Phi();
//     param_R2 = primaryCand2.getR();
//     param_Z2 = primaryCand2.getZ();

//     // Calculate the base and direction vectors of the two candidates
//     TVector3 vtx_base_1, vtx_base_2, vtx_dir_1, vtx_dir_2;

//     // Base vectors
//     vtx_base_1.SetXYZ(param_R1 * std::cos(param_phi1 + TMath::PiOver2()),
//                       param_R1 * std::sin(param_phi1 + TMath::PiOver2()),
//                       param_Z1);

//     vtx_base_2.SetXYZ(param_R2 * std::cos(param_phi2 + TMath::PiOver2()),
//                       param_R2 * std::sin(param_phi2 + TMath::PiOver2()),
//                       param_Z2);

//     // Direction vectors
//     vtx_dir_1.SetXYZ(std::sin(param_theta1) * std::cos(param_theta1),
//                      std::sin(param_theta1) * std::sin(param_phi1),
//                      std::cos(param_theta1));

//     vtx_dir_2.SetXYZ(std::sin(param_theta2) * std::cos(param_theta2),
//                      std::sin(param_theta2) * std::sin(param_phi2),
//                      std::cos(param_theta2));

//     HGeomVector vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_base_1, vtx_geom_base_2;

//     // Direction vectors
//     vtx_geom_dir_1.setX(std::sin(param_theta1) * std::cos(param_theta1));
//     vtx_geom_dir_1.setY(std::sin(param_theta1) * std::sin(param_phi1));
//     vtx_geom_dir_1.setZ(std::cos(param_theta1));
//     vtx_geom_dir_2.setX(std::sin(param_theta2) * std::cos(param_theta2));
//     vtx_geom_dir_2.setY(std::sin(param_theta2) * std::sin(param_phi2));
//     vtx_geom_dir_2.setZ(std::cos(param_theta2));

//     // Base vectors
//     vtx_geom_base_1.setX(param_R1 * std::cos(param_phi1 + TMath::PiOver2()));
//     vtx_geom_base_1.setY(param_R1 * std::sin(param_phi1 + TMath::PiOver2()));
//     vtx_geom_base_1.setZ(param_Z1);
//     vtx_geom_base_2.setX(param_R2 * std::cos(param_phi2 + TMath::PiOver2()));
//     vtx_geom_base_2.setY(param_R2 * std::sin(param_phi2 + TMath::PiOver2()));
//     vtx_geom_base_2.setZ(param_Z2);

//     HGeomVector primaryVertex = HParticleTool::calculatePointOfClosestApproach(vtx_geom_base_1, vtx_geom_dir_1, vtx_geom_dir_2, vtx_geom_dir_2);

//     fPrimaryVertex.SetXYZ(primaryVertex.X(), primaryVertex.Y(), primaryVertex.Z());

//     // if the primary vertex was not found, each coordinate x,y,z is set to -20000
//     if (primaryVertex.X() != -2000)
//     {
//         fPrimaryVertexFound = true;
//     }

//     return fPrimaryVertex;
// }