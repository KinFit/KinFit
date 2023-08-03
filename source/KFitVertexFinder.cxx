#include "KFitVertexFinder.h"

KFitVertexFinder::KFitVertexFinder(std::vector<KFitParticle> &cands) : fVerbose(0), fCands(cands)
{   
    
    fM.ResizeTo(3, 3);
    fSys.ResizeTo(3, 3);
    
    reset();

    for (int i_cand = 0; i_cand < fCands.size(); i_cand++)
    {
        KFitParticle cand = cands[i_cand];

        double param_theta = cand.Theta();
        double param_phi = cand.Phi();
        double param_R = cand.getR();
        double param_Z = cand.getZ();

        // Direction vector
        fDir.SetX(std::sin(param_theta) * std::cos(param_phi));
        fDir.SetY(std::sin(param_theta) * std::sin(param_phi));
        fDir.SetZ(std::cos(param_theta));

        // Base vector
        fBase.SetX(param_R * std::cos(param_phi + TMath::PiOver2()));
        fBase.SetY(param_R * std::sin(param_phi + TMath::PiOver2()));
        fBase.SetZ(param_Z);

        addLinesToVertex(fBase, fDir, 1.0); // Function for adding the lines to the vertex
    }

    findVertex(); //Function for finding the vertex and giving the output

}

void KFitVertexFinder::addLinesToVertex(const TVector3 &r, const TVector3 &alpha,
                             const double w)
{

    if (fVerbose > 0)
    {
        std::cout << " ----------- KFitVertexFinder::addLinesToVertex() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    fM(0, 0) = w * (alpha.Y() * alpha.Y() + alpha.Z() * alpha.Z());
    fM(0, 1) = -w * (alpha.X() * alpha.Y());
    fM(0, 2) = -w * (alpha.X() * alpha.Z());
    //  fM(1,0)=-(alpha.X() * alpha.Y());
    fM(1, 1) = w * (alpha.X() * alpha.X() + alpha.Z() * alpha.Z());
    fM(1, 2) = -w * (alpha.Y() * alpha.Z());
    // fM(2,0)=-(alpha.X() * alpha.Z());
    // fM(2,1)=-(alpha.Y() * alpha.Z());
    fM(2, 2) = w * (alpha.Y() * alpha.Y() + alpha.X() * alpha.X());

    fSys(0, 0) += fM(0, 0);
    fSys(0, 1) += fM(0, 1);
    fSys(0, 2) += fM(0, 2);
    fSys(1, 1) += fM(1, 1);
    fSys(1, 2) += fM(1, 2);
    fSys(2, 2) += fM(2, 2);

    fB.SetX(fB.X()+fM(0, 0) * r.X() + fM(0, 1) * r.Y() + fM(0, 2) * r.Z());
    fB.SetY(fB.Y()+fM(0, 1) * r.X() + fM(1, 1) * r.Y() + fM(1, 2) * r.Z());
    fB.SetZ(fB.Z()+fM(0, 2) * r.X() + fM(1, 2) * r.Y() + fM(2, 2) * r.Z());

}

void KFitVertexFinder::findVertex(){

    if (fVerbose > 0)
    {
        std::cout << " ----------- KFitVertexFinder::findVertex() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    double det = 0;

    //det = fSys.Determinant();
   
    det=fSys(0,0)*fSys(1,1)*fSys(2,2) + fSys(0,1)*fSys(1,2)*fSys(0,2)+ fSys(0,1)*fSys(1,2)*fSys(0,2) - fSys(0,2)*fSys(1,1)*fSys(0,2)-fSys(0,1)*fSys(0,1)*fSys(2,2) - fSys(1,2)*fSys(1,2)*fSys(0,0);

    fM(0, 0) = fSys(1, 1) * fSys(2, 2) - fSys(1, 2) * fSys(1, 2);
    fM(0, 1) = fSys(1, 2) * fSys(0, 2) - fSys(0, 1) * fSys(2, 2);
    fM(0, 2) = fSys(0, 1) * fSys(1, 2) - fSys(1, 1) * fSys(0, 2);
    fM(1, 1) = fSys(0, 0) * fSys(2, 2) - fSys(0, 2) * fSys(0, 2);
    fM(1, 2) = fSys(0, 1) * fSys(0, 2) - fSys(0, 0) * fSys(1, 2);
    fM(2, 2) = fSys(0, 0) * fSys(1, 1) - fSys(0, 1) * fSys(0, 1);
    fM(1, 0) = fM(0, 1);
    fM(2, 0) = fM(0, 2);
    fM(2, 1) = fM(1, 2);

    if (det == 0)
    {
        fVertex.SetXYZ(-1000., -1000., -1000.);

    }
    else
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
                fM(i, j) = fM(i, j) / (det);
        }
        fVertex = fM * fB;
    }
}

void KFitVertexFinder::reset()
{

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
            fM(i, j) = fSys(i, j) = 0.0;
    }
    fB.SetXYZ(0.0, 0.0, 0.0);
}
