//****************************************************************************
//*                    This file is part of KinFit.                          *
//*                                                                          *
//*            	KinFit is distributed under the terms of the                 *
//*              GNU General Public License (GPL) version 3,                 *
//*                 copied verbatim in the file "LICENSE".                   *
//*                                                                          *
//*  				           Copyright 2024                                *
//*		         GSI Helmholtzzentrum für Schwerionenforschung               *
//* 	      This software is distributed under the terms of the            *
//*	          GNU General Public Licence version 3 (GPL Version 3)           *
//*		      			     				                                 *
//*     The copyright holders are listed in the file "COPYRIGHTHOLDERS".     *
//*               The authors are listed in the file "AUTHORS".              *
//****************************************************************************

#include "KFitDecayBuilder.h"

KFitDecayBuilder::KFitDecayBuilder(TString task, std::vector<int> pids, TLorentzVector lv, double mass) : fTask(task),
                                                                                                            fPids(pids),
                                                                                                            fVerbose(0)

{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder() -----------------" << std::endl;
    }

    setIniSys(lv);
    setMass(mass);

}

void KFitDecayBuilder::countCombis()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::countCombis() -----------------" << std::endl;
    }
    fCombiCounter = 0;
    // Determine the number of combinations of the input particles
    fTotalCombis = 1;
    if(particleCounter.size()>0) particleCounter.clear();
    for (size_t i = 0; i < fPids.size(); i++)
    {
        fTotalCombis *= fCands[i].size();
        particleCounter.push_back(0);
    }

}

void KFitDecayBuilder::buildDecay()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::buildDecay() -----------------" << std::endl;
    }
    if (fVerbose > 1)
    {
        std::cout << "Combi counter: " << fCombiCounter << " total Combos: " << fTotalCombis << std::endl;
        std::cout<< "Task is " << fTask <<std::endl;
    }

    // For selection of best event according to probability
    fBestProb = 0;
    fBestChi2 = 1e6;

    while (fCombiCounter < (fTotalCombis ))
    {
        // Take one combination of particles
        fillFitCands();
        
        // Check if correct number of particles
        if (fFitCands.size() != fPids.size())
        {
            if (fVerbose > 1)
                std::cout << "wrong number of particles" << std::endl;
            fCombiCounter++;
            continue;
        }
        // Check for double usage of same particle
        if (doubleParticle)
        {
            if (fVerbose > 1)
                std::cout << "particle used twice" << std::endl;
            continue;
        }
        // Do task that was chosen by user
        doFit();

        if (fVerbose > 1)
        {
            std::cout << "Combi counter: " << fCombiCounter << " total Combos: " << fTotalCombis << std::endl;
            std::cout << "fit finished" <<std::endl;
        }
        
        
        
    }
    
}

Bool_t KFitDecayBuilder::doFit()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::doFit() -----------------" << std::endl;
    }
    KinFitter Fitter(fFitCands);
    if (fTask == "4C") Fitter.add4Constraint(fIniSys);
    else if (fTask == "Vertex") Fitter.addVertexConstraint();
    else if (fTask == "MassVtx") Fitter.addMassVtxConstraint(fMass);
    else if (fTask == "MM" || fTask == "MissingMass" ) Fitter.addMissingMassConstraint(fIniSys, fMass);
    else if (fTask == "Mom" || fTask == "MissingParticle") Fitter.addMissingParticleConstraint(fIniSys, fMass);
    else if (fTask == "Mass"){
        Fitter.addMassConstraint(fMass);
        /*
        cout << "constraint added, Mass = " << fMass << endl;
        cout << "proton momentum = " << fFitCands[0].getMomentum() << endl;
        cout << "pion momentum = " << fFitCands[1].getMomentum() << endl;
        TMatrixD cov_p = fFitCands[0].getCovariance();
        TMatrixD cov_pi = fFitCands[1].getCovariance();
        /*
        cout << "proton momentum error = " << cov_p(0,0) << endl;
        cout << "pion momentum error = " << cov_pi(0,0) << endl;
*/
    } 
    else cout << "Task not available" << endl;

    if (Fitter.fit())
    {
        if (fVerbose > 1)
        {
            std::cout << "fit successful" << std::endl;
        }
        if (Fitter.getProb() > fBestProb)
        {
            fBestProb = Fitter.getProb();
            fBestChi2 = Fitter.getChi2();
            if(fOutputCands.size() != 0) fOutputCands.clear();
            Fitter.getDaughters(fOutputCands);
        }
        return kTRUE;
    }
    else
    {
        if (fVerbose > 1)
        {
            std::cout << "fit not successful" << std::endl;
        }
        return kFALSE;
    }
}

void KFitDecayBuilder::fillFitCands()
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::fillFitCands() -----------------" << std::endl;
    }

    if(fFitCands.size()>0) fFitCands.clear();
    doubleParticle = false;
    
    for (size_t i = 0; i < fPids.size(); i++)
    {
        checkDoubleParticle(i);
        if (fVerbose > 1)
            std::cout << "particle counter " << i << " is " << particleCounter[i] << " now" << std::endl;
        fFitCands.push_back(fCands[i][particleCounter[i]]);
        if (fVerbose > 1)
            std::cout << "fill in pid " << fPids[i] << std::endl;
    }
    int a = fPids.size() - 1;
    if (fVerbose > 1)
        std::cout << "particle counter " << a << " is " << particleCounter[a] << " now" << std::endl;
    while (particleCounter[a] == (fCands[a].size() - 1))
    {
        particleCounter[a] = 0;
        if (fVerbose > 1)
            std::cout << "particle counter " << a << " is " << particleCounter[a] << " now" << std::endl;
        a--;
        if (a < 0)
        {
            if (!(fCombiCounter == fTotalCombis-1))
            if (fVerbose > 1)
                std::cout << "counted wrong: " << fCombiCounter << " != " << fTotalCombis << std::endl;
            break;

            a = 1e6;
        }
    }
    particleCounter[a]++;
    if (fVerbose > 1)
        std::cout << "particle counter " << a << " is " << particleCounter[a] << " now" << std::endl;
    fCombiCounter++;
    

    if (doubleParticle && (fCombiCounter < fTotalCombis))
        fillFitCands(); // If some particle has been filled more than once into fFitCands, repeat the procedure with the next combination
}

void KFitDecayBuilder::checkDoubleParticle(size_t i)
{
    if (fVerbose > 0)
    {
        std::cout << "--------------- KFitDecayBuilder::checkDoubleParticle() -----------------" << std::endl;
    }
    for (size_t j = 0; j < i; j++)
    {
        if ((fPids[j] == fPids[i]) && (particleCounter[j] == particleCounter[i]))
        {
            if (fVerbose > 1)
                std::cout << "double particle" << std::endl;
            doubleParticle = true;
        }
    }
    if (fVerbose > 1)
    {
        std::cout << "no double particle" << std::endl;
    }
}