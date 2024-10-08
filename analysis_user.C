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

#include "KFitAnalyzer.h"

Int_t analysis_user(TString infileList="../KinFit/QA/test_files/input.root", TString outfile = "fitted_test.root", Int_t nEvents=1000){

    // Initialize Analyzer 
    // Arguments are input file(list), output file name and number of events to analyze
    KFitAnalyzer RootAnalyzer(infileList, outfile, nEvents);

    // Create vector that contains the PIDs of the particles that should be analyzed
    std::vector<int> pids;
    pids.push_back(14); pids.push_back(11); pids.push_back(14); pids.push_back(9);    //4C fit
    //pids.push_back(14); pids.push_back(14); pids.push_back(9);    // Missing K fits
    //pids.push_back(14); pids.push_back(9);  // Vertex fit, Mass fit

    // Other input that might be needed, depending on constraint
    TLorentzVector ppSystem(0,0,5.3567,2*0.938272+4.500);
    Double_t mLambda = 1.11568;
    Double_t mK = 0.493677;

    RootAnalyzer.doFitterTask("4C", pids, -1, ppSystem);
    //RootAnalyzer.doFitterTask("Mass", pids, mLambda);  
    //RootAnalyzer.doFitterTask("MassVtx", pids, mLambda); 
    //RootAnalyzer.doFitterTask("MissingMass", pids, mK, ppSystem);  
    //RootAnalyzer.doFitterTask("MissingParticle", pids, mK, ppSystem); 
    //RootAnalyzer.doFitterTask("Vertex", pids, -1);  
    
    return 0;
}
