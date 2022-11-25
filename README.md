KinFit: A library for few-GeV hadronic and leptonic decays

Installation:

1) Clone the repository
    git clone https://github.com/KinFit/KinFit.git
2) Build the library using cmake
    mkdir build
    cd build
    cmake ..
    make
3) Add library to rootlogon.C
    gSystem->Load("pathtobuild/libKinFit.so")



Direkt use of the fitter tools in your own analysis:

1) Create KFitParticles from your particle candidates
    KFitParticle particle(TLorentzVector lv, Double_t R, Double_t Z);
    FillData(particle, particle_covariance)

2) Create a vector of particles which shall be fitted
    std::vector<KFitParticle> cands;
    cands.push_back(particle);

3) Initialize the KinFitter and fit
    KinFitter fitter(cands);
    fitter.addXXConstraint(arguments);
    fitter.fit()

4) Get the fit result
    KFitParticle fcand = fitter.getDaughter(0);

The procedure is illustrated in the example macro fit_toyMC.C.
A more advanced example including the construction of a neutral weakly decaying particle is given in the macro fitLambda3C_toyMC_fromPluto.C.



Automated analysis:

The KFitRootAnalyzer can be used to automatically analyze a root file event by event, performing the selected fit. The best particle combination for each event according to fit probability is written to an output file.

Input: A root file containing a TClonesArray of KFitParticles

Output: A root file containing a TClonesArray of KFitParticles, fit probability and chi2

Usage:
1) User uses the macro analysis_user.C to define input and output and the fit that should be performed.

2) HRootAnalyzer object is created. User adds fitting tasks to it, defines PIDs of particles to be fitted. 

3)HRootAnalyzer creates the output tree. Event loop in wich particles are selected according to their PID. A DecayBuilder object is created if the minimum of necessary particles for this task is found.

4) Decaybuilder performs combinatorics, creates a Fitter object, adds the constraint and calls the fitter for each combination

5) Fitter performes the fit

6) Get the updated particle candidates and write to output tree

