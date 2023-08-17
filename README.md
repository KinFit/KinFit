# KinFit: A Kinematic Fitting Package for Hadron Physics Experiments

**For a more detailed user's guide and description of the package see**

## Installation:

1) Clone the repository
    git clone https://github.com/KinFit/KinFit.git
2) Build the library using cmake   
    mkdir build   
    cd build   
    cmake .. -DCMAKE_INSTALL_PREFIX=/your_preferred_install_location   
    make   
    make install   
3) Add library to rootlogon.C   
    gSystem->Load("pathtoinstall/libKinFit.so")   

KinFit depends on the ROOT libraries Core, Physics and Tree   
When properly installed, KinFit can be included and linked to other CMake projects by   
   
find_package(KinFit REQUIRED)   
target_link_libraries(your_executable other_dependencies KinFit::KinFit ROOT::Core ROOT::Physics ROOT::Tree)   



## Direkt use of the fitter tools in your own analysis:   

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

The procedure is illustrated in the example macro fit_toyMC.C in QA. An example input file is provided in QA/test_files/toy_montecarlo_rz_vtx.root.   



## Automated analysis:   

The KFitAnalyzer can be used to automatically analyze a root file event by event, performing the selected fit. The best particle combination for each event according to fit probability is written to an output file.   

Input: A root file containing a TClonesArray of KFitParticles. Example: QA/test_files/input.root   

Output: A root file containing a TClonesArray of KFitParticles, fit probability and chi2   

### Usage:   
1) User uses the macro analysis_user.C to define input and output and the fit that should be performed.   

2) KFitAnalyzer object is created. User adds fitting tasks to it, defines PIDs of particles to be fitted.    

3) KFitAnalyzer creates the output tree. The particle combination with the best probability is written there for each event.   


