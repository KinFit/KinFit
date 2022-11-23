KinFit: A library for few-GeV hadronic and leptonic decays

Input: A root file containing
- A TLorentzVector / p, theta and phi for every track
- R, Z and PID for every track
- Covariance matrix for every track
- Primary event vertex
To decide: Tree structure: Track parameters summarized in object (TClonesArray of KFitParticle) or flat tree structure

Output: Same structure as input
To decide: Several fitted candidates for same event or just best probability? Orgaization of fitted/unfitted tracks. Fit probablility/chi2 in addition?


Usage:
1) User uses the macro analysis_user.C to define input and output and the fit(s) that should be performed.

2) HRootFitter object is created. User adds fitting tasks to it, defines PIDs of particles to be fitted. Output tree is created. Event loop in wich particles are selected according to their PID. A DecayBuilder object is created if the minimum of necessary particles for this task is found.

3) Decaybuilder performs combinatorics, creates a Fitter object, adds the constraint and calls the fitter for each combination

4) Fitter performes the fit

5) Get the updated particle candidates and write to output tree


Needs to be added:
- Compilation with cmake, CMakeLists.txt
- Multiple fitting steps
- Test on toy MC
- Add more constraints
