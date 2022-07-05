
#include "hcovariancekinfit.h"

//HCovarianceKinFit* HCovarianceKinFit::gCovariance = NULL;

//ClassImp(HCovarianceKinFit) //?
HCovarianceKinFit::HCovarianceKinFit() : fMomDepErrors(false), 
                                         fVerbose(0)

{
    if (fVerbose > 0)
    {
        std::cout << "--------------- HCovarianceKinFit() -----------------" << std::endl;
    }
}


void HCovarianceKinFit::estimateCov(Int_t pid, Double_t mom, double (&covariance)[5])
{
    if (fVerbose > 0)
    {
        std::cout << " ----------- HDecayBuildeer::estimateCovarianceMatrix() -----------" << std::endl;
        std::cout << "" << std::endl;
    }

    if (fMomDepErrors == true && fSetup == "pp45")
    {
        //Momentum dependent uncertainty estimation input for HADES particle candidates
        TFile *momErr = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepMomErrors_allParticles_recoMom_errFunc.root", "read");
        TF1 *momErrP = (TF1 *)momErr->Get("f_pP");
        TF1 *momErrPi = (TF1 *)momErr->Get("f_pPi");
        TF1 *momErrK = (TF1 *)momErr->Get("f_pK");
        TF1 *momErrEm = (TF1 *)momErr->Get("f_pEm");

        TFile *thtErr = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepThtErrors_recoMom_errFunc.root", "read");
        TF1 *thtErrP = (TF1 *)thtErr->Get("f_thtP");
        TF1 *thtErrPi = (TF1 *)thtErr->Get("f_thtPi");
        TF1 *thtErrK = (TF1 *)thtErr->Get("f_thtK");
        TF1 *thtErrEm = (TF1 *)thtErr->Get("f_thtEm");

        TFile *phiErr = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepPhiErrors_recoMom_errFunc.root", "read");
        TF1 *phiErrP = (TF1 *)phiErr->Get("f_phiP");
        TF1 *phiErrPi = (TF1 *)phiErr->Get("f_phiPi");
        TF1 *phiErrK = (TF1 *)phiErr->Get("f_phiK");
        TF1 *phiErrEm = (TF1 *)phiErr->Get("f_phiEm");

        TFile *RZErr_PPi = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepRZErrors_PPi_errFunc.root", "read");
        TFile *RZErr_PK = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepRZErrors_PK_errFunc.root", "read");
        TFile *RZErr_PEm = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepRZErrors_recoMom_PEm_errFunc.root", "read");
        TF1 *RErrP = (TF1 *)RZErr_PEm->Get("f_RP");
        TF1 *RErrPi = (TF1 *)RZErr_PPi->Get("f_RPi");
        TF1 *RErrK = (TF1 *)RZErr_PK->Get("f_RK");
        TF1 *RErrEm = (TF1 *)RZErr_PEm->Get("f_REm");
        TF1 *ZErrP = (TF1 *)RZErr_PEm->Get("f_ZP");
        TF1 *ZErrPi = (TF1 *)RZErr_PPi->Get("f_ZPi");
        TF1 *ZErrK = (TF1 *)RZErr_PK->Get("f_ZK");
        TF1 *ZErrEm = (TF1 *)RZErr_PEm->Get("f_ZEm");


        if (pid == 14) {
			double cov[] = {momErrP->Eval(mom), thtErrP->Eval(mom), phiErrP->Eval(mom),
                               RErrP->Eval(mom), ZErrP->Eval(mom)}; 
            std::copy(cov, cov+5, covariance);
		} else if (pid == 9) {
			double cov[] = {momErrPi->Eval(mom), thtErrPi->Eval(mom), phiErrPi->Eval(mom),
                               RErrPi->Eval(mom), ZErrPi->Eval(mom)}; 
            std::copy(cov, cov+5, covariance);
        } else if (pid == 11) {
			double cov[] = {momErrK->Eval(mom), thtErrK->Eval(mom), phiErrK->Eval(mom),
                               RErrK->Eval(mom), ZErrK->Eval(mom)}; 
            std::copy(cov, cov+5, covariance);
        } else if (pid == 3) {
			double cov[] = {momErrEm->Eval(mom), thtErrEm->Eval(mom), phiErrEm->Eval(mom),
                               RErrEm->Eval(mom), ZErrEm->Eval(mom)};
            std::copy(cov, cov+5, covariance);
        } else cout<<"No momentum dependent error estimate available for this pid"<<endl;
    }

    else if(fMomDepErrors == true && fSetup=="FwDet")
    {
        // FwDet protons
        TFile *momErr_fw = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepMomErrors_fw_recoMom_errFunc.root", "read");
        TF1 *momErrP_fw = (TF1 *)momErr_fw->Get("f_pP");
        TFile *thtErr_fw = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepThtErrors_fw_recoMom_errFunc.root", "read");
        TF1 *thtErrP_fw = (TF1 *)thtErr_fw->Get("f_thtP");
        TFile *phiErr_fw = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepPhiErrors_fw_recoMom_errFunc.root", "read");
        TF1 *phiErrP_fw = (TF1 *)phiErr_fw->Get("f_phiP");
        TFile *RZErr_fw = new TFile("/lustre/hades/user/jrieger/pp_pKSigma0_dalitz/sim/ana/pp_pKSigma0_dalitz_p4500p_momDepRZErrors_fw_recoMom_errFunc.root", "read");
        TF1 *RErrP_fw = (TF1 *)RZErr_fw->Get("f_RP");
        TF1 *ZErrP_fw = (TF1 *)RZErr_fw->Get("f_ZP");

        double cov[] = {momErrP_fw->Eval(mom), thtErrP_fw->Eval(mom), phiErrP_fw->Eval(mom),
                                RErrP_fw->Eval(mom), ZErrP_fw->Eval(mom)};
        std::copy(cov, cov+5, covariance);
    }

    else if (fMomDepErrors == false && fSetup == "pp35") //which momentum/setup is this for?
    {
        if (pid == 14) {
			double cov[] = {1.469 * 1e-5, 2.410 * 1e-3, 5.895 * 1e-3,
                               1.188, 2.652};
            std::copy(cov, cov+5, covariance);
        } else if (pid == 9) {
			double cov[] = {5.959 * 1e-5, 9.316 * 1e-3, 1.991 * 1e-2,
                               4.006, 7.629};
            std::copy(cov, cov+5, covariance);
        } else if (pid == 11) {
			double cov[] = {1.947 * 1e-5, 2.296 * 1e-3, 6.312 * 1e-3,
                               1.404, 2.723};
            std::copy(cov, cov+5, covariance);
        } else cout<<"No error estimate available for this pid"<<endl;
    }
    else {
		cout<<"No suitable covariance estimate found"<<endl;
		double cov[] = {0,0,0,0,0};
        std::copy(cov, cov+5, covariance);
	}
}
