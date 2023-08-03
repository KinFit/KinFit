#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TLorentzVector.h"
#include "TH1F.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "KinFitter.h"
#include "KFitRootAnalyzer.h"
#include "/home/pecarG/Downloads/KinFit-main/include/CoordinateConversion.h"
#include "/home/pecarG/Downloads/KinFit-main/include/KFitParticle.h"

Int_t analysis_user (TString inFileName, TString outFileName, Int_t nEvents){
KFitRootAnalyzer RootAnalyzer (inFileName, outFileName, nEvents);
std::vector<int> pids;
pids.push_back(14);
pids.push_back(9);
Double_t mass = 1.11568;
RootAnalyzer.doFitterTask ("Mass", pids, mass);
return 0;
}

class tokenizer{
    public:
    string val1, val2, val3, val4, val5, val6;

void simple_tokenizer(string s)
{
    stringstream ss(s);
    while (ss >> val1 >> val2 >> val3 >> val4 >> val5 >> val6) {
        //cout << val1 << endl;
        //cout << val2 << endl;
        //cout << val3 << endl;
        //cout << val4 << endl;
        //cout << val5 << endl;
        //cout << val6 << endl;
        }

    }

};


using namespace std;

int Testing(){

CoordinateConversion testconvert;
double ParameterVector[]{0.55, 0.22, -0.58, 0.18, -0.43, 1.29};
double ParameterErrorsVector[]{0.46, 0.19, 1, 0.0002, 0.0006, 0.0017};

testconvert.setparametersCart(ParameterVector);
testconvert.setparametererrorsCart(ParameterErrorsVector);

testconvert.setvarsSpher();
testconvert.seterrorsSpher();

std::cout << "Spherical Vars" << endl;
std::cout << "--------------------------" << endl;
std::cout << testconvert.mom << "+/-" << testconvert.errmom << endl;
std::cout << testconvert.theta << "+/-" << testconvert.errtheta << endl;
std::cout << testconvert.phi << "+/-" << testconvert.errphi << endl;
std::cout << testconvert.R << "+/-" << testconvert.errR << endl;
std::cout << testconvert.Z << "+/-" << testconvert.errZ << endl;
std::cout << "--------------------------" << endl;

CoordinateConversion testconvert2;
double ParameterVector2[]{testconvert.mom, testconvert.theta, testconvert.phi, testconvert.R, testconvert.Z};
double ParameterErrorsVector2[]{testconvert.errmom, testconvert.errtheta, testconvert.errphi, testconvert.errR, testconvert.errZ};

testconvert2.setparametersSpher(ParameterVector2);
testconvert2.setparametererrorsSpher(ParameterErrorsVector2);

testconvert2.setvarsCart();
testconvert2.seterrorsCart();

std::cout << "Cartesian Vars" << endl;
std::cout << "--------------------------" << endl;
std::cout << testconvert2.Px << "+/-" << testconvert2.errPx << endl;
std::cout << testconvert2.Py << "+/-" << testconvert2.errPy << endl;
std::cout << testconvert2.Pz << "+/-" << testconvert2.errPz << endl;
std::cout << testconvert2.X << "+/-" << testconvert2.errX << endl;
std::cout << testconvert2.Y << "+/-" << testconvert2.errY << endl;
std::cout << testconvert2.Z << "+/-" << testconvert2.errZ << endl;

testconvert2.setvarsSpher();
testconvert2.seterrorsSpher();

std::cout << "Spherical Vars 2" << endl;
std::cout << "--------------------------" << endl;
std::cout << testconvert2.mom << "+/-" << testconvert2.errmom << endl;
std::cout << testconvert2.theta << "+/-" << testconvert2.errtheta << endl;
std::cout << testconvert2.phi << "+/-" << testconvert2.errphi << endl;
std::cout << testconvert2.R << "+/-" << testconvert2.errR << endl;
std::cout << testconvert2.Z << "+/-" << testconvert2.errZ << endl;
std::cout << "--------------------------" << endl;

std::vector<double> Px;
std::vector<double> Py;
std::vector<double> Pz;
std::vector<double> X;
std::vector<double> Y;
std::vector<double> Z;

tokenizer obj1;
string line;

ifstream filetest;
filetest.open("cartfiletest.txt");
if (filetest.is_open()){
    while (getline(filetest, line)){
        obj1.simple_tokenizer(line);
        Px.push_back(stod(obj1.val1));
        Py.push_back(stod(obj1.val2));
        Pz.push_back(stod(obj1.val3));
        X.push_back(stod(obj1.val4));
        Y.push_back(stod(obj1.val5));
        Z.push_back(stod(obj1.val6));

    }

//    for (double Pxout: Px) {
//    std::cout << Pxout << endl;}

}

//    analysis_user("cartfiletest2.txt", "outputfiletest.root", 499);




return(0);
}
