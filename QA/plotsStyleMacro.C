#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>

#include "TString.h"
#include "TLatex.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TColor.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFitResult.h"

using std::cout;
using std::endl;

TCanvas* create_canvas(){
	
	TCanvas *can = new TCanvas();
	can->SetCanvasSize(800, 500);
	can->SetWindowSize(801, 501);
	can->SetTopMargin(1.);
	can->SetBottomMargin(0.15);
	can->SetRightMargin(0.03);
	can->SetLeftMargin(0.10);

	return can;
}

TCanvas* create_canvas3(){
	
	TCanvas *can = new TCanvas();
	can->SetCanvasSize(1200, 400);

	can->SetWindowSize(1201, 401);
	can->SetTopMargin(0.01);
	can->SetBottomMargin(0.5);
	can->SetRightMargin(0.01);
	can->SetLeftMargin(0.10);

	can->Divide(3,1,0.00005,0.01);
	

	return can;
}

void plotLayout(TH1F *h){

	gStyle->SetOptStat(0);
		// Fixing the X axis
	//h->SetTitle(" ");
	h->GetXaxis()->SetLabelSize(0.06);
	h->GetXaxis()->SetTitleSize(0.06);
	h->GetXaxis()->SetTitleOffset(1.1);
	h->GetXaxis()->CenterTitle();
	//h->GetXaxis()->SetRangeUser(-40,40);
	//h->GetXaxis()->SetTitle("Distance to wire[mm]");

	// Fixing the Y axis
	h->GetYaxis()->SetLabelSize(0.06);
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetTitleOffset(0.7);
	h->GetYaxis()->CenterTitle();
	if(h->GetMaximum()/1000 >1) h->GetYaxis()->SetMaxDigits(3);
	if(h->GetMaximum()/10000 >1) h->GetYaxis()->SetMaxDigits(4);
	
	// h1->GetYaxis()->SetTitle("Occurrence");
	// h1->GetYaxis()->SetRangeUser(0,1400);

	//h->SetStats(kFALSE);
	h->SetLineWidth(1);
	// h1->SetLineStyle(9);
	h->SetLineColor(kBlack);
	//gStyle->SetTitleSize(0.3);
	//gStyle->SetTitleW(0.5) ;//title width
	//gStyle->SetTitleH(0.1); //title height

}

void plotLayout(TH2F *h){

	gStyle->SetOptStat(0);
		// Fixing the X axis
	//h->SetTitle(" ");
	h->GetXaxis()->SetLabelSize(0.06);
	h->GetXaxis()->SetTitleSize(0.06);
	h->GetXaxis()->SetTitleOffset(1.1);
	h->GetXaxis()->CenterTitle();
	//h->GetXaxis()->SetRangeUser(-40,40);
	//h->GetXaxis()->SetTitle("Distance to wire[mm]");

	// Fixing the Y axis
	h->GetYaxis()->SetLabelSize(0.06);
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetTitleOffset(0.8);
	h->GetYaxis()->CenterTitle();
	if(h->GetYaxis()->GetLast()/1000 >1) h->GetYaxis()->SetMaxDigits(3);
	if(h->GetYaxis()->GetLast()/10000 >1) h->GetYaxis()->SetMaxDigits(4);

	// Fixing the Z axis
	h->GetZaxis()->SetLabelSize(0.06);
	h->GetZaxis()->SetTitleSize(0.06);
	h->GetZaxis()->SetTitleOffset(1.7);
	h->GetZaxis()->CenterTitle();
	
	// h1->GetYaxis()->SetTitle("Occurrence");
	// h1->GetYaxis()->SetRangeUser(0,1400);

	//h->SetStats(kFALSE);
	//h->SetLineWidth(1);
	// h1->SetLineStyle(9);
	//h->SetLineColor(kBlack);
	h->SetOption("COLZ");
	//gStyle->SetTitleW(0.5) ;//title width
	//gStyle->SetTitleH(0.1); //title height

}

void plotsStyleMacro(TH1F *h1, TString outtitle, TString drawopt = "")
{

	plotLayout(h1);

	TCanvas *can = create_canvas();
	can->cd();

	if(drawopt=="bar") h1->GetXaxis()->SetLabelSize(0.08)
;	h1->Draw(drawopt);
	
	TLatex latexText(h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetLast())*0.8, h1->GetMaximum()*0.85, "KinFit");
	latexText.SetTextFont(42);
	latexText.SetTextSize(0.08);
	latexText.SetTextColor(15);
	latexText.DrawClone();

	can->Write();
	can->SaveAs("../pics/"+outtitle+".png");
}

void plotsStyleMacro(TH2F *h1, TString outtitle)
{

	plotLayout(h1);

	TCanvas *can = create_canvas();
	can->cd();

	h1->Draw("");
	
	TLatex latexText(h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetLast())*0.6, h1->GetYaxis()->GetBinUpEdge(h1->GetYaxis()->GetLast())*0.75, "KinFit");
	latexText.SetTextFont(42);
	latexText.SetTextSize(0.08);
	latexText.SetTextColor(15);
	latexText.DrawClone();

	can->Write();
	can->SaveAs("../pics/"+outtitle+".png");
}

void plotsStyleMacro(TH1F *h1, TH1F *h2, TString outtitle, TString h1_title, TString h2_title)
{

	plotLayout(h1);
	h1->SetLineStyle(2);
	plotLayout(h2);

	TCanvas *can = create_canvas();
	can->cd();

	h2->Draw("");
	h1->Draw("same");
	h1->SetMaximum(max(h1->GetMaximum(), h2->GetMaximum()));
	h1->SetMaximum(h1->GetMaximum()*1.05);
	
	TLatex latexText(h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetLast())*0.65, h1->GetMaximum()*0.85, "KinFit");
	latexText.SetTextFont(42);
	latexText.SetTextSize(0.08);
	latexText.SetTextColor(15);
	latexText.DrawClone();

	auto legend = new TLegend(0.15,0.6,0.45,0.75);
	legend->AddEntry(h1, h1_title);
	legend->AddEntry(h2, h2_title);
	legend->Draw();

	can->Write();
	can->SaveAs("../pics/"+outtitle+".png");
}

void plotsStyleMacro_andFit(TH1F *h1, TString outtitle)
{

	plotLayout(h1);

	TCanvas *can = create_canvas();
	can->cd();

	h1->Draw("");

	TF1 *fit = new TF1("fit","gaus",-5,5);
	fit->SetLineColor(kGray);
	fit->SetLineWidth(3);
	TFitResultPtr r = h1->Fit(fit);
	stringstream mean_val;
	stringstream sigma_val;
	stringstream mean_err;
	stringstream sigma_err;
	mean_val << std::fixed << std::setprecision(3) << fit->GetParameter(1);
	sigma_val << std::fixed << std::setprecision(3) << fit->GetParameter(2);
	mean_err << std::fixed << std::setprecision(3) << fit->GetParError(1);
	sigma_err << std::fixed << std::setprecision(3) << fit->GetParError(2);

	//TString mean = "mean = "+to_string(int(std::ceil(fit->GetParameter(1) * 100)) / 100.);
	//TString sigma = "sigma = "+to_string(int(std::ceil(fit->GetParameter(2) * 100)) / 100.);
	TString mean = "#mu = " + mean_val.str() + "#pm" + mean_err.str();
	TString sigma = "#sigma = " + sigma_val.str() + "#pm" + sigma_err.str();
	
	TLatex latexText(h1->GetXaxis()->GetBinUpEdge(h1->GetXaxis()->GetLast())*0.6, h1->GetMaximum()*0.85, "KinFit");
	latexText.SetTextFont(42);
	latexText.SetTextSize(0.08);
	latexText.SetTextColor(15);
	latexText.DrawClone();

	TLatex latexText1(h1->GetXaxis()->GetBinLowEdge(h1->GetXaxis()->GetFirst())*0.9, h1->GetMaximum()*0.9, mean);
	latexText1.SetTextFont(42);
	latexText1.SetTextSize(0.06);
	latexText1.SetTextColor(kBlack);
	latexText1.DrawClone();

	cout<<mean<<"lalala"<<h1->GetXaxis()->GetBinLowEdge(h1->GetXaxis()->GetFirst())*0.9<<"   "<<h1->GetMaximum()*0.9<<endl;

	TLatex latexText2(h1->GetXaxis()->GetBinLowEdge(h1->GetXaxis()->GetFirst())*0.9, h1->GetMaximum()*0.8, sigma);
	latexText2.SetTextFont(42);
	latexText2.SetTextSize(0.06);
	latexText2.SetTextColor(kBlack);
	latexText2.DrawClone();

	can->Write();
	can->SaveAs("../pics/"+outtitle+".png");
}



void plotsStyleMacro_andFit3(TH1F *h1, TH1F *h2, TH1F *h3, TString outtitle)
{

	TCanvas *can = create_canvas3();
	can->cd();
/*
	TPad *smallpad = new TPad("smallpad","",0.0,0.0,1,1);
	smallpad->SetFillStyle(4000);
	smallpad->Draw();
	smallpad->Divide(3,1,0,0);
*/
	std::vector<TH1F*> histos;
	histos.push_back(h1);
	histos.push_back(h2);
	histos.push_back(h3);

	for(int i=0; i<histos.size(); i++){
		TH1F *h = histos[i];
		plotLayout(h);
		can->cd(i+1);

		h->SetMaximum(h->GetMaximum()*1.2);
		h->Draw("");

		TF1 *fit = new TF1("fit","gaus",-5,5);
		fit->SetLineColor(kGray);
		fit->SetLineWidth(3);
		TFitResultPtr r = h->Fit(fit);
		stringstream mean_val;
		stringstream sigma_val;
		stringstream mean_err;
		stringstream sigma_err;
		mean_val << std::fixed << std::setprecision(3) << fit->GetParameter(1);
		sigma_val << std::fixed << std::setprecision(3) << fit->GetParameter(2);
		mean_err << std::fixed << std::setprecision(3) << fit->GetParError(1);
		sigma_err << std::fixed << std::setprecision(3) << fit->GetParError(2);

		//TString mean = "mean = "+to_string(int(std::ceil(fit->GetParameter(1) * 100)) / 100.);
		//TString sigma = "sigma = "+to_string(int(std::ceil(fit->GetParameter(2) * 100)) / 100.);
		TString mean = "#mu = " + mean_val.str() + "#pm" + mean_err.str();
		TString sigma = "#sigma = " + sigma_val.str() + "#pm" + sigma_err.str();
		
		can->cd(i+1);

		TLatex latexText;
		latexText.SetTextFont(42);
		latexText.SetTextSize(0.08);
		latexText.SetTextColor(15);
		latexText.DrawLatex(h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetLast())*0.5, h->GetMaximum()*0.85, "KinFit");

		/*TLatex latexText1(h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst())*0.9, h->GetMaximum()*0.9, mean);
		latexText1.SetTextFont(42);
		latexText1.SetTextSize(0.06);
		latexText1.SetTextColor(kBlack);
		latexText1.DrawClone();*/
		TLatex latexText1;
		latexText1.SetTextFont(42);
		latexText1.SetTextSize(0.06);
		latexText1.SetTextColor(kBlack);
		latexText1.DrawLatex(h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst())*0.9, h->GetMaximum()*0.9, mean);

		cout<<mean<<"lalala"<<h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst())*0.9<<"   "<<h->GetMaximum()*0.9<<endl;

		/*TLatex latexText2(h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst())*0.9, h->GetMaximum()*0.8, sigma);
		latexText2.SetTextFont(42);
		latexText2.SetTextSize(0.06);
		latexText2.SetTextColor(kBlack);
		latexText2.DrawClone();*/
		TLatex latexText2;
		latexText2.SetTextFont(42);
		latexText2.SetTextSize(0.06);
		latexText2.SetTextColor(kBlack);
		latexText2.DrawLatex(h->GetXaxis()->GetBinLowEdge(h->GetXaxis()->GetFirst())*0.9, h->GetMaximum()*0.8,sigma);
	}

	can->Write();
	can->SaveAs("../pics/"+outtitle+".png");
}

void plotsStyleMacro3(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h12, TH1F *h22, TH1F *h32, TString outtitle, TString h1_title, TString h2_title)
{

	TCanvas *can = create_canvas3();
	can->cd();

	std::vector<TH1F*> histos1;
	histos1.push_back(h1);
	histos1.push_back(h2);
	histos1.push_back(h3);

	std::vector<TH1F*> histos2;
	histos2.push_back(h12);
	histos2.push_back(h22);
	histos2.push_back(h32);

	for(int i=0; i<histos1.size(); i++){
		TH1F *h = histos1[i];
		TH1F *h_fit = histos2[i];
		plotLayout(h);
		h->SetLineStyle(2);
		plotLayout(h_fit);
		can->cd(i+1);

		h->Draw("");
		h_fit->Draw("same");
		
		h->SetMaximum(max(h->GetMaximum(), h_fit->GetMaximum()));
		h->SetMaximum(h->GetMaximum()*1.05);
		
		TLatex latexText;
		latexText.SetTextFont(42);
		latexText.SetTextSize(0.08);
		latexText.SetTextColor(15);
		latexText.DrawLatex(h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetLast())*0.5, h->GetMaximum()*0.85, "KinFit");

		auto legend = new TLegend(0.15,0.6,0.45,0.75);
		legend->AddEntry(h, h1_title);
		legend->AddEntry(h_fit, h2_title);
		legend->Draw();
	}

	can->Write();
	can->SaveAs("../pics/"+outtitle+".png");
}
