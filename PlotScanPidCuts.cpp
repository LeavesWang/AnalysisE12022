#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <cmath>

using namespace std;

#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TEllipse.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMath.h"

void ScanPidCuts()
{
	gStyle->SetOptStat("nemri");
	gStyle->SetStatFormat("8.6g");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetFitFormat("8.6g");
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);

	Option_t *opt = "colz";

	const int N = 200;

	int i = 0, j = 0, k = 0, m = 0;

	string strNucl[N];
	int ZPid[N] = {0}, APid[N] = {0}, QPid[N] = {0};
	double ZPidCut[N][2] = {0}, AoQPidCut[N][2] = {0};
	string strRead;
	int iNucl = 0, nNucl = 0;
	ifstream fCutInfo_Man("fCutInfo_Man.dat");
	while (!fCutInfo_Man.eof() && fCutInfo_Man.peek() != EOF)
	{
		if (fCutInfo_Man.peek() != '#')
		{
			fCutInfo_Man >> strRead >> strRead >> strNucl[nNucl] >> ZPid[nNucl] >> APid[nNucl] >> QPid[nNucl] >> ZPidCut[nNucl][0] >> ZPidCut[nNucl][1] >> AoQPidCut[nNucl][0] >> AoQPidCut[nNucl][1];
			getline(fCutInfo_Man, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo_Man, strRead);
	}
	fCutInfo_Man.close();

	string strDir = "/home/kailong/ExpData/Jul2018/Pid/Graphs/PID/PID/";
	if (access(strDir.c_str(), F_OK) != 0)
		system(("mkdir " + strDir).c_str());
	string strGif = strDir + "cPidCuts.gif";
	if (access(strGif.c_str(), F_OK) == 0)
		system(("rm -f " + strGif).c_str());

	TFile *fCutScan = new TFile("fCutScan.root", "UPDATE");
	TTree *tMaxTof = new TTree("tMaxTof", "Most possible ToF among cuts for each isotope");
	string sNucl = "";
	int Z = 0, A = 0, Q = 0;
	double tofMax = 0, tofCorrMax = 0;
	tMaxTof->Branch("sNucl", &sNucl);
	tMaxTof->Branch("Z", &Z, "Z/I");
	tMaxTof->Branch("A", &A, "A/I");
	tMaxTof->Branch("Q", &Q, "Q/I");
	tMaxTof->Branch("tofMax", &tofMax, "tofMax/D");
	tMaxTof->Branch("tofCorrMax", &tofCorrMax, "tofCorrMax/D");
	
	for (iNucl = 0; iNucl < nNucl; iNucl++)
	// for (iNucl = 62; iNucl < 63; iNucl++)
	{
		if (ZPid[iNucl] <= 37)
			continue;

		sNucl = strNucl[iNucl];
		Z = ZPid[iNucl];
		A = APid[iNucl];
		Q = QPid[iNucl];
		tofMax = 0;
		tofCorrMax = 0;

		TTree *tCutScan = fCutScan->Get<TTree>(("tCutScan_" + strNucl[iNucl]).c_str());
		if (!tCutScan)
		{
			cout << "Error read the tree of tCutScan_" << strNucl[iNucl] << "!\n";
			continue;
		}
		if(!(tCutScan->GetBranch("tofCorrCtr")))
			continue;
		tCutScan->SetEstimate(-1);
		int nData = tCutScan->Draw("Entry$:tofCtr:tofCorrCtr", "abs(tofCtr)>0 && abs(tofCorrCtr)>0 ", "goff");
		double *iCut = tCutScan->GetVal(0);
		double *tofCtr = tCutScan->GetVal(1);
		double *tofCorrCtr = tCutScan->GetVal(2);
		if (nData != tCutScan->GetEntries() || nData<10)
			continue;
		int init = nData / 2;

		TGraph *grTofCuts = new TGraph(nData, iCut, tofCtr);
		grTofCuts->SetNameTitle(("grTofCuts_" + strNucl[iNucl]).c_str(), ("ToF centroids of " + strNucl[iNucl] + " VS cuts;#cut;ToF_{c} [ns]").c_str());
		TGraph *grTofCorrCuts = new TGraph(nData, iCut, tofCorrCtr);
		grTofCorrCuts->SetNameTitle(("grTofCorrCuts_" + strNucl[iNucl]).c_str(), ("Corrected ToF centroids of " + strNucl[iNucl] + " VS cuts;#cut;ToF_{c}^{corr} [ns]").c_str());

		double tofInit = tofCtr[init];
		double tofRang[2] = {round((grTofCuts->GetMean(2) - 4 * grTofCuts->GetRMS(2)) * 1000) / 1000, round((grTofCuts->GetMean(2) + 4 * grTofCuts->GetRMS(2)) * 1000) / 1000};
		double tofCorrInit = tofCorrCtr[init];
		double tofCorrRang[2] = {round((grTofCorrCuts->GetMean(2) - 4 * grTofCorrCuts->GetRMS(2)) * 1000) / 1000, round((grTofCorrCuts->GetMean(2) + 4 * grTofCorrCuts->GetRMS(2)) * 1000) / 1000};

		TH2S *hTofCuts = new TH2S(("hTofCuts_" + strNucl[iNucl]).c_str(), ("ToF centroids of " + strNucl[iNucl] + " VS cuts;#cut;ToF_{c} [ns];").c_str(), nData, 0, nData, 1000 * (tofRang[1] - tofRang[0]), tofRang[0], tofRang[1]);
		TH2S *hTofCorrCuts = new TH2S(("hTofCorrCuts_" + strNucl[iNucl]).c_str(), ("Corrected ToF centroids of " + strNucl[iNucl] + " VS cuts;#cut;ToF_{c}^{corr} [ns];").c_str(), nData, 0, nData, 1000 * (tofCorrRang[1] - tofCorrRang[0]), tofCorrRang[0], tofCorrRang[1]);
		for (int iData = 0; iData < nData; iData++)
		{
			hTofCuts->Fill(iCut[iData], tofCtr[iData]);
			hTofCorrCuts->Fill(iCut[iData], tofCorrCtr[iData]);
		}

		string sCanTofCuts = "cTofCuts_" + strNucl[iNucl];
		TCanvas *cTofCuts = new TCanvas(sCanTofCuts.c_str(), sCanTofCuts.c_str());
		cTofCuts->Divide(2, 2);

		cTofCuts->cd(1);
		hTofCuts->Draw(opt);
		hTofCuts->SetStats(0);
		TLine *xInit = new TLine(init, tofRang[0], init, tofRang[1]);
		TLine *yTofInit = new TLine(0, tofInit, nData, tofInit);
		xInit->SetLineWidth(2);
		xInit->SetLineStyle(2);
		yTofInit->SetLineWidth(2);
		yTofInit->SetLineStyle(2);
		xInit->Draw("same");
		yTofInit->Draw("same");

		cTofCuts->cd(2);
		TH1D *hTof = hTofCuts->ProjectionY();
		hTof->SetNameTitle(("hTof_" + strNucl[iNucl]).c_str(), ("ToF centroids dist. of " + strNucl[iNucl] + " with diff. cuts;ToF_{c} [ns];Counts/ps").c_str());
		hTof->Draw();
		hTof->SetAxisRange(hTof->GetMean() - hTof->GetStdDev(), hTof->GetMean() + hTof->GetStdDev(), "X");
		hTof->GetYaxis()->SetMaxDigits(3);
		TLine *xTofInit = new TLine(tofInit, 0, tofInit, hTof->GetMaximum() * 1.1);
		xTofInit->SetLineWidth(4);
		xTofInit->SetLineStyle(2);
		xTofInit->Draw("same");

		cTofCuts->cd(3);
		hTofCorrCuts->Draw(opt);
		hTofCorrCuts->SetStats(0);
		TLine *xCorrInit = new TLine(init, tofCorrRang[0], init, tofCorrRang[1]);
		TLine *yTofCorrInit = new TLine(0, tofCorrInit, nData, tofCorrInit);
		xCorrInit->SetLineWidth(2);
		xCorrInit->SetLineStyle(2);
		yTofCorrInit->SetLineWidth(2);
		yTofCorrInit->SetLineStyle(2);
		xCorrInit->Draw("same");
		yTofCorrInit->Draw("same");

		cTofCuts->cd(4);
		TH1D *hTofCorr = hTofCorrCuts->ProjectionY();
		hTofCorr->SetNameTitle(("hTofCorr_" + strNucl[iNucl]).c_str(), ("Corrected ToF centroids dist. of " + strNucl[iNucl] + " with diff. cuts;ToF_{c}^{corr} [ns];Counts/ps").c_str());
		hTofCorr->Draw();
		hTofCorr->SetAxisRange(hTofCorr->GetMean() - hTofCorr->GetStdDev(), hTofCorr->GetMean() + hTofCorr->GetStdDev(), "X");
		hTofCorr->GetYaxis()->SetMaxDigits(3);
		TLine *xTofCorrInit = new TLine(tofCorrInit, 0, tofCorrInit, hTofCorr->GetMaximum() * 1.1);
		xTofCorrInit->SetLineWidth(4);
		xTofCorrInit->SetLineStyle(2);
		xTofCorrInit->Draw("same");

		tofMax = hTof->GetXaxis()->GetBinCenter(hTof->GetMaximumBin());
		tofCorrMax = hTofCorr->GetXaxis()->GetBinCenter(hTofCorr->GetMaximumBin());
		tMaxTof->Fill();

		cTofCuts->cd(1);
		TLine *yTofMax = new TLine(0, tofMax, nData, tofMax);
		yTofMax->SetLineWidth(2);
		yTofMax->SetLineStyle(2);
		yTofMax->SetLineColor(kRed);
		yTofMax->Draw("same");

		cTofCuts->cd(3);
		TLine *yTofCorrMax = new TLine(0, tofCorrMax, nData, tofCorrMax);
		yTofCorrMax->SetLineWidth(2);
		yTofCorrMax->SetLineStyle(2);
		yTofCorrMax->SetLineColor(kRed);
		yTofCorrMax->Draw("same");

		cTofCuts->Write("", TObject::kOverwrite);
		cTofCuts->SaveAs((strDir + sCanTofCuts + ".png").c_str());
		cTofCuts->SaveAs((strGif + "+50").c_str());

		delete yTofCorrMax;
		delete yTofMax;
		delete xTofCorrInit;
		delete yTofCorrInit;
		delete xCorrInit;
		delete xTofInit;
		delete yTofInit;
		delete xInit;
		delete cTofCuts;
		delete hTofCuts;
		delete grTofCuts;
	}
	tMaxTof->Write("", TObject::kOverwrite);
	fCutScan->Close();
	delete fCutScan;
}

int main(int argc, char **argv)
{
	ScanPidCuts();
	return 0;
}