#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <cmath>

using namespace std;

#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TMultiGraph.h"

void PlotTofQpmt()
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

	const int N = 200, NSGM = 4;

	int i = 0, j = 0, k = 0;

	string strNucl[N];
	string strRead;
	int iNucl = 0, nNucl = 0;
	ifstream fCutInfo("fCutInfo.dat");
	while (!fCutInfo.eof() && fCutInfo.peek() != EOF)
	{
		if (fCutInfo.peek() != '#')
		{
			fCutInfo >> strRead >> strRead >> strNucl[nNucl];
			getline(fCutInfo, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo, strRead);
	}
	fCutInfo.close();

	TFile *fPid = new TFile("tof-run-150--153__270--385.root");
	TTree *tTOFDetail = fPid->Get<TTree>("tTOFDetail");
	tTOFDetail->SetEstimate(-1);

	string *sNucl = new string;
	bool good_CfdTdc = 0;
	double tof_Tdc = 0, qPlaS800 = 0, qPmtS800[4] = {0};
	tTOFDetail->SetBranchAddress("sNucl", &sNucl);
	tTOFDetail->SetBranchAddress("good_CfdTdc", &good_CfdTdc);
	tTOFDetail->SetBranchAddress("tof2XsCor_CfdTdc", &tof_Tdc);
	// tTOFDetail->SetBranchAddress("tof1XmcpCor_CfdTdc", &tof_Tdc);
	// tTOFDetail->SetBranchAddress("tof0_CfdTdc", &tof_Tdc);
	tTOFDetail->SetBranchAddress("qPlaS800", &qPlaS800);
	tTOFDetail->SetBranchAddress("qPmtS800", qPmtS800);

	string strDir = "/home/kailong/ExpData/Jul2018/Pid/Graphs/TimeWalk/";
	if (access(strDir.c_str(), F_OK) != 0)
		system(("mkdir " + strDir).c_str());
	string strGif = strDir + "cTofQpmt.gif";
	if (access(strGif.c_str(), F_OK) == 0)
		system(("rm -f " + strGif).c_str());

	TGraph *grTofQpla[N], *grTofQpmt[N][4];
	for (iNucl = 0; iNucl < nNucl; iNucl++)
	{
		grTofQpla[iNucl] = new TGraph();
		for (k = 0; k < 4; k++)
			grTofQpmt[iNucl][k] = new TGraph();
	}
	long long iEnt = 0, nEnt = tTOFDetail->GetEntries();
	int cntsPla[N] = {0}, cntsPmt[N][4] = {0};
	for (iEnt = 0; iEnt < nEnt; iEnt++)
	{
		tTOFDetail->GetEntry(iEnt);
		if (good_CfdTdc == 1)
			for (iNucl = 0; iNucl < nNucl; iNucl++)
				if (*sNucl == strNucl[iNucl])
				{
					if (qPlaS800 > 0)
					{
						// hTofQpla[iNucl]->Fill(qPlaS800, tof_Tdc);
						grTofQpla[iNucl]->SetPoint(cntsPla[iNucl]++, qPlaS800, tof_Tdc);
					}
					for (k = 0; k < 4; k++)
						if (qPmtS800[k] > 0)
						{
							// hTofQpmt[iNucl][k]->Fill(qPmtS800[k], tof_Tdc);
							grTofQpmt[iNucl][k]->SetPoint(cntsPmt[iNucl][k]++, qPmtS800[k], tof_Tdc);
						}
				}
	}
	TFile *fTofQpmt = new TFile("fTofQpmt.root", "RECREATE");
	double egyRang[2] = {0}, tofRang[2] = {0};
	int nBinEgy = 0, nBinTof = 0;
	for (iNucl = 0; iNucl < nNucl; iNucl++)
	// for (iNucl = 0; iNucl < 5; iNucl++)
		if (cntsPla[iNucl] > 0)
		{
			string sCan = "cTofQpmt_" + strNucl[iNucl];
			TCanvas *cTof = new TCanvas(sCan.c_str(), sCan.c_str(), 1800, 900);
			cTof->Divide(3, 2);

			cTof->cd(3);
			egyRang[0] = grTofQpla[iNucl]->GetMean(1) - NSGM * grTofQpla[iNucl]->GetRMS(1);
			egyRang[1] = grTofQpla[iNucl]->GetMean(1) + NSGM * grTofQpla[iNucl]->GetRMS(1);
			tofRang[0] = grTofQpla[iNucl]->GetMean(2) - NSGM * grTofQpla[iNucl]->GetRMS(2);
			tofRang[1] = grTofQpla[iNucl]->GetMean(2) + NSGM * grTofQpla[iNucl]->GetRMS(2);
			nBinEgy = (int)(0.25*(egyRang[1] - egyRang[0]));
			nBinTof = (int)(100 * (tofRang[1] - tofRang[0]));
			TH2F *hTofQpla = new TH2F(("hTofQpla_" + strNucl[iNucl]).c_str(), ("ToF VS amplitude in S800 plastic for " + strNucl[iNucl] + ";Q_{Pla} [a.u.];ToF [ns]").c_str(), nBinEgy, egyRang[0], egyRang[1], nBinTof, tofRang[0], tofRang[1]);
			for (i = 0; i < cntsPla[iNucl]; i++)
				hTofQpla->Fill(grTofQpla[iNucl]->GetX()[i], grTofQpla[iNucl]->GetY()[i]);
			hTofQpla->Draw(opt);
			TProfile *hTofQpla_pfx = hTofQpla->ProfileX("_pfx");
			// hTofQpla_pfx->SetTitle("Profile of above 2D histogram;Q_{Pla} [a.u.];ToF [ns]");
			hTofQpla_pfx->SetStats(kFALSE);
			hTofQpla_pfx->SetMarkerColor(kRed);
			hTofQpla_pfx->SetMarkerStyle(3);
			hTofQpla_pfx->SetLineColor(kRed);
			hTofQpla_pfx->Draw("same");

			TH2F *hTofQpmt[4];
			TProfile *hTofQpmt_pfx[4];
			for (k = 0; k < 4; k++)
			{
				cTof->cd(k + 1 + k / 2);
				egyRang[0] = grTofQpmt[iNucl][k]->GetMean(1) - NSGM * grTofQpmt[iNucl][k]->GetRMS(1);
				egyRang[1] = grTofQpmt[iNucl][k]->GetMean(1) + NSGM * grTofQpmt[iNucl][k]->GetRMS(1);
				tofRang[0] = grTofQpmt[iNucl][k]->GetMean(2) - NSGM * grTofQpmt[iNucl][k]->GetRMS(2);
				tofRang[1] = grTofQpmt[iNucl][k]->GetMean(2) + NSGM * grTofQpmt[iNucl][k]->GetRMS(2);
				nBinEgy = (int)(0.25*(egyRang[1] - egyRang[0]));
				nBinTof = (int)(100 * (tofRang[1] - tofRang[0]));
				hTofQpmt[k] = new TH2F(("hTofQpmt" + to_string(k + 1) + "_" + strNucl[iNucl]).c_str(), ("ToF VS amplitude of PMT# " + to_string(k + 1) + " in S800 for " + strNucl[iNucl] + ";Q_{PMT#" + to_string(k + 1) + "} [a.u.];ToF [ns]").c_str(), nBinEgy, egyRang[0], egyRang[1], nBinTof, tofRang[0], tofRang[1]);
				for (i = 0; i < cntsPmt[iNucl][k]; i++)
					hTofQpmt[k]->Fill(grTofQpmt[iNucl][k]->GetX()[i], grTofQpmt[iNucl][k]->GetY()[i]);
				hTofQpmt[k]->Draw(opt);
				hTofQpmt_pfx[k] = hTofQpmt[k]->ProfileX("_pfx");
				// hTofQpmt_pfx[k]->SetTitle("Profile of above 2D histogram;Q_{Pla} [a.u.];ToF [ns]");
				hTofQpmt_pfx[k]->SetStats(kFALSE);
				hTofQpmt_pfx[k]->SetMarkerColor(kRed);
				hTofQpmt_pfx[k]->SetMarkerStyle(3);
				hTofQpmt_pfx[k]->SetLineColor(kRed);
				hTofQpmt_pfx[k]->Draw("same");
			}
			cTof->Write("", TObject::kOverwrite);
			cTof->SaveAs((strDir + sCan + ".png").c_str());
			cTof->SaveAs((strGif + "+50").c_str());

			for (k = 0; k < 4; k++)
				delete hTofQpmt[k];
			delete hTofQpla;
			delete cTof;
		}
	fTofQpmt->Close();
	delete fTofQpmt;
	for (iNucl = 0; iNucl < nNucl; iNucl++)
	{
		delete grTofQpla[iNucl];
		for (k = 0; k < 4; k++)
			delete grTofQpmt[iNucl][k];
	}
	delete sNucl;
	fPid->Close();
	delete fPid;
}

int main(int argc, char **argv)
{
	PlotTofQpmt();
	return 0;
}