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

void PlotTofElec()
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
	string sElec[2] = {"TacNoClk", "TacClk"};

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

	TFile *fPid = new TFile("pid-run-150--153__270--385.root");
	TTree *tPid = fPid->Get<TTree>("tPid");
	tPid->SetEstimate(-1);

	int good_CfdTdc = 0, good_TacNoClk = 0, good_TacClk = 0;
	double tof_CfdTdc = 0, tof_TacNoClk = 0, tof_TacClk = 0;
	string *nucl = new string;

	tPid->SetBranchAddress("nucl", &nucl);
	tPid->SetBranchAddress("good_CfdTdc", &good_CfdTdc);
	tPid->SetBranchAddress("tof_CfdTdc", &tof_CfdTdc);
	tPid->SetBranchAddress("good_TacNoClk", &good_TacNoClk);
	tPid->SetBranchAddress("tof_TacNoClk", &tof_TacNoClk);
	tPid->SetBranchAddress("good_TacClk", &good_TacClk);
	tPid->SetBranchAddress("tof_TacClk", &tof_TacClk);

	string strDir = "/home/kailong/ExpData/Jul2018/Pid/Graphs/ElecDif/";
	if (access(strDir.c_str(), F_OK) != 0)
		system(("mkdir " + strDir).c_str());
	string strGif = strDir + "cTofElec.gif";
	if (access(strGif.c_str(), F_OK) == 0)
		system(("rm -f " + strGif).c_str());

	long long iEnt = 0, nEnt = tPid->GetEntries();
	TGraph *grTofDif[2], *grTofDifNucl[N][2];
	TH2F *hTofDif[2];
	for (k = 0; k < 2; k++)
	{
		grTofDif[k] = new TGraph();
		hTofDif[k] = new TH2F(("hTofDif" + to_string(k)).c_str(), ("tof_" + sElec[k] + "-tof_CfdTdc VS tof_CfdTdc;ToF_{CfdTdc} [ns];ToF_{TacNoClk}-ToF_{" + sElec[k] + "} [ps]").c_str(), 4000, 460, 500, 120, -150, 150);
		for (iNucl = 0; iNucl < nNucl; iNucl++)
		{
			grTofDifNucl[iNucl][k] = new TGraph();
		}
	}
	int nData[2] = {0}, cnts[N][2] = {0};
	for (iEnt = 0; iEnt < nEnt; iEnt++)
	{
		tPid->GetEntry(iEnt);
		if (good_CfdTdc == 1 && good_TacNoClk == 1)
		{
			grTofDif[0]->SetPoint(nData[0]++, tof_CfdTdc, 1000 * (tof_TacNoClk - tof_CfdTdc));
			for (iNucl = 0; iNucl < nNucl; iNucl++)
				if (*nucl == strNucl[iNucl])
					grTofDifNucl[iNucl][0]->SetPoint(cnts[iNucl][0]++, tof_CfdTdc, 1000 * (tof_TacNoClk - tof_CfdTdc));
		}
		if (good_CfdTdc == 1 && good_TacClk == 1)
		{
			grTofDif[1]->SetPoint(nData[1]++, tof_CfdTdc, 1000 * (tof_TacClk - tof_CfdTdc));
			for (iNucl = 0; iNucl < nNucl; iNucl++)
				if (*nucl == strNucl[iNucl])
					grTofDifNucl[iNucl][1]->SetPoint(cnts[iNucl][1]++, tof_CfdTdc, 1000 * (tof_TacClk - tof_CfdTdc));
		}
	}
	double meanTofDifNucl[N][2] = {0};
	TH2F *hTofDifNucl[N][2];
	for (k = 0; k < 2; k++)
		for (iNucl = 0; iNucl < nNucl; iNucl++)
			if(cnts[iNucl][k]>0)
			{
				meanTofDifNucl[iNucl][k] = grTofDifNucl[iNucl][k]->GetMean(2);
				double rmsTofDif = grTofDifNucl[iNucl][k]->GetRMS(2);
				double meanTof = grTofDifNucl[iNucl][k]->GetMean(1);
				double rmsTof = grTofDifNucl[iNucl][k]->GetRMS(1);
				hTofDifNucl[iNucl][k] = new TH2F(("hTofDif" + to_string(k) + "_" + strNucl[iNucl]).c_str(), ("tof_" + sElec[k] + "-tof_CfdTdc VS tof_CfdTdc for " + strNucl[iNucl] + ";ToF_{CfdTdc} [ns];ToF_{TacNoClk}-ToF_{" + sElec[k] + "} [ps]").c_str(), 6 * rmsTof * 50, meanTof - 3 * rmsTof, meanTof + 3 * rmsTof, 8 * rmsTofDif, -4 * rmsTofDif, 4 * rmsTofDif);
			}
	double meanTofDif[2] = {0};
	for (k = 0; k < 2; k++)
	{
		meanTofDif[k] = grTofDif[k]->GetMean(2);
		for (i = 0; i < nData[k]; i++)
			hTofDif[k]->Fill(grTofDif[k]->GetX()[i], grTofDif[k]->GetY()[i] - meanTofDif[k]);
		for (iNucl = 0; iNucl < nNucl; iNucl++)
			for (i = 0; i < cnts[iNucl][k]; i++)
				hTofDifNucl[iNucl][k]->Fill(grTofDifNucl[iNucl][k]->GetX()[i], grTofDifNucl[iNucl][k]->GetY()[i] - meanTofDifNucl[iNucl][k]);
	}
	TFile *fTofElecDif = new TFile("fTofElecDif.root", "RECREATE");

	string sCanAll = "cTofElecDif_All";
	TCanvas *cAll = new TCanvas(sCanAll.c_str(), sCanAll.c_str());
	cAll->Divide(2, 2);
	TProfile *hTofDif_pfx[2];
	for (k = 0; k < 2; k++)
	{
		cAll->cd(k + 1);
		hTofDif[k]->Draw(opt);
		cAll->cd(k + 3);
		hTofDif_pfx[k] = hTofDif[k]->ProfileX("_pfx");
		hTofDif_pfx[k]->SetTitle(("Profile of above 2D histogram;ToF_{CfdTdc} [ns];ToF_{" + sElec[k] + "}-ToF_{CfdTdc} [ps]").c_str());
		hTofDif_pfx[k]->SetAxisRange(-70, 70, "Y");
		hTofDif_pfx[k]->SetStats(kFALSE);
		hTofDif_pfx[k]->Draw();
	}
	cAll->Write("", TObject::kOverwrite);
	cAll->SaveAs((strDir + sCanAll + ".png").c_str());
	delete cAll;

	for(iNucl=0; iNucl<nNucl; iNucl++)
	// for (iNucl = 0; iNucl < 5; iNucl++)
		if(cnts[iNucl][0]>0 || cnts[iNucl][1]>0)
		{
			string sCanNucl = "cTofElecDif_" + strNucl[iNucl];
			TCanvas *cNucl = new TCanvas(sCanNucl.c_str(), sCanNucl.c_str());
			cNucl->Divide(2, 2);
			for (k = 0; k < 2; k++)
				if(cnts[iNucl][k]>0)
				{
					cNucl->cd(k + 1);
					hTofDifNucl[iNucl][k]->Draw(opt);
					cNucl->cd(k + 3);
					hTofDif_pfx[k]->Clear();
					hTofDif_pfx[k]->Reset("ICESM");
					hTofDif_pfx[k] = hTofDifNucl[iNucl][k]->ProfileX("_pfx");
					hTofDif_pfx[k]->SetTitle(("Profile of above 2D histogram;ToF_{CfdTdc} [ns];ToF_{" + sElec[k] + "}-ToF_{CfdTdc} [ps]").c_str());
					hTofDif_pfx[k]->SetStats(kFALSE);
					hTofDif_pfx[k]->Draw();
				}
			cNucl->Write("", TObject::kOverwrite);
			cNucl->SaveAs((strDir + sCanNucl + ".png").c_str());
			cNucl->SaveAs((strGif + "+50").c_str());
			delete cNucl;
		}

	fTofElecDif->Close();
	delete fTofElecDif;
	for (k = 0; k < 2; k++)
	{
		for (iNucl = 0; iNucl < nNucl; iNucl++)
		{
			delete grTofDifNucl[iNucl][k];
			if(cnts[iNucl][k]>0)
				delete hTofDifNucl[iNucl][k];
		}
		delete hTofDif[k];
		delete grTofDif[k];
	}
	delete nucl;
	fPid->Close();
	delete fPid;
}

int main(int argc, char **argv)
{
	PlotTofElec();
	return 0;
}