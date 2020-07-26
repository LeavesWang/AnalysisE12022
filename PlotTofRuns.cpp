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

void PlotTofRuns()
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

	string strNucl[N];
	string strRead;
	int iNucl=0, nNucl=0;
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

	int runNum=0, good_CfdTdc=0;
	double tof_CfdTdc=0;
	string *eSet = new string;
	string *nucl = new string;

	tPid->SetBranchAddress("eSet", &eSet);
	tPid->SetBranchAddress("runNum", &runNum);
	tPid->SetBranchAddress("nucl", &nucl);
	tPid->SetBranchAddress("good_CfdTdc", &good_CfdTdc);
	tPid->SetBranchAddress("tof_CfdTdc", &tof_CfdTdc);

	string strDir = "/home/kailong/ExpData/Jul2018/Pid/Graphs/CfdTdc/";
	if (access(strDir.c_str(), F_OK) != 0)
		system(("mkdir " + strDir).c_str());
	string strGif = strDir + "cTofRuns.gif";
	if (access(strGif.c_str(), F_OK) == 0)
		system(("rm -f " + strGif).c_str());

	long long iEnt = 0, nEnt=tPid->GetEntries();
	TH2F *hTofPSrun = new TH2F("hTofPSrun", "ToF from CFD+TDC VS run number;#run;ToF [ns]", 121, 266, 387, 40000, 460, 500);
	TH2F *hTofRSrun = new TH2F("hTofRSrun", "ToF from CFD+TDC VS run number;#run;ToF [ns]", 121, 266, 387, 40000, 460, 500);
	TH2F *hTofPSrunNucl[N], *hTofRSrunNucl[N];
	TGraph *grTofPSrunNucl[N], *grTofRSrunNucl[N];
	for(iNucl=0; iNucl<nNucl; iNucl++)
	{
		hTofPSrunNucl[iNucl] = new TH2F(("hTofPSrun_"+strNucl[iNucl]).c_str(), ("ToF of "+strNucl[iNucl]+" from CFD+TDC VS run number;#run;ToF [ns]").c_str(), 121, 266, 387, 40000, 460, 500);
		hTofRSrunNucl[iNucl] = new TH2F(("hTofRSrun_"+strNucl[iNucl]).c_str(), ("ToF of "+strNucl[iNucl]+" from CFD+TDC VS run number;#run;ToF [ns]").c_str(), 121, 266, 387, 40000, 460, 500);

		grTofPSrunNucl[iNucl] = new TGraph();
		grTofPSrunNucl[iNucl]->SetName(("grTofPSrun_"+strNucl[iNucl]).c_str());
		grTofRSrunNucl[iNucl] = new TGraph();
		grTofRSrunNucl[iNucl]->SetName(("grTofRSrun_"+strNucl[iNucl]).c_str());
	}
	int cnts[N][2]={0};
	for(iEnt=0; iEnt<nEnt; iEnt++)
	{
		tPid->GetEntry(iEnt);
		if(good_CfdTdc==1 && runNum>0 && tof_CfdTdc>0)
		{
			if(*eSet == "PS_270_382")
			{
				hTofPSrun->Fill(runNum, tof_CfdTdc);
				for(iNucl=0; iNucl<nNucl; iNucl++)
					if(*nucl == strNucl[iNucl])
					{
						hTofPSrunNucl[iNucl]->Fill(runNum, tof_CfdTdc);
						grTofPSrunNucl[iNucl]->SetPoint(cnts[iNucl][0]++, runNum, tof_CfdTdc);
					}
			}
			if(*eSet == "RS_270_382")
			{
				hTofRSrun->Fill(runNum, tof_CfdTdc);
				for(iNucl=0; iNucl<nNucl; iNucl++)
					if(*nucl == strNucl[iNucl])
					{
						hTofRSrunNucl[iNucl]->Fill(runNum, tof_CfdTdc);
						grTofRSrunNucl[iNucl]->SetPoint(cnts[iNucl][1]++, runNum, tof_CfdTdc);
					}
			}
		}
	}
	TFile *fTofRuns = new TFile("fTofRuns.root", "RECREATE");
	string sCanAll="cTofVSRun_All";
	TCanvas *cAll = new TCanvas(sCanAll.c_str(), sCanAll.c_str(), 2400, 900);
	cAll->Divide(2,1);

	cAll->cd(1);
	hTofPSrun->Draw(opt);
	// hTofPSrun->SetAxisRange(hTofPSrun->GetMean(2) - 3 * hTofPSrun->GetStdDev(2), hTofPSrun->GetMean(2) + 3 * hTofPSrun->GetStdDev(2), "Y");
	hTofPSrun->SetStats(0);

	hTofRSrun->Draw("colz same");
	// hTofRSrun->SetAxisRange(hTofRSrun->GetMean(2) - 3 * hTofRSrun->GetStdDev(2), hTofRSrun->GetMean(2) + 3 * hTofRSrun->GetStdDev(2), "Y");
	hTofRSrun->SetStats(0);
	
	cAll->cd(2);
	TProfile *hTofPSrun_pfx = hTofPSrun->ProfileX("_pfx");
	hTofPSrun_pfx->SetTitle("ToF centroid VS run number;#run;ToF_{c} [ns]");
	hTofPSrun_pfx->SetStats(0);
	hTofPSrun_pfx->SetMarkerStyle(7);
	hTofPSrun_pfx->Draw();
	hTofPSrun_pfx->SetAxisRange(480.4, 482.2, "Y");

	TProfile *hTofRSrun_pfx = hTofRSrun->ProfileX("_pfx");
	hTofRSrun_pfx->SetTitle("ToF centroid VS run number;#run;ToF_{c} [ns]");
	hTofRSrun_pfx->SetStats(0);
	hTofRSrun_pfx->SetMarkerStyle(3);
	hTofRSrun_pfx->SetLineColor(kRed);
	hTofRSrun_pfx->SetMarkerColor(kRed);
	hTofRSrun_pfx->Draw("same");
	hTofRSrun_pfx->SetAxisRange(480.4, 482.2, "Y");

	cAll->Write("", TObject::kOverwrite);
	cAll->SaveAs((strDir+sCanAll+".png").c_str());

	delete cAll;

	for(iNucl=0; iNucl<nNucl; iNucl++)
	// for(iNucl=0; iNucl<5; iNucl++)
		if(cnts[iNucl][0]>0 || cnts[iNucl][1]>0)
		{
			string sCanNucl="cTofVSRun_"+strNucl[iNucl];
			TCanvas *cNucl = new TCanvas(sCanNucl.c_str(), sCanNucl.c_str(), 2500, 900);
			cNucl->Divide(2,1);

			cNucl->cd(1);
			TMultiGraph *mgrTofRunNucl = new TMultiGraph(("mgrTofRun_"+strNucl[iNucl]).c_str(), ("ToF of "+strNucl[iNucl]+" from CFD+TDC VS run number;#run;ToF [ns]").c_str());

			double tofRang[2]={0};

			grTofPSrunNucl[iNucl]->SetMarkerStyle(7);
			if(cnts[iNucl][0]>0)
			{
				mgrTofRunNucl->Add(grTofPSrunNucl[iNucl]);
				tofRang[0]=grTofPSrunNucl[iNucl]->GetMean(2)-2*grTofPSrunNucl[iNucl]->GetRMS(2);
				tofRang[1]=grTofPSrunNucl[iNucl]->GetMean(2)+2*grTofPSrunNucl[iNucl]->GetRMS(2);
			}
			grTofRSrunNucl[iNucl]->SetMarkerStyle(3);
			grTofRSrunNucl[iNucl]->SetMarkerColor(kRed);
			if(cnts[iNucl][1]>0)
			{
				mgrTofRunNucl->Add(grTofRSrunNucl[iNucl]);
				if(cnts[iNucl][0]==0 || grTofRSrunNucl[iNucl]->GetMean(2)-2*grTofRSrunNucl[iNucl]->GetRMS(2)<tofRang[0])
					tofRang[0]=grTofRSrunNucl[iNucl]->GetMean(2)-2*grTofRSrunNucl[iNucl]->GetRMS(2);
				if(cnts[iNucl][0]==0 || grTofRSrunNucl[iNucl]->GetMean(2)+2*grTofRSrunNucl[iNucl]->GetRMS(2)>tofRang[1])
					tofRang[1]=grTofRSrunNucl[iNucl]->GetMean(2)+2*grTofRSrunNucl[iNucl]->GetRMS(2);
			}
			mgrTofRunNucl->Draw("AP");

			cNucl->cd(2);
			if(cnts[iNucl][0]>0)
			{
				TProfile *hTofPSrunNucl_pfx = hTofPSrunNucl[iNucl]->ProfileX(("hTofPSrun_"+strNucl[iNucl]+"_pfx").c_str());
				hTofPSrunNucl_pfx->SetTitle(("ToF centroid of "+strNucl[iNucl]+" VS run number;#run;ToF_{c} [ns]").c_str());
				hTofPSrunNucl_pfx->SetStats(0);
				hTofPSrunNucl_pfx->SetMarkerStyle(7);
				hTofPSrunNucl_pfx->Draw();
				hTofPSrunNucl_pfx->SetAxisRange(tofRang[0], tofRang[1], "Y");
			}
			if(cnts[iNucl][1]>0)
			{
				TProfile *hTofRSrunNucl_pfx = hTofRSrunNucl[iNucl]->ProfileX(("hTofRSrun_"+strNucl[iNucl]+"_pfx").c_str());
				hTofRSrunNucl_pfx->SetTitle(("ToF centroid of "+strNucl[iNucl]+" VS run number;#run;ToF_{c} [ns]").c_str());
				hTofRSrunNucl_pfx->SetStats(0);
				hTofRSrunNucl_pfx->SetMarkerStyle(3);
				hTofRSrunNucl_pfx->SetLineColor(kRed);
				hTofRSrunNucl_pfx->SetMarkerColor(kRed);
				if(cnts[iNucl][0]>0)
					hTofRSrunNucl_pfx->Draw("same");
				else
					hTofRSrunNucl_pfx->Draw();
				hTofRSrunNucl_pfx->SetAxisRange(tofRang[0], tofRang[1], "Y");
			}		

			cNucl->Write("", TObject::kOverwrite);
			cNucl->SaveAs((strDir+sCanNucl+".png").c_str());
			cNucl->SaveAs((strGif+"+50").c_str());

			delete grTofRSrunNucl[iNucl];
			delete grTofPSrunNucl[iNucl];
			delete mgrTofRunNucl;
			delete cNucl;
		}
	fTofRuns->Close();
	delete fTofRuns;
	for(iNucl=0; iNucl<nNucl; iNucl++)
	{
		delete hTofRSrunNucl[iNucl];
		delete hTofPSrunNucl[iNucl];
	}
	delete hTofRSrun;
	delete hTofPSrun;
	delete nucl;
	delete eSet;
	fPid->Close();
	delete fPid;
}

int main(int argc, char **argv)
{
	PlotTofRuns();
	return 0;
}