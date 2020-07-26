#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <cmath>

using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TTree.h"
#include "TGraph.h"

void PlotTofXmcp()
{
	gStyle->SetOptStat("nemri");
	gStyle->SetStatFormat(".6f");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetFitFormat(".6f");
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);

	Option_t *opt = "colz scat";

	const int NSGM = 3, MINPOIN = 100;
	const int N = 200;
	int i, j, k;

	int iNucl, nNucl, isForCal[N] = {0};
	string strNucl[N];
	int ZOrig[N] = {0}, AOrig[N] = {0}, QOrig[N] = {0};
	string strRead;

	nNucl = 0;
	ifstream fCutInfo("fCutInfo.dat");
	while (!fCutInfo.eof() && fCutInfo.peek() != EOF)
	{
		if (fCutInfo.peek() != '#')
		{
			fCutInfo >> isForCal[nNucl] >> strRead >> strNucl[nNucl] >> ZOrig[nNucl] >> AOrig[nNucl] >> QOrig[nNucl];
			getline(fCutInfo, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo, strRead);
	}
	fCutInfo.close();

	string sRoot = "/home/kailong/ExpData/Jul2018/Pid/pid-run-150--153__270--385.root";
	string sTree = "tPid";
	TFile *fRoot = new TFile(sRoot.c_str());
	TTree *tPid;
	fRoot->GetObject(sTree.c_str(), tPid);
	tPid->SetEstimate(-1);

	int nData=0, iData=0, nPoin=0;
	string sDraw, sCut;

	TFile *fDraw = new TFile("fDraw.root", "UPDATE");

	FILE *fXTof = fopen("fXTof.dat", "w+");
	const int NAna = 3;
	string sAna[3] = {"CfdTdc", "TacNoClk", "TacClk"};

	for (int iAna = 0; iAna < NAna; iAna++)
	{
		fprintf(fXTof, "\n**********%s**********\n", sAna[iAna].c_str());
		string strElec = "./Graphs/" + sAna[iAna];
		if (access(strElec.c_str(), F_OK) != 0)
			system(("mkdir " + strElec).c_str());

		sDraw = "ZPid:APid:QPid:tof_" + sAna[iAna] + ":xMCP";
		sCut = "good_" + sAna[iAna] + "==1";
		nData = tPid->Draw(sDraw.c_str(), sCut.c_str(), "goff");
		double *ZPid = tPid->GetVal(0);
		double *APid = tPid->GetVal(1);
		double *QPid = tPid->GetVal(2);
		double *tof = tPid->GetVal(3);
		double *xMcp = tPid->GetVal(4);
				
		bool *goodData = new bool [nData];

		for (iNucl = 0; iNucl < nNucl; iNucl++)
		{
			printf("\n**********%s with %s**********\n", strNucl[iNucl].c_str(), sAna[iAna].c_str());

			for(iData=0; iData<nData; iData++)
				goodData[iData]=false;

			TGraph *grTofXmcp = new TGraph();
			grTofXmcp->SetNameTitle(("grTofVSXmcp__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("grTofVSXmcp__" + sAna[iAna] + "_" + strNucl[iNucl] + ";X_{MCP} [mm];ToF [ns]").c_str());

			nPoin=0;
			for(iData=0; iData<nData; iData++)
				if (abs(ZPid[iData] - ZOrig[iNucl]) < 1E-3 && abs(APid[iData] - AOrig[iNucl]) < 1E-3 && abs(QPid[iData] - QOrig[iNucl]) < 1E-3 && abs(xMcp[iData]) > 0 && abs(tof[iData]) > 0)
				{
					goodData[iData]=true;
					grTofXmcp->SetPoint(nPoin++, xMcp[iData], tof[iData]);
				}
			if (nPoin < MINPOIN)
				continue;

			double xLow=grTofXmcp->GetMean(1)-NSGM*grTofXmcp->GetRMS(1);
			double xUp=grTofXmcp->GetMean(1)+NSGM*grTofXmcp->GetRMS(1);
			double tofLow=grTofXmcp->GetMean(2)-NSGM*grTofXmcp->GetRMS(2);
			double tofUp=grTofXmcp->GetMean(2)+NSGM*grTofXmcp->GetRMS(2);
			grTofXmcp->Set(0);
			nPoin=0;
			for (iData = 0; iData < nData; iData++)
				if(goodData[iData] && xLow<xMcp[iData] && xMcp[iData] < xUp && tofLow<tof[iData] && tof[iData] < tofUp)
					grTofXmcp->SetPoint(nPoin++, xMcp[iData], tof[iData]);
			if (nPoin < MINPOIN)
				continue;
			double xMean=grTofXmcp->GetMean(1);
			double xRms=grTofXmcp->GetRMS(1);
			double xMeanErr=xRms/sqrt(nPoin);

			double tofMean=grTofXmcp->GetMean(2);
			double tofRms=grTofXmcp->GetRMS(2);
			double tofMeanErr=tofRms/sqrt(nPoin);

			tofLow=tofMean-NSGM*tofRms;
			tofUp=tofMean+NSGM*tofRms;
			int nBinTof=(int)((tofUp - tofLow) * 100);

			TCanvas *cTofXmcp = new TCanvas(("cTofXmcp_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("cTofXmcp_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
			cTofXmcp->Divide(2, 2);

			cTofXmcp->cd(1);
			TH2F *hTofXmcp = new TH2F(("hTofVSXmcp__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTofVSXmcp__" + sAna[iAna] + "_" + strNucl[iNucl] + ";X_{MCP} [mm];ToF [ns]").c_str(), 800, -40, 40, nBinTof, tofLow, tofUp);
			double xMcpVal=0, tofVal=0;
			for (iData = 0; iData < nPoin; iData++)
			{
				grTofXmcp->GetPoint(iData, xMcpVal, tofVal);
				hTofXmcp->Fill(xMcpVal, tofVal);
			}
			hTofXmcp->Draw(opt);
			
			cTofXmcp->cd(2);
			TH1D *hTof = hTofXmcp->ProjectionY(("hTof__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
			hTof->SetTitle(("hTof__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF [ns];Counts/0.01ns").c_str());
			hTof->Draw();

			cTofXmcp->cd(3);
			TH1D *hXmcp= hTofXmcp->ProjectionX(("hXmcp__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
			hXmcp->SetTitle(("hXmcp__" + sAna[iAna] + "_" + strNucl[iNucl] + ";X_{MCP} [mm];Counts/0.1mm").c_str());
			hXmcp->Draw();

			// if(isForCal[iNucl]==1)
			fprintf(fXTof, "%d  %d  %s %d %d %d  %f %f  %f %f\n", isForCal[iNucl], iNucl, strNucl[iNucl].c_str(), ZOrig[iNucl], AOrig[iNucl], QOrig[iNucl], xMean, xMeanErr, tofMean, tofMeanErr);

			cTofXmcp->Write("", TObject::kOverwrite);
			cTofXmcp->SaveAs(("./Graphs/" + sAna[iAna] + "/cTofXmcp_" + sAna[iAna] + "_" + strNucl[iNucl] + ".png").c_str());

			delete hTofXmcp;
			delete grTofXmcp;
			delete cTofXmcp;
		}
		delete [] goodData;
	}
	delete fDraw;
	delete fRoot;
	fclose(fXTof);
}

int main(int argc, char **argv)
{
	PlotTofXmcp();
	return 0;
}