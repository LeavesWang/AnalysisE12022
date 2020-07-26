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
#include "TH1D.h"
#include "TTree.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TBox.h"
#include "TEllipse.h"
#include "TPad.h"

void PlotCutZAoQ()
{
	gStyle->SetOptStat("nemri");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);
	gStyle->SetPalette(1);

	Option_t *opt = "cont colz scat";
	const int N = 200, NSGM = 5, MINPOIN = 100;
	string sCan, sDraw, sCut, sAna;
	int i;

	int iNucl, nNucl;
	string strNucl[N];
	int ZOrig[N] = {0}, AOrig[N] = {0}, QOrig[N] = {0};
	double ZCut[N][2] = {0}, AoQCut[N][2] = {0};
	string strRead;
	ifstream fCutInfo("fCutInfo.dat");
	nNucl = 0;
	while (!fCutInfo.eof() && fCutInfo.peek() != EOF)
	{
		if (fCutInfo.peek() != '#')
		{
			fCutInfo >> strRead >> strRead >> strNucl[nNucl] >> ZOrig[nNucl] >> AOrig[nNucl] >> QOrig[nNucl] >> ZCut[nNucl][0] >> ZCut[nNucl][1] >> AoQCut[nNucl][0] >> AoQCut[nNucl][1];
			getline(fCutInfo, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo, strRead);
	}
	fCutInfo.close();
	double wZ[N] = {0}, wAoQ[N] = {0}, cenZ[N] = {0}, cenAoQ[N] = {0};
	for (iNucl = 0; iNucl < nNucl; iNucl++)
	{
		wZ[iNucl] = (ZCut[iNucl][1] - ZCut[iNucl][0]) / 2;
		wAoQ[iNucl] = (AoQCut[iNucl][1] - AoQCut[iNucl][0]) / 2;
		cenAoQ[iNucl] = (AoQCut[iNucl][1] + AoQCut[iNucl][0]) / 2;
		cenZ[iNucl] = (ZCut[iNucl][1] + ZCut[iNucl][0]) / 2;
	}

	TFile *fCut = new TFile("fCut.root", "UPDATE");

	string strCut = "";
	string strDir = "/home/kailong/ExpData/Jul2018/Pid/Graphs/PID/PID" + strCut;
	if (access(strDir.c_str(), F_OK) != 0)
		system(("mkdir " + strDir).c_str());
	string strGif = strDir + "/" + "cZAoQCut" + strCut + ".gif";
	if (access(strGif.c_str(), F_OK) == 0)
		system(("rm -f " + strGif).c_str());

	TGraph *grZAoQPid;
	fCut->GetObject(("grZAoQPid" + strCut).c_str(), grZAoQPid);
	double *AoQPid = grZAoQPid->GetX();
	double *ZPid = grZAoQPid->GetY();
	int nData = grZAoQPid->GetN();
	TH2F *hZAoQPid;
	fCut->GetObject(("hZAoQPid" + strCut).c_str(), hZAoQPid);

	TCanvas *cZAoQCutFull = new TCanvas("cZAoQCut_Full", "cZAoQCut_Full", 2400, 900);
	cZAoQCutFull->Divide(2, 1);

	cZAoQCutFull->cd(1);
	hZAoQPid->SetAxisRange(2.4, 2.75, "X");
	hZAoQPid->SetAxisRange(29, 49, "Y");
	hZAoQPid->SetStats(kFALSE);
	hZAoQPid->DrawClone("colz scat");

	cZAoQCutFull->cd(2);
	TH2F *hZAoQPidCuts = new TH2F(("hZAoQPidCuts" + strCut).c_str(), ("hZAoQPidCuts" + strCut + ";A/Q;Z").c_str(), 1750, 2.4, 2.75, 1000, 29, 49);

	for (iNucl = 0; iNucl < nNucl; iNucl++)
	//  for(iNucl=0; iNucl<3; iNucl++)
	{
		TH2F *hZAoQCut = new TH2F(("hZAoQ_PidCut_" + strNucl[iNucl]).c_str(), ("hZAoQ_PidCut_" + strNucl[iNucl] + ";A/Q;Z").c_str(), 1750, 2.4, 2.75, 1200, 28, 52);

		for (i = 0; i < nData; i++)
			if (pow((ZPid[i] - cenZ[iNucl]) / wZ[iNucl], 2) + pow((AoQPid[i] - cenAoQ[iNucl]) / wAoQ[iNucl], 2) <= 1) //for circle cut
			// if(ZPid[i]>ZCut[iNucl][0]&&ZPid[i]<ZCut[iNucl][1]&&AoQPid[i]>AoQCut[iNucl][0]&&AoQPid[i]<AoQCut[iNucl][1])  //For box cut
			{
				hZAoQCut->Fill(AoQPid[i], ZPid[i]);
				hZAoQPidCuts->Fill(AoQPid[i], ZPid[i]);
			}

		if (hZAoQCut->Integral() < MINPOIN)
			continue;

		sCan = "cZAoQCut_" + strNucl[iNucl] + strCut;
		TCanvas *cZAoQCut = new TCanvas(sCan.c_str(), sCan.c_str());
		cZAoQCut->Divide(2, 2);

		cZAoQCut->cd(2);
		hZAoQCut->Draw(opt);
		double meanAoQ = hZAoQCut->GetMean(1);
		double rmsAoQ = hZAoQCut->GetRMS(1);
		double meanZ = hZAoQCut->GetMean(2);
		double rmsZ = hZAoQCut->GetRMS(2);
		hZAoQCut->SetAxisRange(meanAoQ - NSGM * rmsAoQ, meanAoQ + NSGM * rmsAoQ, "X");
		hZAoQCut->SetAxisRange(meanZ - NSGM * rmsZ, meanZ + NSGM * rmsZ, "Y");

		cZAoQCut->cd(3);
		TH1D *hZCut = hZAoQCut->ProjectionY();
		hZCut->SetNameTitle(("hZ_" + strNucl[iNucl]).c_str(), ("hZ_" + strNucl[iNucl] + ";Z;counts/0.02").c_str());
		hZCut->Draw();
		hZCut->Fit("gaus", "Q");

		cZAoQCut->cd(4);
		TH1D *hAoQCut = hZAoQCut->ProjectionX();
		hAoQCut->SetNameTitle(("hAoQ_" + strNucl[iNucl]).c_str(), ("hAoQ_" + strNucl[iNucl] + ";A/Q;counts/0.0002").c_str());
		hAoQCut->Draw();
		hAoQCut->Fit("gaus", "Q");

		cZAoQCut->cd(1);
		hZAoQPid->Draw("colz");
		hZAoQPid->SetAxisRange(AoQCut[iNucl][0] - 10 * wAoQ[iNucl], AoQCut[iNucl][1] + 10 * wAoQ[iNucl], "X");
		hZAoQPid->SetAxisRange(ZCut[iNucl][0] - 10 * wZ[iNucl], ZCut[iNucl][1] + 10 * wZ[iNucl], "Y");

		// For box cut
		// TBox *cutBox=new TBox(AoQCut[iNucl][0], ZCut[iNucl][0], AoQCut[iNucl][1], ZCut[iNucl][1]);
		// cutBox->SetFillStyle(0);
		// cutBox->SetLineColor(kRed);
		// cutBox->SetLineWidth(2);
		// cutBox->Draw("same");

		// For circle cut
		TEllipse *cutElps = new TEllipse(cenAoQ[iNucl], cenZ[iNucl], wAoQ[iNucl], wZ[iNucl]);
		cutElps->SetFillStyle(0);
		cutElps->SetLineColor(kRed);
		cutElps->SetLineWidth(2);
		cutElps->DrawClone("same");
		
		cZAoQCutFull->cd(1);
		cutElps->SetLineWidth(1);
		cutElps->DrawClone("same");

		cZAoQCut->Write("", TObject::kOverwrite);
		cZAoQCut->SaveAs((strDir + "/" + sCan + ".png").c_str());
		cZAoQCut->SaveAs((strGif + "+40").c_str());

		delete cutElps;
		delete cZAoQCut;
		delete hZAoQCut;
	}
	cZAoQCutFull->cd(2);
	hZAoQPidCuts->Draw("colz scat");
	hZAoQPidCuts->SetStats(kFALSE);

	cZAoQCutFull->Write("", TObject::kOverwrite);
	cZAoQCutFull->SaveAs((strDir + "/cZAoQCut_Full.png").c_str());

	delete hZAoQPidCuts;
	delete cZAoQCutFull;
	fCut->Close();
	delete fCut;
}

void StandaloneApplication(int argc, char **argv)
{
	PlotCutZAoQ();
}

int main(int argc, char **argv)
{
	TApplication app("ROOT Application", &argc, argv);
	StandaloneApplication(app.Argc(), app.Argv());
	// app.Run();
	return 0;
}