#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include <cstdlib>

using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLinearFitter.h"
#include "TF1.h"
#include "Fit/Fitter.h"
#include "Fit/FitResult.h"
#include "TFormula.h"
#include "Math/QuantFuncMathCore.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

void Pid2Tof_Full()
{
	gStyle->SetOptStat("nemri");
	gStyle->SetStatFormat(".4f");
	gStyle->SetFitFormat(".4f");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	// gStyle->SetOptLogy(1);
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);
	gROOT->ForceStyle();

	Option_t *opt = "colz scat";

	const int NCorr = 6, NSGM = 3, MINPOIN = 100;

	string sRoot, sTree;
	int i, j, k, m;

	const int NAna = 3;
	string sAna[3] = {"CfdTdc", "TacNoClk", "TacClk"};

	const int N = 200;
	string strRead;
	string strNucl[N];
	int isForCal[N] = {0}, ZOrig[N] = {0}, AOrig[N] = {0}, QOrig[N] = {0};
	double MoQOrig[N][2] = {0}, tofPivotAME[N] = {0};

	int iNucl = 0, nNucl = 0;
	ifstream fCutInfo("fCutInfo.dat");
	while (!fCutInfo.eof() && fCutInfo.peek() != EOF)
	{
		if (fCutInfo.peek() != '#')
		{
			fCutInfo >> isForCal[nNucl] >> strRead >> strNucl[nNucl] >> ZOrig[nNucl] >> AOrig[nNucl] >> QOrig[nNucl] >> strRead >> strRead >> strRead >> strRead >> MoQOrig[nNucl][0] >> MoQOrig[nNucl][1] >> tofPivotAME[nNucl];
			getline(fCutInfo, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo, strRead);
	}
	fCutInfo.close();

	TFile *fTof = new TFile("/home/kailong/ExpData/Jul2018/Pid/tof-run-150--153__270--385_Comp.root", "RECREATE");
	TTree *tTOF = new TTree("tTOF", "tree for TOF analysis");
	string sNucl = "";
	int isCal = 0, Z = 0, A = 0, Q = 0, cnts[NAna] = {0};
	bool goodElec[NAna] = {0};
	double AoQ = 0, MoQ[2] = {0};
	double tofRaw[NAna][3] = {0}, tofXmcpCorr[NAna][3], tofXsCorr[NAna][3] = {0};								   //[3]: [0] for mean TOF; [1] for sigma of mean TOF; [2] for sigma of TOF distribution
	double parTofXmcp[NAna][2][2] = {0}, parTofXs[NAna][NCorr + 1][2] = {0};									   //second [2]: [0] value; [1] error
	double tofPmtRaw[NAna][5][9][3] = {0}, tofPmtXmcpCorr[NAna][5][9][3] = {0}, tofPmtXsCorr[NAna][5][9][3] = {0}; //[5][9]: different PMT combinations: [1--4][5--8] for ToFs from PMTs#1--#4 and PMTs#5--#8; [0][1] ToF from PMTs#1#3#5#7; [1][0] Tof from PMTs#2#4#6#8
	double parTofPmtXmcp[NAna][5][9][2][2] = {0}, parTofPmtXs[NAna][5][9][NCorr + 1][2] = {0};

	int iAna = 0, iCorr = 0;
	string sDirTofPmt, sDirElec;
	for (iAna = 0; iAna < NAna; iAna++)
	{
		sDirElec = "./Graphs/" + sAna[iAna];
		if (access(sDirElec.c_str(), F_OK) != 0)
			system(("mkdir " + sDirElec).c_str());

		sDirTofPmt = sDirElec + "/TOF_PMT";
		if (access(sDirTofPmt.c_str(), F_OK) != 0)
			system(("mkdir " + sDirTofPmt).c_str());
	}

	tTOF->Branch("sNucl", &sNucl);
	tTOF->Branch("isCal", &isCal, "isCal/I");
	tTOF->Branch("Z", &Z, "Z/I");
	tTOF->Branch("A", &A, "A/I");
	tTOF->Branch("Q", &Q, "Q/I");
	tTOF->Branch("AoQ", &AoQ, "AoQ/D");
	tTOF->Branch("MoQ", MoQ, "MoQ[2]/D");
	for (iAna = 0; iAna < NAna; iAna++)
	{
		tTOF->Branch(("good_" + sAna[iAna]).c_str(), &goodElec[iAna], ("good_" + sAna[iAna] + "/O").c_str());
		tTOF->Branch(("cnts_" + sAna[iAna]).c_str(), &cnts[iAna], ("cnts_" + sAna[iAna] + "/I").c_str());
		tTOF->Branch(("tofRaw_" + sAna[iAna]).c_str(), tofRaw[iAna], ("tofRaw_" + sAna[iAna] + "[3]/D").c_str());
		tTOF->Branch(("tofXmcpCorr_" + sAna[iAna]).c_str(), tofXmcpCorr[iAna], ("tofXmcpCorr_" + sAna[iAna] + "[3]/D").c_str());
		tTOF->Branch(("parTofXmcp_" + sAna[iAna]).c_str(), parTofXmcp[iAna], ("parTofXmcp_" + sAna[iAna] + "[2][2]/D").c_str());
		tTOF->Branch(("tofXsCorr_" + sAna[iAna]).c_str(), tofXsCorr[iAna], ("tofXsCorr_" + sAna[iAna] + "[3]/D").c_str());
		tTOF->Branch(("parTofXs_" + sAna[iAna]).c_str(), parTofXs[iAna], ("parTofXs_" + sAna[iAna] + "[" + to_string(NCorr + 1) + "][2]/D").c_str());
		tTOF->Branch(("tofPmtRaw_" + sAna[iAna]).c_str(), tofPmtRaw[iAna], ("tofPmtRaw_" + sAna[iAna] + "[5][9][3]/D").c_str());
		tTOF->Branch(("tofPmtXmcpCorr_" + sAna[iAna]).c_str(), tofPmtXmcpCorr[iAna], ("tofPmtXmcpCorr_" + sAna[iAna] + "[5][9][3]/D").c_str());
		tTOF->Branch(("parTofPmtXmcp_" + sAna[iAna]).c_str(), parTofPmtXmcp[iAna], ("parTofPmtXmcp_" + sAna[iAna] + "[5][9][2][2]/D").c_str());
		tTOF->Branch(("tofPmtXsCorr_" + sAna[iAna]).c_str(), tofPmtXsCorr[iAna], ("tofPmtXsCorr_" + sAna[iAna] + "[5][9][3]/D").c_str());
		tTOF->Branch(("parTofPmtXs_" + sAna[iAna]).c_str(), parTofPmtXs[iAna], ("parTofPmtXs_" + sAna[iAna] + "[5][9][" + to_string(NCorr + 1) + "][2]/D").c_str());
	}

	string sCorr[NCorr] = {"Xmcp", "Ymcp", "Qpmt1S800", "Qpmt2S800", "Qpmt3S800", "Qpmt4S800"};
	string sCorrUnit[NCorr] = {" [mm]", " [mm]", "", "", "", ""};
	int scalBinX[NCorr] = {10, 10, 1, 1, 1, 1};
	string sPmt[5][9];
	bool isXForPmt[NAna][5][9][NCorr] = {0}, goodIDPmt[NAna][5][9] = {0};
	int nXsPmt[5][9] = {0};

	sPmt[0][1] = "Pmt1357";
	sPmt[1][0] = "Pmt2468";
	nXsPmt[0][1] = 4;
	nXsPmt[1][0] = 4;
	for (i = 1; i < 5; i++)
		for (j = 5; j < 9; j++)
		{
			sPmt[i][j] = "Pmt" + to_string(i) + to_string(j);
			nXsPmt[i][j] = 3;
		}

	for (iAna = 0; iAna < NAna; iAna++)
		for (j = 0; j < 5; j++)
			for (k = 0; k < 9; k++)
				if ((j == 0 && k == 1) || (j == 1 && k == 0) || (j > 0 && k > 4 && ((iAna == 1 && k == j + 4) || (iAna != 1))))
				{
					goodIDPmt[iAna][j][k] = true;
					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						if (iCorr < 2)
							isXForPmt[iAna][j][k][iCorr] = true;
						else if ((j > 0 && k > 4 && iCorr - 1 == j) || (j == 0 && k == 1 && iCorr % 2 == 0) || (j == 1 && k == 0 && iCorr % 2 == 1))
							isXForPmt[iAna][j][k][iCorr] = true;
					}
				}

	sRoot = "/home/kailong/ExpData/Jul2018/Pid/pid-run-150--153__270--385.root";
	sTree = "tPid";

	TFile *fPid = new TFile(sRoot.c_str());
	TTree *tPid;
	fPid->GetObject(sTree.c_str(), tPid);
	tPid->SetEstimate(-1);

	string sDraw = "ZPid:APid:QPid:xMCP:yMCP:egyPMT[0][0]:egyPMT[0][1]:egyPMT[0][2]:egyPMT[0][3]";
	for (iAna = 0; iAna < NAna; iAna++)
	{
		sDraw += ":good_" + sAna[iAna];
		sDraw += ":tof_" + sAna[iAna];
		sDraw += ":tofPmt_" + sAna[iAna] + "[0][1]:tofPmt_" + sAna[iAna] + "[1][0]";
		for (i = 1; i < 5; i++)
			for (j = 5; j < 9; j++)
				sDraw += ":tofPmt_" + sAna[iAna] + "[" + to_string(i) + "][" + to_string(j) + "]";
	}

	string sCut = "good_" + sAna[0] + "==1";
	for (iAna = 1; iAna < NAna; iAna++)
		sCut += "||good_" + sAna[iAna] + "==1";
	int nData = tPid->Draw(sDraw.c_str(), sCut.c_str(), "goff");
	double *ZVal = tPid->GetVal(0);
	double *AVal = tPid->GetVal(1);
	double *QVal = tPid->GetVal(2);
	double *xVal[NCorr];
	for (iCorr = 0; iCorr < NCorr; iCorr++)
		xVal[iCorr] = tPid->GetVal(3 + iCorr);
	double *goodVal[NAna], *tofVal[NAna], *tofPmtVal[NAna][5][9] = {0};
	for (iAna = 0; iAna < NAna; iAna++)
	{
		goodVal[iAna] = tPid->GetVal(3 + NCorr + 20 * iAna);
		tofVal[iAna] = tPid->GetVal(3 + NCorr + 20 * iAna + 1);
		tofPmtVal[iAna][0][1] = tPid->GetVal(3 + NCorr + 20 * iAna + 2);
		tofPmtVal[iAna][1][0] = tPid->GetVal(3 + NCorr + 20 * iAna + 3);
		k = 3;
		for (i = 1; i < 5; i++)
			for (j = 5; j < 9; j++)
				tofPmtVal[iAna][i][j] = tPid->GetVal(3 + NCorr + 20 * iAna + (++k));
	}

	string strFitTofX = "hyp" + to_string(NCorr);
	TLinearFitter *fitTofXs = new TLinearFitter(NCorr, strFitTofX.c_str());
	fitTofXs->SetName("fitTofXs");
	TLinearFitter *fitTofPmtXs = new TLinearFitter(NCorr, strFitTofX.c_str());
	fitTofPmtXs->SetName("fitTofPmtXs");

	string sFcnTofX = "[0]";
	k = 0;
	for (i = 0; i < NCorr; i++)
		sFcnTofX += "+[" + to_string(++k) + "]*x[" + to_string(i) + "]";
	TFormula *fcnTofXs = new TFormula("fcnTofXs", sFcnTofX.c_str());
	TFormula *fcnTofPmtXs = new TFormula("fcnTofPmtXs", sFcnTofX.c_str());

	int iData = 0, nFilt = 0, nPoin = 0;
	double xRang[NCorr][2] = {0}, tofRang[2] = {0};
	int nBinX[NCorr] = {0};
	bool **goodData = new bool *[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
		goodData[iAna] = new bool[nData];
	double xPoint[NCorr] = {0}, tofCorrVal = 0;
	double xPivot[NCorr] = {0}, tofPivot = 0, parVal[NCorr + 1] = {0};

	TFile *fDraw = new TFile("fDraw.root", "UPDATE");

	for (iNucl = 0; iNucl < nNucl; iNucl++)
	// for (iNucl = 0; iNucl < 2; iNucl++)
	{
		cout << "\n****************************************Now processing " << strNucl[iNucl] << "****************************************\n";
		sNucl = "";
		isCal = 0;
		Z = 0;
		A = 0;
		Q = 0;
		AoQ = 0;
		memset(MoQ, 0, sizeof(MoQ));
		memset(goodElec, 0, sizeof(goodElec));
		memset(cnts, 0, sizeof(cnts));
		memset(tofRaw, 0, sizeof(tofRaw));
		memset(tofXmcpCorr, 0, sizeof(tofXmcpCorr));
		memset(parTofXmcp, 0, sizeof(parTofXmcp));
		memset(tofXsCorr, 0, sizeof(tofXsCorr));
		memset(parTofXs, 0, sizeof(parTofXs));
		memset(tofPmtRaw, 0, sizeof(tofPmtRaw));
		memset(tofPmtXmcpCorr, 0, sizeof(tofPmtXmcpCorr));
		memset(parTofPmtXmcp, 0, sizeof(parTofPmtXmcp));
		memset(tofPmtXsCorr, 0, sizeof(tofPmtXsCorr));
		memset(parTofPmtXs, 0, sizeof(parTofPmtXs));

		for (i = 0; i < NAna; i++)
			for (j = 0; j < nData; j++)
				goodData[i][j] = false;
		bool goodNucl = false;

		TGraph *grTofXs[NCorr];
		for (iCorr = 0; iCorr < NCorr; iCorr++)
			grTofXs[iCorr] = new TGraph();

		for (iAna = 0; iAna < NAna; iAna++)
		{
			memset(xRang, 0, sizeof(xRang));
			memset(tofRang, 0, sizeof(tofRang));
			for (iCorr = 0; iCorr < NCorr; iCorr++)
				grTofXs[iCorr]->Set(0);
			nFilt = 0;
			for (iData = 0; iData < nData; iData++)
			{
				goodData[iAna][iData] = abs(goodVal[iAna][iData] - 1) < 1E-3 && abs(ZVal[iData] - ZOrig[iNucl]) < 1E-3 && abs(AVal[iData] - AOrig[iNucl]) < 1E-3 && abs(QVal[iData] - QOrig[iNucl]) < 1E-3 && abs(tofVal[iAna][iData]) > 0;
				if (!goodData[iAna][iData])
					continue;
				for (iCorr = 0; iCorr < NCorr; iCorr++)
				{
					goodData[iAna][iData] = abs(xVal[iCorr][iData]) > 0;
					if (!goodData[iAna][iData])
						break;
				}
				if (goodData[iAna][iData])
				{
					for (iCorr = 0; iCorr < NCorr; iCorr++)
						grTofXs[iCorr]->SetPoint(nFilt, xVal[iCorr][iData], tofVal[iAna][iData]);
					nFilt++;
				}
			}
			if (nFilt < MINPOIN)
				continue;
			for (iCorr = 0; iCorr < NCorr; iCorr++)
			{
				xRang[iCorr][0] = grTofXs[iCorr]->GetMean(1) - NSGM * grTofXs[iCorr]->GetRMS(1);
				xRang[iCorr][1] = grTofXs[iCorr]->GetMean(1) + NSGM * grTofXs[iCorr]->GetRMS(1);
			}
			tofRang[0] = grTofXs[0]->GetMean(2) - NSGM * grTofXs[0]->GetRMS(2);
			tofRang[1] = grTofXs[0]->GetMean(2) + NSGM * grTofXs[0]->GetRMS(2);

			for (iCorr = 0; iCorr < NCorr; iCorr++)
				grTofXs[iCorr]->Set(0);
			nFilt = 0;
			for (iData = 0; iData < nData; iData++)
				if (goodData[iAna][iData])
				{
					goodData[iAna][iData] = tofVal[iAna][iData] > tofRang[0] && tofVal[iAna][iData] < tofRang[1];
					if (goodData[iAna][iData])
						for (iCorr = 0; iCorr < NCorr; iCorr++)
						{
							goodData[iAna][iData] = xVal[iCorr][iData] > xRang[iCorr][0] && xVal[iCorr][iData] < xRang[iCorr][1];
							if (!goodData[iAna][iData])
								break;
						}
					if (goodData[iAna][iData])
					{
						for (iCorr = 0; iCorr < NCorr; iCorr++)
							grTofXs[iCorr]->SetPoint(nFilt, xVal[iCorr][iData], tofVal[iAna][iData]);
						nFilt++;
					}
				}
			if (nFilt < MINPOIN)
				continue;
			else
				goodElec[iAna] = true;
		}
		goodNucl = goodElec[0];
		for (iAna = 1; iAna < NAna; iAna++)
			goodNucl = goodNucl || goodElec[iAna];
		if (goodNucl)
		{
			sNucl = strNucl[iNucl];
			isCal = isForCal[iNucl];
			Z = ZOrig[iNucl];
			A = AOrig[iNucl];
			Q = QOrig[iNucl];
			AoQ = 1.0 * A / Q;
			MoQ[0] = MoQOrig[iNucl][0];
			MoQ[1] = MoQOrig[iNucl][1];

			for (iAna = 0; iAna < NAna; iAna++)
				if (goodElec[iAna])
				{
					memset(xRang, 0, sizeof(xRang));
					memset(nBinX, 0, sizeof(nBinX));
					memset(tofRang, 0, sizeof(tofRang));
					memset(xPoint, 0, sizeof(xPoint));
					memset(xPivot, 0, sizeof(xPivot));
					memset(parVal, 0, sizeof(parVal));
					tofPivot = 0;
					tofCorrVal = 0;
					TGraph *grTofXmcpCorrXs[NCorr], *grTofXsCorrXs[NCorr];
					TGraph *grTofPmtXs[5][9][NCorr], *grTofPmtXmcpCorrXs[5][9][NCorr], *grTofPmtXsCorrXs[5][9][NCorr];
					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						grTofXs[iCorr]->Set(0);
						grTofXmcpCorrXs[iCorr] = new TGraph();
						grTofXsCorrXs[iCorr] = new TGraph();
						for (j = 0; j < 5; j++)
							for (k = 0; k < 9; k++)
								if (goodIDPmt[iAna][j][k])
								{
									grTofPmtXs[j][k][iCorr] = new TGraph();
									grTofPmtXmcpCorrXs[j][k][iCorr] = new TGraph();
									grTofPmtXsCorrXs[j][k][iCorr] = new TGraph();
								}
					}
					for (iData = 0; iData < nData; iData++)
						if (goodData[iAna][iData])
						{
							for (iCorr = 0; iCorr < NCorr; iCorr++)
							{
								grTofXs[iCorr]->SetPoint(cnts[iAna], xVal[iCorr][iData], tofVal[iAna][iData]);
								for (j = 0; j < 5; j++)
									for (k = 0; k < 9; k++)
										if (isXForPmt[iAna][j][k][iCorr])
											grTofPmtXs[j][k][iCorr]->SetPoint(cnts[iAna], xVal[iCorr][iData], tofPmtVal[iAna][j][k][iData]);
							}
							cnts[iAna]++;
						}
					tofRaw[iAna][0] = grTofXs[0]->GetMean(2);
					tofRaw[iAna][2] = grTofXs[0]->GetRMS(2);
					tofRaw[iAna][1] = tofRaw[iAna][2] / sqrt(cnts[iAna]);

					for (j = 0; j < 5; j++)
						for (k = 0; k < 9; k++)
							if (goodIDPmt[iAna][j][k])
							{
								tofPmtRaw[iAna][j][k][0] = grTofPmtXs[j][k][0]->GetMean(2);
								tofPmtRaw[iAna][j][k][2] = grTofPmtXs[j][k][0]->GetRMS(2);
								tofPmtRaw[iAna][j][k][1] = tofPmtRaw[iAna][j][k][2] / sqrt(cnts[iAna]);
							}

					TFitResultPtr fitTofXmcp = grTofXs[0]->Fit("pol1", "SQ0");
					for (i = 0; i < 2; i++)
					{
						parTofXmcp[iAna][i][0] = fitTofXmcp->Parameter(i);
						parTofXmcp[iAna][i][1] = fitTofXmcp->ParError(i);
					}
					xPivot[0] = 0;
					tofPivot = parTofXmcp[iAna][0][0] + parTofXmcp[iAna][1][0] * xPivot[0];
					nPoin = 0;
					for (iData = 0; iData < nData; iData++)
						if (goodData[iAna][iData])
						{
							tofCorrVal = tofVal[iAna][iData] + tofPivot - parTofXmcp[iAna][0][0] - parTofXmcp[iAna][1][0] * xVal[0][iData];
							for (iCorr = 0; iCorr < NCorr; iCorr++)
								grTofXmcpCorrXs[iCorr]->SetPoint(nPoin, xVal[iCorr][iData], tofCorrVal);
							nPoin++;
						}
					tofXmcpCorr[iAna][0] = grTofXmcpCorrXs[0]->GetMean(2);
					tofXmcpCorr[iAna][2] = grTofXmcpCorrXs[0]->GetRMS(2);
					tofXmcpCorr[iAna][1] = tofXmcpCorr[iAna][2] / sqrt(nPoin);

					for (j = 0; j < 5; j++)
						for (k = 0; k < 9; k++)
							if (isXForPmt[iAna][j][k][0])
							{
								TFitResultPtr fitTofPmtXmcp = grTofPmtXs[j][k][0]->Fit("pol1", "SQ0");
								for (i = 0; i < 2; i++)
								{
									parTofPmtXmcp[iAna][j][k][i][0] = fitTofPmtXmcp->Parameter(i);
									parTofPmtXmcp[iAna][j][k][i][1] = fitTofPmtXmcp->ParError(i);
								}
								xPivot[0] = 0;
								tofPivot = parTofPmtXmcp[iAna][j][k][0][0] + parTofPmtXmcp[iAna][j][k][1][0] * xPivot[0];
								nPoin = 0;
								for (iData = 0; iData < nData; iData++)
									if (goodData[iAna][iData])
									{
										tofCorrVal = tofPmtVal[iAna][j][k][iData] + tofPivot - parTofPmtXmcp[iAna][j][k][0][0] - parTofPmtXmcp[iAna][j][k][1][0] * xVal[0][iData];
										for (iCorr = 0; iCorr < NCorr; iCorr++)
											if (isXForPmt[iAna][j][k][iCorr])
												grTofPmtXmcpCorrXs[j][k][iCorr]->SetPoint(nPoin, xVal[iCorr][iData], tofCorrVal);
										nPoin++;
									}
								for (iCorr = 0; iCorr < NCorr; iCorr++)
									if (isXForPmt[iAna][j][k][iCorr])
									{
										tofPmtXmcpCorr[iAna][j][k][0] = grTofPmtXmcpCorrXs[j][k][0]->GetMean(2);
										tofPmtXmcpCorr[iAna][j][k][2] = grTofPmtXmcpCorrXs[j][k][0]->GetRMS(2);
										tofPmtXmcpCorr[iAna][j][k][1] = tofPmtXmcpCorr[iAna][j][k][2] / sqrt(nPoin);
									}
							}

					fitTofXs->ClearPoints();
					for (iData = 0; iData < nData; iData++)
						if (goodData[iAna][iData])
						{
							memset(xPoint, 0, sizeof(xPoint));
							for (iCorr = 0; iCorr < NCorr; iCorr++)
								xPoint[iCorr] = xVal[iCorr][iData];
							fitTofXs->AddPoint(xPoint, tofVal[iAna][iData]);
						}
					fitTofXs->Eval();
					for (j = 0; j < NCorr + 1; j++)
					{
						parTofXs[iAna][j][0] = fitTofXs->GetParameter(j);
						parTofXs[iAna][j][1] = fitTofXs->GetParError(j);
						parVal[j] = parTofXs[iAna][j][0];
					}
					memset(xPivot, 0, sizeof(xPivot));
					for (iCorr = 2; iCorr < NCorr; iCorr++)
						xPivot[iCorr] = grTofXs[iCorr]->GetMean(1);

					// // Use the same tofPivot for different PMT combinations
					// tofPivot = fcnTofXs->EvalPar(xPivot, parVal);
					// printf("tofPivot=%f, tofPivotAME=%f, Delta=%f\n", tofPivot, tofPivotAME[iNucl], tofPivot-tofPivotAME[iNucl]);
					tofPivot = tofPivotAME[iNucl];

					nPoin = 0;
					for (iData = 0; iData < nData; iData++)
						if (goodData[iAna][iData])
						{
							memset(xPoint, 0, sizeof(xPoint));
							for (iCorr = 0; iCorr < NCorr; iCorr++)
								xPoint[iCorr] = xVal[iCorr][iData];
							tofCorrVal = tofVal[iAna][iData] + tofPivot - fcnTofXs->EvalPar(xPoint, parVal);

							for (iCorr = 0; iCorr < NCorr; iCorr++)
								grTofXsCorrXs[iCorr]->SetPoint(nPoin, xVal[iCorr][iData], tofCorrVal);
							nPoin++;
						}
					tofXsCorr[iAna][0] = grTofXsCorrXs[0]->GetMean(2);
					tofXsCorr[iAna][2] = grTofXsCorrXs[0]->GetRMS(2);
					tofXsCorr[iAna][1] = tofXsCorr[iAna][2] / sqrt(nPoin);

					for (j = 0; j < 5; j++)
						for (k = 0; k < 9; k++)
							if (goodIDPmt[iAna][j][k])
							{
								fitTofPmtXs->ClearPoints();
								for (iData = 0; iData < nData; iData++)
									if (goodData[iAna][iData])
									{
										memset(xPoint, 0, sizeof(xPoint));
										for (iCorr = 0; iCorr < NCorr; iCorr++)
											xPoint[iCorr] = xVal[iCorr][iData];
										fitTofPmtXs->AddPoint(xPoint, tofPmtVal[iAna][j][k][iData]);
									}
								// cout<<"before: fitTofPmtXs->GetNumberFreeParameters() = "<<fitTofPmtXs->GetNumberFreeParameters()<<endl;
								for (iCorr = 0; iCorr < NCorr; iCorr++)
									if (!isXForPmt[iAna][j][k][iCorr])
									{
										fitTofPmtXs->FixParameter(iCorr + 1, 0);
										// cout<<fitTofPmtXs->GetParameter(iCorr+1)<<endl;
									}
								// cout<<"after: fitTofPmtXs->GetNumberFreeParameters() = "<<fitTofPmtXs->GetNumberFreeParameters()<<endl;
								fitTofPmtXs->Eval();
								memset(parVal, 0, sizeof(parVal));
								for (i = 0; i < NCorr + 1; i++)
								{
									parTofPmtXs[iAna][j][k][i][0] = fitTofPmtXs->GetParameter(i);
									parTofPmtXs[iAna][j][k][i][1] = fitTofPmtXs->GetParError(i);
									parVal[i] = parTofPmtXs[iAna][j][k][i][0];
								}
								for (iCorr = 0; iCorr < NCorr; iCorr++)
									if (!isXForPmt[iAna][j][k][iCorr])
									{
										fitTofPmtXs->ReleaseParameter(iCorr + 1);
									}
								memset(xPivot, 0, sizeof(xPivot));
								for (iCorr = 2; iCorr < NCorr; iCorr++)
									xPivot[iCorr] = grTofXs[iCorr]->GetMean(1);

								nPoin = 0;
								for (iData = 0; iData < nData; iData++)
									if (goodData[iAna][iData])
									{
										memset(xPoint, 0, sizeof(xPoint));
										for (iCorr = 0; iCorr < NCorr; iCorr++)
											xPoint[iCorr] = xVal[iCorr][iData];
										tofCorrVal = tofPmtVal[iAna][j][k][iData] + tofPivot - fcnTofPmtXs->EvalPar(xPoint, parVal);

										for (iCorr = 0; iCorr < NCorr; iCorr++)
											if (isXForPmt[iAna][j][k][iCorr])
												grTofPmtXsCorrXs[j][k][iCorr]->SetPoint(nPoin, xVal[iCorr][iData], tofCorrVal);
										nPoin++;
									}
								tofPmtXsCorr[iAna][j][k][0] = grTofPmtXsCorrXs[j][k][0]->GetMean(2);
								tofPmtXsCorr[iAna][j][k][2] = grTofPmtXsCorrXs[j][k][0]->GetRMS(2);
								tofPmtXsCorr[iAna][j][k][1] = tofPmtXsCorr[iAna][j][k][2] / sqrt(nPoin);
							}
					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						xRang[iCorr][0] = grTofXs[iCorr]->GetMean(1) - 2 * NSGM * grTofXs[iCorr]->GetRMS(1);
						xRang[iCorr][1] = grTofXs[iCorr]->GetMean(1) + 2 * NSGM * grTofXs[iCorr]->GetRMS(1);
						nBinX[iCorr] = (int)((xRang[iCorr][1] - xRang[iCorr][0]) * scalBinX[iCorr]);
					}
					tofRang[0] = tofRaw[iAna][0] - 2 * NSGM * tofRaw[iAna][2];
					tofRang[1] = tofRaw[iAna][0] + 2 * NSGM * tofRaw[iAna][2];
					int nBinTof = (int)((tofRang[1] - tofRang[0]) * 100);

					double tofXsCorrRang[2] = {0};
					tofXsCorrRang[0] = tofPivot - 2 * NSGM * tofRaw[iAna][2];
					tofXsCorrRang[1] = tofPivot + 2 * NSGM * tofRaw[iAna][2];

					double tofPmtRang[5][9][2] = {0};
					int nBinTofPmt[5][9] = {0};
					for (j = 0; j < 5; j++)
						for (k = 0; k < 9; k++)
							if (goodIDPmt[iAna][j][k])
							{
								tofPmtRang[j][k][0] = tofPmtRaw[iAna][j][k][0] - 2 * NSGM * tofPmtRaw[iAna][j][k][2];
								tofPmtRang[j][k][1] = tofPmtRaw[iAna][j][k][0] + 2 * NSGM * tofPmtRaw[iAna][j][k][2];
								nBinTofPmt[j][k] = (int)((tofPmtRang[j][k][1] - tofPmtRang[j][k][0]) * 100);
							}
					TH2F *hTofXs[NCorr], *hTofXmcpCorrXs[NCorr], *hTofXsCorrXs[NCorr];
					TH2F *hTofPmtXs[5][9][NCorr], *hTofPmtXmcpCorrXs[5][9][NCorr], *hTofPmtXsCorrXs[5][9][NCorr];
					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						hTofXs[iCorr] = new TH2F(("hTofRawVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTofRawVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{raw} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTof, tofRang[0], tofRang[1]);

						hTofXmcpCorrXs[iCorr] = new TH2F(("hTofXmcpCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTofXmcpCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{XmcpCorr} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTof, tofRang[0], tofRang[1]);

						hTofXsCorrXs[iCorr] = new TH2F(("hTofXsCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTofXsCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{XsCorr} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTof, tofXsCorrRang[0], tofXsCorrRang[1]);

						for (j = 0; j < 5; j++)
							for (k = 0; k < 9; k++)
								if (isXForPmt[iAna][j][k][iCorr])
								{
									hTofPmtXs[j][k][iCorr] = new TH2F(("hTof" + sPmt[j][k] + "RawVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTof" + sPmt[j][k] + "RawVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{raw}^{" + sPmt[j][k] + "} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTofPmt[j][k], tofPmtRang[j][k][0], tofPmtRang[j][k][1]);

									hTofPmtXmcpCorrXs[j][k][iCorr] = new TH2F(("hTof" + sPmt[j][k] + "XmcpCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTof" + sPmt[j][k] + "XmcpCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{XmcpCorr}^{" + sPmt[j][k] + "} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTofPmt[j][k], tofPmtRang[j][k][0], tofPmtRang[j][k][1]);

									hTofPmtXsCorrXs[j][k][iCorr] = new TH2F(("hTof" + sPmt[j][k] + "XsCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTof" + sPmt[j][k] + "XsCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{XsCorr}^{" + sPmt[j][k] + "} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTof, tofXsCorrRang[0], tofXsCorrRang[1]);
								}
					}
					for (iData = 0; iData < cnts[iAna]; iData++)
						for (iCorr = 0; iCorr < NCorr; iCorr++)
						{
							hTofXs[iCorr]->Fill(grTofXs[iCorr]->GetX()[iData], grTofXs[iCorr]->GetY()[iData]);
							hTofXmcpCorrXs[iCorr]->Fill(grTofXmcpCorrXs[iCorr]->GetX()[iData], grTofXmcpCorrXs[iCorr]->GetY()[iData]);
							hTofXsCorrXs[iCorr]->Fill(grTofXsCorrXs[iCorr]->GetX()[iData], grTofXsCorrXs[iCorr]->GetY()[iData]);

							for (j = 0; j < 5; j++)
								for (k = 0; k < 9; k++)
									if (isXForPmt[iAna][j][k][iCorr])
									{
										hTofPmtXs[j][k][iCorr]->Fill(grTofPmtXs[j][k][iCorr]->GetX()[iData], grTofPmtXs[j][k][iCorr]->GetY()[iData]);
										hTofPmtXmcpCorrXs[j][k][iCorr]->Fill(grTofPmtXmcpCorrXs[j][k][iCorr]->GetX()[iData], grTofPmtXmcpCorrXs[j][k][iCorr]->GetY()[iData]);
										hTofPmtXsCorrXs[j][k][iCorr]->Fill(grTofPmtXsCorrXs[j][k][iCorr]->GetX()[iData], grTofPmtXsCorrXs[j][k][iCorr]->GetY()[iData]);
									}
						}
					TCanvas *cTofXs = new TCanvas(("cTofXs_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("cTofXs_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), 400 * NCorr, 900);
					cTofXs->Divide(NCorr, 3);
					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						cTofXs->cd(iCorr + 1);
						hTofXs[iCorr]->Draw(opt);

						cTofXs->cd(iCorr + NCorr + 1);
						hTofXmcpCorrXs[iCorr]->Draw(opt);

						cTofXs->cd(iCorr + 2 * NCorr + 1);
						hTofXsCorrXs[iCorr]->Draw(opt);
					}
					cTofXs->Write("", TObject::kOverwrite);
					cTofXs->SaveAs((sDirElec + "/cTofXs_" + sAna[iAna] + "_" + strNucl[iNucl] + ".png").c_str());
					delete cTofXs;

					TCanvas *cTof = new TCanvas(("cTof_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("cTof_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
					cTof->Divide(2, 2);
					cTof->cd(1);
					TH1D *hTof = hTofXs[0]->ProjectionY(("hTofRaw__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
					hTof->SetTitle(("hTofRaw__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{raw} [ns];Counts/0.01ns").c_str());
					hTof->Draw();
					hTof->Fit("gaus", "Q");

					cTof->cd(2);
					TH1D *hTofXmcpCorr = hTofXmcpCorrXs[0]->ProjectionY(("hTofXmcpCorr__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
					hTofXmcpCorr->SetTitle(("hTofXmcpCorr__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{XmcpCorr} [ns];Counts/0.01ns").c_str());
					hTofXmcpCorr->Draw();
					hTofXmcpCorr->Fit("gaus", "Q");

					cTof->cd(3);
					TH1D *hTofXsCorr = hTofXsCorrXs[0]->ProjectionY(("hTofXsCorr__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
					hTofXsCorr->SetTitle(("hTofXsCorr__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{XsCorr} [ns];Counts/0.01ns").c_str());
					hTofXsCorr->Draw();
					hTofXsCorr->Fit("gaus", "Q");

					cTof->Write("", TObject::kOverwrite);
					cTof->SaveAs((sDirElec + "/cTof_" + sAna[iAna] + "_" + strNucl[iNucl] + ".png").c_str());
					delete cTof;

					TCanvas *cTofPmtXs[5][9], *cTofPmt[5][9];
					for (j = 0; j < 5; j++)
						for (k = 0; k < 9; k++)
							if (goodIDPmt[iAna][j][k])
							{
								cTofPmtXs[j][k] = new TCanvas(("cTof" + sPmt[j][k] + "Xs_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("cTof" + sPmt[j][k] + "Xs_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), 400 * nXsPmt[j][k], 900);
								cTofPmtXs[j][k]->Divide(nXsPmt[j][k], 3);
								m = 1;
								for (iCorr = 0; iCorr < NCorr; iCorr++)
									if (isXForPmt[iAna][j][k][iCorr])
									{
										cTofPmtXs[j][k]->cd(m);
										hTofPmtXs[j][k][iCorr]->Draw(opt);

										cTofPmtXs[j][k]->cd(m + nXsPmt[j][k]);
										hTofPmtXmcpCorrXs[j][k][iCorr]->Draw(opt);

										cTofPmtXs[j][k]->cd(m + 2 * nXsPmt[j][k]);
										hTofPmtXsCorrXs[j][k][iCorr]->Draw(opt);

										m++;
									}
								cTofPmtXs[j][k]->Write("", TObject::kOverwrite);
								cTofPmtXs[j][k]->SaveAs((sDirTofPmt + "/cTof" + sPmt[j][k] + "Xs_" + sAna[iAna] + "_" + strNucl[iNucl] + ".png").c_str());
								delete cTofPmtXs[j][k];

								cTofPmt[j][k] = new TCanvas(("cTof" + sPmt[j][k] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("cTof" + sPmt[j][k] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
								cTofPmt[j][k]->Divide(2, 2);

								cTofPmt[j][k]->cd(1);
								TH1D *hTofPmt = hTofPmtXs[j][k][0]->ProjectionY(("hTof" + sPmt[j][k] + "Raw__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
								hTofPmt->SetTitle(("hTof" + sPmt[j][k] + "Raw__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{raw}^{" + sPmt[j][k] + "} [ns];Counts/0.01ns").c_str());
								hTofPmt->Draw();
								hTofPmt->Fit("gaus", "Q");

								cTofPmt[j][k]->cd(2);
								TH1D *hTofPmtXmcpCorr = hTofPmtXmcpCorrXs[j][k][0]->ProjectionY(("hTof" + sPmt[j][k] + "XmcpCorr__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
								hTofPmtXmcpCorr->SetTitle(("hTof" + sPmt[j][k] + "XmcpCorr__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{XmcpCorr}^{" + sPmt[j][k] + "} [ns];Counts/0.01ns").c_str());
								hTofPmtXmcpCorr->Draw();
								hTofPmtXmcpCorr->Fit("gaus", "Q");

								cTofPmt[j][k]->cd(3);
								TH1D *hTofPmtXsCorr = hTofPmtXsCorrXs[j][k][0]->ProjectionY(("hTof" + sPmt[j][k] + "XsCorr__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
								hTofPmtXsCorr->SetTitle(("hTof" + sPmt[j][k] + "XsCorr__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{XsCorr}^{" + sPmt[j][k] + "} [ns];Counts/0.01ns").c_str());
								hTofPmtXsCorr->Draw();
								hTofPmtXsCorr->Fit("gaus", "Q");

								cTofPmt[j][k]->Write("", TObject::kOverwrite);
								cTofPmt[j][k]->SaveAs((sDirTofPmt + "/cTof" + sPmt[j][k] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ".png").c_str());
								delete cTofPmt[j][k];
							}
					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						delete hTofXs[iCorr];
						delete hTofXmcpCorrXs[iCorr];
						delete hTofXsCorrXs[iCorr];

						for (j = 0; j < 5; j++)
							for (k = 0; k < 9; k++)
								if (isXForPmt[iAna][j][k][iCorr])
								{
									delete hTofPmtXs[j][k][iCorr];
									delete hTofPmtXmcpCorrXs[j][k][iCorr];
									delete hTofPmtXsCorrXs[j][k][iCorr];
								}
					}
					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						delete grTofXmcpCorrXs[iCorr];
						delete grTofXsCorrXs[iCorr];
						for (j = 0; j < 5; j++)
							for (k = 0; k < 9; k++)
								if (goodIDPmt[iAna][j][k])
								{
									delete grTofPmtXs[j][k][iCorr];
									delete grTofPmtXmcpCorrXs[j][k][iCorr];
									delete grTofPmtXsCorrXs[j][k][iCorr];
								}
					}
				}
		}
		tTOF->Fill();
		for (iCorr = 0; iCorr < NCorr; iCorr++)
			delete grTofXs[iCorr];
	}
	fTof->cd();
	tTOF->Write();

	delete fDraw;
	for (iAna = 0; iAna < NAna; iAna++)
		delete[] goodData[iAna];
	delete[] goodData;
	delete fcnTofPmtXs;
	delete fcnTofXs;
	delete fitTofPmtXs;
	delete fitTofXs;
	delete fPid;
	delete tTOF;
	delete fTof;
}

int main()
{
	Pid2Tof_Full();
	return 0;
}