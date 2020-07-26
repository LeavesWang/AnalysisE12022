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
#include "TMultiGraph.h"
#include "TMath.h"
#include "TFormula.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

void Pid2TofGlob()
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
	gStyle->SetPalette(1);
	gROOT->ForceStyle();

	Option_t *opt = "colz scat";

	const int NCorr = 2, NSGM = 3, MINPOIN = 100;
	const int NPar = NCorr + 1;
	string strPar[NCorr] = {"parXmcp", "parXcrdc1"};

	int i, j, k;
	int iPar = 0;

	const int NAna = 1;
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

	int iAna = 0, iCorr = 0;

	TFile *fTof = new TFile("/home/kailong/ExpData/Jul2018/Pid/tofGlob-run-150--153__270--385.root", "RECREATE");
	TTree *tTOF = new TTree("tTOF", "tree for TOF analysis");
	TTree *tTOFDetail = new TTree("tTOFDetail", "tree for TOF analysis with details");
	string sNucl = "";
	int isCal = 0, Z = 0, A = 0, Q = 0, cnts[NAna] = {0};
	bool goodElec[NAna] = {0};
	double AoQ = 0, MoQ[2] = {0};
	double tofRaw[NAna][3] = {0}, tofXsCorr[NAna][3] = {0}; //[3]: [0] for mean TOF; [1] for sigma of mean TOF; [2] for sigma of TOF distribution
	double parTofXs[NAna][NPar][2] = {0};					//second [2]: [0] value; [1] error

	double run = 0;
	double tof0[NAna] = {0}, tof2XsCor[NAna] = {0};
	double xMcp = 0, yMcp = 0;
	double xPlaQS800 = 0, yPlaQS800 = 0, xCrdc[2] = {0}, yCrdc[2] = {0};

	string sDirElec, strGif[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
	{
		sDirElec = "./Graphs/" + sAna[iAna];
		if (access(sDirElec.c_str(), F_OK) != 0)
			system(("mkdir " + sDirElec).c_str());
		strGif[iAna] = sDirElec + "/cTofGlobXs_" + sAna[iAna] + "_XmcpXcrdc1.gif";
		if (access(strGif[iAna].c_str(), F_OK) == 0)
			system(("rm -f " + strGif[iAna]).c_str());
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
		tTOF->Branch(("tofXsCorr_" + sAna[iAna]).c_str(), tofXsCorr[iAna], ("tofXsCorr_" + sAna[iAna] + "[3]/D").c_str());
		tTOF->Branch(("parTofXs_" + sAna[iAna]).c_str(), parTofXs[iAna], ("parTofXs_" + sAna[iAna] + "[" + to_string(NPar) + "][2]/D").c_str());
	}
	tTOFDetail->Branch("run", &run, "run/D");
	tTOFDetail->Branch("sNucl", &sNucl);
	tTOFDetail->Branch("isCal", &isCal, "isCal/I");
	tTOFDetail->Branch("Z", &Z, "Z/I");
	tTOFDetail->Branch("A", &A, "A/I");
	tTOFDetail->Branch("Q", &Q, "Q/I");
	tTOFDetail->Branch("AoQ", &AoQ, "AoQ/D");
	tTOFDetail->Branch("MoQ", MoQ, "MoQ[2]/D");
	for (iAna = 0; iAna < NAna; iAna++)
	{
		tTOFDetail->Branch(("good_" + sAna[iAna]).c_str(), &goodElec[iAna], ("good_" + sAna[iAna] + "/O").c_str());
		tTOFDetail->Branch(("tof0_" + sAna[iAna]).c_str(), &tof0[iAna], ("tof0_" + sAna[iAna] + "/D").c_str());
		tTOFDetail->Branch(("tof2XsCor_" + sAna[iAna]).c_str(), &tof2XsCor[iAna], ("tof2XsCor_" + sAna[iAna] + "/D").c_str());
	}
	tTOFDetail->Branch("xMcp", &xMcp, "xMcp/D");
	tTOFDetail->Branch("yMcp", &yMcp, "yMcp/D");
	tTOFDetail->Branch("xPlaQS800", &xPlaQS800, "xPlaQS800/D");
	tTOFDetail->Branch("yPlaQS800", &yPlaQS800, "yPlaQS800/D");
	tTOFDetail->Branch("xCrdc", xCrdc, "xCrdc[2]/D");
	tTOFDetail->Branch("yCrdc", yCrdc, "yCrdc[2]/D");

	string sCorr[NCorr] = {"Xmcp", "xCrdc1"};
	string sCorrUnit[NCorr] = {" [mm]", " [pad]"};
	int scalBinX[NCorr] = {10, 1};

	TFile *fPid = new TFile("/home/kailong/ExpData/Jul2018/Pid/pid-run-150--153__270--385.root");
	TTree *tPid;
	fPid->GetObject("tPid", tPid);
	tPid->SetEstimate(-1);

	string sDraw = "runNum:ZPid:APid:QPid:xMCP:yMCP:xPlaQ[0]:yPlaQ[0]:xCrdc[0]:yCrdc[0]:xCrdc[1]:yCrdc[1]";
	for (iAna = 0; iAna < NAna; iAna++)
	{
		sDraw += ":good_" + sAna[iAna];
		sDraw += ":tof_" + sAna[iAna];
	}
	string sCut = "ZPid>37&&isForCal>-1 && (good_" + sAna[0] + "==1";
	for (iAna = 1; iAna < NAna; iAna++)
		sCut += "||good_" + sAna[iAna] + "==1";
	sCut += ")";
	int nData = tPid->Draw(sDraw.c_str(), sCut.c_str(), "goff");
	double *runVal = tPid->GetVal(0);
	double *ZVal = tPid->GetVal(1);
	double *AVal = tPid->GetVal(2);
	double *QVal = tPid->GetVal(3);
	double *xVal[NCorr];
	// for (iCorr = 0; iCorr < NCorr; iCorr++)
	// 	xVal[iCorr] = tPid->GetVal(4 + iCorr);
	xVal[0] = tPid->GetVal(4);
	xVal[1] = tPid->GetVal(8);

	double *xPlaVal = tPid->GetVal(6);
	double *yPlaVal = tPid->GetVal(7);
	double *xCrdcVal[2], *yCrdcVal[2];
	for (i = 0; i < 2; i++)
	{
		xCrdcVal[i] = tPid->GetVal(8 + 2 * i);
		yCrdcVal[i] = tPid->GetVal(9 + 2 * i);
	}

	double *goodVal[NAna], *tofVal[NAna];
	int nVs = 2;
	for (iAna = 0; iAna < NAna; iAna++)
	{
		goodVal[iAna] = tPid->GetVal(12 + nVs * iAna);
		tofVal[iAna] = tPid->GetVal(12 + nVs * iAna + 1);
	}

	int iData = 0, nFilt = 0, nPts = 0;
	double xRang[NCorr][2] = {0}, tofRang[2] = {0};
	int nBinX[NCorr] = {0};
	bool **goodData = new bool *[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
		goodData[iAna] = new bool[nData];
	double xPoint[NCorr] = {0}, tofCorrVal = 0;
	double xPivot[NCorr] = {0}, tofPivot = 0, parVal[NPar] = {0};

	TFile *fDraw = new TFile("fDraw.root", "UPDATE");

	TGraph *grTofSgm[NAna], *grTofXsCorrSgm[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
	{
		grTofSgm[iAna] = new TGraph();
		grTofSgm[iAna]->SetNameTitle(("grTofSgm_" + sAna[iAna]).c_str(), "Raw ToF resolutions;Nuclide (=1000#times Z+A);#sigma_{ToF} [ps]");
		grTofSgm[iAna]->SetLineColor(kBlack);
		grTofSgm[iAna]->SetMarkerColor(kBlack);

		grTofXsCorrSgm[iAna] = new TGraph();
		grTofXsCorrSgm[iAna]->SetLineColor(kRed);
		grTofXsCorrSgm[iAna]->SetMarkerColor(kRed);
		grTofXsCorrSgm[iAna]->SetNameTitle(("grTofXsCorrSgm_" + sAna[iAna]).c_str(), "Xs-corrected ToF resolutions;Nuclide (=1000#times Z+A);#sigma_{ToF} [ps]");
	}
	int nuclide[NAna] = {0};
	for (iNucl = 0; iNucl < nNucl; iNucl++)
	// for (iNucl = 0; iNucl < 3; iNucl++)
	{
		if (isForCal[iNucl] < 0 || ZOrig[iNucl]<=37)
			continue;

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
		memset(tofXsCorr, 0, sizeof(tofXsCorr));
		memset(parTofXs, 0, sizeof(parTofXs));

		for (i = 0; i < NAna; i++)
			for (j = 0; j < nData; j++)
				goodData[i][j] = false;
		bool goodNucl = false;

		TGraph *grTofXs[NAna][NCorr], *grTofXsCorrXs[NAna][NCorr];
		for (iAna = 0; iAna < NAna; iAna++)
			for (iCorr = 0; iCorr < NCorr; iCorr++)
			{
				grTofXs[iAna][iCorr] = new TGraph();
				grTofXs[iAna][iCorr]->SetName(("grTofRawVS" + sCorr[iCorr] + "_" + sAna[iAna]).c_str());
				grTofXsCorrXs[iAna][iCorr] = new TGraph();
				grTofXsCorrXs[iAna][iCorr]->SetName(("grTofXsCorrVS" + sCorr[iCorr] + "_" + sAna[iAna]).c_str());
			}
		for (iAna = 0; iAna < NAna; iAna++)
		{
			memset(xRang, 0, sizeof(xRang));
			memset(tofRang, 0, sizeof(tofRang));
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
						grTofXs[iAna][iCorr]->SetPoint(nFilt, xVal[iCorr][iData], tofVal[iAna][iData]);
					nFilt++;
				}
			}
			if (nFilt < MINPOIN)
				continue;
			for (iCorr = 0; iCorr < NCorr; iCorr++)
			{
				xRang[iCorr][0] = grTofXs[iAna][iCorr]->GetMean(1) - NSGM * grTofXs[iAna][iCorr]->GetRMS(1);
				xRang[iCorr][1] = grTofXs[iAna][iCorr]->GetMean(1) + NSGM * grTofXs[iAna][iCorr]->GetRMS(1);
			}
			tofRang[0] = grTofXs[iAna][0]->GetMean(2) - NSGM * grTofXs[iAna][0]->GetRMS(2);
			tofRang[1] = grTofXs[iAna][0]->GetMean(2) + NSGM * grTofXs[iAna][0]->GetRMS(2);

			for (iCorr = 0; iCorr < NCorr; iCorr++)
				grTofXs[iAna][iCorr]->Set(0);
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
							grTofXs[iAna][iCorr]->SetPoint(nFilt, xVal[iCorr][iData], tofVal[iAna][iData]);
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
			memset(parTofXs, 0, sizeof(parTofXs));

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

					cnts[iAna] = grTofXs[iAna][0]->GetN();
					tofRaw[iAna][0] = grTofXs[iAna][0]->GetMean(2);
					tofRaw[iAna][2] = grTofXs[iAna][0]->GetRMS(2);
					tofRaw[iAna][1] = tofRaw[iAna][2] / sqrt(cnts[iAna]);

					TFile *fGlobFitPar = new TFile("fGlobFitPar.root");
					TTree *tGlobFitPar;
					fGlobFitPar->GetObject(("tGlobFitPar_" + sAna[iAna]).c_str(), tGlobFitPar);
					double ZGlob = 0, AGlob = 0, QGlob = 0;
					double globPar[NCorr][2] = {0};
					tGlobFitPar->SetBranchAddress("Z", &ZGlob);
					tGlobFitPar->SetBranchAddress("A", &AGlob);
					tGlobFitPar->SetBranchAddress("Q", &QGlob);
					for (iPar = 0; iPar < NCorr; iPar++)
					{
						tGlobFitPar->SetBranchAddress(("glob_" + strPar[iPar]).c_str(), globPar[iPar]);
					}
					memset(parVal, 0, sizeof(parVal));
					for (int iEnt = 0; iEnt < tGlobFitPar->GetEntries(); iEnt++)
					{
						tGlobFitPar->GetEntry(iEnt);
						if (abs(Z - ZGlob) < 1E-3 && abs(A - AGlob) < 1E-3 && abs(Q - QGlob) < 1E-3)
						{
							for (j = 1; j < NPar; j++)
							{
								parTofXs[iAna][j][0] = globPar[j - 1][0];
								parTofXs[iAna][j][1] = globPar[j - 1][1];
								parVal[j] = parTofXs[iAna][j][0];
							}
							break;
						}
					}
					delete fGlobFitPar;
					TFormula *fcnTofXs = new TFormula("fcnTofXs", "[0] + [1] * x[0] + [2] * x[1]");
					memset(xPivot, 0, sizeof(xPivot));
					xPivot[1] = 111;
					// xPivot[4] = 111;
					// xPivot[5] = 899;
					// xPivot[6] = 114;
					// xPivot[7] = 1055;
					// for (iCorr = 2; iCorr < NCorr; iCorr++)
					// 	xPivot[iCorr] = grTofXs[iAna][iCorr]->GetMean(1);
					tofPivot = fcnTofXs->EvalPar(xPivot, parVal);

					nPts = 0;
					for (iData = 0; iData < nData; iData++)
						if (goodData[iAna][iData])
						{
							memset(xPoint, 0, sizeof(xPoint));
							for (iCorr = 0; iCorr < NCorr; iCorr++)
								xPoint[iCorr] = xVal[iCorr][iData];
							tofCorrVal = tofVal[iAna][iData] + tofPivot - fcnTofXs->EvalPar(xPoint, parVal);

							for (iCorr = 0; iCorr < NCorr; iCorr++)
								grTofXsCorrXs[iAna][iCorr]->SetPoint(nPts, xVal[iCorr][iData], tofCorrVal);
							nPts++;
						}
					delete fcnTofXs;
					tofXsCorr[iAna][0] = grTofXsCorrXs[iAna][0]->GetMean(2);
					tofXsCorr[iAna][2] = grTofXsCorrXs[iAna][0]->GetRMS(2);
					tofXsCorr[iAna][1] = tofXsCorr[iAna][2] / sqrt(nPts);

					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						xRang[iCorr][0] = grTofXs[iAna][iCorr]->GetMean(1) - 2 * NSGM * grTofXs[iAna][iCorr]->GetRMS(1);
						xRang[iCorr][1] = grTofXs[iAna][iCorr]->GetMean(1) + 2 * NSGM * grTofXs[iAna][iCorr]->GetRMS(1);
						nBinX[iCorr] = (int)((xRang[iCorr][1] - xRang[iCorr][0]) * scalBinX[iCorr]);
					}
					tofRang[0] = tofRaw[iAna][0] - 2 * NSGM * tofRaw[iAna][2];
					tofRang[1] = tofRaw[iAna][0] + 2 * NSGM * tofRaw[iAna][2];
					int nBinTof = (int)((tofRang[1] - tofRang[0]) * 100);

					double tofXsCorrRang[2] = {0};
					tofXsCorrRang[0] = tofXsCorr[iAna][0] - 2 * NSGM * tofXsCorr[iAna][2];
					tofXsCorrRang[1] = tofXsCorr[iAna][0] + 2 * NSGM * tofXsCorr[iAna][2];
					int nBinTofXsCorr = (int)((tofXsCorrRang[1] - tofXsCorrRang[0]) * 100);

					TH2F *hTofXs[NCorr], *hTofXsCorrXs[NCorr];
					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						hTofXs[iCorr] = new TH2F(("hTofRawVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTofRawVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{raw} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTof, tofRang[0], tofRang[1]);

						hTofXsCorrXs[iCorr] = new TH2F(("hTofXsCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTofXsCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{XsCorr} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTofXsCorr, tofXsCorrRang[0], tofXsCorrRang[1]);
					}
					for (iData = 0; iData < cnts[iAna]; iData++)
						for (iCorr = 0; iCorr < NCorr; iCorr++)
						{
							hTofXs[iCorr]->Fill(grTofXs[iAna][iCorr]->GetX()[iData], grTofXs[iAna][iCorr]->GetY()[iData]);
							hTofXsCorrXs[iCorr]->Fill(grTofXsCorrXs[iAna][iCorr]->GetX()[iData], grTofXsCorrXs[iAna][iCorr]->GetY()[iData]);
						}
					TCanvas *cTofXs = new TCanvas(("cTofGlobXs_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("cTofGlobXs_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), 600 * (NCorr + 1), 900);
					cTofXs->Divide(NCorr + 1, 2);

					cTofXs->cd(1);
					TH1D *hTof = hTofXs[0]->ProjectionY(("hTofRaw__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
					hTof->SetTitle(("hTofRaw__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{raw} [ns];Counts/0.01ns").c_str());
					hTof->Draw();
					TFitResultPtr fitTof = hTof->Fit("gaus", "SQ");
					// grTofSgm[iAna]->SetPoint(nuclide[iAna], ZOrig[iNucl] * 1000 + AOrig[iNucl], 1000 * fitTof->Parameter(2));
					grTofSgm[iAna]->SetPoint(nuclide[iAna], ZOrig[iNucl] * 1000 + AOrig[iNucl], 1000 * grTofXs[iAna][0]->GetRMS(2));

					cTofXs->cd(NCorr + 2);
					TH1D *hTofXsCorr = hTofXsCorrXs[0]->ProjectionY(("hTofXsCorr__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
					hTofXsCorr->SetTitle(("hTofXsCorr__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{XsCorr} [ns];Counts/0.01ns").c_str());
					hTofXsCorr->Draw();
					TFitResultPtr fitTofXsCorr = hTofXsCorr->Fit("gaus", "SQ");
					grTofXsCorrSgm[iAna]->SetPoint(nuclide[iAna], ZOrig[iNucl] * 1000 + AOrig[iNucl], 1000 * fitTofXsCorr->Parameter(2));

					nuclide[iAna]++;

					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						cTofXs->cd(iCorr + 2);
						hTofXs[iCorr]->Draw(opt);

						cTofXs->cd(iCorr + NCorr + 3);
						hTofXsCorrXs[iCorr]->Draw(opt);
					}
					fDraw->cd();
					cTofXs->Write("", TObject::kOverwrite);
					cTofXs->SaveAs((sDirElec + "/cTofGlobXs_" + sAna[iAna] + "_" + strNucl[iNucl] + ".png").c_str());
					cTofXs->SaveAs((strGif[iAna] + "+20").c_str());
					delete cTofXs;

					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						delete hTofXs[iCorr];
						delete hTofXsCorrXs[iCorr];
					}
				}
			tTOF->Fill();

			memset(cnts, 0, sizeof(cnts));
			for (iData = 0; iData < nData; iData++)
			{
				run = 0;
				memset(tof0, 0, sizeof(tof0));
				memset(tof2XsCor, 0, sizeof(tof2XsCor));
				xMcp = 0;
				yMcp = 0;
				xPlaQS800 = 0;
				yPlaQS800 = 0;
				memset(xCrdc, 0, sizeof(xCrdc));
				memset(yCrdc, 0, sizeof(yCrdc));

				for (iAna = 0; iAna < NAna; iAna++)
					if (goodElec[iAna] && goodData[iAna][iData])
					{
						tof0[iAna] = tofVal[iAna][iData];
						tof2XsCor[iAna] = *(grTofXsCorrXs[iAna][0]->GetY() + cnts[iAna]);
						cnts[iAna]++;
					}
				bool goodDetail = false;
				for (iAna = 0; iAna < NAna; iAna++)
					if (goodElec[iAna] && goodData[iAna][iData])
					{
						run = runVal[iData];
						xMcp = xVal[0][iData];
						yMcp = xVal[1][iData];
						xPlaQS800 = xPlaVal[iData];
						yPlaQS800 = yPlaVal[iData];
						for (k = 0; k < 2; k++)
						{
							xCrdc[k] = xCrdcVal[k][iData];
							yCrdc[k] = yCrdcVal[k][iData];
						}
						goodDetail = true;
						break;
					}
				if (goodDetail)
					tTOFDetail->Fill();
			}
		}
		for (iAna = 0; iAna < NAna; iAna++)
			for (iCorr = 0; iCorr < NCorr; iCorr++)
			{
				delete grTofXs[iAna][iCorr];
				delete grTofXsCorrXs[iAna][iCorr];
			}
	}
	fTof->cd();
	tTOF->Write();
	tTOFDetail->Write();

	fDraw->cd();
	for (iAna = 0; iAna < NAna; iAna++)
	{
		TCanvas *cTofSgm = new TCanvas(("cTofGlobSgm_" + sAna[iAna]).c_str(), ("cTofGlobSgm_" + sAna[iAna]).c_str(), 2400, 900);
		cTofSgm->Divide(2, 1);
		cTofSgm->cd(1);
		grTofSgm[iAna]->Draw("AP*");
		cTofSgm->cd(2);
		grTofXsCorrSgm[iAna]->Draw("AP*");

		cTofSgm->Write("", TObject::kOverwrite);
		cTofSgm->SaveAs((sDirElec + "/cTofGlobSgm_" + sAna[iAna] + ".png").c_str());

		delete cTofSgm;
		delete grTofXsCorrSgm[iAna];
		delete grTofSgm[iAna];
	}

	delete fDraw;
	for (iAna = 0; iAna < NAna; iAna++)
		delete[] goodData[iAna];
	delete[] goodData;
	delete fPid;
	delete tTOFDetail;
	delete tTOF;
	delete fTof;
}

int main()
{
	Pid2TofGlob();
	return 0;
}