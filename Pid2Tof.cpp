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
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Fit/FitResult.h"
#include "Math/WrappedMultiTF1.h"

double myFunc(double *x, double *par)
{
	// return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[1] + par[5] * x[2] + par[6] * x[3] + par[7]*x[4];
	// return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[1];
	// return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[1] + par[5] * x[2] + par[6] * x[3] + par[7] * x[4] + par[8] * x[5];
	// return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[1] + par[5] * x[2] + par[6] * x[3];
	return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[1] + par[5] * x[2];
	// return par[0] + par[1] * x[0] + par[2] * x[1];
}

void Pid2Tof()
{
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "MigradImproved");
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

	const int NCorr = 3, NSGM = 3, MINPOIN = 100;
	const int NPar = NCorr + 3;

	string sRoot, sTree;
	int i, j, k;

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

	TFile *fTof = new TFile("/home/kailong/ExpData/Jul2018/Pid/tof-run-150--153__270--385.root", "RECREATE");
	TTree *tTOF = new TTree("tTOF", "tree for TOF analysis");
	TTree *tTOFDetail = new TTree("tTOFDetail", "tree for TOF analysis with details");
	string sNucl = "";
	int isCal = 0, Z = 0, A = 0, Q = 0, cnts[NAna] = {0};
	bool goodElec[NAna] = {0};
	double AoQ = 0, MoQ[2] = {0};
	double tofRaw[NAna][3] = {0}, tofXmcpCorr[NAna][3], tofXsCorr[NAna][3] = {0};											//[3]: [0] for mean TOF; [1] for sigma of mean TOF; [2] for sigma of TOF distribution
	double parTofXmcp[NAna][2][2] = {0}, parTofXs[NAna][NPar][2] = {0}, tofPivot[NAna] = {0}, meanXs[NAna][NCorr][2] = {0}; //second [2]: [0] value; [1] error
	double tofAME = 0;

	double run = 0;
	double tof0[NAna] = {0}, tof1XmcpCor[NAna] = {0}, tof2XsCor[NAna] = {0};
	double xMcp = 0, yMcp = 0, xPlaQS800 = 0, yPlaQS800 = 0, xCrdc[2] = {0}, yCrdc[2] = {0}, xPlaCrdcS800 = 0, yPlaCrdcS800 = 0, qPmtS800[4] = {0}, qPlaS800 = 0;

	string sDirElec[NAna], strGif[NAna];
	FILE *fParCorrTof[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
	{
		sDirElec[iAna] = "./Graphs/" + sAna[iAna];
		if (access(sDirElec[iAna].c_str(), F_OK) != 0)
			system(("mkdir " + sDirElec[iAna]).c_str());
		strGif[iAna] = sDirElec[iAna] + "/cTofXs_" + sAna[iAna] + ".gif";
		if (access(strGif[iAna].c_str(), F_OK) == 0)
			system(("rm -f " + strGif[iAna]).c_str());

		fParCorrTof[iAna] = fopen(("fParCorrTof_" + sAna[iAna] + ".dat").c_str(), "w+");
		fprintf(fParCorrTof[iAna], "#Nuclide  Z    A   Q     ");
		for (i = 0; i < NPar; i++)
			fprintf(fParCorrTof[iAna], "par_%d        ", i);
		fprintf(fParCorrTof[iAna], "\n");
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
		tTOF->Branch(("parTofXs_" + sAna[iAna]).c_str(), parTofXs[iAna], ("parTofXs_" + sAna[iAna] + "[" + to_string(NPar) + "][2]/D").c_str());
		tTOF->Branch(("tofPivot_" + sAna[iAna]).c_str(), &tofPivot[iAna], ("tofPivot_" + sAna[iAna] + "/D").c_str());
		tTOF->Branch(("meanXs_" + sAna[iAna]).c_str(), meanXs[iAna], ("meanXs_" + sAna[iAna] + "[" + to_string(NCorr) + "][2]/D").c_str());
	}
	tTOF->Branch("tofAME", &tofAME, "tofAME/D");

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
		tTOFDetail->Branch(("tof1XmcpCor_" + sAna[iAna]).c_str(), &tof1XmcpCor[iAna], ("tof1XmcpCor_" + sAna[iAna] + "/D").c_str());
		tTOFDetail->Branch(("tof2XsCor_" + sAna[iAna]).c_str(), &tof2XsCor[iAna], ("tof2XsCor_" + sAna[iAna] + "/D").c_str());
	}
	tTOFDetail->Branch("xMcp", &xMcp, "xMcp/D");
	tTOFDetail->Branch("yMcp", &yMcp, "yMcp/D");
	tTOFDetail->Branch("xPlaQS800", &xPlaQS800, "xPlaQS800/D");
	tTOFDetail->Branch("yPlaQS800", &yPlaQS800, "yPlaQS800/D");
	tTOFDetail->Branch("xCrdc", xCrdc, "xCrdc[2]/D");
	tTOFDetail->Branch("yCrdc", yCrdc, "yCrdc[2]/D");
	tTOFDetail->Branch("xPlaCrdcS800", &xPlaCrdcS800, "xPlaCrdcS800/D");
	tTOFDetail->Branch("yPlaCrdcS800", &yPlaCrdcS800, "yPlaCrdcS800/D");
	tTOFDetail->Branch("qPmtS800", qPmtS800, "qPmtS800[4]/D");
	tTOFDetail->Branch("qPlaS800", &qPlaS800, "qPlaS800/D");

	// string sCorr[NCorr] = {"Xmcp", "Ymcp", "Qpmt1S800", "Qpmt2S800", "Qpmt3S800", "Qpmt4S800"};
	// string sCorrUnit[NCorr] = {" [mm]", " [mm]", "", "", "", ""};
	// int scalBinX[NCorr] = {10, 10, 1, 1, 1, 1};

	// string sCorr[NCorr] = {"Xmcp", "Ymcp", "Qpmt1S800", "Qpmt2S800", "Qpmt3S800", "Qpmt4S800", "xCRDC1", "yCRDC1", "xCRDC2", "yCRDC2", "", ""};
	// string sCorrUnit[NCorr] = {" [mm]", " [mm]", " [ch]", " [ch]", " [ch]", " [ch]", " [pad]", " [ch]", " [pad]", " [ch]", "", ""};
	// int scalBinX[NCorr] = {10, 10, 1, 1, 1, 1, 1, 1, 1, 1, 100, 100};

	// string sCorr[NCorr] = {"Xmcp", "Ymcp", "xPlaQS800", "yPlaQS800", "xCRDC1", "yCRDC1", "xCRDC2", "yCRDC2"};
	// string sCorrUnit[NCorr] = {" [mm]", " [mm]", "", "", " [pad]", " [ch]", " [pad]", " [ch]"};
	// int scalBinX[NCorr] = {10, 10, 100, 100, 1, 1, 1, 1};

	// string sCorr[NCorr] = {"Xmcp", "Ymcp", "xCRDC1", "yCRDC1", "xCRDC2", "yCRDC2"};
	// string sCorrUnit[NCorr] = {" [mm]", " [mm]", " [pad]", " [ch]", " [pad]", " [ch]"};
	// int scalBinX[NCorr] = {10, 10, 1, 1, 1, 1};
	// string sCorr[NCorr] = {"Xmcp", "Ymcp", "xPlaCrdcS800", "yPlaCrdcS800"};
	// string sCorrUnit[NCorr] = {" [mm]", " [mm]", " [pad]", " [ch]"};
	// int scalBinX[NCorr] = {10, 10, 2, 1};
	string sCorr[NCorr] = {"Xmcp", "xCRDC1", "xCRDC2"};
	string sCorrUnit[NCorr] = {" [mm]", " [pad]", " [pad]"};
	int scalBinX[NCorr] = {10, 2, 2};

	// string sCorr[NCorr] = {"Xmcp", "Ymcp", "xPlaQS800", "yPlaQS800"};
	// string sCorrUnit[NCorr] = {" [mm]", " [mm]", "", ""};
	// int scalBinX[NCorr] = {10, 10, 100, 100};

	// string sCorr[NCorr] = {"Xmcp", "Ymcp", "xPlaQS800", "yPlaQS800", "xCrdc1"};
	// string sCorrUnit[NCorr] = {" [mm]", " [mm]", "", "", " [pad]"};
	// int scalBinX[NCorr] = {10, 10, 100, 100, 1};

	// string sCorr[NCorr] = {"Xmcp", "xCrdc1"};
	// string sCorrUnit[NCorr] = {" [mm]", " [pad]"};
	// int scalBinX[NCorr] = {10, 1};

	sRoot = "/home/kailong/ExpData/Jul2018/Pid/pid-run-150--153__270--385.root";
	sTree = "tPid";

	TFile *fPid = new TFile(sRoot.c_str());
	TTree *tPid;
	fPid->GetObject(sTree.c_str(), tPid);
	tPid->SetEstimate(-1);

	// string sDraw = "ZPid:APid:QPid:xMCP:yMCP:egyPMT[0][0]:egyPMT[0][1]:egyPMT[0][2]:egyPMT[0][3]:xCrdc[0]:yCrdc[0]:xCrdc[1]:yCrdc[1]:xPlaQ[0]:yPlaQ[0]";
	string sDraw = "runNum:ZPid:APid:QPid:xMCP:yMCP:xPlaCrdcS800:yPlaCrdcS800:xCrdc[0]:yCrdc[0]:xCrdc[1]:yCrdc[1]:xPlaQ[0]:yPlaQ[0]:egyPMT[0][0]:egyPMT[0][1]:egyPMT[0][2]:egyPMT[0][3]:egyPla[0]";
	for (iAna = 0; iAna < NAna; iAna++)
	{
		sDraw += ":good_" + sAna[iAna];
		sDraw += ":tof_" + sAna[iAna];
	}
	string sCut = "isForCal>-1 && ZPid>37";
	// string sCut = "isForCal>-1";
	int nData = tPid->Draw(sDraw.c_str(), sCut.c_str(), "goff");
	double *runVal = tPid->GetVal(0);
	double *ZVal = tPid->GetVal(1);
	double *AVal = tPid->GetVal(2);
	double *QVal = tPid->GetVal(3);
	double *xVal[NCorr];
	for (iCorr = 0; iCorr < 1; iCorr++)
		xVal[iCorr] = tPid->GetVal(4 + iCorr);
	xVal[1] = tPid->GetVal(8);
	xVal[2] = tPid->GetVal(10);
	// xVal[4] = tPid->GetVal(10);
	// xVal[5] = tPid->GetVal(11);

	double *xPlaCrdcVal = tPid->GetVal(6);
	double *yPlaCrdcVal = tPid->GetVal(7);
	double *xCrdcVal[2], *yCrdcVal[2];
	for (i = 0; i < 2; i++)
	{
		xCrdcVal[i] = tPid->GetVal(8 + 2 * i);
		yCrdcVal[i] = tPid->GetVal(9 + 2 * i);
	}
	double *xPlaQVal = tPid->GetVal(12);
	double *yPlaQVal = tPid->GetVal(13);
	double *qPmtS800Val[4];
	for (i = 0; i < 4; i++)
		qPmtS800Val[i] = tPid->GetVal(14 + i);
	double *qPlaS800Val = tPid->GetVal(18);

	double *goodVal[NAna], *tofVal[NAna];
	int nVs = 2;
	for (iAna = 0; iAna < NAna; iAna++)
	{
		goodVal[iAna] = tPid->GetVal(19 + nVs * iAna);
		tofVal[iAna] = tPid->GetVal(19 + nVs * iAna + 1);
	}

	int iData = 0, nFilt = 0, nPts = 0;
	double xRang[NCorr][2] = {0}, tofRang[2] = {0}, tofXsCorrRang[2] = {0};
	int nBinX[NCorr] = {0};
	bool **goodData = new bool *[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
		goodData[iAna] = new bool[nData];
	double xPt[NCorr] = {0}, tofCorrVal = 0;
	double xPivot[NCorr] = {0}, parVal[NPar] = {0};

	TFile *fDraw = new TFile("fDraw.root", "UPDATE");

	TGraph *grTofSgm[NAna], *grTofXmcpCorrSgm[NAna], *grTofXsCorrSgm[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
	{
		grTofSgm[iAna] = new TGraph();
		grTofSgm[iAna]->SetNameTitle(("grTofSgm_" + sAna[iAna]).c_str(), "Raw ToF resolutions;Nuclide (=1000#times Z+A);#sigma_{ToF} [ps]");
		grTofSgm[iAna]->SetLineColor(kBlack);
		grTofSgm[iAna]->SetMarkerColor(kBlack);

		grTofXmcpCorrSgm[iAna] = new TGraph();
		grTofXmcpCorrSgm[iAna]->SetLineColor(kBlue);
		grTofXmcpCorrSgm[iAna]->SetMarkerColor(kBlue);
		grTofXmcpCorrSgm[iAna]->SetNameTitle(("grTofXmcpCorrSgm_" + sAna[iAna]).c_str(), "Xmcp-corrected ToF resolutions;Nuclide (=1000#times Z+A);#sigma_{ToF} [ps]");

		grTofXsCorrSgm[iAna] = new TGraph();
		grTofXsCorrSgm[iAna]->SetLineColor(kRed);
		grTofXsCorrSgm[iAna]->SetMarkerColor(kRed);
		grTofXsCorrSgm[iAna]->SetNameTitle(("grTofXsCorrSgm_" + sAna[iAna]).c_str(), "Xs-corrected ToF resolutions;Nuclide (=1000#times Z+A);#sigma_{ToF} [ps]");
	}
	int nuclide[NAna] = {0};
	// for (iNucl = 0; iNucl < nNucl; iNucl++)
	for (iNucl = 0; iNucl < 3; iNucl++)
	{
		if (isForCal[iNucl] < 0 || ZOrig[iNucl] <= 37)
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
		memset(tofXmcpCorr, 0, sizeof(tofXmcpCorr));
		memset(parTofXmcp, 0, sizeof(parTofXmcp));
		memset(tofXsCorr, 0, sizeof(tofXsCorr));
		memset(parTofXs, 0, sizeof(parTofXs));
		memset(tofPivot, 0, sizeof(tofPivot));
		memset(meanXs, 0, sizeof(meanXs));
		tofAME = 0;

		for (i = 0; i < NAna; i++)
			for (j = 0; j < nData; j++)
				goodData[i][j] = false;
		bool goodNucl = false;

		TGraph *grTofXs[NAna][NCorr], *grTofXmcpCorrXs[NAna][NCorr], *grTofXsCorrXs[NAna][NCorr];
		for (iAna = 0; iAna < NAna; iAna++)
			for (iCorr = 0; iCorr < NCorr; iCorr++)
			{
				grTofXs[iAna][iCorr] = new TGraph();
				grTofXs[iAna][iCorr]->SetName(("grTofRawVS" + sCorr[iCorr] + "_" + sAna[iAna]).c_str());
				grTofXmcpCorrXs[iAna][iCorr] = new TGraph();
				grTofXmcpCorrXs[iAna][iCorr]->SetName(("grTofXmcpCorrVS" + sCorr[iCorr] + "_" + sAna[iAna]).c_str());
				grTofXsCorrXs[iAna][iCorr] = new TGraph();
				grTofXsCorrXs[iAna][iCorr]->SetName(("grTofXsCorrVS" + sCorr[iCorr] + "_" + sAna[iAna]).c_str());
			}
		for (iAna = 0; iAna < NAna; iAna++)
		{
			memset(xRang, 0, sizeof(xRang));
			memset(tofRang, 0, sizeof(tofRang));
			memset(xPt, 0, sizeof(xPt));
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
			{
				goodElec[iAna] = true;
			}
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
			tofAME = tofPivotAME[iNucl];

			for (iAna = 0; iAna < NAna; iAna++)
				if (goodElec[iAna])
				{
					memset(xRang, 0, sizeof(xRang));
					memset(nBinX, 0, sizeof(nBinX));
					memset(tofRang, 0, sizeof(tofRang));
					memset(tofXsCorrRang, 0, sizeof(tofXsCorrRang));
					memset(xPt, 0, sizeof(xPt));
					memset(xPivot, 0, sizeof(xPivot));
					memset(parVal, 0, sizeof(parVal));
					tofCorrVal = 0;
					for (iCorr = 0; iCorr < NCorr; iCorr++)
						grTofXsCorrXs[iAna][iCorr]->Set(0);

					cnts[iAna] = grTofXs[iAna][0]->GetN();
					tofRaw[iAna][0] = grTofXs[iAna][0]->GetMean(2);
					tofRaw[iAna][2] = grTofXs[iAna][0]->GetRMS(2);
					tofRaw[iAna][1] = tofRaw[iAna][2] / sqrt(cnts[iAna]);

					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						meanXs[iAna][iCorr][0] = grTofXs[iAna][iCorr]->GetMean(1);
						meanXs[iAna][iCorr][1] = grTofXs[iAna][iCorr]->GetRMS(1);
					}

					TFitResultPtr fitTofXmcp = grTofXs[iAna][0]->Fit("pol1", "SQ0");
					for (i = 0; i < 2; i++)
					{
						parTofXmcp[iAna][i][0] = fitTofXmcp->Parameter(i);
						parTofXmcp[iAna][i][1] = fitTofXmcp->ParError(i);
					}
					xPivot[0] = 0;
					nPts = 0;
					for (iData = 0; iData < nData; iData++)
						if (goodData[iAna][iData])
						{
							tofCorrVal = tofVal[iAna][iData] + parTofXmcp[iAna][1][0] * xPivot[0] - parTofXmcp[iAna][1][0] * xVal[0][iData];
							for (iCorr = 0; iCorr < NCorr; iCorr++)
								grTofXmcpCorrXs[iAna][iCorr]->SetPoint(nPts, xVal[iCorr][iData], tofCorrVal);
							nPts++;
						}
					tofXmcpCorr[iAna][0] = grTofXmcpCorrXs[iAna][0]->GetMean(2);
					tofXmcpCorr[iAna][2] = grTofXmcpCorrXs[iAna][0]->GetRMS(2);
					tofXmcpCorr[iAna][1] = tofXmcpCorr[iAna][2] / sqrt(nPts);

					ROOT::Fit::BinData dataTofXs(nData, NCorr, ROOT::Fit::BinData::kNoError);
					for (iData = 0; iData < nData; iData++)
						if (goodData[iAna][iData])
						{
							memset(xPt, 0, sizeof(xPt));
							for (iCorr = 0; iCorr < NCorr; iCorr++)
								xPt[iCorr] = xVal[iCorr][iData];
							dataTofXs.Add(xPt, tofVal[iAna][iData]);
						}
					TF1 fcnTofXs("fcnTofXs", myFunc, 0, 1, NPar);
					ROOT::Math::WrappedMultiTF1 fitFunc(fcnTofXs, NCorr);
					ROOT::Fit::Fitter fitter;
					fitter.SetFunction(fitFunc, true);
					for (i = 0; i < NPar; i++)
						fitter.Config().ParSettings(i).SetValue(1);
					fitter.LinearFit(dataTofXs);
					// fitter.LeastSquareFit(dataTofXs);
					ROOT::Fit::FitResult resultFit = fitter.Result();

					fprintf(fParCorrTof[iAna], "%s   %d   %d  %d    ", sNucl.c_str(), Z, A, Q);
					for (j = 0; j < NPar; j++)
					{
						parTofXs[iAna][j][0] = resultFit.Parameter(j);
						parTofXs[iAna][j][1] = resultFit.ParError(j);
						parVal[j] = parTofXs[iAna][j][0];
						fprintf(fParCorrTof[iAna], "%.7G  ", parVal[j]);
					}
					fprintf(fParCorrTof[iAna], "\n");
					memset(xPivot, 0, sizeof(xPivot));
					// for (iCorr = 0; iCorr < NCorr; iCorr++)
					// 	xPivot[iCorr] = meanXs[iAna][iCorr][0];
					// xPivot[0] = -8.0770565;
					// xPivot[1] = -1.1232244;
					// xPivot[2] = 110.49939;  //X CRDC-1
					// xPivot[3] = 898.22997;  //Y CRDC-1
					// xPivot[4] = 114.42846;  //X CRDC-2
					// xPivot[5] = 1054.6725;  //Y CRDC-2

					// xPivot[2] = 111;  //X Plastic in S800 from CRDCs
					// xPivot[3] = 913;  //Y Plastic in S800 from CRDCs

					// xPivot[0] = -0.899170;
					// xPivot[1] = 110.050144;
					// xPivot[2] = 113.278606;

					xPivot[1] = 110;
					xPivot[2] = 113;

					// xPivot[0] = -10.398019;
					// xPivot[1] = 110.565783;
					// xPivot[2] = 114.595472;

					// xPivot[0] = -1.455118;
					// xPivot[1] = -2.118212;
					// xPivot[2] = 109.676419;
					// xPivot[3] = 902.504500;
					// xPivot[4] = 113.429762;
					// xPivot[5] = 1055.472889;

					// xPivot[0] = -0.899170;
					// xPivot[1] = 110.050144;
					// xPivot[2] = 113.278606;

					tofPivot[iAna] = fcnTofXs(xPivot, parVal);
					// tofPivot[iAna] = tofRaw[iAna][0];
					// tofPivot[iAna] = tofXmcpCorr[iAna][0];
					// tofPivot[iAna] = tofPivotAME[iNucl];

					// if (isCal == 1)
					// tofPivot[iAna] = tofPivotAME[iNucl];
					// else
					// tofPivot[iAna] = 60.763 / 0.299792458 * sqrt(1 + 9.654255907 * pow(1.0 * A / Q / 3.7211, 2));

					// printf("tofPivot=%f, tofPivotAME=%f, Delta=%f\n", tofPivot[iAna], tofPivotAME[iNucl], tofPivot[iAna]-tofPivotAME[iNucl]);
					// tofPivot[iAna] = tofPivotAME[iNucl];

					nPts = 0;
					for (iData = 0; iData < nData; iData++)
						if (goodData[iAna][iData])
						{
							memset(xPt, 0, sizeof(xPt));
							for (iCorr = 0; iCorr < NCorr; iCorr++)
								xPt[iCorr] = xVal[iCorr][iData];
							tofCorrVal = tofVal[iAna][iData] + tofPivot[iAna] - fcnTofXs(xPt, parVal);

							for (iCorr = 0; iCorr < NCorr; iCorr++)
								grTofXsCorrXs[iAna][iCorr]->SetPoint(nPts, xVal[iCorr][iData], tofCorrVal);
							nPts++;
						}
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

					double tofXmcpCorrRang[2] = {0};
					tofXmcpCorrRang[0] = tofXmcpCorr[iAna][0] - 2 * NSGM * tofXmcpCorr[iAna][2];
					tofXmcpCorrRang[1] = tofXmcpCorr[iAna][0] + 2 * NSGM * tofXmcpCorr[iAna][2];
					tofXsCorrRang[0] = tofXsCorr[iAna][0] - 2 * NSGM * tofXsCorr[iAna][2];
					tofXsCorrRang[1] = tofXsCorr[iAna][0] + 2 * NSGM * tofXsCorr[iAna][2];
					int nBinTofXmcpCorr = (int)((tofXmcpCorrRang[1] - tofXmcpCorrRang[0]) * 100);
					int nBinTofXsCorr = (int)((tofXsCorrRang[1] - tofXsCorrRang[0]) * 100);

					TH2F *hTofXs[NCorr], *hTofXmcpCorrXs[NCorr], *hTofXsCorrXs[NCorr];
					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						hTofXs[iCorr] = new TH2F(("hTofRawVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTofRawVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{raw} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTof, tofRang[0], tofRang[1]);

						hTofXmcpCorrXs[iCorr] = new TH2F(("hTofXmcpCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTofXmcpCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{XmcpCorr} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTofXmcpCorr, tofXmcpCorrRang[0], tofXmcpCorrRang[1]);

						hTofXsCorrXs[iCorr] = new TH2F(("hTofXsCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("hTofXsCorrVS" + sCorr[iCorr] + "_" + sAna[iAna] + "_" + strNucl[iNucl] + ";" + sCorr[iCorr] + sCorrUnit[iCorr] + ";ToF_{XsCorr} [ns]").c_str(), nBinX[iCorr], xRang[iCorr][0], xRang[iCorr][1], nBinTofXsCorr, tofXsCorrRang[0], tofXsCorrRang[1]);
					}
					for (iData = 0; iData < cnts[iAna]; iData++)
						for (iCorr = 0; iCorr < NCorr; iCorr++)
						{
							hTofXs[iCorr]->Fill(grTofXs[iAna][iCorr]->GetX()[iData], grTofXs[iAna][iCorr]->GetY()[iData]);
							hTofXmcpCorrXs[iCorr]->Fill(grTofXmcpCorrXs[iAna][iCorr]->GetX()[iData], grTofXmcpCorrXs[iAna][iCorr]->GetY()[iData]);
							hTofXsCorrXs[iCorr]->Fill(grTofXsCorrXs[iAna][iCorr]->GetX()[iData], grTofXsCorrXs[iAna][iCorr]->GetY()[iData]);
						}
					TCanvas *cTofXs = new TCanvas(("cTofXs_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), ("cTofXs_" + sAna[iAna] + "_" + strNucl[iNucl]).c_str(), 1200 * (NCorr + 1), 2700);
					cTofXs->SetCanvasSize(1200 * (NCorr + 1), 2700);
					cTofXs->Divide(NCorr + 1, 3);

					cTofXs->cd(1);
					TH1D *hTof = hTofXs[0]->ProjectionY(("hTofRaw__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
					hTof->SetTitle(("hTofRaw__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{raw} [ns];Counts/0.01ns").c_str());
					hTof->Draw();
					hTof->GetYaxis()->SetMaxDigits(3);
					TFitResultPtr fitTof = hTof->Fit("gaus", "SQ");
					// grTofSgm[iAna]->SetPoint(nuclide[iAna], ZOrig[iNucl] * 1000 + AOrig[iNucl], 1000 * fitTof->Parameter(2));
					grTofSgm[iAna]->SetPoint(nuclide[iAna], ZOrig[iNucl] * 1000 + AOrig[iNucl], 1000 * grTofXs[iAna][0]->GetRMS(2));

					cTofXs->cd(NCorr + 2);
					TH1D *hTofXmcpCorr = hTofXmcpCorrXs[0]->ProjectionY(("hTofXmcpCorr__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
					hTofXmcpCorr->SetTitle(("hTofXmcpCorr__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{XmcpCorr} [ns];Counts/0.01ns").c_str());
					hTofXmcpCorr->Draw();
					hTofXmcpCorr->GetYaxis()->SetMaxDigits(3);
					TFitResultPtr fitTofXmcpCorr = hTofXmcpCorr->Fit("gaus", "SQ");
					grTofXmcpCorrSgm[iAna]->SetPoint(nuclide[iAna], ZOrig[iNucl] * 1000 + AOrig[iNucl], 1000 * fitTofXmcpCorr->Parameter(2));

					cTofXs->cd(2 * NCorr + 3);
					TH1D *hTofXsCorr = hTofXsCorrXs[0]->ProjectionY(("hTofXsCorr__" + sAna[iAna] + "_" + strNucl[iNucl]).c_str());
					hTofXsCorr->SetTitle(("hTofXsCorr__" + sAna[iAna] + "_" + strNucl[iNucl] + ";ToF_{XsCorr} [ns];Counts/0.01ns").c_str());
					hTofXsCorr->Draw();
					hTofXsCorr->GetYaxis()->SetMaxDigits(3);
					TFitResultPtr fitTofXsCorr = hTofXsCorr->Fit("gaus", "SQ");
					grTofXsCorrSgm[iAna]->SetPoint(nuclide[iAna], ZOrig[iNucl] * 1000 + AOrig[iNucl], 1000 * fitTofXsCorr->Parameter(2));

					nuclide[iAna]++;

					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						cTofXs->cd(iCorr + 2);
						hTofXs[iCorr]->Draw(opt);
						hTofXs[iCorr]->GetYaxis()->SetMaxDigits(3);

						cTofXs->cd(iCorr + NCorr + 3);
						hTofXmcpCorrXs[iCorr]->Draw(opt);
						hTofXmcpCorrXs[iCorr]->GetYaxis()->SetMaxDigits(3);

						cTofXs->cd(iCorr + 2 * NCorr + 4);
						hTofXsCorrXs[iCorr]->Draw(opt);
						hTofXsCorrXs[iCorr]->GetYaxis()->SetMaxDigits(3);
					}
					cTofXs->Write("", TObject::kOverwrite);
					cTofXs->SaveAs((sDirElec[iAna] + "/cTofXs_" + sAna[iAna] + "_" + strNucl[iNucl] + ".png").c_str());
					cTofXs->SaveAs((strGif[iAna] + "+20").c_str());
					delete cTofXs;

					for (iCorr = 0; iCorr < NCorr; iCorr++)
					{
						delete hTofXs[iCorr];
						delete hTofXmcpCorrXs[iCorr];
						delete hTofXsCorrXs[iCorr];
					}
				}
			tTOF->Fill();

			memset(cnts, 0, sizeof(cnts));
			for (iData = 0; iData < nData; iData++)
			{
				run = 0;
				memset(tof0, 0, sizeof(tof0));
				memset(tof1XmcpCor, 0, sizeof(tof1XmcpCor));
				memset(tof2XsCor, 0, sizeof(tof2XsCor));
				xMcp = 0;
				yMcp = 0;
				memset(qPmtS800, 0, sizeof(qPmtS800));
				xPlaQS800 = 0;
				yPlaQS800 = 0;
				memset(xCrdc, 0, sizeof(xCrdc));
				memset(yCrdc, 0, sizeof(yCrdc));
				xPlaCrdcS800 = 0;
				yPlaCrdcS800 = 0;
				memset(qPmtS800, 0, sizeof(qPmtS800));
				qPlaS800 = 0;

				for (iAna = 0; iAna < NAna; iAna++)
					if (goodElec[iAna] && goodData[iAna][iData])
					{
						tof0[iAna] = tofVal[iAna][iData];
						tof1XmcpCor[iAna] = *(grTofXmcpCorrXs[iAna][0]->GetY() + cnts[iAna]);
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
						xPlaQS800 = xPlaQVal[iData];
						yPlaQS800 = yPlaQVal[iData];
						for (k = 0; k < 2; k++)
						{
							xCrdc[k] = xCrdcVal[k][iData];
							yCrdc[k] = yCrdcVal[k][iData];
						}
						xPlaCrdcS800 = xPlaCrdcVal[iData];
						yPlaCrdcS800 = yPlaCrdcVal[iData];
						for (k = 0; k < 4; k++)
							qPmtS800[k] = qPmtS800Val[k][iData];
						qPlaS800 = qPlaS800Val[iData];
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
				delete grTofXmcpCorrXs[iAna][iCorr];
				delete grTofXsCorrXs[iAna][iCorr];
			}
	}
	fTof->cd();
	tTOF->Write();
	tTOFDetail->Write();

	fDraw->cd();
	for (iAna = 0; iAna < NAna; iAna++)
	{
		TCanvas *cTofSgm = new TCanvas(("cTofSgm_" + sAna[iAna]).c_str(), ("cTofSgm_" + sAna[iAna]).c_str(), 2400, 1800);
		cTofSgm->SetCanvasSize(2400, 1800);
		cTofSgm->Divide(2, 2);
		cTofSgm->cd(1);
		grTofSgm[iAna]->Draw("AP*");
		cTofSgm->cd(2);
		grTofXmcpCorrSgm[iAna]->Draw("AP*");
		cTofSgm->cd(3);
		grTofXsCorrSgm[iAna]->Draw("AP*");

		cTofSgm->Write("", TObject::kOverwrite);
		cTofSgm->SaveAs((sDirElec[iAna] + "/cTofSgm_" + sAna[iAna] + ".png").c_str());

		delete cTofSgm;
		delete grTofXsCorrSgm[iAna];
		delete grTofXmcpCorrSgm[iAna];
		delete grTofSgm[iAna];
	}

	delete fDraw;
	for (iAna = 0; iAna < NAna; iAna++)
	{
		fclose(fParCorrTof[iAna]);
		delete[] goodData[iAna];
	}
	delete[] goodData;
	delete fPid;
	delete tTOFDetail;
	delete tTOF;
	delete fTof;
}

int main()
{
	Pid2Tof();
	return 0;
}