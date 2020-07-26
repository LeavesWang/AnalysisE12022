#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TF2.h"
#include "TGraph2DErrors.h"
#include "TH2F.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TLegend.h"
#include "TMath.h"
#include "TF1.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Fit/FitResult.h"
#include "Math/WrappedMultiTF1.h"

const int NDIM1 = 1, NDIM2 = 3;
const int NDIM = NDIM1 + NDIM2;
const int PolyOrd1 = 3, PolyOrd2 = 3;

double myFunc(double *x, double *par)
{
	double f = par[0];
	int i, j, k, m, n, iPar = 0;
	
	// if (PolyOrd1 >= 1)
	// 	for (i = 0; i < NDIM1; i++)
	// 		f += par[++iPar] * x[i];
	// if (PolyOrd1 >= 2)
	// 	for (j = 0; j < NDIM1; j++)
	// 		f += par[++iPar] * x[j] * x[j];
	// if (PolyOrd1 >= 3)
	// 	for (k = 0; k < NDIM1; k++)
	// 		f += par[++iPar] * x[k] * x[k] * x[k];
	// if (PolyOrd1 >= 4)
	// 	for (m = 0; m < NDIM1; m++)
	// 		f += par[++iPar] * x[m] * x[m] * x[m] * x[m];
	// if (PolyOrd1 >= 5)
	// 	for (n = 0; n < NDIM1; n++)
	// 		f += par[++iPar] * x[n] * x[n] * x[n] * x[n] * x[n];

	if (PolyOrd1 >= 1)
		for (i = 0; i < NDIM1 + 1; i++)
		{
			f += par[++iPar] * x[i];
			if (PolyOrd1 >= 2)
				for (j = i; j < NDIM1 + 1; j++)
				{
					f += par[++iPar] * x[i] * x[j];
					if (PolyOrd1 >= 3)
						for (k = j; k < NDIM1 + 1; k++)
						{
							f += par[++iPar] * x[i] * x[j] * x[k];
							if (PolyOrd1 >= 4)
								for (m = k; m < NDIM1 + 1; m++)
								{
									f += par[++iPar] * x[i] * x[j] * x[k] * x[m];
									if (PolyOrd1 >= 5)
										for (n = m; n < NDIM1 + 1; n++)
											f += par[++iPar] * x[i] * x[j] * x[k] * x[m] * x[n];
								}
						}
				}
		}
	// f += par[++iPar] * x[1] * x[1] * x[1] * x[1];
	if (PolyOrd2 >= 1)
		for (i = 1; i < NDIM2; i++)
		{
			f += par[++iPar] * x[NDIM1 + i];
			if (PolyOrd2 >= 2)
				for (j = i; j < NDIM2; j++)
				{
					f += par[++iPar] * x[NDIM1 + i] * x[NDIM1 + j];
					if (PolyOrd2 >= 3)
						for (k = j; k < NDIM2; k++)
						{
							f += par[++iPar] * x[NDIM1 + i] * x[NDIM1 + j] * x[NDIM1 + k];
							if (PolyOrd2 >= 4)
								for (m = k; m < NDIM2; m++)
								{
									f += par[++iPar] * x[NDIM1 + i] * x[NDIM1 + j] * x[NDIM1 + k] * x[NDIM1 + m];
									if (PolyOrd2 >= 5)
										for (n = m; n < NDIM2; n++)
											f += par[++iPar] * x[NDIM1 + i] * x[NDIM1 + j] * x[NDIM1 + k] * x[NDIM1 + m] * x[NDIM1 + n];
								}
						}
				}
		}
	return f;
}

void FitMoq()
{
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "MigradImproved");
	gStyle->SetOptStat("nemri");
	// gStyle->SetOptStat(0);
	gStyle->SetStatFormat(".6f");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1111);
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);

	const int MINPOIN = 100, N = 200;
	const double NSGM = 3;

	int i, j, k, m, n;
	double moqErrSys = 0;

	int nPar1 = 1;
	for (i = NDIM1 +1+ PolyOrd1; i > NDIM1 +1; i--)
		nPar1 *= i;
	for (j = PolyOrd1; j > 0; j--)
		nPar1 /= j;
	int nPar2 = 1;
	for (i = NDIM2-1 + PolyOrd2; i > NDIM2-1; i--)
		nPar2 *= i;
	for (j = PolyOrd2; j > 0; j--)
		nPar2 /= j;
	int nPar = nPar1 + nPar2 - 1;
	// int nPar = PolyOrd1 * NDIM1 + 1;
	cout << "nPar = " << nPar << endl;

	const int NAna = 1;
	string sAna[3] = {"CfdTdc", "TacNoClk", "TacClk"};
	int iAna = 0;

	TFile *fMass = new TFile("fMass.root", "UPDATE");
	TFile *fMoq = new TFile("../MachineLearning/fMoq.root", "RECREATE");
	TTree *tMoq[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
		tMoq[iAna] = new TTree(("tMoq_" + sAna[iAna]).c_str(), ("tree of " + sAna[iAna] + " for MoQ fit").c_str());

	int isForCal[N] = {0}, ZOrig[N] = {0}, AOrig[N] = {0}, QOrig[N] = {0};
	double MoQAME[N][2] = {0}, tofPivot[N] = {0};
	int iNucl = 0, nNucl = 0;
	string strRead, strNucl[N];
	for (i = 0; i < N; i++)
		strNucl[iNucl] = "";
	ifstream fCutInfo("fCutInfo.dat");
	while (!fCutInfo.eof() && fCutInfo.peek() != EOF)
	{
		if (fCutInfo.peek() != '#')
		{
			fCutInfo >> isForCal[nNucl] >> strRead >> strNucl[nNucl] >> ZOrig[nNucl] >> AOrig[nNucl] >> QOrig[nNucl] >> strRead >> strRead >> strRead >> strRead >> MoQAME[nNucl][0] >> MoQAME[nNucl][1] >> tofPivot[nNucl];
			getline(fCutInfo, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo, strRead);
	}
	fCutInfo.close();

	FILE *fFitMoq = fopen("fFitMoq.dat", "w+");

	TFile *fRoot = new TFile("tof-run-150--153__270--385.root");
	TTree *tRoot;
	fRoot->GetObject("tTOFDetail", tRoot);

	// TFile *fRoot = new TFile("pid-run-150--153__270--385.root");
	// TTree *tRoot;
	// fRoot->GetObject("tPid", tRoot);

	tRoot->SetEstimate(-1);

	// string strFact[NDIM]={"TOF", "xMcp", "yMcp", "qPMTS800_0", "qPMTS800_1", "qPMTS800_2", "qPMTS800_3", "qPMTA1900_0", "qPMTA1900_1", "qPMTA1900_2", "qPMTA1900_3", "Q"};
	// string strFactUnit[NDIM]={" [ns]", " [mm]", " [mm]", " [ch.]", " [ch.]", " [ch.]", " [ch.]", " [ch.]", " [ch.]", " [ch.]", " [ch.]", ""};

	// string strFact[NDIM]={"TOF", "xMcp", "yMcp", "qPMTS800_1", "qPMTS800_2", "qPMTS800_3", "qPMTS800_4", "Q", "Z", "A"};
	// string strFactUnit[NDIM]={" [ns]", " [mm]", " [mm]", " [ch.]", " [ch.]", " [ch.]", " [ch.]", "", "", ""};

	// string strFact[NDIM]={"TOF", "xMcp", "yMcp", "qPMTS800_1", "qPMTS800_2", "qPMTS800_3", "qPMTS800_4", "AoQ"};
	// string strFactUnit[NDIM]={" [ns]", " [mm]", " [mm]", " [ch.]", " [ch.]", " [ch.]", " [ch.]", ""};

	// string strFact[NDIM] = {"TOF", "xMcp", "yMcp", "xCRDC1", "yCRDC1", "xCRDC2", "yCRDC2", "Q", "Z", "A"};
	// string strFactUnit[NDIM] = {" [ns]", " [mm]", " [mm]", " [pad]", " [ch]", " [pad]", " [ch]", "", "", ""};

	// string strFact[NDIM] = {"TOF", "xMcp", "yMcp", "xCRDC1", "yCRDC1", "xCRDC2", "yCRDC2"};
	// string strFactUnit[NDIM] = {" [ns]", " [mm]", " [mm]", " [pad]", " [ch]", " [pad]", " [ch]"};

	// string strFact[NDIM] = {"TOF", "xMcp", "xCRDC1", "xCRDC2"};
	// string strFactUnit[NDIM] = {" [ns]", " [mm]", " [pad]", " [pad]"};

	// string strFact[NDIM] = {"TOF", "xMcp", "xCRDC1", "xCRDC2", "Q", "Z", "A"};
	// string strFactUnit[NDIM] = {" [ns]", " [mm]", " [pad]", " [pad]", "", "", ""};

	// string strFact[NDIM] = {"TOF", "xMcp", "xCRDC1", "xCRDC2"};
	// string strFactUnit[NDIM] = {" [ns]", " [mm]", " [pad]", " [pad]"};

	// string strFact[NDIM] = {"TOF", "xMcp", "yMcp", "xQS800", "yQS800", "Q", "A", "Z"};
	// string strFactUnit[NDIM] = {" [ns]", " [mm]", " [mm]", "", "", "", "", ""};

	string strFact[NDIM] = {"TOF"};
	string strFactUnit[NDIM] = {" [ns]"};

	for (iAna = 0; iAna < NAna; iAna++)
	{
		printf("\n********************Processing data from electronics setup of %s********************\n", sAna[iAna].c_str());
		fprintf(fFitMoq, "\n********************Processing data from electronics setup of %s********************\n", sAna[iAna].c_str());

		// string sDraw="ZPid:APid:QPid:AoQPid:MoQPid[0]:MoQPid[1]:tof_"+sAna[iAna]+":xMCP:yMCP:egyPMT[0][0]:egyPMT[0][1]:egyPMT[0][2]:egyPMT[0][3]";

		// string sDraw = "ZPid:APid:QPid:AoQPid:MoQPid[0]:MoQPid[1]:tof_" + sAna[iAna] + ":xMCP:yMCP:xCrdc[0]:yCrdc[0]:xCrdc[1]:yCrdc[1]"; //for "pid-run-150--153__270--385.root" file and "tPid" tree
		// string sCut = "isForCal>-1 && good_" + sAna[iAna] + "==1 && ZPid>37";

		// string sDraw = "ZPid:APid:QPid:AoQPid:MoQPid[0]:MoQPid[1]:tof_" + sAna[iAna] + ":xMCP:xCrdc[0]:xCrdc[1]"; //for "pid-run-150--153__270--385.root" file and "tPid" tree
		// string sCut = "isForCal>-1 && good_" + sAna[iAna] + "==1 && ZPid>37";

		// string sDraw = "ZPid:APid:QPid:AoQPid:MoQPid[0]:MoQPid[1]:tof_" + sAna[iAna] + ":xMCP:yMCP:xPlaQ[0]:yPlaQ[0]";  //for "pid-run-150--153__270--385.root" file and "tPid" tree
		// string sCut = "isForCal>-1 && good_" + sAna[iAna] + "==1 && ZPid>37";

		string sDraw = "Z:A:Q:AoQ:MoQ[0]:MoQ[1]:tof2XsCor_" + sAna[iAna];	//for "tof-run-150--153__270--385.root" file and "tTOFDetail" tree
		string sCut = "good_" + sAna[iAna] + "==1 && Z>37";					//for "tof-run-150--153__270--385.root" file and "tTOFDetail" tree

		int nData = tRoot->Draw(sDraw.c_str(), sCut.c_str(), "goff");
		double *ZVal = tRoot->GetVal(0);
		double *AVal = tRoot->GetVal(1);
		double *QVal = tRoot->GetVal(2);
		double *MoQ = tRoot->GetVal(4);
		double *MoQErr = tRoot->GetVal(5);
		double *txyzVal[NDIM];
		for (i = 0; i < NDIM1; i++)
			txyzVal[i] = tRoot->GetVal(6 + i);
		txyzVal[NDIM1 + 0] = tRoot->GetVal(2);
		txyzVal[NDIM1 + 1] = tRoot->GetVal(1);
		txyzVal[NDIM1 + 2] = tRoot->GetVal(0);
		// txyzVal[NDIM1 + 3] = tRoot->GetVal(3);

		//setup the tree of filtered varialbes and MoQ
		double x[NDIM] = {0}, y[2] = {0};
		tMoq[iAna]->Branch("x", x, ("x[" + to_string(NDIM) + "]/D").c_str());
		tMoq[iAna]->Branch("y", y, "y[2]/D");

		int iData = 0, iPar = 0;

		TGraph *grMoqTxyzNucl[NDIM], *grMoqTxyz[NDIM];
		for (i = 0; i < NDIM; i++)
		{
			grMoqTxyzNucl[i] = new TGraph();
			grMoqTxyz[i] = new TGraph();
		}
		TGraph *grMoqErrQ = new TGraph();

		double txyzMean[NDIM1] = {0}, txyzRms[NDIM1] = {0}, txyzRng[NDIM1][2] = {0};
		int nPts = 0, nSel = 0;
		bool good = false;
		nPts = 0;
		bool goodTxyzVal = false;
		for (iNucl = 0; iNucl < nNucl; iNucl++)
			if (isForCal[iNucl] == 1)
			{
				for (i = 0; i < NDIM; i++)
					grMoqTxyzNucl[i]->Set(0);
				memset(txyzMean, 0, sizeof(txyzMean));
				memset(txyzRms, 0, sizeof(txyzRms));
				memset(txyzRng, 0, sizeof(txyzRng));
				nSel = 0;

				for (iData = 0; iData < nData; iData++)
				{
					goodTxyzVal = abs(ZVal[iData] - ZOrig[iNucl]) < 1E-3 && abs(AVal[iData] - AOrig[iNucl]) < 1E-3 && abs(QVal[iData] - QOrig[iNucl]) < 1E-3;
					if (!goodTxyzVal)
						continue;

					for (i = 0; i < NDIM; i++)
					{
						goodTxyzVal = abs(txyzVal[i][iData]) > 0;
						if (!goodTxyzVal)
							break;
					}
					if (goodTxyzVal)
					{
						for (i = 0; i < NDIM; i++)
							grMoqTxyzNucl[i]->SetPoint(nSel, txyzVal[i][iData], MoQ[iData]);
						nSel++;
					}
				}
				if (nSel < MINPOIN)
					continue;

				for (i = 0; i < NDIM1; i++)
				{
					txyzMean[i] = grMoqTxyzNucl[i]->GetMean(1);
					txyzRms[i] = grMoqTxyzNucl[i]->GetRMS(1);
					txyzRng[i][0] = txyzMean[i] - NSGM * txyzRms[i];
					txyzRng[i][1] = txyzMean[i] + NSGM * txyzRms[i];
				}

				for (iData = 0; iData < nData; iData++)
				{
					goodTxyzVal = abs(ZVal[iData] - ZOrig[iNucl]) < 1E-3 && abs(AVal[iData] - AOrig[iNucl]) < 1E-3 && abs(QVal[iData] - QOrig[iNucl]) < 1E-3;
					if (!goodTxyzVal)
						continue;
					for (i = 0; i < NDIM; i++)
					{
						goodTxyzVal = abs(txyzVal[i][iData]) > 0;
						if (!goodTxyzVal)
							break;
					}
					if (goodTxyzVal)
					{
						good = false;
						for (i = 0; i < NDIM1; i++)
						{
							good = txyzVal[i][iData] > txyzRng[i][0] && txyzVal[i][iData] < txyzRng[i][1];
							if (!good)
								break;
						}
						if (good)
						{
							for (i = 0; i < NDIM; i++)
								grMoqTxyz[i]->SetPoint(nPts, txyzVal[i][iData], MoQ[iData]);
							grMoqErrQ->SetPoint(nPts, QVal[iData], MoQErr[iData]);
							nPts++;
						}
					}
				}
			}
		if (nPts < nPar * MINPOIN)
			continue;
		int ndf = nPts - nPar;

		// ROOT::Fit::BinData dataMoq(nData, NDIM, ROOT::Fit::BinData::kNoError);
		ROOT::Fit::BinData dataMoq(nData, NDIM);
		double xTxyz[NDIM] = {0}, yMoq = 0, yMoqErr = 0;
		for (iData = 0; iData < nPts; iData++)
		{
			memset(xTxyz, 0, sizeof(xTxyz));
			memset(x, 0, sizeof(x));
			memset(y, 0, sizeof(y));
			yMoq = 0;
			yMoqErr = 0;
			for (i = 0; i < NDIM; i++)
			{
				xTxyz[i] = *(grMoqTxyz[i]->GetX() + iData);
				x[i] = xTxyz[i];
			}
			yMoq = *(grMoqTxyz[0]->GetY() + iData);
			yMoqErr = *(grMoqErrQ->GetY() + iData);
			y[0] = yMoq;
			y[1] = yMoqErr * sqrt(nPts);
			// dataMoq.Add(xTxyz, yMoq, sqrt(nPts));
			dataMoq.Add(xTxyz, yMoq, yMoqErr * sqrt(nPts));
			// dataMoq.Add(xTxyz, yMoq);
			tMoq[iAna]->Fill();
		}
		TF1 fcnMoq("fcnMoq", myFunc, 0, 10, nPar);
		ROOT::Math::WrappedMultiTF1 fitFunc(fcnMoq, NDIM);
		ROOT::Fit::Fitter fitter;
		fitter.SetFunction(fitFunc, true);
		for (i = 0; i < nPar; i++)
			fitter.Config().ParSettings(i).SetValue(1);
		// fitter.LeastSquareFit(dataMoq);
		fitter.LinearFit(dataMoq);
		ROOT::Fit::FitResult resultFit = fitter.Result();
		resultFit.Print(std::cout);

		double chi2 = resultFit.Chi2();
		double *parVal = new double[nPar];
		for (iPar = 0; iPar < nPar; iPar++)
			parVal[iPar] = resultFit.Parameter(iPar);

		TGraphErrors *grMoQ[NDIM];
		TGraph *grResMoq[NDIM], *grResM[NDIM];
		for (i = 0; i < NDIM; i++)
		{
			grMoQ[i] = new TGraphErrors();
			grResMoq[i] = new TGraph();
			grResM[i] = new TGraph();

			grMoQ[i]->SetNameTitle(("MoQ__" + strFact[i] + "_" + sAna[iAna]).c_str(), ("m/q_vs_" + strFact[i] + ";" + strFact[i] + strFactUnit[i] + ";m/q [keV/e]").c_str());

			grResMoq[i]->SetNameTitle(("Res_Moq__" + strFact[i] + "_" + sAna[iAna]).c_str(), ("#Deltam/q_vs_" + strFact[i] + ";" + strFact[i] + strFactUnit[i] + ";m/q_{fit}-m/q_{AME16} [keV/e]").c_str());

			grResM[i]->SetNameTitle(("Res_M__" + strFact[i] + "_" + sAna[iAna]).c_str(), ("#Deltam_vs_" + strFact[i] + ";" + strFact[i] + strFactUnit[i] + ";m_{fit}-m_{AME16} [keV]").c_str());
		}
		double moqRes = 0, qVal = 0, mRes = 0;
		for (iData = 0; iData < nPts; iData++)
		{
			memset(xTxyz, 0, sizeof(xTxyz));
			yMoq = 0;
			yMoqErr = 0;
			qVal = 0;
			for (j = 0; j < NDIM; j++)
				xTxyz[j] = *(grMoqTxyz[j]->GetX() + iData);
			yMoq = *(grMoqTxyz[0]->GetY() + iData);
			grMoqErrQ->GetPoint(iData, qVal, yMoqErr);

			moqRes = yMoq - fcnMoq(xTxyz, parVal);

			mRes = qVal * moqRes;

			for (j = 0; j < NDIM; j++)
			{
				grMoQ[j]->SetPoint(iData, xTxyz[j], yMoq);
				grMoQ[j]->SetPointError(iData, 0, yMoqErr);
				grResMoq[j]->SetPoint(iData, xTxyz[j], moqRes);
				grResM[j]->SetPoint(iData, xTxyz[j], mRes);
			}
		}
		TH1F *hResMoq = new TH1F(("hRes_Moq__" + sAna[iAna]).c_str(), "#Deltam/q;m/q_{fit}-m/q_{AME16} [keV/e];Counts/(5keV/e)", 1.2 * grResMoq[0]->GetRMS(2), grResMoq[0]->GetMean(2) - 3 * grResMoq[0]->GetRMS(2), grResMoq[0]->GetMean(2) + 3 * grResMoq[0]->GetRMS(2));
		TH1F *hResM = new TH1F(("hRes_M__" + sAna[iAna]).c_str(), "#Deltam;m_{fit}-m_{AME16} [keV];Counts/100keV", 0.06 * grResM[0]->GetRMS(2), grResM[0]->GetMean(2) - 3 * grResM[0]->GetRMS(2), grResM[0]->GetMean(2) + 3 * grResM[0]->GetRMS(2));

		for (iData = 0; iData < nPts; iData++)
		{
			hResMoq->Fill(*(grResMoq[0]->GetY() + iData));
			hResM->Fill(*(grResM[0]->GetY() + iData));
		}
		printf("****************************************************************************************************\n");
		printf("Chi2/NDF=%.0f/%d=%.0f;  Average residual:  m/q: %.1f +- %.1f;  m: %.1f +- %.1f\n", chi2, ndf, chi2 / ndf, grResMoq[0]->GetMean(2), grResMoq[0]->GetRMS(2), grResM[0]->GetMean(2), grResM[0]->GetRMS(2));
		printf("****************************************************************************************************\n");
		fprintf(fFitMoq, "****************************************************************************************************\n");
		fprintf(fFitMoq, "Chi2/NDF=%.0f/%d=%.0f;  Average residual:  m/q: %.1f +- %.1f;  m: %.1f +- %.1f\n", chi2, ndf, chi2 / ndf, grResMoq[0]->GetMean(2), grResMoq[0]->GetRMS(2), grResM[0]->GetMean(2), grResM[0]->GetRMS(2));
		fprintf(fFitMoq, "****************************************************************************************************\n");

		TCanvas *cResMoq = new TCanvas(("c2ResMoq" + sAna[iAna]).c_str(), ("c2ResMoq" + sAna[iAna]).c_str(), 4800, 2700);
		cResMoq->SetCanvasSize(4800, 2700);
		cResMoq->Divide(4, 3);
		for (i = 0; i < NDIM; i++)
		{
			cResMoq->cd(i + 1);
			grResMoq[i]->Draw("AP*");
		}
		cResMoq->SaveAs(".png");

		fMass->cd();
		cResMoq->Write("", TObject::kOverwrite);

		TGraphErrors *grMexCal = new TGraphErrors();
		grMexCal->SetNameTitle(("grMexCal" + sAna[iAna]).c_str(), "grMexCal");
		TGraphErrors *grMexAME = new TGraphErrors();
		grMexAME->SetNameTitle(("grMexAME" + sAna[iAna]).c_str(), "grMexAME");

		TGraphErrors *grResMexRef = new TGraphErrors();
		grResMexRef->SetNameTitle(("grResMexRef" + sAna[iAna]).c_str(), "grResMexRef;Nuclide (=1000#timesZ+A);ME_{cal}-ME_{AME16} [keV]");
		TGraphErrors *grResMexNo = new TGraphErrors();
		grResMexNo->SetNameTitle(("grResMexNo" + sAna[iAna]).c_str(), "grResMexNo;Nuclide (=1000#timesZ+A);ME_{cal}-ME_{AME16} [keV]");

		int iPt = 0, iPtRef = 0, iPtNo = 0;
		printf("\n****************************************************************************************************\n");
		for (iNucl = 0; iNucl < nNucl; iNucl++)
		// if(ZOrig[iNucl]==QOrig[iNucl]&&ZOrig[iNucl]>38)
		{
			for (i = 0; i < NDIM; i++)
			{
				grMoqTxyzNucl[i]->Set(0);
				grMoqTxyz[i]->Set(0);
			}
			grMoqErrQ->Set(0);

			memset(txyzMean, 0, sizeof(txyzMean));
			memset(txyzRms, 0, sizeof(txyzRms));
			memset(txyzRng, 0, sizeof(txyzRng));
			nSel = 0;
			for (iData = 0; iData < nData; iData++)
			{
				goodTxyzVal = abs(ZVal[iData] - ZOrig[iNucl]) < 1E-3 && abs(AVal[iData] - AOrig[iNucl]) < 1E-3 && abs(QVal[iData] - QOrig[iNucl]) < 1E-3;
				if (!goodTxyzVal)
					continue;
				for (i = 0; i < NDIM; i++)
				{
					goodTxyzVal = abs(txyzVal[i][iData]) > 0;
					if (!goodTxyzVal)
						break;
				}
				if (goodTxyzVal)
				{
					for (i = 0; i < NDIM; i++)
						grMoqTxyzNucl[i]->SetPoint(nSel, txyzVal[i][iData], MoQ[iData]);
					nSel++;
				}
			}
			for (i = 0; i < NDIM1; i++)
			{
				txyzMean[i] = grMoqTxyzNucl[i]->GetMean(1);
				txyzRms[i] = grMoqTxyzNucl[i]->GetRMS(1);
				txyzRng[i][0] = txyzMean[i] - NSGM * txyzRms[i];
				txyzRng[i][1] = txyzMean[i] + NSGM * txyzRms[i];
			}
			nPts = 0;
			for (iData = 0; iData < nData; iData++)
			{
				goodTxyzVal = abs(ZVal[iData] - ZOrig[iNucl]) < 1E-3 && abs(AVal[iData] - AOrig[iNucl]) < 1E-3 && abs(QVal[iData] - QOrig[iNucl]) < 1E-3;
				if (!goodTxyzVal)
					continue;
				for (i = 0; i < NDIM; i++)
				{
					goodTxyzVal = abs(txyzVal[i][iData]) > 0;
					if (!goodTxyzVal)
						break;
				}
				if (goodTxyzVal)
				{
					good = false;
					for (i = 0; i < NDIM1; i++)
					{
						good = txyzVal[i][iData] > txyzRng[i][0] && txyzVal[i][iData] < txyzRng[i][1];
						if (!good)
							break;
					}
					if (good)
					{
						for (i = 0; i < NDIM; i++)
							grMoqTxyz[i]->SetPoint(nPts, txyzVal[i][iData], MoQ[iData]);
						grMoqErrQ->SetPoint(nPts, QVal[iData], MoQErr[iData]);
						nPts++;
					}
				}
			}
			if (nPts < MINPOIN)
				continue;

			double *xVal = new double[nPts * NDIM];
			for (i = 0; i < nPts * NDIM; i++)
				xVal[i] = 0;
			for (iData = 0; iData < nPts; iData++)
				for (k = 0; k < NDIM; k++)
					xVal[NDIM * iData + k] = grMoqTxyz[k]->GetX()[iData];
			double *moqErrPar = new double[nPts];
			// resultFit.GetConfidenceIntervals(nPts, NDIM, 1, xVal, moqErrPar, 0.682689492137086, false);
			delete[] xVal;

			TGraph *grMass = new TGraph();
			double mass = 0;
			for (iData = 0; iData < nPts; iData++)
			{
				memset(xTxyz, 0, sizeof(xTxyz));
				for (j = 0; j < NDIM; j++)
					xTxyz[j] = *(grMoqTxyz[j]->GetX() + iData);

				mass = QOrig[iNucl] * fcnMoq(xTxyz, parVal);
				grMass->SetPoint(iData, mass, mass);
			}
			double meanMass = grMass->GetMean(1);
			double rmsMass = grMass->GetRMS(1);
			delete grMass;

			TGraph *grCalMExMass = new TGraph();
			TGraph *grMassWt = new TGraph();
			TGraph *grMassMex2 = new TGraph();
			double mex = 0;
			int nMassPts = 0;
			for (iData = 0; iData < nPts; iData++)
			{
				memset(xTxyz, 0, sizeof(xTxyz));
				for (j = 0; j < NDIM; j++)
					xTxyz[j] = *(grMoqTxyz[j]->GetX() + iData);

				mass = QOrig[iNucl] * fcnMoq(xTxyz, parVal);
				if (abs(mass - meanMass) < NSGM * rmsMass)
				{
					mex = mass - (ZOrig[iNucl] - QOrig[iNucl]) * 511 - 931494.0954 * AOrig[iNucl] + 510.9989461 * ZOrig[iNucl] - 0.0144381 * pow(ZOrig[iNucl], 2.39) - 1.55468E-9 * pow(ZOrig[iNucl], 5.35);
					grCalMExMass->SetPoint(nMassPts, mass, mex);
					// grMassWt->SetPoint(nMassPts, 1 / pow(QOrig[iNucl] * moqErrPar[iData], 2), 1 / pow(QOrig[iNucl] * moqErrPar[iData], 2));
					grMassMex2->SetPoint(nMassPts, mass * mass, mex * mex);
					nMassPts++;
				}
			}
			double massTarg = 0, mExTarg = 0, mExAME[2] = {0};
			double massErr[3] = {0}; //[3]: [0] total error; [1] systematic error; [2] statistical+fit errors;
			double massResolStat = 0, massResol = 0;

			// massTarg = TMath::Mean(nMassPts, grCalMExMass->GetX(), grMassWt->GetX());
			massTarg = grCalMExMass->GetMean(1);

			// mExTarg = TMath::Mean(nMassPts, grCalMExMass->GetY(), grMassWt->GetY());
			mExTarg = grCalMExMass->GetMean(2);

			// massResolStat = sqrt((TMath::Mean(nMassPts, grMassMex2->GetX(), grMassWt->GetX()) - massTarg * massTarg) * 1.0 * nMassPts / (nMassPts - 1));
			massResolStat = grCalMExMass->GetRMS(1) * sqrt(1.0 * nMassPts / (nMassPts - 1));
			// massResolStat = sqrt(nMassPts) / sqrt(grMassWt->GetMean(1));

			mExAME[0] = QOrig[iNucl] * MoQAME[iNucl][0] - (ZOrig[iNucl] - QOrig[iNucl]) * 511 - 931494.0954 * AOrig[iNucl] + 510.9989461 * ZOrig[iNucl] - 0.0144381 * pow(ZOrig[iNucl], 2.39) - 1.55468E-9 * pow(ZOrig[iNucl], 5.35);
			mExAME[1] = QOrig[iNucl] * MoQAME[iNucl][1];

			massErr[1] = QOrig[iNucl] * moqErrSys;

			massErr[2] = massResolStat / sqrt(nMassPts);

			// cout << "\nmassResolStat= " << massResolStat << ",  massErr[2]*sqrt(N)= " << massErr[2] * sqrt(nMassPts) << endl;

			massErr[0] = sqrt(massErr[1] * massErr[1] + massErr[2] * massErr[2]);

			massResol = sqrt(massErr[1] * massErr[1] + massResolStat * massResolStat);

			delete[] moqErrPar;
			delete grCalMExMass;
			delete grMassWt;
			delete grMassMex2;

			grMexCal->SetPoint(iPt, ZOrig[iNucl] * 1000 + AOrig[iNucl], mExTarg);
			grMexCal->SetPointError(iPt, 0, massErr[0]);

			grMexAME->SetPoint(iPt, ZOrig[iNucl] * 1000 + AOrig[iNucl], mExAME[0]);
			grMexAME->SetPointError(iPt, 0, mExAME[1]);

			if (isForCal[iNucl] == 1)
			{
				grResMexRef->SetPoint(iPtRef, ZOrig[iNucl] * 1000 + AOrig[iNucl], mExTarg - mExAME[0]);
				grResMexRef->SetPointError(iPtRef, 0, sqrt(mExAME[1] * mExAME[1] + massErr[0] * massErr[0]));
				iPtRef++;
			}
			if (isForCal[iNucl] == 0)
			{
				grResMexNo->SetPoint(iPtNo, ZOrig[iNucl] * 1000 + AOrig[iNucl], mExTarg - mExAME[0]);
				grResMexNo->SetPointError(iPtNo, 0, sqrt(mExAME[1] * mExAME[1] + massErr[0] * massErr[0]));
				iPtNo++;
			}
			printf("\n%s: isCalib=%d, Z=%d, A=%d, Q=%d,  m_cal=%.1f +- %.1f +- %.1f=%.1f +- %.1f keV, sigma_m=%.1f +- %.1f=%.1f keV (%.1E),  ME_cal=%.1f +- %.1f keV, ME_AME=%.1f +- %.1f keV\n", strNucl[iNucl].c_str(), isForCal[iNucl], ZOrig[iNucl], AOrig[iNucl], QOrig[iNucl], massTarg, massErr[1], massErr[2], massTarg, massErr[0], massErr[1], massResolStat, massResol, massResol / massTarg, mExTarg, massErr[0], mExAME[0], mExAME[1]);

			fprintf(fFitMoq, "\n%s: isCalib=%d, Z=%d, A=%d, Q=%d,  m_cal=%.1f +- %.1f +- %.1f=%.1f +- %.1f keV, sigma_m=%.1f +- %.1f=%.1f keV (%.1E),  ME_cal=%.1f +- %.1f keV, ME_AME=%.1f +- %.1f keV\n", strNucl[iNucl].c_str(), isForCal[iNucl], ZOrig[iNucl], AOrig[iNucl], QOrig[iNucl], massTarg, massErr[1], massErr[2], massTarg, massErr[0], massErr[1], massResolStat, massResol, massResol / massTarg, mExTarg, massErr[0], mExAME[0], mExAME[1]);

			iPt++;
		}
		printf("****************************************************************************************************\n");
		fprintf(fFitMoq, "****************************************************************************************************\n");

		TCanvas *cMex = new TCanvas(("c2Mex" + sAna[iAna]).c_str(), ("c2Mex" + sAna[iAna]).c_str(), 2400, 1800);
		cMex->SetCanvasSize(2400, 1800);
		cMex->Divide(2, 2);

		cMex->cd(1);
		hResMoq->Draw();
		hResMoq->GetYaxis()->SetMaxDigits(3);
		cMex->cd(2);
		hResM->Draw();
		hResM->GetYaxis()->SetMaxDigits(3);

		cMex->cd(3);
		grMexCal->SetLineColor(kRed);
		grMexCal->SetMarkerColor(kRed);

		grMexAME->SetLineColor(kBlue);
		grMexAME->SetMarkerColor(kBlue);

		TMultiGraph *mgrMex = new TMultiGraph(("mgrMex" + sAna[iAna]).c_str(), "mgrMex;Nuclide (=1000#timesZ+A);Mass excess [keV]");
		mgrMex->Add(grMexCal);
		mgrMex->Add(grMexAME);
		mgrMex->GetXaxis()->SetMaxDigits(3);
		mgrMex->GetYaxis()->SetMaxDigits(3);
		mgrMex->Draw("AP*");
		gPad->BuildLegend(0.13, 0.75, 0.4, 0.9)->SetTextSize(0.04);

		cMex->cd(4);
		grResMexNo->SetMarkerColor(kRed);
		grResMexNo->SetLineColor(kRed);
		grResMexNo->SetMarkerSize(2);
		grResMexNo->SetLineWidth(2);
		TMultiGraph *mgrResMex = new TMultiGraph(("mgrResMex" + sAna[iAna]).c_str(), ";Nuclide (=1000#timesZ+A);ME_{cal}-ME_{AME16} [keV]");
		mgrResMex->Add(grResMexRef);
		mgrResMex->Add(grResMexNo);
		mgrResMex->GetXaxis()->SetMaxDigits(3);
		mgrResMex->GetYaxis()->SetMaxDigits(3);
		mgrResMex->DrawClone("AP*");

		TCanvas *cResM = new TCanvas(("c2ResM" + sAna[iAna]).c_str(), ("c2ResM" + sAna[iAna]).c_str());
		mgrResMex->Draw("AP*");

		cResM->SaveAs(".png");
		cMex->SaveAs(".png");
		fMass->cd();
		cMex->Write("", TObject::kOverwrite);
		cResM->Write("", TObject::kOverwrite);

		fMoq->cd();
		tMoq[iAna]->Write();

		delete mgrResMex;
		delete mgrMex;
		delete cResM;
		delete cMex;
		delete cResMoq;
		for (i = 0; i < NDIM; i++)
		{
			delete grMoQ[i];
			delete grResMoq[i];
			delete grResM[i];
			delete grMoqTxyz[i];
			delete grMoqTxyzNucl[i];
		}
		delete grMoqErrQ;
		delete[] parVal;
		delete hResMoq;
	}
	delete fRoot;
	delete fMoq;
	delete fMass;
	fclose(fFitMoq);
}

int main(int argc, char **argv)
{
	FitMoq();
	return 0;
}