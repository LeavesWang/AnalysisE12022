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
#include "TMultiGraph.h"
#include "TPad.h"
#include "TLegend.h"

void FitMoqTof()
{
	// gStyle->SetOptStat("nemri");
	gStyle->SetOptStat(0);
	gStyle->SetStatFormat(".6f");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1111);
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);

	const int MINPOIN = 100;
	const int NORD = 3;
	const int NPAR = (NORD + 2) * (NORD + 1) / 2;
	// const int NORD = 0;
	// const int NPAR = 4;

	const int NAna = 1;
	string sAna[3] = {"CfdTdc", "TacNoClk", "TacClk"};
	int iAna;

	int i, j;

	TFile *fMass = new TFile("fMass.root", "UPDATE");

	FILE *fFitMoqTof = fopen("fFitMoqTof.dat", "w+");

	TFile *fMoqTof = new TFile("../MachineLearning/fMoqTof.root", "RECREATE");
	TTree *tMoqTof[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
		tMoqTof[iAna] = new TTree(("tMoqTof_" + sAna[iAna]).c_str(), ("tree of " + sAna[iAna] + " for MoQ--ToF fit").c_str());

	TFile *fTof = new TFile("/home/kailong/ExpData/Jul2018/Pid/tof-run-150--153__270--385.root");
	// TFile *fTof = new TFile("/home/kailong/ExpData/Jul2018/Pid/tofGlob-run-150--153__270--385.root");
	TTree *tTOF;
	fTof->GetObject("tTOF", tTOF);
	tTOF->SetEstimate(-1);

	for (iAna = 0; iAna < NAna; iAna++)
	{
		printf("\n********************Processing data from electronics setup of %s********************\n", sAna[iAna].c_str());
		fprintf(fFitMoqTof, "\n********************Processing data from electronics setup of %s********************\n", sAna[iAna].c_str());

		string sDraw = "isCal:Z:A:Q:AoQ:MoQ[0]:MoQ[1]:tofXsCorr_" + sAna[iAna] + "[0]:tofXsCorr_" + sAna[iAna] + "[1]:cnts_" + sAna[iAna];
		string sCut = "good_" + sAna[iAna] + " && abs(tofXsCorr_" + sAna[iAna] + "[0])>0 && Z>37";

		int nData = tTOF->Draw(sDraw.c_str(), sCut.c_str(), "goff");

		double *calVal = tTOF->GetVal(0);
		double *ZVal = tTOF->GetVal(1);
		double *AVal = tTOF->GetVal(2);
		double *QVal = tTOF->GetVal(3);
		double *MoQVal = tTOF->GetVal(5);
		double *MoQErrVal = tTOF->GetVal(6);
		double *tofVal = tTOF->GetVal(7);
		double *tofErrVal = tTOF->GetVal(8);
		double *cntsVal = tTOF->GetVal(9);

		double tof[2] = {0}, q = 0, moq[2] = {0};
		int isCalib = 0;
		tMoqTof[iAna]->Branch("tof", tof, "tof[2]/D");
		tMoqTof[iAna]->Branch("q", &q, "q/D");
		tMoqTof[iAna]->Branch("moq", moq, "moq[2]/D");
		tMoqTof[iAna]->Branch("isCalib", &isCalib, "isCalib/I");

		TGraph2DErrors *grMoqTofZVal = new TGraph2DErrors();
		TGraph *grAQVal = new TGraph();
		int iData = 0, nCalib = 0;
		for (iData = 0; iData < nData; iData++)
		{
			memset(tof, 0, sizeof(tof));
			q = 0;
			memset(moq, 0, sizeof(moq));
			isCalib = 0;
			if (abs(calVal[iData] - 1) < 1E-3)
			{
				grMoqTofZVal->SetPoint(nCalib, tofVal[iData], QVal[iData], MoQVal[iData]);
				grMoqTofZVal->SetPointError(nCalib, tofErrVal[iData], 0, MoQErrVal[iData]);
				grAQVal->SetPoint(nCalib, QVal[iData], AVal[iData]);
				nCalib++;
				isCalib = 1;
			}
			tof[0] = tofVal[iData];
			tof[1] = tofErrVal[iData];
			q = QVal[iData];
			moq[0] = MoQVal[iData];
			moq[1] = MoQErrVal[iData];
			tMoqTof[iAna]->Fill();
		}
		double *tofCalib = grMoqTofZVal->GetX();
		double *zCalib = grMoqTofZVal->GetY();
		double *moqCalib = grMoqTofZVal->GetZ();
		double *tofErrCalib = grMoqTofZVal->GetEX();
		double *moqErrCalib = grMoqTofZVal->GetEZ();
		double *qCalib = grAQVal->GetX();
		double *aCalib = grAQVal->GetY();

		TGraph2DErrors *grMoq_Tof_Z = new TGraph2DErrors();
		grMoq_Tof_Z->SetNameTitle(("MoQ_TOF_Z" + sAna[iAna]).c_str(), ("MoQ_TOF_Z" + sAna[iAna] + ";TOF(ns);ZVal;m/q(keV)").c_str());
		TF2 *fcnMoq_Tof_Q, *difFcnMoq;
		if (NORD == 2) // 2nd order polynomial
		{
			fcnMoq_Tof_Q = new TF2("fcnMoq_Tof_Q", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y", 460, 500, 30, 50);
			difFcnMoq = new TF2("difFcnMoq", "0*[0]+[1]+0*[2]+2*[3]*x+[4]*y+0*[5]", 460, 500, 30, 50);
		}
		if (NORD == 3) // 3rd polynomial
		{
			fcnMoq_Tof_Q = new TF2("fcnMoq_Tof_Q", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y+[6]*x*x*x+[7]*x*x*y+[8]*x*y*y+[9]*y*y*y", 460, 500, 30, 50);
			difFcnMoq = new TF2("difFcnMoq", "0*[0]+[1]+0*[2]+2*[3]*x+[4]*y+0*[5]+3*[6]*x*x+2*[7]*x*y+[8]*y*y+0*[9]", 460, 500, 30, 50);
		}
		if (NORD == 4) // 4th polynomial
		{
			fcnMoq_Tof_Q = new TF2("fcnMoq_Tof_Q", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y+[6]*x*x*x+[7]*x*x*y+[8]*x*y*y+[9]*y*y*y+[10]*x*x*x*x+[11]*x*x*x*y+[12]*x*x*y*y+[13]*x*y*y*y+[14]*y*y*y*y", 460, 500, 30, 50);
			difFcnMoq = new TF2("difFcnMoq", "0*[0]+[1]+0*[2]+2*[3]*x+[4]*y+0*[5]+3*[6]*x*x+2*[7]*x*y+[8]*y*y+0*[9]+4*[10]*x*x*x+3*[11]*x*x*y+2*[12]*x*y*y+[13]*y*y*y+0*[14]", 460, 500, 30, 50);
		}
		if (NORD == 5) // 5th polynomial
		{
			fcnMoq_Tof_Q = new TF2("fcnMoq_Tof_Q", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y+[6]*x*x*x+[7]*x*x*y+[8]*x*y*y+[9]*y*y*y+[10]*x*x*x*x+[11]*x*x*x*y+[12]*x*x*y*y+[13]*x*y*y*y+[14]*y*y*y*y+[15]*x*x*x*x*x+[16]*x*x*x*x*y+[17]*x*x*x*y*y+[18]*x*x*y*y*y+[19]*x*y*y*y*y+[20]*y*y*y*y*y", 460, 500, 30, 50);
			difFcnMoq = new TF2("difFcnMoq", "0*[0]+[1]+0*[2]+2*[3]*x+[4]*y+0*[5]+3*[6]*x*x+2*[7]*x*y+[8]*y*y+0*[9]+4*[10]*x*x*x+3*[11]*x*x*y+2*[12]*x*y*y+[13]*y*y*y+0*[14]+5*[15]*x*x*x*x+4*[16]*x*x*x*y+3*[17]*x*x*y*y+2*[18]*x*y*y*y+[19]*y*y*y*y+0*[20]", 460, 500, 30, 50);
		}
		if (NORD == 6) // 6th polynomial
		{
			fcnMoq_Tof_Q = new TF2("fcnMoq_Tof_Q", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y+[6]*x*x*x+[7]*x*x*y+[8]*x*y*y+[9]*y*y*y+[10]*x*x*x*x+[11]*x*x*x*y+[12]*x*x*y*y+[13]*x*y*y*y+[14]*y*y*y*y+[15]*x*x*x*x*x+[16]*x*x*x*x*y+[17]*x*x*x*y*y+[18]*x*x*y*y*y+[19]*x*y*y*y*y+[20]*y*y*y*y*y+[21]*x*x*x*x*x*x+[22]*x*x*x*x*x*y+[23]*x*x*x*x*y*y+[24]*x*x*x*y*y*y+[25]*x*x*y*y*y*y+[26]*x*y*y*y*y*y+[27]*y*y*y*y*y*y", 460, 500, 30, 50);
			difFcnMoq = new TF2("difFcnMoq", "0*[0]+[1]+0*[2]+2*[3]*x+[4]*y+0*[5]+3*[6]*x*x+2*[7]*x*y+[8]*y*y+0*[9]+4*[10]*x*x*x+3*[11]*x*x*y+2*[12]*x*y*y+[13]*y*y*y+0*[14]+5*[15]*x*x*x*x+4*[16]*x*x*x*y+3*[17]*x*x*y*y+2*[18]*x*y*y*y+[19]*y*y*y*y+0*[20]+6*[21]*x*x*x*x*x+5*[22]*x*x*x*x*y+4*[23]*x*x*x*y*y+3*[24]*x*x*y*y*y+2*[25]*x*y*y*y*y+[26]*y*y*y*y*y+0*[27]", 460, 500, 30, 50);
		}
		if (NORD == 0)
		{
			fcnMoq_Tof_Q = new TF2("fcnMoq_Tof_Q", "[0]+sqrt([1]+[2]*x+[3]*x*x)+[4]*y", 460, 500, 30, 50);
			difFcnMoq = new TF2("difFcnMoq", "0*[0]+([2]+2*[3]*x)/(2*sqrt([1]+[2]*x+[3]*x*x))+0*[4]", 460, 500, 30, 50);
		}

		double moqErrSys = 0;
		double *moqErrTot = new double[nCalib];
		double par[NPAR] = {0}, parTem[NPAR] = {0};
		bool goodFit = false;
		TFitResultPtr fitMoq;
		double chi2 = 0;
		int ndf = 0, iter;
		////////using reduced chi-square normalization method to estimate systematic error
		for (iter = 0; iter < 1; iter++)
		// for (iter = 0; iter < 100000; iter++)
		{
			// printf("\n********************iteration=%d********************\n", iter);
			difFcnMoq->SetParameters(par);
			for (j = 0; j < nCalib; j++)
			{
				grMoq_Tof_Z->SetPoint(j, tofCalib[j], zCalib[j], moqCalib[j]);

				moqErrTot[j] = sqrt(moqErrSys * moqErrSys + moqErrCalib[j] * moqErrCalib[j] + pow(difFcnMoq->Eval(tofCalib[j], zCalib[j]) * tofErrCalib[j], 2));

				grMoq_Tof_Z->SetPointError(j, 0, 0, moqErrTot[j]);
			}
			fitMoq = grMoq_Tof_Z->Fit(fcnMoq_Tof_Q, "S");
			chi2 = fitMoq->Chi2();
			ndf = fitMoq->Ndf();
			goodFit = true;
			for (int iPar = 0; iPar < NPAR; iPar++)
			{
				parTem[iPar] = fitMoq->Parameter(iPar);
				goodFit = goodFit && abs((parTem[iPar] - par[iPar]) / parTem[iPar]) < 1E-3;
				// printf("%d %d: par=%f, parTem=%f\n", iter, iPar, par[iPar], parTem[iPar]);
				if (abs((parTem[iPar] - par[iPar]) / parTem[iPar]) > 1E-3)
					par[iPar] = parTem[iPar];
			}
			if (!goodFit)
				continue;
			else
				goodFit = goodFit && abs(chi2 / ndf - 1) < 0.01;
			if (goodFit)
				break;
			moqErrSys = moqErrSys + 0.01;
		}
		for (int iPar = 0; iPar < NPAR; iPar++)
			printf("par_%d=%f\n", iPar, par[iPar]);

		printf("\n**********************************************************************\n");
		printf("After %d iterations: moqErrSys=%.2f, chi2/ndf=%f/%d=%.3f", iter, moqErrSys, chi2, ndf, chi2 / ndf);
		printf("\n**********************************************************************\n");
		fprintf(fFitMoqTof, "\n**********************************************************************\n");
		fprintf(fFitMoqTof, "After %d iterations: moqErrSys=%.2f, chi2/ndf=%f/%d=%.3f", iter, moqErrSys, chi2, ndf, chi2 / ndf);
		fprintf(fFitMoqTof, "\n**********************************************************************\n");

		TCanvas *cMoq = new TCanvas(("c1Moq" + sAna[iAna]).c_str(), ("c1Moq" + sAna[iAna]).c_str());
		cMoq->Divide(2, 2);

		TGraphErrors *resMoq_Tof = new TGraphErrors();
		resMoq_Tof->SetNameTitle(("ResMoq_TOF_" + sAna[iAna]).c_str(), "#Deltam/q_vs_TOF;TOF(ns);m/q_{fit}-m/q_{AME16} [keV/e]");
		TGraphErrors *resMoq_Z = new TGraphErrors();
		resMoq_Z->SetNameTitle(("ResMoq_Z_" + sAna[iAna]).c_str(), "#Deltam/q_vs_Z;Z;m/q_{fit}-m/q_{AME16} [keV/e]");
		TGraph2DErrors *resMoq_Z_AoQ = new TGraph2DErrors();
		resMoq_Z_AoQ->SetNameTitle(("ResMoq_Z_AoQ_" + sAna[iAna]).c_str(), "#Deltam/q_vs_Z_AoQ;Z;AoQ;m/q_{fit}-m/q_{AME16} [keV/e]");

		TGraphErrors *resM_Tof = new TGraphErrors();
		resM_Tof->SetNameTitle(("ResM_TOF_" + sAna[iAna]).c_str(), "#Deltam_vs_TOF;TOF(ns);m_{fit}-m_{AME16} [keV]");
		TGraphErrors *resM_Z = new TGraphErrors();
		resM_Z->SetNameTitle(("ResM_Z_" + sAna[iAna]).c_str(), "#Deltam_vs_Z;Z;m_{fit}-m_{AME16} [keV]");
		TGraph2DErrors *resM_Z_AoQ = new TGraph2DErrors();
		resM_Z_AoQ->SetNameTitle(("ResM_Z_AoQ_" + sAna[iAna]).c_str(), "#Deltam_vs_Z_AoQ;Z;AoQ;m_{fit}-m_{AME16} [keV]");

		double moqRes = 0, moqResErr = 0, mRes = 0, mResErr = 0, aoqRef = 0;
		cMoq->cd(1);
		fcnMoq_Tof_Q->Draw("surf1");
		grMoq_Tof_Z->Draw("same err p0");
		for (j = 0; j < nCalib; j++)
		{
			moqRes = fcnMoq_Tof_Q->Eval(tofCalib[j], zCalib[j]) - moqCalib[j];
			moqResErr = moqErrTot[j];

			mRes = qCalib[j] * moqRes - (zCalib[j] - qCalib[j]) * 511;
			mResErr = qCalib[j] * moqResErr;

			aoqRef = aCalib[j] / qCalib[j];

			resMoq_Tof->SetPoint(j, tofCalib[j], moqRes);
			resMoq_Tof->SetPointError(j, tofErrCalib[j], moqResErr);

			resMoq_Z->SetPoint(j, zCalib[j], moqRes);
			resMoq_Z->SetPointError(j, 0, moqResErr);

			resMoq_Z_AoQ->SetPoint(j, zCalib[j], aoqRef, moqRes);
			resMoq_Z_AoQ->SetPointError(j, 0, 0, moqResErr);

			resM_Tof->SetPoint(j, tofCalib[j], mRes);
			resM_Tof->SetPointError(j, tofErrCalib[j], mResErr);

			resM_Z->SetPoint(j, zCalib[j], mRes);
			resM_Z->SetPointError(j, 0, mResErr);

			resM_Z_AoQ->SetPoint(j, zCalib[j], aoqRef, mRes);
			resM_Z_AoQ->SetPointError(j, 0, 0, mResErr);
		}
		delete[] moqErrTot;

		cMoq->cd(2);
		resMoq_Tof->Draw("AP*");
		cMoq->cd(3);
		resMoq_Z->Draw("AP*");
		// cMoq->cd(4);
		// resM_Z_AoQ->Draw("surf1");
		// resM_Z_AoQ->Draw("same p0 err");

		TGraphErrors *grMexCal = new TGraphErrors();
		grMexCal->SetNameTitle(("grMexCal" + sAna[iAna]).c_str(), "grMexCal");
		TGraphErrors *grMexAME = new TGraphErrors();
		grMexAME->SetNameTitle(("grMexAME" + sAna[iAna]).c_str(), "grMexAME");

		TGraphErrors *grResMexRef = new TGraphErrors();
		grResMexRef->SetNameTitle(("grResMexRef" + sAna[iAna]).c_str(), "grResMexRef;Nuclide (=1000#timesZ+A);ME_{cal}-ME_{AME16} [keV]");
		TGraphErrors *grResMexNo = new TGraphErrors();
		grResMexNo->SetNameTitle(("grResMexNo" + sAna[iAna]).c_str(), "grResMexNo;Nuclide (=1000#timesZ+A);ME_{cal}-ME_{AME16} [keV]");

		TGraphErrors *grResMex = new TGraphErrors();
		grResMex->SetNameTitle(("grResMex" + sAna[iAna]).c_str(), "grResMex;Nuclide (=1000#times ZVal+AVal);ME_{cal}-ME_{ref} (keV)");

		printf("\n****************************************************************************************************\n");
		fprintf(fFitMoqTof, "\n****************************************************************************************************\n");
		string strRead, strNucl;
		ifstream fCutInfo("fCutInfo.dat");
		int iNucl = 0, iNuclRef = 0, iNuclNo = 0;
		while (!fCutInfo.eof() && fCutInfo.peek() != EOF)
		{
			if (fCutInfo.peek() != '#')
			{
				int isForCal = 0, ZTar = 0, ATar = 0, QTar = 0;
				double MoQAME[2] = {0};

				fCutInfo >> isForCal >> strRead >> strNucl >> ZTar >> ATar >> QTar >> strRead >> strRead >> strRead >> strRead >> MoQAME[0] >> MoQAME[1];
				getline(fCutInfo, strRead);

				// if (ZTar != QTar)
				// 	continue;
				if (isForCal == -1)
					continue;

				double tofTar[2] = {0};
				int counts = 0;
				for (i = 0; i < nData; i++)
					if (abs(ZVal[i] - ZTar) < 1E-3 && abs(AVal[i] - ATar) < 1E-3 && abs(QVal[i] - QTar) < 1E-3 && abs(tofVal[i]) > 0)
					{
						tofTar[0] = tofVal[i];
						tofTar[1] = tofErrVal[i];
						counts = cntsVal[i];
					}
				if (counts < MINPOIN)
					continue;

				double massTar = 0, mExTar = 0, mExAME[2] = {0};
				double massErr[4] = {0}; //[4]: [0] total error; [1] systematic error; [2] statistical error; [3] fit error
				double massResol = 0;

				massErr[1] = QTar * moqErrSys;
				difFcnMoq->SetParameters(par);
				massErr[2] = QTar * (difFcnMoq->Eval(tofTar[0], ZTar)) * tofTar[1];

				// double parErr2[NPAR][NPAR] = {0};
				// for (i = 0; i < NPAR; i++)
				// 	for (j = 0; j < NPAR; j++)
				// 		parErr2[i][j] = fitMoq->CovMatrix(i, j);
				// double derPar[28] = {0};
				// derPar[0] = 1;
				// derPar[1] = tofTar[0];
				// derPar[2] = ZTar;
				// derPar[3] = tofTar[0] * tofTar[0];
				// derPar[4] = tofTar[0] * ZTar;
				// derPar[5] = ZTar * ZTar;
				// derPar[6] = tofTar[0] * tofTar[0] * tofTar[0];
				// derPar[7] = tofTar[0] * tofTar[0] * ZTar;
				// derPar[8] = tofTar[0] * ZTar * ZTar;
				// derPar[9] = pow(ZTar, 3);
				// derPar[10]=pow(tofTar[0], 4);
				// derPar[11]=pow(tofTar[0], 3)*ZTar;
				// derPar[12]=pow(tofTar[0]*ZTar, 2);
				// derPar[13]=tofTar[0]*pow(ZTar, 3);
				// derPar[14]=pow(ZTar, 4);
				// derPar[15]=pow(tofTar[0], 5);
				// derPar[16]=pow(tofTar[0], 4)*ZTar;
				// derPar[17]=pow(tofTar[0], 3)*pow(ZTar, 2);
				// derPar[18]=pow(tofTar[0],2)*pow(ZTar, 3);
				// derPar[19]=tofTar[0]*pow(ZTar, 4);
				// derPar[20]=pow(ZTar, 5);
				// derPar[21]=pow(tofTar[0], 6);
				// derPar[22]=pow(tofTar[0], 5)*ZTar;
				// derPar[23]=pow(tofTar[0], 4)*pow(ZTar, 2);
				// derPar[24]=pow(tofTar[0],3)*pow(ZTar, 3);
				// derPar[25]=pow(tofTar[0],2)*pow(ZTar, 4);
				// derPar[26]=tofTar[0]*pow(ZTar, 5);
				// derPar[27]=pow(ZTar, 6);
				// double moqErrFit = 0;
				// for (i = 0; i < NPAR; i++)
				// 	for (j = 0; j < NPAR; j++)
				// 		moqErrFit += derPar[i] * derPar[j] * parErr2[i][j];
				// massErr[3] = QTar * sqrt(moqErrFit);

				// double xVar[2] = {tofTar[0], 1.0 * ZTar}, moqErrFit[1] = {0};
				// fitMoq->GetConfidenceIntervals(1, 2, 1, xVar, moqErrFit, 0.682689492137086, false);
				// massErr[3] = QTar * moqErrFit[0];

				massTar = QTar * fcnMoq_Tof_Q->Eval(tofTar[0], ZTar) - (ZTar - QTar) * 511;
				mExTar = massTar - 931494.0954 * ATar + 510.9989461 * ZTar - 0.0144381 * pow(ZTar, 2.39) - 1.55468E-9 * pow(ZTar, 5.35);

				massErr[0] = sqrt(massErr[1] * massErr[1] + massErr[2] * massErr[2] + massErr[3] * massErr[3]);

				massResol = sqrt(massErr[1] * massErr[1] + counts * massErr[2] * massErr[2] + massErr[3] * massErr[3]);

				mExAME[0] = QTar * MoQAME[0] - (ZTar - QTar) * 511 - 931494.0954 * ATar + 510.9989461 * ZTar - 0.0144381 * pow(ZTar, 2.39) - 1.55468E-9 * pow(ZTar, 5.35);
				mExAME[1] = QTar * MoQAME[1];

				grMexCal->SetPoint(iNucl, ZTar * 1000 + ATar, mExTar);
				grMexCal->SetPointError(iNucl, 0, massErr[0]);

				grMexAME->SetPoint(iNucl, ZTar * 1000 + ATar, mExAME[0]);
				grMexAME->SetPointError(iNucl, 0, mExAME[1]);

				if (isForCal == 1)
				{
					grResMexRef->SetPoint(iNuclRef, ZTar * 1000 + ATar, mExTar - mExAME[0]);
					grResMexRef->SetPointError(iNuclRef, 0, sqrt(mExAME[1] * mExAME[1] + massErr[0] * massErr[0]));
					iNuclRef++;
				}
				if (isForCal == 0)
				{
					grResMexNo->SetPoint(iNuclNo, ZTar * 1000 + ATar, mExTar - mExAME[0]);
					grResMexNo->SetPointError(iNuclNo, 0, sqrt(mExAME[1] * mExAME[1] + massErr[0] * massErr[0]));
					iNuclNo++;
				}

				printf("\n%s: isCalib=%d, Z=%d, A=%d, Q=%d,  ToF=%.3f +- %.3f ns ----> m_cal=%.1f +- %.1f +- %.1f +- %.1f keV=%.1f +- %.1f keV, sigma_m=%.1f +- %.1f +- %.1f keV=%.1f keV (%.1E),  ME_cal=%.1f +- %.1f keV, ME_AME=%.1f +- %.1f keV\n", strNucl.c_str(), isForCal, ZTar, ATar, QTar, tofTar[0], tofTar[1], massTar, massErr[1], massErr[2], massErr[3], massTar, massErr[0], massErr[1], sqrt(counts) * massErr[2], massErr[3], massResol, massResol / massTar, mExTar, massErr[0], mExAME[0], mExAME[1]);

				fprintf(fFitMoqTof, "\n%s: isCalib=%d, Z=%d, A=%d, Q=%d,  ToF=%.3f +- %.3f ns ----> m_cal=%.1f +- %.1f +- %.1f +- %.1f keV=%.1f +- %.1f keV, sigma_m=%.1f +- %.1f +- %.1f keV=%.1f keV (%.1E),  ME_cal=%.1f +- %.1f keV, ME_AME=%.1f +- %.1f keV\n", strNucl.c_str(), isForCal, ZTar, ATar, QTar, tofTar[0], tofTar[1], massTar, massErr[1], massErr[2], massErr[3], massTar, massErr[0], massErr[1], sqrt(counts) * massErr[2], massErr[3], massResol, massResol / massTar, mExTar, massErr[0], mExAME[0], mExAME[1]);

				iNucl++;
			}
			else
				getline(fCutInfo, strRead);
		}
		fCutInfo.close();
		printf("****************************************************************************************************\n");
		fprintf(fFitMoqTof, "****************************************************************************************************\n");

		TCanvas *cMex = new TCanvas(("c1Mex" + sAna[iAna]).c_str(), ("c1Mex" + sAna[iAna]).c_str(), 2400, 900);
		cMex->Divide(2, 1);

		cMex->cd(1);
		grMexCal->SetLineColor(kRed);
		grMexCal->SetMarkerColor(kRed);

		grMexAME->SetLineColor(kBlue);
		grMexAME->SetMarkerColor(kBlue);

		TMultiGraph *mgrMex = new TMultiGraph(("mgrMex" + sAna[iAna]).c_str(), "mgrMex;Nuclide (=1000#times Z+A);Mass excess [keV]");
		mgrMex->Add(grMexCal);
		mgrMex->Add(grMexAME);
		mgrMex->Draw("AP*");
		mgrMex->GetYaxis()->SetMaxDigits(3);
		gPad->BuildLegend(0.15, 0.15, 0.35, 0.3)->SetTextSize(0.04);

		cMex->cd(2);
		grResMexNo->SetMarkerColor(kRed);
		grResMexNo->SetLineColor(kRed);
		grResMexNo->SetMarkerSize(2);
		grResMexNo->SetLineWidth(2);
		TMultiGraph *mgrResMex = new TMultiGraph(("mgrResMex" + sAna[iAna]).c_str(), "Mass residuals;Nuclide (=1000#timesZ+A);ME_{fit}-ME_{AME16} [keV]");
		mgrResMex->Add(grResMexRef);
		mgrResMex->Add(grResMexNo);
		mgrResMex->Draw("AP*");
		mgrResMex->GetYaxis()->SetMaxDigits(3);

		cMoq->cd(4);
		mgrResMex->DrawClone("AP*");

		TCanvas *cResM = new TCanvas(("c1ResM" + sAna[iAna] + "").c_str(), ("c1ResM" + sAna[iAna] + "").c_str());
		mgrResMex->DrawClone("AP*");
		// mgrResMex->SetMinimum(-4500);
		// mgrResMex->SetMaximum(6500);

		cMoq->SaveAs(".png");
		cMex->SaveAs(".png");
		cResM->SaveAs(".png");

		fMass->cd();
		cMoq->Write("", TObject::kOverwrite);
		cMex->Write("", TObject::kOverwrite);
		cResM->Write("", TObject::kOverwrite);

		fMoqTof->cd();
		tMoqTof[iAna]->Write();

		delete mgrMex;
		delete cResM;
		delete cMex;
		delete grResMex;
		delete resMoq_Z_AoQ;
		delete resMoq_Z;
		delete resMoq_Tof;
		delete resM_Z_AoQ;
		delete resM_Z;
		delete resM_Tof;
		delete cMoq;
		delete difFcnMoq;
		delete fcnMoq_Tof_Q;
		delete grMoq_Tof_Z;
		delete grAQVal;
		delete grMoqTofZVal;
	}
	delete fTof;
	delete fMoqTof;
	fclose(fFitMoqTof);
	delete fMass;
}

int main(int argc, char **argv)
{
	FitMoqTof();
	return 0;
}