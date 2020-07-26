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
#include "TTree.h"
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF2.h"
#include "TFormula.h"
#include "TLinearFitter.h"
#include "TAxis.h"

void ScanXpivot()
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

	const int NXs = 3, NParTof = 6, NParMoq = 10;
	const int N = 200;
	const int NX[2] = {-50, 50};
	const double StpX[NXs] = {0.1, 0.1, 0.1};
	const double XINI[NXs]={0, 110, 113};

	int i = 0, j = 0, k = 0;

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

	string strDir = "/home/kailong/ExpData/Jul2018/Pid/Graphs/ResMXpivot/";
	if (access(strDir.c_str(), F_OK) != 0)
		system(("mkdir " + strDir).c_str());
	string strGif = strDir + "cResMXpivot.gif";
	if (access(strGif.c_str(), F_OK) == 0)
		system(("rm -f " + strGif).c_str());

	TFile *fTof = new TFile("tof-run-150--153__270--385.root");
	TTree *tTOF = fTof->Get<TTree>("tTOF");
	tTOF->SetEstimate(-1);

	string sDraw = "isCal:Z:A:Q:AoQ:MoQ[0]:MoQ[1]:tofXsCorr_CfdTdc[1]";
	for (i = 0; i < NParTof; i++)
		sDraw += ":parTofXs_CfdTdc[" + to_string(i) + "][0]";
	string sCut = "good_CfdTdc && abs(tofXsCorr_CfdTdc[0])>0 && Z>37";

	int nData = tTOF->Draw(sDraw.c_str(), sCut.c_str(), "goff");
	double *calVal = tTOF->GetVal(0);
	double *ZVal = tTOF->GetVal(1);
	double *AVal = tTOF->GetVal(2);
	double *QVal = tTOF->GetVal(3);
	double *AoQVal = tTOF->GetVal(4);
	double *MoQVal = tTOF->GetVal(5);
	double *MoQErrVal = tTOF->GetVal(6);
	double *tofErrVal = tTOF->GetVal(7);
	double *parVal[NParTof];
	for (i = 0; i < NParTof; i++)
		parVal[i] = tTOF->GetVal(8 + i);
	
	FILE *fResMScan=fopen("fResMScan.dat", "w+");

	// TFile *fXpivotScan = new TFile("fXpivotScan.root", "RECREATE");
	TFormula *fcnTofXs = new TFormula("fcnTofXs", "[0]+[1]*x[0]+[2]*x[0]*x[0]+[3]*x[0]*x[0]*x[0]+[4]*x[1]+[5]*x[2]");
	TF2 *fcnMoq_Tof_ZQ = new TF2("fcnMoq_Tof_ZQ", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y+[6]*x*x*x+[7]*x*x*y+[8]*x*y*y+[9]*y*y*y");
	TLinearFitter *fitMoq = new TLinearFitter(2, "1++x++y++x*x++x*y++y*y++x*x*x++x*x*y++x*y*y++y*y*y");
	double *tofPivot = new double[nData];
	double xTzq[2] = {0};
	int iter = 0;
	int iX[NXs]={0};
	for (iX[0] = NX[0]; iX[0] <= NX[1]; iX[0]++)
		for (iX[1] = NX[0]; iX[1] <= NX[1]; iX[1]++)
			for (iX[2] = NX[0]; iX[2] <= NX[1]; iX[2]++)
			{
				iter++;
				double xPivot[NXs] = {0};
				for(k=0; k<NXs; k++)
					xPivot[k] = XINI[k] + iX[k] * StpX[k];
				fitMoq->ClearPoints();
				int iData = 0;
				for (iData = 0; iData < nData; iData++)
				{
					tofPivot[iData] = 0;
					double parTof[NParTof] = {0};
					for (j = 0; j < NParTof; j++)
						parTof[j] = parVal[j][iData];
					tofPivot[iData] = fcnTofXs->EvalPar(xPivot, parTof);
					if (abs(calVal[iData] - 1) < 1E-3)
					{
						xTzq[0] = tofPivot[iData];
						xTzq[1] = QVal[iData];
						fitMoq->AddPoint(xTzq, MoQVal[iData], MoQErrVal[iData]);
					}
				}
				fitMoq->Eval();
				// fitMoq->PrintResults(3);
				double parMoq[NParMoq] = {0};
				for (i = 0; i < NParMoq; i++)
					parMoq[i] = fitMoq->GetParameter(i);

				TGraph *grResMRef = new TGraph();
				grResMRef->SetName("grResMRef");
				TGraph *grResMNo = new TGraph();
				grResMNo->SetName("grResMNo");
				int iPtRef = 0, iPtNo = 0;
				iPtRef = 0;
				for (iData = 0; iData < nData; iData++)
				{
					xTzq[0] = tofPivot[iData];
					xTzq[1] = QVal[iData];
					double mass = QVal[iData] * fcnMoq_Tof_ZQ->EvalPar(xTzq, parMoq) - (ZVal[iData] - QVal[iData]) * 511;
					double massAME = 0;
					for (iNucl = 0; iNucl < nNucl; iNucl++)
						if (abs(ZVal[iData] - ZOrig[iNucl]) < 1E-3 && abs(AVal[iData] - AOrig[iNucl]) < 1E-3 && abs(QVal[iData] - QOrig[iNucl]) < 1E-3)
							massAME = QOrig[iNucl] * MoQOrig[iNucl][0] - (ZOrig[iNucl] - QOrig[iNucl]) * 511;
					if (abs(calVal[iData] - 1) < 1E-3)
						grResMRef->SetPoint(iPtRef++, 1000 * ZVal[iData] + AVal[iData], mass - massAME);
					if (abs(calVal[iData]) < 1E-3)
						grResMNo->SetPoint(iPtNo++, 1000 * ZVal[iData] + AVal[iData], mass - massAME);
				}
				if (grResMRef->GetRMS(2) > 1250)
				{
					delete grResMNo;
					delete grResMRef;
					continue;
				}
				string sCan = "cResM_Xpivot_" + to_string(xPivot[0]) + "_" + to_string(xPivot[1]) + "_" + to_string(xPivot[2]);
				TCanvas *cResM = new TCanvas(sCan.c_str(), sCan.c_str());
				TMultiGraph *mgrResM = new TMultiGraph("mgrResM", ("Mass residuals with Xpivots of " + to_string(xPivot[0]) + "__" + to_string(xPivot[1]) + "__" + to_string(xPivot[2]) + ";Nuclide (=1000#timesZ+A);M_{fit}-M_{AME16} [keV]").c_str());
				grResMRef->SetMarkerColor(kBlack);
				grResMNo->SetMarkerColor(kRed);
				mgrResM->Add(grResMRef);
				mgrResM->Add(grResMNo);
				mgrResM->Draw("AP*");
				mgrResM->GetYaxis()->SetMaxDigits(3);

				printf("\nIter#%d:  ", iter);
				for(j=0; j<NXs; j++)
					printf("xPivot_%d = %.2f  ", j, xPivot[j]);
				printf("Average mass residual = %.0f +- %.0f keV\n", grResMRef->GetMean(2), grResMRef->GetRMS(2));
				
				fprintf(fResMScan, "\nIter#%d:  ", iter);
				for(j=0; j<NXs; j++)
					fprintf(fResMScan, "xPivot_%d = %.2f  ", j, xPivot[j]);
				fprintf(fResMScan, "Average mass residual = %.0f +- %.0f keV\n", grResMRef->GetMean(2), grResMRef->GetRMS(2));

				// cResM->Write("", TObject::kOverwrite);
				cResM->SaveAs((strDir + sCan + ".png").c_str());
				cResM->SaveAs((strGif + "+50").c_str());

				delete grResMNo;
				delete grResMRef;
				delete mgrResM;
				delete cResM;
			}
	delete[] tofPivot;
	delete fitMoq;
	delete fcnMoq_Tof_ZQ;
	delete fcnTofXs;
	// fXpivotScan->Close(0);
	// delete fXpivotScan;
	fTof->Close();
	delete fTof;
}

int main()
{
	ScanXpivot();
	return 0;
}