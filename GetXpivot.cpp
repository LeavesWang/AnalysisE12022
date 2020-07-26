#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include <cstdlib>

using namespace std;

#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTree.h"
#include "TLinearFitter.h"

#include "TF1.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Fit/FitResult.h"
#include "Math/WrappedMultiTF1.h"

#include "TMath.h"

double myFunc(double *par, double *x)
{
	// return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[1] + par[5] * x[2] + par[6] * x[3] + par[7] * x[4] + par[8] * x[5];
	return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[1] + par[5] * x[2];
}

void GetXpivot()
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

	const int NXs = 3, NPar = 6;
	const int NAna = 1;
	string sAna[3] = {"CfdTdc", "TacNoClk", "TacClk"};
	int iAna = 0, i, j, k;

	TFile *fTof = new TFile("tof-run-150--153__270--385.root");
	TTree *tTOF, *tTOFDetail;
	fTof->GetObject("tTOF", tTOF);
	fTof->GetObject("tTOFDetail", tTOFDetail);
	tTOF->SetEstimate(-1);
	tTOFDetail->SetEstimate(-1);

	for (iAna = 0; iAna < NAna; iAna++)
	{
		printf("\n********************Getting pivot Xs of %s electronics setup********************\n", sAna[iAna].c_str());

		string sDraw = "xMcp:xCrdc[0]:xCrdc[1]";
		string sCut = "good_" + sAna[iAna] + " && Z>37";

		int nData = tTOFDetail->Draw(sDraw.c_str(), sCut.c_str(), "goff");

		double *xVal[NXs];
		for (i = 0; i < NXs; i++)
			xVal[i] = tTOFDetail->GetVal(i);

		double xPivot[NXs] = {0};
		for (j = 0; j < NXs; j++)
		{
			xPivot[j] = TMath::Mean(nData, xVal[j]);
			printf("xPivot[%d] = %f;\n", j, xPivot[j]);
		}
	}

	// for (iAna = 0; iAna < NAna; iAna++)
	// {
	// 	printf("\n********************Getting pivot Xs of %s electronics setup********************\n", sAna[iAna].c_str());

	// 	string sDraw = "tofAME";
	// 	for(i=0; i<NPar; i++)
	// 		sDraw+=":parTofXs_"+sAna[iAna]+"["+to_string(i)+"][0]";
	// 	string sCut = "good_" + sAna[iAna] + " && abs(tofXsCorr_" + sAna[iAna] + "[0])>0 && Z>37";

	// 	int nData = tTOF->Draw(sDraw.c_str(), sCut.c_str(), "goff");

	// 	double *tofAME=tTOF->GetVal(0);
	// 	double *parVal[NPar];
	// 	for(i=1; i<=NPar; i++)
	// 		parVal[i-1]=tTOF->GetVal(i);

	// 	ROOT::Fit::BinData dataPars(nData, NPar, ROOT::Fit::BinData::kNoError);
	// 	int iData=0;
	// 	double parPt[NPar]={0};
	// 	for (iData = 0; iData < nData; iData++)
	// 	{
	// 		memset(parPt, 0, sizeof(parPt));
	// 		for (k = 0; k < NPar; k++)
	// 			parPt[k] = parVal[k][iData];
	// 		dataPars.Add(parPt, tofAME[iData]);
	// 	}
	// 	TF1 fcnToFamePar("fcnToFamePar", myFunc, 0, 1, NXs);
	// 	ROOT::Math::WrappedMultiTF1 fitFunc(fcnToFamePar, NPar);
	// 	ROOT::Fit::Fitter fitter;
	// 	fitter.SetFunction(fitFunc, true);
	// 	fitter.Config().ParSettings(0).SetLimits(-10,10);
	// 	fitter.Config().ParSettings(1).SetLimits(105,115);
	// 	fitter.Config().ParSettings(2).SetLimits(109,119);

	// 	fitter.LeastSquareFit(dataPars);
	// 	ROOT::Fit::FitResult resultFit = fitter.Result();
	// 	resultFit.Print(std::cout);
	// 	double xPivot[NXs]={0};
	// 	for(j=0; j<NXs; j++)
	// 	{
	// 		xPivot[j]=resultFit.Parameter(j);
	// 		printf("xPivot[%d] = %f;\n", j, xPivot[j]);
	// 	}

	// 	TGraph *grResTofAME = new TGraph(nData);
	// 	grResTofAME->SetNameTitle("grResTofAME", "Residuals of fitting tofAME with correction parameters;tofAME [ns];tofAME_{fit}-tofAME [ps]");
	// 	for (iData = 0; iData < nData; iData++)
	// 	{
	// 		memset(parPt, 0, sizeof(parPt));
	// 		for (k = 0; k < NPar; k++)
	// 			parPt[k] = parVal[k][iData];
	// 		grResTofAME->SetPoint(iData, tofAME[iData], 1000*(fcnToFamePar(parPt, xPivot)-tofAME[iData]));
	// 	}
	// 	TCanvas *cResTofAME = new TCanvas("cResTofAME", "cResTofAME");
	// 	grResTofAME->Draw("AP*");
	// 	cResTofAME->SaveAs("cResTofAME.png");
	// 	delete cResTofAME;
	// 	delete grResTofAME;
	// }
	fTof->Close();
	delete fTof;
}

int main()
{
	GetXpivot();
	return 0;
}