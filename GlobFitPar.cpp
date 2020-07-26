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
#include "TLinearFitter.h"
#include "TFormula.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Fit/FitResult.h"
#include "Math/WrappedMultiTF1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMath.h"

const int Poly = 2, NVar = 3;
const int NPoP = 10;
const int NPar = 6;

double funcParZAQ(double *x, double *par)
{
	// return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[1] + par[4] * x[1] * x[1] + par[5] * x[2];
	// return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[1] + par[4] * x[1] * x[1] + par[5] * x[2] + par[6] * x[2] * x[2];
	// return par[0] + par[1] * x[0] + par[2] * x[1] + par[3] * x[2] + par[4] * x[0] * x[0] + par[5] * x[0] * x[1] + par[6] * x[0] * x[2] + par[7] * x[1] * x[1] + par[8] * x[1] * x[2] + par[9] * x[2] * x[2];

	double f = par[0];
	int i, j, k, m, n, iPar = 0;
	if (Poly >= 1)
		for (i = 0; i < NVar; i++)
		{
			f += par[++iPar] * x[i];
			if (Poly >= 2)
				for (j = i; j < NVar; j++)
				{
					f += par[++iPar] * x[i] * x[j];
					if (Poly >= 3)
						for (k = j; k < NVar; k++)
						{
							f += par[++iPar] * x[i] * x[j] * x[k];
							if (Poly >= 4)
								for (m = k; m < NVar; m++)
								{
									f += par[++iPar] * x[i] * x[j] * x[k] * x[m];
									if (Poly >= 5)
										for (n = m; n < NVar; n++)
											f += par[++iPar] * x[i] * x[j] * x[k] * x[m] * x[n];
								}
						}
				}
		}
	return f;
}

void GlobFitPar()
{
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "MigradImproved");
	// gStyle->SetOptStat("nemri");
	gStyle->SetOptStat(0);
	gStyle->SetStatFormat(".6f");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1111);
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);

	int i, j, k, m;
	int iPar, iData;

	// string strPar[NPar] = {"parXmcp", "parXcrdc1"};
	// string strParUnit[NPar] = {" [ns/mm]", " [ns/pad]"};

	// string strPar[NPar] = {"parXmcp", "parXmcp2", "parXmcp3", "parYmcp", "parXQPlaS800", "parYQPlaS800", "parXcrdc1", "parYcrdc1", "parXcrdc2", "parYcrdc2"};
	// string strParUnit[NPar] = {" [ns/mm]", " [ns/mm^2]", " [ns/mm^3]", " [ns/mm]", " [ns/a.u.]", " [ns/a.u.]", " [ns/pad]", " [ns/a.u.]", " [ns/pad]", " [ns/a.u.]"};
	string strPar[NPar] = {"parXmcp", "parXmcp2", "parXmcp3", "parYmcp", "parXPlaCrdcS800", "parYPlaCrdcS800"};
	string strParUnit[NPar] = {" [ns/mm]", " [ns/mm^2]", " [ns/mm^3]", " [ns/mm]", " [ns/pad]", " [ns/a.u.]"};

	string sZAQ[NVar] = {"AoQ", "Z", "A"};

	const int NAna = 1;
	string sAna[3] = {"CfdTdc", "TacNoClk", "TacClk"};
	int iAna = 0;

	TFile *fGlobFitPar = new TFile("fGlobFitPar.root", "RECREATE");
	TTree *tGlobFitPar[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
		tGlobFitPar[iAna] = new TTree(("tGlobFitPar_" + sAna[iAna]).c_str(), ("tree of " + sAna[iAna] + " global fitting parameters").c_str());

	TFile *fPar = new TFile("../MachineLearning/fPar.root", "RECREATE");
	TTree *tPar[NAna];
	for (iAna = 0; iAna < NAna; iAna++)
		tPar[iAna] = new TTree(("tPar_" + sAna[iAna]).c_str(), ("tree of " + sAna[iAna] + " to fit global parameters").c_str());

	TFile *fTof = new TFile("/home/kailong/ExpData/Jul2018/Pid/tof-run-150--153__270--385.root");
	TTree *tTOF;
	fTof->GetObject("tTOF", tTOF);
	tTOF->SetEstimate(-1);

	for (iAna = 0; iAna < NAna; iAna++)
	{
		string sDraw = "Z:A:Q:AoQ";
		for (iPar = 0; iPar < NPar; iPar++)
			for (j = 0; j < 2; j++)
				sDraw += ":parTofXs_" + sAna[iAna] + "[" + to_string(iPar + 1) + "][" + to_string(j) + "]";
		string sCut = "good_" + sAna[iAna] + " && abs(tofXsCorr_" + sAna[iAna] + "[0])>0 && Z>37";
		int nData = tTOF->Draw(sDraw.c_str(), sCut.c_str(), "goff");
		double *ZVal = tTOF->GetVal(0);
		double *AVal = tTOF->GetVal(1);
		double *QVal = tTOF->GetVal(2);
		double *AoQVal = tTOF->GetVal(3);
		double *par[NPar], *parErr[NPar];
		for (iPar = 0; iPar < NPar; iPar++)
		{
			par[iPar] = tTOF->GetVal(4 + 2 * iPar);
			parErr[iPar] = tTOF->GetVal(5 + 2 * iPar);
		}
		double *xVal = new double[nData * NVar];
		double **errFit = new double *[NPar];
		for (iPar = 0; iPar < NPar; iPar++)
			errFit[iPar] = new double[nData];

		double Z = 0, A = 0, Q = 0, AoQ = 0;
		double globPar[NPar][2] = {0}, parOfPar[NPar][NPoP] = {0}, parErrOfPar[NPar][NPoP] = {0};
		tGlobFitPar[iAna]->Branch("Z", &Z, "Z/D");
		tGlobFitPar[iAna]->Branch("A", &A, "A/D");
		tGlobFitPar[iAna]->Branch("Q", &Q, "Q/D");
		tGlobFitPar[iAna]->Branch("AoQ", &AoQ, "AoQ/D");
		for (iPar = 0; iPar < NPar; iPar++)
		{
			tGlobFitPar[iAna]->Branch(("glob_" + strPar[iPar]).c_str(), globPar[iPar], ("glob_" + strPar[iPar] + "[2]/D").c_str());
			tGlobFitPar[iAna]->Branch(("parOf_" + strPar[iPar]).c_str(), parOfPar[iPar], ("parOf_" + strPar[iPar] + "[" + to_string(NPoP) + "]/D").c_str());
			tGlobFitPar[iAna]->Branch(("parErrOf_" + strPar[iPar]).c_str(), parErrOfPar[iPar], ("parErrOf_" + strPar[iPar] + "[" + to_string(NPoP) + "]/D").c_str());
		}
		double myPar[NPar][2] = {0};
		tPar[iAna]->Branch("myPar", myPar, ("myPar[" + to_string(NPar) + "][2]/D").c_str());
		tPar[iAna]->Branch("Z", &Z, "Z/D");
		tPar[iAna]->Branch("A", &A, "A/D");
		tPar[iAna]->Branch("Q", &Q, "Q/D");
		tPar[iAna]->Branch("AoQ", &AoQ, "AoQ/D");

		for (iData = 0; iData < nData; iData++)
		{
			Z = ZVal[iData];
			A = AVal[iData];
			Q = QVal[iData];
			AoQ = AoQVal[iData];
			for (iPar = 0; iPar < NPar; iPar++)
			{
				myPar[iPar][0] = par[iPar][iData];
				myPar[iPar][1] = parErr[iPar][iData];
			}
			tPar[iAna]->Fill();
		}
		fPar->cd();
		tPar[iAna]->Write();

		TF1 fcnParZAQ("fcnParZAQ", funcParZAQ, 0, 1, NPoP); //define a TF1 function using "general C function with parameters"; xmin, xmax (their values don't matter) and parameter number must be given
		double xPt[NVar] = {0};
		for (iPar = 0; iPar < NPar; iPar++)
		{
			cout << "\n\n********************Global fitting parameter of " << strPar[iPar] << "********************\n";
			for (i = 0; i < NVar * nData; i++)
				xVal[i] = 0;
			for (i = 0; i < nData; i++)
				errFit[iPar][i] = 0;

			// ROOT::Fit::BinData dataParAoQ(nData, NVar); // set the data set to be fitted; nData: number of data points (necessary and >= atcual data points number); NVar: dimension (necessary == actual dimensional number); also ErrorType can be given as ROOT::Fit::BinData::{kNoError, kValueError, kCoordError, kAsymError} "kValueError" is the default
			ROOT::Fit::BinData dataParAoQ(nData, NVar, ROOT::Fit::BinData::kNoError);
			for (iData = 0; iData < nData; iData++)
			{
				xPt[0] = AoQVal[iData];
				xPt[1] = ZVal[iData];
				xPt[2] = AVal[iData];
				// dataParAoQ.Add(xPt, par[iPar][iData], parErr[iPar][iData]);
				dataParAoQ.Add(xPt, par[iPar][iData]);
			}
			ROOT::Math::WrappedMultiTF1 fitFunc(fcnParZAQ, NVar); //generate the general parametric function by wrapping the TF1 function; NVarmeaning dimension should be given
			ROOT::Fit::Fitter fitter;							  //define fitter object
			fitter.SetFunction(fitFunc, true);					  //set fit function; false (usually) means NOT use gradient of fucntion
			for (i = 0; i < NPoP; i++)
				fitter.Config().ParSettings(i).SetValue(1); //set initial parameter values; Not very necessary in linear case, but better to set
			// fitter.LeastSquareFit(dataParAoQ);				  //perform fit
			fitter.LinearFit(dataParAoQ);					  //perform fit using linear
			ROOT::Fit::FitResult resultFit = fitter.Result(); //get fit results
			resultFit.Print(std::cout);						  //print fit results

			for (i = 0; i < NPoP; i++)
			{
				parOfPar[iPar][i] = resultFit.Parameter(i);
				parErrOfPar[iPar][i] = resultFit.ParError(i);
			}
			TGraphErrors *grParZAQ[NVar], *grResParZAQ[NVar];
			for (k = 0; k < NVar; k++)
			{
				grParZAQ[k] = new TGraphErrors(nData);
				grParZAQ[k]->SetNameTitle(("gr_" + strPar[iPar] + "_" + sZAQ[k] + "_" + sAna[iAna]).c_str(), (strPar[iPar] + " VS " + sZAQ[k] + ";" + sZAQ[k] + ";" + strPar[iPar] + strParUnit[iPar]).c_str());

				grResParZAQ[k] = new TGraphErrors(nData);
				grResParZAQ[k]->SetNameTitle(("gr_Res_" + strPar[iPar] + "_" + sZAQ[k] + "_" + sAna[iAna]).c_str(), (strPar[iPar] + " residuals VS " + sZAQ[k] + ";" + sZAQ[k] + ";#Delta" + strPar[iPar] + strParUnit[iPar]).c_str());
			}
			for (iData = 0; iData < nData; iData++)
			{
				xPt[0] = AoQVal[iData];
				xPt[1] = ZVal[iData];
				xPt[2] = AVal[iData];

				for (k = 0; k < NVar; k++)
				{
					grParZAQ[k]->SetPoint(iData, xPt[k], par[iPar][iData]);
					// grParZAQ[k]->SetPointError(iData, 0, parErr[iPar][iData]);
					grParZAQ[k]->SetPointError(iData, 0, 0);
					xVal[NVar * iData + k] = xPt[k];
				}
				for (k = 0; k < NVar; k++)
					grResParZAQ[k]->SetPoint(iData, xPt[k], par[iPar][iData] - fcnParZAQ(xPt, parOfPar[iPar]));
			}
			// resultFit.GetConfidenceIntervals(nData, NVar, 1, xVal, errFit[iPar], 0.682689492137086, true);
			for (iData = 0; iData < nData; iData++)
				for (k = 0; k < NVar; k++)
					grResParZAQ[k]->SetPointError(iData, 0, errFit[iPar][iData]);

			TCanvas *cPar = new TCanvas(("c" + strPar[iPar] + "_" + sAna[iAna]).c_str(), ("c" + strPar[iPar] + "_" + sAna[iAna]).c_str(), 600 * NVar, 900);
			cPar->Divide(NVar, 2);
			for (k = 0; k < NVar; k++)
			{
				cPar->cd(k + 1);
				grParZAQ[k]->Draw("AP*");
				// grParZAQ[k]->GetYaxis()->SetRangeUser(grParZAQ[k]->GetMean(2) - 6 * grParZAQ[k]->GetRMS(2), grParZAQ[k]->GetMean(2) + 6 * grParZAQ[k]->GetRMS(2));
				// grParZAQ[k]->SetMarkerSize(2);
				grParZAQ[k]->GetYaxis()->SetMaxDigits(3);

				cPar->cd(k + NVar + 1);
				grResParZAQ[k]->Draw("AP*");
				// grResParZAQ[k]->SetMarkerSize(2);
				// grResParZAQ[k]->GetYaxis()->SetRangeUser(grResParZAQ[k]->GetMean(2)-5*grResParZAQ[k]->GetRMS(2), grResParZAQ[k]->GetMean(2)+5*grResParZAQ[k]->GetRMS(2));
				grResParZAQ[k]->GetYaxis()->SetMaxDigits(3);
			}
			printf("Mean Residual: %f, Sigma Residual: %f\n", grResParZAQ[0]->GetMean(2), grResParZAQ[0]->GetRMS(2));
			fGlobFitPar->cd();
			cPar->Write("", TObject::kOverwrite);
			cPar->SaveAs(("./Graphs/CfdTdc/c" + strPar[iPar] + "_" + sAna[iAna] + ".png").c_str());

			for (k = 0; k < NVar; k++)
			{
				delete grResParZAQ[k];
				delete grParZAQ[k];
			}
			delete cPar;
		}
		for (iData = 0; iData < nData; iData++)
		{
			Z = ZVal[iData];
			A = AVal[iData];
			Q = QVal[iData];
			AoQ = AoQVal[iData];
			for (iPar = 0; iPar < NPar; iPar++)
			{
				xPt[0] = AoQVal[iData];
				xPt[1] = ZVal[iData];
				xPt[2] = AVal[iData];
				globPar[iPar][0] = fcnParZAQ(xPt, parOfPar[iPar]);
				globPar[iPar][1] = errFit[iPar][iData];
			}
			tGlobFitPar[iAna]->Fill();
		}
		for (iPar = 0; iPar < NPar; iPar++)
			delete[] errFit[iPar];
		delete[] errFit;
		delete[] xVal;
	}
	fGlobFitPar->cd();
	for (iAna = 0; iAna < NAna; iAna++)
		tGlobFitPar[iAna]->Write();

	delete fTof;
	delete fPar;
	delete fGlobFitPar;
}

int main(int argc, char **argv)
{
	GlobFitPar();
	return 0;
}