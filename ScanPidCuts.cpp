#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <cmath>

using namespace std;

#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGraph2D.h"

void ScanPidCuts()
{
	gStyle->SetOptStat("nemri");
	gStyle->SetStatFormat("8.6g");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetFitFormat("8.6g");
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);

	const int N = 200, NMIN = 10;
	const int LZ[2] = {-10, 10}, UZ[2] = {-10, 10}, LAoQ[2] = {-10, 10}, UAoQ[2] = {-10, 10};
	const double StZ = 0.02, StAoQ = 0.0002;

	int i = 0, j = 0, k = 0, m = 0;

	int iNucl = 0, nNucl = 0;
	string strNucl[N];
	int isCal[N] = {0}, idN[N] = {0}, ZPid[N] = {0}, APid[N] = {0}, QPid[N] = {0};
	double ZPidCut[N][2] = {0}, AoQPidCut[N][2] = {0}, MoQ[N][2] = {0}, tofPivotAME[N] = {0};
	string strRead;
	TFile *fTof = new TFile("tof-run-150--153__270--385.root");
	TTree *tTOF = fTof->Get<TTree>("tTOF");
	bool good_CfdTdc = 0;
	int Z = 0, A = 0, Q = 0;
	double parTofXs_CfdTdc[6][2] = {0};
	tTOF->SetBranchAddress("good_CfdTdc", &good_CfdTdc);
	tTOF->SetBranchAddress("Z", &Z);
	tTOF->SetBranchAddress("A", &A);
	tTOF->SetBranchAddress("Q", &Q);
	tTOF->SetBranchAddress("parTofXs_CfdTdc", parTofXs_CfdTdc);
	double parTofXs[N][6] = {0};
	bool goodParTofXs[N] = {0};
	nNucl = 0;
	ifstream fCutInfo_Man("fCutInfo_Man.dat");
	while (!fCutInfo_Man.eof() && fCutInfo_Man.peek() != EOF)
	{
		if (fCutInfo_Man.peek() != '#')
		{
			fCutInfo_Man >> isCal[nNucl] >> idN[nNucl] >> strNucl[nNucl] >> ZPid[nNucl] >> APid[nNucl] >> QPid[nNucl] >> ZPidCut[nNucl][0] >> ZPidCut[nNucl][1] >> AoQPidCut[nNucl][0] >> AoQPidCut[nNucl][1] >> MoQ[nNucl][0] >> MoQ[nNucl][1] >> tofPivotAME[nNucl];
			getline(fCutInfo_Man, strRead);
			for (i = 0; i < tTOF->GetEntries(); i++)
			{
				good_CfdTdc = 0;
				Z = 0;
				A = 0;
				Q = 0;
				memset(parTofXs_CfdTdc, 0, sizeof(parTofXs_CfdTdc));
				tTOF->GetEntry(i);
				if (good_CfdTdc && Z == ZPid[nNucl] && A == APid[nNucl] && Q == QPid[nNucl])
				{
					goodParTofXs[nNucl] = true;
					for (j = 0; j < 6; j++)
					{
						parTofXs[nNucl][j] = parTofXs_CfdTdc[j][0];
						goodParTofXs[nNucl] = goodParTofXs[nNucl] && abs(parTofXs_CfdTdc[j][0]) > 0;
					}
					break;
				}
			}
			nNucl++;
		}
		else
			getline(fCutInfo_Man, strRead);
	}
	fCutInfo_Man.close();

	TFile *fCluPid = new TFile("fCluPid.root");

	TTree *tCluPid = fCluPid->Get<TTree>("tCluPid");
	int nData = tCluPid->Draw("AoQPid:ZPid:tof:xMCP:xCrdc[0]:xCrdc[1]", "", "goff");
	double *AoQVal = tCluPid->GetVal(0);
	double *ZVal = tCluPid->GetVal(1);
	double *tofVal = tCluPid->GetVal(2);
	double *xMCP = tCluPid->GetVal(3);
	double *xCrdc1 = tCluPid->GetVal(4);
	double *xCrdc2 = tCluPid->GetVal(5);

	double cenZ = 0, wZ = 0;
	double cenAoQ = 0, wAoQ = 0;
	TFile *fCutScan = new TFile("fCutScan.root", "UPDATE");

	for (iNucl = 0; iNucl < nNucl; iNucl++)
	// for (iNucl = 62; iNucl < 63; iNucl++)
	{
		if (ZPid[iNucl] <= 37 || !goodParTofXs[iNucl])
			continue;

		double AoQPid = 1.0 * APid[iNucl] / QPid[iNucl];
		printf("\n\n********************Now scaning PID cuts of %s: Z=%d, A=%d, Q=%d, A/Q=%f********************\n", strNucl[iNucl].c_str(), ZPid[iNucl], APid[iNucl], QPid[iNucl], AoQPid);
		// printf("\n********************ToF correction parameters********************\n");
		// for(i=0; i<6; i++)
		// 	printf("parTofXs[%d] = %.6G\n", i, parTofXs[iNucl][i]);

		TTree *tCutScan = new TTree(("tCutScan_" + strNucl[iNucl]).c_str(), ("Scan Z-AoQ PID cuts of " + strNucl[iNucl]).c_str());
		double ZCut[2] = {0}, AoQCut[2] = {0}, tofCtr = 0, tofCorrCtr = 0;
		tCutScan->Branch("ZCut", ZCut, "ZCut[2]/D");
		tCutScan->Branch("AoQCut", AoQCut, "AoQCut[2]/D");
		tCutScan->Branch("tofCtr", &tofCtr, "tofCtr/D");
		tCutScan->Branch("tofCorrCtr", &tofCorrCtr, "tofCorrCtr/D");

		double tofPivot = parTofXs[iNucl][0] + parTofXs[iNucl][4] * 110 + parTofXs[iNucl][5] * 113; //xCrdc1's pivot is 110; xCrdc2's pivot is 113

		k = 0;
		for (int iZL = LZ[0]; iZL <= LZ[1]; iZL++)
			for (int iZU = UZ[0]; iZU <= UZ[1]; iZU++)
				for (int iAoQL = LAoQ[0]; iAoQL <= LAoQ[1]; iAoQL++)
					for (int iAoQU = UAoQ[0]; iAoQU <= UAoQ[1]; iAoQU++)
					{
						ZCut[0] = ZPidCut[iNucl][0] + iZL * StZ;
						ZCut[1] = ZPidCut[iNucl][1] + iZU * StZ;
						AoQCut[0] = AoQPidCut[iNucl][0] + iAoQL * StAoQ;
						AoQCut[1] = AoQPidCut[iNucl][1] + iAoQU * StAoQ;
						if (ZCut[1] > ZCut[0] && AoQCut[1] > AoQCut[0])
						{
							cenZ = (ZCut[0] + ZCut[1]) / 2;
							wZ = (ZCut[1] - ZCut[0]) / 2;
							cenAoQ = (AoQCut[0] + AoQCut[1]) / 2;
							wAoQ = (AoQCut[1] - AoQCut[0]) / 2;

							j = 0;
							tofCtr = 0;
							tofCorrCtr = 0;
							for (i = 0; i < nData; i++)
								if (pow((ZVal[i] - cenZ) / wZ, 2) + pow((AoQVal[i] - cenAoQ) / wAoQ, 2) <= 1)
								{
									tofCtr += tofVal[i];
									tofCorrCtr += tofVal[i] + tofPivot - (parTofXs[iNucl][0] + parTofXs[iNucl][1] * xMCP[i] + parTofXs[iNucl][2] * xMCP[i] * xMCP[i] + parTofXs[iNucl][3] * xMCP[i] * xMCP[i] * xMCP[i] + parTofXs[iNucl][4] * xCrdc1[i] + parTofXs[iNucl][5] * xCrdc2[i]);
									j++;
								}
							if (j < NMIN)
								continue;
							tofCtr /= j;
							tofCorrCtr /= j;

							tCutScan->Fill();

							k++;
						}
					}

		fCutScan->cd();
		tCutScan->Write("", TObject::kOverwrite);
		delete tCutScan;
	}
	fCutScan->Close();
	delete fCutScan;
	fCluPid->Close();
	delete fCluPid;
}

int main(int argc, char **argv)
{
	ScanPidCuts();
	return 0;
}