#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <cmath>

using namespace std;

#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"

void AdjPidCut()
{
	gStyle->SetOptStat("nemri");
	gStyle->SetStatFormat("8.6g");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetFitFormat("8.6g");
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);

	const int N = 200;
	const double ZRang = 0.04, AoQRang = 0.0004;

	int i = 0, j = 0, k = 0;

	string strNucl[N];
	int isCal[N] = {0}, idN[N] = {0}, ZPid[N] = {0}, APid[N] = {0}, QPid[N] = {0};
	double ZPidCut[N][2] = {0}, AoQPidCut[N][2] = {0}, MoQ[N][2] = {0}, tofPivotAME[N] = {0};
	string strRead;
	int iNucl = 0, nNucl = 0;
	ifstream fCutInfo_Man("fCutInfo_Man.dat");
	while (!fCutInfo_Man.eof() && fCutInfo_Man.peek() != EOF)
	{
		if (fCutInfo_Man.peek() != '#')
		{
			fCutInfo_Man >> isCal[nNucl] >> idN[nNucl] >> strNucl[nNucl] >> ZPid[nNucl] >> APid[nNucl] >> QPid[nNucl] >> ZPidCut[nNucl][0] >> ZPidCut[nNucl][1] >> AoQPidCut[nNucl][0] >> AoQPidCut[nNucl][1] >> MoQ[nNucl][0] >> MoQ[nNucl][1] >> tofPivotAME[nNucl];
			getline(fCutInfo_Man, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo_Man, strRead);
	}
	fCutInfo_Man.close();

	FILE *fCutInfo = fopen("fCutInfo.dat", "w+");

	TFile *fCutScan = new TFile("fCutScan.root");
	TTree *tMaxTof = fCutScan->Get<TTree>("tMaxTof");
	int Z = 0, A = 0, Q = 0;
	double tofMax = 0, tofCorrMax = 0;
	tMaxTof->SetBranchAddress("Z", &Z);
	tMaxTof->SetBranchAddress("A", &A);
	tMaxTof->SetBranchAddress("Q", &Q);
	tMaxTof->SetBranchAddress("tofMax", &tofMax);
	tMaxTof->SetBranchAddress("tofCorrMax", &tofCorrMax);

	for (iNucl = 0; iNucl < nNucl; iNucl++)
	// for (iNucl = 62; iNucl < 64; iNucl++)
	{
		TTree *tCutScan = fCutScan->Get<TTree>(("tCutScan_" + strNucl[iNucl]).c_str());
		if(!(tCutScan->GetBranch("tofCorrCtr")))
			continue;
		double ZCut[2] = {0}, AoQCut[2] = {0}, tofCtr = 0, tofCorrCtr = 0;
		tCutScan->SetBranchAddress("ZCut", ZCut);
		tCutScan->SetBranchAddress("AoQCut", AoQCut);
		tCutScan->SetBranchAddress("tofCtr", &tofCtr);
		tCutScan->SetBranchAddress("tofCorrCtr", &tofCorrCtr);
		int nData = tCutScan->GetEntries();

		double maxTof = 0, maxTofCorr=0;
		for (i = 0; i < tMaxTof->GetEntries(); i++)
		{
			tMaxTof->GetEntry(i);
			if (Z == ZPid[iNucl] && A == APid[iNucl] && Q == QPid[iNucl])
			{
				maxTof = tofMax;
				maxTofCorr = tofCorrMax;
				break;
			}
		}
		bool isAdj = false;
		for (int iData = 0; iData < nData; iData++)
		{
			tCutScan->GetEntry(iData);
			// if (abs(tofCtr - maxTof) < 0.001 && abs(ZCut[0] - ZPidCut[iNucl][0]) < ZRang && abs(ZCut[1] - ZPidCut[iNucl][1]) < ZRang && abs(AoQCut[0] - AoQPidCut[iNucl][0]) < AoQRang && abs(AoQCut[1] - AoQPidCut[iNucl][1]) < AoQRang)
			if (abs(tofCorrCtr - maxTofCorr) < 0.001 && abs(ZCut[0] - ZPidCut[iNucl][0]) < ZRang && abs(ZCut[1] - ZPidCut[iNucl][1]) < ZRang && abs(AoQCut[0] - AoQPidCut[iNucl][0]) < AoQRang && abs(AoQCut[1] - AoQPidCut[iNucl][1]) < AoQRang)
			{
				fprintf(fCutInfo, "%d  %d  %s  %d  %d  %d  %.2f %.2f  %.4f %.4f  %.3f %.3f  %.6f\n", isCal[iNucl], idN[iNucl], strNucl[iNucl].c_str(), ZPid[iNucl], APid[iNucl], QPid[iNucl], ZCut[0], ZCut[1], AoQCut[0], AoQCut[1], MoQ[iNucl][0], MoQ[iNucl][1], tofPivotAME[iNucl]);
				isAdj = true;
				break;
			}
		}
		if (!isAdj)
		{
			printf("Initial PID cut for %s is NOT chanegd!\n", strNucl[iNucl].c_str());
			fprintf(fCutInfo, "%d  %d  %s  %d  %d  %d  %.2f %.2f  %.4f %.4f  %.3f %.3f  %.6f\n", isCal[iNucl], idN[iNucl], strNucl[iNucl].c_str(), ZPid[iNucl], APid[iNucl], QPid[iNucl], ZPidCut[iNucl][0], ZPidCut[iNucl][1], AoQPidCut[iNucl][0], AoQPidCut[iNucl][1], MoQ[iNucl][0], MoQ[iNucl][1], tofPivotAME[iNucl]);
		}
	}
	fclose(fCutInfo);
	fCutScan->Close();
	delete fCutScan;
}

int main(int argc, char **argv)
{
	AdjPidCut();
	return 0;
}