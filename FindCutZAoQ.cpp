#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>

using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TCutG.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TPad.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TGraph2D.h"

void FindCutZAoQ()
{	
	gStyle->SetOptStat("nemri");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetFitFormat(".6f");
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);
	
	const int N=200, NMIN=100;
	
	int i=0, j=0, k=0;
	
	int iNucl=0, nNucl=0;
	string strNucl[N];
	int isCal[N]={0}, idN[N]={0}, ZPid[N]={0}, APid[N]={0}, QPid[N]={0};
	double ZPidCut[N][2]={0}, AoQPidCut[N][2]={0}, MoQ[N][2]={0}, tofPivot[N]={0};
	string strRead;
	nNucl=0;
	ifstream fCutInfo_Man("fCutInfo_Man.dat");
	while(!fCutInfo_Man.eof()&&fCutInfo_Man.peek()!=EOF)
	{
		if(fCutInfo_Man.peek()!='#')
		{
			fCutInfo_Man>>isCal[nNucl]>>idN[nNucl]>>strNucl[nNucl]>>ZPid[nNucl]>>APid[nNucl]>>QPid[nNucl]>>ZPidCut[nNucl][0]>>ZPidCut[nNucl][1]>>AoQPidCut[nNucl][0]>>AoQPidCut[nNucl][1]>>MoQ[nNucl][0]>>MoQ[nNucl][1]>>tofPivot[nNucl];
			getline(fCutInfo_Man, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo_Man, strRead);
	}
	fCutInfo_Man.close();
	
	TFile *fCut=new TFile("fCut.root");
	TGraph2D *grTofZAoQPid;
	fCut->GetObject("grTofZAoQPid", grTofZAoQPid);
	double *Z=grTofZAoQPid->GetY();
	double *AoQ=grTofZAoQPid->GetX();
	int nData=grTofZAoQPid->GetN();
	
	FILE *fCutInfo=fopen("fCutInfo.dat", "w+");
	for(iNucl=0; iNucl<nNucl; iNucl++)
	{
		printf("\n********************Now drawing PID information of: Z=%d, A=%d, Q=%d, A/Q=%f\n", ZPid[iNucl], APid[iNucl], QPid[iNucl], 1.0*APid[iNucl]/QPid[iNucl]);
		
		double ZCut[2]={ZPidCut[iNucl][0], ZPidCut[iNucl][1]};
		double cenZ=(ZCut[0]+ZCut[1])/2;
		double wZ=(ZCut[1]-ZCut[0])/2;

		double AoQCut[2]={AoQPidCut[iNucl][0], AoQPidCut[iNucl][1]};
		double cenAoQ=(AoQCut[0]+AoQCut[1])/2;
		double wAoQ=(AoQCut[1]-AoQCut[0])/2;

		
		printf("Initially: Z cut: %.2f-->%.2f,  A/Q cut: %.4f-->%.4f\n", ZCut[0], ZCut[1], AoQCut[0], AoQCut[1]);
		// printf("Initially: Z cut: %.2f +- %.2f;  A/Q cut: %.4f +- %.4f\n", cenZ, wZ, cenAoQ, wAoQ);
		
		TGraph *grZAoQ=new TGraph();
		k=0;
		for(i=0; i<nData; i++)
			if(pow((Z[i]-cenZ)/wZ,2)+pow((AoQ[i]-cenAoQ)/wAoQ,2)<=1)
			// if(Z[i]>ZCut[0]&&Z[i]<ZCut[1] && AoQ[i]>AoQCut[0]&&AoQ[i]<AoQCut[1])
			{
				grZAoQ->SetPoint(k++, AoQ[i], Z[i]);
			}
		if(k<NMIN)
		{
			delete grZAoQ;
			continue;
		}
		// int nPoin;
		// double ZCutTem[2]={0}, AoQCutTem[2]={0};
		// double wZTem=0, wAoQTem=0;
		// double cenZTem=0, cenAoQTem=0;
		// for(j=0; j<100; j++)
		// {
		// 	cenAoQTem=grZAoQ->GetMean(1);
		// 	wAoQTem=NSgmAoQ*grZAoQ->GetRMS(1);
		// 	AoQCutTem[0]=cenAoQTem-wAoQTem;
		// 	AoQCutTem[1]=cenAoQTem+wAoQTem;
			
		// 	cenZTem=grZAoQ->GetMean(2);
		// 	wZTem=NSgmZ*grZAoQ->GetRMS(2);
		// 	ZCutTem[0]=cenZTem-wZTem;
		// 	ZCutTem[1]=cenZTem+wZTem;
			
		// 	if(abs(ZCutTem[0]-ZCut[0])<0.01&&abs(ZCutTem[1]-ZCut[1])<0.01&&abs(AoQCutTem[0]-AoQCut[0])<0.0001&&abs(AoQCutTem[1]-AoQCut[1])<0.0001 && abs(wZTem-wZ)<0.01&&abs(wAoQTem-wAoQ)<0.0001&&abs(cenZTem-cenZ)<0.01&&abs(cenAoQTem-cenAoQ)<0.0001)
		// 	// if(abs(AoQCutTem[0]-AoQCut[0])<0.0001&&abs(AoQCutTem[1]-AoQCut[1])<0.0001 &&abs(wAoQTem-wAoQ)<0.0001&&abs(cenAoQTem-cenAoQ)<0.0001)
		// 		break;
		// 	ZCut[0]=ZCutTem[0];
		// 	ZCut[1]=ZCutTem[1];
		// 	wZ=wZTem;
		// 	cenZ=cenZTem;
		// 	AoQCut[0]=AoQCutTem[0];
		// 	AoQCut[1]=AoQCutTem[1];
		// 	wAoQ=wAoQTem;
		// 	cenAoQ=cenAoQTem;
		// printf("Iteration#%d: Z cut: %.2f-->%.2f,  A/Q cut: %.4f-->%.4f\n", j, ZCut[0], ZCut[1], AoQCut[0], AoQCut[1]);
		// // printf("Iteration#%d: Z cut: %.2f +- %.2f;  A/Q cut: %.4f +- %.4f\n", j, cenZ, wZ, cenAoQ, wAoQ);
				
		// 	nPoin=grZAoQ->GetN();
		// 	double *zVal=grZAoQ->GetY();
		// 	double *aoqVal=grZAoQ->GetX();
		// 	grZAoQ->Set(0);
		// 	k=0;
		// 	for(i=0; i<nPoin; i++)
		// 		if(pow((zVal[i]-cenZ)/wZ,2)+pow((aoqVal[i]-cenAoQ)/wAoQ,2)<=1)
		// 		// if(zVal[i]>ZCut[0]&&zVal[i]<ZCut[1]&&aoqVal[i]>AoQCut[0]&&aoqVal[i]<AoQCut[1])
		// 		{
		// 			grZAoQ->SetPoint(k++, aoqVal[i], zVal[i]);
		// 		}
		// 	if(k<NMIN)
		// 		break;
		// }
		// if(grZAoQ->GetN()<NMIN)
		// {
		// 	delete grZAoQ;
		// 	continue;
		// }
		delete grZAoQ;
		
		fprintf(fCutInfo, "%d  %d  %s  %d  %d  %d  %.2f %.2f  %.4f %.4f  %.3f %.3f  %.6f\n", isCal[iNucl], idN[iNucl], strNucl[iNucl].c_str(), ZPid[iNucl], APid[iNucl], QPid[iNucl], ZCut[0], ZCut[1], AoQCut[0], AoQCut[1], MoQ[iNucl][0], MoQ[iNucl][1], tofPivot[iNucl]);	
	}
	fclose(fCutInfo);
	// delete fCluPid;
	delete fCut;
}

int main(int argc, char** argv)
{
	FindCutZAoQ();
	return 0;
}