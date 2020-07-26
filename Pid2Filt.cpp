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
#include "TLinearFitter.h"
#include "TFormula.h"
#include "TH2F.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TLegend.h"
#include "TMath.h"

void Pid2Filt()
{	
	gStyle->SetOptStat("nemri");
	// gStyle->SetOptStat(0);
	gStyle->SetStatFormat(".6f");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1111);
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);
	
	const int MINPOIN=100, N=200;
	const double NSGM=3;
	const int NDIM=19, NAna=2;
	string sAna[3]={"CfdTdc", "TacNoClk", "TacClk"};
	int i=0, j=0, k=0, iAna=0, iData=0;
	
	TFile *fFilt=new TFile("fFilt.root", "RECREATE");
	TTree *tFilt=new TTree("tFilt", "tree for fit analysis");
	int isCal=0, Z=0, A=0, Q=0;
	bool good[NAna]={0};
	double MoQ[2]={0}, tof[NAna]={0}, xMCP=0, yMCP=0, egyPMT[2][4]={0}, xPlaQ[2]={0}, yPlaQ[2]={0}, xPlaT[2]={0}, yPlaT[2]={0}, delE0=0;
	tFilt->Branch("isCal", &isCal, "isCal/I");
	tFilt->Branch("Z", &Z, "Z/I");
	tFilt->Branch("A", &A, "A/I");
	tFilt->Branch("Q", &Q, "Q/I");
	tFilt->Branch("MoQ", MoQ, "MoQ[2]/D");
	for(iAna=0; iAna<NAna; iAna++)
	{
		tFilt->Branch(("good_"+sAna[iAna]).c_str(), &good[iAna], ("good_"+sAna[iAna]+"/O").c_str());
		tFilt->Branch(("tof_"+sAna[iAna]).c_str(), &tof[iAna], ("tof_"+sAna[iAna]+"/D").c_str());
	}
	tFilt->Branch("xMCP", &xMCP, "xMCP/D");
	tFilt->Branch("yMCP", &yMCP, "yMCP/D");
	tFilt->Branch("egyPMT", egyPMT, "egyPMT[2][4]/D");
	tFilt->Branch("xPlaQ", xPlaQ, "xPlaQ[2]/D");
	tFilt->Branch("yPlaQ", yPlaQ, "yPlaQ[2]/D");
	tFilt->Branch("xPlaT", xPlaT, "xPlaT[2]/D");
	tFilt->Branch("yPlaT", yPlaT, "yPlaT[2]/D");
	tFilt->Branch("delE0", &delE0, "delE0/D");
	
	int isForCal[N]={0}, ZOrig[N]={0}, AOrig[N]={0}, QOrig[N]={0};
	double MoQRef[N][2]={0};
	int iNucl=0, nNucl=0;
	string strRead;
	ifstream fCutInfo("fCutInfo.dat");
	while(!fCutInfo.eof()&&fCutInfo.peek()!=EOF)
	{
		if(fCutInfo.peek()!='#')
		{
			fCutInfo>>isForCal[nNucl]>>strRead>>strRead>>ZOrig[nNucl]>>AOrig[nNucl]>>QOrig[nNucl]>>strRead>>strRead>>strRead>>strRead>>MoQRef[nNucl][0]>>MoQRef[nNucl][1];				
			getline(fCutInfo, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo, strRead);
	}
	fCutInfo.close();
	
	TFile *fPid=new TFile("/home/kailong/ExpData/Jul2018/Pid/pid-run-150--153__270--385.root");
	TTree *tPid;
	fPid->GetObject("tPid",tPid);
	tPid->SetEstimate(-1);
	
	string sDraw="ZPid:APid:QPid";
	for(iAna=0; iAna<NAna; iAna++)
	{
		sDraw+=":good_"+sAna[iAna];
		sDraw+=":tof_"+sAna[iAna];
	}
	sDraw+=":xMCP:yMCP:egyPMT[0][0]:egyPMT[0][1]:egyPMT[0][2]:egyPMT[0][3]:egyPMT[1][0]:egyPMT[1][1]:egyPMT[1][2]:egyPMT[1][3]:xPlaQ[0]:yPlaQ[0]:xPlaQ[1]:yPlaQ[1]:xPlaT_CfdTdc[0]:yPlaT_CfdTdc[0]:xPlaT_CfdTdc[1]:yPlaT_CfdTdc[1]:delE[0]";
	string sCut="good_"+sAna[0]+"==1";
	for(iAna=1; iAna<NAna; iAna++)
		sCut+="||good_"+sAna[iAna]+"==1";
	int nData=tPid->Draw(sDraw.c_str(), sCut.c_str(), "goff");
	double *ZVal=tPid->GetVal(0);
	double *AVal=tPid->GetVal(1);
	double *QVal=tPid->GetVal(2);
	double *goodVal[NAna], *tofVal[NAna];
	for(iAna=0; iAna<NAna; iAna++)
	{
		goodVal[iAna]=tPid->GetVal(3+2*iAna);
		tofVal[iAna]=tPid->GetVal(3+2*iAna+1);
	}
	double *xVal[NDIM];
	for(i=0; i<NDIM; i++)
		xVal[i]=tPid->GetVal(3+2*NAna+i);
	
	TGraph *grTofX[NDIM];
	for(i=0; i<NDIM; i++)
		grTofX[i]=new TGraph();
	double xRng[NAna][NDIM][2]={0}, tofRng[NAna][2]={0};
	bool goodTofX=0, goodRng=0;
	int nSel=0;
	for(iNucl=0; iNucl<nNucl; iNucl++)
	{		
		memset(xRng, 0, sizeof(xRng));
		memset(tofRng, 0, sizeof(tofRng));
		for(iAna=0; iAna<NAna; iAna++)
		{
			double xMean[NDIM]={0}, xRms[NDIM]={0}, tofMean=0, tofRms=0;
			for(i=0; i<NDIM; i++)
				grTofX[i]->Set(0);
			nSel=0;
			goodTofX=0;
			for(iData=0; iData<nData; iData++)
			{
				goodTofX=abs(goodVal[iAna][iData]-1)<1E-3 && abs(ZVal[iData]-ZOrig[iNucl])<1E-3 && abs(AVal[iData]-AOrig[iNucl])<1E-3 && abs(QVal[iData]-QOrig[iNucl])<1E-3 && abs(tofVal[iAna][iData])>0;
				
				if(!goodTofX)
					continue;
				for(i=0; i<NDIM; i++)
				{
					goodTofX=abs(xVal[i][iData])>0;
					if(!goodTofX)
						break;
				}
				if(goodTofX)
				{
					for(i=0; i<NDIM; i++)
						grTofX[i]->SetPoint(nSel, xVal[i][iData], tofVal[iAna][iData]);
					nSel++;
				}				
			}
			if(nSel<MINPOIN)
				continue;
			
			for(i=0; i<NDIM; i++)
			{
				xMean[i]=grTofX[i]->GetMean(1);
				xRms[i]=grTofX[i]->GetRMS(1);
				xRng[iAna][i][0]=xMean[i]-NSGM*xRms[i];
				xRng[iAna][i][1]=xMean[i]+NSGM*xRms[i];
			}
			tofMean=grTofX[0]->GetMean(2);
			tofRms=grTofX[0]->GetRMS(2);
			tofRng[iAna][0]=tofMean-NSGM*tofRms;
			tofRng[iAna][1]=tofMean+NSGM*tofRms;
		}
		for(iData=0; iData<nData; iData++)
		{
			isCal=0;
			Z=0;
			A=0;
			Q=0;
			MoQ[0]=0;  MoQ[1]=0;
			for(iAna=0; iAna<NAna; iAna++)
			{
				good[iAna]=0;
				tof[iAna]=0;
			}
			xMCP=0;
			yMCP=0;
			for(i=0; i<2; i++)
			{
				for(j=0; j<4; j++)
					egyPMT[i][j]=0;
				xPlaQ[i]=0;
				yPlaQ[i]=0;
				xPlaT[i]=0;
				yPlaT[i]=0;
			}
			delE0=0;
			goodTofX=0;
			for(iAna=0; iAna<NAna; iAna++)
			{
				goodTofX=abs(goodVal[iAna][iData]-1)<1E-3 && abs(ZVal[iData]-ZOrig[iNucl])<1E-3 && abs(AVal[iData]-AOrig[iNucl])<1E-3 && abs(QVal[iData]-QOrig[iNucl])<1E-3 && abs(tofVal[iAna][iData])>0;
				if(!goodTofX)
					continue;
				for(i=0; i<NDIM; i++)
				{
					goodTofX=abs(xVal[i][iData])>0;
					if(!goodTofX)
						break;
				}
				if(goodTofX && tofVal[iAna][iData]>tofRng[iAna][0] && tofVal[iAna][iData]<tofRng[iAna][1])
				{
					good[iAna]=1;
					for(i=0; i<NDIM; i++)
					{
						good[iAna]=xVal[i][iData]>xRng[iAna][i][0] && xVal[i][iData]<xRng[iAna][i][1];
						if(!good[iAna])
							break;
					}
				}
			}
			goodRng=good[0];
			for(iAna=1; iAna<NAna; iAna++)
				goodRng=goodRng || good[iAna];
			k=0;
			if(goodRng)
			{
				isCal=isForCal[iNucl];
				Z=ZOrig[iNucl];
				A=AOrig[iNucl];
				Q=QOrig[iNucl];
				MoQ[0]=MoQRef[iNucl][0];
				MoQ[1]=MoQRef[iNucl][1];
				for(iAna=0; iAna<NAna; iAna++)
					tof[iAna]=tofVal[iAna][iData];
				xMCP=xVal[k++][iData];
				yMCP=xVal[k++][iData];
				for(i=0; i<2; i++)
					for(j=0; j<4; j++)
						egyPMT[i][j]=xVal[k++][iData];
				for(i=0; i<2; i++)
				{
					xPlaQ[i]=xVal[k++][iData];
					yPlaQ[i]=xVal[k++][iData];
				}
				for(i=0; i<2; i++)
				{
					xPlaT[i]=xVal[k++][iData];
					yPlaT[i]=xVal[k++][iData];
				}
				delE0=xVal[k++][iData];
				tFilt->Fill();
			}
		}
	}
	for(i=0; i<NDIM; i++)
		delete grTofX[i];
	fFilt->cd();
	tFilt->Write();
	fFilt->Close();
	fPid->Close();
	delete fPid;
	delete fFilt;
}

int main(int argc, char** argv)
{
	Pid2Filt();
	return 0;
}