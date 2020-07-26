#include <iostream>
#include <string>
#include <cstdio>
#include <fstream>
#include <cmath>

using namespace std;

#include "TApplication.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TPad.h"

void PlotCutZAoQChange()
{	
	gStyle->SetOptStat("nemri");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);
	gStyle->SetPalette(1);
	
	const int N=200;
	const int NCUT=8;
	
	string sCan, sDraw, sCut, sAna;
	int i;	
	
	int iNucl, nNucl;
	string strNucl[N];
	string strRead;	
	ifstream fCutInfo("fCutInfo.dat");	
	nNucl=0;	
	while(!fCutInfo.eof()&&fCutInfo.peek()!=EOF)
	{
		if(fCutInfo.peek()!='#')
		{
			fCutInfo>>strRead>>strRead>>strNucl[nNucl];
			getline(fCutInfo, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo, strRead);
	}
	fCutInfo.close();
	
	TFile *fCut=new TFile("fCut.root", "UPDATE");
	
	string strCut[NCUT]={"_0", "_PlusCutMcp", "_PlusCutDelE012", "_PlusCutZQ", "_PlusCutHodo", "_PlusCutDelE3", "_PlusCutDelE4", "_PlusCutZdAm3Z"};
	double iCut[NCUT]={0}, counts[NCUT]={0}, meanZ[NCUT]={0}, sgmZ[NCUT]={0}, meanAoQ[NCUT]={0}, sgmAoQ[NCUT]={0};
	for(i=0; i<NCUT; i++)
		iCut[i]=i;
	
	string strDir="/home/kailong/ExpData/Jul2018/Pid/Graphs/PID/PIDCutChange";
	if(access(strDir.c_str(), F_OK)!=0)
		system(("mkdir "+strDir).c_str());

	for(iNucl=0; iNucl<nNucl; iNucl++)
	// for(iNucl=0; iNucl<10; iNucl++)
	{
		memset(counts, 0, sizeof(counts));
		memset(meanZ, 0, sizeof(meanZ));
		memset(sgmZ, 0, sizeof(sgmZ));
		memset(meanAoQ, 0, sizeof(meanAoQ));
		memset(sgmAoQ, 0, sizeof(sgmAoQ));
		
		for(i=0; i<NCUT; i++)
		{
			TCanvas *cZAoQCut;
			fCut->GetObject(("cZAoQCut_"+strNucl[iNucl]+strCut[i]).c_str(), cZAoQCut);
			cZAoQCut->Draw();
			cZAoQCut->cd(2);
			TH2D *hZAoQ_PidCut=(TH2D*)(gPad->GetPrimitive(("hZAoQ_PidCut_"+strNucl[iNucl]).c_str()));
			
			counts[i]=hZAoQ_PidCut->GetEntries();
			meanZ[i]=hZAoQ_PidCut->GetMean(2);
			sgmZ[i]=hZAoQ_PidCut->GetStdDev(2);
			meanAoQ[i]=hZAoQ_PidCut->GetMean(1);
			sgmAoQ[i]=hZAoQ_PidCut->GetStdDev(1);
			
			cZAoQCut->Close();
		}
		
		TGraph *grCounts=new TGraph(NCUT, iCut, counts);
		grCounts->SetNameTitle(("grCounts_"+strNucl[iNucl]).c_str(), ("grCounts_"+strNucl[iNucl]+";cut# ;counts").c_str());
		
		TGraph *grMeanZ=new TGraph(NCUT, iCut, meanZ);
		grMeanZ->SetNameTitle(("grMeanZ_"+strNucl[iNucl]).c_str(), ("grMeanZ_"+strNucl[iNucl]+";cut# ;Z_c").c_str());
		
		TGraph *grSgmZ=new TGraph(NCUT, iCut, sgmZ);
		grSgmZ->SetNameTitle(("grSgmZ_"+strNucl[iNucl]).c_str(), ("grSgmZ_"+strNucl[iNucl]+";cut# ;#sigma_Z").c_str());
		
		TGraph *grMeanAoQ=new TGraph(NCUT, iCut, meanAoQ);
		grMeanAoQ->SetNameTitle(("grMeanAoQ_"+strNucl[iNucl]).c_str(), ("grMeanAoQ_"+strNucl[iNucl]+";cut# ;A/Q_c").c_str());
		
		TGraph *grSgmAoQ=new TGraph(NCUT, iCut, sgmAoQ);
		grSgmAoQ->SetNameTitle(("grSgmAoQ_"+strNucl[iNucl]).c_str(), ("grSgmAoQ_"+strNucl[iNucl]+";cut# ;#sigma_{A/Q}").c_str());
		
		sCan="cCutChange_"+strNucl[iNucl];
		TCanvas *cCutChange=new TCanvas(sCan.c_str(), sCan.c_str(), 1800, 900);
		cCutChange->Divide(3,2);
		cCutChange->cd(1);
		grMeanZ->Draw("APL*");
		cCutChange->cd(4);
		grSgmZ->Draw("APL*");
		cCutChange->cd(2);
		grMeanAoQ->Draw("APL*");
		cCutChange->cd(5);
		grSgmAoQ->Draw("APL*");
		cCutChange->cd(3);
		grCounts->Draw("APL*");
		
		cCutChange->Write("", TObject::kOverwrite);
		cCutChange->SaveAs((strDir+"/"+sCan+".png").c_str());
		
		delete cCutChange;
		delete grSgmAoQ;
		delete grMeanAoQ;
		delete grSgmZ;
		delete grMeanZ;
		delete grCounts;
	}
	fCut->Close();
	delete fCut;
}

void StandaloneApplication(int argc, char** argv)
{
	PlotCutZAoQChange();
}

int main(int argc, char** argv)
{
	TApplication app("ROOT Application", &argc, argv);
	StandaloneApplication(app.Argc(), app.Argv());
	// app.Run();
	return 0;
}