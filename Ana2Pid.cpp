#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <unistd.h>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TStyle.h"
#include "TH2F.h"

using namespace std;

struct StrtAna
{
	double tof[4]; //[4]: [0] TAC+ADC+clock; [1] TAC+ADC; [2]: regular CFD+TDC; [3] MCFD16+TDC
	double tD[4][8][8];
	double egyPMT[2][4]; //energies in 4 PMTs and each plastic
	double egyPla[2];	//[2]: [0] is energy in plastic at S800, [1] is energy in plastic at A1900
	double xPlaT[4][2];  //[2]: [0] for S800 plastic from time info, [1] for A1900 plastic from time info
	double yPlaT[4][2];
	double xPlaQ[2]; //[2]: [0] for S800 plastic from amp info, [1] for A1900 plastic from amp info
	double yPlaQ[2];
	double xMCP[2]; //two gain settings
	double yMCP[2]; //two gain settings
	double delE[5];
	double tke;
	double beta[4];
	double gamma[4];
	double Z[4];
	double dZ[4];
	double brho[2];
	double AoQ[4][2];
	double Q[4][2];
	double ZmQ[4][2];
	double ZImQ[4][2];
	double A[4][2];
	double Araw[4][2];
	double Am2Q[4][2];
	double Am3Q[4][2];
	double Am2Z[4][2];
	double Am3Z[4][2];
	double dAm2Z[4][2];
	double dAm3Z[4][2];
	int numTime[4][2];
	int sigTime[4][2];
	int Zi[4];
	int Am2Zi[4][2];
	int Am3Zi[4][2];
	double xCrdc[2][2]; //first [2]: two CRDCs; second [2]: [0] is from gravity center; [1] is from Gaussian fit
	double xEgyCrdc[2][256]; //256 energies for 256 cathode pads
	double yCrdc[2];	//y is only considered from electron's drift time
	int mulHodo;
	double egyHodo[32];
};

struct StrtPid
{
	double tof;
	double beta;
	double gamma;
	double Z;
	double dZ;
	double AoQ;
	double Q;
	double ZmQ;
	double ZImQ;
	double A;
	double Araw;
	double Am2Q;
	double Am3Q;
	double Am2Z;
	double Am3Z;
	double dAm2Z;
	double dAm3Z;
	int Zi;
	int Am2Zi;
	int Am3Zi;
};

void Ana2Pid()
{
	gStyle->SetOptStat("nemri");
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptFit(1);
	gStyle->SetCanvasDefH(900);
	gStyle->SetCanvasDefW(1200);

	const int N = 200;
	const int NAna = 3;
	string sAna[NAna] = {"CfdTdc", "TacNoClk", "TacClk"};
	int i, j, k;
	string strRead;

	int iNucl, nNucl;
	string strNucl[N];
	int isCal[N]={0}, ZOrig[N] = {0}, AOrig[N] = {0}, QOrig[N] = {0};
	double ZCut[N][2] = {0}, AoQCut[N][2] = {0}, MoQOrig[N][2] = {0};

	nNucl = 0;
	ifstream fCutInfo("fCutInfo.dat");
	// printf("\n********************information from fCutInfo.dat********************\n");
	while (!fCutInfo.eof() && fCutInfo.peek() != EOF)
	{
		if (fCutInfo.peek() != '#')
		{
			fCutInfo >> isCal[nNucl] >> strRead >> strNucl[nNucl] >> ZOrig[nNucl] >> AOrig[nNucl] >> QOrig[nNucl] >> ZCut[nNucl][0] >> ZCut[nNucl][1] >> AoQCut[nNucl][0] >> AoQCut[nNucl][1] >> MoQOrig[nNucl][0] >> MoQOrig[nNucl][1];
			// printf("%s %d %d %d  %f %f  %f %f  %f %f\n", strNucl[nNucl].c_str(), ZOrig[nNucl], AOrig[nNucl], QOrig[nNucl], ZCut[nNucl][0], ZCut[nNucl][1], AoQCut[nNucl][0], AoQCut[nNucl][1], MoQOrig[nNucl][0], MoQOrig[nNucl][1]);
			getline(fCutInfo, strRead);
			nNucl++;
		}
		else
			getline(fCutInfo, strRead);
	}
	fCutInfo.close();
	
	TFile *fCut = new TFile("fCut.root");
	TH2F *hMCP_y_x_cut, *h_dE1_dE0_cut, *h_dE2_dE1_cut, *h_dE2_dE0_cut, *hS800PlaXQ_T_cut, *hS800PlaYQ_T_cut, *hZ_ZmQ_Pid_cut, *hZ_Q_Pid_cut;
	fCut->GetObject("hMCP_y_x_cut", hMCP_y_x_cut);
	fCut->GetObject("h_dE1_dE0_cut", h_dE1_dE0_cut);
	fCut->GetObject("h_dE2_dE1_cut", h_dE2_dE1_cut);
	fCut->GetObject("h_dE2_dE0_cut", h_dE2_dE0_cut);
	fCut->GetObject("hS800PlaXQ_T_cut", hS800PlaXQ_T_cut);
	fCut->GetObject("hS800PlaYQ_T_cut", hS800PlaYQ_T_cut);
	fCut->GetObject("hZ_ZmQ_Pid_cut", hZ_ZmQ_Pid_cut);
	fCut->GetObject("hZ_Q_Pid_cut", hZ_Q_Pid_cut);

	TFile *fPid = new TFile("/home/kailong/ExpData/Jul2018/Pid/pid-run-150--153__270--385.root", "RECREATE");
	TTree *tPid = new TTree("tPid", "tree for pid analysis");
	int runNum=0, ZPid=0, APid=0, QPid=0, isForCal=0;
	double AoQPid=0, MoQPid[2]={0}, xMCP=0, yMCP=0, delE[5]={0}, egyPla[2]={0}, egyPMT[2][4]={0}, xPlaQ[2]={0}, yPlaQ[2]={0}, xCrdc[2]={0}, yCrdc[2]={0}, xPlaCrdcS800=0, yPlaCrdcS800=0;
	string eSet, nucl;

	tPid->Branch("eSet", &eSet);
	tPid->Branch("runNum", &runNum, "runNum/I");
	tPid->Branch("nucl", &nucl);
	tPid->Branch("ZPid", &ZPid, "ZPid/I");
	tPid->Branch("APid", &APid, "APid/I");
	tPid->Branch("QPid", &QPid, "QPid/I");
	tPid->Branch("isForCal", &isForCal, "isForCal/I");
	tPid->Branch("AoQPid", &AoQPid, "AoQPid/D");
	tPid->Branch("MoQPid", MoQPid, "MoQPid[2]/D");
	tPid->Branch("xMCP", &xMCP, "xMCP/D");
	tPid->Branch("yMCP", &yMCP, "yMCP/D");
	tPid->Branch("delE", delE, "delE[5]/D");
	tPid->Branch("egyPla", egyPla, "egyPla[2]/D");
	tPid->Branch("egyPMT", egyPMT, "egyPMT[2][4]/D");
	tPid->Branch("xPlaQ", xPlaQ, "xPlaQ[2]/D");
	tPid->Branch("yPlaQ", yPlaQ, "yPlaQ[2]/D");
	tPid->Branch("xCrdc", xCrdc, "xCrdc[2]/D");
	tPid->Branch("yCrdc", yCrdc, "yCrdc[2]/D");
	tPid->Branch("xPlaCrdcS800", &xPlaCrdcS800, "xPlaCrdcS800/D");
	tPid->Branch("yPlaCrdcS800", &yPlaCrdcS800, "yPlaCrdcS800/D");

	double Z[NAna]={0}, AoQ[NAna]={0};
	double tof[NAna]={0}, tofPmt[NAna][5][9]={0};  //[5][9]: different PMT combinations: [1--4][5--8] for ToFs from PMTs#1--#4 and PMTs#5--#8; [0][1] ToF from PMTs#1#3#5#7; [1][0] Tof from PMTs#2#4#6#8
	double xPlaT[NAna][2]={0}, yPlaT[NAna][2]={0};
	int good[NAna]={0};

	int iAna;
	for (iAna = 0; iAna < NAna; iAna++)
	{
		tPid->Branch(("good_" + sAna[iAna]).c_str(), &good[iAna], ("good_" + sAna[iAna] + "/I").c_str());
		tPid->Branch(("Z_" + sAna[iAna]).c_str(), &Z[iAna], ("Z_" + sAna[iAna] + "/D").c_str());
		tPid->Branch(("AoQ_" + sAna[iAna]).c_str(), &AoQ[iAna], ("AoQ_" + sAna[iAna] + "/D").c_str());
		tPid->Branch(("tof_" + sAna[iAna]).c_str(), &tof[iAna], ("tof_" + sAna[iAna] + "/D").c_str());
		tPid->Branch(("tofPmt_" + sAna[iAna]).c_str(), tofPmt[iAna], ("tofPmt_" + sAna[iAna] + "[5][9]/D").c_str());
		tPid->Branch(("xPlaT_" + sAna[iAna]).c_str(), xPlaT[iAna], ("xPlaT_" + sAna[iAna] + "[2]/D").c_str());
		tPid->Branch(("yPlaT_" + sAna[iAna]).c_str(), yPlaT[iAna], ("yPlaT_" + sAna[iAna] + "[2]/D").c_str());
	}

	TFile *fAna = new TFile("/home/kailong/ExpData/Jul2018/AnaData/ana-run-150--153__270--385.root");
	TTree *tAna;
	fAna->GetObject("tAna", tAna);
	int run = 0;
	string *setting = new string;
	StrtAna ana;
	StrtPid pid;
	memset(&ana, 0, sizeof(ana));
	memset(&pid, 0, sizeof(pid));
	tAna->SetBranchAddress("run", &run);
	tAna->SetBranchAddress("setting", &setting);
	tAna->SetBranchAddress("ana", &ana);
	tAna->SetBranchAddress("pid", &pid);

	long long iEnt = 0, nEnt=tAna->GetEntries();
	for (iEnt = 0; iEnt < nEnt; iEnt++)
	{
		tAna->GetEntry(iEnt);
		eSet = *setting;
		runNum = run;
		nucl="";
		ZPid = 0;
		APid = 0;
		QPid = 0;
		isForCal = 0;
		AoQPid = 0;
		xMCP = 0;
		yMCP = 0;
		memset(MoQPid, 0, sizeof(MoQPid));
		memset(delE, 0, sizeof(delE));
		memset(egyPla, 0, sizeof(egyPla));
		memset(egyPMT, 0, sizeof(egyPMT));
		memset(xPlaQ, 0, sizeof(xPlaQ));
		memset(yPlaQ, 0, sizeof(yPlaQ));
		memset(good, 0, sizeof(good));
		memset(xCrdc, 0, sizeof(xCrdc));
		memset(yCrdc, 0, sizeof(yCrdc));
		xPlaCrdcS800=0;
		yPlaCrdcS800=0;
		memset(Z, 0, sizeof(Z));
		memset(AoQ, 0, sizeof(AoQ));
		memset(tof, 0, sizeof(tof));
		memset(tofPmt, 0, sizeof(tofPmt));
		memset(xPlaT, 0, sizeof(xPlaT));
		memset(yPlaT, 0, sizeof(yPlaT));

		for (iNucl = 0; iNucl < nNucl; iNucl++)
			// if(run==150||run==152||run==153||run==384||run==385||eSet=="PS_270_382"||eSet=="RS_270_382")
			if (eSet == "PS_270_382" || eSet == "RS_270_382")
			// if (eSet == "PS_270_382")
			// if (eSet == "RS_270_382")
			{
				double wZ = (ZCut[iNucl][1] - ZCut[iNucl][0]) / 2;
				double wAoQ = (AoQCut[iNucl][1] - AoQCut[iNucl][0]) / 2;
				double cenAoQ = (AoQCut[iNucl][1] + AoQCut[iNucl][0]) / 2;
				double cenZ = (ZCut[iNucl][1] + ZCut[iNucl][0]) / 2;

				if (pow((pid.Z - cenZ) / wZ, 2) + pow((pid.AoQ - cenAoQ) / wAoQ, 2) <= 1 && ana.mulHodo == 0 && hMCP_y_x_cut->GetBinContent(hMCP_y_x_cut->FindFixBin(ana.xMCP[1], ana.yMCP[1])) > 0 && h_dE1_dE0_cut->GetBinContent(h_dE1_dE0_cut->FindFixBin(ana.delE[0], ana.delE[1])) > 0 && h_dE2_dE1_cut->GetBinContent(h_dE2_dE1_cut->FindFixBin(ana.delE[1], ana.delE[2])) > 0 && h_dE2_dE0_cut->GetBinContent(h_dE2_dE0_cut->FindFixBin(ana.delE[0], ana.delE[2])) > 0 && ana.delE[3] > 10 && ana.delE[4] < 10 && hS800PlaXQ_T_cut->GetBinContent(hS800PlaXQ_T_cut->FindFixBin(ana.xPlaT[2][0], ana.xPlaQ[0])) > 0 && hS800PlaYQ_T_cut->GetBinContent(hS800PlaYQ_T_cut->FindFixBin(ana.yPlaT[2][0], ana.yPlaQ[0])) > 0)
				{
					nucl = strNucl[iNucl];
					ZPid = ZOrig[iNucl];
					APid = AOrig[iNucl];
					QPid = QOrig[iNucl];
					isForCal = isCal[iNucl];
					AoQPid = 1.0 * AOrig[iNucl] / QOrig[iNucl];
					MoQPid[0] = MoQOrig[iNucl][0];
					MoQPid[1] = MoQOrig[iNucl][1];

					if (abs(ana.xMCP[1]) > 0 && abs(ana.yMCP[1]) > 0)
					{
						xMCP = ana.xMCP[1];
						yMCP = ana.yMCP[1];
					}

					for (i = 0; i < 5; i++)
						delE[i] = ana.delE[i];

					for (i = 0; i < 2; i++)
					{
						egyPla[i] = ana.egyPla[i];
						for (j = 0; j < 4; j++)
							egyPMT[i][j] = ana.egyPMT[i][j];
					}

					if (abs(ana.xPlaQ[0]) > 0 && abs(ana.yPlaQ[0]) > 0)
					{
						xPlaQ[0] = ana.xPlaQ[0];
						yPlaQ[0] = ana.yPlaQ[0];
					}

					if (abs(ana.xPlaQ[1]) > 0 && abs(ana.yPlaQ[1]) > 0)
					{
						xPlaQ[1] = ana.xPlaQ[1];
						yPlaQ[1] = ana.yPlaQ[1];
					}

					for (i = 0; i < 2; i++)
						if (ana.xCrdc[i][0] > 0 && ana.yCrdc[i] > 0)
						{
							xCrdc[i] = ana.xCrdc[i][0];
							yCrdc[i] = ana.yCrdc[i];
						}
					if(xCrdc[0]>0 && xCrdc[1]>0 && yCrdc[0]>0 && yCrdc[1]>0)
					{
						xPlaCrdcS800=(97*xCrdc[0]+10*xCrdc[1])/107;  //unit: pad;  //Distance from CRDC1 to CRDC2 is 107 cm, and assume distance from CRDC1 to plastic is 10 cm
						yPlaCrdcS800=(97*yCrdc[0]+10*yCrdc[1])/107;  //unit: channel
					}

					for (iAna = 0; iAna < NAna; iAna++)
					{
						j = 2 - iAna;
						if (ana.numTime[j][0] == 4 && ana.numTime[j][1] == 4 && abs(ana.Z[j]) > 0 && abs(ana.AoQ[j][1]) > 0 && abs(ana.tof[j]) > 0)
						{
							good[iAna] = 1;
							Z[iAna] = ana.Z[j];
							AoQ[iAna] = ana.AoQ[j][1];
							tof[iAna] = ana.tof[j];

							if(iAna==1)
							{
								for(i=0; i<4; i++)
									tofPmt[iAna][i+1][i+5] = ana.tD[j][i][i + 4] / 1000 + 550;
								tofPmt[iAna][0][1]=(ana.tD[j][0][4]+ana.tD[j][2][6])/2.0/1000 + 550;
								tofPmt[iAna][1][0]=(ana.tD[j][1][5]+ana.tD[j][3][7])/2.0/1000 + 550;
							}
							else
							{
								for(i=0; i<4; i++)
									for(k=4; k<8; k++)
										tofPmt[iAna][i+1][k+1] = ana.tD[j][i][k] / 1000 + 550;
								tofPmt[iAna][0][1]=(ana.tD[j][0][4]+ana.tD[j][2][6]+ana.tD[j][0][6]+ana.tD[j][2][4])/4.0/1000 + 550;
								tofPmt[iAna][1][0]=(ana.tD[j][1][5]+ana.tD[j][3][7]+ana.tD[j][1][7]+ana.tD[j][3][5])/4.0/1000 + 550;
							}
							if (abs(ana.xPlaT[j][0]) > 0 && abs(ana.yPlaT[j][0]) > 0)
							{
								xPlaT[iAna][0] = ana.xPlaT[j][0];
								yPlaT[iAna][0] = ana.yPlaT[j][0];
							}
							if (abs(ana.xPlaT[j][1]) > 0 && abs(ana.yPlaT[j][1]) > 0)
							{
								xPlaT[iAna][1] = ana.xPlaT[j][1];
								yPlaT[iAna][1] = ana.yPlaT[j][1];
							}
						}
					}
					if (good[0]==1 || good[1]==1 || good[2]==1)
						tPid->Fill();
				}
			}
	}
	fPid->cd();
	tPid->Write();
	fPid->Close();
	fAna->Close();
	delete fAna;
	delete setting;
	delete fPid;
	delete fCut;
} //end of whole function

int main(int argc, char **argv)
{
	Ana2Pid();
	return 0;
}