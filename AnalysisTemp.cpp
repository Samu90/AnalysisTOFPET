#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TSystem.h"

#include "ReadFileList.cc"
#include "FunctionsTemp.cc"


void Analysis(){

  gROOT->Reset();
  gROOT->SetBatch(kTRUE);

  gSystem->Exec("mkdir Plot");
  gSystem->Exec("mkdir Plot/EnergyTemp");
  gSystem->Exec("mkdir Plot/EnergyTemp/Projections/");

  vector<string> FileList;
  FileList=ReadData("TestStability3/fileov7.txt");

  Int_t Nch=2;
  Int_t EnergyBin=100,TempBin=15;
  Double_t EMin=0,EMax=100,TempMin=26,TempMax=35;

  bool printOpenedFileNumber=true;
  bool PrintVectorTGraph=true;
  
  TFile* DataFile[3]; //File dei dati
  TTree* DataTree[3]; //Tree dei dati
  TH1D* Histo[(int)FileList.size()][Nch]; // [2] canali 59 e 315

  TH2D* HistoEnergyVsTemp[2]; 
  HistoEnergyVsTemp[0] = new TH2D("HistoEnergyVsTempCh59","HistoEnergyVsTempCh59",TempBin,TempMin,TempMax,EnergyBin,EMin,EMax);
  HistoEnergyVsTemp[1] = new TH2D("HistoEnergyVsTempCh315","HistoEnergyVsTempCh315",TempBin,TempMin,TempMax,EnergyBin,EMin,EMax);

  Double_t MeanPedCh59,MeanPedCh315; // mean value of pedestal taken before and after Phys for both channel
  
  for(int i=0;i<(int)FileList.size()-2;i+=3){
    
    for(int j=0; j<3;j++){
      DataFile[j]=TFile::Open(("TestStability3/"+FileList.at(i+j)).c_str());
      if(printOpenedFileNumber) std::cout << "FileOpened: " << i+j << std::endl;
      DataTree[j] = (TTree*)DataFile[j]->Get("data");
    }//chiudo for j
    
    Histo[i][0]  =  new TH1D(("HistoPedCh59BNum"+to_string(i+1)).c_str(),("HistoPedCh59BNum"+to_string(i+1)).c_str(),600,0,600);
    Histo[i][1]  =  new TH1D(("HistoPedCh315BNum"+to_string(i+1)).c_str(),("HistoPedCh315BNum"+to_string(i+1)).c_str(),600,0,600);
    GetPedestal(DataTree[0], Histo[i][0], Histo[i][1]);
    
    
    Histo[i+2][0]  =  new TH1D(("HistoPedCh59ANum"+to_string(i+1)).c_str(),("HistoPedCh59ANum"+to_string(i+1)).c_str(),600,0,600);
    Histo[i+2][1]  =  new TH1D(("HistoPedCh315ANum"+to_string(i+1)).c_str(),("HistoPedCh315ANum"+to_string(i+1)).c_str(),600,0,600);
    GetPedestal(DataTree[2],Histo[i+2][0],Histo[i+2][1]);
    
    MeanPedCh59=(Histo[i][0]->GetMean()+Histo[i+2][0]->GetMean())/2;
    MeanPedCh315=(Histo[i][1]->GetMean()+Histo[i+2][1]->GetMean())/2;
    
    FillETemperature(DataTree[1],HistoEnergyVsTemp[0],HistoEnergyVsTemp[1],MeanPedCh59,MeanPedCh315);
  }
  
  Double_t XTemp[TempBin];
  std::vector<std::vector<Double_t> > YTemp(8, std::vector<Double_t>(TempBin));
  //Y1->[0]:Y511ch59,[1]:EY511ch59,[2]Y1275ch59,[3]EY1275ch59
  //    [4]:Y511ch315,[5]:EY511ch315,[6]Y1275ch315,[7]EY1275ch315
  
  GetProfiles(HistoEnergyVsTemp[0],HistoEnergyVsTemp[1],XTemp,YTemp,"EnergyTemp");
  
  TGraphErrors* graphTemp1Ch59 = new TGraphErrors(TempBin,XTemp,&YTemp[0][0],0,&YTemp[1][0]);
  TGraphErrors* graphTemp2Ch59 = new TGraphErrors(TempBin,XTemp,&YTemp[2][0],0,&YTemp[3][0]);
  TGraphErrors* graphTemp1Ch315 = new TGraphErrors(TempBin,XTemp,&YTemp[4][0],0,&YTemp[5][0]);
  TGraphErrors* graphTemp2Ch315 = new TGraphErrors(TempBin,XTemp,&YTemp[6][0],0,&YTemp[7][0]);
  
  if(PrintVectorTGraph){
    for(int i=0;i<8;i++){
      for(int j=0;j<TempBin;j++){
	std::cout<<XTemp[j] <<"    " << YTemp[i][j]<< std::endl;
      }//chiudo for j
    }//chiudo for i
  }//chiudo if
  
  
  TCanvas* canvETemp = new TCanvas("EnergySpectrumVsTemp","EnergySpectrumVsTemp",1200,600);
  canvETemp->Divide(2,1);
  
  canvETemp->cd(1);
  HistoEnergyVsTemp[0]->Draw("COLZ");
  graphTemp1Ch59->Draw("SAMEP");
  graphTemp2Ch59->Draw("SAMEP");
  
  canvETemp->cd(2);
  HistoEnergyVsTemp[1]->Draw("COLZ");
  graphTemp1Ch315->Draw("SAMEP");
  graphTemp2Ch315->Draw("SAMEP");
  
  
  canvETemp->SaveAs("Plot/EnergyTemp/EnergySpectrumVsTemp.png");
  
  
}
