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

#include "ReadFileList.cc"
#include "Functions.cc"

void Analysis(){
  
  gROOT->SetBatch(kTRUE);

  vector<string> FileList;
  FileList=ReadData("TestStability3/fileov7.txt");

  Int_t Nch=2;
  Int_t TimeBin=30,EnergyBin=100;
  Float_t TMin=0,TMax=20000,EMin=0,EMax=100;
  

  bool printFileList= true;
  bool printOpenedFileNumber= true;

  if(printFileList){
    
    for(int i=0;i<(int)FileList.size();i++){
      std::cout<<"File " << i << "   " << FileList.at(i) << std::endl;
    }//chiudo for
    
  }//chiudo if
  



  TFile* DataFile[3]; //File dei dati
  TTree* DataTree[3]; //Tree dei dati
  TH1F* Histo[(int)FileList.size()][Nch]; // [2] canali 59 e 315

  TH2F* HistoEnergyVsTime[2]; 
  HistoEnergyVsTime[0] = new TH2F("HistoEnergyVsTimeCh59","HistoEnergyVsTimeCh59",TimeBin,TMin,TMax,EnergyBin,EMin,EMax);
  HistoEnergyVsTime[1] = new TH2F("HistoEnergyVsTimeCh315","HistoEnergyVsTimeCh315",TimeBin,TMin,TMax,EnergyBin,EMin,EMax);
  
  TCanvas* canvino = new TCanvas("HistoCanvas","HistoCanvas",1200,800);
  canvino->Divide(3,Nch);
  
  Double_t MeanPedCh59,MeanPedCh315; // mean value of pedestal taken before and after Phys for both channel

  for(int i=0;i<(int)FileList.size()-2;i+=3){
    
    for(int j=0; j<3;j++){
      DataFile[j]=TFile::Open(("TestStability3/"+FileList.at(i+j)).c_str());
      if(printOpenedFileNumber) std::cout << "FileOpened: " << i+j << std::endl;
      DataTree[j] = (TTree*)DataFile[j]->Get("data");
    }//chiudo for j

    Histo[i][0]  =  new TH1F(("HistoPedCh59BNum"+to_string(i+1)).c_str(),("HistoPedCh59BNum"+to_string(i+1)).c_str(),600,0,600);
    Histo[i][1]  =  new TH1F(("HistoPedCh315BNum"+to_string(i+1)).c_str(),("HistoPedCh315BNum"+to_string(i+1)).c_str(),600,0,600);
    GetPedestal(DataTree[0], Histo[i][0], Histo[i][1]);
 
    
    Histo[i+2][0]  =  new TH1F(("HistoPedCh59ANum"+to_string(i+1)).c_str(),("HistoPedCh59ANum"+to_string(i+1)).c_str(),600,0,600);
    Histo[i+2][1]  =  new TH1F(("HistoPedCh315ANum"+to_string(i+1)).c_str(),("HistoPedCh315ANum"+to_string(i+1)).c_str(),600,0,600);
    GetPedestal(DataTree[2],Histo[i+2][0],Histo[i+2][1]);
   
    MeanPedCh59=(Histo[i][0]->GetMean()+Histo[i+2][0]->GetMean())/2;
    MeanPedCh315=(Histo[i][1]->GetMean()+Histo[i+2][1]->GetMean())/2;

    Histo[i+1][0]  =  new TH1F(("HistoPhysCh59Num"+to_string(i+1)).c_str(),("HistoPhysCh59Num"+to_string(i+1)).c_str(),600,0,600);
    Histo[i+1][1]  =  new TH1F(("HistoPhysCh315Num"+to_string(i+1)).c_str(),("HistoPhysCh315Num"+to_string(i+1)).c_str(),600,0,600);
    GetSpectrum(DataTree[1],Histo[i+1][0],Histo[i+1][1],MeanPedCh59,MeanPedCh315);

    FillETime(DataTree[1],HistoEnergyVsTime[0],HistoEnergyVsTime[1],MeanPedCh59,MeanPedCh315);

    canvino->cd(1);
    //Histo[i][0]->GetXaxis()->SetRangeUser(Histo[i][0]->GetMean()-7*Histo[i][0]->GetRMS(),Histo[i][0]->GetMean()+7*Histo[i][0]->GetRMS());
    Histo[i][0]->GetXaxis()->SetRangeUser(100,300);
    Histo[i][0]->Draw();
    
    canvino->cd(2);
    //Histo[i][1]->GetXaxis()->SetRangeUser(Histo[i][1]->GetMean()-7*Histo[i][1]->GetRMS(),Histo[i][1]->GetMean()+7*Histo[i][1]->GetRMS());
    Histo[i+2][0]->GetXaxis()->SetRangeUser(100,300);
    Histo[i+2][0]->Draw();
    
    canvino->cd(3);
    Histo[i+1][0]->GetXaxis()->SetRangeUser(0,100);
    Histo[i+1][0]->Draw();
    
    canvino->cd(4);
    //Histo[i+2][0]->GetXaxis()->SetRangeUser(Histo[i+2][0]->GetMean()-7*Histo[i+2][0]->GetRMS(),Histo[i+2][0]->GetMean()+7*Histo[i+2][0]->GetRMS());
    Histo[i][1]->GetXaxis()->SetRangeUser(100,300);
    Histo[i][1]->Draw();
    
    canvino->cd(5);
    //Histo[i+2][1]->GetXaxis()->SetRangeUser(Histo[i+2][1]->GetMean()-7*Histo[i+2][1]->GetRMS(),Histo[i+2][1]->GetMean()+7*Histo[i+2][1]->GetRMS());
    Histo[i+2][1]->GetXaxis()->SetRangeUser(100,300);
    Histo[i+2][1]->Draw();
    
    canvino->cd(6);
    Histo[i+1][1]->GetXaxis()->SetRangeUser(0,100);
    Histo[i+1][1]->Draw();
    
    canvino->SaveAs(("Plot/Plot"+to_string(i+1)+".png").c_str());
  
  }//chiudo for i+=3

  //gROOT->SetBatch(kTRUE);

  TCanvas* canvET = new TCanvas("EnergySpectrumVsTime","EnergySpectrumVsTime",1200,600);
  canvET->Divide(2,1);
  canvET->cd(1);
  HistoEnergyVsTime[0]->Draw("COLZ");
  canvET->cd(2);
  HistoEnergyVsTime[1]->Draw("COLZ");
  
  canvET->SaveAs("Plot/EnergyPlot/EnergySpectrumVsTime.png");

  //gROOT->SetBatch(kTRUE);
  
  Double_t X[TimeBin],Y1[TimeBin],EY1[TimeBin],Y2[TimeBin],EY2[TimeBin];
  
  GetProfiles(HistoEnergyVsTime[0],HistoEnergyVsTime[1],X,Y1,EY1,Y2,EY2);
  
}
