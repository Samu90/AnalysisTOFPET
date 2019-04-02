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
#include "Functions.cc"

void Analysis(string DirData){
  
  gROOT->Reset();
  gROOT->SetBatch(kTRUE);

  gSystem->Exec(("mkdir "+DirData+"/Plot").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/EnergyTime").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/EnergyTime/Projection").c_str());
  //  gSystem->Exec("mkdir Plot/EnergyTemp");
  //  gSystem->Exec("mkdir Plot/EnergyTemp/Projections/");

  //gStyle->SetOptStat("000001000");

  vector<string> FileList;
  FileList=ReadData(DirData+"/fileov7.txt");

  Int_t Nch=2;
  Int_t TimeBin=30,EnergyBin=100,TempBin=40;
  Double_t TMin=0,TMax=20000,EMin=0,EMax=100,TempMin=26,TempMax=35;
  

  bool printFileList= false;
  bool printOpenedFileNumber= true;
  bool PrintVectorTGraph=false;

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
  
  
  /*TH2D* HistoEnergyVsTemp[2]; 
  HistoEnergyVsTemp[0] = new TH2D("HistoEnergyVsTempCh59","HistoEnergyVsTempCh59",TempBin,26,33,EnergyBin,EMin,EMax);
  HistoEnergyVsTemp[1] = new TH2D("HistoEnergyVsTempCh315","HistoEnergyVsTempCh315",TempBin,28,35,EnergyBin,EMin,EMax);*/
  
  
  TCanvas* canvino = new TCanvas("HistoCanvas","HistoCanvas",1200,800);
  canvino->Divide(3,Nch);
  
  Double_t MeanPedCh59,MeanPedCh315; // mean value of pedestal taken before and after Phys for both channel

  
  for(int i=0;i<(int)FileList.size()-2;i+=3){
    
    for(int j=0; j<3;j++){
      DataFile[j]=TFile::Open((DirData+"/"+FileList.at(i+j)).c_str());
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
    //FillETemperature(DataTree[1],HistoEnergyVsTemp[0],HistoEnergyVsTemp[1],MeanPedCh59Bis,MeanPedCh315Bis);
    
    canvino->cd(1);
    Histo[i][0]->GetXaxis()->SetRangeUser(100,300);
    Histo[i][0]->Draw();
    
    canvino->cd(2);
    Histo[i+2][0]->GetXaxis()->SetRangeUser(100,300);
    Histo[i+2][0]->Draw();
    
    canvino->cd(3);
    Histo[i+1][0]->GetXaxis()->SetRangeUser(0,100);
    Histo[i+1][0]->Draw();
    
    canvino->cd(4);
    Histo[i][1]->GetXaxis()->SetRangeUser(100,300);
    Histo[i][1]->Draw();
    
    canvino->cd(5);
    Histo[i+2][1]->GetXaxis()->SetRangeUser(100,300);
    Histo[i+2][1]->Draw();
    
    canvino->cd(6);
    Histo[i+1][1]->GetXaxis()->SetRangeUser(0,100);
    Histo[i+1][1]->Draw();
    
    canvino->SaveAs((DirData+"/Plot/Plot"+to_string(i+1)+".png").c_str());
  
  }//chiudo for i+=3

 
  Double_t XTime[TimeBin];
  std::vector<std::vector<Double_t> > YTime(8, std::vector<Double_t>(TimeBin));
  //Y1->[0]:Y511ch59,[1]:EY511ch59,[2]Y1275ch59,[3]EY1275ch59
  //    [4]:Y511ch315,[5]:EY511ch315,[6]Y1275ch315,[7]EY1275ch315

  GetProfiles(HistoEnergyVsTime[0],HistoEnergyVsTime[1],XTime,YTime,"EnergyTime");
  
   
  TGraphErrors* graphTime1Ch59 = new TGraphErrors(TimeBin,XTime,&YTime[0][0],0,&YTime[1][0]);
  TGraphErrors* graphTime2Ch59 = new TGraphErrors(TimeBin,XTime,&YTime[2][0],0,&YTime[3][0]);
  TGraphErrors* graphTime1Ch315 = new TGraphErrors(TimeBin,XTime,&YTime[4][0],0,&YTime[5][0]);
  TGraphErrors* graphTime2Ch315 = new TGraphErrors(TimeBin,XTime,&YTime[6][0],0,&YTime[7][0]);
      
  if(PrintVectorTGraph){
    for(int i=0;i<8;i++){
      for(int j=0;j<TimeBin;j++){
	std::cout<<XTime[j] <<"    " << YTime[i][j]<< std::endl;
      }//chiudo for j
    }//chiudo for i
  }//chiudo if
  
  
  TCanvas* canvETime = new TCanvas("EnergySpectrumVsTime","EnergySpectrumVsTime",1200,600);
  canvETime->Divide(2,1);

  canvETime->cd(1);
  HistoEnergyVsTime[0]->Draw("COLZ");
  graphTime1Ch59->Draw("SAMEP");
  graphTime2Ch59->Draw("SAMEP");

  canvETime->cd(2);
  HistoEnergyVsTime[1]->Draw("COLZ");
  graphTime1Ch315->Draw("SAMEP");
  graphTime2Ch315->Draw("SAMEP");


  canvETime->SaveAs((DirData+"/Plot/EnergyTime/EnergySpectrumVsTime.png").c_str());

  /*
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
  graphTime1Ch59->Draw("SAMEP");
  graphTime2Ch59->Draw("SAMEP");
  
  canvETemp->cd(2);
  HistoEnergyVsTemp[1]->Draw("COLZ");
  graphTime1Ch315->Draw("SAMEP");
  graphTime2Ch315->Draw("SAMEP");
  
  canvETemp->SaveAs("Plot/EnergyTemp/EnergySpectrumVsTemp.png");
  */
  

}
