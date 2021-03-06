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
  
  gSystem->Exec("mkdir Plot/EnergyTempV2");
  gSystem->Exec("mkdir Plot/EnergyTempV2/Partials/");

  //gStyle->SetOptStat("000001000");

  vector<string> FileListPedestal;
  FileListPedestal=ReadData("TestStability3/PedFile.txt");

  vector<string> FileListPhysics;
  FileListPhysics=ReadData("TestStability3/PhysFile.txt");

  
  //Get Pedestal
  Double_t Pedestal[(int)FileListPhysics.size()][2];
  Double_t PedestalCh1[2];
  Double_t PedestalCh2[2];

  int k=0;

  for(int i=0;i < (int)FileListPedestal.size()-1;i+=2){
    
    TFile* f0= TFile::Open(("TestStability3/"+FileListPedestal.at(i)).c_str());
    TFile* f1= TFile::Open(("TestStability3/"+FileListPedestal.at(i)).c_str());

    TTree* tree0 = (TTree*)f0->Get("data"); //Before
    TTree* tree1 = (TTree*)f1->Get("data"); //After
    
    GetPedestal(tree0,&PedestalCh1[0],&PedestalCh2[0]);
    GetPedestal(tree1,&PedestalCh1[1],&PedestalCh2[1]);
    
    Pedestal[k][0]=(PedestalCh1[0]+PedestalCh1[1])/2;
    Pedestal[k][1]=(PedestalCh2[0]+PedestalCh2[1])/2;
    
    k++;
  }// chiudo for

  
  for(int i=0;i<(int)FileListPhysics.size();i++){
    std::cout <<(int)FileListPhysics.size() <<"   "<< (int)FileListPedestal.size()<< Pedestal[i][0] << "    " << Pedestal[i][1] << std::endl;
  }

  

  TH1D* HistoCh59[(int)FileListPhysics.size()];
  TH1D* HistoCh315[(int)FileListPhysics.size()];
    
  Double_t MeanTCh59[(int)FileListPhysics.size()], MeanTCh315[(int)FileListPhysics.size()];
  Double_t SigmaTCh59[(int)FileListPhysics.size()], SigmaTCh315[(int)FileListPhysics.size()];
  
  Double_t Peak1Ch59[(int)FileListPhysics.size()], Peak2Ch59[(int)FileListPhysics.size()];
  Double_t SigmaPeak1Ch59[(int)FileListPhysics.size()], SigmaPeak2Ch59[(int)FileListPhysics.size()];

  Double_t Peak1Ch315[(int)FileListPhysics.size()], Peak2Ch315[(int)FileListPhysics.size()];
  Double_t SigmaPeak1Ch315[(int)FileListPhysics.size()], SigmaPeak2Ch315[(int)FileListPhysics.size()];
  

  TFile* f0;
  TTree* tree0;

  TF1* FitSpectrum[(int)FileListPhysics.size()][2];
  
  gStyle->SetOptFit(0111);

  for(int i=0;i < (int)FileListPhysics.size();i++){
    
    f0= TFile::Open(("TestStability3/"+FileListPhysics.at(i)).c_str());
    tree0 = (TTree*)f0->Get("data");
    
    HistoCh59[i] = new TH1D(("HistoCh59N"+to_string(i)).c_str(),("HistoCh59N"+to_string(i)).c_str(),100,0,100);
    HistoCh315[i] = new TH1D(("HistoCh315N"+to_string(i)).c_str(),("HistoCh315N"+to_string(i)).c_str(),100,0,100);
    
    GetSpectrum(tree0,HistoCh59[i],HistoCh315[i],Pedestal[i][0],Pedestal[i][1]);
    GetMeanTemperature(tree0,&MeanTCh59[i],&SigmaTCh59[i],&MeanTCh315[i],&SigmaTCh315[i]);
    
    TCanvas* canvino = new TCanvas("Canvino","Canvino",1200,600);
    canvino->Divide(2,1);


    FitSpectrum[i][0]=FitNaSpectrum(HistoCh59[i]);
    FitSpectrum[i][1]=FitNaSpectrum(HistoCh315[i]);
    
    canvino->cd(1);
    HistoCh59[i]->GetXaxis()->SetTitle("E [DU]");
    HistoCh59[i]->GetYaxis()->SetTitle("Counts");
    HistoCh59[i]->Draw();
    FitSpectrum[i][0]->Draw("SAME");

    canvino->cd(2);
    HistoCh315[i]->GetXaxis()->SetTitle("E [DU]");
    HistoCh315[i]->GetYaxis()->SetTitle("Counts");
    HistoCh315[i]->Draw();
    FitSpectrum[i][1]->Draw("SAME");    

    canvino->SaveAs(("Plot/EnergyTempV2/Partials/Canvas"+to_string(i)+".png").c_str());
    delete canvino;
  }

  for(int i=0;i<(int)FileListPhysics.size();i++){
    std::cout << MeanTCh59[i] << " " <<SigmaTCh59[i]<< " " << MeanTCh315[i] << " " << SigmaTCh315[i] << std::endl;

    Peak1Ch59[i]=FitSpectrum[i][0]->GetParameter(1);
    SigmaPeak1Ch59[i]=FitSpectrum[i][0]->GetParError(1);
    Peak2Ch59[i]=FitSpectrum[i][0]->GetParameter(6);
    SigmaPeak2Ch59[i]=FitSpectrum[i][0]->GetParError(6);

    Peak1Ch315[i]=FitSpectrum[i][1]->GetParameter(1);
    SigmaPeak1Ch315[i]=FitSpectrum[i][1]->GetParError(1);
    Peak2Ch315[i]=FitSpectrum[i][1]->GetParameter(6);
    SigmaPeak2Ch315[i]=FitSpectrum[i][1]->GetParError(6);
    
    std::cout << Peak1Ch59[i]<< " " << SigmaPeak1Ch59[i]<< " " << Peak2Ch59[i]<< " " << SigmaPeak2Ch59[i]<< " "<< Peak1Ch315[i]<< " " <<SigmaPeak1Ch315[i] <<" "<< Peak2Ch315[i]<< " " <<SigmaPeak2Ch315[i]<< std::endl; 
    
  }
  
  TCanvas* PlotEVsT = new TCanvas("PlotEVsT","PlotEVsT",1200,600);

  TGraphErrors* Graph1Ch59 = new TGraphErrors((int)FileListPhysics.size(),MeanTCh59,Peak1Ch59,SigmaTCh59,SigmaPeak1Ch59);
  TGraphErrors* Graph2Ch59 = new TGraphErrors((int)FileListPhysics.size(),MeanTCh59,Peak2Ch59,SigmaTCh59,SigmaPeak2Ch59);
  
  TGraphErrors* Graph1Ch315 = new TGraphErrors((int)FileListPhysics.size(),MeanTCh315,Peak1Ch315,SigmaTCh315,SigmaPeak1Ch315);
  TGraphErrors* Graph2Ch315 = new TGraphErrors((int)FileListPhysics.size(),MeanTCh315,Peak2Ch315,SigmaTCh315,SigmaPeak2Ch315);

  PlotEVsT->Divide(2,1);
  PlotEVsT->cd(1);
  Graph1Ch59->SetMaximum(90);
  Graph1Ch59->SetMinimum(35);
  Graph1Ch59->GetXaxis()->SetTitle("TMeanRun [°C]");
  Graph1Ch59->GetYaxis()->SetTitle("E [DU]");
  Graph1Ch59->Draw("AP");
  Graph2Ch59->Draw("SAMEP");

  PlotEVsT->cd(2);
  Graph1Ch315->SetMaximum(90);
  Graph1Ch315->SetMinimum(35);
  Graph1Ch315->GetXaxis()->SetTitle("TMeanRun [°C]");
  Graph1Ch315->GetYaxis()->SetTitle("E [DU]");
  Graph1Ch315->Draw("AP");
  Graph2Ch315->Draw("SAMEP");


  PlotEVsT->SaveAs("Plot/EnergyTempV2/PlotEVsT.png");
}
