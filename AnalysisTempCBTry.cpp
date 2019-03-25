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
#include "TPad.h"
#include "TMath.h"

#include "ReadFileList.cc"
#include "FunctionsTempCB.cc"


void Analysis(){
  
  gROOT->Reset();
  gROOT->SetBatch(kTRUE);

  gSystem->Exec("mkdir Plot");
  
  gSystem->Exec("mkdir Plot/EnergyTempCB");
  gSystem->Exec("mkdir Plot/EnergyTempCB/Partials/");

  //gStyle->SetOptStat("000001000");

  vector<string> FileListPedestal;
  FileListPedestal=ReadData("TestStability3/PedFile.txt");

  vector<string> FileListPhysics;
  FileListPhysics=ReadData("TestStability3/PhysFile.txt");

  int NFilePhys=(int)FileListPhysics.size();
  
  //Get Pedestal
  Double_t Pedestal[NFilePhys][2];
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

  
  for(int i=0;i<NFilePhys;i++){
    std::cout <<NFilePhys <<"   "<< (int)FileListPedestal.size()<< Pedestal[i][0] << "    " << Pedestal[i][1] << std::endl;
  }


  TH1D* HistoCh59[NFilePhys];
  TH1D* HistoCh315[NFilePhys];
    
  Double_t MeanTCh59[NFilePhys], MeanTCh315[NFilePhys];
  Double_t SigmaTCh59[NFilePhys], SigmaTCh315[NFilePhys];
  Double_t MeanTGlobal[NFilePhys], SigmaTGlobal[NFilePhys];

  Double_t Peak1Ch59[NFilePhys], Peak2Ch59[NFilePhys];
  Double_t SigmaPeak1Ch59[NFilePhys], SigmaPeak2Ch59[NFilePhys];

  Double_t Peak1Ch315[NFilePhys], Peak2Ch315[NFilePhys];
  Double_t SigmaPeak1Ch315[NFilePhys], SigmaPeak2Ch315[NFilePhys];
  

  TFile* f0;
  TTree* tree0;

  TF1* FitSpectrum[NFilePhys][2];
  
  gStyle->SetOptFit(0111);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.15);

  for(int i=0;i < NFilePhys;i++){
    
    f0= TFile::Open(("TestStability3/"+FileListPhysics.at(i)).c_str());
    tree0 = (TTree*)f0->Get("data");
    
    HistoCh59[i] = new TH1D(("HistoCh59N"+to_string(i)).c_str(),("HistoCh59N"+to_string(i)).c_str(),100,0,100);
    HistoCh315[i] = new TH1D(("HistoCh315N"+to_string(i)).c_str(),("HistoCh315N"+to_string(i)).c_str(),100,0,100);
    
    GetSpectrum(tree0,HistoCh59[i],HistoCh315[i],Pedestal[i][0],Pedestal[i][1]);
    GetMeanTemperature(tree0,MeanTCh59,SigmaTCh59,MeanTCh315,SigmaTCh315,MeanTGlobal,SigmaTGlobal,i);
    
    TCanvas* canvino = new TCanvas("Canvino","Canvino",1200,600);
    canvino->Divide(2,1);


    FitSpectrum[i][0]=FitNaSpectrumCB(HistoCh59[i]);
    FitSpectrum[i][1]=FitNaSpectrumCB(HistoCh315[i]);
    
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

    canvino->SaveAs(("Plot/EnergyTempCB/Partials/Canvas"+to_string(i)+".png").c_str());
    delete canvino;
  }

  for(int i=0;i<NFilePhys;i++){
    std::cout << MeanTCh59[i] << " " <<SigmaTCh59[i]<< " " << MeanTCh315[i] << " " << SigmaTCh315[i] <<" " << MeanTGlobal[i]<< " " << SigmaTGlobal[i] <<std::endl;

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

  TGraphErrors* Graph1Ch59 = new TGraphErrors(NFilePhys,MeanTCh59,Peak1Ch59,SigmaTCh59,SigmaPeak1Ch59);
  TGraphErrors* Graph2Ch59 = new TGraphErrors(NFilePhys,MeanTCh59,Peak2Ch59,SigmaTCh59,SigmaPeak2Ch59);
  
  TGraphErrors* Graph1Ch315 = new TGraphErrors(NFilePhys,MeanTCh315,Peak1Ch315,SigmaTCh315,SigmaPeak1Ch315);
  TGraphErrors* Graph2Ch315 = new TGraphErrors(NFilePhys,MeanTCh315,Peak2Ch315,SigmaTCh315,SigmaPeak2Ch315);

  
  PlotEVsT->Divide(2,1);
  PlotEVsT->cd(1)->SetGridx();
  PlotEVsT->cd(1)->SetGridy();
  
  Graph1Ch59->SetMaximum(90);
  Graph1Ch59->SetMinimum(35);
  Graph1Ch59->SetTitle("EnergyVsTempCh59");
  Graph1Ch59->GetXaxis()->SetTitle("TMeanRun [°C]");
  Graph1Ch59->GetYaxis()->SetTitle("E [DU]");
  Graph1Ch59->Draw("AP");
  Graph2Ch59->Draw("SAMEP");
  

  PlotEVsT->cd(2)->SetGridx();
  PlotEVsT->cd(2)->SetGridy();
  
  Graph1Ch315->SetMaximum(90);
  Graph1Ch315->SetMinimum(35);
  Graph1Ch315->SetTitle("EnergyVsTempCh315");
  Graph1Ch315->GetXaxis()->SetTitle("TMeanRun [°C]");
  Graph1Ch315->GetYaxis()->SetTitle("E [DU]");
  Graph1Ch315->Draw("AP");
  Graph2Ch315->Draw("SAMEP");
  
  
  
  
  PlotEVsT->SaveAs("Plot/EnergyTempCB/PlotEVsT.png");

  //gROOT->SetBatch(kFALSE);
  //gStyle->SetOptFit(1111);
  
  TCanvas* CanvGlobalTemp = new TCanvas("CanvasGlobalTemp", "CanvasGlobalTemp",1200,600);
  
  
  TGraphErrors* Graph1GlobalTempCh59= new TGraphErrors(NFilePhys,MeanTGlobal,Peak1Ch59,SigmaTGlobal,SigmaPeak1Ch59);
  TGraphErrors* Graph2GlobalTempCh59= new TGraphErrors(NFilePhys,MeanTGlobal,Peak2Ch59,SigmaTGlobal,SigmaPeak2Ch59);
  TGraphErrors* Graph1GlobalTempCh315= new TGraphErrors(NFilePhys,MeanTGlobal,Peak1Ch315,SigmaTGlobal,SigmaPeak1Ch315);
  TGraphErrors* Graph2GlobalTempCh315= new TGraphErrors(NFilePhys,MeanTGlobal,Peak2Ch315,SigmaTGlobal,SigmaPeak2Ch315);
  ////////////////////////////////////////////////////////////////////////////////////
  Double_t RatioPeakCh59[NFilePhys], RatioPeakCh315[NFilePhys],SigmaRPeakCh59[NFilePhys],SigmaRPeakCh315[NFilePhys];
  
  RatioWithError(Peak1Ch59,Peak2Ch59,SigmaPeak1Ch59,SigmaPeak2Ch59,RatioPeakCh59,SigmaRPeakCh59,NFilePhys);
  RatioWithError(Peak1Ch315,Peak2Ch315,SigmaPeak1Ch315,SigmaPeak2Ch315,RatioPeakCh315,SigmaRPeakCh315,NFilePhys);
  
  TGraphErrors* GraphRatioCh59 = new TGraphErrors(NFilePhys,MeanTGlobal,RatioPeakCh59,SigmaTGlobal,SigmaRPeakCh59);
  TGraphErrors* GraphRatioCh315 = new TGraphErrors(NFilePhys,MeanTGlobal,RatioPeakCh315,SigmaTGlobal,SigmaRPeakCh315);

  /////////////////////////////////////////////////////////////////////////////////////
  //RatioEnergyCh1 (TGlobal)
  CanvGlobalTemp->Divide(2,1);
  CanvGlobalTemp->cd(1);

  TF1* fitRatioCh59 = new TF1("fitRatioCh59","[0]+[1]*x");
  TF1* fitRatioCh315 = new TF1("fitRatioCh315","[0]+[1]*x");
  
  TPad *p2Ch59 = new TPad("p2","p3",0.,0.,1.,0.3); p2Ch59->Draw();
  TPad *p1Ch59 = new TPad("p1","p1",0.,0.3,1.,1.); p1Ch59->Draw();
  p1Ch59->SetBottomMargin(0.001);
  p2Ch59->SetTopMargin(0.001);
  p2Ch59->SetBottomMargin(0.3);

  p1Ch59->cd()->SetGridx();
  p1Ch59->cd()->SetGridy();
  
  Graph1GlobalTempCh59->SetMaximum(90);
  Graph1GlobalTempCh59->SetMinimum(35);
  Graph1GlobalTempCh59->SetTitle("EnergyVsBoxTempCh59");
  Graph1GlobalTempCh59->GetYaxis()->SetTitle("E [DU]");
  Graph1GlobalTempCh59->Draw("AP");
  Graph2GlobalTempCh59->Draw("SAMEP");

  p2Ch59->cd()->SetGridx();
  p2Ch59->cd()->SetGridy();
  SetStyleRatioPlot(GraphRatioCh59,0.35,0.65);
  GraphRatioCh59->Draw("AP");
  GraphRatioCh59->Fit("fitRatioCh59");
  //////////////////////////////////////////////////////////////////////////////////////
  //RatioEnergyCh1 (TGlobal)
  CanvGlobalTemp->cd(2);
  
 
  TPad *p2Ch315 = new TPad("p2","p3",0.,0.,1.,0.3); p2Ch315->Draw();
  TPad *p1Ch315 = new TPad("p1","p1",0.,0.3,1.,1.); p1Ch315->Draw();
  p1Ch315->SetBottomMargin(0.001);
  p2Ch315->SetTopMargin(0.001);
  p2Ch315->SetBottomMargin(0.3);

  p1Ch315->cd()->SetGridx();
  p1Ch315->cd()->SetGridy();
  
  Graph1GlobalTempCh315->SetMaximum(90);
  Graph1GlobalTempCh315->SetMinimum(35);
  Graph1GlobalTempCh315->SetTitle("EnergyVsBoxTempCh315");
  Graph1GlobalTempCh315->GetYaxis()->SetTitle("E [DU]");
  Graph1GlobalTempCh315->Draw("AP");
  Graph2GlobalTempCh315->Draw("SAMEP");

  p2Ch315->cd()->SetGridx();
  p2Ch315->cd()->SetGridy();
  SetStyleRatioPlot(GraphRatioCh315,0.35,0.65);
  GraphRatioCh315->Draw("AP");
  GraphRatioCh315->Fit("fitRatioCh315");
  ////////////////////////////////////////////////////////////////////////////////
  CanvGlobalTemp->SaveAs("Plot/EnergyTempCB/PlotEVsTglobal.png");
  ////////////////////////////////////////////////////////////////////////////////
  //RatioPeak1(511)PeackCh59/315 (TGlobal)

  Double_t ratioPeak1[NFilePhys], sigmaRatioPeak1[NFilePhys], ratioPeak2[NFilePhys], sigmaRatioPeak2[NFilePhys];
  
  RatioWithError(Peak1Ch59,Peak1Ch315,SigmaPeak1Ch59,SigmaPeak1Ch315,ratioPeak1,sigmaRatioPeak1,NFilePhys);
  RatioWithError(Peak2Ch59,Peak2Ch315,SigmaPeak2Ch59,SigmaPeak2Ch315,ratioPeak2,sigmaRatioPeak2,NFilePhys);

  TGraphErrors* GraphRatioPeak1 = new TGraphErrors(NFilePhys,MeanTGlobal,ratioPeak1,SigmaTGlobal,sigmaRatioPeak1);
  TGraphErrors* GraphRatioPeak2 = new TGraphErrors(NFilePhys,MeanTGlobal,ratioPeak2,SigmaTGlobal,sigmaRatioPeak2);
  
    
  TCanvas* CanvasComparisonPeak= new TCanvas("CanvasComparisonPeak","CanvasComparisonPeak",1200,600);
  CanvasComparisonPeak->Divide(2,1);
  CanvasComparisonPeak->cd(1);

  TF1* fitRatioPeak1 = new TF1("fitRatioPeak1","[0]+[1]*x");
  TF1* fitRatioPeak2 = new TF1("fitRatioPeak2","[0]+[1]*x");
  
  TPad *pad2Peak1 = new TPad("p2","p3",0.,0.,1.,0.3); pad2Peak1->Draw();
  TPad *pad1Peak1 = new TPad("p1","p1",0.,0.3,1.,1.); pad1Peak1->Draw();
  pad1Peak1->SetBottomMargin(0.001);
  pad2Peak1->SetTopMargin(0.001);
  pad2Peak1->SetBottomMargin(0.3);

  pad1Peak1->cd()->SetGridx();
  pad1Peak1->cd()->SetGridy();
  
  Graph1GlobalTempCh59->SetMaximum(45);
  Graph1GlobalTempCh59->SetMinimum(35);
  Graph1GlobalTempCh315->SetMarkerStyle(4);
  Graph1GlobalTempCh59->SetMarkerStyle(8);
  //Graph1GlobalTempCh315->SetMarkerSize(.7);
  //Graph1GlobalTempCh59->SetMarkerSize(.7);
  Graph1GlobalTempCh59->SetTitle("EnergyVsBoxTempPeak511Kev");
  Graph1GlobalTempCh59->GetYaxis()->SetTitle("E [DU]");
  Graph1GlobalTempCh59->Draw("AP");
  Graph1GlobalTempCh315->Draw("SAMEP");

  pad2Peak1->cd()->SetGridx();
  pad2Peak1->cd()->SetGridy();
  SetStyleRatioPlot(GraphRatioPeak1,0.9,1.1);
  GraphRatioPeak1->Draw("AP");
  GraphRatioPeak1->Fit("fitRatioPeak1");
  
  //////////////////////////////////////////////////////////////
  CanvasComparisonPeak->cd(2);
  TPad *pad2Peak2 = new TPad("p2","p3",0.,0.,1.,0.3); pad2Peak2->Draw();
  TPad *pad1Peak2 = new TPad("p1","p1",0.,0.3,1.,1.); pad1Peak2->Draw();
  pad1Peak2->SetBottomMargin(0.001);
  pad2Peak2->SetTopMargin(0.001);
  pad2Peak2->SetBottomMargin(0.3);

  pad1Peak2->cd()->SetGridx();
  pad1Peak2->cd()->SetGridy();
  
  Graph2GlobalTempCh315->SetMaximum(90);
  Graph2GlobalTempCh315->SetMinimum(75);
  Graph2GlobalTempCh315->SetMarkerStyle(4);
  Graph2GlobalTempCh59->SetMarkerStyle(8);
  //Graph2GlobalTempCh315->SetMarkerSize(.7);
  //Graph2GlobalTempCh59->SetMarkerSize(.7);
  Graph2GlobalTempCh315->SetTitle("EnergyVsBoxTempPeak1275KeV");
  Graph2GlobalTempCh315->GetYaxis()->SetTitle("E [DU]");
  Graph2GlobalTempCh315->Draw("AP");
  Graph2GlobalTempCh59->Draw("SAMEP");

  pad2Peak2->cd()->SetGridx();
  pad2Peak2->cd()->SetGridy();
  SetStyleRatioPlot(GraphRatioPeak2,0.9,1.1);
  GraphRatioPeak2->Draw("AP");
  GraphRatioPeak2->Fit("fitRatioPeak2");
  
  CanvasComparisonPeak->SaveAs("Plot/EnergyTempCB/PeakComparisonVsTemp.png");

  ////////////////////////////////////////////////////////////////////////////////
  TH1D* HistoSumCh59 = new TH1D("HistoSumCh59","HistoSumCh59",100,0,100);
  TH1D* HistoSumCh315 = new TH1D("HistoSumCh315","HistoSumCh315",100,0,100);
  
  for(int i=0;i<NFilePhys;i++){
    
    HistoSumCh59->Add(HistoCh59[i]);
    HistoSumCh315->Add(HistoCh315[i]);

  }

  TF1* fitHistoSum[2];
  
  TCanvas* canvasSum = new TCanvas("CanvasSumHisto","CanvasSumHisto",1200,600);
  canvasSum->Divide(2,1);

  canvasSum->cd(1);
  HistoSumCh59->Draw();
  fitHistoSum[0]=FitNaSpectrumCB(HistoSumCh59);
  fitHistoSum[0]->Draw("SAME");

  canvasSum->cd(2);
  HistoSumCh315->Draw();
  fitHistoSum[1]=FitNaSpectrumCB(HistoSumCh315);
  fitHistoSum[1]->Draw("SAME");

  canvasSum->SaveAs("Plot/EnergyTempCB/HistoSum.png");

  ///////////////////////////////////////////////////////////////////////////////
  Double_t MeanCh59P1[NFilePhys],MeanCh315P1[NFilePhys],SigmaCh59P1[NFilePhys],SigmaCh315P1[NFilePhys];
  Double_t ErrMeanCh59P1[NFilePhys],ErrMeanCh315P1[NFilePhys],ErrSigmaCh59P1[NFilePhys],ErrSigmaCh315P1[NFilePhys];

  Double_t MeanCh59P2[NFilePhys],MeanCh315P2[NFilePhys],SigmaCh59P2[NFilePhys],SigmaCh315P2[NFilePhys];
  Double_t ErrMeanCh59P2[NFilePhys],ErrMeanCh315P2[NFilePhys],ErrSigmaCh59P2[NFilePhys],ErrSigmaCh315P2[NFilePhys];

  Double_t ResCh59P1[NFilePhys],ResCh315P1[NFilePhys],SigmaResCh59P1[NFilePhys],SigmaResCh315P1[NFilePhys];
  Double_t ResCh59P2[NFilePhys],ResCh315P2[NFilePhys],SigmaResCh59P2[NFilePhys],SigmaResCh315P2[NFilePhys];
    
  for(int i=0;i<NFilePhys;i++){

    MeanCh59P1[i] = FitSpectrum[i][0]->GetParameter(1);
    MeanCh315P1[i] = FitSpectrum[i][1]->GetParameter(1);
    SigmaCh59P1[i] = FitSpectrum[i][0]->GetParameter(2);
    SigmaCh315P1[i] = FitSpectrum[i][1]->GetParameter(2);

    ErrMeanCh59P1[i] = FitSpectrum[i][0]->GetParError(1);
    ErrMeanCh315P1[i] = FitSpectrum[i][1]->GetParError(1);
    ErrSigmaCh59P1[i] = FitSpectrum[i][0]->GetParError(2);
    ErrSigmaCh315P1[i] = FitSpectrum[i][1]->GetParError(2);
    
    MeanCh59P2[i] = FitSpectrum[i][0]->GetParameter(6);
    MeanCh315P2[i] = FitSpectrum[i][1]->GetParameter(6);
    SigmaCh59P2[i] = FitSpectrum[i][0]->GetParameter(7);
    SigmaCh315P2[i] = FitSpectrum[i][1]->GetParameter(7);

    ErrMeanCh59P2[i] = FitSpectrum[i][0]->GetParError(6);
    ErrMeanCh315P2[i] = FitSpectrum[i][1]->GetParError(6);
    ErrSigmaCh59P2[i] = FitSpectrum[i][0]->GetParError(7);
    ErrSigmaCh315P2[i] = FitSpectrum[i][1]->GetParError(7);

    std::cout<< MeanCh59P2[i] << " " << SigmaCh59P2[i] << std::endl;
  }

  RatioWithError(SigmaCh59P1,MeanCh59P1,ErrSigmaCh59P1,ErrMeanCh59P1,ResCh59P1,SigmaResCh59P1,NFilePhys);
  RatioWithError(SigmaCh315P1,MeanCh315P1,ErrSigmaCh315P1,ErrMeanCh315P1,ResCh315P1,SigmaResCh315P1,NFilePhys);
  RatioWithError(SigmaCh59P2,MeanCh59P2,ErrSigmaCh59P2,ErrMeanCh59P2,ResCh59P2,SigmaResCh59P2,NFilePhys);
  RatioWithError(SigmaCh315P2,MeanCh315P2,ErrSigmaCh315P2,ErrMeanCh315P2,ResCh315P2,SigmaResCh315P2,NFilePhys);

  //////////////////////////////////////////////////////////////////////////////////

  TCanvas* ResolutionVsTemp= new TCanvas("ResolutionVsTemp","ResolutionVsTem",1200,600);
  ResolutionVsTemp->Divide(2,2);
  
  TGraphErrors* PlotResCh59P1 = new TGraphErrors(NFilePhys,MeanTGlobal,ResCh59P1,SigmaTGlobal,SigmaResCh59P1);
  TGraphErrors* PlotResCh315P1 = new TGraphErrors(NFilePhys,MeanTGlobal,ResCh315P1,SigmaTGlobal,SigmaResCh315P1);
  TGraphErrors* PlotResCh59P2 = new TGraphErrors(NFilePhys,MeanTGlobal,ResCh59P2,SigmaTGlobal,SigmaResCh59P2);
  TGraphErrors* PlotResCh315P2 = new TGraphErrors(NFilePhys,MeanTGlobal,ResCh315P2,SigmaTGlobal,SigmaResCh315P2);

  PlotResCh59P1->SetTitle("ResCh59P 511KeV");
  PlotResCh315P1->SetTitle("ResCh315P 511KeV");

  PlotResCh59P2->SetTitle("ResCh59P 1275KeV");
  PlotResCh315P2->SetTitle("ResCh315P 1275KeV");

  PlotResCh59P1->GetXaxis()->SetTitle("TMeanBox [°C]");
  PlotResCh59P1->GetYaxis()->SetTitle("EnergyResolution");
  //PlotResCh59P1->GetYaxis()->SetLimits(0.5,0.9);
  //PlotResCh59P1->GetYaxis()->SetRangeUser(0.5,0.9);
  PlotResCh315P1->GetXaxis()->SetTitle("TMeanBox [°C]");
  PlotResCh315P1->GetYaxis()->SetTitle("EnergyResolution");
  //PlotResCh315P1->GetYaxis()->SetLimits(0.5,0.9);
  //PlotResCh315P1->GetYaxis()->SetRangeUser(0.5,0.9);
  PlotResCh59P2->GetXaxis()->SetTitle("TMeanBox [°C]");
  PlotResCh59P2->GetYaxis()->SetTitle("EnergyResolution");
  //PlotResCh59P2->GetYaxis()->SetLimits(0.5,0.9);
  //PlotResCh59P2->GetYaxis()->SetRangeUser(0.5,0.9);
  PlotResCh315P2->GetXaxis()->SetTitle("TMeanBox [°C]");
  PlotResCh315P2->GetYaxis()->SetTitle("EnergyResolution");
  //PlotResCh315P2->GetYaxis()->SetLimits(0.5,0.9);
  //PlotResCh315P2->GetYaxis()->SetRangeUser(0.5,0.9);
  
  ResolutionVsTemp->cd(1);
  PlotResCh59P1->Draw("AP");
  ResolutionVsTemp->cd(2);
  PlotResCh315P1->Draw("AP");
  ResolutionVsTemp->cd(3);
  PlotResCh59P2->Draw("AP");
  ResolutionVsTemp->cd(4);
  PlotResCh315P2->Draw("AP");

  ResolutionVsTemp->SaveAs("Plot/EnergyTempCB/ResolutionVsTemp.png");
  
}
