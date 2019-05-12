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
#include "TLine.h"

void GetPedestal(TTree* tree,Double_t* pedch1,Double_t* pedch2){
  
  TH1F* histoCh1 = new TH1F("pippoch1","pippoch1",600,0,600);
  TH1F* histoCh2 = new TH1F("pippoch2","pippoch2",600,0,600);

  Float_t energy;
  UShort_t chID;

  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("channelID",&chID);
   
  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(chID==59) { histoCh1->Fill(energy); }//chiudo if
    else if(chID==288){ histoCh2->Fill(energy); }//chiudo else
  }//chiudo for
  
  std::cout << "Filled"<< std::endl;

  *pedch1 = histoCh1->GetMean();
  *pedch2 = histoCh2->GetMean();
  
  delete histoCh1;
  delete histoCh2;
}

TF1* CalibrationCurve(TGraphErrors* graph,int i){
  
  TF1* fit = new TF1(("fitCalib"+std::to_string(i)).c_str()," [0] * (1-exp(-[1]*x) )",0,1300);

  fit->SetParameter(0,144);
  fit->SetParameter(1,4e-4);

  graph->Fit(fit,"Q");

  return fit;
}


void GetSpectrum(TTree* tree, TH1D* histoCh1, TH1D* histoCh2, Double_t MeanPedCh59, Double_t MeanPedCh288){
  
  Float_t energy;
  UShort_t chID;

  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("channelID",&chID);
   
  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(chID==59) { histoCh1->Fill(energy-MeanPedCh59); }//chiudo if
    else if(chID==288){ histoCh2->Fill(energy-MeanPedCh288); }//chiudo else
  }//chiudo for
  
  std::cout << "Filled"<< std::endl;
  
  
}


void GetMeanTemperature(TTree* tree, Double_t* MeanTempCh59, Double_t* SigmaTempCh59,Double_t* MeanTempCh288, Double_t* SigmaTempCh288, Double_t* MeanTGlobal,Double_t* SigmaTGlobal,int j){

  UShort_t chID;
  Double_t temp1,temp2,temp3;

  TH1D* histo1 = new TH1D("pippoch1","pippoch1",20,10,35);
  TH1D* histoCh59 = new TH1D("pippoch59","pippoch59",20,20,40);
  TH1D* histoCh288 = new TH1D("pippoch288","pippoch288",20,20,40);

  tree->SetBranchAddress("temp1",&temp1);
  tree->SetBranchAddress("temp2",&temp2);
  tree->SetBranchAddress("temp3",&temp3);
  tree->SetBranchAddress("channelID",&chID);
   
  for(int i=0; i<tree->GetEntries(); i++){

    tree->GetEntry(i);
    histo1->Fill(temp1);
    if(chID==59) { histoCh59->Fill(temp2); }//chiudo if
    else if(chID==288){ histoCh288->Fill(temp3); }//chiudo else

  }//chiudo for
  
  MeanTGlobal[j]= histo1->GetMean();
  SigmaTGlobal[j]= histo1->GetMeanError();
  MeanTempCh59[j]= histoCh59->GetMean();
  SigmaTempCh59[j]= histoCh59->GetMeanError();
  SigmaTempCh288[j]= histoCh288->GetMeanError();
  MeanTempCh288[j]= histoCh288->GetMean();

  std::cout << "TMeanDone"<< std::endl;

}



TF1* FitNaSpectrumCB(TH1D* Profile){

  
  Double_t min;
  Double_t peak1,peak2;
  Profile->GetXaxis()->SetRangeUser(20,40);
  min = Profile->GetBinCenter(Profile->GetMinimumBin());
  Profile->GetXaxis()->UnZoom();

  Profile->GetXaxis()->SetRangeUser(30,55);
  peak1 = Profile->GetBinCenter(Profile->GetMaximumBin());
  Profile->GetXaxis()->UnZoom();

  Profile->GetXaxis()->SetRangeUser(75,98);
  peak2 = Profile->GetBinCenter(Profile->GetMaximumBin());
  Profile->GetXaxis()->UnZoom();

  std::cout << "p1 , p2  ->  " << peak1 << " " << peak2 << std::endl;

  TF1* spectrum = new TF1(Form("SpectrumFit_%s", Profile->GetName()),"[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) + [3] / (exp( (x*[4]-(2*[1]*[1]/([1]+2*[1])))) + 1)+ [5] * exp(-( x-[6] )*( x-[6] )/( 2* [7]* [7])) +crystalball([8],[9],[10],[11],[12])",min,92);

  Double_t max;
  max= Profile->GetMaximum();
  

  spectrum->SetParameter(0,max);
  spectrum->SetParameter(1,peak1);
  spectrum->SetParameter(2,3);
  spectrum->SetParameter(3,max/5);
  spectrum->SetParameter(4,0.82);
  spectrum->SetParameter(5,max/20);
  spectrum->SetParameter(6,peak2);
  spectrum->SetParameter(7,3.2);
  spectrum->SetParameter(8,max/12);
  spectrum->SetParameter(9,peak2/1.21);
  spectrum->SetParameter(10,4.3);
  spectrum->SetParameter(11,0.04);
  spectrum->SetParameter(12,-1.4);

  //std::cout << "kk6" << std::endl;
  spectrum->SetParLimits(10,3.8,6);
  spectrum->SetParLimits(3,4,2700);
  spectrum->SetParLimits(11,0.02,1);
  
  //spectrum->SetParLimits(5,800,2000);
  //spectrum->SetParLimits(8,1000,10000);
  
  //std::cout << "kk7" << std::endl;
  Profile->Fit(Form("SpectrumFit_%s", Profile->GetName()),"R0");

  //std::cout << "kk8" << std::endl;
  return spectrum;
}


TF1* FitNaSpectrumCBFull(TH1D* Profile){

  
  Double_t min;
  Double_t peak1,peak2, peak0;

  Profile->GetXaxis()->SetRangeUser(20,40);
  min = Profile->GetBinCenter(Profile->GetMinimumBin());
  Profile->GetXaxis()->UnZoom();

  Profile->GetXaxis()->SetRangeUser(30,55);
  peak1 = Profile->GetBinCenter(Profile->GetMaximumBin());
  Profile->GetXaxis()->UnZoom();

  Profile->GetXaxis()->SetRangeUser(75,98);
  peak2 = Profile->GetBinCenter(Profile->GetMaximumBin());
  Profile->GetXaxis()->UnZoom();

  Profile->GetXaxis()->SetRangeUser(0,30);
  peak0 = Profile->GetBinCenter(Profile->GetMaximumBin());
  Profile->GetXaxis()->UnZoom();

  std::cout << "p1 , p2  ->  " << peak1 << " " << peak2 << std::endl;

  TF1* spectrum = new TF1(Form("SpectrumFit_%s", Profile->GetName()),"[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) + [3] / (exp( (x*[4]-(2*[1]*[1]/([1]+2*[1])))) + 1)+ [5] * exp(-( x-[6] )*( x-[6] )/( 2* [7]* [7])) +crystalball([8],[9],[10],[11],[12])+crystalball([13],[14],[15],[16],[17])",10,92);

  Double_t max;
  max= Profile->GetMaximum();
  
  spectrum->SetParameter(0,max);
  spectrum->SetParameter(1,peak1);
  spectrum->SetParameter(2,3);
  spectrum->SetParameter(3,max/5);
  spectrum->SetParameter(4,0.82);
  spectrum->SetParameter(5,max/20);
  spectrum->SetParameter(6,peak2);
  spectrum->SetParameter(7,3.2);

  spectrum->SetParameter(8,max/12);
  spectrum->SetParameter(9,peak2/1.21);
  spectrum->SetParameter(10,4.3);
  spectrum->SetParameter(11,0.04);
  spectrum->SetParameter(12,-1.4);

  spectrum->SetParameter(13,4500);
  spectrum->SetParameter(14,peak0);
  spectrum->SetParameter(15,9);
  spectrum->SetParameter(16,1);
  spectrum->SetParameter(17,0.05);
  

 
  spectrum->SetParLimits(10,3.8,6);
  spectrum->SetParLimits(3,4,2700);
  spectrum->SetParLimits(11,0.02,1);
  spectrum->SetParLimits(17,0.002,3);
  
 
  Profile->Fit(Form("SpectrumFit_%s", Profile->GetName()),"R0");
 
  return spectrum;
}


void RatioWithError(Double_t* A,Double_t* B,Double_t* sA,Double_t* sB,Double_t* ratio, Double_t* sigmaRatio,int Ndata){

  for(int i=0;i<Ndata;i++){
    ratio[i]=A[i]/B[i];
    sigmaRatio[i]=TMath::Sqrt((sA[i]*sA[i])/(B[i]*B[i])+(A[i]*A[i])/(B[i]*B[i]*B[i]*B[i])*(sB[i]*sB[i]));
  }
  
}

void SetStyleRatioPlot(TGraphErrors* ratioPlot,Double_t minRange,Double_t maxRange){

  ratioPlot->SetTitle("");
  ratioPlot->GetXaxis()->SetTitle("TMeanBoxRun [°C]");
  ratioPlot->GetYaxis()->SetTitle("PeakRatio");
  ratioPlot->GetXaxis()->SetTitleSize(0.08);
  ratioPlot->GetYaxis()->SetTitleSize(0.08);
  ratioPlot->GetYaxis()->SetTitleOffset(0.5);
  ratioPlot->GetXaxis()->SetLabelSize(0.075);
  ratioPlot->GetYaxis()->SetLabelSize(0.075);
  ratioPlot->GetYaxis()->SetLimits(minRange,maxRange);
  ratioPlot->GetYaxis()->SetRangeUser(minRange,maxRange);
  
}



std::vector<std::string> ReadData(std::string FileName){
 
  std::vector<std::string> FileList;
  std::ifstream file(FileName);
  std::string str; 
  
  if(!file.is_open()){std::cout << "file non trovato" << std::endl;}

  while (std::getline(file, str))
    {
      FileList.push_back(str);
    }
  //  if(debug){
  //  for(int i=0; i< (int)FileList.size();i++){
  //   std::cout << FileList.at(i) << std::endl;
  // }
  //}
  return FileList;
}


int main(int argc, char* argv[] ){

  std::cout << "Analisys" << std::endl;
  
  std::string DirData(argv[1]);
  std::string RootFileName=DirData;
  RootFileName.erase(0,2);
  std::string OV(argv[2]);

  gSystem->Exec(("ls "+DirData+"/*PED*"+OV+"_singles.root > "+DirData+"/PedFile.txt").c_str());
  gSystem->Exec(("ls "+DirData+"/*PHYS*"+OV+"_singles.root > "+DirData+"/PhysFile.txt").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot").c_str());
  
  gSystem->Exec(("mkdir "+DirData+"/Plot/EnergyTempCB").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/EnergyTempCB/Partials/").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/EnergyTempCB/CalibPlot/").c_str());
  gSystem->Exec(("mkdir "+DirData+"/../RootFileGraph").c_str());
  
  
  //gStyle->SetOptStat("000001000");

  TFile* f = new TFile(("../RootFileGraph/"+RootFileName+".root").c_str(),"RECREATE");

  if(f) std::cout << "Root File Graph opened->"<<RootFileName << std::endl;

  std::vector<std::string> FileListPedestal;
  std::string ListFilePed = DirData+"/PedFile.txt";
  std::cout << "Lista File Pedestal: "<< ListFilePed << std::endl;    
  FileListPedestal=ReadData(ListFilePed);
  
  std::vector<std::string> FileListPhysics;
  std::string ListFilePhys = DirData+"/PhysFile.txt";
  std::cout << "Lista File Pedestal: "<< ListFilePed << std::endl;
  FileListPhysics=ReadData(ListFilePhys);

  int NFilePhys=(int)FileListPhysics.size();
  
  //Get Pedestal
  Double_t Pedestal[NFilePhys][2];
  Double_t PedestalCh1[2];
  Double_t PedestalCh2[2];

  int k=0;

  for(int i=0;i < (int)FileListPedestal.size()-1;i+=2){
    
    std::cout<< "open file:  " << (DirData+"/"+FileListPedestal.at(i)).c_str() << std::endl;
    TFile* f0= TFile::Open((DirData+"/"+FileListPedestal.at(i)).c_str());
    TFile* f1= TFile::Open((DirData+"/"+FileListPedestal.at(i)).c_str());

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
  TH1D* HistoCh288[NFilePhys];
    
  Double_t MeanTCh59[NFilePhys], MeanTCh288[NFilePhys];
  Double_t SigmaTCh59[NFilePhys], SigmaTCh288[NFilePhys];
  Double_t MeanTGlobal[NFilePhys], SigmaTGlobal[NFilePhys];

  Double_t Peak1Ch59[NFilePhys], Peak2Ch59[NFilePhys];
  Double_t SigmaPeak1Ch59[NFilePhys], SigmaPeak2Ch59[NFilePhys];

  Double_t Peak1Ch288[NFilePhys], Peak2Ch288[NFilePhys];
  Double_t SigmaPeak1Ch288[NFilePhys], SigmaPeak2Ch288[NFilePhys];
  

  TFile* f0;
  TTree* tree0;

  TF1* FitSpectrum[NFilePhys][2];
  
  gStyle->SetOptFit(1111);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.15);

  for(int i=0;i < NFilePhys;i++){
    
    std::cout << "i: " << i << std::endl;
    f0= TFile::Open((DirData+"/"+FileListPhysics.at(i)).c_str());
    tree0 = (TTree*)f0->Get("data");

        
    HistoCh59[i]  = new TH1D(Form("HistoCh59N%d", i),Form("HistoCh59N%d", i), 100,0,100);
    HistoCh288[i] = new TH1D(Form("HistoCh288N%d",i),Form("HistoCh288N%d",i), 100,0,100);

        
    GetSpectrum(tree0,HistoCh59[i],HistoCh288[i],Pedestal[i][0],Pedestal[i][1]);
    GetMeanTemperature(tree0,MeanTCh59,SigmaTCh59,MeanTCh288,SigmaTCh288,MeanTGlobal,SigmaTGlobal,i);

    //HistoCh59[i]->SetTitle(("HistoCh59Temp"+std::to_string(MeanTGlobal[i])).c_str());
    //HistoCh288[i]->SetTitle(("HistoCh288Temp"+std::to_string(MeanTGlobal[i])).c_str());
        
    TCanvas* canvino = new TCanvas("Canvino","Canvino",1200,600);
    canvino->Divide(2,1);
    
    //do{
    FitSpectrum[i][0]=FitNaSpectrumCB(HistoCh59[i]);
    //}while(isnan(FitSpectrum[i][0]->GetParError(5)));
   
    //do{
    FitSpectrum[i][1]=FitNaSpectrumCB(HistoCh288[i]);
    //}while(isnan(FitSpectrum[i][1]->GetParError(5)));
        
    canvino->cd(1);
    HistoCh59[i]->GetXaxis()->SetTitle("E [DU]");
    HistoCh59[i]->GetYaxis()->SetTitle("Counts");
    HistoCh59[i]->Draw();
    FitSpectrum[i][0]->Draw("SAME");

    canvino->cd(2);
    HistoCh288[i]->GetXaxis()->SetTitle("E [DU]");
    HistoCh288[i]->GetYaxis()->SetTitle("Counts");
    HistoCh288[i]->Draw();
    FitSpectrum[i][1]->Draw("SAME");    

    canvino->SaveAs((DirData+"/Plot/EnergyTempCB/Partials/Canvas"+std::to_string(i)+".png").c_str());
    delete canvino;
  }

  for(int i=0;i<NFilePhys;i++){
    std::cout << MeanTCh59[i] << " " <<SigmaTCh59[i]<< " " << MeanTCh288[i] << " " << SigmaTCh288[i] <<" " << MeanTGlobal[i]<< " " << SigmaTGlobal[i] <<std::endl;

    Peak1Ch59[i]=FitSpectrum[i][0]->GetParameter(1);
    SigmaPeak1Ch59[i]=FitSpectrum[i][0]->GetParError(1);
    Peak2Ch59[i]=FitSpectrum[i][0]->GetParameter(6);
    SigmaPeak2Ch59[i]=FitSpectrum[i][0]->GetParError(6);

    Peak1Ch288[i]=FitSpectrum[i][1]->GetParameter(1);
    SigmaPeak1Ch288[i]=FitSpectrum[i][1]->GetParError(1);
    Peak2Ch288[i]=FitSpectrum[i][1]->GetParameter(6);
    SigmaPeak2Ch288[i]=FitSpectrum[i][1]->GetParError(6);
    
    std::cout << Peak1Ch59[i]<< " " << SigmaPeak1Ch59[i]<< " " << Peak2Ch59[i]<< " " << SigmaPeak2Ch59[i]<< " "<< Peak1Ch288[i]<< " " <<SigmaPeak1Ch288[i] <<" "<< Peak2Ch288[i]<< " " <<SigmaPeak2Ch288[i]<< std::endl; 
    
  }
  
  TCanvas* PlotEVsT = new TCanvas("PlotEVsT","PlotEVsT",1200,600);

  TGraphErrors* Graph1Ch59 = new TGraphErrors(NFilePhys,MeanTCh59,Peak1Ch59,SigmaTCh59,SigmaPeak1Ch59);
  TGraphErrors* Graph2Ch59 = new TGraphErrors(NFilePhys,MeanTCh59,Peak2Ch59,SigmaTCh59,SigmaPeak2Ch59);
  
  TGraphErrors* Graph1Ch288 = new TGraphErrors(NFilePhys,MeanTCh288,Peak1Ch288,SigmaTCh288,SigmaPeak1Ch288);
  TGraphErrors* Graph2Ch288 = new TGraphErrors(NFilePhys,MeanTCh288,Peak2Ch288,SigmaTCh288,SigmaPeak2Ch288);

  
  PlotEVsT->Divide(2,1);
  PlotEVsT->cd(1)->SetGridx();
  PlotEVsT->cd(1)->SetGridy();
  
  Graph1Ch59->SetMaximum(90);
  Graph1Ch59->SetMinimum(30);
  Graph1Ch59->SetTitle("EnergyVsTempCh59P1");
  Graph2Ch59->SetTitle("EnergyVsTempCh59P2");
  Graph1Ch59->GetXaxis()->SetTitle("TMeanRun [°C]");
  Graph1Ch59->GetYaxis()->SetTitle("E [DU]");
  Graph1Ch59->Draw("AP");
  Graph2Ch59->Draw("SAMEP");
  

  PlotEVsT->cd(2)->SetGridx();
  PlotEVsT->cd(2)->SetGridy();
  
  Graph1Ch288->SetMaximum(90);
  Graph1Ch288->SetMinimum(30);
  Graph1Ch288->SetTitle("EnergyVsTempCh288P1");
  Graph2Ch288->SetTitle("EnergyVsTempCh288P2");
  Graph1Ch288->GetXaxis()->SetTitle("TMeanRun [°C]");
  Graph1Ch288->GetYaxis()->SetTitle("E [DU]");
  Graph1Ch288->Draw("AP");
  Graph2Ch288->Draw("SAMEP");
  
  f->cd();
  Graph1Ch59->SetName(Graph1Ch59->GetTitle());
  Graph1Ch288->SetName(Graph1Ch288->GetTitle());
  Graph2Ch59->SetName(Graph2Ch59->GetTitle());
  Graph2Ch288->SetName(Graph2Ch288->GetTitle());

  Graph1Ch59->Write();
  Graph1Ch288->Write();
  Graph2Ch59->Write();
  Graph2Ch288->Write();

  PlotEVsT->SaveAs((DirData+"/Plot/EnergyTempCB/PlotEVsT.png").c_str());

  //gROOT->SetBatch(kFALSE);
  //gStyle->SetOptFit(1111);
  
  TCanvas* CanvGlobalTemp = new TCanvas("CanvasGlobalTemp", "CanvasGlobalTemp",1200,600);
  
  
  TGraphErrors* Graph1GlobalTempCh59= new TGraphErrors(NFilePhys,MeanTGlobal,Peak1Ch59,SigmaTGlobal,SigmaPeak1Ch59);
  TGraphErrors* Graph2GlobalTempCh59= new TGraphErrors(NFilePhys,MeanTGlobal,Peak2Ch59,SigmaTGlobal,SigmaPeak2Ch59);
  TGraphErrors* Graph1GlobalTempCh288= new TGraphErrors(NFilePhys,MeanTGlobal,Peak1Ch288,SigmaTGlobal,SigmaPeak1Ch288);
  TGraphErrors* Graph2GlobalTempCh288= new TGraphErrors(NFilePhys,MeanTGlobal,Peak2Ch288,SigmaTGlobal,SigmaPeak2Ch288);
  ////////////////////////////////////////////////////////////////////////////////////
  Double_t RatioPeakCh59[NFilePhys], RatioPeakCh288[NFilePhys],SigmaRPeakCh59[NFilePhys],SigmaRPeakCh288[NFilePhys];
  
  RatioWithError(Peak1Ch59,Peak2Ch59,SigmaPeak1Ch59,SigmaPeak2Ch59,RatioPeakCh59,SigmaRPeakCh59,NFilePhys);
  RatioWithError(Peak1Ch288,Peak2Ch288,SigmaPeak1Ch288,SigmaPeak2Ch288,RatioPeakCh288,SigmaRPeakCh288,NFilePhys);
  
  TGraphErrors* GraphRatioCh59 = new TGraphErrors(NFilePhys,MeanTGlobal,RatioPeakCh59,SigmaTGlobal,SigmaRPeakCh59);
  TGraphErrors* GraphRatioCh288 = new TGraphErrors(NFilePhys,MeanTGlobal,RatioPeakCh288,SigmaTGlobal,SigmaRPeakCh288);

  /////////////////////////////////////////////////////////////////////////////////////
  //RatioEnergyCh1 (TGlobal)
  CanvGlobalTemp->Divide(2,1);
  CanvGlobalTemp->cd(1);

  TF1* fitRatioCh59 = new TF1("fitRatioCh59","[0]+[1]*x");
  TF1* fitRatioCh288 = new TF1("fitRatioCh288","[0]+[1]*x");

  fitRatioCh59->SetParameter(0,0.5);
  fitRatioCh59->SetParameter(1,0.00001);
  fitRatioCh288->SetParameter(0,0.5);
  fitRatioCh288->SetParameter(1,0.00001);
  
  TPad *p2Ch59 = new TPad("p2","p3",0.,0.,1.,0.3); p2Ch59->Draw();
  TPad *p1Ch59 = new TPad("p1","p1",0.,0.3,1.,1.); p1Ch59->Draw();
  p1Ch59->SetBottomMargin(0.001);
  p2Ch59->SetTopMargin(0.001);
  p2Ch59->SetBottomMargin(0.3);

  p1Ch59->cd()->SetGridx();
  p1Ch59->cd()->SetGridy();
  
  Graph1GlobalTempCh59->SetMaximum(90);
  Graph1GlobalTempCh59->SetMinimum(30);
  Graph1GlobalTempCh59->SetTitle("EnergyVsBoxTempCh59P1");
  Graph2GlobalTempCh59->SetTitle("EnergyVsBoxTempCh59P2");
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
  
 
  TPad *p2Ch288 = new TPad("p2","p3",0.,0.,1.,0.3); p2Ch288->Draw();
  TPad *p1Ch288 = new TPad("p1","p1",0.,0.3,1.,1.); p1Ch288->Draw();
  p1Ch288->SetBottomMargin(0.001);
  p2Ch288->SetTopMargin(0.001);
  p2Ch288->SetBottomMargin(0.3);

  p1Ch288->cd()->SetGridx();
  p1Ch288->cd()->SetGridy();
  
  Graph1GlobalTempCh288->SetMaximum(90);
  Graph1GlobalTempCh288->SetMinimum(30);
  Graph1GlobalTempCh288->SetTitle("EnergyVsBoxTempCh288P1");
  Graph2GlobalTempCh288->SetTitle("EnergyVsBoxTempCh288P2");
  Graph1GlobalTempCh288->GetYaxis()->SetTitle("E [DU]");
  Graph1GlobalTempCh288->Draw("AP");
  Graph2GlobalTempCh288->Draw("SAMEP");

  p2Ch288->cd()->SetGridx();
  p2Ch288->cd()->SetGridy();
  SetStyleRatioPlot(GraphRatioCh288,0.35,0.65);
  GraphRatioCh288->Draw("AP");
  GraphRatioCh288->Fit("fitRatioCh288");
  ////////////////////////////////////////////////////////////////////////////////
  f->cd();

  Graph1GlobalTempCh59->SetName(Graph1GlobalTempCh59->GetTitle());
  Graph1GlobalTempCh288->SetName(Graph1GlobalTempCh288->GetTitle());
  Graph2GlobalTempCh59->SetName(Graph2GlobalTempCh59->GetTitle());
  Graph2GlobalTempCh288->SetName(Graph2GlobalTempCh288->GetTitle());

  Graph1GlobalTempCh59->Write();
  Graph2GlobalTempCh59->Write();
  Graph1GlobalTempCh288->Write();
  Graph2GlobalTempCh288->Write();

  CanvGlobalTemp->SaveAs((DirData+"/Plot/EnergyTempCB/PlotEVsTglobal.png").c_str());
  
  ////////////////////////////////////////////////////////////////////////////////
  //RatioPeak1(511)PeackCh59/288 (TGlobal)

  Double_t ratioPeak1[NFilePhys], sigmaRatioPeak1[NFilePhys], ratioPeak2[NFilePhys], sigmaRatioPeak2[NFilePhys];
  
  RatioWithError(Peak1Ch59,Peak1Ch288,SigmaPeak1Ch59,SigmaPeak1Ch288,ratioPeak1,sigmaRatioPeak1,NFilePhys);
  RatioWithError(Peak2Ch59,Peak2Ch288,SigmaPeak2Ch59,SigmaPeak2Ch288,ratioPeak2,sigmaRatioPeak2,NFilePhys);

  TGraphErrors* GraphRatioPeak1 = new TGraphErrors(NFilePhys,MeanTGlobal,ratioPeak1,SigmaTGlobal,sigmaRatioPeak1);
  TGraphErrors* GraphRatioPeak2 = new TGraphErrors(NFilePhys,MeanTGlobal,ratioPeak2,SigmaTGlobal,sigmaRatioPeak2);
  
  TF1* fitRatioPeak1 = new TF1("fitRatioPeak1","[0]+[1]*x");
  TF1* fitRatioPeak2 = new TF1("fitRatioPeak2","[0]+[1]*x");

  fitRatioPeak1->SetParameter(0,1);
  fitRatioPeak1->SetParameter(1,0.00001);

  fitRatioPeak2->SetParameter(0,1);
  fitRatioPeak2->SetParameter(1,0.00001);
  
  TCanvas* CanvasComparisonPeak= new TCanvas("CanvasComparisonPeak","CanvasComparisonPeak",1200,600);
  CanvasComparisonPeak->Divide(2,1);
  CanvasComparisonPeak->cd(1);
    
  TPad *pad2Peak1 = new TPad("p2","p3",0.,0.,1.,0.3); pad2Peak1->Draw();
  TPad *pad1Peak1 = new TPad("p1","p1",0.,0.3,1.,1.); pad1Peak1->Draw();
  pad1Peak1->SetBottomMargin(0.001);
  pad2Peak1->SetTopMargin(0.001);
  pad2Peak1->SetBottomMargin(0.3);

  pad1Peak1->cd()->SetGridx();
  pad1Peak1->cd()->SetGridy();
  
  Graph1GlobalTempCh59->SetMaximum(45);
  Graph1GlobalTempCh59->SetMinimum(35);
  Graph1GlobalTempCh288->SetMarkerStyle(4);
  Graph1GlobalTempCh59->SetMarkerStyle(8);
  //Graph1GlobalTempCh288->SetMarkerSize(.7);
  //Graph1GlobalTempCh59->SetMarkerSize(.7);
  Graph1GlobalTempCh59->SetTitle("EnergyVsBoxTempPeak511Kev");
  Graph1GlobalTempCh59->GetYaxis()->SetTitle("E [DU]");
  Graph1GlobalTempCh59->Draw("AP");
  Graph1GlobalTempCh288->Draw("SAMEP");

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
  
  Graph2GlobalTempCh288->SetMaximum(90);
  Graph2GlobalTempCh288->SetMinimum(75);
  Graph2GlobalTempCh288->SetMarkerStyle(4);
  Graph2GlobalTempCh59->SetMarkerStyle(8);
  //Graph2GlobalTempCh288->SetMarkerSize(.7);
  //Graph2GlobalTempCh59->SetMarkerSize(.7);
  Graph2GlobalTempCh288->SetTitle("EnergyVsBoxTempPeak1275KeV");
  Graph2GlobalTempCh288->GetYaxis()->SetTitle("E [DU]");
  Graph2GlobalTempCh288->Draw("AP");
  Graph2GlobalTempCh59->Draw("SAMEP");

  pad2Peak2->cd()->SetGridx();
  pad2Peak2->cd()->SetGridy();
  SetStyleRatioPlot(GraphRatioPeak2,0.9,1.1);
  GraphRatioPeak2->Draw("AP");
  GraphRatioPeak2->Fit("fitRatioPeak2");
    
  CanvasComparisonPeak->SaveAs((DirData+"/Plot/EnergyTempCB/PeakComparisonVsTemp.png").c_str());
  ////////////////////////////////////////////////////////////////////////////////
  TH1D* HistoSumCh59 = new TH1D("HistoSumCh59","HistoSumCh59",100,0,100);
  TH1D* HistoSumCh288 = new TH1D("HistoSumCh288","HistoSumCh288",100,0,100);
  
  for(int i=0;i<NFilePhys;i++){
    
    HistoSumCh59->Add(HistoCh59[i]);
    HistoSumCh288->Add(HistoCh288[i]);

  }

  TF1* fitHistoSum[2];
  
  TCanvas* canvasSum = new TCanvas("CanvasSumHisto","CanvasSumHisto",1200,600);
  canvasSum->Divide(2,1);

  canvasSum->cd(1);
  HistoSumCh59->Draw();
  fitHistoSum[0]=FitNaSpectrumCB(HistoSumCh59);
  fitHistoSum[0]->Draw("SAME");

  canvasSum->cd(2);
  HistoSumCh288->Draw();
  fitHistoSum[1]=FitNaSpectrumCB(HistoSumCh288);
  fitHistoSum[1]->Draw("SAME");

  canvasSum->SaveAs((DirData+"/Plot/EnergyTempCB/HistoSum.png").c_str());

  ///////////////////////////////////////////////////////////////////////////////
  Double_t MeanCh59P1[NFilePhys],MeanCh288P1[NFilePhys],SigmaCh59P1[NFilePhys],SigmaCh288P1[NFilePhys];
  Double_t ErrMeanCh59P1[NFilePhys],ErrMeanCh288P1[NFilePhys],ErrSigmaCh59P1[NFilePhys],ErrSigmaCh288P1[NFilePhys];

  Double_t MeanCh59P2[NFilePhys],MeanCh288P2[NFilePhys],SigmaCh59P2[NFilePhys],SigmaCh288P2[NFilePhys];
  Double_t ErrMeanCh59P2[NFilePhys],ErrMeanCh288P2[NFilePhys],ErrSigmaCh59P2[NFilePhys],ErrSigmaCh288P2[NFilePhys];

  Double_t ResCh59P1[NFilePhys],ResCh288P1[NFilePhys],SigmaResCh59P1[NFilePhys],SigmaResCh288P1[NFilePhys];
  Double_t ResCh59P2[NFilePhys],ResCh288P2[NFilePhys],SigmaResCh59P2[NFilePhys],SigmaResCh288P2[NFilePhys];
    
  for(int i=0;i<NFilePhys;i++){

    MeanCh59P1[i] = FitSpectrum[i][0]->GetParameter(1);
    MeanCh288P1[i] = FitSpectrum[i][1]->GetParameter(1);
    SigmaCh59P1[i] = FitSpectrum[i][0]->GetParameter(2);
    SigmaCh288P1[i] = FitSpectrum[i][1]->GetParameter(2);

    ErrMeanCh59P1[i] = FitSpectrum[i][0]->GetParError(1);
    ErrMeanCh288P1[i] = FitSpectrum[i][1]->GetParError(1);
    ErrSigmaCh59P1[i] = FitSpectrum[i][0]->GetParError(2);
    ErrSigmaCh288P1[i] = FitSpectrum[i][1]->GetParError(2);
    
    MeanCh59P2[i] = FitSpectrum[i][0]->GetParameter(6);
    MeanCh288P2[i] = FitSpectrum[i][1]->GetParameter(6);
    SigmaCh59P2[i] = FitSpectrum[i][0]->GetParameter(7);
    SigmaCh288P2[i] = FitSpectrum[i][1]->GetParameter(7);

    ErrMeanCh59P2[i] = FitSpectrum[i][0]->GetParError(6);
    ErrMeanCh288P2[i] = FitSpectrum[i][1]->GetParError(6);
    ErrSigmaCh59P2[i] = FitSpectrum[i][0]->GetParError(7);
    ErrSigmaCh288P2[i] = FitSpectrum[i][1]->GetParError(7);

    //  std::cout<< MeanCh59P2[i] << " " << SigmaCh59P2[i] << std::endl;
  }

  RatioWithError(SigmaCh59P1,MeanCh59P1,ErrSigmaCh59P1,ErrMeanCh59P1,ResCh59P1,SigmaResCh59P1,NFilePhys);
  RatioWithError(SigmaCh288P1,MeanCh288P1,ErrSigmaCh288P1,ErrMeanCh288P1,ResCh288P1,SigmaResCh288P1,NFilePhys);
  RatioWithError(SigmaCh59P2,MeanCh59P2,ErrSigmaCh59P2,ErrMeanCh59P2,ResCh59P2,SigmaResCh59P2,NFilePhys);
  RatioWithError(SigmaCh288P2,MeanCh288P2,ErrSigmaCh288P2,ErrMeanCh288P2,ResCh288P2,SigmaResCh288P2,NFilePhys);

  //////////////////////////////////////////////////////////////////////////////////

  TCanvas* ResolutionVsTemp= new TCanvas("ResolutionVsTemp","ResolutionVsTem",1200,600);
  ResolutionVsTemp->Divide(2,2);
  
  TGraphErrors* PlotResCh59P1 = new TGraphErrors(NFilePhys,MeanTGlobal,ResCh59P1,SigmaTGlobal,SigmaResCh59P1);
  TGraphErrors* PlotResCh288P1 = new TGraphErrors(NFilePhys,MeanTGlobal,ResCh288P1,SigmaTGlobal,SigmaResCh288P1);
  TGraphErrors* PlotResCh59P2 = new TGraphErrors(NFilePhys,MeanTGlobal,ResCh59P2,SigmaTGlobal,SigmaResCh59P2);
  TGraphErrors* PlotResCh288P2 = new TGraphErrors(NFilePhys,MeanTGlobal,ResCh288P2,SigmaTGlobal,SigmaResCh288P2);

  Double_t YCh59P1Line= fitHistoSum[0]->GetParameter(2)/fitHistoSum[0]->GetParameter(1);
  TLine* ResSumCh59P1Line = new TLine(PlotResCh59P1->GetXaxis()->GetXmin(),YCh59P1Line,PlotResCh59P1->GetXaxis()->GetXmax(), YCh59P1Line);

  Double_t YCh288P1Line= fitHistoSum[1]->GetParameter(2)/fitHistoSum[1]->GetParameter(1);
  TLine* ResSumCh288P1Line = new TLine(PlotResCh288P1->GetXaxis()->GetXmin(),YCh288P1Line,PlotResCh59P1->GetXaxis()->GetXmax(), YCh288P1Line);

  Double_t YCh59P2Line= fitHistoSum[0]->GetParameter(7)/fitHistoSum[0]->GetParameter(6);
  TLine* ResSumCh59P2Line = new TLine(PlotResCh59P2->GetXaxis()->GetXmin(),YCh59P2Line,PlotResCh59P1->GetXaxis()->GetXmax(), YCh59P2Line);

  Double_t YCh288P2Line= fitHistoSum[1]->GetParameter(7)/fitHistoSum[1]->GetParameter(6);
  TLine* ResSumCh288P2Line = new TLine(PlotResCh288P2->GetXaxis()->GetXmin(),YCh288P2Line,PlotResCh59P1->GetXaxis()->GetXmax(), YCh288P2Line);

  PlotResCh59P1->SetTitle("ResCh59P 511KeV");
  PlotResCh288P1->SetTitle("ResCh288P 511KeV");

  PlotResCh59P2->SetTitle("ResCh59P 1275KeV");
  PlotResCh288P2->SetTitle("ResCh288P 1275KeV");

  PlotResCh59P1->GetXaxis()->SetTitle("TMeanBox [°C]");
  PlotResCh59P1->GetYaxis()->SetTitle("EnergyResolution");
  PlotResCh59P1->SetMinimum(0.06);
  PlotResCh59P1->SetMaximum(0.08);
  PlotResCh288P1->GetXaxis()->SetTitle("TMeanBox [°C]");
  PlotResCh288P1->GetYaxis()->SetTitle("EnergyResolution");
  PlotResCh288P1->SetMinimum(0.06);
  PlotResCh288P1->SetMaximum(0.08);
  PlotResCh59P2->GetXaxis()->SetTitle("TMeanBox [°C]");
  PlotResCh59P2->GetYaxis()->SetTitle("EnergyResolution");
  PlotResCh59P2->SetMinimum(0.03);
  PlotResCh59P2->SetMaximum(0.05);
  PlotResCh288P2->GetXaxis()->SetTitle("TMeanBox [°C]");
  PlotResCh288P2->GetYaxis()->SetTitle("EnergyResolution");
  PlotResCh288P2->SetMinimum(0.03);
  PlotResCh288P2->SetMaximum(0.05);
  
  ResolutionVsTemp->cd(1);
  PlotResCh59P1->Draw("AP");
  ResSumCh59P1Line->Draw("SAME");
  
  ResolutionVsTemp->cd(2);
  PlotResCh288P1->Draw("AP");
  ResSumCh288P1Line->Draw("SAME");

  ResolutionVsTemp->cd(3);
  PlotResCh59P2->Draw("AP");
  ResSumCh59P2Line->Draw("SAME");

  ResolutionVsTemp->cd(4);
  PlotResCh288P2->Draw("AP");
  ResSumCh288P2Line->Draw("SAME");

  ResolutionVsTemp->SaveAs((DirData+"/Plot/EnergyTempCB/ResolutionVsTemp.png").c_str());
  
    ///////////////////////////////////////////////////////////////////////////////////////////

  TGraphErrors* CalibPlot[NFilePhys][2];
  TCanvas* CanvasCalib;
  TF1* FitCalib[2];

  Double_t ACh59[NFilePhys],sACh59[NFilePhys],BCh59[NFilePhys],sBCh59[NFilePhys];
  Double_t ACh288[NFilePhys],sACh288[NFilePhys],BCh288[NFilePhys],sBCh288[NFilePhys];
  
  for(int i=0;i<NFilePhys;i++){

    CanvasCalib= new TCanvas("CanvasCalib","CanvasCalib",1200,600);
    CanvasCalib->Divide(2,1);

    CanvasCalib->cd(1);
    CalibPlot[i][0] = new TGraphErrors();
    CalibPlot[i][0]->SetTitle(("CalibPlotCh59"+std::to_string(i)).c_str());
    CalibPlot[i][0]->SetPoint(0,0,0);
    CalibPlot[i][0]->SetPoint(1,511,MeanCh59P1[i]);
    CalibPlot[i][0]->SetPointError(1,0,ErrMeanCh59P1[i]);
    CalibPlot[i][0]->SetPoint(2,1275,MeanCh59P2[i]);
    CalibPlot[i][0]->SetPointError(2,0,ErrMeanCh59P2[i]);
    FitCalib[0]= CalibrationCurve(CalibPlot[i][0],i);
    CalibPlot[i][0]->Draw("AP");
    FitCalib[0]->Draw("SAME");
    
    CanvasCalib->cd(2);
    CalibPlot[i][1] = new TGraphErrors();
    CalibPlot[i][1]->SetTitle(("CalibPlotCh288"+std::to_string(i)).c_str());
    CalibPlot[i][1]->SetPoint(0,0,0);
    CalibPlot[i][1]->SetPoint(1,511,MeanCh288P1[i]);
    CalibPlot[i][1]->SetPointError(1,0,ErrMeanCh288P1[i]);
    CalibPlot[i][1]->SetPoint(2,1275,MeanCh288P2[i]);
    CalibPlot[i][1]->SetPointError(2,0,ErrMeanCh288P2[i]);
    FitCalib[1]= CalibrationCurve(CalibPlot[i][1],i);
    CalibPlot[i][1]->Draw("AP");
    FitCalib[1]->Draw("SAME");
    
    CanvasCalib->SaveAs((DirData+"/Plot/EnergyTempCB/CalibPlot/CalibPlot"+std::to_string(i)+".png").c_str());

    ACh59[i]=FitCalib[0]->GetParameter(0);
    sACh59[i]=FitCalib[0]->GetParError(0);
    BCh59[i]=FitCalib[0]->GetParameter(1);
    sBCh59[i]=FitCalib[0]->GetParError(1);

    ACh288[i]=FitCalib[1]->GetParameter(0);
    sACh288[i]=FitCalib[1]->GetParError(0);
    BCh288[i]=FitCalib[1]->GetParameter(1);
    sBCh288[i]=FitCalib[1]->GetParError(1);
    
    delete CanvasCalib;
    delete FitCalib[0];
    delete FitCalib[1];
  }

  TGraphErrors* SatValVsMeanTemp[2];
  TGraphErrors* LYValVsMeanTemp[2];
  
  SatValVsMeanTemp[0]= new TGraphErrors(NFilePhys,MeanTGlobal,ACh59,SigmaTGlobal,sACh59);
  SatValVsMeanTemp[0]->SetTitle("SatValVsMeanTempCh59");
  SatValVsMeanTemp[0]->GetYaxis()->SetRangeUser(136,170);
  SatValVsMeanTemp[1]= new TGraphErrors(NFilePhys,MeanTGlobal,ACh288,SigmaTGlobal,sACh288);
  SatValVsMeanTemp[1]->SetTitle("SatValVsMeanTempCh288");
  SatValVsMeanTemp[1]->GetYaxis()->SetRangeUser(136,170);  

  LYValVsMeanTemp[0]= new TGraphErrors(NFilePhys,MeanTGlobal,BCh59,SigmaTGlobal,sBCh59);
  LYValVsMeanTemp[0]->SetTitle("LYValVsMeanTempCh59");
  LYValVsMeanTemp[0]->GetYaxis()->SetRangeUser(0.4e-3,0.7e-3);
  LYValVsMeanTemp[1]= new TGraphErrors(NFilePhys,MeanTGlobal,BCh288,SigmaTGlobal,sBCh288);
  LYValVsMeanTemp[1]->SetTitle("LYValVsMeanTempCh288");
  LYValVsMeanTemp[1]->GetYaxis()->SetRangeUser(0.4e-3,0.7e-3);
  
  TCanvas* CanvasValSat = new TCanvas("CanvasValSat","CanvasValSat",1200,600);
  CanvasValSat->Divide(2,1);
  CanvasValSat->cd(1);
  SatValVsMeanTemp[0]->Draw("AP");
  CanvasValSat->cd(2);
  SatValVsMeanTemp[1]->Draw("AP");
  CanvasValSat->SaveAs((DirData+"/Plot/EnergyTempCB/SatValueVsTemp.png").c_str());
  
  TCanvas* CanvasValLY = new TCanvas("CanvasValLY","CanvasValLY",1200,600);
  CanvasValLY->Divide(2,1);
  CanvasValLY->cd(1);
  LYValVsMeanTemp[0]->Draw("AP");
  CanvasValLY->cd(2);
  LYValVsMeanTemp[1]->Draw("AP");
  CanvasValLY->SaveAs((DirData+"/Plot/EnergyTempCB/LYValueVsTemp.png").c_str());

 
  f->Save();
  f->Close();
}
