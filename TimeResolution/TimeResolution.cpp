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

///////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////
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
    else if(chID==315){ histoCh2->Fill(energy); }//chiudo else
  }//chiudo for
  
  std::cout << "Filled"<< std::endl;

  *pedch1 = histoCh1->GetMean();
  *pedch2 = histoCh2->GetMean();
  
  delete histoCh1;
  delete histoCh2;
}



////////////////////////////////////////////////////////////////////////
void GetSpectrum(TTree* tree, TH1D* histoCh1, TH1D* histoCh2, Double_t MeanPedCh59, Double_t MeanPedCh315){
  
  Double_t energy[2];
  
  tree->SetBranchAddress("energy",energy);
     
  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
     histoCh1->Fill(energy[0]-MeanPedCh59); 
     histoCh2->Fill(energy[1]-MeanPedCh315);
  }//chiudo for
  
  std::cout << "Filled"<< std::endl;
  
  
}

////////////////////////////////////////////////////////////////////////////////////

void GetMeanTemperature(TTree* tree, Double_t* MeanTempCh59, Double_t* SigmaTempCh59,Double_t* MeanTempCh315, Double_t* SigmaTempCh315, Double_t* MeanTGlobal,Double_t* SigmaTGlobal,int j){

  Double_t temp1,temp2,temp3;

  TH1D* histo1 = new TH1D(Form("pippoch1%i",j),Form("pippoch1%i",j),20,10,35);
  TH1D* histoCh59 = new TH1D(Form("pippoch59%i",j),Form("pippoch59%i",j),20,20,40);
  TH1D* histoCh315 = new TH1D(Form("pippoch315%i",j),Form("pippoch315%i",j),20,20,40);

  tree->SetBranchAddress("temp1",&temp1);
  tree->SetBranchAddress("temp2",&temp2);
  tree->SetBranchAddress("temp3",&temp3);
  
  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    histo1->Fill(temp1);
    histoCh59->Fill(temp2);
    histoCh315->Fill(temp3); 
  }//chiudo for
  

  MeanTGlobal[j]= histo1->GetMean();
  SigmaTGlobal[j]= histo1->GetMeanError();
  MeanTempCh59[j]= histoCh59->GetMean();
  SigmaTempCh59[j]= histoCh59->GetMeanError();
  SigmaTempCh315[j]= histoCh315->GetMeanError();
  MeanTempCh315[j]= histoCh315->GetMean();

  std::cout << "TMeanDone"<< std::endl;

  delete histo1;
  delete histoCh59;
  delete histoCh315;
}


//////////////////////////////////////////////////////////////////////////////////////////////


TF1* FitNaSpectrumCB(TH1D* Profile){

  
  Double_t max;
  Double_t peak1;
  

  Profile->GetXaxis()->SetRangeUser(30,55);
  peak1 = Profile->GetBinCenter(Profile->GetMaximumBin());
  max = Profile->GetMaximum();
  Profile->GetXaxis()->UnZoom();


  TF1* spectrum = new TF1(Form("SpectrumFit_%s", Profile->GetName()),"[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) + [3] / (exp( (x*[4]-(2*[1]*[1]/([1]+2*[1])))) + 1)",12,60);


  spectrum->SetParameter(0,max);
  spectrum->SetParameter(1,peak1);
  spectrum->SetParameter(2,3);
  spectrum->SetParameter(3,max/5);
  spectrum->SetParameter(4,0.82);
  
  Profile->Fit(Form("SpectrumFit_%s", Profile->GetName()),"R0Q");

  //std::cout << "kk8" << std::endl;
  return spectrum;
}


//////////////////////////////////////////////////////////////////////////

void GetTimeRes(TTree* tree, TF1* FitCh59, TF1* FitCh315, Double_t PedCh59,Double_t PedCh315,TH1D* Tdiff,TF1* FitTdiff){
  
  Double_t time[2];
  Double_t energy[2];
  Double_t mean59=FitCh59->GetParameter(1);
  Double_t mean315=FitCh315->GetParameter(1);
  Double_t RMS59=FitCh59->GetParameter(2);
  Double_t RMS315=FitCh315->GetParameter(2);
  Float_t NS=2.5;

  tree->SetBranchAddress("time",time);
  tree->SetBranchAddress("energy",energy);
  
  std::cout << "_____CutInterval59____  " << mean59-NS*RMS59 << " " << mean59+NS*RMS59 << "  _____CutInterval315____  " << mean315-NS*RMS315 << " " << mean315+NS*RMS315 << std::endl;

  for(int i=0; i<tree->GetEntries();i++){
    tree->GetEntry(i);
    
    if((energy[0]-PedCh59)>(mean59-NS*RMS59) && (energy[0]-PedCh59)<(mean59+NS*RMS59) && (energy[1]-PedCh315)>(mean315-NS*RMS315) && (energy[1]-PedCh315)<(mean315+NS*RMS315)){ 
      Tdiff->Fill(time[0]-time[1]);
    }//chiudo if
    
  }//chiudo for
  TF1* preFit = new TF1("TempFit","gaus");
  

  preFit->SetParameter(0,Tdiff->GetMaximum());
  preFit->SetParameter(1,Tdiff->GetBinCenter(Tdiff->GetMaximumBin()));
  
  Tdiff->Fit(preFit,"Q");

  FitTdiff->SetParameter(0,preFit->GetParameter(0));
  FitTdiff->SetParameter(1,preFit->GetParameter(1));
  FitTdiff->SetParameter(2,preFit->GetParameter(2));

  Tdiff->Fit(FitTdiff,"Q");
  
}
///////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[] ){
  
  std::string DirData(argv[1]);
  std::string RootFileName=DirData;
  RootFileName.erase(0,2);
  std::string OV(argv[2]);

  gSystem->Exec(("ls "+DirData+"/*PED*"+OV+"_singles.root > "+DirData+"/PedFile.txt").c_str());
  gSystem->Exec(("ls "+DirData+"/*PHYS*"+OV+"_coincidences.root > "+DirData+"/PhysCoincFile.txt").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot").c_str());
  
  gSystem->Exec(("mkdir "+DirData+"/Plot/TimeRes").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TimeRes/Partials/").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TimeRes/Tdiff/").c_str());
       
  std::vector<std::string> FileListPedestal;
  std::string ListFilePed = DirData+"/PedFile.txt";
  std::cout << "Lista File Pedestal: "<< ListFilePed << std::endl;    
  FileListPedestal=ReadData(ListFilePed);
  
  std::vector<std::string> FileListPhysicsCoinc;
  std::string ListFilePhys = DirData+"/PhysCoincFile.txt";
  std::cout << "Lista File Data: "<< ListFilePhys << std::endl;
  FileListPhysicsCoinc=ReadData(ListFilePhys);

  int NFilePhys=(int)FileListPhysicsCoinc.size();
  
  //Get Pedestal
  Double_t Pedestal[NFilePhys][2];
  Double_t PedestalCh1[2];
  Double_t PedestalCh2[2];

  int k=0;
  ////////////////////////////////////////////////////CALCOLO PIEDISTALLI
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

  //////////////////////////////////////////////////////SPETTRI COINCIDENZA
  
  for(int i=0;i<NFilePhys;i++){
    std::cout <<NFilePhys <<"   "<< (int)FileListPedestal.size()<< Pedestal[i][0] << "    " << Pedestal[i][1] << std::endl;
  }


  TH1D* HistoCh59[NFilePhys];
  TH1D* HistoCh315[NFilePhys];
    
  Double_t MeanTCh59[NFilePhys], MeanTCh315[NFilePhys];
  Double_t SigmaTCh59[NFilePhys], SigmaTCh315[NFilePhys];
  Double_t MeanTGlobal[NFilePhys], SigmaTGlobal[NFilePhys];
  
  TFile* f0;
  TTree* tree0;

  TF1* FitSpectrum[NFilePhys][2];
  
  TH1D* HistoTdiff[NFilePhys];
  TF1* FitTRes[NFilePhys];
  
  Double_t TimeRes[NFilePhys];
  Double_t SigmaTimeRes[NFilePhys];
  
  gStyle->SetOptFit(111);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.15);

  ///////////////////////////////////////////////////////////////////////////

  for(int i=0;i < NFilePhys;i++){
    
    std::cout << "i : " << i <<" File: " << DirData+"/"+FileListPhysicsCoinc.at(i) << std::endl;
    f0= TFile::Open((DirData+"/"+FileListPhysicsCoinc.at(i)).c_str());
    
    if(f0) std::cout << "File Correctly Opend" << std::endl;
    
    tree0 = (TTree*)f0->Get("data");

        
    HistoCh59[i]  = new TH1D(Form("HistoCh59N%d", i),Form("HistoCh59N%d", i), 100,0,100);
    HistoCh315[i] = new TH1D(Form("HistoCh315N%d",i),Form("HistoCh315N%d",i), 100,0,100);

        
    GetSpectrum(tree0,HistoCh59[i],HistoCh315[i],Pedestal[i][0],Pedestal[i][1]);
    
  
    GetMeanTemperature(tree0,MeanTCh59,SigmaTCh59,MeanTCh315,SigmaTCh315,MeanTGlobal,SigmaTGlobal,i);
  
    HistoTdiff[i] = new TH1D(Form("Tdiff%i",i),Form("Tdiff%i",i),200,-2000,2000);
    FitTRes[i] = new TF1(Form("FitTRes%i",i),"gaus");

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

    canvino->SaveAs((DirData+"/Plot/TimeRes/Partials/Canvas"+std::to_string(i)+".png").c_str());
    delete canvino;

    //    GetTimeRes(TTree* tree, TF1* FitCh59, TF1* FitCh315, Double_t PedCh59,Double_t PedCh315,TH1D* Tdiff,TF1* FitTdiff,int j)
    
    GetTimeRes(tree0,FitSpectrum[i][0],FitSpectrum[i][1],Pedestal[i][0],Pedestal[i][1],HistoTdiff[i],FitTRes[i]);
    
    TimeRes[i]= FitTRes[i]->GetParameter(2);
    SigmaTimeRes[i] = FitTRes[i]->GetParError(2);
    
    TCanvas* CanvTimeRes = new TCanvas("TResCanv","TResCanv",1600,800);
    HistoTdiff[i]->GetXaxis()->SetTitle("T_{ch59}-T_{ch315}");
    HistoTdiff[i]->GetYaxis()->SetTitle("Counts");
    HistoTdiff[i]->Draw();
    FitTRes[i]->Draw("SAME");

    CanvTimeRes->SaveAs((DirData+"/Plot/TimeRes/Tdiff/TdiffN"+std::to_string(i)+".png").c_str());
    delete CanvTimeRes;
  }

  TCanvas* TResVsTemp = new TCanvas("TResVsTemp","TResVsTemp",1600,800);
  TGraphErrors* GraphTResVsTemp = new TGraphErrors(NFilePhys,MeanTGlobal,TimeRes,SigmaTGlobal,SigmaTimeRes);
  GraphTResVsTemp->SetTitle("GraphTimeResolutionVsTemp");
  GraphTResVsTemp->GetXaxis()->SetTitle("TempMeanBox [°C]");
  GraphTResVsTemp->GetYaxis()->SetTitle("TimeResolution [ps]");

  GraphTResVsTemp->Draw();

  TResVsTemp->SaveAs((DirData+"/Plot/TimeRes/TimeResVsTemp.png").c_str());

  return 0;
}
