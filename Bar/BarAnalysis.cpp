#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>

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
#include "TCut.h"

std::vector<std::string> ReadData(std::string FileName){
  
  std::vector<std::string> FileList;
  std::ifstream file(FileName);
  std::string str; 
  
  if(!file.is_open()){std::cout << "file non trovato" << std::endl;}
  
  while (std::getline(file, str))
    {
      FileList.push_back(str);
    }

  return FileList;
}

void GetPedestal(TTree*,Double_t*, std::vector<int>,std::string,int);
void GetSpectrum(TTree* ,TH1D**,Double_t* ,std::vector<int>);
TF1* FitCoincSpectrum(TH1D*);
void GetTdiff(TTree* , TH1D** , TF1** , Double_t* ped, std::vector<int>, std::vector<std::string>);

int main(int argc, char* argv[] ){
  
  std::cout << "Analisys" << std::endl;
  
  gStyle->SetOptFit(1111);

  std::string DirData(argv[1]);
  std::string RootFileName=DirData;
  RootFileName.erase(0,2);
  std::string OV(argv[2]);
  
  gSystem->Exec(("ls "+DirData+"/*PED*"+OV+"_singles.root > "+DirData+"/PedFile.txt").c_str());
  gSystem->Exec(("ls "+DirData+"/*PHYS*"+OV+"_coincidences.root > "+DirData+"/PhysCoicFile.txt").c_str());
  gSystem->Exec(("ls "+DirData+"/*PHYS*"+OV+".txt > "+DirData+"/ConfigFile.txt").c_str());
  
  gSystem->Exec(("mkdir "+DirData+"/Plot").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/Pedestal").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/Spectra").c_str());    
  
  
  //gStyle->SetOptStat("000001000");
  
  TFile* f = new TFile(("../RootFileGraph/"+RootFileName+".root").c_str(),"RECREATE");

  if(f) std::cout << "Root File Graph opened->"<<RootFileName << std::endl;

  std::vector<std::string> FileListPedestal;
  std::string ListFilePed = DirData+"/PedFile.txt";
  std::cout << "Lista File Pedestal: "<< ListFilePed << std::endl;    
  FileListPedestal=ReadData(ListFilePed);

  std::vector<std::string> FileListPhysics;
  std::string ListFilePhys = DirData+"/PhysCoicFile.txt";
  std::cout << "Lista File Pedestal: "<< ListFilePed << std::endl;
  FileListPhysics=ReadData(ListFilePhys);

  std::vector<std::string> FileListConfig;
  std::string ListFileConfig = DirData+"/ConfigFile.txt";
  std::cout << "Lista File Configuration: "<< ListFileConfig << std::endl;    
  FileListConfig=ReadData(ListFileConfig);
  
  for(int i=0; i< (int )FileListConfig.size();i++){
    std::cout << FileListConfig.at(i) << std::endl;
  }
  
  std::vector<int> Channels;
  std::vector<std::string> ChannelType;
  std::vector<std::string> ConfigData;
  ConfigData=ReadData(FileListConfig[0]);

  const int size =30 ;
  
  for(int i=0; i<(int)ConfigData.size();i++){
    std::string arr[size];
    int k = 0;
    std::stringstream ssin(ConfigData.at(i));
    while (ssin.good() && k < size){
      ssin >> arr[k];
      std::cout << "instram    " << arr[k] << std::endl; 
      ++k;
    }
    if(arr[0]=="CH"){ 
      Channels.push_back(std::atoi((arr[5]).c_str())*64+std::atoi((arr[6]).c_str()));
      ChannelType.push_back(arr[15]);
    }
    
  }//chiudo for

  for(int i=0;i< (int)Channels.size(); i++){
    std::cout<< Channels[i] <<"  " <<ChannelType[i]<< std::endl;
  }
  
  std::map<int,int> mapCh;
  
  for(int i=0;i<(int)Channels.size();i++){
    mapCh[Channels[i]]=i;
    std::cout << "originalmap  " << Channels[i] << " " << mapCh[Channels[i]] << std::endl;
  }
  

  int NFilePhys=(int)FileListPhysics.size();
  
  Double_t Pedestal[NFilePhys][(int)Channels.size()];
  Double_t PedestalBefore[(int)Channels.size()];
  Double_t PedestalAfter[(int)Channels.size()];

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  int k=0;
  for(int i=0;i < (int)FileListPedestal.size()-1;i+=2){
    
    std::cout<< "open file ped:  " << (DirData+"/"+FileListPedestal.at(i)).c_str() << std::endl;
    TFile* f0= TFile::Open((DirData+"/"+FileListPedestal.at(i)).c_str());
    TFile* f1= TFile::Open((DirData+"/"+FileListPedestal.at(i)).c_str());
    
    TTree* tree0 = (TTree*)f0->Get("data"); //Before
    TTree* tree1 = (TTree*)f1->Get("data"); //After
    
    GetPedestal(tree0,PedestalBefore,Channels,DirData,i);
    GetPedestal(tree1,PedestalAfter,Channels,DirData,i+1);
  
    for(int l=0; l<(int)Channels.size();l++){
      Pedestal[k][mapCh[Channels[l]]]=(PedestalBefore[mapCh[Channels[l]]]+PedestalAfter[mapCh[Channels[l]]])/2;
      
    }
      
    std::cout << Pedestal[k][mapCh[Channels[0]]] << " " << Pedestal[k][mapCh[Channels[1]]] << " " << Pedestal[k][mapCh[Channels[2]]] << std::endl;
    k++;
  }// chiudo for
  ///////////////////////////////////////////////////////////////////////////////
  
  TH1D* Spectrum[NFilePhys][(int)Channels.size()];
  TFile* f0;
  TTree* tree0;
  TF1* FitSpectrum[NFilePhys][(int)Channels.size()];
  TCanvas* canvino;
  
  TH1D* tdiff[3];

  tdiff[0] = new TH1D("tref-tave","tref-tave",2000,-2000,2000);
  tdiff[1] = new TH1D("tref-t1","tref-t1",2000,-2000,2000);
  tdiff[2] = new TH1D("tref-t2","tref-t2",2000,-2000,2000);
  
  for(int i=0;i<NFilePhys;i++){
   
    
    f0= TFile::Open((DirData+"/"+FileListPhysics.at(i)).c_str());
    tree0 = (TTree*)f0->Get("data");
    
    for(int j=0; j< (int)Channels.size();j++){
      
      Spectrum[i][mapCh[Channels[j]]] = new TH1D(Form("Ch%i_%i_%s",Channels[j],i,(ChannelType[mapCh[Channels[j]]]).c_str()),Form("Ch%i_%i_%s",Channels[j],i,(ChannelType[mapCh[Channels[j]]]).c_str()),100,0,100);
      
    }//CHIUDO FOR J
    
    GetSpectrum(tree0,Spectrum[i],Pedestal[i],Channels);
    
    for(int k=0; k<(int)Channels.size(); k++){
      FitSpectrum[i][mapCh[Channels[k]]]=FitCoincSpectrum(Spectrum[i][mapCh[Channels[k]]]);
    }
    
    canvino = new TCanvas("canvino","canvino",1500,700);
    canvino->Divide((int)Channels.size(),1);
    
    for(int k=0; k<(int)Channels.size(); k++){
      canvino->cd(k+1);
      Spectrum[i][mapCh[Channels[k]]]->GetXaxis()->SetRange(1,100);
      Spectrum[i][mapCh[Channels[k]]]->Draw();
      FitSpectrum[i][mapCh[Channels[k]]]->Draw("SAME");
    }
    
    canvino->SaveAs((DirData+"/Plot/Spectra/"+"Spectrum"+std::to_string(i)+".png").c_str());
      
    delete canvino;

    //  GetTdiff(TTree* tree, TH1D** histo, TF1** fitspectrum, Double_t* ped, std::vector<int> NCH, std::vector<std::string> ChType) 
    GetTdiff(tree0,tdiff,FitSpectrum[i],Pedestal[i],Channels,ChannelType);


   }//CHIUDO FOR
  
 
}//////////////////////////////CHIUDO MAIN

///////////////////////////////////////////////////////////////////////////////////////////////
void GetPedestal(TTree* tree,Double_t* ped,std::vector<int> NCH,std::string DirData,int index){
  
  TH1D* histo[(int)NCH.size()]; 
  TCanvas* CavasPedestal = new TCanvas("canvasPedestal","canvasPedestal",1400,600);
  CavasPedestal->Divide((int)NCH.size(),1);
  
  std::map<int,int> mapCh;
  
  for(int i=0;i<(int)NCH.size();i++){
    mapCh[NCH[i]]=i;
    histo[i]= new TH1D(Form("pedestal%i_%i",NCH[i],index),Form("pedestal%i_%i",NCH[i],index),600,0,600);
  }
  
  Float_t energy;
  UShort_t chID;

  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("channelID",&chID);
   
  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    histo[mapCh[chID]]->Fill(energy);
    
  }//chiudo for
  
  std::cout << "Filled"<< std::endl;

  for(int i=0;i<(int)NCH.size();i++){
    
    ped[mapCh[NCH[i]]]= histo[mapCh[NCH[i]]]->GetMean();
    CavasPedestal->cd(1+i);
    histo[mapCh[NCH[i]]]->Draw();
    //std::cout << "PED NCH" << NCH[i] << "->" << mapCh[NCH[i]] << std::endl;
  }
  
  CavasPedestal->SaveAs((DirData+"/Plot/Pedestal/"+"pedestal"+std::to_string(index)+".png").c_str());
  
  delete CavasPedestal;
  
}

///////////////////////////////////////////////////////////////////////////////////////////////
void GetSpectrum(TTree* tree, TH1D** histo,Double_t* ped, std::vector<int> NCH){

  
  Double_t energy[(int)NCH.size()];
  Double_t chID[(int)NCH.size()];
  Int_t check;
  
  std::map<int,int> mapCh;

  tree->SetBranchAddress("energy",energy);
  tree->SetBranchAddress("chId",chID);
  
  for(int i=0;i<(int)NCH.size();i++){
    mapCh[NCH[i]]=i;
    //  std::cout << "PHYS NCH" << NCH[i] << "->" << mapCh[NCH[i]] << std::endl; 
  }
  
  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    check=0;
    
    for(int j=0;j<(int)NCH.size();j++){
      if(chID[mapCh[NCH[j]]]==-9) {check++;}
    }//CHIUDO FOR
    
    if(check==0){
      
      for(int j=0;j<(int)NCH.size();j++){
	histo[mapCh[chID[j]]]->Fill(energy[mapCh[chID[j]]]-ped[mapCh[chID[j]]]);
      }//CHIUDO FOR J
      
    }//CHIUDO IF
    
  }//chiudo for
  
  
}
///////////////////////////////////////////////////////////////////////////////////////////////

void GetTdiff(TTree* tree, TH1D** histo, TF1** fitspectrum, Double_t* ped, std::vector<int> NCH, std::vector<std::string> ChType){
  //Histo0 tref-tave, histo1 tref-t1, histo2 tref-t2
  Double_t energy[(int)NCH.size()];
  Double_t chID[(int)NCH.size()];
  Double_t time[(int)NCH.size()];
  
  std::map<int,int> mapCh;
  std::map<std::string,int> mapTime;
  
  int NTime=0;
  int NTimeRef=0;
  int h=1;

  tree->SetBranchAddress("energy",energy);
  tree->SetBranchAddress("chId",chID);
  tree->SetBranchAddress("time",time);  

  for(int i=0;i<(int)NCH.size();i++){
    mapCh[NCH[i]]=i;
    if(ChType[mapCh[NCH[i]]]=="LYSO_bar"){
      mapTime["time"+std::to_string(h)]=i;
      h++;
    }//CHIUDO IF
    else{ 
      mapTime["timeRef"]=i; 
    }//CHIUDO ELSE
  }//CHIUDO FOR
  
  TCut cut[(int)NCH.size()];
  TCut TotalCut("1");
  for(int i=0;i< (int)NCH.size();i++){
    cut[mapCh[NCH[i]]]=Form("energy[mapCh[NCH[%i]]]-ped[mapCh[NCH[%i]]]<fitspectrum[mapCh[NCH[%i]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[%i]]]->GetParameter(2) && energy[mapCh[NCH[%i]]]-ped[mapCh[NCH[%i]]])>fitspectrum[mapCh[NCH[%i]]]->GetParameter(1)-3*fitspectrum[mapCh[NCH[%i]]]->GetParameter(2)",i,i,i,i,i,i,i,i);
    std::cout << cut[mapCh[NCH[i]]] << std::endl;
    TotalCut+=cut[mapCh[NCH[i]]];
  }
  
  for(int i=0;i< tree->GetEntries();i++){
    tree->GetEntry(i);
    //Histo0 tref-tave, histo1 tref-t1, histo2 tref-t2
    if( TotalCut){
      histo[0]->Fill(time[mapCh[NCH[mapTime["timeRef"]]]]-(time[mapCh[NCH[mapTime["time1"]]]]+time[mapCh[NCH[mapTime["time2"]]]])/2);
      histo[1]->Fill(time[mapCh[NCH[mapTime["timeRef"]]]]-time[mapCh[NCH[mapTime["time1"]]]]);
      histo[2]->Fill(time[mapCh[NCH[mapTime["timeRef"]]]]-time[mapCh[NCH[mapTime["time2"]]]]);
    }
    
  }


  
  std::cout << Form("time1 %i ",  NCH[mapTime["time1"]]) << Form("time2 %i ",     NCH[mapTime["time2"]]) << Form("timeRef %i ", NCH[mapTime["timeRef"]]) << std::endl;
  std::cout << Form("time1 %s ",  ChType[mapTime["time1"]].c_str()) << Form("time2 %s ",    ChType[mapTime["time2"]].c_str()) << Form("timeRef %s ", ChType[mapTime["timeRef"]].c_str()) << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////
TF1* FitCoincSpectrum(TH1D* Profile){

  Double_t max;
  Double_t peak1;
  
  peak1 = Profile->GetBinCenter(Profile->GetMaximumBin());
  max = Profile->GetMaximum();
  Profile->GetXaxis()->UnZoom();

  TF1* spectrum = new TF1(Form("SpectrumFit_%s", Profile->GetName()),"[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) + [3] / (exp( (x*[4]-(2*[1]*[1]/([1]+2*[1])))) + 1)",10,60);
  
  spectrum->SetParameter(0,max);
  spectrum->SetParameter(1,peak1);
  spectrum->SetParameter(2,3);
  spectrum->SetParameter(3,700);
  spectrum->SetParameter(4,0.82);
  
  Profile->Fit(Form("SpectrumFit_%s", Profile->GetName()),"R0");

  return spectrum;
}

