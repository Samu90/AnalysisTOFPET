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

void GetPedestal(TTree*,Double_t*,Double_t*, std::vector<int>,std::string,int);
void GetSpectrum(TTree* ,TH1D**,Double_t* ,std::vector<int>);
TF1* FitCoincSpectrum(TH1D*);
void GetTdiff(TTree* , TH1D**,TH2D** ,TH2D**,TH2D*,Double_t*, TF1** , Double_t*, std::vector<int>, std::vector<std::string>);
void GetTdiffCorr(TTree* ,TH2D**,TH2D* , TF1** ,std::vector<TF1*>, Double_t*,Double_t* , std::vector<int> , std::vector<std::string> );
void ProjectionTemp(TH2D**,std::string);


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
  gSystem->Exec(("mkdir "+DirData+"/Plot/TR").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TR/Pedestal").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TR/Spectra").c_str());    
  gSystem->Exec(("mkdir "+DirData+"/Plot/TR/ProjectionAEff").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TR/ProjectionTemp").c_str());
  
  //gStyle->SetOptStat("000001000");
  
  TFile* f = new TFile(("../RootFileGraphBar/"+RootFileName+".root").c_str(),"RECREATE");

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
  int NFilePed = (int)FileListPedestal.size()-1;

  /*  
   NFilePhys=1;
   NFilePed= NFilePhys*2;
  */

  Double_t Pedestal[NFilePhys][(int)Channels.size()];
  Double_t RMSPedestal[NFilePhys][(int)Channels.size()];
  Double_t PedestalBefore[(int)Channels.size()];
  Double_t PedestalAfter[(int)Channels.size()];
  Double_t RMSPedestalBefore[(int)Channels.size()];
  Double_t RMSPedestalAfter[(int)Channels.size()];
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  int k=0;

  for(int i=0; i < NFilePed ; i+=2){
    
    std::cout<< "open file ped:  " << (DirData+"/"+FileListPedestal.at(i)).c_str() << std::endl;
    TFile* f0= TFile::Open((DirData+"/"+FileListPedestal.at(i)).c_str());
    TFile* f1= TFile::Open((DirData+"/"+FileListPedestal.at(i)).c_str());
    
    TTree* tree0 = (TTree*)f0->Get("data"); //Before
    TTree* tree1 = (TTree*)f1->Get("data"); //After
    
    GetPedestal(tree0,PedestalBefore,RMSPedestalBefore,Channels,DirData,i);
    GetPedestal(tree1,PedestalAfter,RMSPedestalAfter,Channels,DirData,i+1);
    
    for(int l=0; l<(int)Channels.size();l++){
      Pedestal[k][mapCh[Channels[l]]]=(PedestalBefore[mapCh[Channels[l]]]+PedestalAfter[mapCh[Channels[l]]])/2;
      RMSPedestal[k][mapCh[Channels[l]]]=(RMSPedestalBefore[mapCh[Channels[l]]]+RMSPedestalAfter[mapCh[Channels[l]]])/2;
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

  tdiff[0] = new TH1D("tref-tave","tref-tave",200,-2000,2000);
  tdiff[1] = new TH1D("tref-t1","tref-t1",200,-2000,2000);
  tdiff[2] = new TH1D("tref-t2","tref-t2",200,-2000,2000);

  TH2D* tdiffVsE[2];

  tdiffVsE[0] = new TH2D("tdiffVsE_time1","tdiffVsE_time1",70,0,70,200,-3000,3000);
  tdiffVsE[1] = new TH2D("tdiffVsE_time2","tdiffVsE_time2",70,0,70,200,-3000,3000);

  TH2D* tdiffVsT[3];

  tdiffVsT[0] = new TH2D("tdiffVsT_time1","tdiffVsT_time1",40,21,29,200,-2000,2000);
  tdiffVsT[1] = new TH2D("tdiffVsT_time2","tdiffVsT_time2",40,21,29,200,-2000,2000);
  tdiffVsT[2] = new TH2D("tdiffVsT_tave","tdiffVsT_time2",40,21,29,200,-2000,2000);
  
  TH2D* tdiffVsAmpEff = new TH2D("tdiffVsAmpEff","tdiffVsAmpEff",50,-10,40,100,-3500,3000);
  
  TLine* line[(int)Channels.size()][2];
    
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
      
      Double_t min=FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(1)-2.5*FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(2);
      Double_t max=FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(1)+2.5*FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(2);
      
      line[mapCh[Channels[k]]][0]= new TLine(min,0,min,Spectrum[i][mapCh[Channels[k]]]->GetMaximum());
      line[mapCh[Channels[k]]][1]= new TLine(max,0,max,Spectrum[i][mapCh[Channels[k]]]->GetMaximum());
      
      canvino->cd(k+1);
      Spectrum[i][mapCh[Channels[k]]]->GetXaxis()->SetRange(1,100);
      Spectrum[i][mapCh[Channels[k]]]->GetXaxis()->SetTitle("Energy [D.U.]");
      Spectrum[i][mapCh[Channels[k]]]->GetYaxis()->SetTitle("Counts");
      Spectrum[i][mapCh[Channels[k]]]->Draw();
      FitSpectrum[i][mapCh[Channels[k]]]->Draw("SAME");
      line[mapCh[Channels[k]]][0]->Draw("SAME");
      line[mapCh[Channels[k]]][1]->Draw("SAME");
      
    }
    
    canvino->SaveAs((DirData+"/Plot/TR/Spectra/"+"Spectrum"+std::to_string(i)+".png").c_str());
    
    delete canvino;
    
    GetTdiff(tree0,tdiff,tdiffVsE,tdiffVsT,tdiffVsAmpEff,RMSPedestal[i],FitSpectrum[i],Pedestal[i],Channels,ChannelType);
    
    
   }//CHIUDO FOR
  
  TF1* FitTreso[3];
  for(int i=0; i<(int)Channels.size();i++){
    FitTreso[i] = new TF1(Form("t_reso%i",i),"gaus",tdiff[i]->GetMean()-1.5*tdiff[i]->GetRMS(),tdiff[i]->GetMean()+1.5*tdiff[i]->GetRMS());
  }

  TCanvas* TimeResoCavas = new TCanvas("TimeResoCavas","TimeResoCavas", 1700,800);
  TimeResoCavas->Divide((int)Channels.size(),1);
  for(int i=0; i<(int)Channels.size();i++){
    TimeResoCavas->cd(1+i);
    tdiff[i]->GetXaxis()->SetTitle("tdiff [ps]");
    tdiff[i]->GetYaxis()->SetTitle("counts");
    tdiff[i]->Draw();
    tdiff[i]->Fit(FitTreso[i],"R");
  }
  
  TimeResoCavas->SaveAs((DirData+"/Plot/TR/TimeReso.png").c_str());

  
  TCanvas* TimeResoLogCavas = new TCanvas("TimeResoLogCavas","TimeResoLogCavas", 1700,800);
  TimeResoLogCavas->Divide((int)Channels.size(),1);
  for(int i=0; i<(int)Channels.size();i++){
    TimeResoLogCavas->cd(1+i)->SetLogy();
    tdiff[i]->Draw();
    tdiff[i]->Fit(FitTreso[i]);
  }
  
  TimeResoLogCavas->SaveAs((DirData+"/Plot/TR/TimeResoLog.png").c_str());
  
  std::vector<TF1*> fitTdiffVsE;
  TF1* pippo1 = new TF1("fitTdiffVsE0pippo","pol5");
  TF1* pippo2 = new TF1("fitTdiffVsE1pippo","pol5");

  fitTdiffVsE.push_back(pippo1);
  fitTdiffVsE.push_back(pippo2);

  TCanvas* CanvasTdiffVsE = new TCanvas("CanvasTdiffVsE","CanvasTdiffVsE",1500,700);
  CanvasTdiffVsE->Divide((int)Channels.size()-1,1);
  CanvasTdiffVsE->cd(1);
  
  tdiffVsE[0]->GetXaxis()->SetTitle("Energy [D.U.]");
  tdiffVsE[0]->GetYaxis()->SetTitle("time1-timeRef [ps]");
  tdiffVsE[0]->GetYaxis()->SetTitleOffset(0.6);
  tdiffVsE[0]->Draw("COLZ");
  tdiffVsE[0]->Fit( fitTdiffVsE[0] );
  
  std::cout <<"FitAddress1:" <<fitTdiffVsE[0] << std::endl;
  
  CanvasTdiffVsE->cd(2);
  tdiffVsE[1]->GetXaxis()->SetTitle("Energy [D.U.]");
  tdiffVsE[1]->GetYaxis()->SetTitle("time2-timeRef [ps]");
  tdiffVsE[1]->GetYaxis()->SetTitleOffset(0.6);
  tdiffVsE[1]->Draw("COLZ");
  tdiffVsE[1]->Fit( fitTdiffVsE[1] );
  
  std::cout <<"FitAddress2:" <<fitTdiffVsE[1] << std::endl;
  
  std::cout << "Address1_:" << fitTdiffVsE[0] << " "<< fitTdiffVsE[0]->GetParameter(0) << " "<< fitTdiffVsE[0]->GetParameter(1)<< " "<< fitTdiffVsE[0]->GetParameter(2)<< " "<< fitTdiffVsE[0]->GetParameter(3)<< " "<< fitTdiffVsE[0]->GetParameter(4)<< " "<< fitTdiffVsE[0]->GetParameter(5)<< std::endl;
  
  std::cout << "Address2_:" << fitTdiffVsE[1] << " "<< fitTdiffVsE[1]->GetParameter(0) << " "<< fitTdiffVsE[1]->GetParameter(1)<< " "<< fitTdiffVsE[1]->GetParameter(2)<< " "<< fitTdiffVsE[1]->GetParameter(3)<< " "<< fitTdiffVsE[1]->GetParameter(4)<< " "<< fitTdiffVsE[1]->GetParameter(5) << std::endl;

  

  CanvasTdiffVsE->SaveAs((DirData+"/Plot/TR/TdiffVsE.png").c_str());

  TF1* fitTdiffVsT[2];
  fitTdiffVsT[0]= new TF1("fitTdiffVsT1","pol1");
  fitTdiffVsT[1]= new TF1("fitTdiffVsT2","pol1");
  fitTdiffVsT[2]= new TF1("fitTdiffVsTave","pol1");
  
  
  TCanvas* CanvasTdiffVsT = new TCanvas("CanvasTdiffVsT","CanvasTdiffVsT",1500,700);
  CanvasTdiffVsT->Divide((int)Channels.size(),1);
  
  CanvasTdiffVsT->cd(1);
  tdiffVsT[0]->GetXaxis()->SetTitle("Temperature [°C]");
  tdiffVsT[0]->GetYaxis()->SetTitle("time1-timeRef [ps]");
  tdiffVsT[0]->Draw("COLZ");
  tdiffVsT[0]->Fit(fitTdiffVsT[0]);

  CanvasTdiffVsT->cd(2);
  tdiffVsT[1]->GetXaxis()->SetTitle("Temperature [°C]");
  tdiffVsT[1]->GetYaxis()->SetTitle("time2-timeRef [ps]");
  tdiffVsT[1]->Draw("COLZ");
  tdiffVsT[1]->Fit(fitTdiffVsT[1]);

  CanvasTdiffVsT->cd(3);
  tdiffVsT[2]->GetXaxis()->SetTitle("Temperature [°C]");
  tdiffVsT[2]->GetYaxis()->SetTitle("tave-timeRef [ps]");
  tdiffVsT[2]->Draw("COLZ");
  tdiffVsT[2]->Fit(fitTdiffVsT[2]);

  CanvasTdiffVsT->SaveAs((DirData+"/Plot/TR/TdiffVsT.png").c_str());
  
  ProjectionTemp(tdiffVsT, DirData);
  
  /////////////////CORREZIONE//////////////////////
  
  TH1D* tdiffCorr[3];
  
  tdiffCorr[0] = new TH1D("tref-tave_corr","tref-tave_corr",200,-2000,2000);
  tdiffCorr[1] = new TH1D("tref-t1_corr","tref-t1_corr",200,-2000,2000);
  tdiffCorr[2] = new TH1D("tref-t2_corr","tref-t2_corr",200,-2000,2000);

  TH2D* tdiffCorrVsE[2];

  tdiffCorrVsE[0] = new TH2D("tdiffCorrVsE_time1","tdiffCorrVsE_time1",70,0,70,200,-3000,3000);
  tdiffCorrVsE[1] = new TH2D("tdiffCorrVsE_time2","tdiffCorrVsE_time2",70,0,70,200,-3000,3000);
  
  TH2D* tdiffCorrVsAmpEff = new TH2D("tdiffCorrVsAmpEff","tdiffCorrVsAmpEff",50,-10,40,150,-3500,3000);


  for(int i=0;i<NFilePhys;i++){
    
    f0= TFile::Open((DirData+"/"+FileListPhysics.at(i)).c_str());
    tree0 = (TTree*)f0->Get("data"); 
    
    //std::cout<<"before function " << fitTdiffVsE[0] << " " << fitTdiffVsE[1] << std::endl;
    
    GetTdiffCorr(tree0,tdiffCorrVsE,tdiffCorrVsAmpEff,FitSpectrum[i],fitTdiffVsE,Pedestal[i],RMSPedestal[i],Channels,ChannelType);
    // GetTdiffCorr(TTree* tree, TH2D** tdiffVsE ,TF1** fitspectrum,TF1** fitcorr, Double_t* ped, std::vector<int> NCH, std::vector<std::string> ChType)
  }

   

  TF1* FitTresoCorr[3];
  for(int i=0; i<(int)Channels.size();i++){
    FitTresoCorr[i] = new TF1(Form("t_reso_corr%i",i),"gaus",tdiffCorr[i]->GetMean()-1.5*tdiffCorr[i]->GetRMS(),tdiffCorr[i]->GetMean()+1.5*tdiffCorr[i]->GetRMS());
  }

  TCanvas* TimeResoCorrCavas = new TCanvas("TimeResoCorrCavas","TimeResoCorrCavas", 1700,800);
  TimeResoCorrCavas->Divide((int)Channels.size(),1);
  for(int i=0; i<(int)Channels.size();i++){
    TimeResoCorrCavas->cd(1+i);
    tdiffCorr[i]->GetXaxis()->SetTitle("tdiffCorr [ps]");
    tdiffCorr[i]->GetYaxis()->SetTitle("counts");
    tdiffCorr[i]->Draw();
    tdiffCorr[i]->Fit(FitTresoCorr[i],"R");
  }
  
  TimeResoCorrCavas->SaveAs((DirData+"/Plot/TR/TimeResoCorr.png").c_str());
  
 TCanvas* TimeResoCorrLogCavas = new TCanvas("TimeResoCorrLogCavas","TimeResoCorrLogCavas", 1700,800);
  TimeResoCorrLogCavas->Divide((int)Channels.size(),1);
  for(int i=0; i<(int)Channels.size();i++){
    TimeResoCorrLogCavas->cd(1+i)->SetLogy();
    tdiffCorr[i]->Draw();
    tdiffCorr[i]->Fit(FitTreso[i]);
  }
  
  TimeResoCorrLogCavas->SaveAs((DirData+"/Plot/TR/TimeResoCorrLog.png").c_str());

  TF1* fitTdiffCorrVsE[2];
  fitTdiffCorrVsE[0] = new TF1("fitTdiffCorrVsE0","pol5");
  fitTdiffCorrVsE[1] = new TF1("fitTdiffCorrVsE1","pol5");

  TCanvas* CanvasTdiffCorrVsE = new TCanvas("CanvasTdiffCorrVsE","CanvasTdiffCorrVsE",1500,700);
  
  CanvasTdiffCorrVsE->Divide((int)Channels.size()-1,1);
  
  CanvasTdiffCorrVsE->cd(1);  
  tdiffCorrVsE[0]->GetXaxis()->SetTitle("Energy [D.U.]");
  tdiffCorrVsE[0]->GetYaxis()->SetTitle("time1-timeRef [ps]");
  tdiffCorrVsE[0]->Draw("COLZ");
  tdiffCorrVsE[0]->Fit(fitTdiffCorrVsE[0]);

  CanvasTdiffCorrVsE->cd(2);
  tdiffCorrVsE[1]->GetXaxis()->SetTitle("Energy [D.U.]");
  tdiffCorrVsE[1]->GetYaxis()->SetTitle("time2-timeRef [ps]");
  tdiffCorrVsE[1]->Draw("COLZ");
  tdiffCorrVsE[1]->Fit(fitTdiffCorrVsE[1]);

  
  CanvasTdiffCorrVsE->SaveAs((DirData+"/Plot/TR/TdiffCorrVsE.png").c_str());
  
  TCanvas* CanvasAEff = new TCanvas("CanvasAEff","CanvasAEff",900,1000);
  tdiffVsAmpEff->GetXaxis()->SetTitle("A_Eff/#sigma_{n} [D.U]");
  tdiffVsAmpEff->GetYaxis()->SetTitle("T_diff [ps]");
  tdiffVsAmpEff->GetYaxis()->SetTitleOffset(0.7);
  tdiffVsAmpEff->Draw("COLZ");

  CanvasAEff->SaveAs((DirData+"/Plot/TR/TdiffVsAEff.png").c_str());  
  
  TCanvas* CanvasTdiffCorrVsAmpEff = new TCanvas("CanvasTdiffCorrVsAmpEff","CanvasTdiffCorrVsAmpEff",900,1000);
  tdiffCorrVsAmpEff->GetXaxis()->SetTitle("AmpEff/#sigma_{n}[D.U.]");
  tdiffCorrVsAmpEff->GetYaxis()->SetTitle("t_ave-t_ref[ps]");
  tdiffCorrVsAmpEff->GetYaxis()->SetTitleOffset(0.7);
  tdiffCorrVsAmpEff->Draw("COLZ");
  CanvasTdiffCorrVsAmpEff->SaveAs((DirData+"/Plot/TR/TdiffCorrVsAmpEff.png").c_str());

 
  Int_t nbins = tdiffVsAmpEff->GetXaxis()->GetNbins();
  
  TH1D* Projection;
  TF1* FitProjection;
  TCanvas* CanvasProjection;

  std::vector<Double_t> AEffValue;
  std::vector<Double_t> TimeResValue;
  std::vector<Double_t> ErrTimeResValue;

  Int_t fitStatus;

  for(int i=0; i<nbins-1; i+=2){


    CanvasProjection = new TCanvas("CanvasProjection","CanvasProjection",900,1000);
    
    Projection = tdiffCorrVsAmpEff->ProjectionY(Form("ProjectionCorr%i_Amp%lf",i,tdiffCorrVsAmpEff->GetXaxis()->GetBinCenter(i)),i,i);
    Projection->GetXaxis()->SetTitle("t_diff [ps]");
    Projection->GetYaxis()->SetTitle("Counts");
    Projection->SetTitle(Projection->GetName());
    Projection->Draw();
    
    FitProjection = new TF1(Form("fitProjectionCorr%i",i) , "gaus" , Projection->GetBinCenter(Projection->GetMaximumBin())-1.1*Projection->GetRMS(),Projection->GetBinCenter(Projection->GetMaximumBin())+1.1*Projection->GetRMS());
    
    std::cout << "Range_______:    " << Projection->GetBinCenter(Projection->GetMaximumBin())-1.1*Projection->GetRMS() << "   " << Projection->GetBinCenter(Projection->GetMaximumBin())+1.1*Projection->GetRMS() << std::endl;
          
    FitProjection->SetParameter(0 , Projection->GetMaximum());
    FitProjection->SetParameter(1 , Projection->GetMaximumBin());
    
    fitStatus = Projection->Fit(FitProjection,"RQ");
    
    CanvasProjection->SaveAs((DirData+"/Plot/TR/ProjectionAEff/Projection"+std::to_string(i)+".png").c_str());
    
    if(fitStatus==0 && tdiffCorrVsAmpEff->GetXaxis()->GetBinCenter(i)>0){
      AEffValue.push_back(tdiffCorrVsAmpEff->GetXaxis()->GetBinCenter(i));
      TimeResValue.push_back(sqrt(FitProjection->GetParameter(2)*FitProjection->GetParameter(2)-93*93));
      ErrTimeResValue.push_back(FitProjection->GetParError(2));
    }
    
    
    delete Projection;
    delete CanvasProjection;

    std::cout << tdiffCorrVsAmpEff->GetXaxis()->GetBinCenter(i) << " " << FitProjection->GetParameter(2) << " " << FitProjection->GetParError(2) << std::endl;
  }
  
  TCanvas* CanvasTimeResVsAEff = new TCanvas("CanvasTimeResVsAEff","CanvasTimeResVsAEff",900,1000);
  
  
  TGraphErrors* TimeResVsAEff = new TGraphErrors((int)AEffValue.size(),&AEffValue[0],&TimeResValue[0],0,&ErrTimeResValue[0]);
  
  TF1* fitAEff = new TF1("fitAEff", "sqrt([0]/(x*x)+sqrt(2)*[1]*[1])",1.8,22);
  
  fitAEff->SetParameter(0,2.58e6);
  fitAEff->SetParameter(1,70);
  fitAEff->SetParLimits(1,0,500);
  
  TimeResVsAEff->SetTitle("TimeResVsAEff");
  TimeResVsAEff->SetName(TimeResVsAEff->GetTitle());
  TimeResVsAEff->GetXaxis()->SetTitle("A_Eff/#sigma_{n} [D.U]");
  TimeResVsAEff->GetYaxis()->SetTitle("TimeRes [ps]");
  
  TimeResVsAEff->SetMarkerStyle(8);
  
  TimeResVsAEff->Draw("AP");
  TimeResVsAEff->Fit(fitAEff,"RW");
  
  CanvasTimeResVsAEff->SaveAs((DirData+"/Plot/TR/TimeResVsAEff.png").c_str()); 
  
  
  TFile* f1 = new TFile((DirData+"/Plot/TResVsAEff.root").c_str(),"RECREATE");
  f1->cd();
  
  TimeResVsAEff->Write();
  fitAEff->Write();

  f1->Save();
  f1->Close();
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}//////////////////////////////CHIUDO MAIN



///////////////////////////////////////////////////////////////////////////////////////////////


void GetPedestal(TTree* tree,Double_t* ped,Double_t* RMSPed,std::vector<int> NCH,std::string DirData,int index){
  
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
    
    if( NCH[i]==288 ){ ped[mapCh[NCH[i]]]= histo[mapCh[NCH[i]]]->GetMean()-6; }  //rimuovi quando il canale funziona
    else{ ped[mapCh[NCH[i]]]= histo[mapCh[NCH[i]]]->GetMean(); }   //rimuovi quando il canale funziona
    
    //ped[mapCh[NCH[i]]]= histo[mapCh[NCH[i]]]->GetMean();
    RMSPed[mapCh[NCH[i]]]= histo[mapCh[NCH[i]]]->GetRMS();
    CavasPedestal->cd(1+i);
    histo[mapCh[NCH[i]]]->GetXaxis()->SetTitle("Energy [D.U]");
    histo[mapCh[NCH[i]]]->GetYaxis()->SetTitle("Counts");
    histo[mapCh[NCH[i]]]->Draw();
    //std::cout << "PED NCH" << NCH[i] << "->" << mapCh[NCH[i]] << std::endl;
  }
  
  CavasPedestal->SaveAs((DirData+"/Plot/TR/Pedestal/"+"pedestal"+std::to_string(index)+".png").c_str());
  
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




void GetTdiff(TTree* tree, TH1D** histo,TH2D** tdiffVsE,TH2D** tdiffVsT,TH2D* tdiffVsAmpEff,Double_t* RMSPed, TF1** fitspectrum, Double_t* ped, std::vector<int> NCH, std::vector<std::string> ChType){
  //Histo0 tref-tave, histo1 tref-t1, histo2 tref-t2
  Double_t energy[(int)NCH.size()];
  Double_t chID[(int)NCH.size()];
  Double_t time[(int)NCH.size()];
  Double_t temp1;

  std::map<int,int> mapCh;
  std::map<std::string,int> mapTime;
  
  int h=1;

  tree->SetBranchAddress("energy",energy);
  tree->SetBranchAddress("chId",chID);
  tree->SetBranchAddress("time",time);  
  tree->SetBranchAddress("temp1",&temp1);  

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
    cut[mapCh[NCH[i]]]=Form("energy[mapCh[NCH[%i]]]-ped[mapCh[NCH[%i]]]<fitspectrum[mapCh[NCH[%i]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[%i]]]->GetParameter(2) && energy[mapCh[NCH[%i]]]-ped[mapCh[NCH[%i]]]>fitspectrum[mapCh[NCH[%i]]]->GetParameter(1)-3*fitspectrum[mapCh[NCH[%i]]]->GetParameter(2)",i,i,i,i,i,i,i,i);
    //std::cout << cut[mapCh[NCH[i]]] << std::endl;
    TotalCut+=cut[mapCh[NCH[i]]];
  }
  
  Double_t NSigma=2;
  
  Double_t AmpEff;
  Double_t A1,A2;
  Double_t MeanRMSPed = (RMSPed[mapCh[NCH[mapTime["time1"]]]] + RMSPed[mapCh[NCH[mapTime["time2"]]]])/2;


  for(int i=0;i< tree->GetEntries();i++){
    tree->GetEntry(i);
    //Histo0 tref-tave, histo1 tref-t1, histo2 tref-t2
   
    if( energy[mapCh[NCH[0]]]-ped[mapCh[NCH[0]]]<fitspectrum[mapCh[NCH[0]]]->GetParameter(1)+NSigma*fitspectrum[mapCh[NCH[0]]]->GetParameter(2) && 
	energy[mapCh[NCH[0]]]-ped[mapCh[NCH[0]]]>fitspectrum[mapCh[NCH[0]]]->GetParameter(1)-NSigma*fitspectrum[mapCh[NCH[0]]]->GetParameter(2) && 
        energy[mapCh[NCH[1]]]-ped[mapCh[NCH[1]]]<fitspectrum[mapCh[NCH[1]]]->GetParameter(1)+NSigma*fitspectrum[mapCh[NCH[1]]]->GetParameter(2) && 
        energy[mapCh[NCH[1]]]-ped[mapCh[NCH[1]]]>fitspectrum[mapCh[NCH[1]]]->GetParameter(1)-NSigma*fitspectrum[mapCh[NCH[1]]]->GetParameter(2) && 
        energy[mapCh[NCH[2]]]-ped[mapCh[NCH[2]]]<fitspectrum[mapCh[NCH[2]]]->GetParameter(1)+NSigma*fitspectrum[mapCh[NCH[2]]]->GetParameter(2) && 
        energy[mapCh[NCH[2]]]-ped[mapCh[NCH[2]]]>fitspectrum[mapCh[NCH[2]]]->GetParameter(1)-NSigma*fitspectrum[mapCh[NCH[2]]]->GetParameter(2) ){
      
      histo[0]->Fill(time[mapCh[NCH[mapTime["timeRef"]]]]-(time[mapCh[NCH[mapTime["time1"]]]]+time[mapCh[NCH[mapTime["time2"]]]])/2);
      histo[1]->Fill(time[mapCh[NCH[mapTime["timeRef"]]]]-time[mapCh[NCH[mapTime["time1"]]]]);
      histo[2]->Fill(time[mapCh[NCH[mapTime["timeRef"]]]]-time[mapCh[NCH[mapTime["time2"]]]]);

      tdiffVsT[0]->Fill(temp1,time[mapCh[NCH[mapTime["time1"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]]);
      tdiffVsT[1]->Fill(temp1,time[mapCh[NCH[mapTime["time2"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]]);
      tdiffVsT[2]->Fill(temp1,(time[mapCh[NCH[mapTime["time2"]]]]+time[mapCh[NCH[mapTime["time1"]]]])/2-time[mapCh[NCH[mapTime["timeRef"]]]]);
    
    }//chiudo if grande
    

    if(energy[mapCh[NCH[mapTime["timeRef"]]]]-ped[mapCh[NCH[mapTime["timeRef"]]]] < fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(2) && 
       energy[mapCh[NCH[mapTime["timeRef"]]]]-ped[mapCh[NCH[mapTime["timeRef"]]]] > fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(1)-3*fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]] < fitspectrum[mapCh[NCH[mapTime["time1"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["time1"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]] < fitspectrum[mapCh[NCH[mapTime["time2"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["time2"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]]>10 && energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]]>10 &&
       chID[mapCh[NCH[mapTime["timeRef"]]]]!=-9 && chID[mapCh[NCH[mapTime["time1"]]]]!=-9 && chID[mapCh[NCH[mapTime["time2"]]]]!=-9){
      
      A1=energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]]; 
      A2=energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]];

      AmpEff=A1*A2/(sqrt(A1*A1+A2*A2));
      //AmpEff=A1*A2/(A1+A2);

      tdiffVsAmpEff->Fill(AmpEff/MeanRMSPed,time[mapCh[NCH[mapTime["timeRef"]]]]-(time[mapCh[NCH[mapTime["time1"]]]]+time[mapCh[NCH[mapTime["time2"]]]])/2);
      //tdiffVsAmpEff->Fill(AmpEff,time[mapCh[NCH[mapTime["timeRef"]]]]-(time[mapCh[NCH[mapTime["time1"]]]]+time[mapCh[NCH[mapTime["time2"]]]])/2);

      tdiffVsE[0]->Fill(energy[mapCh[NCH[mapTime["time1"]]]] - ped[mapCh[NCH[mapTime["time1"]]]] , time[mapCh[NCH[mapTime["time1"]]]] - time[mapCh[NCH[mapTime["timeRef"]]]]);
      tdiffVsE[1]->Fill(energy[mapCh[NCH[mapTime["time2"]]]] - ped[mapCh[NCH[mapTime["time2"]]]] , time[mapCh[NCH[mapTime["time2"]]]] - time[mapCh[NCH[mapTime["timeRef"]]]]);

    }//chiudo altro if


  }


  
  std::cout << Form("time1 %i ",  NCH[mapTime["time1"]]) << Form("time2 %i ",     NCH[mapTime["time2"]]) << Form("timeRef %i ", NCH[mapTime["timeRef"]]) << std::endl;
  std::cout << Form("time1 %s ",  ChType[mapTime["time1"]].c_str()) << Form("time2 %s ",    ChType[mapTime["time2"]].c_str()) << Form("timeRef %s ", ChType[mapTime["timeRef"]].c_str()) << std::endl;
}




///////////////////////////////////////////////////////////////////////////////////////////////




void GetTdiffCorr(TTree* tree, TH2D** tdiffVsE,TH2D* tdiffVsAmpEffCorr ,TF1** fitspectrum, std::vector<TF1*> fitcorr, Double_t* ped,Double_t* RMSPed, std::vector<int> NCH, std::vector<std::string> ChType){
  //Histo0 tref-tave, histo1 tref-t1, histo2 tref-t2
  
  Double_t energy[(int)NCH.size()];
  Double_t chID[(int)NCH.size()];
  Double_t time[(int)NCH.size()];
  
 
  std::map<int,int> mapCh;
  std::map<std::string,int> mapTime;
  
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
    cut[mapCh[NCH[i]]]=Form("energy[mapCh[NCH[%i]]]-ped[mapCh[NCH[%i]]]<fitspectrum[mapCh[NCH[%i]]]->GetParameter(1)+NSigma*fitspectrum[mapCh[NCH[%i]]]->GetParameter(2) && energy[mapCh[NCH[%i]]]-ped[mapCh[NCH[%i]]]>fitspectrum[mapCh[NCH[%i]]]->GetParameter(1)-2.5*fitspectrum[mapCh[NCH[%i]]]->GetParameter(2)",i,i,i,i,i,i,i,i);
    //std::cout << cut[mapCh[NCH[i]]] << std::endl;
    TotalCut+=cut[mapCh[NCH[i]]];
  }
  
 
  Double_t AmpEff;
  Double_t A1,A2;
  Double_t MeanRMSPed = (RMSPed[mapCh[NCH[mapTime["time1"]]]] + RMSPed[mapCh[NCH[mapTime["time2"]]]])/2;


  std::cout << fitcorr[0] << "____" << fitcorr[1] << std::endl;
  
  std::cout <<"Address:"<< fitcorr[0] << " "<<" Param:" << fitcorr[0]->GetParameter(0) << " "<< fitcorr[0]->GetParameter(1) << " "<< fitcorr[0]->GetParameter(2) << " "<< fitcorr[0]->GetParameter(3) << " "<< fitcorr[0]->GetParameter(4) << " "<< fitcorr[0]->GetParameter(5) << std::endl;
  
  std::cout <<"Address:"<< fitcorr[1] << " "<<" Param:" << fitcorr[1]->GetParameter(0) << " "<< fitcorr[1]->GetParameter(1) << " "<< fitcorr[1]->GetParameter(2) << " "<< fitcorr[1]->GetParameter(3) << " "<< fitcorr[1]->GetParameter(4) << " "<< fitcorr[1]->GetParameter(5) << std::endl;
  
  for(int i=0;i< tree->GetEntries();i++){
    tree->GetEntry(i);
    
     if(energy[mapCh[NCH[mapTime["timeRef"]]]]-ped[mapCh[NCH[mapTime["timeRef"]]]] < fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(2) && 
       energy[mapCh[NCH[mapTime["timeRef"]]]]-ped[mapCh[NCH[mapTime["timeRef"]]]] > fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(1)-3*fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]] < fitspectrum[mapCh[NCH[mapTime["time1"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["time1"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]] < fitspectrum[mapCh[NCH[mapTime["time2"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["time2"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]] >8 &&
       energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]] >8 &&
       chID[mapCh[NCH[mapTime["timeRef"]]]]!=-9 && chID[mapCh[NCH[mapTime["time1"]]]]!=-9 && chID[mapCh[NCH[mapTime["time2"]]]]!=-9){
      

       A1=energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]]; 
       A2=energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]];
       

       AmpEff=A1*A2/(sqrt(A1*A1+A2*A2));
       

       tdiffVsAmpEffCorr->Fill(AmpEff/MeanRMSPed,time[mapCh[NCH[mapTime["timeRef"]]]]-(time[mapCh[NCH[mapTime["time1"]]]]- fitcorr[0]->Eval(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]])+time[mapCh[NCH[mapTime["time2"]]]]- fitcorr[1]->Eval(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]]))/2);
      

       tdiffVsE[0]->Fill(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]], time[mapCh[NCH[mapTime["time1"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]] - fitcorr[0]->Eval( energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]] ) );

       //std:: cout << "Ch1___:"<< energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]]<<" ->("<< time[mapCh[NCH[mapTime["time1"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]] << ")-(" << fitcorr[0]->Eval( energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]] ) << ")=" << time[mapCh[NCH[mapTime["time1"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]] - fitcorr[0]->Eval( energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]] ) << std::endl;

       
       tdiffVsE[1]->Fill(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]], time[mapCh[NCH[mapTime["time2"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]] - fitcorr[1]->Eval( energy[mapCh[NCH[mapTime["time2"]]]] - ped[mapCh[NCH[mapTime["time2"]]]] ) );

       //std:: cout << "Ch2___"<< energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]]<<" ->("<< time[mapCh[NCH[mapTime["time2"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]] << ")-(" << fitcorr[1]->Eval( energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]] ) << ")=" << time[mapCh[NCH[mapTime["time2"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]] - fitcorr[1]->Eval( energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]] ) << std::endl;
 
     }//chiudo altro if
     
     
  }
  
  
  
  std::cout << Form("time1 %i ",  NCH[mapTime["time1"]]) << Form("time2 %i ",     NCH[mapTime["time2"]]) << Form("timeRef %i ", NCH[mapTime["timeRef"]]) << std::endl;
  std::cout << Form("time1 %s ",  ChType[mapTime["time1"]].c_str()) << Form("time2 %s ",    ChType[mapTime["time2"]].c_str()) << Form("timeRef %s ", ChType[mapTime["timeRef"]].c_str()) << std::endl;
  
}




///////////////////////////////////////////////////////////////////////////////////////////////





TF1* FitCoincSpectrum(TH1D* Profile){
  
  Double_t max;
  Double_t peak1;
  Int_t minRange=0;
  
  for(int i=1 ; i<Profile->GetNbinsX();i++){
    
    if( Profile->GetBinContent(i+1)>Profile->GetBinContent(i) ){  minRange=i; 
      std::cout << Profile->GetBinContent(i+1)<< " > " << Profile->GetBinContent(i) << "  " << minRange << std::endl;} 
    else { break; }
    
  }

  peak1 = Profile->GetBinCenter(Profile->GetMaximumBin());
  max = Profile->GetMaximum();
  Profile->GetXaxis()->UnZoom();

  TF1* spectrum = new TF1(Form("SpectrumFit_%s", Profile->GetName()),"[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) + [3] / (exp( (x*[4]-(2*[1]*[1]/([1]+2*[1])))) + 1)",minRange,peak1+10);

  spectrum->SetParameter(0,max);
  spectrum->SetParameter(1,peak1);
  spectrum->SetParameter(2,3);
  spectrum->SetParameter(3,max/2.1);
  spectrum->SetParameter(4,0.82);
  
  Profile->Fit(Form("SpectrumFit_%s", Profile->GetName()),"R0Q");

  return spectrum;
}



void ProjectionTemp(TH2D** histoTemp,std::string DirData){

  TH1D* projection[3];
  TF1* fitProjection[3];
  Int_t fitStatus[3]={0,0,0};
  TGraphErrors* Plot[3]={new TGraphErrors(),new TGraphErrors(),new TGraphErrors()};
  TCanvas* canvino;
  Int_t h=0,l=0,m=0;

  for(int i=0;i< histoTemp[0]->GetNbinsX();i++){
    
    canvino= new TCanvas(Form("CanvasProjection_%i",i),Form("CanvasProjection_%i",i),1400,700);
    canvino->Divide(3,1);
    
    fitProjection[0]= new TF1(Form("FitProjectionTemp1_%i",i),"gaus");
    fitProjection[1]= new TF1(Form("FitProjectionTemp2_%i",i),"gaus");
    fitProjection[2]= new TF1(Form("FitProjectionTempAve_%i",i),"gaus");
    
    canvino->cd(1);
    projection[0] = histoTemp[0]->ProjectionY(Form("ProjectionTemp1_%i",i),i,i);
    projection[0]->SetTitle(Form("ProjectionTemp1_%i",i));
    projection[0]->GetXaxis()->SetTitle("Time1-time_ref [ps]");
    projection[0]->GetYaxis()->SetTitle("Counts");
    fitProjection[0]->SetRange(projection[0]->GetMean()-1.5*projection[0]->GetRMS(),projection[0]->GetMean()+1.5*projection[0]->GetRMS());
    //fitStatus[0]=projection[0]->Fit(fitProjection[0],"R");
    fitStatus[0]=projection[0]->Fit(fitProjection[0],"RM");
    
    canvino->cd(2);
    projection[1] = histoTemp[1]->ProjectionY(Form("ProjectionTemp2_%i",i),i,i);
    projection[1]->SetTitle(Form("ProjectionTemp2_%i",i));
    projection[1]->GetXaxis()->SetTitle("Time2-time_ref [ps]");
    projection[1]->GetYaxis()->SetTitle("Counts");
    fitProjection[1]->SetRange(projection[1]->GetMean()-1.5*projection[1]->GetRMS(),projection[1]->GetMean()+1.5*projection[1]->GetRMS());
    //fitStatus[1]=projection[1]->Fit(fitProjection[1],"R");
    fitStatus[1]=projection[1]->Fit(fitProjection[1],"RM");

    canvino->cd(3);
    projection[2] = histoTemp[2]->ProjectionY(Form("ProjectionTempAve_%i",i),i,i);
    projection[2]->SetTitle(Form("ProjectionTempAve_%i",i));
    projection[2]->GetXaxis()->SetTitle("T_{ave}-time_ref [ps]");
    projection[2]->GetYaxis()->SetTitle("Counts");
    fitProjection[2]->SetRange(projection[2]->GetMean()-1.5*projection[2]->GetRMS(),projection[2]->GetMean()+1.5*projection[2]->GetRMS());
    //fitStatus[2]=projection[2]->Fit(fitProjection[2],"R");
    fitStatus[2]=projection[2]->Fit(fitProjection[2],"R");
    
    std::cout << " fitstaus1: "<< fitStatus[0] << std::endl;
    if( fitStatus[0]!=-1 && fitProjection[0]->GetParError(2)<40 ){       
      Plot[0]->SetPoint(h,histoTemp[0]->GetXaxis()->GetBinCenter(i),fitProjection[0]->GetParameter(2));
      Plot[0]->SetPointError(h,0,fitProjection[0]->GetParError(2));
      h++;
      std::cout << "filled1" << std::endl;
    }
    
    std::cout << " fitstaus2: "<< fitStatus[1] << std::endl;
    if( fitStatus[1]!=-1 && fitProjection[1]->GetParError(2)<40 ){      
      Plot[1]->SetPoint(l,histoTemp[1]->GetXaxis()->GetBinCenter(i),fitProjection[1]->GetParameter(2));
      Plot[1]->SetPointError(l,0,fitProjection[1]->GetParError(2));
      l++;
      std::cout << "filled2" << std::endl;
    }
    
    std::cout << " fitstausave: "<< fitStatus[3] << std::endl;
    if( fitStatus[1]!=-1 ){ 
      Plot[2]->SetPoint(m,histoTemp[2]->GetXaxis()->GetBinCenter(i),fitProjection[2]->GetParameter(2));
      Plot[2]->SetPointError(m,0,fitProjection[2]->GetParError(2));
      m++;
      std::cout << "filledAve" << std::endl;
    }
    
    
    canvino->SaveAs((DirData+"/Plot/TR/ProjectionTemp/"+"Projection"+std::to_string(i)+".png").c_str());
    
 }

  std::cout << "PlotAddress :" << Plot[0] << " " << Plot[1] << " " << Plot[2] << std::endl;

  TCanvas* canvasPlot= new TCanvas("canvasPlot","canvasPlot",1400,700);
  canvasPlot->Divide(3,1);
  
  canvasPlot->cd(1);
  Plot[0]->SetTitle("TimeResolution1VsTempBox");
  Plot[0]->SetName(Plot[0]->GetTitle());
  Plot[0]->GetXaxis()->SetTitle("Temp [°C]");
  Plot[0]->GetYaxis()->SetTitle("Time Resolution [ps]");
  Plot[0]->Draw("AP");
  
  canvasPlot->cd(2);
  Plot[1]->SetTitle("TimeResolution2VsTempBox");
  Plot[1]->SetName(Plot[1]->GetTitle());
  Plot[1]->GetXaxis()->SetTitle("Temp [°C]");
  Plot[1]->GetYaxis()->SetTitle("Time Resolution [ps]");
  Plot[1]->Draw("AP");

  canvasPlot->cd(3);
  Plot[2]->SetTitle("TimeResolutionAveVsTempBox");
  Plot[2]->SetName(Plot[2]->GetTitle());
  Plot[2]->GetXaxis()->SetTitle("Temp [°C]");
  Plot[2]->GetYaxis()->SetTitle("Time Resolution [ps]");
  Plot[2]->Draw("AP");

  canvasPlot->SaveAs((DirData+"/Plot/TR/TResVsTemp.png").c_str()); 


  TFile* f1 = new TFile((DirData+"/Plot/TReVsTemp.root").c_str(),"RECREATE");
  f1->cd();
  
  Plot[0]->Write();
  Plot[1]->Write();
  Plot[2]->Write();

  f1->Save();
  f1->Close();


}
