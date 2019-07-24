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
#include "TLegend.h"
       
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
Double_t GetTdiffCorr(TTree* ,TH2D**,TH2D* ,TH2D*, TF1** ,std::vector<TF1*>, Double_t*,Double_t* , std::vector<int> , std::vector<std::string> );
void ProjectionTemp(TH2D**,std::string);
void GetTresVsAEff(TH2D* ,Double_t, std::string, Double_t*, Double_t*, Double_t*, Double_t*);
TGraphErrors* GetTimeResVsE(TH2D* , std::string ,std::string );

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
  gSystem->Exec(("mkdir "+DirData+"/Plot/TR/TdiffVsAmpEffTemp").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TR/TimeResVsAEffTemp").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TR/ProjTResE1").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TR/ProjTResE2").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TR/ProjTResEMean").c_str());
  

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
   NFilePhys=2;
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

  tdiffVsE[0] = new TH2D("tdiffVsE_time1","tdiffVsE_time1",150,0,150,200,-3000,3000);
  tdiffVsE[1] = new TH2D("tdiffVsE_time2","tdiffVsE_time2",150,0,150,200,-3000,3000);

  TH2D* tdiffVsT[3];

  tdiffVsT[0] = new TH2D("tdiffVsT_time1","tdiffVsT_time1",40,21,29,200,-2000,2000);
  tdiffVsT[1] = new TH2D("tdiffVsT_time2","tdiffVsT_time2",40,21,29,200,-2000,2000);
  tdiffVsT[2] = new TH2D("tdiffVsT_tave","tdiffVsT_time_ave",40,21,29,200,-2000,2000);
  
  TH2D* tdiffVsAmpEff = new TH2D("tdiffVsAmpEff","tdiffVsAmpEff",100,-10,90,100,-3500,3000);
  
  TLine* line[(int)Channels.size()][2];
    
  for(int i=0;i<NFilePhys;i++){
   
    
    f0= TFile::Open((DirData+"/"+FileListPhysics.at(i)).c_str());
    tree0 = (TTree*)f0->Get("data");
    
    for(int j=0; j< (int)Channels.size();j++){
      
      Spectrum[i][mapCh[Channels[j]]] = new TH1D(Form("Ch%i_%i_%s",Channels[j],i,(ChannelType[mapCh[Channels[j]]]).c_str()),Form("Ch%i_%i_%s",Channels[j],i,(ChannelType[mapCh[Channels[j]]]).c_str()),150,0,150);
      
    }//CHIUDO FOR J
    
    GetSpectrum(tree0,Spectrum[i],Pedestal[i],Channels);
    
    for(int k=0; k<(int)Channels.size(); k++){
      FitSpectrum[i][mapCh[Channels[k]]]=FitCoincSpectrum(Spectrum[i][mapCh[Channels[k]]]);
    }
    
    canvino = new TCanvas("canvino","canvino",1500,700);
    canvino->Divide((int)Channels.size(),1);
    
    for(int k=0; k<(int)Channels.size(); k++){
      
      Double_t min=FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(1)-2*FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(2);
      Double_t max=FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(1)+2*FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(2);
      
      line[mapCh[Channels[k]]][0]= new TLine(min,0,min,Spectrum[i][mapCh[Channels[k]]]->GetMaximum());
      line[mapCh[Channels[k]]][1]= new TLine(max,0,max,Spectrum[i][mapCh[Channels[k]]]->GetMaximum());
      
      canvino->cd(k+1);
      //Spectrum[i][mapCh[Channels[k]]]->GetXaxis()->SetRange(1,100);
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

  TH2D* tdiffCorrVsE[3];

  tdiffCorrVsE[0] = new TH2D("tdiffCorrVsE_time1","tdiffCorrVsE_time1",150,0,150,200,-3000,3000);
  tdiffCorrVsE[1] = new TH2D("tdiffCorrVsE_time2","tdiffCorrVsE_time2",150,0,150,200,-3000,3000);
  tdiffCorrVsE[2] = new TH2D("tdiffCorrVsEMean_time_ave","tdiffCorrVsEMean_time_ave",150,0,150,200,-3000,3000);

  TH2D* tdiffCorrVsAmpEff = new TH2D("tdiffCorrVsAmpEff","tdiffCorrVsAmpEff",100,-10,90,150,-3500,3000);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  std::vector<TH2D*> tdiffCorrVsAEffTemp(NFilePhys);
  Double_t tempRun[NFilePhys];
  TCanvas* CanvasTdiffAEffTemp = new TCanvas("CanvasTdiffAEffTemp","CanvasTdiffAEffTemp",700,700);

  Double_t N[NFilePhys], sN[NFilePhys], C[NFilePhys], sC[NFilePhys];
  
  for(int i=0;i<NFilePhys;i++){
    
    f0= TFile::Open((DirData+"/"+FileListPhysics.at(i)).c_str());
    tree0 = (TTree*)f0->Get("data"); 
    
    //std::cout<<"before function " << fitTdiffVsE[0] << " " << fitTdiffVsE[1] << std::endl;
    tdiffCorrVsAEffTemp[i]= new TH2D(Form("tdiffCorrVsAmpEff%i",i),Form("tdiffCorrVsAmpEff%i",i),100,-10,90,150,-3500,3000);
    
    tempRun[i] = GetTdiffCorr(tree0,tdiffCorrVsE,tdiffCorrVsAmpEff,tdiffCorrVsAEffTemp[i],FitSpectrum[i],fitTdiffVsE,Pedestal[i],RMSPedestal[i],Channels,ChannelType);
        
    CanvasTdiffAEffTemp = new TCanvas("CanvasTdiffAEffTemp","CanvasTdiffAEffTemp",700,700);
    tdiffCorrVsAEffTemp[i]->GetXaxis()->SetTitle("#A_{eff}/#\sigma_{n}");
    tdiffCorrVsAEffTemp[i]->GetYaxis()->SetTitle("t_ave-t_ref [ps]");
    tdiffCorrVsAEffTemp[i]->SetTitle((((std::string)tdiffCorrVsAEffTemp[i]->GetTitle())+"_"+std::to_string(tempRun[i])).c_str());
    tdiffCorrVsAEffTemp[i]->Draw("COLZ");
    
    CanvasTdiffAEffTemp->SaveAs((DirData+"/Plot/TR/TdiffVsAmpEffTemp/Plot"+std::to_string(i)+".png").c_str());
    
    GetTresVsAEff(tdiffCorrVsAEffTemp[i],tempRun[i],DirData,&N[i],&sN[i],&C[i],&sC[i]);
    
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
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

  TF1* fitTdiffCorrVsE[3];
  fitTdiffCorrVsE[0] = new TF1("fitTdiffCorrVsE0","pol5");
  fitTdiffCorrVsE[1] = new TF1("fitTdiffCorrVsE1","pol5");
  fitTdiffCorrVsE[2] = new TF1("fitTdiffCorrVsEMean","pol5");

  TCanvas* CanvasTdiffCorrVsE = new TCanvas("CanvasTdiffCorrVsE","CanvasTdiffCorrVsE",1500,700);
  
  CanvasTdiffCorrVsE->Divide(3,1);
  
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

  CanvasTdiffCorrVsE->cd(3);
  tdiffCorrVsE[2]->GetXaxis()->SetTitle("MeanEnergy [D.U.]"); 
  tdiffCorrVsE[2]->GetYaxis()->SetTitle("timeave-timeRef [ps]"); 
  tdiffCorrVsE[2]->Draw("COLZ");  
  tdiffCorrVsE[2]->Fit(fitTdiffCorrVsE[2]);

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

    if(fitStatus==0 && tdiffCorrVsAmpEff->GetXaxis()->GetBinCenter(i)>0 && Projection->GetEntries()>150){    
     
      CanvasProjection->SaveAs((DirData+"/Plot/TR/ProjectionAEff/Projection"+std::to_string(i)+".png").c_str());
     
      AEffValue.push_back(tdiffCorrVsAmpEff->GetXaxis()->GetBinCenter(i));
      TimeResValue.push_back(sqrt(FitProjection->GetParameter(2)*FitProjection->GetParameter(2)-91.9*91.9));
      ErrTimeResValue.push_back(FitProjection->GetParError(2));
    
    }
    
    
    delete Projection;
    delete CanvasProjection;

    std::cout << tdiffCorrVsAmpEff->GetXaxis()->GetBinCenter(i) << " " << FitProjection->GetParameter(2) << " " << FitProjection->GetParError(2) << std::endl;
  }
  
  TCanvas* CanvasTimeResVsAEff = new TCanvas("CanvasTimeResVsAEff","CanvasTimeResVsAEff",900,1000);
  
  
  TGraphErrors* TimeResVsAEff = new TGraphErrors((int)AEffValue.size(),&AEffValue[0],&TimeResValue[0],0,&ErrTimeResValue[0]);
  
  TF1* fitAEff = new TF1("fitAEff", "sqrt([0]*[0]/(x*x)+2*[1]*[1])");
  
  fitAEff->SetParameter(0,1600);
  fitAEff->SetParameter(1,38);
  fitAEff->SetParLimits(0,500,2500);
  fitAEff->SetParLimits(1,0,500);
  
  TimeResVsAEff->SetTitle("TimeResVsAEff");
  TimeResVsAEff->SetName(TimeResVsAEff->GetTitle());
  TimeResVsAEff->GetXaxis()->SetTitle("A_Eff/#sigma_{n} [D.U]");
  TimeResVsAEff->GetYaxis()->SetTitle("TimeRes [ps]");
  
  TimeResVsAEff->SetMarkerStyle(8);
  
  TimeResVsAEff->GetXaxis()->SetLimits(0,30);
  TimeResVsAEff->GetXaxis()->SetRangeUser(0,30);
  TimeResVsAEff->GetYaxis()->SetRangeUser(0,400);

  TimeResVsAEff->Draw("AP");
  TimeResVsAEff->Fit(fitAEff,"WEM");
  
  
  CanvasTimeResVsAEff->SaveAs((DirData+"/Plot/TR/TimeResVsAEff.png").c_str()); 
  
  TGraphErrors* NoiseTermVsTemp = new TGraphErrors(NFilePhys,tempRun,N,0,sN);
  NoiseTermVsTemp->SetTitle("NoiseTermVsTemp");
  NoiseTermVsTemp->SetName(NoiseTermVsTemp->GetTitle());
  NoiseTermVsTemp->GetXaxis()->SetTitle("TempRun [°C]");
  NoiseTermVsTemp->GetXaxis()->SetTitle("N");

  Double_t MeanNoise,RMSNoise;
  MeanNoise=TMath::Mean(NFilePhys,N);
  RMSNoise = TMath::RMS(NFilePhys,N);
  
  TLegend* legNoise = new TLegend();
  legNoise->SetHeader(Form("Mean: %.1lf \n RMS: %.1lf",MeanNoise,RMSNoise));
  
  TGraphErrors* ConstTermVsTemp = new TGraphErrors(NFilePhys,tempRun,C,0,sC);
  ConstTermVsTemp->SetTitle("ConstTermVsTemp");
  ConstTermVsTemp->SetName(ConstTermVsTemp->GetTitle());
  ConstTermVsTemp->GetXaxis()->SetTitle("TempRun [°C]");
  ConstTermVsTemp->GetXaxis()->SetTitle("C");
  
  Double_t MeanConst,RMSConst;
  MeanConst=TMath::Mean(NFilePhys,C);
  RMSConst = TMath::RMS(NFilePhys,C);
  
  TLegend* legConst = new TLegend();
  legConst->SetHeader(Form("Mean: %.1lf \n RMS: %.1lf",MeanConst,RMSConst));


  TCanvas* CanvasConstNoise = new TCanvas("CanvasConstNoise","CanvasConstNoise",1400,700);
  CanvasConstNoise->Divide(2,1);
  
  CanvasConstNoise->cd(1);
  NoiseTermVsTemp->Draw("AP");
  legNoise->Draw("SAME");
  
  CanvasConstNoise->cd(2);
  ConstTermVsTemp->Draw("AP");
  legConst->Draw("SAME");  

  CanvasConstNoise->SaveAs((DirData+"/Plot/TR/NCVsTemp.png").c_str());

  TH1D* CProjection = new TH1D("CProjection","CProjection",40,0,100);
  TH1D* NProjection = new TH1D("NProjection","NProjection",40,0,2000);
  
  for(int i=0;i< ConstTermVsTemp->GetN();i++){
    CProjection->Fill(*(ConstTermVsTemp->GetY()+i));
  }
  
  for(int i=0;i< NoiseTermVsTemp->GetN();i++){
    NProjection->Fill(*(NoiseTermVsTemp->GetY()+i));
  }
  
  TGraphErrors* TimeResVsE1;
  TGraphErrors* TimeResVsE2;
  TGraphErrors* TimeResVsEMean;
  
  TimeResVsE1= GetTimeResVsE(tdiffCorrVsE[0],"1",DirData);
  TimeResVsE2= GetTimeResVsE(tdiffCorrVsE[1],"2",DirData);
  TimeResVsEMean= GetTimeResVsE(tdiffCorrVsE[1],"Mean",DirData);
  
  TF1* fitTResVsECh1 = new TF1("fitTResVsECh1","sqrt( [0]*[0]/(x*x) + [1]*[1]/(x) + [2]*[2] )");
  TF1* fitTResVsECh2 = new TF1("fitTResVsECh2","sqrt( [0]*[0]/(x*x) + [1]*[1]/(x) + [2]*[2] )");
  TF1* fitTResVsEMean = new TF1("fitTResVsEMean","sqrt( [0]*[0]/(x*x) + [1]*[1]/(x) + [2]*[2] )");

  TCanvas* CanvTRVsE = new TCanvas("CanvTRVsE","CanvTRVsE",1400,700);
  CanvTRVsE->Divide(2,1);
  
  CanvTRVsE->cd(1);
  TimeResVsE1->Draw("AP");
  fitTResVsECh1->SetParameter(0,1000);
  fitTResVsECh1->SetParameter(1,0.1);
  fitTResVsECh1->SetParameter(2,40);
  TimeResVsE1->Fit(fitTResVsECh1);
  
  CanvTRVsE->cd(2);
  TimeResVsE2->Draw("AP");
  fitTResVsECh1->SetParameter(0,1000);
  fitTResVsECh1->SetParameter(1,0.1);
  fitTResVsECh1->SetParameter(2,40);
  TimeResVsE2->Fit(fitTResVsECh2);
  
  CanvTRVsE->SaveAs((DirData+"/Plot/TR/TResVsE.png").c_str());
  

  TCanvas* CanvasTRVsEMean = new TCanvas("CanvasTRVsEMean","CanvasTRVsEMean",800,800);
  TimeResVsEMean->Draw("AP");
  fitTResVsEMean->SetParameter(0,1000);
  fitTResVsEMean->SetParameter(1,400);
  fitTResVsEMean->SetParameter(2,40);
  TimeResVsEMean->Fit(fitTResVsEMean); 
  
  CanvasTRVsEMean->SaveAs((DirData+"/Plot/TR/TResVsEMean.png").c_str()); 
  
  TFile* f1 = new TFile((DirData+"/Plot/TResVsAEff.root").c_str(),"RECREATE");
  f1->cd();
  
  TimeResVsAEff->Write();
  fitAEff->Write();
  NoiseTermVsTemp->Write();
  ConstTermVsTemp->Write();
  tdiffCorrVsE[0]->Write();
  fitTdiffCorrVsE[0]->Write();
  tdiffCorrVsE[1]->Write();
  fitTdiffCorrVsE[1]->Write();
  tdiffVsE[0]->Write();
  fitTdiffVsE[0]->Write();
  tdiffVsE[1]->Write();
  fitTdiffVsE[1]->Write();
  tdiffVsAmpEff->Write();
  tdiffCorrVsAmpEff->Write();
  NoiseTermVsTemp->Write();
  ConstTermVsTemp->Write();
  CProjection->Write();
  NProjection->Write();
  TimeResVsE1->Write();
  TimeResVsE2->Write();
  tdiff[0]->SetName("TResAve");
  tdiff[0]->Write();
  tdiff[1]->SetName("TRes1"); 
  tdiff[1]->Write();
  tdiff[2]->SetName("TRes1"); 
  tdiff[2]->Write();
  tdiffCorrVsE[2]->Write();
  TimeResVsEMean->Write();
   

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
    
    
    ped[mapCh[NCH[i]]]= histo[mapCh[NCH[i]]]->GetMean(); //rimuovi quando il canale funziona
    
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
  Double_t tempSiPMTest;

  std::map<int,int> mapCh;
  std::map<std::string,int> mapTime;
  
  int h=1;

  tree->SetBranchAddress("energy",energy);
  tree->SetBranchAddress("chId",chID);
  tree->SetBranchAddress("time",time);  
  tree->SetBranchAddress("tempSiPMTest",&tempSiPMTest);  

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

      tdiffVsT[0]->Fill(tempSiPMTest,time[mapCh[NCH[mapTime["time1"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]]);
      tdiffVsT[1]->Fill(tempSiPMTest,time[mapCh[NCH[mapTime["time2"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]]);
      tdiffVsT[2]->Fill(tempSiPMTest,(time[mapCh[NCH[mapTime["time2"]]]]+time[mapCh[NCH[mapTime["time1"]]]])/2-time[mapCh[NCH[mapTime["timeRef"]]]]);
    
    }//chiudo if grande

    

    if(energy[mapCh[NCH[mapTime["timeRef"]]]]-ped[mapCh[NCH[mapTime["timeRef"]]]] < fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(1)+1*fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(2) && 
       energy[mapCh[NCH[mapTime["timeRef"]]]]-ped[mapCh[NCH[mapTime["timeRef"]]]] > fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(1)-1*fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(2) &&
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




Double_t GetTdiffCorr(TTree* tree, TH2D** tdiffVsE,TH2D* tdiffVsAmpEffCorr,TH2D* tdiffVsAmpEffTemp ,TF1** fitspectrum, std::vector<TF1*> fitcorr, Double_t* ped,Double_t* RMSPed, std::vector<int> NCH, std::vector<std::string> ChType){
  //Histo0 tref-tave, histo1 tref-t1, histo2 tref-t2
  
  Double_t energy[(int)NCH.size()];
  Double_t chID[(int)NCH.size()];
  Double_t time[(int)NCH.size()];
  Double_t tempSiPMTest;
 
  Double_t TMean=0;
  
  std::map<int,int> mapCh;
  std::map<std::string,int> mapTime;
  
  int h=1;

  tree->SetBranchAddress("energy",energy);
  tree->SetBranchAddress("chId",chID);
  tree->SetBranchAddress("time",time);  
  tree->SetBranchAddress("tempSiPMTest",&tempSiPMTest);  
  
  
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
  Double_t EMean;

  std::cout << fitcorr[0] << "____" << fitcorr[1] << std::endl;
  
  std::cout <<"Address:"<< fitcorr[0] << " "<<" Param:" << fitcorr[0]->GetParameter(0) << " "<< fitcorr[0]->GetParameter(1) << " "<< fitcorr[0]->GetParameter(2) << " "<< fitcorr[0]->GetParameter(3) << " "<< fitcorr[0]->GetParameter(4) << " "<< fitcorr[0]->GetParameter(5) << std::endl;
  
  std::cout <<"Address:"<< fitcorr[1] << " "<<" Param:" << fitcorr[1]->GetParameter(0) << " "<< fitcorr[1]->GetParameter(1) << " "<< fitcorr[1]->GetParameter(2) << " "<< fitcorr[1]->GetParameter(3) << " "<< fitcorr[1]->GetParameter(4) << " "<< fitcorr[1]->GetParameter(5) << std::endl;
  
  for(int i=0;i< tree->GetEntries();i++){
    tree->GetEntry(i);
    
    TMean+=tempSiPMTest;

     if(energy[mapCh[NCH[mapTime["timeRef"]]]]-ped[mapCh[NCH[mapTime["timeRef"]]]] < fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(1)+1*fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(2) && 
       energy[mapCh[NCH[mapTime["timeRef"]]]]-ped[mapCh[NCH[mapTime["timeRef"]]]] > fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(1)-1*fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]] < fitspectrum[mapCh[NCH[mapTime["time1"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["time1"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]] < fitspectrum[mapCh[NCH[mapTime["time2"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["time2"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]] >8 &&
       energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]] >8 &&
       chID[mapCh[NCH[mapTime["timeRef"]]]]!=-9 && chID[mapCh[NCH[mapTime["time1"]]]]!=-9 && chID[mapCh[NCH[mapTime["time2"]]]]!=-9){
      

       A1=energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]]; 
       A2=energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]];
       

       AmpEff=A1*A2/(sqrt(A1*A1+A2*A2));
       
       EMean= (A1+A2)/2;

       tdiffVsAmpEffCorr->Fill(AmpEff/MeanRMSPed,time[mapCh[NCH[mapTime["timeRef"]]]]-(time[mapCh[NCH[mapTime["time1"]]]]- fitcorr[0]->Eval(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]])+time[mapCh[NCH[mapTime["time2"]]]]- fitcorr[1]->Eval(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]]))/2);
      
       tdiffVsAmpEffTemp->Fill(AmpEff/MeanRMSPed,time[mapCh[NCH[mapTime["timeRef"]]]]-(time[mapCh[NCH[mapTime["time1"]]]]- fitcorr[0]->Eval(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]])+time[mapCh[NCH[mapTime["time2"]]]]- fitcorr[1]->Eval(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]]))/2);




       tdiffVsE[0]->Fill(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]], time[mapCh[NCH[mapTime["time1"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]] - fitcorr[0]->Eval( energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]] ) );
              
       tdiffVsE[1]->Fill(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]], time[mapCh[NCH[mapTime["time2"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]] - fitcorr[1]->Eval( energy[mapCh[NCH[mapTime["time2"]]]] - ped[mapCh[NCH[mapTime["time2"]]]] ) );

       tdiffVsE[2]->Fill(EMean,time[mapCh[NCH[mapTime["timeRef"]]]]-(time[mapCh[NCH[mapTime["time1"]]]]- fitcorr[0]->Eval(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]])+time[mapCh[NCH[mapTime["time2"]]]]- fitcorr[1]->Eval(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]]))/2);
        
     }//chiudo altro if
       
  }// chiudo for
    
  std::cout << Form("time1 %i ",  NCH[mapTime["time1"]]) << Form("time2 %i ",     NCH[mapTime["time2"]]) << Form("timeRef %i ", NCH[mapTime["timeRef"]]) << std::endl;
  std::cout << Form("time1 %s ",  ChType[mapTime["time1"]].c_str()) << Form("time2 %s ",    ChType[mapTime["time2"]].c_str()) << Form("timeRef %s ", ChType[mapTime["timeRef"]].c_str()) << std::endl;
  
  TMean/=tree->GetEntries();

  return TMean;
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
    
    fitProjection[0]= new TF1(Form("FitProjectionTempSiPMRef_%i",i),"gaus");
    fitProjection[1]= new TF1(Form("FitProjectionTemp2_%i",i),"gaus");
    fitProjection[2]= new TF1(Form("FitProjectionTempAve_%i",i),"gaus");
    
    canvino->cd(1);
    projection[0] = histoTemp[0]->ProjectionY(Form("ProjectionTempSiPMRef_%i",i),i,i);
    projection[0]->SetTitle(Form("ProjectionTempSiPMRef_%i",i));
    projection[0]->GetXaxis()->SetTitle("Time1-time_ref [ps]");
    projection[0]->GetYaxis()->SetTitle("Counts");
    fitProjection[0]->SetRange(projection[0]->GetMean()-2*projection[0]->GetRMS(),projection[0]->GetMean()+2*projection[0]->GetRMS());
    //fitStatus[0]=projection[0]->Fit(fitProjection[0],"R");
    fitStatus[0]=projection[0]->Fit(fitProjection[0],"RM");
    
    canvino->cd(2);
    projection[1] = histoTemp[1]->ProjectionY(Form("ProjectionTemp2_%i",i),i,i);
    projection[1]->SetTitle(Form("ProjectionTemp2_%i",i));
    projection[1]->GetXaxis()->SetTitle("Time2-time_ref [ps]");
    projection[1]->GetYaxis()->SetTitle("Counts");
    fitProjection[1]->SetRange(projection[1]->GetMean()-2*projection[1]->GetRMS(),projection[1]->GetMean()+2*projection[1]->GetRMS());
    //fitStatus[1]=projection[1]->Fit(fitProjection[1],"R");
    fitStatus[1]=projection[1]->Fit(fitProjection[1],"RM");

    canvino->cd(3);
    projection[2] = histoTemp[2]->ProjectionY(Form("ProjectionTempAve_%i",i),i,i);
    projection[2]->SetTitle(Form("ProjectionTempAve_%i",i));
    projection[2]->GetXaxis()->SetTitle("T_{ave}-time_ref [ps]");
    projection[2]->GetYaxis()->SetTitle("Counts");
    fitProjection[2]->SetRange(projection[2]->GetMean()-2*projection[2]->GetRMS(),projection[2]->GetMean()+2*projection[2]->GetRMS());
    //fitStatus[2]=projection[2]->Fit(fitProjection[2],"R");
    fitStatus[2]=projection[2]->Fit(fitProjection[2],"R");
    
   

    Double_t StDev1,StDev2,StDevAve;
    Double_t EStDev1,EStDev2,EStDevAve;
    Double_t PixelRes = 91.9;
    Double_t SPiexelRes = 0.1;

    std::cout << " fitstaus1: "<< fitStatus[0] << std::endl;
    if( fitStatus[0]!=-1 && fitProjection[0]->GetParError(2)<40 && projection[0]->GetEntries()>150 ){       
      StDev1=fitProjection[0]->GetParameter(2);
      EStDev1=fitProjection[0]->GetParError(2);

      Plot[0]->SetPoint(h,histoTemp[0]->GetXaxis()->GetBinCenter(i),sqrt(StDev1*StDev1-PixelRes*PixelRes));
      Plot[0]->SetPointError(h,0,sqrt( StDev2*StDev1/(StDev1*StDev1-PixelRes*PixelRes)*EStDev1*EStDev1 + PixelRes*PixelRes/(StDev1*StDev1-PixelRes*PixelRes)*SPiexelRes*SPiexelRes ));
      h++;
      std::cout << "filled1" << std::endl;
    }
    
    std::cout << " fitstaus2: "<< fitStatus[1] << std::endl;
    if( fitStatus[1]!=-1 && fitProjection[1]->GetParError(2)<40 && projection[1]->GetEntries()>150){
      StDev2 = fitProjection[1]->GetParameter(2);
      EStDev2 = fitProjection[1]->GetParError(2);
 
      Plot[1]->SetPoint(l,histoTemp[1]->GetXaxis()->GetBinCenter(i),sqrt(StDev2*StDev2-PixelRes*PixelRes));
      Plot[1]->SetPointError(l,0, sqrt( StDev2*StDev2/(StDev2*StDev2-PixelRes*PixelRes)*EStDev2*EStDev2 + PixelRes*PixelRes/(StDev2*StDev2-PixelRes*PixelRes)*SPiexelRes*SPiexelRes ) );
      l++;
      std::cout << "filled2" << std::endl;
    }
    
    std::cout << " fitstausave: "<< fitStatus[3] << std::endl;
    if( fitStatus[1]!=-1 && projection[2]->GetEntries()>150){ 
      StDevAve=fitProjection[2]->GetParameter(2);
      EStDevAve=fitProjection[2]->GetParError(2);
      
     
      Plot[2]->SetPoint(m,histoTemp[2]->GetXaxis()->GetBinCenter(i),sqrt(StDevAve*StDevAve-PixelRes*PixelRes));
      Plot[2]->SetPointError(m,0,sqrt( StDevAve*StDevAve/(StDevAve*StDevAve-PixelRes*PixelRes)*EStDevAve*EStDevAve + PixelRes*PixelRes/(StDevAve*StDevAve-PixelRes*PixelRes)*SPiexelRes*SPiexelRes ) );
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




void GetTresVsAEff(TH2D* Histo,Double_t temp,std::string DirData,Double_t* N,Double_t* sN,Double_t* C,Double_t* sC){

  gSystem->Exec(("mkdir "+DirData+"/Plot/TR/TimeResVsAEffTemp/TimeResVsAEffTemp"+std::to_string(temp*1000).erase(12-7,12)).c_str());

  Int_t nbins = Histo->GetXaxis()->GetNbins();
  
  TH1D* Projection;
  TF1* FitProjection;
  TCanvas* CanvasProjection;

  std::vector<Double_t> AEffValue;
  std::vector<Double_t> TimeResValue;
  std::vector<Double_t> ErrTimeResValue;

  Int_t fitStatus;

  for(int i=0; i<nbins-1; i+=2){


    CanvasProjection = new TCanvas("CanvasProjection","CanvasProjection",900,1000);
    
    Projection = Histo->ProjectionY(Form("ProjectionCorr%i_Amp%lf",i,Histo->GetXaxis()->GetBinCenter(i)),i,i);
    Projection->GetXaxis()->SetTitle("t_diff [ps]");
    Projection->GetYaxis()->SetTitle("Counts");
    Projection->SetTitle(Projection->GetName());
    Projection->Draw();
    
    FitProjection = new TF1(Form("fitProjectionCorr%i",i) , "gaus" , Projection->GetBinCenter(Projection->GetMaximumBin())-1.1*Projection->GetRMS(),Projection->GetBinCenter(Projection->GetMaximumBin())+1.1*Projection->GetRMS());
    
    std::cout << "Range_______:    " << Projection->GetBinCenter(Projection->GetMaximumBin())-1.1*Projection->GetRMS() << "   " << Projection->GetBinCenter(Projection->GetMaximumBin())+1.1*Projection->GetRMS() << std::endl;
          
    FitProjection->SetParameter(0 , Projection->GetMaximum());
    FitProjection->SetParameter(1 , Projection->GetMaximumBin());
    
    fitStatus = Projection->Fit(FitProjection,"RQ");
    
    
    
    if(fitStatus==0 && Histo->GetXaxis()->GetBinCenter(i)>0 && Projection->GetEntries()>150){
      
      CanvasProjection->SaveAs((DirData+"/Plot/TR/TimeResVsAEffTemp/TimeResVsAEffTemp"+std::to_string(temp*1000).erase(12-7,12)+"/Proj"+std::to_string(i)+".png").c_str());
      
      AEffValue.push_back(Histo->GetXaxis()->GetBinUpEdge(i));
      TimeResValue.push_back(sqrt(FitProjection->GetParameter(2)*FitProjection->GetParameter(2)-91.9*91.9));
      ErrTimeResValue.push_back(FitProjection->GetParError(2));
    }
    
    
    delete Projection;
    delete CanvasProjection;

    std::cout << Histo->GetXaxis()->GetBinCenter(i) << " " << FitProjection->GetParameter(2) << " " << FitProjection->GetParError(2) << std::endl;
  }
  
  TCanvas* CanvasTimeResVsAEff = new TCanvas("CanvasTimeResVsAEff","CanvasTimeResVsAEff",900,1000);
  
  
  TGraphErrors* TimeResVsAEff = new TGraphErrors((int)AEffValue.size(),&AEffValue[0],&TimeResValue[0],0,&ErrTimeResValue[0]);
  
  TF1* fitAEff = new TF1("fitAEff", "sqrt([0]*[0]/(x*x)+2*[1]*[1])");
  
  fitAEff->SetParameter(0,1600);
  fitAEff->SetParameter(1,60);
  fitAEff->SetParLimits(1,1000,2500);
  fitAEff->SetParLimits(1,0,500);
  
  TimeResVsAEff->SetTitle("TimeResVsAEff");
  TimeResVsAEff->SetName(TimeResVsAEff->GetTitle());
  TimeResVsAEff->GetXaxis()->SetTitle("A_Eff/#sigma_{n} [D.U]");
  TimeResVsAEff->GetYaxis()->SetTitle("TimeRes [ps]");
  
  TimeResVsAEff->SetMarkerStyle(8);
  
  TimeResVsAEff->GetXaxis()->SetLimits(0,TMath::MaxElement(TimeResVsAEff->GetN(),TimeResVsAEff->GetX()));
  TimeResVsAEff->GetXaxis()->SetRangeUser(0,TMath::MaxElement(TimeResVsAEff->GetN(),TimeResVsAEff->GetX()));
  TimeResVsAEff->GetYaxis()->SetRangeUser(0,400);
  
  TimeResVsAEff->Draw("AP");
  Int_t FitStatus =  TimeResVsAEff->Fit(fitAEff,"WEM");
  
  
  CanvasTimeResVsAEff->SaveAs((DirData+"/Plot/TR/TimeResVsAEffTemp/TimeResVsAEff"+std::to_string(temp)+".png").c_str()); 

  // if( fitAEff->GetParError(0) < 2000 && fitAEff->GetParError(1) < 2000){

    *N=fitAEff->GetParameter(0);
    *sN=fitAEff->GetParError(0);
    *C=fitAEff->GetParameter(1);
    *sC=fitAEff->GetParError(1);

    std::cout << "up" << std::endl;
    
    /*} else {
    
    *N=-10;
    *sN=0;
    *C=-10;
    *sC=0;
    
    }*/
  
  
  delete CanvasTimeResVsAEff;
}



TGraphErrors* GetTimeResVsE(TH2D* Histo, std::string Ch,std::string DirData){
  
  TH1D* proj;
  TF1* fitProj = new TF1("fitProj","gaus");
  TCanvas* canv = new TCanvas("canv","canv",700,700);
  
  Double_t ResPixel=91.9,SResPixel=0.1;
  Double_t Res,SRes;
  

  std::vector<Double_t> X,Y,sY;


  Double_t Value[2], Weight[2];

  for(int i=0; i<Histo->GetXaxis()->GetNbins()-2; i+=2){ 
    
    proj = Histo->ProjectionY(Form("Histo%i_%i",i,(int)Histo->GetXaxis()->GetBinCenter(i+1)),i,i+1); 
    proj->Draw();
    proj->Fit(fitProj);

    fitProj->SetRange(fitProj->GetParameter(1)-2*fitProj->GetParameter(2), fitProj->GetParameter(1)+2*fitProj->GetParameter(2));
    
    proj->Fit(fitProj,"R");

    Res=fitProj->GetParameter(2);
    SRes=fitProj->GetParError(2);

    Value[0]=Histo->GetXaxis()->GetBinCenter(i);
    Value[1]=Histo->GetXaxis()->GetBinCenter(i+1);
    
    Weight[0]=Histo->GetBinContent(i);
    Weight[1]=Histo->GetBinContent(i+1);

    if(proj->GetEntries()>1000 && fitProj->GetParError(2)<20){

      X.push_back( TMath::Mean(2,Value,Weight) );
      Y.push_back( sqrt(Res*Res - ResPixel*ResPixel) );
      sY.push_back(  sqrt( Res*Res/(Res*Res+ResPixel*ResPixel)*SRes*SRes + ResPixel*ResPixel/(Res*Res+ResPixel*ResPixel)*SResPixel*SResPixel ) ); 

      canv->SaveAs((DirData+"/Plot/TR/ProjTResE"+Ch+"/TimeResVsAEff"+std::to_string(i)+".png").c_str());  
    }//chiudo if
    
 
  }//chiudo for
  
  TGraphErrors* graph = new TGraphErrors(X.size(),&X[0],&Y[0],0,&sY[0]);
  graph->SetTitle(Form("TimeResVsECh%s",Ch.c_str()));
  graph->SetName(graph->GetTitle());
  graph->GetXaxis()->SetTitle("ADC [D.U]");
  graph->GetYaxis()->SetTitle("TimeRes [ps]");
  
  return graph;
  
}
