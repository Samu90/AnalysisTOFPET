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
void GetTdiffCorr(TTree*, TH1D** ,TH2D**,  TF1** ,TF1**, Double_t* , std::vector<int> , std::vector<std::string> );

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
  Double_t RMSPedestal[NFilePhys][(int)Channels.size()];
  Double_t PedestalBefore[(int)Channels.size()];
  Double_t PedestalAfter[(int)Channels.size()];
  Double_t RMSPedestalBefore[(int)Channels.size()];
  Double_t RMSPedestalAfter[(int)Channels.size()];
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  int k=0;

  for(int i=0;i < (int)FileListPedestal.size()-1;i+=2){
    
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

  tdiffVsE[0] = new TH2D("tdiffVsE_time1","tdiffVsE_time1",24,24,48,200,-2000,2000);
  tdiffVsE[1] = new TH2D("tdiffVsE_time2","tdiffVsE_time2",24,24,48,200,-2000,2000);

  TH2D* tdiffVsT[2];

  tdiffVsT[0] = new TH2D("tdiffVsT_time1","tdiffVsT_time1",40,21,25,200,-2000,2000);
  tdiffVsT[1] = new TH2D("tdiffVsT_time2","tdiffVsT_time2",40,21,25,200,-2000,2000);
  
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
      
      Double_t min=FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(1)-3*FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(2);
      Double_t max=FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(1)+3*FitSpectrum[i][mapCh[Channels[k]]]->GetParameter(2);
      
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
  
  TF1* fitTdiffVsE[2];
  fitTdiffVsE[0] = new TF1("fitTdiffVsE0","pol1");
  fitTdiffVsE[1] = new TF1("fitTdiffVsE1","pol1");

  TCanvas* CanvasTdiffVsE = new TCanvas("CanvasTdiffVsE","CanvasTdiffVsE",1500,700);
  CanvasTdiffVsE->Divide((int)Channels.size()-1,1);
  CanvasTdiffVsE->cd(1);
  
  tdiffVsE[0]->GetXaxis()->SetTitle("Energy [D.U.]");
  tdiffVsE[0]->GetYaxis()->SetTitle("time1-timeRef [ps]");
  tdiffVsE[0]->Draw("COLZ");
  tdiffVsE[0]->Fit(fitTdiffVsE[0]);

  CanvasTdiffVsE->cd(2);
  tdiffVsE[1]->GetXaxis()->SetTitle("Energy [D.U.]");
  tdiffVsE[1]->GetYaxis()->SetTitle("time2-timeRef [ps]");
  tdiffVsE[1]->Draw("COLZ");
  tdiffVsE[1]->Fit(fitTdiffVsE[1]);

  CanvasTdiffVsE->SaveAs((DirData+"/Plot/TR/TdiffVsE.png").c_str());
  
  
  TCanvas* CanvasTdiffVsT = new TCanvas("CanvasTdiffVsT","CanvasTdiffVsT",1500,700);
  CanvasTdiffVsT->Divide((int)Channels.size()-1,1);
  
  CanvasTdiffVsT->cd(1);
  tdiffVsT[0]->GetXaxis()->SetTitle("Energy [D.U.]");
  tdiffVsT[0]->GetYaxis()->SetTitle("time1-timeRef [ps]");
  tdiffVsT[0]->Draw("COLZ");

  CanvasTdiffVsT->cd(2);
  tdiffVsT[1]->GetXaxis()->SetTitle("Energy [D.U.]");
  tdiffVsT[1]->GetYaxis()->SetTitle("time2-timeRef [ps]");
  tdiffVsT[1]->Draw("COLZ");

  CanvasTdiffVsT->SaveAs((DirData+"/Plot/TR/TdiffVsT.png").c_str());
  
  /////////////////CORREZIONE//////////////////////
  
  TH1D* tdiffCorr[3];
  
  tdiffCorr[0] = new TH1D("tref-tave_corr","tref-tave_corr",200,-2000,2000);
  tdiffCorr[1] = new TH1D("tref-t1_corr","tref-t1_corr",200,-2000,2000);
  tdiffCorr[2] = new TH1D("tref-t2_corr","tref-t2_corr",200,-2000,2000);

  TH2D* tdiffCorrVsE[2];

  tdiffCorrVsE[0] = new TH2D("tdiffCorrVsE_time1","tdiffCorrVsE_time1",24,24,48,200,-2000,2000);
  tdiffCorrVsE[1] = new TH2D("tdiffCorrVsE_time2","tdiffCorrVsE_time2",24,24,48,200,-2000,2000);
  
  /*for(int i=0;i<NFilePhys;i++){
    
    f0= TFile::Open((DirData+"/"+FileListPhysics.at(i)).c_str());
    tree0 = (TTree*)f0->Get("data"); 

    GetTdiffCorr(tree0,tdiffCorr,tdiffCorrVsE,FitSpectrum[i],fitTdiffVsE,Pedestal[i],Channels,ChannelType);
    
    }*/

  

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
  fitTdiffCorrVsE[0] = new TF1("fitTdiffCorrVsE0","pol1");
  fitTdiffCorrVsE[1] = new TF1("fitTdiffCorrVsE1","pol1");

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
  
  TCanvas* CanvasAEff = new TCanvas("CanvasAEff","CanvasAEff",800,800);
  tdiffVsAmpEff->GetXaxis()->SetTitle("A_Eff/#sigma_{n} [D.U]");
  tdiffVsAmpEff->GetYaxis()->SetTitle("T_diff [ps]");
  tdiffVsAmpEff->Draw("COLZ");

  CanvasAEff->SaveAs((DirData+"/Plot/TR/TdiffVsAEff.png").c_str());  
  
  
  Int_t nbins = tdiffVsAmpEff->GetXaxis()->GetNbins();
  
  TH1D* Projection;
  TF1* FitProjection;
  TCanvas* CanvasProjection;

  std::vector<Double_t> AEffValue;
  std::vector<Double_t> TimeResValue;
  std::vector<Double_t> ErrTimeResValue;

  Int_t fitStatus;

  for(int i=0; i<nbins-1; i+=2){

    //ProjectionAEff
    CanvasProjection = new TCanvas("CanvasProjection","CanvasProjection",800,800);
    
    Projection = tdiffVsAmpEff->ProjectionY(Form("Projection%i_Amp%lf",i,tdiffVsAmpEff->GetXaxis()->GetBinCenter(i)),i,i);
    Projection->GetXaxis()->SetTitle("t_diff [ps]");
    Projection->GetYaxis()->SetTitle("Counts");
    Projection->GetXaxis()->SetRangeUser(-3500,3000);  
    Projection->SetTitle(Projection->GetName());
    Projection->Draw();
    
    FitProjection = new TF1(Form("fitProjection%i",i) , "gaus" , Projection->GetBinCenter(Projection->GetMaximumBin())-1.1*Projection->GetRMS(),Projection->GetBinCenter(Projection->GetMaximumBin())+1.1*Projection->GetRMS());
    
    std::cout << "Range_______:    " << Projection->GetBinCenter(Projection->GetMaximumBin())-1.1*Projection->GetRMS() << "   " << Projection->GetBinCenter(Projection->GetMaximumBin())+1.1*Projection->GetRMS() << std::endl;
          
    FitProjection->SetParameter(0 , Projection->GetMaximum());
    FitProjection->SetParameter(1 , Projection->GetMaximumBin());
    
    fitStatus = Projection->Fit(FitProjection,"RQ");
    
    CanvasProjection->SaveAs((DirData+"/Plot/TR/ProjectionAEff/Projection"+std::to_string(i)+".png").c_str());
    
    if(fitStatus==0 && tdiffVsAmpEff->GetXaxis()->GetBinCenter(i)>0){
      AEffValue.push_back(tdiffVsAmpEff->GetXaxis()->GetBinCenter(i));
      TimeResValue.push_back(sqrt(FitProjection->GetParameter(2)*FitProjection->GetParameter(2)-93*93));
      ErrTimeResValue.push_back(FitProjection->GetParError(2));
    }
    
    
    delete Projection;
    delete CanvasProjection;

    std::cout << tdiffVsAmpEff->GetXaxis()->GetBinCenter(i) << " " << FitProjection->GetParameter(2) << " " << FitProjection->GetParError(2) << std::endl;
  }
  
  TCanvas* CanvasTimeResVsAEff = new TCanvas("CanvasTimeResVsAEff","CanvasTimeResVsAEff",800,800);
  
  
  TGraphErrors* TimeResVsAEff = new TGraphErrors((int)AEffValue.size(),&AEffValue[0],&TimeResValue[0],0,&ErrTimeResValue[0]);
  
  TF1* fitAEff = new TF1("fitAEff", "[0]/x + [1]",0,23);
  
  fitAEff->SetParameter(0,1000);
  fitAEff->SetParameter(1,40);

  TimeResVsAEff->SetTitle("TimeResVsAEff");
  TimeResVsAEff->SetName(TimeResVsAEff->GetTitle());
  TimeResVsAEff->GetXaxis()->SetTitle("A_Eff/#sigma_{n} [D.U]");
  TimeResVsAEff->GetYaxis()->SetTitle("TimeRes [ps]");
  
  TimeResVsAEff->SetMarkerStyle(8);
  
  TimeResVsAEff->Draw("AP");
  TimeResVsAEff->Fit(fitAEff,"RW");
  
  CanvasTimeResVsAEff->SaveAs((DirData+"/Plot/TR/TimeResVsAEff.png").c_str()); 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
    
    ped[mapCh[NCH[i]]]= histo[mapCh[NCH[i]]]->GetMean();
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

      tdiffVsE[0]->Fill(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]],time[mapCh[NCH[mapTime["time1"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]]);
      tdiffVsE[1]->Fill(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]],time[mapCh[NCH[mapTime["time2"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]]);
      
      tdiffVsT[0]->Fill(temp1,time[mapCh[NCH[mapTime["time2"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]]);
      tdiffVsT[1]->Fill(temp1,time[mapCh[NCH[mapTime["time2"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]]);
    
    }//chiudo if grande
    

    if(energy[mapCh[NCH[mapTime["timeRef"]]]]-ped[mapCh[NCH[mapTime["timeRef"]]]] < fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(2) && 
       energy[mapCh[NCH[mapTime["timeRef"]]]]-ped[mapCh[NCH[mapTime["timeRef"]]]] > fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(1)-3*fitspectrum[mapCh[NCH[mapTime["timeRef"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]] < fitspectrum[mapCh[NCH[mapTime["time1"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["time1"]]]]->GetParameter(2) &&
       energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]] < fitspectrum[mapCh[NCH[mapTime["time2"]]]]->GetParameter(1)+3*fitspectrum[mapCh[NCH[mapTime["time2"]]]]->GetParameter(2) &&
       chID[mapCh[NCH[mapTime["timeRef"]]]]!=-9 && chID[mapCh[NCH[mapTime["time1"]]]]!=-9 && chID[mapCh[NCH[mapTime["time2"]]]]!=-9){
      
      A1=energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]]+6; //Più 6 solo per far tornare il problema con l'elettroica e le eergie effettive nnegative//
      A2=energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]];

      AmpEff=A1*A2/(sqrt(A1*A1+A2*A2));

      tdiffVsAmpEff->Fill(AmpEff/MeanRMSPed,time[mapCh[NCH[mapTime["timeRef"]]]]-(time[mapCh[NCH[mapTime["time1"]]]]+time[mapCh[NCH[mapTime["time2"]]]])/2);
      
    }//chiudo altro if


  }


  
  std::cout << Form("time1 %i ",  NCH[mapTime["time1"]]) << Form("time2 %i ",     NCH[mapTime["time2"]]) << Form("timeRef %i ", NCH[mapTime["timeRef"]]) << std::endl;
  std::cout << Form("time1 %s ",  ChType[mapTime["time1"]].c_str()) << Form("time2 %s ",    ChType[mapTime["time2"]].c_str()) << Form("timeRef %s ", ChType[mapTime["timeRef"]].c_str()) << std::endl;
}




///////////////////////////////////////////////////////////////////////////////////////////////




void GetTdiffCorr(TTree* tree, TH1D** histo,TH2D** tdiffVsE, TF1** fitspectrum,TF1** fitcorr, Double_t* ped, std::vector<int> NCH, std::vector<std::string> ChType){
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
  
  Double_t NSigma=2.5;
  
  for(int i=0;i< tree->GetEntries();i++){
    tree->GetEntry(i);
    //Histo0 tref-tave, histo1 tref-t1, histo2 tref-t2
    
    if( energy[mapCh[NCH[0]]]-ped[mapCh[NCH[0]]]<fitspectrum[mapCh[NCH[0]]]->GetParameter(1)+NSigma*fitspectrum[mapCh[NCH[0]]]->GetParameter(2) && 
	energy[mapCh[NCH[0]]]-ped[mapCh[NCH[0]]]>fitspectrum[mapCh[NCH[0]]]->GetParameter(1)-NSigma*fitspectrum[mapCh[NCH[0]]]->GetParameter(2) && 
        energy[mapCh[NCH[1]]]-ped[mapCh[NCH[1]]]<fitspectrum[mapCh[NCH[1]]]->GetParameter(1)+NSigma*fitspectrum[mapCh[NCH[1]]]->GetParameter(2) && 
        energy[mapCh[NCH[1]]]-ped[mapCh[NCH[1]]]>fitspectrum[mapCh[NCH[1]]]->GetParameter(1)-NSigma*fitspectrum[mapCh[NCH[1]]]->GetParameter(2) && 
        energy[mapCh[NCH[2]]]-ped[mapCh[NCH[2]]]<fitspectrum[mapCh[NCH[2]]]->GetParameter(1)+NSigma*fitspectrum[mapCh[NCH[2]]]->GetParameter(2) && 
        energy[mapCh[NCH[2]]]-ped[mapCh[NCH[2]]]>fitspectrum[mapCh[NCH[2]]]->GetParameter(1)-NSigma*fitspectrum[mapCh[NCH[2]]]->GetParameter(2)){
      
      histo[0]->Fill(time[mapCh[NCH[mapTime["timeRef"]]]]-(time[mapCh[NCH[mapTime["time1"]]]]-fitcorr[0]->Eval(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]])+time[mapCh[NCH[mapTime["time2"]]]]-fitcorr[1]->Eval(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]]))/2);
      
      histo[1]->Fill(time[mapCh[NCH[mapTime["timeRef"]]]]-time[mapCh[NCH[mapTime["time1"]]]]-fitcorr[0]->Eval(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]]));
      
      histo[2]->Fill(time[mapCh[NCH[mapTime["timeRef"]]]]-time[mapCh[NCH[mapTime["time2"]]]]-fitcorr[1]->Eval(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]]));
      
      tdiffVsE[0]->Fill(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]],time[mapCh[NCH[mapTime["time1"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]]-fitcorr[0]->Eval(energy[mapCh[NCH[mapTime["time1"]]]]-ped[mapCh[NCH[mapTime["time1"]]]]));
      
      tdiffVsE[1]->Fill(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]],time[mapCh[NCH[mapTime["time2"]]]]-time[mapCh[NCH[mapTime["timeRef"]]]]-fitcorr[1]->Eval(energy[mapCh[NCH[mapTime["time2"]]]]-ped[mapCh[NCH[mapTime["time2"]]]]));
        
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

  TF1* spectrum = new TF1(Form("SpectrumFit_%s", Profile->GetName()),"[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) + [3] / (exp( (x*[4]-(2*[1]*[1]/([1]+2*[1])))) + 1)",10,48);
  
  spectrum->SetParameter(0,max);
  spectrum->SetParameter(1,peak1);
  spectrum->SetParameter(2,3);
  spectrum->SetParameter(3,700);
  spectrum->SetParameter(4,0.82);
  
  Profile->Fit(Form("SpectrumFit_%s", Profile->GetName()),"R0");

  return spectrum;
}

