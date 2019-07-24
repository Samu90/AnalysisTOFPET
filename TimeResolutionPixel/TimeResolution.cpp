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
void GetPedestal(TTree* tree,Double_t* pedch1,Double_t* pedch2, Double_t* RMSPed1, Double_t* RMSPed2){
  
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
  *RMSPed1 = histoCh1->GetRMS();
  
  *pedch2 = histoCh2->GetMean();
  *RMSPed2 = histoCh2->GetRMS();  

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
  Float_t NS=3;

  tree->SetBranchAddress("time",time);
  tree->SetBranchAddress("energy",energy);
  
  std::cout << "_____CutInterval59____  " << mean59-NS*RMS59 << " " << mean59+NS*RMS59 << "  _____CutInterval315____  " << mean315-NS*RMS315 << " " << mean315+NS*RMS315 << std::endl;

  for(int i=0; i<tree->GetEntries();i++){
    tree->GetEntry(i);
    
    if((energy[0]-PedCh59)>(mean59-NS*RMS59) && (energy[0]-PedCh59)<(mean59+NS*RMS59) && (energy[1]-PedCh315)>(mean315-NS*RMS315) && (energy[1]-PedCh315)<(mean315+NS*RMS315)){ 
      Tdiff->Fill(time[0]-time[1]);
    }//chiudo if
    
  }//chiudo for
  TF1* preFit = new TF1("TempFit","gaus", (Tdiff->GetMean())-2*(Tdiff->GetRMS()), (Tdiff->GetMean())+2*(Tdiff->GetRMS()));

  preFit->SetParameter(0,Tdiff->GetMaximum());
  preFit->SetParameter(1,Tdiff->GetBinCenter(Tdiff->GetMaximumBin()));
  
  Tdiff->Fit(preFit,"RQ");

  FitTdiff->SetParameter(0,preFit->GetParameter(0));
  FitTdiff->SetParameter(1,preFit->GetParameter(1));
  FitTdiff->SetParameter(2,preFit->GetParameter(2));

  FitTdiff->SetRange((Tdiff->GetMean())-2*(Tdiff->GetRMS()), (Tdiff->GetMean())+2*(Tdiff->GetRMS()));

  Tdiff->Fit(FitTdiff,"RQ");
  
}
///////////////////////////////////////////////////////////////////////////////////////

void AmpCorrection(TTree*, TH2D*, TH2D*, TF1*, TF1*, Double_t, Double_t,TF1*,TF1*);
void GetTDiffVsAEff(TTree* , TH2D* , TF1*, TF1*, Double_t, Double_t , TF1*, TF1*, Double_t, Double_t);
TGraphErrors* GetTimeResVsE(TH2D*, Int_t, std::string,TF1*  );

//////////////////////////////////////////////////////////////////////////////////////

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
  gSystem->Exec(("mkdir "+DirData+"/Plot/TimeRes/ProjAEff").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TimeRes/ProjTResE59").c_str());
  gSystem->Exec(("mkdir "+DirData+"/Plot/TimeRes/ProjTResE315").c_str());


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

  Double_t RMSPedestal[NFilePhys][2];
  Double_t RMSPedestalCh1[2];
  Double_t RMSPedestalCh2[2];

  int k=0;
  ////////////////////////////////////////////////////CALCOLO PIEDISTALLI
  for(int i=0;i < (int)FileListPedestal.size()-1;i+=2){
    
    std::cout<< "open file:  " << (DirData+"/"+FileListPedestal.at(i)).c_str() << std::endl;
    TFile* f0= TFile::Open((DirData+"/"+FileListPedestal.at(i)).c_str());
    TFile* f1= TFile::Open((DirData+"/"+FileListPedestal.at(i)).c_str());

    TTree* tree0 = (TTree*)f0->Get("data"); //Before
    TTree* tree1 = (TTree*)f1->Get("data"); //After
    
    GetPedestal(tree0,&PedestalCh1[0],&PedestalCh2[0],&RMSPedestalCh1[0],&RMSPedestalCh2[0]);
    GetPedestal(tree1,&PedestalCh1[1],&PedestalCh2[1],&RMSPedestalCh1[1],&RMSPedestalCh2[1]);
    
    Pedestal[k][0]=(PedestalCh1[0]+PedestalCh1[1])/2;
    Pedestal[k][1]=(PedestalCh2[0]+PedestalCh2[1])/2;

    RMSPedestal[k][0] = (RMSPedestalCh1[0]+RMSPedestalCh1[1])/2;
    RMSPedestal[k][1] = (RMSPedestalCh2[0]+RMSPedestalCh2[1])/2;
    
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

  std::string FileName =DirData;
  FileName.erase(0,3);

  TFile* f = new TFile(("../RootFileGraphPixel/"+FileName+"TRes.root").c_str(),"RECREATE");
  
  TCanvas* TResVsTemp = new TCanvas("TResVsTemp","TResVsTemp",1600,800);
  TGraphErrors* GraphTResVsTemp = new TGraphErrors(NFilePhys,MeanTGlobal,TimeRes,SigmaTGlobal,SigmaTimeRes);
  GraphTResVsTemp->SetTitle("GraphTimeResolutionVsTemp");
  GraphTResVsTemp->SetName(GraphTResVsTemp->GetTitle());
  GraphTResVsTemp->GetXaxis()->SetTitle("TempMeanBox [°C]");
  GraphTResVsTemp->GetYaxis()->SetTitle("TimeResolution [ps]");

  GraphTResVsTemp->Draw("AP");
  
  f->cd();
  GraphTResVsTemp->Write();

  TResVsTemp->SaveAs((DirData+"/Plot/TimeRes/TimeResVsTemp.png").c_str());

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas* CanvasTotalTimeRes= new TCanvas("CanvasTotalTimeRes","CanvasTotalTimeRes",1600,800);
  TH1D* TimeResTot = new TH1D("TimeResTot","TimeResTot",200,-2000,2000);
  TF1* FitTimeResTot = new TF1("FitTimeResTot","gaus");

  for(int i=0; i< NFilePhys; i++){
    TimeResTot->Add(HistoTdiff[i]);
  }
  
  FitTimeResTot->SetParameter(0,TimeResTot->GetMaximum());
  FitTimeResTot->SetParameter(1,TimeResTot->GetBinCenter(TimeResTot->GetMaximumBin()));
  TimeResTot->Draw();    
  FitTimeResTot->SetRange( (TimeResTot->GetMean())-2*(TimeResTot->GetRMS()) , (TimeResTot->GetMean())+2*(TimeResTot->GetRMS()) );
  TimeResTot->Fit(FitTimeResTot,"R");
  
  CanvasTotalTimeRes->SaveAs((DirData+"/Plot/TimeRes/TimeResTotal.png").c_str());

  //void AmpCorrection(TTree* tree, TH2D* Histo59, TH2D* Histo315, TF1* Spectrum59, TF1* Spectrum315, Double_t PedCh59, Double_t PedCh315)

  TH2D* TDiffVsE59 = new TH2D("TDiffVsE59","TDiffVsE59",100,0,100,200,-2000,2000);
  TH2D* TDiffVsE315 = new TH2D("TDiffVsE315","TDiffVsE315",100,0,100,200,-2000,2000);
  
  TDiffVsE59->GetXaxis()->SetTitle("ADC [D.U]");
  TDiffVsE59->GetXaxis()->SetTitle("Tdiff[ps]");
  
  TDiffVsE315->GetXaxis()->SetTitle("ADC [D.U]");
  TDiffVsE315->GetXaxis()->SetTitle("Tdiff[ps]");
  
  for(int i=0;i < NFilePhys;i++){
    
    std::cout << "i : " << i <<" File: " << DirData+"/"+FileListPhysicsCoinc.at(i) << std::endl;
    f0= TFile::Open((DirData+"/"+FileListPhysicsCoinc.at(i)).c_str());
    
    if(f0) std::cout << "File Correctly Opend" << std::endl;
    
    tree0 = (TTree*)f0->Get("data");
    
    AmpCorrection(tree0,TDiffVsE59,TDiffVsE315,FitSpectrum[i][0],FitSpectrum[i][1],Pedestal[i][0],Pedestal[i][1],0,0);
   
  }

  TF1* FitTdiffVsE59 = new TF1("FitTdiffVsE59","pol3");
  TF1* FitTdiffVsE315 = new TF1("FitTdiffVsE315","pol3");

  TCanvas* CanvasTDiffVsE = new TCanvas("CanvasTDiffVsE","CanvasTDiffVsE",1400,700);
  CanvasTDiffVsE->Divide(2,1);
  
  CanvasTDiffVsE->cd(1);
  TDiffVsE59->Draw("COLZ");
  TDiffVsE59->Fit(FitTdiffVsE59);

  CanvasTDiffVsE->cd(2);
  TDiffVsE315->Draw("COLZ");
  TDiffVsE315->Fit(FitTdiffVsE315);  

  CanvasTDiffVsE->SaveAs((DirData+"/Plot/TimeRes/TdiffVsE.png").c_str());
  
  TH2D* TDiffCorrVsE59 = new TH2D("TDiffCorrVsE59","TDiffCorrVsE59",100,0,100,200,-2000,2000);
  TH2D* TDiffCorrVsE315 = new TH2D("TDiffCorrVsE315","TDiffCorrVsE315",100,0,100,200,-2000,2000);

  TDiffCorrVsE59->GetXaxis()->SetTitle("ADC [D.U]");
  TDiffCorrVsE59->GetXaxis()->SetTitle("Tdiff[ps]");
  
  TDiffCorrVsE315->GetXaxis()->SetTitle("ADC [D.U]");
  TDiffCorrVsE315->GetXaxis()->SetTitle("Tdiff[ps]");


  for(int i=0;i < NFilePhys;i++){
    
    std::cout << "i : " << i <<" File: " << DirData+"/"+FileListPhysicsCoinc.at(i) << std::endl;
    f0= TFile::Open((DirData+"/"+FileListPhysicsCoinc.at(i)).c_str());
    
    if(f0) std::cout << "File Correctly Opend" << std::endl;
    
    tree0 = (TTree*)f0->Get("data");
    
    AmpCorrection(tree0,TDiffCorrVsE59,TDiffCorrVsE315,FitSpectrum[i][0],FitSpectrum[i][1],Pedestal[i][0],Pedestal[i][1],FitTdiffVsE59,FitTdiffVsE315);
    
  }

  
  
  TF1* FitTdiffCorrVsE59 = new TF1("FitTdiffCorrVsE59","pol3");
  TF1* FitTdiffCorrVsE315 = new TF1("FitTdiffCorrVsE315","pol3");
  
  TCanvas* CanvasTdiffCorrVsE = new TCanvas("CanvasTdiffCorrVsE","CanvasTdiffCorrVsE",1400,700);
  CanvasTdiffCorrVsE->Divide(2,1);
  
  CanvasTdiffCorrVsE->cd(1);
  TDiffCorrVsE59->Draw("COLZ");
  TDiffCorrVsE59->Fit(FitTdiffCorrVsE59);

  CanvasTdiffCorrVsE->cd(2);
  TDiffCorrVsE315->Draw("COLZ");
  TDiffCorrVsE315->Fit(FitTdiffCorrVsE315);
  
  CanvasTdiffCorrVsE->SaveAs((DirData+"/Plot/TimeRes/TdiffVsECorr.png").c_str());

  std::cout <<FitTdiffCorrVsE59->GetParameter(0) << " "<<FitTdiffCorrVsE59->GetParameter(1) << " "<<FitTdiffCorrVsE59->GetParameter(2) << " "<<FitTdiffCorrVsE59->GetParameter(3) << std::endl;
  std::cout <<FitTdiffCorrVsE315->GetParameter(0) << " "<<FitTdiffCorrVsE315->GetParameter(1) << " "<<FitTdiffCorrVsE315->GetParameter(2) << " "<<FitTdiffCorrVsE315->GetParameter(3) << std::endl;

  TH2D* TDiffVsAEff = new TH2D("TDiffVsAEff","TDiffVsAEff",100,0,100,200,-2000,2000);
  TH2D* TDiffCorrVsAEff = new TH2D("TDiffCorrVsAEff","TDiffCorrVsAEff",100,0,100,200,-2000,2000);

  TDiffVsAEff->GetXaxis()->SetTitle("AEff/#sigma_n");
  TDiffVsAEff->GetYaxis()->SetTitle("Tdiff[ps]");

  TDiffCorrVsAEff->GetXaxis()->SetTitle("AEff/#sigma_n");
  TDiffCorrVsAEff->GetYaxis()->SetTitle("Tdiff[ps]");
  

  for(int i=0;i < NFilePhys;i++){
    
    std::cout << "i : " << i <<" File: " << DirData+"/"+FileListPhysicsCoinc.at(i) << std::endl;
    f0= TFile::Open((DirData+"/"+FileListPhysicsCoinc.at(i)).c_str());
    
    if(f0) std::cout << "File Correctly Opend" << std::endl;
    
    tree0 = (TTree*)f0->Get("data");
       
    GetTDiffVsAEff(tree0,TDiffVsAEff,FitSpectrum[i][0],FitSpectrum[i][1],Pedestal[i][0],Pedestal[i][1],0,0,RMSPedestal[i][0],RMSPedestal[i][1]); //uncorrecter
    GetTDiffVsAEff(tree0,TDiffCorrVsAEff,FitSpectrum[i][0],FitSpectrum[i][1],Pedestal[i][0],Pedestal[i][1],FitTdiffVsE59,FitTdiffVsE315,RMSPedestal[i][0],RMSPedestal[i][1]); //corrected

  }
  
    

  TCanvas* CanvasTDiffVsAEff = new TCanvas("CanvasTDiffVsAEff","CanvasTDiffVsAEff",1400,700);
  CanvasTDiffVsAEff->Divide(2,1);
  
  CanvasTDiffVsAEff->cd(1);
  TDiffVsAEff->GetXaxis()->SetRangeUser(0,40);
  TDiffVsAEff->Draw("COLZ");
  
  CanvasTDiffVsAEff->cd(2);
  TDiffCorrVsAEff->GetXaxis()->SetRangeUser(0,40); 
  TDiffCorrVsAEff->Draw("COLZ");

  CanvasTDiffVsAEff->SaveAs((DirData+"/Plot/TimeRes/TdiffVsAEff.png").c_str()); 


  
  TH1D* proj;
  TF1* FitProj = new TF1("FitProj","gaus");
  TCanvas* CanvasPlot = new TCanvas("CanvasPlot","CanvasPlot",700,700);
  
  std::vector<Double_t> X,Y,sY;
  
  Double_t ResPixel,SResPixel;
  Double_t Res,SRes;

  ResPixel= FitTimeResTot->GetParameter(2)/sqrt(2);
  SResPixel= FitTimeResTot->GetParError(2)/sqrt(2);

  for(int i=0; i<(TDiffCorrVsAEff->GetXaxis()->GetLast()-1); i+=2){
    
    proj=TDiffCorrVsAEff->ProjectionY(Form("Proj%i",i),i,i+1);
    
    CanvasPlot->cd(); 
    proj->Draw();
    proj->Fit(FitProj);
    
    FitProj->SetRange( FitProj->GetParameter(1)-1.7*FitProj->GetParameter(2), FitProj->GetParameter(1)+1.7*FitProj->GetParameter(2) );

    proj->Fit(FitProj,"R");
    
    Res = FitProj->GetParameter(2); 
    SRes = FitProj->GetParError(2);

    if(FitProj->GetParError(2)<50 && proj->GetEntries()>100){
      

      X.push_back( TDiffCorrVsAEff->GetXaxis()->GetBinCenter(i+1) );
      Y.push_back( sqrt(Res*Res-ResPixel*ResPixel) );
      sY.push_back( sqrt( Res*Res/(Res*Res+ResPixel*ResPixel)*SRes*SRes + ResPixel*ResPixel/(Res*Res+ResPixel*ResPixel)*SResPixel*SResPixel ) );
                            
      CanvasPlot->SaveAs((DirData+"/Plot/TimeRes/ProjAEff/Proj"+std::to_string(i)+".png").c_str());
    }
    
  }
  
  
  
  TGraphErrors* TResVsAEff = new TGraphErrors(X.size(),&X[0],&Y[0],0,&sY[0]);

  TF1* fitTresVsAEff = new TF1("fitTresVsAEff","sqrt( [0]*[0]/(x*x)+2*[1]*[1] )",X.at(0)+1, X.at(X.size()-1)-1);
  // TF1* fitTresVsAEff = new TF1("fitTresVsAEff","sqrt( [0]*[0]/(x)+2*[1]*[1] )",X.at(0)+1, X.at(X.size()-1)-1);
  
  TCanvas* CanvasTResVsAEff = new TCanvas("CanvasTResVsAEff","CanvasTResVsAEff",700,700);
  TResVsAEff->SetMarkerStyle(8);
  TResVsAEff->SetTitle("TResVsAEff");
  TResVsAEff->SetName(TResVsAEff->GetTitle());
  TResVsAEff->Draw("AP");
  f->cd();
  TResVsAEff->Write();
  
  fitTresVsAEff->SetParameter(0,1800);
  fitTresVsAEff->SetParameter(1,100);
  
  TResVsAEff->Fit(fitTresVsAEff,"RWE");

  CanvasTResVsAEff->SaveAs((DirData+"/Plot/TimeRes/TimeResVsAEff.png").c_str()); 
  
  TGraphErrors* TimeResVsEnergyCh59;
  TGraphErrors* TimeResVsEnergyCh315;

  TimeResVsEnergyCh59 = GetTimeResVsE(TDiffCorrVsE59, 59, DirData, FitTimeResTot);
  TimeResVsEnergyCh315 = GetTimeResVsE(TDiffCorrVsE315, 315, DirData, FitTimeResTot);

  TF1* fitTResVsECh59 = new TF1("fitTResVsECh59","sqrt( [0]*[0]/(x*x) + [1]*[1]/(x) + [2]*[2] )");
  TF1* fitTResVsECh315 = new TF1("fitTResVsECh315","sqrt( [0]*[0]/(x*x) + [1]*[1]/(x) + [2]*[2] )");

  TCanvas* CanvTResVsE = new TCanvas("CanvTResVsE","CanvTResVsE",1400,700);
  CanvTResVsE->Divide(2,1);

  CanvTResVsE->cd(1);
  TimeResVsEnergyCh59->Draw("AP");
  fitTResVsECh59->SetParameter(0,1000);
  fitTResVsECh59->SetParameter(1,40);
  TimeResVsEnergyCh59->Fit(fitTResVsECh59);

  CanvTResVsE->cd(2);
  TimeResVsEnergyCh315->Draw("AP");
  fitTResVsECh315->SetParameter(0,1000);
  fitTResVsECh315->SetParameter(1,40);
  TimeResVsEnergyCh315->Fit(fitTResVsECh315);
  
  
  CanvTResVsE->SaveAs((DirData+"/Plot/TimeRes/TimeResVsE.png").c_str());
  
  
  f->cd();
  TimeResTot->Write();
  FitTimeResTot->Write();
  TDiffVsE59->Write();
  TDiffVsE315->Write();
  TDiffCorrVsE59->Write();
  TDiffCorrVsE315->Write();
  TDiffVsAEff->Write();
  TDiffCorrVsAEff->Write();
  TimeResVsEnergyCh59->Write();
  TimeResVsEnergyCh315->Write();

  f->Save();
  f->Close();

  std::cout << "comand->  " <<("cp ../RootFileGraphPixel/"+FileName+"TRes.root "+DirData+"/Plot/.").c_str() << std::endl; 
  gSystem->Exec(("cp ../RootFileGraphPixel/"+FileName+"TRes.root "+DirData+"/Plot/.").c_str());

  return 0;














}



void AmpCorrection(TTree* tree, TH2D* Histo59, TH2D* Histo315, TF1* Spectrum59, TF1* Spectrum315, Double_t PedCh59, Double_t PedCh315,TF1* corr59,TF1* corr315){
  
  Double_t energy[2];
  Double_t time[2];
  
  tree->SetBranchAddress("time",time); 
  tree->SetBranchAddress("energy",energy);
  
  for(int i=0; i< tree->GetEntries();i++){
    tree->GetEntry(i);
    
    if(energy[0]-PedCh59 > 8 && energy[1]-PedCh315 > (Spectrum315->GetParameter(1)-3*Spectrum315->GetParameter(2)) 
       && energy[1]-PedCh315 < (Spectrum315->GetParameter(1)+1*Spectrum315->GetParameter(2)) 
       && energy[0]-PedCh59 < (Spectrum59->GetParameter(1)+1*Spectrum59->GetParameter(2))){

      if(corr59) Histo59->Fill(energy[0]-PedCh59, time[0]-time[1] - corr59->Eval(energy[0]-PedCh59) );
      else Histo59->Fill(energy[0]-PedCh59, time[0]-time[1]);
      }//chiudo if

    if(energy[1]-PedCh315 > 8 && energy[0]-PedCh59 > (Spectrum59->GetParameter(1)-3*Spectrum59->GetParameter(2)) 
       && energy[0]-PedCh59 < (Spectrum59->GetParameter(1)+1*Spectrum59->GetParameter(2)) 
       && energy[1]-PedCh315 < (Spectrum315->GetParameter(1)+1*Spectrum315->GetParameter(2))){
      
      if(corr315) Histo315->Fill(energy[1]-PedCh315, time[0]-time[1] - corr315->Eval(energy[1]-PedCh315) );
      else Histo315->Fill(energy[1]-PedCh315, time[0]-time[1] );
    }//chiudo if
    
    
  }//chiudo for
  
  
}



void GetTDiffVsAEff(TTree* tree, TH2D* Histo, TF1* Spectrum59, TF1* Spectrum315, Double_t PedCh59, Double_t PedCh315, TF1* corr59, TF1* corr315, Double_t RMSPedCh59, Double_t RMSPedCh315){
  
  Double_t energy[2];
  Double_t time[2];
  
  Double_t A1,A2;
  Double_t AmpEff;
  
  tree->SetBranchAddress("time",time); 
  tree->SetBranchAddress("energy",energy);
  
  
  
  for(int i=0; i< tree->GetEntries();i++){
    tree->GetEntry(i);
    
    A1=energy[0]-PedCh59;
    A2=energy[1]-PedCh315;

    if(A1 > 8 && A1 < ( Spectrum59->GetParameter(1)+2.5*Spectrum59->GetParameter(3) ) && A2 > 8 && A2 < ( Spectrum315->GetParameter(1)+2.5*Spectrum315->GetParameter(3) )){
     
      AmpEff=A1*A2/(sqrt(A1*A1+A2*A2))/((RMSPedCh59+RMSPedCh315)/2);
      //    AmpEff=A1*A2/(A1+A2);
      
      if(corr59 && corr315) {
	
	Histo->Fill(AmpEff, ( time[0]-corr59->Eval(A1) ) - ( time[1]+corr315->Eval(A2) ) );
	
      }
      else Histo->Fill(AmpEff, time[0]-time[1]);
      
    }//chiudo if
    
  }//chiudo for
  
}


TGraphErrors* GetTimeResVsE(TH2D* Histo, Int_t Ch,std::string DirData,TF1* FitTResTot){
  
  TH1D* proj;
  TF1* fitProj = new TF1("fitProj","gaus");
  TCanvas* canv = new TCanvas("canv","canv",700,700);
  
  std::vector<Double_t> X,Y,sY;

  Double_t ResPixel,SResPixel;
  Double_t Res,SRes;
  
  ResPixel= FitTResTot->GetParameter(2)/sqrt(2);
  SResPixel= FitTResTot->GetParError(2)/sqrt(2);

  Int_t step=2;

  for(int i=0; i<Histo->GetXaxis()->GetNbins()-step; i+=step){ 
    
    proj = Histo->ProjectionY(Form("Histo%i_%i",i,(int)Histo->GetXaxis()->GetBinCenter(i+1)),i,i+1); 
    proj->Draw();
    proj->Fit(fitProj);

    fitProj->SetRange(fitProj->GetParameter(1)-1.7*fitProj->GetParameter(2), fitProj->GetParameter(1)+1.7*fitProj->GetParameter(2));
    
    proj->Fit(fitProj,"R");

    Res = fitProj->GetParameter(2); 
    SRes = fitProj->GetParError(2);

    Double_t Value[step],Weight[step];

    if(proj->GetEntries()>1000 && fitProj->GetParError(2)<20){
      Value[0]=Histo->GetXaxis()->GetBinCenter(i);
      Value[1]=Histo->GetXaxis()->GetBinCenter(i+1);
      Weight[0]=Histo->GetBinContent(i);
      Weight[1]=Histo->GetBinContent(i+1);

      X.push_back( TMath::Mean(step,Value,Weight) );
      Y.push_back( sqrt(Res*Res-ResPixel*ResPixel) );
      sY.push_back( sqrt( Res*Res/(Res*Res+ResPixel*ResPixel)*SRes*SRes + ResPixel*ResPixel/(Res*Res+ResPixel*ResPixel)*SResPixel*SResPixel ) );
      

      canv->SaveAs((DirData+"/Plot/TimeRes/ProjTResE"+std::to_string(Ch)+"/TimeRes"+std::to_string(i)+".png").c_str());  
    }//chiudo if
    
 
  }//chiudo for
  
  TGraphErrors* graph = new TGraphErrors(X.size(),&X[0],&Y[0],0,&sY[0]);
  graph->SetTitle(Form("TimeResVsECh%i",Ch));
  graph->SetName(graph->GetTitle());
  graph->GetXaxis()->SetTitle("ADC [D.U]");
  graph->GetYaxis()->SetTitle("TimeRes [ps]");
  
  return graph;
  
}
