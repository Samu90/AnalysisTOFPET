#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <utility> 

#include "TSystem.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TPaveStats.h"
#include "TVirtualFitter.h"
#include "TFormula.h"
#include "TList.h" 

struct dataPlot {
  Double_t A[2]={0,0};
  Double_t SA[2]={0,0};
  Double_t B[2]={0,0};
  Double_t SB[2]={0,0};
};

void ResidualPlot(TMultiGraph* ,TF1* , std::string );
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


TF1* CalibrationCurve(TGraphErrors* graph,std::string channel){
  
  TF1* fit = new TF1(("fitCalib"+channel).c_str()," [0] * (1-exp(-[1]*x) )",0,1300);

  fit->SetParameter(0,144);
  fit->SetParameter(1,4e-4);

  graph->Fit(fit,"QS");

  return fit;
}


void LYAtFixedTemp(Double_t TempEval, const std::pair<TF1*, Double_t>  FitTotCh59P1, const std::pair<TF1*, Double_t>  FitTotCh59P2, const std::pair<TF1*, Double_t> FitTotCh315P1,const  std::pair<TF1*, Double_t> FitTotCh315P2, dataPlot* data){

  
  Double_t YCh59[2], SYCh59[2];
  Double_t YCh315[2], SYCh315[2];
  Double_t XVector[2]={511,1275};
  
  Double_t sm59[2]={FitTotCh59P1.first->GetParError(0),FitTotCh59P2.first->GetParError(0)};
  Double_t sq59[2]={FitTotCh59P1.first->GetParError(1),FitTotCh59P2.first->GetParError(1)};
  
  Double_t sm315[2]={FitTotCh315P1.first->GetParError(0),FitTotCh315P2.first->GetParError(0)};
  Double_t sq315[2]={FitTotCh315P1.first->GetParError(1),FitTotCh315P2.first->GetParError(1)};

  
  YCh59[0]=FitTotCh59P1.first->Eval(TempEval);
  SYCh59[0]=sqrt(TempEval*TempEval*sm59[0]*sm59[0]+sq59[0]*sq59[0]+2*TempEval*FitTotCh59P1.second);
  
  YCh59[1]=FitTotCh59P2.first->Eval(TempEval);
  SYCh59[1]=sqrt(TempEval*TempEval*sm59[1]*sm59[1]+sq59[1]*sq59[1]+2*TempEval*FitTotCh59P2.second);
  
  YCh315[0]=FitTotCh315P1.first->Eval(TempEval);
  SYCh315[0]=sqrt(TempEval*TempEval*sm315[0]*sm315[0]+sq315[0]*sq315[0]+2*TempEval*FitTotCh315P1.second);
  
  YCh315[1]=FitTotCh315P2.first->Eval(TempEval);  
  SYCh315[1]=sqrt(TempEval*TempEval*sm315[1]*sm315[1]+sq315[1]*sq315[1]+2*TempEval*FitTotCh315P2.second);
  

  TGraphErrors* CalibrationCurve59 = new TGraphErrors(2,XVector,YCh59,0,SYCh59);
  TGraphErrors* CalibrationCurve315 = new TGraphErrors(2,XVector,YCh315,0,SYCh315);
  TF1* CalibFitCh59;
  TF1* CalibFitCh315;

  CalibrationCurve59->SetTitle(Form("CalibPlotCh59Temp%.1lf",TempEval));
  CalibrationCurve59->SetName(CalibrationCurve59->GetTitle());
  CalibrationCurve59->GetXaxis()->SetTitle("Energy [KeV]"); 
  CalibrationCurve59->GetYaxis()->SetTitle("Counts [D.U.]"); 

  CalibrationCurve315->SetTitle(Form("CalibPlotCh315Temp%.1lf",TempEval));
  CalibrationCurve315->SetName(CalibrationCurve315->GetTitle());
  CalibrationCurve315->GetXaxis()->SetTitle("Energy [KeV]"); 
  CalibrationCurve315->GetYaxis()->SetTitle("Counts [D.U.]"); 

  TCanvas* CanvasCalib = new TCanvas(Form("CanvasCalib%.0lf",TempEval),Form("CanvasCalib%.0lf",TempEval),1600,800);
  CanvasCalib->Divide(2,1);

  CanvasCalib->cd(1);
  CalibFitCh59 = CalibrationCurve(CalibrationCurve59,"Ch59");
  CalibrationCurve59->SetMarkerStyle(8);
  CalibrationCurve59->GetXaxis()->SetLimits(0,2000);
  CalibrationCurve59->GetYaxis()->SetLimits(0,100);
  CalibrationCurve59->GetYaxis()->SetRangeUser(0,100);
  CalibrationCurve59->Draw("AP");
  CalibFitCh59->Draw("SAME");

  CanvasCalib->cd(2);
  CalibFitCh315 = CalibrationCurve(CalibrationCurve315,"Ch315");
  CalibrationCurve315->SetMarkerStyle(8);
  CalibrationCurve315->GetXaxis()->SetLimits(0,2000);
  CalibrationCurve315->GetYaxis()->SetLimits(0,100);
  CalibrationCurve315->GetYaxis()->SetRangeUser(0,100);
  CalibrationCurve315->Draw("AP");
  CalibFitCh315->Draw("SAME");
  
  CanvasCalib->SaveAs(Form("Plot/CalibVsTempFit/CalibPlotT%.0lf.png",TempEval));

  data->A[0]=CalibFitCh59->GetParameter(0);
  data->SA[0]=CalibFitCh59->GetParError(0);
  data->A[1]=CalibFitCh315->GetParameter(0);
  data->SA[1]=CalibFitCh315->GetParError(0);

  data->B[0]=CalibFitCh59->GetParameter(1);
  data->SB[0]=CalibFitCh59->GetParError(1);
  data->B[1]=CalibFitCh315->GetParameter(1);
  data->SB[1]=CalibFitCh315->GetParError(1);
  
}


int main(int argc, char* argv[]){
  bool PlotDataset = true;
  Int_t hadd=std::atoi(argv[1]);

  gSystem->Exec("ls Test*.root > listfiletemp.txt");
  gSystem->Exec("mkdir Plot");
  gSystem->Exec("mkdir Plot/CalibVsTempFit");
  gSystem->Exec("mkdir Plot/ControlPlotCorrected");

  std::vector<std::string> FileListTemp;
  FileListTemp=ReadData("listfiletemp.txt");
  int NFileTemp = FileListTemp.size();
  
  std::vector<std::string> FileListOrigin;
  
  for(int i=0; i<NFileTemp-1;i+=2){
    if(hadd==1)gSystem->Exec(("hadd -f File"+std::to_string(i)+".root "+FileListTemp[i]+" "+FileListTemp[i+1]).c_str());
    FileListOrigin.push_back(FileListTemp[i]); 
  }
  
  gSystem->Exec("ls File*.root > filelist.txt");

  std::vector<std::string> FileList;
  FileList=ReadData("filelist.txt");
  int NFile = FileList.size();

  if(PlotDataset){
    
    for(int i=0; i<NFile;i++){
      std::cout << FileList.at(i) << std::endl;
    }//chiudo for
    
  }//chiudo if
  
  TFile* PlotFile[NFile];
  
  TGraphErrors* GraphCh59LTP1[NFile];
  TGraphErrors* GraphCh59LTP2[NFile];
  TGraphErrors* GraphCh315LTP1[NFile];
  TGraphErrors* GraphCh315LTP2[NFile];
  
  TGraphErrors* GraphCh59GTP1[NFile];
  TGraphErrors* GraphCh59GTP2[NFile];
  TGraphErrors* GraphCh315GTP1[NFile];
  TGraphErrors* GraphCh315GTP2[NFile];

  TGraphErrors* GraphBetaLYCh59[NFile];
  TGraphErrors* GraphBetaLYCh315[NFile];

  TGraphErrors* GraphSatValCh59[NFile];
  TGraphErrors* GraphSatValCh315[NFile];
 
  TGraphErrors* GraphTResVsTemp[NFile];

  TMultiGraph* Ch59P1 = new TMultiGraph();
  TMultiGraph* Ch59P2 = new TMultiGraph();
  TMultiGraph* Ch315P1 = new TMultiGraph();
  TMultiGraph* Ch315P2 = new TMultiGraph();
  
  TMultiGraph* LYVsTempCh59 = new TMultiGraph();
  TMultiGraph* LYVsTempCh315 = new TMultiGraph();

  TMultiGraph* SatValVsTempCh59 = new TMultiGraph();
  TMultiGraph* SatValVsTempCh315 = new TMultiGraph();

  for(int i=0;i<NFile;i++){
    
    PlotFile[i]= TFile::Open( (FileList.at(i)).c_str() );
    
    GraphCh59LTP1[i] = (TGraphErrors*)PlotFile[i]->Get("EnergyVsTempCh59P1"); 
    GraphCh59LTP2[i] = (TGraphErrors*)PlotFile[i]->Get("EnergyVsTempCh59P2"); 
    GraphCh315LTP1[i] = (TGraphErrors*)PlotFile[i]->Get("EnergyVsTempCh315P1"); 
    GraphCh315LTP2[i] = (TGraphErrors*)PlotFile[i]->Get("EnergyVsTempCh315P2"); 
    
    GraphCh59GTP1[i] = (TGraphErrors*)PlotFile[i]->Get("EnergyVsBoxTempCh59P1"); 
    GraphCh59GTP2[i] = (TGraphErrors*)PlotFile[i]->Get("EnergyVsBoxTempCh59P2"); 
    GraphCh315GTP1[i] = (TGraphErrors*)PlotFile[i]->Get("EnergyVsBoxTempCh315P1"); 
    GraphCh315GTP2[i] = (TGraphErrors*)PlotFile[i]->Get("EnergyVsBoxTempCh315P2"); 
    
    GraphTResVsTemp[i] = (TGraphErrors*)PlotFile[i]->Get("GraphTimeResolutionVsTemp");
       
    GraphBetaLYCh59[i]=(TGraphErrors*)PlotFile[i]->Get("LYValVsMeanTempCh59");
    GraphBetaLYCh315[i]=(TGraphErrors*)PlotFile[i]->Get("LYValVsMeanTempCh315");
  
    GraphSatValCh59[i]=(TGraphErrors*)PlotFile[i]->Get("SatValVsMeanTempCh59");
    GraphSatValCh315[i]=(TGraphErrors*)PlotFile[i]->Get("SatValVsMeanTempCh315");
  }
  
  GraphCh59LTP1[0]->SetTitle("Ch59LocalTempP1");
  GraphCh59LTP2[0]->SetTitle("Ch59LocalTempP2");
  GraphCh315LTP1[0]->SetTitle("Ch315LocalTempP1");
  GraphCh315LTP2[0]->SetTitle("Ch315LocalTempP2");

  GraphCh59GTP1[0]->SetTitle("Ch59GlobalTempP1");
  GraphCh59GTP2[0]->SetTitle("Ch59GlobalTempP2");
  GraphCh315GTP1[0]->SetTitle("Ch315GlobalTempP1");
  GraphCh315GTP2[0]->SetTitle("Ch315GlobalTempP2");
  
  /////////////////////////////////////////////CH59//////////////////////////////////////////
  TCanvas* CanvasCh59 = new TCanvas("CanvasCh59LT","CanvasCh59LT",1500,800);
  CanvasCh59->Divide(2,1);
  TLegend* LegendCh59 = new TLegend(0.2, 0.2, .2, .2);
  
  CanvasCh59->cd(1);
  GraphCh59LTP1[0]->GetXaxis()->SetLimits(26,36);
  GraphCh59LTP1[0]->Draw("AP");

  CanvasCh59->cd(2);
  GraphCh59GTP1[0]->GetXaxis()->SetLimits(17,27);
  GraphCh59GTP1[0]->Draw("AP");

  for(int i=0;i<NFile;i++){
    CanvasCh59->cd(1);
    GraphCh59LTP1[i]->SetMarkerStyle(24+i);
    GraphCh59LTP2[i]->SetMarkerStyle(24+i);
    LegendCh59->AddEntry(GraphCh59LTP1[i],(FileListOrigin.at(i).erase(FileListOrigin.at(i).size()-3,FileListOrigin.at(i).size())).c_str());
    GraphCh59LTP1[i]->Draw("SAMEP");
    GraphCh59LTP2[i]->Draw("SAMEP");

    CanvasCh59->cd(2);
    GraphCh59GTP1[i]->SetMarkerStyle(24+i);
    GraphCh59GTP2[i]->SetMarkerStyle(24+i);
    GraphCh59GTP1[i]->Draw("SAMEP");
    Ch59P1->Add(GraphCh59GTP1[i]);
    GraphCh59GTP2[i]->Draw("SAMEP");
    Ch59P2->Add(GraphCh59GTP2[i]);
    
    GraphBetaLYCh59[i]->SetMarkerStyle(24+i);
    LYVsTempCh59->Add(GraphBetaLYCh59[i]);
    
    GraphSatValCh59[i]->SetMarkerStyle(24+i);
    SatValVsTempCh59->Add(GraphSatValCh59[i]);
  }
  
  CanvasCh59->cd(1);
  LegendCh59->Draw("SAME");
  CanvasCh59->cd(2);
  LegendCh59->Draw("SAME");
  
  CanvasCh59->SaveAs("Plot/GlobalPlotCh59.png");
  delete LegendCh59;
  delete CanvasCh59;

  ///////////////////////////////////CH315////////////////////////////////////
  TCanvas* CanvasCh315 = new TCanvas("CanvasCh315LT","CanvasCh315LT",1500,800);
  CanvasCh315->Divide(2,1);
  TLegend* LegendCh315 = new TLegend(0.2, 0.2, .2, .2);
  
  CanvasCh315->cd(1);
  GraphCh315LTP1[0]->GetXaxis()->SetLimits(26,36);
  GraphCh315LTP1[0]->Draw("AP");
  
  CanvasCh315->cd(2);
  GraphCh315GTP1[0]->GetXaxis()->SetLimits(17,27);
  GraphCh315GTP1[0]->Draw("AP");
 
  for(int i=0;i<NFile;i++){
  
    CanvasCh315->cd(1);
    GraphCh315LTP1[i]->SetMarkerStyle(24+i);
    GraphCh315LTP2[i]->SetMarkerStyle(24+i);
    LegendCh315->AddEntry(GraphCh315LTP1[i],(FileListOrigin.at(i).erase(FileListOrigin.at(i).size()-3,FileListOrigin.at(i).size())).c_str());
    GraphCh315LTP1[i]->Draw("SAMEP");
    GraphCh315LTP2[i]->Draw("SAMEP");

    CanvasCh315->cd(2);
    GraphCh315GTP1[i]->SetMarkerStyle(24+i);
    GraphCh315GTP2[i]->SetMarkerStyle(24+i);
    GraphCh315GTP1[i]->Draw("SAMEP");
    Ch315P1->Add(GraphCh315GTP1[i]);
    GraphCh315GTP2[i]->Draw("SAMEP");
    Ch315P2->Add(GraphCh315GTP2[i]);
 
    GraphBetaLYCh315[i]->SetMarkerStyle(24+i);
    LYVsTempCh315->Add(GraphBetaLYCh315[i]);

    GraphSatValCh315[i]->SetMarkerStyle(24+i);
    SatValVsTempCh315->Add(GraphSatValCh315[i]);
  }
  
  CanvasCh315->cd(1);
  LegendCh315->Draw("SAME");
  CanvasCh315->cd(2);
  LegendCh315->Draw("SAME");
  
  CanvasCh315->SaveAs("Plot/GlobalPlotCh315.png");
  delete LegendCh315;
  delete CanvasCh315;

  TCanvas* CanvasTimeResVsTemp = new TCanvas("TimeResTot","TimeResTot",1600,800);
  TLegend* LegendTResVsTemp = new TLegend(0.2, 0.2, .2, .2);
  Double_t Mean=0;
  Double_t RMSVal=0;
  int count=0;
  std::vector<Double_t> YVal;
  
  Mean+=GraphTResVsTemp[0]->GetMean(2);
  for(int j=0;j< GraphTResVsTemp[0]->GetN();j++){
    YVal.push_back(GraphTResVsTemp[0]->GetY()[j]);
    count++;
  }
  
  GraphTResVsTemp[0]->GetXaxis()->SetLimits(21,26.5);
  GraphTResVsTemp[0]->SetMarkerStyle(24);
  GraphTResVsTemp[0]->SetName(GraphTResVsTemp[0]->GetTitle());
  GraphTResVsTemp[0]->Draw("AP");
  LegendTResVsTemp->AddEntry(GraphTResVsTemp[0],(FileListOrigin.at(0)).c_str());
			    
  for(int i=1; i< NFile;i++){
    for(int j=0;j< GraphTResVsTemp[i]->GetN();j++){
      YVal.push_back(GraphTResVsTemp[i]->GetY()[j]);
      count++;
    }
    Mean+=GraphTResVsTemp[i]->GetMean(2);
    GraphTResVsTemp[i]->SetMarkerStyle(24+i);
    GraphTResVsTemp[i]->SetName(GraphTResVsTemp[i]->GetTitle());
    GraphTResVsTemp[i]->Draw("SAMEP");
    LegendTResVsTemp->AddEntry(GraphTResVsTemp[i],(FileListOrigin.at(i)).c_str());
  }
  
  RMSVal=TMath::RMS(YVal.begin(),YVal.end());
  LegendTResVsTemp->SetHeader(Form("Mean: %lf  MeanError: %lf",Mean/NFile,RMSVal/sqrt(count)));
  LegendTResVsTemp->Draw("SAME");
    
  CanvasTimeResVsTemp->SaveAs("Plot/TimeResVsTempTot.png");
  delete CanvasTimeResVsTemp;
  delete LegendTResVsTemp;
  
  //Light yeld at fixed Temperature with fit
  ////////////////////////////////////////////////////////////////////////////////FIT
  gStyle->SetOptFit(1111);

  TF1* FitTotCh59P1 = new TF1("FitTotCh59P1","[0]*x+[1]", 20., 30.);
  TF1* FitTotCh315P1 = new TF1("FitTotCh315P1","[0]*x+[1]", 20., 30.);
  TF1* FitTotCh59P2 = new TF1("FitTotCh59P2","[0]*x+[1]", 20., 30.);
  TF1* FitTotCh315P2 = new TF1("FitTotCh315P2","[0]*x+[1]", 20., 30.);
 

  TCanvas* CanvasFitTotal = new TCanvas("CanvasFitTotal","CanvasFitTotal",1600,800);
  CanvasFitTotal->Divide(2,2);
  
  CanvasFitTotal->cd(1);
  Ch59P1->SetMinimum(36);
  Ch59P1->SetMaximum(44);
  Ch59P1->SetTitle("Ch59P511KeV "+FitTotCh59P1->GetExpFormula());
  Ch59P1->Draw("AP");
  Ch59P1->Fit(FitTotCh59P1, "R");
  
  TVirtualFitter* fitterCh59P1 = TVirtualFitter::GetFitter();
  Double_t* CMCh59P1 = fitterCh59P1->GetCovarianceMatrix();
  std::cout <<"CorrelationElements " << CMCh59P1[0] << " " <<CMCh59P1[1] << " " << CMCh59P1[2] << std::endl;
  
  std::pair<TF1*,Double_t> FullCh59P1;
  FullCh59P1.first=FitTotCh59P1;
  FullCh59P1.second=CMCh59P1[1];


  CanvasFitTotal->cd(2);
  Ch315P1->SetMinimum(36);
  Ch315P1->SetMaximum(44);
  Ch315P1->SetTitle("Ch315P511KeV");
  Ch315P1->Draw("AP");
  Ch315P1->Fit(FitTotCh315P1);
  
  TVirtualFitter* fitterCh315P1 = TVirtualFitter::GetFitter();
  Double_t* CMCh315P1 = fitterCh315P1->GetCovarianceMatrix();  
  std::cout <<"CorrelationElements "  << CMCh315P1[0] << " " <<CMCh315P1[1] << " " << CMCh315P1[2] << std::endl;

  std::pair<TF1*,Double_t> FullCh315P1;
  FullCh315P1.first=FitTotCh315P1;
  FullCh315P1.second=CMCh315P1[1];


  CanvasFitTotal->cd(3);
  Ch59P2->SetMinimum(74);
  Ch59P2->SetMaximum(84);
  Ch59P2->SetTitle("Ch59P1275KeV");
  Ch59P2->Draw("AP");
  Ch59P2->Fit(FitTotCh59P2);
  TVirtualFitter* fitterCh59P2 = TVirtualFitter::GetFitter();
  Double_t* CMCh59P2 = fitterCh59P2->GetCovarianceMatrix();
  std::cout <<"CorrelationElements "  << CMCh59P2[0] << " " <<CMCh59P2[1] << " " << CMCh59P2[2] << std::endl;

  std::pair<TF1*,Double_t> FullCh59P2;
  FullCh59P2.first=FitTotCh59P2;
  FullCh59P2.second=CMCh59P2[1];

  
  CanvasFitTotal->cd(4);
  Ch315P2->SetMinimum(74);
  Ch315P2->SetMaximum(84);
  Ch315P2->SetTitle("Ch315P1275KeV");
  Ch315P2->Draw("AP");
  Ch315P2->Fit(FitTotCh315P2);
  TVirtualFitter* fitterCh315P2 = TVirtualFitter::GetFitter();
  Double_t* CMCh315P2 = fitterCh315P2->GetCovarianceMatrix();
  std::cout <<"CorrelationElements "  << CMCh315P2[0] << " " <<CMCh315P2[1] << " " << CMCh315P2[2] << std::endl;

  std::pair<TF1*,Double_t> FullCh315P2;
  FullCh315P2.first=FitTotCh315P2;
  FullCh315P2.second=CMCh315P2[1];

  ResidualPlot(Ch59P1,FitTotCh59P1,"FitTotalShiftPeak59P1");
  ResidualPlot(Ch59P2,FitTotCh59P2,"FitTotalShiftPeak59P2");
  ResidualPlot(Ch315P1,FitTotCh315P1,"FitTotalShiftPeak315P1");
  ResidualPlot(Ch315P2,FitTotCh315P2,"FitTotalShiftPeak315P2");

  CanvasFitTotal->SaveAs("Plot/FitTotalShiftPeak.png");
  
  /////////////////////////////////////////////////////////////////////////////////CalibCurve
  
  const Double_t Tmin=21;
  const Double_t Tmax=27;
  const int NPlot = (int)(Tmax-Tmin);
  
  Double_t GTemp[NPlot];
  Double_t AVal[2][NPlot];
  Double_t BVal[2][NPlot];
  Double_t SAVal[2][NPlot];
  Double_t SBVal[2][NPlot];

  dataPlot data;

  for(int i=0; i<NPlot; i++){

    LYAtFixedTemp(i+Tmin, FullCh59P1, FullCh59P2, FullCh315P1, FullCh315P2, &data);         
    GTemp[i]=i+Tmin;
   
    AVal[0][i]=data.A[0];
    AVal[1][i]=data.A[1];
    
    SAVal[0][i]=data.SA[0];
    SAVal[1][i]=data.SA[1];
    
    BVal[0][i]=data.B[0];
    BVal[1][i]=data.B[1];

    SBVal[0][i]=data.SB[0];
    SBVal[1][i]=data.SB[1];
    
    std::cout << GTemp[i] <<"\t"<< AVal[0][i]<<"\t" <<SAVal[0][i]<< "\t" <<AVal[1][i]<< "\t" <<SAVal[1][i]<< "\t" <<BVal[0][i] << "\t" << SBVal[0][i]     << "\t" << BVal[1][i] << "\t" << SBVal[1][i] << std::endl;
  }

  
  TGraphErrors* AValue[2];
  AValue[0]= new TGraphErrors(NPlot,GTemp,AVal[0],0,SAVal[0]);
  AValue[0]->SetTitle("SaturationValueCh59");
  AValue[0]->SetMarkerStyle(8);
  AValue[0]->SetName(AValue[0]->GetTitle());

  AValue[1]= new TGraphErrors(NPlot,GTemp,AVal[1],0,SAVal[1]);
  AValue[1]->SetTitle("SaturationValueCh315");
  AValue[1]->SetMarkerStyle(8);  
  AValue[1]->SetName(AValue[1]->GetTitle());

  TGraphErrors* BValue[2];
  BValue[0]= new TGraphErrors(NPlot,GTemp,BVal[0],0,SBVal[0]);
  BValue[0]->SetTitle("LYValueCh59");
  BValue[0]->SetMarkerStyle(8);
  BValue[0]->SetName(BValue[0]->GetTitle());

  BValue[1]= new TGraphErrors(NPlot,GTemp,BVal[1],0,SBVal[1]);
  BValue[1]->SetTitle("LYValueCh315");
  BValue[1]->SetMarkerStyle(8);
  BValue[1]->SetName(BValue[1]->GetTitle());
  
  TCanvas* AandBValues=new TCanvas("AandBValues","AandBValues",1600,800);
  AandBValues->Divide(2,2);
  AandBValues->cd(1);
  AValue[0]->Draw("AP");
  AandBValues->cd(2);
  AValue[1]->Draw("AP");
  AandBValues->cd(3);
  BValue[0]->Draw("AP");
  AandBValues->cd(4);
  BValue[1]->Draw("AP");
  
  AandBValues->SaveAs("Plot/ABValuesWithFit.png");
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  Double_t T0,TRef;
  T0=24;
  TRef= *(GraphCh59GTP2[0]->GetX()+100%GraphCh59GTP2[0]->GetN()); //Take a quasi-random value
  

  TFormula* CorrCh59[2]; //So called correction functionn
  TFormula* CorrCh315[2]; 

  CorrCh59[0]= new TFormula("CorrCh59Peak1",("1/("+(std::string)(FitTotCh59P1->GetExpFormula())+")*"+std::to_string(FitTotCh59P1->Eval(T0))).c_str());
  CorrCh59[0]->SetParameter(0,FitTotCh59P1->GetParameter(0));
  CorrCh59[0]->SetParameter(1,FitTotCh59P1->GetParameter(1));
  std::cout << CorrCh59[0]->GetExpFormula()<< " "<< CorrCh59[0]->GetParameter(0) << " "<<CorrCh59[0]->GetParameter(1) <<std::endl;
    
  CorrCh59[1]= new TFormula("CorrCh59Peak2",("1/("+(std::string)(FitTotCh59P2->GetExpFormula())+")*"+std::to_string(FitTotCh59P2->Eval(T0))).c_str());
  CorrCh59[1]->SetParameter(0,FitTotCh59P2->GetParameter(0));
  CorrCh59[1]->SetParameter(1,FitTotCh59P2->GetParameter(1));
  std::cout << CorrCh59[1]->GetExpFormula()<< " "<< CorrCh59[1]->GetParameter(0) << " "<<CorrCh59[1]->GetParameter(1) <<std::endl;
  
  CorrCh315[0] = new TFormula("CorrCh315Peak1",("1/("+(std::string)(FitTotCh315P1->GetExpFormula())+")*"+std::to_string(FitTotCh315P1->Eval(T0))).c_str());
  CorrCh315[0]->SetParameter(0,FitTotCh315P1->GetParameter(0));
  CorrCh315[0]->SetParameter(1,FitTotCh315P1->GetParameter(1));
  std::cout << CorrCh315[0]->GetExpFormula()<< " "<< CorrCh315[0]->GetParameter(0) << " "<<CorrCh315[0]->GetParameter(1) <<std::endl;

  CorrCh315[1] = new TFormula("CorrCh315Peak2",("1/("+(std::string)(FitTotCh315P2->GetExpFormula())+")*"+std::to_string(FitTotCh315P2->Eval(T0))).c_str());
  CorrCh315[1]->SetParameter(0,FitTotCh315P2->GetParameter(0));
  CorrCh315[1]->SetParameter(1,FitTotCh315P2->GetParameter(1));
  std::cout << CorrCh315[1]->GetExpFormula()<< " "<< CorrCh315[1]->GetParameter(0) << " "<<CorrCh315[1]->GetParameter(1) <<std::endl;

  std::cout << TRef << std::endl;
  //////////////////////////////////////////////////// AT T_REF
  Double_t ADCRefT0Ch59[2]={*(GraphCh59GTP1[0]->GetY()+100%GraphCh59GTP1[0]->GetN()) * CorrCh59[0]->Eval(TRef) , *(GraphCh59GTP2[0]->GetY()+100%GraphCh59GTP2[0]->GetN()) * CorrCh59[1]->Eval(TRef)};

  Double_t ADCRefT0Ch315[2]={*(GraphCh315GTP1[0]->GetY()+100%GraphCh315GTP1[0]->GetN()) * CorrCh315[0]->Eval(TRef), *(GraphCh315GTP2[0]->GetY()+100%GraphCh315GTP2[0]->GetN()) * CorrCh315[0]->Eval(TRef)};
  
  std::vector<Double_t> LYRelT0Ch59P1;
  std::vector<Double_t> LYRelT0Ch59P2;
  std::vector<Double_t> LYRelT0Ch315P1;
  std::vector<Double_t> LYRelT0Ch315P2;
    
  std::vector<Double_t> BValCh59;
  std::vector<Double_t> BValCh315;  
  std::vector<Double_t> SigmaBValCh59;
  std::vector<Double_t> SigmaBValCh315;  

  std::vector<Double_t> GlobalTemperature;

  TGraphErrors* CalibCh59Corr;
  TGraphErrors* CalibCh315Corr;
  
  TF1* FitCalibCh59Corr;
  TF1* FitCalibCh315Corr;
  
  TCanvas* CalibCanvasPlotCorr;

  for(int i=0; i<NFile ; i++){
    for(int j=0; j<GraphCh59GTP1[i]->GetN();j++){

      CalibCanvasPlotCorr = new TCanvas("CalibCanvasPlotCorr","CalibCanvasPlotCorr",1500,700);
      CalibCanvasPlotCorr->Divide(2,1);

      CalibCanvasPlotCorr->cd(1);
      CalibCh59Corr = new TGraphErrors();
      CalibCh59Corr->SetName(("CalibPlotCh59"+std::to_string(i)+std::to_string(j)).c_str());
      CalibCh59Corr->SetPoint(0,511,*(GraphCh59GTP1[i]->GetY()+j)  *  CorrCh59[0]->Eval(  *(GraphCh59GTP1[i]->GetX()+j) ) );
      CalibCh59Corr->SetPoint(1,1275, *(GraphCh59GTP2[i]->GetY()+j)  *  CorrCh59[1]->Eval(  *(GraphCh59GTP2[i]->GetX()+j)  ) );
      //std::cout <<"CHECK_____________"<< *(GraphCh59GTP2[i]->GetY()+j)  *  CorrCh59[1]->Eval(  *(GraphCh59GTP2[i]->GetX()+j)  ) << std::endl;
      FitCalibCh59Corr=CalibrationCurve(CalibCh59Corr,"Ch59");
      BValCh59.push_back(FitCalibCh59Corr->GetParameter(1));
      SigmaBValCh59.push_back(FitCalibCh59Corr->GetParError(1));
      //std::cout << FitCalibCh59Corr->GetParError(1) << std::endl;
      CalibCh59Corr->SetMarkerStyle(8);
      CalibCh59Corr->GetXaxis()->SetLimits(0,1700);
      CalibCh59Corr->GetXaxis()->SetRangeUser(0,1700);
      CalibCh59Corr->GetXaxis()->SetTitle("E [KeV]");
      CalibCh59Corr->GetYaxis()->SetTitle("ADC [D.U.]");
      CalibCh59Corr->GetYaxis()->SetLimits(0,110);
      CalibCh59Corr->GetYaxis()->SetRangeUser(0,110);
      CalibCh59Corr->Draw("AP");
      FitCalibCh59Corr->Draw("SAME");

      GlobalTemperature.push_back(*(GraphCh59GTP1[i]->GetX()+j));

      LYRelT0Ch59P1.push_back(  *(GraphCh59GTP1[i]->GetY()+j)  *  CorrCh59[0]->Eval(  *(GraphCh59GTP1[i]->GetX()+j)  ) /  (  ADCRefT0Ch59[0] )        );
      LYRelT0Ch59P2.push_back(  *(GraphCh59GTP2[i]->GetY()+j)  *  CorrCh59[1]->Eval(  *(GraphCh59GTP2[i]->GetX()+j)  ) /  (  ADCRefT0Ch59[1] )        );

      CalibCanvasPlotCorr->cd(2);
      CalibCh315Corr = new TGraphErrors();
      CalibCh59Corr->SetName(("CalibPlotCh315"+std::to_string(i)+std::to_string(j)).c_str());
      CalibCh315Corr->SetPoint(0,511,*(GraphCh315GTP1[i]->GetY()+j)  *  CorrCh315[0]->Eval(  *(GraphCh315GTP1[i]->GetX()+j) ) );
      CalibCh315Corr->SetPoint(1,1275, *(GraphCh315GTP2[i]->GetY()+j)  *  CorrCh315[1]->Eval(  *(GraphCh315GTP2[i]->GetX()+j)  ) );
      //std::cout << "CHECK____________"<< *(GraphCh315GTP2[i]->GetY()+j)  *  CorrCh315[1]->Eval(  *(GraphCh315GTP2[i]->GetX()+j)  ) << std::endl;
      FitCalibCh315Corr=CalibrationCurve(CalibCh315Corr,"Ch315");
      BValCh315.push_back(FitCalibCh315Corr->GetParameter(1)); 
      SigmaBValCh315.push_back(FitCalibCh315Corr->GetParError(1));
      //std::cout << FitCalibCh315Corr->GetParError(1) << std::endl;     
      CalibCh315Corr->SetMarkerStyle(8);
      CalibCh315Corr->GetXaxis()->SetLimits(0,1700);
      CalibCh315Corr->GetXaxis()->SetRangeUser(0,1700);
      CalibCh315Corr->GetXaxis()->SetTitle("E [KeV]");
      CalibCh315Corr->GetYaxis()->SetTitle("ADC [D.U.]");
      CalibCh315Corr->GetYaxis()->SetLimits(0,110);
      CalibCh315Corr->GetYaxis()->SetRangeUser(0,110);
      CalibCh315Corr->Draw("AP");
      FitCalibCh315Corr->Draw("SAME");

      LYRelT0Ch315P1.push_back(  *(GraphCh315GTP1[i]->GetY()+j)  *  CorrCh315[0]->Eval(  *(GraphCh315GTP1[i]->GetX()+j)  ) /  (  ADCRefT0Ch59[0] )        );
      LYRelT0Ch315P2.push_back(  *(GraphCh315GTP2[i]->GetY()+j)  *  CorrCh315[1]->Eval(  *(GraphCh315GTP2[i]->GetX()+j)  ) /  (  ADCRefT0Ch59[1] )        );
      
      std::cout << LYRelT0Ch59P1.at(j) <<" " << LYRelT0Ch59P2.at(j) <<" "<< LYRelT0Ch315P1.at(j) <<" "<<LYRelT0Ch315P2.at(j) <<" "<< std::endl;
      
      CalibCanvasPlotCorr->SaveAs(("Plot/ControlPlotCorrected/CalibCurve"+std::to_string(i)+std::to_string(j)+".png").c_str());
      
      delete CalibCanvasPlotCorr;
      delete FitCalibCh59Corr;
      delete FitCalibCh315Corr;   
      delete CalibCh59Corr;
      delete CalibCh315Corr;
    }//chiudo for j
  }//chiudo for i
  
  
  TH1D* HistoLYRelT0Ch59[2];
  TH1D* HistoLYRelT0Ch315[2];

  Double_t minTH=0.95;
  Double_t maxTH= 1.05;
  Int_t NBin=70;
  
  HistoLYRelT0Ch59[0]= new TH1D("HistoLYRelT0Ch59P1","HistoLYRelT0Ch59P1",NBin,minTH,maxTH);
  HistoLYRelT0Ch59[1]= new TH1D("HistoLYRelT0Ch59P2","HistoLYRelT0Ch59P2",NBin,minTH,maxTH);

  HistoLYRelT0Ch315[0]= new TH1D("HistoLYRelT0Ch315P1","HistoLYRelT0Ch315P1",NBin,minTH,maxTH);
  HistoLYRelT0Ch315[1]= new TH1D("HistoLYRelT0Ch315P2","HistoLYRelT0Ch315P2",NBin,minTH,maxTH);

  for(int i=0; i<(int)LYRelT0Ch315P1.size();i++){
    
    HistoLYRelT0Ch59[0]->Fill(LYRelT0Ch59P1.at(i));
    HistoLYRelT0Ch59[1]->Fill(LYRelT0Ch59P2.at(i));
    HistoLYRelT0Ch315[0]->Fill(LYRelT0Ch315P1.at(i));
    HistoLYRelT0Ch315[1]->Fill(LYRelT0Ch315P2.at(i));
  }
  
  TF1* FitLYResCh59P1 = new TF1("FitLYResCh59P1","gaus", 0.5, 1.5);
  TF1* FitLYResCh315P1 = new TF1("FitLYResCh315P1","gaus", 0.5,1.5);
  TF1* FitLYResCh59P2 = new TF1("FitLYResCh59P2","gaus", 0.5,1.5);
  TF1* FitLYResCh315P2 = new TF1("FitLYResCh315P2","gaus", 0.5,1.5);
 

  TCanvas* CanvasLYRelT0 = new TCanvas("CanvasLYRelT0","CanvasLYRelT0",1400,700);
  CanvasLYRelT0->Divide(2,2);
  CanvasLYRelT0->cd(1);
  HistoLYRelT0Ch59[0]->Draw();
  HistoLYRelT0Ch59[0]->Fit(FitLYResCh59P1,"Q0");
  FitLYResCh59P1->SetLineWidth(2);
  FitLYResCh59P1->Draw("SAME");
  
  CanvasLYRelT0->cd(2);
  HistoLYRelT0Ch315[0]->Draw();
  HistoLYRelT0Ch315[0]->Fit(FitLYResCh315P1,"Q0");
  FitLYResCh315P1->SetLineWidth(2);
  FitLYResCh315P1->Draw("SAME");
  
  CanvasLYRelT0->cd(3);
  HistoLYRelT0Ch59[1]->Draw();
  HistoLYRelT0Ch59[1]->Fit(FitLYResCh59P2,"Q0");
  FitLYResCh59P2->SetLineWidth(2);
  FitLYResCh59P2->Draw("SAME");

  CanvasLYRelT0->cd(4);
  HistoLYRelT0Ch315[1]->Draw();
  HistoLYRelT0Ch315[1]->Fit(FitLYResCh315P2,"Q0");
  FitLYResCh315P2->SetLineWidth(2);
  FitLYResCh59P2->Draw("SAME");

  CanvasLYRelT0->SaveAs("Plot/LYRelT0AllChannel.png");

  TCanvas* CanvasBValuCorrected = new TCanvas("CanvasBValuCorrected", "CanvasBValuCorrected",1600,800);
  CanvasBValuCorrected->Divide(2,1);
  
  CanvasBValuCorrected->cd(1);
  
  
  TGraphErrors* GraphBValueCorrCh59 = new TGraphErrors((int)BValCh59.size(),&GlobalTemperature[0],&BValCh59[0],0,&SigmaBValCh59[0]);
  TLegend* legCh59 = new TLegend();
  legCh59->AddEntry(GraphBValueCorrCh59,Form("Mean: %lf",TMath::Mean(BValCh59.size(),&BValCh59[0])));
  legCh59->AddEntry(GraphBValueCorrCh59,Form("Mean: %lf",TMath::RMS(BValCh59.size(),&BValCh59[0])));
 
  GraphBValueCorrCh59->SetMarkerStyle(8);
  GraphBValueCorrCh59->GetXaxis()->SetTitle("Temp [°C]");
  GraphBValueCorrCh59->GetYaxis()->SetTitle("FitValue");  
  GraphBValueCorrCh59->Draw("AP");
  legCh59->Draw("SAME");

  CanvasBValuCorrected->cd(2);
  TGraphErrors* GraphBValueCorrCh315 = new TGraphErrors((int)BValCh315.size(),&GlobalTemperature[0],&BValCh315[0],0,&SigmaBValCh315[0]);
  TLegend* legCh315 = new TLegend();
  legCh315->AddEntry(GraphBValueCorrCh315,Form("Mean: %lf",TMath::Mean(BValCh315.size(),&BValCh315[0])));
  legCh315->AddEntry(GraphBValueCorrCh315,Form("Mean: %lf",TMath::RMS(BValCh315.size(),&BValCh315[0])));

  GraphBValueCorrCh315->SetMarkerStyle(8);    
  GraphBValueCorrCh315->GetXaxis()->SetTitle("Temp [°C]");
  GraphBValueCorrCh315->GetYaxis()->SetTitle("FitValue");  
  GraphBValueCorrCh315->Draw("AP");
  legCh315->Draw("SAME");
 
  CanvasBValuCorrected->SaveAs("Plot/CavasBValuesCorrected.png");
  



  for(int i=0;i<NFile;i++){
    for(int j=0; j<GraphBetaLYCh315[i]->GetN();j++ ){

      std::cout<<j << " " << *(GraphBetaLYCh59[i]->GetY()+j) << std::endl;
      if(*(GraphBetaLYCh59[i]->GetY()+j)==(Double_t)0){
	GraphBetaLYCh59[i]->RemovePoint(j);
	std::cout << "Ch59:LY-> point "<< j<< " removed: value "<< *(GraphBetaLYCh59[i]->GetY()+j) << std::endl;   
      }//chiudo if

      std::cout<< j<< " " << *(GraphBetaLYCh315[i]->GetY()+j) << std::endl; 
      if(*(GraphBetaLYCh315[i]->GetY()+j)==(Double_t)0){
	std::cout << "Ch315:LY-> point "<< j << " removed: value "<< *(GraphBetaLYCh315[i]->GetY()+j) << std::endl;
	GraphBetaLYCh315[i]->RemovePoint(j);
      }//chiudo if
      
      std::cout<< j<< " " << *(GraphSatValCh59[i]->GetY()+j) << std::endl; 
      if(*(GraphSatValCh59[i]->GetY()+j)==(Double_t)0){
	std::cout << "Ch59:SV-> point "<< j << " removed: value "<< *(GraphSatValCh59[i]->GetY()+j) << std::endl;
	GraphSatValCh59[i]->RemovePoint(j);
      }//chiudo if

      std::cout<< j<< " " << *(GraphSatValCh315[i]->GetY()+j) << std::endl; 
      if(*(GraphSatValCh315[i]->GetY()+j)==(Double_t)0){
	std::cout << "Ch315:SV-> point "<< j << " removed: value "<< *(GraphSatValCh315[i]->GetY()+j) << std::endl;
	GraphSatValCh315[i]->RemovePoint(j);
      }//chiudo if
      
      
      
    }//chiudo for j
  }//chiudo for i
  
  //GraphSatValCh59[i]

  TF1* FitLYCh59= new TF1("FitLYCh59","pol1");
  TF1* FitLYCh315= new TF1("FitLYCh315","pol1");
  
  FitLYCh59->SetParameter(1,-1.4e-5);
  FitLYCh315->SetParameter(1,-1.4e-5);

  TCanvas* CanvasLYUncorr = new TCanvas("CanvasLYUncorr","CanvasLYUncorr",1400,700);
  CanvasLYUncorr->Divide(2,1);
  
  CanvasLYUncorr->cd(1);
  LYVsTempCh59->SetTitle("LYVsTempCh59Total");
  LYVsTempCh59->SetMinimum(0.3e-3);
  LYVsTempCh59->SetMaximum(0.9e-3);
  LYVsTempCh59->Draw("AP");
  LYVsTempCh59->Fit(FitLYCh59,"FM");

  CanvasLYUncorr->cd(2);
  LYVsTempCh315->SetTitle("LYVsTempCh315Total");
  LYVsTempCh315->SetMinimum(0.3e-3);
  LYVsTempCh315->SetMaximum(0.9e-3);
  LYVsTempCh315->Draw("AP");        
  LYVsTempCh315->Fit(FitLYCh315,"FM");

  ResidualPlot(LYVsTempCh59,FitLYCh59,"LYUncorrTot59");
  ResidualPlot(LYVsTempCh315,FitLYCh315,"LYUncorrTot315");
  CanvasLYUncorr->SaveAs("Plot/LYUncorrTot.png");
  
  
  TF1* FitSatCh59= new TF1("FitSatCh59","pol1");
  TF1* FitSatCh315= new TF1("FitSatCh315","pol1");
  
  FitSatCh59->SetParameter(0,250);
  FitSatCh59->SetParameter(1,-4.19);
  FitSatCh315->SetParameter(0,250);
  FitSatCh315->SetParameter(1,-4.19);


  TCanvas* CavasSatVal = new TCanvas("CavasSatVal","CavasSatVal",1400,700);
  CavasSatVal->Divide(2,1);
  
  CavasSatVal->cd(1);
  SatValVsTempCh59->Draw("AP");
  SatValVsTempCh59->Fit(FitSatCh59);
  SatValVsTempCh59->SetMinimum(130);
  SatValVsTempCh59->SetMaximum(180);
  SatValVsTempCh59->Draw("AP");
  
  CavasSatVal->cd(2);
  SatValVsTempCh315->Draw("AP");
  SatValVsTempCh315->Fit(FitSatCh315);
  SatValVsTempCh315->SetMinimum(130);
  SatValVsTempCh315->SetMaximum(180);
  SatValVsTempCh315->Draw("AP");

  CavasSatVal->SaveAs("Plot/SatValUncorrTot.png");    

  ResidualPlot(SatValVsTempCh59,FitSatCh59,"SatValUncorrTot59");
  ResidualPlot(SatValVsTempCh315,FitSatCh315,"SatValUncorrTot315");
 

}//chiudo main

void ResidualPlot(TMultiGraph* graph,TF1* fit, std::string name){
  
  TGraphErrors* temp;
  TList* lista = graph->GetListOfGraphs();
  std::vector<Double_t> vec;
  
  for(int i=0; i< lista->GetSize();i++){
    temp = (TGraphErrors*)(lista->At(i));
    
    for(int j=0; j<temp->GetN(); j++){
      vec.push_back( *(temp->GetY()+j)-fit->Eval(*(temp->GetX()+j)) );
    }//chiudo for j
    
  }//chiudo for i
  
  TH1D* residual = new TH1D(("Res"+name).c_str(),("Res"+name).c_str(),80, TMath::Mean(vec.begin(),vec.end())-5*TMath::RMS(vec.begin(),vec.end()),TMath::Mean(vec.begin(),vec.end())+5*TMath::RMS(vec.begin(),vec.end()));
  
  for(int i=0; i < (int)vec.size();i++){
    
    residual->Fill(vec.at(i));
    
  }//chiudo for i

  TCanvas* CanvasRes = new TCanvas("CanvasRes","CanvasRes",700,700);
  residual->Draw();

  CanvasRes->SaveAs(("Plot/"+name+"Res.png").c_str());

  delete CanvasRes;
  delete residual;
}
