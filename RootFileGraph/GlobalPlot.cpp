#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TSystem.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TF1.h"

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

int main(int argc, char* argv[]){
  bool PlotDataset = true;

  gSystem->Exec("ls Test*.root > listfiletemp.txt");
  gSystem->Exec("mkdir Plot");

   
  std::vector<std::string> FileListTemp;
  FileListTemp=ReadData("listfiletemp.txt");
  int NFileTemp = FileListTemp.size();
  
  std::vector<std::string> FileListOrigin;
  
  for(int i=0; i<NFileTemp-1;i+=2){
    gSystem->Exec(("hadd -f File"+std::to_string(i)+".root "+FileListTemp[i]+" "+FileListTemp[i+1]).c_str());
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

  TGraphErrors* GraphTResVsTemp[NFile];
  TH1D* HistoTimeRes[NFile];
  TF1* FitHistoTimeRes[NFile];

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
    HistoTimeRes[i] = (TH1D*)PlotFile[i]->Get("TimeResTot");
    FitHistoTimeRes[i] = (TF1*)PlotFile[i]->Get("FitTimeResTot");
      
  }
  
  GraphCh59LTP1[1]->SetTitle("Ch59LocalTempP1");
  GraphCh59LTP2[1]->SetTitle("Ch59LocalTempP2");
  GraphCh315LTP1[1]->SetTitle("Ch315LocalTempP1");
  GraphCh315LTP2[1]->SetTitle("Ch315LocalTempP2");

  GraphCh59GTP1[1]->SetTitle("Ch59GlobalTempP1");
  GraphCh59GTP2[1]->SetTitle("Ch59GlobalTempP2");
  GraphCh315GTP1[1]->SetTitle("Ch315GlobalTempP1");
  GraphCh315GTP2[1]->SetTitle("Ch315GlobalTempP2");

  /////////////////////////////////////////////CH59//////////////////////////////////////////
  TCanvas* CanvasCh59 = new TCanvas("CanvasCh59LT","CanvasCh59LT",1500,800);
  CanvasCh59->Divide(2,1);
  TLegend* LegendCh59 = new TLegend(0.2, 0.2, .2, .2);
  
  CanvasCh59->cd(1);
  GraphCh59LTP1[1]->GetXaxis()->SetLimits(26,36);
  GraphCh59LTP1[1]->Draw("AP");
  CanvasCh59->cd(2);
  GraphCh59GTP1[1]->GetXaxis()->SetLimits(17,27);
  GraphCh59GTP1[1]->Draw("AP");

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
    GraphCh59GTP2[i]->Draw("SAMEP");
  }
  
  CanvasCh59->cd(1);
  LegendCh59->Draw("SAME");
  CanvasCh59->cd(2);
  LegendCh59->Draw("SAME");
  
  CanvasCh59->SaveAs("Plot/GlobalPlotCh59.png");

  ///////////////////////////////////CH315////////////////////////////////////
  TCanvas* CanvasCh315 = new TCanvas("CanvasCh315LT","CanvasCh315LT",1500,800);
  CanvasCh315->Divide(2,1);
  TLegend* LegendCh315 = new TLegend(0.2, 0.2, .2, .2);
  
  CanvasCh315->cd(1);
  GraphCh315LTP1[1]->GetXaxis()->SetLimits(26,36);
  GraphCh315LTP1[1]->Draw("AP");
  CanvasCh315->cd(2);
  GraphCh315GTP1[1]->GetXaxis()->SetLimits(17,27);
  GraphCh315GTP1[1]->Draw("AP");
  
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
    GraphCh315GTP2[i]->Draw("SAMEP");
  }
  
  CanvasCh315->cd(1);
  LegendCh315->Draw("SAME");
  CanvasCh315->cd(2);
  LegendCh315->Draw("SAME");
  
  CanvasCh315->SaveAs("Plot/GlobalPlotCh315.png");
  
  TCanvas* CanvasTimeResVsTemp = new TCanvas("TimeResTot","TimeResTot",1600,800);
  TLegend* LegendTResVsTemp = new TLegend(0.2, 0.2, .2, .2);
  
  GraphTResVsTemp[0]->GetXaxis()->SetLimits(17,28);
  GraphTResVsTemp[0]->SetMarkerStyle(24);
  GraphTResVsTemp[0]->SetName(GraphTResVsTemp[0]->GetTitle());
  GraphTResVsTemp[0]->Draw("AP");
  LegendTResVsTemp->AddEntry(GraphTResVsTemp[0],(FileListOrigin.at(0)).c_str());
			     //.erase(FileList.at(0).size()-3,FileList.at(0).size())).c_str()
  for(int i=1; i< NFile;i++){
    GraphTResVsTemp[i]->SetMarkerStyle(24+i);
    GraphTResVsTemp[i]->SetName(GraphTResVsTemp[i]->GetTitle());
    GraphTResVsTemp[i]->Draw("SAMEP");
    LegendTResVsTemp->AddEntry(GraphTResVsTemp[i],(FileListOrigin.at(i)).c_str());
  }
  LegendTResVsTemp->Draw("SAME");
  CanvasTimeResVsTemp->SaveAs("Plot/TimeResVsTempTot.png");
  
  

  /*
  TH1D* TimeResolutionAllRun = new TH1D("TimeResolutionAllRun","TimeResolutionTot",200,-2000,2000);

  for(int i=0; i< NFile;i++){
    TimeResolutionTot->Add(HistoTimeRes[i])
  }
  
  TF1* pippo = new TF1("FitTimeResAllRun")
  */
  
  return 0;  
}
