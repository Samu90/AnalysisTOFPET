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

  gSystem->Exec("ls *.root > listfile.txt");
  gSystem->Exec("mkdir Plot");
 
  std::vector<std::string> FileList;
  FileList=ReadData("listfile.txt");
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
    LegendCh59->AddEntry(GraphCh59LTP1[i],(FileList.at(i).erase(FileList.size()-3,FileList.size())).c_str());
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
    LegendCh315->AddEntry(GraphCh315LTP1[i],(FileList.at(i).erase(FileList.size()-3,FileList.size())).c_str());
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
  
  

}
