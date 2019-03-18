
void GetPedestal(TTree* tree, TH1F* histoCh1, TH1F* histoCh2){
  
  Float_t energy;
  UShort_t chID;

  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("channelID",&chID);
   
  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(chID==59) { histoCh1->Fill(energy); }//chiudo if
    else{ histoCh2->Fill(energy); }//chiudo else
  }//chiudo for
  
  std::cout << "Filled"<< std::endl;

  
}

void GetSpectrum(TTree* tree, TH1F* histoCh1, TH1F* histoCh2, Double_t MeanPedCh59, Double_t MeanPedCh315){
  
  Float_t energy;
  UShort_t chID;

  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("channelID",&chID);
   
  for(int i=0; i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(chID==59) { histoCh1->Fill(energy-MeanPedCh59); }//chiudo if
    else{ histoCh2->Fill(energy-MeanPedCh315); }//chiudo else
  }//chiudo for
  
  std::cout << "Filled"<< std::endl;
  
  
}

void FillETime(TTree* tree, TH2F* HistoCh59,TH2F* HistoCh315, Double_t MeanPedCh59, Double_t MeanPedCh315){

  Float_t energy;
  UShort_t chID;
  Long64_t time;

  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("channelID",&chID);
  tree->SetBranchAddress("unixTime",&time);

  for(int i=0; i<tree->GetEntries();i++){
    tree->GetEntry(i);
    if(chID==59){HistoCh59->Fill(time-1552505309,energy-MeanPedCh59);}
    else if(chID==315){HistoCh315->Fill(time-1552505309,energy-MeanPedCh315);}
    if(i%50000==0) std::cout<< time-1552505322 << std::endl;
  }//chiudo for

  std::cout<<"Filled EnergyTime" << std::endl;

  
  
}

TF1* FitNaSpectrum(TH1D* Profile){

  TF1* spectrum = new TF1("SpectrumFit","[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) + [3] / (exp( (x*[4]-(2*[1]*[1]/([1]+2*[1])))) + 1)+ [5] * exp(-( x-[6] )*( x-[6] )/( 2* [7]* [7])) + [8] / (exp( (x*[9]-(2*[6]*[6]/([1]+2*[6]))))+ 1)",30,90);
  
  spectrum->SetParameter(0,6000);
  spectrum->SetParameter(1,40);
  spectrum->SetParameter(2,3);
  spectrum->SetParameter(3,1700);
  spectrum->SetParameter(4,0.8);
  spectrum->SetParameter(5,350);
  spectrum->SetParameter(6,78);
  spectrum->SetParameter(7,3.2);
  spectrum->SetParameter(8,550);
  spectrum->SetParameter(9,0.8);

  spectrum->SetParLimits(5,250,500);
  spectrum->SetParLimits(7,3,4);
  spectrum->SetParLimits(8,300,700);
  spectrum->SetParLimits(9,0.8,1);
  



  Profile->Fit("SpectrumFit","R0");

  return spectrum;
}

void GetProfiles(TH2F* HistoCh59, TH2F* HistoCh315,Double_t* X,std::vector<std::vector<Double_t> >& Y){

  //gStyle->SetOptFit(0111);
   
  TH1D* HistoTemp59;
  TH1D* HistoTemp315;

  TCanvas* canvino;

  TF1* TempFitCh59;
  TF1* TempFitCh315;
  
  std::cout << "Getting Profiles" << std::endl;
  
  Double_t XMax=HistoCh59->GetXaxis()->GetXmax();
  Int_t Nbins=HistoCh59->GetXaxis()->GetNbins();

  for(int k=0;k < Nbins;k++){
    
    canvino= new TCanvas(("Projection"+to_string(k)).c_str(),("Projection"+to_string(k)).c_str(),1200,700);
    canvino->Divide(2,1);
    
    HistoTemp59=HistoCh59->ProjectionY(("ProjectionCh59"+to_string(k)).c_str(),k,k);
    HistoTemp315=HistoCh315->ProjectionY(("ProjectionCh315"+to_string(k)).c_str(),k,k); 

    TempFitCh59=FitNaSpectrum(HistoTemp59);
    TempFitCh315=FitNaSpectrum(HistoTemp315);
    
    HistoTemp59->SetTitle(("ProjectionCh59N"+to_string(k)).c_str());
    HistoTemp315->SetTitle(("ProjectionCh315N"+to_string(k)).c_str());
    
    canvino->cd(1);
    HistoTemp59->Draw();
    TempFitCh59->Draw("SAME");
    
    canvino->cd(2);
    HistoTemp315->Draw();
    TempFitCh315->Draw("SAME");
    
    canvino->SaveAs(("Plot/EnergyTime/Projections/Projection"+to_string(k)+".png").c_str());


    X[k]=XMax/Nbins*(k+1)-(XMax/Nbins)/2;
    
    Y[0][k]=TempFitCh59->GetParameter(1);
    Y[1][k]=TempFitCh59->GetParError(1);
    Y[2][k]=TempFitCh59->GetParameter(6);
    Y[3][k]=TempFitCh59->GetParError(6);

    Y[4][k]=TempFitCh315->GetParameter(1);
    Y[5][k]=TempFitCh315->GetParError(1);
    Y[6][k]=TempFitCh315->GetParameter(6);
    Y[7][k]=TempFitCh315->GetParError(6);
    
    
  }
  
}


void FillETemperature(TTree* tree, TH2D* HistoCh59,TH2D* HistoCh315, Double_t MeanPedCh59, Double_t MeanPedCh315){

  Float_t energy;
  UShort_t chID;
  Double_t temp2;
  Double_t temp3;

  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("channelID",&chID);
  tree->SetBranchAddress("temp2",&temp2);
  tree->SetBranchAddress("temp3",&temp3);

  for(int i=0; i<tree->GetEntries();i++){
    tree->GetEntry(i);
    if(chID==59){HistoCh59->Fill(temp2,energy-MeanPedCh59);}
    else if(chID==315){
      HistoCh315->Fill(temp3,energy-MeanPedCh315);
      std::cout << temp3 << "   " << energy-MeanPedCh315 << std::endl;
    }
    
  }//chiudo for
  
  std::cout<<"Filled EnergyTemperature" << std::endl;
  
  
  
}
