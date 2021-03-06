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
    else{ histoCh2->Fill(energy); }//chiudo else
  }//chiudo for
  
  std::cout << "Filled"<< std::endl;

  *pedch1 = histoCh1->GetMean();
  *pedch2 = histoCh2->GetMean();
  
  delete histoCh1;
  delete histoCh2;
}


void GetSpectrum(TTree* tree, TH1D* histoCh1, TH1D* histoCh2, Double_t MeanPedCh59, Double_t MeanPedCh315){
  
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



void GetMeanTemperature(TTree* tree, Double_t* MeanTempCh1, Double_t* SigmaTempCh1,Double_t* MeanTempCh2, Double_t* SigmaTempCh2){

  UShort_t chID;
  Double_t temp2,temp3;

  TH1D* histoCh1 = new TH1D("pippoch1","pippoch1",20,20,40);
  TH1D* histoCh2 = new TH1D("pippoch2","pippoch2",20,20,40);

  tree->SetBranchAddress("temp2",&temp2);
  tree->SetBranchAddress("temp3",&temp3);
  tree->SetBranchAddress("channelID",&chID);
   
  for(int i=0; i<tree->GetEntries(); i++){

    tree->GetEntry(i);
    if(chID==59) { histoCh1->Fill(temp2); }//chiudo if
    else{ histoCh2->Fill(temp3); }//chiudo else

  }//chiudo for
  
  *MeanTempCh1= histoCh1->GetMean();
  *SigmaTempCh1=histoCh1->GetMeanError();
  *MeanTempCh2= histoCh2->GetMean();
  *SigmaTempCh2=histoCh2->GetMeanError();
  
  std::cout << "TMeanDone"<< std::endl;

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

TF1* FitNaSpectrumCB(TH1D* Profile){

  
  TF1* spectrum = new TF1("SpectrumFit","[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) + [3] / (exp( (x*[4]-(2*[1]*[1]/([1]+2*[1])))) + 1)+ [5] * exp(-( x-[6] )*( x-[6] )/( 2* [7]* [7])) +crystalball([8],[9],[10],[11],[12])",30,88);
  
  spectrum->SetParameter(0,6000);
  spectrum->SetParameter(1,41);
  spectrum->SetParameter(2,3);
  spectrum->SetParameter(3,1700);
  spectrum->SetParameter(4,0.82);
  spectrum->SetParameter(5,350);
  spectrum->SetParameter(6,80);
  spectrum->SetParameter(7,3.2);
  spectrum->SetParameter(8,572);
  spectrum->SetParameter(9,68);
  spectrum->SetParameter(10,4.3);
  spectrum->SetParameter(11,0.04);
  spectrum->SetParameter(12,-1.4);

  Profile->Fit("SpectrumFit","R0");

  return spectrum;
}
