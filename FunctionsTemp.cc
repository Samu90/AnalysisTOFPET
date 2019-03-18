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
