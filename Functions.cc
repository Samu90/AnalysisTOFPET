
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
    if(chID==59){HistoCh59->Fill(time-540000102715,energy-MeanPedCh59);}
    else if(chID==315){HistoCh315->Fill(time-540000102715,energy-MeanPedCh315);}
    
  }//chiudo for

  std::cout<<"Filled EnergyTime" << std::endl;

}

