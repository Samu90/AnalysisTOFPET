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


void GetMeanTemperature(TTree* tree, Double_t* MeanTempCh59, Double_t* SigmaTempCh59,Double_t* MeanTempCh315, Double_t* SigmaTempCh315, Double_t* MeanTGlobal,Double_t* SigmaTGlobal,int j){

  UShort_t chID;
  Double_t temp1,temp2,temp3;

  TH1D* histo1 = new TH1D("pippoch1","pippoch1",20,10,35);
  TH1D* histoCh59 = new TH1D("pippoch59","pippoch59",20,20,40);
  TH1D* histoCh315 = new TH1D("pippoch315","pippoch315",20,20,40);

  tree->SetBranchAddress("temp1",&temp1);
  tree->SetBranchAddress("temp2",&temp2);
  tree->SetBranchAddress("temp3",&temp3);
  tree->SetBranchAddress("channelID",&chID);
   
  for(int i=0; i<tree->GetEntries(); i++){

    tree->GetEntry(i);
    histo1->Fill(temp1);
    if(chID==59) { histoCh59->Fill(temp2); }//chiudo if
    else{ histoCh315->Fill(temp3); }//chiudo else

  }//chiudo for
  
  MeanTGlobal[j]= histo1->GetMean();
  SigmaTGlobal[j]= histo1->GetMeanError();
  MeanTempCh59[j]= histoCh59->GetMean();
  SigmaTempCh59[j]= histoCh59->GetMeanError();
  SigmaTempCh315[j]= histoCh315->GetMeanError();
  MeanTempCh315[j]= histoCh315->GetMean();

  std::cout << "TMeanDone"<< std::endl;

}



TF1* FitNaSpectrumCB(TH1D* Profile){

std::cout << "kk1" << std::endl;
  Double_t min;
std::cout << Profile << std::endl;
std::cout << Profile->GetNbinsX() << std::endl;
std::cout << Profile->GetXaxis() << std::endl;
std::cout << Profile->GetMinimum() << std::endl;
std::cout << Profile->GetName() << std::endl;
  Profile->GetXaxis()->SetRangeUser(20,40);
std::cout << "kk2" << std::endl;
  min = Profile->GetBinCenter(Profile->GetMinimumBin());
std::cout << "kk3" << std::endl;
  Profile->GetXaxis()->UnZoom();

std::cout << "kk4" << std::endl;
  std::cout<<"minimum: "<< min<< std::endl;
  
  
  TF1* spectrum = new TF1(Form("SpectrumFit_%s", Profile->GetName()),"[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) + [3] / (exp( (x*[4]-(2*[1]*[1]/([1]+2*[1])))) + 1)+ [5] * exp(-( x-[6] )*( x-[6] )/( 2* [7]* [7])) +crystalball([8],[9],[10],[11],[12])",32,92);

  Double_t max;
  max= Profile->GetMaximum();
  
std::cout << "kk5" << std::endl;
  spectrum->SetParameter(0,max);
  spectrum->SetParameter(1,41);
  spectrum->SetParameter(2,3);
  spectrum->SetParameter(3,max/5);
  spectrum->SetParameter(4,0.82);
  spectrum->SetParameter(5,max/20);
  spectrum->SetParameter(6,80);
  spectrum->SetParameter(7,3.2);
  spectrum->SetParameter(8,max/12);
  spectrum->SetParameter(9,68);
  spectrum->SetParameter(10,4.3);
  spectrum->SetParameter(11,0.04);
  spectrum->SetParameter(12,-1.4);

std::cout << "kk6" << std::endl;
  spectrum->SetParLimits(10,3.8,6);
  spectrum->SetParLimits(3,4,2700);
  
std::cout << "kk7" << std::endl;
  Profile->Fit("SpectrumFit","R0");

std::cout << "kk8" << std::endl;
  return spectrum;
}


void RatioWithError(Double_t* A,Double_t* B,Double_t* sA,Double_t* sB,Double_t* ratio, Double_t* sigmaRatio,int Ndata){

  for(int i=0;i<Ndata;i++){
    ratio[i]=A[i]/B[i];
    sigmaRatio[i]=TMath::Sqrt((sA[i]*sA[i])/(B[i]*B[i])+(A[i]*A[i])/(B[i]*B[i]*B[i]*B[i])*(sB[i]*sB[i]));
  }
  
}

void SetStyleRatioPlot(TGraphErrors* ratioPlot,Double_t minRange,Double_t maxRange){

  ratioPlot->SetTitle("");
  ratioPlot->GetXaxis()->SetTitle("TMeanBoxRun [°C]");
  ratioPlot->GetYaxis()->SetTitle("PeakRatio");
  ratioPlot->GetXaxis()->SetTitleSize(0.08);
  ratioPlot->GetYaxis()->SetTitleSize(0.08);
  ratioPlot->GetYaxis()->SetTitleOffset(0.5);
  ratioPlot->GetXaxis()->SetLabelSize(0.075);
  ratioPlot->GetYaxis()->SetLabelSize(0.075);
  ratioPlot->GetYaxis()->SetLimits(minRange,maxRange);
  ratioPlot->GetYaxis()->SetRangeUser(minRange,maxRange);
  
}



