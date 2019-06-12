
TF1* FitNaSpectrumCBBar(TH1D* Profile,Int_t* fitStatus ,Int_t chID){

  
  Double_t min;
  Double_t peak1,peak2;
  Double_t RMS1,RMS2;
  Int_t EndPlot;
  
  peak1 = Profile->GetBinCenter(Profile->GetMaximumBin());
  
  Profile->GetXaxis()->SetRangeUser(peak1/2,peak1);
  min = Profile->GetBinCenter(Profile->GetMinimumBin());
  Profile->GetXaxis()->UnZoom();

  Profile->GetXaxis()->SetRangeUser( peak1 - (peak1 - min) , peak1+ (peak1 - min) );
  RMS1 = Profile->GetRMS();
  Profile->GetXaxis()->UnZoom();

  Int_t Last= Profile->GetXaxis()->GetLast();
  std::cout << Last << std::endl;

  for(int i=Last;i>0;i--){
    if( Profile->GetBinContent(i) > 100) {
      EndPlot=i;
      break;
    }
  }
  
  std::cout << EndPlot << std::endl;

  peak2=EndPlot/1.09;

  std::cout << peak2 << std::endl;
  
  Profile->GetXaxis()->SetRangeUser( peak2 - ( EndPlot - peak2 ) , EndPlot );
  peak2 = Profile->GetBinCenter(Profile->GetMaximumBin());
  RMS2 = Profile->GetRMS();
  Profile->GetXaxis()->UnZoom();

  std::cout << "p1 , RMS1 , p2 , RMS2 ->  " << peak1<< " " << RMS1 << " " << peak2 << " " <<RMS2<< std::endl;

  TF1* spectrum = new TF1(Form("SpectrumFit_%s", Profile->GetName()),"[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) + [3] / (exp( (x*[4]-(2*[1]*[1]/([1]+2*[1])))) + 1)+ [5] * exp(-( x-[6] )*( x-[6] )/( 2* [7]* [7])) +crystalball([8],[9],[10],[11],[12])",min+1,EndPlot);
  
  //TF1* spectrumBar = new TF1(Form("SpectrumBarFit_%s", Profile->GetName()),"[0] * exp(-( x-[1] )*( x-[1] )/( 2* [2]* [2])) +[3] * exp(-( x-[6] )*( x-[6] )/( 2* [5]* [5])) +[4] / (exp( (x*[7]-(2*[6]*[6]/([6]+2*[1])))) + 1)");
  

  Double_t max;
  max= Profile->GetMaximum();
  

  if(chID==59 ){
    
    spectrum->SetParameter(0,max);
    spectrum->SetParameter(1,peak1);
    spectrum->SetParameter(2,3);
    spectrum->SetParameter(3,max/5);
    spectrum->SetParameter(4,0.82);
    spectrum->SetParameter(5,max/20);
    spectrum->SetParameter(6,peak2);
    //spectrum->SetParameter(7,3.2);
    spectrum->SetParameter(7,RMS2);
    spectrum->SetParameter(8,max/12);
    spectrum->SetParameter(9,peak2/1.21);
    spectrum->SetParameter(10,4.3);
    spectrum->SetParameter(11,0.04);
    spectrum->SetParameter(12,-1.4);
    
    spectrum->SetParLimits(10,3.8,4.6);
    spectrum->SetParLimits(3,4,2700);
    spectrum->SetParLimits(11,0.02,1);

    *fitStatus=Profile->Fit(Form("SpectrumFit_%s", Profile->GetName()),"RM0");
    
    return spectrum;

  } else {
      
    spectrum->SetParameter(0,max);
    spectrum->SetParameter(1,peak1);
    spectrum->SetParameter(2,3);
    spectrum->SetParameter(3,max/5);
    spectrum->SetParameter(4,0.82);
    spectrum->SetParameter(5,max/20);
    spectrum->SetParameter(6,peak2);
    //spectrum->SetParameter(7,3.2);
    spectrum->SetParameter(7,RMS2);
    spectrum->SetParameter(8,max/12);
    spectrum->SetParameter(9,peak2/1.21);
    spectrum->SetParameter(10,3.8);
    spectrum->SetParameter(11,0.02);
    spectrum->SetParameter(12,-4e5);
    
    //spectrum->SetParLimits(10,3.6,4.6);
    spectrum->SetParLimits(3,4,2700);
    spectrum->SetParLimits(11,0.02,1);
    spectrum->SetParLimits(12,-5e-7,-2000);

    
    *fitStatus=Profile->Fit(Form("SpectrumFit_%s", Profile->GetName()),"R0M");
      
    return spectrum;
  }  


}



