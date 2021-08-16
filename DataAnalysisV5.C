/*
 To run this macro:

 .L DataAnalysisV5.C
 DataAnalysisV5(ch,run,filename)

 where "ch" is the channel number
  "run" is the run number
  "filename" is a string containing the file name with processed data.

 for example:

 DataAnalysisV5(19,448,"SignalIntegralV5_2014_FitRange_270.root")

 you can analyze another channel and/or another run
 in the same ROOT session.

 DataAnalysisV5()
 will use as default arguments
 ch=16
 run=420
 filename="SignalIntegralV5_2014_FitRange_270.root"

 */

Double_t xMaxValI = 0.;
Double_t xMaxValA = 0.;

void DataAnalysisV5(Int_t ch=16, Int_t run=420, TString SignalFile="SignalIntegralV5_2014_FitRange_270.root"){

//   gROOT->Reset();

  gSystem->Load("functions_C.so");

  gStyle->SetOptFit(1111);

  if( gROOT->FindObject(SignalFile.Data()) )
    ((TFile *)gROOT->FindObject(SignalFile.Data()))->Close(); //file location


  TFile *file = TFile::Open(SignalFile.Data());

  file->cd(Form("%d",run)); //go into folder

  TTree *treeParam = (TTree *) gROOT->FindObject("ADRIANOParam"); //find first TTree

  Float_t BiasVoltage, Energy, cnt, PartType, Angle, x, z, pressure, NEvents, NTrigMB1, Gain1, Gain2;

  treeParam->SetBranchAddress("BiasVoltage", &BiasVoltage); //finds branches in TTree
  treeParam->SetBranchAddress("Energy", &Energy);
  treeParam->SetBranchAddress("cnt", &cnt);
  treeParam->SetBranchAddress("PartType", &PartType);
  treeParam->SetBranchAddress("Angle", &Angle);
  treeParam->SetBranchAddress("x", &x);
  treeParam->SetBranchAddress("z", &z);
  treeParam->SetBranchAddress("pressure", &pressure);
  treeParam->SetBranchAddress("NTrigMB1", &NTrigMB1);
  treeParam->SetBranchAddress("NTrigMB2", &NEvents);
  treeParam->SetBranchAddress("Gain1", &Gain1);
  treeParam->SetBranchAddress("Gain2", &Gain2);

  treeParam->GetEvent(0);

  //     if(NEvents != NTrigMB2) cout << Form("***Warning!!! NTrigMB1=%d while NTrigMB2=%d",(Int_t)NEvents, (Int_t)NTrigMB2);

  Char_t *PartName[4]={"pion", "proton", "muon","n.a."};
  if(PartType<1 || PartType>3) PartType = 4;


  cout << "Run | beam |Energy|   X  |   Z  |Angle| NEvents| psia| BiasV| counts| gain1| gain2\n";
  cout << Form("%3d |%6s|%sGeV| %4d | %4d |%4d | %6d | %4.1f| %2.1fV| %4dK| %3d| %3d \n", run, PartName[(Int_t)PartType-1], (Energy>0.9)?Form("%3d",(Int_t)Energy):Form("%1.1f",Energy), (Int_t)x, (Int_t)z, (Int_t)Angle, (Int_t)NEvents, pressure, BiasVoltage, (Int_t)cnt, (Int_t)Gain1/1e5, (Int_t)Gain2/1e5);


  TTree *ADRIANOAmp1 = (TTree *) gROOT->FindObject("ADRIANOAmp1"); //find second TTree

  ADRIANOAmp1->AddFriend("ADRIANOFitCheck"); //adds different branches to tree
  ADRIANOAmp1->AddFriend("ADRIANOIntegral1");
  ADRIANOAmp1->AddFriend("ADRIANOAmp2");
  ADRIANOAmp1->AddFriend("ADRIANOIntegral2");

  if( gROOT->FindObject("c1") )
    ((TCanvas *)gROOT->FindObject("c1"))->Close();

  TCanvas *c1 = new TCanvas("c1","c1",1400,500);
  c1->Divide(2);

  c1->cd(1);
  ADRIANOAmp1->Draw(Form("ADRIANOIntegral1.ch%d>>htmp(1000,0.,0.)", ch, ch), Form("(ADRIANOFitCheck.ch%d&1) && (ADRIANOFitCheck.ch%d&2)", ch, ch));

  Double_t xmin = htmp->GetBinLowEdge(1);
  Int_t xBinMax = htmp->FindLastBinAbove(10.);
  Double_t xmax = htmp->GetBinLowEdge(xBinMax+1);
  htmp->Reset();

  ADRIANOAmp1->Draw(Form("ADRIANOIntegral1.ch%d>>hInteg(150,%f,%f)", ch, xmin, xmax), Form("(ADRIANOFitCheck.ch%d&1) && (ADRIANOFitCheck.ch%d&2)", ch, ch));

  hInteg->SetTitle(Form("Integral ch %d run %d", ch, run));
  xBinMax = hInteg->GetMaximumBin();
  xmax = hInteg->GetBinLowEdge(xBinMax);
  Int_t nbins = hInteg->GetNbinsX();
  Double_t yExpo = hInteg->GetBinContent(3);
  Double_t Area = hInteg->Integral()*nbins;
  Double_t FitRangeI[2]   = {hInteg->GetBinLowEdge(5), hInteg->GetBinLowEdge(xBinMax*2)};
  Double_t StartValueI[6]    = {9.e2,     xmax,       Area, 1.e1, yExpo, -1.e-5};
  Double_t LowLimitValueI[6] = {5.e2, xmax/10.,   Area/10.,   0., TMath::Min(1., yExpo/10.), -1.e-1};
  Double_t HiLimitValueI[6]  = {9.e4, xmax*10.,  Area*1.e3, 1.e5, TMath::Max(1.e4, yExpo*10.), 0.};

//   for(Int_t idx=0; idx<6; idx++)
//     cout << Form("%e\t%e\t%e\n", LowLimitValueI[idx], StartValueI[idx], HiLimitValueI[idx]);


  TF1 *fitsnrI = langaufit3(hInteg, FitRangeI, StartValueI, LowLimitValueI, HiLimitValueI, 3, "QRBWI", "RBWI");

  xMaxValI = fitsnrI->GetMaximumX();

  cout << Form("\n*** Integral: The maximum is at %e\n\n", xMaxValI);

  //=================================================

  c1->cd(2);
  ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp(1000,0.,0.)", ch, ch), Form("(ADRIANOFitCheck.ch%d&1) && (ADRIANOFitCheck.ch%d&2)", ch, ch));

  xmin = htmp->GetBinLowEdge(1);
  xBinMax = htmp->FindLastBinAbove(10.);
  xmax = htmp->GetBinLowEdge(xBinMax+1);
  htmp->Reset();

  ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>hAmp(80,%f,%f)", ch, xmin, xmax), Form("(ADRIANOFitCheck.ch%d&1) && (ADRIANOFitCheck.ch%d&2)", ch, ch));

  hAmp->SetTitle(Form("Amplitude ch %d run %d", ch, run));
  xBinMax = hAmp->GetMaximumBin();
  xmax = hAmp->GetBinLowEdge(xBinMax);
  nbins = hAmp->GetNbinsX();
  yExpo = hAmp->GetBinContent(5);
  Area = hAmp->Integral()*nbins/10.;
  Double_t FitRangeA[2]   = {hAmp->GetBinLowEdge(6), hAmp->GetBinLowEdge(xBinMax*2)};
  Double_t StartValueA[6]    = {10.,      xmax,       Area,   2., yExpo, -1.e-3};
  Double_t LowLimitValueA[6] = { 1.,  xmax/10.,   Area/10.,   0., TMath::Min(1., yExpo/10.), -1.e-1};
  Double_t HiLimitValueA[6]  = {1.e2, xmax*10.,  Area*1.e3, 1.e2, TMath::Max(1.e4, yExpo*10.), 0.};

//   for(Int_t idx=0; idx<6; idx++)
//     cout << Form("%e\t%e\t%e\n", LowLimitValueA[idx], StartValueA[idx], HiLimitValueA[idx]);


  TF1 *fitsnrA = langaufit3(hAmp, FitRangeA, StartValueA, LowLimitValueA, HiLimitValueA, 3, "QRBWI", "RBWI");

  xMaxValA = fitsnrA->GetMaximumX();

  Double_t xMaxValRatio = xMaxValI/xMaxValA;

  cout << Form("\n*** Amplitude: The maximum is at %e\n", xMaxValA);
  cout << Form("*** Integral:  The maximum is at %e\n", xMaxValI);
  cout << Form("*** Maximum ratio: %e\n\n", xMaxValRatio);


}
