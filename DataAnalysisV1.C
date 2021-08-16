/*
 To run this macro:

 .L DataAnalysisV1.C
 DataAnalysisV1(ch,run,fit,filename)

 where "ch" is the channel number
  "run" is the run number
  "filename" is a string containing the file name with processed data.

 for example:
 DataAnalysisV1(19,448,1,"SignalIntegralV5_2014_FitRange_270.root")

 you can analyze another channel and/or another run
 in the same ROOT session.
 */

R__LOAD_LIBRARY(functions_C)
R__LOAD_LIBRARY(DRCEventV5_C)

void DataAnalysisV1(Int_t ch, Int_t run, int fit, TString SignalFile)
{

  gStyle->SetOptFit(1111);

  if( gROOT->FindObject(SignalFile.Data()) )
    ((TFile *)gROOT->FindObject(SignalFile.Data()))->Close(); //file location


  TFile *file = TFile::Open(SignalFile.Data());


  file->cd(Form("%d",run)); //go into folder

  TTree *treeParam = (TTree *) gROOT->FindObject("ADRIANOParam"); //find first TTree

  Float_t BiasVoltage, Energy, cnt, PartType, Angle, x, z, pressure, NEvents, NTrigMB1;

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

  treeParam->GetEvent(0);

  //     if(NEvents != NTrigMB2) cout << Form("***Warning!!! NTrigMB1=%d while NTrigMB2=%d",(Int_t)NEvents, (Int_t)NTrigMB2);

  const char *PartName[4]={"pion", "proton", "muon","n.a."};
  if(PartType<1 || PartType>3) PartType = 4;


  cout << "Run | beam |Energy|   X  |   Z  |Angle| NEvents| psia| BiasV| counts\n";
  cout << Form("%3d |%6s|%sGeV| %4d | %4d |%4d | %6d | %4.1f| %2.1fV| %4dK  \n", run, PartName[(Int_t)PartType-1], (Energy>0.9)?Form("%3d",(Int_t)Energy):Form("%1.1f",Energy), (Int_t)x, (Int_t)z, (Int_t)Angle, (Int_t)NEvents, pressure, BiasVoltage, (Int_t)cnt);


  TTree *ADRIANOAmp1 = (TTree *) gROOT->FindObject("ADRIANOAmp1"); //finds amplitude TTree
    ADRIANOAmp1->AddFriend("ADRIANOPed1");
    ADRIANOAmp1->AddFriend("ADRIANOPedSigma1");
    ADRIANOAmp1->AddFriend("ADRIANOTimeAmp1");

  //ADRIANOAmp1->AddFriend("ADRIANOFitCheck"); //adds branches from different trees to relevant tree
  //ADRIANOAmp1->AddFriend("ADRIANOIntegral1");


  if( gROOT->FindObject("c1") )
    ((TCanvas *)gROOT->FindObject("c1"))->Close();


  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);
  c1->cd(1);

  //smoothing + fitting
  if (ch > 15){ //plastic
    ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(80,0.,10.)", ch, ch), Form("ADRIANOTimeAmp1.ch%d > 21 && ADRIANOTimeAmp1.ch%d < 27 && ADRIANOPedSigma1.ch%d < 3", ch, ch, ch));
  }
  else { //glass
      ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(600,0.,0.)", ch, ch), Form("ADRIANOAmp1.ch%d < 400", ch));
  }

  TH1F *htmp = ((TH1F *)gROOT->FindObject(Form("htmp%d", ch)));
    if (ch > 15){
        htmp->Smooth(3);
    }

    Double_t xmin = htmp->GetBinLowEdge(1);
    Double_t xBinMax = htmp->FindLastBinAbove(1.);
    Double_t xmax = htmp->GetBinLowEdge(xBinMax+1);
    //htmp->Reset();
    htmp->Draw();
    htmp->SetTitle(Form("Amplitude ch %d run %d", ch, run));

    if(fit){
      xBinMax = htmp->GetMaximumBin(); //bin of maximum value
      xmax = htmp->GetBinLowEdge(xBinMax); //x value of maximum bin
      Double_t nbins = htmp->GetNbinsX(); //# of bins
      //Double_t yExpo = htmp->GetBinContent(5); //?
      Double_t Area = htmp->Integral();//*nbins/10.;
      Double_t FitRangeA[2]   = {htmp->GetBinLowEdge(1), htmp->GetBinLowEdge(nbins+1)}; //x-axis fit range
      Double_t StartValueA[4]    = {10.,      xmax,       Area,   2.};
      Double_t LowLimitValueA[4] = { 0.5,  xmax-10.,   Area/10.,   0.};
      Double_t HiLimitValueA[4]  = {1.e2, xmax+10.,  Area*1.e3, 100};
      Double_t fitparams[4];
      Double_t fiterrors[4];
      Double_t chi;
      Int_t ndf;

    //   for(Int_t idx=0; idx<6; idx++)
    //     cout << Form("%e\t%e\t%e\n", LowLimitValueA[idx], StartValueA[idx], HiLimitValueA[idx]);

    TF1 *fitsnrA = langaufit(htmp, FitRangeA, StartValueA, LowLimitValueA, HiLimitValueA, fitparams, fiterrors, &chi, &ndf);
    //, 3, "QRBWI", "RBWI");
    cout << endl;
    cout << Form("%d, %f, %f, %f, %f, %f, %f, %f, %d", (Int_t)htmp->GetEntries(), htmp->GetMean(), htmp->GetStdDev(), fitparams[0], fiterrors[0], fitparams[1], fiterrors[1], chi, ndf) << endl;
}

  //==============================
  c1->cd(2);
    ADRIANOAmp1->Draw(Form("ADRIANOTimeAmp1.ch%d>>htmp1(300,0.,0.)", ch));
    TH1F *htmp1 = ((TH1F *)gROOT->FindObject("htmp1"));
    htmp1->SetTitle(Form("CH %d Time Amp", ch));

  //==============================
  c1->cd(3);
    ADRIANOAmp1->Draw(Form("ADRIANOPed1.ch%d>>htmp2(300,0.,0.)", ch));
    TH1F *htmp2 = ((TH1F *)gROOT->FindObject("htmp2"));
    htmp2->SetTitle(Form("CH %d Ped", ch));


  //==============================
  c1->cd(4);
    ADRIANOAmp1->Draw(Form("ADRIANOPedSigma1.ch%d>>htmp3(300,0.,0.)", ch));//, Form("ADRIANOAmp1.ch%d < 1000", ch[i]));
    TH1F *htmp3 = ((TH1F *)gROOT->FindObject("htmp3"));
    htmp3->SetTitle(Form("CH %d Ped Sigma", ch));

    c1->Draw();
    //c1->SaveAs(Form("./RUN_%d_ch%d_ogfitfn2.png",run, ch));
}
