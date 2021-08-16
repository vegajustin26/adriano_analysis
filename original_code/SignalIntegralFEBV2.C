/****
 * ROOT macro to process raw waveforms
 *
 * To run this macro it is needed to compile DRCEventV5.C and functions.C first
 * At ROOT prompt run
 * root [0] .L DRCEventV5.C++
 * root [1] .L functions.C++
 * root [2] .L SignalIntegralFEBV2.C++
 *
 * If those macros have already been compiled, it is enough to load its library
 *
 * root [0] gSystem->Load("SignalIntegralFEBV2_C.so");
 *
 * This macro is usually executed through AnalysisIntegralFEBV2.C macro that helps to setup all arguments.
 * As it is now this macro supports up to 2 FEB boards.
 * Usually channels 16 and 17 have special meaning, but his can change for different TestBeams.
 * Some TB setups uses SiPM and/or PMT that usually have opposite polarity,
 * also MWPC data can be present or not.
 *
 ****/

#ifndef ROOT_SignalIntegralFEBV2
#define ROOT_SignalIntegralFEBV2

#include <Riostream.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TF1.h>
#include <TSystemDirectory.h>
#include <TClassTable.h>
#include <TString.h>
#include <TObjString.h>
#include <TPaveStats.h>
#include <TMinuit.h>
#include <TROOT.h>
#include <TStyle.h>

#include "DRCEventV5.h"
#include "functions.C"

UShort_t GetMWPCxy(Double_t &TB4TimeEpoc, TTree *trMWPC, Int_t &MWPCTrig, Short_t &MWPCxwire, Short_t &MWPCywire, Short_t &MWPCxtdc, Short_t &MWPCytdc, Double_t &dTime, Short_t TB42MPWCTimeOffset);

void SignalIntegralFEBV2(Int_t RunNo, const Int_t Nch, Int_t *chlist, Bool_t *IsSiPMChList, Int_t *XMaxList, Int_t FirstEvent=0, Int_t EventsToProcess=-1, Int_t dataset=2019, Int_t subdataset=1, TString path="../data/", Bool_t debug=0, Int_t wfsize = 127, Int_t nfit = 1, Bool_t bShift=1, Bool_t bFast=0, Bool_t bUseMWPC=0)
{

  gROOT->Reset();
  gStyle->SetOptFit(1111);
  new TMinuit;
  Int_t iconverge=0;
//   TCanvas c1("c1","c1");


  if( gClassTable->GetID("DRCEvent") < 0 )
  {
    cout << "Loading DRCEvent library\n";

    #if defined(__linux__)
    gSystem->Load("DRCEventV5_C.so");
    gSystem->Load("functions_C.so");
    #else
    gSystem->Load("DRCEventV5_C.dll");
    gSystem->Load("functions_C.dll");
    #endif


  }



  TString sAnalysisDir(gSystem->pwd());
  TFile *pfiled = new TFile(Form("SignalIntegralFEBV1_%03d_%02d.root",dataset,subdataset),"UPDATE", "ADRIANO", 9);

  TString sParam("Energy:x:z:Angle:BiasVoltage:cnt:PartType:pressure:NTrigMB1:NTrigMB2:DataWfSize:FitWfSize:Nch:RunBeginTime:RunEndTime:Gain1:Gain2");
  for(Int_t j=0; j<Nch; j++)
    sParam.Append(Form(":ch%d",j));

  TNtuple *ntADRIANOParam = new TNtuple("ADRIANOParam","ADRIANOParam",sParam.Data(),320000);

  TString sNtuples("");
  for(Int_t j=0; j<Nch; j++)
    sNtuples.Append(Form("ch%d:",chlist[j]));

  sNtuples.Append("evt");

  TNtuple *ntADRIANOAmp1 = new TNtuple("ADRIANOAmp1","ADRIANOAmp1",sNtuples.Data(),320000);
  TNtuple *ntADRIANOTimeAmp1 = new TNtuple("ADRIANOTimeAmp1","ADRIANOTimeAmp1",sNtuples.Data(),320000);
  TNtuple *ntADRIANOIntegral1 = new TNtuple("ADRIANOIntegral1","ADRIANOIntegral1",sNtuples.Data(),320000);
  TNtuple *ntADRIANOPed1 = new TNtuple("ADRIANOPed1","ADRIANOPed1",sNtuples.Data(),320000);
  TNtuple *ntADRIANOPedSigma1 = new TNtuple("ADRIANOPedSigma1","ADRIANOPedSigma1",sNtuples.Data(),320000);


  TNtuple *ntADRIANOAmp2 = new TNtuple("ADRIANOAmp2","ADRIANOAmp2",sNtuples.Data(),320000);
  TNtuple *ntADRIANOTimeAmp2 = new TNtuple("ADRIANOTimeAmp2","ADRIANOTimeAmp2",sNtuples.Data(),320000);
  TNtuple *ntADRIANOIntegral2 = new TNtuple("ADRIANOIntegral2","ADRIANOIntegral2",sNtuples.Data(),320000);
  TNtuple *ntADRIANOPed2Left = new TNtuple("ADRIANOPed2Left","ADRIANOPed2Left",sNtuples.Data(),320000);
  TNtuple *ntADRIANOPed2Right = new TNtuple("ADRIANOPed2Right","ADRIANOPed2Right",sNtuples.Data(),320000);

  TNtuple *ntADRIANOFitCheck = new TNtuple("ADRIANOFitCheck","ADRIANOFitCheck",sNtuples.Data(),320000);

  TNtuple *ntADRIANOFitParams16 = new TNtuple("ADRIANOFitParams16","ADRIANOFitParams16","param0:param1:param2:param3:param4:param5:param6:param7:chi2:edm:evt",320000);
  TNtuple *ntADRIANOFitParams17 = new TNtuple("ADRIANOFitParams17","ADRIANOFitParams17","param0:param1:param2:param3:param4:param5:param6:param7:chi2:edm:evt",320000);
  TNtuple *ntADRIANOFitParams18 = new TNtuple("ADRIANOFitParams18","ADRIANOFitParams18","param0:param1:param2:param3:param4:param5:param6:param7:chi2:edm:evt",320000);
  TNtuple *ntADRIANOFitParams19 = new TNtuple("ADRIANOFitParams19","ADRIANOFitParams19","param0:param1:param2:param3:param4:param5:param6:param7:chi2:edm:evt",320000);

  TNtuple *ntMWPC = 0x0;
  if(bUseMWPC){
    cout << "MWPC data will be used\n";

    ntMWPC = new TNtuple("MWPC","MWPC","timeEPOC:xwire:ywire:xtdc:ytdc:dTime:MWPCStatus:evt",320000);
  }

  TH1F *htmp = new TH1F("htmp","htmp",wfsize,0.,wfsize);


  const char *ext=".root";



  TTree *trMWPC = 0x0;

  if(bUseMWPC)
  {

    //     const char *ext=".root";
    TString fnameMWPC;
    string inputMWPC = "";

    TString pathMWPC="../data/mwpc/data/";
    TSystemDirectory dirMWPC("dataMWPC", pathMWPC.Data());

    TList*filesMWPC = dirMWPC.GetListOfFiles();


    if (filesMWPC)
    {
      TSystemFile *fileMWPC;
      TIter nextMWPC(filesMWPC);
      while ((fileMWPC=(TSystemFile*)nextMWPC())) {
	fnameMWPC = fileMWPC->GetName();
	if (!fileMWPC->IsDirectory() && fnameMWPC.EndsWith(ext)) {
	  if(fnameMWPC.BeginsWith(Form("RUN_%d_",RunNo))){
	    cout << fnameMWPC.Data() << endl;

	    inputMWPC = fnameMWPC.Data();

	    fnameMWPC.Prepend(pathMWPC);
	    string RootFileNameMWPC = string(fnameMWPC);

	    TFile *fMWPC = new TFile(Form("%s",RootFileNameMWPC.data()),"READ");

	    trMWPC = (TTree *)fMWPC->Get("trMWPC");

	    break;
	  }
	}
      }
    }
  }







  TString fname;
  string input = "";

  TSystemDirectory dir("data", path.Data());

  TList*files = dir.GetListOfFiles();


  Float_t BiasVoltage=-999., Energy=-999., cnt=-999., PartType=-999., Angle=-999., x=-999., z=-999., pressure=-999., DataWfSize=-999., RunBeginTime=-999., RunEndTime=-999., Gain1=-1., Gain2=-1.;

  if (files)
  {
    TSystemFile *file;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext)) {
	if(fname.BeginsWith(Form("RUN_%d_",RunNo)) || fname.BeginsWith(Form("RUN_FEB_%d_",RunNo)) || fname.BeginsWith(Form("RUN_TB4_%d_",RunNo)) ){
	  cout << fname.Data() << endl;

	  input = fname.Data();


	  //================================================================
	  const string delimiters = "_";
	  // Skip delimiters at beginning.
	  string::size_type lastPos = input.find_first_not_of(delimiters, 0);
	  // Find first "non-delimiter".
	  string::size_type pos     = input.find_first_of(delimiters, lastPos);

	  while (string::npos != pos || string::npos != lastPos){

	    TString subinput = input.substr(lastPos, pos - lastPos);

	    if(subinput.Contains("BV")){
	      BiasVoltage = subinput.Atof();
	    }
	    else if(subinput.Contains("GeV")){
	      Energy = subinput.Atof();
	    }
	    else if(subinput.Contains("Kcnt")){
	      cnt = subinput.Atoi();
	    }
	    else if(subinput.Contains("deg")){
	      Angle = subinput.Atoi();
	    }
	    else if(subinput.Contains("x")){
	      x = TString(subinput.Strip(TString::kLeading,'x')).Atoi();
	    }
	    else if(subinput.Contains("z")){
	      z = TString(subinput.Strip(TString::kLeading,'z')).Atoi();
	    }
	    else if(subinput.Contains("psi")){
	      pressure = subinput.Atof();
	    }
	    else if(subinput.Contains("Gain")){
	      subinput.ReplaceAll("Gain","");

	      TObjArray *tokenArr = subinput.Tokenize("-");
	      Int_t nGains=tokenArr->GetEntries();

	      Gain1 = ((TObjString *)(tokenArr->At(0)))->String().Atof()*10000.;
	      Gain1 += ((TObjString *)(tokenArr->At(1)))->String().Atof();

	      if(nGains==4){
		Gain2 = ((TObjString *)(tokenArr->At(2)))->String().Atof()*10000.;
		Gain2 += ((TObjString *)(tokenArr->At(3)))->String().Atof();
	      }

              // cout << Gain1 << " " << Gain2 << endl;
              //
              // cout << Form("%8.0f  %8.0f  \n", Gain1, Gain2);
              //
              // cout << Form("%d  %d  \n", (Int_t)(Gain1/10000), ((Int_t)Gain1)%10000);
              // cout << Form("%d  %d  \n", (Int_t)(Gain2/10000), ((Int_t)Gain2)%10000);

	      if(tokenArr){
		tokenArr->Clear();
		tokenArr->Delete();
		delete tokenArr;
	      }

	    }

	    // Skip delimiters.  Note the "not_of"
	    lastPos = input.find_first_not_of(delimiters, pos);
	    // Find next "non-delimiter"
	    pos = input.find_first_of(delimiters, lastPos);

	  }

	  if(input.find("electron")!=string::npos){
	    PartType = 0;
	  }
	  else if(input.find("pion")!=string::npos){
	    PartType = 1;
	  }
	  else if(input.find("proton")!=string::npos){
	    PartType = 2;
	  }
	  else if(input.find("muon")!=string::npos){
	    PartType = 3;
	  }
	  else if(input.find("LED")!=string::npos){
	    PartType = 4;
	  }

	  break;
	}
      }
    }
  }

  if((input.empty()))
  {
    cout << Form("\n\n $$ No Run %d in your data $$\n", RunNo);
    pfiled->Close();
    return;
  }

  Bool_t bad=0;
  if (!(fname.EndsWith(ext)) )
    bad=1;

  if(bad){
    cout << Form("RunNo: %d don't exist in data\n", RunNo);
    pfiled->Close();
    return;
  }

  fname.Prepend(path);
  string RootFileName = string(fname);




  ntADRIANOAmp1->Reset();
  ntADRIANOTimeAmp1->Reset();
  ntADRIANOIntegral1->Reset();
  ntADRIANOPed1->Reset();
  ntADRIANOPedSigma1->Reset();
  ntADRIANOAmp2->Reset();
  ntADRIANOTimeAmp2->Reset();
  ntADRIANOIntegral2->Reset();
  ntADRIANOPed2Left->Reset();
  ntADRIANOPed2Right->Reset();
  ntADRIANOFitCheck->Reset();
  ntADRIANOParam->Reset();
  ntADRIANOFitParams16->Reset();
  ntADRIANOFitParams17->Reset();
  ntADRIANOFitParams18->Reset();
  ntADRIANOFitParams19->Reset();


  cout << Form("filename: %s RunNo: %d\n", RootFileName.data(), RunNo);
  TString sCurDir(gSystem->pwd());

  gSystem->cd(sAnalysisDir.Data());

  pfiled->rmdir(Form("%d",RunNo));
  pfiled->mkdir(Form("%d",RunNo));

  gSystem->cd(sCurDir.Data());



  cout << gSystem->pwd() << endl;
  TFile f(Form("%s",RootFileName.data()));

  if(f.IsZombie()){
    f.Close();
    cout << Form("filename: %s RunNo: %d is bad\n", RootFileName.data(), RunNo);
    bad=1;
  }

  if(bad){
    pfiled->rmdir(Form("%d",RunNo));
    pfiled->Close();
    return;
  }

  TTree *T1 = (TTree*)f.Get("DRCMB1");
  DRCEvent *event1 = 0;
  T1->SetBranchAddress("DRCEvent", &event1);

  TTree *T2 = (TTree*)f.Get("DRCMB2");
  DRCEvent *event2 = 0;
  T2->SetBranchAddress("DRCEvent", &event2);


  Int_t nEntriesMB1 = T1->GetEntries();
  Int_t nEntriesMB2 = T2->GetEntries();

  Int_t nEntries = TMath::Max(nEntriesMB1, nEntriesMB2);

  cout << Form("\n### Run: %d\n", RunNo);
  cout << Form("\nMB1 has %d events   ###   ", nEntriesMB1);
  cout << Form("MB2 has %d events\n", nEntriesMB2);

  cout << "\nchannels selected: ";
  for(Int_t j=0; j<Nch; j++)
    cout << Form("%d ",  chlist[j]);

  cout << "\n\nParameters:\nEnergy\tx\tz\tAngle\tBiasVoltage\tcounts\tPartType\tpressure\tGain1\tGain2\n";
  cout << Form("%.0f\t%.0f\t%.0f\t%.0f\t%.1f\t\t%.0f\t%.0f\t\t%.1f\t\t%.0f\t\t%.0f", Energy, x, z, Angle, BiasVoltage, cnt, PartType, pressure, Gain1, Gain2);


  Int_t ChMaskMB1[32]={ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
  Int_t ChMaskMB2[32]={ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};


  Int_t LastEvent = nEntries;

  if(EventsToProcess>0)
    LastEvent = FirstEvent+EventsToProcess;

  if( LastEvent > nEntries ) LastEvent = nEntries;

  cout << Form("\n\nEvents to be processed: from %d to %d", FirstEvent, LastEvent-1);

  cout << "\n\n";

  Int_t TB4MWPCTrigOffSet = 0;

  Float_t *ArrNtADRIANOAmp1 = new Float_t[Nch+1];
  Float_t *ArrNtADRIANOTimeAmp1 = new Float_t[Nch+1];
  Float_t *ArrNtADRIANOIntegral1 = new Float_t[Nch+1];
  Float_t *ArrNtADRIANOPed1 = new Float_t[Nch+1];
  Float_t *ArrNtADRIANOPedSigma1 = new Float_t[Nch+1];
  Float_t *ArrNtADRIANOAmp2 = new Float_t[Nch+1];
  Float_t *ArrNtADRIANOTimeAmp2 = new Float_t[Nch+1];
  Float_t *ArrNtADRIANOIntegral2 = new Float_t[Nch+1];
  Float_t *ArrNtADRIANOPed2Left = new Float_t[Nch+1];
  Float_t *ArrNtADRIANOPed2Right = new Float_t[Nch+1];
  Float_t *ArrNtADRIANOFitCheck = new Float_t[Nch+1];
  Float_t ArrNtADRIANOFitParams16[11];
  Float_t ArrNtADRIANOFitParams17[11];
  Float_t ArrNtADRIANOFitParams18[11];
  Float_t ArrNtADRIANOFitParams19[11];

  for(Int_t iEvt=FirstEvent;iEvt<LastEvent;iEvt++){

    for(Int_t idval=0; idval<Nch; idval++){
      ArrNtADRIANOAmp1[idval] = -999.;
      ArrNtADRIANOTimeAmp1[idval] = -999.;
      ArrNtADRIANOIntegral1[idval] = -999.;
      ArrNtADRIANOPed1[idval] = -999.;
      ArrNtADRIANOPedSigma1[idval] = -999.;
      ArrNtADRIANOAmp2[idval] = -999.;
      ArrNtADRIANOTimeAmp2[idval] = -999.;
      ArrNtADRIANOIntegral2[idval] = -999.;
      ArrNtADRIANOPed2Left[idval] = -999.;
      ArrNtADRIANOPed2Right[idval] = -999.;
      ArrNtADRIANOFitCheck[idval] = -999.;
    }

    for (Int_t idval=0; idval<11; idval++){
      ArrNtADRIANOFitParams16[idval]=0.;
      ArrNtADRIANOFitParams17[idval]=0.;
      ArrNtADRIANOFitParams18[idval]=0.;
      ArrNtADRIANOFitParams19[idval]=0.;
    }

    ArrNtADRIANOAmp1[Nch] = iEvt;
    ArrNtADRIANOTimeAmp1[Nch] = iEvt;
    ArrNtADRIANOIntegral1[Nch] = iEvt;
    ArrNtADRIANOPed1[Nch] = iEvt;
    ArrNtADRIANOPedSigma1[Nch] = iEvt;
    ArrNtADRIANOAmp2[Nch] = iEvt;
    ArrNtADRIANOTimeAmp2[Nch] = iEvt;
    ArrNtADRIANOIntegral2[Nch] = iEvt;
    ArrNtADRIANOPed2Left[Nch] = iEvt;
    ArrNtADRIANOPed2Right[Nch] = iEvt;
    ArrNtADRIANOFitCheck[Nch] = iEvt;
    ArrNtADRIANOFitParams16[10]=iEvt;
    ArrNtADRIANOFitParams17[10]=iEvt;
    ArrNtADRIANOFitParams18[10]=iEvt;
    ArrNtADRIANOFitParams19[10]=iEvt;


    Bool_t bKeep = 0;
    UShort_t uMWPCTrigStatus = 0;

    //===============================
    //       UShort_t GetMWPCxy(Double_t &TB4TimeEpoc, TTree *trMWPC, Int_t &MWPCTrig, Short_t &MWPCxwire, Short_t &MWPCywire, Short_t &MWPCxtdc, Short_t &MWPCytdc, Double_t &dTime, Short_t TB42MPWCTimeOffset=345);

    Double_t TB4TimeEpoc=0.;

    Int_t MWPCTrig=iEvt;
    Short_t MWPCxwire=-1, MWPCywire=-1, MWPCxtdc=-1, MWPCytdc=-1, TB42MPWCTimeOffset=345;
    Double_t dTime=-999.;

    Bool_t bIsNewEvent=1;

//     GetMWPCxy(TB4TimeEpoc, trMWPC, MWPCTrig, MWPCxwire, MWPCywire, MWPCxtdc, MWPCytdc, dTime, TB42MPWCTimeOffset);
//     TB4TimeEpoc, MWPCxwire, MWPCywire, MWPCxtdc, MWPCytdc, dTime

    //===============================


    if(!(((Int_t)((iEvt-FirstEvent)*LastEvent/(LastEvent - FirstEvent)))%(LastEvent/10))){
      if(LastEvent>100)
	cout << Form("progress... %d%%   %d/%d   %d\r", (Int_t)((iEvt-FirstEvent)*LastEvent/(LastEvent - FirstEvent)/(LastEvent/100)), (iEvt-FirstEvent), (LastEvent - FirstEvent), iEvt) << flush;
      else
	cout << Form("progress... %d%%   %d/%d   %d\r", (Int_t)((iEvt-FirstEvent)*LastEvent/(LastEvent - FirstEvent)/(LastEvent/100.)), (iEvt-FirstEvent), (LastEvent - FirstEvent), iEvt) << flush;
    }

    for(Int_t iich=0;iich<Nch;iich++){


      Int_t ich = chlist[iich];

      //       cout << Form("progress... %d%%   %d/%d   %d ch: %d\r", (Int_t)((iEvt-FirstEvent+1)*100/(LastEvent - FirstEvent)), (iEvt-FirstEvent+1), (LastEvent - FirstEvent), iEvt, ich) << flush;

      htmp->Reset();

      Int_t XMaxHisto = -999;
      Float_t XMaxFit = -999.;
      Double_t MaxVal1 = -999., MaxVal2 = -999.;
      Double_t Integral1 = -999., Integral2 = -999.;
      Double_t ped1=-999., ped2Left=-999., ped2Right=-999., sigma=-999.;
      Int_t ChCheck=0;


      DRCHit *drchit = 0x0;

      //***************** MB1 **********************//
      if(ich<16 && iEvt<nEntriesMB1){

	T1->GetEntry(iEvt);
	TClonesArray *DRCHits1=event1->GetDRCHits();
	if(!(DRCHits1->GetEntries())) continue;

	if(iEvt==FirstEvent)
	{
	  for(Int_t idxch=0; idxch<DRCHits1->GetEntries(); idxch++)
	    ChMaskMB1[((DRCHit*)DRCHits1->At(idxch))->GetCh()] = idxch;
	}

	TB4TimeEpoc = event1->GetEventTime();

	drchit = (DRCHit*)DRCHits1->At(ChMaskMB1[ich]);
	if(!drchit) continue;

	if(bUseMWPC && bIsNewEvent)
	{
	  MWPCTrig+=TB4MWPCTrigOffSet;
	  uMWPCTrigStatus = GetMWPCxy(TB4TimeEpoc, trMWPC, MWPCTrig, MWPCxwire, MWPCywire, MWPCxtdc, MWPCytdc, dTime, TB42MPWCTimeOffset);
	  if(uMWPCTrigStatus==2)
	    TB4MWPCTrigOffSet -= 1;
	  bIsNewEvent=0;
	}

      }


      //***************** MB2 **********************//
      else if(ich>15 && iEvt<nEntriesMB2){

	T2->GetEntry(iEvt);
	TClonesArray *DRCHits2=event2->GetDRCHits();
	if(!(DRCHits2->GetEntries())) continue;

	if(iEvt==FirstEvent)
	{
	  for(Int_t idxch=0; idxch<DRCHits2->GetEntries(); idxch++){
	    ChMaskMB2[((DRCHit*)DRCHits2->At(idxch))->GetCh()-16] = idxch;
	  }
	}

	TB4TimeEpoc = event2->GetEventTime();

	drchit = (DRCHit*)DRCHits2->At(ChMaskMB2[ich-16]);

	if(!drchit) continue;

	if(bUseMWPC && bIsNewEvent)
	{
	  MWPCTrig+=TB4MWPCTrigOffSet;
	  uMWPCTrigStatus = GetMWPCxy(TB4TimeEpoc, trMWPC, MWPCTrig, MWPCxwire, MWPCywire, MWPCxtdc, MWPCytdc, dTime, TB42MPWCTimeOffset);
	  if(uMWPCTrigStatus==2)
	    TB4MWPCTrigOffSet -= 1;
	  bIsNewEvent=0;
	}

	if(iEvt==FirstEvent)
	{
	  DataWfSize = (Float_t)drchit->GetSize();
	  RunBeginTime = event2->GetEventTime();
	}
	if(iEvt==LastEvent-1)
	{
	  RunEndTime = event2->GetEventTime();
	}

      }

      for(Int_t i1=0;i1<wfsize;i1++){
	htmp->Fill(i1,(drchit->GetW().At(i1)-2048.*bShift));
      }


      if( ich == 16 || ich == 17 )
      {
	XMaxHisto = htmp->GetMaximumBin();
	MaxVal1 = htmp->GetBinContent(XMaxHisto);

	Double_t AmpHisto = htmp->GetBinContent(htmp->GetMaximumBin()) - htmp->GetBinContent(htmp->GetMinimumBin());
	if(AmpHisto>12) ChCheck |= 1;

	if(TMath::Abs(XMaxHisto-XMaxList[iich])<20.) ChCheck |= 1<<1;

	GetPedestalPADE(drchit->GetW(),ped1, sigma, bShift);
	MaxVal1 = TMath::Abs(htmp->GetBinContent(XMaxHisto) - ped1);


      }
      else
      {
	if(IsSiPMChList[iich])
	  GetPedestalFEB(drchit->GetW(),ped1, sigma, bShift);
	else
	  GetPedestalFEB(drchit->GetW(),ped1, sigma, bShift);

	if(bFast)
	{
	  if(IsSiPMChList[iich])
	    XMaxHisto = htmp->GetMaximumBin();
	  else
	    XMaxHisto = htmp->GetMinimumBin();

	  MaxVal1 = TMath::Abs(htmp->GetBinContent(XMaxHisto) - ped1);
	  Integral1 = TMath::Abs(htmp->Integral(1,wfsize)-wfsize*ped1);

	  Double_t AmpHisto = htmp->GetBinContent(htmp->GetMaximumBin()) - htmp->GetBinContent(htmp->GetMinimumBin());

	  if(AmpHisto>12) ChCheck |= 1;

	  if(TMath::Abs(XMaxHisto-XMaxList[iich])<20.) ChCheck |= 1<<1;
	}
	else
	{

	  TF1 *fitsnr = new TF1("fitsnr",fitfn2,0.,wfsize,8);
	  fitsnr->SetParameters(60., 2000.-2048.*bShift, -20000.,-0.05,-5000., -0.02,  -0.001, 2010.-2048.*bShift);
	  for(int ifit=0;ifit<nfit;ifit++){
	    htmp->Fit("fitsnr", "B", "", 0., wfsize);
	    if(gMinuit->fEDM < 1.e-5) break;
	  }
	  for(int ipar = 0;ipar<8;ipar++){
	    if(ich == 16)
	      ArrNtADRIANOFitParams16[ipar] = fitsnr->GetParameter(ipar);
	    if(ich == 17)
	      ArrNtADRIANOFitParams17[ipar] = fitsnr->GetParameter(ipar);
	    if(ich ==18)
	      ArrNtADRIANOFitParams18[ipar] = fitsnr->GetParameter(ipar);
	    if(ich ==19)
	      ArrNtADRIANOFitParams19[ipar] = fitsnr->GetParameter(ipar);
	  }

	  if(ich == 16){
	    ArrNtADRIANOFitParams16[8] = fitsnr->GetChisquare()/fitsnr->GetNDF();
	    ArrNtADRIANOFitParams16[9] = gMinuit->fEDM;
	  }
	  if(ich == 17){
	    ArrNtADRIANOFitParams17[8] = fitsnr->GetChisquare()/fitsnr->GetNDF();
	    ArrNtADRIANOFitParams17[9] = gMinuit->fEDM;
	  }
	  if(ich ==18){
	    ArrNtADRIANOFitParams18[8] = fitsnr->GetChisquare()/fitsnr->GetNDF();
	    ArrNtADRIANOFitParams18[9] = gMinuit->fEDM;
	  }
	  if(ich ==19){
	    ArrNtADRIANOFitParams19[8] = fitsnr->GetChisquare()/fitsnr->GetNDF();
	    ArrNtADRIANOFitParams19[9] = gMinuit->fEDM;
	  }

	  // 	c1.Modified();
	  // 	c1.Update();
	  // 	TPaveStats *stats = (TPaveStats*)c1.GetPrimitive("stats");
	  // 	stats->SetY1NDC(.2);
	  // 	stats->SetY2NDC(.6);
	  // 	c1.Modified();
	  // 	c1.Update();
	  // 	c1.WaitPrimitive();


	  Bool_t status = 0;
	  if(gMinuit->fCstatu.EqualTo("CONVERGED ")){
	    status = 1;
	    iconverge++;
	  }
	  //cout<<status<<endl;

	  if(IsSiPMChList[iich])
	    XMaxHisto = htmp->GetMaximumBin();
	  else
	    XMaxHisto = htmp->GetMinimumBin();

	  MaxVal1 = TMath::Abs(htmp->GetBinContent(XMaxHisto) - ped1);
	  Integral1 = TMath::Abs(htmp->Integral(1,wfsize)-wfsize*ped1);

	  Double_t AmpHisto = htmp->GetBinContent(htmp->GetMaximumBin()) - htmp->GetBinContent(htmp->GetMinimumBin());

	  if(AmpHisto>12) ChCheck |= 1;

	  if(TMath::Abs(XMaxHisto-XMaxList[iich])<20.) ChCheck |= 1<<1;

	  ped2Left = fitsnr->GetParameter(1);
	  ped2Right = fitsnr->GetParameter(7);

	  if((ped2Right-ped2Left<50.)&&(ped2Right>ped2Left)) ChCheck |= 1<<2; //If the difference between the right and left pedestals is larger then 50.

	  if(IsSiPMChList[iich])
	    XMaxFit = fitsnr->GetMaximumX();
	  else
	    XMaxFit = fitsnr->GetMinimumX();

	  if(TMath::Abs(XMaxFit-XMaxList[iich])<20.) ChCheck |= 1<<3;

	  if(status) ChCheck |= 1<<4;

	  if(IsSiPMChList[iich])
	    MaxVal2 = TMath::Abs(fitsnr->GetMaximum() - ped2Left);
	  else
	    MaxVal2 = TMath::Abs(fitsnr->GetMinimum() - ped2Left);

	  Float_t finalPoint = fitsnr->GetX(ped2Left, fitsnr->GetParameter(0)+10, wfsize);
	  if(ped2Left<ped2Right){
	    Integral2 = -(fitsnr->Integral(0.,finalPoint)-ped2Left*finalPoint);
	    ChCheck |= 1<<5;
	  }else{
	    Integral2 = -999.;
	  }

	}
      }

      ArrNtADRIANOAmp1[iich] = MaxVal1;
      ArrNtADRIANOTimeAmp1[iich] = XMaxHisto;
      ArrNtADRIANOIntegral1[iich] = Integral1;
      ArrNtADRIANOPed1[iich] = ped1;
      ArrNtADRIANOPedSigma1[iich] = sigma;
      ArrNtADRIANOFitCheck[iich] = ChCheck;

      ArrNtADRIANOAmp2[iich] = MaxVal2;
      ArrNtADRIANOTimeAmp2[iich] = XMaxFit;
      ArrNtADRIANOIntegral2[iich] = Integral2;
      ArrNtADRIANOPed2Left[iich] = ped2Left;
      ArrNtADRIANOPed2Right[iich] = ped2Right;

      bKeep = 1;

    }

    Float_t ArrNtMWPC[8] = {(Float_t)TB4TimeEpoc, (Float_t)MWPCxwire, (Float_t)MWPCywire, (Float_t)MWPCxtdc, (Float_t)MWPCytdc, (Float_t)dTime, (Float_t)(!uMWPCTrigStatus), (Float_t)iEvt};

    if(bKeep){
      ntADRIANOAmp1->Fill(ArrNtADRIANOAmp1);
      ntADRIANOTimeAmp1->Fill(ArrNtADRIANOTimeAmp1);
      ntADRIANOIntegral1->Fill(ArrNtADRIANOIntegral1);
      ntADRIANOPed1->Fill(ArrNtADRIANOPed1);
      ntADRIANOPedSigma1->Fill(ArrNtADRIANOPedSigma1);
      ntADRIANOAmp2->Fill(ArrNtADRIANOAmp2);
      ntADRIANOTimeAmp2->Fill(ArrNtADRIANOTimeAmp2);
      ntADRIANOIntegral2->Fill(ArrNtADRIANOIntegral2);
      ntADRIANOPed2Left->Fill(ArrNtADRIANOPed2Left);
      ntADRIANOPed2Right->Fill(ArrNtADRIANOPed2Right);
      ntADRIANOFitCheck->Fill(ArrNtADRIANOFitCheck);
      ntADRIANOFitParams16->Fill(ArrNtADRIANOFitParams16);
      ntADRIANOFitParams17->Fill(ArrNtADRIANOFitParams17);
      ntADRIANOFitParams18->Fill(ArrNtADRIANOFitParams18);
      ntADRIANOFitParams19->Fill(ArrNtADRIANOFitParams19);

      if(bUseMWPC)
	ntMWPC->Fill(ArrNtMWPC);

    }

    if(debug) cout << endl;


  }

  delete[] ArrNtADRIANOAmp1;
  delete[] ArrNtADRIANOTimeAmp1;
  delete[] ArrNtADRIANOIntegral1;
  delete[] ArrNtADRIANOPed1;
  delete[] ArrNtADRIANOPedSigma1;
  delete[] ArrNtADRIANOAmp2;
  delete[] ArrNtADRIANOTimeAmp2;
  delete[] ArrNtADRIANOIntegral2;
  delete[] ArrNtADRIANOPed2Left;
  delete[] ArrNtADRIANOPed2Right;
  delete[] ArrNtADRIANOFitCheck;

  cout << endl;

  pfiled->cd(Form("%d",RunNo));

  Int_t NchOffSet=17;
  Float_t *ArrNtADRIANOParam = new Float_t[NchOffSet+Nch] ;
  ArrNtADRIANOParam[0] = Energy;
  ArrNtADRIANOParam[1] = x;
  ArrNtADRIANOParam[2] = z;
  ArrNtADRIANOParam[3] = Angle;
  ArrNtADRIANOParam[4] = BiasVoltage;
  ArrNtADRIANOParam[5] = cnt;
  ArrNtADRIANOParam[6] = PartType;
  ArrNtADRIANOParam[7] = pressure;
  ArrNtADRIANOParam[8] = nEntriesMB1;
  ArrNtADRIANOParam[9] = nEntriesMB2;
  ArrNtADRIANOParam[10] = DataWfSize;
  ArrNtADRIANOParam[11] = wfsize;
  ArrNtADRIANOParam[12] = Nch;
  ArrNtADRIANOParam[13] = RunBeginTime;
  ArrNtADRIANOParam[14] = RunEndTime;
  ArrNtADRIANOParam[15] = Gain1;
  ArrNtADRIANOParam[16] = Gain2;
  for(Int_t j=0; j<Nch; j++)
    ArrNtADRIANOParam[NchOffSet+j] = chlist[j];

  ntADRIANOParam->Fill(ArrNtADRIANOParam);


  ntADRIANOAmp1->Write(0,TObject::kOverwrite);
  ntADRIANOTimeAmp1->Write(0,TObject::kOverwrite);
  ntADRIANOIntegral1->Write(0,TObject::kOverwrite);
  ntADRIANOPed1->Write(0,TObject::kOverwrite);
  ntADRIANOPedSigma1->Write(0,TObject::kOverwrite);
  ntADRIANOFitCheck->Write(0,TObject::kOverwrite);
  if(!bFast){
    ntADRIANOAmp2->Write(0,TObject::kOverwrite);
    ntADRIANOTimeAmp2->Write(0,TObject::kOverwrite);
    ntADRIANOIntegral2->Write(0,TObject::kOverwrite);
    ntADRIANOPed2Left->Write(0,TObject::kOverwrite);
    ntADRIANOPed2Right->Write(0,TObject::kOverwrite);
    ntADRIANOFitParams16->Write(0,TObject::kOverwrite);
    ntADRIANOFitParams17->Write(0,TObject::kOverwrite);
    ntADRIANOFitParams18->Write(0,TObject::kOverwrite);
    ntADRIANOFitParams19->Write(0,TObject::kOverwrite);
  }

  if(bUseMWPC){
    ntMWPC->Write(0,TObject::kOverwrite);
    cout << Form("TB4 got %d more triggers than MWPC\n", -TB4MWPCTrigOffSet);
  }

  if(FirstEvent==0)
    ntADRIANOParam->Write(0,TObject::kOverwrite);

  gSystem->cd(sAnalysisDir.Data());

  pfiled->Close();
  cout<<iconverge<<endl;

}

UShort_t GetMWPCxy(Double_t &TB4TimeEpoc, TTree *trMWPC, Int_t &MWPCTrig, Short_t &MWPCxwire, Short_t &MWPCywire, Short_t &MWPCxtdc, Short_t &MWPCytdc, Double_t &dTime, Short_t TB42MPWCTimeOffset){

  if(!trMWPC) return 1;

  Double_t MWPCTimeEpoc=0.;
//   ULong64_t MWPCuiNanoSecond=0;

//   TFile *fileMWPC = TFile::Open(MWPCRunFile.Data());

//   TTree *trMWPC = (TTree *)gROOT->FindObject("trMWPC");

//   int MWPCxwire_local, MWPCywire_local, MWPCxtdc_local, MWPCytdc_local;

  // x-wire hits
  Int_t ichamber = 1;
  TString str1 = "xwire"; str1 += ichamber;
  trMWPC->SetBranchAddress(str1,&MWPCxwire);
  str1 = "xtime"; str1 += ichamber;
  trMWPC->SetBranchAddress(str1,&MWPCxtdc);

  // y-wire hits
  str1 = "ywire"; str1 += ichamber;
  trMWPC->SetBranchAddress(str1,&MWPCywire);
  str1 = "ytime"; str1 += ichamber;
  trMWPC->SetBranchAddress(str1,&MWPCytdc);

  // time hits
  str1 = "timeEPOC";
  trMWPC->SetBranchAddress(str1,&MWPCTimeEpoc);
//   str1 = "ns";
//   trMWPC->SetBranchAddress(str1,&MWPCuiNanoSecond);


//   Double_t TB4TimeEpoc=0.;

//   TFile *fileTB4 = TFile::Open(TB4RunFile.Data());
//   TTree *trTB4 = (TTree *)fileTB4->Get("DRCMB2");

//   DRCEvent *TB4Event = 0;
//   str1 = "DRCEvent";
//   trTB4->SetBranchAddress(str1, &TB4Event);



//   Long64_t TB4nEvents = trTB4->GetEntries();
//   Long64_t MWPCnEvents = trMWPC->GetEntries();
//   //   Long64_t nEvents = TMath::Max(TB4nEvents, MWPCnEvents);
//   Long64_t nEvents = MWPCnEvents;
//
//   cout << Form(" Events: TB4 %lld   MWPC %lld\n", TB4nEvents, MWPCnEvents);

//   Int_t TB4Trig=0;
//   Int_t TB4Trig_OK=0;
//   Int_t MWPCTrig=0;
//   for(MWPCTrig=0; MWPCTrig<nEvents; MWPCTrig++){
//     trTB4->GetEvent(TB4Trig);

//     TB4TimeEpoc = TB4Event->GetEventTime();

    trMWPC->GetEvent(MWPCTrig);

//     MWPCxwire = MWPCxwire_local;
//     MWPCywire = MWPCywire_local;
//     MWPCxtdc = MWPCxtdc_local;
//     MWPCytdc = MWPCytdc_local;

//     dTime = (MWPCTimeEpoc+MWPCuiNanoSecond/1.e9)-TB4TimeEpoc;
    dTime = MWPCTimeEpoc-TB4TimeEpoc;

    if(dTime > TB42MPWCTimeOffset)
      return 2;

//     while(dTime > 345.)
//     {
//       cout << Form("Event %d skipped (out of time %f)\n", TB4Trig, dTime);
//       TB4Trig++;
//       trTB4->GetEvent(TB4Trig);
//       TB4TimeEpoc = TB4Event->GetEventTime();
//       dTime = (MWPCTimeEpoc+MWPCuiNanoSecond/1.e9)-TB4TimeEpoc;
//     }


//     TB4Trig++;
//     TB4Trig_OK++;
//   }

//   cout << Form(" Events: TB4 %d (%d)  MWPC %d\n", TB4Trig-1, TB4Trig_OK, MWPCTrig);


    return 0;

}

#endif
