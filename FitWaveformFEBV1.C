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

R__LOAD_LIBRARY(DRCEventV5_C)
R__LOAD_LIBRARY(functions_C)

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

void FitWaveformFEBV1(Int_t RunNo, Int_t FirstEvent, TString path, Bool_t debug=0, Int_t wfsize = 127, Int_t nfit = 1)
{

    const Int_t Nch=1;
    Int_t chlist[Nch] = {3};

  gROOT->Reset();
  gStyle->SetOptFit(1111);
  new TMinuit;
  Int_t iconverge=0;

  TString sAnalysisDir(gSystem->pwd());

  TH1F *htmp = new TH1F("htmp","htmp",wfsize,0.,wfsize);

  const char *ext=".root";

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

	  break;
	}
      }
    }
  }

  if((input.empty()))
  {
    cout << Form("\n\n $$ No Run %d in your data $$\n", RunNo);
    return;
  }

  Bool_t bad=0;
  if (!(fname.EndsWith(ext)) )
    bad=1;

  if(bad){
    cout << Form("RunNo: %d don't exist in data\n", RunNo);
    return;
  }

  fname.Prepend(path);
  string RootFileName = string(fname);


  cout << Form("filename: %s RunNo: %d\n", RootFileName.data(), RunNo);
  TString sCurDir(gSystem->pwd());

  gSystem->cd(sAnalysisDir.Data());

  gSystem->cd(sCurDir.Data());



  cout << gSystem->pwd() << endl;
  TFile f(Form("%s",RootFileName.data()));

  if(f.IsZombie()){
    f.Close();
    cout << Form("filename: %s RunNo: %d is bad\n", RootFileName.data(), RunNo);
    bad=1;
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


  if( FirstEvent > nEntries ) FirstEvent = nEntries;

  cout << Form("\n\nEvents to be processed: %d", FirstEvent);

  cout << "\n\n";



  Int_t iEvt=FirstEvent;
  {

    Double_t dTime=-999.;

    Bool_t bIsNewEvent=1;

    //===============================

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

	drchit = (DRCHit*)DRCHits1->At(ChMaskMB1[ich]);
	if(!drchit) continue;

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

	drchit = (DRCHit*)DRCHits2->At(ChMaskMB2[ich-16]);

	if(!drchit) continue;

      }

      for(Int_t i1=0;i1<wfsize;i1++){
	htmp->Fill(i1,(drchit->GetW().At(i1)));
      }

	  TF1 *fitsnr = new TF1("fitsnr",fitfn2,0.,wfsize,8);
	  // fitsnr->SetParameters(60., 2000.-2048.*bShift, -20000.,-0.05,-5000., -0.02,  -0.001, 2010.-2048.*bShift);
          fitsnr->SetParameters(15., 5, 40000., -0.05, 20000., -0.02, -0.0001, -20);
	  for(int ifit=0;ifit<nfit;ifit++){
	    htmp->Fit("fitsnr", "B", "", 0., wfsize);
	    if(gMinuit->fEDM < 1.e-5) break;
	  }

          htmp->Draw("HIST");
          fitsnr->Draw("SAME");

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
	  cout<<status<<endl;


    }

    if(debug) cout << endl;

  }

  cout << endl;

  gSystem->cd(sAnalysisDir.Data());

  cout<<iconverge<<endl;

}
