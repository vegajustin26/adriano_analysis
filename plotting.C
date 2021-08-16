//master plotting script
/* This script was used by Justin Vega in Summer 2021 as part of the SULI Program.
You can use it in a ROOT notebook by importing #include "plotting.C", or use these as individual macros within another script

Be sure to modify the path to where your data is at the start of each function

These functions take ROOT files that are in the format of the SignalIntegralFebV3 script and plots their various ntuples

Here are short descriptions of what each one does:
DataAnalysisV1 - shows amplitude, time amplitude, pedestal, and pedestal sigma of a particular channel of a run, with optional fit
ch_amp - makes plots of amplitudes for each specified channel for up to 8 channels, optional fit
time_amp - same thing as ch_amp, but for time amplitude
ch_ped - '', but for pedestal
ch_ped_sigma - '', but for pedestal sigma
all_ch_amp - shows plot of 32 channels for a given run (might be broken? might have something to do with the fitting)
run_stats - shows the beam parameters of a particular run
merge_runs - merges similar runs together (but you have to figure out which ones are similar)
*/
// #ifndef ROOT_functions
// #define ROOT_functions

R__LOAD_LIBRARY(DRCEventV5_C)
R__LOAD_LIBRARY(functions_C)

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TArrayI.h>
#include <TClassTable.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"

#include "DRCEventV5.h"
#include "functions.C"

//DATA ANALYSIS
void DataAnalysisV1(Int_t ch, Int_t run, int fit, int merged, int *bin_params, int *timeamp_cuts, TString SignalFile)
{
//
// ./sigint/SignalIntegralFEBV3_2019_ 1_proton_fitfn2.root
  gSystem->cd("/Users/justinvega/Documents/Fermilab/src/adriano_analysis/");

  gStyle->SetOptFit(1111);

  if( gROOT->FindObject(SignalFile.Data()) )
    ((TFile *)gROOT->FindObject(SignalFile.Data()))->Close(); //file location


  TFile *file = TFile::Open(SignalFile.Data());

  if (!merged){
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
}

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
  if (ch > 15){ //plastic 125, 0-800
    ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(%d,0.,%d)", ch, ch, bin_params[2], bin_params[3]), Form("ADRIANOTimeAmp1.ch%d > %d && ADRIANOTimeAmp1.ch%d < %d && ADRIANOPedSigma1.ch%d < 3", ch, timeamp_cuts[2]-1, ch, timeamp_cuts[3]+1, ch));
    //ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch16>>hInteg3(125,0.,800.)"),Form("ADRIANOTimeAmp1.ch16>22&&ADRIANOTimeAmp1.ch16<25&&ADRIANOPedSigma1.ch16<3."),"")
  }
  else if (ch < 15){ //glass 200, 0-400
      ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(%d,0.,%d)", ch, ch, bin_params[0], bin_params[1]), Form("ADRIANOTimeAmp1.ch%d > %d && ADRIANOTimeAmp1.ch%d < %d && ADRIANOPedSigma1.ch%d < 10", ch, timeamp_cuts[0]-1, ch, timeamp_cuts[1]+1, ch));
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
      xmax = htmp->GetBinLowEdge(xBinMax+1); //x value of maximum bin
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
    //cout << Form("%f, %f, %f, %d", fitparams[1], fiterrors[1], chi, ndf) << endl;
}

  //==============================
  c1->cd(2);
    ADRIANOAmp1->Draw(Form("ADRIANOTimeAmp1.ch%d>>htmp1(300,0.,0.)", ch));
    TH1F *htmp1 = ((TH1F *)gROOT->FindObject("htmp1"));
    htmp1->SetTitle(Form("CH %d Time Amp", ch));

  //==============================
  c1->cd(3);
    ADRIANOAmp1->Draw(Form("ADRIANOPed1.ch%d>>htmp2(100,0.,0.)", ch), Form("ADRIANOPed1.ch%d > 0 && ADRIANOPed1.ch%d < 10", ch, ch));
    TH1F *htmp2 = ((TH1F *)gROOT->FindObject("htmp2"));
    htmp2->SetTitle(Form("CH %d Ped", ch));


  //==============================
  c1->cd(4);
    ADRIANOAmp1->Draw(Form("ADRIANOPedSigma1.ch%d>>htmp3(100,0.,0.)", ch), Form("ADRIANOPedSigma1.ch%d < 5 && ADRIANOPedSigma1.ch%d > 0", ch, ch));
    TH1F *htmp3 = ((TH1F *)gROOT->FindObject("htmp3"));
    htmp3->SetTitle(Form("CH %d Ped Sigma", ch));

    c1->Draw();
    //c1->SaveAs(Form("./RUN_%d_ch%d_ogfitfn2.png",run, ch));
}


//Plot of CHANNEL AMPLITUDES
void ch_amp(vector<int> ch, Int_t run, int fit, int smooth, int merged, int *bin_params, int *timeamp_cuts, TString SignalFile){

  gSystem->cd("/Users/justinvega/Documents/Fermilab/src/adriano_analysis/");

  gStyle->SetOptFit(1111);

  if( gROOT->FindObject(SignalFile.Data()) )
    ((TFile *)gROOT->FindObject(SignalFile.Data()))->Close(); //file location


  TFile *file = TFile::Open(SignalFile.Data());

  if (!merged){
  file->cd(Form("%d",run)); //go into folder
  }

  TTree *ADRIANOAmp1 = (TTree *) gROOT->FindObject("ADRIANOAmp1"); //finds amplitude TTree
    ADRIANOAmp1->AddFriend("ADRIANOFitCheck");
    ADRIANOAmp1->AddFriend("ADRIANOPedSigma1");
    ADRIANOAmp1->AddFriend("ADRIANOTimeAmp1");

  if( gROOT->FindObject("c1") )
    ((TCanvas *)gROOT->FindObject("c1"))->Close();

  TCanvas *c1 = new TCanvas("c1","c1",1000,600);

  c1->Divide(4,2);
   // TLatex *   tex = new TLatex(0.23,0.97, Form("Run %d", run));
   // tex->SetNDC();
   // tex->SetTextSize(0.025);
   // tex->Draw();

   int glass_ch[] = {4, 5, 6, 7, 13, 14, 21, 22, 29, 30};
   int pstic_ch[] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32};

    for (int i = 0; i < ch.size(); i++){
        c1->cd(i+1);
//         original
//         //ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(80,0.,0.)", ch[i], ch[i]), Form("(ADRIANOAmp1.ch%d < 100) && (ADRIANOFitCheck.ch%d&16) != 0", ch[i], ch[i]));
//         ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(300,0.,0.)", ch[i], ch[i]), Form("ADRIANOAmp1.ch%d < 600", ch[i]));
//         TH1F *htmp = ((TH1F *)gROOT->FindObject(Form("htmp%d", ch[i])));
//         htmp->SetTitle(Form("CH %d Amp", ch[i]));

        int glass = 0;
        int plastic = 0;

        if (std::find(std::begin(glass_ch), std::end(glass_ch), ch[i]) != std::end(glass_ch)){
          glass = 1;
        }
        else if (std::find(std::begin(pstic_ch), std::end(pstic_ch), ch[i]) != std::end(pstic_ch)){
          plastic = 1;
        }
        else {
          cout << "no compatibility for your channel(s) of interest" << endl;
        }

        //smoothing + fitting
        if (plastic){ //plastic 125, 0-800
          ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(%d,0.,%d)", ch[i], ch[i], bin_params[2], bin_params[3]), Form("ADRIANOTimeAmp1.ch%d > %d && ADRIANOTimeAmp1.ch%d < %d && ADRIANOPedSigma1.ch%d < 3", ch[i], timeamp_cuts[2]-1, ch[i], timeamp_cuts[3]+1, ch[i]));
          //ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch16>>hInteg3(125,0.,800.)"),Form("ADRIANOTimeAmp1.ch16>22&&ADRIANOTimeAmp1.ch16<25&&ADRIANOPedSigma1.ch16<3."),"")
        }
        else if (glass){ //glass 200, 0-400
            ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(%d,0.,%d)", ch[i], ch[i], bin_params[0], bin_params[1]), Form("ADRIANOTimeAmp1.ch%d > %d && ADRIANOTimeAmp1.ch%d < %d && ADRIANOPedSigma1.ch%d < 10", ch[i], timeamp_cuts[0]-1, ch[i], timeamp_cuts[1]+1, ch[i]));
        }

        TH1F *htmp = ((TH1F *)gROOT->FindObject(Form("htmp%d", ch[i])));
        if(smooth){
          if (ch[i] > 15){ //plastic
              htmp->Smooth(3);
          }
        }
        Double_t xmin = htmp->GetBinLowEdge(1);
        Double_t xBinMax = htmp->FindLastBinAbove(1.);
        Double_t xmax = htmp->GetBinLowEdge(xBinMax+1);

        htmp->Draw();
        if (!merged){
        htmp->SetTitle(Form("Amplitude ch %d run %d", ch[i], run));
        }
        else if (merged){
        htmp->SetTitle(Form("Amplitude ch %d merged", ch[i]));
        }

//         if(!fit3){
//         xBinMax = htmp->GetMaximumBin(); //bin of maximum value
//         xmax = htmp->GetBinLowEdge(xBinMax); //x value of maximum bin
//         Double_t nbins = htmp->GetNbinsX(); //# of bins
//         Double_t yExpo = htmp->GetBinContent(5); //?
//         Double_t Area = htmp->Integral()*nbins/10.;
//         Double_t FitRangeA[2]   = {htmp->GetBinLowEdge(1), htmp->GetBinLowEdge(htmp->FindLastBinAbove(1.))}; //x-axis fit range
//         Double_t StartValueA[4]    = {10.,      xmax,       Area,   2.};// yExpo, -1.e-3};
//         Double_t LowLimitValueA[4] = { 1.,  xmax/10.,   Area/10.,   0.};// TMath::Min(1., yExpo/10.), -1.e-1};
//         Double_t HiLimitValueA[4]  = {1.e2, xmax*10.,  Area*1.e3, 1.e2};// TMath::Max(1.e4, yExpo*10.), 0.};

//         TF1 *fitsnrA = langaufit(htmp, FitRangeA, StartValueA, LowLimitValueA, HiLimitValueA);
//             }
        if (fit){
            xBinMax = htmp->GetMaximumBin(); //bin of maximum value
            xmax = htmp->GetBinLowEdge(xBinMax); //x value of maximum bin
            Double_t nbins = htmp->GetNbinsX(); //# of bins
            Double_t Area = htmp->Integral();
            Double_t FitRangeA[2]   = {htmp->GetBinLowEdge(1), htmp->GetBinLowEdge(nbins+1)}; //x-axis fit range
            Double_t StartValueA[4]    = {1.,      xmax,       Area,   2.};
            Double_t LowLimitValueA[4] = {0.1,  xmax-10.,   Area/10.,   0.};
            Double_t HiLimitValueA[4]  = {100, xmax+10.,  Area*1.e3, 100};
            Double_t fitparams[4];
            Double_t fiterrors[4];
            Double_t chi;
            Int_t ndf;


        TF1 *fitsnrA = langaufit(htmp, FitRangeA, StartValueA, LowLimitValueA, HiLimitValueA, fitparams, fiterrors, &chi, &ndf);
        // cout << htmp->GetBinLowEdge(1) << endl;
        // cout << xBinMax << endl;
        // cout << xmax << endl;
        // cout << htmp->GetBinLowEdge(nbins+1) << endl;
        cout << Form("%d, %f, %f, %f, %f, %f, %f, %f, %d", (Int_t)htmp->GetEntries(), htmp->GetMean(), htmp->GetStdDev(), fitparams[0], fiterrors[0], fitparams[1], fiterrors[1], chi, ndf) << endl;

        }
    c1->Modified();
    c1->Draw();
  }
}

void time_amp(vector<int> ch, Int_t run, TString SignalFile="./sigint/histo/SignalIntegralFEBV3_1025_data.root"){
//./sigint/SignalIntegralFEBV3_2019_ 1_fitfn2.root
  gSystem->cd("/Users/justinvega/Documents/Fermilab/src/adriano_analysis/");

  gStyle->SetOptFit(1111);

  if( gROOT->FindObject(SignalFile.Data()) )
    ((TFile *)gROOT->FindObject(SignalFile.Data()))->Close(); //file location


  TFile *file = TFile::Open(SignalFile.Data());
  file->cd(Form("%d",run)); //go into folder

  TTree *ADRIANOTimeAmp1 = (TTree *) gROOT->FindObject("ADRIANOTimeAmp1"); //finds amplitude TTree
    //TTree *ADRIANOFitCheck = (TTree *) gROOT->FindObject("ADRIANOFitCheck");
    ADRIANOTimeAmp1->AddFriend("ADRIANOFitCheck");

if( gROOT->FindObject("c2") )
    ((TCanvas *)gROOT->FindObject("c2"))->Close();

  TCanvas *c2 = new TCanvas("c2","c2",1000,600);

  c2->Divide(4,2);
    TLatex *   tex = new TLatex(0.23,0.97, Form("Run %d", run));
   tex->SetNDC();
   tex->SetTextSize(0.025);
   tex->Draw();
    for (int i = 0; i < ch.size(); i++){
        c2->cd(i+1);
        ADRIANOTimeAmp1->Draw(Form("ADRIANOTimeAmp1.ch%d>>htmp%d(300,0.,0.)", ch[i], ch[i]));
        TH1F *htmp = ((TH1F *)gROOT->FindObject(Form("htmp%d", ch[i])));
        htmp->SetTitle(Form("CH %d Time Amp", ch[i]));
        c2->Modified();
    c2->Draw();
}

}

void ch_ped(vector<int> ch, Int_t run, TString SignalFile="./sigint/histo/SignalIntegralFEBV3_1025_data.root"){

  gSystem->cd("/Users/justinvega/Documents/Fermilab/src/adriano_analysis/");

  gStyle->SetOptFit(1111);

  if( gROOT->FindObject(SignalFile.Data()) )
    ((TFile *)gROOT->FindObject(SignalFile.Data()))->Close(); //file location


  TFile *file = TFile::Open(SignalFile.Data());
  file->cd(Form("%d",run)); //go into folder

  TTree *ADRIANOPed1 = (TTree *) gROOT->FindObject("ADRIANOPed1"); //finds amplitude TTree
    //ADRIANOAmp1->AddFriend("ADRIANOFitCheck");

  if( gROOT->FindObject("c3") )
    ((TCanvas *)gROOT->FindObject("c3"))->Close();

  TCanvas *c3 = new TCanvas("c3","c3",1000,600);

  c3->Divide(4,2);
   TLatex *   tex = new TLatex(0.23,0.97, Form("Run %d", run));
   tex->SetNDC();
   tex->SetTextSize(0.025);
   tex->Draw();
    for (int i = 0; i < ch.size(); i++){
        c3->cd(i+1);
        //ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(80,0.,0.)", ch[i], ch[i]), Form("(ADRIANOAmp1.ch%d < 100) && (ADRIANOFitCheck.ch%d&16) != 0", ch[i], ch[i]));
        ADRIANOPed1->Draw(Form("ADRIANOPed1.ch%d>>htmp%d(300,0.,0.)", ch[i], ch[i]), Form("ADRIANOPed1.ch%d > -25 && ADRIANOPed1.ch%d < 15", ch[i], ch[i]));
        TH1F *htmp = ((TH1F *)gROOT->FindObject(Form("htmp%d", ch[i])));
        htmp->SetTitle(Form("CH %d Ped", ch[i]));
        c3->Modified();
    c3->Draw();
}

}

void ch_ped_sigma(vector<int> ch, Int_t run, TString SignalFile="./sigint/histo/SignalIntegralFEBV3_1025_data.root")
{

  gSystem->cd("/Users/justinvega/Documents/Fermilab/src/adriano_analysis/");

  gStyle->SetOptFit(1111);

  if( gROOT->FindObject(SignalFile.Data()) )
    ((TFile *)gROOT->FindObject(SignalFile.Data()))->Close(); //file location


  TFile *file = TFile::Open(SignalFile.Data());
  file->cd(Form("%d",run)); //go into folder

  TTree *ADRIANOPedSigma1 = (TTree *) gROOT->FindObject("ADRIANOPedSigma1"); //finds amplitude TTree
    //ADRIANOAmp1->AddFriend("ADRIANOFitCheck");

  if( gROOT->FindObject("c4") )
    ((TCanvas *)gROOT->FindObject("c4"))->Close();

  TCanvas *c4 = new TCanvas("c4","c4",1000,600);

  c4->Divide(4,2);
   TLatex *   tex = new TLatex(0.23,0.97, Form("Run %d", run));
   tex->SetNDC();
   tex->SetTextSize(0.025);
   tex->Draw();
    for (int i = 0; i < ch.size(); i++){
        c4->cd(i+1);
        //ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(80,0.,0.)", ch[i], ch[i]), Form("(ADRIANOAmp1.ch%d < 100) && (ADRIANOFitCheck.ch%d&16) != 0", ch[i], ch[i]));
        ADRIANOPedSigma1->Draw(Form("ADRIANOPedSigma1.ch%d>>htmp%d(300,0.,0.)", ch[i], ch[i]));//, Form("ADRIANOAmp1.ch%d < 1000", ch[i]));
        TH1F *htmp = ((TH1F *)gROOT->FindObject(Form("htmp%d", ch[i])));
        htmp->SetTitle(Form("CH %d Ped Sigma", ch[i]));
        c4->Modified();
    c4->Draw();
  }
}

//Plot of all CHANNEL AMPLITUDES for protons
void all_ch_amp(vector<int> ch, vector<int> runslist, int fit, int smooth, TString SignalFile){

  gSystem->cd("/Users/justinvega/Documents/Fermilab/src/adriano_analysis/");
    gStyle->SetLineScalePS(0.5);
  gStyle->SetOptFit(1111);
    //gStyle->SetLineWidth(-1);
  if( gROOT->FindObject(SignalFile.Data()) )
    ((TFile *)gROOT->FindObject(SignalFile.Data()))->Close(); //file location

    if( gROOT->FindObject("allplots") )
        ((TCanvas *)gROOT->FindObject("allplots"))->Close();

        TCanvas *allplots = new TCanvas("allplots","allplots", 2000, 6250);
        allplots->SetLineWidth(1);
        allplots->Divide(8,34, 0.01, 0.000001);

  TFile *file = TFile::Open(SignalFile.Data());
    int idxrun = 0;

    for (int run = 0; run < runslist.size(); run++){ //for all runs
        file->cd(Form("%d",runslist[run])); //go into folder

        TTree *ADRIANOAmp1 = (TTree *) gROOT->FindObject("ADRIANOAmp1"); //finds amplitude TTree
        ADRIANOAmp1->AddFriend("ADRIANOFitCheck");
        ADRIANOAmp1->AddFriend("ADRIANOPedSigma1");
        ADRIANOAmp1->AddFriend("ADRIANOTimeAmp1");

        for (int i = 0; i < ch.size(); i++){ //for all channels
            allplots->cd(idxrun+1);

            //smoothing + fitting
            if (ch[i] > 15){ //plastic
              ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(250,0.,0.)", ch[i], ch[i]));
            }
            else { //glass
                ADRIANOAmp1->Draw(Form("ADRIANOAmp1.ch%d>>htmp%d(1000,0.,0.)", ch[i], ch[i]));//, Form("ADRIANOAmp1.ch%d < 600", ch[i]));
            }

            TH1F *htmp = ((TH1F *)gROOT->FindObject(Form("htmp%d", ch[i])));
            if(smooth){
              if (ch[i] > 15){ //for plastic
                  htmp->Smooth(3);
              }
            }
            Double_t xmin = htmp->GetBinLowEdge(1);
            Double_t xBinMax = htmp->FindLastBinAbove(1.);
            Double_t xmax = htmp->GetBinLowEdge(xBinMax+1);

            htmp->Draw();
            htmp->SetTitle(Form("Amplitude ch %d run %d", ch[i], runslist[run]));

            // if (fit){
            //     xBinMax = htmp->GetMaximumBin(); //bin of maximum value
            //     xmax = htmp->GetBinLowEdge(xBinMax); //x value of maximum bin
            //     Double_t nbins = htmp->GetNbinsX(); //# of bins
            //     Double_t yExpo = htmp->GetBinContent(5); //?
            //     Double_t Area = htmp->Integral()*nbins/10.;
            //     Double_t FitRangeA[2]   = {htmp->GetBinLowEdge(1), htmp->GetBinLowEdge(htmp->FindLastBinAbove(1.))}; //x-axis fit range
            //     Double_t StartValueA[4]    = {10.,      xmax,       Area,   2.};
            //     Double_t LowLimitValueA[4] = { 1.,  xmax/10.,   Area/10.,   0.};
            //     Double_t HiLimitValueA[4]  = {1.e2, xmax*10.,  Area*1.e3, 1.e2};
            //
            // TF1 *fitsnrA = langaufit(htmp, FitRangeA, StartValueA, LowLimitValueA, HiLimitValueA);
            // }
        idxrun++;
        }
    }
    allplots->Modified();
    allplots->Draw();
}

void run_stats(int run, TString SignalFile){
    gSystem->cd("/Users/justinvega/Documents/Fermilab/src/adriano_analysis/");

    gStyle->SetOptFit(1111);

    std::vector<int> protonlist;

    if( gROOT->FindObject(SignalFile.Data()) )
        ((TFile *)gROOT->FindObject(SignalFile.Data()))->Close(); //file location

    TFile *file = TFile::Open(SignalFile.Data());

        TDirectory *dir = file->GetDirectory(Form("%d", run));
        if (dir){
            file->cd(Form("%d", run)); //go into folder
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

            const char *PartName[4]={"pion", "proton", "muon","n.a."};
            if(PartType<1 || PartType>3) PartType = 4;

              cout << "Run | beam |Energy|   X  |   Z  |Angle| NEvents| psia| BiasV| counts| gain1| gain2\n";
              cout << Form("%3d |%6s|%sGeV| %4d | %4d |%4d | %6d | %4.1f| %2.1fV| %4dK| %3d  | %3d \n", run, PartName[(Int_t)PartType-1], (Energy>0.9)?Form("%3d",(Int_t)Energy):Form("%1.1f",Energy), (Int_t)x, (Int_t)z, (Int_t)Angle, (Int_t)NEvents, pressure, BiasVoltage, (Int_t)cnt, (Int_t)(Gain1/1e4), (Int_t)(Gain2/1e4));

//             // identifying runs with protons
//             TString particle = PartName[(Int_t)PartType-1];
//             if (particle.Contains("proton")){
//                 protonlist.push_back(run);
//             }
//             }
//         else {
//             continue;
//         }
    }
//             for (int i = 0; i < protonlist.size(); i++){
//                 cout << Form("%d, ", protonlist[i]);
//                    }
    }

void merge_runs(vector<int> runslist, TString input = "./sigint/histo/SignalIntegralFEBV3_protonruns_data.root"){
//similar runs #1: x: 1302, z: 274, gain 1&2: 900

  gSystem->cd("/Users/justinvega/Documents/Fermilab/src/adriano_analysis/");

  //grabs beam parameters from first run in runslist(for filename)
  if( gROOT->FindObject(input.Data()) )
          ((TFile *)gROOT->FindObject(input.Data()))->Close(); //file location

      TFile *file = TFile::Open(input.Data());

      file->cd(Form("%d", runslist[0])); //go into folder
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

      const char *PartName[4]={"pion", "proton", "muon","n.a."};
      if(PartType<1 || PartType>3) PartType = 4;

  //         cout << "Run | beam |Energy|   X  |   Z  |Angle| NEvents| psia| BiasV| counts| gain1| gain2\n";
  //         cout << Form("%3d |%6s|%sGeV| %4d | %4d |%4d | %6d | %4.1f| %2.1fV| %4dK| %3d  | %3d \n", runslist[0], PartName[(Int_t)PartType-1], (Energy>0.9)?Form("%3d",(Int_t)Energy):Form("%1.1f",Energy), (Int_t)x, (Int_t)z, (Int_t)Angle, (Int_t)NEvents, pressure, BiasVoltage, (Int_t)cnt, (Int_t)(Gain1/1e4), (Int_t)(Gain2/1e4));

  file->Close();

  gSystem->cd("/Users/justinvega/Documents/Fermilab/src/adriano_analysis/");

  //make chains for each tree
  TChain run1("ADRIANOAmp1");
  TChain runtime("ADRIANOTimeAmp1");
  TChain runped("ADRIANOPed1");
  TChain runpedsigma("ADRIANOPedSigma1");
  TChain runfit("ADRIANOFitCheck");

  for (int runs=0; runs<runslist.size(); runs++){
      run1.Add(Form("%s/%d/ADRIANOAmp1", input.Data(), runslist[runs]));
      runtime.Add(Form("%s/%d/ADRIANOTimeAmp1", input.Data(), runslist[runs]));
      runped.Add(Form("%s/%d/ADRIANOPed1", input.Data(), runslist[runs]));
      runpedsigma.Add(Form("%s/%d/ADRIANOPedSigma1", input.Data(), runslist[runs]));
      runfit.Add(Form("%s/%d/ADRIANOFitCheck", input.Data(), runslist[runs]));
  }

  TString newrun = Form("RUN_FEB_merged_%d-%d_%6s_%sGeV_%2.1fV_%3d-%3d_%dDEG_x%4d_z%3d.root", runslist[0], runslist.back(), PartName[(Int_t)PartType-1], (Energy>0.9)?Form("%3d",(Int_t)Energy):Form("%1.1f",Energy), BiasVoltage, (Int_t)(Gain1/1e4), (Int_t)(Gain2/1e4), (Int_t)Angle, (Int_t)x, (Int_t)z);

  cout << Form("Creating new file: %s", newrun.Data()) << endl;

  TFile *outfile = TFile::Open(Form("./merged/%s", newrun.Data()),"RECREATE");

  //write each tree to outfile
  run1.CloneTree(-1,"fast");
  runtime.CloneTree(-1, "fast");
  runped.CloneTree(-1, "fast");
  runpedsigma.CloneTree(-1, "fast");
  runfit.CloneTree(-1, "fast");

  outfile->Write();
  cout << "complete" << endl;

  }
