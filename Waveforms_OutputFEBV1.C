/****
 * ROOT macro to visualize data waveform from ROOT file
 * this opens 2 TCanvas' the first with channels from 0-15
 * the second with channels from 16-31
 * i.e. this supports up to 2 FEB boards
 *
 * To run this macro it is needed to compile DRCEventV5.C
 * At ROOT prompt run
 * root [0] .L DRCEventV5.C++
 * root [1] .x WaveformFEBV1.C
 *
 * If DRCEventV5.C has already been compiled, it is enough to run
 * root [0] .x WaveformFEBV1.C
 *
 * This macro can be also executed directly from the terminal
 * root WaveformFEBV1.C
 *

 Justin Note:
  This is mostly the same as WaveformFEBV1, except some differences:
  1: this generates plots for all 64 channels of a SiPM
  2: you can see which channels are above a certain threshold by enabling 'cherrypick'. This will evaluate if each event is above a certain threshold for each channel.
  If there are waveforms above that threshold in a channel, then it will log the run and the channels of interest in a log.txt file

 ****/

#include "TFile.h"
#include "TTree.h"
#include <Riostream.h>
#include <TSystemDirectory.h>
#include <TString.h>
#include <TError.h>
#include <TObjString.h>
#include "TSystem.h"
#include "TROOT.h"
#include "TRootCanvas.h"
#include <TMath.h>
#include<vector>

#include "DRCEventV5.h"
#include <TStyle.h>
#include "TCanvas.h"
#include "TH1.h"

R__LOAD_LIBRARY(DRCEventV5_C)
using namespace std;


void Waveforms_Output(int RunNo, TString Path, int plot, int cherrypick)
{
  gROOT->Delete("fname");
  gROOT->Reset();

  TSystemDirectory dir("data",Path.Data());//TSystemDirectory dir("data","/data/ADRIANO2/TB_2019_dec/")
  int start=20;        //First event
  int nevt=2;       //how many events
  double hmin=-500.;   //min val in vertical axis
  double hmax=2500.; //max val in vertical axis
  bool persist = 1;   //0 = single waveforms; 1 = superimposed waveforms
  RunNo =  RunNo; //1009; //1196;//810;   //RUN to select
  int ncheck = 10;   //print event number every ncheck events


  TString sCurDir(gSystem->pwd());


  gStyle->SetOptStat(0);

  #if defined(__linux__)
  gSystem->Load("DRCEventV5_C.so");
  #else
  gSystem->Load("DRCEventV5_C.dll");
  #endif

  const char *ext=".root";
  TString fname;

  TList*files = dir.GetListOfFiles();

  // finds relevant .data file of interest and stores it
  if (files) {
    TSystemFile *file;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext)) {
        if(fname.Contains(Form("RUN_FEB_%d",RunNo))){
          cout << Form("Filename: %s", fname.Data()) << endl;
          break;
        }
      }
    }
  }

  fname.Prepend(dir.GetTitle());
  cout << Form("File Location: %s", fname.Data()) << endl;

  TFile f(fname.Data(),"","fname");

  // added MB3 and MB4
  TTree *T1 = (TTree*)f.Get("DRCMB1");
  DRCEvent *event1 = 0;
  T1->SetBranchAddress("DRCEvent", &event1);

  TTree *T2 = (TTree*)f.Get("DRCMB2");
  DRCEvent *event2 = 0;
  T2->SetBranchAddress("DRCEvent", &event2);

  TTree *T3 = (TTree*)f.Get("DRCMB3");
  DRCEvent *event3 = 0;
  T3->SetBranchAddress("DRCEvent", &event3);

  TTree *T4 = (TTree*)f.Get("DRCMB4");
  DRCEvent *event4 = 0;
  T4->SetBranchAddress("DRCEvent", &event4);

  // added MB3 and MB4
  Int_t nEntriesMB1 = T1->GetEntries();
  Int_t nEntriesMB2 = T2->GetEntries();
  Int_t nEntriesMB3 = T3->GetEntries();
  Int_t nEntriesMB4 = T4->GetEntries();

  // variables for cherrypick
  Double_t ch_max = 0; //initialize
  std::vector<double> ch_list; //records channels of interest
  Double_t threshold = 250; //threshold for interesting channel

  // cout << nEntriesMB3 << endl;
  // cout << nEntriesMB4 << endl;
  // //cout << (!(nEntriesMB3 && nEntriesMB4)) << endl;

  //return(Form("Run %d has %d MB3 entries, %d MB4 entries", RunNo, nEntriesMB3, nEntriesMB4));

  if (!(nEntriesMB3 && nEntriesMB4)) // if the data does not have MB3 and MB4
  {
    // code for MB1 and MB2
    Int_t nEntries = TMath::Max(nEntriesMB1, nEntriesMB2);

    TCanvas c1(Form("WaveformPersistentFEB1_RUN_%d",RunNo),Form("Waveform %s",fname.Data()),1500,900);
    c1.Divide(4,4);
    TCanvas c2(Form("WaveformPersistentFEB2_RUN_%d",RunNo),Form("Waveform %s",fname.Data()),1500,900);
    c2.Divide(4,4);

    Int_t ChToArr[32]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};

    TH1D *htmpM[33];
    for(Int_t ch=0; ch<33;ch++){
      htmpM[ch] = new TH1D(Form("htmp%.2d%.6d", ch, start),Form("ch %.2d",ch+1),127,0.,127.);
    }
    gStyle->SetTitleFontSize(0.1);

    Float_t MinY[32];
    Float_t MaxY[32];

    for(Int_t ch=0; ch<32;ch++){
      MinY[ch]=1.e5; //initialize
      MaxY[ch]=-1.;
    }

    if (plot){ //for plotting capabilities
      for(int evt=start; evt<start+nevt; evt++){ //for each event
        if(evt>nEntries){
          cout << Form("event %d not available, there are only %d events in RUN %d\n", evt, nEntries, RunNo);
          continue;
        }
        if(!(evt%ncheck)) cout <<Form("evt: %d\n", evt);


      //************************* MB1 *************************//
      T1->GetEntry(evt);
      TClonesArray *DRCHits1=event1->GetDRCHits();
      if(evt<nEntriesMB1){
        for(Int_t ch=2; ch<3;ch++){ //for each channel Int_t ch=0; ch<16;ch++)
           DRCHit *drchit1 = (DRCHit*)DRCHits1->At(ch); //assigns channel number to GetCh()
  	       if(!drchit1) continue; //if no event, skip loop once
           if(evt==0){
              cout << Form("ch: %d %d\n", ChToArr[ch], drchit1->GetCh());
            }
           if(persist){ //if you want waveform persistent/superimposed in same canvas
          	  htmpM[drchit1->GetCh()] = new TH1D(Form("htmp%.2d%.6d", drchit1->GetCh(), evt),Form("ch %.2d",drchit1->GetCh()+1),127,0.,127.);
          	}

          	htmpM[drchit1->GetCh()]->Reset(); //clearing histogram
          	for(int i1=0;i1<drchit1->GetSize();i1++){ //for each element of waveform, drchit1->GetSize() = # of elements (bins?) in waveform ~127
          	  htmpM[drchit1->GetCh()]->Fill(i1,(drchit1->GetW().At(i1))); //Fill histogram with (x, y)
            }
          	c1.cd(ch+1); //activate next element of c1.Divide matrix
          	MaxY[drchit1->GetCh()] = TMath::Max((Double_t)MaxY[drchit1->GetCh()],(Double_t)htmpM[drchit1->GetCh()]->GetBinContent(htmpM[drchit1->GetCh()]->GetMaximumBin())+50.);
          	MinY[drchit1->GetCh()] = TMath::Min((Double_t)MinY[drchit1->GetCh()],(Double_t)htmpM[drchit1->GetCh()]->GetBinContent(htmpM[drchit1->GetCh()]->GetMinimumBin())-50.);

            //for extra scale-setting histogram
            htmpM[32] = (TH1D *) gROOT->FindObject(Form("htmp%.2d%.6d", ChToArr[drchit1->GetCh()], start));
          	//         htmpM[28]->GetYaxis()->SetRangeUser(MinY[drchit1->GetCh()],MaxY[drchit1->GetCh()]);
          	htmpM[32]->GetYaxis()->SetRangeUser(hmin,hmax);
          	//       cout << Form("#%d Range %f \t %f \n", drchit1->GetCh(), MinY[drchit1->GetCh()],MaxY[drchit1->GetCh()]);
          	htmpM[drchit1->GetCh()]->Draw((evt==start)?"hist":"hist same"); //if evt==start, "hist", else "hist same" (draw in same canvas)

                      c1.Modified();
                      c1.Update();


                  //************************* MB2 *************************//
                  T2->GetEntry(evt);
                  TClonesArray *DRCHits2=event2->GetDRCHits();
                  if(evt<nEntriesMB2){
                    for(Int_t ch=0; ch<16;ch++){
                      DRCHit *drchit2 = (DRCHit*)DRCHits2->At(ch);
                      if(!drchit2) continue;
              	//       if(evt==0)
              	//      cout << Form("ch: %d %d\n", ChToArr[ch], drchit2->GetCh());
                      if(persist){
              	  htmpM[drchit2->GetCh()] = new TH1D(Form("htmp%.2d%.6d", drchit2->GetCh(), evt),Form("ch %.2d",drchit2->GetCh()+1),127,0.,127.);
                      }


                      htmpM[drchit2->GetCh()]->Reset();
                      for(int i1=0;i1<drchit2->GetSize();i1++){
              	  htmpM[drchit2->GetCh()]->Fill(i1,(drchit2->GetW().At(i1)));
                      }
                      c2.cd(ch+1);
              	MaxY[drchit2->GetCh()] = TMath::Max((Double_t)MaxY[drchit2->GetCh()],(Double_t)htmpM[drchit2->GetCh()]->GetBinContent(htmpM[drchit2->GetCh()]->GetMaximumBin())+50.);
              	MinY[drchit2->GetCh()] = TMath::Min((Double_t)MinY[drchit2->GetCh()],(Double_t)htmpM[drchit2->GetCh()]->GetBinContent(htmpM[drchit2->GetCh()]->GetMinimumBin())-50.);
              	htmpM[32] = (TH1D *) gROOT->FindObject(Form("htmp%.2d%.6d", ChToArr[drchit2->GetCh()], start));
              	//         htmpM[28]->GetYaxis()->SetRangeUser(MinY[drchit2->GetCh()],MaxY[drchit2->GetCh()]);
              	htmpM[32]->GetYaxis()->SetRangeUser(hmin,hmax);
              	//       cout << Form("#%d Range %f \t %f \n", drchit2->GetCh(), MinY[drchit2->GetCh()],MaxY[drchit2->GetCh()]);
              	htmpM[drchit2->GetCh()]->Draw((evt==start)?"hist":"hist same");
                    }
                  }

                  c2.Modified();
                  c2.Update();

                  //gSystem->ProcessEvents();

                  if(!persist){
                    TDirectory *curdir = gDirectory;
                    curdir->cd();
                  }

                }
                //c1.SaveAs(Form("./wfoutput/RUN_%d_MB1_WaveformPersistentFEB.png",RunNo));
                //c2.SaveAs(Form("./wfoutput/RUN_%d_MB2_WaveformPersistentFEB.png",RunNo));

                if(persist){
                  TDirectory *curdir = gDirectory;
                  curdir->cd();
                }

                gSystem->cd(sCurDir.Data());



        }
      }
    }
          if (cherrypick){   //for extracting channels of interest
            //************************* MB1 *************************//
              TClonesArray *DRCHits1=event1->GetDRCHits(); //get DRChits

              for(Int_t ch=0; ch<16;ch++){ //for each channel
                DRCHit *drchit = 0x0; //initialize

                for(int evt=start; evt<start+nevt; evt++){ // for each event
                  T1->GetEntry(evt); // get event
                  drchit = (DRCHit*)DRCHits1->At(ch);
                  if(!drchit) continue;
                if(persist){ //if you want waveform persistent/superimposed in same canvas
                   htmpM[((DRCHit*)DRCHits1->At(ch))->GetCh()] = new TH1D(Form("htmp%.2d%.6d", ((DRCHit*)DRCHits1->At(ch))->GetCh(), evt),Form("ch %.2d",((DRCHit*)DRCHits1->At(ch))->GetCh()+1),127,0.,127.);
                 }

                htmpM[ch]->Reset(); //clearing histogram
               	for(int i1=0;i1<drchit->GetSize();i1++){ //for each element of waveform, drchit1->GetSize() = # of elements (bins?) in waveform ~127
               	  htmpM[((DRCHit*)DRCHits1->At(ch))->GetCh()]->Fill(i1,((DRCHit*)DRCHits1->At(ch))->GetW().At(i1)); //Fill histogram with (x, y)
                 }
                 Double_t event_max = (Double_t)htmpM[drchit->GetCh()]->GetMaximum();

                 if (event_max > ch_max) {
                   ch_max = event_max; //stores as channel maximum
                 }
               }
               if (ch_max > threshold){
                 ch_list.push_back(ch+1); //add to channel list
               }
                ch_max = 0; //reset
              }

              //************************* MB2 *************************//
                TClonesArray *DRCHits2=event2->GetDRCHits(); //get DRChits

                for(Int_t ch=0; ch<16;ch++){ //for each channel
                  DRCHit *drchit = 0x0; //initialize

                  for(int evt=start; evt<start+nevt; evt++){ // for each event
                    T2->GetEntry(evt); // get event
                    drchit = (DRCHit*)DRCHits2->At(ch);
                    if(!drchit) continue;
                  if(persist){ //if you want waveform persistent/superimposed in same canvas
                     htmpM[((DRCHit*)DRCHits2->At(ch))->GetCh()] = new TH1D(Form("htmp%.2d%.6d", ((DRCHit*)DRCHits2->At(ch))->GetCh(), evt),Form("ch %.2d",((DRCHit*)DRCHits2->At(ch))->GetCh()+1),127,0.,127.);
                   }

                  htmpM[ch]->Reset(); //clearing histogram
                  for(int i1=0;i1<drchit->GetSize();i1++){ //for each element of waveform, drchit1->GetSize() = # of elements (bins?) in waveform ~127
                    htmpM[((DRCHit*)DRCHits2->At(ch))->GetCh()]->Fill(i1,((DRCHit*)DRCHits2->At(ch))->GetW().At(i1)); //Fill histogram with (x, y)
                   }
                   Double_t event_max = (Double_t)htmpM[drchit->GetCh()]->GetMaximum();

                   if (event_max > ch_max) {
                     ch_max = event_max; //stores as channel maximum
                   }
                 }
                 if (ch_max > threshold){
                   ch_list.push_back(ch+17); //add to channel list
                 }
                  ch_max = 0; //reset
                }
        }
  }

  if (nEntriesMB3 && nEntriesMB4) { //if the data has entries for MB3 and MB4
    // code for MB1/2/3/4
    cout << "it's got 4 motherboards" << endl;
    Int_t MB12max = TMath::Max(nEntriesMB1, nEntriesMB2);
    Int_t MB34max = TMath::Max(nEntriesMB3, nEntriesMB4);

    Int_t nEntries = TMath::Max(MB12max, MB34max);

    TCanvas c1(Form("WaveformPersistentFEB1_RUN_%d",RunNo),Form("Waveform %s",fname.Data()),1500,900);
    c1.Divide(4,4);
    TCanvas c2(Form("WaveformPersistentFEB2_RUN_%d",RunNo),Form("Waveform %s",fname.Data()),1500,900);
    c2.Divide(4,4);
    TCanvas c3(Form("WaveformPersistentFEB3_RUN_%d",RunNo),Form("Waveform %s",fname.Data()),1500,900);
    c3.Divide(4,4);
    TCanvas c4(Form("WaveformPersistentFEB4_RUN_%d",RunNo),Form("Waveform %s",fname.Data()),1500,900);
    c4.Divide(4,4);

    Int_t ChToArr[64]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63};

    TH1D *htmpM[65];
    for(Int_t ch=0; ch<65;ch++)
      htmpM[ch] = new TH1D(Form("htmp%.2d%.6d", ch, start),Form("ch %.2d",ch+1),127,0.,127.);

    gStyle->SetTitleFontSize(0.1);

    Float_t MinY[64];
    Float_t MaxY[64];


    for(Int_t ch=0; ch<64;ch++){
      MinY[ch]=1.e5;
      MaxY[ch]=-1.;
    }


    if (plot){ //for plotting capabilities
      for(int evt=start; evt<start+nevt; evt++){
        if(evt>nEntries){
          cout << Form("event %d not available, there are only %d events in RUN %d\n", evt, nEntries, RunNo);
          continue;
        }
        if(!(evt%ncheck)) cout <<Form("evt: %d\n", evt);


        //************************* MB1 *************************//
        T1->GetEntry(evt);
        TClonesArray *DRCHits1=event1->GetDRCHits();
        if(evt<nEntriesMB1){
          for(Int_t ch=0; ch<16;ch++){
    	DRCHit *drchit1 = (DRCHit*)DRCHits1->At(ch);
    	if(!drchit1) continue;
            if(evt==0)
                cout << Form("ch: %d %d\n", ChToArr[ch], drchit1->GetCh());
    	if(persist){
    	  htmpM[drchit1->GetCh()] = new TH1D(Form("htmp%.2d%.6d", drchit1->GetCh(), evt),Form("ch %.2d",drchit1->GetCh()+1),127,0.,127.);
    	}


    	htmpM[drchit1->GetCh()]->Reset();
    	for(int i1=0;i1<drchit1->GetSize();i1++){
    	  htmpM[drchit1->GetCh()]->Fill(i1,(drchit1->GetW().At(i1)));
    	}
    	c1.cd(ch+1);
    	MaxY[drchit1->GetCh()] = TMath::Max((Double_t)MaxY[drchit1->GetCh()],(Double_t)htmpM[drchit1->GetCh()]->GetBinContent(htmpM[drchit1->GetCh()]->GetMaximumBin())+50.);
    	MinY[drchit1->GetCh()] = TMath::Min((Double_t)MinY[drchit1->GetCh()],(Double_t)htmpM[drchit1->GetCh()]->GetBinContent(htmpM[drchit1->GetCh()]->GetMinimumBin())-50.);
    	htmpM[32] = (TH1D *) gROOT->FindObject(Form("htmp%.2d%.6d", ChToArr[drchit1->GetCh()], start));
    	//         htmpM[28]->GetYaxis()->SetRangeUser(MinY[drchit1->GetCh()],MaxY[drchit1->GetCh()]);
    	htmpM[32]->GetYaxis()->SetRangeUser(hmin,hmax);
    	//       cout << Form("#%d Range %f \t %f \n", drchit1->GetCh(), MinY[drchit1->GetCh()],MaxY[drchit1->GetCh()]);
    	htmpM[drchit1->GetCh()]->Draw((evt==start)?"hist":"hist same");
          }
        }


        c1.Modified();
        c1.Update();


        //************************* MB2 *************************//
        T2->GetEntry(evt);
        TClonesArray *DRCHits2=event2->GetDRCHits();
        if(evt<nEntriesMB2){
          for(Int_t ch=0; ch<16;ch++){
            DRCHit *drchit2 = (DRCHit*)DRCHits2->At(ch);
            if(!drchit2) continue;
    	//       if(evt==0)
    	//      cout << Form("ch: %d %d\n", ChToArr[ch], drchit2->GetCh());
            if(persist){
    	  htmpM[drchit2->GetCh()] = new TH1D(Form("htmp%.2d%.6d", drchit2->GetCh(), evt),Form("ch %.2d",drchit2->GetCh()+1),127,0.,127.);
            }


            htmpM[drchit2->GetCh()]->Reset();
            for(int i1=0;i1<drchit2->GetSize();i1++){
    	  htmpM[drchit2->GetCh()]->Fill(i1,(drchit2->GetW().At(i1)));
            }
            c2.cd(ch+1);
    	MaxY[drchit2->GetCh()] = TMath::Max((Double_t)MaxY[drchit2->GetCh()],(Double_t)htmpM[drchit2->GetCh()]->GetBinContent(htmpM[drchit2->GetCh()]->GetMaximumBin())+50.);
    	MinY[drchit2->GetCh()] = TMath::Min((Double_t)MinY[drchit2->GetCh()],(Double_t)htmpM[drchit2->GetCh()]->GetBinContent(htmpM[drchit2->GetCh()]->GetMinimumBin())-50.);
    	htmpM[32] = (TH1D *) gROOT->FindObject(Form("htmp%.2d%.6d", ChToArr[drchit2->GetCh()], start));
    	//         htmpM[28]->GetYaxis()->SetRangeUser(MinY[drchit2->GetCh()],MaxY[drchit2->GetCh()]);
    	htmpM[32]->GetYaxis()->SetRangeUser(hmin,hmax);
    	//       cout << Form("#%d Range %f \t %f \n", drchit2->GetCh(), MinY[drchit2->GetCh()],MaxY[drchit2->GetCh()]);
    	htmpM[drchit2->GetCh()]->Draw((evt==start)?"hist":"hist same");
          }
        }


        c2.Modified();
        c2.Update();

        //************************* MB3 *************************//
        T3->GetEntry(evt);
        TClonesArray *DRCHits3=event3->GetDRCHits();
        if(evt<nEntriesMB3){
          for(Int_t ch=0; ch<16;ch++){
        DRCHit *drchit3 = (DRCHit*)DRCHits3->At(ch);
        if(!drchit3) continue;
            if(evt==0)
                cout << Form("ch: %d %d\n", ChToArr[ch], drchit3->GetCh());
        if(persist){
        htmpM[drchit3->GetCh()] = new TH1D(Form("htmp%.2d%.6d", drchit3->GetCh(), evt),Form("ch %.2d",drchit3->GetCh()+1),127,0.,127.);
        }

        htmpM[drchit3->GetCh()]->Reset();
        for(int i1=0;i1<drchit3->GetSize();i1++){
          htmpM[drchit3->GetCh()]->Fill(i1,(drchit3->GetW().At(i1)));
        }
        c3.cd(ch+1);
        MaxY[drchit3->GetCh()] = TMath::Max((Double_t)MaxY[drchit3->GetCh()],(Double_t)htmpM[drchit3->GetCh()]->GetBinContent(htmpM[drchit3->GetCh()]->GetMaximumBin())+50.);
        MinY[drchit3->GetCh()] = TMath::Min((Double_t)MinY[drchit3->GetCh()],(Double_t)htmpM[drchit3->GetCh()]->GetBinContent(htmpM[drchit3->GetCh()]->GetMinimumBin())-50.);
        htmpM[64] = (TH1D *) gROOT->FindObject(Form("htmp%.2d%.6d", ChToArr[drchit3->GetCh()], start));
        //         htmpM[28]->GetYaxis()->SetRangeUser(MinY[drchit1->GetCh()],MaxY[drchit1->GetCh()]);
        htmpM[64]->GetYaxis()->SetRangeUser(hmin,hmax);
        //       cout << Form("#%d Range %f \t %f \n", drchit1->GetCh(), MinY[drchit1->GetCh()],MaxY[drchit1->GetCh()]);
        htmpM[drchit3->GetCh()]->Draw((evt==start)?"hist":"hist same");
          }
        }


        c3.Modified();
        c3.Update();


        //************************* MB4 *************************//
        T4->GetEntry(evt);
        TClonesArray *DRCHits4=event4->GetDRCHits();
        if(evt<nEntriesMB4){
          for(Int_t ch=0; ch<16;ch++){
            DRCHit *drchit4 = (DRCHit*)DRCHits4->At(ch);
            if(!drchit4) continue;
        //       if(evt==0)
        //      cout << Form("ch: %d %d\n", ChToArr[ch], drchit2->GetCh());
            if(persist){
        htmpM[drchit4->GetCh()] = new TH1D(Form("htmp%.2d%.6d", drchit4->GetCh(), evt),Form("ch %.2d",drchit4->GetCh()+1),127,0.,127.);
            }

            htmpM[drchit4->GetCh()]->Reset();
            for(int i1=0;i1<drchit4->GetSize();i1++){
        htmpM[drchit4->GetCh()]->Fill(i1,(drchit4->GetW().At(i1)));
            }
            c4.cd(ch+1);
        MaxY[drchit4->GetCh()] = TMath::Max((Double_t)MaxY[drchit4->GetCh()],(Double_t)htmpM[drchit4->GetCh()]->GetBinContent(htmpM[drchit4->GetCh()]->GetMaximumBin())+50.);
        MinY[drchit4->GetCh()] = TMath::Min((Double_t)MinY[drchit4->GetCh()],(Double_t)htmpM[drchit4->GetCh()]->GetBinContent(htmpM[drchit4->GetCh()]->GetMinimumBin())-50.);
        htmpM[64] = (TH1D *) gROOT->FindObject(Form("htmp%.2d%.6d", ChToArr[drchit4->GetCh()], start));
        //         htmpM[28]->GetYaxis()->SetRangeUser(MinY[drchit2->GetCh()],MaxY[drchit2->GetCh()]);
        htmpM[64]->GetYaxis()->SetRangeUser(hmin,hmax);
        //       cout << Form("#%d Range %f \t %f \n", drchit2->GetCh(), MinY[drchit2->GetCh()],MaxY[drchit2->GetCh()]);
        htmpM[drchit4->GetCh()]->Draw((evt==start)?"hist":"hist same");
          }
        }


        c4.Modified();
        c4.Update();

        gSystem->ProcessEvents();

        if(!persist){
          TDirectory *curdir = gDirectory;
          curdir->cd();
        }

      }
      // c1.SaveAs(Form("./wfoutput/RUN_%d_MB1_WaveformPersistentFEB.png",RunNo));
      // c2.SaveAs(Form("./wfoutput/RUN_%d_MB2_WaveformPersistentFEB.png",RunNo));
      // c3.SaveAs(Form("./wfoutput/RUN_%d_MB3_WaveformPersistentFEB.png",RunNo));
      // c4.SaveAs(Form("./wfoutput/RUN_%d_MB4_WaveformPersistentFEB.png",RunNo));

      if(persist){
        TDirectory *curdir = gDirectory;
        curdir->cd();
      }

      gSystem->cd(sCurDir.Data());
  }

  if(cherrypick){
    cout << "MB1 test" << endl;
    //************************* MB1 *************************//
      TClonesArray *DRCHits1=event1->GetDRCHits(); //get DRChits

      for(Int_t ch=0; ch<16;ch++){ //for each channel
        DRCHit *drchit = 0x0; //initialize

        for(int evt=start; evt<start+nevt; evt++){ // for each event
          T1->GetEntry(evt); // get event
          drchit = (DRCHit*)DRCHits1->At(ch);
          if(!drchit) continue;
        if(persist){ //if you want waveform persistent/superimposed in same canvas
           htmpM[((DRCHit*)DRCHits1->At(ch))->GetCh()] = new TH1D(Form("htmp%.2d%.6d", ((DRCHit*)DRCHits1->At(ch))->GetCh(), evt),Form("ch %.2d",((DRCHit*)DRCHits1->At(ch))->GetCh()+1),127,0.,127.);
         }

        htmpM[ch]->Reset(); //clearing histogram
        for(int i1=0;i1<drchit->GetSize();i1++){ //for each element of waveform, drchit1->GetSize() = # of elements (bins?) in waveform ~127
          htmpM[((DRCHit*)DRCHits1->At(ch))->GetCh()]->Fill(i1,((DRCHit*)DRCHits1->At(ch))->GetW().At(i1)); //Fill histogram with (x, y)
         }
         Double_t event_max = (Double_t)htmpM[drchit->GetCh()]->GetMaximum();

         if (event_max > ch_max) {
           ch_max = event_max; //stores as channel maximum
         }
       }
       if (ch_max > threshold){
         ch_list.push_back(ch+1); //add to channel list
       }
        ch_max = 0; //reset
      }
      cout << "MB2 test" << endl;
      //************************* MB2 *************************//
        TClonesArray *DRCHits2=event2->GetDRCHits(); //get DRChits

        for(Int_t ch=0; ch<16;ch++){ //for each channel
          DRCHit *drchit = 0x0; //initialize

          for(int evt=start; evt<start+nevt; evt++){ // for each event
            T2->GetEntry(evt); // get event
            drchit = (DRCHit*)DRCHits2->At(ch);
            if(!drchit) continue;
          if(persist){ //if you want waveform persistent/superimposed in same canvas
             htmpM[((DRCHit*)DRCHits2->At(ch))->GetCh()] = new TH1D(Form("htmp%.2d%.6d", ((DRCHit*)DRCHits2->At(ch))->GetCh(), evt),Form("ch %.2d",((DRCHit*)DRCHits2->At(ch))->GetCh()+1),127,0.,127.);
           }

          htmpM[ch]->Reset(); //clearing histogram
          for(int i1=0;i1<drchit->GetSize();i1++){ //for each element of waveform, drchit1->GetSize() = # of elements (bins?) in waveform ~127
            htmpM[((DRCHit*)DRCHits2->At(ch))->GetCh()]->Fill(i1,((DRCHit*)DRCHits2->At(ch))->GetW().At(i1)); //Fill histogram with (x, y)
           }
           Double_t event_max = (Double_t)htmpM[drchit->GetCh()]->GetMaximum();

           if (event_max > ch_max) {
             ch_max = event_max; //stores as channel maximum
           }
         }
         if (ch_max > threshold){
           ch_list.push_back(ch+17); //add to channel list
         }
          ch_max = 0; //reset
        }
        cout << "MB3 test" << endl;
        //************************* MB3 *************************//
          TClonesArray *DRCHits3=event3->GetDRCHits(); //get DRChits

          for(Int_t ch=0; ch<16;ch++){ //for each channel
            DRCHit *drchit = 0x0; //initialize

            for(int evt=start; evt<start+nevt; evt++){ // for each event
              T3->GetEntry(evt); // get event
              drchit = (DRCHit*)DRCHits3->At(ch);
              if(!drchit) continue;
            if(persist){ //if you want waveform persistent/superimposed in same canvas
               htmpM[((DRCHit*)DRCHits3->At(ch))->GetCh()] = new TH1D(Form("htmp%.2d%.6d", ((DRCHit*)DRCHits3->At(ch))->GetCh(), evt),Form("ch %.2d",((DRCHit*)DRCHits3->At(ch))->GetCh()+1),127,0.,127.);
             }

            htmpM[ch]->Reset(); //clearing histogram
            for(int i1=0;i1<drchit->GetSize();i1++){ //for each element of waveform, drchit1->GetSize() = # of elements (bins?) in waveform ~127
              htmpM[((DRCHit*)DRCHits3->At(ch))->GetCh()]->Fill(i1,((DRCHit*)DRCHits3->At(ch))->GetW().At(i1)); //Fill histogram with (x, y)
             }
             Double_t event_max = (Double_t)htmpM[drchit->GetCh()]->GetMaximum();

             if (event_max > ch_max) {
               ch_max = event_max; //stores as channel maximum
             }
           }
           if (ch_max > threshold){
             ch_list.push_back(ch+33); //add to channel list
           }
            ch_max = 0; //reset
          }
          cout << "MB4 test" << endl;
          //************************* MB4 *************************//
            TClonesArray *DRCHits4=event4->GetDRCHits(); //get DRChits

            for(Int_t ch=0; ch<16;ch++){ //for each channel
              DRCHit *drchit = 0x0; //initialize

              for(int evt=start; evt<start+nevt; evt++){ // for each event
                T4->GetEntry(evt); // get event
                drchit = (DRCHit*)DRCHits4->At(ch);
                if(!drchit) continue;
              if(persist){ //if you want waveform persistent/superimposed in same canvas
                 htmpM[((DRCHit*)DRCHits4->At(ch))->GetCh()] = new TH1D(Form("htmp%.2d%.6d", ((DRCHit*)DRCHits4->At(ch))->GetCh(), evt),Form("ch %.2d",((DRCHit*)DRCHits4->At(ch))->GetCh()+1),127,0.,127.);
               }

              htmpM[ch]->Reset(); //clearing histogram
              for(int i1=0;i1<drchit->GetSize();i1++){ //for each element of waveform, drchit1->GetSize() = # of elements (bins?) in waveform ~127
                htmpM[((DRCHit*)DRCHits4->At(ch))->GetCh()]->Fill(i1,((DRCHit*)DRCHits4->At(ch))->GetW().At(i1)); //Fill histogram with (x, y)
               }
               Double_t event_max = (Double_t)htmpM[drchit->GetCh()]->GetMaximum();

               if (event_max > ch_max) {
                 ch_max = event_max; //stores as channel maximum
               }
             }
             if (ch_max > threshold){
               ch_list.push_back(ch+49); //add to channel list
             }
              ch_max = 0; //reset
            }
  }
}

// opens log.txt for readout
TString datadir = "/Users/justinvega/Documents/Fermilab/src/adriano_analysis/wfoutput/";
ofstream log;
log.open(datadir + "ch_list.txt", ios::app);

log << Form("RUN %d", RunNo) << endl;
if (!(ch_list.empty())){
  for (int i=0; i < (int)ch_list.size(); i++){
    if (i == (int)ch_list.size()-1){
      log << Form("%d", (int)ch_list[i]) << endl;
}
    else {
      log << Form("%d, ", (int)ch_list[i]);
    }
}
}
else {
  log << "null" << endl;
}


log.close();
}
