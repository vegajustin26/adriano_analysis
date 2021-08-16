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
 ****/


{
  R__LOAD_LIBRARY(DRCEventV5_C);

  gROOT->Delete("fname");
  gROOT->Reset();

  TString datadir = "./ntup/";

  TSystemDirectory dir("data",datadir.Data());//TSystemDirectory dir("data","/data/ADRIANO2/TB_2019_dec/")
  int start=20;        //First event
  int nevt=200;       //how many events
  double hmin=-500.;   //min val in vertical axis
  double hmax=2500.; //max val in vertical axis
  bool persist = 1;   //0 = single waveforms; 1 = superimposed waveforms
  int RunNo =  1037; //1009; //1196;//810;   //RUN to select

  int ncheck = 10;   //print event number every ncheck events



  TString sCurDir(gSystem->pwd());


  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(0.7);
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
          cout << fname.Data() << endl;
          break;
        }
      }
    }
  }

  fname.Prepend(dir.GetTitle());
  cout << fname.Data() << endl;

  TFile f(fname.Data(),"","fname");

  TTree *T1 = (TTree*)f.Get("DRCMB1");
  DRCEvent *event1 = 0;
  T1->SetBranchAddress("DRCEvent", &event1);

  TTree *T2 = (TTree*)f.Get("DRCMB2");
  DRCEvent *event2 = 0;
  T2->SetBranchAddress("DRCEvent", &event2);

  Int_t nEntriesMB1 = T1->GetEntries();
  Int_t nEntriesMB2 = T2->GetEntries();

  cout << Form("MB1 has %d events\n", nEntriesMB1);
  cout << Form("MB2 has %d events\n", nEntriesMB2);

  Int_t nEntries = TMath::Max(nEntriesMB1, nEntriesMB2);

  TCanvas c1(Form("WaveformPersistentFEB1_RUN_%d",RunNo),Form("Waveform %s",fname.Data()),1500,900);
  c1.Divide(4,4);
  TCanvas c2(Form("WaveformPersistentFEB2_RUN_%d",RunNo),Form("Waveform %s",fname.Data()),1500,900);
  c2.Divide(4,4);

  Int_t ChToArr[32]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};

  TH1D *htmpM[33];
  for(Int_t ch=0; ch<33;ch++)
    htmpM[ch] = new TH1D(Form("htmp%.2d%.6d", ch, start),Form("ch %.2d",ch+1),127,0.,127.);

  gStyle->SetTitleFontSize(0.1);

  Float_t MinY[32];
  Float_t MaxY[32];


  for(Int_t ch=0; ch<32;ch++){
    MinY[ch]=1.e5;
    MaxY[ch]=-1.;
  }

  for(int evt=start; evt<start+nevt; evt++){
    if(evt>nEntries){
      cout << Form("event %d not available, there are only %d events in RUN %d\n", evt, nEntries, RunNo);
      continue;
    }
    if(!(evt%ncheck)) cout <<Form("evt: %d\n", evt);


    //************************* MB1 *************************//
    T1->GetEntry(evt); // for each event
    TClonesArray *DRCHits1=event1->GetDRCHits();
    if(evt<nEntriesMB1){
      for(Int_t ch=0; ch<16;ch++){
      	DRCHit *drchit1 = (DRCHit*)DRCHits1->At(ch);
      	if(!drchit1) continue; //if no event
              if(evt==0)
                  cout << Form("ch: %d %d\n", ChToArr[ch], drchit1->GetCh());
      	if(persist){ //if you want waveform persistent/superimposed in same canvas
      	  htmpM[drchit1->GetCh()] = new TH1D(Form("htmp%.2d%.6d", drchit1->GetCh(), evt),Form("ch %.2d",drchit1->GetCh()+1),127,0.,127.);
      	}


      	htmpM[drchit1->GetCh()]->Reset(); //clearing histogram
      	for(int i1=0;i1<drchit1->GetSize();i1++){ //one waveform, drchit1->GetSize() = # of elements in waveform ~127
      	  htmpM[drchit1->GetCh()]->Fill(i1,(drchit1->GetW().At(i1)));
      	}
      	c1.cd(ch+1);
      	MaxY[drchit1->GetCh()] = TMath::Max((Double_t)MaxY[drchit1->GetCh()],(Double_t)htmpM[drchit1->GetCh()]->GetBinContent(htmpM[drchit1->GetCh()]->GetMaximumBin())+50.);
      	MinY[drchit1->GetCh()] = TMath::Min((Double_t)MinY[drchit1->GetCh()],(Double_t)htmpM[drchit1->GetCh()]->GetBinContent(htmpM[drchit1->GetCh()]->GetMinimumBin())-50.);
      	htmpM[32] = (TH1D *) gROOT->FindObject(Form("htmp%.2d%.6d", ChToArr[drchit1->GetCh()], start)); //extra histogram
      	//         htmpM[28]->GetYaxis()->SetRangeUser(MinY[drchit1->GetCh()],MaxY[drchit1->GetCh()]);
      	htmpM[32]->GetYaxis()->SetRangeUser(hmin,hmax);
      	//       cout << Form("#%d Range %f \t %f \n", drchit1->GetCh(), MinY[drchit1->GetCh()],MaxY[drchit1->GetCh()]);
      	htmpM[drchit1->GetCh()]->Draw((evt==start)?"hist":"hist same"); //if evt==start, "hist", else "hist same" (draw in same canvas)
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


    if(!persist){
      TDirectory *curdir = gDirectory;
      curdir->cd();
    }

  }

c1.SaveAs(Form("./wfoutput/RUN_%d_MB1_WaveformPersistentFEB_scale.pdf",RunNo));
  if(persist){
    TDirectory *curdir = gDirectory;
    curdir->cd();
  }

  gSystem->cd(sCurDir.Data());
}
