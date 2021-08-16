#ifndef ROOT_functions
#define ROOT_functions

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

#include "DRCEventV5.h"


void GetPedestalFEB(TArrayI w, double &ped, double &sigma, bool bShift=0)
{

  gROOT->Delete("hped");
  TH1D *hped = new TH1D("hped","hped",25,-100.,100.);

  for(int i1=0;i1<w.GetSize()*0+14;i1++){
//     if(w.At(i1)>0.)
      hped->Fill(w.At(i1)-100.*bShift);
  }

  sigma = hped->GetRMS();


  hped->BufferEmpty();
  TH1D *hped2 = new TH1D("hped2","hped2",25,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());

  for(int i1=0;i1<25;i1++){
    if(hped->GetBinContent(i1) > 4 /*&& hped->GetBinCenter(i1) < 200.-100.*bShift && hped->GetBinCenter(i1) > 50.-100.*bShift*/ ){
      hped2->SetBinContent(i1,hped->GetBinContent(i1));
    }
  }

  ped = hped2->GetMean();

  delete hped2;
  delete hped;

}

void GetPedestalPADE(TArrayI w, double &ped, double &sigma, bool bShift=0)
{

  gROOT->Delete("hped");
  TH1D *hped = new TH1D("hped","hped",25,0.,200.);

  for(int i1=0;i1<w.GetSize()*0+22;i1++){
    if(w.At(i1)>0.)
      hped->Fill(w.At(i1)-100.*bShift);
  }

  sigma = hped->GetRMS();


  hped->BufferEmpty();
  TH1D *hped2 = new TH1D("hped2","hped2",25,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());

  for(int i1=0;i1<25;i1++){
    if(hped->GetBinContent(i1) > 4 && hped->GetBinCenter(i1) < 200.-100.*bShift && hped->GetBinCenter(i1) > 50.-100.*bShift ){
      hped2->SetBinContent(i1,hped->GetBinContent(i1));
    }
  }


//   {
//     TH1D *hw = new TH1D("hw","hw",105,0.,105.);
//
//     for(int i1=0;i1<w.GetSize();i1++)
//       hw->Fill(i1,w.At(i1));
//
//
//     TCanvas c("c","c",1200,800);
//     c.Divide(2,2);
//     c.cd(1);
//     hw->Draw();
//     c.cd(2);
//     hped->Draw();
//     c.cd(3);
//     hped2->Draw();
//     c.Modified();
//     c.Update();
//     c.WaitPrimitive();
//     hw->Delete();
//   }

  ped = hped2->GetMean();
//   sigma = hped2->GetRMS();

  delete hped2;
  delete hped;

}

double GetPedestal(vector <int> w)
{

  TH1D *htmp = new TH1D("htmp","htmp",500,0.,8600.);

  for(int i1=0;i1<(int)w.size();i1++)
    htmp->Fill(w[i1]);

  TF1 *fitgaus = new TF1("fitgaus", "gaus",0.,10000.);

  htmp->Fit("fitgaus","NQ","",8000.,8600.);

  double ped = fitgaus->GetParameter(1);

  delete htmp;
  delete fitgaus;

  return ped;

}

void GetPedestal2(vector <int> w, double &ped, double &sigma)
{

  gROOT->Delete("hped");
  TH1D *hped = new TH1D("hped","hped",500,0.,0.);

  for(int i1=0;i1<(int)w.size();i1++){
    if(i1>50 && i1 <250) continue;
    hped->Fill(w[i1]);
  }

  hped->GetEntries();

  TF1 *fitgaus = new TF1("fitgaus", "gaus",0.,10000.);

  hped->Fit("fitgaus","NQ","",8000.,8600.);

  ped = fitgaus->GetParameter(1);
  sigma = fitgaus->GetParameter(2);

  delete hped;
  delete fitgaus;

//   return hped;

}

void GetPedestal2(TArrayI w, double &ped, double &sigma)
{

  gROOT->Delete("hped");
  TH1D *hped = new TH1D("hped","hped",500,0.,0.);

  for(int i1=0;i1<w.GetSize();i1++){
    if(i1>50 && i1 <250) continue;
    hped->Fill(w.At(i1));
  }

  hped->BufferEmpty();
  TF1 *fitgaus = new TF1("fitgaus", "gaus",0.,10000.);

  hped->Fit("fitgaus","NQ","",8000.,8600.);

  ped = fitgaus->GetParameter(1);
  sigma = fitgaus->GetParameter(2);

  delete hped;
  delete fitgaus;

//   return hped;

}

void GetPedestalSiPM(TArrayI *w, double &ped, double &sigma)
{

  gROOT->Delete("hped");
  TH1D *hped = new TH1D("hped","hped",100,0.,0.);

  for(int i1=0;i1<w->GetSize();i1++){
    hped->Fill(w->At(i1));
  }

  hped->BufferEmpty();
  TH1D *hped2 = new TH1D("hped2","hped2",100,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());

  for(int i1=0;i1<100;i1++){
    if(hped->GetBinContent(i1) > 7 && hped->GetBinCenter(i1) < 9000. && hped->GetBinCenter(i1) > 7000. ){
      hped2->SetBinContent(i1,hped->GetBinContent(i1));
    }
  }


  ped = hped2->GetMean();
  sigma = hped2->GetRMS();

  delete hped2;
  delete hped;

}


void GetPedestalSiPM(TArrayI w, double &ped, double &sigma, bool bShift=0)
{

  gROOT->Delete("hped");
  TH1D *hped = new TH1D("hped","hped",100,0.,0.);
//   TH1D *hw = new TH1D("hw","hw",500,0.,500.);

  for(int i1=0;i1<w.GetSize();i1++){
//     if(i1>50 && i1 <250) continue;
//     if(i1<80 || i1 > 400)
      hped->Fill(w.At(i1)-8192.*bShift);

//     hw->Fill(i1,w.At(i1));
  }

  hped->BufferEmpty();
  TH1D *hped2 = new TH1D("hped2","hped2",100,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());
//   TH1D *hped3 = new TH1D("hped3","hped3",100,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());
//   cout << Form("%f %f\n", hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());

  for(int i1=0;i1<100;i1++){
    if(hped->GetBinContent(i1) > 40 && hped->GetBinCenter(i1) < 9000.-8192.*bShift && hped->GetBinCenter(i1) > 7000.-8192.*bShift ){
      hped2->SetBinContent(i1,hped->GetBinContent(i1));
//       for(int i2=0;i2<hped->GetBinContent(i1);i2++){
//         hped3->Fill(hped->GetBinCenter(i1));
//       }
    }
  }

//   hped2->Sumw2();

//   double xMax = hped->GetBinCenter(hped->GetMaximumBin());
//   cout << Form("params: %f %f %f\n", hped->GetMaximum(), hped->GetMean(), hped->GetRMS());
//   hped->GetXaxis()->SetRangeUser(xMax-200.,xMax+200.);
//   cout << Form("params: %f %f %f\n", hped->GetMaximum(), hped->GetMean(), hped->GetRMS());
//   TF1 *fitgaus = new TF1("fitgaus", "gaus",0.,10000.);
//   TF1 *fitgaus2 = new TF1("fitgaus2", "gaus",0.,10000.);

//   hped->Fit("fitgaus","NQ","",8000.,8600.);
//   fitgaus->SetParameter(0,hped->GetMaximum());
//   fitgaus->SetParLimits(0,30.,1000.);
//   fitgaus->SetParameter(1,hped->GetMean());
//   fitgaus->SetParLimits(1,hped->GetMean()-3.*hped->GetRMS(),hped->GetMean()+3.*hped->GetRMS());
//   fitgaus->SetParameter(2,/*hped->GetRMS()*/1.5);
//   fitgaus->SetParLimits(2,.9,1.);

//   hped2->Fit("fitgaus2","","",hped2->GetMean()-3.*hped2->GetRMS(),hped2->GetMean()+3.*hped2->GetRMS());
//   hped3->Fit("gaus","","",hped3->GetMean()-3.*hped3->GetRMS(),hped3->GetMean()+3.*hped3->GetRMS());
//   cout << Form("params: %f %f %f\n", hped->GetMaximum(), hped->GetMean(), hped->GetRMS());
//   cout << Form("params: %f %f %f\n", fitgaus->GetParameter(0), fitgaus->GetParameter(1), fitgaus->GetParameter(2));
//   hped->Fit("fitgaus","","",hped->GetMean()-3.*hped->GetRMS(),hped->GetMean()+3.*hped->GetRMS());
//   hped->Fit("fitgaus","Q","",8000.,8300.);

//   ped = fitgaus->GetParameter(1);
//   sigma = fitgaus->GetParameter(2);

  ped = hped2->GetMean();
  sigma = hped2->GetRMS();

//   cout << Form("params: %f %f %f\n", fitgaus->GetParameter(0), fitgaus->GetParameter(1), fitgaus->GetParameter(2));

    //     gSystem->Sleep(1);
//   TCanvas c("c","c",1200,800);
//   c.Divide(2,2);
//   c.cd(1);
//   hped->Draw();
//   c.cd(2);
//   hw->Draw();
//   TLine loriz(0., fitgaus->GetParameter(1), 500., fitgaus->GetParameter(1));
//   loriz.SetLineColor(1);
//   loriz.SetLineWidth(4);
//   loriz.Draw("same");
//   TLine loriz1(0., hped2->GetMean(), 500., hped2->GetMean());
//   loriz1.SetLineColor(4);
//   loriz1.SetLineWidth(2);
//   loriz1.Draw("same");
//   TLine loriz2(0., fitgaus2->GetParameter(1), 500., fitgaus2->GetParameter(1));
//   loriz2.SetLineColor(3);
//   loriz2.SetLineWidth(1);
//   loriz2.Draw("same");
//   c.cd(3);
//   hped2->Draw();
//   c.cd(4);
//   hped3->Draw();

//   c.Modified();
//   c.Update();
//   c.WaitPrimitive();


//   delete hw;
  delete hped2;
  delete hped;
//   delete fitgaus;

//   return hped;

}

void GetPedestalSiPM2(TArrayI w, double &ped, double &sigma)
{

  gROOT->Delete("hped");
  TH1D *hped = new TH1D("hped","hped",100,0.,0.);

  for(int i1=0;i1<w.GetSize();i1++){
    if(i1<80 || i1 > 400)
      hped->Fill(w.At(i1)/4.);
  }

  hped->BufferEmpty();
  TH1D *hped2 = new TH1D("hped2","hped2",100,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());

  for(int i1=0;i1<100;i1++){
    if(hped->GetBinContent(i1) > 7 && hped->GetBinCenter(i1) < 9000./4. && hped->GetBinCenter(i1) > 7000./4. ){
      hped2->SetBinContent(i1,hped->GetBinContent(i1));
    }
  }


  ped = hped2->GetMean();
  sigma = hped2->GetRMS();


  delete hped2;
  delete hped;

}

void GetPedestalSiPMDebug(TArrayI w, double &ped, double &sigma)
{

  gROOT->Delete("hped");
  TH1D *hped = new TH1D("hped","hped",100,0.,0.);
  TH1D *hw = new TH1D("hw","hw",w.GetSize(),0.,w.GetSize());

  for(int i1=0;i1<w.GetSize();i1++){
//     if(i1>50 && i1 <250) continue;
//     if(i1<80 || i1 > 400)
    hped->Fill(w.At(i1));
    hw->Fill(i1,w.At(i1));
  }

  hped->BufferEmpty();
  TH1D *hped2 = new TH1D("hped2","hped2",100,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());
//   TH1D *hped3 = new TH1D("hped3","hped3",100,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());
  cout << Form("%f %f\n", hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());

  for(int i1=0;i1<100;i1++){
    if(hped->GetBinContent(i1) > 5 && hped->GetBinCenter(i1) < 9000. && hped->GetBinCenter(i1) > 7000. ){
      hped2->SetBinContent(i1,hped->GetBinContent(i1));
//       for(int i2=0;i2<hped->GetBinContent(i1);i2++){
//         hped3->Fill(hped->GetBinCenter(i1));
//       }
    }
  }

//   hped2->Sumw2();

//   double xMax = hped->GetBinCenter(hped->GetMaximumBin());
//   cout << Form("params: %f %f %f\n", hped->GetMaximum(), hped->GetMean(), hped->GetRMS());
//   hped->GetXaxis()->SetRangeUser(xMax-200.,xMax+200.);
//   cout << Form("params: %f %f %f\n", hped->GetMaximum(), hped->GetMean(), hped->GetRMS());
  TF1 *fitgaus = new TF1("fitgaus", "gaus",0.,10000.);
  TF1 *fitgaus2 = new TF1("fitgaus2", "gaus",0.,10000.);

//   hped->Fit("fitgaus","NQ","",8000.,8600.);
//   fitgaus->SetParameter(0,hped->GetMaximum());
//   fitgaus->SetParLimits(0,30.,1000.);
//   fitgaus->SetParameter(1,hped->GetMean());
//   fitgaus->SetParLimits(1,hped->GetMean()-3.*hped->GetRMS(),hped->GetMean()+3.*hped->GetRMS());
//   fitgaus->SetParameter(2,/*hped->GetRMS()*/1.5);
//   fitgaus->SetParLimits(2,.9,1.);

  hped2->Fit("fitgaus2","","",hped2->GetMean()-3.*hped2->GetRMS(),hped2->GetMean()+3.*hped2->GetRMS());
//   hped3->Fit("gaus","","",hped3->GetMean()-3.*hped3->GetRMS(),hped3->GetMean()+3.*hped3->GetRMS());
//   cout << Form("params: %f %f %f\n", hped->GetMaximum(), hped->GetMean(), hped->GetRMS());
//   cout << Form("params: %f %f %f\n", fitgaus->GetParameter(0), fitgaus->GetParameter(1), fitgaus->GetParameter(2));
  hped->Fit("fitgaus","","",hped->GetMean()-3.*hped->GetRMS(),hped->GetMean()+3.*hped->GetRMS());
//   hped->Fit("fitgaus","Q","",8000.,8300.);

  ped = fitgaus->GetParameter(1);
  sigma = fitgaus->GetParameter(2);
  cout << Form("params: %f %f %f\n", fitgaus->GetParameter(0), fitgaus->GetParameter(1), fitgaus->GetParameter(2));

    //     gSystem->Sleep(1);
  TCanvas c("c","c",1200,800);
  c.Divide(2,2);
  c.cd(1);
  hped->Draw();
  c.cd(2);
  hw->Draw();
  TLine loriz(0., fitgaus->GetParameter(1), w.GetSize(), fitgaus->GetParameter(1));
  loriz.SetLineColor(1);
  loriz.SetLineWidth(4);
  loriz.Draw("same");
  TLine loriz1(0., hped2->GetMean(), w.GetSize(), hped2->GetMean());
  loriz1.SetLineColor(4);
  loriz1.SetLineWidth(2);
  loriz1.Draw("same");
  TLine loriz2(0., fitgaus2->GetParameter(1), w.GetSize(), fitgaus2->GetParameter(1));
  loriz2.SetLineColor(3);
  loriz2.SetLineWidth(1);
  loriz2.Draw("same");
  c.cd(3);
  hped2->Draw();
//   c.cd(4);
//   hped3->Draw();

  c.Modified();
  c.Update();
  c.WaitPrimitive();


  delete hw;
  delete hped2;
  delete hped;
  delete fitgaus;

//   return hped;

}

void GetPedestalPMT(TArrayI *w, double &ped, double &sigma)
{

  gROOT->Delete("hped");
  gROOT->Delete("hped2");
  TH1D *hped = new TH1D("hped","hped",100,0.,0.);

  for(int i1=0;i1<w->GetSize();i1++){
    hped->Fill(w->At(i1));
  }

  hped->BufferEmpty();
  TH1D *hped2 = new TH1D("hped2","hped2",100,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());

  for(int i1=0;i1<100;i1++){
    if(hped->GetBinContent(i1) > 5 && hped->GetBinCenter(i1) < 9000. && hped->GetBinCenter(i1) > 7000. ){
      hped2->SetBinContent(i1,hped->GetBinContent(i1));
    }
  }

  ped = hped2->GetMean();
  sigma = hped2->GetRMS();

  delete hped2;
  delete hped;

}


void GetPedestalPMT(TArrayI w, double &ped, double &sigma, bool bShift=0)
{

  gROOT->Delete("hped");
  gROOT->Delete("hped2");
  TH1D *hped = new TH1D("hped","hped",100,0.,0.);
//   TH1D *hw = new TH1D("hw","hw",500,0.,500.);

  for(int i1=0;i1<w.GetSize();i1++){
//     if(i1<150)
    hped->Fill((w.At(i1)-8192.*bShift));

//     hw->Fill(i1,w.At(i1));
  }

  hped->BufferEmpty();
  TH1D *hped2 = new TH1D("hped2","hped2",100,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());
//   cout << Form("%f %f\n", hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());

  for(int i1=0;i1<100;i1++){
    if(hped->GetBinContent(i1) > 40 && hped->GetBinCenter(i1) < 9000.-8192.*bShift && hped->GetBinCenter(i1) > 7000.-8192.*bShift ){
      hped2->SetBinContent(i1,hped->GetBinContent(i1));
    }
  }

//   TF1 *fitgaus = new TF1("fitgaus", "gaus",0.,10000.);
//   TF1 *fitgaus2 = new TF1("fitgaus2", "gaus",0.,10000.);
//
//
//   hped2->Fit("fitgaus2","","",hped2->GetMean()-3.*hped2->GetRMS(),hped2->GetMean()+3.*hped2->GetRMS());
//
//   hped->Fit("fitgaus","","",hped->GetMean()-3.*hped->GetRMS(),hped->GetMean()+3.*hped->GetRMS());

  ped = hped2->GetMean();
  sigma = hped2->GetRMS();

//   ped = fitgaus->GetParameter(1);
//   sigma = fitgaus->GetParameter(2);
//   cout << Form("params: %f %f %f\n", fitgaus->GetParameter(0), fitgaus->GetParameter(1), fitgaus->GetParameter(2));

    //     gSystem->Sleep(1);
//   TCanvas c("c","c",1200,800);
//   c.Divide(2,2);
//   c.cd(1);
//   hped->Draw();
//   c.cd(2);
//   hw->Draw();
//   TLine loriz(0., fitgaus->GetParameter(1), 500., fitgaus->GetParameter(1));
//   loriz.SetLineColor(1);
//   loriz.SetLineWidth(4);
//   loriz.Draw("same");
//   TLine loriz1(0., hped2->GetMean(), 500., hped2->GetMean());
//   loriz1.SetLineColor(4);
//   loriz1.SetLineWidth(2);
//   loriz1.Draw("same");
//   TLine loriz2(0., fitgaus2->GetParameter(1), 500., fitgaus2->GetParameter(1));
//   loriz2.SetLineColor(3);
//   loriz2.SetLineWidth(1);
//   loriz2.Draw("same");
//   c.cd(3);
//   hped2->Draw();
//   c.cd(4);
//   hped3->Draw();

//   c.Modified();
//   c.Update();
//   c.WaitPrimitive();


//   delete hw;
  delete hped2;
  delete hped;
//   delete fitgaus;

//   return hped;

}

void GetPedestalPMT2(TArrayI w, double &ped, double &sigma)
{

  gROOT->Delete("hped");
  gROOT->Delete("hped2");
  TH1D *hped = new TH1D("hped","hped",100,0.,0.);

  for(int i1=0;i1<w.GetSize();i1++){
//     if(i1<150)
      hped->Fill(w.At(i1)/4.);
  }

  hped->BufferEmpty();
  TH1D *hped2 = new TH1D("hped2","hped2",100,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());
//   cout << Form("%f %f\n", hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());

  for(int i1=0;i1<100;i1++){
    if(hped->GetBinContent(i1) > 5 && hped->GetBinCenter(i1) < 9000./4. && hped->GetBinCenter(i1) > 7000./4. ){
      hped2->SetBinContent(i1,hped->GetBinContent(i1));
    }
  }


  ped = hped2->GetMean();
  sigma = hped2->GetRMS();


  delete hped2;
  delete hped;

}

void GetPedestalPMTDebud(TArrayI w, double &ped, double &sigma)
{

  gROOT->Delete("hped");
  gROOT->Delete("hped2");
  gROOT->Delete("hw");
  TH1D *hped = new TH1D("hped","hped",100,0.,0.);
  TH1D *hw = new TH1D("hw","hw",w.GetSize(),0.,w.GetSize());

  for(int i1=0;i1<w.GetSize();i1++){
//     if(i1<150)
      hped->Fill(w.At(i1));
    hw->Fill(i1,w.At(i1));
  }

  hped->BufferEmpty();
  TH1D *hped2 = new TH1D("hped2","hped2",100,hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());
  cout << Form("%f %f\n", hped->GetXaxis()->GetXmin(),hped->GetXaxis()->GetXmax());

  for(int i1=0;i1<100;i1++){
    if(hped->GetBinContent(i1) > 5 && hped->GetBinCenter(i1) < 9000. && hped->GetBinCenter(i1) > 7000. ){
      hped2->SetBinContent(i1,hped->GetBinContent(i1));
    }
  }

  TF1 *fitgaus = new TF1("fitgaus", "gaus",0.,10000.);
  TF1 *fitgaus2 = new TF1("fitgaus2", "gaus",0.,10000.);


  hped2->Fit("fitgaus2","","",hped2->GetMean()-3.*hped2->GetRMS(),hped2->GetMean()+3.*hped2->GetRMS());

  hped->Fit("fitgaus","","",hped->GetMean()-3.*hped->GetRMS(),hped->GetMean()+3.*hped->GetRMS());

  ped = hped2->GetMean();
  sigma = hped2->GetRMS();

//   ped = fitgaus->GetParameter(1);
//   sigma = fitgaus->GetParameter(2);
//   cout << Form("params: %f %f %f\n", fitgaus->GetParameter(0), fitgaus->GetParameter(1), fitgaus->GetParameter(2));

    //     gSystem->Sleep(1);
  TCanvas c("c","c",1200,800);
  c.Divide(2,2);
  c.cd(1);
  hped->Draw();
  c.cd(2);
  hw->Draw();
  TLine loriz(0., fitgaus->GetParameter(1), w.GetSize(), fitgaus->GetParameter(1));
  loriz.SetLineColor(1);
  loriz.SetLineWidth(4);
  loriz.Draw("same");
  TLine loriz1(0., hped2->GetMean(), w.GetSize(), hped2->GetMean());
  loriz1.SetLineColor(4);
  loriz1.SetLineWidth(2);
  loriz1.Draw("same");
  TLine loriz2(0., fitgaus2->GetParameter(1), w.GetSize(), fitgaus2->GetParameter(1));
  loriz2.SetLineColor(3);
  loriz2.SetLineWidth(1);
  loriz2.Draw("same");
  c.cd(3);
  hped2->Draw();

  c.Modified();
  c.Update();
  c.WaitPrimitive();


  delete hw;
  delete hped2;
  delete hped;
  delete fitgaus;
  delete fitgaus2;

//   return hped;

}


void GetPedestalPMT(TArrayI w, double &ped, double &sigma, int lBin, int hBin)
{
//[lBin, hBin] is the exclusion interval (where signal is expected)
  gROOT->Delete("hped");
  TH1D *hped = new TH1D("hped","hped",500,0.,0.);

  for(int i1=0;i1<w.GetSize();i1++){
    if(i1>lBin && i1 <hBin) continue;
    hped->Fill(w.At(i1));
  }

  hped->BufferEmpty();

  TF1 *fitgaus = new TF1("fitgaus", "gaus",0.,10000.);

  hped->Fit("fitgaus","NQ","",8150.,8250.);

  ped = fitgaus->GetParameter(1);
  sigma = fitgaus->GetParameter(2);

  delete hped;
  delete fitgaus;

//   return hped;

}


Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
  //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


      // MP shift correction
  mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}



TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
  //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

  Int_t i;
  Char_t FunName[100];

  sprintf(FunName,"Fitfcn_%s",his->GetName());

  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MP","Area","GSigma");

  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

  his->Fit(FunName,"QRB0");   // fit within specified range, use ParLimits, do not plot

  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf

  return (ffit);              // return fit function

}

Double_t langaufun2(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3] + (par[4] + par[5]*x[0]) );
}



TF1 *langaufit2(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi/*, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF*/, Int_t nloops=1, Option_t *opt1="QRB0WI", Option_t *opt2="QRB0WI")
{
  // Once again, here are the Landau * Gaussian parameters:
  //   par[0]=Width (scale) parameter of Landau density
  //   par[1]=Most Probable (MP, location) parameter of Landau density
  //   par[2]=Total area (integral -inf to inf, normalization constant)
  //   par[3]=Width (sigma) of convoluted Gaussian function
  //
  // Variables for langaufit call:
  //   his             histogram to fit
  //   fitrange[2]     lo and hi boundaries of fit range
  //   startvalues[6]  reasonable start values for the fit
  //   parlimitslo[6]  lower parameter limits
  //   parlimitshi[6]  upper parameter limits
  //   fitparams[6]    returns the final fit parameters
  //   fiterrors[6]    returns the final fit errors
  //   ChiSqr          returns the chi square
  //   NDF             returns ndf

  Int_t i;
  Char_t FunName[100];

  sprintf(FunName,"Fitfcn_%s",his->GetName());

  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName,langaufun2,fitrange[0],fitrange[1],6);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MP","Area","GSigma", "p0", "p1");

  for (i=0; i<6; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

  //   his->Fit(FunName,"QRB0WI");   // fit within specified range, use ParLimits, do not plot

  for(Int_t ii=0; ii<nloops-1; ii++)
    his->Fit(FunName,opt1);


  his->Fit(FunName,opt2);

//   ffit->GetParameters(fitparams);    // obtain fit parameters
//   for (i=0; i<6; i++) {
//     fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
//   }
//   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
//   NDF[0] = ffit->GetNDF();           // obtain ndf

  return (ffit);              // return fit function

}

Double_t langaufun3(Double_t *x, Double_t *par) {

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;


  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);

    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3] + (par[4]*TMath::Exp(par[5]*x[0]) ) );
}



TF1 *langaufit3(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi/*, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF*/, Int_t nloops=1, Option_t *opt1="QRB0WI", Option_t *opt2="QRB0WI")
{
  // Once again, here are the Landau * Gaussian parameters:
  //   par[0]=Width (scale) parameter of Landau density
  //   par[1]=Most Probable (MP, location) parameter of Landau density
  //   par[2]=Total area (integral -inf to inf, normalization constant)
  //   par[3]=Width (sigma) of convoluted Gaussian function
  //
  // Variables for langaufit call:
  //   his             histogram to fit
  //   fitrange[2]     lo and hi boundaries of fit range
  //   startvalues[6]  reasonable start values for the fit
  //   parlimitslo[6]  lower parameter limits
  //   parlimitshi[6]  upper parameter limits
  //   fitparams[6]    returns the final fit parameters
  //   fiterrors[6]    returns the final fit errors
  //   ChiSqr          returns the chi square
  //   NDF             returns ndf

  Int_t i;
  Char_t FunName[100];

  sprintf(FunName,"Fitfcn_%s",his->GetName());

  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName,langaufun3,fitrange[0],fitrange[1],6);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MP","Area","GSigma", "Const", "Slope");

  for (i=0; i<6; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

  //   his->Fit(FunName,"QRB0WI");   // fit within specified range, use ParLimits, do not plot

  for(Int_t ii=0; ii<nloops-1; ii++)
    his->Fit(FunName,opt1);


  his->Fit(FunName,opt2);

  //   ffit->GetParameters(fitparams);    // obtain fit parameters
  //   for (i=0; i<6; i++) {
  //     fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  //   }
  //   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  //   NDF[0] = ffit->GetNDF();           // obtain ndf

  return (ffit);              // return fit function

}


Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   // Seaches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
  //
   // The search is probably not very efficient, but it's a first try.

  Double_t p,x,fy,fxr,fxl;
  Double_t step;
  Double_t l,lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;


   // Search for maximum

  p = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l    = -1.0;


  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;

    lold = l;
    x = p + step;
    l = langaufun(&x,params);

    if (l < lold)
      step = -step/10;

    p += step;
  }

  if (i == MAXCALLS)
    return (-1);

  maxx = x;

  fy = l/2;


   // Search for right x location of fy

  p = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;


  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;

    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);

    if (l > lold)
      step = -step/10;

    p += step;
  }

  if (i == MAXCALLS)
    return (-2);

  fxr = x;


   // Search for left x location of fy

  p = maxx - 0.5 * params[0];
  step = -params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;

  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;

    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);

    if (l > lold)
      step = -step/10;

    p += step;
  }

  if (i == MAXCALLS)
    return (-3);


  fxl = x;

  FWHM = fxr - fxl;
  return (0);
}

Double_t fitfn(Double_t *x, Double_t *par){
  Double_t fn = -999.;
  if (x[0] < par[4]){
    fn = par[0];
  }else{
    fn = par[1]*(TMath::Exp(par[2]*(x[0]-par[4])))*(1-TMath::Exp(par[3]*(x[0]-par[4])))+par[0];
  }
  return fn;
}

Double_t fitfn2(Double_t *x, Double_t *par){

  Double_t fn = -999.;

  if (x[0]< par[0]){
    fn = par[1];
  }else{
    fn = (par[2]*TMath::Exp(par[3]*(x[0]-par[0]))+par[4]*TMath::Exp(par[5]*(x[0]-par[0])))*(1-TMath::Exp(par[6]*(x[0]-par[0])))+par[7];
  }
  return fn;
}

#endif
