/****
 * ROOT macro to process raw waveforms
 * Arguments are:
 * runset * to select the case in the switch construct
 * subdataset * is used in the output filename
 * dataset * is used in the output filename
 *
 * To run this macro to compile depends on DRCEventV5.C, functions.C and SignalIntegralFEBV2.C
 * so those macros need to be compiled
 * At ROOT prompt run
 * root [0] .L DRCEventV5.C++
 * root [1] .L functions.C++
 * root [2] .L SignalIntegralFEBV2.C++
 *
 * If those macros have already been compiled, it is enough to load SignalIntegralFEBV2_C.so library
 * root [0] gSystem->Load("SignalIntegralFEBV2_C.so");
 * then this macro can be execute
 * root [1] .x AnalysisIntegralFEBV2.C(100,1,2019)
 *
 * The switch 'runset' helps to organize the similar runs that uses same parameters, then loop over run indexes, few examples are below.
 * The actual work is done by SignalIntegralFEBV2 macro.
 *
Justin note:
  Here are some helpful pointers on what the arguments to the SignalIntegralScript mean:
  


 ****/

#include <TSystem.h>
#include <TStopwatch.h>
#include "TFile.h"
#include "TTree.h"
#include <Riostream.h>
#include <TSystemDirectory.h>
#include <TString.h>
#include <TError.h>
#include <TObjString.h>
#include "TSystem.h"
#include "TROOT.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include<algorithm>

using std::stoi;

void AnalysisIntegralFEBV2(int runset = 100, int subdataset = 12, int dataset = 2019){
  TStopwatch s;
  s.Start();

  gROOT->Reset();

//   if( gClassTable->GetID("DRCEvent") < 0 )
//   {
//     cout << Form("Loading DRCEvent library\n");
//
//     gSystem->Load("DRCEventV5_C.so");
//     gSystem->Load("functions_C.so");
//     gSystem->Load("SignalIntegralFEBV2_C.so");
//
//   }


  TString path = "../data/";
  //  Int_t angle = 0;

  switch (runset){

    case 100: {

      path = "./ntup/"; //look for root files here
      subdataset = 1;
      const int Nch=32; //How many channels have information you are looking at in them
      //int chlist[]={4, 5, 6, 7, 16, 17, 18, 19}; //channel list of interest
      int chlist[]={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32}; //channel list of interest
      Bool_t IsSiPMChList[]={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      //Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
      //Int_t XMaxList[]={20, 20, 20, 20, 20, 20};
      Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
      Bool_t debug = 0;
      Int_t wfsize = 127;
      Int_t FirstEvent=0;
      Int_t EventsToProcess=-1;
      Int_t nfit=1;
      Bool_t bShift=0; //baseline shift
      Bool_t bFast=1; //fitting/histo, fitting is 0
      Bool_t bUseMWPC=0;

      // grab list of RunNo from success txt file
      //TString datadir = "./ntup/";
      //std::vector<int> successlist;
      //int num;
      int idxRun;
      //int runslist[] = {1063, 1066, 1067, 1068, 1069, 1070, 1071, 1089, 1090, 1091, 1092, 1093};
      int runslist[] = {988, 991, 992, 993, 996, 997, 999, 1003, 1004, 1007, 1009, 1014, 1015};
      //int runslist[] = {1063};
      //
      // ifstream myfile(datadir + "success.txt");
      //
      // if (!myfile) {
      //   cout<<"Error opening output file"<<endl;
      // }
      // while(myfile >> num){ //adds RunNo from each line of file to runslist
      //   successlist.push_back(num);
      // }
      // myfile.close();

      //parse ch_list.txt for chlist
      // for(Int_t run=1190; run<1191; idxrun++){

  //       //finding index of successlist
  //       auto it = find(successlist.begin(), successlist.end(), run);
  //       if(it != successlist.end()){
  //         idxRun = it - successlist.begin();
  //         cout << idxRun << endl;
  //       }
  //       else {
  //         cout << Form("Run %d not in successlist!", run) << endl;
  //       }
  //
  //       TString datadir = "/Users/justinvega/Documents/Fermilab/src/adriano_analysis/wfoutput/";
  //
  //       ifstream ch_input(datadir + "ch_list.txt");
  //       string input;
  //       string ch;
  //       string str_check;
  //       std::vector<int> newchlist;
  //
  //       if (ch_input.is_open()){
  //         while(getline(ch_input, input, '\n')){
  //           TString sInput = input;
  //           if (sInput.Contains(Form("RUN %d", successlist[idxRun]))){ //finds RunNo from txtfile
  //             cout << Form("RUN %d", successlist[idxRun]) << endl;
  //             streampos run = ch_input.tellg(); // sets position
  //             getline(ch_input, str_check, '\n');
  //             TString null = str_check;
  //
  //             if (null.Contains("null")){ //checks if null
  //                 break;
  //             }
  //             else {
  //               ch_input.seekg(run);
  //               while(getline(ch_input, ch, ',')){
  //                 TString test = ch;
  //                 if (!test.Contains('\n')){
  //                   newchlist.push_back(stoi(ch)); //add to vector
  //               }
  //               else {
  //                 newchlist.push_back(stoi(ch.substr(0, ch.find("\n"))));
  //                 break;
  //                 }
  //               }
  //             }
  //           }
  //           else {
  //             continue;
  //           }
  //         }
  //       }
  //
  // ch_input.close();
  //
  // //setting variables
  // const int Nch = newchlist.size();
  // int chlist[Nch];
  // copy(newchlist.begin(), newchlist.end(), chlist);
  // Bool_t IsSiPMChList[Nch];
  // Int_t XMaxList[Nch];
  //
  // fill_n(IsSiPMChList, Nch, 1);
  // fill_n(XMaxList, Nch, 20);
for(Int_t idxrun=0; idxrun<sizeof(runslist)/sizeof(runslist[0]); idxrun++){
	Int_t iRun = runslist[idxrun];

	TStopwatch timer;
	timer.Start();

  //cout << "is this working?" << endl;
	SignalIntegralFEBV3(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);

	timer.Stop();
	timer.Print();
}
      }

    break;
  //
  //
  //   case 101: {
  //
  //     path = "/data/T-1015/TB112015/data/";
  //     subdataset = 1;
  //     const int Nch=32; //How many channels have information you are looking at in them
  //     int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31}; //channel list of interest
  //     Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  //     Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
  //     Bool_t debug = 0;
  //     Int_t wfsize = 127;
  //     Int_t FirstEvent=0;
  //     Int_t EventsToProcess=-1;
  //     Int_t nfit=1;
  //     Bool_t bShift=0;
  //     Bool_t bFast=1;
  //     Bool_t bUseMWPC=0;
  //
  //     for(Int_t idxRun=125; idxRun<=203; idxRun++){
	// if(idxRun==125) continue;
	// Int_t iRun = idxRun;
  //
	// TStopwatch timer;
	// timer.Start();
  //
	// SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);
  //
	// timer.Stop();
	// timer.Print();
  //
  //     }
  //   }
  //   break;
  //
  //
  //   case 102: {
  //
  //     path = "/data/T-1015/TB112015/data/";
  //     subdataset = 1;
  //     const int Nch=16; //How many channels have information you are looking at in them
  //     int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}; //channel list of interest
  //     Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  //     Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
  //     Bool_t debug = 0;
  //     Int_t wfsize = 127;
  //     Int_t FirstEvent=0;
  //     Int_t EventsToProcess=-1;
  //     Int_t nfit=1;
  //     Bool_t bShift=0;
  //     Bool_t bFast=1;
  //     Bool_t bUseMWPC=0;
  //
  //     for(Int_t idxRun=1; idxRun<=123; idxRun++){
	// //if(idxRun==125) continue;
	// Int_t iRun = idxRun;
  //
	// TStopwatch timer;
	// timer.Start();
  //
	// SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);
  //
	// timer.Stop();
	// timer.Print();
  //
  //     }
  //   }
  //   break;
  //
  //
  //   case 103: {
  //
  //     path = "/run/media/vito/dati/T-1015/TB112015/";
  //     subdataset = 1;
  //     const int Nch=19; //How many channels have information you are looking at in them
  //     int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}; //channel list of interest
  //     Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  //     Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 40, 40, 40};
  //     Bool_t debug = 0;
  //     Int_t wfsize = 127;
  //     Int_t FirstEvent=0;
  //     Int_t EventsToProcess=-1;
  //     Int_t nfit=1;
  //     Bool_t bShift=0;
  //     Bool_t bFast=1;
  //     Bool_t bUseMWPC=0;
  //
  //     for(Int_t idxRun=176; idxRun<=343; idxRun++){
	// if(idxRun<282) continue;
	// Int_t iRun = idxRun;
  //
	// TStopwatch timer;
	// timer.Start();
  //
	// SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);
  //
	// timer.Stop();
	// timer.Print();
  //
  //     }
  //   }
  //   break;
  //
  //
  //
  //   case 200: {
  //
  //     path = "/run/media/vito/dati/T-1015/TB112015/LED/";
  //     subdataset = 2;
  //     const int Nch=32; //How many channels have information you are looking at in them
  //     int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31}; //channel list of interest
  //     Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  //     Int_t XMaxList[]={40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
  //     Bool_t debug = 0;
  //     Int_t wfsize = 127;
  //     Int_t FirstEvent=0;
  //     Int_t EventsToProcess=-1;
  //     Int_t nfit=1;
  //     Bool_t bShift=0;
  //     Bool_t bFast=1;
  //     Bool_t bUseMWPC=0;
  //
  //     for(Int_t idxRun=183; idxRun<=345; idxRun++){
	// Int_t iRun = idxRun;
  //
	// TStopwatch timer;
	// timer.Start();
  //
	// SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);
  //
	// timer.Stop();
	// timer.Print();
  //
  //     }
  //   }
  //   break;
  //
  //
  //   case 300: {
  //
  //     path = "/run/media/vito/dati/T-1015/TB112015/";
  //     subdataset = 3;
  //     const int Nch=20; //How many channels have information you are looking at in them
  //     int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}; //channel list of interest
  //     Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  //     Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 40, 40, 40, 40};
  //     Bool_t debug = 0;
  //     Int_t wfsize = 127;
  //     Int_t FirstEvent=0;
  //     Int_t EventsToProcess=-1;
  //     Int_t nfit=1;
  //     Bool_t bShift=0;
  //     Bool_t bFast=1;
  //     Bool_t bUseMWPC=0;
  //
  //     for(Int_t idxRun=801; idxRun<=818; idxRun++){
	// Int_t iRun = idxRun;
  //
	// TStopwatch timer;
	// timer.Start();
  //
	// SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);
  //
	// timer.Stop();
	// timer.Print();
  //
  //     }
  //   }
  //   break;



    default:
      return;

  }
  s.Stop();
  s.Print();
}
