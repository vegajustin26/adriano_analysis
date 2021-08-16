/****
 * ROOT macro to process raw waveforms
 * Arguments are:
 * runset * to select the case in the swithc construct
 * subdataset * is used in the output filename
 * dataset * is used in the output filename
 *
 * To run this macro it is needed to compile DRCEventV5.C, functions.C and SignalIntegralFEBV2.C first
 * At ROOT prompt run
 * root [0] .L DRCEventV5.C++
 * root [1] .L functions.C++
 * root [2] .L SignalIntegralFEBV2.C++
 * root [3] .L AnalysisIntegralFEBV2.C++
 *
 * If those macros have already been compiled, it is enough to load its library
 *
 * root [0] gSystem->Load("SignalIntegralFEBV2_C.so");
 * root [1] gSystem->Load("AnalysisIntegralFEBV2_C.so");
 *
 * The switch 'runset' helps to organize the similar runs that uses same parameters, then loop over run indexes, few examples are below.
 * The actual work is done by SignalIntegralFEBV2 macro.
 *
 * This macro can be executes as:
 * AnalysisIntegralFEBV2(100,1,2019)
 *
 ****/

#include <TSystem.h>
#include <TStopwatch.h>

#include "functions.C"
#include "SignalIntegralFEBV2.C"


void AnalysisIntegralFEBV2(int runset = 100, int subdataset = 1, int dataset = 2019){
  TStopwatch s;
  s.Start();

  gROOT->Reset();

  if( gClassTable->GetID("DRCEvent") < 0 )
  {
    cout << Form("Loading DRCEvent library\n");

    gSystem->Load("DRCEventV5_C.so");
    gSystem->Load("functions_C.so");
    gSystem->Load("SignalIntegralFEBV2_C.so");

  }


  TString path = "../data/";
  //  Int_t angle = 0;

  switch (runset){

    case 100: {

      path = "/data/ADRIANO2/TB_2019_dec/";
      subdataset = 1;
      const int Nch=32; //How many channels have information you are looking at in them
      int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31}; //channel list of interest
      Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
      Bool_t debug = 0;
      Int_t wfsize = 127;
      Int_t FirstEvent=0;
      Int_t EventsToProcess=-1;
      Int_t nfit=1;
      Bool_t bShift=0;
      Bool_t bFast=1;
      Bool_t bUseMWPC=0;

      for(Int_t idxRun=1190; idxRun<=1190; idxRun++){
	// 	if(idxRun==5) continue;
	Int_t iRun = idxRun;

	TStopwatch timer;
	timer.Start();

	SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);

	timer.Stop();
	timer.Print();

      }
    }
    break;


    case 101: {

      path = "/data/T-1015/TB112015/data/";
      subdataset = 1;
      const int Nch=32; //How many channels have information you are looking at in them
      int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31}; //channel list of interest
      Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
      Bool_t debug = 0;
      Int_t wfsize = 127;
      Int_t FirstEvent=0;
      Int_t EventsToProcess=-1;
      Int_t nfit=1;
      Bool_t bShift=0;
      Bool_t bFast=1;
      Bool_t bUseMWPC=0;

      for(Int_t idxRun=125; idxRun<=203; idxRun++){
	if(idxRun==125) continue;
	Int_t iRun = idxRun;

	TStopwatch timer;
	timer.Start();

	SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);

	timer.Stop();
	timer.Print();

      }
    }
    break;


    case 102: {

      path = "/data/T-1015/TB112015/data/";
      subdataset = 1;
      const int Nch=16; //How many channels have information you are looking at in them
      int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}; //channel list of interest
      Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
      Bool_t debug = 0;
      Int_t wfsize = 127;
      Int_t FirstEvent=0;
      Int_t EventsToProcess=-1;
      Int_t nfit=1;
      Bool_t bShift=0;
      Bool_t bFast=1;
      Bool_t bUseMWPC=0;

      for(Int_t idxRun=1; idxRun<=123; idxRun++){
	//if(idxRun==125) continue;
	Int_t iRun = idxRun;

	TStopwatch timer;
	timer.Start();

	SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);

	timer.Stop();
	timer.Print();

      }
    }
    break;


    case 103: {

      path = "/run/media/vito/dati/T-1015/TB112015/";
      subdataset = 1;
      const int Nch=19; //How many channels have information you are looking at in them
      int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}; //channel list of interest
      Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 40, 40, 40};
      Bool_t debug = 0;
      Int_t wfsize = 127;
      Int_t FirstEvent=0;
      Int_t EventsToProcess=-1;
      Int_t nfit=1;
      Bool_t bShift=0;
      Bool_t bFast=1;
      Bool_t bUseMWPC=0;

      for(Int_t idxRun=176; idxRun<=343; idxRun++){
	if(idxRun<282) continue;
	Int_t iRun = idxRun;

	TStopwatch timer;
	timer.Start();

	SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);

	timer.Stop();
	timer.Print();

      }
    }
    break;



    case 200: {

      path = "/run/media/vito/dati/T-1015/TB112015/LED/";
      subdataset = 2;
      const int Nch=32; //How many channels have information you are looking at in them
      int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31}; //channel list of interest
      Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Int_t XMaxList[]={40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};
      Bool_t debug = 0;
      Int_t wfsize = 127;
      Int_t FirstEvent=0;
      Int_t EventsToProcess=-1;
      Int_t nfit=1;
      Bool_t bShift=0;
      Bool_t bFast=1;
      Bool_t bUseMWPC=0;

      for(Int_t idxRun=183; idxRun<=345; idxRun++){
	Int_t iRun = idxRun;

	TStopwatch timer;
	timer.Start();

	SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);

	timer.Stop();
	timer.Print();

      }
    }
    break;


    case 300: {

      path = "/run/media/vito/dati/T-1015/TB112015/";
      subdataset = 3;
      const int Nch=20; //How many channels have information you are looking at in them
      int chlist[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}; //channel list of interest
      Bool_t IsSiPMChList[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      Int_t XMaxList[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 40, 40, 40, 40};
      Bool_t debug = 0;
      Int_t wfsize = 127;
      Int_t FirstEvent=0;
      Int_t EventsToProcess=-1;
      Int_t nfit=1;
      Bool_t bShift=0;
      Bool_t bFast=1;
      Bool_t bUseMWPC=0;

      for(Int_t idxRun=801; idxRun<=818; idxRun++){
	Int_t iRun = idxRun;

	TStopwatch timer;
	timer.Start();

	SignalIntegralFEBV2(iRun, Nch, chlist, IsSiPMChList, XMaxList, FirstEvent, EventsToProcess, dataset, subdataset, path, debug, wfsize, nfit, bShift, bFast, bUseMWPC);

	timer.Stop();
	timer.Print();

      }
    }
    break;



    default:
      return;

  }
  s.Stop();
  s.Print();
}


