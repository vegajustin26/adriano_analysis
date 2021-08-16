/****
 * ROOT macro to convert .data raw files into ROOT files
 * this supports up to 4 FEB boards
 * Arguments are:
 * RunNumber * run number to process
 * Path * path where raw data are located
 * debug * set on/off the debug flag
 *
 * To run this macro it is needed to compile DRCEventV5.C
 * At ROOT prompt run
 * root [0] .L DRCEventV5.C++
 * root [1] .L FEBDataConverterV4.C++
 *
 * If the macro has already been compiled, it is enough to load its library.
 *
 * root [0] gSystem->Load("FEBDataConverterV4_C.so");
 *
 * then the macro can be executes as:
 * FEBDataConverter(1190,"/data/TB_2019_dec/")
 *
 * that is equivalent to:
 * FEBDataConverter(1190,"/data/TB_2019_dec/",0)
 *
 ****/

#include "TFile.h"
#include "TTree.h"
#include <Riostream.h>
#include <TSystemDirectory.h>
#include <TString.h>
#include <TError.h>
#include <TObjString.h>

#include "DRCEventV5.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>


using namespace std;


void FEBDataConverter(Int_t RunNumber, TString Path = "../data/", UInt_t debug=0)
{


  const char *ext=".data";
  TString txtfilename("");
  TString datfilename("");

  TSystemDirectory dir("data",Path.Data());

  TList*files = dir.GetListOfFiles();

  if (files) {
    TSystemFile *file;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      txtfilename = file->GetName();
      if (!file->IsDirectory() && txtfilename.EndsWith(ext)) {
        if(txtfilename.BeginsWith(Form("RUN_FEB_%d_",RunNumber))){ //RUN_FEB_%d_
            datfilename = txtfilename;
            break;
        }
      }
    }
  }

  if((datfilename.IsNull()))
  {
    cout << Form("\n\n ## No Run %d in your data ##\n", RunNumber);
    return;
  }

  TString dirpath = "ntup/";
  datfilename.Prepend(dir.GetTitle()); //+dirpath);

  TString rootfilename = datfilename;
  rootfilename.ReplaceAll(".data",".root");

  TFile f;
  f.Open(rootfilename.Data(),"RECREATE");

  DRCEvent *event = new DRCEvent;
  const Int_t nFEBS=4;
  TTree *T[nFEBS];
  for (Int_t idx=0; idx<nFEBS; idx++){
      T[idx] = new TTree(Form("DRCMB%d",idx+1),Form("DRCMB%d",idx+1));
      T[idx]->Branch("DRCEvent","DRCEvent",&event);
  }
  vector <Hit> drchits[nFEBS];

  //datfilename = txtfilename;
  ifstream InputDataStream(datfilename.Data());

  cout << Form("\n\n===>\nConverting Run %d\n", RunNumber);
  cout << Form("filename: %s \n", txtfilename.Data());

  Int_t ch=-1;              //channel number
  vector <Int_t> wdat;     //waveform data
  Int_t TimeSinceLastTrig=-1; //Time (trigger leap)
  string input;
  Int_t ich=0;
  Int_t dataval=0;
  Bool_t bChDat=0;
  Int_t Nch_x_FEB=16;

  Int_t previousChFEB[nFEBS]={-1,15,31,47};
  Int_t previousCh=-1;
  Int_t iFEB=1;
  Bool_t bNewFEB=0;

  Int_t ThisSpillWordCounts=0;
  Int_t ThisSpillTriggers=0;
  Int_t GlobalSpillIdx=0;
  Int_t ChannelMask=0;
  Int_t BoardID=0;
  Int_t SpillStatus=0;
  Int_t PrevSpillIdx=-1;
  Bool_t bDuplicatedSpill=0;
  const Int_t MaxSpillBytes=20971519;
  Int_t SpillBytesCount=0;
  Bool_t bForceSpillEnd=0;

  Int_t ThisEventWordCounts=0;
  Int_t ThisEventTimestamp=0;
  Int_t PrevTriggerIdx=0;
  Int_t TriggersIdx=0;
  Int_t WFSamples=0;
  Int_t TriggetType=0;
  Int_t EventStatus=0;


  while(true){

    bForceSpillEnd=0;
    SpillBytesCount=0;
    //Header
    getline(InputDataStream, input, '\n');
    TString sInput = input;
    //--Begin of spill
    if (sInput.Contains("--Begin of spill")){
        cout << "Begin of spill\n";
    }
    else if(sInput.Contains("-- END OF RUN ")){
        cout << "END OF RUN \n";
        break;
    }
    else{
        cout << "Run header corrupted *** ABORTING***\n";
        cout << sInput << endl;
        exit(1);
    }
    //--** SOURCE = FEB1
    getline(InputDataStream, input, '\n');
    //Spill header;
    getline(InputDataStream, input, '\n'); SpillBytesCount+=16;
    sInput=input;
    TObjArray *arrThisSpillHeader = sInput.Tokenize(" ");

    ThisSpillWordCounts =((TObjString *)(arrThisSpillHeader->At(0)))->String().Atoi()*0x1000000;
    ThisSpillWordCounts+=((TObjString *)(arrThisSpillHeader->At(1)))->String().Atoi()*0x10000;
    ThisSpillWordCounts+=((TObjString *)(arrThisSpillHeader->At(2)))->String().Atoi()*0x100;
    ThisSpillWordCounts+=((TObjString *)(arrThisSpillHeader->At(3)))->String().Atoi();

    ThisSpillTriggers =((TObjString *)(arrThisSpillHeader->At(4)))->String().Atoi()*0x1000000;
    ThisSpillTriggers+=((TObjString *)(arrThisSpillHeader->At(5)))->String().Atoi()*0x10000;
    ThisSpillTriggers+=((TObjString *)(arrThisSpillHeader->At(6)))->String().Atoi()*0x100;
    ThisSpillTriggers+=((TObjString *)(arrThisSpillHeader->At(7)))->String().Atoi();

    GlobalSpillIdx =((TObjString *)(arrThisSpillHeader->At(8)))->String().Atoi()*0x100;
    GlobalSpillIdx+=((TObjString *)(arrThisSpillHeader->At(9)))->String().Atoi();

    ChannelMask =((TObjString *)(arrThisSpillHeader->At(10)))->String().Atoi()*0x100;
    ChannelMask+=((TObjString *)(arrThisSpillHeader->At(11)))->String().Atoi();

    BoardID =((TObjString *)(arrThisSpillHeader->At(12)))->String().Atoi()*0x100;
    BoardID+=((TObjString *)(arrThisSpillHeader->At(13)))->String().Atoi();

    SpillStatus =((TObjString *)(arrThisSpillHeader->At(15)))->String().Atoi();

    cout << "ThisSpillWordCounts: " << ThisSpillWordCounts << endl;
    cout << "ThisSpillTriggers: " << ThisSpillTriggers << endl;
    cout << "GlobalSpillIdx: " << GlobalSpillIdx << endl;
    cout << "ChannelMask: " << ChannelMask << endl;
    cout << "BoardID: " << BoardID << endl;
    cout << "SpillStatus: " << SpillStatus << endl;

    if(ThisSpillWordCounts*2>MaxSpillBytes){
        cout << Form("\n*** This spill has %d bytes, but we can read only %d bytes\n", ThisSpillWordCounts*2, MaxSpillBytes);
        cout << "*** Most likely this spill is trunkated\n\n";
    }

    if(PrevSpillIdx==GlobalSpillIdx){
        bDuplicatedSpill=1;
        cout << Form("\n***Spill #%d is buplicated, it will be skipped ***", GlobalSpillIdx);
    }

    if (ThisSpillTriggers==0){
        cout << "Empty spill... skipping event section\n";
        //--wrote
        getline(InputDataStream, input, '\n');
        //--Read took
        getline(InputDataStream, input, '\n');
        //--Save took
        getline(InputDataStream, input, '\n');
        continue;
    }



    iFEB=0;
    previousCh=previousChFEB[iFEB];
    bNewFEB=0;
    PrevTriggerIdx=0;
    while(true){
        //Event header;
        streampos oldpos = InputDataStream.tellg();
        getline(InputDataStream, input, '\n'); SpillBytesCount+=16;
        sInput=input;
        if(sInput.BeginsWith("--wrote")){
            InputDataStream.seekg(oldpos);
            bForceSpillEnd=1;
            break;
        }
        TObjArray *arrThisEventHeader = sInput.Tokenize(" ");

        ThisEventWordCounts =((TObjString *)(arrThisEventHeader->At(0)))->String().Atoi()*0x100;
        ThisEventWordCounts+=((TObjString *)(arrThisEventHeader->At(1)))->String().Atoi();

        ThisEventTimestamp =((TObjString *)(arrThisEventHeader->At(2)))->String().Atoi()*0x1000000;
        ThisEventTimestamp+=((TObjString *)(arrThisEventHeader->At(3)))->String().Atoi()*0x10000;
        ThisEventTimestamp+=((TObjString *)(arrThisEventHeader->At(4)))->String().Atoi()*0x100;
        ThisEventTimestamp+=((TObjString *)(arrThisEventHeader->At(5)))->String().Atoi();

        TriggersIdx =((TObjString *)(arrThisEventHeader->At(6)))->String().Atoi()*0x1000000;
        TriggersIdx+=((TObjString *)(arrThisEventHeader->At(7)))->String().Atoi()*0x10000;
        TriggersIdx+=((TObjString *)(arrThisEventHeader->At(8)))->String().Atoi()*0x100;
        TriggersIdx+=((TObjString *)(arrThisEventHeader->At(9)))->String().Atoi();

        WFSamples =((TObjString *)(arrThisEventHeader->At(11)))->String().Atoi();

        TriggetType =((TObjString *)(arrThisEventHeader->At(13)))->String().Atoi();

        EventStatus =((TObjString *)(arrThisEventHeader->At(15)))->String().Atoi();

        if(debug){
            cout << "Event header:\n" << sInput << endl;
            cout << "ThisEventWordCounts: " << ThisEventWordCounts << endl;
            cout << "ThisEventTimestamp: " << ThisEventTimestamp*6.28e-9 << endl;
            cout << "TriggersIdx: " << TriggersIdx << endl;
            cout << "WFSamples: " << WFSamples << endl;
            cout << "TriggetType: " << TriggetType << endl;
            cout << "EventStatus: " << EventStatus << endl;
            cout << "iFEB: " << iFEB << endl;
        }

        if(TriggersIdx<PrevTriggerIdx){
            cout << "*** TriggersIdx: " << TriggersIdx << " < PrevTriggerIdx: " << PrevTriggerIdx << endl;
            iFEB++;
            PrevTriggerIdx=0;
            bNewFEB=1;
        }

        if(bNewFEB || iFEB==0){
            if(iFEB==0) iFEB=1;
            cout << "iFEB: " << iFEB << endl;
        }

        if(TriggersIdx-PrevTriggerIdx != 1){
            cout << endl << sInput << endl;
            cout << Form("PrevTriggerIdx: %d TriggersIdx: %d ThisSpillTriggers: %d SpillBytesCount: %d\n", PrevTriggerIdx, TriggersIdx, ThisSpillTriggers, SpillBytesCount);
            cout << "Missing trigger data\n\n";
            //exit(1);
        }
        PrevTriggerIdx=TriggersIdx;

        //Event section
        Int_t NChToProcess=Nch_x_FEB;
        while(NChToProcess>0){

            //channel data block
            Int_t iSample=WFSamples*2/16;
            bChDat=0;
            while(iSample>0){
                streampos oldpos = InputDataStream.tellg();
                getline(InputDataStream, input, '\n'); SpillBytesCount+=16;
                if(SpillBytesCount>MaxSpillBytes){
                    cout << input << endl;
                    cout << Form("\n*** We are trying to read %d bytes from this spill\n", SpillBytesCount);
                    cout << Form("*** but the buffer limit is %d\n", MaxSpillBytes);
                    InputDataStream.seekg(oldpos);
                    bForceSpillEnd=1;
                }

                if(bForceSpillEnd) break;

                sInput=input;

                if(0&&debug)
                    cout << sInput << endl;

                TObjArray *chData = sInput.Tokenize(" ");

                // cout << "chData->GetEntries(): " << chData->GetEntries() << endl;
                for (Int_t i = 0; i < chData->GetEntries()-1; i++){

                    if(!bChDat){
                        Int_t chCheck = ((TObjString *)(chData->At(i++)))->String().Atoi();
                        if( chCheck != 128){
                            cout << "SpillBytesCount: " << SpillBytesCount << endl;
                            cout << "Possible corruption in data file...\n";
                            cout << "expecting to read 128, but read " << chCheck << endl;
                            cout << sInput << endl;
                            cout << "*** exiting\n";
                            exit(1);
                        }
                        ich=((TObjString *)(chData->At(i++)))->String().Atoi();
                        ch=ich;
                        // cout << Form("ich: %d - previousCh: %d iFEB: %d\n", ich, previousCh, iFEB);
                        // cout << Form("## previousCh: %d ich: %d iFEB: %d\n", previousCh, ich, iFEB);
                        if (bNewFEB){
                            previousCh=previousChFEB[iFEB-1];
                            bNewFEB=0;
                        }
                        if(0&&ich-previousCh != 1){
                            cout << sInput << endl;
                            cout << Form("previousCh: %d ich: %d - %d\n", previousCh, ich, previousChFEB[iFEB-1]);
                            cout << "Possible corruption in data file... exiting\n";
                            exit(1);
                        }
                        previousCh=ich;
                        if ( ich == (16*iFEB-1) ) previousCh=previousChFEB[iFEB-1];

                        bChDat=1;
                    }

                    if(!bDuplicatedSpill){
                        dataval=((TObjString *)(chData->At(i++)))->String().Atoi()*0x100;
                        dataval+=((TObjString *)(chData->At(i)))->String().Atoi();
                        if(dataval > 0x7ff) dataval -= 0xfff;
                        wdat.push_back(dataval);
                    }

                }

                if(chData){
                    chData->Clear();
                    chData->Delete();
                    delete chData;
                }

                iSample--;

            }

            if(!bDuplicatedSpill){
                Hit h;
                h.fCh=ch;
                h.fTimeSinceLastTrig=TimeSinceLastTrig;

                for(UInt_t i5=0;i5<wdat.size();i5++)
                    h.fW.push_back(wdat[i5]);

                drchits[iFEB-1].push_back(h);

                wdat.clear();
            }

            if(bForceSpillEnd) break;

            NChToProcess--;
            // cout << "NChToProcess: " << NChToProcess << endl;
        }

        if(drchits[iFEB-1].size() != 0){
            event->Build(GlobalSpillIdx, TriggersIdx, ThisEventTimestamp, drchits[iFEB-1]);
            T[iFEB-1]->Fill();
            drchits[iFEB-1].clear();
            while(!drchits[iFEB-1].empty()){}
        }

        if(0&&TriggersIdx==ThisSpillTriggers ){
            iFEB++;
            PrevTriggerIdx=0;
            bNewFEB=1;
            SpillBytesCount=0;
        }


        if (iFEB>nFEBS){
            iFEB=1;
            PrevTriggerIdx=0;
            bNewFEB=1;
            break;
        }

        if (0&&SpillBytesCount>ThisSpillWordCounts){
            cout << " **** SpillBytesCount>ThisSpillWordCounts **** " << SpillBytesCount << " " << ThisSpillWordCounts << endl;
            iFEB++;
//             PrevTriggerIdx=0;
            bNewFEB=1;
            SpillBytesCount=0;
        }

        if(bForceSpillEnd) break;

    }

    if(bForceSpillEnd) bForceSpillEnd=0;

    PrevSpillIdx=GlobalSpillIdx;
    //--wrote
    getline(InputDataStream, input, '\n');
    cout << "wrote: " << input << endl;
    //--Read took
    getline(InputDataStream, input, '\n');
    //--Save took
    getline(InputDataStream, input, '\n');
    cout << "SpillBytesCount: " << SpillBytesCount << endl;

    cout << endl;

  }

  cout << "<===\n";

  InputDataStream.close();

  for (Int_t idx=0; idx<nFEBS; idx++){
    T[idx]->Write(0,TObject::kOverwrite);
    T[idx]->Delete();
  }
  f.Close();

}
