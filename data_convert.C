/****
 * ROOT macro to convert .data binary files into ROOT files

first load FEBDataConverterV4.C++ (.L FEBDataConverterV4.C++)
then execute data_convert.C (.x data_convert.C)

FEBDataConverterV4 does all the work, this just loops the process until it crashes
optional save to a txt file, so that you can see which ones converted and which ones failed

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
#include<vector>


void data_convert() {

  TString datadir = "/Volumes/G-Drive Mobile USB 2/T1604_data_2019_Dec/";
  std::vector<int> runslist;
  int num;

  // grab list of RunNo from txt file
  ifstream myfile(datadir + "runs.txt");

  if (!myfile) {
    cout<<"Error opening output file"<<endl;
  }
  while(myfile >> num){ //adds RunNo from each line of file to runslist
    runslist.push_back(num);
  }
  myfile.close();

  // data convert + error handling
  std::vector<int> datacrpt;
  std::vector<int> hdrcrpt;
  std::vector<int> success;

  // opens log.txt for readout
  ofstream log;
  log.open(datadir + "ntup/log.txt", ios::app);

  // cycles through all the run numbers
  for (int i = 169; i < runslist.size(); i++){
    int result = FEBDataConverter(runslist[i], datadir);
    switch(result){
      case 100: //no run
        cout << "there's no run here";
        continue;
      case 101: //header corrupt
        log << Form("Run %d header corrupted!", runslist[i]) << endl;
        //hdrcrpt.push_back(runslist[i]);
        //cout << hdrcrpt[i] << endl;
        continue;
      case 102: //data corrupt
        log << Form("Run %d data corrupted!", runslist[i]) << endl;
        //datacrpt.push_back(runslist[i]);
        //cout << datacrpt[i] << endl;
        continue;
      case 103: //success!
        log << Form("Run %d was successful!", runslist[i]) << endl;
        //success.push_back(runslist[i]);
        continue;
    }
  }
log.close();
}
