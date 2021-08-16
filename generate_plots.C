/**** macro that generates plots of up to 32 channels per run and saves it
This is a good way to visualize all of the waveforms in a test beam in a pretty quick way
you can run ROOT in batch mode ('root -b') so that it doesn't plot histograms and the macro works faster for you

*optional grab all of the successfully converted runs from a txt file and convert them into ints of a vector

To run:
first load Waveforms_OutputFEBV1 (.L Waveforms_OutputFEBV1.C++)
then execute generate_plots (.x generate_plots.C)

****/
void generate_plots(){

TString datadir = "./ntup/";
// std::vector<int> successlist;
// int num;
//
// // grab list of RunNo from success txt file
// ifstream myfile(datadir + "success.txt");
//
// if (!myfile) {
//   cout<<"Error opening output file"<<endl;
// }
// while(myfile >> num){ //adds RunNo from each line of file to runslist
//   successlist.push_back(num);
// }
// myfile.close();

// opens MB4exist.txt for readout
//ofstream mother;
//mother.open(datadir + "MB4exist.txt", ios::trunc);

for (int i = 0; i < 1; i++){
  Waveforms_Output(1025, datadir, 1, 0);
}
}
