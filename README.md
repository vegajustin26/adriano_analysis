This repo is the data analysis pipeline needed to process, visualize, and extract waveform data from the ADRIANO2 calorimeter prototype of REDTOP. All scripts should have some documentation to them, but here's a quick rundown of what everything does:

### Libraries
`DRCEventV5.C` has all of the functions needed to convert waveforms from binary into DRCEvent format, it's mainly used in FEBDataConverterV4.C

`functions.C` has fitting functions, functions to get pedestal values, and other useful things

## Analysis
`FEBDataConverterV4.C` - convert binary output from DAQ into a root file in DRCEvent format

`data_convert.C` - basically runs FEBDataConverter in a loop and has some form of error handling, just in case FEBDataConverter crashes

`WaveformFEBV1.C` - visualize waveforms of up to 32 channels of the DAQ given a DRCEvent .root file

`generate_plots.C` - loops WaveformFebV1

`Waveforms_OutputFEBV1.C` - basically does the same thing as WaveformFEBV1, but with additional channel support and some useful event filtering

`SignalIntegralFEBV3.C` - processes the waveforms and extracts observables into ROOT n-tuples

`AnalysisIntegralFEBV2` - loops SignalIntegral for specific cases

`plotting.C` - my little nifty script that has useful functions for analysis

`FitWaveformFEBV1.C` - visualize what the fit of a waveform is doing for a single event

This work was done under the DOE SULI Program at Fermilab in Summer 2021 by Justin Vega.
