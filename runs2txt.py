#prints out run number for 'FEBDataConverterV4.C' macro and stores it in .txt

from os import path
import numpy as np
import shutil

#datadir = '/Volumes/G-Drive Mobile USB 2/T1604_data_2019_Dec/' #location where data is stored
datadir = '/Users/justinvega/Documents/Fermilab/src/adriano_analysis/ntup/'

def get_runs():
    events = open(datadir+ 'run_list3.txt')
    lines = events.readlines()
    runs =[]

    for i in lines:
        runs.append(int(i[8:12].strip('_')))

    runs = np.array(runs)

def write_txt():
    txtfile = open(datadir + "runs.txt", "x")

    for d in lines:
        txtfile.write("{run} \n".format(run = d[8:12].strip('_')))

    txtfile.close()

def get_value(RunNo):
    value = np.where(runs == RunNo)
    print(value)

#print(len(runs)) 358

success = 0
datcorrupt = 0
hdrcorrupt = 0

def datareadout():
    log = open(datadir+"ntup/log.txt")
    lines = log.readlines()
    for i in lines:
        if ("success" in i):
            success +=1
        elif ("header corrupted" in i):
            datcorrupt +=1
        elif ("data corrupted" in i):
            hdrcorrupt +=1

    print("{success} were successful, {data} had corrupted data, and {hdr} had corrupted headers".format(success=success, data=datcorrupt, hdr=hdrcorrupt))

def get_success():
    log = open(datadir+"log.txt")
    lines = log.readlines()
    txtfile = open(datadir + "success.txt", "w")

    for d in lines:
        if ("success" in d):
            txtfile.write("{run} \n".format(run = d[4:8].strip(' ')))
        else:
            continue

    txtfile.close()

def null_waveforms():
    datadir = "/Users/justinvega/Documents/Fermilab/src/adriano_analysis/wfoutput/"
    wflist = open(datadir+"null_waveforms.txt")
    runslist = []

    lines = wflist.readlines()
    for i in lines:
        RunNo = i[71:75].strip("_")
        runslist.append(int(RunNo))

    print(list(np.unique(runslist)))

# def copy_success():
#     log = open(datadir+"log.txt")
#     lines = log.readlines()
#     txtfile = open(datadir + "success.txt", "r")
#
#     lines = txtfile.readlines()
#
#     for i in lines:
#         if i in ______:
#
