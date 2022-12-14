# GOAL: to extract important instrument configuration that determines the 'sampling' restrictions

import ifcb # follow instructions to install on : https://github.com/joefutrelle/pyifcb. Make sure to clone the git repository, and run the commands from the directory where the repo files are downloaded
url = 'https://ifcb-data.whoi.edu/EXPORTS/D20210511T025400_IFCB125.hdr'
#sample_bin = ifcb.open_raw(url)



#with sample_bin:
    #all_data = sample_bin.read()

with ifcb.open_url(url, images=False) as sample_bin:
    #lid = sample_bin.lid
    minBlobArea = sample_bin.headers['minimumBlobArea']
    PMTAhighVoltage = sample_bin.headers['PMTAhighVoltage']
    PMTBhighVoltage= sample_bin.headers['PMTBhighVoltage']
    PMTChighVoltage= sample_bin.headers['PMTChighVoltage']
    SyringeSampleVolume = sample_bin.headers['SyringeSampleVolume']
    PMTAtriggerThreshold_DAQ_MCConly= sample_bin.headers['PMTAtriggerThreshold_DAQ_MCConly']
    PMTBtriggerThreshold_DAQ_MCConly= sample_bin.headers['PMTBtriggerThreshold_DAQ_MCConly']
    PMTCtriggerThreshold_DAQ_MCConly= sample_bin.headers['PMTCtriggerThreshold_DAQ_MCConly']
    PMTDtriggerThreshold_DAQ_MCConly= sample_bin.headers['PMTDtriggerThreshold_DAQ_MCConly']

    #print('{} has {} image(s)'.format(lid, len(sample_bin.images)))

# parameters to get:



#  PMTAhighVoltage= sample_bin.headers['PMTAhighVoltage']


#

