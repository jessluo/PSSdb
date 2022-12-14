# GOAL: to extract important instrument configuration that determines the 'sampling' restrictions

import pyifcb as ifcb
url = 'http://ifcb-data.whoi.edu/IFCB101_BigelowJun2016/D20160709T210837_IFCB101.html'

#with ifcb.open_url(url) as sample_bin:
    #lid = sample_bin.lid
    print('{} has {} image(s)'.format(lid, len(sample_bin.images)))

# parameters to get:

#  minBlobArea= sample_bin.headers['minimumBlobArea']

#  PMTAhighVoltage= sample_bin.headers['PMTAhighVoltage']
#  PMTBhighVoltage= sample_bin.headers['PMTBhighVoltage']
#  PMTChighVoltage= sample_bin.headers['PMTChighVoltage']

#  SyringeSampleVolume= sample_bin.headers['SyringeSampleVolume']

#  PMTAtriggerThreshold_DAQ_MCConly= sample_bin.headers['PMTAtriggerThreshold_DAQ_MCConly']
#  PMTBtriggerThreshold_DAQ_MCConly= sample_bin.headers['PMTBtriggerThreshold_DAQ_MCConly']
#  PMTCtriggerThreshold_DAQ_MCConly= sample_bin.headers['PMTCtriggerThreshold_DAQ_MCConly']
#  PMTDtriggerThreshold_DAQ_MCConly= sample_bin.headers['PMTDtriggerThreshold_DAQ_MCConly']