#Goal: create script that will iterate through downloaded projects and create a dataframe
# with the project name, the instrument used to collect data and lat & long

# We'll eventually have to convert this into a function I assume
# def LatLong (pathname):

import numpy as np
import pandas as pd
import os

pathname = '/Users/mc4214/GIT/PSSdb_sandbox/raw_data/' #replace with the local directory that contains
# the tsv data files
os.chdir(pathname)

# create empty lists for each of the metadata we want to get
lat = []
lon = []
files = []
instrument = []
cruise = []
per_validated = []
for file in os.listdir(pathname):
    if file.endswith(".tsv"):

        file_data = os.path.join(pathname, file)  # creating the element pathname allows to
        # combine the filenames with the directory
        data = pd.read_csv(file_data, sep='\t', encoding='latin_1', encoding_errors='ignore', header=0)
        files.append(file)
        instrument.append(data.at[0, 'acq_instrument'])
        # cruise.append(data.at[0, 'sample_cruise']) had to take this one out, some data sheets dont have it
        lat.append(data.at[0, 'object_lat'])
        lon.append(data.at[0, 'object_lon'])
        t = data['object_annotation_status'].value_counts()
        per_validated.append((t[0]/data.shape[0])*100)

df = pd.DataFrame (list(zip(files, instrument, lat, lon, per_validated)),
                           columns = ['File_name', 'instrument', 'lat', 'lon', '%_validation'])
df.to_csv(os.path.join(pathname, 'Ecotaxa_metadata.tsv'))
