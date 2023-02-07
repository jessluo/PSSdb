# gridding standardized files generated in 2.

## load packages
import numpy as np
import pandas as pd
import os
from funcs_read import * # functions to list standardized and binned files
from funcs_gridding import *
from pathlib import Path
from glob import glob
import yaml  # requires installation of PyYAML package
import os
import shutil
from tqdm import tqdm

## set up queries :
instrument = input ('for which instrument do you want to grid the standardized dataset? \n Enter IFCB, Zooscan or UVP ')
cats = input ('the gridding process by default removes ROIs that have been identified as \n '
                    'BEADS, DETRITUS, ARTEFACTS or BUBBLES, would you like to keep any of these categories?'
                    '  \n Enter Y/N ')
st_increment = float(input('select the size of the spatial grid to group the data i.e for a 1x1 degree bin type 1 \n' ))

if cats == 'Y':
    keep_cat=[]
    n = int(input("how many categories do you want to keep? : "))
    # iterating till the range
    for i in range(0, n):
        ele = input('type category to keep')
        keep_cat.append(ele)  # adding the element
elif cats == 'N':
    keep_cat = 'none'

date_group = input ('input how will the data be grouped by date \n  yyyymm for month and year \n yyyy for year \n mm by month \n week for weekly average for each month or \n None to keep all time information as is')
depth_binning = input ('Would you like to bin the data by depth? \n Enter Y/N')
if depth_binning == 'Y':
    custom_depth_bin = input ('default depth binning is done at depths: 0, 25, 50, 100, 200, 500, 1000, 3000, 8000 \n '
                              'Would you like to bin data using different depths \n Enter Y/N ')
    if custom_depth_bin == 'Y':
        depth_bins = []
        n = int(input("how many depths would you like to use for binning? : "))
        # iterating till the range
        for i in range(0, n):
            ele = int(input('type depth #' + str(i+1) ))
            depth_bins .append(ele)  # adding the element
    elif custom_depth_bin == 'N':
        depth_bins = [0, 25, 50, 100, 200, 500, 1000, 3000, 8000]


## processing starts here
path_to_data, df_list = proj_id_list_func(instrument, data_status='standardized')  # generate path and project ID's
#if instrument == 'IFCB':  # this removal is only for testing,
        #since the standardizer should be able to deal with projects with empty rows
    #IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
    #for element in IFCB_empty_rows:
        #if element in df_list:
            #df_list.remove(element)
dirpath = str(path_to_data) + '/gridded_data/'

if os.path.isdir(dirpath) and len(os.listdir(dirpath)) != 0:  # and  os.path.exists(path_download)
    replace = input('There is already gridded data in ' + dirpath + ' do you want to replace the files? \n Y/N')
    if replace == 'Y':
        print('Overwriting gridded file(s), please wait')
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
        for i in tqdm(df_list):
            print('gridding and binning ' + i)
            df = read_func(path_to_data, i, data_status='standardized')
            df = df.dropna(subset=['Latitude']).reset_index(drop=True)
            df = df[df['Area'] != 0].reset_index(drop=True)
            df = biovol_func(df, instrument, keep_cat='none')
            df['date_bin'], df['light_cond'] = date_binning_func(df['Sampling_date'], df['Sampling_time'], df['Latitude'], df['Longitude'], group_by=date_group)
            df['date_bin'] = df['date_bin'].astype(str)
            df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(st_increment, df['Latitude'], df['Longitude'])
            if depth_binning == 'N':
                metadata_bins = pd.DataFrame(
                    {'Variables': ['Biovolume', 'date_bin', 'light_cond', 'Station_location', 'midLatBin', 'midLonBin'],
                     'Variable_types': ['float64', 'int64', 'str', 'object', 'float64', 'float64'],
                     'Units/Values/Timezone': ['cubic_micrometer', date_group, '', 'lat_lon', 'degree', 'degree'],
                     'Description': ['Biovolume calculated as a spherical projection of Area in cubic micrometers',
                                     'binned date information',
                                     'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                                     'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                     'latitude of the center point of the 1x1 degree cell',
                                     'longitude of the center point of the 1x1 degree cell']})
            elif depth_binning == 'Y':
                df['midDepthBin'] = depth_binning_func(df['Depth_max'], depth_bins=depth_bins)
                metadata_bins = pd.DataFrame({'Variables': ['Biovolume', 'date_bin', 'light_cond', 'Station_location', 'midLatBin',
                                                            'midLonBin', 'midDepthBin'],
                                              'Variable_types': ['float64', 'int64', 'str', 'object', 'float64', 'float64',
                                                                 'float64'],
                                              'Units/Values/Timezone': ['cubic_micrometer', date_group, '', 'lat_lon',
                                                                        'degree', 'degree', 'meters'],
                                              'Description': [
                                                  'Biovolume calculated as a spherical projection of of Area in cubic micrometers',
                                                  'binned date information',
                                                  'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                                                  'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                                  'latitude of the center point of the 1x1 degree cell',
                                                  'longitude of the center point of the 1x1 degree cell',
                                                  'middle point within a depth bin']})

            i = i.replace("standardized", "gridded")
            i = i.replace(".tsv", ".csv")
            df.to_csv(str(path_to_data) + '/gridded_data/' + instrument + '_' + str(i), index=False)
            metadata_bins.to_csv(str(path_to_data) + '/gridded_data/' + instrument + '_' + 'metadata_' + str(i), index=False)
    elif replace == 'N':
        print('previously gridded files will be kept')
elif not os.path.exists(dirpath):
    os.mkdir(dirpath)
    for i in tqdm(df_list):
        print('gridding and binning ' + i)
        df = read_func(path_to_data, i)
        df = df.dropna(subset=['Latitude']).reset_index(drop=True)
        df = df[df['Area'] != 0].reset_index(drop=True)
        df = biovol_func(df, instrument, keep_cat='none')
        df['date_bin'], df['light_cond'] = date_binning_func(df['Sampling_date'], df['Sampling_time'], df['Latitude'], df['Longitude'], group_by=date_group)
        df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(st_increment, df['Latitude'], df['Longitude'])
        if depth_binning == 'N':
            metadata_bins = pd.DataFrame(
                {'Variables': ['Biovolume', 'date_bin', 'light_cond', 'Station_location', 'midLatBin', 'midLonBin'],
                 'Variable_types': ['float64', 'int64', 'str', 'object', 'float64', 'float64'],
                 'Units/Values/Timezone': ['cubic_micrometer', date_group, '', 'lat_lon', 'degree', 'degree'],
                 'Description': ['Biovolume calculated as a spherical projection of Area in cubic micrometers',
                                 'binned date information',
                                 'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                 'latitude of the center point of the 1x1 degree cell',
                                 'longitude of the center point of the 1x1 degree cell']})
        elif depth_binning == 'Y':
            df['midDepthBin'] = depth_binning_func(df['Depth_max'], depth_bins=depth_bins)
            metadata_bins = pd.DataFrame(
                {'Variables': ['Biovolume', 'date_bin', 'light_cond', 'Station_location', 'midLatBin',
                               'midLonBin', 'midDepthBin'],
                 'Variable_types': ['float64', 'int64', 'str', 'object', 'float64', 'float64',
                                    'float64'],
                 'Units/Values/Timezone': ['cubic_micrometer', date_group, '', 'lat_lon',
                                           'degree', 'degree', 'meters'],
                 'Description': [
                     'Biovolume calculated as a spherical projection of of Area in cubic micrometers',
                     'binned date information',
                     'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                     'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                     'latitude of the center point of the 1x1 degree cell',
                     'longitude of the center point of the 1x1 degree cell',
                     'middle point within a depth bin']})
        i = i.replace("standardized", "gridded")
        i = i.replace('.tsv', '.csv')
        df.to_csv(str(path_to_data) + '/gridded_data/' + instrument + '_' + str(i), index=False)
        metadata_bins.to_csv(str(path_to_data) + '/gridded_data/' + instrument + '_' +'metadata_'+ str(i),sep='\t', index=False)


    #return df, metadata_bins
