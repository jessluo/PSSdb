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

## define parameters for gridding:
instrument = 'IFCB'
removecat = 'none'
area_type = 'object_area'
date_group = 'yyyymm'
ignore_depth = 'yes'

## processing starts here
path_to_data, df_list = proj_id_list_func(instrument, data_status='standardized')  # generate path and project ID's
if instrument == 'IFCB':  # this removal is only for testing,
        #since the standardizer should be able to deal with projects with empty rows
    IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
    for element in IFCB_empty_rows:
        if element in df_list:
            df_list.remove(element)
dirpath = str(path_to_data) + '/gridded_data/'
if os.path.exists(dirpath) and os.path.isdir(dirpath):
    shutil.rmtree(dirpath)
os.mkdir(dirpath)

for i in df_list:
    df = read_func(path_to_data, i)
    df = df.dropna(subset = ['Latitude']).reset_index(drop = True)
    df = df[df['Area'] != 0].reset_index(drop=True)
    df = biovol_func(df, instrument= instrument, area_type= area_type, remove_cat=removecat)
    df['date_bin'] = date_binning_func(df['Sampling_date'], df['Sampling_time'], group_by =  date_group)
    df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(df['Latitude'], df['Longitude'])
    if ignore_depth == 'yes':
        metadata_bins = pd.DataFrame({'Variables':['Biovolume', 'date_bin', 'Station_location', 'midLatBin', 'midLonBin'],
                                 'Variable_types': ['float64',  'int64', 'object', 'float64', 'float64'],
                                 'Units/Values/Timezone': ['cubic_micrometer',  date_group, 'lat_lon', 'degree', 'degree'],
                                 'Description': ['Biovolume calculated as a spherical projection of ' + area_type ,
                                                 'binned date information',
                                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                                 'latitude of the center point of the 1x1 degree cell',
                                                 'longitude of the center point of the 1x1 degree cell']})
    elif ignore_depth == 'no':
        df['midDepthBin'] = depth_binning_func(df['Depth_max'])
        metadata_bins = pd.DataFrame({'Variables':['Biovolume',  'date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'midDepthBin'],
                                 'Variable_types': ['float64',  'int64', 'object', 'float64', 'float64', 'float64'],
                                 'Units/Values/Timezone': ['cubic_micrometer',  date_group, 'lat_lon', 'degree', 'degree', 'meters'],
                                 'Description': ['Biovolume calculated as a spherical projection of ' + area_type ,
                                                 'binned date information',
                                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                                 'latitude of the center point of the 1x1 degree cell',
                                                 'longitude of the center point of the 1x1 degree cell',
                                                 'middle point within a depth bin']})

    df.to_csv(str(path_to_data) + '/gridded_data/' + instrument + '_' + str(i) + '_gridded.tsv', sep='\t')
    metadata_bins.to_csv(str(path_to_data) + '/gridded_data/' + instrument + '_' + str(i) + '_gridded_metadata.tsv',
                           sep='\t')
    #return df, metadata_bins
