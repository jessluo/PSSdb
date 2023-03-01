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

day_night = input('Would you like to group data by day/night? Enter Y/N \n' )

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
files_data = proj_id_list_func(instrument, data_status='standardized')  # generate path and project ID's
test_1 = int(input('Would you like to test step3 and step4 on a small number of standardized files? \n Enter how many files you want to run through the pipeline, 0 for all files'))
if test_1 == 0:
    files_data = files_data
else:
    files_data = files_data[0:test_1]
#if instrument == 'IFCB':  # this removal is only for testing,
        #since the standardizer should be able to deal with projects with empty rows
    #IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
    #for element in IFCB_empty_rows:
        #if element in df_list:
            #df_list.remove(element)

path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

gridpath = Path(cfg['raw_dir']).expanduser() / 'gridded_data'
if not os.path.exists(gridpath):
    os.mkdir(gridpath)

dirpath = Path(cfg['raw_dir']).expanduser() / 'gridded_data' / instrument

if os.path.isdir(dirpath) and len(os.listdir(dirpath)) != 0:  # and  os.path.exists(path_download)
    replace = input('There is already gridded data in ' + str(dirpath) + ' do you want to replace the files? \n Y/N')
    if replace == 'Y':
        print('Overwriting gridded file(s), please wait')
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
        for i in tqdm(files_data):
            filename= i.split('/')[-1]
            print('gridding and binning ' + i)
            df = pd.read_csv(i, header=0)
            df = df.dropna(subset=['Latitude']).reset_index(drop=True)
            if len(df) == 0:
                print('no data left after removing ROIS with no Lat-lon information in' + filename)
                continue
            df = df[df['Area'] != 0].reset_index(drop=True)
            if instrument == 'UVP':
                df = remove_UVP_noVal(df)
                if len(df) == 0:
                    print ('no data left after assessing validation status for ' + filename)
                    continue
            if (instrument == 'Zooscan') or (instrument == 'UVP'):
                df = depth_parsing_func(df, instrument)
                if len(df) == 0:
                    print('no data left after restricting depths to less than 200 meters in ' + filename)
                    continue
            elif instrument == 'IFCB':
                df = depth_parsing_func(df, instrument)
                df = df.replace(np.nan, '')
                sampling_type_to_remove = ['test', 'exp', 'junk', 'culture']
                df = df[~df.Sampling_type.str.contains('|'.join(sampling_type_to_remove))]
                df = df.replace('', np.nan)
                if len(df) == 0:
                    print('no data left after restricting depths to less than 200 meters in ' + filename)
                    continue
            df = biovol_func(df, instrument, keep_cat='none')
            if day_night == 'Y':
                df = date_binning_func(df, group_by=date_group, day_night=True, ignore_high_lat=True)
            elif day_night == 'N':
                df = date_binning_func(df, group_by=date_group, day_night=False, ignore_high_lat=True)
            if len(df)==0:
                print('no data left after removing data without proper time stamp in ' + filename) # necessary because some projects don't have time info: '/Users/mc4214/GIT/PSSdb/raw/raw_standardized/ecotaxa/Zooscan/standardized_project_5785_20221103_1928.csv'
                continue
            df['date_bin'] = df['date_bin'].astype(str)
            df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(st_increment, df['Latitude'], df['Longitude'])
            if depth_binning == 'N':
                metadata_bins = pd.DataFrame(
                    {'Variables': ['Biovolume', 'date_bin',  'Station_location', 'midLatBin', 'midLonBin'], # 'light_cond',
                     'Variable_types': ['float64', 'int64',  'object', 'float64', 'float64'], # 'str',
                     'Units/Values/Timezone': ['cubic_micrometer', date_group,  'lat_lon', 'degree', 'degree'], #'',
                     'Description': ['Biovolume calculated as a spherical projection of Area in cubic micrometers',
                                     'binned date information',
                                     #'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                                     'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                     'latitude of the center point of the 1x1 degree cell',
                                     'longitude of the center point of the 1x1 degree cell']})
            elif depth_binning == 'Y':
                df['midDepthBin'] = depth_binning_func(df['Depth_max'], depth_bins=depth_bins)
                metadata_bins = pd.DataFrame({'Variables': ['Biovolume', 'date_bin',  'Station_location', 'midLatBin', # 'light_cond',
                                                            'midLonBin', 'midDepthBin'],
                                              'Variable_types': ['float64', 'int64',  'object', 'float64', 'float64', #'str',
                                                                 'float64'],
                                              'Units/Values/Timezone': ['cubic_micrometer', date_group, 'lat_lon', # '',
                                                                        'degree', 'degree', 'meters'],
                                              'Description': [
                                                  'Biovolume calculated as a spherical projection of of Area in cubic micrometers',
                                                  'binned date information',
                                                  #'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                                                  'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                                  'latitude of the center point of the 1x1 degree cell',
                                                  'longitude of the center point of the 1x1 degree cell',
                                                  'middle point within a depth bin']})

            filename = filename.replace("standardized", "gridded")
            df.to_csv(str(dirpath) +'/'+ instrument + '_' + filename, index=False)
            metadata_bins.to_csv(str(dirpath) +'/'+ instrument + '_' + 'metadata_' + filename, index=False)
    elif replace == 'N':
        print('previously gridded files will be kept')
elif not os.path.exists(dirpath):
    os.mkdir(dirpath)
    for i in tqdm(files_data):
        filename = i.split('/')[-1]
        print('gridding and binning ' + i)
        df = pd.read_csv(i, header=0)
        if instrument == 'UVP':
            df = remove_UVP_noVal(df)
            if len(df) == 0:
                print('no data left after assessing validation status for ' + filename)
                continue
        df = df.dropna(subset=['Latitude']).reset_index(drop=True)
        df = df[df['Area'] != 0].reset_index(drop=True)
        if (instrument == 'Zooscan') or (instrument == 'UVP'):
            df = depth_parsing_func(df, instrument)
            if len(df) == 0:
                print('no data left after restricting depths to less than 200 meters in ' + filename)
                continue
        elif instrument == 'IFCB':
            df = depth_parsing_func(df, instrument)
            df = df.replace(np.nan, '')
            sampling_type_to_remove = ['test', 'exp', 'junk', 'culture']
            df = df[~df.Sampling_type.str.contains('|'.join(sampling_type_to_remove))]
            df = df.replace('', np.nan)
            if len(df) == 0:
                print('no data left after restricting depths to less than 200 meters in ' + filename)
                continue
        df = biovol_func(df, instrument, keep_cat='none')
        if day_night == 'Y':
            df = date_binning_func(df, group_by=date_group, day_night=True, ignore_high_lat=True)
        elif day_night == 'N':
            df = date_binning_func(df, group_by=date_group, day_night=False, ignore_high_lat=True)
        df['date_bin'] = df['date_bin'].astype(str)
        df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(st_increment, df['Latitude'],df['Longitude'])
        if depth_binning == 'N':
            metadata_bins = pd.DataFrame(
                {'Variables': ['Biovolume', 'date_bin',  'Station_location', 'midLatBin', 'midLonBin'], # 'light_cond',
                 'Variable_types': ['float64', 'int64',  'object', 'float64', 'float64'], # 'str',
                 'Units/Values/Timezone': ['cubic_micrometer', date_group,  'lat_lon', 'degree', 'degree'], # '',
                 'Description': ['Biovolume calculated as a spherical projection of Area in cubic micrometers',
                                 'binned date information',
                                 #'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                 'latitude of the center point of the 1x1 degree cell',
                                 'longitude of the center point of the 1x1 degree cell']})
        elif depth_binning == 'Y':
            df['midDepthBin'] = depth_binning_func(df['Depth_max'], depth_bins=depth_bins)
            metadata_bins = pd.DataFrame(
                {'Variables': ['Biovolume', 'date_bin',  'Station_location', 'midLatBin', # 'light_cond',
                               'midLonBin', 'midDepthBin'],
                 'Variable_types': ['float64', 'int64',  'object', 'float64', 'float64', # 'str',
                                    'float64'],
                 'Units/Values/Timezone': ['cubic_micrometer', date_group,  'lat_lon', # '',
                                           'degree', 'degree', 'meters'],
                 'Description': [
                     'Biovolume calculated as a spherical projection of of Area in cubic micrometers',
                     'binned date information',
                     #'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                     'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                     'latitude of the center point of the 1x1 degree cell',
                     'longitude of the center point of the 1x1 degree cell',
                     'middle point within a depth bin']})
        filename = filename.replace("standardized", "gridded")
        df.to_csv(str(dirpath)  + '/'+ instrument + '_' + filename, index=False)
        metadata_bins.to_csv(str(dirpath)  + '/'+ instrument + '_' +'metadata_'+ filename, index=False)


    #return df, metadata_bins
