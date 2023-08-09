# gridding standardized files generated in 2.

## load packages
import numpy as np
import pandas as pd
import os
try:
    from funcs_read import * # functions to list standardized and binned files
    from funcs_gridding import *
except:
    from scripts.funcs_read import * # functions to list standardized and binned files
    from scripts.funcs_gridding import *
from pathlib import Path
from glob import glob
import yaml  # requires installation of PyYAML package
import os
import shutil
from tqdm import tqdm

path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

## set up parameters for script:
cats = cfg['keep_cats']
if cats == 'Y':
    keep_cat= cfg['cats_to_keep']
elif cats == 'N':
    keep_cat = 'none'

day_night = cfg['day_night']
st_increment = float(cfg['st_increment'])
date_group = cfg['date_group']
depth_binning = cfg['depth_binning']
if depth_binning == 'Y':
    depth_bins = cfg['depth_bins']

## processing starts here
for instrument in ['Scanner', 'UVP', 'IFCB']:
    files_data = proj_id_list_func(instrument, data_status='standardized')  # generate path and project ID's
    test_1 = cfg['N_test']
    if test_1 == 0:
        files_data = files_data
    else:
        files_data = files_data[0:test_1]
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
            n_del_files = 0
            for i in tqdm(files_data):
                filename= i.split('/')[-1]
                print('gridding and binning ' + i)
                df = pd.read_csv(i, header=0)
                df = df.dropna(subset=['Latitude']).reset_index(drop=True)
                if len(df) == 0:
                    print('no data left after removing ROIS with no Lat-lon information in' + filename)
                    n_del_files += 1
                    continue
                df = df[df['Area'] != 0].reset_index(drop=True)
            #if instrument == 'UVP':
                #df = remove_UVP_noVal(df)
                #if len(df) == 0:
                    #print ('no data left after assessing validation status for ' + filename)
                    #n_del_files += 1
                    #continue
                if (instrument == 'Scanner') or (instrument == 'UVP'):
                    df = depth_parsing_func(df, instrument)
                    if len(df) == 0:
                        print('no data left after restricting depths to less than 200 meters in ' + filename)
                        n_del_files += 1
                        continue
                elif instrument == 'IFCB':
                    df = depth_parsing_func(df, instrument)
                    df = df.replace(np.nan, '')
                    sampling_type_to_remove = ['test', 'exp', 'junk', 'culture']
                    df = df[~df.Sampling_type.str.contains('|'.join(sampling_type_to_remove))]
                    df = df.replace('', np.nan)
                    if len(df) == 0:
                        print('no data left after restricting depths to less than 200 meters in ' + filename)
                        n_del_files += 1
                        continue
                df = biovol_func(df, instrument, keep_cat='none')
                if day_night == 'Y':
                    df = date_binning_func(df, group_by=date_group, day_night=True)
                elif day_night == 'N':
                    df = date_binning_func(df, group_by=date_group, day_night=False)
                if len(df)==0:
                    print('no data left after removing data without proper time stamp in ' + filename) # necessary because some projects don't have time info: '/Users/mc4214/GIT/PSSdb/raw/raw_standardized/ecotaxa/Zooscan/standardized_project_5785_20221103_1928.csv'
                    n_del_files += 1
                    continue
                df['date_bin'] = df['date_bin'].astype(str)
                df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(st_increment, df['Latitude'], df['Longitude'])
                filename = filename.replace("standardized", "gridded")
                df.to_csv(str(dirpath) +'/'+ instrument + '_' + filename, index=False)
            print(str(n_del_files) + ' files were not gridded due to missing Lat/Lon, deep water sample, time information or validation status < 95% for UVP')
            if depth_binning == 'N':
                metadata_bins = pd.DataFrame(
                    {'Variables': ['Biovolume', 'date_bin', 'Station_location', 'midLatBin', 'midLonBin'],
                     # 'light_cond',
                     'Variable_types': ['float64', 'int64', 'object', 'float64', 'float64'],  # 'str',
                     'Units/Values/Timezone': ['cubic_micrometer', date_group, 'lat_lon', 'degree', 'degree'],  # '',
                     'Description': ['Biovolume calculated as a spherical projection of Area in cubic micrometers',
                                     'binned date information',
                                     # 'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                                     'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                     'latitude of the center point of the 1x1 degree cell',
                                     'longitude of the center point of the 1x1 degree cell']})
            elif depth_binning == 'Y':
                df['midDepthBin'] = depth_binning_func(df['Depth_max'], depth_bins=depth_bins)
                metadata_bins = pd.DataFrame(
                    {'Variables': ['Biovolume', 'date_bin', 'Station_location', 'midLatBin','midLonBin', 'midDepthBin'],  # 'light_cond',
                     'Variable_types': ['float64', 'int64', 'object', 'float64', 'float64', 'float64'], # 'str',
                     'Units/Values/Timezone': ['cubic_micrometer', date_group, 'lat_lon', 'degree', 'degree', 'meters'], # '',
                     'Description': ['Biovolume calculated as a spherical projection of of Area in cubic micrometers',
                                     'binned date information',
                                     # 'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                                    'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                    'latitude of the center point of the 1x1 degree cell',
                                    'longitude of the center point of the 1x1 degree cell',
                                    'middle point within a depth bin']})
            metadata_bins.to_csv(str(dirpath) + '/' + instrument + '_' + 'metadata_gridded.csv', index=False)
            grid_list = group_gridded_files_func(instrument, already_gridded='N')
        elif replace == 'N':
            print('previously gridded files will be kept')

    elif not os.path.exists(dirpath):
        os.mkdir(dirpath)
        n_del_files = 0
        for i in tqdm(files_data):
            filename = i.split('/')[-1]
            print('gridding and binning ' + i)
            df = pd.read_csv(i, header=0)
        #if instrument == 'UVP':
            #df = remove_UVP_noVal(df)
            #if len(df) == 0:
                #print('no data left after assessing validation status for ' + filename)
                #n_del_files += 1
                #continue
            df = df.dropna(subset=['Latitude']).reset_index(drop=True)
            df = df[df['Area'] != 0].reset_index(drop=True)
            if (instrument == 'Scanner') or (instrument == 'UVP'):
                df = depth_parsing_func(df, instrument)
                if len(df) == 0:
                    print('no data left after restricting depths to less than 200 meters in ' + filename)
                    n_del_files += 1
                    continue
            elif instrument == 'IFCB':
                df = depth_parsing_func(df, instrument)
                df = df.replace(np.nan, '')
                sampling_type_to_remove = ['test', 'exp', 'junk', 'culture']
                df = df[~df.Sampling_type.str.contains('|'.join(sampling_type_to_remove))]
                df = df.replace('', np.nan)
                if len(df) == 0:
                    print('no data left after restricting depths to less than 200 meters in ' + filename)
                    n_del_files += 1
                    continue
            df = biovol_func(df, instrument, keep_cat='none')
            if day_night == 'Y':
                df = date_binning_func(df, group_by=date_group, day_night=True)
            elif day_night == 'N':
                df = date_binning_func(df, group_by=date_group, day_night=False)
            if len(df) == 0:
                print('no data left after removing data without proper time stamp in ' + filename)  # necessary because some projects don't have time info: '/Users/mc4214/GIT/PSSdb/raw/raw_standardized/ecotaxa/Zooscan/standardized_project_5785_20221103_1928.csv'
                n_del_files += 1
                continue
            df['date_bin'] = df['date_bin'].astype(str)
            df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(st_increment, df['Latitude'],df['Longitude'])
            filename = filename.replace("standardized", "gridded")
            df.to_csv(str(dirpath)  + '/'+ instrument + '_' + filename, index=False)
        print(str(n_del_files) + ' files were not gridded due to missing Lat/Lon, deep water sample, time information or validation status < 95% for UVP')
        if depth_binning == 'N':
            metadata_bins = pd.DataFrame(
                {'Variables': ['Biovolume', 'date_bin', 'Station_location', 'midLatBin', 'midLonBin'],
                 # 'light_cond',
                 'Variable_types': ['float64', 'int64', 'object', 'float64', 'float64'],  # 'str',
                 'Units/Values/Timezone': ['cubic_micrometer', date_group, 'lat_lon', 'degree', 'degree'],  # '',
                 'Description': ['Biovolume calculated as a spherical projection of Area in cubic micrometers',
                                 'binned date information',
                                 # 'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                 'latitude of the center point of the 1x1 degree cell',
                                 'longitude of the center point of the 1x1 degree cell']})
        elif depth_binning == 'Y':
            df['midDepthBin'] = depth_binning_func(df['Depth_max'], depth_bins=depth_bins)
            metadata_bins = pd.DataFrame(
                {'Variables': ['Biovolume', 'date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'midDepthBin'],
                 # 'light_cond',
                 'Variable_types': ['float64', 'int64', 'object', 'float64', 'float64', 'float64'],
                 # 'str',
                 'Units/Values/Timezone': ['cubic_micrometer', date_group, 'lat_lon', 'degree', 'degree', 'meters'],
                 # '',
                 'Description': ['Biovolume calculated as a spherical projection of of Area in cubic micrometers',
                                 'binned date information',
                                 # 'sampling occured during daylight or nightime. Considers time zone and geopgraphical information',
                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                 'latitude of the center point of the 1x1 degree cell',
                                 'longitude of the center point of the 1x1 degree cell',
                                 'middle point within a depth bin']})
        metadata_bins.to_csv(str(dirpath) + '/' + instrument + '_' + 'metadata_gridded.csv', index=False)

        grid_list = group_gridded_files_func(instrument, already_gridded='N')
