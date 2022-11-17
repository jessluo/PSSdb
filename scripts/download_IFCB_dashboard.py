# Process to download raw data from IFCB dashboards



import pandas as pd
import os
import shutil
import numpy as np
import yaml
from tqdm import tqdm
import time
import datetime as dt
from pathlib import Path
from funcs_IFCB_dashboard import *

# set paths:
path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
# open the metadata of the standardized files
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
all_files_data = str(Path(cfg['raw_dir']).expanduser()) + '/IFCB_dashboard_projects_list.xlsx'
IFCB_dashboard_data_path = str(Path(cfg['IFCB_dir']).expanduser())


# generate a dataframe with the timeseries list for each dashboard:
timeseries_data =  pd.read_excel(all_files_data)

# subset timeseries based on downloading the ones for testing or all

testing = input (' Do you want to download all projects or just the tests? \n Enter all or tests ')
if testing == 'tests':
    timeseries_data = timeseries_data[timeseries_data['Project_test'] == True].reset_index(drop=True)

subset = input (' Do you want to get data for all the time series or a date range or the test dates? \n Enter all, range or tests ')
if subset == 'range':
    start_date= input(' enter starting date of the data as YYYYMMDD')
    start_date= int(start_date)
    end_date = input(' enter final date of the data as YYYYMMDD')
    end_date = int(end_date)
elif subset == 'all':
    start_date = 20000101
    end_date = 21000101
elif subset == 'tests':
    start_date = 20180825
    end_date = 20180826



# download starts here
for n in range (0, len(timeseries_data)):
    #start = time.time()
    Project_source = timeseries_data.loc[n, 'dashboard_url']
    Project_ID = timeseries_data.loc[n, 'Project_ID']
    # define path for download
    path_download = IFCB_dashboard_data_path + '/' + Project_ID  ## NOTE: the first part of the directory needs to be changed for the GIT PSSdb project
    pathname = Project_source + Project_ID + '/'

    if os.path.exists(path_download) and os.path.isdir(path_download):
        shutil.rmtree(path_download)
    os.mkdir(path_download)

    # PROCESSING HERE
    METRIC = 'temperature'  # needs to be defined based on Karina's script
    # use the file list to download data and store it locally
    file_info = get_df_list_IFCB(base_url=Project_source, dataset=Project_ID, startdate=start_date, enddate=end_date) # remember this fuction has the option to define an interval of dates

    # the loop uses the file_info dataframe to get the 'pid' (file name) and download features and class scores files
    print ('extracting features files, metadata and top 5 class scores')
    for i in tqdm(file_info .loc[:, 'pid']): # timeit and progressbar can be used here
        try:
            # generate features files
            pid_id = str(i)
            features_filename = pathname + pid_id + '_features.csv'
            features_df = df_from_url(features_filename)
            # obtain metadata
            met_dict = metadata_dict(dashboard=Project_source, pid_id=pid_id)
            for c in ['texture', 'Wedge', 'Ring', 'HOG', 'moment_invariant']:
                features_df = features_df[features_df.columns.drop(list(features_df.filter(regex=c)))]
            # url = metadata_link + dataset_name +'&bin='+ str(subset_file_info.loc[i, 'pid'])
            features_df['sample_id'] = Project_ID
            features_df['datetime'] = met_dict['datetime']
            features_df['date'] = str(pd.to_datetime(met_dict['datetime']).date().strftime('%Y%m%d'))
            features_df['time'] = str(pd.to_datetime(met_dict['datetime']).time().strftime('%H%M%S'))
            features_df['pixel_to_micron'] = met_dict['scale']
            features_df['Latitude'] = met_dict['latitude']
            features_df['Longitude'] = met_dict['longitude']
            features_df['depth'] = met_dict['depth']
            features_df['vol_analyzed'] = met_dict['ml_analyzed']
            features_df['sample_type'] = met_dict['sample_type']
            features_df['sample_cruise'] = met_dict['cruise']
            features_df['number_of_rois'] = met_dict['number_of_rois']
            features_df['concentration'] = met_dict['concentration']
            # now generate dataframe for the class scores
            class_filename = pathname + str(i) + '_class_scores.csv'
            class_df = df_from_url(class_filename)
            # class_df = class_df.drop(columns = 'pid')
            # create column with taxonomic category based on class score
            try:
                for r in class_df.index.values:
                    # class_df.iloc[r, 2:len(class_df.columns)]= class_df.iloc[r, 2:len(class_df.columns)].astype(float)
                    scores = class_df.iloc[r - 1, 2:len(class_df.columns)].astype(float)
                    for col in class_df.columns:
                        if class_df.loc[r, col] == str(max(scores)):
                            features_df.loc[r, 'class'] = col
            except:
                print (i + ' does not have a class score')

            features_df.to_csv(path_download + '/' + str(i) + '_features.csv', sep='\t')
            #class_df.to_csv(path_download + '/' + str(i) + '_class_scores.csv', sep='\t') 11/16/2022 decided not to keep class scores
            print(str(i) + ' download done ')

        except:
            print('there is no features or class_scores files for ' + str(file_info.loc[i, 'pid']))

   # elapsed_time_fl = (time.time() - start)
    # test = 'D20150910T212605_IFCB104'


