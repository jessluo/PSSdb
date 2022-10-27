# Process to download raw data from IFCB dashboards

import pandas as pd
import os
import shutil
import numpy as np
import yaml
import datetime as dt
from pathlib import Path
from funcs_IFCB_dashboard import *

# set paths:
path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
# open the metadata of the standardized files
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
IFCB_dashboard_data_path = str(Path(cfg['IFCB_dir']).expanduser())


# Set names of variables (date range, dataset, dashboard name) NOTE: Project_ID and Project_source should be obtained from the
#IFCB standardizer
which_dashboard = input('from which IFCB dashboard do you want to download (WHOI or CALOOS)?')
if which_dashboard == 'WHOI':
    Project_source = 'https://ifcb-data.whoi.edu/'
elif which_dashboard == 'CALOOS':
    Project_source = 'https://ifcb.caloos.org/'

Project_ID = input('what is the name of the project you want to download? name has to be as in the dashboard')

start = int(input('what is the start date you want to set? ( format: YYYYMMDD)'))

end = int(input('what is the end date you want to set? ( format: YYYYMMDD)'))

#define path for download
path_download = IFCB_dashboard_data_path + '/' + Project_ID ## NOTE: the first part of the directory needs to be changed for the GIT PSSdb project
pathname = Project_source + Project_ID + '/'

if os.path.exists(path_download) and os.path.isdir(path_download):
    shutil.rmtree(path_download)
os.mkdir(path_download)

#PROCESSING HERE
METRIC = 'temperature' # needs to be defined based on Karina's script
# use the file list to download data and store it locally
file_info = get_df_list_IFCB(startdate= start ,enddate= end , base_url = Project_source, dataset = Project_ID)

#the loop uses the file_info dataframe to get the 'pid' (file name) and download features and class scores files
for i in range (len(file_info)):
    try:
        # generate features files
        pid_id = str(file_info.loc[i, 'pid'])
        features_filename = pathname + pid_id + '_features.csv'
        features_df = df_from_url(features_filename)
        #features_df.head()
        localpath = path_download + '/' + pid_id# create local directory to store the files
        if os.path.exists(localpath) and os.path.isdir(localpath):
            shutil.rmtree(localpath)
        os.mkdir(localpath)
        # obtain metadata
        met_dict = metadata_dict(dashboard= Project_source, pid_id=pid_id)
        for c in ['texture', 'Wedge', 'Ring', 'HOG', 'moment_invariant']:
            features_df = features_df[features_df.columns.drop(list(features_df.filter(regex=c)))]
        #url = metadata_link + dataset_name +'&bin='+ str(subset_file_info.loc[i, 'pid'])
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
        #now generate dataframe for the class scores
        class_filename = pathname + str(file_info.loc[i, 'pid']) + '_class_scores.csv'
        class_df = df_from_url(class_filename)
        #class_df = class_df.drop(columns = 'pid')
        # create column with taxonomic category based on class score
        for r in class_df.index.values:
            #class_df.iloc[r, 2:len(class_df.columns)]= class_df.iloc[r, 2:len(class_df.columns)].astype(float)
            scores = class_df.iloc[r-1, 2:len(class_df.columns)].astype(float)
            for col in class_df.columns:
                if class_df.loc[r, col] == str(max(scores)):
                    features_df.loc[r, 'class'] = col

            ## continue trying to get the class for each ROI


        features_df.to_csv(localpath + '/' +str(file_info.loc[i, 'pid']) + '_features.csv', sep='\t')
        class_df.to_csv(localpath + '/' +str(file_info.loc[i, 'pid']) + '_class_scores.csv', sep='\t')
        print(str(file_info.loc[i, 'pid']) + ' download done ')

    except:
        print('there is no features or class_scores files for ' + str(file_info.loc[i, 'pid']))
#test = 'D20150910T212605_IFCB104'