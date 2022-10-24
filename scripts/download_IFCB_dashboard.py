import datetime as dt
import re
import requests
import urllib.request
import pandas as pd
import json
import os
import shutil
import numpy as np
import yaml
from pathlib import Path
import re

# set paths:
path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
# open the metadata of the standardized files
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
IFCB_dashboard_data_path = str(Path(cfg['IFCB_dir']).expanduser())


# Set names of variables (date range, dataset, dashboard name)
y = {}
y['startdate'] = dt.datetime(2015,8,21);
y['enddate'] = dt.datetime(2015,9,10);
y['size_cutoff'] = 50;
Project_ID = 'EXPORTS'
Project_source = 'https://ifcb-data.whoi.edu/'
BASE_URL=Project_source # needs to be defined based on Karina's script
METRIC = 'temperature' # needs to be defined based on Karina's script
y['path_download'] = IFCB_dashboard_data_path + '/' + Project_ID ## NOTE: the first part of the directory needs to be changed for the GIT PSSdb project
pathname = Project_source + Project_ID  + '/'

if os.path.exists(y['path_download']) and os.path.isdir(y['path_download']):
    shutil.rmtree(y['path_download'])
os.mkdir(y['path_download'])

## FUNCTIONS HERE

## create list of dates available for the project of interest
def get_df_list_IFCB ( base_url, dataset, startdate=20000101,enddate=21000101):
    """
    Objective: generate list of files for a dataset in the IFCB dashboard
    :param startdate:
    :param enddate:
    :param base_url:
    :param dataset:
    :return:
    """
    link = base_url + 'api/list_bins?dataset=' + dataset # this path is for the data LIST only
    #metadata_link = "https://ifcb.caloos.org/timeline?dataset="
    r = requests.get(link)
    r_content = r.content
    r_content = json.loads(r_content.decode('utf-8'))
    all_file_info = pd.DataFrame.from_dict(r_content['data'], orient = 'columns')
    all_file_info['date_info'] = [int(sample[1:9]) for sample in all_file_info['pid']]
    # constrain date ranges for file download.
    #start_stamp=int(startdate.strftime('%Y%m%d'))
    #end_stamp=int(y['enddate'].strftime('%Y%m%d'))
    file_info = all_file_info.loc[all_file_info['date_info']>= startdate].reset_index(drop =True) # | all_file_info['date_info']<= end_stamp]
    file_info = file_info.loc[file_info['date_info']<= enddate].reset_index(drop =True)
    return file_info

# We adopted some code contributed by Joe Futrelle (WHOI): https://github.com/joefutrelle/pyifcb

# make function to get datasets from urls
def df_from_url(url):
    """
    Objective: use requests to build dataframes with ROI information from contents obtained from IFCB dashboards
    :param url: url for the dataset
    :return: a dataframe with the desired information (see outputs of IFCB dashboards)
    """
    data = requests.get(url)
    data = data.content
    data = data.decode('utf-8')
    data = data.split('\n')
    data_dict = {}
    for n, r in enumerate(data):
        roi_info = r.split(',')
        data_dict[n] = roi_info
    data_dict.popitem()
    df = pd.DataFrame.from_dict(data_dict, orient='index')
    df.columns = df.iloc[0]
    df.drop(df.index[0], inplace = True)
    return df


def metadata_dict(pid_id):
    """
    :param pid: the identifier for the data in each time bin displayed on the timeline in the dashboard
    :return: dictionary with the metadata
    """
    url = f'{BASE_URL}api/bin/{pid_id}?include_coordinates=False'
    r = requests.get(url)
    assert r.ok
    record = r.json()
    return {
        'datetime': record['timestamp_iso'],
        'previous_bin': record['previous_bin_id'],
        'next_bin': record['next_bin_id'],
        'latitude': record['lat'],
        'longitude': record['lng'],
        'depth': record['depth'],
        'instrument': record['instrument'],
        'number_of_rois': record['num_images'],
        'ml_analyzed': float(re.sub(' .*', '', record['ml_analyzed'])),
        'concentration': record['concentration'],
        'sample_type': record['sample_type'],
        'cruise': record['cruise'],
        'cast': record['cast'],
        'niskin': record['niskin'],
    }

#PROCESSING HERE



# use the file list to download data and store it locally
file_info = get_df_list_IFCB (startdate=20180811,enddate=20180812, base_url = Project_source, dataset = Project_ID)
for i in range (len(file_info)):
    try:
        # generate features files
        features_filename = pathname + str(file_info.loc[i, 'pid']) + '_features.csv'
        features_df = df_from_url(features_filename)
        features_df.head()
        localpath = y['path_download'] + '/' + str(file_info.loc[i, 'pid'])# create local directory to store the files
        if os.path.exists(localpath) and os.path.isdir(localpath):
            shutil.rmtree(localpath)
        os.mkdir(localpath)
        # obtain metadata
        pid_id = str(file_info.loc[i, 'pid'])
        met_dict = metadata_dict(pid_id)
        for c in ['texture', 'Wedge', 'Ring', 'HOG', 'moment_invariant']:
            features_df = features_df[features_df.columns.drop(list(features_df.filter(regex=c)))]
        #url = metadata_link + dataset_name +'&bin='+ str(subset_file_info.loc[i, 'pid'])
        features_df['datetime'] = met_dict['datetime']
        features_df['Latitude'] = met_dict['latitude']
        features_df['Longitude'] = met_dict['longitude']
        features_df['depth'] = met_dict['depth']
        features_df['vol_analyzed'] = met_dict['ml_analyzed']
        features_df['sample_type'] = met_dict['sample_type']
        features_df['cruise'] = met_dict['cruise']
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
        print('there is no features or class_scores files for' + str(file_info.loc[i, 'pid']))
#test = 'D20150910T212605_IFCB104'