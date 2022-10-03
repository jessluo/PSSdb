import datetime as dt
import requests
import urllib.request
import pandas as pd
import json
import os
import shutil
import numpy as np
import yaml
from pathlib import Path
# set paths, dates for query, and dataset name:
path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
# open the metadata of the standardized files
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
IFCB_dashboard_data_path = str(Path(cfg['IFCB_dir']).expanduser())

y = {}
y['startdate']=dt.datetime(2015,8,21);
y['enddate']=dt.datetime(2015,9,10);
y['size_cutoff']=50;
dataset_name = 'santa-cruz-municipal-wharf'
y['path_download']= IFCB_dashboard_data_path + '/' + dataset_name ## NOTE: the first part of the directory needs to be changed for the GIT PSSdb project
pathname = 'https://ifcb.caloos.org/' + dataset_name + '/'

if os.path.exists(y['path_download']) and os.path.isdir(y['path_download']):
    shutil.rmtree(y['path_download'])
os.mkdir(y['path_download'])

## create list of dates available for the project of interest
link = 'https://ifcb.caloos.org/api/list_bins?dataset=' + dataset_name # this path is for the data LIST only
r = requests.get(link)
r_content = r.content
r_content = json.loads(r_content.decode('utf-8'))
all_file_info = pd.DataFrame.from_dict(r_content['data'], orient = 'columns')
all_file_info['date_info'] = [int(sample[1:9]) for sample in all_file_info['pid']]
# constrain date ranges for file download
start_stamp=int(y['startdate'].strftime('%Y%m%d'))
end_stamp=int(y['enddate'].strftime('%Y%m%d'))
subset_file_info = all_file_info.loc[all_file_info['date_info']>= start_stamp].reset_index(drop =True) # | all_file_info['date_info']<= end_stamp]
subset_file_info = all_file_info.loc[all_file_info['date_info']<= end_stamp].reset_index(drop =True)

# use the file list to download data and store it locally
pathname = 'https://ifcb.caloos.org/' + dataset_name + '/' # notice the different url for the data
for i in range (len(subset_file_info)):
    localpath = y['path_download'] + '/' + str(subset_file_info.loc[i, 'pid'])# create local directory to store the files
    if os.path.exists(localpath) and os.path.isdir(localpath):
        shutil.rmtree(localpath)
    os.mkdir(localpath)

    filename = pathname + str(subset_file_info.loc[i, 'pid']) + '_features.csv'
    features = requests.get(filename)
    features = features.content
    features = features.decode('utf-8')
    features = features.split('\n')
    features_df = pd.DataFrame(features)
    features_dict = {}
    for n, i in enumerate(features):
        roi_info=i.split(',')
        features_dict[n] = roi_info
    features_dict.popitem()
    features_df = pd.DataFrame.from_dict(features_dict, orient='index')
    features_df.columns = features_df.iloc[0]
    features_df.drop(features_df.index[0], inplace = True)
    features_df.to_csv(localpath + '/' +str(subset_file_info.loc[0, 'pid']) + '_features.csv', sep='\t')

