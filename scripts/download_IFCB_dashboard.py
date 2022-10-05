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
import re
# Module for web crawling with javascript inputs
# Selenium for web scraping: https://stackoverflow.com/questions/29385156/invoking-onclick-event-with-beautifulsoup-python
# Download chromedriver at:https://chromedriver.chromium.org/downloads
# Go in System Preferences > Security and Privacy > General > Allow
from selenium import webdriver # Use: pip3 install -U selenium
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager # Use: pip3 install webdriver-manager
from selenium.webdriver.chrome.options import Options
import warnings
warnings.filterwarnings("ignore")

import time

chromedriver = '{}/Downloads/chromedriver'.format(str(Path.home())) # Local path to the chromedriver. Note: Deprecated
options=Options()
options.add_argument("--headless") # Use this to keep the browser session closed

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


## create list of dates available for the project of interest
link = 'https://ifcb.caloos.org/api/list_bins?dataset=' + dataset_name # this path is for the data LIST only
metadata_link = "https://ifcb.caloos.org/timeline?dataset="
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
    #generate features files
    features_filename = pathname + str(subset_file_info.loc[i, 'pid']) + '_features.csv'
    features_df = df_from_url(features_filename)
    for c in ['texture', 'Wedge', 'Ring', 'HOG']:
        features_df = features_df[features_df.columns.drop(list(features_df.filter(regex=c)))]
    url = metadata_link + dataset_name +'&bin='+ str(subset_file_info.loc[i, 'pid'])
    with webdriver.Chrome(chromedriver, options=options, keep_alive=False) as driver:
        driver.get(url)
        time.sleep(60)  # Time needed for the webpage to be fully loaded, in seconds
        vol_analyzed = float(re.sub("[^\d\.]", "", driver.find_element(by='id', value='stat-ml-analyzed').text)) if len(driver.find_element(by='id', value='stat-ml-analyzed').text) > 0 else pd.NA
        depth = float(re.sub("[^\d\.]", "", driver.find_element(by='id', value='stat-depth').text)) if len(driver.find_element(by='id', value='stat-depth').text) > 0 else pd.NA
        Latitude = float(driver.find_element(by='id', value='stat-lat').text) if len(driver.find_element(by='id', value='stat-lat').text) > 0 else pd.NA
        Longitude = float(driver.find_element(by='id', value='stat-lon').text) if len(driver.find_element(by='id', value='stat-lon').text) > 0 else pd.NA
        driver.quit()
    features_df['Latitude'] = Latitude
    features_df['Longitude'] = Longitude
    features_df['depth'] = depth
    features_df['vol_analyzed'] = vol_analyzed

    #now generate dataframe for the class scores
    class_filename = pathname + str(subset_file_info.loc[i, 'pid']) + '_class_scores.csv'
    class_df = df_from_url(class_filename)

    class_df.iloc[1, 2:len(class_df.columns)]= class_df.iloc[1, 2:len(class_df.columns)].astype(float)
    scores = class_df.iloc[1, 2:len(class_df.columns)].astype(float)
    for col in class_df.columns:
        if class_df.loc[1, col] == max (class_df.iloc[1, 2:len(class_df.columns)]):
            print('yes')
        else:
            print('no')

            ## continue trying to get the class for each ROI




    df.to_csv(localpath + '/' +str(subset_file_info.loc[0, 'pid']) + t, sep='\t')

