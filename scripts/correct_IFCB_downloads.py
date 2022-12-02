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

# create list of directories with the data files:
dir_list = []
for l, i in enumerate(timeseries_data['Project_ID']):
    dataset = IFCB_dashboard_data_path + '/' + i  + '/'
    if timeseries_data.loc[l, 'dashboard_url'] == 'https://ifcb.caloos.org/':
        dashboard_id = 'CALOOS'
    elif timeseries_data.loc[l, 'dashboard_url'] == 'https://ifcb-data.whoi.edu/':
        dashboard_id = 'WHOI'
    try:
        file_list = os.listdir(dataset)
        for i in file_list:
            filename = dataset + i
            print(filename)
            df = pd.read_csv(filename, sep='\t', header=0, index_col=[0])
            df_clean = pd.DataFrame()
            df_clean['roi_number'] = df['roi_number']
            df_clean['Area'] = df['Area']
            df_clean['Biovolume'] = df['Biovolume']
            df_clean['EquivDiameter'] = df['EquivDiameter']
            df_clean['MajorAxisLength'] = df['MajorAxisLength']
            df_clean['MinorAxisLength'] = df['MinorAxisLength']
            df_clean['Solidity'] = df['Solidity']
            df_clean['summedArea'] = df['summedArea']
            df_clean['summedBiovolume'] = df['summedBiovolume']
            df_clean['summedMajorAxisLength'] = df['summedMajorAxisLength']
            df_clean['summedMinorAxisLength'] = df['summedMinorAxisLength']
            df_clean['Project_ID'] = df['Project_ID']
            df_clean['datetime'] = df['datetime']
            df_clean['pixel_to_micron'] = df['pixel_to_micron']
            df_clean['Latitude'] = df['Longitude']
            df_clean['depth'] = df['depth']
            df_clean['vol_analyzed'] = df['vol_analyzed']
            df_clean['sample_type'] = df['sample_type']
            df_clean['sample_cruise'] = df['sample_cruise']
            df_clean['number_of_rois'] = df['number_of_rois']
            df_clean['concentration'] = df['concentration']
            df_clean['roi_id'] = df['roi_id']
            df_clean['class_1'] = df['class_1']
            df_clean['class_1_score'] = df['class_1_score']
            df_clean['class_2'] = df['class_2']
            df_clean['class_2_score'] = df['class_2_score']
            df_clean['class_3'] = df['class_3']
            df_clean['class_3_score'] = df['class_3_score']
            df_clean['class_4'] = df['class_4']
            df_clean['class_4_score'] = df['class_4_score']
            df_clean['class_5'] = df['class_5']
            df_clean['class_5_score'] = df['class_5_score']
            # replace class and class  scores of 0 to Nan
            cols = ["class_1", "class_2", "class_3", "class_4", "class_5"]
            cols_scores = ["class_1_score", "class_2_score", "class_3_score", "class_4_score", "class_5_score"]
            for n, c in enumerate(cols_scores):
                df_clean[c] = df[c].mask(df[c] < 0.0001, np.nan)
                df_clean[cols[n]] = df[cols[n]].mask(df[c] < 0.0001, np.nan)
            os.remove(filename)

            # CONTINUE HERE, include in name dashboard info,  dataset ID, dataset name, year, 'features' and file number  :
            df_clean.to_csv(path_download + '/' + Project_ID + '_features_' + str(file_numbers) + '.tsv',
                                   sep='\t')



    except:
        pass






