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
        print(dataset)
        try:
            for num, id in enumerate(file_list):
                filename = dataset + id
                df = pd.read_csv(filename, sep='\t', header=0, index_col=[0])
                df_clean = pd.DataFrame()
                #df_clean['roi_number'] = df['roi_number']
                df_clean['area'] = df['Area']
                df_clean['biovolume'] = df['Biovolume']
                df_clean['equiv_diameter'] = df['EquivDiameter']
                df_clean['major_axis_length'] = df['MajorAxisLength']
                df_clean['minor_axis_length'] = df['MinorAxisLength']
                df_clean['solidity'] = df['Solidity']
                df_clean['summed_area'] = df['summedArea']
                df_clean['summed_biovolume'] = df['summedBiovolume']
                df_clean['summed_major_axis_length'] = df['summedMajorAxisLength']
                df_clean['summed_minor_axis_length'] = df['summedMinorAxisLength']
                df_clean['project_ID'] = df['Project_ID']
                df_clean['datetime'] = df['datetime']
                df_clean['pixel_to_micron'] = df['pixel_to_micron']
                df_clean['latitude'] = df['Latitude']
                df_clean['longitude'] = df['Longitude']
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
                df['datetime'] = pd.to_datetime(df['datetime'])
                year = df.loc[0, 'datetime'].strftime('%Y')
                dataset_id = df_clean.loc[0, 'roi_id'].split('_')
                dataset_id = dataset_id[1]
                os.remove(filename)
                #CONTINUE HERE, include in name dashboard info,  dataset ID, dataset name, year, 'features' and file number  :
                print(dataset + dashboard_id +'_'+ dataset_id +'_'+ df_clean.loc[0, 'project_ID'] +'_'+ year + '_features_' + str(num) + '.tsv')
                df_clean.to_csv(dataset + dashboard_id +'_'+ dataset_id +'_'+ df_clean.loc[0, 'project_ID'] +'_'+ year + '_features_' + str(num) + '.tsv', sep='\t')
        except:
            pass
    except:
        pass






