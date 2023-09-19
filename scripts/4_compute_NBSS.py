# compute NBS from step3 gridded files
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import statistics as st
import os
from datetime import datetime
from pathlib import Path
from glob import glob
import shutil
# Config modules
import yaml  # requires installation of PyYAML package
try:
    from funcs_read import *
    from funcs_NBS import *
except:
    from scripts.funcs_read import *
    from scripts.funcs_NBS import *

from tqdm import tqdm
import matplotlib.pyplot as plt



#open config file and extract parameters:
path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

day_night = cfg['day_night']
depth_binning = cfg['depth_binning']


bin_loc = cfg['st_increment_avg']
group_by= cfg['date_group_avg']
big_grid_status = cfg['big_grid']
# processing starts here
for instrument in ['Scanner', 'UVP', 'IFCB']:
    #first, break apart datasets by big global grids, to avoid making one HUGE file of all gridded datasets per instrument
    grid_list = group_gridded_files_func(instrument, already_gridded='Y')

    #get paths to gridded files
    file_list = glob(str(Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir']) + '*/**/' + instrument +'*_temp_binned_*.csv', recursive=True) #generate path and project ID's but ONLY for parsed data

    NBSSpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data'
    if not os.path.exists(NBSSpath):
        os.mkdir(NBSSpath)

    currentMonth = str(datetime.now().month).rjust(2, '0')
    currentYear = str(datetime.now().year)

    NBSS_full_var_full = pd.DataFrame()
    NBSS_1a_raw_full = pd.DataFrame()
    NBSS_1a_full = pd.DataFrame()
    lin_fit_1b_full = pd.DataFrame()
    Sample_NB_ID = pd.Series()
    for i in tqdm(grid_list):
        #print('grouping data and calculating NBSS for cell number ' + i)
        file_subset = [file for file in file_list if i  in file]
        df = pd.concat(map((lambda path: (pd.read_csv(path))), file_subset)).reset_index(drop=True)
        #rows 66-74 contain routine to detect duplicates in the data
        df['NB_Sample_ID'] = df.date_bin.astype(str) + df.Station_location.astype(str)
        contains_sample= pd.Series(df.NB_Sample_ID.unique()).isin(Sample_NB_ID)
        Sample_NB_ID = pd.concat([Sample_NB_ID, pd.Series(df.NB_Sample_ID.unique())])
        if contains_sample.any()==True:
            print ('Warning! during 15 degree gridding, duplicate datasets were generated')
            print('grid' +i + 'files in more than one 15 degree cell')
            print([f for f in file_subset])
            print(Sample_NB_ID[Sample_NB_ID.duplicated()==True])
            break
        df = df.dropna(subset=['Biovolume']).reset_index(drop=True) # project 5785 from Zooscan is completely empty
        if len(df) == 0:
            print('no data left after removing empty biovolumes in grid cell' + i)
            continue
        #df = df.filter(['Sample','date_bin', 'Station_location','light_cond', 'midDepthBin', 'Sampling_lower_size', 'Sampling_upper_size', 'Min_obs_depth', 'Max_obs_depth', 'ROI_number','Biovolume', 'midLatBin', 'midLonBin', 'Volume_imaged'], axis=1)
        if depth_binning == 'Y':
            NBSS_full_var, NBSS_1a_raw, NBSS_1a,  lin_fit_1b = parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], light_parsing=True, depth_parsing= True, bin_loc = bin_loc, group_by = group_by)
        if day_night == 'Y':
            NBSS_full_var, NBSS_1a_raw, NBSS_1a,  lin_fit_1b= parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], light_parsing=True, bin_loc = bin_loc, group_by = group_by)
        else:
            NBSS_full_var, NBSS_1a_raw, NBSS_1a,  lin_fit_1b = parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], bin_loc = bin_loc, group_by = group_by)

        NBSS_full_var_full = pd.concat([NBSS_full_var_full, NBSS_full_var])
        NBSS_1a_raw_full = pd.concat([NBSS_1a_raw_full, NBSS_1a_raw])
        NBSS_1a_full = pd.concat([NBSS_1a_full, NBSS_1a])
        lin_fit_1b_full = pd.concat([lin_fit_1b_full, lin_fit_1b])
    # Save NBSS results and sorting

    NBSS_full_var_full.to_csv(str(NBSSpath) + '/' + instrument +'_Biovolume-by-Size_all_var'+currentYear+'-'+currentMonth+'.csv', index=False)

    NBSS_1a_raw_full.to_csv(str(NBSSpath) + '/' + instrument +'_Biovolume-by-Size_raw_v'+currentYear+'-'+currentMonth+'.csv', index=False)

    NBSS_1a_full['month_int'] = NBSS_1a_full['month'].astype(int)
    NBSS_1a_full['year_int'] = NBSS_1a_full['year'].astype(int)
    NBSS_1a_full = NBSS_1a_full.sort_values(by=['year_int', 'month_int'])
    NBSS_1a_full = NBSS_1a_full.drop(['year_int', 'month_int'], axis=1)
    NBSS_1a_full.to_csv(str(NBSSpath) + '/' + instrument +'_1a_Biovolume-by-Size_v' +currentYear+'-'+currentMonth+'.csv', index=False)

    lin_fit_1b_full['month_int'] = lin_fit_1b_full['month'].astype(int)
    lin_fit_1b_full['year_int'] = lin_fit_1b_full['year'].astype(int)
    lin_fit_1b_full = lin_fit_1b_full.sort_values(by=['year_int', 'month_int'])
    lin_fit_1b_full = lin_fit_1b_full.drop(['year_int', 'month_int'], axis=1)
    lin_fit_1b_full.to_csv(str(NBSSpath) + '/' + instrument + '_1b_NBSS-fit_v' +currentYear+'-'+currentMonth+'.csv', index=False)
