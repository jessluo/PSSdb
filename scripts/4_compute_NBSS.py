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

light_parsing = cfg['day_night']
depth_parsing = cfg['depth_binning']


bin_loc = cfg['st_increment_avg']
group_by= cfg['date_group_avg']
sensitivity = cfg['sensitivity']
# processing starts here
for instrument in [ 'UVP', 'IFCB']:
    #first, break apart datasets by big global grids, to avoid making one HUGE file of all gridded datasets per instrument
    #get paths to gridded files
    file_list = glob(str(Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir']) + '*/**/' + instrument +'*_temp_binned_*.csv', recursive=True) #generate path and project ID's but ONLY for parsed data
    grid_list = [re.search('N_(.+?).csv', file) for file in file_list]#group_gridded_files_func(instrument, already_gridded='Y')
    NBSSpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data'
    if not os.path.exists(NBSSpath):
        os.mkdir(NBSSpath)
    if sensitivity == True:
        if not os.path.exists((NBSSpath / 'Sensitivity_analysis')):
            os.mkdir((NBSSpath / 'Sensitivity_analysis'))

    currentMonth = str(datetime.datetime.now().month).rjust(2, '0')
    currentYear = str(datetime.datetime.now().year)
    #set the list of biovolume estimates when doing the sensitivity analysis
    if sensitivity == True:
        print ('NOTE: products for Sensitivity analyses will be generated, this process will take longer')
        biovol_list = ['Biovolume_area', 'Biovolume_ellipsoid', 'Biovolume_orig'] if instrument == 'IFCB' else ['Biovolume_area', 'Biovolume_ellipsoid']
    else:
        biovol_list = ['Biovolume_area']


    for biovol in biovol_list:
        print ('generating PSSdb products for ' +instrument+' based on ' + biovol)
        NBSS_binned_all = pd.DataFrame()  # NBSS dataset
        lin_fit_data = pd.DataFrame()
        for i in tqdm(grid_list):
            i = i.group(1)
            print ('grouping data and calculating NBSS for cell number ' + i)
            file_subset = [file for file in file_list if i in file]
            NBS_biovol_df= pd.concat(map((lambda path: (pd.read_csv(path))), file_subset)).reset_index(drop=True)
            NBS_biovol_df, df_bins = size_binning_func(NBS_biovol_df, biovol)  # create size bins
            NBS_biovol_df = NB_SS_func(NBS_biovol_df, df_bins, biovol_estimate=biovol, sensitivity=sensitivity)
            lin_fit = linear_fit_func(NBS_biovol_df)
            NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])
            lin_fit_data = pd.concat([lin_fit_data, lin_fit])
            print(lin_fit_data[lin_fit_data.columns[0]].count())
            print(NBSS_binned_all[NBSS_binned_all.columns[0]].count())
        NBSS_raw = NBSS_binned_all.filter(['date_bin', 'midLatBin', 'midLonBin', 'light_cond', 'size_class_mid', 'ECD_mean', 'NB', 'PSD','Min_obs_depth', 'Max_obs_depth'], axis=1)
        if light_parsing == True:
            NBSS_1a_full = NBSS_stats_func(NBSS_raw, light_parsing=True, bin_loc=bin_loc, group_by=group_by)
            lin_fit_1b_full = stats_linfit_func(lin_fit_data, light_parsing=True, bin_loc=bin_loc, group_by=group_by)
        else:
            NBSS_1a_full = NBSS_stats_func(NBSS_raw, bin_loc=bin_loc, group_by=group_by)
            lin_fit_1b_full= stats_linfit_func(lin_fit_data, bin_loc=bin_loc, group_by=group_by)

        NBSS_1a_full['month_int'] = NBSS_1a_full['month'].astype(int)
        NBSS_1a_full['year_int'] = NBSS_1a_full['year'].astype(int)
        NBSS_1a_full = NBSS_1a_full.sort_values(by=['year_int', 'month_int'])
        NBSS_1a_full = NBSS_1a_full.drop(['year_int', 'month_int'], axis=1)
        NBSS_1a_full = ocean_label_func(NBSS_1a_full, 'longitude', 'latitude')

        lin_fit_1b_full['month_int'] = lin_fit_1b_full['month'].astype(int)
        lin_fit_1b_full['year_int'] = lin_fit_1b_full['year'].astype(int)
        lin_fit_1b_full = lin_fit_1b_full.sort_values(by=['year_int', 'month_int'])
        lin_fit_1b_full = lin_fit_1b_full.drop(['year_int', 'month_int'], axis=1)
        lin_fit_1b_full = ocean_label_func(lin_fit_1b_full, 'longitude', 'latitude')
        biovol = biovol.replace('_', '-')
        if sensitivity == True:
            NBSS_binned_all.to_csv(str(NBSSpath) + '/Sensitivity_analysis/' + instrument + '_' + biovol + '_by-Size_all_var' + currentYear + '-' + currentMonth + '.csv',index=False)
            NBSS_1a_full.to_csv(str(NBSSpath) + '/Sensitivity_analysis/' + instrument + '_1a_' + biovol + '_by-Size_v' + currentYear + '-' + currentMonth + '.csv',index=False)
            lin_fit_1b_full.to_csv(str(NBSSpath) + '/Sensitivity_analysis/' + instrument + '_1b_' + biovol + '_NBSS-fit_v' + currentYear + '-' + currentMonth + '.csv',index=False)
        else:
            NBSS_binned_all.to_csv(str(NBSSpath) + '/' + instrument + '_Biovolume-by-Size_all_var' + currentYear + '-' + currentMonth + '.csv',index=False)

            NBSS_1a_full.to_csv(str(NBSSpath) + '/' + instrument + '_1a_Biovolume-by-Size_v' + currentYear + '-' + currentMonth + '.csv',index=False)
            lin_fit_1b_full.to_csv(str(NBSSpath) + '/' + instrument + '_1b_NBSS-fit_v' + currentYear + '-' + currentMonth + '.csv',index=False)


