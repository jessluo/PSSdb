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
for instrument in ['Scanner', 'UVP', 'IFCB']:
    #first, break apart datasets by big global grids, to avoid making one HUGE file of all gridded datasets per instrument
    grid_list = [re.search('N_(.+?).csv', file) for file in file_list]#group_gridded_files_func(instrument, already_gridded='Y')

    #get paths to gridded files
    file_list = glob(str(Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir']) + '*/**/' + instrument +'*_temp_binned_*.csv', recursive=True) #generate path and project ID's but ONLY for parsed data

    NBSSpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data'
    if not os.path.exists(NBSSpath):
        os.mkdir(NBSSpath)
        if sensitivity == True:
            os.mkdir((NBSSpath / 'Sensitivity_analysis'))

    currentMonth = str(datetime.datetime.now().month).rjust(2, '0')
    currentYear = str(datetime.datetime.now().year)
    #set the list of biovolume estimates when doing the sensitivity analysis
    if sensitivity == True:
        print ('NOTE: products for Sensitivity analyses will be generated, this process will take longer')
        biovol_list = ['Biovolume_area', 'Biovolume_ellipsoid', 'Biovolume_orig'] if instrument == 'IFCB' else ['Biovolume_area', 'Biovolume_ellipsoid']
    else:
        biovol_list = ['Biovolume_area']

        NBSS_full_var_full = pd.DataFrame()
        NBSS_1a_raw_full = pd.DataFrame()
        NBSS_1a_full = pd.DataFrame()
        lin_fit_1b_full = pd.DataFrame()
        Sample_NB_ID = pd.Series()
        
    for biovol in biovol_list:
        print ('generating PSSdb products for ' +instrument+' based on ' + biovol)
        NBSS_binned_all = pd.DataFrame()  # NBSS dataset
        lin_fit_data = pd.DataFrame()
        for i in tqdm(file_list):
            print ('grouping data and calculating NBSS for cell number ' + re.search('N_(.+?).csv', i).group(1))
            df_binned, df_bins = size_binning_func(i, biovol)  # create size bins
            NBS_biovol_df, NBSS_1a_raw = NB_SS_func(df_binned, df_bins, biovol_estimate=biovol, sensitivity=sensitivity)
            lin_fit = linear_fit_func(NBS_biovol_df)
            NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])
            lin_fit_data = pd.concat([lin_fit_data, lin_fit])
 # CONTINUE REESTRUCTURING STEP 4 from here

        for i in tqdm(grid_list):
            #print('grouping data and calculating NBSS for cell number ' + i)
            file_subset = [file for file in file_list if i in file]
            df = pd.concat(map((lambda path: (pd.read_csv(path))), file_subset)).reset_index(drop=True)
            #rows 66-74 contain routine to detect duplicates in the data
            df['NB_Sample_ID'] = df.date_bin.astype(str) + df.Station_location.astype(str)
            contains_sample= pd.Series(df.NB_Sample_ID.unique()).isin(Sample_NB_ID)
            Sample_NB_ID = pd.concat([Sample_NB_ID, pd.Series(df.NB_Sample_ID.unique())])
            #if contains_sample.any()==True:
                #print ('Warning! during 15 degree gridding, duplicate datasets were generated')
                #print('grid' +i + 'files in more than one 15 degree cell')
                #print([f for f in file_subset])
                #print(Sample_NB_ID[Sample_NB_ID.duplicated()==True])
                #break
            #df = df.dropna(subset=['Biovolume']).reset_index(drop=True) # project 5785 from Zooscan is completely empty
            if len(df) == 0:
                print('no data left after removing empty biovolumes in grid cell' + i)
                continue
            #df = df.filter(['Sample','date_bin', 'Station_location','light_cond', 'midDepthBin', 'Sampling_lower_size', 'Sampling_upper_size', 'Min_obs_depth', 'Max_obs_depth', 'ROI_number','Biovolume', 'midLatBin', 'midLonBin', 'Volume_imaged'], axis=1)

            NBSS_full_var, NBSS_1a_raw, NBSS_1a,  lin_fit_1b = parse_NBS_linfit_func(df,  biovol, sensitivity, light_parsing, depth_parsing, bin_loc = bin_loc, group_by = group_by)


            NBSS_full_var_full = pd.concat([NBSS_full_var_full, NBSS_full_var])
            NBSS_1a_raw_full = pd.concat([NBSS_1a_raw_full, NBSS_1a_raw])
            NBSS_1a_full = pd.concat([NBSS_1a_full, NBSS_1a])
            lin_fit_1b_full = pd.concat([lin_fit_1b_full, lin_fit_1b])

        NBSS_full_var_full = ocean_label_func(NBSS_full_var_full, 'midLonBin', 'midLatBin')
        NBSS_1a_raw_full = ocean_label_func(NBSS_1a_raw_full, 'longitude', 'latitude')

        # Save NBSS results and sorting
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

        if sensitivity == True:
            NBSS_full_var_full.to_csv(str(NBSSpath) + '/Sensitivity_analysis/' + instrument + '_'+biovol+'-by-Size_all_var'+currentYear+'-'+currentMonth+'.csv', index=False)
            NBSS_1a_raw_full.to_csv(str(NBSSpath) + '/Sensitivity_analysis/' + instrument +'_'+biovol+'-by-Size_raw_v'+currentYear+'-'+currentMonth+'.csv', index=False)
            NBSS_1a_full.to_csv(str(NBSSpath) + '/Sensitivity_analysis/' + instrument +'_1a_'+biovol+'-by-Size_v' +currentYear+'-'+currentMonth+'.csv', index=False)
            lin_fit_1b_full.to_csv(str(NBSSpath) + '/Sensitivity_analysis/' + instrument + '_1b_'+biovol+'-NBSS-fit_v' +currentYear+'-'+currentMonth+'.csv', index=False)
        else:
            NBSS_full_var_full.to_csv(str(NBSSpath) + '/' + instrument + '_Biovolume-by-Size_all_var' + currentYear + '-' + currentMonth + '.csv',index=False)
            NBSS_1a_raw_full.to_csv(str(NBSSpath) + '/' + instrument + '_Biovolume-by-Size_raw_v' + currentYear + '-' + currentMonth + '.csv',index=False)
            NBSS_1a_full.to_csv(str(NBSSpath) + '/' + instrument + '_1a_Biovolume-by-Size_v' + currentYear + '-' + currentMonth + '.csv',index=False)
            lin_fit_1b_full.to_csv(str(NBSSpath) + '/' + instrument + '_1b_NBSS-fit_v' + currentYear + '-' + currentMonth + '.csv',index=False)
