
import numpy as np
import pandas as pd
import statistics as st
import os

from pathlib import Path
from glob import glob
import shutil
# Config modules
import yaml  # requires installation of PyYAML package
from funcs_read import *
from funcs_gridding import *
from funcs_NBS import *
from tqdm import tqdm
import matplotlib.pyplot as plt
#look for coastline layer and add it to that map, UNABLE to import cartopy, see issue
import cartopy.crs as ccrs


start_time = input('type start time to create the dummy dataset example 2009-01-01 00:00:00 \n')
end_time = input('type end time to create the dummy dataset example 2009-01-01 00:00:00 \n')

timestamp_df = pd.date_range(start=pd.Timestamp(start_time), end=pd.Timestamp(end_time), freq='15T')

date_df = timestamp_df.strftime('%Y-%m-%d')
date_df = [date.replace('-', '') for date in date_df]

time_df = timestamp_df.strftime('%H-%M-%S')
time_df = [time.replace('-', '') for time in time_df]

full_df = pd.DataFrame(data = {'Datetime': timestamp_df, 'Sampling_date': date_df, 'Sampling_time': time_df, 'Latitude': 37.80305 , 'Longitude':-122.39723, 'Volume_imaged': 5,
                               'Min_obs_depth':0, 'Max_obs_depth': 5 })

full_df = pd.concat([full_df]*30).reset_index(drop=True)
full_df = full_df.sort_values(by = 'Datetime').reset_index(drop=True)
full_df = date_binning_func(full_df, group_by= 'week',day_night=False, ignore_high_lat=True)
full_df['Station_location'],full_df['midLatBin'], full_df['midLonBin'] = gridding_func(0.5, lat= full_df['Latitude'], lon= full_df['Longitude'])
#assign sample id
sample_id = list(range(0, len(timestamp_df), 1))
sample_id = sample_id*30
sample_id.sort()
full_df['Sample'] = sample_id
full_df['Sample'] = 'id_' + full_df['Sample'].astype(str)

#assign biovolume
full_df['Biovolume'] = np.random.exponential(50000, size=len(full_df))
full_df = full_df[full_df['Biovolume'] > 0.522]  # had to do this since some values are below our size classes


path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)


save_path = Path(cfg['raw_dir']).expanduser() / 'dummy_data'
if not os.path.exists(save_path):
    os.mkdir(save_path)

full_df.to_csv(str(save_path) + '/dummy_data_timeseries.csv', index=False)


NBSS_1a, lin_fit_1b = parse_NBS_linfit_func(full_df, parse_by=['Station_location', 'date_bin'], light_parsing=False, depth_parsing=False, bin_loc = 1, group_by = 'yyyymm')


NBSS_1a['month_int'] = NBSS_1a['month'].astype(int)
NBSS_1a['year_int'] = NBSS_1a['year'].astype(int)
NBSS_1a= NBSS_1a.sort_values(by = ['year_int', 'month_int'])
NBSS_1a = NBSS_1a.drop(['year_int', 'month_int'], axis=1)

lin_fit_1b['month_int'] = lin_fit_1b['month'].astype(int)
lin_fit_1b['year_int'] = lin_fit_1b['year'].astype(int)
lin_fit_1b = lin_fit_1b.sort_values(by=['year_int', 'month_int'])
lin_fit_1b = lin_fit_1b.drop(['year_int', 'month_int'], axis=1)

NBSS_1a.to_csv(str(save_path) + '/NBSS_1a.csv', index=False)
lin_fit_1b.to_csv(str(save_path) + '/lin_fit_1b.csv', index=False)






