# compute NBS from step3 gridded files

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
#look for coastline layer and add it to that map, UNABLE to import cartopy, see issue
import cartopy.crs as ccrs


#define the instrument to calculate NBS:

path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

instrument = cfg['instrument_data_source']
day_night = cfg['day_night']
depth_binning = cfg['depth_binning']


bin_loc = cfg['st_increment_avg']
group_by= cfg['date_group_avg']
big_grid_status = cfg['big_grid']
# processing starts here

#first, break apart datasets by big global grids, to avoid making one HUGE file of all gridded datasets per instrument
grid_list = group_gridded_files_func(instrument, already_gridded=big_grid_status)

#get paths to gridded files
file_list = proj_id_list_func(instrument, data_status ='gridded', big_grid = True)#generate path and project ID's but ONLY for parsed data



NBSSpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data'
if not os.path.exists(NBSSpath):
    os.mkdir(NBSSpath)

currentMonth = str(datetime.now().month).rjust(2, '0')
currentYear = str(datetime.now().year)

NBSS_1a_binned_full = pd.DataFrame()
NBSS_1a_raw_full = pd.DataFrame()
NBSS_1a_full = pd.DataFrame()
lin_fit_1b_full = pd.DataFrame()
for i in tqdm(grid_list):
    print('grouping data and calculating NBSS for cell number ' + i)
    file_subset = [file for file in file_list if i in file]
    df = pd.concat(map((lambda path: (pd.read_csv(path))), file_subset)).reset_index(drop=True)
    df = df.dropna(subset=['Biovolume']).reset_index(drop=True) # project 5785 from Zooscan is completely empty
    if len(df) == 0:
        print('no data left after removing empty biovolumes in grid cell' + i)
        continue
     #df = df.filter(['Sample','date_bin', 'Station_location','light_cond', 'midDepthBin', 'Sampling_lower_size', 'Sampling_upper_size', 'Min_obs_depth', 'Max_obs_depth', 'ROI_number','Biovolume', 'midLatBin', 'midLonBin', 'Volume_imaged'], axis=1)
    if depth_binning == 'Y':
        NBSS_binned_all, NBSS_1a_raw, NBSS_1a,  lin_fit_1b = parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], light_parsing=True, depth_parsing= True, bin_loc = bin_loc, group_by = group_by)
    if day_night == 'Y':
        NBSS_binned_all, NBSS_1a_raw, NBSS_1a,  lin_fit_1b= parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], light_parsing=True, bin_loc = bin_loc, group_by = group_by)
    else:
        NBSS_binned_all, NBSS_1a_raw, NBSS_1a,  lin_fit_1b = parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], bin_loc = bin_loc, group_by = group_by)

    NBSS_1a_binned_full= pd.concat([NBSS_1a_binned_full, NBSS_binned_all])
    NBSS_1a_raw_full = pd.concat([NBSS_1a_raw_full, NBSS_1a_raw])
    NBSS_1a_full = pd.concat([NBSS_1a_full, NBSS_1a])
    lin_fit_1b_full = pd.concat([lin_fit_1b_full, lin_fit_1b])
# Save NBSS results and sorting
NBSS_1a_binned_full.to_csv(str(NBSSpath) + '/' + instrument +'_Biovolume-by-Size_FULL_var_v'+currentYear+'-'+currentMonth+'.csv', index=False)
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



#lat = lin_fit_1b['latitude']
#lon = lin_fit_1b['longitude']
#slope = lin_fit_1b['slope_mean']
#intercept_t = lin_fit_1b['intercept_mean']

#intercept_plot = [x*3 for x in intercept_t]
#ax = plt.axes(projection=ccrs.PlateCarree())
#plt.gca().coastlines('50m')
#g1=ax.gridlines(draw_labels=True)
#g1.xlines = False
#g1.ylines = False
#plt.scatter(lon, lat, label=None, c=slope, cmap='viridis', s=intercept_plot, linewidth=0, alpha=0.5, transform=ccrs.PlateCarree())
#plt.gca().coastlines('50m')
#ax.set_extent([-180, 180, -90, 90])
#ax.xlabel('longitude')
#ax.ylabel('latitude')
#plt.colorbar(label='slope', orientation='horizontal', anchor=(0.5, 1))
#ax.clim(min(slope), max(slope))


#labels = [str(np.round(min(intercept_t), decimals=2)),str(np.round(st.median(intercept_t), decimals=2)), str(np.round(max(intercept_t), decimals=2))]

#for n, area in enumerate([18, 28.5, 73.5]):
    #plt.scatter([], [], c='k', alpha=0.3, s=area, label=labels[n], transform=ccrs.PlateCarree())

#plt.legend(bbox_to_anchor=(0.75, 0), ncol = 3, scatterpoints=1, frameon=False,
           #labelspacing=1, title='intercept')

#figname = 'step4_slopes_intercept_updated_firstrelease' + instrument + '.pdf'

#path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
#with open(path_to_config, 'r') as config_file:
    #cfg = yaml.safe_load(config_file)

#savepath = Path(cfg['git_dir']).expanduser() / 'figures' / figname

#plt.savefig(savepath)
