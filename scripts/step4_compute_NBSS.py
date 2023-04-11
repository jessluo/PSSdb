# compute NBS from step3 gridded files

import numpy as np
import pandas as pd
import statistics as st
import os

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

instrument = input ('for which instrument do you want to calculate the size spectra? \n Enter IFCB, Zooscan or UVP ')
day_night = input('Would you like to group data by day/night? Enter Y/N \n' )
depth_binning  = input ('Would you like to bin the data by depth? \n Enter Y/N \n')
bin_loc = float(input('select the size of the spatial grid to group the data i.e for a 1x1 degree bin type 1 \n' ))
group_by= input('input how will the data be grouped by date \n  yyyymm for month and year \n yyyy for year \n ')
big_grid_status = input('have the gridded files already been parsed and labelled according to their large global grid? \n type Y or N \n' )
# processing starts here

#first, break apart datasets by big global grids, to avoid making one HUGE file of all gridded datasets per instrument
grid_list = group_gridded_files_func(instrument, already_gridded=big_grid_status)

#get paths to gridded files
file_list = proj_id_list_func(instrument, data_status ='gridded', big_grid = True)#generate path and project ID's but ONLY for parsed data

path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

NBSSpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data'
if not os.path.exists(NBSSpath):
    os.mkdir(NBSSpath)

dirpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data' / instrument


if os.path.isdir(dirpath) and len(os.listdir(dirpath)) != 0:  # and  os.path.exists(path_download)
    replace = input('There is already NBS data in ' + str(dirpath) + ' do you want to replace the files? \n Y/N')
    if replace == 'Y':
        print('Overwriting normalized biomass data file(s), please wait')
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
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
            df = df.filter(['Sample','date_bin', 'Station_location','light_cond', 'midDepthBin', 'Min_obs_depth', 'Max_obs_depth', 'ROI_number','Biovolume', 'midLatBin', 'midLonBin', 'Volume_imaged'], axis=1)
            if depth_binning == 'Y':
                NBSS_1a_raw, NBSS_1a,  lin_fit_1b = parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], light_parsing=True, depth_parsing= True, bin_loc = bin_loc, group_by = group_by)
            if day_night == 'Y':
                NBSS_1a_raw, NBSS_1a,  lin_fit_1b= parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], light_parsing=True, bin_loc = bin_loc, group_by = group_by)
            else:
                NBSS_1a_raw, NBSS_1a,  lin_fit_1b = parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], bin_loc = bin_loc, group_by = group_by)

            NBSS_1a_raw_full = pd.concat([NBSS_1a_raw_full, NBSS_1a_raw])
            NBSS_1a_full = pd.concat([NBSS_1a_full, NBSS_1a])
            lin_fit_1b_full = pd.concat([lin_fit_1b_full, lin_fit_1b])
        # Save NBSS results and sorting
        NBSS_1a_raw_full.to_csv(str(dirpath) + '/' + instrument +'_NBSS_raw.csv', index=False)

        NBSS_1a_full['month_int'] = NBSS_1a_full['month'].astype(int)
        NBSS_1a_full['year_int'] = NBSS_1a_full['year'].astype(int)
        NBSS_1a_full = NBSS_1a_full.sort_values(by=['year_int', 'month_int'])
        NBSS_1a_full = NBSS_1a_full.drop(['year_int', 'month_int'], axis=1)
        NBSS_1a_full.to_csv(str(dirpath) + '/' + instrument +'_NBSS_1a.csv', index=False)

        lin_fit_1b_full['month_int'] = lin_fit_1b_full['month'].astype(int)
        lin_fit_1b_full['year_int'] = lin_fit_1b_full['year'].astype(int)
        lin_fit_1b_full = lin_fit_1b_full.sort_values(by=['year_int', 'month_int'])
        lin_fit_1b_full = lin_fit_1b_full.drop(['year_int', 'month_int'], axis=1)
        lin_fit_1b_full.to_csv(str(dirpath) + '/' + instrument + '_lin_fit_1b.csv', index=False)

    elif replace == 'N':
        print('previous products 1a and 1b will be kept')

elif not os.path.exists(dirpath):
    os.mkdir(dirpath)
    NBSS_1a_raw_full = pd.DataFrame()
    NBSS_1a_full = pd.DataFrame()
    lin_fit_1b_full = pd.DataFrame()
    for i in tqdm(grid_list):
        print('grouping data and calculating NBSS for cell number ' + i)
        file_subset = [file for file in file_list if i in file]
        df = pd.concat(map((lambda path: (pd.read_csv(path))), file_subset)).reset_index(drop=True)
        df = df.dropna(subset=['Biovolume']).reset_index(drop=True)  # project 5785 from Zooscan is completely empty
        if len(df) == 0:
            print('no data left after removing empty biovolumes in grid cell ' + i)
            continue
        df = df.filter(
            ['Sample', 'date_bin', 'Station_location', 'light_cond', 'midDepthBin', 'Min_obs_depth', 'Max_obs_depth', 'ROI_number',
             'Biovolume', 'midLatBin', 'midLonBin', 'Volume_imaged'], axis=1)
        if depth_binning == 'Y':
            NBSS_1a_raw, NBSS_1a,  lin_fit_1b = parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'],
                                                        light_parsing=True, depth_parsing=True, bin_loc = bin_loc, group_by = group_by)
        if day_night == 'Y':
            NBSS_1a_raw, NBSS_1a,  lin_fit_1b = parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'],
                                                        light_parsing=True, bin_loc = bin_loc, group_by = group_by)
        else:
            NBSS_1a_raw, NBSS_1a,  lin_fit_1b = parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], bin_loc = bin_loc, group_by = group_by)

        NBSS_1a_raw_full = pd.concat([NBSS_1a_raw_full, NBSS_1a_raw])
        NBSS_1a_full = pd.concat([NBSS_1a_full, NBSS_1a])
        lin_fit_1b_full = pd.concat([lin_fit_1b_full, lin_fit_1b])
    # Save NBSS results
    NBSS_1a_raw_full.to_csv(str(dirpath) + '/' + instrument + '_NBSS_raw.csv', index=False)

    NBSS_1a_full['month_int'] = NBSS_1a_full['month'].astype(int)
    NBSS_1a_full['year_int'] = NBSS_1a_full['year'].astype(int)
    NBSS_1a_full = NBSS_1a_full.sort_values(by=['year_int', 'month_int'])
    NBSS_1a_full = NBSS_1a_full.drop(['year_int', 'month_int'], axis=1)
    NBSS_1a_full.to_csv(str(dirpath) + '/' + instrument + '_NBSS_1a.csv', index=False)

    lin_fit_1b_full['month_int'] = lin_fit_1b_full['month'].astype(int)
    lin_fit_1b_full['year_int'] = lin_fit_1b_full['year'].astype(int)
    lin_fit_1b_full = lin_fit_1b_full.sort_values(by=['year_int', 'month_int'])
    lin_fit_1b_full = lin_fit_1b_full.drop(['year_int', 'month_int'], axis=1)
    lin_fit_1b_full.to_csv(str(dirpath) + '/' + instrument + '_lin_fit_1b.csv', index=False)

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
