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
from funcs_read import *
from funcs_NBS import *
from tqdm import tqdm
import matplotlib.pyplot as plt
#look for coastline layer and add it to that map, UNABLE to import cartopy, see issue
import cartopy.crs as ccrs


#define the instrument to calculate NBS:

instrument = input ('for which instrument do you want to calculate the size spectra? \n Enter IFCB, Zooscan or UVP ')
depth_binning  = input ('Would you like to bin the data by depth? \n Enter Y/N')

# processing starts here

file_list = proj_id_list_func(instrument, data_status ='gridded')#generate path and project ID's

path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

NBSSpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data'
if not os.path.exists(NBSSpath):
    os.mkdir(NBSSpath)

dirpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data' / instrument


if os.path.isdir(dirpath) and len(os.listdir(dirpath)) != 0:  # and  os.path.exists(path_download)
    replace = input('There is already NBS data in ' + dirpath + ' do you want to replace the files? \n Y/N')
    if replace == 'Y':
        print('Overwriting normalized biomass data file(s), please wait')
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
        df = pd.concat(map((lambda path: (pd.read_csv(path))), file_list)).reset_index(drop=True)
        df = df.filter(['Sample','date_bin', 'Station_location', 'midDepthBin', 'Min_obs_depth', 'Max_obs_depth', 'range_size_bin', 'sizeClasses', 'Biovolume', 'midLatBin', 'midLonBin', 'Volume_imaged'], axis=1)

        if depth_binning == 'Y':
            NBSS_1a, lin_fit_1b = parse_NBS_linfit_func(df, parse_by=['Station_location', 'date_bin'],
                                                                  depth_bin=True)
        else:
            NBSS_1a, lin_fit_1b = parse_NBS_linfit_func(df, parse_by=['Station_location', 'date_bin'],
                                                                  depth_bin=False)

        # Save NBSS results
        NBSS_1a.to_csv(str(dirpath) + '/' + instrument +'_NBSS_1a.csv', index=False)
        lin_fit_1b.to_csv(str(dirpath) + '/' + instrument + '_lin_fit_1b.csv', index=False)

    elif replace == 'N':
        print('previous products 1a and 1b will be kept')

elif not os.path.exists(dirpath):
    os.mkdir(dirpath)
    df = pd.concat(map((lambda path: (pd.read_csv(path))), file_list)).reset_index(drop=True)
    df = df.filter(
        ['Sample', 'date_bin', 'Station_location', 'midDepthBin', 'Min_obs_depth', 'Max_obs_depth', 'range_size_bin',
         'sizeClasses', 'Biovolume', 'midLatBin', 'midLonBin', 'Volume_imaged'], axis=1)

    if depth_binning == 'Y':
        NBSS_1a, lin_fit_1b = parse_NBS_linfit_func(df, parse_by=['Station_location', 'date_bin'],
                                                    depth_bin=True)
    else:
        NBSS_1a, lin_fit_1b = parse_NBS_linfit_func(df, parse_by=['Station_location', 'date_bin'],
                                                    depth_bin=False)

    # Save NBSS results
    NBSS_1a.to_csv(str(dirpath) + '/' + instrument + '_NBSS_1a.csv', index=False)
    lin_fit_1b.to_csv(str(dirpath) + '/' + instrument + '_lin_fit_1b.csv', index=False)


lat = lin_fit_1b['latitude']
lon = lin_fit_1b['longitude']
slope = lin_fit_1b['slope']
intercept_t = lin_fit_1b['intercept']

intercept_plot = [x*3 for x in intercept_t]
ax = plt.axes(projection=ccrs.PlateCarree())
plt.gca().coastlines('50m')
g1=ax.gridlines(draw_labels=True)
g1.xlines = False
g1.ylines = False
plt.scatter(lon, lat, label=None, c=slope, cmap='viridis', s=intercept_plot, linewidth=0, alpha=0.5, transform=ccrs.PlateCarree())
plt.gca().coastlines('50m')
ax.set_extent([-180, 180, -90, 90])
#ax.xlabel('longitude')
#ax.ylabel('latitude')
plt.colorbar(label='slope', orientation='horizontal', anchor=(0.5, 1))
#ax.clim(min(slope), max(slope))


labels = [str(np.round(min(intercept_t), decimals=2)),
          str(np.round(st.median(intercept_t), decimals=2)),
          str(np.round(max(intercept_t), decimals=2))]

for n, area in enumerate([18, 28.5, 73.5]):
    plt.scatter([], [], c='k', alpha=0.3, s=area, label=labels[n], transform=ccrs.PlateCarree())

plt.legend(bbox_to_anchor=(0.75, 0), ncol = 3, scatterpoints=1, frameon=False,
           labelspacing=1, title='intercept')

figname = 'step4_slopes_intercept_updated_firstrelease' + instrument + '.pdf'

path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

savepath = Path(cfg['git_dir']).expanduser() / 'figures' / figname

#plt.savefig(savepath)
