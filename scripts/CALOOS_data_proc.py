# GOAL: to read the files from the IFCB dashboard download, and append date, volume filtered, latitude, longitude
# and depth information to the features files.

#1) read the .hdr files and extract the valuable information (sampling date and volume filtered)

# focus on getting NBSS from caloos

import pandas as pd
import yaml
from pathlib import Path
import glob
import os
import ntpath
import pandas as pd
import numpy as np
import math as m

# returns path to project data and list of projects based on intrument (IFCB, Zooscan, UVP)
# read git-tracked config file (text file) with inputs:  project ID, output directory
path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_to_data = Path(cfg['IFCB_dir']).expanduser()

features_filelist = glob.glob(str(path_to_data) + '/**/**/*.csv')


df_SC_wharf = pd.DataFrame()

for i in features_filelist:
    df = pd.read_csv(i, sep=',',  header=0, index_col=[0])
    df = df.loc[(df['Area'] != 0)]
    head, tail = ntpath.split(i)
    df.loc[:, 'Sampling_date'] = ''.join(list(tail)[1:9]) # the variables included below
    df.loc[:, 'Sampling_time'] = ''.join(list(tail)[10:16])
    df.loc[:, 'Latitude'] = 36.96149
    df.loc[:, 'Longitude'] = -122.02187
    df.loc[:, 'Volume_analyzed'] = 0.005 # THIS IS NOT REAL, we definitely need to extract this from the dashboard
    df.loc[:, 'Project_ID'] = 'SC_wharf'
    pixel_to_micron = (3.14*3.14*3.14) # got this scale from ecotaxa projects, cubed to get biovolume
    df.loc[:, 'Area_um'] = df.loc[:, 'Area'] * pixel_to_micron
    df_SC_wharf = pd.concat([df_SC_wharf, df]).reset_index(drop=True)

for i in ['texture', 'Wedge', 'Ring', 'HOG']:
    df_SC_wharf = df_SC_wharf[df_SC_wharf.columns.drop(list(df_SC_wharf.filter(regex=i)))]

# at this point, the dataset cannot be incorporated into the pipeline because there is no metadata associated to it,
#biovolume will be calculated here and other binning functions will be used here:
for n, i in enumerate(df_SC_wharf['Area_um']):
    r = np.sqrt((i/m.pi))
    df_SC_wharf.loc[n, 'Biovolume'] = (4/3) * m.pi * (r**3)

df = df_SC_wharf

df['date_bin'] = date_binning_func(df['Sampling_date'], df['Sampling_time'], group_by =  'yyyymm')
df['Station_location'], df['midLatBin'], df['midLonBin'] = station_binning_func(df['Latitude'], df['Longitude'])

NBSS_binned_all = parse_NBS_func(df, parse_by = ['date_bin'])

NBSS_binned_all['date_bin'] = NBSS_binned_all['date_bin'].astype(str)

fig, ax1 = plt.subplots(1,1)
for i in set(NBSS_binned_all['date_bin']):
    NBS_df = pd.DataFrame()
    NBS_df = NBSS_binned_all.loc[(NBSS_binned_all['date_bin'] == i)]
    ax1.plot(NBS_df['logSize'], NBS_df['logNBSS'], marker='.', label=i)
ax1.set_ylabel('log (NBS)')
ax1.set_xlabel('log (Biovolume)')
ax1.legend()
ax1.title.set_text('NBSS from Santa Cruz Municipal Wharf, separated by month and year (yyyymm)')
fig.savefig('NBSS_prelim_results_CALOOS.png', bbox_inches='tight')
