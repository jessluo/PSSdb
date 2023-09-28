


import os
from pathlib import Path
import re
from glob import glob
import yaml
import pandas as pd
import calendar
import numpy as np
import datetime as dt
from tqdm import tqdm
import pytz # to define the timezone of the time info
from datetime import datetime
from astral.sun import sun #run pip install astral
# to assign day/night to an ROI, first create a LocationInfo element. See instructions: https://astral.readthedocs.io/en/latest/
from astral import LocationInfo
# an alternative to astral
from suncalc import get_position, get_times # pip install suncalc
#from dtki import d_n_calc # R. Kiko's module, see PSSdb documentation to install
try:
    from funcs_read import *
except:
    from scripts.funcs_read import *

# Module astral for day/night determination:
import astral
from astral.sun import sun # Use pip3 install astral
from timezonefinder import TimezoneFinder
tf = TimezoneFinder()
import datetime
from suntime import Sun
import ephem
from calendar import month_abbr

def daynight_M(df_standardized):
    """
    Objective: Determine the period of sampling (day or night) according to sampling time and sampling position
    :param df_standardized: a STANDARDIZED dataframe that contains sample date (Sampling_date), time (Sampling_time). latitude (Latitude), longtiude (Longitude)
    :return: a string Day or Night, or Empty if required columns are missing
    """
    if ('Sampling_date' in df_standardized.index) and ('Sampling_time' in df_standardized.index):
        df_standardized['Sampling_date'] = str(df_standardized['Sampling_date']).zfill(6)
        df_standardized['Sampling_time'] = str(df_standardized['Sampling_time']).zfill(6)
        df_standardized['Sampling_time'] = pd.to_datetime(df_standardized['Sampling_date']+" "+df_standardized['Sampling_time'],format="%Y%m%d %H%M%S", utc=True)

        try:
             return 'Day' if (pd.to_datetime(df_standardized.Sampling_time, format="%Y%m%d %H%M%S", utc=True).to_pydatetime() >= Sun(df_standardized.Latitude, df_standardized.Longitude).get_sunrise_time(pd.to_datetime(df_standardized.Sampling_time, format="%Y%m%d %H%M%S", utc=True).to_pydatetime())) &  (pd.to_datetime(df_standardized.Sampling_time, format="%Y%m%d %H%M%S", utc=True).to_pydatetime() <= Sun(df_standardized.Latitude, df_standardized.Longitude).get_sunset_time(pd.to_datetime(df_standardized.Sampling_time, format="%Y%m%d %H%M%S", utc=True).to_pydatetime())) else 'Night'
        except: # Required as high latitudes result in error
            o=ephem.Observer()
            o.lat, o.long, o.date = str(df_standardized.Latitude), str(df_standardized.Longitude),pd.to_datetime(df_standardized.Sampling_time, format="%Y%m%d %H%M%S", utc=True)
            sun = ephem.Sun(o)
            try :
                next_sunset = o.next_setting(sun, start=o.date)
            except ephem.NeverUpError as error:
                return 'Night'
            except ephem.AlwaysUpError as error:
                return 'Day'
    else:
        return ''

def daynight(date_time, Latitude, Longitude):
    """
    Objective: Determine the period of sampling (day or night) according to sampling time and sampling position
    :param df_standardized: a STANDARDIZED dataframe that contains sample date (Sampling_date), time (Sampling_time). latitude (Latitude), longtiude (Longitude)
    :return: a string Day or Night, or Empty if required columns are missing
    """

    try:
        if (pd.to_datetime(date_time, format="%Y%m%d %H%M%S", utc=True).to_pydatetime() >= Sun(Latitude,Longitude).get_sunrise_time(pd.to_datetime(date_time, format="%Y%m%d %H%M%S", utc=True).to_pydatetime())) \
                & (pd.to_datetime(date_time, format="%Y%m%d %H%M%S", utc=True).to_pydatetime() <= Sun(Latitude,Longitude).get_sunset_time(pd.to_datetime(date_time, format="%Y%m%d %H%M%S", utc=True).to_pydatetime())):
            light_cond = 'Day'
        else:
            light_cond = 'Night'
    except: # Required as high latitudes result in error
        o=ephem.Observer()
        o.lat, o.long, o.date = str(Latitude), str(Longitude),pd.to_datetime(date_time, format="%Y%m%d %H%M%S", utc=True)
        sun = ephem.Sun(o)
        try :
            next_sunset = o.next_setting(sun, start=o.date)
        except ephem.NeverUpError as error:
            light_cond =  'Night'
        except ephem.AlwaysUpError as error:
            light_cond = 'Day'

    return light_cond

# NOTE:  biovol_func_old deprecated as of 1/11/2023. This function allowed the user to choose which type of area estimation to use to
#get biovolume
def biovol_func_old(df, instrument, area_type= 'object_area', remove_cat='none'):
    """
    Objective: calculate biovolume (in cubic micrometers) of each object, only for UVP and Zooscan projects, following
    the volume of  a sphere. Also, determine which area will be used to calculate the biovolume.
    This function also removes Biovolumes= 0 in the IFCB files
    :param df: a STANDARDIZED dataframe that contains object's Area (should be in micrometers^2)
     and object's Minor_axis (should be in micrometers)
    :param area_type: name of the column that contains the area to be used to calculate biovolume
    :param geom_shape: geometric shape to be used to calculate biovolume, either 'sphere' or 'ellipse'
    :return: Biovolume in cubic micrometers
    """
    import math as m
    from pathlib import Path
    from glob import glob
    import yaml
    import pandas as pd
    # find the Area column in the metadata file and use its position to select the area
    ID = df.loc[0, 'Project_ID']
    path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    path_to_metadata = glob(str(Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'])+ '/**/*'+ str(ID)+ '*metadata*')
    metadata = pd.read_csv(path_to_metadata[0], sep='\t')
    if instrument == 'IFCB':# for IFCB projects, it needs to remove the column that is repeated. Also, remove bubbles
        if 'summed' in area_type:
            df = df.drop(df.columns[[19]], axis=1)
            df = df[df.Biovolume != 0]
            df = df[df.Category != 'bubble']
            df = df.reset_index(drop=True)
        else:
            df = df.drop(df.columns[[20]], axis=1)
            df = df[df.Biovolume != 0]
            df = df[df.Category != 'bubble']
            df = df.reset_index(drop=True)
            if 'summed' in area_type:
                ind = metadata.loc[metadata['Description'].str.contains('summed')].index[0]
            else:
                ind = metadata.loc[(metadata['Description'].str.contains('object_area') == True) &
                                   (metadata['Description'].str.contains('summed') == False)].index[0]
    if instrument == 'IFCB':
        if remove_cat == 'beads':
            df = df[df['Category'].str.contains("bead") == False]
            df = df.reset_index(drop=True)
        if remove_cat == 'artefact':
            df = df[df['Category'].str.contains('artefact') == False]
            df = df.reset_index(drop=True)
        elif remove_cat == 'none':
            df = df
            if 'summed' in area_type:
                ind = metadata.loc[metadata['Description'].str.contains('summed')].index[0]
            else:
                ind = metadata.loc[(metadata['Description'].str.contains('object_area') == True) &
                                   (metadata['Description'].str.contains('summed') == False)].index[0]

    elif instrument == 'Scanner':
        df = df[df['Category'].str.lower().apply(lambda annotation: len(re.findall(r'bead|bubble|artefact|artifact|glue', annotation)) == 0 if str(annotation) != 'nan' else True)]

        df = df.reset_index(drop=True)
        ind = metadata.loc[(metadata['Description'].str.contains('object_area')==True) &
                           (metadata['Description'].str.contains('exc')== False)].index[0]
    for i in range(0, len(df)):
        r = m.sqrt((df.iloc[i, [ind]]/m.pi))
        #df.loc[i, 'ESD'] = r*2
        df.loc[i, 'Biovolume'] = (4/3) * m.pi * (r**3)

    return df

def remove_UVP_noVal(df):
    '''
    Objective: first, remove all UVP ROI's with status = unclassified. For the remaining ROI's with annotation,
    calculate the percentage of validated ROI's. if the sample has <95% validation, discard
    '''
    df['Category'].replace('', np.nan, inplace=True)
    df = df.dropna(subset = ['Category']).reset_index(drop=True)

    if len(df.dropna(subset=['Category'])) == 0:
        df = pd.DataFrame()
        return df
    else:
        val = df.groupby(['Sample']).apply(lambda x: pd.Series({'Validation_percentage':len(x.dropna(subset=['Category'])[x.dropna(subset=['Category'])['Annotation'].isin(['validated'])].ROI) /
                                                                            len(x.dropna(subset=['Category']).ROI)}))
        for s in df.Sample.unique():
            if val.loc[s, 'Validation_percentage'] < 0.95:
                index = df[df['Sample'] == s].index
                df.drop(index, inplace=True)
            else:
                df.loc[df['Sample']==s, 'Validation_percentage'] = val.loc[s, 'Validation_percentage']

        return df

def biovol_func(df, instrument, keep_cat='none'):
    """
    Objective: calculate biovolume (in cubic micrometers) of each object  following
    the volume of  a sphere.
        :param df: a STANDARDIZED dataframe (output of step2_standardize_projects)
    This function by default removes objects based on 'Category':
       instrument = 'IFCB': removes 'beads' and 'artefact'
       instrument = 'Zooscan': removes 'detritus' and 'artefact'
       instrument = 'UVP': removes 'detritus' and 'artefact'
    :param keep_cat: list of strings that determines if categories that should be removed are kept
    :return: the original dataframe Biovolume in cubic micrometers
    """
    import math as m

    path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    path_to_taxonomy = str(Path(cfg['git_dir']).expanduser()) + '/ancillary/plankton_annotated_taxonomy.xlsx'
    taxonomy_df = pd.read_excel(path_to_taxonomy)
    df['Category'] = df['Category'].fillna('')
    df['Cat_remove'] = pd.merge(df, taxonomy_df, how='left', on=['Category'])['EcoTaxa_hierarchy']
    df['Cat_remove'] = df['Cat_remove'].fillna('')
    for n, i in enumerate(df['Cat_remove']):
        if i == '':
            df.loc[n, 'Cat_remove'] = df.loc[n, 'Category']
        else:
            pass

    if instrument == 'IFCB':
        cat_remove = ['beads','detritus', 'artefact', 'bubble']
        if keep_cat == 'none':
            cat_remove = cat_remove
        else:
            cat_remove = [x for x in cat_remove if x not in keep_cat]
        for i in cat_remove:
            try:
                df = df[df['Cat_remove'].str.contains(i) == False]
                df = df.reset_index(drop=True)
            except:
                pass
    elif (instrument == 'Zooscan'):
        cat_remove = ['artefact','detritus', 'plastic']
        if keep_cat == 'none':
            cat_remove = cat_remove
        else:
            cat_remove = [x for x in cat_remove if x not in keep_cat]
        for i in cat_remove:
            try:
                df = df[df['Cat_remove'].str.contains(i) == False]#df[~df.Category.str.contains('|'.join(i))]
                df = df.reset_index(drop=True)
            except:
                pass
    elif (instrument == 'UVP'):
        cat_remove = ['artefact','detritus', 'plastic']
        if keep_cat == 'none':
            cat_remove = cat_remove
        else:
            cat_remove = [x for x in cat_remove if x not in keep_cat]
        for i in cat_remove:
            try:
                df = df[df['Cat_remove'].str.contains(i) == False]#df[~df.Category.str.contains('|'.join(i))]
                df = df.reset_index(drop=True)
            except:
                pass
    df['Biovolume_area'] = (4/3)* m.pi * (((df['Area']/m.pi)**0.5) **3) # convert area to ESD, then calculate biovolume
    if 'Minor_axis' in df.columns:
        df['Biovolume_ellipsoid'] = (4/3)*df['Area']*df['Minor_axis']
    else:
        None
    if 'Biovolume' in df.columns:
        df = df.rename(columns={'Biovolume': 'Biovolume_orig'}).reset_index()
    else:
        None

    df = df.drop(columns=['Cat_remove']).reset_index(drop=True)


    return df




def depth_binning_func(depth='Depth_min', depth_bins=[0, 25, 50, 100, 200, 500, 1000, 3000, 8000]):
    # create depth bins based on Jessica's suggestion (Jessica's suggestions as default)
    import pandas as pd
    depth_bin = pd.cut(x=depth, bins=depth_bins, include_lowest=True)
    # create a categorical variable with the middle depth of each bin
    midDepth_bin = []
    for i in range(0, len(depth)):
        midDepth_bin.append(depth_bin.iloc[i].mid)
    return midDepth_bin



    #6.3 by station (geographical location)
def gridding_func(st_increment, lat= 'Latitude', lon= 'Longitude'):
    import numpy as np
    import pandas as pd
    # create a categorical variable that will serve as the station ID based on binned lat and lon
    # 1) create arrays that will serve to bin lat and long with 1 degree increments
    lat_increments = np.arange(np.floor(min(lat) - st_increment),
                               np.ceil(max(lat) + st_increment), st_increment)
    lon_increments = np.arange(np.floor(min(lon) - st_increment),
                               np.ceil(max(lon) + st_increment), st_increment)
    # 2) create bins
    lat_bin = pd.cut(x=lat, bins=lat_increments, include_lowest=True)

    lon_bin = pd.cut(x=lon, bins=lon_increments, include_lowest=True)
    # 3) get middle lat after assigning each sample to a lat/lon bin
    midLat_bin = []
    for i in range(0, len(lat_bin)):
        midLat_bin.append(lat_bin.iloc[i].mid)
    midLon_bin = []
    for i in range(0, len(lon_bin)):
        midLon_bin.append(lon_bin.iloc[i].mid)
    # create a station ID
    Station_ID = []
    for i in range(0, len(midLon_bin)):
        Station_ID.append(str(midLat_bin[i]) + '_' + str(midLon_bin[i]))

    return Station_ID, midLat_bin, midLon_bin

def group_gridded_files_func(instrument, already_gridded= 'N'):
    """
    Objective: assign a label to each gridded dataset, so when the compilation of files occurs, it is only done with a group of data and prevents
    the creation of a huge file
    """
    #first, create series of numbers that break the globe into 15x15 degree cells:
    #lat_left = np.arange(-90, 90+15, 1, dtype=int)
    #lat_right = np.arange(-75,105 + 15, 1, dtype=int)
    #lat_int = pd.IntervalIndex.from_arrays(lat_left, lat_right, closed='both')
    #lon_left = np.arange(-180, 180+15, 1, dtype=int)
    #lon_right = np.arange(-165, 195+15, 1, dtype=int)
    #lon_int = pd.IntervalIndex.from_arrays(lon_left, lon_right, closed='both')
    lat_int = pd.interval_range(start=-90, end=90, closed='both')
    lon_int = pd.interval_range(start=-180, end=180, closed='both')
    grid_list = []
    # create a dictionary that has unique keys for each grid interval
    #grid_dict = {}
    #for x in range(0, len(lon_int)):
        #for y in range(0, len(lat_int)):
            #grid_dict[(str(x) + '_' + str(y))] = {'lon': lon_int[x], 'lat': lat_int[y]} ## however, this might not be necessary
    # try now to assing lat & lon grid numbers to the dataframe directly
    file_list = proj_id_list_func(instrument, data_status='gridded')
    for i in tqdm(file_list):
        print('generating subsets of data for 1 degree grids from ' + i)
        filename = i.split('/')[-1]
        dirpath = i.replace(filename, '')
        df = pd.read_csv(i, header = 0)
        df['lat_grid'] = None
        df['lon_grid'] = None
        for la in df['midLatBin'].unique():
            df.loc[df['midLatBin'] == la, 'lat_grid']= str(lat_int[lat_int.contains(la) == True][0].left) + '-' + str(lat_int[lat_int.contains(la) == True][0].right)
        for lo in df['midLonBin'].unique():
            df.loc[df['midLonBin'] == lo, 'lon_grid'] = str(lon_int[lon_int.contains(lo) == True][0].left) + '-' + str(lon_int[lon_int.contains(lo) == True][0].right)
        # combine these numbers to form a grid identifier
        df['grid_id'] = df.lat_grid.astype(str).str.cat(df.lon_grid.astype(str), sep='_')
        df_subgrouped = [x for _, x in df.groupby('grid_id')]
        for s in df_subgrouped:
            grid_list.append(str(s['grid_id'].unique()[0]))
            if already_gridded == 'N':
                s.to_csv(str(dirpath) + 'grid_N_'+str(s['grid_id'].unique()[0]) + '_'+ filename, index=False)
            #os.remove(i)
    grid_list_unique = [*set(grid_list)]
    grid_list_unique = ['N_' + s for s in grid_list_unique]
    return grid_list_unique

def date_binning_func(df, group_by= 'yyyymm'): # , ignore_high_lat=True, consider adding a day/night column, this will have to consider latitude and month
    """
    Objective: reduce the date information, so that the data can be binned by month, year, or month and year. Also create a column that assigns a 'day' or 'night' category to each ROI
    :param date: column of a  standardized dataframe containing date information ('Sampling_date' in standardized ecotaxa projects)
    with format yyyymmdd (UTC)
    :param time: column of a  standardized dataframe containing time ('Sampling_time' in standardized ecotaxa projects)
    with format HHMMSS  (UTC)
    :param lat:
    :param lon:
    :param group_by: range of time to bin the data, options are: 'year', 'month' and 'year_month'
    :return:
    """
    df['Sampling_time']=df['Sampling_time'].astype(str)
    df['Sampling_date'] = df['Sampling_date'].astype(str)
    df['Sampling_time'] = df['Sampling_time'].apply(lambda x: x.zfill(6)) #x.ljust(6, '0')

    date = df['Sampling_date'].astype(str)
    date_bin = pd.to_datetime(date, format='%Y%m%d')
    df['year'] = pd.DatetimeIndex(date_bin).year.values
    df['month'] = pd.DatetimeIndex(date_bin).month.values
    df['week'] = pd.DatetimeIndex(date_bin).isocalendar().week.values
    dict_week = df.groupby(['week']).apply(lambda x: pd.Series({ 'week': x.week.unique()[0],
                                                               'month': np.round(x.month.mean())})).set_index('week')['month'].to_dict()
    df['month'] = df['week'].map(dict_week)

    if group_by == 'yyyy':
        df['date_bin'] = df['year'].astype(str)
    elif group_by == 'yyyymm':
        df['date_bin'] = df['year'].astype(str) + '_' + df['month'].astype(int).astype(str).str.zfill(2)
    elif group_by == 'mm':
        df['date_bin'] = df['month'].astype(int).astype(str).str.zfill(2)
    elif group_by == 'week':
        df['date_bin'] = df['year'].astype(str) + '_' + df['month'].astype(int).astype(str).str.zfill(2) + '_wk' +df['week'].astype(int).astype(str).str.zfill(2)
    elif group_by == 'None':
        df['date_bin'] == str(date_bin)

    #df['date_grouping'] = df['year'].astype(str) + '-' + df['month'].astype(int).astype(str).str.zfill(2)

    df['time_full'] = pd.to_datetime(df['Sampling_date'].astype(str) + " " + df['Sampling_time'].astype(str), format="%Y%m%d %H%M%S", utc=True)
    df['light_cond'] = df.apply(lambda x: daynight(x['time_full'], x['Latitude'], x['Longitude']), axis = 1)
    df = df.drop(columns=['year', 'month', 'week', 'time_full']).reset_index(drop=True)
    return df


#daynight(date_time, Latitude, Longitude):