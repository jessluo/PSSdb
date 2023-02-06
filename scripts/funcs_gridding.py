

import numpy as np
import pandas as pd
import os
from pathlib import Path
from glob import glob
import yaml

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
    path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
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

    elif instrument == 'Zooscan':
        df = df[df['Category'].str.contains('artefact') == False]
        df = df[df['Category'].str.contains('detritus') == False]
        df = df.reset_index(drop=True)
        ind = metadata.loc[(metadata['Description'].str.contains('object_area')==True) &
                           (metadata['Description'].str.contains('exc')== False)].index[0]
    for i in range(0, len(df)):
        r = m.sqrt((df.iloc[i, [ind]]/m.pi))
        #df.loc[i, 'ESD'] = r*2
        df.loc[i, 'Biovolume'] = (4/3) * m.pi * (r**3)

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
    if instrument == 'IFCB':
        cat_remove = ['beads','detritus', 'artefact', 'bubble']
        if keep_cat == 'none':
            cat_remove = cat_remove
        else:
            cat_remove = [x for x in cat_remove if x not in keep_cat]
        for i in cat_remove:
            try:
                df = df[df['Category'].str.contains(i) == False]
                df = df.reset_index(drop=True)
            except:
                pass
    elif (instrument == 'Zooscan' or instrument == 'UVP'):
        cat_remove = ['artefacts','detritus']
        if keep_cat == 'none':
            cat_remove = cat_remove
        else:
            cat_remove = [x for x in cat_remove if x not in keep_cat]
        for i in cat_remove:
            try:
                df = df[df['Category'].str.contains(i) == False]
                df = df.reset_index(drop=True)
            except:
                pass
    for i in range(0, len(df)):
        r = m.sqrt((df.loc[i, ['Area']]/m.pi))
        #df.loc[i, 'ESD'] = r*2
        df.loc[i, 'Biovolume'] = (4/3) * m.pi * (r**3)

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

def date_binning_func(date, group_by= 'yyyymm'): # consider adding a day/night column, this will have to consider latitude and month
    """
    Objective: reduce the date information, so that the data can be binned by month, year, or month and year
    :param date: column of a  standardized dataframe containing date information ('Sampling_date' in standardized ecotaxa projects)
    with format yyyymmdd (GMT)
    :param group_by: range of time to bin the data, options are: 'year', 'month' and 'year_month'
    :return:
    """
    import pandas as pd
    import calendar
    import numpy as np
    calendar.setfirstweekday(6)
    date = date.astype(str)
    date_bin = pd.to_datetime(date, format='%Y%m%d')
    #time = time.astype(str)

    # time of day is not considered due to discrepancies of the time format
    #if time.any() == 0:
        #date_bin = date
        #date_bin = date_bin.astype(int)
        #date_bin = pd.to_datetime(date_bin, format='%Y%m%d')
    #else:
        #date_bin = (date + time)
        #date_bin = date_bin.astype(int)
        #date_bin = pd.to_datetime(date_bin, format='%Y%m%d%H%M%S')

    if group_by == 'yyyy':
        date_bin = date_bin.dt.strftime('%Y')
    elif group_by == 'yyyymm':
        date_bin = date_bin.dt.strftime('%Y%m')
    elif group_by == 'mm':
        date_bin = date_bin.dt.strftime('%m')
    elif group_by == 'week':
        def get_week_of_month(year, month, day):
            x = np.array(calendar.monthcalendar(year, month))
            week_of_month = np.where(x == day)[0][0]
            if week_of_month == 0:
                week_of_month = week_of_month +1
            elif week_of_month == 5:
                week_of_month = week_of_month - 1
            return (week_of_month)
        day = date_bin.dt.day
        month = date_bin.dt.month
        year = date_bin.dt.year
        for i in range(0, len(date_bin)):
            week=get_week_of_month(year[i], month[i], day[i])
            date_bin[i] = str(year[i]) + str(month[i]) + str(week)
    elif group_by == 'None':
        date_bin == date_bin

    return date_bin

