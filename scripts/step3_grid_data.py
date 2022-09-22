# functions related to gridding standardized files generated in 2. bin_func performs all spatiotemporal gridding of
# the data in one step

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import read_funcs as rd # functions to list standardized and binned files

from pathlib import Path
from glob import glob
# Config modules
import yaml  # requires installation of PyYAML package

def biovol_func(df, instrument, area_type= 'object_area', remove_cat='none'):
    """
    Objective: calculate biovolume (in cubic micrometers) of each object, only for UVP and Zooscan projects, following
    the volume of an ellipsoid OR a sphere. Also, determine which area will be used to calculate the biovolume.
    This function also removes Biovolumes= 0 in the IFCB files
    :param df: a STANDARDIZED dataframe that contains object's Area (should be in micrometers^2)
     and object's Minor_axis (should be in micrometers)
    :param area_type: name of the column that contains the area to be used to calculate biovolume
    :param geom_shape: geometric shape to be used to calculate biovolume, either 'sphere' or 'ellipse'
    :return: Biovolume in cubic micrometers
    """
    import math as m
    from glob import glob
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



def depth_binning_func(depth='Depth_min', depth_bins=[0, 25, 50, 100, 200, 500, 1000, 3000, 8000]):
    # create depth bins based on Jessica's suggestion (Jessica's suggestions as default)
    depth_bin = pd.cut(x=depth, bins=depth_bins, include_lowest=True)
    # create a categorical variable with the middle depth of each bin
    midDepth_bin = []
    for i in range(0, len(depth)):
        midDepth_bin.append(depth_bin.iloc[i].mid)
    return midDepth_bin

    #6.3 by station (geographical location)
def gridding_func(lat= 'Latitude', lon= 'Longitude', st_increment=1):
    import numpy as np
    # create a categorical variable that will serve as the station ID based on binned lat and lon
    # 1) create arrays that will serve to bin lat and long with 1 degree increments
    lat_increments = np.arange(np.floor(min(lat) - st_increment),
                               np.ceil(max(lat) + st_increment), 1)
    lon_increments = np.arange(np.floor(min(lon) - st_increment),
                               np.ceil(max(lon) + st_increment), 1)
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

def date_binning_func(date, time, group_by= 'yyyymm'):
    """
    Objective: reduce the date information, so that the data can be binned by month, year, or month and year
    :param date: column of a  standardized dataframe containing date information ('Sampling_date' in standardized ecotaxa projects)
    with format yyyymmdd (GMT)
    :param group_by: range of time to bin the data, options are: 'year', 'month' and 'year_month'
    :return:
    """
    date = date.astype(str)
    time = time.astype(str)
    date_bin = (date + time)
    date_bin = date_bin.astype(int)
    date_bin = pd.to_datetime(date_bin, format='%Y%m%d%H%M%S')
    if group_by == 'yyyy':
        date_bin = date_bin.dt.strftime('%Y')
    elif group_by == 'yyyymm':
        date_bin = date_bin.dt.strftime('%Y%m')
    elif group_by == 'mm':
        date_bin = date_bin.dt.strftime('%m')
    elif group_by == 'None':
        date_bin == date_bin
    return date_bin


def bin_func(instrument, removecat, area_type = 'object_area', date_group = 'yyyymm',   ignore_depth = 'no'):
    """
    Objective: same as binning_NBS_func but for a single dataframe
    :param df: a standardized dataframe
    :param removecat: remove a specific category, currently restricted to IFCB bubbles, artefacts or beads
    :return: a binned dataframe with NBS
    """
    import os
    import shutil
    path_to_data, id_list = rd.proj_id_list_func(instrument, data_status='standardized')  # generate path and project ID's
    if instrument == 'IFCB':  # this removal is only for testing,
        # since the standardizer should be able to deal with projects with empty rows
        IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
        for element in IFCB_empty_rows:
            if element in id_list:
                id_list.remove(element)
    dirpath = str(path_to_data) + '/gridded_data/'
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    for i in id_list:
        df = rd.read_func(path_to_data, i)
        df = biovol_func(df, instrument= instrument, area_type= area_type, remove_cat=removecat)
        df['date_bin'] = date_binning_func(df['Sampling_date'], df['Sampling_time'], group_by =  date_group)
        df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(df['Latitude'], df['Longitude'])
        metadata_bins = pd.DataFrame({'Variables':['Biovolume', 'date_bin', 'Station_location', 'midLatBin', 'midLonBin'],
                                 'Variable_types': ['float64',  'int64', 'object', 'float64', 'float64'],
                                 'Units/Values/Timezone': ['cubic_micrometer',  date_group, 'lat_lon', 'degree', 'degree'],
                                 'Description': ['Biovolume calculated as a spherical projection of ' + area_type ,
                                                 'binned date information',
                                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                                 'latitude of the center point of the 1x1 degree cell',
                                                 'longitude of the center point of the 1x1 degree cell']})
        if ignore_depth == 'no':
            df['midDepthBin'] = depth_binning_func(df['Depth_max'])
            metadata_bins = pd.DataFrame({'Variables':['Biovolume',  'date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'midDepthBin'],
                                 'Variable_types': ['float64',  'int64', 'object', 'float64', 'float64', 'float64'],
                                 'Units/Values/Timezone': ['cubic_micrometer',  date_group, 'lat_lon', 'degree', 'degree', 'meters'],
                                 'Description': ['Biovolume calculated as a spherical projection of ' + area_type ,
                                                 'binned date information',
                                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                                 'latitude of the center point of the 1x1 degree cell',
                                                 'longitude of the center point of the 1x1 degree cell',
                                                 'middle point within a depth bin']})

        df.to_csv(str(path_to_data) + '/gridded_data/' + instrument + '_' + str(i) + '_gridded.tsv', sep='\t')
        metadata_bins.to_csv(str(path_to_data) + '/gridded_data/' + instrument + '_' + str(i) + '_gridded_metadata.tsv',
                           sep='\t')
    #return df, metadata_bins
