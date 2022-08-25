# create a python function that receives a folder as an imput (such folder contains csv files)
# another imput should be bin count
# and returns: slope, intercept and R^2 of the regression (do we want plots?)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

from pathlib import Path
from glob import glob
# Config modules
import yaml  # requires installation of PyYAML package
from scipy.stats import linregress

# functions here




# 2) function to create list of projects. Inputs: instrument
def proj_id_list_func(instrument, standardized='no',  testing=False):
    """
    Objective:
    generate a list of the downloaded Ecotaxa projects, and return the path of the files
    :param instrument: the device used for image adcquisition. important since the path will change according to it
    :param standardized: will change the paths whether the data of interest is standardized or not
    :param testing: restrict the id list to only the ecotaxa projects that are for testing functions
    :return path_to_data: string where the tsv files are stores
    :return id_list: list with the ecotaxa project identifiers (might not be relevant once we expand
    to get data from other platforms)
    """
    import pandas as pd
    import yaml
    from pathlib import Path
    # returns path to project data and list of projects based on intrument (IFCB, Zooscan, UVP)
    # read git-tracked config file (text file) with inputs:  project ID, output directory
    path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    # read config file (text file) with inputs:  project ID, output directory
    # prepare storage based on path stored in the yaml config file and instrument type
    if standardized == 'no':
        path_to_data = Path(cfg['git_dir']).expanduser() / cfg['dataset_subdir'] / instrument
    elif standardized == 'yes':
        path_to_data = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / instrument
    elif standardized == 'binned':
        path_to_data = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / instrument / cfg['binned_subdir']
    # create a list  projects that we have access to, based on project_list_all.xlsx
    path_to_proj_id_list = Path(cfg['git_dir']).expanduser() / cfg['proj_list']
    proj_list = pd.read_excel(path_to_proj_id_list)
    if testing == False:
        id_list = proj_list['Project_ID'].loc[ (proj_list['Instrument'] == str(instrument)) & (proj_list['PSSdb_access'] == True)].tolist() # instrument type filtered here
    else:
        if instrument == 'IFCB':
            id_list = [3315, 3318, 3326]  # THIS IS THE TEST LIST, USE THIS WHEN DEBUGGING for IFCB
        elif instrument == 'Zooscan':
            id_list = [1]  # ask mathilde about her test projects
        elif instrument == 'UVP':
            id_list = [2]  # ask mathilde about her test projects

    return path_to_data, id_list

# 3) Standardize all files:
def mass_st_func(instrument):
    """"
    Objective:
    Complete standardization of all projects for a given instrument
    """
    #import functions to 1. get project ID lists and
    from funcs_size_spectra import proj_id_list_func as id_list
    #get the standardizer function
    from funcs_export_standardizer import standardization_func as standardize
    from pathlib import Path
    import yaml
    import glob
    instrument = str(instrument)
    #define the paths to the config file and the standardizer using the yaml file
    path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    # read config file (text file) with inputs:  project ID, output directory
    standardizer_path = glob.glob(str(Path(cfg['raw_dir']).expanduser()) + '/' + '*' + instrument+ '*xlsx')[0]
    config_path = str(Path('~/GIT/PSSdb/scripts/Ecotaxa_API_pw.yaml').expanduser())
    path_to_data, ID_list = id_list(instrument, testing=False)
    #for i in ID_list:
        #standardize(config_path, standardizer_path, i)
    if instrument == 'IFCB':# this removal is only for testing,
        #since the standardizer should be able to deal with projects with empty rows
        IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
        for element in IFCB_empty_rows:
            if element in ID_list:
                ID_list.remove(element)
        for i in ID_list:
            standardize(config_path, standardizer_path, i)

    else:
        for i in ID_list:
            standardize(config_path, standardizer_path, i)



# 4) Create clean dataframes in python, input: a path and an ID. modify this to prevent removing rows
def read_func(path_to_data, ID):
    """
    Objective: returns a dataframe for each STANDARDIZED project
    :param path_to_data: path where the standardized files are stored
    :param ID: id number in ecotaxa of eatch project
    """
    import pandas as pd
    from glob import glob
    search_pattern = path_to_data / ("*" + str(ID) + ".tsv")
    filename = glob(str(search_pattern))
    df = pd.read_csv(filename[0], sep='\t',  header=0, index_col=[0]) #
    df = df.loc[:, ~df.columns.str.match("Unnamed")]
    # convert taxonomic Id column into categorical.
    # NOTE all of these column names might need to be modified after standardization of datasets
    #df.loc[:, cat].astype('category')
    #df['proj_ID'] = [ID] * (len(df))
    return df

# 5) function to remove rows without data
def clean_df_func(df, esd='ESD', lat= 'Latitude', lon= 'Longitude', cat='Category'):
 # FILTERING DATA: filter taxonomic categories that are artefacts
    data_clean = df[df[cat].str.contains("artefact") == False]
    # remove data points where size metrics = 0
    data_clean = data_clean[data_clean.esd != 0]
    # finally, remove any row with coordinates =nan
    data_clean = data_clean[data_clean[lat].notna()]
    data_clean = data_clean[data_clean[lon].notna()]
    data_clean = data_clean.reset_index()
    return data_clean

# include here a creation of a spreadsheet that says which information is missing

#5 biovol_standardizer
def biovol_func(df, instrument, area_type= 'object_area', geom_shape = 'sphere', remove_cat='none'):
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
    if instrument == 'IFCB':
        if remove_cat == 'beads':
            df = df[df['Category'].str.contains("bead") == False]
            df = df.reset_index(drop=True)
        if remove_cat == 'artefact':
            df = df[df['Category'].str.contains('artefact') == False]
            df = df.reset_index(drop=True)
        elif remove_cat == 'none':
            df = df

    elif instrument == 'Zooscan':
        df = df[df['Category'].str.contains('artefact') == False]
        df = df[df['Category'].str.contains('detritus') == False]
        df = df.reset_index(drop=True)
        if 'exc' in area_type:# this condition might be exclusive for Zooscan only
            ind = metadata.loc[metadata['Description'].str.contains('exc')].index[0]
        else:
            ind = metadata.loc[(metadata['Description'].str.contains('object_area')==True) &
                           (metadata['Description'].str.contains('exc')== False)].index[0]
        for i in range(0, len(df)):
            r = m.sqrt((df.iloc[i, [ind]]/m.pi))
            df.loc[i, 'ESD'] = r*2
            if geom_shape =='sphere':
                df.loc[i, 'Biovolume'] = (4/3) * m.pi * (r**3)
            elif geom_shape == 'ellipse':
                df.loc[i, 'Biovolume'] = (4/3) * m.pi * ((r/2) * (df.loc[i, 'Minor_axis']/2) * (df.loc[i, 'Minor_axis']/2))
    return df
#6) data binning
    # by depth, size, lat/lon and size, input: a dataframe, lat/lon increments that define the stations, and a list
# of depth bins. NOTE: separate this into three functions

#1) function to create size bins
def log_bins_func(start, bin_number, power=10):
    import numpy as np
    # start: lower size limit of the entire series
    # bin_number: a high enough value to contain all objects, maybe should choose something to be able to integrate blue whales == largest marine living object
    iter = range(0, bin_number)
    limits_list = [start]
    for i in iter:
        y = limits_list[i]
        x = pow(power, (np.log10(y) + np.log10(pow(2, 0.25))))
        limits_list.append(x)
    return limits_list

    #6.1 by size bins. Inputs: a column of the dataframe with biovolume, and the number of log bins to use.
    #returns two lists : the sizeClasses bins (categorical) and range_size_bins (float)
def size_binning_func(biovolume, N_log_bins=200):
    import numpy as np
    PD_SizeBins_um3 = log_bins_func((min(biovolume)), N_log_bins)  # number of bins work for this data
    SizeBins_um3 = np.array(PD_SizeBins_um3)
    # create a categorical variable, which are bins defined by the numbers on the sizeClasses list
    # Assign bin to each data point based on biovolume
    sizeClasses= pd.cut(x=biovolume, bins=SizeBins_um3, include_lowest=True)# size classes defined by biovolume
    # obtain the size of each size class bin and append it to the stats dataframe
    # unfortunately the Interval type produced by pd.cut is hard to work with. So they need to be modified
    range_size_bin = []
    for i in range(0, len(biovolume)):
        range_size_bin.append(sizeClasses.iloc[i].length)
    return sizeClasses, range_size_bin

   # 6.2 by depth. Inputs: a column of the dataframe with depth data, and a list of depth bins
def depth_binning_func(depth='Depth_min', depth_bins=[0, 25, 50, 100, 200, 500, 1000, 3000, 8000]):
    # create depth bins based on Jessica's suggestion (Jessica's suggestions as default)
    depth_bin = pd.cut(x=depth, bins=depth_bins, include_lowest=True)
    # create a categorical variable with the middle depth of each bin
    midDepth_bin = []
    for i in range(0, len(depth)):
        midDepth_bin.append(depth_bin.iloc[i].mid)
    return midDepth_bin

    #6.3 by station (geographical location)
def station_binning_func(lat= 'Latitude', lon= 'Longitude', st_increment=1):
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


def bin_func(df, instrument, removecat, area_type = 'object_area', geom_shape = 'ellipse', date_group = 'yyyymm',   ignore_depth = 'no'):
    """
    Objective: same as binning_NBS_func but for a single dataframe
    :param df: a standardized dataframe
    :param removecat: remove a specific category, currently restricted to IFCB bubbles, artefacts or beads
    :return: a binned dataframe with NBS
    """
    df = biovol_func(df, instrument= instrument, area_type= area_type, geom_shape = geom_shape, remove_cat=removecat)
    df['sizeClasses'], df['range_size_bin'] = size_binning_func(df['Biovolume'])
    df['date_bin'] = date_binning_func(df['Sampling_date'], df['Sampling_time'], group_by =  date_group)
    df['Station_location'], df['midLatBin'], df['midLonBin'] = station_binning_func(df['Latitude'], df['Longitude'])
    metadata_bins = pd.DataFrame({'Variables':['Biovolume', 'sizeClasses', 'range_size_bin', 'date_bin', 'Station_location', 'midLatBin', 'midLonBin'],
                                 'Variable_types': ['float64', 'object', 'float64', 'int64', 'object', 'float64', 'float64'],
                                 'Units/Values/Timezone': ['cubic_micrometer', 'cubic_micrometer','cubic_micrometer', date_group, 'lat_lon', 'degree', 'degree'],
                                 'Description': ['Biovolume calculated following' + geom_shape + 'projection based on' + area_type ,
                                                 'minimum and maximum value of the size bin',
                                                 'difference between max and min value of a size bin',
                                                 'binned date information',
                                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                                 'latitude of the center point of the 1x1 degree cell',
                                                 'longitude of the center point of the 1x1 degree cell']})
    if ignore_depth == 'no':
        df['midDepthBin'] = depth_binning_func(df['Depth_max'])
        metadata_bins = pd.DataFrame({'Variables':['Biovolume', 'sizeClasses', 'range_size_bin', 'date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'midDepthBin'],
                                 'Variable_types': ['float64', 'object', 'float64', 'int64', 'object', 'float64', 'float64', 'float64'],
                                 'Units/Values/Timezone': ['cubic_micrometer', 'cubic_micrometer','cubic_micrometer', date_group, 'lat_lon', 'degree', 'degree', 'meters'],
                                 'Description': ['Biovolume calculated following' + geom_shape + 'projection based on' + area_type ,
                                                 'minimum and maximum value of the size bin',
                                                 'difference between max and min value of a size bin',
                                                 'binned date information',
                                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                                 'latitude of the center point of the 1x1 degree cell',
                                                 'longitude of the center point of the 1x1 degree cell',
                                                 'middle point within a depth bin']})
    return df, metadata_bins

#6) open files as dataframes and bin them by location and depth, incorporates all binning functions (6.1-6.3) :
# do we want to do the size class limits here? 7/27/2022: silenced the binning since we want to do that with the
# merged files of IFCB and zooscan
def proj_NBS_func(instrument, removecat, time_grouping= 'yyyymm', area_type = 'object_area', geom_shape = 'ellipse', ignore_depth = 'no'):
    """
    Objective: read into the standardized tsv files and bin the data by size (biovolume), station and depth
    :param instrument: the device used for image adcquisition. important since the path will change
    :return: tsv files with the same data but binned (would we like to delete the standardized data?)
    """
    import os
    path_to_data, id_list = proj_id_list_func(instrument, standardized= 'yes', testing=False)#generate path and project ID's
    if instrument == 'IFCB':# this removal is only for testing,
        #since the standardizer should be able to deal with projects with empty rows
        IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
        for element in IFCB_empty_rows:
            if element in id_list:
                id_list.remove(element)
    os.mkdir(str(path_to_data) + '/binned_data/')


    for n, i in enumerate(id_list):
        df = read_func(path_to_data, i)# get a dataframe for each project
        # bin standardized data and obtain the metadata of the binned columns
        df, metadata_bins = bin_func(df, instrument, removecat, area_type = area_type, geom_shape = geom_shape, date_group = time_grouping,  ignore_depth = ignore_depth)
        path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
        # open the metadata of the standardized files
        with open(path_to_config, 'r') as config_file:
            cfg = yaml.safe_load(config_file)
        path_to_metadata = glob(str(Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir']) + '/**/*' + str(i) + '*metadata*')
        metadata_std = pd.read_csv(path_to_metadata[0], sep='\t', header=0, index_col=[0])
        #remove old biovolume descriptions in the metadata file
        metadata_std = metadata_std[ metadata_std['Variables'].str.contains('Biovolume')== False]
        #concatenate metadata files to generate binned metadata
        metadata_binned = pd.concat([metadata_std, metadata_bins], axis=0)
        NBS_data_binned = NB_SS_func(df, ignore_depth= ignore_depth )
        #generate metadata for the NBS files
        Variables = NBS_data_binned.columns.to_list()
        Variable_types = NBS_data_binned.dtypes.to_list()
        Units_Values = [time_grouping, 'lat_lon', 'cubic micrometers', 'cubic micrometers', 'degree', 'degree', '', 'cubic decimeters', 'cubic micrometers', '', 'cubic micrometers','counts/ cubic decimeters', 'log(counts/ cubic decimeters)', 'log (cubic micrometers)']
        Description = ['binned date information','string that serves as an identifier of a single cell of a  1x1 degree spatial grid','minimum and maximum value of the size bin, calculated from biovolume obtained using ' + area_type + ' and using a projection of ' + geom_shape,'difference between max and min value of a size bin','latitude of the center point of the 1x1 degree cell','longitude of the center point of the 1x1 degree cell','Project ID in Ecotaxa','Volume analyzed (not accounting for sample dilution and/or fractionation)','Sum of the biovolume of individual objects classified into a biovolume based size bin','number of objects assigned into the size bin','mean biovolume for each size bin','Normalized biomass size spectra based on biovolume','logarithmic transformation of the NBSS, use as y axis when performing size spectra analysis','logarithmic transformation of the median of a size bin, use as x axis when performing size spectra analysis']
        NBS_metadata = pd.DataFrame({'Variables': Variables, 'Variable_types': Variable_types,'Units_Values': Units_Values,'Description': Description})

        df.to_csv(str(path_to_data) + '/binned_data/' + str(i) + '_'+ instrument + '_binned.csv',  sep = '\t')
        metadata_binned.to_csv(str(path_to_data) + '/binned_data/' + str(i) + '_'+ instrument + '_binned_metadata.csv',  sep = '\t')
        NBS_data_binned.to_csv(str(path_to_data) + '/binned_data/' + str(i) + '_' + instrument + '_NBSS.csv', sep='\t')
        NBS_metadata.to_csv(str(path_to_data) + '/binned_data/' + str(i) + '_' + instrument + '_NBSS_metadata.csv', sep='\t')



# 7)calculate x and y for NBSS, this includes adding total biovolume per size bin for each station and depth bin,
# inputs: a dataframe and the volume sampled of each (in cubic meters). other inputs are the column names
# from the original dataset that want to be retained (i.e;  proj ID, station ID, depth, size classes,
# range of size classes, biovolume, lat & lon ) plus the variables that will be used to group the data by
# stations, depths and size classes
# using 5 ml for IFCB
def NB_SS_func(df, ignore_depth = 'no', dates='date_bin', station= 'Station_location', depths = 'midDepthBin',\
               size_range= 'range_size_bin', sizeClasses= 'sizeClasses', biovolume='Biovolume',\
               lat='midLatBin', lon='midLonBin', project_ID= 'Project_ID', vol_filtered='Volume_analyzed'):
    import numpy as np
    import pandas as pd
    if ignore_depth == 'no':
        # create a dataframe with summary statistics for each station, depth bin and size class
        # these column names should all be the same, since the input is a dataframe from the 'binning' and 'biovol' functions
        # group data by bins
        NBS_biovol_df = df.groupby([dates, station, sizeClasses, depths]).agg(
            {depths:'first', size_range:'first', lat: 'first', lon: 'first', project_ID: 'first', vol_filtered: 'first',
             biovolume:['sum', 'count', 'mean'] }).reset_index()

  # remove bins that have zero values
        # standardize by volume sample
    elif ignore_depth == 'yes':
        NBS_biovol_df = df.groupby([dates, station, sizeClasses]).agg(
            {size_range:'first', lat: 'first', lon: 'first', project_ID: 'first', vol_filtered: 'first',
             biovolume:['sum', 'count' , 'mean'] }).reset_index()

    NBS_biovol_df.columns = NBS_biovol_df.columns.map('_'.join).str.removesuffix("first")
    NBS_biovol_df.columns = NBS_biovol_df.columns.str.removesuffix("_")
    NBS_biovol_df = NBS_biovol_df[NBS_biovol_df['Biovolume_count'] != 0]
    NBS_biovol_df['NBSS'] = (NBS_biovol_df['Biovolume_sum'] / (NBS_biovol_df[vol_filtered])) / NBS_biovol_df[
        size_range]
    # create two more columns with the parameters of normalized size spectra,:
    NBS_biovol_df['logNBSS'] = np.log(NBS_biovol_df['NBSS'])
    NBS_biovol_df['logSize'] = np.log(NBS_biovol_df[size_range])

    return NBS_biovol_df



# 7)  create a function to remove data based on JO's comments and Haentjens THIS IS WORK IN PROGRESS:
# low threshold will be filtered because they are smaller than the next standardized size spectrum
# high threshold: obtain the highest size class, then iterate from highest to lowest and remove the large size bins
# that have less than an (instrument specific) threshold
# imputs: the data frame, and the column that contains the NBSS
def clean_lin_fit(binned_data, instrument, method = 'MAX'):
    """
    Objective: remove size bins based on JO's comments and Haentjens et al (2022)
    :param df: a binned dataframe with NBSS already calculated
    :return: a dataset without the bins that contain inaccurate NBSS
    """
    # step 1: 'We set the lower boundary of the spectra (bo) to be the operational lower size measured by the
    # flow cytometer (bo = 0.6 um in ESD) and the upper boundary (bn) as the largest cytobot size bins accross all
    # NAAMES samples that contained at least 10 particles (bn = 65 um in ESD)'
    #reverse_list = list(range(0, len(df))), reverse_list.sort(reverse=True)

    #upper limit for IFCB: the largest size bin with at least 10 cells
    if instrument  == 'IFCB':
        index_10 = binned_data.index[binned_data['Biovolume_count'] >= 10].tolist()
        binned_data_filt = binned_data.loc[0:max(index_10)]
        binned_data_filt = binned_data_filt.reset_index(drop=True)
        if method == 'MAX':
            list_max_NBSS = binned_data.NBSS.tolist()
            binned_data_filt = binned_data.loc[list_max_NBSS.index(max(list_max_NBSS)):len(binned_data)]
            binned_data_filt = binned_data_filt.reset_index(drop=True)
        elif method == 'instrument_specific': # lower limit for IFCB following Haentjens: size bin with the highest cell abundance among the <524 um3 size bins.
            index_524_sizes = binned_data_filt.index[binned_data_filt['Biovolume_mean'] < 524].tolist()
            count_maxSmBins = max(binned_data_filt.loc[index_524_sizes]['Biovolume_count'])
            lg524_binList = binned_data_filt.index[binned_data_filt['Biovolume_count'] == count_maxSmBins].tolist()
            binned_data_filt = binned_data_filt.loc[lg524_binList[0]:len(binned_data_filt)]
            binned_data_filt = binned_data_filt.reset_index(drop=True)
    #step 2: if a value is lower than the next one, remove that row, THIS DOES NOT WORK
    #low_bin_index = []
    #for i in range (1, len(binned_data_filt)):
        #if binned_data_filt.loc[i-1, 'NBSS'] < binned_data_filt.loc[i, 'NBSS']:
            #low_bin_index.append(i)
    #min_bin_index = []
    #for i in range(1, len(low_bin_index)):
        #if (low_bin_index[i] - low_bin_index[i-1]) > 1:
            #min_bin_index.append(low_bin_index[i])
    #binned_data_filt = binned_data_filt.drop(binned_data_filt.index[0:(min(min_bin_index)+1)])
    #binned_data_filt = binned_data_filt.reset_index(drop=True)
    # for zooscan:
    elif instrument == 'Zooscan':
        if method == 'MAX': # select the highest NBSS value following da Rocha Marcolin et al.
            list_max_NBSS = binned_data.NBSS.tolist()
            binned_data_filt = binned_data.loc[list_max_NBSS.index(max(list_max_NBSS)):len(binned_data)]
            binned_data_filt = binned_data_filt.reset_index(drop=True)
        elif method == 'instrument_specific': #adapted from Haentjens: theoretically, minimum sampled volume of an organism
            #with a WP2 is 4188790 um3 (volume of a sphere with 200 um diameter). So, size bin with the highest cell abundance
            #among the < 4188790
            index_4188790_sizes = binned_data.index[binned_data['Biovolume_mean'] < 4188790].tolist()
            count_maxSmBins = max(binned_data.loc[index_4188790_sizes]['Biovolume_count'])
            lg4188790_binList = binned_data.index[binned_data['Biovolume_count'] == count_maxSmBins].tolist()
            binned_data_filt = binned_data.loc[lg4188790_binList[0]:len(binned_data)]
            binned_data_filt = binned_data_filt.reset_index(drop=True)
        #elif binned_data_filt.loc[i, 'NBSS'] < binned_data_filt.loc[i+1, 'NBSS']:
             #subset_stDepth.drop(subset_stDepth.index[i:len(subset_stDepth[logNBSS]])
    #subset_stDepth = subset_stDepth.reset_index(drop=True)
    return binned_data_filt



def stitching_func(IFCB_proj_type= 'large_sample', method1= 'MAX', method2 = 'MAX'):
    path_to_binned_IFCB, ID_list_IFCB = proj_id_list_func('IFCB', standardized='binned', testing=False)
    path_to_binned_Zooscan, ID_list_Zooscan = proj_id_list_func('Zooscan', standardized='binned', testing=False)

    if IFCB_proj_type == 'large_sample':
        IFCB_list= [3302, 3309, 3321, 3325, 3334] # this should also handle proj 3324 (that project contains two stations)
    elif IFCB_proj_type == 'small_sample':
        IFCB_list= [3147, 2248, 3289, 3296]
    IFCB_list= [str(x) for x in IFCB_list]
    path_to_Zooscan = glob(str(path_to_binned_Zooscan) + '/' + str(378) + '*NBSS*') # for 3309, there is an issue in a single point where comparing NBSS between size bins does not work for filtering the shorter side of the spectrum
    IFCB_list= [str(x) for x in IFCB_list]
    IFCB_dict = {}
    for i in IFCB_list:
        path_to_IFCB = glob(str(path_to_binned_IFCB) + '/' + i + '*NBSS*')
        IFCB_dict[i] = pd.read_csv(path_to_IFCB[0], sep='\t', header=0, index_col=[0])
        IFCB_dict[i] = clean_lin_fit(IFCB_dict[i], instrument='IFCB', method= method1)
        binned_data_Zooscan = pd.read_csv(path_to_Zooscan[0], sep='\t', header=0, index_col=[0])
        subset_binned_Zooscan = binned_data_Zooscan.loc[(binned_data_Zooscan['Station_location'] == IFCB_dict[i].loc[0, 'Station_location'])]
        subset_binned_Zooscan = subset_binned_Zooscan.sort_values(by='range_size_bin')
        subset_binned_Zooscan = subset_binned_Zooscan.reset_index(drop=True)
        data_filt_Zooscan = clean_lin_fit(subset_binned_Zooscan, instrument='Zooscan', method= method2)
        IFCB_dict[i] = pd.concat([IFCB_dict[i], data_filt_Zooscan], axis=0)
        IFCB_dict[i] = IFCB_dict[i].reset_index(drop=True)
    Station_list = []
    for i in IFCB_list:
        Station_list.append(IFCB_dict[i].loc[0, 'Station_location'])
        station = (IFCB_dict[i].loc[0, 'Station_location'])
        IFCB_dict[station] = IFCB_dict.pop(i)

    fig, ax1 = plt.subplots(1,1)
    for st in Station_list:
        ax1.plot(IFCB_dict[st]['logSize'], IFCB_dict[st]['logNBSS'], marker='.', label=st)
    ax1.set_ylabel('log (NBS)')
    ax1.set_xlabel('log (Biovolume)')
    ax1.legend()
    ax1.title.set_text('NBSS for IFCB ' + method1 + ' and Zooscan ' + method2  + ' for ' + IFCB_proj_type)

# 8) perform linear regressions. Imput: a dataframe with biovolume and NB_SS
# (this can be modified based on Mathilde's comment
# to remove size bins based on  non linear least squares). In this example, the clean_lin_fit will be used for now
def NB_SS_regress_func(stats_biovol_SC, station, depths, lat, lon,  ID):
    import numpy as np
    import pandas as pd
    from scipy.stats import linregress
    # create empty lists that will be populated with the data of interest
    Middepth_range = []
    station_lat = []
    station_lon = []
    slopes = []
    intercept = []
    r_val = []
    p_val = []
    project_number = []
    # parse data by station
    for n, s in enumerate(stats_biovol_SC[station].unique()):
        subset_station = stats_biovol_SC[(stats_biovol_SC[station] == s)]
        subset_station = subset_station.reset_index(drop=True)
        # parse data by depth
        for o, d in enumerate(subset_station[depths].unique()):
            subset_stDepth = subset_station[(subset_station[depths] == d)]
            subset_stDepth = subset_station.reset_index(drop=True)
            # make the linear fit between normalized biovolume and size
            model = linregress(subset_stDepth['logSize'], subset_stDepth['logNBSS'])
            # populate the lists of lat long, and depth range
            Middepth_range.append(subset_stDepth[depths][0])
            station_lat.append(subset_stDepth[lat][0])
            station_lon.append(subset_stDepth[lon][0])
            # add the statistics to the empty lists
            slopes.append(model.slope)
            intercept.append(model.intercept)
            r_val.append(model.rvalue)
            p_val.append(model.pvalue)
        proj_number = [ID] * (n + o + 1)
    project_number = project_number + proj_number
    # compile all lists into a dataframe
    results_SS = pd.DataFrame(
        list(zip(project_number, station_lat, station_lon, Middepth_range, slopes, intercept, r_val, p_val)),
        columns=['project_id', 'station_lat', 'station_lon', 'Middepth_m', 'slope', 'intercept', 'r_val', 'p_val'])
    return results_SS




# 9) compile all regression results for all the projects, input: a list of projects, and use functions 2-7 to compile
# everything. This should be the only call, from telling it what instrument and then doing the rest.
# It needs to loop through the project ID and merge all the resulting dataframes at the end
def PSS(instrument, testing=False):  # highly redundant to function 2,
    # but I think it's still better to keep it separate
    import pandas as pd
    if testing == True:
        path_to_data, ID = proj_id_list(instrument, testing=True)
    else:
        path_to_data, ID = proj_id_list(instrument, testing=False)
    results_PSS_all = pd.DataFrame()
    for i in ID:
        data_clean = read_clean(path_to_data, i)
        data_clean = biovol(data_clean, instrument)
        data_clean = binning(data_clean)
        stats_biovol_SC = NB_SS(data_clean)
        results_PSS = NB_SS_regress(stats_biovol_SC, i)
        results_PSS_all = pd.concat([results_PSS_all, results_PSS], axis=0)
    results_PSS_all.to_csv(os.path.join(path_to_data, 'NPSS_IFCB.tsv'), sep='\t')
    return results_PSS_all
