# create a python function that receives a folder as an imput (such folder contains csv files)
# another imput should be bin count
# and returns: slope, intercept and R^2 of the regression (do we want plots?)
# def SS_fit (pathname, n_bins):

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

# 1) function to create size bins
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


# 2) function to create list of projects. Inputs: instrument

def proj_id_list(instrument, testing= False ):
    import pandas as pd
    import yaml
    from pathlib import Path
    # returns path to project data and list of projects based on intrument (IFCB, Zooscan, UVP)
    # read git-tracked config file (text file) with inputs:  project ID, output directory
    path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    # read config file (text file) with inputs:  project ID, output directory
    path_to_git = Path(cfg['git_dir']).expanduser()
    # prepare storage based on path stored in the yaml config file and instrument type
    path_to_data = Path(cfg['git_dir']).expanduser() / cfg['dataset_subdir'] / instrument
    # create a list  projects that we have access to, based on project_list_all.xlsx
    path_to_proj_id_list = Path(cfg['git_dir']).expanduser() / cfg['proj_list']
    proj_list = pd.read_excel(path_to_proj_id_list)
    if testing == False:
        id_list = proj_list['Project_ID'].loc[
            (proj_list['Instrument'] == str(instrument)) &  # instrument type filtered here
            (proj_list['PSSdb_access'] == True)].tolist()
    else:
        if instrument == 'IFCB':
            id_list = [3315, 3318, 3326]  # THIS IS THE TEST LIST, USE THIS WHEN DEBUGGING for IFCB
        elif instrument == 'Zooscan':
            id_list = [1] # ask mathilde about her test projects
        elif instrument == 'UVP':
            id_list = [2] # ask mathilde about her test projects
    return path_to_data, id_list

#3) Create clean dataframes, input: a path and an ID
def read_clean(path_to_data, ID):
    #returns a dataframe for each project, excluding datapoints with errors
    import pandas as pd
    from glob import glob
    search_pattern = path_to_data / ("*" + str(ID) + "*")  # hard coded to find the ICFB folder
    filename = glob(str(search_pattern))
    df = pd.read_csv(filename[0], sep='\t', encoding='latin_1', encoding_errors='ignore', header=0)
    # convert taxonomic Id column into categorical.
    #NOTE all of these column names might need to be modified after standardization of datasets
    df['object_annotation_category'] = df['object_annotation_category'].astype('category')#t
    # FILTERING DATA: filter taxonomic categories that are artefacts
    data_clean = df[df['object_annotation_category'].str.contains("artefact") == False]
    # remove data points where size metrics = 0
    data_clean = data_clean[data_clean.object_equivdiameter != 0]
    # finally, remove any row with coordinates =nan
    data_clean = data_clean[data_clean['object_lat'].notna()]
    data_clean = data_clean[data_clean['object_lon'].notna()]
    return data_clean


#4) calculate biovolume of each object. input: a dataframe
def biovol(data_clean, instrument):
    import pandas as pd
    # these calculations should be done differently for each instrument, perhaps?
    if instrument == 'IFCB':
        data_clean['ESD_um'] = data_clean['object_equivdiameter'] / data_clean['process_pixel_to_micron']
        data_clean['biovol_um3'] = (4 / 3) * 3.14 * (
                (data_clean['ESD_um'] / 2) ** 2)  # based on sphere volume, likely better ways
        # to do this
    elif instrument == 'UVP':
        # a set of calculations here, repeated calculations as a placeholder
        data_clean['ESD_um'] = data_clean['object_equivdiameter'] / data_clean['process_pixel_to_micron']
        data_clean['biovol_um3'] = (4 / 3) * 3.14 * (
                (data_clean['ESD_um'] / 2) ** 2)  # based on sphere volume, likely better ways
    elif instrument == 'Zooscan':
        # a set of calculations here, repeated calculations as a placeholder
        data_clean['ESD_um'] = data_clean['object_equivdiameter'] / data_clean['process_pixel_to_micron']
        data_clean['biovol_um3'] = (4 / 3) * 3.14 * (
                (data_clean['ESD_um'] / 2) ** 2)  # based on sphere volume, likely better ways
    return data_clean

#5) bin by depth, size, lat/lon and size, input: a dataframe, lat/lon increments that define the stations, and a list
#of depth bins
def binning (data_clean, st_increment=1, depths = [0, 25, 50, 100, 200, 500, 1000, 3000, 8000]):
    import numpy as np
    import pandas as pd
    # create size bins, dependent of log_bins_func
    PD_SizeBins_um3 = log_bins_func((min(data_clean['biovol_um3'])), 50)  # number of bins work for this data
    SizeBins_um3 = np.array(PD_SizeBins_um3)
    # create a categorical variable, which are bins defined by the numbers on the sizeClasses list
    # Assign bin to each data point based on biovolume
    data_clean['sizeClasses'] = pd.cut(x=data_clean[('biovol_um3')],
                                       bins=PD_SizeBins_um3, include_lowest=True)
    # the dataframe now needs to be parsed, first by station and second, by depth ranges
    # create a categorical variable that will serve as the station ID based on binned lat and lon
    # 1) create arrays that will serve to bin lat and long with 1 degree increments
    lat_increments = np.arange(np.floor(min(data_clean['object_lat']) - st_increment),
                              np.ceil(max(data_clean['object_lat']) + st_increment), 1)
    lon_increments = np.arange(np.floor(min(data_clean['object_lon']) - st_increment),
                              np.ceil(max(data_clean['object_lon']) + st_increment), 1)
    # 2) create bins
    data_clean['lat_bin'] = pd.cut(x=data_clean[('object_lat')],
                                   bins=lat_increments, include_lowest=True)

    data_clean['lon_bin'] = pd.cut(x=data_clean[('object_lon')],
                                   bins=lon_increments, include_lowest=True)
    # 3) get middle lat after assigning each sample to a lat/lon bin
    midLat_bin = []
    for i in range(0, len(data_clean.index)):
        midLat_bin.append(data_clean['lat_bin'].iloc[i].mid)
    data_clean['midLat_bin'] = midLat_bin
    midLon_bin = []
    for i in range(0, len(data_clean.index)):
        midLon_bin.append(data_clean['lon_bin'].iloc[i].mid)
    data_clean['midLon_bin'] = midLon_bin
    # create categorical variable for station ID
    data_clean['station_ID'] = data_clean['midLat_bin'].astype(str) + '_' + data_clean['midLon_bin'].astype(str)
    data_clean['station_ID'] = data_clean['station_ID'].astype('category')

    # create depth bins based on Jessica's suggestion (Jessica's suggestions as default)
    data_clean['depth_bin'] = pd.cut(x=data_clean[('object_depth_min')],
                                     # min and max depths are the same in the IFCB data
                                     bins=depths, include_lowest=True)
    # create a categorical variable with the middle depth of each bin
    midDepth_bin = []
    for i in range(0, len(data_clean.index)):
        midDepth_bin.append(data_clean['depth_bin'].iloc[i].mid)
    data_clean['midDepth_bin'] = midDepth_bin
    data_clean['midDepth_bin'] = data_clean['midDepth_bin'].astype('category')

    return data_clean


#6)calculate NBSS,this includes adding total biovolume per size bin for each station and depth bin

#7) perform linear regression, input: a dataframe
#8) compile all regression results for all the projects, input: a list of projects, and use functions 3-6 to compile
# everything


# create lists for the statistics of the Size spectrum analysis
project_number = []
slopes = []
intercept = []
r_val = []
p_val = []
# create lists for sampling sites and depths
Middepth_range = []
station_lat = []
station_lon = []

# create list of size bins based on Ecopart, converted to um
EcopartSizeBins_mm = np.array(
    [0.001000, 0.001260, 0.001590, 0.002000, 0.002520, 0.003170, 0.004000, 0.005040, 0.006350, 0.008000
        , 0.0101, 0.0127, 0.0160, 0.0202, 0.0254, 0.0320, 0.0403, 0.0508, 0.0640, 0.0806
        , 0.1020, 0.1280, 0.1610, 0.2030, 0.2560, 0.3230, 0.4060, 0.5120, 0.6450, 0.8130
        , 1.0200, 1.2900, 1.6300, 2.0500, 2.5800, 3.2500, 4.1000, 5.1600, 6.5000, 8.1900
        , 10.3000, 13.0000, 16.4000, 20.6000, 26.0000, 10E7])

EcopartSizeBins_um = np.array([bin * 1000 for bin in EcopartSizeBins_mm])



# Size spectra calculations for each project starts here
for id in id_list:

    # scale the dimensions of the particles to metric units and append to the dataframe
    # THIS SECTION NEEDS TO BE COMPARED TO THE IFCB LITERATURE and/or Mathilde could advise us on this
    data_clean['ESD_um'] = data_clean['object_equivdiameter'] / data_clean['process_pixel_to_micron']
    data_clean['biovol_um3'] = (4 / 3) * 3.14 * (
            (data_clean['ESD_um'] / 2) ** 2)  # based on sphere volume, likely better ways
    # to do this

    # create the "characteristic size  of each size class"
    # three ways of doing this:  (1) by applying a log scale to the size range of the particles
    # sizeClassesLog = np.linspace(np.log(data_clean.ESD_um.min() - 0.0000001),
    # np.log(data_clean.ESD_um.max() + 0.0000001),
    # endpoint=True, num=12)  # num parameter defines the number of bins,
    # how many bins do we want with this method?
    # remove log scale to end up with the final array of size classes
    # sizeBins_um= np.exp(sizeClassesLog)

    # (2) Use the fixed Ecopart size classes defined on lines 57-62
    # sizeBins_um= EcopartSizeBins_um
    # (3) use the harmonic scaling proposed by Platt and Denman (1977). Method used here


    # create a dataframe with summary statistics for each station, depth bin and size class
    stats_biovol_SC = data_clean[
        ['station_ID', 'midDepth_bin', 'sizeClasses', 'biovol_um3', 'midLat_bin', 'midLon_bin']] \
        .groupby(['station_ID', 'midDepth_bin', 'sizeClasses']).describe()
    stats_biovol_SC = stats_biovol_SC.reset_index()
    stats_biovol_SC.columns = stats_biovol_SC.columns.map('_'.join).str.strip('_')

    # add column of summed biovolume per size class to the stats_biovol_SC dataframe (needs to be different lat long)
    sum_biovol_SC = data_clean.groupby(['station_ID', 'midDepth_bin', 'sizeClasses']).agg({'biovol_um3': ['sum']})
    sum_biovol_SC = sum_biovol_SC.reset_index()
    sum_biovol_SC.columns = sum_biovol_SC.columns.map('_'.join).str.strip('_')
    sum_biovol_SC = sum_biovol_SC[sum_biovol_SC['biovol_um3_sum'] != 0]
    sum_biovol_SC = sum_biovol_SC.reset_index(drop=True)
    stats_biovol_SC['sum_biovol'] = sum_biovol_SC['biovol_um3_sum']

    # obtain the size of each size class bin and append it to the stats dataframe
    # unfortunately the Interval type produced by pd.cut is hard to work with. So they need to be modified
    range_size_bin = []
    for i in range(0, len(stats_biovol_SC.index)):
        range_size_bin.append(stats_biovol_SC['sizeClasses'].iloc[i].length)
    stats_biovol_SC['range_size_bin'] = range_size_bin

    # create two more columns with the parameters of normalized size spectra:
    stats_biovol_SC['logBm/size'] = np.log(stats_biovol_SC['sum_biovol'] / stats_biovol_SC['range_size_bin'])
    stats_biovol_SC['logSize'] = np.log(stats_biovol_SC['range_size_bin'])

    # change depth and station ID to strings so that they are easily indexed
    # stats_biovol_SC['midDepth_bin'] = stats_biovol_SC.midDepth_bin.astype(str)
    # stats_biovol_SC['station_ID'] = stats_biovol_SC.station_ID.astype(str)

    # do fitting for each station and for each depth bin. I could only figure out how to parse the data by stations
    # even then, an issue to resolve is to determine the cutoff between stations of lat and long
    # any help in also parsing data by 'midDepth_bin' is appreciated

    for n, s in enumerate(stats_biovol_SC.station_ID.unique()):
        subset_station = stats_biovol_SC[(stats_biovol_SC['station_ID'] == s)]
        subset_station = subset_station.reset_index(drop=True)
        for o, d in enumerate(subset_station.midDepth_bin.unique()):
            subset_stDepth = subset_station[(subset_station['midDepth_bin'] == d)]
            subset_stDepth = subset_station.reset_index(drop=True)
            # populate the lists of lat long, and depth range. Lots of work to be done here
            Middepth_range.append(subset_stDepth['midDepth_bin'][0])
            station_lat.append(subset_stDepth['midLat_bin_mean'][0])
            station_lon.append(subset_stDepth['midLon_bin_mean'][0])
            # make the linear fit between normalized biovolume and size
            model = linregress(subset_stDepth['logSize'], subset_stDepth['logBm/size'])
            # add the statistics to the empty lists
            slopes.append(model.slope)
            intercept.append(model.intercept)
            r_val.append(model.rvalue)
            p_val.append(model.pvalue)
        proj_number = [id] * (n + o + 1)
    project_number = project_number + proj_number
# compile all the lists of values into a dataframe. This is still missing depth bins,
results_SS = pd.DataFrame(
    list(zip(project_number, station_lat, station_lon, Middepth_range, slopes, intercept, r_val, p_val)),
    columns=['project_id', 'station_lat', 'station_lon', 'Middepth_m', 'slope', 'intercept', 'r_val', 'p_val'])

results_SS.to_csv(os.path.join(path_to_data, 'NPSS_IFCB.tsv'), sep='\t')
