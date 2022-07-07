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

#function to create size bins
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


# 1) function to create list of projects. Inputs: instrument
def proj_id_list_func(instrument, testing=False):
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
            id_list = [1]  # ask mathilde about her test projects
        elif instrument == 'UVP':
            id_list = [2]  # ask mathilde about her test projects
    return path_to_data, id_list

# 2) Standardize all files:
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
    if instrument == 'IFCB':# this removal is only for testing,
        # since the standardizer should be able to deal with projects with empty rows
        IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
        for element in IFCB_empty_rows:
            if element in ID_list:
                ID_list.remove(element)
        for i in ID_list:
            standardize(config_path, standardizer_path, i)

    else:
        for i in ID_list:
            standardize(config_path, standardizer_path, i)



# 3) Create clean dataframes in python, input: a path and an ID. modify this to prevent removing rows
def read_func(path_to_data, ID, cat='Category'):
    # returns a dataframe for each project, excluding datapoints with errors
    import pandas as pd
    from glob import glob
    search_pattern = path_to_data / ("*" + str(ID) + "*")
    filename = glob(str(search_pattern))
    df = pd.read_csv(filename[0], sep='\t', encoding='latin_1', encoding_errors='ignore', header=0)
    # convert taxonomic Id column into categorical.
    # NOTE all of these column names might need to be modified after standardization of datasets
    df[cat] = df[cat].astype('category')  # t
    df['proj_ID'] = [ID] * (len(df))
    return df

#function to remove rows without data
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

#3.1standardizer function should go here in the work flow (see outputs from functions below)

#3.2biovol_standardizer (should scaling happen here?). Imputs: column with size measurement (area or biovolume)
def biovol_standardizer_func(measurement, scale, instrument):
    if instrument == 'IFCB':
        biovol_um3 = measurement / (scale**3)
        # what should we do here for UVP and Zooscan
    return biovol_um3

#3.3 area_to_ESD: should scaling happen here?

#3.4 ESD_to_biovol: note: ESD should be in micrometers already
def ESD_to_biovol_func(ESD):
    biovol_um3 = (4 / 3) * 3.14 * ((ESD / 2) ** 2)

#5) data binning
    # by depth, size, lat/lon and size, input: a dataframe, lat/lon increments that define the stations, and a list
# of depth bins. NOTE: separate this into three functions
    #5.1 by size bins. Inputs: a column of the dataframe with biovolume, and the number of log bins to use.
    #returns two lists : the sizeClasses bins (categorical) and range_size_bins (float)
def size_binning_func(biovolume, N_log_bins=200):
    import numpy as np
    PD_SizeBins_um3 = log_bins_func((min(biovolume)), N_log_bins)  # number of bins work for this data
    SizeBins_um3 = np.array(PD_SizeBins_um3)
    # create a categorical variable, which are bins defined by the numbers on the sizeClasses list
    # Assign bin to each data point based on biovolume
    sizeClasses= pd.cut(x=biovolume,  # size classes defined by biovolume
                                       bins=SizeBins_um3, include_lowest=True)
    # obtain the size of each size class bin and append it to the stats dataframe
    # unfortunately the Interval type produced by pd.cut is hard to work with. So they need to be modified
    range_size_bin = []
    for i in range(0, len(biovolume)):
        range_size_bin.append(sizeClasses.iloc[i].length)
    return sizeClasses, range_size_bin

   # 5.3 by depth. Inputs: a column of the dataframe with depth data, and a list of depth bins
def depth_binning_func(depth='Depth_min', depth_bins=[0, 25, 50, 100, 200, 500, 1000, 3000, 8000]):
    # create depth bins based on Jessica's suggestion (Jessica's suggestions as default)
    depth_bin = pd.cut(x=depth, bins=depth_bins, include_lowest=True)
    # create a categorical variable with the middle depth of each bin
    midDepth_bin = []
    for i in range(0, len(depth)):
        midDepth_bin.append(depth_bin.iloc[i].mid)
    return midDepth_bin

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

# 6)calculate x and y for NBSS, this includes adding total biovolume per size bin for each station and depth bin,
# inputs: a dataframe and the volume sampled of each (in cubic meters). other inputs are the column names
# from the original dataset that want to be retained (i.e;  proj ID, station ID, depth, size classes,
# range of size classes, biovolume, lat & lon ) plus the variables that will be used to group the data by
# stations, depths and size classes
# using 5 ml for IFCB
def NB_SS_func(data_clean, station= 'Station_ID', depths = 'midDepth_bin',\
               size_range= 'range_size_bin', sizeClasses= 'sizeClasses', biovolume='Biovolume',\
               lat='midLat_bin', lon= 'midLon_bin', project_ID= 'proj_ID', vol_filtered=5e-6):
    import numpy as np
    import pandas as pd
    # create a dataframe with summary statistics for each station, depth bin and size class
    # these column names should all be the same, since the input is a dataframe from the 'binning' and 'biovol' functions
    # group data by bins
    stats_biovol_SC = data_clean[
        [station, depths, sizeClasses, size_range, biovolume, lat, lon, project_ID]] \
        .groupby([station, depths, sizeClasses, size_range, lat, lon, project_ID]).describe()
    # reset index and rename columns to facilitate further calculations
    stats_biovol_SC = stats_biovol_SC.reset_index()
    stats_biovol_SC.columns = stats_biovol_SC.columns.map('_'.join).str.strip('_')

    # add column of summed biovolume per size class to the stats_biovol_SC dataframe
    sum_biovol_SC = data_clean.groupby([station, depths, sizeClasses]).agg({biovolume: ['sum']})
    # reset index and rename columns to facilitate further calculations
    sum_biovol_SC = sum_biovol_SC.reset_index()
    sum_biovol_SC.columns = sum_biovol_SC.columns.map('_'.join).str.strip('_')
    sum_biovol_SC = sum_biovol_SC[sum_biovol_SC['biovol_um3_sum'] != 0]  # remove bins that have zero values
    sum_biovol_SC = sum_biovol_SC.reset_index(drop=True)
    stats_biovol_SC['sum_biovol'] = sum_biovol_SC['biovol_um3_sum']
    # standardize by volume sample
    stats_biovol_SC['NBSS'] = (stats_biovol_SC['sum_biovol'].div(vol_filtered)) / stats_biovol_SC[size_range]
    # ( based on Buonassisi and Diersen's (2010) cut of undersampled size bins with 10 particles per L= 10 000 particles/m3)
    stats_biovol_SC = stats_biovol_SC[stats_biovol_SC['NBSS'] > 10000]
    stats_biovol_SC = stats_biovol_SC.reset_index(drop=True)

    # create two more columns with the parameters of normalized size spectra,:
    stats_biovol_SC['logNBSS'] = np.log(stats_biovol_SC['NBSS'])
    stats_biovol_SC['logSize'] = np.log(stats_biovol_SC['range_size_bin'])

    return stats_biovol_SC

# 7)  create a function to remove data based on JO's comments THIS IS WORK IN PROGRESS:
# low threshold will be filtered because they are smaller than the next standardized size specttrum
# high threshold: still unsure. same concept? remove points that are higher than the previous one?
# imputs: the data frame, and the column that contains the NBSS
#def clean_lin_fit(subset_stDepth, logNBSS='logNBSS'):
    # step 1: if a value is lower than the next one, remove that row
    #for  i in  range  (1, len(subset_stDepth[logNBSS]):
        #if subset_stDepth[logNBSS][i-1] < subset_stDepth[logNBSS][i]:
            #subset_stDepth.drop([i-1])
        #elif subset_stDepth[logNBSS][i] < subset_stDepth[logNBSS][i+1]:
            #subset_stDepth.drop(subset_stDepth.index[i:len(subset_stDepth[logNBSS]])
    #subset_stDepth = subset_stDepth.reset_index(drop=True)
    #return subset_stDepth


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
