# functions to calculate normalized biomass size spectra, based on ecopart size bins and a logarithmic scale
# data needs to be gridded and temporally and spatially binned. proj_nbs_func summarizes all the steps

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

from pathlib import Path
from glob import glob
# Config modules
import yaml  # requires installation of PyYAML package

import read_funcs as rd #functions to list standardized and binned files

  #6.1 by size bins. Inputs: a column of the dataframe with biovolume, and the number of log bins to use.
    #returns two lists : the sizeClasses bins (categorical) and range_size_bins (float)
def size_binning_func(df_subset):
    """
    Objective: bin the dataframes (parsed by stations and/or dates and/or depths)
    :param biovolume: column of a dataframe that contains the biovolume (in cubic micrometers)
    :return: a binned dataframe by sizes, containing the range of the size bin for each
    """
    import numpy as np
    # open dataframe with the size bins
    path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
    # open the metadata of the standardized files
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    path_to_bins = str(Path(cfg['raw_dir']).expanduser() / 'ecopart_size_bins.txt')
    bins_df = pd.read_csv(path_to_bins, sep='\t')
    # create a categorical variable, which are bins defined by the numbers on the sizeClasses list
    # Assign bin to each data point based on biovolume, append bin range and bin class to the dataframe
    df_subset.loc[:, 'sizeClasses']= pd.cut(x=df_subset['Biovolume'], bins=bins_df['biovol_um3'], include_lowest=True)# size classes defined by biovolume
    range_size_bin = []
    for i in range(0, len(df_subset['sizeClasses'])):
        range_size_bin.append(df_subset.loc[i, 'sizeClasses'].length)
    df_subset.loc[:, 'range_size_bin'] = range_size_bin
    bins = df_subset['sizeClasses'].unique()
    bin_width = []# list that will contain the width of each size_bin
    for i in range(0, len(bins.categories)):
        bin_width.append(bins.categories[i].length)

    d = {'sizeClasses':bins.categories, 'range_size_bin': bin_width}
    df_bins = pd.DataFrame(d)
    # obtain the size of each size class bin and append it to the stats dataframe
    # unfortunately the Interval type produced by pd.cut is hard to work with. So they need to be modified

    return df_subset, df_bins


def parse_NBS_func(df, parse_by = ['Station_location', 'date_bin', 'midDepthBin']):
    """
    Objective: parse out any df that come from files that include multiple stations, dates and depths.
    dataframe needs to be already binned.
    :param df: a standardized, binned dataframe
    :param parse_by: list of columns that will be use][-=d to parse the dataframe
    :return: a binned, tresholded NBS dataframe
    """
    NBSS_binned_all = pd.DataFrame()
    for st in set(df[parse_by[0]]):
        df_subset = df[df[parse_by[0]] == st].reset_index(drop=True)
        df_binned, df_bins = size_binning_func(df_subset)  # create size bins
        NBS_biovol_df = NB_SS_func(df_binned, df_bins, ignore_depth='yes', dates='date_bin', station='Station_location',
                               depths='midDepthBin', size_range='range_size_bin', sizeClasses='sizeClasses',
                               biovolume='Biovolume', lat='midLatBin', lon='midLonBin',
                               project_ID='Project_ID', vol_filtered='Volume_analyzed')  # calculate NBSS
        NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])

    return NBSS_binned_all

# 7)calculate x and y for NBSS, this includes adding total biovolume per size bin for each station and depth bin,
# inputs: a dataframe and the volume sampled of each (in cubic meters). other inputs are the column names
# from the original dataset that want to be retained (i.e;  proj ID, station ID, depth, size classes,
# range of size classes, biovolume, lat & lon ) plus the variables that will be used to group the data by
# stations, depths and size classes
# using 5 ml for IFCB
def NB_SS_func(df, df_bins, ignore_depth = 'no', dates='date_bin', station= 'Station_location', depths = 'midDepthBin',\
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
    NBS_biovol_df['NBSS'] = (NBS_biovol_df['Biovolume_sum']) / (np.nansum(list(set(NBS_biovol_df[vol_filtered])))) / (NBS_biovol_df[size_range])
    NBS_biovol_df[vol_filtered] = np.nansum(list(set(NBS_biovol_df[vol_filtered]))) # this is in the case there are different volumes analyzed based on binning, they have to be summed
    # create two more columns with the parameters of normalized size spectra,:
    NBS_biovol_df['logNBSS'] = np.log(NBS_biovol_df['NBSS'])
    NBS_biovol_df['logSize'] = np.log(NBS_biovol_df[size_range])

    # now calculate the log biovolume for the df_bins dataframe
    df_bins['logSize'] = np.log(df_bins[size_range])
    # merge the two dataframes
    binned_NBS= pd.merge(df_bins, NBS_biovol_df, how='left', on=[sizeClasses, size_range, 'logSize'])
    # now fill the columns of date, station, lat/lon, project ID and volume
    for i in [dates, station, lat, lon, project_ID, vol_filtered]:
        binned_NBS[i] = binned_NBS[i].fillna(NBS_biovol_df[i].unique()[0])
    # let's do the thresholding here:
    NBS_binned_thres = threshold_func(binned_NBS)

    return NBS_binned_thres

# 7)  create a function to remove data based on JO's comments and Haentjens THIS IS WORK IN PROGRESS:
# low threshold will be filtered because they are smaller than the next standardized size spectrum
# high threshold: obtain the highest size class, then iterate from highest to lowest and remove the large size bins
# that have less than an (instrument specific) threshold
# imputs: the data frame, and the column that contains the NBSS
def threshold_func(binned_data):
    """
    Objective: remove size bins based on Fabien's approach: (lower limit: max NBSS, upper limit: 3 NBSS with Nans)
    :param df: a binned dataframe with NBSS already calculated
    :return: a dataset without the bins that contain inaccurate NBSS
    """
    import numpy as np
    import pandas as pd

    #upper threshold based on max NBSS
    list_max_NBSS = binned_data.NBSS.tolist()
    binned_data_filt = binned_data.loc[list_max_NBSS.index(np.nanmax(list_max_NBSS)):len(binned_data)] #a ask if its better to just remove all data
    binned_data_filt = binned_data_filt.reset_index(drop=True)

    #lower threshold based on three consecutive size bins with nans
    for n, i in enumerate(np.isnan(binned_data_filt['NBSS'])):
        if (i == True) and (np.isnan(binned_data_filt.loc[n+1, 'NBSS']) == True) and (np.isnan(binned_data_filt.loc[n+2, 'NBSS']) == True):
            binned_data_filt = binned_data_filt.loc[0:n-1]
            break

    return binned_data_filt



def proj_NBS_func(instrument,  ignore_depth = 'yes'):
    """
    Objective: read into the standardized, gridded tsv files (output from ) and
    :param instrument: the device used for image adcquisition. important since the path will change
    :return: tsv files with the same data but binned (would we like to delete the standardized data?)
    """
    import os
    import shutil
    path_to_data, id_list = rd.proj_id_list_func(instrument, data_status ='gridded')#generate path and project ID's
    if instrument == 'IFCB':# this removal is only for testing,
        #since the standardizer should be able to deal with projects with empty rows
        IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
        for element in IFCB_empty_rows:
            if element in id_list:
                id_list.remove(element)
    dirpath = str(path_to_data) + '/NBSS_data/'
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    for i in id_list:
        df = rd.read_func(path_to_data, i)# get a dataframe for each project
        # bin standardized data and obtain the metadata of the binned columns
        path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
        # open the metadata of the standardized files
        with open(path_to_config, 'r') as config_file:
            cfg = yaml.safe_load(config_file)
        path_to_metadata = glob(str(Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir']) + '/**/*' + str(i) + '*metadata*')
        metadata_std = pd.read_csv(path_to_metadata[0], sep='\t', header=0, index_col=[0])
        #remove old biovolume descriptions in the metadata file
        metadata_std = metadata_std[ metadata_std['Variables'].str.contains('Biovolume')== False]
        #concatenate metadata files to generate binned metadata
        #metadata_binned = pd.concat([metadata_std, metadata_bins], axis=0)
        # and generate the NBSS WITH THRESHOLDING INCLUDED
        if ignore_depth == 'no':
            NBS_data_binned = parse_NBS_func(df, parse_by=['Station_location', 'date_bin', 'midDepthBin'])
        else:
            NBS_data_binned = parse_NBS_func(df, parse_by=['Station_location', 'date_bin'])

        #generate metadata for the NBS files
        Variables = NBS_data_binned.columns.to_list()
        Variable_types = NBS_data_binned.dtypes.to_list()
        Units_Values = ['refer to gridded metadata', 'lat_lon', 'cubic micrometers', 'cubic micrometers', 'degree', 'degree', '', 'cubic decimeters', 'cubic micrometers', '', 'cubic micrometers','counts/ cubic decimeters', 'log(counts/ cubic decimeters)', 'log (cubic micrometers)'] # notice that we dont have date unit info here, fix this eventually
        Description = ['binned date information','string that serves as an identifier of a single cell of a  1x1 degree spatial grid','minimum and maximum value of the size bin, calculated from biovolume obtained using ' + 'refer to gridded metadata' + ' and using a projection of a sphere','difference between max and min value of a size bin','latitude of the center point of the 1x1 degree cell','longitude of the center point of the 1x1 degree cell','Project ID in Ecotaxa','Volume analyzed (not accounting for sample dilution and/or fractionation)','Sum of the biovolume of individual objects classified into a biovolume based size bin','number of objects assigned into the size bin','mean biovolume for each size bin','Normalized biomass size spectra based on biovolume','logarithmic transformation of the NBSS, use as y axis when performing size spectra analysis','logarithmic transformation of the median of a size bin, use as x axis when performing size spectra analysis']
        NBS_metadata = pd.DataFrame({'Variables': Variables, 'Variable_types': Variable_types,'Units_Values': Units_Values,'Description': Description})

        NBS_data_binned.to_csv(str(path_to_data) + '/NBSS_data/' + str(i) + '_' + instrument + '_NBSS.csv', sep='\t')
        NBS_metadata.to_csv(str(path_to_data) + '/NBSS_data/' + str(i) + '_' + instrument + '_NBSS_metadata.csv', sep='\t')

### NOTE: 'refer to gridded metadata' units in the metadata needs to be replaced, perhaps it might be worth keeping everythin as it was before (gridded and NBS files produced in one call)

