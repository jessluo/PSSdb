# functions to calculate normalized biomass size spectra, based on ecopart size bins and a logarithmic scale
# data needs to be gridded and temporally and spatially binned. proj_nbs_func summarizes all the steps
import numpy as np
from pathlib import Path
import pandas as pd
import yaml
try:
    from funcs_read import *
    from funcs_gridding import *
except:
    from scripts.funcs_read import *
    from scripts.funcs_gridding import *

import os
from tqdm import tqdm

def group_gridded_files_func(instrument, already_gridded= 'N'):
    """
    Objective: assign a label to each gridded dataset, so when the compilation of files occurs, it is only done with a group of data and prevents
    the creation of a huge file
    """
    #first, create series of numbers that break the globe into 15x15 degree cells:
    lat_left = np.arange(-90, 90+15, 15, dtype=int)
    lat_right = np.arange(-75,105 + 15, 15, dtype=int)
    lat_int = pd.IntervalIndex.from_arrays(lat_left, lat_right)
    lon_left = np.arange(-180, 180+15, 15, dtype=int)
    lon_right = np.arange(-165, 195+15, 15, dtype=int)
    lon_int = pd.IntervalIndex.from_arrays(lon_left, lon_right)
    grid_list = []
    # create a dictionary that has unique keys for each grid interval
    #grid_dict = {}
    #for x in range(0, len(lon_int)):
        #for y in range(0, len(lat_int)):
            #grid_dict[(str(x) + '_' + str(y))] = {'lon': lon_int[x], 'lat': lat_int[y]} ## however, this might not be necessary
    # try now to assing lat & lon grid numbers to the dataframe directly
    file_list = proj_id_list_func(instrument, data_status='gridded')
    for i in tqdm(file_list):
        if already_gridded == 'N':
            print('generating subsets of data for 15 degree grids from ' + i)
        else:
            print('extracting big grid labels for subsetting from ' + i)
        filename = i.split('/')[-1]
        dirpath = i.replace(filename, '')
        df = pd.read_csv(i, header = 0)
        df['lat_grid'] = None
        df['lon_grid'] = None
        for la in df['midLatBin'].unique():
            df.loc[df['midLatBin'] == la, 'lat_grid']= int(np.where(lat_int.contains(la) == True)[0])
        for lo in df['midLonBin'].unique():
            df.loc[df['midLonBin'] == lo, 'lon_grid'] = int(np.where(lon_int.contains(lo) == True)[0])
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



  #6.1 by size bins. Inputs: a column of the dataframe with biovolume, and the number of log bins to use.
    #returns two lists : the sizeClasses bins (categorical) and range_size_bins (float)
def size_binning_func(df_subset):
    """
    Objective: bin the dataframes (parsed by stations and/or dates and/or depths)
    :param biovolume: column of a dataframe that contains the biovolume (in cubic micrometers)
    :return: a binned dataframe by sizes, containing the range of the size bin for each
    """

    # open dataframe with the size bins
    path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
    # open the metadata of the standardized files
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    path_to_bins = str(Path(cfg['raw_dir']).expanduser() / 'ecopart_size_bins.tsv')
    bins_df = pd.read_csv(path_to_bins)
    # create a categorical variable, which are bins defined by the numbers on the sizeClasses list
    # Assign bin to each data point based on biovolume, append bin range and bin class to the dataframe
    df_subset.loc[:, 'sizeClasses']= pd.cut(x=df_subset['Biovolume'], bins=bins_df['biovol_um3'], include_lowest=True)# size classes defined by biovolume
    df_subset = df_subset.dropna(subset=['sizeClasses']).reset_index(drop=True) # hard pass, debug this
    range_size_bin = []
    size_class_mid = []
    for i in range(0, len(df_subset['sizeClasses'])):
        range_size_bin.append(df_subset.loc[i, 'sizeClasses'].length)
        size_class_mid.append(df_subset.loc[i, 'sizeClasses'].mid)
    df_subset.loc[:, 'range_size_bin'] = range_size_bin
    df_subset.loc[:, 'size_class_mid'] = size_class_mid
    bins = df_subset['sizeClasses'].unique()
    bin_width = []# list that will contain the width of each size_bin
    bin_mid = []
    for i in range(0, len(bins.categories)):
        bin_width.append(bins.categories[i].length)
        bin_mid.append(bins.categories[i].mid)
    d = {'sizeClasses':bins.categories, 'range_size_bin': bin_width,'size_class_mid': bin_mid}
    df_bins = pd.DataFrame(d)
    # obtain the size of each size class bin and append it to the stats dataframe
    # unfortunately the Interval type produced by pd.cut is hard to work with. So they need to be modified

    return df_subset, df_bins


# 7)calculate x and y for NBSS, this includes adding total biovolume per size bin for each station and depth bin,
# inputs: a dataframe and the volume sampled of each (in cubic meters). other inputs are the column names
# from the original dataset that want to be retained (i.e;  proj ID, station ID, depth, size classes,
# range of size classes, biovolume, lat & lon ) plus the variables that will be used to group the data by
# stations, depths and size classes
# using 5 ml for IFCB
def NB_SS_func(df_binned, df_bins, light_parsing = False, depth_parsing = False,thresholding=True):
    """
    """
    import numpy as np
    import pandas as pd
    if depth_parsing == True:
        # create a dataframe with summary statistics for each station, depth bin and size class
        # these column names should all be the same, since the input is a dataframe from the 'binning' and 'biovol' functions
        # group data by bins
        #test = df_binned.groupby([sizeClasses])

        df_binned['Biovolume'] = df_binned['Biovolume'] * df_binned['ROI_number']
        NBS_biovol_df = df_binned.groupby(['date_bin', 'Station_location', 'light_cond', 'midDepthBin', 'sizeClasses']).agg(
            {'Sample':'first', 'midDepthBin':'first', 'range_size_bin':'first', 'size_class_mid':'first', 'midLatBin': 'first', 'midLonBin': 'first', 'Volume_imaged': 'first', 'Min_obs_depth': 'first', 'Max_obs_depth':'first',
             'Biovolume':['sum'], 'ROI_number':['sum']}).reset_index() #roi_n:['sum']

        #average_biovol = df_binned.groupby([dates, station, light, depths]).apply(lambda x: pd.Series().reset_index(drop=True))

        df_vol = df_binned.groupby(['date_bin', 'Station_location', 'light_cond', 'midDepthBin']).apply(lambda x: pd.Series({'cumulative_volume': x[['Sample', 'Min_obs_depth', 'Max_obs_depth','Volume_imaged']].drop_duplicates().Volume_imaged.sum()})).reset_index(drop=True)

        # remove bins that have zero values
        # standardize by volume sample
    else:
        if light_parsing == True:
            df_binned['Biovolume'] = df_binned['Biovolume'] * df_binned['ROI_number']
            NBS_biovol_df = df_binned.groupby(['date_bin', 'Station_location', 'light_cond', 'sizeClasses']).agg(
                {'Sample':'first', 'range_size_bin':'first', 'size_class_mid':'first', 'midLatBin': 'first', 'midLonBin': 'first', 'Volume_imaged': 'first', min_depth: 'first', max_depth:'first',
                 'Biovolume':['sum'], 'ROI_number':['sum']}).reset_index()
        # cumulative volume for 1x1 degree bin. For UVP: volume imaged is the volume of a picture. Consider sample ID and depths, to be abble to get the cumulative volume imaged for a grid
            df_vol = df_binned.groupby(['date_bin', 'Station_location', 'light_cond']).apply(lambda x: pd.Series({ 'cumulative_volume': x[['Sample', 'Min_obs_depth', 'Max_obs_depth', 'Volume_imaged']].drop_duplicates().Volume_imaged.sum()})).reset_index(drop=True) # don't understand why sampl
        else:
            df_binned['Biovolume'] = df_binned['Biovolume'] * df_binned['ROI_number']
            NBS_biovol_df = df_binned.groupby(['date_bin', 'Station_location', 'sizeClasses']).agg(
                {'Sample': 'first', 'range_size_bin': 'first', 'size_class_mid':'first', 'midLatBin': 'first', 'midLonBin': 'first', 'Volume_imaged': 'first', 'Min_obs_depth': 'first', 'Max_obs_depth': 'first',
                 'Biovolume': ['sum'], 'ROI_number': ['sum']}).reset_index()

            # cumulative volume for 1x1 degree bin. For UVP: volume imaged is the volume of a picture. Consider sample ID and depths, to be abble to get the cumulative volume imaged for a grid
            df_vol = df_binned.groupby(['date_bin', 'Station_location']).apply(lambda x: pd.Series({'cumulative_volume': x[['Sample', 'Min_obs_depth', 'Max_obs_depth','Volume_imaged']].drop_duplicates().Volume_imaged.sum()})).reset_index(drop=True)  # don't understand why sample:


    NBS_biovol_df['cumulative_vol'] =pd.merge(NBS_biovol_df, df_vol,how='left',on=['date_bin', 'Station_location'])['cumulative_vol']# df_vol.loc[0, 'cumulative_volume']


    NBS_biovol_df.columns = NBS_biovol_df.columns.map('_'.join).str.removesuffix("first")
    NBS_biovol_df.columns = NBS_biovol_df.columns.str.removesuffix("_")
    NBS_biovol_df['Biovolume_mean'] = NBS_biovol_df['Biovolume_sum']/ NBS_biovol_df['ROI_number_sum']
    #NBS_biovol_df = NBS_biovol_df[NBS_biovol_df['Biovolume_count'] != 0]
    NBS_biovol_df['NB'] = (NBS_biovol_df['Biovolume_sum']) / (NBS_biovol_df['cumulative_vol']) / (NBS_biovol_df['range_size_bin'])

    #two lines below deprecated as of 1/31/2023, use new biovolume:
    #NBS_biovol_df['NBSS'] = (NBS_biovol_df['Biovolume_sum']) / (np.nansum(list(set(NBS_biovol_df[vol_filtered])))) / (NBS_biovol_df[size_range])
    #NBS_biovol_df[vol_filtered] = np.nansum(list(set(NBS_biovol_df[vol_filtered]))) # this is in the case there are different volumes analyzed based on binning, they have to be summed

    # create two more columns with the parameters of normalized size spectra,:
    NBS_biovol_df['logNB'] = np.log(NBS_biovol_df['NB'])
    NBS_biovol_df['logSize'] = np.log(NBS_biovol_df['size_class_mid'])

    # now calculate the log biovolume for the df_bins dataframe. NOTE updated on 2/23/2023. middle point of the size class
    df_bins['logSize'] = np.log(df_bins['size_class_mid'])
    # merge the two dataframes
    binned_NBS= pd.merge(df_bins, NBS_biovol_df, how='left', on=['sizeClasses', 'size_class_mid', 'logSize'])
    # now fill the columns of date, station, lat/lon, project ID and volume
    for i in ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Volume_imaged']:
        try:
            if (np.isnan(NBS_biovol_df[i].unique()).any()) == True:
                val =np.unique(NBS_biovol_df[i][~np.isnan(NBS_biovol_df[i])])[0]
        except TypeError:
            val = NBS_biovol_df[i].unique()[0]
        binned_NBS[i] = binned_NBS[i].fillna(val)
    # let's do the thresholding here:
    if thresholding:
        NBS_binned_thres = threshold_func(binned_NBS)
    else:
        NBS_binned_thres = binned_NBS
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
    list_max_NBSS = binned_data.NB.tolist()
    binned_data_filt = binned_data.loc[list_max_NBSS.index(np.nanmax(list_max_NBSS)):len(binned_data)] #a ask if its better to just remove all data
    binned_data_filt = binned_data_filt.reset_index(drop=True)

    #lower threshold based on three consecutive size bins with nans
    for n, i in enumerate(np.isnan(binned_data_filt['NB'])):
        if (i == True) and (np.isnan(binned_data_filt.loc[n+1, 'NB']) == True) and (np.isnan(binned_data_filt.loc[n+2, 'NB']) == True):
            binned_data_filt = binned_data_filt.loc[0:n-1]
            break

    return binned_data_filt

def linear_fit_func(df1, light_parsing = False, depth_parsing = False):
    """
    Objective: perform linear fit of NBSS. method of linear fit might change, current one is Least squares regression
    :param df: a dataframe with size bins and NBSS belonging to only one sample (one unique station, date and depth)
    from https://www.edureka.co/blog/least-square-regression/
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    #remove Nans
    df1 = df1.dropna()

    # extract x and y
    X = df1['logSize'].values#.reshape(-1, 1)
    Y = df1['logNB'].values#.reshape(-1, 1)
    # Mean X and Y
    mean_x = np.mean(X)
    mean_y = np.mean(Y)
    # Total number of values
    n = len(X)
    # Model 1 -linear model : Y = m*X + b
    numer = 0 # numerator
    denom = 0 # denominator
    for i in range(n):
        numer += (X[i] - mean_x) * (Y[i] - mean_y)
        denom += (X[i] - mean_x) ** 2
    m = numer / denom # slope
    b = mean_y - (m * mean_x) # intercept
    #print("Coefficients")
    #print(m, c)
    # Calculating Root Mean Squares Error and R2 score
    rmse = 0 # root mean square error
    ss_tot = 0 # total sum of squares
    ss_res = 0 # total sum of squares of residuals
    for i in range(n):
        y_pred = b + m * X[i]
        rmse += (Y[i] - y_pred) ** 2
        ss_tot += (Y[i] - mean_y) ** 2
        ss_res += (Y[i] - y_pred) ** 2
    rmse = np.sqrt(rmse / n)
    #print("RMSE")
    #print(rmse)
    R2 = 1 - (ss_res / ss_tot)
    #print("R2 Score")
    #print(R2)

    # generate dataframe and append the results in the corresponding variable
    lin_fit = pd.DataFrame()
    if depth_parsing == True:
        lin_fit.loc[0, 'depth'] = df1.loc[0, 'midDepthBin']

    lin_fit.loc[0, 'Station_location'] = str(df1.loc[0, 'Station_location'])
    lin_fit.loc[0, 'Date'] = str(df1.loc[0, 'date_bin'])

    #lin_fit.loc[0, 'year'] = str(df1.loc[0, 'year'])
    #lin_fit.loc[0, 'month'] = str(df1.loc[0, 'month'])

    lin_fit.loc[0, 'latitude'] = df1.loc[0, 'midLatBin']
    lin_fit.loc[0, 'longitude'] = df1.loc[0, 'midLonBin']
    if light_parsing==True:
        lin_fit.loc[0, 'light_condition'] = df1.loc[0, 'light_cond']
    lin_fit.loc[0, 'min_depth'] = df1.loc[0, 'Min_obs_depth']
    lin_fit.loc[0, 'max_depth'] = df1.loc[0, 'Max_obs_depth']

    lin_fit.loc[0, 'slope'] = m
    lin_fit.loc[0, 'intercept'] = b
    lin_fit.loc[0, 'rmse'] = rmse
    lin_fit.loc[0, 'r2'] = R2

    return (lin_fit)

def parse_NBS_linfit_func(df, instrument, parse_by=['Station_location', 'date_bin'], light_parsing=False, depth_parsing=False, bin_loc = 1, group_by = 'yyyymm'):
    """
    Objective: parse out any df that come from files that include multiple stations, dates and depths.
    dataframe needs to be already binned. Steps are: 1) parse the data 2) bin ROIs by size  3) calculate the NBS
    (threshold included) and 4) create a dataframe with the Least squares regression result
    :param df: a standardized, binned dataframe
    :param parse_by: list of columns that will be used to parse the dataframe
    :return: a binned, tresholded NBS dataframe
    """
    import pandas as pd
    if light_parsing == True:
        parse_by.append('light_cond')
    if depth_parsing == True:
        parse_by.append('midDepthBin')
        df['midDepthBin'] = df['midDepthBin'].astype(str)
    df['date_bin'] = df['date_bin'].astype(str)
    NBSS_binned_all = pd.DataFrame()  # NBSS dataset
    lin_fit_data = pd.DataFrame()
    for st in list(set(df[parse_by[0]])):
        # print(st)
        df_st_subset = df[df[parse_by[0]] == st].reset_index(drop=True)
        for t in list(set(df_st_subset[parse_by[1]])):
            # print(t)
            df_st_t_subset = df_st_subset[df_st_subset[parse_by[1]] == t].reset_index(drop=True)
            if light_parsing == True:
                for l in list(set(df_st_t_subset[parse_by[2]])):
                    df_st_t_l_subset = df_st_subset[df_st_subset[parse_by[2]] == l].reset_index(drop=True)
                    if depth_parsing == True:
                        for d in list(set(df_st_t_l_subset[parse_by[3]])):
                            df_st_t_l_d_subset = df_st_t_l_subset[df_st_t_l_subset[parse_by[3]] == d].reset_index(drop=True)
                            df_binned, df_bins = size_binning_func(df_st_t_l_d_subset)  # create size bins
                            NBS_biovol_df = NB_SS_func(df_binned, df_bins, light_parsing = True, depth_parsing = True)  # calculate NBSS
                            if instrument == 'IFCB':
                                NBS_biovol_df= NBS_biovol_df[NBS_biovol_df['logSize'] !=17.778833578677386].reset_index(drop=True) # this and previous line are temporary, to remove faulty dataset from san francisco without having to go over all the pipeline
                            lin_fit = linear_fit_func(NBS_biovol_df, light_parsing = True, depth_parsing = True)
                            NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])
                            lin_fit_data = pd.concat([lin_fit_data, lin_fit])
                    else:
                        # print ('getting NBS for station ' + st + 'and date' + t)
                        df_binned, df_bins = size_binning_func(df_st_t_l_subset)  # create size bins
                        NBS_biovol_df = NB_SS_func(df_binned, df_bins, light_parsing = True)  # calculate NBSS
                        if instrument == 'IFCB':
                            NBS_biovol_df = NBS_biovol_df[NBS_biovol_df['logSize'] != 17.778833578677386].reset_index(drop=True)
                        lin_fit = linear_fit_func(NBS_biovol_df, light_parsing = True)
                        NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])
                        lin_fit_data = pd.concat([lin_fit_data, lin_fit])
            else:
                df_binned, df_bins = size_binning_func(df_st_t_subset)  # create size bins
                NBS_biovol_df = NB_SS_func(df_binned, df_bins)  # calculate NBSS
                if instrument == 'IFCB':
                    NBS_biovol_df = NBS_biovol_df[NBS_biovol_df['logSize'] <= 17.5].reset_index(drop=True)
                lin_fit = linear_fit_func(NBS_biovol_df)
                NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])
                lin_fit_data = pd.concat([lin_fit_data, lin_fit])

    #NBSS_binned_all['Station_location'], NBSS_binned_all['midLatBin'], NBSS_binned_all['midLonBin'] = gridding_func(bin_loc, NBSS_binned_all['midLatBin'], NBSS_binned_all['midLonBin'])
    NBSS_raw = NBSS_binned_all.filter(['date_bin', 'midLatBin', 'midLonBin', 'light_cond','size_class_mid', 'NB', 'Min_obs_depth', 'Max_obs_depth'], axis=1)
    NBSS_1a_raw = NBSS_raw.rename(columns={'date_bin': 'date_year_month_week', 'midLatBin': 'latitude', 'midLonBin': 'longitude', 'light_cond': 'light_condition', 'size_class_mid': 'biovolume_size_class',
                                      'NB': 'normalized_biovolume', 'Min_obs_depth': 'min_depth', 'Max_obs_depth': 'max_depth'})
    if light_parsing==True:
        NBSS_1a = NBSS_stats_func(NBSS_raw, light_parsing = True, bin_loc = bin_loc, group_by = group_by)
        lin_fit_1b = stats_linfit_func(lin_fit_data, light_parsing = True, bin_loc = bin_loc, group_by = group_by)
    else:
        NBSS_1a = NBSS_stats_func(NBSS_raw,  bin_loc=bin_loc, group_by=group_by)
        lin_fit_1b = stats_linfit_func(lin_fit_data, bin_loc = bin_loc, group_by = group_by)


    return NBSS_1a_raw, NBSS_1a,  lin_fit_1b


def reshape_date_func(df, group_by ='yyyymm'):
    import numpy as np
    date_df = df['date_bin'].str.split("_", expand=True)
    if group_by == 'yyyy':
        df['year'] = date_df[0]
    elif group_by == 'yyyymm':
        df['year'] = date_df[0]
        df['month'] = date_df[1]
    return df



def NBSS_stats_func(df, light_parsing = False, bin_loc = 1, group_by = 'yyyymm'):
    """

    """
    from funcs_gridding import gridding_func
    import numpy as np
    df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(bin_loc, df['midLatBin'], df['midLonBin'])
    date_df = df['date_bin'].str.split("_", expand=True)
    if group_by == 'yyyy':
        df['year'] = date_df[0]
        df['date_bin'] = df['year']
    elif group_by == 'yyyymm':
        df['year'] = date_df[0]
        df['month'] = date_df[1]
        df['date_bin'] = df['year'] + '_' + df['month']

    if light_parsing == True:
        NBSS_avg = df.groupby(['date_bin', 'Station_location', 'light_cond', 'size_class_mid']).agg(
            {'year': 'first', 'month': 'first', 'midLatBin':'first', 'midLonBin': 'first', 'NB':['count', 'mean', 'std'], 'Min_obs_depth':'first', 'Max_obs_depth':'first'}).reset_index()
    else:
        NBSS_avg = df.groupby(['date_bin', 'Station_location', 'size_class_mid']).agg(
            {'year': 'first', 'month': 'first', 'midLatBin':'first', 'midLonBin': 'first', 'NB':['count', 'mean', 'std'], 'Min_obs_depth':'first', 'Max_obs_depth':'first'}).reset_index()



    # generate a cleaner dataset (without grouping variables) to finally present 1a
    NBSS_avg.columns = NBSS_avg.columns.map('_'.join).str.removesuffix("first")
    NBSS_avg.columns = NBSS_avg.columns.str.removesuffix("_")
    NBSS_avg= NBSS_avg.filter(['year', 'month', 'midLatBin', 'midLonBin', 'light_cond', 'Min_obs_depth', 'Max_obs_depth', 'size_class_mid','NB_mean', 'NB_std', 'NB_count'], axis=1)
    NBSS_avg= NBSS_avg.rename(columns={'midLatBin':'latitude', 'midLonBin': 'longitude', 'light_cond':'light_condition',
                                       'Min_obs_depth':'min_depth', 'Max_obs_depth':'max_depth', 'size_class_mid': 'biovolume_size_class',
                                       'NB_count':'N', 'NB_mean':'normalized_biovolume_mean', 'NB_std': 'normalized_biovolume_std' })
    NBSS_avg= NBSS_avg[NBSS_avg.N !=0].reset_index(drop = True)
    return NBSS_avg

def stats_linfit_func(df, light_parsing = False, bin_loc = 1, group_by = 'yyyymm'):
    """
    Objective: create summary statistics for the linear fit dataframes
    param df: a dataframe containing the results of the linear fits for each station, time and depth for a given imaging instrument
    """
    from funcs_gridding import gridding_func
    import numpy as np
    #bin_loc =input('Group data by location? \n Enter Y/N ') commented out query since new method of parsing data added 3/2/2023 makes this a nuisance, will be asked in step4
    #if bin_loc == 'Y':
        #st_increment = float(input('select the size of the spatial grid to group the data i.e for a 1x1 degree bin type 1 \n' ))
    df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(bin_loc, df['latitude'], df['longitude'])
    #bin_date =input('Group data by date? \n Enter Y/N ')
    #if bin_date == 'Y':
        #group_by= input('input how will the data be grouped by date \n  yyyymm for month and year \n yyyy for year \n ')
    date_df= df['Date'].str.split("_", expand=True)
    if group_by == 'yyyy':
        df['year'] = date_df[0]
        df['date_bin'] = df['year']
    elif group_by == 'yyyymm':
        df['year'] = date_df[0]
        df['month'] = date_df[1]
        df['date_bin'] = df['year'] + '_' + df['month']
    #if group_by == 'yyyy':
        #df['date_bin'] = date_df[0]
    #elif group_by == 'yyyymm':
        #df['year_week'] = date_df[0] + '_' + date_df[2]
        #df['month'] = date_df[1]
        #week_dict = {key:None for key in df['year_week'].unique()}
        #for i in df['year_week'].unique():
            #week_dict[i] =list(df['month'].where(df['year_week'] == i).dropna().unique()) #.astype(int)
            #if len(week_dict[i]) == 1:
                #week_dict[i] = week_dict[i][0]
            #else:
                #week_dict[i] = list(map(int, week_dict[i]))
                #week_dict[i] = str(int(np.round(np.mean(week_dict[i])))).zfill(2)
        #df['year'] = date_df[0]
        #df['month'] = df['year_week'].map(week_dict)
        #df['date_bin'] = date_df[0] + '_' + df['month']
        #df = df.drop(columns=['Date', 'year_week'])

    if light_parsing ==True:
        lin_fit_stats = df.groupby(['Station_location', 'date_bin', 'light_condition']).agg(
            {'year': 'first', 'month': 'first', 'midLatBin': 'first', 'midLonBin':'first', 'min_depth': 'first', 'max_depth':'first',
             'slope': ['count', 'mean', 'std'], 'intercept': ['count', 'mean', 'std'], 'r2': ['count', 'mean', 'std']}).reset_index()
    else:
        lin_fit_stats = df.groupby(['Station_location', 'date_bin']).agg(
            {'year': 'first', 'month': 'first', 'midLatBin': 'first', 'midLonBin': 'first', 'min_depth': 'first','max_depth': 'first',
             'slope': ['count', 'mean', 'std'], 'intercept': ['count', 'mean', 'std'], 'r2': ['count', 'mean', 'std']}).reset_index()

    lin_fit_stats.columns = lin_fit_stats.columns.map('_'.join).str.removesuffix("first")
    lin_fit_stats.columns = lin_fit_stats.columns.str.removesuffix("_")
    lin_fit_stats= lin_fit_stats.rename(columns={'slope_count': 'N', 'midLatBin': 'latitude', 'midLonBin': 'longitude'})
    lin_fit_stats= lin_fit_stats[lin_fit_stats.columns.drop(list(lin_fit_stats.filter(regex='count')))]
    del lin_fit_stats['Station_location']
    del lin_fit_stats['date_bin']

    return lin_fit_stats




# Plotting Values and Regression Line
#max_x = np.max(X) + 100
#min_x = np.min(X) - 100
# Calculating line values x and y
#x = np.linspace(min_x, max_x, 1000)
#y = c + m * x
# Ploting Line
#plt.plot(x, y, color='#58b970', label='Regression Line')
# Ploting Scatter Points
#plt.scatter(X, Y, c='#ef5423', label='Scatter Plot')
#plt.xlabel('logSize')
#plt.ylabel('logNBSS')
#plt.legend()
#plt.show()