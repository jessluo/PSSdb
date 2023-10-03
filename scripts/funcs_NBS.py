# functions to calculate normalized biomass size spectra, based on ecopart size bins and a logarithmic scale
# data needs to be gridded and temporally and spatially binned. proj_nbs_func summarizes all the steps
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from wcmatch.pathlib import Path # Handling of path object
import pandas as pd
import yaml
import math as m
import shutil
import tarfile
try:
    from funcs_read import *
    from funcs_gridding import *
except:
    from scripts.funcs_read import *
    from scripts.funcs_gridding import *

import os
from tqdm import tqdm
from operator import attrgetter
import geopandas as gpd
from scipy.stats import poisson,norm,lognorm # Estimate uncertainties assuming count detection follows a Poisson distribution and size a normal distribution
from itertools import compress

# line 23-32: adds the oceans shape file needed to assing ocean basin to the dataset. See https://www.marineregions.org/. https://doi.org/
if len(list(Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent.rglob('goas_v01.shp'))) == 0:
    path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
    # open the metadata of the standardized files
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    with tarfile.open(str(Path(cfg['git_dir']).expanduser() / cfg['ancillary_subdir'] /  cfg['oceans_subdir_tar']), 'r') as tar:
        tar.extractall(str(Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent))

oceans = gpd.read_file(list(Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent.rglob('goas_v01.shp'))[0])
  #6.1 by size bins. Inputs: a column of the dataframe with biovolume, and the number of log bins to use.
    #returns two lists : the sizeClasses bins (categorical) and range_size_bins (float)
def size_binning_func(df_subset, biovol_estimate):
    """
    Objective: bin the dataframes (parsed by stations and/or dates and/or depths)
    :param biovolume: column of a dataframe that contains the biovolume (in cubic micrometers)
    :return: a binned dataframe by sizes, containing the range of the size bin for each
    """

    # open dataframe with the size bins
    path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
    # open the metadata of the standardized files
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    path_to_bins = str(Path(cfg['size_bins_path']).expanduser())
    bins_df = pd.read_csv(path_to_bins)
    # create a categorical variable, which are bins defined by the numbers on the sizeClasses list
    # Assign bin to each data point based on biovolume, append bin range and bin class to the dataframe
    df_subset.loc[:, 'sizeClasses']= pd.cut(x=df_subset[biovol_estimate], bins=bins_df['biovol_um3'], include_lowest=True)# size classes defined by biovolume
    df_subset = df_subset.dropna(subset=['sizeClasses']).reset_index(drop=True) # hard pass, debug this
    range_size_bin = []
    size_class_mid = []
    for i in range(0, len(df_subset['sizeClasses'])):
        range_size_bin.append(df_subset.loc[i, 'sizeClasses'].length) # consider replacing this loop with df_binned['sizeClasses'].map(attrgetter('length'))
        size_class_mid.append(df_subset.loc[i, 'sizeClasses'].mid)
    df_subset.loc[:, 'range_size_bin'] = range_size_bin
    df_subset.loc[:, 'size_class_mid'] = size_class_mid
    df_subset.loc[:, 'size_range_ECD']=((np.array(range_size_bin) * 6) / m.pi) ** (1. / 3.)
    df_subset.loc[:, 'ECD_mean'] = ((np.array(size_class_mid) * 6) / m.pi) ** (1. / 3.)

    bins = df_subset['sizeClasses'].unique()
    bin_width = []# list that will contain the width of each size_bin
    bin_mid = []
    for i in range(0, len(bins.categories)):
        bin_width.append(bins.categories[i].length)
        bin_mid.append(bins.categories[i].mid)
    bin_width_ECD=((np.array(bin_width) * 6) / m.pi) ** (1. / 3.)
    bin_mean_ECD = ((np.array(bin_mid) * 6) / m.pi) ** (1. / 3.)
    d = {'sizeClasses':bins.categories, 'range_size_bin': bin_width,'size_class_mid': bin_mid, 'size_range_ECD': bin_width_ECD, 'ECD_mean': bin_mean_ECD}
    df_bins = pd.DataFrame(d).astype(dict(zip(['range_size_bin', 'size_class_mid', 'size_range_ECD', 'ECD_mean'], [str]*4)))

    # obtain the size of each size class bin and append it to the stats dataframe
    # unfortunately the Interval type produced by pd.cut is hard to work with. So they need to be modified

    return df_subset, df_bins


## Different ways of calculating biovolume for sensitivity analyses
def biovolume_metric(df_binned):
    ###
    ###
    df_binned= df_binned.drop(columns=['Biovolume']).reset_index()

# 7)calculate x and y for NBSS, this includes adding total biovolume per size bin for each station and depth bin,
# inputs: a dataframe and the volume sampled of each (in cubic meters). other inputs are the column names
# from the original dataset that want to be retained (i.e;  proj ID, station ID, depth, size classes,
# range of size classes, biovolume, lat & lon ) plus the variables that will be used to group the data by
# stations, depths and size classes
# using 5 ml for IFCB
def NB_SS_func(NBS_biovol_df, df_bins, biovol_estimate = 'Biovolume_area',sensitivity = False, light_parsing = False, depth_parsing = False,thresholding=True): #, niter=0
    """
    """
    import numpy as np
    import pandas as pd
    if sensitivity == True:
        NBS_biovol_df = NBS_biovol_df.dropna(subset=['Minor_axis']).reset_index(drop =True)
    if depth_parsing == True:
        # create a dataframe with summary statistics for each station, depth bin and size class
        # these column names should all be the same, since the input is a dataframe from the 'binning' and 'biovol' functions
        # group data by bins
        NBS_biovol_df[biovol_estimate] = NBS_biovol_df[biovol_estimate] * NBS_biovol_df['ROI_number']
        df_vol = NBS_biovol_df.groupby(['date_bin', 'Station_location']).apply(lambda x: pd.Series({'cumulative_vol': x[['Sample', 'Min_obs_depth', 'Max_obs_depth','Volume_imaged']].drop_duplicates().Volume_imaged.sum()})).reset_index()
        NBS_biovol_df['cumulative_vol'] = pd.merge(NBS_biovol_df, df_vol, how='left', on=['date_bin', 'Station_location'])['cumulative_vol']  # df_vol.loc[0, 'cumulative_volume']

        NBS_biovol_df = NBS_biovol_df.groupby(['date_bin', 'Station_location', 'sizeClasses', 'light_cond', 'midDepthBin']).apply(lambda x: pd.Series({'Sample': x.Sample.unique()[0],
                                                                                                                                    'size_class_mid': x.size_class_mid.unique()[0],
                                                                                                                                    'midLatBin': x.midLatBin.unique()[0],
                                                                                                                                    'midLonBin': x.midLonBin.unique()[0],
                                                                                                                                    'Min_obs_depth': x.Min_obs_depth.unique()[0],
                                                                                                                                    'Max_obs_depth': x.Max_obs_depth.unique()[0],
                                                                                                                                    'Biovolume_mean': x[biovol_estimate].sum() / x.ROI_number.sum(),
                                                                                                                                    'NB': (x[biovol_estimate].sum() / x.cumulative_vol / x.range_size_bin).unique()[0]})).reset_index()  # , 'ROI_number':x['ROI_number'].sum()}))

    else:
        if light_parsing == True:
            NBS_biovol_df[biovol_estimate] = NBS_biovol_df[biovol_estimate] * NBS_biovol_df['ROI_number']
            df_vol = NBS_biovol_df.groupby(['date_bin', 'Station_location']).apply(lambda x: pd.Series({'cumulative_vol': x[
                ['Sample', 'Min_obs_depth', 'Max_obs_depth',
                 'Volume_imaged']].drop_duplicates().Volume_imaged.sum()})).reset_index()
            NBS_biovol_df['cumulative_vol'] = pd.merge(NBS_biovol_df, df_vol, how='left', on=['date_bin', 'Station_location'])[
                'cumulative_vol']  # df_vol.loc[0, 'cumulative_volume']

            NBS_biovol_df = NBS_biovol_df.groupby(['date_bin', 'Station_location', 'sizeClasses', 'light_cond']).apply(lambda x: pd.Series({'Sample': x.Sample.unique()[0],
                                                                                                                                        'size_class_mid': x.size_class_mid.unique()[0],
                                                                                                                                        'midLatBin': x.midLatBin.unique()[0],
                                                                                                                                        'midLonBin': x.midLonBin.unique()[0],
                                                                                                                                        'Min_obs_depth': x.Min_obs_depth.unique()[0],
                                                                                                                                        'Max_obs_depth': x.Max_obs_depth.unique()[0],
                                                                                                                                        'Biovolume_mean': x[biovol_estimate].sum() / x.ROI_number.sum(),
                                                                                                                                        'NB': (x[biovol_estimate].sum() / x.cumulative_vol / x.range_size_bin).unique()[0]})).reset_index()  # , 'ROI_number':x['ROI_number'].sum()}))

        else:
            NBS_biovol_df[biovol_estimate] = NBS_biovol_df[biovol_estimate] * NBS_biovol_df['ROI_number']
            grouping = ['Sample', 'Sampling_lower_size', 'Sampling_upper_size', 'date_bin', 'Station_location', 'midLatBin', 'midLonBin','Depth_min', 'Depth_max', 'Min_obs_depth', 'Max_obs_depth'] if NBS_biovol_df.Instrument.unique()[0] == 'Scanner' else ['date_bin', 'Station_location','midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth'] # , 'size_class_mid', 'range_size_bin', 'ECD_mean', 'size_range_ECD'
            df_vol = NBS_biovol_df.groupby(grouping, dropna=False).apply(lambda x: pd.Series({'cumulative_vol': x[['Sample', 'Depth_min', 'Depth_max','Volume_imaged']].drop_duplicates().Volume_imaged.sum(),
                                                                                          'Depth_range_min':x.Depth_min.astype(float).min(),
                                                                                          'Depth_range_max':x.Depth_max.astype(float).max()})).reset_index()
            NBS_biovol_df= pd.merge(NBS_biovol_df, df_vol, how='left', on=grouping )  # df_vol.loc[0, 'cumulative_volume']
            #computing NBSS for individual size fractions in each samples. Required for zooscan projects, since a single net tow is fractionated by sieves
            #for index_group in list(NBS_biovol_df.astype(dict(zip(grouping+['Depth_range_min', 'Depth_range_max'],[str]*(2+len(grouping))))).groupby(grouping+['sizeClasses','Depth_range_min', 'Depth_range_max']).groups.values()):
            #    x=NBS_biovol_df.loc[index_group]
            #    x.Station_location.unique()[0] if 'Sample' not in grouping else x.Sample.unique()[0]
            #    x[biovol_estimate].sum() / x.ROI_number.sum()
            #    x[biovol_estimate].sum()
            #    x.ROI_number.sum()
            #    x.cumulative_vol.unique()[0]
            #    (x.ROI_number.sum() / x.cumulative_vol.unique()[0] / (((x.range_size_bin.astype(float).unique()[0] * 6) / m.pi) ** (1. / 3.)))
            #    (x[biovol_estimate].sum() / x.cumulative_vol.unique()[0] / x.range_size_bin.astype(float).unique()[0])
            NBS_biovol_df = NBS_biovol_df.astype(dict(zip(grouping+['sizeClasses','Depth_range_min', 'Depth_range_max'],[str]*(3+len(grouping))))).groupby(grouping+['sizeClasses','Depth_range_min', 'Depth_range_max'],observed=True).apply(lambda x: pd.Series({'Sample_id': x.Station_location.unique()[0] if 'Sample' not in grouping else x.Sample.unique()[0],  # 'fake' sample identifier, so that grouping by Sample in the next step is the same as grouing it by Station_location
                                                                                                   'pixel_size':np.nanmean(x.Pixel),
                                                                                                   'Biovolume_mean': x[biovol_estimate].sum() / x.ROI_number.sum(),
                                                                                                   'Biovolume_sum': x[biovol_estimate].sum(),
                                                                                                   'ROI_number_sum': x.ROI_number.sum(),
                                                                                                   'Total_volume': x.cumulative_vol.unique()[0],
                                                                                                   'PSD': (x.ROI_number.sum() /x.cumulative_vol.unique()[0] / (((x.range_size_bin.astype(float).unique()[0]*6)/m.pi)**(1./3.))),
                                                                                                   #'NB_std': np.std(list(map(lambda count: np.random.normal( loc=np.random.choice(np.repeat(x[biovol_estimate],x.ROI_number),size=count,replace=True),
                                                                                                                                                             #scale=(1 / 3) * np.pi * (np.sqrt((x.Pixel.unique()[0]) * 1e-03 * x.size_class_mid.astype(float).unique()[0]) ** 3) / ( (1e-03 * x.Pixel.unique()[0]) ** 3),
                                                                                                                                                             #size=None).sum() / x.cumulative_volume.unique()[0] / (x.range_size_bin.astype(float).unique())[0],
                                                                                                                                                             #np.random.poisson(lam=sum(x.ROI_number), size=niter)))),

                                                                                                   'NB': (x[biovol_estimate].sum() /x.cumulative_vol.unique()[0] / x.range_size_bin.astype(float).unique()[0])})).reset_index()#, 'ROI_number':x['ROI_number'].sum()}))

            #df_bins['sizeClasses']=df_bins.sizeClasses.astype(str)
            #NBS_biovol_df = pd.merge(NBS_biovol_df, df_bins[['sizeClasses', 'size_class_mid', 'range_size_bin', 'ECD_mean', 'size_range_ECD']].drop_duplicates(),how='left', on='sizeClasses')
            NBS_biovol_df.sizeClasses = pd.Categorical(NBS_biovol_df.sizeClasses, df_bins.sizeClasses.astype(str),ordered=True) # this generates a categorical index for the column, this index can have a different lenght that the actual number of unique values, especially here since the categories come from df_bins
            NBS_biovol_df = pd.merge(NBS_biovol_df, NBS_biovol_df.drop_duplicates(subset=grouping, ignore_index=True)[grouping].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left', on=grouping)
            multiindex = pd.MultiIndex.from_product([list(NBS_biovol_df.astype({column: 'category' for column in ['Group_index', 'sizeClasses']})[column].cat.categories) for column in ['Group_index', 'sizeClasses']],names=['Group_index', 'sizeClasses'])
            NBS_biovol_df =  pd.merge(NBS_biovol_df.drop_duplicates(['Group_index'])[grouping + ['Sample_id','Depth_range_min', 'Depth_range_max','Group_index','Total_volume']],NBS_biovol_df.set_index(['Group_index', 'sizeClasses']).reindex(multiindex,fill_value=pd.NA).reset_index().drop(columns=grouping+['Sample_id','Depth_range_min', 'Depth_range_max','Total_volume']), how='right', on=['Group_index']).sort_values(['Group_index']).reset_index(drop=True)


            NBS_biovol_df= pd.merge(NBS_biovol_df, df_bins[['sizeClasses','size_class_mid', 'range_size_bin','ECD_mean', 'size_range_ECD']].astype({'sizeClasses':str}).drop_duplicates(), how='left', on='sizeClasses')
            NBS_biovol_df =NBS_biovol_df.astype({'ECD_mean':float}).sort_values(['Group_index','ECD_mean']).reset_index(drop=True).astype({'ECD_mean':str})

            # summing the NB for each size class that come from separate size fractions(because of sieving) to obtain the overall NB for a size class in a net
            NBS_biovol_df= NBS_biovol_df.astype(dict(zip(['midLatBin', 'midLonBin','size_class_mid', 'range_size_bin','ECD_mean', 'size_range_ECD'],[str]*6))).groupby(['Sample_id',  'date_bin', 'Station_location','midLatBin', 'midLonBin','sizeClasses','size_class_mid', 'range_size_bin','ECD_mean', 'size_range_ECD']).apply(lambda x: pd.Series({ #.astype(dict(zip(['Depth_range_min','Depth_range_max'],[str]*2)))
                                                                                                                     'Depth_range_min': x.Depth_range_min.unique()[0],
                                                                                                                     'Depth_range_max': x.Depth_range_max.unique()[0],
                                                                                                                     'Biovolume_mean': x.Biovolume_sum.sum() / x.ROI_number_sum.sum(),
                                                                                                                     'Biovolume_sum': x.Biovolume_sum.sum(),
                                                                                                                     'Pixel_mean':np.nanmean(x.pixel_size.astype(float) * 1e-03 * x.ECD_mean.astype(float)),
                                                                                                                     'ROI_number_sum': x.ROI_number_sum.sum(),
                                                                                                                     'ROI_abundance': np.nansum(x.ROI_number_sum/ x.Total_volume),
                                                                                                                     #'NB_std': np.sqrt(sum(x.NB_std ** 2)),
                                                                                                                     'PSD': x.PSD.sum(),
                                                                                                                     'NB': x.NB.sum()})).reset_index()#, 'ROI_nu
            #average NB for each date and station bin, from NB that come from different net tows. NOTE: here, things need to be done differently depending on


            nbss_depths = NBS_biovol_df.groupby(['date_bin','midLatBin', 'midLonBin','Station_location']).apply(lambda x: pd.Series({ 'Min_obs_depth':np.nanmin(x.Depth_range_min.astype(float)),
                                                                                                                                            'Max_obs_depth': np.nanmax(x.Depth_range_max.astype(float))})).reset_index()


            NBS_biovol_df = pd.merge(NBS_biovol_df, nbss_depths, how='left', on=['date_bin','midLatBin', 'midLonBin','Station_location'])

            # depth integration OR weighted sum
            NBS_biovol_df = NBS_biovol_df.astype(dict(zip(['Min_obs_depth','Max_obs_depth'],[str]*2))).groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin','Min_obs_depth','Max_obs_depth', 'sizeClasses', 'size_class_mid', 'range_size_bin','ECD_mean', 'size_range_ECD']).apply(lambda x: pd.Series({
                                                                                                                        'Biovolume_mean': x.Biovolume_mean.mean(),# in cubic micrometers
                                                                                                                        'size_class_pixel': x.Pixel_mean.mean(),# Number of pixels
                                                                                                                        'ROI_number_sum': np.nansum(x.ROI_number_sum),
                                                                                                                        'ROI_abundance_mean': np.nansum(x.ROI_abundance * np.where((x.Depth_range_max.astype(float) - x.Depth_range_min.astype(float)) > 1,(x.Depth_range_max.astype(float) - x.Depth_range_min.astype(float)), 1)) / np.where( (x.Max_obs_depth.astype(float).unique()[0] - x.Min_obs_depth.astype(float).unique()[0]) > 1,(x.Max_obs_depth.astype(float).unique()[0] - x.Min_obs_depth.astype(float).unique()[0]), 1), # also need to integrate

                                                                                                                        #'NB_std': np.sqrt(sum((((x.NB_std) ** 2) * ((np.where((x.Max_obs_depth.astype(float) - x.Min_obs_depth.astype(float)).values > 1,(x.Max_obs_depth.astype(float) - x.Min_obs_depth.astype(float)).values, 1) / np.where((x.Depth_range_max.astype(float).unique()[0] - x.Depth_range_min.astype(float).unique()[0]) > 1,(x.Depth_range_max.astype(float).unique()[0] - x.Depth_range_min.astype(float).unique()[0]),1)) ** 2)))),
                                                                                                                        'NB': np.nansum(x.NB * np.where((x.Depth_range_max.astype(float) - x.Depth_range_min.astype(float)) > 1,(x.Depth_range_max.astype(float) - x.Depth_range_min.astype(float)), 1)) / np.where( (x.Max_obs_depth.astype(float).unique()[0] - x.Min_obs_depth.astype(float).unique()[0]) > 1,(x.Max_obs_depth.astype(float).unique()[0] - x.Min_obs_depth.astype(float).unique()[0]), 1),
                                                                                                                        'PSD': np.nansum(x.PSD * np.where((x.Depth_range_max.astype(float) - x.Depth_range_min.astype(float)) > 1,(x.Depth_range_max.astype(float) - x.Depth_range_min.astype(float)), 1)) / np.where( (x.Max_obs_depth.astype(float).unique()[0] - x.Min_obs_depth.astype(float).unique()[0]) > 1,(x.Max_obs_depth.astype(float).unique()[0] - x.Min_obs_depth.astype(float).unique()[0]), 1)})).reset_index()

            NBS_biovol_df= NBS_biovol_df.astype({'ECD_mean':float}).sort_values(['date_bin','midLatBin', 'midLonBin','Station_location','ECD_mean']).reset_index(drop=True)
            NBS_biovol_df=NBS_biovol_df.assign(count_uncertainty=poisson.pmf(k=NBS_biovol_df['ROI_number_sum'], mu=NBS_biovol_df['ROI_number_sum']),
                                 size_uncertainty=NBS_biovol_df.astype({'size_class_pixel':float}).groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin','Min_obs_depth','Max_obs_depth']).size_class_pixel.transform(lambda x: norm.pdf((1/3)*np.pi*x**3,loc=(1/3)*np.pi*x.min()**3,scale=(1/3)*1*np.pi)).values.tolist())

    # create three more columns with the parameters of particle size distribution and normalized size spectra,:
    NBS_biovol_df['logNB'] = np.log10(NBS_biovol_df['NB'])
    NBS_biovol_df['logSize'] = np.log10(NBS_biovol_df['size_class_mid'].astype(float))

    NBS_biovol_df['logPSD'] = np.log10(NBS_biovol_df['PSD'])
    NBS_biovol_df['logECD'] = np.log10(NBS_biovol_df['ECD_mean'].astype(float))

    # now calculate the log biovolume for the df_bins dataframe. NOTE updated on 2/23/2023. middle point of the size class
    df_bins['logSize'] = np.log10(df_bins['size_class_mid'].astype(float))
    # merge the two dataframes
    df_bins=df_bins.astype(dict(zip(['range_size_bin', 'size_class_mid', 'size_range_ECD', 'ECD_mean','logSize'], [float] * 5)))
    NBS_biovol_df = NBS_biovol_df.astype(dict(zip([ 'range_size_bin', 'size_class_mid', 'size_range_ECD', 'ECD_mean', 'logSize'], [float] * 5)))

    #NBS_biovol_df2= pd.merge(df_bins, NBS_biovol_df, how='left', on=['sizeClasses','range_size_bin', 'size_class_mid', 'size_range_ECD', 'ECD_mean','logSize'])
    # now fill the columns of date, station, lat/lon, project ID and volume
    #for i in ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth']:
        #NBS_biovol_df[i] = NBS_biovol_df[i].unique()[1]
    # let's do the thresholding here:
    if thresholding==True:
        NBS_biovol_df2 = NBS_biovol_df.groupby(['date_bin', 'Station_location']).apply(lambda x: threshold_func(x))
    else:
        NBS_biovol_df= NBS_biovol_df

    NBS_biovol_df = NBS_biovol_df.astype(dict(zip(['midLatBin', 'midLonBin','size_class_mid', 'range_size_bin','ECD_mean', 'size_range_ECD'], [float] * 6)))


    return NBS_biovol_df

## assing ocean label using gdp
def ocean_label_func(df, Lon, Lat):
    df.columns = [c.replace(' ', '_') for c in df.columns]
    df = pd.merge(df.astype(dict(zip([Lon, Lat], [str] * 2))),df.astype(
                                    dict(zip([Lon, Lat], [str] * 2))).drop_duplicates(subset=[Lon, Lat], ignore_index=True)[
                                    [Lon, Lat]].reset_index().rename({'index': 'Group_index'},axis='columns'),
                                how='left', on=[Lon, Lat]).astype(dict(zip([Lon, Lat], [float] * 2)))
    gdf = gpd.GeoDataFrame(df[['Group_index', Lon, Lat]].drop_duplicates().dropna(),
                           geometry=gpd.points_from_xy(df[['Group_index', Lon, Lat]].drop_duplicates().dropna()[Lon],
                               df[['Group_index', Lon, Lat]].drop_duplicates().dropna()[Lat]))
    df['ocean'] = pd.merge(df, gpd.tools.sjoin_nearest(gdf, oceans, how='left')[['Group_index', 'name']], how='left',on='Group_index')['name'].astype(str)
    ocean = df.pop('ocean')
    df.insert(4, ocean.name, ocean)
    df = df.drop(columns=['Group_index'])
    return df

# 7)  create a function to remove data based on JO's comments and Haentjens THIS IS WORK IN PROGRESS:
# low threshold will be filtered because they are smaller than the next standardized size spectrum
# high threshold: obtain the highest size class, then iterate from highest to lowest and remove the large size bins
# that have less than an (instrument specific) threshold
# imputs: the data frame, and the column that contains the NBSS
def threshold_func(binned_data, empty_bins = 3,threshold_count=0.2,threshold_size=0.2):
    """
    Objective: remove size bins based on Fabien's approach: (lower limit: max NBSS, upper limit: 3 NBSS with Nans)
    :param df: a binned dataframe with NBSS already calculated
    :return: a dataset without the bins that contain inaccurate NBSS
    """
    binned_data= binned_data.reset_index(drop=True)
    #upper threshold based on max NB
    binned_data_filt = binned_data.iloc[binned_data.NB.idxmax():len(binned_data), :].reset_index(drop=True)#a ask if its better to just remove all data

    #binned_data_filt = binned_data_filt.loc[0: max(binned_data_filt['NB'].isnull()[binned_data_filt['NB'].isnull() == False].index.to_numpy())] # cut dataframe to include only size classes up to the largest observed particle

    empty_SC = binned_data_filt['NB'].isnull()  # get row indices with empty size bins
    sum_empty = empty_SC.ne(empty_SC.shift()).cumsum() # add a value each time there is a shift in true/false
    cum_empty = sum_empty.map(sum_empty.value_counts()).where(empty_SC) # add occurrences of a value if the row is nan (based on empty SC) and find if the empty_bins (consecutive empty bins) is in the array
    if empty_bins in cum_empty.to_numpy():
        binned_data_filt = binned_data_filt.loc[0:(min(cum_empty[cum_empty.isnull() == False].index.to_numpy())-1)]
    else:
        binned_data_filt = binned_data_filt

    #lower threshold based on three consecutive size bins with nans
    #for n, i in enumerate(np.isnan(binned_data_filt['NB'])):
        #if (i == True) and (np.isnan(binned_data_filt.loc[n+1, 'NB']) == True) and (np.isnan(binned_data_filt.loc[n+2, 'NB']) == True): #edited to test thresholding
            #binned_data_filt = binned_data_filt.loc[0:n-1]
            #break

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
    df1 = df1[df1['NB'].notna()].reset_index()

    # extract x and y
    def regression(df1, X_var, Y_var):
        '''
        :param Y_var: string that identifies the column header of the df1 dataframe to use as Y variable
        '''
        X = df1[X_var].values#.reshape(-1, 1)
        Y = df1[Y_var].values#.reshape(-1, 1)
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
        return m,b,rmse,R2
    m_NB, b_NB, rmse_NB, R2_NB = regression(df1, 'logSize','logNB')
    m_PSD, b_PSD, rmse_PSD, R2_PSD = regression(df1, 'logECD', 'logPSD')
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

    lin_fit.loc[0, 'NBSS_slope'] = m_NB
    lin_fit.loc[0, 'NBSS_intercept'] = b_NB
    lin_fit.loc[0, 'NBSS_rmse'] = rmse_NB
    lin_fit.loc[0, 'NBSS_r2'] = R2_NB

    lin_fit.loc[0, 'PSD_slope'] = m_PSD
    lin_fit.loc[0, 'PSD_intercept'] = b_PSD
    lin_fit.loc[0, 'PSD_rmse'] = rmse_PSD
    lin_fit.loc[0, 'PSD_r2'] = R2_PSD


    return (lin_fit)

def parse_NBS_linfit_func(df, biovol_estimate, sensitivity, light_parsing, depth_parsing, bin_loc = 1, group_by = 'yyyymm'):
    """
    Objective: parse out any df that come from files that include multiple stations, dates and depths.
    dataframe needs to be already binned. Steps are: 1) parse the data 2) bin ROIs by size  3) calculate the NBS
    (threshold included) and 4) create a dataframe with the Least squares regression result
    :param df: a standardized, binned dataframe
    :param parse_by: list of columns that will be used to parse the dataframe
    :return: a binned, tresholded NBS dataframe
    """
    import pandas as pd
    parse_by=['Station_location', 'date_bin']
    if light_parsing == True:
        parse_by.append('light_cond')
    if depth_parsing == True:
        parse_by.append('midDepthBin')
        df['midDepthBin'] = df['midDepthBin'].astype(str)
    df['date_bin'] = df['date_bin'].astype(str)
    NBSS_binned_all = pd.DataFrame()  # NBSS dataset
    lin_fit_data = pd.DataFrame()
    #df_binned_test, df_bins = df.groupby(['Station_location', 'date_bin']).apply(lambda x: size_binning_func(x)).reset_index()
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
                            df_binned, df_bins = size_binning_func(df_st_t_l_d_subset, biovol_estimate)  # create size bins
                            NBS_biovol_df = NB_SS_func(df_binned, df_bins, biovol_estimate = biovol_estimate,sensitivity = sensitivity, light_parsing = True, depth_parsing = True)  # calculate NBSS
                            #if instrument == 'IFCB':
                                #NBS_biovol_df= NBS_biovol_df[NBS_biovol_df['logSize'] !=17.778833578677386].reset_index(drop=True) # this and previous line are temporary, to remove faulty dataset from san francisco without having to go over all the pipeline
                            lin_fit = linear_fit_func(NBS_biovol_df,  light_parsing = True, depth_parsing = True)
                            NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])
                            lin_fit_data = pd.concat([lin_fit_data, lin_fit])
                    else:
                        # print ('getting NBS for station ' + st + 'and date' + t)
                        df_binned, df_bins = size_binning_func(df_st_t_l_subset, biovol_estimate)  # create size bins
                        NBS_biovol_df = NB_SS_func(df_binned, df_bins, biovol_estimate = biovol_estimate,sensitivity = sensitivity, light_parsing = True)  # calculate NBSS
                        #if instrument == 'IFCB':
                            #NBS_biovol_df = NBS_biovol_df[NBS_biovol_df['logSize'] != 17.778833578677386].reset_index(drop=True)
                        lin_fit = linear_fit_func(NBS_biovol_df, light_parsing = True)
                        NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])
                        lin_fit_data = pd.concat([lin_fit_data, lin_fit])
            else:
                df_binned, df_bins = size_binning_func(df_st_t_subset, biovol_estimate)  # create size bins
                NBS_biovol_df = NB_SS_func(df_binned, df_bins, biovol_estimate = biovol_estimate,sensitivity = sensitivity)  # calculate NBSS
                #if instrument == 'IFCB':
                    #NBS_biovol_df = NBS_biovol_df[NBS_biovol_df['logSize'] <= 17.5].reset_index(drop=True)
                lin_fit = linear_fit_func(NBS_biovol_df)
                NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])
                lin_fit_data = pd.concat([lin_fit_data, lin_fit])

    #NBSS_binned_all['Station_location'], NBSS_binned_all['midLatBin'], NBSS_binned_all['midLonBin'] = gridding_func(bin_loc, NBSS_binned_all['midLatBin'], NBSS_binned_all['midLonBin'])
    NBSS_raw = NBSS_binned_all.filter(['date_bin', 'midLatBin', 'midLonBin', 'light_cond','size_class_mid', 'ECD_mean', 'NB', 'PSD', 'Min_obs_depth', 'Max_obs_depth'], axis=1)
    NBSS_1a_raw = NBSS_raw.rename(columns={'date_bin': 'date_year_month_week', 'midLatBin': 'latitude', 'midLonBin': 'longitude', 'light_cond': 'light_condition', 'size_class_mid': 'biovolume_size_class',
                                      'NB': 'normalized_biovolume', 'ECD_mean': 'equivalent_circular_diameter_mean', 'PSD': 'normalized_abundance',  'Min_obs_depth': 'min_depth', 'Max_obs_depth': 'max_depth'})
    if light_parsing==True:
        NBSS_1a = NBSS_stats_func(NBSS_raw, light_parsing = True, bin_loc = bin_loc, group_by = group_by)
        lin_fit_1b = stats_linfit_func(lin_fit_data, light_parsing = True, bin_loc = bin_loc, group_by = group_by)
    else:
        NBSS_1a = NBSS_stats_func(NBSS_raw,  bin_loc=bin_loc, group_by=group_by)
        lin_fit_1b = stats_linfit_func(lin_fit_data, bin_loc = bin_loc, group_by = group_by)


    return NBSS_binned_all, NBSS_1a_raw, NBSS_1a,  lin_fit_1b


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
            {'year': 'first', 'month': 'first', 'midLatBin':'first', 'midLonBin': 'first', 'NB':['count', 'mean', 'std'],'ECD_mean':['mean'], 'PSD':['count', 'mean', 'std'] , 'Min_obs_depth':'first', 'Max_obs_depth':'first'}).reset_index()
    else:
        NBSS_avg = df.groupby(['date_bin', 'Station_location', 'size_class_mid']).agg(
            {'year': 'first', 'month': 'first', 'midLatBin':'first', 'midLonBin': 'first', 'NB':['count', 'mean', 'std'], 'ECD_mean':['mean'], 'PSD':['count', 'mean', 'std'],'Min_obs_depth':'first', 'Max_obs_depth':'first'}).reset_index()



    # generate a cleaner dataset (without grouping variables) to finally present 1a
    NBSS_avg.columns = NBSS_avg.columns.map('_'.join).str.removesuffix("first")
    NBSS_avg.columns = NBSS_avg.columns.str.removesuffix("_")
    NBSS_avg= NBSS_avg.filter(['year', 'month', 'midLatBin', 'midLonBin', 'light_cond', 'Min_obs_depth', 'Max_obs_depth', 'size_class_mid','NB_mean', 'NB_std', 'ECD_mean_mean','PSD_mean', 'PSD_std', 'NB_count'], axis=1)
    NBSS_avg= NBSS_avg.rename(columns={'midLatBin':'latitude', 'midLonBin': 'longitude', 'light_cond':'light_condition',
                                       'Min_obs_depth':'min_depth', 'Max_obs_depth':'max_depth', 'size_class_mid': 'biovolume_size_class',
                                       'NB_mean':'normalized_biovolume_mean', 'NB_std': 'normalized_biovolume_std', 'ECD_mean_mean': 'equivalent_circular_diameter_mean',
                                       'PSD_mean':'normalized_abundance_mean', 'PSD_std': 'normalized_abundance_std', 'NB_count':'N'})
    NBSS_avg= NBSS_avg[NBSS_avg.N !=0].reset_index(drop = True)
    return NBSS_avg

def stats_linfit_func(df, light_parsing = False, bin_loc = 1, group_by = 'yyyymm'):
    """
    Objective: create summary statistics for the linear fit dataframes
    param df: a dataframe containing the results of the linear fits for each station, time and depth for a given imaging instrument
    """
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
             'NBSS_slope': ['count', 'mean', 'std'], 'NBSS_intercept': ['count', 'mean', 'std'], 'NBSS_r2': ['count', 'mean', 'std'],
             'PSD_slope': ['count', 'mean', 'std'], 'PSD_intercept': ['count', 'mean', 'std'], 'PSD_r2': ['count', 'mean', 'std']}).reset_index()
    else:
        lin_fit_stats = df.groupby(['Station_location', 'date_bin']).agg(
            {'year': 'first', 'month': 'first', 'midLatBin': 'first', 'midLonBin': 'first', 'min_depth': 'first','max_depth': 'first',
             'NBSS_slope': ['count', 'mean', 'std'], 'NBSS_intercept': ['count', 'mean', 'std'], 'NBSS_r2': ['count', 'mean', 'std'],
             'PSD_slope': ['count', 'mean', 'std'], 'PSD_intercept': ['count', 'mean', 'std'], 'PSD_r2': ['count', 'mean', 'std']}).reset_index()

    lin_fit_stats.columns = lin_fit_stats.columns.map('_'.join).str.removesuffix("first")
    lin_fit_stats.columns = lin_fit_stats.columns.str.removesuffix("_")
    lin_fit_stats= lin_fit_stats.rename(columns={'slope_NB_count': 'N', 'midLatBin': 'latitude', 'midLonBin': 'longitude'})
    lin_fit_stats= lin_fit_stats[lin_fit_stats.columns.drop(list(lin_fit_stats.filter(regex='count')))]
    del lin_fit_stats['Station_location']
    del lin_fit_stats['date_bin']

    return lin_fit_stats




