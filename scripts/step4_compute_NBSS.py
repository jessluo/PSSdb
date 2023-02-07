# compute NBS from step3 gridded files

import numpy as np
import pandas as pd
import statistics as st
import os

from pathlib import Path
from glob import glob
import shutil
# Config modules
import yaml  # requires installation of PyYAML package
from funcs_read import *
from funcs_NBS import *
from tqdm import tqdm
import matplotlib.pyplot as plt
#look for coastline layer and add it to that map, UNABLE to import cartopy, see issue
import cartopy.crs as ccrs


#define the instrument to calculate NBS:

instrument = input ('for which instrument do you want to calculate the size spectra? \n Enter IFCB, Zooscan or UVP ')
depth_binning  = input ('Would you like to bin the data by depth? \n Enter Y/N')
summary_stats = input ('Would you like the summary statistics of the linear regression results? \n Enter Y/N')

# processing starts here

path_to_data, id_list = proj_id_list_func(instrument, data_status ='gridded')#generate path and project ID's

dirpath = str(path_to_data) + '/NBSS_data/'
if os.path.isdir(dirpath) and len(os.listdir(dirpath)) != 0:  # and  os.path.exists(path_download)
    replace = input('There is already NBS data in ' + dirpath + ' do you want to replace the files? \n Y/N')
    if replace == 'Y':
        print('Overwriting normalized biomass data file(s), please wait')
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
        os.mkdir(dirpath + 'Lin_regress/')
        lin_fit_all = pd.DataFrame()
        for i in tqdm(id_list):
            print('calculating normalized biomass size spectra for ' + i)
            df = read_func(path_to_data, i)  # get a dataframe for each project
            # bin standardized data and obtain the metadata of the binned columns

            # 1/24/2023 commented out section below was for obtaining info from metadata about biovolume calculations. might be unnecesary now
            #path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
            # open the metadata of the standardized files
            #with open(path_to_config, 'r') as config_file:
                #cfg = yaml.safe_load(config_file)
            #path_to_metadata = glob(
                #str(Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir']) + '/**/*' + str(i) + '*metadata*')
            #metadata_std = pd.read_csv(path_to_metadata[0], sep='\t', header=0, index_col=[0])
            # remove old biovolume descriptions in the metadata file
            #metadata_std = metadata_std[metadata_std['Variables'].str.contains('Biovolume') == False]
            # concatenate metadata files to generate binned metadata
            # metadata_binned = pd.concat([metadata_std, metadata_bins], axis=0)
            # and generate the NBSS WITH THRESHOLDING INCLUDED
            if depth_binning == 'Y':
                NBS_data_binned, lin_fit_data = parse_NBS_linfit_func(df, parse_by = ['Station_location','date_bin'], depth_bin= True)
            else:
                NBS_data_binned, lin_fit_data = parse_NBS_linfit_func(df, parse_by = [ 'Station_location', 'date_bin'], depth_bin= False)

            # generate metadata for the NBS files
            Variables = NBS_data_binned.columns.to_list()
            Variable_types = NBS_data_binned.dtypes.to_list()
            Units_Values = ['cubic micrometers (range)', 'cubic micrometers',
                            'log(cubic micrometers)','yyyy or yyyymm (user defined)',
                            'lat_lon (string)', 'degree','degree',
                            '','cubic decimeters',
                            'cubic micrometers', '',
                            'cubic micrometers','cubic decimeters',
                            'counts/ cubic decimeters', 'log(counts/ cubic decimeters)']  # notice that we dont have date unit info here, fix this eventually
            Description = ['size bins in which particles are classified',
                           'minimum and maximum biovolume value of the size bin, calculated from biovolume using a projection of a sphere',
                           'logarithmic transformation of the size bin range',
                            'binned date information',
                           'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                           'latitude of the center point of the 1x1 degree cell',
                           'longitude of the center point of the 1x1 degree cell',
                           'Project identifier',
                           'Volume analyzed (not accounting for sample dilution and/or fractionation)',
                           'Sum of the biovolume of individual objects classified into a biovolume based size bin',
                           'number of objects assigned into the size bin',
                           'mean biovolume for each size bin',
                           'sum of the volumes analyzed within a 1x1 degree space, used to calculate NBSS',
                           'Normalized biomass size spectra based on biovolume',
                           'logarithmic transformation of the NBSS, use as y axis when performing size spectra analysis']
            NBS_metadata = pd.DataFrame(
                {'Variables': Variables, 'Variable_types': Variable_types, 'Units_Values': Units_Values,
                 'Description': Description})
            i = i.replace( "gridded", "NBSS")
            NBS_data_binned.to_csv(str(path_to_data) + '/NBSS_data/' + str(i), index=False)
            i = i.replace("NBSS", "metadata_NBSS")
            NBS_metadata.to_csv(str(path_to_data) + '/NBSS_data/'  + str(i), index=False)
            # append linear regression results to the dataset of slopes per instrument:
            lin_fit_all = pd.concat([lin_fit_all, lin_fit_data], ignore_index=True)


        # save the Linear regression results and metadata
        lin_fit_Variables = lin_fit_all.columns.to_list()
        lin_fit_Variable_types = lin_fit_all.dtypes.to_list()
        lin_fit_Units_Values = ['', 'degree', 'degree', 'refer to gridded metadata', '', '', '', ''] # NOTE for the future, binning information should be cinluded here
        lin_fit_Description = ['Project Identifier',
                        'latitude of the center point of the 1x1 degree cell',
                        'longitude of the center point of the 1x1 degree cell',
                        'binned date information',
                        'slope of the least-squares linear regression performed on the datased defined by the gridding and binning configuration',
                        'intercept of the least-squares linear regression performed on the datased defined by the gridding and binning configuration',
                        'Root mean square error of the least-squares linear regression performed on the datased defined by the gridding and binning configuration',
                        'Coefficiento of determination of the least-squares linear regression performed on the datased defined by the gridding and binning configuration',]
        lin_fit_metadata = pd.DataFrame(
            {'Variables': lin_fit_Variables, 'Variable_types': lin_fit_Variable_types, 'Units_Values': lin_fit_Units_Values,
                'Description': lin_fit_Description})
        # saving only needs to happen once for each instrument, yes?
        lin_fit_all.to_csv(str(path_to_data) + '/NBSS_data/Lin_regress/' + instrument + '_Linear_regression_results.csv', index=False)
        lin_fit_metadata.to_csv(str(path_to_data) + '/NBSS_data/Lin_regress/' + instrument + '_Linear_regression_metadata.csv', index=False)

        # create summary statistics for the linear regressions:
        if summary_stats == 'Y':
            lin_fit_stats = stats_linfit_func(lin_fit_all)
            lin_fit_stats.to_csv(str(path_to_data) + '/NBSS_data/Lin_regress/' + instrument + '_Linear_regression_summary.csv', index=False)


    elif replace == 'N':
        print('previously calculated normalized biomass files will be kept')

elif not os.path.exists(dirpath):
    os.mkdir(dirpath)
    os.mkdir(dirpath + 'Lin_regress/')
    lin_fit_all = pd.DataFrame()
    for i in tqdm(id_list):
        print('calculating normalized biomass size spectra for ' + i)
        df = read_func(path_to_data, i)  # get a dataframe for each project
        # bin standardized data and obtain the metadata of the binned columns

        # 1/24/2023 commented out section below was for obtaining info from metadata about biovolume calculations. might be unnecesary now
        # path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
        # open the metadata of the standardized files
        # with open(path_to_config, 'r') as config_file:
        # cfg = yaml.safe_load(config_file)
        # path_to_metadata = glob(
        # str(Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir']) + '/**/*' + str(i) + '*metadata*')
        # metadata_std = pd.read_csv(path_to_metadata[0], sep='\t', header=0, index_col=[0])
        # remove old biovolume descriptions in the metadata file
        # metadata_std = metadata_std[metadata_std['Variables'].str.contains('Biovolume') == False]
        # concatenate metadata files to generate binned metadata
        # metadata_binned = pd.concat([metadata_std, metadata_bins], axis=0)
        # and generate the NBSS WITH THRESHOLDING INCLUDED
        if depth_binning == 'Y':
            NBS_data_binned, lin_fit_data = parse_NBS_linfit_func(df, parse_by=['Station_location', 'date_bin', 'midDepthBin'])
        else:
            NBS_data_binned, lin_fit_data = parse_NBS_linfit_func(df, parse_by=['Station_location', 'date_bin'])

        # generate metadata for the NBS files
        Variables = NBS_data_binned.columns.to_list()
        Variable_types = NBS_data_binned.dtypes.to_list()
        Units_Values = ['refer to gridded metadata', 'lat_lon', 'cubic micrometers', 'cubic micrometers', 'degree',
                        'degree', '', 'cubic decimeters', 'cubic micrometers', '', 'cubic micrometers',
                        'counts/ cubic decimeters', 'log(counts/ cubic decimeters)',
                        'log (cubic micrometers)']  # notice that we dont have date unit info here, fix this eventually
        Description = ['binned date information',
                       'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                       'minimum and maximum biovolume value of the size bin, calculated from biovolume using a projection of a sphere',
                       'difference between max and min value of a size bin',
                       'latitude of the center point of the 1x1 degree cell',
                       'longitude of the center point of the 1x1 degree cell', 'Project ID in Ecotaxa',
                       'Volume analyzed (not accounting for sample dilution and/or fractionation)',
                       'Sum of the biovolume of individual objects classified into a biovolume based size bin',
                       'number of objects assigned into the size bin', 'mean biovolume for each size bin',
                       'Normalized biomass size spectra based on biovolume',
                       'logarithmic transformation of the NBSS, use as y axis when performing size spectra analysis',
                       'logarithmic transformation of the median of a size bin, use as x axis when performing size spectra analysis']
        NBS_metadata = pd.DataFrame(
            {'Variables': Variables, 'Variable_types': Variable_types, 'Units_Values': Units_Values,
             'Description': Description})

        i = i.replace("gridded", "NBSS")
        NBS_data_binned.to_csv(str(path_to_data) + '/NBSS_data/' + str(i), index=False)
        i = i.replace("NBSS", "metadata_NBSS")
        NBS_metadata.to_csv(str(path_to_data) + '/NBSS_data/' + str(i), index=False)
        # append linear regression results to the dataset of slopes per instrument:
        lin_fit_all = pd.concat([lin_fit_all, lin_fit_data], ignore_index=True)

    # save the Linear regression results
    lin_fit_Variables = lin_fit_all.columns.to_list()
    lin_fit_Variable_types = lin_fit_all.dtypes.to_list()
    lin_fit_Units_Values = ['', 'degree', 'degree', 'refer to gridded metadata', '', '', '',
                                '']  # NOTE for the future, binning information should be cinluded here
    lin_fit_Description = ['Project Identifier',
                               'latitude of the center point of the 1x1 degree cell',
                               'longitude of the center point of the 1x1 degree cell',
                               'binned date information',
                               'slope of the least-squares linear regression performed on the datased defined by the gridding and binning configuration',
                               'intercept of the least-squares linear regression performed on the datased defined by the gridding and binning configuration',
                               'Root mean square error of the least-squares linear regression performed on the datased defined by the gridding and binning configuration',
                               'Coefficiento of determination of the least-squares linear regression performed on the datased defined by the gridding and binning configuration', ]
    lin_fit_metadata = pd.DataFrame(
            {'Variables': lin_fit_Variables, 'Variable_types': lin_fit_Variable_types,
             'Units_Values': lin_fit_Units_Values,
             'Description': lin_fit_Description})
        # saving only needs to happen once for each instrument, yes?
    lin_fit_all.to_csv(str(path_to_data) + '/NBSS_data/Lin_regress/' + instrument + '_Linear_regression_results.csv', index=False)
    lin_fit_metadata.to_csv(str(path_to_data) + '/NBSS_data/Lin_regress/' + instrument + '_Linear_regression_metadata.csv', index=False)
    # create summary statistics for the linear regressions:
    if summary_stats == 'Y':
        lin_fit_stats = stats_linfit_func(lin_fit_all)
        lin_fit_stats.to_csv(str(path_to_data) + '/NBSS_data/Lin_regress/' + instrument + '_Linear_regression_summary.csv', index=False)

### NOTE: 'refer to gridded metadata' units in the metadata needs to be replaced, perhaps it might be worth keeping everythin as it was before (gridded and NBS files produced in one call)


lat = lin_fit_all['Latitude']
lon = lin_fit_all['Longitude']
slope = lin_fit_all['Slope']
intercept_t = lin_fit_all['Intercept']

intercept_plot = [x*3 for x in intercept_t]
ax = plt.axes(projection=ccrs.PlateCarree())
plt.gca().coastlines('50m')
g1=ax.gridlines(draw_labels=True)
g1.xlines = False
g1.ylines = False
plt.scatter(lon, lat, label=None, c=slope, cmap='viridis', s=intercept_plot, linewidth=0, alpha=0.5, transform=ccrs.PlateCarree())
plt.gca().coastlines('50m')
ax.set_extent([-180, 180, -90, 90])
#ax.xlabel('longitude')
#ax.ylabel('latitude')
plt.colorbar(label='slope', orientation='horizontal', anchor=(0.5, 1))
#ax.clim(min(slope), max(slope))


labels = [str(np.round(min(intercept_t), decimals=2)),
          str(np.round(st.median(intercept_t), decimals=2)),
          str(np.round(max(intercept_t), decimals=2))]

for n, area in enumerate([18, 28.5, 73.5]):
    plt.scatter([], [], c='k', alpha=0.3, s=area, label=labels[n], transform=ccrs.PlateCarree())

plt.legend(bbox_to_anchor=(0.75, 0), ncol = 3, scatterpoints=1, frameon=False,
           labelspacing=1, title='intercept')

figname = 'step4_slopes_intercept_' + instrument + '.pdf'

path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

savepath = Path(cfg['git_dir']).expanduser() / 'figures' / figname

plt.savefig(savepath)
