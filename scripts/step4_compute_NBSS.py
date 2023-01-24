# compute NBS from step3 gridded files

import numpy as np
import pandas as pd
import os

from pathlib import Path
from glob import glob
import shutil
# Config modules
import yaml  # requires installation of PyYAML package
from funcs_read import *
from funcs_NBS import *
from tqdm import tqdm

#define the instrument to calculate NBS:

instrument = input ('for which instrument do you want to grid the standardized dataset? \n Enter IFCB, Zooscan or UVP ')
depth_binning  = input ('Would you like to bin the data by depth? \n Enter Y/N')

# processing starts here

path_to_data, id_list = proj_id_list_func(instrument, data_status ='gridded')#generate path and project ID's

dirpath = str(path_to_data) + '/NBSS_data/'
if os.path.isdir(dirpath) and len(os.listdir(dirpath)) != 0:  # and  os.path.exists(path_download)
    replace = input('There is already gridded data in ' + dirpath + ' do you want to replace the files? \n Y/N')
    if replace == 'Y':
        print('Overwriting normalized biomass data file(s), please wait')
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
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
                NBS_data_binned = parse_NBS_func(df, parse_by=['Station_location', 'date_bin', 'midDepthBin'])
            else:
                NBS_data_binned = parse_NBS_func(df, parse_by=['Station_location', 'date_bin'])

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
            i = i.replace( "gridded", "NBSS")
            NBS_data_binned.to_csv(str(path_to_data) + '/NBSS_data/' + str(i), sep='\t')
            i = i.replace("NBSS", "metadata_NBSS")
            NBS_metadata.to_csv(str(path_to_data) + '/NBSS_data/'  + str(i), sep='\t')

    elif replace == 'N':
        print('previously calculated normalized biomass files will be kept')

elif not os.path.exists(dirpath):
    os.mkdir(dirpath)
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
            NBS_data_binned = parse_NBS_func(df, parse_by=['Station_location', 'date_bin', 'midDepthBin'])
        else:
            NBS_data_binned = parse_NBS_func(df, parse_by=['Station_location', 'date_bin'])

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
        NBS_data_binned.to_csv(str(path_to_data) + '/NBSS_data/' + str(i), sep='\t')
        i = i.replace("NBSS", "metadata_NBSS")
        NBS_metadata.to_csv(str(path_to_data) + '/NBSS_data/' + str(i), sep='\t')

### NOTE: 'refer to gridded metadata' units in the metadata needs to be replaced, perhaps it might be worth keeping everythin as it was before (gridded and NBS files produced in one call)

