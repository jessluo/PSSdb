# compute NBS from step3 gridded files

import numpy as np
import pandas as pd
import os

from pathlib import Path
from glob import glob
import shutil
# Config modules
import yaml  # requires installation of PyYAML package
from read_funcs import *
from NBS_funcs import *

#define the instrument to calculate NBS:

instrument = 'IFCB'
ignore_depth = 'yes'

# processing starts here

path_to_data, id_list = proj_id_list_func(instrument, data_status ='gridded')#generate path and project ID's
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
    df = read_func(path_to_data, i)# get a dataframe for each project
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

