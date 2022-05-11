



from pathlib import Path
from glob import glob

# Config modules
import yaml # requires installation of PyYAML package

import pandas as pd


## Workflow starts here

# read git-tracked config file (text file) with inputs:  project ID, output directory
path_to_config=Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
# read config file (text file) with inputs:  project ID, output directory

path_to_git=Path(cfg['git_dir']).expanduser()

# prepare storage based on path stored in the yaml config file and instrument type

path_to_data = Path(cfg['git_dir']).expanduser() / cfg['dataset_subdir']


### need to create an excel sheet to provide input for below calculation,
#              e.g. needs to contain pixel to size conversion, volume information
#               which categories to exclude (artefact???)

#we should have an excel sheet that has a test set of projects, and one that has all projects for production


test = 1

if test == 1:
    #read in the test excel sheet to only work on few projects
if test == 0:
    print("This is the whole bunch, you really want to do it???")

id_dict = {3657 : ["ml"]} #just an example, this info should be in excel sheet that is read in and converted to dict


for id in id_dict:
    search_pattern = path_to_data / ("*" + str(id) + "*")
    volume_info = id_dict[id][0]

    filename = glob(str(search_pattern))
    df = pd.read_csv(filename) #consider throwing out info we do not need anymore

    #need to convert pixel size to area size

    #consider that a size spectrum for zooscan needs to be calculated per sample

    #use ESD, area or biovolume to calculate the size spectrum???

    ##size spectrum per sample....


    #concatenate the resulting df for the different projects in one big df_sizespectrum, which we safe outside the loop





#df_sizespectrum.to_csv("")









