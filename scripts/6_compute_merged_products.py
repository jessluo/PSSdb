# goal: to produce a merged product that uses Soviadan et al. (in review) method of selecting the maximum NBSS for a given sample
# Python modules and path specifications:
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
try:
    from funcs_merged_products import *
    from funcs_NBS import *
except:
    from scripts.funcs_merged_products import *
    from scripts.funcs_NBS import *
import yaml# requires installation of PyYAML package
from pathlib import Path
#open config file and extract parameters:
path_to_git=Path('~/GIT/PSSdb').expanduser()
path_to_config = path_to_git /'scripts'/'configuration_masterfile.yaml'
import statistics as st
import os
import datetime

from glob import glob
import shutil

from tqdm import tqdm

# Processing here:

#1) calculate NBSS products no adjustments, and without considering plankton functional groups
merged_raw_list = merge_taxa_products(grouping_factors= ['date_bin', 'Station_location', 'midLatBin', 'midLonBin'])





#2) merge PFT products calculated with different size metrics and calculate
merged_adjusted_raw_list = merge_adjust_taxa_products(grouping_factors= ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'PFT'])


#2) Calculate NBSS for each of these products


