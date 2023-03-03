
# 2) function to create list of projects. Inputs: instrument
import os
import pandas as pd
import yaml
from pathlib import Path

def proj_id_list_func(instrument, data_status, big_grid = False):
    """
    Objective:
    generate a list and create the path  of the standardized Ecotaxa projects
    :param instrument: the device used for image adcquisition. important since the path will change according to it
    :return path_to_data: string where the tsv files are stores
    :return file_list: list with the filenames that have data
    """
    # returns path to project data and list of projects based on intrument (IFCB, Zooscan, UVP)
    # read git-tracked config file (text file) with inputs:  project ID, output directory
    path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    # read config file (text file) with inputs:  project ID, output directory
    # prepare access based on path stored in the yaml config file and instrument type
    path_to_ecotaxa = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / cfg['Ecotaxa_subdir']
    path_to_consolidated_UVP= Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / cfg['UVP_consolidation_subdir']/ cfg['UVP_consolidation_subdir'] #tiny glitch here (repeated subdirectories)
    path_to_IFCBdashboard = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / cfg['IFCB_dir']
    if data_status == 'standardized':
        if instrument  == 'Zooscan':
            path_to_data =  path_to_ecotaxa / instrument
            file_list = os.listdir(path_to_data)
            files_data = [(str(path_to_data / x)) for x in file_list if not 'metadata' in x and '.csv' in x]
        elif instrument == 'UVP':
            file_list = os.listdir(path_to_consolidated_UVP)
            files_data = [(str(path_to_consolidated_UVP / x)) for x in file_list if not 'metadata' in x and '.csv' in x]
        elif instrument  == 'IFCB':
            #get data from ecotaxa
            path_to_ecotaxa_data = path_to_ecotaxa / instrument / instrument # little glitch here
            ecotaxa_list = os.listdir(path_to_ecotaxa_data)
            ecotaxa_data = [(str(path_to_ecotaxa_data / x)) for x in ecotaxa_list  if not 'metadata' in x and '.csv' in x]
            #get data from dashboards
            dashboard_list = [f for f in os.listdir(path_to_IFCBdashboard) if not f.startswith('.')]
            dashboard_data = []
            for i in dashboard_list:
                proj_list = os.listdir(path_to_IFCBdashboard / i)
                proj_data = [(str(path_to_IFCBdashboard / i / x )) for x in proj_list if not 'metadata' in x and '.csv' in x]
                dashboard_data = dashboard_data + proj_data
            files_data = ecotaxa_data + dashboard_data
    elif data_status == 'gridded':
        path_to_gridded = Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir'] / instrument
        file_list = os.listdir(path_to_gridded)
        if big_grid == False:
            files_data = [(str(path_to_gridded / x)) for x in file_list if not 'metadata' in x and not 'grid_N_' in x and  '.csv' in x]
        elif big_grid == True:
            files_data = [(str(path_to_gridded / x)) for x in file_list if 'grid_N_' in x]
    return files_data


def depth_parsing_func(df, instrument):
    '''
    Objective: restrict samples based on the agreed parameters for first PSSdb release.
    param df: a standardized dataframe
    param instrument: the imaging system used to collect the data
    '''
    if instrument == 'UVP':
        df_d = df.drop(df[(df['Depth_min'] > 200)].index)
        df_d['Min_obs_depth'] = 0
        df_d['Max_obs_depth'] = 200
    elif instrument == 'Zooscan':
        df_d = df.drop(df[(df['Depth_min'] > 250) | (df['Depth_max'] > 250)].index)
        df_d['Min_obs_depth'] = 0
        df_d['Max_obs_depth'] = 250
    elif instrument == 'IFCB': # IFCB depth specification occurs in step4 due to the need of getting a depth interval of ALL the projects
        df_d = df.drop(df[(df['Depth_min'] > 200)].index)
        df_d['Min_obs_depth'] = min(df_d['Depth_min'])
        df_d['Max_obs_depth'] = max(df_d['Depth_max'])
    return df_d

# 4) Create clean dataframes in python, input: a path and an ID. modify this to prevent removing rows
def read_func(path_to_data, ID, data_status='gridded'):
    """
    Objective: returns a dataframe for each STANDARDIZED project
    :param path_to_data: path where the standardized files are stored
    :param ID: id number in ecotaxa of eatch project
    """
    import pandas as pd
    from glob import glob
    filename = path_to_data / ID #  DEPRECATED: 1/12/2023 ("*" + str(ID) + "*.tsv") #'[^metadata]'
    if data_status == 'standardized':
        df = pd.read_csv(filename, sep='\t',  header=0) #
    else:
        df = pd.read_csv(filename, header=0)
    df = df.loc[:, ~df.columns.str.match("Unnamed")]
    # convert taxonomic Id column into categorical.
    # NOTE all of these column names might need to be modified after standardization of datasets
    #df.loc[:, cat].astype('category')
    #df['proj_ID'] = [ID] * (len(df))
    return df


def df_list_func(instrument, data_status):
    """
    Objective:
    generate a list and create the path  of the standardized Ecotaxa projects
    :param instrument: the device used for image adcquisition. important since the path will change according to it
    :return path_to_data: string where the tsv files are stores
    :return id_list: list with the ecotaxa project identifiers (might not be relevant once we expand
    to get data from other platforms)
    """
    import pandas as pd
    import yaml
    from pathlib import Path
    # returns path to project data and list of projects based on intrument (IFCB, Zooscan, UVP)
    # read git-tracked config file (text file) with inputs:  project ID, output directory
    path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    # read config file (text file) with inputs:  project ID, output directory
    # prepare storage based on path stored in the yaml config file and instrument type
    if data_status == 'standardized':
        path_to_data = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / instrument
        df_list = os.listdir(path_to_data)
        df_list = [x for x in df_list if 'metadata' not in x]
        df_list = [x for x in df_list if '.tsv' in x]
        df_list = [str(path_to_data) +'/'+ x for x in df_list]
    elif data_status == 'gridded':
        path_to_data = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / instrument / cfg['gridded_subdir']
        df_list = os.listdir(path_to_data)
        df_list = [x for x in df_list if 'metadata' not in x]
        df_list = [x for x in df_list if '.tsv' in x]
        df_list = [str(path_to_data) +'/'+ x for x in df_list]
    elif data_status == 'NBSS':
        path_to_data = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / instrument / cfg['gridded_subdir'] / 'NBSS_data'
        df_list = os.listdir(path_to_data)
        df_list = [x for x in df_list if 'metadata' not in x]
        df_list = [x for x in df_list if '.tsv' in x]
        df_list = [str(path_to_data) +'/'+ x for x in df_list]

    return path_to_data, df_list
