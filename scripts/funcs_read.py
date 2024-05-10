
# 2) function to create list of projects. Inputs: instrument
import os
import pandas as pd
import yaml
from pathlib import Path
from glob import glob

Instrument_dict = {'IFCB': ['IFCB'], 'UVP': ['UVP'], 'Scanner': ['Scanner', 'ZooCam', 'Zooscan'],'PlanktoScope': ['PlanktoScope'], 'FlowCam': ['FlowCam']}


def proj_id_list_func(instrument, data_status):
    """
    Objective:
    read standaroizer files, and generate an instrument specific list of either standardized or gridded files
    :param instrument: the device used for image adcquisition. important since the path will change according to it
    :return path_to_data: string where the tsv files are stores
    :return file_list: list with the filenames that have data
    """
    # read config file
    path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    # generate paths to standardized, gridded files, and the standardizer spreadsheets
    path_to_standardized = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir']
    path_to_gridded = Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir']
    #11/27/2024 changes: flag files will be used instead of the standardizer to generate the files to be gridded
    #path_to_standardizer=Path('~/GIT/PSSdb/raw').expanduser()
    #make a dataframe with all the projects and then subset based on the instrument:
    #standardizer_files=list(path_to_standardizer.glob('project_*_standardizer.xlsx'))
    #standardizer_df = pd.concat(map((lambda path: (pd.read_excel(path))), standardizer_files))

    path_flags = Path(cfg['raw_dir']).expanduser() / cfg['flag_subdir']
    flags_file_list = glob(str(path_flags) + '*/**/*flags.csv', recursive=True)
    flags_df = pd.concat(map((lambda path: (pd.read_csv(path, sep = ','))), flags_file_list))
    flags_df_subset = pd.concat([flags_df.loc[(((flags_df.Flag == 0) & (flags_df.Overrule == False)) | ((flags_df.Flag == 1) & (flags_df.Overrule == True))) & (flags_df.Instrument==i)]  for i in Instrument_dict[instrument]]) # subsetting files that passed through the QC and data owner check
    # generate the file list (standardized or gridded files)

    #standardizer_df_subset = standardizer_df[standardizer_df.Instrument == i].reset_index(drop=True)
    if data_status == 'standardized':
        #file_list = [file for project in [glob(str(path_to_standardized) + '*/**/*' + i + '*/**/*' +str(proj_id) + '*.csv', recursive=True) for proj_id in standardizer_df_subset.Project_ID] for file in project] # https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
        file_list = [os.path.expanduser(path) for path in flags_df_subset.Sample_localpath.tolist()]
        file_list_local = list(set([file for project in [glob(str(path_to_standardized) + '*/**/*standardized_project_' + str(proj_id) + '_*.csv', recursive=True) for proj_id in flags_df_subset.Project_ID] for file in project]))
        file_list = list(set([file for file in file_list if file in file_list_local])) # local file search is necessary for testing, since a personal computer won't have all the standardized files
        #print(file_list)
    elif data_status == 'gridded':
        file_list = list(set([file for project in [glob(str(path_to_gridded) + '*/**/*_gridded_project_' + str(proj_id) + '_*.csv', recursive=True) for proj_id in flags_df_subset.Project_ID] for file in project]))
    return file_list

def proj_id_list_func_old(instrument, data_status, big_grid = False):
    """
    Objective:
    generate a list and create the path  of the standardized Ecotaxa projects
    :param instrument: the device used for image adcquisition. important since the path will change according to it
    :return path_to_data: string where the tsv files are stores
    :return file_list: list with the filenames that have data
    """
    # returns path to project data and list of projects based on intrument (IFCB, Zooscan, UVP)
    # read git-tracked config file (text file) with inputs:  project ID, output directory
    path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    # read config file (text file) with inputs:  project ID, output directory
    # prepare access based on path stored in the yaml config file and instrument type
    path_standardized = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir']
    #path_to_ecotaxa = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / cfg['Ecotaxa_subdir']
    #path_to_consolidated_UVP= Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / cfg['UVP_consolidation_subdir']/ cfg['UVP_consolidation_subdir'] #tiny glitch here (repeated subdirectories)
    #path_to_IFCBdashboard = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / cfg['IFCB_dir']
    if instrument == 'UVP':
        id = 'ecopart'
    else:
        id = instrument
    if data_status == 'standardized':
        file_list = glob(str(path_standardized) + '*/**/*' + id +'*/**/*.csv', recursive=True)
        files_data = [x for x in file_list if not 'metadata' in x]
        files_data_clean = []
        [files_data_clean.append(x) for x in files_data if x not in files_data_clean]
    elif data_status == 'gridded':
        path_to_gridded = Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir']
        file_list = glob(str(path_to_gridded) + '/**/*' + instrument +'*/**/*.csv', recursive=True)
        if big_grid == False:
            files_data = [(str(path_to_gridded / x)) for x in file_list if not 'metadata' in x and not 'grid_N_' in x and  '.csv' in x]
            files_data_clean = []
            [files_data_clean.append(x) for x in files_data if x not in files_data_clean]
        elif big_grid == True:
            files_data = [(str(path_to_gridded / x)) for x in file_list if 'grid_N_' in x]
            files_data_clean = []
            [files_data_clean.append(x) for x in files_data if x not in files_data_clean]
    return files_data_clean


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
    elif instrument == 'Scanner':
        df_d = df.drop(df[(df['Depth_min'] > 250) | (df['Depth_max'] > 250)].index)
        df_d['Min_obs_depth'] = 0
        df_d['Max_obs_depth'] = 250
    else: # IFCB depth specification occurs in step4 due to the need of getting a depth interval of ALL the projects
        df_d = df.drop(df[(df['Depth_min'] > 200)].index)
        df_d['Min_obs_depth'] = min(df_d['Depth_min'])
        df_d['Max_obs_depth'] = max(df_d['Depth_max'])
    df_d = df_d.reset_index(drop=True)
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
    path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
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
