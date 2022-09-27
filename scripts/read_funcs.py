
# 2) function to create list of projects. Inputs: instrument
import os


def proj_id_list_func(instrument, data_status):
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
    elif data_status == 'gridded':
        path_to_data = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / instrument / cfg['gridded_subdir']
    elif data_status == 'NBSS':
        path_to_data = Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'] / instrument / cfg['gridded_subdir'] / 'NBSS_data'
    # create a list  projects that we have access to, based on project_list_all.xlsx
    path_to_proj_id_list = Path(cfg['git_dir']).expanduser() / cfg['proj_list']
    proj_list = pd.read_excel(path_to_proj_id_list)
    id_list = proj_list['Project_ID'].loc[ (proj_list['Instrument'] == str(instrument)) & (proj_list['PSSdb_access'] == True)].tolist() # instrument type filtered here

    return path_to_data, id_list


# 4) Create clean dataframes in python, input: a path and an ID. modify this to prevent removing rows
def read_func(path_to_data, ID):
    """
    Objective: returns a dataframe for each STANDARDIZED project
    :param path_to_data: path where the standardized files are stored
    :param ID: id number in ecotaxa of eatch project
    """
    import pandas as pd
    from glob import glob
    search_pattern = path_to_data / ("*" + str(ID) + "*.tsv") #'[^metadata]'
    filename = glob(str(search_pattern))
    filename = [x for x in filename if 'metadata' not in x]
    df = pd.read_csv(filename[0], sep='\t',  header=0, index_col=[0]) #
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
