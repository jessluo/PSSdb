
 3) Standardize all files:
def mass_st_func(instrument):
    """"
    Objective:
    Complete standardization of all projects for a given instrument
    """
    #import functions to 1. get project ID lists and
    from step2_standardize_data import proj_id_list_func as id_list
    #get the standardizer function
    from funcs_export_standardizer import standardization_func as standardize
    from pathlib import Path
    import yaml
    import glob
    instrument = str(instrument)
    #define the paths to the config file and the standardizer using the yaml file
    path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    # read config file (text file) with inputs:  project ID, output directory
    standardizer_path = glob.glob(str(Path(cfg['raw_dir']).expanduser()) + '/' + '*' + instrument+ '*xlsx')[0]
    config_path = str(Path('~/GIT/PSSdb/scripts/Ecotaxa_API_pw.yaml').expanduser())
    path_to_data, ID_list = id_list(instrument, testing=False)
    #for i in ID_list:
        #standardize(config_path, standardizer_path, i)
    if instrument == 'IFCB':# this removal is only for testing,
        #since the standardizer should be able to deal with projects with empty rows
        IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
        for element in IFCB_empty_rows:
            if element in ID_list:
                ID_list.remove(element)
        for i in ID_list:
            standardize(config_path, standardizer_path, i)

    else:
        for i in ID_list:
            standardize(config_path, standardizer_path, i)

