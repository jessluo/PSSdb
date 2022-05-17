#GOAL: function to unzip folders and obtain the tsv files of a given instrument

def extract_tsv(instrument):
    # inst: string that corresponds to the instrument to use: IFCB, Zooscan,
    import os
    import zipfile
    from pathlib import Path # Handling of path object
    import shutil # Delete uncompressed export zip folder
    # Config modules
    import yaml # requires installation of PyYAML package

    ## Workflow starts here

    # read git-tracked config file (text file) with inputs:  project ID, output directory
    path_to_config=Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
    with open(path_to_config ,'r') as config_file:
        cfg = yaml.safe_load(config_file)
    # read config file (text file) with inputs:  project ID, output directory

    # prepare storage based on path stored in the yaml config file and instrument type

    path_to_data = Path(cfg['git_dir']).expanduser() / cfg['dataset_subdir'] / instrument

    for file in os.listdir(path_to_data):
        if file.endswith('.zip'):
            with zipfile.ZipFile((path_to_data / file), 'r') as zip_ref:
                zip_ref.extractall(path_to_data)
