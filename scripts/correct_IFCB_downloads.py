import shutil
import yaml
from tqdm import tqdm
from pathlib import Path
from glob import glob
import pandas as pd
import os


# set paths:
path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
# open the metadata of the standardized files
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
all_files_data = str(Path(cfg['raw_dir']).expanduser()) + '/IFCB_dashboard_projects_list.xlsx'
IFCB_dashboard_data_path = str(Path(cfg['IFCB_dir']).expanduser())

path_IFCB_dashboard_data = Path(cfg['raw_dir']).expanduser() / cfg['IFCB_dir']

file_list = glob(str(path_IFCB_dashboard_data) + '*/**/*.tsv', recursive=True)

projects_v3= ['arctic', 'harpswell', 'EXPORTS', 'NESLTER_broadscale', 'NESLTER_transect', 'SPIROPA']

projects_v2= ['santa-cruz-municipal-wharf', 'san-francisco-pier-17']

print('deleting data from IFCB014 and setting appropiate pixel_to_micron_values')
for i in tqdm(file_list):
    try:
        df = pd.read_csv(i, sep= '\t', index_col=0)
        print ('checking file ' + i)
        df = df[~df.sample_id.str.contains('|'.join(['IFCB014']))].reset_index(drop=True)
        print ('this file has a scale of ')
        print(3.400 if any(proj in i for proj in projects_v2) else 2.77 if any(proj in i for proj in projects_v3) else df.loc[0, 'pixel_to_micron'])
        df['pixel_to_micron'] = 3.4 if any(proj in i for proj in projects_v2) else 2.77 if any(proj in i for proj in projects_v3) else df.loc[0,'pixel_to_micron']
        df.to_csv(i, sep='\t', index=True)
    except:
        pass


