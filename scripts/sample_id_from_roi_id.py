##Objective: script to correct ALREADY DOWNLOADED  IFCB dashboard files, to generate sample ID

import pandas as pd
from pathlib import Path
import yaml
import os
from tqdm import tqdm

#setup path to downloaded data
path_to_config=Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_to_git=Path(cfg['git_dir']).expanduser()
path_to_IFCB_downloads = path_to_git / cfg['dataset_subdir'] / cfg['IFCB_dir']

files_list = []
for path, subdirs, files in os.walk(path_to_IFCB_downloads):
    for name in files:
        files_list.append(os.path.join(path, name))

files_list = [s for s in files_list if '.tsv' in s]

for d in tqdm(files_list):
    print('looking for sample_id column in' + d)
    df = pd.read_csv(d, sep='\t',  header=0)
    if 'sample_id' not in df.columns.values:
        print('generating sample_id')
        df['sample_id'] = ''
        for i in range(0, len(df)):
            roi_id = df.loc[i, 'roi_id']
            id_split = roi_id.split('_')
            sample_id = id_split [0] + '_' + id_split[1]
            df.loc[i, 'sample_id'] = sample_id
        df.to_csv(d, sep='\t')
    else:
        print(d + ' already has sample_id')
