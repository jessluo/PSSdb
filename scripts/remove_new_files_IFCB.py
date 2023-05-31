
import pandas as pd
import os
import yaml
from tqdm import tqdm
from pathlib import Path
from glob import glob


project_ID = input('type IFCB dashboard project ID')

path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

path_to_standardized_data = str(Path(cfg['raw_dir']).expanduser() / Path(cfg['standardized_raw_subdir']) / Path(cfg['IFCB_dir']))
data_list= glob(path_to_standardized_data+ '/' + project_ID +'/*', recursive=True)

path_to_flags = str(Path(cfg['raw_dir']).expanduser() / Path(cfg['flag_subdir']) / Path(cfg['IFCB_dir']))
path_to_flag = path_to_flags + '/project_' + project_ID + '_flags.csv'
flag_file = pd.read_csv(path_to_flag)
file_list = flag_file.Sample.values.tolist()
for i in tqdm(data_list):
    print('checking data in ' + i)
    df = pd.read_csv(i)
    df = df[df['Sample'].isin(file_list)].reset_index(drop=True)
    if len(df) == 0:
        os.remove(i)
    else:
        df.to_csv(i, index=False)

