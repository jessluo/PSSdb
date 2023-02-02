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
    df = pd.read_csv(d, sep='\t',  header=0, index_col=[0])
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    if 'sample_id' not in df.columns.values:
        print('generating sample_id')
        df['sample_id'] = ''
        for i in range(0, len(df)):
            roi_id = df.loc[i, 'roi_id']
            id_split = roi_id.split('_')
            sample_id = id_split [0] + '_' + id_split[1]
            df.loc[i, 'sample_id'] = sample_id
        col = df.pop('sample_id')
        df.insert(12, col.name, col)
        df.to_csv(d, sep='\t', index = False)
    else:
        print(d + ' already has sample_id, checking and relocating sample_id column')
        col_index_id = df.columns.get_loc("project_ID")
        print('relocating sample_id for ' +d)
        sample_id_col = df.pop('sample_id')
        df.insert(col_index_id+1, sample_id_col.name, sample_id_col)

        print('relocating sampling columns info for ' + d)
        minArea_col = df.pop('minBlobArea')
        sampling_index_id = df.columns.get_loc("class_5_score")
        df.insert(sampling_index_id + 1, minArea_col.name, minArea_col)
        sampling_desc = ['minBlobArea','PMTAhighVoltage', 'PMTBhighVoltage', 'PMTChighVoltage', 'SyringeSampleVolume',
                             'PMTAtriggerThreshold_DAQ_MCConly', 'PMTBtriggerThreshold_DAQ_MCConly', 'PMTCtriggerThreshold_DAQ_MCConly']
        for i, s in enumerate(sampling_desc):
            if s == 'minBlobArea':
                pass
            else:
                col = df.pop(s)
                n = df.columns.get_loc(sampling_desc[i-1])
                df.insert(n+1, col.name, col)
        df.to_csv(d, sep='\t', index = False)

