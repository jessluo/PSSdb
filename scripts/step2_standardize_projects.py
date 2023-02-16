##  Objective: Flag and standardize projects listed and exported on steps 0 and 1

## Requirements: Completion of standardizer spreadsheets (step0.list_projects.py) and projects export (step1.export_projects.py)

## Python modules

## Path modules
from pathlib import Path
import os
import yaml
## Data handling
import pandas as pd
from natsort import natsorted

## standardization functions
from funcs_consolidate_UVP_files import *


## standardization functions
from funcs_standardize_projects import *

path_to_config=Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_to_git=Path(cfg['git_dir']).expanduser()

# Workflow starts here:
path_to_standardizer=Path('~/GIT/PSSdb/raw').expanduser()
standardizer_files=list(path_to_standardizer.glob('project_*_standardizer.xlsx'))
standardizer_files=[path for path in standardizer_files if 'dashboard' not in str(path)]

#1) Consolidate UVP projects before control quality check and standardization
print('Consolidating UVP particle sizes (Ecotaxa: Large particles, Ecopart: Small particles), please wait')
standardizer=[standardizer for standardizer in standardizer_files if 'UVP' in str(standardizer.stem)][0]
df_standardizer = pd.read_excel(standardizer, index_col=0)
for project in list(df_standardizer.index):
    try:
        print('Project: {}'.format(str(project)))
        consolidate_ecotaxa_project(project_id=project, standardizer=df_standardizer, sheetname='ecotaxa', run_macro=True, upload_metadata=False,localpath_metadata=Path(cfg['raw_dir']).expanduser() / cfg['UVP_consolidation_subdir'])
    except Exception as e:
        print('Skipping consolidation of project ', str(project), '\n', e, sep='')

#2) Flag samples and generate a report using the filling_standardizer_flag_func function in funcs_standardize_projects.py
print('Performing project control quality check based on the following criteria, please wait:\nFlag_missing: Missing data/metadata\nFlag_GPScoordinatesonland: GPS coordinates on land\nFlag_dubiousGPScoordinates: Dubious GPS coordinates\nFlag_count: Low ROI counts (yielding uncertainties>5%)\nFlag_artefacts: High percentage of artefacts (>20%)\nFlag_size: Multiple size calibration factors\n(0:flag, 1:no flag)')
for standardizer in natsorted(standardizer_files)[::-1]:
    df_standardizer=pd.read_excel(standardizer,index_col=0)
    for project in list(df_standardizer.index):
        # Flagging project
        report_file='report_project_'+str(project)+'.html'
        report_path=path_to_standardizer.parent / 'reports'
        try:
            print('Project: {}'.format(str(project)))
            filling_standardizer_flag_func(standardizer_path=standardizer, project_id=project,report_path=report_path)
        except Exception as e:
            print('Skipping flagging of project ',str(project),'\n',e,sep='')

#3) Standardize project export files using the standardization_func function in funcs_standardize_projects.py
print('Performing project standardization')
for standardizer in natsorted(standardizer_files)[::-1]:
    df_standardizer = pd.read_excel(standardizer, index_col=0)
    
    for project in list(df_standardizer.index):
        try:
            print('Project: {}'.format(str(project)))
            standardization_func(standardizer_path=standardizer, project_id=project,plot='diversity')
        except Exception as e:
            print('Skipping standardization of project ', str(project),'\n',e, sep='')