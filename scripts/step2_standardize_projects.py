##  Objective: Flag and standardize projects listed and exported on steps 0 and 1

## Requirements: Completion of standardizer spreadsheets (step0.list_projects.py) and projects export (step1.export_projects.py)

## Python modules

## Path modules
from pathlib import Path
import os

## Data handling
import pandas as pd
from natsort import natsorted

## standardization functions
from funcs_standardize_projects import *

# Workflow starts here:
path_to_standardizer=Path('~/GIT/PSSdb/raw').expanduser()
standardizer_files=list(path_to_standardizer.glob('project_*_standardizer.xlsx'))
standardizer_files=[path for path in standardizer_files if 'dashboard' not in str(path)]
#1) Flag samples and generate a report using the filling_standardizer_flag_func function in funcs_standardize_projects.py
print('Performing project control quality check based on the following criteria, please wait:\nFlag_missing: Missing data/metadata\nFlag_GPScoordinatesonland: GPS coordinates on land\nFlag_dubiousGPScoordinates: Dubious GPS coordinates\nFlag_count: Low ROI counts (yielding uncertainties>5%)\nFlag_artefacts: High percentage of artefacts (>20%)\nFlag_size: Multiple size calibration factors\n(0:flag, 1:no flag)')
for standardizer in natsorted(standardizer_files)[::-1]:
    df_standardizer=pd.read_excel(standardizer,index_col=0)
    for project in list(df_standardizer.index):
        # Flagging project
        report_file='report_project_'+str(project)+'.html'
        report_path=path_to_standardizer.parent / 'reports'
        try:
            print('Project: {}'.format(str(project)))
            filling_standardizer_flag_func(standardizer_path=standardizer, project_id=int(project),report_path=report_path)
        except Exception as e:
            print('Skipping flagging of project ',str(project),'\n',e,sep='')

#2) Standardize project export files using the standardization_func function in funcs_standardize_projects.py
print('Performing project standardization')
for standardizer in natsorted(standardizer_files)[::-1]:
    df_standardizer = pd.read_excel(standardizer, index_col=0)
    
    for project in list(df_standardizer.index):
        try:
            print('Project: {}'.format(str(project)))
            standardization_func(standardizer_path=standardizer, project_id=project,plot='diversity')
        except Exception as e:
            print('Skipping standardization of project ', str(project),'\n',e, sep='')