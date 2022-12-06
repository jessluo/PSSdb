##  Objective: Flag and standardize projects in batch

## Requirements: Completion of standardizer spreadsheets (step0.list_projects.py) and projects export (step1.export_projects.py)

## Python modules

## Path modules
from pathlib import Path
import os

## Data handling
import pandas as pd

## standardization functions
from funcs_standardize_projects import *

# Workflow starts here:
path_to_standardizer=Path('~/GIT/PSSdb/raw').expanduser()
standardizer_files=list(path_to_standardizer.glob('project_*_standardizer.xlsx'))
standardizer_files=[path for path in standardizer_files if 'dashboard' not in str(path)]

print('Performing project control quality check based on the following criteria, please wait:\nFlag_missing: Missing data/metadata\nFlag_GPSonland: GPS coordinates on land\nFlag_count: Low ROI counts (yielding uncertainties>5%)\nFlag_artefacts: High percentage of artefacts (>20%)\nFlag_size: Multiple size calibration factors\n(0:flag, 1:no flag)')
for standardizer in standardizer_files:
    df_standardizer=pd.read_excel(standardizer,index_col=0)
    for project in list(df_standardizer.index):
        # Flagging project
        report_file='report_project_'+str(project)+'.html'
        report_path=path_to_standardizer.parent / 'reports' / str(df_standardizer.loc[project]['Instrument']) / report_file
        if report_path.is_file()==False:
           try:
              print('Project: {}'.format(str(project)))
              filling_standardizer_flag_func(standardizer_path=standardizer, project_id=int(project),report_path=report_path.parent)
           except Exception as e:
               print('Skipping flagging of project ',str(project),'\n',e,sep='')
        else:
            pass

print('Performing project standardization')
for standardizer in standardizer_files:
    df_standardizer = pd.read_excel(standardizer, index_col=0)
    
    for project in list(df_standardizer.index):
        # Standardizing project
        standard_file = 'standardized_project_{}.tsv'.format(str(project))
        standard_path=path_to_standardizer/ 'raw_standardized' /str(df_standardizer.loc[project]['Instrument']) /standard_file
        if standard_path.is_file() == False:
            try:
                print('Project: {}'.format(str(project)))
                standardization_func(standardizer_path=standardizer, project_id=int(project),plot='diversity')
            except Exception as e:
                print('Skipping standardization of project ', str(project),'\n',e, sep='')
        else:
            pass
