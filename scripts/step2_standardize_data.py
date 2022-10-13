##  Objective: Flag and standardize projects in batch

## Requirements: Completion of standardizer spreadsheets (step0.list_projects.py) and projects export (step1.export_projects.py)

## Python modules

## Path modules
from pathlib import Path
import os

## Data handling
import pandas as pd

## standardization functions
from funcs_project_standardizer import *

# Workflow starts here:
path_to_standardizer=Path('~/GIT/PSSdb/raw').expanduser()
standardizer_files=list(path_to_standardizer.glob('project_*_standardizer.xlsx'))

for standardizer in standardizer_files:
    df_standardizer=pd.read_excel(standardizer,index_col=0)
    for project in list(df_standardizer.index):
        # Flagging project
        report_file='report_project_'+str(project)+'.html'
        report_path=path_to_standardizer.parent / 'reports' / report_file
        if report_path.is_file()==False:
           try:
              filling_standardizer_flag_func(standardizer_path=standardizer, project_id=int(project),report_path=report_path.parent)
           except:
               print('Skipping flagging of project ',str(project),sep='')
        else:
            pass

        # Standardizing project
        standard_file = 'standardized_export_{}.tsv'.format(str(project))
        standard_path=path_to_standardizer/ 'raw_standardized' /str(df_standardizer.loc[project]['Instrument']) /standard_file
        if standard_path.is_file() == False:
            try:
                standardization_func(standardizer_path=standardizer, project_id=int(project),plot='diversity')
            except:
                print('Skipping standardization of project ', str(project), sep='')
        else:
            pass