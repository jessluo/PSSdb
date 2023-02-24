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

#1) Consolidate UVP projects before control quality check and standardization, see README file for a description of consolidation steps and objectives
# There are 4 options to consolidate a UVP project:
#a- The data owner runs the ImageJ macro on their computer and upload the metadata on Ecotaxa following the instructions of the README file and macro dialog window
#b- Metadata are missing and the script is run outside the Complex server Marie (home/plankton/uvp5/6_missions): for testing purposes, the consolidation will be performed using vignettes object_area
# Attention: Your local home path should be added to the funcs_consolidate_UVP_files script on l.259 and 300. Allowed local users:/Users/dugennem (Mathilde Dugenne) and /Users/mc4214 (Marco Corrales)
#c- Set run_macro=True. Metadata are missing and the script is run on the Complex server Marie (home/plankton/uvp5/6_missions): the consolidation will be performed by re-computing vignettes object_bru_area using a custom ImageJ macro
# Attention: The pyImageJ macro is great but it requires a lot of RAM. Do not use unless you're the only one working on the server!
#d- Ser run_macro=False. Metadata are missing and the script is run on the Complex server Marie (home/plankton/uvp5/6_missions): the consolidation will be performed by matching-up the vignettes to the native project files id (image index, blob index) and extracting the corresponding area saved in the bru files
# Attention: This is a workaround to pyimageJ that requires less memory since it is based on python code only. However, you still need an access to the server and the native UVP project folders should be clean (no compressed sub-directories, no profiles renamed with _old versions saved)!
# For disorganized UVP native projects, non-matching vignettes will be assign a size equivalent to object_area

print('Consolidating UVP particle sizes (Ecotaxa: Large particles + Ecopart: Small particles), please wait')
standardizer=[standardizer for standardizer in standardizer_files if 'UVP' in str(standardizer.stem)][0]
df_standardizer = pd.read_excel(standardizer, index_col=0)
for project in list(df_standardizer.index):
    if len(list(Path(df_standardizer.at[project,'Project_localpath']).expanduser().glob('ecotaxa_export_{}_*.tsv'.format(str(project))))):
        print('Consolidated files found. Skipping project {}'.format(str(project)))
    else:
        try:
            print('Consolidating project: {}'.format(str(project)))
            consolidate_ecotaxa_project(project_id=project, standardizer=df_standardizer, sheetname='ecotaxa', run_macro=False, upload_metadata=False,localpath_metadata=Path(cfg['raw_dir']).expanduser() / cfg['UVP_consolidation_subdir'])
        except Exception as e:
            print('Skipping consolidation of project ', str(project), '\n', e, sep='')
"""
#2) Flag samples and generate a report using the filling_standardizer_flag_func function in funcs_standardize_projects.py
print('Performing project control quality check based on the following criteria, please wait:\nFlag_missing: Missing data/metadata\nFlag_GPScoordinatesonland: GPS coordinates on land\nFlag_dubiousGPScoordinates: Dubious GPS coordinates\nFlag_count: Low ROI counts (yielding uncertainties>5%)\nFlag_artefacts: High percentage of artefacts (>20%)\nFlag_size: Multiple size calibration factors\n(0:flag, 1:no flag)')
for standardizer in natsorted(standardizer_files)[::-1]:
    for project in list(df_standardizer.index):
        # Flagging project
        report_file='report_project_'+str(project)+'.html'
        report_path=path_to_standardizer.parent / 'reports'
        try:
            print('Flagging project: {}'.format(str(project)))
            filling_standardizer_flag_func(standardizer_path=standardizer, project_id=project,report_path=report_path)
        except Exception as e:
            print('Skipping flagging of project ',str(project),'\n',e,sep='')

#3) Standardize project export files using the standardization_func function in funcs_standardize_projects.py
print('Performing project standardization')
for standardizer in natsorted(standardizer_files)[::-1]:
    for project in list(df_standardizer.index):
        try:
            print('Standardizing project: {}'.format(str(project)))
            standardization_func(standardizer_path=standardizer, project_id=project,plot='diversity')
        except Exception as e:
            print('Skipping standardization of project ', str(project),'\n',e, sep='')
"""