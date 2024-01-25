##  Objective: Flag and standardize projects listed and exported on steps 0 and 1

## Requirements: Completion of standardizer spreadsheets (step0.list_projects.py) and projects export (step1.export_projects.py)

## Python modules

## Path modules
from pathlib import Path
import os
import yaml
import sys
# Progress bar modules
#import progressbar # Attention: Use pip install progressbar2 to install
from tqdm import tqdm

## Data handling
import pandas as pd
from natsort import natsorted

## standardization functions
try:
    from funcs_consolidate_UVP_files import *
except:
    from scripts.funcs_consolidate_UVP_files import *

## standardization functions
try:
    from funcs_standardize_projects import *
except:
    from scripts.funcs_standardize_projects import *

## Updating annotations functions
try:
    from funcs_export_projects import *
except:
    from scripts.funcs_export_projects import *


path_to_config=Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_to_git=Path(cfg['git_dir']).expanduser()

# Workflow starts here:
path_to_standardizer=Path(cfg['raw_dir']).expanduser()
standardizer_files=list(path_to_standardizer.glob('project_*_standardizer.xlsx'))
standardizer_files=[path for path in standardizer_files if 'dashboard' not in str(path)]

#1) Update taxonomic annotations for Ecotaxa projects to ensure good validation control quality check
confirmation=input("Do you wish to update the taxonomic annotations of Ecotaxa projects before proceeding to project control quality check? Enter Y or N\n")
if confirmation=='Y':
    print('Updating taxonomic annotations of Ecotaxa projects before control quality check, please wait')
    for standardizer in natsorted(standardizer_files)[::-1]:
        df_standardizer = pd.read_excel(standardizer, index_col=0)
        project_ecotaxa=[project for project in list(df_standardizer.index) if str(Path(df_standardizer.at[project,'Project_localpath']).stem).lower()=='ecotaxa']
        for project in project_ecotaxa:
            # Updating annotations of project
            try:
                Ecotaxa_update_annotations(project_id=project, localpath=Path(df_standardizer.at[project,'Project_localpath']) / df_standardizer.at[project,'Instrument'])
            except Exception as e:
                print('Skipping update of project ',str(project),'\n',e,sep='')
else:
    print('Skipping Ecotaxa projects update')

#2) Consolidate UVP projects before control quality check and standardization, see README file for a description of consolidation steps and objectives
# There are 4 options to consolidate a UVP project:
#a- The data owner runs the ImageJ macro on their computer and upload the metadata on Ecotaxa following the instructions of the README file and macro dialog window
#b- Metadata are missing and the script is run outside the Complex server Marie (home/plankton/uvp5/6_missions): for testing purposes, the consolidation will be performed using vignettes object_area
# Attention: Your local home path should be added to the funcs_consolidate_UVP_files script on l.259 and 300. Allowed local users:/Users/dugennem (Mathilde Dugenne) and /Users/mc4214 (Marco Corrales)
#c- Set run_macro=True. Metadata are missing and the script is run on the Complex server Marie (home/plankton/uvp5/6_missions): the consolidation will be performed by re-computing vignettes object_bru_area using a custom ImageJ macro
# Attention: The pyImageJ macro is great but it requires a lot of RAM. Do not use unless you're the only one working on the server!
#d- Ser run_macro=False. Metadata are missing and the script is run on the Complex server Marie (home/plankton/uvp5/6_missions): the consolidation will be performed by matching-up the vignettes to the native project files id (image index, blob index) and extracting the corresponding area saved in the bru files
# Attention: This is a workaround to pyimageJ that requires less memory since it is based on python code only. However, you still need an access to the server and the native UVP project folders should be clean (no compressed sub-directories, no profiles renamed with _old versions saved)!
# For disorganized UVP native projects, non-matching vignettes will be assign a size equivalent to object_area

print('Consolidating UVP particle sizes (Ecotaxa: Large particles + Ecopart: Small particles), please wait\n')
standardizer=[standardizer for standardizer in standardizer_files if 'UVP' in str(standardizer.stem)][0]
df_standardizer = pd.read_excel(standardizer, index_col=0,sheet_name='ecotaxa')
with tqdm(desc='{} dataset'.format(standardizer.name.split('_')[1]), total=len(list(df_standardizer.index)), bar_format='{desc}{bar}', position=0, leave=True) as bar:
    for project in list(df_standardizer.index):
        df_standardizer = pd.read_excel(standardizer, index_col=0,sheet_name='ecotaxa') # Standardizer should be updated for each loop iteration
        percent = np.round(100 * (bar.n / len(list(df_standardizer.index))), 1)
        bar.set_description('{} dataset ID {} (%s%%)'.format(standardizer.name.split('_')[1], project) % percent, refresh=True)

        if len(list(Path(df_standardizer.at[project,'Project_localpath']).expanduser().glob('ecotaxa_export_{}_*.tsv'.format(str(project))))):
            continue
        else:

            try:
                consolidate_ecotaxa_project(project_id=project, standardizer=df_standardizer, sheetname='ecotaxa',run_matchup=False, run_macro=False, upload_metadata=False,localpath_metadata=Path(cfg['raw_dir']).expanduser() / cfg['UVP_consolidation_subdir'])

            except Exception as e:
                print('\nSkipping consolidation of project ', str(project), '\n', e, sep='')
        ok = bar.update(n=1)

#3) Standardize project export files using the standardization_func function in funcs_standardize_projects.py
print('Performing project standardization')
for standardizer in natsorted(standardizer_files)[::-1]:
    df_standardizer = pd.read_excel(standardizer, index_col=0)
    with tqdm(desc='{} dataset'.format(standardizer.name.split('_')[1]), total=len(list(df_standardizer.index)), bar_format='{desc}{bar}', position=0, leave=True) as bar:
        for project in list(df_standardizer.index):
            percent = np.round(100 * (bar.n / len(list(df_standardizer.index))), 1)
            bar.set_description('{} dataset ID {} (%s%%)'.format(standardizer.name.split('_')[1], project) % percent, refresh=True)

            try:
                standardization_func(standardizer_path=standardizer, project_id=project,plot='nbss')

            except Exception as e:
                print('\nSkipping standardization of project ', str(project),'\n',e, sep='')
            ok = bar.update(n=1)

#4) Flag samples and generate a report using the quality_control_func function in funcs_standardize_projects.py
print('Performing project control quality check based on the following criteria, please wait:\nFlag_missing: Missing data/metadata\nFlag_GPScoordinatesonland: GPS coordinates on land\nFlag_dubiousGPScoordinates: Dubious GPS coordinates\nFlag_count: Low ROI counts (yielding uncertainties>5%)\nFlag_validation: Low percentage of taxonomic annotations validation (<95%). Only for Zooscan and UVP\nFlag_artefacts: High percentage of artefacts (>20%)\nFlag_size: Multiple size calibration factors\n(1:flagged, 0:no flag)\n')
for standardizer in natsorted(standardizer_files)[::-1]:
    df_standardizer = pd.read_excel(standardizer, index_col=0)
    with tqdm(desc='{} dataset'.format(standardizer.name.split('_')[1]), total=len(list(df_standardizer.index)), bar_format='{desc}{bar}', position=0, leave=True) as bar:
        for project in list(df_standardizer.index):
            percent = np.round(100 * (bar.n / len(list(df_standardizer.index))), 1)
            bar.set_description('{} dataset ID {} (%s%%)'.format(standardizer.name.split('_')[1], project) % percent, refresh=True)

            # Flagging project
            report_file='report_project_'+str(project)+'.html'
            report_path=path_to_standardizer.parent / cfg['report_subdir']
            try:
                quality_control_func(standardizer_path=standardizer, project_id=project,report_path=report_path,validation_threshold=0.95)
            except Exception as e:
                print('\nSkipping flagging of project ',str(project),'\n',e,sep='')
            ok = bar.update(n=1)

print('\nStandardization and flagging processes finished.\n Please check project standardization (under figures/standardizer) to verify column names and units were correctly specified in the corresponding standardizer spreadsheet.\n Please send the flagging report (under reports) and flag tables (optional, under raw/flags) to individual PI before releasing the products. This is to ensure that samples have been correctly located and flagged, or alternatively, that certain flags should be overruled (use the overrule column in the flag table).')