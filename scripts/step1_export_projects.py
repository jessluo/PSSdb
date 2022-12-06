## Objective: Exporting projects listed in project_list_all spreadsheet based on portal-specific functions

## TO DO: Export projects from additional data sources (e.g. IFCB dashboards)


## Python modules

# Path modules
from pathlib import Path # Handling of path object
import shutil # Delete uncompressed export zip folder
import re

# Config modules
import yaml # requires installation of PyYAML package
# read git-tracked config file (text file) with inputs:  project ID, output directory
path_to_config=Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_to_git=Path(cfg['git_dir']).expanduser()
path_to_data=path_to_git / cfg['dataset_subdir']

path_to_config_usr=Path('~/GIT/PSSdb/scripts/Ecotaxa_API_pw.yaml').expanduser()
with open(path_to_config_usr ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

# Progress bar modules
#import progressbar # Attention: Use pip install progressbar2 to install
from tqdm import tqdm
import datetime, time # Time module for break and current date

# Prompt for export confirmation
import sys

# Ecotaxa API modules. Install module via git first
import ecotaxa_py_client
from funcs_export_projects import *
from sys import argv
import requests

# Dataframe modules
import pandas as pd
import numpy as np

# Workflow starts here

# 1) Exporting Ecotaxa projects:
print('Exporting projects hosted on Ecotaxa (https://ecotaxa.obs-vlfr.fr/):')
# prepare storage based on project list stored in the yaml config file and instrument type
project_list=pd.read_excel(path_to_data.parent / cfg['proj_list'],sheet_name="ecotaxa",usecols=['Project_ID','Instrument','PSSdb_access','Project_test','Project_localpath'])
project_ids = np.array(project_list[project_list['PSSdb_access']==True].Project_ID.astype(str))#str(cfg['proj_id'])
path_to_projects=Path(project_list.at[0,'Project_localpath']).expanduser()

dict_instruments={'IFCB':'IFCB','UVP':['UVP5HD','UVP5SD','UVP5Z','UVP6'],'Zooscan':'Zooscan','Unknown':'?','AMNIS':'AMNIS','CPICS':'CPICS','CytoSense':'CytoSense','FastCam':'FastCam','FlowCam':'FlowCam','ISIIS':'ISIIS','LISST':['LISST','LISST-Holo'],'Loki':'Loki','Other':['Other camera','Other flowcytometer','Other microscope','Other scanner'],'PlanktoScope':'PlanktoScope','VPR':'VPR','ZooCam':'ZooCam','eHFCM':'eHFCM'}
project_inst=dict(zip(project_ids,list(map(lambda inst: [key for key,value in dict_instruments.items() if inst in value]  ,project_list[project_list['PSSdb_access']==True].Instrument.astype(str)))))

# prompting a warning to export all accessible projects or test set only
test_confirmation=input("Do you wish to export all accessible project(s) or the test subset? Enter all or test\n")
if test_confirmation=='test':
    project_ids=project_list[(project_list['PSSdb_access']==True) & (project_list['Project_test']==True)].Project_ID.astype(str)

# prompting a warning if export files for the list of projects already exist
# asking for confirmation to overwrite. Attention: all outdated project export files will be erased!
existing_project_path=list(path_to_projects.rglob('*export_*'))
existing_project_ids=list(map(lambda x: str(x.stem)[(str(x.stem).find('export_')+7):([i.start() for i in re.finditer('_',str(x.stem))])[2]],existing_project_path))
if len(existing_project_path)!=0:

  confirmation=input("{} project(s) already downloaded. Do you wish to overwrite the export file(s)? Enter Y or N\n If Y, outdated project export(s) will be erased\n".format(str(len(existing_project_path))))
  if confirmation!='Y':
      project_ids=np.array([x for x in project_ids if x not in existing_project_ids])
      print("Skipping existing project(s). Creating export file for other projects")

  if confirmation == 'Y':
     print("Overwriting project export file(s), please wait")
     for file in existing_project_path:
         shutil.rmtree(file, ignore_errors=True)
         file.unlink(missing_ok=True)
else:
  print("Creating project export file using",str(path_to_data.parent / cfg['proj_list']) ,"file, please wait",sep=" ")

# Loop through project using Ecotaxa export function (see funcs_export_projects.py)
for i in range(len(project_ids)):
    proj_id = int(project_ids[i])
    Ecotaxa_export(project=proj_id,username=cfg_pw['ecotaxa_user'], password=cfg_pw['ecotaxa_pass'],localpath=path_to_projects / project_inst[project_ids[i]][0])

# 2) Exporting Ecopart projects:
print("Exporting projects hosted on EcoPart (https://ecopart.obs-vlfr.fr/):")
# prepare storage based on project list stored in the yaml config file and instrument type
project_list=pd.read_excel(path_to_data.parent / cfg['proj_list'],sheet_name="ecopart",usecols=['Project_ID','Project_title','Project_localpath','Instrument','PSSdb_access','Project_test'])
project_ids = np.array(project_list[project_list['PSSdb_access']==True].Project_ID.astype(str))#str(cfg['proj_id'])
path_to_projects=Path(project_list.at[0,'Project_localpath']).expanduser()

# prompting a warning to export all accessible projects or test set only
test_confirmation=input("Do you wish to export all accessible project(s) or the test subset? Enter all or test\n")
if test_confirmation=='test':
    project_ids=project_list[(project_list['PSSdb_access']==True) & (project_list['Project_test']==True)].Project_ID.astype(str)

# prompting a warning if export files for the list of projects already exist
# asking for confirmation to overwrite. Attention: all outdated project export files will be erased!
existing_project_path=list(path_to_projects.rglob('*export_*'))
existing_project_ids=list(set(list(map(lambda x: str(x.stem)[1+([i.start() for i in re.finditer('_',str(x.stem))])[2]:([i.start() for i in re.finditer('_',str(x.stem))])[3]],existing_project_path))))
if len(existing_project_path)!=0:

  confirmation=input("{} project(s) already downloaded. Do you wish to overwrite the export file(s)? Enter Y or N\n If Y, outdated project export(s) will be erased\n".format(str(len(existing_project_path))))
  if confirmation!='Y':
      project_ids=np.array([x for x in project_ids if x not in existing_project_ids])
      print("Skipping existing project(s). Creating export file for other projects")

  if confirmation == 'Y':
     print("Overwriting project export file(s), please wait")
     for file in existing_project_path:
         shutil.rmtree(file, ignore_errors=True)
         file.unlink(missing_ok=True)
else:
  print("Creating project export file using",str(path_to_data / cfg['proj_list']) ,"file, please wait",sep=" ")

# Loop through project using Ecopart export function (see funcs_export_projects.py)
for i in 0: #range(len(project_ids))
    proj_id = project_ids[i]
    Ecopart_export(project={proj_id:project_list.loc[project_list.Project_ID==proj_id,'Project_title'].values[0]},username=cfg_pw['ecotaxa_user'], password=cfg_pw['ecotaxa_pass'],localpath=path_to_projects)

# 3) Exporting IFCB dashboard projects:

quit()