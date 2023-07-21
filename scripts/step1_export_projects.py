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
path_to_config=Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_to_git=Path(cfg['git_dir']).expanduser()
path_to_data=path_to_git / cfg['dataset_subdir']

path_to_config_usr=Path('~/GIT/PSSdb/scripts/configuration_masterfile_pw.yaml').expanduser()
with open(path_to_config_usr ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

# Progress bar modules
#import progressbar # Attention: Use pip install progressbar2 to install
from tqdm import tqdm
import datetime, time # Time module for break and current date

# Prompt for export confirmation
import sys

# Dataframe/Array modules
import pandas as pd
import numpy as np
import os

# Ecotaxa API modules. Install module via git first
import ecotaxa_py_client

# Data portal-specific web scraping functions to export projects
from sys import argv
import requests
from funcs_export_projects import *

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
    project_ids=np.array(project_list[(project_list['PSSdb_access']==True) & (project_list['Project_test']==True)].Project_ID.astype(str))

# prompting a warning if export files for the list of projects already exist
# asking for confirmation to overwrite. Attention: all outdated project export files will be erased!
existing_project_path=list(path_to_projects.rglob('*export_*'))
existing_project_ids=list(map(lambda x: str(x.stem)[(str(x.stem).find('export_')+7):([i.start() for i in re.finditer('_',str(x.stem))])[2]],existing_project_path))
if len(existing_project_path)!=0:

  confirmation=input("{} project(s) already downloaded. Do you wish to overwrite the export file(s)? Enter Y or N\n If Y, outdated project export(s) will be erased\n".format(str(len(existing_project_ids))))
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
    project_ids=np.array(project_list[(project_list['PSSdb_access']==True) & (project_list['Project_test']==True)].Project_ID.astype(str))

# prompting a warning if export files for the list of projects already exist
# asking for confirmation to overwrite. Attention: all outdated project export files will be erased!
existing_project_path=list(path_to_projects.rglob('*export_*'))
existing_project_ids=list(set(list(map(lambda x: str(x.stem)[1+([i.start() for i in re.finditer('_',str(x.stem))])[2]:([i.start() for i in re.finditer('_',str(x.stem))])[3]],existing_project_path))))
if len(existing_project_path)!=0:

  confirmation=input("{} project(s) already downloaded. Do you wish to overwrite the export file(s)? Enter Y or N\n If Y, outdated project export(s) will be erased\n".format(str(len(existing_project_ids))))
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

# Loop through project using Ecopart export function (see funcs_export_projects.py)
for i in range(len(project_ids)):
    proj_id = project_ids[i]
    Ecopart_export(project={proj_id:project_list.loc[project_list.Project_ID==int(proj_id),'Project_title'].values[0]},username=cfg_pw['ecotaxa_user'], password=cfg_pw['ecotaxa_pass'],localpath=path_to_projects)

# 3) Exporting IFCB dashboard projects:
print("Exporting projects hosted on CALOOS (https://ifcb.caloos.org/) and WHOI (https://ifcb-data.whoi.edu/) IFCB dashboards:")
project_list=pd.read_excel(path_to_data.parent / cfg['proj_list'],sheet_name="ifcb",usecols=[ 'dashboard_url', 'Project_source', 'Project_ID','Instrument','PSSdb_access','Project_test','Project_localpath', 'Latest_update'])
project_ids = np.array(project_list[project_list['PSSdb_access']==True].Project_ID.astype(str))#str(cfg['proj_id'])
path_to_projects=Path(project_list.at[0,'Project_localpath']).expanduser()

# query pattern to restrict file download is different due to the structure of IFCB dashboard datasets,

testing = input (' Do you want to download all projects or just the test? \n Enter all or test ')
if testing == 'test':
    timeseries_data = project_list[project_list['Project_test'] == True].reset_index(drop=True)
elif testing == 'all':
    timeseries_data = project_list


subset = input (' Do you want to get data for all the time series or a date range or the test dates? \n Enter all, range or test ')
if subset == 'range':
    start_date= input(' enter starting date of the data as YYYYMMDD')
    start_date= int(start_date)
    end_date = input(' enter final date of the data as YYYYMMDD')
    end_date = int(end_date)
elif subset == 'all':
    start_date = 20000101
    end_date = 21000101
elif subset == 'test':
    start_date = 20180825
    end_date = 20180825

dashboard = input (' Do you want to get data from WHOI or CALOOS dashboards? \n Enter WHOI or CALOOS ')
if dashboard == 'CALOOS':
    timeseries_data = timeseries_data.loc[timeseries_data['dashboard_url'] == 'https://ifcb.caloos.org/'].reset_index(drop=True)
    timeseries_data
if dashboard == 'WHOI':
    timeseries_data = timeseries_data.loc[timeseries_data['dashboard_url'] == 'https://ifcb-data.whoi.edu/'].reset_index(drop=True)

specify_proj = input (' Would you like to automatically download all data or data of just one project? \n Enter ALL or ONE ')

update_project = input ('Are the full dataset being downloaded or are the projects being updated with the latest data? \n Enter FULL or UPDATE')

if specify_proj == 'ONE':
    Project_ID= input(' Would you like to automatically download all data or data of just one project? \n Enter project ID ')
    Project_source = timeseries_data.loc[timeseries_data['Project_ID'] == Project_ID, 'Project_source'].iloc[0]
    Project_localpath = timeseries_data.loc[timeseries_data['Project_ID'] == Project_ID, 'Project_localpath'].iloc[0]
    dashboard_url = timeseries_data.loc[timeseries_data['Project_ID'] == Project_ID, 'dashboard_url'].iloc[0]
    path_download = str(path_to_projects) + '/' + Project_ID  ## NOTE: the first part of the directory needs to be changed for the GIT PSSdb project
    # pathname = Project_source + Project_ID + '/'
    if update_project =='FULL':
        if os.path.isdir(path_download) and len(os.listdir(path_download)) != 0:  # and  os.path.exists(path_download)
            replace = input('There is already downloaded data in ' + path_download + ' do you want to replace the files? \n Y/N')
            if replace == 'Y':
                print('Overwriting ' + Project_ID + ' file(s), please wait')
                shutil.rmtree(path_download)
                os.mkdir(path_download)
            elif replace == 'N':
                print('Skipping ' + Project_ID)
        elif not os.path.exists(path_download):
            os.mkdir(path_download)
        IFCB_dashboard_export(dashboard_url, Project_source, Project_ID, path_download, start_date, end_date,update_project)
    else:
        IFCB_dashboard_export(dashboard_url, Project_source, Project_ID, path_download, start_date, end_date, update_project)

else:
    for n in range (0, len(timeseries_data)):
        Project_source = timeseries_data.loc[n, 'Project_source']
        Project_ID = timeseries_data.loc[n, 'Project_ID']
        Project_localpath = timeseries_data.loc[n, 'Project_localpath']
        dashboard_url = timeseries_data.loc[n, 'dashboard_url']
        # define path for download
        path_download = str(path_to_projects) + '/' + Project_ID  ## NOTE: the first part of the directory needs to be changed for the GIT PSSdb project
        #pathname = Project_source + Project_ID + '/'
        if update_project =='FULL':
            if os.path.isdir(path_download) and len(os.listdir(path_download) ) != 0: # and  os.path.exists(path_download)
                replace = input('There is already downloaded data in ' +path_download+ ' do you want to replace the files? \n Y/N')
                if replace == 'Y':
                    print('Overwriting '+ Project_ID  +' file(s), please wait')
                    shutil.rmtree(path_download)
                    os.mkdir(path_download)
                elif replace == 'N':
                    print('Skipping ' + Project_ID)
                    continue
            elif not os.path.exists(path_download):
                os.mkdir(path_download)
            IFCB_dashboard_export(dashboard_url, Project_source, Project_ID, path_download, start_date, end_date, update_project)
        else:
            IFCB_dashboard_export(dashboard_url, Project_source, Project_ID, path_download, start_date, end_date,update_project)





quit()