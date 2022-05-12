## Objective: Automate project search and listing using Ecotaxa API by generating a project_list_all and instrument-specific project files

## Requirements:
# Note that you need to install ecotaxa API using git first:
# in terminal: pip install git+https://github.com/ecotaxa/ecotaxa_py_client.git
# Note that you need to authenticate on EcoTaxa to search projects based on title/instrument filter, regardless of access rights
# Note that you need two text files in your git repository that contain the info relative to:
# authentication (Ecotaxa_API_pw.yaml, in .gitignore), output path (Ecotaxa_API.yaml, tracked)
# Note that you need to change the path of the config files on lines 51 & 55

## Workflow:
# This script uses 4 API functions:
# Entry 1: Search projects ID given instrument filter.
# Entry 2: Generate an exhaustive list of accessible/inaccessible project(s) based on your Ecotaxa login/pw (info stored in untracked Ecotaxa_API_pw.yaml). Projects list stored as project_list_all
# Entry 3: Query projects/objects summary to inform project uploading/updating/sampling date, spatial coverage, images predictions/annotations status. Create instrument-specific project lists
# Entry 4: Bridge EcoPart-EcoTaxa projects (UVP-specific)

## Useful documentations:
# See documentation for entry 1: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectsApi.md#search_projects
# See documentation for entry 2: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectsApi.md#search_projects
# See documentation for entry 3: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectsApi.md#project_set_get_column_stats,https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ObjectsApi.md#get_object_set_summary
# See documentation for entry 4: upcoming

## TO DO: Automate search for person to contact, Bridge EcoPart and EcoTaxa projects


## Python modules

# Path modules
from pathlib import Path # Handling of path object

# Config modules
import yaml # requires installation of PyYAML package

# Prompt for export confirmation
import sys
import time

# Panda module (used for excel files)
import pandas as pd

# Ecotaxa API modules. Install module via git first
import ecotaxa_py_client
from ecotaxa_py_client.api import projects_api
from ecotaxa_py_client.model.project_model import ProjectModel
from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.model.login_req import LoginReq

# Workflow starts here
# read git-tracked config file (text file) with inputs:  output directory
path_to_config=Path('~/GIT/PSSdb/Scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
# read config file (text file) with inputs:  EcoTaxa username, password. Make sure to secure this file if needed
path_to_config_usr=Path('~/GIT/PSSdb/Scripts/Ecotaxa_API_pw.yaml').expanduser()
with open(path_to_config_usr ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

path_to_git=Path(cfg['git_dir']).expanduser()

# prepare storage based on path stored in the yaml config file
path_to_data = path_to_git / Path(cfg['dataset_subdir']).parent

# Step 1: Search all visible projects based on instrument filter
# This step does require your authentication in EcoTaxa
with ecotaxa_py_client.ApiClient() as client:
    api = authentification_api.AuthentificationApi(client)
    token = api.login(LoginReq(
                               username=cfg_pw['ecotaxa_user'],
                               password=cfg_pw['ecotaxa_pass']
                               ))

configuration = ecotaxa_py_client.Configuration(host = "https://ecotaxa.obs-vlfr.fr/api",access_token=token, discard_unknown_keys=True)

with ecotaxa_py_client.ApiClient(configuration) as api_client:
    api_instance = projects_api.ProjectsApi(api_client)
    try:
        # Search Projects
        api_response_project_search_visible = api_instance.search_projects(also_others=True, # Return visible projects with no access rights
                                                                   title_filter='', # Optional title project filter
                                                                   instrument_filter='', # All instruments
                                                                   order_field='projid' # Sorting variable. Use instrument or projid
                                                                   )
        api_response_project_search_accessible = api_instance.search_projects(also_others=False,# Return visible projects with no access rights
                                                                           title_filter='',# Optional title project filter
                                                                           instrument_filter='',  # All instruments
                                                                           order_field='projid'# Sorting variable. Use instrument or projid
                                                                           )

    except ecotaxa_py_client.ApiException as e:
        print("Exception when calling ProjectsApi->search_projects: %s\n" % e)


print("Searching for projects:",len(api_response_project_search_visible)+len(api_response_project_search_accessible),"projects found", sep=' ')

# Generate a dataframe with variables of interest
df=pd.concat([pd.DataFrame({'Project_ID':list(map(lambda x: x.projid,api_response_project_search_accessible)),
'Instrument':list(map(lambda x: x.instrument,api_response_project_search_accessible)),
'Contact_name':list(map(lambda x: x.managers[0].name if (x.contact is None) else x.contact['name'],api_response_project_search_accessible)),
'Contact_email':list(map(lambda x: x.managers[0].email if (x.contact is None) else x.contact['email'],api_response_project_search_accessible)),
'PSSdb_access':list(map(lambda x: str(cfg_pw['ecotaxa_user'] in list(map(lambda y: (y['email']),x.annotators))),api_response_project_search_accessible)),
'Percentage_classified':list(map(lambda x: x.pctclassified,api_response_project_search_accessible)),
'Percentage_validated':list(map(lambda x: x.pctvalidated, api_response_project_search_accessible))}),
pd.DataFrame({'Project_ID':list(map(lambda x: x.projid,api_response_project_search_visible)),
'Instrument':list(map(lambda x: x.instrument,api_response_project_search_visible)),
'Contact_name':list(map(lambda x: x.managers[0].name if (x.contact is None) else x.contact['name'],api_response_project_search_visible)),
'Contact_email':list(map(lambda x: x.managers[0].email if (x.contact is None) else x.contact['email'],api_response_project_search_visible)),
'PSSdb_access':list(map(lambda x: str(cfg_pw['ecotaxa_user'] in list(map(lambda y: (y['email']),x.annotators))),api_response_project_search_visible)),
'Percentage_classified':list(map(lambda x: x.pctclassified,api_response_project_search_visible)),
'Percentage_validated':list(map(lambda x: x.pctvalidated, api_response_project_search_visible))})
])

df_metadata=pd.DataFrame({'Variables':df.columns,'Variable_types':df.dtypes,
'Units/Values':['','','','','','%','%'],
'Description':['Project ID in EcoTaxa','Project instrument','Name of the project contact','Email of the project contact','Project accessibility. If True, export is possible','Percentage of predicted images','Percentage of validated images']})

# Prompting a warning if a project_list_all file already exists
# asking for confirmation to overwrite. Attention: the outdated project list will be erased!

existing_project_path = list(path_to_data.glob('project_list_all*'))

if len(existing_project_path) != 0:

    confirmation = input(
        "Projects list already created. Do you wish to overwrite the existing list? Enter Y or N\n")
    if confirmation != 'Y':
        quit()

    print("Overwriting project_list_all file, please wait")
    with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'), engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        df.to_excel(writer, sheet_name='Data', index=False)

else:
    print("Creating project_list_all file, please wait")
    # Generating project_list_all file
    with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'),engine="xlsxwriter") as writer:
        df.to_excel(writer, sheet_name='Data', index=False)
        df_metadata.to_excel(writer, sheet_name='Metadata', index=False)
