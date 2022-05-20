## Objective: Automate project search and listing using Ecotaxa API by generating a project_list_all and instrument-specific project files

## Requirements:
# Note that you need to install ecotaxa API using git first:
# in terminal: pip install git+https://github.com/ecotaxa/ecotaxa_py_client.git
# Note that you need to authenticate on EcoTaxa to search projects based on title/instrument filter, regardless of access rights
# Note that you need two text files in your git repository that contain the info relative to:
# authentication (Ecotaxa_API_pw.yaml, in .gitignore), output path (Ecotaxa_API.yaml, tracked)
# Note that you need to change the path of the config files on lines 58 & 62

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

# Dictionary/interpreter modules (used to create instrument-specific directory of variables of interest)
from collections import ChainMap
import difflib
from funcy import join_with # use pip install funcyfrom Scripts.funcsdictmatch import dict_match
from scripts.funcs_dictmatch import dict_match
import re

# Ecotaxa API modules. Install module via git first
import ecotaxa_py_client
from ecotaxa_py_client.api import projects_api
from ecotaxa_py_client.model.project_model import ProjectModel
from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.model.login_req import LoginReq

# Workflow starts here
# read git-tracked config file (text file) with inputs:  output directory
path_to_config=Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
# read config file (text file) with inputs:  EcoTaxa username, password. Make sure to secure this file if needed
path_to_config_usr=Path('~/GIT/PSSdb/scripts/Ecotaxa_API_pw.yaml').expanduser()
with open(path_to_config_usr ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

path_to_git=Path(cfg['git_dir']).expanduser()

# prepare storage based on path stored in the yaml config file
path_to_data = path_to_git / Path(cfg['dataset_subdir']).parent

# Step 1: Search visible and accessible projects based on instrument filter
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
        api_response_project_search_accessible = api_instance.search_projects(also_others=False,# Return accessible projects given account rights
                                                                           title_filter='',# Optional title project filter
                                                                           instrument_filter='',  # All instruments
                                                                           order_field='projid'# Sorting variable. Use instrument or projid
                                                                           )

    except ecotaxa_py_client.ApiException as e:
        print("Exception when calling ProjectsApi->search_projects: %s\n" % e)


print("Searching for projects on EcoTaxa:",len(api_response_project_search_visible)+len(api_response_project_search_accessible),"projects found.",len(api_response_project_search_accessible),"projects accessible", sep=' ')

# Step 2: Generate a dataframe with variables of interest (project ID, instrument, contact infos, access rights, and image annotation statistics) and save/overwrite in project_list_all.xslx
df=pd.concat([pd.DataFrame({'Project_ID':list(map(lambda x: x.projid,api_response_project_search_accessible)),
'Instrument':list(map(lambda x: x.instrument,api_response_project_search_accessible)),
'Contact_name':list(map(lambda x: x.managers[0].name if (x.contact is None) else x.contact['name'],api_response_project_search_accessible)),
'Contact_email':list(map(lambda x: x.managers[0].email if (x.contact is None) else x.contact['email'],api_response_project_search_accessible)),
'PSSdb_access':list(map(lambda x: str(cfg_pw['ecotaxa_user'] in list(map(lambda y: (y['email']),x.annotators))),api_response_project_search_accessible)),
'Objects_number':list(map(lambda x: x.objcount,api_response_project_search_accessible)),
'Percentage_classified':list(map(lambda x: x.pctclassified,api_response_project_search_accessible)),
'Percentage_validated':list(map(lambda x: x.pctvalidated, api_response_project_search_accessible))}),
pd.DataFrame({'Project_ID':list(map(lambda x: x.projid,api_response_project_search_visible)),
'Instrument':list(map(lambda x: x.instrument,api_response_project_search_visible)),
'Contact_name':list(map(lambda x: x.managers[0].name if (x.contact is None) else x.contact['name'],api_response_project_search_visible)),
'Contact_email':list(map(lambda x: x.managers[0].email if (x.contact is None) else x.contact['email'],api_response_project_search_visible)),
'PSSdb_access':list(map(lambda x: str(cfg_pw['ecotaxa_user'] in list(map(lambda y: (y['email']),x.annotators))),api_response_project_search_visible)),
'Objects_number':list(map(lambda x: x.objcount,api_response_project_search_visible)),
'Percentage_classified':list(map(lambda x: x.pctclassified,api_response_project_search_visible)),
'Percentage_validated':list(map(lambda x: x.pctvalidated, api_response_project_search_visible))})
])

df_metadata=pd.DataFrame({'Variables':df.columns,'Variable_types':df.dtypes,
'Units/Values':['','','','','','#','%','%'],
'Description':['Project ID in EcoTaxa','Project instrument','Name of the project contact','Email of the project contact','Project accessibility. If True, export is possible with current authentication','Total number of objects (images) in the project','Percentage of predicted images','Percentage of validated images']})

# Prompting a warning if a project_list_all file already exists
# asking for confirmation to overwrite. Attention: the outdated project list will be erased!

existing_project_path = list(path_to_data.glob('project_list_all*'))

if len(existing_project_path) != 0:

    confirmation = input(
        "Projects list already created. Do you wish to overwrite the existing list? Enter Y or N\n")
    if confirmation == 'Y':
        print("Overwriting project_list_all file, please wait")
        with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'), engine="openpyxl", mode="a",
                            if_sheet_exists="replace") as writer:
            df.to_excel(writer, sheet_name='Data', index=False)



else:
    print("Creating project_list_all file, please wait")
    # Generating project_list_all file
    with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'),engine="xlsxwriter") as writer:
        df.to_excel(writer, sheet_name='Data', index=False)
        df_metadata.to_excel(writer, sheet_name='Metadata', index=False)
# Step 3: Generate an instrument-specific overview of accessible projects (project ID, Entry/latest update date, sampling ranges, image annotation statistics, and EcoPart ID) and save/overwrite in instrument_project_list.xslx
# Step 3a: Build instrument-specific interpreters(dictionaries) to search specific fields in export files/EcoTaxa database
dict_instruments={'IFCB':['IFCB'],'UVP':['UVP5HD','UVP5SD','UVP5Z','UVP6'],'Zooscan':['Zooscan'],'Unknown':['?'],'AMNIS':['AMNIS'],'CPICS':['CPICS'],'CytoSense':['CytoSense'],'FastCam':['FastCam'],'FlowCam':['FlowCam'],'ISIIS':['ISIIS'],'LISST':['LISST','LISST-Holo'],'Loki':['Loki'],'Other':['Other camera','Other flowcytometer','Other microscope','Other scanner'],'PlanktoScope':['PlanktoScope'],'VPR':['VPR'],'ZooCam':['ZooCam'],'eHFCM':['eHFCM']}
for i, n in enumerate(api_response_project_search_accessible):
   api_response_project_search_accessible[i]['instrument']=''.join([key for key,value in dict_instruments.items() if api_response_project_search_accessible[i]['instrument'] in value])

inst=list(set(list(map(lambda x: x.instrument, api_response_project_search_accessible))))
inst_all=list(map(lambda x: x.instrument, api_response_project_search_accessible))

print("Building standardizer for",'/'.join(inst),"projects",sep=" ")
# Step 3b: Use instrument-specific interpreters to provide summary/custom fields of interest and store in instrument_project_list files
subset_df=df[df['PSSdb_access']=='True']
#inst=[key for key, value in dict_instruments.items() if any(pd.Series(value).isin(subset_df.Instrument.unique().tolist()))]#subset_df.Instrument.unique()
for instrument in inst:
    if len(list(path_to_data.glob('project_{}_standardizer*'.format(str(instrument)))))==0:
        subset_df_instrument = pd.DataFrame(
            {'Project_ID': subset_df[subset_df['Instrument'].isin(dict_instruments[instrument])]['Project_ID'],
             'Instrument': subset_df[subset_df['Instrument'].isin(dict_instruments[instrument])]['Instrument']})

        df_sub = pd.concat(list(map(lambda y: pd.DataFrame.from_dict({
            'Cruise.field': ([""] + y)[
                dict_match("program|cruise", [re.subn("object_|acq_|process_|sample_", "", string)[0] for string in y])[
                    0]['Index'] + 1],
            'Station.field': ([""] + y)[
                dict_match("station", [re.subn("object_|acq_|process_|sample_", "", string)[0] for string in y])[0][
                    'Index'] + 1],
            'Profile.field': ([""] + y)[
                dict_match("profile|cast", [re.subn("object_|acq_|process_|sample_", "", string)[0] for string in y])[
                    0]['Index'] + 1],
            'Latitude.field': 'object_lat',
            'Latitude.unit': '',
            'Longitude.field': 'object_lon',
            'Longitude.unit': '',
            'Depth_min.field': 'object_depth_min',
            'Depth_min.unit': '',
            'Depth_max.field': 'object_depth_max',
            'Depth_max.unit': '',
            'Volume_analyzed.field': ([""] + y)[
                dict_match("vol", [re.subn("object_|acq_|process_|sample_", "", string)[0] for string in y])[0][
                    'Index'] + 1],
            'Volume_analyzed.unit': '',
            'Dilution.field': ([""] + y)[
                dict_match("dil", [re.subn("object_|acq_|process_|sample_", "", string)[0] for string in y])[0][
                    'Index'] + 1],
            'ESD.field': ([""] + y)[
                dict_match("esd", [re.subn("object_|acq_|process_|sample_", "", string)[0] for string in y])[0][
                    'Index'] + 1],
            'ESD.unit': '',
            'Biovolume.field': ([""] + y)[
                dict_match("biovolume", [re.subn("object_|acq_|process_|sample_", "", string)[0] for string in y])[0][
                    'Index'] + 1],
            'Biovolume.unit': '',
            'Area.field': ([""] + y)[
                dict_match("area", [re.subn("object_|acq_|process_|sample_", "", string)[0] for string in y])[0][
                    'Index'] + 1],
            'Area.unit': '',
            'Pixel.field': ([""] + y)[
                dict_match("pixel", [re.subn("object_|acq_|process_|sample_", "", string)[0] for string in y])[0][
                    'Index'] + 1],
            'Pixel.unit': '',
            'Category.field': ''
        }, orient="index").T,
                                    list(
                                        map(lambda x: ["acq_" + str for str in list(x.acquisition_free_cols.keys())] + [
                                            "object_" + str for str in list(x.obj_free_cols.keys())] + ["sample_" + str
                                                                                                        for str in list(
                                                x.sample_free_cols.keys())] + ["process_" + str for str in
                                                                               list(x.process_free_cols.keys())], list(
                                            filter(lambda w: any(subset_df_instrument['Project_ID'].isin([w.projid])),
                                                   api_response_project_search_accessible))))
                                    )), ignore_index=True)
        df_sub.index = subset_df_instrument.index.tolist()
        subset_df_instrument = pd.concat([subset_df_instrument, df_sub], axis=1)
        subset_df_instrument_metadata = pd.DataFrame(
            {'Variables': subset_df_instrument.columns, 'Variable_types': subset_df_instrument.dtypes,
             'Description': ['Project ID in EcoTaxa',
                             'Project instrument',
                             'Name of the cruise ID column in project export file',
                             'Name of the station ID column in project export file',
                             'Name of the profile ID column in project export file',
                             'Name of the latitude column in project export file',
                             'Unit of the latitude column in project export file',
                             'Name of the longitude column in project export file',
                             'Unit of the longitude column in project export file',
                             'Name of the minimum depth column in project export file',
                             'Unit of the minimum depth column in project export file',
                             'Name of the maximum depth column in project export file',
                             'Unit of the maximum depth column in project export file',
                             'Name of the sample volume analyzed column in project export file',
                             'Unit of the sample volume analyzed column in project export file',
                             'Name of the dilution factor column in project export file',
                             'Name of the equivalent spherical diameter column in project export file',
                             'Unit of the equivalent spherical diameter column in project export file',
                             'Name of the biovolume column in project export file',
                             'Unit of the biovolume column in project export file',
                             'Name of the area column in project export file',
                             'Unit of the area column in project export file',
                             'Name of the pixel conversion factor column in project export file',
                             'Unit of the pixel conversion factor column in project export file',
                             'Name of the category ID column in project export file']})
        print("Creating project_", instrument, "_standardizer file, please wait", sep="")
        # Generating project_list_all file
        with pd.ExcelWriter(str(path_to_data / 'project_{}_standardizer.xlsx'.format(str(instrument))),
                            engine="xlsxwriter") as writer:
            subset_df_instrument.to_excel(writer, sheet_name='Data', index=False)
            subset_df_instrument_metadata.to_excel(writer, sheet_name='Metadata', index=False)