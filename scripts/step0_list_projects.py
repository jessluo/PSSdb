## Objective: Automate project search and listing using Ecotaxa API by generating a project_list_all and instrument-specific project standardizer

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
# Entry 2: Generate an exhaustive list of visible and accessible project(s) based on your Ecotaxa login/pw (info stored in untracked Ecotaxa_API_pw.yaml). Projects list stored as project_list_all
# Entry 3: Generate instrument-specific project spreadsheets used to standardize and harmonize fields of interest (i.e. needed to calculate NBSS)

## Useful documentations:
# See documentation for entry 1: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectsApi.md#search_projects
# See documentation for entry 2: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectsApi.md#search_projects
# See documentation for entry 3: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectsApi.md#project_set_get_column_stats,https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ObjectsApi.md#get_object_set_summary

## TO DO: List projects downloaded from other sources (e.g. IFCB dashboard)


## Python modules

# Path modules
from pathlib import Path # Handling of path object
import os

# Config modules
import yaml # requires installation of PyYAML package

# Prompt for export confirmation
import sys
import time

# Panda module (used for excel files)
import pandas as pd

# Unit convention to fill standardizer spreadsheets
from pint import UnitRegistry # Use pip install pint to install
ureg=UnitRegistry()
ureg.load_definitions(Path.home()/'GIT'/'PSSdb'/'scripts' /'units_def.txt')  # This text file is used to define custom units from standard units  (e.g. square_pixel etc.)
full_list_units = list(dict.fromkeys(sorted(dir(ureg))))  # list(dict.fromkeys(sorted(list(np.concatenate([dir(getattr(ureg.sys,system)) for system in dir(ureg.sys)]).flat))))


# Ecotaxa API modules. Install module via git first
import ecotaxa_py_client
from ecotaxa_py_client.api import projects_api
from ecotaxa_py_client.api import samples_api
from ecotaxa_py_client.model.project_model import ProjectModel
from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.model.login_req import LoginReq
from ecotaxa_py_client.api import objects_api
from ecotaxa_py_client.model.project_filters import ProjectFilters

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
path_to_projects=path_to_git / Path(cfg['dataset_subdir'])
# Step 1: Search visible and accessible projects based on instrument filter
# This step does require your authentication in EcoTaxa
configuration = ecotaxa_py_client.Configuration( host = "https://ecotaxa.obs-vlfr.fr/api")
configuration.verify_ssl=False
with ecotaxa_py_client.ApiClient() as client:
    api = authentification_api.AuthentificationApi(client)
    token = api.login(LoginReq(
                               username=cfg_pw['ecotaxa_user'],
                               password=cfg_pw['ecotaxa_pass']
                               ))

configuration = ecotaxa_py_client.Configuration(host = "https://ecotaxa.obs-vlfr.fr/api",access_token=token, discard_unknown_keys=True)
configuration.verify_ssl=False
print("Searching for projects on EcoTaxa (please wait):")
with ecotaxa_py_client.ApiClient(configuration) as api_client:
    api_instance = projects_api.ProjectsApi(api_client)
    api_object = objects_api.ObjectsApi(api_client)
    try:
        # Search visible projects
        api_response_project_search_visible = api_instance.search_projects(also_others=True, # Return visible projects with no access rights
                                                                   title_filter='', # Optional title project filter
                                                                   instrument_filter='', # All instruments
                                                                   order_field='projid' # Sorting variable. Use instrument or projid
                                                                   )
        # Search accessible projects based on config file authentication
        api_response_project_search_accessible = api_instance.search_projects(also_others=False,# Return accessible projects given account rights
                                                                           title_filter='',# Optional title project filter
                                                                           instrument_filter='',  # All instruments
                                                                           order_field='projid'# Sorting variable. Use instrument or projid
                                                                           )
        # Retrieve timestamp of the most recent annotation for accessible projects (Note this step takes a long time so it is currently limited to accessible projects, but can be applied to visible projects as well)

        project_list_accessible = list(map(lambda x: str(x.projid), api_response_project_search_accessible))
        api_response_project_get_accessible = list(map(lambda x: api_object.get_object_set(project_id=int(x), project_filters=ProjectFilters(statusfilter="PVD"),fields='obj.classif_when').details, project_list_accessible))
        api_response_project_df_accessible = list(map(lambda x:pd.DataFrame(x,columns=['Annotation_timestamp']), api_response_project_get_accessible))
        api_response_project_latestupdate_accessible = list(map(lambda x: max([date for date in x.Annotation_timestamp if date is not None and date is not pd.NaT]) if len([date for date in x.Annotation_timestamp if date is not None and date is not pd.NaT])>0 else pd.NaT,api_response_project_df_accessible))


    except ecotaxa_py_client.ApiException as e:
        print("Exception when calling EcoTaca API function search_projects: %s\n" % e)

print(len(api_response_project_search_visible)+len(api_response_project_search_accessible),"projects found.",len(api_response_project_search_accessible),"projects accessible", sep=' ')

# Step 2: Generate a dataframe with variables of interest (project ID, instrument, contact infos, access rights, and image annotation statistics) and save/overwrite in project_list_all.xslx
df=pd.concat([pd.DataFrame({'Project_ID':list(map(lambda x: x.projid,api_response_project_search_accessible)),
'Project_title':list(map(lambda x: x.title,api_response_project_search_accessible)),
'Instrument':list(map(lambda x: x.instrument,api_response_project_search_accessible)),
'Contact_name':list(map(lambda x: x.managers[0].name if (x.contact is None) else x.contact['name'],api_response_project_search_accessible)),
'Contact_email':list(map(lambda x: x.managers[0].email if (x.contact is None) else x.contact['email'],api_response_project_search_accessible)),
'PSSdb_access':list(map(lambda x: str(cfg_pw['ecotaxa_user'] in list(map(lambda y: (y['email']),x.annotators))),api_response_project_search_accessible)),
'Objects_number':list(map(lambda x: x.objcount,api_response_project_search_accessible)),
'Percentage_classified':list(map(lambda x: x.pctclassified,api_response_project_search_accessible)),
'Percentage_validated':list(map(lambda x: x.pctvalidated, api_response_project_search_accessible)),
'Latest_update':api_response_project_latestupdate_accessible}),
pd.DataFrame({'Project_ID':list(map(lambda x: x.projid,api_response_project_search_visible)),
'Project_title':list(map(lambda x: x.title,api_response_project_search_visible)),
'Instrument':list(map(lambda x: x.instrument,api_response_project_search_visible)),
'Contact_name':list(map(lambda x: x.managers[0].name if (x.contact is None) else x.contact['name'],api_response_project_search_visible)),
'Contact_email':list(map(lambda x: x.managers[0].email if (x.contact is None) else x.contact['email'],api_response_project_search_visible)),
'PSSdb_access':list(map(lambda x: str(cfg_pw['ecotaxa_user'] in list(map(lambda y: (y['email']),x.annotators))),api_response_project_search_visible)),
'Objects_number':list(map(lambda x: x.objcount,api_response_project_search_visible)),
'Percentage_classified':list(map(lambda x: x.pctclassified,api_response_project_search_visible)),
'Percentage_validated':list(map(lambda x: x.pctvalidated, api_response_project_search_visible)),
'Latest_update': pd.NaT})
])
df=df.assign(Project_test=df['Project_ID'].isin([3315,3318,3326,377,378,714,560,579,5693]).astype(str)) # Add boolean variable to identify fixed test set (3 projects per instrument)

df_metadata=pd.DataFrame({'Variables':df.columns,'Variable_types':df.dtypes,
'Units/Values':['','','','','','','#','%','%','',''],
'Description':['Project ID','Project title','Project instrument','Name of the project contact','Email of the project contact','Project accessibility. If True, export is possible with current authentication','Total number of objects (images) in the project','Percentage of predicted images','Percentage of validated images','Test set boolean identifier. If True, project is one of the test set projects','Timestamp of the latest annotation']})

# Prompting a warning if a project_list_all file already exists
# asking for confirmation to overwrite. Attention: the outdated project list will be erased!

existing_project_path = list(path_to_data.glob('project_list_all*'))

if len(existing_project_path) != 0:

    confirmation = input("Projects list already created. Do you wish to overwrite the existing list? Enter Y or N\n")
    if confirmation == 'Y':
        print("Overwriting project_list_all file, please wait")
        with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
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

print("Building/Updating standardizer for",'/'.join(inst),"projects",sep=" ")
# Step 3b: Use instrument-specific interpreters to provide summary/custom fields of interest and store in instrument_project_list files
subset_df=df[df['PSSdb_access']=='True']
#inst=[key for key, value in dict_instruments.items() if any(pd.Series(value).isin(subset_df.Instrument.unique().tolist()))]#subset_df.Instrument.unique()
for instrument in inst:
    if len(list(path_to_data.glob('project_{}_standardizer*'.format(str(instrument)))))==0:
        subset_df_instrument = pd.DataFrame(
            {'Project_ID': subset_df[subset_df['Instrument'].isin(dict_instruments[instrument])]['Project_ID'],
             'Project_source': ['https://ecotaxa.obs-vlfr.fr/prj/'+str(project_id) for project_id in subset_df[subset_df['Instrument'].isin(dict_instruments[instrument])]['Project_ID']],
             'Project_localpath':str(Path(path_to_projects / instrument)).replace(str(Path.home()),'~'),
             'Instrument':instrument,
             'Cruise_field': '',
             'Station_field': '',
             'Profile_field': '',
             'Sample_field': '',
             'Latitude_field': '',
             'Latitude_unit': '',
             'Longitude_field': '',
             'Longitude_unit': '',
             'Depth_min_field': '',
             'Depth_min_unit': '',
             'Depth_max_field': '',
             'Depth_max_unit': '',
             'Sampling_date_field': '',
             'Sampling_date_format': '',
             'Sampling_time_field': '',
             'Sampling_time_format': '',
             'Volume_analyzed_field': '',
             'Volume_analyzed_unit': '',
             'Dilution_field': '',
             'ROI_field': '',
             'ESD_field': '',
             'ESD_unit': '',
             'Biovolume_field': '',
             'Biovolume_unit': '',
             'Area_field': '',
             'Area_unit': '',
             'Minor_axis_field': '',
             'Minor_axis_unit': '',
             'Pixel_field': '',
             'Pixel_unit': '',
             'Category_field': '',
             'Annotation_field': '',
             'Flag_path': '',
             'NA_value': '',
             'EcoPart_project': '',
             'Sampling_type_field': '',
             'Sampling_lower_size_field': '',
             'Sampling_lower_size_unit': '',
             'Sampling_upper_size_field': '',
             'Sampling_upper_size_unit': '',
             'Sampling_description': ''
             })

        subset_df_instrument_metadata = pd.DataFrame(
            {'Variables': subset_df_instrument.columns, 'Variable_types': subset_df_instrument.dtypes,
             'Description': ['Unique Project ID (e.g EcoTaxa ID or cruise ID)',
                             'URL of project data source',
                             'Local path of directory where project native datafile is stored',
                             'Project instrument',
                             'Name of the cruise ID column in project file',
                             'Name of the station ID column in project file',
                             'Name of the profile ID column in project file',
                             'Name of the sample ID column in project file',
                             'Name of the latitude column in project file',
                             'Unit of the latitude column in project file',
                             'Name of the longitude column in project file',
                             'Unit of the longitude column in project file',
                             'Name of the minimum depth column in project file',
                             'Unit of the minimum depth column in project file',
                             'Name of the maximum depth column in project file',
                             'Unit of the maximum depth column in project file',
                             'Name of the sampling date (UTC) column in project file',
                             'Format of the sampling date (UTC) column in project file  (e.g. yyyymmdd)',
                             'Name of the sampling time (UTC) column in project file',
                             'Format of the sampling time(UTC) column in project file(e.g.hhmmss)',
                             'Name of the sample volume analyzed column in project file',
                             'Unit of the sample volume analyzed column in project file',
                             'Name of the dilution factor column in project file',
                             'Name of the region of interest (ROI) ID column in project file',
                             'Name of the equivalent spherical diameter column in project file',
                             'Unit of the equivalent spherical diameter column in project file',
                             'Name of the biovolume column in project file',
                             'Unit of the biovolume column in project file',
                             'Name of the area column in project file',
                             'Unit of the area column in project file',
                             'Name of the minor ellipsoidal axis column in project file',
                             'Unit of the minor ellipsoidal axis column in project file',
                             'Name of the pixel conversion factor column in project file',
                             'Unit of the pixel conversion factor column in project file',
                             'Name of the category ID column in project file',
                             'Name of the annotation status column in project file (e.g predicted, validated)',
                             'Path of the file containing the ID of flagged samples to be removed from standardized project. Automatic flagging available in funcs_standardize_projects.py',
                             'Character for missing values',
                             'ID of corresponding EcoPart project(s) (*UVP-only). Use comma to separate multiple projects',
                             'Name of the sampling type column in project file (e.g. platform, gear, strategy). Used to describe the sampling method',
                             'Name of the smallest sampled size column in project file (e.g. mesh size). Used to describe the sampling collection threshold',
                             'Unit of the smallest sampled size column in project file (e.g. mesh size). Used to describe the sampling collection threshold',
                             'Name of the largest sampled size column in project file (e.g. mesh size). Used to describe the sampling collection threshold',
                             'Unit of the largest sampled size column in project file (e.g. mesh size). Used to describe the sampling collection threshold',
                             'Additoinal description of the sampling method or protocol (e.g. net_aperture:{field:net_surf,unit:square_meter}, fixative:{chemical:glutaraldehyde,concentration:0.1%}, reference: {https://doi.org/...,https://www.protocols.io/view/zooscan-protocol-yxmvmk8j9g3p/v1})']})
        print("Creating project_", instrument, "_standardizer file, please wait", sep="")
        # Generating instrument-specific standardizer file
        with pd.ExcelWriter(str(path_to_data / 'project_{}_standardizer.xlsx'.format(str(instrument))), engine="xlsxwriter") as writer:
            subset_df_instrument.to_excel(writer, sheet_name='Data', index=False)
            subset_df_instrument_metadata.to_excel(writer, sheet_name='Metadata', index=False)

    else:
            subset_df_instrument=pd.read_excel(list(path_to_data.glob('project_{}_standardizer*'.format(str(instrument))))[0], sheet_name='Data')
            subset_df_instrument_metadata=pd.read_excel(list(path_to_data.glob('project_{}_standardizer*'.format(str(instrument))))[0], sheet_name='Metadata')
            print("Updating project_", instrument, "_standardizer file, please wait", sep="")
            if any((subset_df['Instrument'].isin(dict_instruments[instrument])) & (subset_df['Project_ID'].isin(subset_df_instrument.Project_ID.to_list())==False)):
               df=pd.DataFrame({'Project_ID': subset_df[(subset_df['Instrument'].isin(dict_instruments[instrument])) & (subset_df['Project_ID'].isin(subset_df_instrument.Project_ID.to_list())==False)]['Project_ID'],
                              'Project_source': ['https://ecotaxa.obs-vlfr.fr/prj/'+str(project_id) for project_id in subset_df[(subset_df['Instrument'].isin(dict_instruments[instrument])) & (subset_df['Project_ID'].isin(subset_df_instrument.Project_ID.to_list())==False)]['Project_ID']],
                              'Project_localpath':str(Path(path_to_projects / instrument)).replace(str(Path.home()),'~'),
                              'Instrument': instrument,
                              'Cruise_field': '',
                              'Station_field': '',
                              'Profile_field': '',
                              'Sample_field': '',
                              'Latitude_field': '',
                              'Latitude_unit': '',
                              'Longitude_field': '',
                              'Longitude_unit': '',
                              'Depth_min_field': '',
                              'Depth_min_unit': '',
                              'Depth_max_field': '',
                              'Depth_max_unit': '',
                              'Sampling_date_field': '',
                              'Sampling_date_format': '',
                              'Sampling_time_field': '',
                              'Sampling_time_format': '',
                              'Volume_analyzed_field': '',
                              'Volume_analyzed_unit': '',
                              'Dilution_field': '',
                              'ROI_field': '',
                              'ESD_field': '',
                              'ESD_unit': '',
                              'Biovolume_field': '',
                              'Biovolume_unit': '',
                              'Area_field': '',
                              'Area_unit': '',
                              'Minor_axis_field': '',
                              'Minor_axis_unit': '',
                              'Pixel_field': '',
                              'Pixel_unit': '',
                              'Category_field': '',
                              'Annotation_field': '',
                              'Flag_path': '',
                              'NA_value': '',
                              'EcoPart_project': '',
                              'Sampling_type_field': '',
                              'Sampling_lower_size_field': '',
                              'Sampling_lower_size_unit': '',
                              'Sampling_upper_size_field': '',
                              'Sampling_upper_size_unit': '',
                              'Sampling_description': ''
                              })

               subset_df_instrument=pd.concat([subset_df_instrument, df], axis=0)
            else:
               df=subset_df_instrument
         # Generating instrument-specific standardizer file
            with pd.ExcelWriter(str(path_to_data / 'project_{}_standardizer.xlsx'.format(str(instrument))),engine="xlsxwriter") as writer:
                       subset_df_instrument.to_excel(writer, sheet_name='Data', index=False)
                       subset_df_instrument_metadata.to_excel(writer, sheet_name='Metadata', index=False)
print('Please fill out standardizer for each project. Read the metadata spreadsheet to get variable description.\n Use unit from the list below or define custom units in {}:\n'.format(path_to_git/'scripts'/'units_def.txt'),full_list_units,sep='')