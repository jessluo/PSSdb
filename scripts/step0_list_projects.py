## Objective: Generate a project_list_all and instrument-specific project standardizer using Ecotaxa, IFCB API

## TO DO: Add projects downloaded from other platforms (e.g. IFCB dashboard)


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

# Portal-specifc functions to list available projects
from funcs_list_projects import *

# Workflow starts here

# Generate an empty dataframe with list of Project url, ID, title from Ecotaxa and IFCB dashboard for projects standardizers.
# Not necessary for Ecopart since the output are already formatted (unique format across projects) and already binned
df_projects_list=pd.DataFrame({})

# read git-tracked config file (text file) with inputs:  git directory and data subdirectory, username/password for Ecotaxa/Ecopart authentication
path_to_config=Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
# read config file (text file) with inputs:  EcoTaxa username, password. Make sure to secure this file if needed
path_to_config_usr=Path('~/GIT/PSSdb/scripts/Ecotaxa_API_pw.yaml').expanduser()
with open(path_to_config_usr ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

path_to_git=Path(cfg['git_dir']).expanduser()

# Set storage path based on path stored in the yaml config file
path_to_data = path_to_git / cfg['dataset_subdir']
path_to_ecotaxa_projects=path_to_data / cfg['Ecotaxa_subdir']
path_to_ecopart_projects=path_to_data / cfg['Ecopart_subdir']
path_to_ifcb_projects=path_to_data / cfg['IFCB_dir']

#1) List Ecotaxa projects using Scripts.funcs_list_projects.py function: Ecotaxa_list
# Prompting a warning if a project_list_all file already exists
# asking for confirmation to overwrite. Attention: the outdated project list will be erased!
existing_project_path = list(path_to_data.glob('project_list_all.xlsx'))
if len(existing_project_path) != 0:
    sheets = pd.ExcelFile(existing_project_path[0]).sheet_names
    if ('ecotaxa' in sheets) & ('metadata' in sheets):
        df_metadata=pd.read_excel(existing_project_path[0],sheet_name='metadata')
        df_ecotaxa = pd.read_excel(existing_project_path[0], sheet_name='ecotaxa',dtype=dict(zip(df_metadata.Variables, df_metadata.Variable_types)))

        confirmation = input("List of Ecotaxa projects already created. Do you wish to overwrite the existing list? Enter Y or N\n")
        if confirmation == 'Y':
            print("Overwriting list of Ecotaxa projects. Please wait")
            # get projects list
            df_Ecotaxa_list = Ecotaxa_list(username=cfg_pw['ecotaxa_user'], password=cfg_pw['ecotaxa_pass'],localpath=str(path_to_ecotaxa_projects).replace(str(Path.home()),'~'))
            df_metadata=df_Ecotaxa_list['metadata']
            with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                df_Ecotaxa_list['data'].to_excel(writer, sheet_name='ecotaxa', index=False)

        elif confirmation=='N':
            df_Ecotaxa_list={'data':df_ecotaxa}
    else:
        print("Creating list of Ecotaxa projects. Please wait")
        # get projects list
        df_Ecotaxa_list = Ecotaxa_list(username=cfg_pw['ecotaxa_user'], password=cfg_pw['ecotaxa_pass'],localpath=str(path_to_ecotaxa_projects).replace(str(Path.home()),'~'))
        df_metadata = df_Ecotaxa_list['metadata']
        # Generating project_list_all file
        with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
            df_Ecotaxa_list['data'].to_excel(writer, sheet_name='ecotaxa', index=False)
            df_metadata.to_excel(writer, sheet_name='metadata', index=False)

else:
    print("Creating list of projects hosted on Ecotaxa, Ecopart, or IFCB dashboard. Please wait")
    # get projects list
    df_Ecotaxa_list = Ecotaxa_list(username=cfg_pw['ecotaxa_user'], password=cfg_pw['ecotaxa_pass'],localpath=str(path_to_ecotaxa_projects).replace(str(Path.home()),'~'))
    df_metadata = df_Ecotaxa_list['metadata']

    # Generating project_list_all file
    with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'),engine="xlsxwriter") as writer:
        df_Ecotaxa_list['data'].to_excel(writer, sheet_name='ecotaxa', index=False)
        df_metadata.to_excel(writer, sheet_name='metadata', index=False)

df_projects_list=pd.concat([df_projects_list,df_Ecotaxa_list['data'][['Project_source','Project_ID','Instrument','Project_localpath','PSSdb_access']]],axis=0)
#2) List Ecopart projects using Scripts.funcs_list_projects.py function: Ecopart_list
existing_project_path = list(path_to_data.glob('project_list_all.xlsx'))
if len(existing_project_path) != 0:
    sheets = pd.ExcelFile(existing_project_path[0]).sheet_names
    if ('ecopart' in sheets) & ('metadata' in sheets):
        df_metadata=pd.read_excel(existing_project_path[0],sheet_name='metadata')
        df_ecopart = pd.read_excel(existing_project_path[0], sheet_name='ecopart',dtype=dict(zip(df_metadata.Variables,df_metadata.Variable_types)))

        confirmation = input("List of Ecopart projects already created. Do you wish to overwrite the existing list? Enter Y or N\n")
        if confirmation == 'Y':
            print("Overwriting list of Ecopart projects. Please wait")
            # get projects list
            df_Ecopart_list = Ecopart_list(username=cfg_pw['ecotaxa_user'], password=cfg_pw['ecotaxa_pass'],localpath=str(path_to_ecopart_projects).replace(str(Path.home()),'~'))
            df_metadata=pd.concat([df_metadata,df_Ecopart_list['metadata'].reset_index(drop=True)],axis=0).drop_duplicates(['Variables'])
            with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                df_Ecopart_list['data'].to_excel(writer, sheet_name='ecopart', index=False)
                df_metadata.to_excel(writer, sheet_name='metadata', index=False)
        elif confirmation=='N':
            df_Ecopart_list={'data':df_ecopart}
    else:
        print("Creating list of Ecopart projects. Please wait")
        # get projects list
        df_Ecopart_list = Ecopart_list(username=cfg_pw['ecotaxa_user'], password=cfg_pw['ecotaxa_pass'],localpath=str(path_to_ecopart_projects).replace(str(Path.home()),'~'))
        df_metadata = pd.concat([df_metadata, df_Ecopart_list['metadata'].reset_index(drop=True)], axis=0).drop_duplicates(['Variables'])

        # Generating project_list_all file
        with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
            df_Ecopart_list['data'].to_excel(writer, sheet_name='ecopart', index=False)
            df_metadata.to_excel(writer, sheet_name='metadata', index=False)
else:
    print("Creating list of projects hosted on Ecotaxa, Ecopart, or IFCB dashboard. Please wait")
    # get projects list
    df_Ecopart_list = Ecotaxa_list(username=cfg_pw['ecotaxa_user'], password=cfg_pw['ecotaxa_pass'],localpath=str(path_to_ecopart_projects).replace(str(Path.home()),'~'))
    df_metadata = pd.concat([df_metadata, df_Ecopart_list['metadata'].reset_index(drop=True)], axis=0).drop_duplicates(['Variables'])

    # Generating project_list_all file
    with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'), engine="xlsxwriter") as writer:
        df_Ecopart_list['data'].to_excel(writer, sheet_name='ecopart', index=False)
        df_metadata.to_excel(writer, sheet_name='metadata', index=False)

#Appending Ecopart/Ecotaxa project ID to standardizer
df_projects_list['External_project']=df_projects_list.Project_ID.astype(str).apply(lambda id: ';'.join(df_Ecopart_list['data'][df_Ecopart_list['data'].Ecotaxa_ID.astype(str)==str(id)].Project_ID.astype(str).values.tolist()))
df_projects_list=pd.concat([df_projects_list,df_Ecopart_list['data'].rename(columns={'Ecotaxa_ID':'External_project'})[['Project_source','Project_ID','Instrument','Project_localpath','PSSdb_access','External_project']]],axis=0)

#3) List IFCB dashboard projects using Scripts.funcs_list_projects.py function: IFCB_dashboard_list
existing_project_path = list(path_to_data.glob('project_list_all.xlsx'))
if len(existing_project_path) != 0:
    sheets = pd.ExcelFile(existing_project_path[0]).sheet_names
    if ('ifcb' in sheets) & ('metadata' in sheets):
        df_metadata=pd.read_excel(existing_project_path[0],sheet_name='metadata')
        df_ifcb = pd.read_excel(existing_project_path[0], sheet_name='ifcb')

        confirmation = input("List of IFCB dashboard projects already created. Do you wish to overwrite the existing list? Enter Y or N\n")
        if confirmation == 'Y':
            print("Overwriting list of IFCB projects. Please wait")
            # get projects list
            df_ifcb_list = IFCB_dashboard_list(localpath=str(path_to_ifcb_projects).replace(str(Path.home()),'~'))
            df_metadata = pd.concat([df_metadata, df_ifcb_list['metadata'].reset_index(drop=True)], axis=0).drop_duplicates(['Variables'])
            with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                df_ifcb_list['data'].to_excel(writer, sheet_name='ifcb', index=False)
                df_metadata.to_excel(writer, sheet_name='metadata', index=False)

        elif confirmation=='N':
            df_ifcb_list={'data':df_ifcb}
    else:
        print("Creating list of IFCB projects. Please wait")
        # get projects list
        df_ifcb_list = IFCB_dashboard_list(localpath=str(path_to_ifcb_projects).replace(str(Path.home()),'~'))
        df_metadata = pd.concat([df_metadata,  df_ifcb_list['metadata'].reset_index(drop=True)], axis=0).drop_duplicates(['Variables'])
        # Generating project_list_all file
        with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
            df_ifcb_list['data'].to_excel(writer, sheet_name='ifcb', index=False)
            df_metadata.to_excel(writer, sheet_name='metadata', index=False)

else:
    print("Creating list of projects hosted on Ecotaxa, Ecopart, or IFCB dashboard. Please wait")
    # get projects list
    df_ifcb_list = IFCB_dashboard_list(localpath=str(path_to_ifcb_projects).replace(str(Path.home()),'~'))
    df_metadata = pd.concat([df_metadata, df_ifcb_list['metadata'].reset_index(drop=True)], axis=0).drop_duplicates(['Variables'])

    # Generating project_list_all file
    with pd.ExcelWriter(str(path_to_data / 'project_list_all.xlsx'),engine="xlsxwriter") as writer:
        df_ifcb_list['data'].to_excel(writer, sheet_name='ifcb', index=False)
        df_metadata.to_excel(writer, sheet_name='metadata', index=False)

df_projects_list=pd.concat([df_projects_list,df_ifcb_list['data'][['Project_source','Project_ID','Instrument','Project_localpath','PSSdb_access']]],axis=0)

#df_projects_list=pd.concat([df_projects_list,df_IFCB_list['data'][['Project_source','Project_ID','Instrument','Project_localpath','PSSdb_access']]],axis=0)

#4) Generate an instrument-specific standardizer for  accessible projects ()
dict_instruments={'IFCB':['IFCB'],'UVP':['UVP5','UVP5HD','UVP5SD','UVP5Z','UVP6','UVP6REMOTE'],'Zooscan':['Zooscan'],'Unknown':['?'],'AMNIS':['AMNIS'],'CPICS':['CPICS'],'CytoSense':['CytoSense'],'FastCam':['FastCam'],'FlowCam':['FlowCam'],'ISIIS':['ISIIS'],'LISST':['LISST','LISST-Holo'],'Loki':['Loki'],'Other':['Other camera','Other flowcytometer','Other microscope','Other scanner'],'PlanktoScope':['PlanktoScope'],'VPR':['VPR'],'ZooCam':['ZooCam'],'eHFCM':['eHFCM']}
df_projects_list['instrument']=df_projects_list.Instrument.apply(lambda instrument: [key for key,values in dict_instruments.items() if instrument in values][0])
if df_projects_list.PSSdb_access.dtypes!='bool':
    df_projects_list.PSSdb_access=df_projects_list.PSSdb_access=='True'
inst=list(set(df_projects_list[df_projects_list.PSSdb_access==True].instrument.tolist()))
inst_all=df_projects_list.instrument.tolist()

print("Building/Updating standardizer for",'/'.join(inst),"projects",sep=" ")
subset_df=df_projects_list[df_projects_list.PSSdb_access==True]
for instrument in inst:
    if len(list(path_to_data.glob('project_{}_standardizer*'.format(str(instrument)))))==0:
        subset_df_instrument = pd.DataFrame(
            {'Project_ID': subset_df[subset_df['Instrument'].isin(dict_instruments[instrument])]['Project_ID'],
             'Project_source':subset_df[subset_df['Instrument'].isin(dict_instruments[instrument])]['Project_source'],
             'Project_localpath':subset_df[subset_df['Instrument'].isin(dict_instruments[instrument])]['Project_localpath'],
             'Instrument':instrument})

        df_sub = pd.concat(map(lambda y: pd.DataFrame.from_dict({
            'Cruise_field': '',
            'Station_field': '',
            'Profile_field': '',
            'Sample_field':'',
            'Latitude_field': '',
            'Latitude_unit': '',
            'Longitude_field': '',
            'Longitude_unit': '',
            'Depth_min_field': '',
            'Depth_min_unit': '',
            'Depth_max_field': '',
            'Depth_max_unit': '',
            'Sampling_date_field':'',
            'Sampling_date_format': '',
            'Sampling_time_field':'',
            'Sampling_time_format': '',
            'Volume_analyzed_field': '',
            'Volume_analyzed_unit': '',
            'Dilution_field':'',
            'ROI_field':'',
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
            'Annotation_field':'',
            'Flag_path':'',
            'NA_value':'',
            'EcoPart_project': subset_df[(subset_df['Instrument'].isin(dict_instruments[instrument])) & (subset_df['Project_ID'].isin([y]))]['Ecopart_project'].values.tolist()[0],
            'Sampling_type_field':'',
            'Sampling_lower_size_field': '',
            'Sampling_lower_size_unit':'',
            'Sampling_upper_size_field': '',
            'Sampling_upper_size_unit': '',
            'Sampling_description':''
        }, orient="index").T,subset_df[subset_df['Instrument'].isin(dict_instruments[instrument])]['Project_ID']), ignore_index=True)
        df_sub.index = subset_df_instrument.index.tolist()
        subset_df_instrument = pd.concat([subset_df_instrument, df_sub], axis=1)
        subset_df_instrument_metadata = pd.DataFrame(
            {'Variables': subset_df_instrument.columns, 'Variable_types': subset_df_instrument.dtypes,
             'Description': ['Unique Project ID (e.g EcoTaxa ID or cruise ID)',
                             'URL of project data source',
                             'Local path of directory where projects will be exported',
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
                             'Additional description of the sampling method or protocol (e.g. net_aperture:{field:net_surf,unit:square_meter}, fixative:{chemical:glutaraldehyde,concentration:0.1%}, reference: {https://doi.org/...,https://www.protocols.io/view/zooscan-protocol-yxmvmk8j9g3p/v1})']})
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
                              'Project_source':subset_df[(subset_df['Instrument'].isin(dict_instruments[instrument])) & (subset_df['Project_ID'].isin(subset_df_instrument.Project_ID.to_list())==False)]['Project_source'],
                              'Project_localpath':subset_df[(subset_df['Instrument'].isin(dict_instruments[instrument])) & (subset_df['Project_ID'].isin(subset_df_instrument.Project_ID.to_list())==False)]['Project_localpath'],
                              'Instrument': instrument})

               df_sub = pd.concat(map(lambda y: pd.DataFrame.from_dict({
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
                              'Volume_analyzed_field':'',
                              'Volume_analyzed_unit': '',
                              'Dilution_field': '',
                              'ROI_field': '',
                              'ESD_field': '',
                              'ESD_unit': '',
                              'Biovolume_field':'',
                              'Biovolume_unit': '',
                              'Area_field': '',
                              'Area_unit': '',
                              'Minor_axis_field':'',
                              'Minor_axis_unit': '',
                              'Pixel_field': '',
                              'Pixel_unit': '',
                              'Category_field': '',
                              'Annotation_field': '',
                              'Flag_path': '',
                              'NA_value':'',
                              'EcoPart_project': subset_df[(subset_df['Instrument'].isin(dict_instruments[instrument])) & (subset_df['Project_ID'].isin([y]))]['Ecopart_project'].values.tolist()[0],
                              'Sampling_type_field':'',
                              'Sampling_lower_size_field': '',
                              'Sampling_lower_size_unit':'',
                              'Sampling_upper_size_field': '',
                              'Sampling_upper_size_unit': '',
                              'Sampling_description':''
                              }, orient="index").T,subset_df[(subset_df['Instrument'].isin(dict_instruments[instrument])) & (subset_df['Project_ID'].isin(subset_df_instrument.Project_ID.to_list())==False)]['Project_ID']), ignore_index=True)
               df_sub.index = df.index.tolist()
               df = pd.concat([df, df_sub], axis=1)
               subset_df_instrument=pd.concat([subset_df_instrument, df], axis=0)
            else:
               df=subset_df_instrument
         # Generating instrument-specific standardizer file
            with pd.ExcelWriter(str(path_to_data / 'project_{}_standardizer.xlsx'.format(str(instrument))),engine="xlsxwriter") as writer:
                       subset_df_instrument.to_excel(writer, sheet_name='Data', index=False)
                       subset_df_instrument_metadata.to_excel(writer, sheet_name='Metadata', index=False)