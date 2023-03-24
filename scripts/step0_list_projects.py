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
try:
    from funcs_list_projects import *
except:
    from scripts.funcs_list_projects import *

columns_for_flag_summary=['Number_samples', 'Number_flagged_samples', 'Percentage_flagged_missing', 'Percentage_flagged_GPScoordinatesonland','Percentage_flagged_dubiousGPScoordinates', 'Percentage_flagged_count','Percentage_flagged_artefact','Percentage_flagged_validation', 'Percentage_flagged_size']

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

df_projects_list=pd.concat([df_projects_list,df_Ecotaxa_list['data'][['Project_source','Project_ID','Instrument','Project_localpath','PSSdb_access','Contact_name']+[column for column in columns_for_flag_summary if column in df_Ecotaxa_list['data'].columns.tolist()]]],axis=0)
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
df_projects_list=pd.concat([df_projects_list,df_Ecopart_list['data'].rename(columns={'Ecotaxa_ID':'External_project'})[['Project_source','Project_ID','Instrument','Project_localpath','PSSdb_access','External_project','Contact_name']+[column for column in columns_for_flag_summary if column in df_Ecopart_list['data'].columns.tolist()]]],axis=0).reset_index(drop=True)

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

df_projects_list=pd.concat([df_projects_list,df_ifcb_list['data'][['Project_source','Project_ID','Instrument','Project_localpath','PSSdb_access','Contact_name']+[column for column in columns_for_flag_summary if column in df_ifcb_list['data'].columns.tolist()]]],axis=0).reset_index(drop=True)

#df_projects_list=pd.concat([df_projects_list,df_IFCB_list['data'][['Project_source','Project_ID','Instrument','Project_localpath','PSSdb_access']]],axis=0)

#4) Generate an instrument-specific standardizer for  accessible projects ()
dict_instruments={'IFCB':['IFCB'],'UVP':['UVP5','UVP5HD','UVP5SD','UVP5Z','UVP6','UVP6REMOTE'],'Zooscan':['Zooscan'],'Unknown':['?'],'AMNIS':['AMNIS'],'CPICS':['CPICS'],'CytoSense':['CytoSense'],'FastCam':['FastCam'],'FlowCam':['FlowCam'],'ISIIS':['ISIIS'],'LISST':['LISST','LISST-Holo'],'Loki':['Loki'],'Other':['Other camera','Other flowcytometer','Other microscope','Other scanner'],'PlanktoScope':['PlanktoScope'],'VPR':['VPR'],'ZooCam':['ZooCam'],'eHFCM':['eHFCM']}
df_projects_list['instrument']=df_projects_list.Instrument.apply(lambda instrument: [key for key,values in dict_instruments.items() if instrument in values][0] if str(instrument)!='nan' else '')
if df_projects_list.PSSdb_access.dtypes!='bool':
    df_projects_list.PSSdb_access=df_projects_list.PSSdb_access.astype(str).str.capitalize()=='True'
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
            'External_project': subset_df[(subset_df['Instrument'].isin(dict_instruments[instrument])) & (subset_df['Project_ID'].isin([y]))]['External_project'].values.tolist()[0],
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
                             'ID of corresponding EcoPart/Ecotaxa project(s) (*UVP-only). Use semicolon to separate multiple projects',
                             'Name of the sampling type column in project file (e.g. platform, gear, strategy). Used to describe the sampling method',
                             'Name of the smallest sampled size column in project file (e.g. mesh size). Used to describe the sampling collection threshold',
                             'Unit of the smallest sampled size column in project file (e.g. mesh size). Used to describe the sampling collection threshold',
                             'Name of the largest sampled size column in project file (e.g. mesh size). Used to describe the sampling collection threshold',
                             'Unit of the largest sampled size column in project file (e.g. mesh size). Used to describe the sampling collection threshold',
                             'Additional description of the sampling method or protocol (e.g. net_aperture:{field:net_surf,unit:square_meter}, fixative:{chemical:glutaraldehyde,concentration:0.1%}, reference: {https://doi.org/...,https://www.protocols.io/view/zooscan-protocol-yxmvmk8j9g3p/v1})']})
        print("Creating project_", instrument, "_standardizer file, please wait", sep="")
        if instrument=='UVP':
            # Generating instrument-specific and portal standardizer file
            with pd.ExcelWriter(str(path_to_data / 'project_{}_standardizer.xlsx'.format(str(instrument))), engine="xlsxwriter") as writer:
                for portal in subset_df_instrument.Project_localpath.apply(lambda path: Path(path).name.lower()).unique():
                    subset_df_instrument_portal=subset_df_instrument[subset_df_instrument.Project_localpath.apply(lambda path: Path(path).name.lower())==portal]
                    subset_df_instrument_portal.to_excel(writer, sheet_name=portal, index=False)
                subset_df_instrument_metadata.to_excel(writer, sheet_name='Metadata', index=False)
        else:
            # Generating instrument-specific standardizer file
            with pd.ExcelWriter(str(path_to_data / 'project_{}_standardizer.xlsx'.format(str(instrument))), engine="xlsxwriter") as writer:
                subset_df_instrument.to_excel(writer, sheet_name='Data', index=False)
                subset_df_instrument_metadata.to_excel(writer, sheet_name='Metadata', index=False)

    else:
            if instrument != 'UVP':
                subset_df_instrument=pd.read_excel(list(path_to_data.glob('project_{}_standardizer*'.format(str(instrument))))[0], sheet_name='Data')
            else:
                subset_df_instrument = pd.concat([pd.read_excel(list(path_to_data.glob('project_{}_standardizer*'.format(str(instrument))))[0],sheet_name='ecotaxa'), pd.read_excel(list(path_to_data.glob('project_{}_standardizer*'.format(str(instrument))))[0],sheet_name='ecopart')], axis=0)

            subset_df_instrument_metadata=pd.read_excel(list(path_to_data.glob('project_{}_standardizer*'.format(str(instrument))))[0], sheet_name='Metadata')
            print("Updating project_", instrument, "_standardizer file, please wait", sep="")
            # Check for project we lost access to and remove from standardizer, export/standardized/flag files directories
            for portal in subset_df_instrument.Project_localpath.unique():
                df_projects_list_sub=df_projects_list[df_projects_list.Project_localpath==portal]
                if any(df_projects_list_sub[(df_projects_list_sub.PSSdb_access==False)  & (df_projects_list_sub.Instrument.isin(dict_instruments[instrument]))].Project_ID.isin(subset_df_instrument[subset_df_instrument.Project_localpath==portal].Project_ID.to_list())):
                    project_revoked=df_projects_list_sub[(df_projects_list_sub.PSSdb_access == False) & (df_projects_list_sub.Instrument.isin(dict_instruments[instrument]))].loc[df_projects_list_sub[(df_projects_list_sub.PSSdb_access == False) & (df_projects_list_sub.Instrument.isin(dict_instruments[instrument]))].Project_ID.isin(subset_df_instrument[subset_df_instrument.Project_localpath==portal].Project_ID.to_list()),'Project_ID']
                    contact_revoked=df_projects_list_sub[(df_projects_list_sub.PSSdb_access == False) & (df_projects_list_sub.Instrument.isin(dict_instruments[instrument]))].loc[df_projects_list_sub[(df_projects_list_sub.PSSdb_access == False) & (df_projects_list_sub.Instrument.isin(dict_instruments[instrument]))].Project_ID.isin(subset_df_instrument[subset_df_instrument.Project_localpath==portal].Project_ID.to_list()),'Contact_name']
                    print("Project with revoked access rights found in standardizer spreadsheets. Removing project from standardizer, export/standardized/flag files directories before continuing:\n{}".format('\n'.join(list(map(lambda proj:'Project:{}, Contact: {}'.format(proj[0],proj[1]),dict(zip(project_revoked.values,contact_revoked.values)).items())))))
                    subset_df_instrument=subset_df_instrument.drop(index=np.argwhere(((subset_df_instrument.Project_localpath==portal) & (subset_df_instrument.Project_ID.isin(project_revoked.values.tolist()).values)).values)[0]).reset_index(drop=True)

                    for id in df_projects_list_sub[(df_projects_list_sub.PSSdb_access == False) & (df_projects_list_sub.Instrument.isin(dict_instruments[instrument]))].loc[df_projects_list_sub[(df_projects_list_sub.PSSdb_access == False) & (df_projects_list_sub.Instrument.isin(dict_instruments[instrument]))].Project_ID.isin(subset_df_instrument[subset_df_instrument.Project_localpath==portal].Project_ID.to_list())].index:
                        exportfiles_path = list((Path(df_projects_list_sub.at[id,'Project_localpath']).expanduser() /[ind for ind in dict_instruments.keys() if df_projects_list_sub.at[id,'Project_localpath'].at[id, 'Instrument'] in dict_instruments.get(ind)][0] ).glob('ecotaxa_export_{}_*.tsv'.format(df_projects_list_sub.at[id,'Project_localpath'].at[id,'Project_ID']))) if (Path(portal).stem.lower()=='ecotaxa') else list((Path(df_projects_list_sub.at[id,'Project_localpath']).expanduser() ).rglob('ecopart_export_raw_{}_*.tsv'.format(df_projects_list_sub.at[id,'Project_localpath'].at[id,'Project_ID'])))
                        if len(exportfiles_path):
                            [file.unlink(missing_ok=True) for file in exportfiles_path]
                        flagfiles_path = list((Path(df_projects_list_sub.at[id,'Project_localpath']).expanduser().parent /'flags'/'ecotaxa' ).glob('project_{}_flags.csv'.format(df_projects_list.at[id,'Project_ID'])))+list((Path(df_projects_list.at[id,'Project_localpath']).expanduser().parent /'flags'/'ecotaxa_ecopart_consolidation' ).glob('project_{}_flags.csv'.format(df_projects_list_sub.at[id,'Project_localpath'].at[id,'Project_ID']))) if (Path(portal).stem.lower()=='ecotaxa') else []
                        if len(flagfiles_path ):
                            [file.unlink(missing_ok=True) for file in flagfiles_path]
                        standardizedfiles_path = list((Path(df_projects_list_sub.at[id,'Project_localpath']).expanduser().parent /'raw_standardized'/'ecotaxa' /[ind for ind in dict_instruments.keys() if df_projects_list.at[id, 'Instrument'] in dict_instruments.get(ind)][0] ).rglob('standardized_project_{}_*'.format(df_projects_list.at[id,'Project_ID'])))+list((Path(df_projects_list_sub.at[id,'Project_localpath'].at[id,'Project_localpath']).expanduser().parent /'raw_standardized'/'ecotaxa_ecopart_consolidation' ).rglob('standardized_project_{}_*'.format(df_projects_list_sub.at[id,'Project_localpath'].at[id,'Project_ID']))) if (Path(portal).stem.lower()=='ecotaxa') else []
                        if len(standardizedfiles_path ):
                            [file.unlink(missing_ok=True) for file in standardizedfiles_path]

            # Check for additional projects
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
                              'External_project': subset_df[(subset_df['Instrument'].isin(dict_instruments[instrument])) & (subset_df['Project_ID'].isin([y]))]['External_project'].values.tolist()[0],
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
            if instrument != 'UVP':
                with pd.ExcelWriter(str(path_to_data / 'project_{}_standardizer.xlsx'.format(str(instrument))), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                       subset_df_instrument.to_excel(writer, sheet_name='Data', index=False)
                       subset_df_instrument_metadata.to_excel(writer, sheet_name='Metadata', index=False)
            else:
                portal_dict = {'ecotaxa': ['ecotaxa', 'ecotaxa_ecopart_consolidation'], 'ecopart': ['ecopart']}
                with pd.ExcelWriter(str(path_to_data / 'project_{}_standardizer.xlsx'.format(str(instrument))), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                    for portal in portal_dict.items():#subset_df_instrument.Project_localpath.apply(lambda path: Path(path).name.lower()).unique():
                        subset_df_instrument_portal = subset_df_instrument[subset_df_instrument.Project_localpath.apply(lambda path: Path(path).name.lower() in portal[1])]
                        subset_df_instrument_portal.to_excel(writer, sheet_name=portal[0], index=False)
                    subset_df_instrument_metadata.to_excel(writer, sheet_name='Metadata', index=False)

