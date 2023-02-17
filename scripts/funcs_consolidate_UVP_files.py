## Objective: This script provide functions to match-up Ecotaxa vignettes with bru files and generate extended Ecotaxa table with small particles based on the bru files

## Documentation: This script is based on Zooprocess version 8.12. See https://sites.google.com/view/piqv/softwares/uvp5?authuser=0
## and EcoPart py module. See https://github.com/ecotaxa/ecopart/blob/master/py/part_app/funcs/uvp_sample_import.py

## Python modules
import os
import io
import itertools
from itertools import compress
from natsort import natsorted
import pandas as pd
import numpy as np # Arrays
from pathlib import Path # Handling of path object
import re
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore",category=FutureWarning)
# PyImageJ package: Python wrapper for ImageJ2 (https://pypi.org/project/pyimagej/)
import sys
sys.path.extend(['~/opt/apache-maven-3.8.7/bin']) # Download Maven from : https://maven.apache.org/
import imagej
ij = imagej.init('net.imagej:imagej:2.5.0')  # Attention: mode='gui' not working on Mac
ij.getVersion()  # To check if pyimagej environment had been activated during initialization

# To install:
# conda install mamba -n base -c conda-forge
# mamba create -n pyimagej -c conda-forge pyimagej openjdk=8
# conda activate pyimagej
# The pyimageJ environment is created in ~/opt.anaconda3/envs/pyimagej
# Alternatively, use pip install pyimagej
# Download Maven (https://maven.apache.org/) and add to path in your terminal: export M2_HOME="~/opt/apache-maven-3.8.7"
# export PATH=${PATH}:${M2_HOME}/bin
# Modules for sftp server connection to create UVP metadata:
import yaml  # pip install pysftp
import shutil
# Ecotaxa API for project update:
import ecotaxa_py_client
from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.model.login_req import LoginReq
from ecotaxa_py_client.model.import_req import ImportReq
from ecotaxa_py_client.api import projects_api
from ecotaxa_py_client.api import files_api
from ecotaxa_py_client.api import jobs_api

# Modules for webpage handling/scraping:
import urllib3
import requests
from bs4 import BeautifulSoup # Use pip install bs4 and pip install lxml

from tqdm import tqdm
import datetime, time # Time module for break and current date

## Path definition:
path_to_git=Path('~/GIT/PSSdb').expanduser()

path_to_config=path_to_git/ 'scripts'/ 'Ecotaxa_API.yaml'
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)

path_to_config_pw=path_to_git/ 'scripts'/ 'Ecotaxa_API_pw.yaml'
with open(path_to_config_pw ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)
password_ecotaxa=cfg_pw['ecotaxa_pass']
user_ecotaxa=cfg_pw['ecotaxa_user']

path_to_macro=path_to_git/ 'scripts' / 'PyImageJ_ecotaxa_append_metadata.txt'
# Consolidate Ecotaxa/Ecopart projects using the UVP standardizer spreadsheet
path_to_standardizer=path_to_git/ 'raw' / 'project_UVP_standardizer.xlsx'
if path_to_standardizer.exists():
     df_standardizer_ecotaxa = pd.read_excel(path_to_standardizer,sheet_name='ecotaxa', index_col=0)
     columns_for_description = df_standardizer_ecotaxa['Sampling_description'].apply(lambda description: dict(zip([item.split(':', 1)[0].capitalize() for item in description.split(';')],[eval(re.sub(r'(\w+)', r'"\1"', item.split(':', 1)[1])) for item in description.split(';')])) if str(description) != 'nan' else description)
     fields_for_description = columns_for_description.apply(lambda description: pd.Series([description[key][index] for key, value in description.items() if (type(value) == dict) for index, fields in description[key].items() if ('field' in index)],[key.capitalize() for key, value in description.items() if (type(value) == dict) for index, fields in description[key].items() if ('field' in index)]) if description != 'nan' else pd.Series({}))
     df_standardizer_ecotaxa = pd.concat([df_standardizer_ecotaxa, fields_for_description], axis=1)

     df_standardizer_ecopart = pd.read_excel(path_to_standardizer, sheet_name='ecopart', index_col=0)
     columns_for_description = df_standardizer_ecopart['Sampling_description'].apply(lambda description: dict(zip([item.split(':', 1)[0].capitalize() for item in description.split(';')],[eval(re.sub(r'(\w+)', r'"\1"', item.split(':', 1)[1])) for item in description.split(';')])) if str(description) != 'nan' else description)
     fields_for_description = columns_for_description.apply(lambda description: pd.Series([description[key][index] for key, value in description.items() if (type(value) == dict) for index, fields in description[key].items() if ('field' in index)],[key.capitalize() for key, value in description.items() if (type(value) == dict) for index, fields in description[key].items() if ('field' in index)]) if description != 'nan' else pd.Series({}))
     df_standardizer_ecopart =pd.concat([df_standardizer_ecopart,fields_for_description],axis=1)
# Set headers for bru and dat UVP files
dat_headers=['image_index','image','pressure','tilt_1','tilt_2','board_temperature','voltage','unknown_1','unknown_2','unknown_3','unknown_4','cooler_temperature','camera_temperature','internal_time','nb_blobs_SMbase','mean_area_SMbase','mean_grey_SMbase','nb_blobs_SMzoo','mean_grey_SMzoo']
bru_headers=['image_index','blob','area','mean_grey','x_center','y_center']

# Select native Ecotaxa file variables to be dropped during consoldation (otherwise consolidated files are big)
ecotaxa_fields_to_drop=['object_link','object_annotation_person_name','object_annotation_person_email','object_annotation_date','object_annotation_time','object_mean','object_stdev','object_mode','object_min','object_max','object_intden','object_median','object_skew','object_kurt','object_fractal','object_skelarea','object_slope','object_histcum1','object_histcum2','object_histcum3','object_xmg5','object_ymg5','object_nb1','object_nb2','object_nb3','object_compentropy','object_compmean','object_compsloope','object_compm1','object_compm2','object_compm3','object_symetrieh','object_symetriev','object_symetrievc','object_convperim','object_convarea','object_fcons','object_thickr','object_tag','object_range','object_meanpos','object_centroids','object_cv','object_sr','object_circex','object_cdexc','object_kurt_mean','object_skew_mean','object_symetrieh_area','object_smetriev_area','object_nb1_area','object_nb2_area','object_nb3_area','object_nb1_range','object_nb2_range','object_nb3_range','object_median_mean_range','object_skeleton_area','sample_dataportal_descriptor','sample_winddir','sample_windspeed','sample_seastate','sample_nebulousness','acq_file_description','acq_disktype','acq_shutterspeed','acq_gain','acq_exposure']

## Functions start here:

def update_ecotaxa_project(project_id,path_to_metadata,login=user_ecotaxa,password=password_ecotaxa):
    """
    Objective: This function will merge ecotaxa metadata to existing project using custom web crawling codes.\n
    Metadata are generated by running ImageJ reprocessing macro (pick option according to your need): \n
      -PyImageJ_ecotaxa_append_metadata.txt: to be run from python using function run_imageJ_macro below \n
      -ImageJ_ecotaxa_append_metadata.txt: to be run from pan ImageJ local application \n
    :param project_ID: Project ID on Ecotaxa.
    :param path_to_metadata: Full path to metadata zip file.
    :param login: Ecotaxa account email
    :param password: Ecotaxa account password
    :return: Returns NULL or 200 upon success, 422 for error
    """
    with requests.Session() as session:
        # Authenticate with email/password
        data_login = {"email": login, "password": password}
        session.post('https://ecotaxa.obs-vlfr.fr/login', data=data_login, headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
        # Check that access is granted and export is allowed
        url_import="https://ecotaxa.obs-vlfr.fr/Job/Create/ImportUpdate?p={}".format(project_id)
        response = session.get(url=url_import)
        if response.ok:
            # Start import task with metadata.zip file
            r = session.post(url_import,files={'uploadfile':open(path_to_metadata,'rb')})
            task = r.text[(r.text.find("jqxhr.statusText);\n          });\n          StartMonitor(") + 56):(r.text.find(");\n      });\n  </script>\n\n  ") + 0)]
            # Check status of task
            status = session.get('https://ecotaxa.obs-vlfr.fr/Job/GetStatus/{}'.format(task)).json()
            with tqdm(desc='Working on updated project', total=1000, bar_format='{desc}{bar}', position=0, leave=True) as bar:
                while ('IsComplete' not in status['d'].keys()):
                    time.sleep(3)
                    status = session.get('https://ecotaxa.obs-vlfr.fr/Job/GetStatus/{}'.format(task)).json()
                    percent = status['d']['PercentComplete']
                    bar.set_description("Working on updated project %s%%" % percent, refresh=True)
                    # and update progress bar
                    ok = bar.update(n=10)
                    if ('IsComplete' in status['d'].keys()):
                        if (status['d']['IsComplete'] == 'Y'):
                            break
            print("Update of project {} is complete. You may check and export the project on https://ecotaxa.obs-vlfr.fr/prj/{}".format(project_id,project_id))

def update_ecotaxa_project_API(project_id,path_to_metadata,login=user_ecotaxa,password=password_ecotaxa):
    """
    Objective: This function will merge ecotaxa metadata to existing project using Ecotaxa API.\n
    Metadata are generated by running ImageJ reprocessing macro (pick option according to your need): \n
      -PyImageJ_ecotaxa_append_metadata.txt: to be run from python using function run_imageJ_macro below \n
      -ImageJ_ecotaxa_append_metadata.txt: to be run from pan ImageJ local application \n
    :param project_ID: Project ID on Ecotaxa.
    :param path_to_metadata: Full path to metadata zip file.
    :param login: Ecotaxa account email
    :param password: Ecotaxa account password
    :return: Returns NULL or 200 upon success, 422 for error
    """
    # Generate a token with authentication info
    with ecotaxa_py_client.ApiClient() as client:
        api = authentification_api.AuthentificationApi(client)
        token = api.login(LoginReq(username=login, password=password))
    configuration = ecotaxa_py_client.Configuration(host="https://ecotaxa.obs-vlfr.fr/api", access_token=token, discard_unknown_keys=True)
    configuration.verify_ssl = False
    # Upload and import zip file to existing project
    with ecotaxa_py_client.ApiClient(configuration) as api_client:
        # Create a file instance of the API class for zip file upload
        api_file = files_api.FilesApi(api_client)
        file = open(str(path_to_metadata), 'rb')  # file_type |
        file_response = api_file.post_user_file(file)

        # Create a project instance of the API class
        api_instance = projects_api.ProjectsApi(api_client)
        import_req = ImportReq(source_path=file_response,update_mode="Yes")
        api_response = api_instance.import_file(project_id, import_req)
        api_jobinstance = jobs_api.JobsApi(api_client)
    job_id=api_response['job_id']
    with tqdm(desc='Working on project metadata import file', total=1000, bar_format='{desc}{bar}', position=0,leave=True) as bar:
        job_status = 'R'  # percent = 0
        while job_status not in ('F', 'E'):  # percent!=100:
            time.sleep(2)  # Check the job status every 2 seconds. Modify as needed
            thread = api_jobinstance.get_job(job_id, async_req=True)
            result = thread.get()
            job_status = result.state
            percent = result.progress_pct
            bar.set_description("Working on project metadata import file %s%%" % percent, refresh=True)
            # and update progress bar
            ok = bar.update(n=10)
            if job_status == 'F':
                break
            if job_status == 'E':
                print("Error creating job. Please check your connection or EcoTaxa server")

    print("Update of project {} is complete. You may check and export the project on https://ecotaxa.obs-vlfr.fr/prj/{}".format(project_id,project_id))

def run_imageJ_macro(macro_path,arguments=None):
    """
    Objective: This function will run ImageJ macro from a python script using the PyImageJ package.\n
    The PyImageJ package can be initialize with online or local versions of ImageJ/Fiji software \n
    :param macro_path: ImageJ macro (.txt or .ijm) to be run.
    :param arguments: Macro arguments. Default is None.
    :return: Return depends on the macro
    """
    with open(macro_path, 'r') as file:
        imageJ_macro = file.read()
    macro=ij.py.run_macro(imageJ_macro, args=arguments)
    return macro


def consolidate_ecotaxa_project(project_id,standardizer=df_standardizer_ecotaxa,sheetname='ecotaxa',run_macro=True,upload_metadata=False,localpath_metadata=Path(cfg['raw_dir']).expanduser()/cfg['UVP_consolidation_subdir']):
   """
   Objective: This function consolidate Ecotaxa/Ecopart raw export datafiles.\n
   Additional rows will include particles detected and segmented by Zooprocess/UVPapp but too small to generate vignettes and vice-versa.\n
   Additional variable will include:\n
    -object_bru_id: unique ROI identifier based on image name (yyyymmddhhmmss_ms) and the blob index included in the original bru files\n
    -object_bru_area: area (in square pixels) measured after application of first threshold\n
    -sample_corrected_depth_bin: 1m depth bins (m) defined in Ecopart raw export file\n
    -sample_image_number_depth_bin: Number of images in individual depth bins. This will be used to re-calculate bins sampled volume\n
   :param project_id: project ID.
   :param standardizer: standardizer spreadsheet containing the mapping fields to consolidate the project datafiles (ecotaxa/ecopart depending on accessibility).
   :param sheetname: the name of the corresponding sheet of the standardizer (only works for ecotaxa).
   :param run_macro: boolean factor indicating whether the ImageJ macro to generate missing metadata should be run automatically. If False, return will be null until metadata have been created and uploaded to Ecotaxa manually.
   :param upload_metadata: boolean factor indicating whether the project metadata should be uploaded on Ecotaxa.
   :param localpath_metadata: local path of the project metadata storage (must have writing access, canno't be plankton server).
   :return: Extended dataframe ecotaxa_ecopart/ecotaxa_export_projectID_yyyymmdd_hhmm_extended.tsv
   """
   if project_id  in standardizer.index and str(standardizer.at[project_id, "External_project"])!='nan':
       path_to_standardizer=path_to_git / 'raw' / 'project_UVP_standardizer.xlsx'
       path_to_extended_file = Path(localpath_metadata).expanduser()
       path_to_extended_file.mkdir(parents=True, exist_ok=True)
       path_to_project = Path(standardizer.at[project_id, "Project_localpath"]).expanduser()
       path_to_project =path_to_project if path_to_project.stem!=path_to_extended_file.stem else Path(path_to_extended_file.parent).expanduser() / sheetname.capitalize()
       path_to_export = list(path_to_project.rglob('**/*_{}_*'.format(str(project_id))))[0]
       path_to_external_project = path_to_project.parent / {'ecotaxa':'ecopart','ecopart':'ecotaxa',path_to_extended_file.stem:"ecopart"}[path_to_project.name]
       path_external_export=list(itertools.chain.from_iterable([path_to_external_project.rglob('**/{}_export_*_{}_*.tsv'.format({'ecotaxa': 'ecopart', 'ecopart': 'ecotaxa',path_to_extended_file.stem:"ecopart"}[path_to_project.name].lower(),project)) for project in str(standardizer.at[project_id, "External_project"]).split(';')]))
       # Add small particles (area<=sm_zoo) to Ecotaxa/Ecopart and save to export_*_extended.tsv:
       if path_to_external_project.name=='ecopart':
           path_ecopart_par_export =[path for path in path_external_export if ('_particles.tsv' in path.name) & ('ecopart_export_raw' in path.name)]
           path_ecopart_zoo_export = [path for path in path_external_export if ('_zooplankton.tsv' in path.name) & ('ecopart_export_raw' in path.name)]
           path_ecopart_metadata = [path for path in path_external_export if ('_metadata.tsv' in path.name) & ('ecopart_export_raw' in path.name)]
           dict_metadata_columns={'depth':'object_corrected_min_depth_bin','imgcount':'object_imgcount_bin'}
           df_ecopart = pd.concat(map(lambda path:pd.read_table(path, encoding='latin-1').assign(Project_ID=int(re.findall(r'\d+',path.name)[0])).sort_values( ['profileid', 'depth']),path_ecopart_par_export)).rename(columns=dict_metadata_columns) if len(path_ecopart_par_export) else pd.DataFrame({})
           df_ecopart = df_ecopart.reset_index(drop=True)
           df_ecopart['object_corrected_max_depth_bin']=df_ecopart['object_corrected_min_depth_bin']+1
           df_metadata = pd.concat(map(lambda path:pd.read_table(path, encoding='latin-1').assign(Project_ID=int(re.findall(r'\d+',path.name)[0])),path_ecopart_metadata)) if len(path_ecopart_metadata) else pd.DataFrame({})
           df_metadata = df_metadata.reset_index(drop=True)
           df_metadata['Localpath_project']=df_metadata.rawfolder.apply(lambda path: Path(path).expanduser().stem)
           df_ecopart=pd.merge(df_ecopart,df_metadata[['profileid','acq_volimage']],how='left',on='profileid')
           df_ecopart['object_volume_bin']=df_ecopart.acq_volimage*df_ecopart.object_imgcount_bin
           if path_to_export.exists() and len(df_ecopart) and len(df_metadata):
              df_ecotaxa_columns=pd.read_table(path_to_export,nrows=0).columns
              df_ecotaxa=pd.read_table(path_to_export,usecols=[column for column in df_ecotaxa_columns if column not in ecotaxa_fields_to_drop])
              depth_offset= df_metadata[df_standardizer_ecopart[df_standardizer_ecopart.index.astype(str).isin(standardizer.at[project_id, "External_project"].split(';'))]["Depth_offset"].values[0]].values[0]
              sm_zoo=df_metadata.at[0,'acq_smzoo']
              nativefolder_dict=df_metadata[['Project_ID','Localpath_project']].drop_duplicates().to_dict('r')

              #1) Generate extended ecopart dataframe with all particles number/area using df_ecopart nbr (=number of ROI with corresponding area)
              df_ecopart_extended=df_ecopart[df_ecopart.area<=sm_zoo].groupby(by=['Project_ID','profileid','object_corrected_min_depth_bin','object_corrected_max_depth_bin','object_volume_bin','object_imgcount_bin'],observed=True).apply(lambda df:pd.DataFrame({'object_id':['particle_unknown_id_{}'.format(str(index)) for index in pd.Series(np.repeat(df.area.values,repeats=df.nbr.values)).index],'object_bru_id':['particle_unknown_id_{}'.format(str(index)) for index in pd.Series(np.repeat(df.area.values,repeats=df.nbr.values)).index],'object_bru_area':np.repeat(df.area.values,repeats=df.nbr.values)})).reset_index()
              df_ecopart_extended=df_ecopart_extended.drop(['level_6'],axis=1)
              #2) Merge df_ecopart_extended to df_ecotaxa based on object_areai (raw ecopart), area (plankton ecopart), or object_bru_area (updated ecotaxa)
              if all(pd.Series(['object_bru_id','object_bru_area','object_image_name','object_image_index','object_corrected_depth','object_corrected_min_depth_bin','object_corrected_max_depth_bin','object_volume_bin','object_imgcount_bin']).isin(list(df_ecotaxa.columns))):
                  df_ecopart_extended = pd.concat([pd.merge(df_ecopart_extended.assign(sample_id=df_ecopart_extended.profileid, object_annotation_category='', object_annotation_status='unclassified'), df_ecotaxa[ ['object_lat', 'object_lon', 'object_date', 'object_time'] + [column for column in df_ecotaxa.columns if len(re.findall( r'sample_|acq_|process_',column)) > 0]].drop_duplicates(), how='left', on='sample_id'), df_ecotaxa], axis=0)
                  df_ecopart_extended = df_ecopart_extended[df_ecotaxa.columns].sort_values(['sample_id', 'object_corrected_min_depth_bin', 'object_bru_area'], ascending=[True, True, False]).reset_index(drop=True)
                  df_ecopart_extended['acq_smbase'] = df_ecopart_extended['acq_smbase'] + 1
                  df_ecopart_extended.to_csv(path_to_extended_file / str(path_to_export.name).replace('.tsv', '_ecopart_{}_consolidated.tsv').format( '_'.join(path_to_projects.Project_ID.astype(str))),sep='\t', index=False)
                  print('Consolidated project created at {}'.format(path_to_extended_file / str( path_to_export.name).replace('.tsv', '_ecopart_{}_consolidated.tsv').format('_'.join(standardizer.at[project_id, "External_project"].split(';')))))
                  print( 'Updating standardizer spreadsheet: Project_localpath ({}), Depth_min (object_corrected_min_depth_bin), Depth_max (object_corrected_max_depth_bin), Area (object_bru_area), Volume_analyzed (object_volume_bin), Sampling_lower_size (acq_smbase) '.format(str(path_to_extended_file).replace(str(Path.home()), '~')))
                  index_columns_to_replace = [index for index, column in enumerate(pd.Series(['Project_localpath', 'Flag_path', 'Depth_min_field', 'Depth_max_field', 'Area_field','Minor_axis_field', 'Minor_axis_unit', 'Volume_analyzed_field', 'Sampling_lower_size_field','Sampling_lower_size_unit'])) if str(standardizer.loc[project_id, column]) != 'nan']
                  standardizer.loc[project_id, list(pd.Series(['Project_localpath', 'Flag_path', 'Depth_min_field', 'Depth_max_field', 'Area_field','Minor_axis_field', 'Minor_axis_unit', 'Volume_analyzed_field', 'Sampling_lower_size_field','Sampling_lower_size_unit'])[index_columns_to_replace])] = list(pd.Series([str(path_to_extended_file).replace(str(Path.home()), '~'), pd.NA,'object_corrected_min_depth_bin', 'object_corrected_max_depth_bin', 'object_bru_area', pd.NA,pd.NA, 'object_volume_bin', 'acq_smbase', 'pixel'])[index_columns_to_replace])

                  with pd.ExcelWriter(str(path_to_standardizer), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                      standardizer.rename_axis('Project_ID').reset_index().to_excel(writer, sheet_name=sheetname, index=False)

              else:
                  decision='Please run the ImageJ macro on your computer and update the metadata of your project by following the macro instructions' if not run_macro else 'Running the ImageJ macro to generate metadata automatically.'
                  print('Additional metadata (e.g. object_bru_area) not found in Ecotaxa export datafile.\n{}'.format(decision))
                  if run_macro:
                      path_to_metadata =Path(localpath_metadata).expanduser()/cfg['UVP_consolidation_metadata_subdir'] / "ecotaxa_{}_metadata".format(project_id)
                      path_to_metadata.mkdir(parents=True,exist_ok=True)
                      # Create readme file:
                      if not Path(path_to_metadata.parent.parent / 'README.txt').is_file():
                          with open(str(Path(path_to_metadata.parent.parent / 'README.txt')), 'w') as file:
                              file.write( "README file for UVP project consolidated files (First created on February 10, 2023):\n\nThis directory contains the consolidated files for Ecotaxa/Ecopart UVP projects with the format:ecotaxa_export_projectID(Ecotaxa)_Exportdate_Exporttime_ecopart_projectID(Ecopart).\nThe consolidation include the following steps:\n\n1) Ecotaxa default exported file and Ecopart raw particles exported file are loaded \n\n2) An extended Ecopart datatable is created using (1) the nbr/area columns to generate a row for each particle of the observed sizes, (2) the imgcount column to calculate the volume imaged in each 1-m depth bins.\nEach particle is assigned an id (object_id:particle_unknown_id_#), a null taxonomic annotation and the 'unclassified' status.\nAttention: this table contains small particles only (particles>sm_zoo are filtered out since they correspond to Ecotaxa vignettes)  \n\n3) Additional variables are appended to the native Ecotaxa table, including the object_corrected_depth (accounting for the spatial offset between the imaging frame and the depth sensor), object_min/max_depth_bin (1m-depth bins matching those of Ecopart), object_volume_bin (matching the cumulative volume in 1m-depth bins, in liters). The 'object_bru_area' variable is created (if run_macro=True) using a custom PyimageJ macro that reproduce the initial processing step of all particles segmentation, from the native project folder located in ~/plankton/uvpx_missions\nAttention: If the project folder does not include the required sub-directories (raw with all bmp images, config, and work directories), the object_bru_area will correspond to 'object_bru'  \n\n4) Both data tables are concatenated and saved in {}\n\nContact us: nmfs.pssdb@noaa.gov".format(str(path_to_metadata.parent.parent).replace(str(Path.home()),'~')))

                      #1) Step 1: Generate ecotaxa metadata using the ImageJ macro PyImageJ_ecotaxa_append_metadata.txt

                      if str(Path.home())=='/home/rkiko':
                          localpath_to_project = Path('/home/rkiko/plankton/uvp5_missions').expanduser()
                      elif (str(Path.home())=='/Users/dugennem') or (str(Path.home())=='/Users/mc4214'):
                          localpath_to_project = Path.home().expanduser()
                      else:
                          print("Consolidation only available on the server. Please run the script on Marie")
                          return

                      path_to_projects=df_metadata[['Project_ID','ptitle','Localpath_project']].drop_duplicates().reset_index().drop('index',axis=1)
                      for path in path_to_projects.Localpath_project:
                          arg = {'project_path':str(localpath_to_project/ path) + os.sep ,'depth_offset': str(depth_offset),'metadata_path':str(path_to_metadata)+ os.sep}
                          #arg = {'project_path': '/home/rkiko/plankton/uvp5_missions/uvp5_sn010_2014_m108/', 'depth_offset': '1.2','metadata_path':'~/GIT/PSSdb/raw/ecotaxa_ecopart_consolidation/ecotaxa_metadata/ecotaxa_586_metadata/'}
                          #arg = {'project_path': '/Volumes/upv5_missions/uvp5_sn000_tara2009/', 'depth_offset': '1.2','metadata_path':'~/GIT/PSSdb/raw/ecotaxa_ecopart_consolidation/ecotaxa_metadata/ecotaxa_586_metadata/'}

                          ## Note: The PyimageJ macro save results WITH row numbers (PyimageJ bug): tables with row number cannot be uploaded on Ecotaxa
                          ## I modified PyImageJ_ecotaxa_append_metadata so that it returns the table in string that can be imported and saved on Python
                          try:
                              macro=run_imageJ_macro(macro_path=path_to_macro, arguments=arg)
                              result = str(macro.getOutput('result'))
                          except Exception as e:
                              print('Error running imageJ macro on {}. Please check the folder, known issues include the absence of required file(s) in the project folder (i.e raw subdirectory with all bmp vignettes, config and work subdirectories)'.format(str(localpath_to_project/ path) + os.sep))
                              result='None'

                          ## Re-save metadata files without row numbers (pyImageJ bug)
                          if len(result) and (result!='None'):
                                table=pd.read_csv(io.StringIO(result,newline='\n'), sep=",")
                                table.columns=table.columns.str.strip()
                                #table=table.drop(table.columns[0],1) # Drop row numbers
                                table_headings=table.iloc[0].apply(lambda x:x.strip()).values
                                table=table.drop(0,0) # Drop Ecotaxa sub-headers used to describe the variable type
                                table=table[(table.apply(lambda x: x.isin(table_headings).any()==False, axis=1)) &  (table.apply(lambda x: x.isin(list(table.columns)).any()==False, axis=1))]
                                table['sample_id']=table.object_id.apply(lambda id: id[0:(id.rfind("_"))])
                                table.groupby(by='sample_id').apply(lambda x: pd.concat([pd.DataFrame(table_headings,index=table.columns[0:len(table.columns)-1]).T,x.drop('sample_id',axis=1)],ignore_index=True,axis=0).to_csv(path_to_metadata/"ecotaxa_{}_metadata_v1.tsv".format(x.sample_id.unique()[0]), sep="\t",index=False))
                                print('Saving project metadata to {}.'.format(str(path_to_metadata)))
                          else:
                                print('Consolidated file will be done with vignettes object_area')
                                df_ecotaxa['object_bru_area'] =df_ecotaxa['object_area']
                                df_ecotaxa['object_bru_id']=df_ecotaxa['object_id']
                                df_ecotaxa['object_corrected_depth'] = df_ecotaxa['object_depth_min']+depth_offset
                                df_ecotaxa['object_corrected_min_depth_bin'] =pd.cut(df_ecotaxa['object_corrected_depth'],bins=list(np.arange(0,6001,step=1)),right=False,labels=np.arange(0,6000,step=1)).astype(str).astype(int)
                                df_ecotaxa['object_corrected_max_depth_bin'] =df_ecotaxa['object_corrected_min_depth_bin'] +1
                                df_ecotaxa=pd.merge(df_ecotaxa,df_ecopart_extended[['profileid','object_corrected_min_depth_bin','object_volume_bin']].drop_duplicates(),how='left',left_on=['sample_id','object_corrected_min_depth_bin'],right_on=['profileid','object_corrected_min_depth_bin'])

                                df_ecopart_extended = pd.concat([pd.merge(df_ecopart_extended.drop(columns=['Project_ID']).assign(sample_id=df_ecopart_extended.profileid, object_annotation_category='',object_annotation_status='unclassified'), df_ecotaxa[ ['object_lat', 'object_lon', 'object_date', 'object_time'] + [column for column in df_ecotaxa.columns if len(re.findall( r'sample_|acq_|process_', column)) > 0]].drop_duplicates(),how='left', on='sample_id'), df_ecotaxa], axis=0)
                                df_ecopart_extended = df_ecopart_extended[df_ecotaxa.columns].sort_values(['sample_id', 'object_corrected_min_depth_bin', 'object_bru_area'],ascending=[True, True, False]).reset_index(drop=True)
                                df_ecopart_extended['acq_smbase'] = df_ecopart_extended['acq_smbase'] + 1
                                df_ecopart_extended.to_csv(path_to_extended_file / str(path_to_export.name).replace('.tsv', '_ecopart_{}_consolidated.tsv').format('_'.join(path_to_projects.Project_ID.astype(str))), sep="\t", index=False)
                                print('Consolidated project created at {}'.format(path_to_extended_file / str(path_to_export.name).replace('.tsv','_ecopart_{}_consolidated.tsv').format('_'.join(path_to_projects.Project_ID.astype(str)))))
                                # Update standardizer Depth_min (object_corrected_min_depth_bin), Depth_max (object_corrected_max_depth_bin), Area (object_bru_area), Profile (sample_id), Sample (object_corrected_min_depth_bin), Volume_analyzed (object_volume_bin)
                                print('Updating standardizer spreadsheet: Project_localpath ({}), Depth_min (object_corrected_min_depth_bin), Depth_max (object_corrected_max_depth_bin), Area (object_bru_area), Volume_analyzed (object_volume_bin), Sampling_lower_size (acq_smbase) '.format( str(path_to_extended_file).replace(str(Path.home()), '~')))
                                index_columns_to_replace = [index for index,column in enumerate(pd.Series(['Project_localpath', 'Flag_path', 'Depth_min_field','Depth_max_field', 'Area_field', 'Minor_axis_field','Minor_axis_unit', 'Volume_analyzed_field','Sampling_lower_size_field', 'Sampling_lower_size_unit'])) if str(standardizer.loc[project_id, column]) != 'nan']
                                standardizer.loc[ project_id, list(pd.Series(['Project_localpath', 'Flag_path', 'Depth_min_field', 'Depth_max_field', 'Area_field', 'Minor_axis_field', 'Minor_axis_unit', 'Volume_analyzed_field', 'Sampling_lower_size_field', 'Sampling_lower_size_unit'])[index_columns_to_replace])] = list(pd.Series([str(path_to_extended_file).replace(str(Path.home()), '~'), pd.NA,'object_corrected_min_depth_bin', 'object_corrected_max_depth_bin', 'object_bru_area', pd.NA, pd.NA, 'object_volume_bin', 'acq_smbase', 'pixel'])[index_columns_to_replace])
                                with pd.ExcelWriter(str(path_to_standardizer), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                                    standardizer.rename_axis('Project_ID').reset_index().to_excel(writer, sheet_name=sheetname,index=False)
                      if len(list(path_to_metadata.glob("*.tsv"))):
                          # 2) Step 2: Zip all ecotaxa_profileID_metadata files and merge to existing Ecotaxa project using the API
                          try: # Fast compression
                                os.system('zip -r %s %s' % (str(path_to_metadata / (path_to_metadata.stem + '.zip')), str(path_to_metadata) + os.sep))
                          except Exception as e: # Slow compression
                                shutil.make_archive(base_name=str(path_to_metadata/ 'ecotaxa_{}_metadata'.format(project_id)), format='zip',base_dir=path_to_metadata)

                          # 3) Update Ecotaxa project using the API
                          if upload_metadata:
                              print('Updating Ecotaxa project with metadata zipfile {}. Please wait'.format(str(path_to_metadata/ 'ecotaxa_{}_metadata.zip'.format(project_id))))
                              update_ecotaxa_project_API(project_id, path_to_metadata=path_to_metadata/ 'ecotaxa_{}_metadata.zip'.format(project_id))

                          # 4) Merge df_ecotaxa to all metadata_v1_tables
                          path_to_ecotaxa_metadata=natsorted(list(path_to_metadata.parent.rglob('*_metadata_v1.tsv')))
                          df_ecotaxa_metadata=pd.concat(map(lambda path:pd.read_table(path,skiprows=0),path_to_ecotaxa_metadata)).drop(index=0)
                          df_ecotaxa_metadata=df_ecotaxa_metadata.astype(dict(zip(['object_id','object_bru_id','object_bru_area','object_image_name','object_image_index','object_corrected_depth','object_corrected_min_depth_bin','object_corrected_max_depth_bin','object_volume_bin','object_imgcount_bin'],[str,str,float,str,int,float,float,float,float,int])))
                          df_ecotaxa_metadata =df_ecotaxa_metadata.sort_values(['object_bru_id'])
                          df_ecotaxa=pd.merge(df_ecotaxa,df_ecotaxa_metadata,how="left",on=['object_id'])
                          #df_ecotaxa = df_ecotaxa[['object_id', 'object_bru_id', 'object_area','object_bru_area','profileid','sample_corrected_min_depth_bin','sample_corrected_max_depth_bin','sample_volume_bin','sample_imgcount_bin']]
                          df_ecotaxa = df_ecotaxa.sort_values(['sample_id', 'object_depth_min','object_id'])
                          ## Append sample, process, and acquistion fields to df_ecopart_extended and concatenate to ecotaxa
                          df_ecopart_extended = pd.concat([pd.merge(df_ecopart_extended.assign(sample_id=df_ecopart_extended.profileid,object_annotation_category='',object_annotation_status='unclassified'),df_ecotaxa[['object_lat','object_lon','object_date','object_time']+[column for column in df_ecotaxa.columns if len(re.findall( r'sample_|acq_|process_', column))>0]].drop_duplicates(),how='left',on='sample_id'), df_ecotaxa], axis=0)
                          df_ecopart_extended =df_ecopart_extended [df_ecotaxa.columns].sort_values(['sample_id', 'object_corrected_min_depth_bin','object_bru_area'],ascending=[True,True,False]).reset_index(drop=True)
                          df_ecopart_extended['acq_smbase']=df_ecopart_extended['acq_smbase']+1
                          df_ecopart_extended.to_csv(path_to_extended_file / str(path_to_export.name).replace('.tsv', '_ecopart_{}_consolidated.tsv').format('_'.join(path_to_projects.Project_ID.astype(str))),sep="\t" ,index=False)
                          print('Consolidated project created at {}'.format(path_to_extended_file / str(path_to_export.name).replace('.tsv', '_ecopart_{}_consolidated.tsv').format('_'.join(path_to_projects.Project_ID.astype(str)))))
                          # 5) Update standardizer Depth_min (object_corrected_min_depth_bin), Depth_max (object_corrected_max_depth_bin), Area (object_bru_area), Profile (sample_id), Sample (object_corrected_min_depth_bin), Volume_analyzed (object_volume_bin)
                          print('Updating standardizer spreadsheet: Project_localpath ({}), Depth_min (object_corrected_min_depth_bin), Depth_max (object_corrected_max_depth_bin), Area (object_bru_area), Volume_analyzed (object_volume_bin), Sampling_lower_size (acq_smbase) '.format( str(path_to_extended_file).replace(str(Path.home()), '~')))
                          index_columns_to_replace = [index for index, column in enumerate(pd.Series(['Project_localpath', 'Flag_path', 'Depth_min_field', 'Depth_max_field', 'Area_field','Minor_axis_field', 'Minor_axis_unit', 'Volume_analyzed_field', 'Sampling_lower_size_field','Sampling_lower_size_unit'])) if str(standardizer.loc[project_id, column]) != 'nan']
                          standardizer.loc[project_id, list(pd.Series(['Project_localpath', 'Flag_path', 'Depth_min_field', 'Depth_max_field', 'Area_field','Minor_axis_field', 'Minor_axis_unit', 'Volume_analyzed_field', 'Sampling_lower_size_field','Sampling_lower_size_unit'])[index_columns_to_replace])] = list(pd.Series([str(path_to_extended_file).replace(str(Path.home()), '~'), pd.NA,'object_corrected_min_depth_bin', 'object_corrected_max_depth_bin', 'object_bru_area', pd.NA,pd.NA, 'object_volume_bin', 'acq_smbase', 'pixel'])[index_columns_to_replace])
                          with pd.ExcelWriter(str(path_to_standardizer), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                              standardizer.rename_axis('Project_ID').reset_index().to_excel(writer, sheet_name=sheetname, index=False)
                      else:
                          print('Consolidation done')
                          return

                  else:
                      print('Please run the ImageJ macro ImageJ_ecotaxa_append_metadata and merge metadata files to your existing Ecotaxa project before continuing.\nTo merge the metadata files (project folder/ecotaxa_metadata):' '\n1) Please compress the files in a single .zip file\n2) Upload the zip file on Ecotaxa (https://ecotaxa.obs-vlfr.fr/prj/{}):\n   To add metadata, go to your existing project:\n   Project>Import images and metadata>Start re-import TSV files to update metadata and data>Upload tsv file(s) compressed in a zip file.\nNew metadata will be mapped and added to existing data based on object_id (Safe for annotations)  '.format(project_id))
                      return
           else:
                print("Ecotaxa/Ecopart export files not found. Run script 1.export_projects to get exported files")
                #os.system("python Scripts/1.export_projects.py")
                return
       else:
           print('Current consolidation only works for Ecotaxa projects')
           return
   else:
       print('Project not found in {} or corresponding Ecopart project not found. Skipping project, please check on Ecotaxa/Ecopart',format(str(path_to_standardizer)))
       return

def intersect_size_conversion(aa,exp,pixel):
    """
       Objective: This function returns the intersect (as equivalent circular diameter, micrometer) of UVP size conversions\n
       Reference: Picheral et al. (2010): The Underwater Vision Profiler 5: An advanced instrument for high spatial resolution studies of particle size spectra and zooplankton.  Limnol. Oceanogr.: Methods 8, 2010, 462–473\n
       :param aa: scalar coefficient derived from the intersect of the log-log regression between ROI surface area (square pixels) and microscopic measurements (square millimeters). (square_millimeters_per_square_pixels)
       :param exp: exponent coefficient derived from the slope of the log-log regression between ROI surface area (square pixels) and microscopic measurements (square millimeters). (unitless)
       :param pixel: fixed pixel size coefficient (millimeter_per_pixel)
       :return: equivalent circular diameter (micrometer)
    """
    intercept=(pixel * 2e+03 * ((np.exp(np.log((pixel ** 2) / aa) / (exp - 1))) / np.pi) ** 0.5)
    return intercept


def read_bru(path_bru):
   if path_bru.exists():
      with open(path_bru) as input_text:
         headers=input_text.read().split('\n')[0]
      df_bru=pd.read_csv( path_bru,encoding='latin_1', sep=r'\;|\t',header=0 if headers.split(";")[0]=='index' else None)
      df_bru=df_bru[df_bru.columns[df_bru.apply(lambda col:all(col.isna()),axis=0)==False]]#[[column for column in df_bru.columns if 'Unnamed' not in column]]
      df_bru.columns=['image_index','image','blob','area','mean_grey','x_center','y_center'] if headers.split(";")[0]=='index' else bru_headers
      return df_bru

def read_dat(path_dat):
   if path_dat.exists():
      df_dat=pd.read_table(path_dat,sep=r'\;|\t',encoding="latin-1",header=None,index_col=0)
      df_dat.insert(0, 'image_index', df_dat.index)
      df_dat = df_dat.loc[:, df_dat.apply(lambda x: all(pd.isna(x)), axis=0).index[df_dat.apply(lambda x: all(pd.isna(x)),axis=0) == False]]  # df_dat = df_dat[[column for column in df_dat.columns if 'Unnamed' not in column]].dropna().reset_index(drop=True)
      df_dat.columns = dat_headers
      df_dat = df_dat[['image_index', 'image', 'pressure', 'nb_blobs_SMbase', 'nb_blobs_SMzoo']]
   return df_dat

def read_pid(path_pid):
   pid= {}
   if path_pid.exists():
      # Read pid file
      with open(path_pid) as input_text:
         df_pid = pd.Series(input_text.read().split('\n'))
         df_pid_config = df_pid[df_pid.str.contains("=")]
         df_pid_config=  pd.concat(map(lambda x: pd.DataFrame({x[0]:x[1].replace('\t','')},index=[Path(path_pid).parent.stem]),df_pid_config.str.split(r'= ').tolist()),axis=1)
         pid['metadata_process']=df_pid_config
         df_pid_data= pd.DataFrame(df_pid.iloc[(df_pid.where(df_pid.str.contains("[Data]",regex=False)).first_valid_index()+2):-1].apply(lambda row: row.split(';')).tolist(),columns=df_pid.iloc[(df_pid.where(df_pid.str.contains("[Data]",regex=False)).first_valid_index()+1)].split(';'))
         pid['dat'] = df_pid_data
   return pid

def profile_first_image(project_path,portal='Ecopart'):
   """
   This function returns the first UVP image used to process raw datafiles on Ecopart (firstimage in project / meta / uvp5_header_project.txt) or Ecotaxa (derived from First_image in project / work/ profile / profiledat1.pid)
   Ecotaxa (stable depth allowed during descent)
       :param project_path: path of the raw project folder.
       :param portal: Ecotaxa or Ecopart.
       :return: Dataframe with first image index/name used to process UVP raw datafiles on Ecotaxa or Ecopart per profile
   """
   df=pd.DataFrame({})
   if portal=='Ecopart':
      path_meta=Path(project_path / "meta" /"{}_header{}.txt".format(project_path.stem.split("_")[0],project_path.stem.replace(project_path.stem.split("_")[0],'')))
      df_meta=pd.read_csv(path_meta,encoding='latin-1',sep=';',index_col='profileid')
      df=pd.DataFrame({'first_image':df_meta.loc[:, 'firstimage'].values},index=df_meta.index)
   if portal=='Ecotaxa':
       for profile in next(os.walk(Path(project_path) / "work"))[1]:
            # Set path to configurations, dat, and bru files
            ## Attention: first image used to generate vignettes UVP profiles may differ from the indices recorded in project_path/meta/uvp5_header_project.txt
            path_pid = list(Path(project_path / "work" / profile).glob('*dat1.pid'))[0]  # Use pid table for processed first images
            df_pid = read_pid(path_pid)['metadata_process']
            first_image = int(df_pid.at[profile, 'First_image'])
            first_image0 = first_image
            path_dat=path_dat=list(Path(project_path / "work" / profile).glob('*_datfile.txt'))[0]
            df_dat = read_dat(path_dat)
            first_image_date = (df_dat.at[int(first_image), 'image'])[0:13]  ## Attention: the date used in Zooprocess is truncated at the 10 second
            firstimage = df_dat.loc[(df_dat['image'].apply(lambda date: date[0:13]).values == first_image_date), 'image_index'].index[0]
            df=pd.concat([df,pd.DataFrame({'first_image':firstimage},index=[profile])],axis=0)
   return df


def profile_last_image(image_series,depth_series,portal='Ecopart',first_image=None,depth_min=0):
   """
   This function returns the last UVP image used to process raw datafiles on Ecopart (stable depth not allowed during descent) or
   Ecotaxa (stable depth allowed during descent)
       :param image_series: series of images index/name.
       :param depth_series: series of images depth (meters).
       :param portal: Ecotaxa or Ecopart.
       :param first_image: first image index/name used to filter images. If None, first_image is extracted from the image_series index 0
       :param depth_min: minimum depth allowed (meters). Default is 0.5 for portal Ecopart, 0 for portal Ecotaxa
       :return: Last image index/name used to process UVP raw datafiles on Ecotaxa or Ecopart
   """
   index_images=pd.Series(image_series)
   depth_series= list(depth_series) if not isinstance(depth_series,list) else depth_series
   image_series = list(image_series) if not isinstance(image_series,list) else image_series
   first_image=image_series[0] if not first_image else first_image
   last_image =image_series[len(image_series)-1]
   depth_series = list(compress(depth_series, ((index_images.index>=index_images.index[index_images.values==first_image].values[0]) & (index_images.index<=index_images.index[index_images.values==last_image].values[0]))))
   image_series=list(compress(image_series, ((index_images.index>=index_images.index[index_images.values==first_image].values[0]) & (index_images.index<=index_images.index[index_images.values==last_image].values[0]))))

   depth_min=depth_min if not (depth_min==0) | (depth_min==0.5) else 0 if portal=='Ecopart' else 0.5
   prev_depth_image = depth_min
   condition='(depth >= prev_depth_image)' if portal=="Ecotaxa" else '(depth > prev_depth_image)'
   for index,image in enumerate(image_series):
      depth =depth_series[index] if depth_series[index] > 0 else 0
      if eval(condition) & (depth >= depth_min): # if (depth > prev_depth_image) & (depth >= depth_min): Attention the filter to extract the last image for Ecopart does not allow for stable depth (see l.481 of function GenerateRawHistogram https://github.com/ecotaxa/ecopart/blob/master/py/part_app/funcs/uvp_sample_import.py)
         prev_depth_image = depth
         lastimage = image
   return lastimage


def descending_images(image_series,depth_series,portal='Ecopart',first_image=None,last_image=None,depth_min=0):
   """
   This function filters out UVP images saved while the rosette was ascending from a series of images index/name and associated depth
       :param image_series: series of images index/name.
       :param depth_series: series of images depth (meters).
       :param portal: Ecotaxa or Ecopart.
       :param first_image: first image index/name used to filter images. If None, first_image is extracted from the image_series index 0
       :param last_image: last image index/name used to filter images. If None, last_image corresponds to the deepest image index/name
       :param depth_min: minimum depth allowed (meters). Default is 0.5 for portal Ecopart, 0 for portal Ecotaxa
       :return: list of images saved during descent
   """
   index_images=pd.Series(image_series)
   depth_series= list(depth_series) if not isinstance(depth_series,list) else depth_series
   image_series = list(image_series) if not isinstance(image_series,list) else image_series
   first_image=image_series[0] if not first_image else first_image
   last_image =image_series[len(image_series)-1] if not last_image else last_image
   depth_series = list(compress(depth_series, ((index_images.index>=index_images.index[index_images.values==first_image].values[0]) & (index_images.index<=index_images.index[index_images.values==last_image].values[0]))))
   image_series=list(compress(image_series, ((index_images.index>=index_images.index[index_images.values==first_image].values[0]) & (index_images.index<=index_images.index[index_images.values==last_image].values[0]))))
   images_in_profile = []
   depth_min=depth_min if not (depth_min==0) | (depth_min==0.5) else 0 if portal=='Ecopart' else 0.5
   prev_depth_image = depth_min
   for index,image in enumerate(image_series):
      depth =depth_series[index] if depth_series[index] > 0 else 0
      if (depth >= prev_depth_image) & (depth >= depth_min): # if (depth > prev_depth_image) & (depth >= depth_min): Attention the filter to extract the last image for Ecopart does not allow for stable depth (see l.481 of function GenerateRawHistogram https://github.com/ecotaxa/ecopart/blob/master/py/part_app/funcs/uvp_sample_import.py)
         images_in_profile = images_in_profile+[image_series[index]]
         prev_depth_image = depth
   return images_in_profile


def UVP_vignettes_matchup(project_path,pressure_gain=10,depth_offset=1.2):
   """
   This function matches the blob index from UVP bru files to each vignette (with unique id: imagename_blobindex) and return dataframes for vignettes and images selected on the descent following previous versions of Zooprocess (up to 8.12)
         :param project_path: path of the UVP project.
         :param pressure_gain: gain to correct the pressure sensor reading from centibars to decibars (meters). Default is 10.
         :param depth_offset: spatial offset between the imaging field and the pressure sensor (meters). Default is 1.2 meters.
         :return: dictionary of dataframes vignettes/images
   """
   df_images=pd.DataFrame({})
   df_vignettes=pd.DataFrame({})
   df_process=pd.read_table(Path(project_path / "config" / "process_install_config.txt" ),encoding='latin_1', sep=';',header=None)
   depth_min=float(df_process[0].loc[df_process[0].str.contains("profmin=")].values[0].split("= ")[1])
   for profile in next(os.walk(Path(project_path) / "work"))[1]:
      # Set path to configurations, dat, and bru files
      ## Attention: first image used to generate vignettes UVP profiles may differ from the indices recorded in project_path/meta/uvp5_header_project.txt
      path_config=list(Path(project_path / "work" / profile).glob('*dat1.pid'))[0] # Use pid table for processed first images
      path_dat=list(Path(project_path / "work" / profile).glob('*_datfile.txt'))[0]
      #path_bru = Path(project_path) / "work"/ profile / (profile+'.bru') if  Path(project_path / "work"/ profile / (profile+'.bru')).exists() else list(Path(project_path / "work" / profile).glob('HDR*.bru'))[0]
      if path_config.exists() & path_dat.exists():
         # Read config file to get first image index
         df_config = read_pid(path_config)['metadata_process']
         first_image=int(df_config.at[profile,'First_image'])
         first_image0=first_image
         df_dat=read_dat(path_dat)
         first_image_date=(df_dat.at[int(first_image),'image'])[0:13] ## Attention: the date used in Zooprocess is truncated at the 10 second
         first_image=df_dat.loc[(df_dat['image'].apply(lambda date: date[0:13]).values==first_image_date),'image_index'].index[0]
         # Filter out ascending depths
         images_in_profile={}
         vignettes=pd.DataFrame({})
         prev_depth =0#df_dat.at[int(first_image),"pressure"]/pressure_gain if (df_dat.at[int(first_image),"pressure"]/pressure_gain)>0 else 0
         prev_depth_image=prev_depth
         for i in df_dat.image_index.loc[df_dat.image_index>=first_image].index:
             depth=df_dat.at[i,"pressure"]/pressure_gain if (df_dat.at[i,"pressure"]/pressure_gain)>0 else 0
             if depth>=prev_depth_image:
                if depth>=depth_min:
                    images_in_profile[i] = depth
                    prev_depth_image=depth
                    last_image = i
             if depth>=prev_depth:
                if (df_dat.at[i,"nb_blobs_SMzoo"]>0) & (depth>=depth_min):
                    vignettes=pd.concat([vignettes,pd.DataFrame({'Profile':profile,'Label_ecotaxa':[profile+'_'+str(len(vignettes)+blob) for blob in np.arange(1,df_dat.at[i,"nb_blobs_SMzoo"]+1,step=1)],'Label_bru':[df_dat.at[i, "image"]+'_'+str(blob).zfill(4) for blob in np.arange(0,df_dat.at[i,"nb_blobs_SMzoo"],step=1)],'image_index':df_dat.at[i, "image_index"],'image':df_dat.at[i,"image"],'depth':depth,'depth_corrected':depth+depth_offset},index=[df_dat.at[i, "image"]+'_'+str(blob).zfill(4) for blob in np.arange(0,df_dat.at[i,"nb_blobs_SMzoo"],step=1)])],axis=0)
                    prev_depth=depth # Update previous depth
         df_images=pd.concat([df_images,pd.DataFrame({'Profile': profile, 'image_depth': images_in_profile.values(), 'image_index': images_in_profile.keys(),'image':(df_dat.query('image_index.isin({})'.format(list(images_in_profile.keys()))).image.values)})])
         df_vignettes =pd.concat([df_vignettes,vignettes])

      else:
          print("File(s) {} not found for profile {}".format(list(compress([str(path_config),str(path_dat),str(path_bru)], [not path_config.exists(),not path_dat.exists(),not path_bru.exists()])),profile))
          continue
   return {'vignettes':df_vignettes.reset_index(drop=True),'images':df_images.reset_index(drop=True)}
