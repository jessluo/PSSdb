## Objective: Exporting Ecotaxa project tsv file using Ecotaxa API in batch

## Requirements:
# Note that you need annotator right on the project you wish to download
# Note that you need to install ecotaxa API using git first:
# in terminal: pip install git+https://github.com/ecotaxa/ecotaxa_py_client.git
# Note that you need two text files in your git repository that contain the info relative to:
# authentication (Ecotaxa_API_pw.yaml, in .gitignore), project ID/output path (Ecotaxa_API.yaml, tracked)
# Note that exported zipfiles may be large (>50 MB), hence set the output path accordingly (info stored in Ecotaxa_API.yaml)
# Note that you need to change the path of the config files on lines 72 & 76

## Workflow:
# This script uses 4 API steps:
# Step 1: Configure the API using your Ecotaxa login/pw (info stored in untracked Ecotaxa_API_pw.yaml)
# Step 2: Generate a request to export a project tsv based on project ID (info stored in git-tracked Ecotaxa_API.yaml)
# Step 3: Generate a job/task(=export) based on the request (step 2)
# Step 4: Get the job file (ecotaxa_export_projectID_date_time.zip)

## Useful documentations:
# See documentation for step 1: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/AuthentificationApi.md
# See documentation for step 2: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ExportReq.md
# See documentation for step 3: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ObjectsApi.md#export_object_set
# See documentation for step 4: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/JobsApi.md#get_job_file

## TO DO: Export projects from additional data sources (e.g. IFCB dashboards)


## Python modules

# Path modules
from pathlib import Path # Handling of path object
import shutil # Delete uncompressed export zip folder
import re

# Config modules
import yaml # requires installation of PyYAML package

# Progress bar modules
#import progressbar # Attention: Use pip install progressbar2 to install
from tqdm import tqdm
import datetime, time # Time module for break and current date

# Prompt for export confirmation
import sys

# Ecotaxa API modules. Install module via git first
import ecotaxa_py_client
### Authentification on Ecotaxa
from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.model.login_req import LoginReq
### Object, export, task API on Ecotaxa
from ecotaxa_py_client.api import objects_api
from ecotaxa_py_client.model.http_validation_error import HTTPValidationError
from ecotaxa_py_client.model.body_export_object_set_object_set_export_post import BodyExportObjectSetObjectSetExportPost
from ecotaxa_py_client.model.export_rsp import ExportRsp
from ecotaxa_py_client.model.export_req import ExportReq
from ecotaxa_py_client.api import jobs_api
from ecotaxa_py_client.model.job_model import JobModel
from ecotaxa_py_client.api import projects_api
from ecotaxa_py_client.model.project_filters import ProjectFilters
from ecotaxa_py_client.api import instruments_api
from sys import argv
import requests

# Dataframe modules
import pandas as pd
import numpy as np

# Workflow starts here

# read git-tracked config file (text file) with inputs:  project ID, output directory
path_to_config=Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
# read config file (text file) with inputs:  EcoTaxa username, password. Make sure to secure this file if needed
path_to_config_usr=Path('~/GIT/PSSdb/scripts/Ecotaxa_API_pw.yaml').expanduser()
with open(path_to_config_usr ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

path_to_git=Path(cfg['git_dir']).expanduser()

# Step 1: Authentication in EcoTaxa
with ecotaxa_py_client.ApiClient() as client:
    api = authentification_api.AuthentificationApi(client)
    token = api.login(LoginReq(
                               username=cfg_pw['ecotaxa_user'],
                               password=cfg_pw['ecotaxa_pass']
                               ))

configuration = ecotaxa_py_client.Configuration(host = "https://ecotaxa.obs-vlfr.fr/api",access_token=token, discard_unknown_keys=True)

# Step 1b: Authentication for JSON
BASE_URL = "https://ecotaxa.obs-vlfr.fr"
with requests.Session() as sess:
    rsp = sess.post(BASE_URL + "/api/login", json={"username": cfg_pw['ecotaxa_user'],
                                                   "password": cfg_pw['ecotaxa_pass']})
    json_token =rsp.json()


# prepare storage based on project list  stored in the yaml config file and instrument type
path_to_data = Path(cfg['git_dir']).expanduser() / cfg['dataset_subdir']
project_list=pd.read_excel(path_to_data.parent / Path(cfg['proj_list']).name,sheet_name="Data",usecols=['Project_ID','Instrument','PSSdb_access','Project_test'])
project_ids = project_list[project_list['PSSdb_access']==True].Project_ID.unique().astype(str)#str(cfg['proj_id'])
project_inst=project_list[project_list['PSSdb_access']==True].Instrument.unique().astype(str)

dict_instruments={'IFCB':'IFCB','UVP':['UVP5HD','UVP5SD','UVP5Z','UVP6'],'Zooscan':'Zooscan','Unknown':'?','AMNIS':'AMNIS','CPICS':'CPICS','CytoSense':'CytoSense','FastCam':'FastCam','FlowCam':'FlowCam','ISIIS':'ISIIS','LISST':['LISST','LISST-Holo'],'Loki':'Loki','Other':['Other camera','Other flowcytometer','Other microscope','Other scanner'],'PlanktoScope':'PlanktoScope','VPR':'VPR','ZooCam':'ZooCam','eHFCM':'eHFCM'}

#Prompting a warning to export all accessible projects or test set only
test_confirmation=input("Do you wish to export all accessible project(s) or the test subset? Enter all or test\n")
if test_confirmation=='test':
    project_ids=project_list[(project_list['PSSdb_access']==True) & (project_list['Project_test']==True)].Project_ID.unique().astype(str)


# Prompting a warning if export files for the list of projects already exist
# asking for confirmation to overwrite. Attention: all outdated project export files will be erased!
existing_project_path=list(path_to_data.rglob('*export_*'))
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


for i in range(len(project_ids)):
    proj_id = int(project_ids[i])
    path_to_export = path_to_data / [key for key, value in dict_instruments.items() if project_list[(project_list['PSSdb_access']==True) & (project_list['Project_ID']==proj_id)].Instrument.unique().astype(str)[0] in value][0]
    path_to_export.mkdir(parents=True, exist_ok=True)

    # Enter a context with an instance of the API client
    with ecotaxa_py_client.ApiClient(configuration) as api_client:
        # Create an instance of the API class
        api_instance = objects_api.ObjectsApi(api_client)
        # Step 2: Create a request to export project tsv file
        body_export_object_set_object_set_export_post = BodyExportObjectSetObjectSetExportPost(
            filters=ProjectFilters(),#ProjectFilters(statusfilter="PVD"),
            # get P(redicted), V(alidated), D(ubious) images. Check other options for filter here: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectFilters.md
            request=ExportReq(
                project_id=proj_id,  # the unique project ID of interest (integer)
                exp_type="TSV",
                use_latin1=False,
                tsv_entities="OPASH",
                # entities to be exported: O(bjects), P(rocess), A(cquisition), S(ample), classification H(istory)
                split_by="",  # no split=single table with all images/objects
                coma_as_separator=False,  # set decimal separator to point
                format_dates_times=False,
                with_images=False,  # exporting images
                with_internal_ids=False,
                only_first_image=False,
                sum_subtotal="A",
                out_to_ftp=False)
        )

        # Step 3: Generating a job/task (=Export Object Set)
        try:
            api_jobresponse = api_instance.export_object_set(body_export_object_set_object_set_export_post)

        except ecotaxa_py_client.ApiException as e:
            print("Exception when calling ObjectsApi->export_object_set: %s\n" % e)

        # Report job ID. You may check the job created here: https://ecotaxa.obs-vlfr.fr/Jobs/listall
        job_id = api_jobresponse['job_id']
        api_jobinstance = jobs_api.JobsApi(api_client)
        print("Creating export file with job ID: ", job_id, sep=' ')

        # Necessary break between step 3 and 4
        # Insert a progress bar to allow for the job to be done based on get_job status.
        # Attention, the break cannot be timed with job progress=(percentage) due to a small temporal offset between job progress and status

        with  tqdm(desc='Working on Job', total=1000, bar_format='{desc}{bar}', position=0, leave=True) as bar:
            job_status = 'R'  # percent = 0
            while job_status not in ('F', 'E'):  # percent!=100:
                time.sleep(2)  # Check the job status every 2 seconds. Modify as needed
                thread = api_jobinstance.get_job(job_id, async_req=True)
                result = thread.get()
                job_status = result.state
                percent = result.progress_pct
                bar.set_description("Working on Job %s%%" % percent, refresh=True)
                # and update progress bar
                ok = bar.update(n=10)
                if job_status == 'F':
                    break
                if job_status == 'E':
                    print("Error creating job. Please check your connection or EcoTaxa server")

    # Lines 195:201 are not working. Laurent is updating get_job_file.py.
    # Uncomment line 197 when get_job_file.py is working to proceed to step 4
    # api_downloadresponse = api_jobinstance.get_job_file(job_id)
    # api_downloadresponse = api_jobinstance.get_job_file(job_id, async_req=True, _preload_content=False,_return_http_data_only=False)
    # urlresponse = api_downloadresponse.get()
    # zip_file = urlresponse.headers['content-disposition'][22:-1]
    # url = urlresponse.geturl()

    # Step 4: Get/Download export zipfile
    zip_file = "ecotaxa_export_{}_{}Z.zip".format(str(proj_id), datetime.datetime.utcnow().strftime(
        "%Y%m%d_%H%M"))  # "ecotaxa_{}".format(str(result.result['out_file']))
    path_to_zip = path_to_export / zip_file
    path_to_log=path_to_export / "job_{}.log".format(str(job_id))
    print("\nExporting file", zip_file, "to", path_to_export, ", please wait", sep=' ')

    with requests.Session() as sess:
        url = BASE_URL + "/api/jobs/%d/file" % job_id
        rsp = sess.get(url, headers={"Authorization": "Bearer " + json_token}, stream=True)
        with open(path_to_zip, "wb") as fd:
            for a_chunk in rsp.iter_content():  # Loop over content, i.e. eventual HTTP chunks
                # rsp.raise_for_status()
                fd.write(a_chunk)


    print("Download completed: ", path_to_zip,"\nUnpacking zip file", sep='')
    shutil.unpack_archive(path_to_zip, path_to_export) # Unzip export file
    path_to_zip.unlink(missing_ok=True) # Delete zip file
    path_to_log.unlink(missing_ok=True) # Delete job log

quit()