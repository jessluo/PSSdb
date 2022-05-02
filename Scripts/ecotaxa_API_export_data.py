## Objective: Exporting Ecotaxa project tsv file using Ecotaxa API

## Requirements:
# Note that you need annotator right on the project you wish to download
# Note that you need to install ecotaxa API using git first:
# in terminal: pip install git+https://github.com/ecotaxa/ecotaxa_py_client.git
# Note that you need to change the local path to exported zipfile (l. 173)
# Note that you need two text files in your git repository that contain the info relative to:
# authentication (Ecotaxa_API_pw.yaml, in .gitignore), project ID/output path (Ecotaxa_API.yaml, tracked)
# Note that you need to change the path to the config files on lines 65 & 69

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

## TO DO: Explore ways to make the script faster, Loop through all project ID of interest


## Python modules

# Path modules
from pathlib import Path # Handling of path object

# Config modules
import yaml # requires installation of PyYAML package

# Progress bar modules
import progressbar # Creating progress bar of unknown length
import datetime

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

# Printing module
from pprint import pprint

## Workflow starts here

# read git-tracked config file (text file) with inputs:  project ID, output directory
path_to_config=Path('~/GIT/PSSdb/Scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
# read config file (text file) with inputs:  project ID, output directory
path_to_config_usr=Path('~/GIT/PSSdb/Scripts/Ecotaxa_API_pw.yaml').expanduser()
with open(path_to_config_usr ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

path_to_git=Path(cfg['git_dir']).expanduser()


# prepare storage based on path stored in the yaml config file and instrument type

path_to_data = Path(cfg['git_dir']).expanduser() / cfg['dataset_subdir']


# Query list of instruments
with ecotaxa_py_client.ApiClient(ecotaxa_py_client.Configuration(host = "https://ecotaxa.obs-vlfr.fr/api")) as api_client:
    # Create an instance of the API class
    api_instrumentinstance = instruments_api.InstrumentsApi(api_client)
    project_ids =  str(cfg['proj_id'])

    # example passing only required values which don't have defaults set
    try:
        # Instrument Query
        api_instrumentresponse = api_instrumentinstance.instrument_query(project_ids)

    except ecotaxa_py_client.ApiException as e:
        print("Exception when calling InstrumentsApi->instrument_query: %s\n" % e)


path_to_export = path_to_data / api_instrumentresponse[0]
#path_to_export.mkdir(parents=True,exist_ok=True)

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


# Enter a context with an instance of the API client
with ecotaxa_py_client.ApiClient(configuration) as api_client:
   # Create an instance of the API class
    api_instance = objects_api.ObjectsApi(api_client)
   # Step 2: Create a request to export project tsv file
    body_export_object_set_object_set_export_post = BodyExportObjectSetObjectSetExportPost(
        filters=ProjectFilters(statusfilter="PVD"),
        request=ExportReq(
              project_id=cfg['proj_id'],
              exp_type="TSV",
              use_latin1=False,
              tsv_entities="OPASH",  # entities to be exported: O(bjects), P(rocess), A(cquisition), S(ample), classification H(istory)
              split_by="", # No split
              coma_as_separator=False, # Set decimal separator to point
              format_dates_times=False,
              with_images=False,
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

    # Report job ID
    job_id=api_jobresponse['job_id']
    api_jobinstance = jobs_api.JobsApi(api_client)
    print("Creating export file with job ID: ", job_id, sep=' ')

    # Necessary break between step 3 and 4
    # Insert a progress bar to allow for the job to be done based on get_job state. Use time.sleep instead
    bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
    i = 0
    while True:
        i += 1
        thread = api_jobinstance.get_job(job_id, async_req=True)
        result = thread.get()
        if result.state=='F':
            break
        bar.update(i)
    # Lines 163:169 are not working. Laurent is updating get_job_file.py.
    # Uncomment line 165 when get_job_file.py is working to proceed to step 4
    # api_downloadresponse = api_jobinstance.get_job_file(job_id)
    # api_downloadresponse = api_jobinstance.get_job_file(job_id, async_req=True, _preload_content=False,_return_http_data_only=False)
    # urlresponse = api_downloadresponse.get()
    # zip_file = urlresponse.headers['content-disposition'][22:-1]
    # url = urlresponse.geturl()


# Step 4: Get/Download export zipfile
zip_file = "ecotaxa_export_{}_{}.zip".format(str(cfg['proj_id']), datetime.datetime.now().strftime("%Y%m%d_%H%M"))
path_to_zip=Path('/Volumes/WD/Pssdb/Data/API/').expanduser() / zip_file
print("\nExporting file", zip_file, "to", path_to_export, ", please wait", sep=' ')

with requests.Session() as sess:
    url = BASE_URL + "/api/jobs/%d/file" % job_id
    rsp = sess.get(url, headers={"Authorization": "Bearer " + json_token}, stream=True)
    with open(path_to_zip, "wb") as fd:
        for a_chunk in rsp.iter_content():  # Loop over content, i.e. eventual HTTP chunks
            fd.write(a_chunk)

print("Download completed: ",path_to_zip, sep='')

quit()