## Objective: This script defines functions to export projects from various data portal (i.e. Ecotaxa, Ecopart, IFCB dashboard)

# Workflow for Ecotaxa_export function:
# This script uses 4 API steps:
# Step 1: Configure the API using your Ecotaxa login/pw (info stored in untracked Ecotaxa_API_pw.yaml)
# Step 2: Generate a request to export a project tsv based on project ID (info stored in git-tracked Ecotaxa_API.yaml)
# Step 3: Generate a job/task(=export) based on the request (step 2)
# Step 4: Get the job file (ecotaxa_export_projectID_date_timeZ.zip)

## Useful documentations:
# See documentation for step 1 of Ecotaxa_export: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/AuthentificationApi.md
# See documentation for step 2 of Ecotaxa_export: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ExportReq.md
# See documentation for step 3 of Ecotaxa_export: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ObjectsApi.md#export_object_set
# See documentation for step 4 of Ecotaxa_export: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/JobsApi.md#get_job_file
# https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests

## Python modules:

# Modules for data and path handling:
from pathlib import Path  # Handling of path object
import shutil # Delete uncompressed export zip folder
import pandas as pd

# Modules for webpage handling/scraping:
import urllib3
import requests

# Modules for Ecotaxa API interface:
import warnings
warnings.filterwarnings('ignore', module='urllib3')

import ecotaxa_py_client
from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.model.login_req import LoginReq
from ecotaxa_py_client.api import objects_api
from ecotaxa_py_client.model.body_export_object_set_object_set_export_post import BodyExportObjectSetObjectSetExportPost
from ecotaxa_py_client.model.project_filters import ProjectFilters
from ecotaxa_py_client.model.export_req import ExportReq
from ecotaxa_py_client.api import jobs_api

# Progress bar modules
#import progressbar # Attention: Use pip install progressbar2 to install
from tqdm import tqdm
import datetime, time # Time module for break and current date

## Functions start here:

def Ecopart_export(project,localpath,username,password):
    """
    This function uses basic web scraping python modules to export data hosted on Ecopart (https://ecopart.obs-vlfr.fr/). Exported files include:\n
      -metadata: profile ID, cruise ID, station ID, data owner, raw profile folder, instrument, CTD name, datetime (yyyy-mm-dd hh:mm:ss), latitude (decimal degree), longitude (decimal degree), size calibration coefficients (aa,exp), pixel size (millimeter_per_pixel), particle filename, plankton filename, project name\n
      -particles: profile ID, raw profile folder, datetime (yyyy-mm-dd hh:mm:ss), project name, mid depth bin (meters), cumulative volume in depth bin (number of images within depth bin x image volume), abundances in size bins (particles_per_liter), summed biovolume in size bins (cubic_millimeter_per_liter), matching CTD data\n
      -zooplankton: profile ID, raw profile folder, datetime (yyyy-mm-dd hh:mm:ss), mid depth bin (meters), cumulative volume in depth bin (number of images within depth bin x image volume), project name, group-specific abundances (particles_per_cubic_meter), group-specific summed biovolumes (cubic_millimeter_per_liter), group-specific average equivalent circular diameter (millimeter)
             :param project: dictionary with key ID and name value for Ecopart project.
             :param localpath: Path locating the directory where project export will be stored
             :param username: Email for Ecopart acount authentication
             :param password: Password for Ecopart account authentication
             :return: Ecopart export files in localpath/project/ecopart_export_detailed_uvp5_sn010_2014_m108_metadata/particles/zooplankton.tsv.
    """
    path_to_datafile = localpath
    path_to_datafile.mkdir(parents=True, exist_ok=True)
    with requests.Session() as session:
        # Get project ID given project name
        print('Exporting Ecopart project {}. Please wait'.format(list(project.keys())[0]))
        url_export = 'https://ecopart.obs-vlfr.fr/Task/Create/TaskPartExport?MapN=&MapW=&MapE=&MapS=&filt_fromdate=&filt_todate=&filt_instrum=&filt_proftype=&filt_uproj={}&filt_depthmin=&filt_depthmax=&XScale=I&TimeScale=R'.format(list(project.keys())[0])
        data_login = {"email": username, "password": password}
        # Authenticate with email/password
        session.post('https://ecopart.obs-vlfr.fr/login', data=data_login, headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
        # Check that access is granted and export is allowed
        response = session.get(url=url_export)
        if response.ok:
            # Start export task for detailed data
            r = session.post(url_export,data={'what': 'DET', 'fileformatd': 'TSV', 'aggregatefilesd': 'Y', 'starttask': 'Y'})
            task = r.text[(r.text.find("show_url = ") + 11):(r.text.find("show_url = ") + 50)].split('\n')[0].replace('"/Task/Show/', '').replace('"', '')
            status = session.get('https://ecopart.obs-vlfr.fr/Task/GetStatus/{}'.format(task)).json()
            with tqdm(desc='Working on Job', total=1000, bar_format='{desc}{bar}', position=0, leave=True) as bar:
                while ('IsComplete' not in status['d'].keys()):
                    time.sleep(3)
                    status = session.get('https://ecopart.obs-vlfr.fr/Task/GetStatus/{}'.format(task)).json()
                    percent = status['d']['PercentComplete']
                    bar.set_description("Working on Job %s%%" % percent, refresh=True)
                    # and update progress bar
                    ok = bar.update(n=10)
                    if ('IsComplete' in status['d'].keys()):
                        if (status['d']['IsComplete'] == 'Y'):
                            break

            url_to_zip = 'https://ecopart.obs-vlfr.fr/{}'.format(status['d']['ExtraAction'][(status['d']['ExtraAction'].find('href') + 6):(status['d']['ExtraAction'].find('.zip') + 4)])
    with requests.Session() as session:
            session.post('https://ecopart.obs-vlfr.fr/login', data=data_login, headers={ 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
            rsp = session.get(url_to_zip, stream=True)
            path_to_zip = path_to_datafile / list(project.values())[0] / "ecopart_detailed_export_{}_{}Z.zip".format(list(project.keys())[0], datetime.datetime.utcnow().strftime("%Y%m%d_%H%M"))
            path_to_zip.parent.mkdir(parents=True, exist_ok=True)
            with open(path_to_zip, 'wb') as f:
                shutil.copyfileobj(rsp.raw, f,length=16*1024*1024)
            rsp.close()
    shutil.unpack_archive(path_to_zip, path_to_zip.parent)  # Unzip export file
    # Rename export files
    if len(list(path_to_zip.parent.glob('export_detailed_*_Export_metadata_summary.tsv'))):
        Path(list(path_to_zip.parent.glob('export_detailed_*_Export_metadata_summary.tsv'))[0]).rename(Path(path_to_zip.parent / 'ecopart_export_detailed_{}_metadata.tsv'.format(list(project.keys())[0])))
    if len(list(path_to_zip.parent.glob('export_detailed_*_PAR_Aggregated.tsv'))):
        Path(list(path_to_zip.parent.glob('export_detailed_*_PAR_Aggregated.tsv'))[0]).rename(Path(path_to_zip.parent / 'ecopart_export_detailed_{}_particles.tsv'.format(list(project.keys())[0])))
    if len(list(path_to_zip.parent.glob('export_detailed_*_ZOO_Aggregated.tsv'))):
        Path(list(path_to_zip.parent.glob('export_detailed_*_ZOO_Aggregated.tsv'))[0]).rename(Path(path_to_zip.parent / 'ecopart_export_detailed_{}_zooplankton.tsv'.format(list(project.keys())[0])))
    # Replace particles/plankton filename in metadata
    df_meta=pd.read_table(Path(path_to_zip.parent / 'ecopart_export_detailed_{}_metadata.tsv'.format(list(project.keys())[0])))
    df_meta['Particle filename']=Path(path_to_zip.parent / 'ecopart_export_detailed_{}_particles.tsv'.format(list(project.keys())[0])).name
    df_meta['Plankton filename'] = Path(path_to_zip.parent / 'ecopart_export_detailed_{}_zooplankton.tsv'.format(list(project.keys())[0])).name
    df_meta.to_csv(Path(path_to_zip.parent / 'ecopart_export_detailed_{}_metadata.tsv'.format(list(project.keys())[0])), sep="\t",index=False)
    # Delete original zip file
    path_to_zip.unlink(missing_ok=True)
    with requests.Session() as session:
        # Authenticate with email/password
        session.post('https://ecopart.obs-vlfr.fr/login', data=data_login, headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
        # Check that access is granted and export is allowed
        response = session.get(url=url_export)
        if response.ok:
            # Start export task for raw data
            r = session.post(url_export,data={'what': 'RAW', 'includenotvalidatedr': 'Y', 'starttask': 'Y'})
            task = r.text[(r.text.find("show_url = ") + 11):(r.text.find("show_url = ") + 50)].split('\n')[0].replace('"/Task/Show/', '').replace('"', '')
            status = session.get('https://ecopart.obs-vlfr.fr/Task/GetStatus/{}'.format(task)).json()
            with tqdm(desc='Working on Job', total=1000, bar_format='{desc}{bar}', position=0, leave=True) as bar:
                while ('IsComplete' not in status['d'].keys()):
                    time.sleep(3)
                    status = session.get('https://ecopart.obs-vlfr.fr/Task/GetStatus/{}'.format(task)).json()
                    percent = status['d']['PercentComplete']
                    bar.set_description("Working on Job %s%%" % percent, refresh=True)
                    # and update progress bar
                    ok = bar.update(n=10)
                    if ('IsComplete' in status['d'].keys()):
                        if (status['d']['IsComplete'] == 'Y'):
                            break
            url_to_zip = 'https://ecopart.obs-vlfr.fr/{}'.format(status['d']['ExtraAction'][(status['d']['ExtraAction'].find('href') + 6):(status['d']['ExtraAction'].find('.zip') + 4)])
    with requests.Session() as session:
         # Authenticate with email/password
        session.post('https://ecopart.obs-vlfr.fr/login', data=data_login, headers={ 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
        rsp = session.get(url_to_zip, stream=True)
        path_to_zip = path_to_datafile / list(project.values())[0] / "ecopart_raw_export_{}_{}Z.zip".format(list(project.keys())[0], datetime.datetime.utcnow().strftime("%Y%m%d_%H%M"))
        path_to_zip.parent.mkdir(parents=True, exist_ok=True)
        with open(path_to_zip, 'wb') as f:
            shutil.copyfileobj(rsp.raw, f,length=16*1024*1024)
        rsp.close()
    print("Download completed: ", path_to_zip, "\nUnpacking zip file", sep='')
    shutil.unpack_archive(path_to_zip, path_to_zip.parent)  # Unzip export file
    # Concatenate and rename export files
    Path(list(path_to_zip.parent.glob('export_raw_*_Export_metadata_summary.tsv'))[0]).rename(Path(path_to_zip.parent / 'ecopart_export_raw_{}_metadata.tsv'.format(list(project.keys())[0])))
    if len(list(path_to_zip.parent.glob('*_PAR_raw*.tsv'))):
        pd.concat(map(lambda path: pd.read_csv(path, sep='\t',encoding='latin1'),list(path_to_zip.parent.glob('*_PAR_raw*.tsv'))),axis=0).to_csv(Path(path_to_zip.parent / 'ecopart_export_raw_{}_particles.tsv'.format(list(project.keys())[0])),sep="\t",index=False)
        [file.unlink(missing_ok=True) for file in  list(path_to_zip.parent.glob('*_PAR_raw*.tsv'))]
    if len(list(path_to_zip.parent.glob('*_ZOO_raw*.tsv'))):
         pd.concat(map(lambda path: pd.read_csv(path, sep='\t', encoding='latin1'),list(path_to_zip.parent.glob('*_ZOO_raw*.tsv'))), axis=0).to_csv(Path(path_to_zip.parent / 'ecopart_export_raw_{}_zooplankton.tsv'.format(list(project.keys())[0])),sep="\t", index=False)
         [file.unlink(missing_ok=True) for file in list(path_to_zip.parent.glob('*_ZOO_raw*.tsv'))]
    if len(list(path_to_zip.parent.glob('*_CTD_raw*.tsv'))):
         pd.concat(map(lambda path: pd.read_csv(path, sep='\t', encoding='latin1'),list(path_to_zip.parent.glob('*_CTD_raw*.tsv'))), axis=0).to_csv(Path(path_to_zip.parent / 'ecopart_export_raw_{}_CTD.tsv'.format(list(project.keys())[0])),sep="\t", index=False)
         [file.unlink(missing_ok=True) for file in list(path_to_zip.parent.glob('*_CTD_raw*.tsv'))]
    # Delete original zip file
    path_to_zip.unlink(missing_ok=True)


def Ecotaxa_export(project,localpath,username,password):
    """
    This function uses Ecotaxa API module (https://github.com/ecotaxa/ecotaxa_py_client.git) to export data hosted on Ecotaxa (https://ecotaxa.obs-vlfr.fr/). \n
    Exported file includes:\n
      -object (Region of interest) fields: ROI ID, latitude, longitude, date (yyyymmd), time (hmmdd), url, minimum depth, maximum depth, annotation status (e.g. unclassified, predicted, validated, dubious), annotator name, annotator email, annotation date, annotation group, morphometric measurements (project/processing-specific)
      \n-sample fields: Sample ID, dataportal, scan operator (Zooscan), cruise ID, station ID, bottom depth, net characteristics (e.g. number of tows, mesh size, aperture surface area), volume analyzed, sampling platform and characteristics (e.g. tow duration, speed), latitude, longitude
      \n-processing fields: Processing ID (e.g. software), date (yyyymmdd), time (hhmmss), pixel size
      \n-acquisition fields: Acquisition ID, instrument, size fractionation (Zooscan)
             :param project: Project ID on Ecotaxa (i.e. unique integer).
             :param localpath: Path locating the folder where project export will be stored
             :param username: Email for Ecopart acount authentication
             :param password: Password for Ecopart account authentication
             :return: Ecotaxa export files in localpath/instrument/ecotaxa_export_projectID_yyyymmdd_hhmm.tsv. Export date is reported in UTC for version tracking
    """
    path_to_export = localpath
    path_to_export.mkdir(parents=True, exist_ok=True)
    # Generate a token based on authentication infos to login on Ecotaxa
    configuration = ecotaxa_py_client.Configuration(host="https://ecotaxa.obs-vlfr.fr/api")
    configuration.verify_ssl = False
    with ecotaxa_py_client.ApiClient(configuration) as client:
        api = authentification_api.AuthentificationApi(client)
        token = api.login(LoginReq( username=username,password=password))

    configuration = ecotaxa_py_client.Configuration(host="https://ecotaxa.obs-vlfr.fr/api", access_token=token,discard_unknown_keys=True)
    configuration.verify_ssl = False
    # Step 1b: Authentication for JSON export file
    BASE_URL = "https://ecotaxa.obs-vlfr.fr"
    with requests.Session() as sess:
        rsp = sess.post(BASE_URL + "/api/login", json={"username": username, "password": password})
        json_token = rsp.json()
    # Step 2: Create a request to export project tsv file
    with ecotaxa_py_client.ApiClient(configuration) as api_client:
         # Create an instance of the API class
        api_instance = objects_api.ObjectsApi(api_client)
        body_export_object_set_object_set_export_post = BodyExportObjectSetObjectSetExportPost(
                filters=ProjectFilters(),  # ProjectFilters(statusfilter="PVD") for Predicted, Validated, dubioud. Leave empty to get unclassified ROI
                # get P(redicted), V(alidated), D(ubious) images. Check other options for filter here: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectFilters.md
                request=ExportReq(project_id=project,  # the unique project ID of interest (integer)
                    exp_type="TSV",
                    use_latin1=False,
                    tsv_entities="OPASH", # entities to be exported: O(bjects), P(rocess), A(cquisition), S(ample), classification H(istory)
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
    with tqdm(desc='Working on Job', total=1000, bar_format='{desc}{bar}', position=0, leave=True) as bar:
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
    # Step 4: Get/Download export zipfile
    zip_file = "ecotaxa_export_{}_{}Z.zip".format(str(project), datetime.datetime.utcnow().strftime("%Y%m%d_%H%M"))  # "ecotaxa_{}".format(str(result.result['out_file']))
    path_to_zip = path_to_export / zip_file
    path_to_log = path_to_export / "job_{}.log".format(str(job_id))
    print("\nExporting file ", zip_file, " to ", path_to_export, ", please wait", sep='')

    with requests.Session() as sess:
        url = BASE_URL + "/api/jobs/%d/file" % job_id
        rsp = sess.get(url, headers={"Authorization": "Bearer " + json_token}, stream=True)
        with open(path_to_zip, "wb") as fd:
            for a_chunk in rsp.iter_content():  # Loop over content, i.e. eventual HTTP chunks
                # rsp.raise_for_status()
                fd.write(a_chunk)

    print("Download completed: ", path_to_zip, "\nUnpacking zip file", sep='')
    shutil.unpack_archive(path_to_zip, path_to_export)  # Unzip export file
    path_to_zip.unlink(missing_ok=True)  # Delete zip file
    path_to_log.unlink(missing_ok=True)  # Delete job log