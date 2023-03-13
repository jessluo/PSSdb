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
import yaml # requires installation of PyYAML package
# read git-tracked config file (text file) with inputs:  project ID, output directory
path_to_config=Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_to_git=Path(cfg['git_dir']).expanduser()
path_to_data=path_to_git / cfg['dataset_subdir']
path_to_config_pw = path_to_git / 'scripts' / 'Ecotaxa_API_pw.yaml'
with open(path_to_config_pw, 'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

# Modules for webpage handling/scraping:
import urllib3
import requests
from bs4 import BeautifulSoup

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
from natsort import natsorted

## Functions start here:

def Ecopart_export(project,localpath,username,password):
    """
    This function uses basic web scraping python modules to export data hosted on Ecopart (https://ecopart.obs-vlfr.fr/). Exported files include:\n
      -metadata: profile ID, cruise ID, station ID, data owner, raw profile folder, instrument, CTD name, datetime (yyyy-mm-dd hh:mm:ss), latitude (decimal degree), longitude (decimal degree), size calibration coefficients (aa,exp), pixel size (millimeter_per_pixel), particle filename, plankton filename, project name\n
      -particles: profile ID, raw profile folder, datetime (yyyy-mm-dd hh:mm:ss), project name, mid depth bin (meters), cumulative volume in depth bin (number of images within depth bin x image volume), abundances in size bins (particles_per_liter), summed biovolume in size bins (cubic_millimeter_per_liter), matching CTD data\n
      -zooplankton: profile ID, raw profile folder, datetime (yyyy-mm-dd hh:mm:ss), mid depth bin (meters), cumulative volume in depth bin (number of images within depth bin x image volume), project name, group-specific abundances (particles_per_cubic_meter), group-specific summed biovolumes (cubic_millimeter_per_liter), group-specific average equivalent circular diameter (millimeter)
             :param project: dictionary with key ID and name value for Ecopart project.
             :param localpath: Path locating the directory where project export will be stored
             :param username: Email for Ecopart account authentication
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
        soup = BeautifulSoup(response.text, 'lxml')
        if len(soup.find_all('table')[-1].find_all('span')): # If export allowed, return button to start Task, else return span stating No data to export
            if (soup.find_all('table')[-1].find_all('span')[0].contents[0]=='No data to export, change your criteria'):
                print('Project not exportable. Updating project_list_all/ecopart and skipping export')
                sheets = pd.ExcelFile(path_to_data.parent / cfg['proj_list']).sheet_names
                if ('ecopart' in sheets):
                    project_list = pd.read_excel(path_to_data.parent / cfg['proj_list'], sheet_name="ecopart")
                    project_list.loc[project_list.Project_ID.astype(str)==str(list(project.keys())[0]),'PSSdb_access']=False
                    project_list=project_list.sort_values(['PSSdb_access','Project_ID'],ascending=False)
                    with pd.ExcelWriter(str(path_to_data.parent / cfg['proj_list']), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                        project_list.to_excel(writer, sheet_name='ecopart', index=False)
                return
        if response.ok:
            # Start export task for detailed data
            r = session.post(url_export,data={'what': 'DET', 'fileformatd': 'TSV', 'aggregatefilesd': 'Y', 'starttask': 'Y'})
            task = r.text[(r.text.find("show_url = ") + 11):(r.text.find("show_url = ") + 50)].split('\n')[0].replace('"/Task/Show/', '').replace('"', '')
            status = session.get('https://ecopart.obs-vlfr.fr/Task/GetStatus/{}'.format(task)).json()
            with tqdm(desc='Working on detailed export', total=1000, bar_format='{desc}{bar}', position=0, leave=True) as bar:
                while ('IsComplete' not in status['d'].keys()):
                    time.sleep(3)
                    status = session.get('https://ecopart.obs-vlfr.fr/Task/GetStatus/{}'.format(task)).json()
                    percent = status['d']['PercentComplete']
                    bar.set_description("Working on detailed export %s%%" % percent, refresh=True)
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
    df_meta=pd.read_table(Path(path_to_zip.parent / 'ecopart_export_detailed_{}_metadata.tsv'.format(list(project.keys())[0])),encoding='latin-1')
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
            with tqdm(desc='Working on raw export', total=1000, bar_format='{desc}{bar}', position=0, leave=True) as bar:
                while ('IsComplete' not in status['d'].keys()):
                    time.sleep(3)
                    status = session.get('https://ecopart.obs-vlfr.fr/Task/GetStatus/{}'.format(task)).json()
                    percent = status['d']['PercentComplete']
                    bar.set_description("Working on raw export %s%%" % percent, refresh=True)
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
    # Replace particles/plankton filename in metadata
    df_meta = pd.read_table(Path(path_to_zip.parent / 'ecopart_export_raw_{}_metadata.tsv'.format(list(project.keys())[0])),encoding='latin-1')
    df_meta['Particle filename'] = df_meta['Particle filename'].apply(lambda file: Path(path_to_zip.parent / 'ecopart_export_raw_{}_particles.tsv'.format(list(project.keys())[0])).name if file != 'no data available' else '')
    df_meta['Plankton filename'] = df_meta['Plankton filename'].apply(lambda file: Path(path_to_zip.parent / 'ecopart_export_raw_{}_zooplankton.tsv'.format(list(project.keys())[0])).name if file != 'no data available' else '')
    df_meta['CTD filename'] = df_meta['CTD filename'].apply(lambda file: Path(path_to_zip.parent / 'ecopart_export_raw_{}_CTD.tsv'.format(list(project.keys())[0])).name if file != 'no data available' else '')
    # Correcting newline character created in ctd desc column for some files
    df_meta['ctd_desc'] = df_meta.ctd_desc.str.replace('\r|\n', '', regex=True)
    df_meta.to_csv(Path(path_to_zip.parent / 'ecopart_export_raw_{}_metadata.tsv'.format(list(project.keys())[0])),sep="\t", index=False)
    rawfiles_path = list(path_to_zip.parent.glob('*_PAR_raw*.tsv'))
    if len([path for path in rawfiles_path if 'black.tsv' in str(path.name)]):
        [file.unlink(missing_ok=True) for file in rawfiles_path if 'black.tsv' in str(file.name)]

    # Filter out dark-field raw particles files used for UVP6 data calibration
    rawfiles_path = [path for path in rawfiles_path if 'black.tsv' not in str(path.name)]

    if len(rawfiles_path): # Attention: Using index_col=False since some UVP6 raw particle datafiles have an extra column with datetime (uvp6_sn000130hf_2021_amazomix)
         pd.concat(map(lambda path: pd.read_csv(path, sep='\t', encoding='latin1', index_col=False).assign(profileid=path.name[1 + path.name.find('_'):path.name.find('_PAR_raw_{}'.format(status['d']['ExtraAction'][ 11 + status['d']['ExtraAction'].find('export_raw_'): status['d'][ 'ExtraAction'].find('.zip')]))]),natsorted(rawfiles_path)), axis=0).to_csv(Path(path_to_zip.parent / 'ecopart_export_raw_{}_particles.tsv'.format(list(project.keys())[0])), sep="\t", index=False)
         [file.unlink(missing_ok=True) for file in  list(path_to_zip.parent.glob('*_PAR_raw*.tsv'))]
    if len(list(path_to_zip.parent.glob('*_ZOO_raw*.tsv'))):
         pd.concat(map(lambda path: pd.read_csv(path, sep='\t', encoding='latin1').assign(profileid=path.name[1+path.name.find('_'):path.name.find('_ZOO_raw_{}'.format(status['d']['ExtraAction'][11+status['d']['ExtraAction'].find('export_raw_'):status['d']['ExtraAction'].find('.zip')] ))]),natsorted(list(path_to_zip.parent.glob('*_ZOO_raw*.tsv')))), axis=0).to_csv(Path(path_to_zip.parent / 'ecopart_export_raw_{}_zooplankton.tsv'.format(list(project.keys())[0])),sep="\t", index=False)
         [file.unlink(missing_ok=True) for file in list(path_to_zip.parent.glob('*_ZOO_raw*.tsv'))]
    if len(list(path_to_zip.parent.glob('*_CTD_raw*.tsv'))):
         pd.concat(map(lambda path: pd.read_csv(path, sep='\t', encoding='latin1').assign(profileid=path.name[1+path.name.find('_'):path.name.find('_CTD_raw_{}'.format(status['d']['ExtraAction'][11+status['d']['ExtraAction'].find('export_raw_'):status['d']['ExtraAction'].find('.zip')] ))]),natsorted(list(path_to_zip.parent.glob('*_CTD_raw*.tsv')))), axis=0).to_csv(Path(path_to_zip.parent / 'ecopart_export_raw_{}_CTD.tsv'.format(list(project.keys())[0])),sep="\t", index=False)
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
                filters=ProjectFilters(),  # ProjectFilters(statusfilter="PVD") for Predicted, Validated, dubious. Leave empty to get unclassified ROI
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
    with tqdm(desc='Working on project export', total=1000, bar_format='{desc}{bar}', position=0, leave=True) as bar:
            job_status = 'R'  # percent = 0
            while job_status not in ('F', 'E'):  # percent!=100:
                time.sleep(2)  # Check the job status every 2 seconds. Modify as needed
                thread = api_jobinstance.get_job(job_id, async_req=True)
                result = thread.get()
                job_status = result.state
                percent = result.progress_pct
                bar.set_description("Working on project export %s%%" % percent, refresh=True)
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

# IFCB dashboard functions here:
## create list of dates available for the project of interest
def get_df_list_IFCB (base_url, Project_ID,  startdate=20000101, enddate=21000101):
    """
    Objective: generate list of files for a dataset in the IFCB dashboard
    :param startdate: set start date to  select data from the project, in format YYYYMMDD
    :param enddate: set end  date to  select data from the project, in format YYYYMMDD
    :param base_url: WHOI or CALOOS url
    :param Project_ID: data set contained within the WHOI or CALOOS dashboard
    :return: data frame that contains the dataset names under the 'pid' column
    """
    import requests
    import json
    import pandas as pd
    link = base_url + 'api/list_bins?dataset=' + Project_ID # this path is for the data LIST only, modified 9/12/2022, Project_source in the project list obtained in step0  has the complete url
    #metadata_link = "https://ifcb.caloos.org/timeline?dataset="
    r = requests.get(link)
    r_content = r.content
    r_content = json.loads(r_content.decode('utf-8'))
    all_file_info = pd.DataFrame.from_dict(r_content['data'], orient = 'columns')
    # generate a new column to index by dates:
    date_info = []
    for i in all_file_info.loc[:, 'sample_time']:
        date_info.append(int(pd.to_datetime(i).date().strftime('%Y%m%d')))
    all_file_info['date_info'] = date_info # [int(sample[1:9]) for sample in all_file_info['pid']] changed
    # constrain date ranges for file download.
    #start_stamp=int(startdate.strftime('%Y%m%d'))
    #end_stamp=int(y['enddate'].strftime('%Y%m%d'))
    file_info = all_file_info.loc[all_file_info['date_info']>= startdate].reset_index(drop =True) # | all_file_info['date_info']<= end_stamp]
    file_info = file_info.loc[file_info['date_info']<= enddate].reset_index(drop =True)
    return file_info

# We adopted some code contributed by Joe Futrelle (WHOI): https://github.com/joefutrelle/pyifcb

# make function to get datasets from urls
def df_from_url(url):
    """
    Objective: use requests to build dataframes with ROI information from contents obtained from IFCB dashboards
    :param url: url for the dataset
    :return: a dataframe with the desired information (see outputs of IFCB dashboards)
    """
    import requests
    import pandas as pd
    data = requests.get(url)
    data = data.content
    data = data.decode('utf-8')
    data = data.split('\n')
    data_dict = {}
    for n, r in enumerate(data):
        roi_info = r.split(',')
        data_dict[n] = roi_info
    data_dict.popitem()
    df = pd.DataFrame.from_dict(data_dict, orient='index')
    df.columns = df.iloc[0]
    df.drop(df.index[0], inplace = True)
    df.columns = df.columns.str.replace('\r', '')
    if df.columns.__contains__(None):
        print( url + ' file does not exist')
    else:
        return df


def metadata_dict(dashboard, pid_id):
    """
    :param dashboard: url to either the WHOI or CALOOS dashboard
    :param pid: the identifier for the data in each time bin displayed on the timeline in the dashboard
    :return: dictionary with the metadata
    """
    import requests
    import re
    url = f'{dashboard}api/bin/{pid_id}?include_coordinates=False'
    r = requests.get(url)
    assert r.ok
    record = r.json()
    return {
        'scale': record['scale'],
        'datetime': record['timestamp_iso'],
        'previous_bin': record['previous_bin_id'],
        'next_bin': record['next_bin_id'],
        'latitude': record['lat'],
        'longitude': record['lng'],
        'depth': record['depth'],
        'instrument': record['instrument'],
        'number_of_rois': record['num_images'],
        'ml_analyzed': float(re.sub(' .*', '', record['ml_analyzed'])),
        'concentration': record['concentration'],
        'sample_type': record['sample_type'],
        'cruise': record['cruise'],
        'cast': record['cast'],
        'niskin': record['niskin'],
    }

def ifcb_config(dashboard_url, Project_ID, pid_id):
    """
    :param dashboard: url to either the WHOI or CALOOS dashboard
    :param Project_ID: name of the dataset
    :param pid_id: the identifier for the data in each time bin displayed on the timeline in the dashboard
    :return: dictionary with the metadata
    NOTE: # follow instructions to install on : https://github.com/joefutrelle/pyifcb.
    Make sure to clone the git repository, and run the commands from the directory where the repo files are downloaded
    """
    import ifcb
    url = f'{dashboard_url}{Project_ID}/{pid_id}.hdr'
    with ifcb.open_url(url, images=False) as sample_bin:
        minBlobArea = sample_bin.headers['minimumBlobArea']
        imageResizeFactor = sample_bin.headers['imageResizeFactor']
        blobXgrowAmount = sample_bin.headers['blobXgrowAmount']
        blobYgrowAmount = sample_bin.headers['blobYgrowAmount']
        PMTAhighVoltage = sample_bin.headers['PMTAhighVoltage']
        PMTBhighVoltage = sample_bin.headers['PMTBhighVoltage']
        PMTChighVoltage = sample_bin.headers['PMTChighVoltage']
        SyringeSampleVolume = sample_bin.headers['SyringeSampleVolume']
        PMTAtriggerThreshold_DAQ_MCConly = sample_bin.headers['PMTAtriggerThreshold_DAQ_MCConly']
        PMTBtriggerThreshold_DAQ_MCConly = sample_bin.headers['PMTBtriggerThreshold_DAQ_MCConly']
        PMTCtriggerThreshold_DAQ_MCConly = sample_bin.headers['PMTCtriggerThreshold_DAQ_MCConly']
    return{
        'minBlobArea': minBlobArea,
        'imageResizeFactor': imageResizeFactor,
        'blobXgrowAmount': blobXgrowAmount,
        'blobYgrowAmount': blobYgrowAmount,
        'PMTAhighVoltage': PMTAhighVoltage,
        'PMTBhighVoltage': PMTBhighVoltage,
        'PMTChighVoltage': PMTChighVoltage,
        'SyringeSampleVolume': SyringeSampleVolume,
        'PMTAtriggerThreshold_DAQ_MCConly': PMTAtriggerThreshold_DAQ_MCConly,
        'PMTBtriggerThreshold_DAQ_MCConly': PMTBtriggerThreshold_DAQ_MCConly,
        'PMTCtriggerThreshold_DAQ_MCConly': PMTCtriggerThreshold_DAQ_MCConly
    }

def roi_number(dashboard, pid_id):
    """
    objective: simplified version of the function above. ideal to just get the number of ROIs
    :param dashboard: url to either the WHOI or CALOOS dashboard
    :param pid: the identifier for the data in each time bin displayed on the timeline in the dashboard
    :return: dictionary with the metadata
    """
    import requests
    import re
    url = f'{dashboard}api/bin/{pid_id}?include_coordinates=False'
    r = requests.get(url)
    assert r.ok
    record = r.json()
    return record['num_images']

def IFCB_dashboard_export(dashboard_url, Project_source, Project_ID, path_download, start_date, end_date):

    """
    objective: read the project_list dataframe and export the CALOOS and WHOI IFCB dashboard datasets
    :param df: dataframe with project source, local path to export
    """
    import numpy as np
    import time
    start = time.time()
    METRIC = 'temperature'  # needs to be defined based on Karina's script
    # use the file list to download data and store it locally
    file_info = get_df_list_IFCB(dashboard_url, Project_ID,  startdate=start_date, enddate=end_date) # remember this fuction has the option to define an interval of dates
    last_file = file_info.loc[len(file_info.index) - 1, 'pid']
    # the loop uses the file_info dataframe to get the 'pid' (file name) and download features and class scores files
    print ('extracting features files, metadata and top 5 class scores for timeseries ' + Project_ID +' stored in ' + Project_source)
    file_numbers = 0
    df_concatenated = pd.DataFrame()
    for i in tqdm(file_info.loc[:, 'pid']): # timeit and progressbar can be used here
        try:
            # generate features files
            pid_id = str(i)
            #print(pid_id)
            features_filename = dashboard_url + Project_ID + '/' + pid_id + '_features.csv'
            try:
                features_df = df_from_url(features_filename)
                # obtain metadata
                met_dict = metadata_dict(dashboard=dashboard_url, pid_id=pid_id)
                sample_desc = ifcb_config(dashboard_url=dashboard_url, Project_ID =Project_ID, pid_id=pid_id)
                features_clean = pd.DataFrame()

                # extract data of interest from features file:
                features_clean['area'] = features_df['Area']
                features_clean['biovolume'] = features_df['Biovolume']
                features_clean['equiv_diameter'] = features_df['EquivDiameter']
                features_clean['major_axis_length'] = features_df['MajorAxisLength']
                features_clean['minor_axis_length'] = features_df['MinorAxisLength']
                features_clean['solidity'] = features_df['Solidity']
                features_clean['summed_area'] = features_df['summedArea']
                features_clean['summed_biovolume'] = features_df['summedBiovolume']
                features_clean['summed_major_axis_length'] = features_df['summedMajorAxisLength']
                features_clean['summed_minor_axis_length'] = features_df['summedMinorAxisLength']

                #add metadata
                features_clean['project_ID'] = Project_ID
                features_clean['sample_id'] = pid_id
                features_clean['datetime'] = met_dict['datetime']
                #features_df['date'] = str(pd.to_datetime(met_dict['datetime']).date().strftime('%Y%m%d'))
                #features_df['time'] = str(pd.to_datetime(met_dict['datetime']).time().strftime('%H%M%S'))
                features_clean['pixel_to_micron'] = met_dict['scale']
                features_clean['latitude'] = met_dict['latitude']
                #features_clean=features_clean.replace({'latitude': {0: float('nan'), -9999999: float('nan')}})# replace 0 lat values with Nan
                features_clean['longitude'] = met_dict['longitude']
                #features_clean=features_clean.replace({'longitude': {0: float('nan'), -9999999: float('nan')}})  # replace 0 lat values with Nan
                features_clean['depth'] = met_dict['depth']
                features_clean['vol_analyzed'] = met_dict['ml_analyzed']
                features_clean['sample_type'] = met_dict['sample_type']
                features_clean['sample_cruise'] = met_dict['cruise']
                features_clean['number_of_rois'] = met_dict['number_of_rois']
                features_clean['concentration'] = met_dict['concentration']

                #add sampling description:
                features_clean['minBlobArea'] = sample_desc['minBlobArea']
                features_clean['imageResizeFactor'] = sample_desc['imageResizeFactor']
                features_clean['blobXgrowAmount'] = sample_desc['blobXgrowAmount']
                features_clean['blobYgrowAmount'] = sample_desc['blobYgrowAmount']
                features_clean['PMTAhighVoltage'] = sample_desc['PMTAhighVoltage']
                features_clean['PMTBhighVoltage'] = sample_desc['PMTBhighVoltage']
                features_clean['PMTChighVoltage'] = sample_desc['PMTChighVoltage']
                features_clean['SyringeSampleVolume'] = sample_desc['SyringeSampleVolume']
                features_clean['PMTAtriggerThreshold_DAQ_MCConly'] = sample_desc['PMTAtriggerThreshold_DAQ_MCConly']
                features_clean['PMTBtriggerThreshold_DAQ_MCConly'] = sample_desc['PMTBtriggerThreshold_DAQ_MCConly']
                features_clean['PMTCtriggerThreshold_DAQ_MCConly'] = sample_desc['PMTCtriggerThreshold_DAQ_MCConly']



                # now generate dataframe for the class scores
                class_filename = dashboard_url + Project_ID + '/' + pid_id + '_class_scores.csv'
                class_df = df_from_url(class_filename)
                features_clean['roi_id'] = class_df['pid']
                class_df = class_df.set_index('pid')
                class_df = class_df.astype(float)
                # create lists for top five classes and their scores
                class_1 = []
                class_1_score = []
                class_2 = []
                class_2_score = []
                class_3 = []
                class_3_score = []
                class_4 = []
                class_4_score = []
                class_5 = []
                class_5_score = []
                try:
                    for roi in class_df.index:
                        #print(roi)
                        index_top5 = np.argsort(-class_df.loc[roi].values)[:5]
                        #print(index_top5)
                        class_1.append(class_df.columns[index_top5[0]])
                        class_1_score.append(class_df.loc[roi, class_df.columns[index_top5[0]]])

                        class_2.append(class_df.columns[index_top5[1]])
                        class_2_score.append(class_df.loc[roi, class_df.columns[index_top5[1]]])

                        class_3.append(class_df.columns[index_top5[2]])
                        class_3_score.append(class_df.loc[roi, class_df.columns[index_top5[2]]])

                        class_4.append(class_df.columns[index_top5[3]])
                        class_4_score.append(class_df.loc[roi, class_df.columns[index_top5[3]]])

                        class_5.append(class_df.columns[index_top5[4]])
                        class_5_score.append(class_df.loc[roi, class_df.columns[index_top5[4]]])
                except:
                    print(i + ' does not have a class score')
                features_clean['class_1'] = class_1
                features_clean['class_1_score'] = class_1_score

                features_clean['class_2'] = class_2
                features_clean['class_2_score'] = class_2_score

                features_clean['class_3'] = class_3
                features_clean['class_3_score'] = class_3_score

                features_clean['class_4'] = class_4
                features_clean['class_4_score'] = class_4_score

                features_clean['class_5'] = class_5
                features_clean['class_5_score'] = class_5_score


                #replace with nan the class and class scores below 0.0001
                cols = ["class_1", "class_2", "class_3", "class_4", "class_5"]
                cols_scores = ["class_1_score", "class_2_score", "class_3_score", "class_4_score", "class_5_score"]
                for no, c in enumerate(cols_scores):
                    features_clean[c] = features_clean[c].mask(features_clean[c] < 0.0001, np.nan)
                    features_clean[cols[no]] = features_clean[cols[no]].mask(features_clean[c] < 0.0001, np.nan)

                features_clean['datetime'] = pd.to_datetime(features_clean['datetime'])
                year = features_clean.loc[1, 'datetime'].strftime('%Y')
                dataset_id = features_clean.loc[1, 'roi_id'].split('_')
                dataset_id = dataset_id[1]
                if dashboard_url == 'https://ifcb.caloos.org/':
                    dashboard_id = 'CALOOS'
                elif dashboard_url== 'https://ifcb-data.whoi.edu/':
                    dashboard_id = 'WHOI'

                # add sampling description:
                features_clean['minBlobArea'] = sample_desc['minBlobArea']
                features_clean['PMTAhighVoltage'] = sample_desc['PMTAhighVoltage']
                features_clean['PMTBhighVoltage'] = sample_desc['PMTBhighVoltage']
                features_clean['PMTChighVoltage'] = sample_desc['PMTChighVoltage']
                features_clean['SyringeSampleVolume'] = sample_desc['SyringeSampleVolume']
                features_clean['PMTAtriggerThreshold_DAQ_MCConly'] = sample_desc['PMTAtriggerThreshold_DAQ_MCConly']
                features_clean['PMTBtriggerThreshold_DAQ_MCConly'] = sample_desc['PMTBtriggerThreshold_DAQ_MCConly']
                features_clean['PMTCtriggerThreshold_DAQ_MCConly'] = sample_desc['PMTCtriggerThreshold_DAQ_MCConly']

                df_concatenated = pd.concat([df_concatenated, features_clean], ignore_index=True)

                #if len(df_concatenated.index) <= 100000:
                    #df_concatenated = pd.concat([df_concatenated, features_df], ignore_index=True)
                    #pass
                if pid_id == last_file :
                    file_numbers = file_numbers + 1
                    print ('saving file # ' + str(file_numbers) + ' of '+ Project_ID)
                    df_concatenated.to_csv(path_download + '/' + dashboard_id +'_'+ dataset_id +'_'+ df_concatenated.loc[1, 'project_ID'] +'_'+ year + '_features_' + str(file_numbers) +'.tsv', sep='\t')
                    df_concatenated = pd.DataFrame()
                    #pass
                elif len(df_concatenated.index) > 500000:
                    file_numbers = file_numbers + 1
                    print ('saving file # ' + str(file_numbers) + ' of '+ Project_ID)
                    df_concatenated.to_csv(path_download + '/' + dashboard_id +'_'+ dataset_id +'_'+ df_concatenated.loc[1, 'project_ID'] +'_'+ year + '_features_' + str(file_numbers) +'.tsv', sep='\t')
                    df_concatenated = pd.DataFrame()


                # class_df.to_csv(path_download + '/' + str(i) + '_class_scores.csv', sep='\t') 11/16/2022 decided not to keep class scores
                # print(str(i) + ' download done ')
            except:
                print('there is no features or class_scores files for ' + str(file_info.loc[i, 'pid']))
        except:
            pass
    elapsed_time_fl = (time.time() - start)
    print(Project_ID + ' download took ' + str((elapsed_time_fl/3600)) + ' hours')
    print (Project_ID + ' download done')
