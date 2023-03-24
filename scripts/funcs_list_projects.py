## Objective: This script defines functions to list projects from various data portal (i.e. Ecotaxa, Ecopart, IFCB dashboard)

## Requirements:
# Note that you need to install ecotaxa API to list Ecotaxa projects using git:
# in terminal: pip install git+https://github.com/ecotaxa/ecotaxa_py_client.git
# Note that you need to authenticate on EcoTaxa to search projects based on title/instrument filter, regardless of access rights
# Note that you need two text files in your git repository that contain the info relative to:
# authentication (Ecotaxa_API_pw.yaml, in .gitignore), output path (Ecotaxa_API.yaml, tracked)
# Note that you need to change the path of the config files on lines 58 & 62

## Workflow for Ecotaxa project list:
# This script uses 4 API functions:
# Entry 1: Search projects ID given instrument filter.
# Entry 2: Generate an exhaustive list of visible and accessible project(s) based on your Ecotaxa login/pw (info stored in untracked Ecotaxa_API_pw.yaml). Projects list stored as project_list_all
# Entry 3: Generate instrument-specific project spreadsheets used to standardize and harmonize fields of interest (i.e. needed to calculate NBSS)
# Entry 4: Bridge EcoPart-EcoTaxa projects (UVP-specific)

## Useful documentations:
# See documentation for entry 1: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectsApi.md#search_projects
# See documentation for entry 2: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectsApi.md#search_projects
# See documentation for entry 3: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ProjectsApi.md#project_set_get_column_stats,https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ObjectsApi.md#get_object_set_summary
# See documentation for entry 4: upcoming

## Usage: Run step0.list_projects to use these functions in batch and create list of all available projects


## Python modules:

# Path modules
from pathlib import Path # Handling of path object
import os

# Config modules
import yaml # requires installation of PyYAML package

# Prompt for export confirmation
import sys
import time

# dataframe/array module (used for excel files)
import pandas as pd
import numpy as np

# Ecotaxa API modules. Install module via git first
import ecotaxa_py_client
from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.model.login_req import LoginReq
from ecotaxa_py_client.api import projects_api
from ecotaxa_py_client.api import objects_api
from ecotaxa_py_client.model.project_filters import ProjectFilters

#Web scraping for Ecopart interactions
import requests
try:
    from funcs_standardize_annotations import get_all_forms
except:
    from scripts.funcs_standardize_annotations import get_all_forms
from concurrent.futures import as_completed
from requests_futures import sessions # Use pip install requests-futures
from bs4 import BeautifulSoup
import re

# IFCB dashboard functions:
try:
    from funcs_IFCB_dashboard import *
except:
    from scripts.funcs_IFCB_dashboard import *

import warnings
warnings.filterwarnings('ignore', module='urllib3')



## Functions start here:

def Ecotaxa_list(username,password,localpath):
    """
    This function uses Ecotaxa API module (https://github.com/ecotaxa/ecotaxa_py_client.git) to list projects hosted on Ecotaxa (https://ecotaxa.obs-vlfr.fr/). \n
    :param username: Username used for Ecotaxa account authentication.
    :param password: Password used for Ecotaxa account authentication.
    :param localpath: Local path locating the directory where projects will be exported.
    :return: dictionary with 'data', a dataframe including project(s) data portal (i.e https://ecotaxa.obs-vlfr.fr/), ID, title, instrument, contact name/email, accessibility under username/password, objects number, percentages of classification/validation, test set
    and 'metadata'
    """

    ## Step 1: Search for all visible and accessible projects based on authentication infos.
    with ecotaxa_py_client.ApiClient() as client:
        api = authentification_api.AuthentificationApi(client)
        token = api.login(LoginReq( username=username, password=password))
    configuration = ecotaxa_py_client.Configuration(host = "https://ecotaxa.obs-vlfr.fr/api",access_token=token, discard_unknown_keys=True)
    configuration.verify_ssl=False
    print("Searching for projects on EcoTaxa (https://ecotaxa.obs-vlfr.fr/):")
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
            api_response_project_search_accessible = api_instance.search_projects(also_others=False,title_filter='',instrument_filter='', order_field='projid')
        except ecotaxa_py_client.ApiException as e:
            print("Exception when calling EcoTaca API function search_projects: %s\n" % e)
    # Retrieve timestamp of the most recent annotation for accessible projects (Note this step takes a long time so it is currently limited to accessible projects, but can be applied to visible projects as well)
    #project_list_accessible = list(map(lambda x: str(x.projid), api_response_project_search_accessible))
    #api_response_project_get_accessible = list(map(lambda x: api_object.get_object_set(project_id=int(x), project_filters=ProjectFilters(statusfilter="PVD"),fields='obj.classif_when').details, project_list_accessible))
    #api_response_project_df_accessible = list(map(lambda x:pd.DataFrame(x,columns=['Annotation_timestamp']), api_response_project_get_accessible))
    #api_response_project_latestupdate_accessible = list(map(lambda x: max([date for date in x.Annotation_timestamp if date is not None and date is not pd.NaT]) if len([date for date in x.Annotation_timestamp if date is not None and date is not pd.NaT])>0 else pd.NaT,api_response_project_df_accessible))

    print(len(api_response_project_search_visible)+len(api_response_project_search_accessible),"projects found.",len(api_response_project_search_accessible),"projects accessible", sep=' ')
    # Step 2: Generate a dataframe with variables of interest (project ID, instrument, contact infos, access rights, and image annotation statistics) and save/overwrite in project_list_all.xslx
    df=pd.concat([pd.DataFrame({'Project_source':list(map(lambda x: 'https://ecotaxa.obs-vlfr.fr/prj/{}'.format(str(x.projid)),api_response_project_search_accessible)),
           'Project_localpath':list(map(lambda x: localpath,api_response_project_search_accessible)),
           'Project_ID':list(map(lambda x: x.projid,api_response_project_search_accessible)),
           'Project_title':list(map(lambda x: x.title,api_response_project_search_accessible)),
           'Instrument':list(map(lambda x: x.instrument,api_response_project_search_accessible)),
           'Contact_name':list(map(lambda x: x.managers[0].name if (x.contact is None) else x.contact['name'],api_response_project_search_accessible)),
           'Contact_email':list(map(lambda x: x.managers[0].email if (x.contact is None) else x.contact['email'],api_response_project_search_accessible)),
           'PSSdb_access':list(map(lambda x: str(username in list(map(lambda y: (y['email']),x.annotators))),api_response_project_search_accessible)),
           'Objects_number':list(map(lambda x: x.objcount,api_response_project_search_accessible)),
           'Percentage_classified':list(map(lambda x: x.pctclassified,api_response_project_search_accessible)),
           'Percentage_validated':list(map(lambda x: x.pctvalidated, api_response_project_search_accessible))}),
           #'Latest_update':api_response_project_latestupdate_accessible}),
                  pd.DataFrame({'Project_source':list(map(lambda x: 'https://ecotaxa.obs-vlfr.fr/prj/{}'.format(str(x.projid)),api_response_project_search_visible)),
           'Project_localpath': list( map(lambda x: localpath, api_response_project_search_visible)),
           'Project_ID':list(map(lambda x: x.projid,api_response_project_search_visible)),
           'Project_title':list(map(lambda x: x.title,api_response_project_search_visible)),
           'Instrument':list(map(lambda x: x.instrument,api_response_project_search_visible)),
           'Contact_name':list(map(lambda x: x.managers[0].name if (x.contact is None) else x.contact['name'],api_response_project_search_visible)),
           'Contact_email':list(map(lambda x: x.managers[0].email if (x.contact is None) else x.contact['email'],api_response_project_search_visible)),
           'PSSdb_access':list(map(lambda x: str(username in list(map(lambda y: (y['email']),x.annotators))),api_response_project_search_visible)),
           'Objects_number':list(map(lambda x: x.objcount,api_response_project_search_visible)),
           'Percentage_classified':list(map(lambda x: x.pctclassified,api_response_project_search_visible)),
           'Percentage_validated':list(map(lambda x: x.pctvalidated, api_response_project_search_visible))})
           #'Latest_update': pd.NaT})
                  ])
    df=df.assign(Project_test=df['Project_ID'].isin([3315,3318,3326,377,378,714,560,579,5693]).astype(str)) # Add boolean variable to identify fixed test set (3 projects per instrument)

    df_metadata=pd.DataFrame({'Variables':df.columns,'Variable_types':df.dtypes,
        'Units/Values':['','','','','','','','','#','%','%',''],
        'Description':['Project_source (URL where original project files can be exported)','Local path indicating the location of exported files storage','Project ID','Project title','Project instrument','Name of the project contact','Email of the project contact','Project accessibility. If True, export is possible with current authentication','Total number of objects (images) in the project','Percentage of predicted images','Percentage of validated images','Test set boolean identifier. If True, project is one of the test set projects']})
    return {'data':df,'metadata':df_metadata}


def Ecopart_get_sample_popover(sample_url):
    resp= requests.get(url=sample_url)
    info= {re.sub(r'\s+', '', info.split(' : ')[0]): re.sub(r'\s+', '', info.split(' : ')[1]) for info in resp.text.split('<br>\n ')} if resp.status_code == 200 else {'ID':None,'ProfileID':None,'Project':None,'Ship':None,'Cruise':None,'EcotaxaProject': None,'Lat/Lon':None,'Date/Time':None}
    if resp.status_code == 200:
        return info
    else:
        return None

def Ecopart_format_stat_table(xlml_table):
    """
        This function convert a xlml table in pandas dataframe. \n
        See documentation from:https://medium.com/analytics-vidhya/how-to-scrape-a-table-from-website-using-python-ce90d0cfb607
        :param xlml_table: A xml table obtained from BeautifulSoup.
        :return: original table reformatted with pandas
    """
    headers = []
    for column in xlml_table.find_all('th'):
        title = re.sub(r'\s+', '_',column.text)
        headers.append(title)
    df=pd.DataFrame(columns=headers)
    for row in xlml_table.find_all('tr')[1:]:
        row_data = row.find_all('td')
        items = [i.text.strip(' ')+' ({})'.format(i.find_all('a')[0].get('href').replace('mailto:','')) if len(i.find_all('a')) & len(i.text.strip(' ')) else i.text.strip(' ') for i in row_data]
        df.loc[len(df)] = items
    return df.mask(df=='')

def Ecopart_list(username,password,localpath):
    """
    This function uses basic web requests to list projects hosted on Ecopart (https://ecopart.obs-vlfr.fr/). \n
    :param username: Username used for Ecotaxa/Ecopart account authentication.
    :param password: Password used for Ecotaxa/Ecopart account authentication.
    :param localpath: Local path locating the directory where projects will be exported
    :return: dictionary with 'data', a dataframe including project(s) data portal (i.e https://ecopart.obs-vlfr.fr/), ID, title
             and 'metadata', a dataframe with variables name, type, unit , and description
    """
    print("Searching for projects on EcoPart (https://ecopart.obs-vlfr.fr/):")
    with requests.Session() as session_loggedin:
        # List projects on Ecopart
        post = session_loggedin.post('https://ecopart.obs-vlfr.fr/#')
        form = get_all_forms(post.url)[1]
        form_instrument_filter=[form for form in form.find_all('select') if form.attrs['name']=='filt_instrum']
        instruments_dict= {content.getText().split(' (')[0]: content.get('value') for content in form_instrument_filter[0].contents}
        # Retrieve form's inputs= project id/name for each instrument
        form_project_filter = [form for form in form.find_all('select') if form.attrs['name'] == 'filt_uproj']
        projects_dict = {content.getText().split(' (')[0]: content.get('value') for content in form_project_filter[0].contents}
        df_projects=pd.DataFrame.from_dict(projects_dict,'index',columns=['Project_ID']).reset_index().rename(columns={'index':'Project_title'})
        df_projects['Project_source']=df_projects.Project_ID.apply(lambda id:'https://ecopart.obs-vlfr.fr/?MapN=&MapW=&MapE=&MapS=&filt_fromdate=&filt_todate=&filt_instrum=&filt_proftype=&filt_uproj={}&filt_depthmin=&filt_depthmax=&XScale=I&TimeScale=R'.format(id))
        df_projects['Export_url'] = df_projects.Project_ID.apply(lambda id: 'https://ecopart.obs-vlfr.fr/Task/Create/TaskPartExport?MapN=&MapW=&MapE=&MapS=&filt_fromdate=&filt_todate=&filt_instrum=&filt_proftype=&filt_uproj={}&filt_depthmin=&filt_depthmax=&XScale=I&TimeScale=R'.format(id))

     # Authenticate with email/password
    form_login=get_all_forms(post.url)[0]
    login_url='https://ecopart.obs-vlfr.fr'+form_login.get('action')
    login_method=form_login.get('method').lower()
    login_info=list(map(lambda input: input.get('id'),form_login.find_all('input')))
    data_login = {"email": username, "password": password}
    getattr(session_loggedin,login_method)(url=login_url, data=data_login, headers={ 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
    #session_loggedin.post(login_url, data=data_login, headers={ 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/66.0.3359.181 Safari/537.36'})
    # Get redirection urls from the various buttons
    #json=requests.get('https://ecopart.obs-vlfr.fr/#').json()
    json=BeautifulSoup(requests.get('https://ecopart.obs-vlfr.fr/#').text,'html.parser')
    script=json.find_all('script')[6]
    functions_dict=dict(filter(lambda tuple: not tuple[0] is None,map(lambda script: (script[(script.find('function ')+9):script.find(' {')],script[(script.find('url')+6):script.find('" + $')]) if script.find('function ')!=-1 & script.find('url')!=-1 else (None,None),script.text.split('\n\n'))))

    # Assign an instrument to each project using ShowStat()
    table=pd.DataFrame()
    for instrument in list(filter(None, instruments_dict.values())):
        response = session_loggedin.get('https://ecopart.obs-vlfr.fr/statsample?MapN=&MapW=&MapE=&MapS=&filt_fromdate=&filt_todate=&filt_instrum={}&filt_proftype=&filt_depthmin=&filt_depthmax=&XScale=I&TimeScale=R'.format(instrument))
        while not response.ok:
            time.sleep(5)
            response = session_loggedin.get('https://ecopart.obs-vlfr.fr/statsample?MapN=&MapW=&MapE=&MapS=&filt_fromdate=&filt_todate=&filt_instrum={}&filt_proftype=&filt_depthmin=&filt_depthmax=&XScale=I&TimeScale=R'.format(instrument))
        soup = BeautifulSoup(response.text, 'lxml')
        stats=Ecopart_format_stat_table(soup.find_all('table')[4]) #table 1: Sample statistics, table 2: Sample count, table 3: Pressure slice, table 4-6: project stats, table 7: list of taxonomic categories
        stats.columns=['Project_title','Instrument','Number_depth_cast','Number_timeseries','Number_plankton_cast','Percentage_validated','Contact_name','Particles_manager_name','Plankton_manager_name','Export_status','Plankton_export']
        stats.mask(stats.isna(),other='',inplace=True)
        stats['Contact_name']=stats.apply(lambda x: x.Particles_manager_name if x.Particles_manager_name==x.Particles_manager_name else x.Plankton_manager_name if x.Particles_manager_name==x.Particles_manager_name else x.Contact_name if x.Contact_name==x.Contact_name else '(',axis=1 )
        stats['Contact_email']=stats.apply(lambda x: x.Contact_name[x.Contact_name.find('(')+1:x.Contact_name.find(')')] if x.Contact_name.find(')')!=-1 else x.Particles_manager_name[x.Particles_manager_name.find('(')+1:x.Particles_manager_name.find(')')] if x.Particles_manager_name.find(')')!=-1 else x.Plankton_manager_name[x.Plankton_manager_name.find('(')+1:x.Plankton_manager_name.find(')')] if x.Plankton_manager_name.find(')')!=-1 else '',axis=1)#stats.Contact_name.apply(lambda name: name[name.find('(')+1:name.find(')')] if name.find(')')!=-1 else np.NaN)
        stats['Contact_name'] = stats.Contact_name.apply(lambda name: name.split('(')[0].capitalize())
        stats['Instrument']=stats.Instrument.str.upper()
        table = pd.concat([table, stats], axis=0)
        # The table 'Sample per taxonomy project' contains link to Ecotaxa project, but the project title may differ between Ecotaxa/Ecopart.
        #ecotaxa=Ecopart_format_stat_table(soup.find_all('table')[5])
        #ecotaxa['Ecotaxa_ID']=ecotaxa.Project.apply(lambda name: name[name.rfind('(')+1:name.rfind(')')].replace('https://ecotaxa.obs-vlfr.fr/prj/','') if name.find(')')!=-1 else np.NaN)
        #ecotaxa['Project_title'] = ecotaxa.Project.apply(lambda name: name[name.rfind('(') + 1:name.rfind(')')].replace('https://ecotaxa.obs-vlfr.fr/prj/','') if name.find(')') != -1 else np.NaN)
    df_projects=pd.merge(df_projects,table[['Project_title','Instrument','Contact_name','Contact_email','Percentage_validated','Export_status','Plankton_export']],how='left',on=['Project_title'])
    df_projects.mask(df_projects.isna(), other='', inplace=True)

    # Check for access right and get corresponding ecotaxa project ID using DoSearch(). Looping is too long so using request threads module requests_futures
    # https://github.com/ross/requests-futures
    df_projects['Sample_search']=df_projects['Project_source'].str.replace(r"https://ecopart.obs-vlfr.fr/",r"https://ecopart.obs-vlfr.fr/searchsample", regex=True)
    #samples=list(map(lambda url: requests.get(url).json()[0]['id'] if len(requests.get(url).json()) else None,df_projects.Sample_search.tolist()))
    rights={'VY':'True','YY':'True','YN':'True','VN':'False','NN':'False','NY':'False'} # Y/Nx: login necessary/unnecessary Vx: login unnecessary  xY: exportable xN: not exportable
    # YY: exportable if logged in (green), VY: exportable (login not necessary), VN: not exportable (red), YN: exportable if logged in but no plankton data (orange)
    # NN: not exportable (black), NY: not exportable (green)
    access={}
    access_code = {}
    samples = {}
    with sessions.FuturesSession(session=session_loggedin) as session: #replace by sessions.FuturesSession(session_loggedin) if necessary
        futures = [session.get(url) for url in df_projects.Sample_search.tolist()]
        for future in as_completed(futures):
            resp = future.result()
            project_code=list(set(list(map(lambda info: info['visibility'],resp.json()))))[0] if len(resp.json()) else 'NN'
            project_access=rights[list(set(list(map(lambda info: info['visibility'],resp.json()))))[0]] if len(resp.json()) else False
            first_sample=resp.json()[0]['id'] if len(resp.json()) else ''
            samples[df_projects.query('Sample_search=="{}"'.format(resp.url)).Project_ID.values[0]]=first_sample
            access[df_projects.query('Sample_search=="{}"'.format(resp.url)).Project_ID.values[0]] = project_access
            access_code[df_projects.query('Sample_search=="{}"'.format(resp.url)).Project_ID.values[0]] = project_code
    #print(time.time() - start)
    df_projects=pd.merge(df_projects,pd.DataFrame.from_dict(samples,columns=['first_sample'],orient='index'),how='left',left_on='Project_ID',right_index=True)
    df_projects = pd.merge(df_projects, pd.DataFrame.from_dict(access, columns=['PSSdb_access'], orient='index'),how='left', left_on='Project_ID', right_index=True)
    df_projects = pd.merge(df_projects, pd.DataFrame.from_dict(access_code, columns=['Access_code'], orient='index'),how='left', left_on='Project_ID', right_index=True)
    print(len(df_projects) , "projects found.",len(df_projects[df_projects.PSSdb_access.astype(str).str.capitalize()=='True']), "projects accessible", sep=' ')
    df_projects['PSSdb_access']=df_projects['PSSdb_access']=='True'
    # Get corresponding Ecotaxa project ID from first sample popover infos
    df_projects['Sample_url']=df_projects.first_sample.apply(lambda sample: 'https://ecopart.obs-vlfr.fr/getsamplepopover/{}'.format(sample))
    project_id = dict(filter(None,map(lambda item: (info:=Ecopart_get_sample_popover(sample_url=item[1]),(df_projects.at[item[0],'Project_ID'],info['EcotaxaProject'][(info['EcotaxaProject'].rfind('(') + 1):info['EcotaxaProject'].rfind(')')]) if Ecopart_get_sample_popover(sample_url=item[1]) else None)[-1],enumerate(df_projects.Sample_url))))
    df_projects = pd.merge(df_projects, pd.DataFrame.from_dict(project_id, columns=['Ecotaxa_ID'], orient='index'),how='left', left_on='Project_ID', right_index=True)
    df_projects.mask((df_projects=='None') | (df_projects==''),inplace=True,other=np.NaN)
    #test=df_projects.dropna(subset=['Project_title','Project_ID','Ecotaxa_ID'])
    df_projects['Project_localpath']=localpath
    df_projects['Project_test']=df_projects['Ecotaxa_ID'].isin(['3315', '3318', '3326', '377', '378', '714', '560', '579', '5693']).astype(str)  # Add boolean variable to identify fixed test set (3 projects per instrument)
    df_projects['Project_ID'] =df_projects['Project_ID'].astype(int)
    df_projects=df_projects.sort_values(['PSSdb_access','Project_ID'],ascending=[False,True])
    df_projects=df_projects[['Project_source','Project_localpath','Project_ID','Project_title','Instrument','Contact_name','Contact_email','Access_code','Plankton_export','PSSdb_access','Project_test','Ecotaxa_ID','Percentage_validated']]
    df_metadata = pd.DataFrame({'Variables': df_projects.columns, 'Variable_types': df_projects.dtypes,
                                'Units/Values': ['','', '', '', '','','','YY: particles/plankton exportable if logged in (green), VY: particles/plankton exportable (login not necessary), VN: not exportable (red), YN: particles exportable if logged in (orange), NN: not exportable (black), NY: not exportable (green)','','','','',''],
                                'Description': ['Project_source (URL where original project files can be exported)','Local path indicating the location of exported files storage','Project ID', 'Project title','Project instrument','Name of the project contact','Email of the project contact','Internal code used to define project access rights','Plankton annotations accessibility on Ecopart','Project accessibility. If True, export is possible with current authentication','Test set boolean identifier. If True, project is one of the test set projects','Project ID on Ecotaxa','Percentage of validated images']})
    return {'data': df_projects, 'metadata': df_metadata}


def IFCB_dashboard_list(localpath):
    print('Searching for projects on CALOOS (https://ifcb.caloos.org/) and WHOI-IFCB (https://ifcb-data.whoi.edu/)')
    contacts = ['Clarissa Anderson', 'Heidi Sozik']
    emails = ['cra002@ucsd.edu', 'hsosik@whoi.edu']
    urls = ['https://ifcb.caloos.org/', 'https://ifcb-data.whoi.edu/']

    # lists to be filled with the project
    dashboard_url = [] # necessary to obtain metadata for each file in IFCB_dashboard_export function on the funcs_export_projects module
    Project_ID = []
    Project_title = []
    Project_source = []
    Project_localpath = []
    Instrument = []
    Contact_name = []
    Contact_email = []
    Pssdb_access = []
    Objects_number = [] # object count will not be consdiered here since it can be done together with the downloads
    Percentage_classified = []
    Percentage_validated = []
    Latest_update = []
    Project_test = []
    # n_roi = 0

    for ds, url in enumerate(urls):
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
        keys_proj_numbers = [key for key, val in data_dict.items() if
                                     any('timeline?dataset=' in s for s in val)]
        for n in keys_proj_numbers:
            for i in data_dict[n]:
                proj_info = re.search('timeline\?dataset=(.*)', i)
                try:
                    l = proj_info.group(1).split('>')
                    l[0] = l[0].replace('\'', '').replace('\"', '')
                    l[1] = l[1].replace('</a', '')
                    if l[0] not in Project_ID:
                        dashboard_url.append(urls[ds])
                        Project_ID.append(l[0])
                        Project_title.append(l[1])
                        Project_localpath.append(localpath)
                        Project_source.append(urls[ds] + 'timeline?dataset=' +l[0])
                        Instrument.append('IFCB')
                        Contact_name.append(contacts[ds])
                        Contact_email.append(emails[ds])
                        Pssdb_access.append('True')
                        Objects_number.append(float('nan'))
                        Percentage_classified.append(100.000)
                        Percentage_validated.append(float('nan'))
                        if (l[0] == 'santa-cruz-municipal-wharf' or l[
                            0] == 'EXPORTS'):  # test projects defined here
                            Project_test.append('True')
                        else:
                            Project_test.append('False')
                        # step here to get ROI count for each projects:
                        try:
                            pid_df = get_df_list_IFCB(urls[ds], l[0])
                            # print('obtaining data for project ' + l[0])
                            # print('counting the number of ROIs for '+ l[0])
                            # for i in tqdm(pid_df.loc[:, 'pid']):
                            # n_roi = n_roi + roi_number(urls[ds], i)
                            # Objects_number.append(n_roi)
                            Latest_update.append(pid_df.loc[len(pid_df) - 1, 'sample_time'])
                        except:
                            # Objects_number.append(float('nan'))
                            Latest_update.append(float('nan'))

                except:
                    pass

    df_projects = pd.DataFrame()
    df_projects['dashboard_url'] = dashboard_url
    df_projects['Project_source']= Project_source
    df_projects['Project_localpath'] = Project_localpath
    df_projects['Project_ID'] = Project_ID
    df_projects['Project_title'] = Project_title
    df_projects['Instrument'] = Instrument
    df_projects['Contact_name'] = Contact_name
    df_projects['Contact_email'] = Contact_email
    df_projects['PSSdb_access'] = Pssdb_access
    df_projects['Objects_number'] = Objects_number
    df_projects['Percentage_classified'] = Percentage_classified
    df_projects['Percentage_validated'] = Percentage_validated
    df_projects['Latest_update'] = Latest_update
    df_projects['Project_test'] = Project_test

    df_projects = df_projects[df_projects['Latest_update'].notna()]
    #os.remove(localpath)
    df_metadata = pd.DataFrame({'Variables': df_projects.columns, 'Variable_types': df_projects.dtypes,
                                'Units/Values': ['', '', '', '', '', '', '', '', '', '#', '%', '%', '', ''],
                                'Description': ['URL of IFCB dashboard, necessary to obtain metadata for the IFCB exported files',
                                                'Project_source (URL where original project files can be exported)',
                                                'Local path indicating the location of exported files storage',
                                                'Project ID', 'Project title', 'Project instrument',
                                                'Name of the project contact', 'Email of the project contact',
                                                'Project accessibility. If True, export is possible with current authentication',
                                                'Total number of objects (images) in the project. IFCB dashboard downloads count objects in the export step ',
                                                'Percentage of predicted images', 'Percentage of validated images',
                                                'Latest date the dataset has been updated',
                                                'Test set boolean identifier. If True, project is one of the test set projects']})

    return {'data': df_projects, 'metadata': df_metadata}
