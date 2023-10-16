##  Objective: This file contains one function to:
# (1) print the fields included in Ecotaxa projects using the API (downloading project is thus not necessary). This should facilitate the completion of the standardizer spreadsheets
# (1bis) print the list of standard units from pint module
# (2) save interactive plot of raw ecotaxa export file to identify outliers, missing values, and units. This should facilitate the completion of the standardizer spreadsheets
# (3) flag samples to be left out of standardized projects based on 5 criteria (missing necessary data/metadata, location on land, low ROI count, high number of bubbles and other artefacts, multiple size calibration factors) and generate corresponding interactive report
# (4) perform EcoTaxa export files standardization (standard column names and units) and harmonization (across instruments) based on standardizer spreadsheets
#  Running function (4) is required after project export in order to reduce (i.e. select only the raw columns necessary for further processing), harmonize (all the PSSdb projects now have similar columns, regardless of specific instrumentation and image acquisition/processing routines), and standardize (All columns follow standardized notations and units according to international database, including taxonomic annotations) project export files

## Requirements: Change path on
# l.46: configuration file with username/password information (configuration_masterfile_pw.yaml)
# l.93: custom units definition file (units_def.txt)
# l.106: plankton taxonomic information (plankton_annotated_taxonomy.xlsx). File created/updated below

## TO DO: Perform standardization on other imaging sensor projects (at the moment function has been tested on ZooCAM, generic scanner, Zooscan, IFCB, and UVP)


## Python modules
# Path modules
from pathlib import Path # Handling of path object

# Config modules
import wget
import yaml # requires installation of PyYAML package
# Multiprocessing mapping
from multiprocessing.pool import ThreadPool # Use: pip install multiprocess

# Prompt for export confirmation
import sys
# Ecotaxa API modules. Install module via git first: pip install git+https://github.com/ecotaxa/ecotaxa_py_client.git
import ecotaxa_py_client
from ecotaxa_py_client.api import projects_api
from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.api import samples_api
from ecotaxa_py_client.model.login_req import LoginReq
from ecotaxa_py_client.api import objects_api
from ecotaxa_py_client.model.project_filters import ProjectFilters
from ecotaxa_py_client.api import taxonomy_tree_api
from ecotaxa_py_client.model.taxon_model import TaxonModel
from collections import ChainMap
path_to_config=Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_to_git=Path(cfg['git_dir']).expanduser()
path_to_project_list=path_to_git / cfg['proj_list']

path_to_config_usr=Path('~/GIT/PSSdb/scripts/configuration_masterfile_pw.yaml').expanduser()
with open(path_to_config_usr ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

with ecotaxa_py_client.ApiClient() as client:
    api = authentification_api.AuthentificationApi(client)
    token = api.login(LoginReq(username=cfg_pw['ecotaxa_user'],password=cfg_pw['ecotaxa_pass']))
configuration = ecotaxa_py_client.Configuration(host = "https://ecotaxa.obs-vlfr.fr/api",access_token=token, discard_unknown_keys=True)
configuration.verify_ssl=False

# Panda and numpy modules (used for dataframes and arays)
import pandas as pd
import numpy as np
from natsort import natsorted
from itertools import compress
from datetime import datetime

# Globe data modules
import geopandas as gpd
# This dataset includes EEZ, using GOAS shapefile instead
#path_to_zip = Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent / 'ne_10m_admin_0_countries.zip'
path_to_zip = Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent /'GOaS_v1_20211214.zip'# 'ne_10m_land.zip'
path_to_longhurst=Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent /'longhurst_v4_2010.zip'

# Download high resolution countries/oceans shapefile
import requests
import shutil

if path_to_longhurst.with_suffix('').is_dir()==False:
     with requests.Session() as sess:
        url ='https://www.marineregions.org/download_file.php?name=longhurst_v4_2010.zip'
        params={'name':'NOAA','organisation':'NOAA','email':cfg_pw['ecotaxa_user'],'country':'United States','user_category':'academia','purpose_category':'Research','agree':"1",'submit':'Download'}
        rsp=sess.post(url,data=params, headers={'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.0.0 Safari/537.36'})
        with open(path_to_longhurst, "wb") as fd:
            for a_chunk in rsp.iter_content():  # Loop over content, i.e. eventual HTTP chunks
                # rsp.raise_for_status()
                fd.write(a_chunk)
        shutil.unpack_archive(path_to_longhurst, path_to_longhurst.with_suffix(''))  # Unzip export file
        path_to_longhurst.unlink(missing_ok=True)

if path_to_zip.with_suffix('').is_dir()==False:
  with requests.Session() as sess:
    url = 'https://www.marineregions.org/download_file.php?name=GOaS_v1_20211214.zip'#'https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip'
    #rsp = sess.get(url, headers={'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.0.0 Safari/537.36',}, stream=True)
    params={'name':'NOAA','organisation':'NOAA','email':cfg_pw['ecotaxa_user'],'country':'United States','user_category':'academia','purpose_category':'Research','agree':"1",'submit':'Download'}
    rsp=sess.post(url,data=params, headers={'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.0.0 Safari/537.36'})
    with open(path_to_zip, "wb") as fd:
        for a_chunk in rsp.iter_content():  # Loop over content, i.e. eventual HTTP chunks
            # rsp.raise_for_status()
            fd.write(a_chunk)
    shutil.unpack_archive(path_to_zip, path_to_zip.with_suffix(''))  # Unzip export file
    path_to_zip.unlink(missing_ok=True)

oceans = gpd.read_file(list(Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent.rglob('goas_v01.shp'))[0])
Longhurst = gpd.read_file(list(Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent.rglob('Longhurst_world_v4_2010.shp'))[0])

from shapely.geometry import Polygon, mapping, Point, MultiPolygon
from shapely import wkt
from shapely.ops import unary_union

# Plot module
import plotly.express as px # Use pip install plotly==5.8.0
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotnine import * # Python equivalent to ggplot2. Use pip install plotnine. Do not execute in Pycharm (bug- no fix yet): https://youtrack.jetbrains.com/issue/PY-42390
from colorspace import sequential_hcl # Use: pip install git+https://github.com/retostauffer/python-colorspace


import warnings
warnings.filterwarnings('ignore', module='urllib3')

# Unit conversion
import os
import pint # Use pip install pint
warnings.filterwarnings("ignore")
warnings.filterwarnings('ignore', module='pint')
import pint_pandas # Use pip install pint-pandas
PA=pint_pandas.PintArray
from pint import UnitRegistry # Use pip install pint to install
ureg=UnitRegistry()
PQ=ureg.Quantity
#ureg.default_format = '~' # Add this if you want to use abbreviated unit names.
ureg.load_definitions(Path(cfg['git_dir']).expanduser()/cfg['units_file'])  # This text file is used to define custom units from standard units  (e.g. square_pixel etc.)
full_list_units = list(dict.fromkeys(sorted(dir(ureg))))  # list(dict.fromkeys(sorted(list(np.concatenate([dir(getattr(ureg.sys,system)) for system in dir(ureg.sys)]).flat))))
import re
import json
# Poisson distribution
from scipy.stats import poisson # Estimate counts uncertainties assuming detection follows a Poisson distribution
import math

# NBSS computation
try:
    from script_NBSS_datacheck import *
except:
    from scripts.script_NBSS_datacheck import *

import ast
from tqdm import tqdm

# Standardization of taxonomic annotations
try:
    from funcs_standardize_annotations import taxon_color_palette,annotation_in_WORMS
except:
    from scripts.funcs_standardize_annotations import taxon_color_palette, annotation_in_WORMS
import ast
from tqdm import tqdm

path_to_taxonomy=Path(Path.home()/'GIT'/'PSSdb'/cfg['annotations_lookup']).expanduser()
if not path_to_taxonomy.exists():
    print ('Creating a taxonomic annotation standardization spreadsheet using the World Register of Marine Species (https://www.marinespecies.org/). Please wait')
    # 1) Retrieve fields of interest (including taxonomic hierarchies) in EcoTaxa projects
    # see: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ObjectsApi.md#get_object_set
    fields_of_interest = "txo.id,txo.display_name";
    colnames = ['ROI_Annotation_ID', 'Category']
    with ecotaxa_py_client.ApiClient(configuration) as api_client:
        api_project_instance = projects_api.ProjectsApi(api_client)
        project = api_project_instance.search_projects(also_others=False,  # Return visible projects with access rights
                                                       title_filter='',  # Optional title project filter
                                                       instrument_filter='',  # All instruments
                                                       order_field='projid'# Sorting variable. Use instrument or projid
                                                       )
        project_list = list(map(lambda x: str(x.projid), project))
        api_object_instance = objects_api.ObjectsApi(api_client)
        objects = api_object_instance.get_object_set(project_list[0], ProjectFilters(statusfilter="PVD"),fields=fields_of_interest)
        df = pd.DataFrame(objects.details, columns=colnames).drop_duplicates()
        for project in project_list[1:]:
            objects = api_object_instance.get_object_set(project, ProjectFilters(statusfilter="PVD"),fields=fields_of_interest)
            df = pd.concat([df, pd.DataFrame(objects.details, columns=colnames)], axis=0).drop_duplicates()
        # df['ROI_ID']=df['ROI_ID'].apply(lambda x: ''.join(re.split('([0-9]+)', x)[0:-2]) + re.split('([0-9]+)', x)[-2].zfill(6))
        api_taxo_instance = taxonomy_tree_api.TaxonomyTreeApi(api_client)
        df['EcoTaxa_hierarchy'] = ['>'.join(api_taxo_instance.query_taxa(int(id)).lineage[::-1]) for id in df.ROI_Annotation_ID]
        df = df.sort_values(by=['EcoTaxa_hierarchy'])
        df = df.reset_index()
    # 2) Post/get annotation on World Register of Marine Species (https://www.marinespecies.org/)
    # Loop through EcoTaxa hierarchy and merge to annotated taxonomy
    df_taxonomy = pd.DataFrame({})
    with tqdm(desc='', total=len(df), bar_format='{desc}{bar}', position=0, leave=True) as bar:
        print('Searching for taxon in World Register of Marine Species ({}):\n'.format('https://www.marinespecies.org/'))
        try:
            for hierarchy in df.EcoTaxa_hierarchy:
                # print(hierarchy)
                # Update annotated taxonomy and progress bar
                bar.set_description(hierarchy, refresh=True)
                # and update progress bar
                ok = bar.update(n=1)
                data = annotation_in_WORMS(hierarchy)
                df_taxonomy = pd.concat([df_taxonomy, data],axis=0)
        except KeyboardInterrupt:  # Press escape to exit loop safely
            print('{} not found. Skipping annotation'.format(hierarchy))
    # 3) Assign additional variables
    df_taxonomy['URL'] = df_taxonomy.WORMS_ID.apply(lambda id: 'https://www.marinespecies.org/aphia.php?p=taxdetails&id={}'.format(id.replace('urn:lsid:marinespecies.org:taxname:', '')) if len(id) else '')
    df_taxonomy['Life_stage'] = df_taxonomy.Functional_group.apply(lambda group: ';'.join([ast.literal_eval(dict)['Life stage'] for dict in group.split(';') if len(group) > 0 and 'Life stage' in ast.literal_eval(dict).keys()]))
    df_taxonomy['functional_group'] = df_taxonomy.Functional_group.apply(lambda group: ';'.join([dict.replace('{','').replace(', ', ' (').replace('}', ')').replace("'", "") if len(group) > 0 and len(ast.literal_eval(dict)) > 1 else dict.replace('{', '').replace(', ', ' (').replace('}', '').replace("'","") if len(group) > 0 and len(ast.literal_eval(dict)) == 1 else '' for dict in group.split( ';')]))
    df_taxonomy = pd.merge(df[['Category', 'EcoTaxa_hierarchy']], df_taxonomy, on='EcoTaxa_hierarchy', how='right')
    df_taxonomy = df_taxonomy[['Category', 'EcoTaxa_hierarchy', 'Full_hierarchy', 'Rank', 'Type', 'Domain', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Functional_group', 'functional_group', 'Life_stage', 'WORMS_ID', 'Reference', 'Citation','URL']]
    df_taxonomy = df_taxonomy.sort_values(['Type', 'Category'], ascending=[False, True]).reset_index(drop=True)
    # 4) Save data and metadata
    df_taxonomy_metadata = pd.DataFrame({'Variables': df_taxonomy.columns, 'Variable_types': df_taxonomy.dtypes,
                                          'Description': ['Region of interest (object) annotation category',
                                                          'Full hierarchy of annotation category',
                                                          'Full taxonomic hierarchy in World Register of Marine Species (WORMS)',
                                                          'Taxonomic rank of the annotation category in WORMS',
                                                          'Type of particle for the annotation category (e.g. Living)',
                                                          'Taxonomic domain/kingdom of the annotation category',
                                                          'Taxonomic phylum of the annotation category',
                                                          'Taxonomic class of the annotation category',
                                                          'Taxonomic order of the annotation category',
                                                          'Taxonomic family of the annotation category',
                                                          'Taxonomic genus of the annotation category',
                                                          'Functional group of the annotation category (dictionary format)',
                                                          'Functional group of the annotation category (string format)',
                                                          'Life stage of the annotation category in WORMS',
                                                          'Unique ID of the annotation category in the WORMS database',
                                                          'Reference for the annotation category description',
                                                          'Citation for the annotation category in WORMS',
                                                          'URL of the annotation category in WORMS']})
    with pd.ExcelWriter(str(path_to_taxonomy), engine="xlsxwriter") as writer:
        df_taxonomy.to_excel(writer, sheet_name='Data', index=False)
        df_taxonomy_metadata.to_excel(writer, sheet_name='Metadata', index=False)
    print('Taxonomic annotations lookup table saved: {}'.format(str(path_to_taxonomy)))
metadata=pd.read_excel(path_to_taxonomy,sheet_name='Metadata')
columns_types={metadata.Variables[i]:metadata.Variable_types[i] for i in metadata.index}
df_taxonomy=pd.read_excel(path_to_taxonomy,dtype=columns_types )
rank_categories=[rank.partition("(")[2].replace(')','') for hierarchy in pd.Series(df_taxonomy.loc[df_taxonomy['Full_hierarchy'].astype(str).str.len().nlargest(1, keep="all").index[0],'Full_hierarchy']) for rank in hierarchy.split('>')]
df_taxonomy['Rank']=pd.Categorical(df_taxonomy.Rank,ordered=True,categories=rank_categories)

# Functions start here:

# Function (1): Use EcoTaxa API to search for project fields automatically without downloading export files
def filling_standardizer_field_func(config_path,project_id):
    """
    Objective: This function uses EcoTaxa API to print the free fields in accessible projects.
    This should facilitate the completion of projects standardizer spreadsheets.
    :param config_path: Full path of the configuration file containing Ecotaxa authentication info
    :param project_id: integer for specific project ID
    :return: free fields for all accessible projects or specific project ID
    """
    kwargs = project_id
    project_id = kwargs.get('project_id')

    # Open configuration file to get: EcoTaxa authentication info
    path_to_config_usr = Path(config_path).expanduser()
    with open(path_to_config_usr, 'r') as config_file:
        cfg_pw = yaml.safe_load(config_file)
    # Search visible and accessible projects based on instrument filter (Authentication in EcoTaxa required)
    with ecotaxa_py_client.ApiClient() as client:
        api = authentification_api.AuthentificationApi(client)
        token = api.login(LoginReq(
            username=cfg_pw['ecotaxa_user'],
            password=cfg_pw['ecotaxa_pass']
        ))

    configuration = ecotaxa_py_client.Configuration(host="https://ecotaxa.obs-vlfr.fr/api", access_token=token,discard_unknown_keys=True)
    with ecotaxa_py_client.ApiClient(configuration) as api_client:
        api_instance = projects_api.ProjectsApi(api_client)
        try:
            # Search Projects and return fields
           api_response_project_search_accessible = api_instance.search_projects(also_others=False,# Return accessible projects given account rights
                                                                                  title_filter='',# Optional title project filter
                                                                                  instrument_filter='',# All instruments
                                                                                  order_field='projid'# Sorting variable. Use instrument or projid
                                                                                  )

        except ecotaxa_py_client.ApiException as e:
            print("Exception when calling ProjectsApi->search_projects: %s\n" % e)
    # Step 2: Generate a dictionary for individual project with variables of interest (i.e. free fields)
    subset_df_standardizer=pd.read_excel(standardizer_path, usecols=["Project_ID"])
    project_fields=dict(map(lambda x: (x.projid,["acq_" + str for str in list(x.acquisition_free_cols.keys())] +
                                      ["object_" + str for str in list(x.obj_free_cols.keys())] +
                                      ["sample_" + str for str in list(x.sample_free_cols.keys())] +
                                      ["process_" + str for str in list(x.process_free_cols.keys())]),
                     list(filter(lambda w: any(subset_df_standardizer['Project_ID'].isin([w.projid])),api_response_project_search_accessible))))

    if not project_id is None:
        if int(project_id) in list(project_fields.keys()):
            project_fields = {int(project_id): project_fields[int(project_id)]}


    print(yaml.dump(project_fields))
    return project_fields

# Function (1b): Print possible units
def filling_standardizer_unit_func(custom_units_path=Path(cfg['git_dir']).expanduser()/cfg['units_file']):
    """
    Objective: This function print the units available to fill out the standardizer spreadsheets.
    :param custom_units_path: Full path of the custom unit definition file
    :return: units
    """
    ureg = UnitRegistry()
    # ureg.default_format = '~' # Add this if you want to use abbreviated unit names.
    ureg.load_definitions(Path(custom_units_path).expanduser())  # This text file is used to define custom units from standard units  (e.g. square_pixel etc.)
    full_list_units = list(dict.fromkeys(sorted( dir(ureg))))  # list(dict.fromkeys(sorted(list(np.concatenate([dir(getattr(ureg.sys,system)) for system in dir(ureg.sys)]).flat))))
    print('Use unit from the list below or define custom units in {}:\n'.format(custom_units_path),full_list_units,sep='')


# Function (2): Generate interactive plot of variable whose units need to be checked automatically
def diagnostic_plots_func(standardizer_path,project_id):
    """
       Objective: This function uses standardizer spreadsheets to plot variable whose units need to be checked interactively for a unique project.
       This should facilitate (1) datapoints flagging and filtering (2) the completion of projects standardizer spreadsheets.
       :param standardizer_path: Full path of the standardizer spreadsheet containing project ID of interest
       :param project_id: unique integer for project ID to be plotted
       :return: print interactive plots
       """

    path_to_git = Path('~/GIT/PSSdb').expanduser()

    df_standardizer = pd.read_excel(standardizer_path,index_col=0)
    path_to_plotdir = path_to_git / cfg['figures_subdir'] / 'Standardizer'
    path_to_plotdir.mkdir(parents=True, exist_ok=True)

    path_to_data =Path(df_standardizer.at[project_id,"Project_localpath"]).expanduser()# path_to_git / Path(cfg['dataset_subdir']) / df_standardizer.at[project_id,"Instrument"]
    columns_of_interest =   ["Station_ID","Sample_ID","Cruise_field"]+[string.replace("_unit", "_field") for string in list(df_standardizer.columns) if "_unit" in string]
    fields_of_interest =dict(zip(columns_of_interest,["sample_stationid","sample_id"]+df_standardizer.loc[project_id][columns_of_interest[2:]].values.flatten().tolist()))
    fields_of_interest = dict((key.replace("_field",""),value) for key,value in fields_of_interest.items() if isinstance(value, str))
    path_to_datafile_list=list(path_to_data.glob('ecotaxa_export_{}*'.format(str(project_id))))
    path_to_datafile_list=[path for  path in path_to_datafile_list if 'flag' not in str(path)]
    df=pd.read_table(path_to_datafile_list[0],usecols=fields_of_interest.values())
    old_columns=df.columns
    df.columns = [(list(fields_of_interest.keys())[list(fields_of_interest.values()).index(value)]) for value in old_columns] #pd.MultiIndex.from_arrays([[list(fields_of_interest.keys())[list(fields_of_interest.values()).index(value)] for value in df.columns],df.columns])
    spatial_bins=1 # Fixed spatial bining 1x1 degrees
    vertical_bins=5 # Fixed vertical bin 5m
    df=df.assign(Latitude_group=pd.cut(df.Latitude,np.arange(-90,90,step=spatial_bins)),
                 Longitude_group=pd.cut(df.Longitude,np.arange(-180,180,step=spatial_bins)),
                 Depth_med=np.add(df.Depth_min,df.Depth_max)/2,
                 Depth_group=pd.cut(np.add(df.Depth_min,df.Depth_max)/2,np.arange(0,6000,step=vertical_bins),labels=np.arange(0,6000,step=vertical_bins)[:-1]))
    df=pd.merge(df,df.groupby(['Sample_ID']).agg(objects=(list(fields_of_interest.keys())[9],'count')).reset_index(),on=['Sample_ID'])
    sample_df=df.groupby(['Latitude_group','Longitude_group','Depth_med'],observed=True).apply(lambda x:pd.Series({'n_objects':np.nanmean(x.objects),
                                                                                                                   'n_objects_upper':np.nanmean(x.objects)+np.nanstd(x.objects),
                                                                                                                   'n_objects_lower': np.nanmean(x.objects)-np.nanstd(x.objects),
                                                                                                                   'std_objects':np.nanstd(x.objects)
                                                                                                       })).reset_index()
    profiles_df=df.groupby(['Latitude','Longitude','Station_ID','Depth_group'],observed=True).apply(lambda x:pd.Series({'n_objects':np.nanmean(x.objects)})).reset_index()
    df=pd.merge(df,sample_df,on=['Latitude_group','Longitude_group','Depth_med'])

    plot=(ggplot(sample_df) +
                geom_point(mapping=aes(x='n_objects', y='std_objects'),color='#222222') +
                geom_abline(mapping=aes(intercept=0,slope=0.5))+
                labs(title=r'Count uncertainties in spatial bins ({}x{}$^\circ$)'.format(str(spatial_bins),str(spatial_bins)), x='Number of objects(images) per sample in spatial bins',y='Std of the number of objects per sample in spatial bins') +
                scale_y_log10(labels=lambda y: list(map(lambda x: (("{0:.{1:d}e}".format(x, 0)).split("e"))[0] + "x10$^{" +(("{0:.{1:d}e}".format(x, 2)).split("e"))[1] + "}$", y))) +
                scale_x_log10(labels=lambda y: list(map(lambda x: (("{0:.{1:d}e}".format(x, 0)).split("e"))[0] + "x10$^{" +(("{0:.{1:d}e}".format(x, 2)).split("e"))[1] + "}$", y))) +
                theme(panel_grid=element_blank(), panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'), axis_title=element_text(size=10),
                      axis_text_x=element_text(size=10), axis_text_y=element_text(size=10, rotation=90),
                      plot_background=element_blank()))
    #print(plot)
    subset_df=df.drop_duplicates(subset=list(fields_of_interest.keys())[0:8])
    subset_df=subset_df.assign(Area=subset_df.Area if 'Area' in subset_df.columns else None,
                     Biovolume=subset_df.Biovolume if 'Biovolume' in subset_df.columns else None,
                     Minor_axis=subset_df.Minor_axis if 'Minor_axis' in subset_df.columns else None)

    subset_df =subset_df.sort_values('Station_ID')
    #subset_df = df.groupby(by=list(fields_of_interest.keys())[0:8],dropna=False).apply(lambda x:pd.Series({'n_objects':len(x[list(fields_of_interest.keys())[9]]),'Volume_analyzed':x.Volume_analyzed.unique()[0] if 'Volume_analyzed' in x.columns else None,
    #                                                                                                       'Area':x.Area.unique()[0] if 'Area' in x.columns else None,'Minor':x.Minor_axis.unique()[0] if 'Minor_axis' in x.columns else None,'ESD':x.ESD.unique()[0] if 'ESD' in x.columns else None,'Biovolume':x.Biovolume.unique()[0] if 'Biovolume' in x.columns else None,'Pixel':x.Pixel.unique()[0] if 'Pixel' in x.columns else None})).reset_index()
    cols = subset_df.columns
    fig = make_subplots(rows=2, cols=2,#shared_xaxes='columns', # Use shared x/yaxes to fix axis
                       # subplot_titles=("","Number of objects", "Volume analyzed", "Pixel:size ratio", "Area","Minor ellipsoidal axis","Equivalent spherical diameter","Biovolume"),
                      #  specs=[[{"type": "scattergeo", "rowspan": 7}, {"type": "scatter"}],[None, {"type": "scatter"}],[None, {"type": "scatter"}],[None, {"type": "scatter"}],[None, {"type": "scatter"}],[None, {"type": "scatter"}],[None, {"type": "scatter"}]],
                        specs=[[{"type": "scattergeo", "rowspan": 2}, {"type": "scatter"}],[None, {"type": "scatter"}]],
                       column_widths=[0.5,0.5],vertical_spacing=0.07)
    # Map, left panel
    fig.add_trace(go.Scattergeo(lon=subset_df.Longitude, lat=subset_df.Latitude,
                                  hovertext='Sample ID: '+subset_df.Sample_ID, hoverinfo="lon+lat+text",
                                marker=dict(color='black'),geojson="natural earth",showlegend=False),
                   row=1,col=1).update_geos(projection_rotation=dict(lon=np.nanmean(subset_df.Longitude), lat=np.nanmean(subset_df.Latitude), roll=0),center=dict(lon=np.nanmean(subset_df.Longitude), lat=np.nanmean(subset_df.Latitude)),projection_type='orthographic')
    # Scatterplot 1, top-right panel
    fig.add_trace(go.Scatter(x=subset_df.Station_ID, y=subset_df.n_objects,
                             hovertext="Sample ID: "+subset_df.Sample_ID,hoverinfo="text",marker=dict(color='black'),
                             showlegend=False,visible=True),row=1,col=2) # Use mode='markers' for points only
    # Scatterplot 2, bottom-right panel
    profiles_df=profiles_df.sort_values(['Depth_group','Station_ID'])
    colors=sequential_hcl(h = [-220,20], c = [0, 90, 0], l = [95, 30], power = [1.5,1.5], rev = True)(len(subset_df['Station_ID'].unique()))

    fig.add_trace(go.Box(y=df.Depth_group, x=df.n_objects,
                         fillcolor='black',line=dict(color='#222222'),
                         hoverinfo='skip',
                         # error_y=dict(type="data",array=subset_df.std_objects),
                         # hovertext="Sample ID: " + subset_df.Sample_ID, hoverinfo="text",
                         notched=True,
                         showlegend=False, visible=True), row=2, col=2)
    fig.update_traces(orientation='h',row=2,col=2)

    for i,station in enumerate(np.sort(profiles_df.Station_ID.unique())):
        subset_profiles= profiles_df[profiles_df['Station_ID']==station].sort_values(['Depth_group'])
        fig.add_trace(go.Scatter(y=subset_profiles.Depth_group, x=subset_profiles.n_objects,
                                 marker=dict(color=colors[i]), line=dict(color=colors[i]),
                                 name="Station "+str(station), hovertext=subset_profiles.Station_ID, hoverinfo="text",
                                 legendrank=i,
                                 # line_shape='spline',
                                 # error_y=dict(type="data",array=subset_df.std_objects),
                                 # hovertext="Sample ID: " + subset_df.Sample_ID, hoverinfo="text",
                                 mode='markers',showlegend=True, visible="legendonly"), row=2, col=2)





    fig.update_yaxes(type="linear", autorange="reversed", row=2, col=2)

    button_scatter1 = [dict(method="restyle",  # argument to change data
                            args=[{'x': [subset_df['Station_ID'], 'undefined'],
                                   'y': [subset_df[cols[k]], 'undefined'],
                                   'visible': [True, True, True]}, [1]], # visible should have the same length as number of subplots, second argument specifies which subplot is updated
                            label=cols[k]) for k in np.where(
        np.isin(list(cols), ['n_objects', 'Volume_analyzed', 'Area', 'Minor', 'ESD', 'Biovolume', 'Pixel']))[0]]


    fig.update_layout(updatemenus=list([dict(active=0, buttons=button_scatter1, x=0.87, y=1.1, xanchor='left', yanchor='top')]),
        title={'text': 'Cruise:' + subset_df["Cruise"][0] + "<br> Hover on datapoint to see sample ID",'xanchor': 'center', 'yanchor': 'top', 'x': 0.5},
        xaxis={'title': 'Station ID'},xaxis2={'title': 'Number of objects imaged'},
        yaxis={'title': 'Variable selected in dropdown menu <br>(fill unit in standardizer)'},yaxis2={'title':'Depth (m)'},
        hovermode="x",
       boxmode='group')

    fig.for_each_xaxis(lambda x: x.update(showgrid=False))
    fig.for_each_yaxis(lambda x: x.update(showgrid=False))
    fig.update_yaxes(type="log", row=1, col=2)


    fig.show()
    path_to_plot=path_to_plotdir / "Plot_standardizer_{}_{}.html".format(df_standardizer.at[project_id,"Instrument"],str(project_id))
    print("Saving plots to: ",path_to_plot,sep=" ")
    fig.write_html(path_to_plot)
    """  
    # Currently, API function get_object_set only accepts object related fields, not acquisition/process/sample fields 
    # Search visible and accessible projects based on instrument filter (Authentication in EcoTaxa required)
    with ecotaxa_py_client.ApiClient() as client:
        api = authentification_api.AuthentificationApi(client)
        token = api.login(LoginReq(
            username=cfg_pw['ecotaxa_user'],
            password=cfg_pw['ecotaxa_pass']
        ))
    configuration = ecotaxa_py_client.Configuration(host="https://ecotaxa.obs-vlfr.fr/api", access_token=token,
                                                        discard_unknown_keys=True)
    fields_of_interest = "obj.latitude,obj.longitude,obj.depth_min,obj.depth_max," + ','.join(['fre.' + item[item.find("_") + 1:] for item in fields_of_interest if isinstance(item, str) and 'object_' in item])

   # with ecotaxa_py_client.ApiClient(configuration) as api_client:
   #     api_instance =  objects_api.ObjectsApi(api_client)
   #     objects= api_instance.get_object_set(int(project_id),ProjectFilters(statusfilter="PVD"),fields=fields_of_interest)
   df=pd.DataFrame(objects['details'], columns=fields_of_interest.split(','))

   """

# Function (3): Flag project samples based on 5 criteria, update standardizer spreadsheet flagged sample ID, and generate automatic report
# converter function: Reformat the na values to match the type of individual columns
def is_float(na) -> bool:
    try:
        float(na)
        return True if na != 'nan' else False
    except ValueError:
        return False


def convert(na, type):
    try:
        return str(pd.DataFrame({'NA': na}, index=[0]).NA.astype(type)[0])
    except ValueError:
        str(na)


def flag_func(dataframe,count_uncertainty_threshold=cfg['count_uncertainty_threshold'],artefacts_threshold=cfg['artefacts_threshold'],validation_threshold=cfg['validation_threshold']):
    """
        Objective: This function assign flags to samples contained in a imaging dataset based on multiple criteria.
        Criteria include missing data/metadata, anomalous GPS location, low ROI counts/validation percentage, presence of artefacts, multiple size calibration factors
        :param dataframe: dataframe whose samples need to be quality controlled
        :param count_uncertainty_threshold: threshold used to flag samples with low ROI counts [0-1]. Default is 5% (0.05), change in comfiguration masterfile if needed
        :param artefacts_threshold: threshold used to flag samples with high percentage of artefacts [0-1]. Default is 20% (0.2), change in comfiguration masterfile if needed
        :param validation_threshold: threshold used to flag samples with low percentage of taxonomic annotation [0-1]. Default is 95% (0.95), change in comfiguration masterfile if needed
        :return: initial dataframe with flags, sample-specific summary
    """

    # Replace sample if missing
    if 'Sample' not in dataframe.columns and 'Profile' in dataframe.columns:
        dataframe['Sample']=dataframe['Profile']
    if 'Sample' in dataframe.columns:
        # Flag #1: Missing required data/metadata
        dataframe['Flag_missing'] = 0
        missing_flag = dataframe.groupby(['Sample']).apply(lambda x :pd.Series({'Sample':x.Sample.unique()[0],'Sample_missing':True if len(pd.isnull(x).any(1).to_numpy().nonzero()[0]) else False,'Missing_field':';'.join(list(pd.isnull(x).any(0).index[pd.isnull(x).any(0).to_numpy().nonzero()[0]])) if len(pd.isnull(x).any(1).to_numpy().nonzero()[0]) else ''})).reset_index(drop=True)# pd.isnull(dataframe).any(1).to_numpy().nonzero()[0]
        dataframe.loc[dataframe.Sample.isin(list(compress(missing_flag.Sample, missing_flag.Sample_missing))), 'Flag_missing'] = 1
        dataframe=pd.merge(dataframe, missing_flag[['Sample','Missing_field']],how='left',on='Sample')
        df_summary=dataframe[['Sample','Missing_field']].drop_duplicates().reset_index(drop=True)
        # Flag #2/3: anomalous GPS location
        dataframe['Flag_GPScoordinatesonland'] = 0
        dataframe[['Longitude', 'Latitude']]=np.round(dataframe[['Longitude', 'Latitude']],4)
        gdf = gpd.GeoDataFrame(dataframe[['Sample','Longitude', 'Latitude']].drop_duplicates().dropna(), geometry=gpd.points_from_xy(dataframe[['Sample','Longitude', 'Latitude']].drop_duplicates().dropna().Longitude, dataframe[['Sample','Longitude', 'Latitude']].drop_duplicates().dropna().Latitude))
        point_in_polygons = gpd.tools.sjoin(gdf, oceans, predicate="within", how='left')
        summary=pd.merge(point_in_polygons[['Sample', 'name']].rename(columns={'name':'Study_area'}),gpd.tools.sjoin(gdf, Longhurst, predicate="within", how='left')[['Sample', 'ProvDescr']].rename(columns={'ProvDescr':'Longhurst_province'}),how='left',on='Sample').drop_duplicates(subset=['Sample','Study_area','Longhurst_province'])
        df_summary=pd.merge(df_summary,summary,how='left',on=['Sample'])
        if any(point_in_polygons.index_right.astype(str) == 'nan'):
            gdf['Flag_GPScoordinatesonland'] =point_in_polygons.index_right.astype(str) == 'nan'# point_in_polygons.index_right.astype(str) != 'nan'
        else:
            gdf['Flag_GPScoordinatesonland'] = False
        dataframe['Flag_GPScoordinatesonland'] = np.where(dataframe['Sample'].isin(gdf[gdf['Flag_GPScoordinatesonland'] == True]['Sample'].tolist()), 1, 0)
        dataframe['Flag_dubiousGPScoordinates'] = 0
        dataframe['Flag_dubiousGPScoordinates'] = np.where(dataframe['Sample'].isin(gdf.get('Sample')[list(map(lambda x: x.contains(Point(0, 0)), gdf.geometry.tolist()))]), 1,dataframe['Flag_dubiousGPScoordinates'])
         # Flag #4 (count), #5 (artefacts), #6(validation, UVP/Zooscan-only)
        # Calculate count uncertainties assuming Poisson distribution. Upper-lower count limit are based on 5% uncertainty
        df_summary=pd.merge(df_summary,dataframe.groupby(['Sample']).apply(lambda x :pd.Series({'ROI_count':len(x.ROI) if 'ROI_number' not in dataframe.columns else x['ROI_number'].sum(),'Count_uncertainty':poisson.pmf(k=len(x.ROI),mu=len(x.ROI)) if 'ROI_number' not in dataframe.columns else poisson.pmf(k=x['ROI_number'].sum(),mu=x['ROI_number'].sum()),'Count_lower':poisson.ppf((0.05/2), mu=len(x.ROI)) if 'ROI_number' not in dataframe.columns else poisson.ppf((0.05/2), mu=x['ROI_number'].sum()),'Count_upper':poisson.ppf(1-(0.05/2), mu=len(x.ROI)) if 'ROI_number' not in dataframe.columns else poisson.ppf(1-(0.05/2), mu=x['ROI_number'].sum())})).reset_index(),how='left',on='Sample')
        # kwargs={'mapping':'x:Sample,y:ROI_count,ymin:Count_lower,ymax:Count_upper','scale_y_log10()':'','theme(axis_text_x=element_blank())':''}
        # ggplot_funcs(dataframe=summary,output_path='~/GIT/PSSdb_LOV_orga/Plots/EcoTaxa_3315_count_uncertainty.png',xlab='Sample ID',ylab='Number of pictures per sample',geom='geom_pointrange',**kwargs)
        dataframe['Flag_count'] = 0
        dataframe['Flag_count'] =np.where(dataframe['Sample'].isin(df_summary[df_summary.Count_uncertainty>count_uncertainty_threshold].Sample.tolist()),1,dataframe['Flag_count'])
        # dataframe'0' in dataframe['Flag_count'].astype('str')].Sample.unique()
        # Flag #5: Presence of artefacts (>20% of the total sample ROIs)
        dataframe['Flag_artefacts'] = 0
        dataframe['Flag_validation'] = 1 # Flagged by default to allow correct flagging of unclassified projects

        if ('Category' in dataframe.columns):
            dataframe.Category = np.where(dataframe.Category == '', pd.NA, dataframe.Category)
            if len(dataframe.dropna(subset=['Category'])):
                # Extract number of artefacts per samples
                df_summary =pd.merge(df_summary,dataframe.dropna(subset=['Category']).groupby(['Sample']).apply(lambda x: pd.Series({'Artefacts_count': len(x[x.Category.str.lower().apply(lambda annotation:len(re.findall(r'bead|bubble|artefact|artifact',annotation))>0)].ROI) if len(x.ROI) else 0,'Artefacts_percentage': len(x[x.Category.str.lower().apply(lambda annotation:len(re.findall(r'bead|bubble|artefact|artifact',annotation))>0)].ROI) / len(x.ROI) if  len(x.ROI) else 0,'Validation_percentage': len(x[x.Annotation.str.lower().apply(lambda annotation:len(re.findall(r'validated',annotation))>0)].ROI) / len(x.ROI) if  len(x.ROI) else 0})).reset_index(),how='left',on=['Sample'])
                dataframe['Flag_artefacts'] = np.where(dataframe['Sample'].isin(df_summary[df_summary.Artefacts_percentage > artefacts_threshold].Sample.tolist()),1, dataframe['Flag_artefacts'])
                dataframe['Flag_validation'] =np.where(dataframe['Sample'].isin(df_summary[df_summary.Validation_percentage >=validation_threshold].Sample.tolist()),0, dataframe['Flag_validation'])

        if validation_threshold==0:
            dataframe['Flag_validation'] = 0

        # dataframe['0' in dataframe['Flag_artefacts'].astype('str')].Sample.unique()
        # Flag #7: Multiple size calibration factors
        if 'Pixel' in dataframe.columns:
            dataframe['Flag_size']=0
            # Extract pixel size per sample/profile
            summary_pixel=dataframe.astype({'Cruise':str,'Sample':str,'Pixel':float})[['Cruise','Sample','Pixel']].drop_duplicates().groupby(['Cruise','Sample','Pixel']).apply(lambda x: pd.Series({'Sample_percentage':len(x.Sample.unique())//len(dataframe.Sample.unique()),'Profile_percentage':len(x.Profile.unique())/len(dataframe.Profile.unique()) if 'Profile' in x.columns else 1})).reset_index()
            summary_pixel=summary_pixel.sort_values(by=['Cruise','Profile_percentage'],ascending=[True,False]).reset_index()
            summary_subset =summary_pixel.groupby(by=['Cruise','Sample']).apply(lambda x: pd.Series({'Flag_size':1 if (len(summary_pixel[(summary_pixel.Cruise==x.Cruise.unique()[0])].Pixel.unique())>1) and (x.Pixel.astype('str').isin((summary_pixel[(summary_pixel.Cruise==x.Cruise.unique()[0])].groupby(['Cruise','Pixel']).agg({'Profile_percentage':'sum'}).reset_index().sort_values(by=['Cruise','Profile_percentage'],ascending=[True,False])).reset_index(drop=True).loc[1:,'Pixel'].astype('str').to_list()).values[0]) else 0})).reset_index()
            #pixel_list=summary.loc[1:,'Pixel'].astype('str').to_list() if len(summary.Pixel.values)>1 else ['inf']
            dataframe['Flag_size'] = np.where(dataframe['Sample'].astype('str').isin(list(summary_subset[summary_subset.Flag_size==1].Sample.unique())),1, dataframe['Flag_size']) #np.where(dataframe['Pixel'].astype('str').isin(pixel_list),1, dataframe['Flag_size'])
            # dataframe['0' in dataframe['Flag_size'].astype('str')].Sample.unique()

    dataframe['Flag']=np.where(dataframe[[column for column in dataframe.columns if 'Flag_' in column]].sum(axis=1)==0,0,1) # Default is 0, aka no anomaly
    df_summary = pd.merge(df_summary, dataframe[['Project_ID','Instrument','Sample','Sample_localpath','Sample_URL','datetime','Longitude','Latitude','Depth_min','Depth_max']+[column for column in dataframe.columns if 'Flag' in column]].drop_duplicates(subset=['Project_ID','Instrument','Sample','Sample_localpath','Sample_URL']).rename(columns={ 'datetime': 'Sample_datetime'}).astype({'Project_ID': str, 'Instrument': str, 'Sample': str, 'Sample_localpath': str,'Flag': str, 'Longitude': str, 'Latitude': str,'Sample_datetime': str}).groupby(['Project_ID', 'Instrument', 'Sample', 'Sample_localpath','Sample_URL', 'Flag', 'Longitude', 'Latitude', 'Sample_datetime']+[column for column in dataframe.columns if 'Flag_' in column]).apply(lambda x: pd.Series({'Sampling_depth_min_range': x.Depth_min.min(), 'Sampling_depth_max_range': x.Depth_max.max()})).reset_index(), how='right',on='Sample')
    df_summary=df_summary if 'Artefacts_count' in df_summary.columns else  df_summary.assign(Artefacts_count= 0,Artefacts_percentage=0, Validation_percentage= 0)
    return dataframe,df_summary[['Project_ID','Instrument','Sample','Sample_localpath','Sample_URL','Sample_datetime','Sampling_depth_min_range','Sampling_depth_max_range','Longitude','Latitude', 'Study_area', 'Longhurst_province','ROI_count', 'Count_uncertainty', 'Count_lower', 'Count_upper','Artefacts_count', 'Artefacts_percentage', 'Validation_percentage','Missing_field','Flag_missing', 'Flag_GPScoordinatesonland','Flag_dubiousGPScoordinates', 'Flag_count', 'Flag_artefacts','Flag_validation','Flag_size','Flag']].drop_duplicates(subset=['Project_ID','Instrument','Sample','Sample_localpath','Sample_URL'])

def save_flag_summary(df_flagged,df_standardizer,project_id,path_to_summary):
    gdf = gpd.GeoDataFrame(df_flagged[['Sample', 'Longitude', 'Latitude']].drop_duplicates().dropna(),geometry=gpd.points_from_xy(df_flagged[['Sample', 'Longitude', 'Latitude']].drop_duplicates().dropna().Longitude,df_flagged[['Sample', 'Longitude', 'Latitude']].drop_duplicates().dropna().Latitude))
    df_flagged['Study_area'] = pd.merge(df_flagged, gpd.tools.sjoin(gdf, oceans, predicate="within", how='left')[['Sample', 'name']], how='left', on='Sample')['name'].astype(str)
    df_flagged['Longhurst_province'] = pd.merge(df_flagged, gpd.tools.sjoin(gdf, Longhurst, predicate="within", how='left')[['Sample', 'ProvDescr']], how='left', on='Sample')['ProvDescr'].astype(str)
    instrument = df_standardizer['Instrument'][project_id]
    df_flagged= df_flagged.assign(Project_ID=project_id, Instrument=instrument)
    df_summary = df_flagged.rename(
        columns={'Flag': 'Sample_flag', 'Sample_URL': 'Sample_url', 'datetime': 'Sample_datetime'}).astype(
        {'Project_ID': str, 'Instrument': str, 'Sample': str, 'Sample_url': str, 'Sample_localpath': str,
         'Sample_flag': str, 'Longitude': str, 'Latitude': str, 'Study_area': str, 'Longhurst_province': str,
         'Sample_datetime': str}).groupby(
        ['Project_ID', 'Instrument', 'Sample', 'Sample_url', 'Sample_localpath', 'Sample_flag', 'Longitude', 'Latitude',
         'Study_area', 'Longhurst_province', 'Sample_datetime']).apply(lambda x: pd.Series(
        {'Sampling_depth_min_range': x.Depth_min.min(), 'Sampling_depth_max_range': x.Depth_max.max(),
         'object_number': len(x) if 'ROI_number' not in df_flagged.columns else x.ROI_number.sum()})).reset_index()
    path_to_summary_file = str(path_to_summary / instrument / 'Summary_project_{}.csv'.format(project_id)) if path_to_summary.stem.lower() == 'ecotaxa' else str(path_to_summary / 'Summary_project_{}.csv'.format(project_id))
    Path(path_to_summary_file).parent.mkdir(parents=True, exist_ok=True)
    df_summary.to_csv(path_to_summary_file, sep=',', index=False)


#quality_control_func(standardizer_path='~/GIT/PSSdb/raw/project_Zooscan_standardizer.xlsx',project_id=4023,report_path='~/GIT/PSSdb/reports')
def quality_control_func(standardizer_path,project_id,report_path,validation_threshold=0.95,project_list=[]):
    """
    Objective: This function adds flags to image samples based on multiple criteria and automatically generate a report for a given project
        :param standardizer_path: Full path of the standardizer spreadsheet containing project ID of interest
        :param project_id: ID of the project to be flagged in the standardizer spreadsheet
        :param report_path: Full path for the storage of the project control quality check interactive report
        :param validation_threshold: Threshold for sample/profile taxonomic annotation validation in percentage [0-1]. UVP and Zooscan projects only. Default is 95%
        :param project_list: List of project path used to process multiple files using parallel computation. Do not change
    :return: interactive report showing project's flags and update standardizer spreadsheet with flagged samples ID
    """
    sheets = [sheet for sheet in pd.ExcelFile(path_to_project_list).sheet_names if 'metadata' not in sheet]
    df_projects = pd.concat(map(lambda sheet: pd.read_excel(path_to_project_list, sheet_name=sheet), sheets))
    df_projects['Project_portal'] = df_projects.Project_localpath.apply(lambda path: Path(path).expanduser().name)
    df_projects = df_projects[df_projects['PSSdb_access'] == True]
    sheet = 0 if df_projects[df_projects.Project_ID.astype(str).isin([str(project_id)])].Project_portal.unique()[0] != 'ecopart' else df_projects[df_projects.Project_ID.astype(str).isin([str(project_id)])].Project_portal.unique()[ 0]  # Default is first sheet unless project is from Ecopart
    sheetname = 'Data' if (sheet == 0) and len(re.findall(r'UVP', df_projects[ df_projects.Project_ID.astype(str).isin([str(project_id)])].Instrument.unique()[0])) == 0 else 'ecotaxa' if (len(re.findall( r'UVP', df_projects[df_projects.Project_ID.astype(str).isin([str(project_id)])].Instrument.unique()[ 0])) > 0) and (sheet == 0) else sheet

    df_standardizer = pd.read_excel(standardizer_path, sheet_name=sheet, index_col=0)
    df_standardizer['Flag_path'][df_standardizer['Flag_path'].isna()] = ''
    flag_overrule_path = df_standardizer['Flag_path'][project_id]
    df_standardizer_metadata = pd.read_excel(standardizer_path, sheet_name='Metadata')

    # Read project datafile, subsetting variables of interest
    if len(project_list)==0:
        project_path_list = list(Path(df_standardizer['Project_localpath'][project_id].replace(cfg['raw_dir'],cfg['raw_dir']+cfg['standardized_raw_subdir']+'/')).expanduser().rglob("standardized_project_"+str(project_id)+'_*')) #list(Path(df_standardizer['Project_localpath'][project_id]).expanduser().rglob("*_"+str(project_id)+'_*')) if Path(df_standardizer['Project_localpath'][project_id]).expanduser().stem!=cfg['UVP_consolidation_subdir'] else list(Path(df_standardizer['Project_localpath'][project_id]).expanduser().glob("ecotaxa_export_"+str(project_id)+'_*'))
        project_path = natsorted([path for path in project_path_list if ('flag' not in str(path)) & ("._" not in str(path))])
    else:
        project_path=natsorted([path for path in project_list if ('flag' not in str(path)) & ("._" not in str(path))])
    # Check if samples have been flagged already and generate flags for new samples only.
    if (Path(flag_overrule_path).expanduser().is_file()) and len(project_list) == 0:
        df_flags = pd.read_csv(flag_overrule_path, sep=',')
        project_path=[path for path in project_path if str(path).replace(str(Path.home()),'~') not in df_flags.Sample_localpath.unique()]
        if len(project_path):
            print('\nExisting flags found at {}. Updating the file with new samples'.format(flag_overrule_path))
        else:
            print('\nExisting flags found at {}, but no additional samples.\nSaving flags summary to {} (if missing) and skipping project'.format(flag_overrule_path, path_to_project_list))
            sheets_project_list = [sheet for sheet in pd.ExcelFile(path_to_project_list).sheet_names if 'metadata' not in sheet]
            columns_project_list = dict( map(lambda sheet: (sheet, pd.read_excel(path_to_project_list, sheet_name=sheet).columns.tolist()), sheets_project_list))
            df_projects_metadata = pd.read_excel(path_to_project_list, sheet_name='metadata')
            df_projects = pd.concat( map(lambda sheet: pd.read_excel(path_to_project_list, sheet_name=sheet).assign(Portal=sheet), sheets_project_list)).reset_index(drop=True)
            df_project = df_projects[(df_projects.Project_ID.astype(str) == str(project_id)) & ( df_projects.Portal.astype(str).isin(['ecotaxa', 'ifcb']))]
            idx_project = df_project.index
            if (len(df_project) == 1):
                df_flags_summary = df_flags.groupby(['Project_ID']).apply(lambda x: pd.Series({'Number_samples': int(len(df_flags)), 'Number_flagged_samples': int(len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)'))),'Percentage_flagged_missing': len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[ x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_missing'] == 1]) / len(df_flags),'Percentage_flagged_GPScoordinatesonland': len( x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[ x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_GPScoordinatesonland'] == 1]) / len(df_flags),'Percentage_flagged_dubiousGPScoordinates': len( x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[ x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_dubiousGPScoordinates'] == 1]) / len(df_flags),'Percentage_flagged_count': len( x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[ x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_count'] == 1]) / len(df_flags), 'Percentage_flagged_artefact': len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[ x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_artefacts'] == 1]) / len(df_flags),'Percentage_flagged_validation': len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_validation'] == 1]) / len(df_flags) if 'flag_validation' in x.columns else 0,'Percentage_flagged_size': len( x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[ x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[ 'Flag_size'] == 1]) / len(df_flags)})).reset_index()
                df_project = pd.merge(df_project[['Project_ID'] + [column for column in df_project.columns if column not in df_flags_summary.columns]], df_flags_summary, how='left', on='Project_ID').set_index(idx_project)
                df_projects = pd.concat([df_projects[df_projects.Portal == df_project['Portal'].values[0]].drop(index=df_project.index),df_project], axis=0).sort_values(['PSSdb_access', 'Instrument', 'Project_ID'],ascending=[False, False, True])
                print('\nAdding flag summary to project list')
                with pd.ExcelWriter(str(path_to_project_list), engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
                    df_projects[[column for column in columns_project_list.get(df_project['Portal'].values[0]) if column not in df_flags_summary.columns.tolist()[1:]] + df_flags_summary.columns.tolist()[1:]].to_excel(writer,sheet_name=df_project['Portal'].values[0],index=False)
                    if any(df_flags_summary.columns.isin(df_projects_metadata.Variables) == False):
                        pd.concat([df_projects_metadata, pd.DataFrame({'Variables': df_flags_summary.columns[1:], 'Variable_types': df_flags_summary.dtypes[1:],'Units/Values': ['#', '#'] + ['[0-1]'] * 7,'Description': ['Number of samples/nets/profiles in project','Number of samples flagged during the project control quality check','Percentage of samples flagged due to missing data/metadata','Percentage of samples flagged due to GPS coordinates on land (according to a 10m resolution)','Percentage of samples flagged due to dubious GPS coordinates (0x0 degrees)','Percentage of samples flagged due to low ROI count per sample', 'Percentage of samples flagged due to high percentage of artefacts', 'Percentage of samples flagged due to low validation of the taxonomic annotations (* UVP and Zooscan only)','Percentage of samples flagged due to multiple pixel size conversion factor']})],axis=0).to_excel(writer, sheet_name='metadata', index=False)
            return
    else:
        df_flags = pd.DataFrame()
        df_flagged_existing = pd.DataFrame()

    if len(project_path)>0:

         # Set missing values to NA
         na = str(df_standardizer.loc[project_id]['NA_value']).split(';')
         columns_dict=dict(zip(list(df_standardizer.loc[ project_id, [column for column in df_standardizer.columns if 'field' in column]].dropna().values),list(df_standardizer.loc[ project_id, [column for column in df_standardizer.columns if 'field' in column]].dropna().index.str.replace('_field',''))))
         #columns_dict['object_number']='object_number' # Adding this column for UVP consolidated projects

         dtypes_dict_all = dict(zip(['Cruise', 'Station', 'Profile', 'Sample', 'Latitude', 'Longitude', 'Depth_min', 'Depth_max','Sampling_date', 'Sampling_time', 'Volume_analyzed', 'Volume_imaged', 'ROI', 'Area', 'Pixel', 'Minor_axis','Major_axis', 'Biovolume','ESD', 'Category', 'Annotation', 'Sampling_type', 'Sampling_lower_size', 'Sampling_upper_size','ROI_number', 'Sampling_description'],[str, str, str, str, float, float, float, float, str, str, float,float, str, float,float, float, float, float, float, str, str, str, float, float, int, str]))
         #df_samples=list(map(lambda path:{str(path).replace(str(Path.home()),'~'):pd.read_table(path,sep=',',usecols=['Sample']).Sample.unique()},project_path))
         flagged_df = pd.DataFrame()
         overruled_df = pd.DataFrame()
         if (df_standardizer.loc[project_id]['Instrument']=='IFCB') & (len(project_path)>1):
             # Read all datafiles
             #chunk = 1000
             with tqdm(desc='', total=len(project_path), bar_format='{desc}{bar}', position=0, leave=True) as bar,ThreadPool() as pool:
                 for path in project_path:
             #with ThreadPool() as pool:
                 #for result in pool.map(lambda path: control_quality_func(standardizer_path,project_id,report_path,validation_threshold=0.95,project_list=[path]), project_path, chunksize=chunk):  # pool.map(lambda path: (columns := pd.read_table(path,sep=",", nrows=0).columns, pd.read_table(path,sep=",",usecols=[column for column in ['Project_ID','Instrument','Longitude','Latitude','Sample','Sampling_date','Sampling_time','Depth_min','Depth_max','Volume_imaged','ROI','ROI_number','Area','Category'] if column in columns],dtype=dtypes_dict_all))[-1],standardized_files, chunksize=chunk):
                     percent = np.round(100 * ((1+bar.n )/ len(project_path)), 1)
                     bar.set_description('Flagging datafile {} (%s%%)'.format(str(1+bar.n) + "/" + str(len(project_path))) % percent, refresh=True)
                     result=quality_control_func(standardizer_path,project_id,report_path,validation_threshold=0.95,project_list=[path])
                     overruled_df = pd.concat([overruled_df, result], axis=0)
                     ok = bar.update(n=1)
             overruled_df=overruled_df.reset_index(drop=True)
             flagged_df=overruled_df

         else:
             df = pd.concat(map(lambda path:  pd.read_table(path,sep=',',usecols=list(set(list(df_standardizer.loc[ project_id, [column for column in df_standardizer.columns if ('field' in column) & (column!='Dilution_field')]].dropna().index.str.replace('_field',''))+['Sampling_lower_size','Sampling_upper_size','Sampling_date','Sampling_time','Depth_max','Depth_min','ROI_number'])),dtype=dtypes_dict_all,parse_dates={'datetime':['Sampling_date','Sampling_time']}).assign(Sample_localpath=str(path).replace(str(Path.home()),'~')), project_path)).reset_index(drop=True) #, na_values=na, keep_default_na=True#pd.concat(map(lambda path: (columns := pd.read_table(path, nrows=0).columns, pd.read_table(path, usecols=[header for  header in ['object_number']+list(df_standardizer.loc[project_id, [  column  for column in df_standardizer.columns if 'field' in column]].dropna().values) if columns.isin([header]).any()]).assign(Sample_localpath=str(path).replace(str(Path.home()),'~')).rename(columns=columns_dict))[ -1], project_path)).reset_index(drop=True) #, na_values=na, keep_default_na=True
             if 'ROI_number' not in df.columns:
                 df['ROI_number']=1
         if len(flagged_df) == 0:
             df['datetime']=df.datetime.dt.strftime('%Y-%m-%dT%H%M%SZ').astype(str)
             # Replace NA ROI/mino/annotation by '' to avoid flagging project on missing ROI id (e.g small particles in UVP consolidated files )/ annotations (e.g. mixture of predicted,unclassified project)
             df.ROI = np.where(df.ROI.isna(), '', df.ROI)
             if 'Minor_axis' in df.columns:
                 df.Minor_axis = np.where(df.Minor_axis.isna(), '', df.Minor_axis)
             if 'Major_axis' in df.columns:
                 df.Major_axis = np.where(df.Major_axis.isna(), '', df.Major_axis)
             if 'Sampling_type'  in df.columns:
                 df.Sampling_type=np.where(df.Sampling_type.isna(), '', df.Sampling_type)

             if 'Category' not in df.columns:
                 df['Category'] = ''
                 df['Annotation'] = 'unclassified'
             if ('Category' in df.columns) and ('Annotation' not in df.columns):
                 df['Annotation'] = 'predicted'
             df.Category = np.where((df.Category.isna()) | (df.Category.astype(str)=='nan'),'',df.Category) # Avoid flagging of small particles without annotations for UVP projects
             df.Annotation = np.where((df.Annotation.isna()) | (df.Annotation.astype(str) == 'nan'), 'unclassified',df.Annotation)
             df = df.astype({key:value for key,value in dict(zip(['Sample', 'Latitude', 'Longitude', 'Volume_analyzed', 'Pixel'], [str, float, float, float, float])).items() if key in df.columns})
             if 'Cruise' not in df.columns:
                 df['Cruise']='nan'

             # Adding URL to check flagged samples
             # Query EcoTaxa API to retrieve samples ID/URL
             if (df_standardizer['Project_source'][project_id] == 'https://ecotaxa.obs-vlfr.fr/prj/' + str(project_id)):
                 with ecotaxa_py_client.ApiClient(configuration) as api_client:
                     api_instance = samples_api.SamplesApi(api_client)
                     samples = api_instance.samples_search(project_ids=int(project_id), id_pattern='')  #
                 df = pd.merge(df, pd.DataFrame({'Sample': pd.Series(map(lambda x: x['orig_id'], samples)).astype(str(df.dtypes['Sample'])), 'Sample_ID': list(map(lambda x: x['sampleid'], samples))}), how='left', on='Sample')
             df = pd.merge(df, df.groupby(['Sample']).apply(lambda x: pd.Series({'Sample_URL': r'{}?taxo=&taxochild=&ipp=100&zoom=100&sortby=&magenabled=0&popupenabled=0&statusfilter=&samples={}&sortorder=asc&dispfield=&projid={}&pageoffset=0"'.format(df_standardizer['Project_source'][project_id],str(x.Sample_ID.unique()[0]),project_id) if df_standardizer['Project_source'][project_id] == 'https://ecotaxa.obs-vlfr.fr/prj/' + str(project_id) else r'{}&bin={}'.format(df_standardizer['Project_source'][ project_id],str(x.Sample.unique()[0])) if 'ifcb' in df_standardizer['Project_source'][project_id] else ''})),how='left', on='Sample')
             if 'Sample_ID' in df.columns:
                 df=df.drop(columns=['Sample_ID'])


         if len(flagged_df) == 0:
             # Append flags, Avoid flagging on low validation for IFCB projects by setting validation threshold to 0
             validation_threshold = validation_threshold if df_standardizer.loc[project_id]['Instrument'] != 'IFCB' else 0
             flagged_df,overruled_df=flag_func(df.assign(Project_ID=project_id,Instrument=df_standardizer['Instrument'][project_id]),validation_threshold=validation_threshold)
             # Override flag_missing for optional variables: Sampling_size, Dilution, Annotation, Category
             flagged_df = pd.merge(flagged_df.drop(columns=['Missing_field']),flagged_df.groupby(['Sample']).apply(lambda x:pd.Series({'Missing_field':'' if all(pd.Series((x['Missing_field'].astype(str).unique()[0]).split(';')).isin(['Sampling_upper_size','Sampling_lower_size','Category','Annotation','Dilution'])) else x['Missing_field'].astype(str).unique()[0]})).reset_index(),how='left',on='Sample')
             #flagged_df['Missing_field'] = flagged_df['Missing_field'].apply(lambda x: '' if pd.Series((field in ['Sampling_upper_size', 'Sampling_lower_size', 'Category', 'Annotation', 'Dilution'] for field in str(x).split(';'))).all() else x)
             flagged_df.loc[flagged_df['Missing_field'] == '', 'Flag_missing'] = 0
             flagged_df['Flag'] = np.where(flagged_df[[column for column in flagged_df.columns if 'Flag_' in column]].sum(axis=1)==0,0,1)
             overruled_df=pd.merge(overruled_df.drop(columns=['Missing_field','Flag']),flagged_df[['Sample','Missing_field','Flag']].drop_duplicates(),how='left',on='Sample')[['Project_ID','Instrument','Sample','Sample_localpath','Sample_URL','Sample_datetime','Sampling_depth_min_range','Sampling_depth_max_range','Longitude','Latitude', 'Study_area', 'Longhurst_province','ROI_count', 'Count_uncertainty', 'Count_lower', 'Count_upper','Artefacts_count', 'Artefacts_percentage', 'Validation_percentage','Missing_field','Flag_missing', 'Flag_GPScoordinatesonland','Flag_dubiousGPScoordinates', 'Flag_count', 'Flag_artefacts','Flag_validation','Flag_size','Flag']]




         if len(project_list) != 0:
             return overruled_df#flagged_df

         # Generating flags overruling datafile that will be read to filter samples out during standardization
         path_to_datafile = Path(cfg['raw_dir']).expanduser() / cfg['flag_subdir'] / Path(df_standardizer['Project_localpath'][project_id]).stem
         path_to_datafile.mkdir(parents=True, exist_ok=True)
         report_path = Path(report_path).expanduser() / Path(df_standardizer['Project_localpath'][project_id]).stem
         report_path.mkdir(parents=True, exist_ok=True)
         report_filename = 'Report_project_' + str(project_id) + '.html'
         path_to_report = Path(report_path).expanduser() / report_filename

         # Create readme file for the report and flag folders:
         if not Path(path_to_datafile.parent / 'README.txt').is_file():
             with open(str(Path(path_to_datafile.parent / 'README.txt')), 'w') as file:
                 file.write("README file for project quality control (First created on February 10, 2023):\n\nThis directory contains a summary table of project control quality flags\nEach table includes individual sample info, such as:\n-Project_ID (native format): Unique ID of the project/cruise/program\n-Instrument (Zooscan, IFCB, UVP, Epson scanner): Project imaging device\n-Sample (native format): Unique sample ID\n-Sample_URL: Link to the images/project files online repository (e.g. Ecotaxa, IFCB dashboard)\n-Sample_localpath: Local path of project datafile storage (identical for Ecotaxa Zooscan/IFCB projects, different for UVP consolidated and IFCB dashboard projects)\n-Sampling datetime (yyyy-mm-ddThh:mm:ssZ, UTC): Sampling time\n-Longitude (decimal degree): Longitudinal coordinate [-180:180]\n-Latitude (decimal degree): Latitudinal coordinate [-90:90]\n-Study_area: Study area according to the Global Oceans and Seas of the standard georeferenced marine regions website (https://www.marineregions.org/)\n-Longhurst_province: Longhurst province of the sample location according to the standard georeferenced marine regions website (https://www.marineregions.org/)\n-Sampling_min_depth_range (meters): Minimum/Shallowest sampling depth\n-Sampling_max_depth_range (meters): Maximum/Deepest sampling depth\n-ROI_count (#): Number of Region Of Interest imaged in the sample.\n-Count_uncertainty (probability [0-1]): Uncertainty of the count according to the Poisson distribution.\n-Artefacts_count (#): Number of artefacts imaged in the sample.\n-Artafacts_percentage (% [0-1]): Percentage of artefacts in the sample.\n-Validation_percentage (% [0-1]): Percentage of 'validated' annotations in the sample.\n\n\nCurrent flagging is done based on 7 critera:\n\n-GPS coordinates on land (Falg_GPScoordinatesonland): takes 1 if GPS coordinates correspond to a location on land (according to a 10m spatial resolution)\n\n-Dubious GPS coordinates (Flag_dubiousGPScoordinates): takes 1 if location is at 0 degrees latitude and longitude\n\n-Missing variables (Flag_missing): takes 1 if a variable required for project standardization is missing (check required variables on the project standardizer spreadsheets)\n\n-Low particles count per sample (Flag_count): takes 1 if the total number of particles detected in a given sample yields high uncertainty (>5%) assuming counts follow a Poisson distribution\n\n-High percentage of artefacts (Flag_artefacts): takes 1 if the percentage of artefacts in a given sample is superior to 20%\n\n-Low percentage of taxonomic annotations validation (Flag_validation, *UVP and Zooscan only): takes 1 if the percentage of validation in a given sample/profile is less than 95%\n\n-Multiple pixel-to-size calibration factors (Flag_size): takes 1 if the project include multiple size calibration factors per cruise\n\nEach flag is assigned a boolean factor which takes 0 if the given sample has not been flagged or 1 if it has been flagged.\nThe overall flag (Flag) is predicted based on the result of 1 or more flagged critera.\nThe overrule boolean factor allows to overrule this flag (reset False to True) if the data owner deemed the sample should be kept for further processing.\nFlagged samples will be discarded during standardization.\nThe flag table will be updated if a project have been updated with additional samples.\nInteractive reports of project control quality check can be found at: {}\n\nContact us: nmfs.pssdb@noaa.gov".format(report_path.parent))
         if not Path(report_path.parent / 'README.txt').is_file():
             with open(str(Path(report_path.parent / 'README.txt')), 'w') as file:
                 file.write("README file for project quality control report (First created on February 10, 2023):\n\nThis directory contains interactive reports of project control quality flags.\nThe bottom table include flagged samples that can be checked by clicking on the URL of the first column.\nCurrent flagging is done based on 7 critera:\n\n-GPS coordinates on land (Falg_GPScoordinatesonland): takes 1 if GPS coordinates correspond to a location on land (according to a 10m spatial resolution)\n\n-Dubious GPS coordinates (Flag_dubiousGPScoordinates): takes 1 if location is at 0 degrees latitude and longitude\n\n-Missing variables (Flag_missing): takes 1 if a variable required for project standardization is missing (check required variables on the project standardizer spreadsheets)\n\n-Low particles count per sample (Flag_count): takes 1 if the total number of particles detected in a given sample yields high uncertainty (>5%) assuming counts follow a Poisson distribution\n\n-High percentage of artefacts (Flag_artefacts): takes 1 if the percentage of artefacts in a given sample is superior to 20%\n\n-Low percentage of taxonomic annotations validation (Flag_validation, *UVP and Zooscan only): takes 1 if the percentage of validation in a given sample/profile is less than 95%\n\n-Multiple pixel-to-size calibration factors (Flag_size): takes 1 if the project include multiple size calibration factors per cruise\n\nEach flag is assigned a boolean factor which takes 0 if the given sample has not been flagged or 1 if it has been flagged.\nThe overall flag (Flag) is predicted based on the result of 1 or more flagged critera.\nFlagged samples will be discarded during standardization.\nThe flag report will be updated if a project have been updated with additional samples.\nTables of project control quality check can be found at: {}\n\nContact us: nmfs.pssdb@noaa.gov".format(path_to_datafile.parent))

         #if len(flag_overrule_path) == 0 or Path(flag_overrule_path).expanduser().is_file() == False:
         overrule_name ='project_{}_flags.csv'.format(str(project_id))
         print('\nSaving flags to', str(path_to_datafile / overrule_name),'\nSet Overrule to True if you wish to keep samples for further processing', sep=" ")
         #overruled_df = flagged_df[['Sample', 'Sample_URL','Sample_localpath', 'Flag'] + [column for column in flagged_df.columns if('Flag_' in column) or ('Missing_field' in column)]].drop_duplicates()  # flagged_df[flagged_df['Flag']==0][['Sample','Flag']].drop_duplicates()
         overruled_df['Overrule'] = False
         overruled_df = pd.concat([df_flags.reset_index(drop=True), overruled_df.reset_index(drop=True)], axis=0).sort_values( by=['Flag', 'Sample'], ascending=[False, True])
         overruled_df.to_csv(str(path_to_datafile / overrule_name), sep=',',index=False)

         # Update standardizer spreadsheet with flagged samples path and save
         df_standardizer['Flag_path'] = df_standardizer['Flag_path'].astype(str)
         df_standardizer['Flag_path'][project_id] = str(path_to_datafile / overrule_name).replace(str(Path.home()), '~')
         print('\nUpdating standardizer spreadsheet with path of flagged samples/profiles ID datafile')
         df_standardizer['Project_ID'] = df_standardizer.index
         df_standardizer = df_standardizer[['Project_ID'] + df_standardizer.columns[:-1].tolist()]
         with pd.ExcelWriter(str(standardizer_path), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                df_standardizer.to_excel(writer, sheet_name=sheetname, index=False)


         # Add flag summary to project list and save
         sheets_project_list = [sheet for sheet in pd.ExcelFile(path_to_project_list).sheet_names if 'metadata' not in sheet]
         columns_project_list = dict( map(lambda sheet: (sheet, pd.read_excel(path_to_project_list, sheet_name=sheet).columns.tolist()), sheets_project_list))
         df_projects_metadata = pd.read_excel(path_to_project_list, sheet_name='metadata')
         df_projects = pd.concat( map(lambda sheet: pd.read_excel(path_to_project_list, sheet_name=sheet).assign(Portal=sheet), sheets_project_list)).reset_index(drop=True)
         df_project = df_projects[(df_projects.Project_ID.astype(str) == str(project_id)) & (df_projects.Portal.astype(str).isin(['ecotaxa', 'ifcb']))]
         idx_project = df_project.index
         if (len(df_project) == 1):
             df_flags_summary = overruled_df.groupby(['Project_ID']).apply(lambda x: pd.Series({'Number_samples': int(len(overruled_df)), 'Number_flagged_samples': int(len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)'))), 'Percentage_flagged_missing': len(x[x['Flag_missing'] == 1]) / len(overruled_df), 'Percentage_flagged_GPScoordinatesonland': len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_GPScoordinatesonland'] == 1]) / len( overruled_df),'Percentage_flagged_dubiousGPScoordinates': len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_dubiousGPScoordinates'] == 1]) / len(overruled_df), 'Percentage_flagged_count': len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_count'] == 1]) / len(overruled_df),'Percentage_flagged_artefact': len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_artefacts'] == 1]) / len(overruled_df),'Percentage_flagged_validation': len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_validation'] == 1]) / len(overruled_df), 'Percentage_flagged_size': len(x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')[x.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Flag_size'] == 1]) / len(overruled_df)})).reset_index()
             df_project = pd.merge(df_project[['Project_ID']+[column for column in df_project.columns if column not in df_flags_summary.columns]], df_flags_summary, how='left', on='Project_ID').set_index(idx_project)
             df_projects = pd.concat([df_projects[df_projects.Portal == df_project['Portal'].values[0]].drop(index=df_project.index),df_project], axis=0).sort_values(['PSSdb_access', 'Instrument', 'Project_ID'], ascending=[False, False, True])
             print('\nAdding flag summary to project list')
             with pd.ExcelWriter(str(path_to_project_list), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                 df_projects[[column for column in columns_project_list.get(df_project['Portal'].values[0]) if column not in df_flags_summary.columns.tolist()[1:]] + df_flags_summary.columns.tolist()[1:]].to_excel(writer, sheet_name=df_project['Portal'].values[0], index=False)
                 if any(df_flags_summary.columns.isin(df_projects_metadata.Variables) == False):
                     pd.concat([df_projects_metadata, pd.DataFrame({'Variables': df_flags_summary.columns[1:], 'Variable_types': df_flags_summary.dtypes[1:],'Units/Values': ['#', '#'] + ['[0-1]'] * 7,'Description': ['Number of samples/nets/profiles in project','Number of samples flagged during the project control quality check','Percentage of samples flagged due to missing data/metadata','Percentage of samples flagged due to GPS coordinates on land (according to a 10m resolution)','Percentage of samples flagged due to dubious GPS coordinates (0x0 degrees)','Percentage of samples flagged due to low ROI count per sample','Percentage of samples flagged due to high percentage of artefacts', 'Percentage of samples flagged due to low validation of the taxonomic annotations (* UVP and Zooscan only)','Percentage of samples flagged due to multiple pixel size conversion factor']})],axis=0).to_excel(writer, sheet_name='metadata', index=False)

         # Generating interactive project report
         summary_df = overruled_df.sort_values(['Sample_datetime'])#flagged_df_subset.groupby(['Sample','Sample_URL',  'Latitude', 'Longitude']+[column for column in flagged_df.columns if ('Flag' in column) or ('Missing_field' in column)], dropna=False).apply(lambda x: pd.Series({'ROI_count': x.ROI_number.sum(),'Count_error':np.diff(poisson.ppf([0.05/2,1-(0.05/2)], mu=sum(x.ROI_number)))[0] if len(x.ROI) else 0,'Validation_percentage':len(x[x['Annotation'].isin(['validated'])].ROI) / len(x.ROI) if len(x.ROI) else 0,'Artefacts_percentage':len(x[x.Category.astype(str).str.lower().apply(lambda annotation:len(re.findall(r'bead|bubble|artefact|artifact',annotation))>0)].ROI) / len(x.ROI) if len(x.ROI) else 0})).reset_index()
         summary_df['Sample_URL'] =  summary_df[['Sample', 'Sample_URL']].apply(lambda x: pd.Series({'Sample_URL': r'<a href="{}">{}</a>'.format( x.Sample_URL,x.Sample)}),axis=1)
         summary_df['Count_error']=summary_df.ROI_count.apply(lambda count:np.diff(poisson.ppf([0.05/2,1-(0.05/2)], mu=count))[0])
         summary_df=summary_df.astype(dict(zip([column for column in summary_df.columns if 'Flag' in column],[int]*len([column for column in summary_df.columns if 'Flag' in column]))))
         subset_summary_df=summary_df.dropna(subset=[ 'Sample','Sample_URL'])

         if len(subset_summary_df):
             fig = make_subplots(subplot_titles=['', '', 'Flag: 0 (no flag), 1 (flagged)', 'Percentage flagged samples/profiles: {}%'.format(str(np.round( 100 * len(subset_summary_df[subset_summary_df.Flag == 1]) / len(subset_summary_df), 2)))], rows=3, cols=2,
                                 specs=[[{"type": "scattergeo", "rowspan": 2}, {"type": "scatter", 't': 0.05}],
                                        [None, {"type": "scatter", 'b': 0.05}],
                                        [{"type": 'table', "colspan": 2}, None]], row_heights=[0.3, 0.3, 0.4],vertical_spacing=0.02)

             fig.layout.annotations[0].update(xanchor='left', yanchor='top', x=0.0, y=1.0, font={'size': 13})
             fig.layout.annotations[1].update(xanchor='left', yanchor='top', x=0.0, y=0.98, font={'size': 13})

             # Map, left column
             subset_summary_df['colors']=np.where(subset_summary_df.Flag_GPScoordinatesonland==1, 'red','black')
             subset_summary_df.Sample = pd.Categorical(subset_summary_df.Sample, categories=subset_summary_df[['Sample','Sample_datetime']].drop_duplicates().sort_values(['Sample_datetime']).Sample.values,ordered=True)
             #gdf = gpd.GeoDataFrame(subset_summary_df,geometry=gpd.points_from_xy(subset_summary_df.Longitude,subset_summary_df.Latitude)).set_index('Sample').set_crs(epsg=4326, inplace=True)
             if len(subset_summary_df[subset_summary_df.Flag_GPScoordinatesonland==0])>0:
                 sub_data = subset_summary_df[subset_summary_df.Flag_GPScoordinatesonland == 0]
                 data_geo =dict(type='scattergeo',
                             name='No flag',
                             lon=sub_data.Longitude.astype(float),
                             lat=sub_data.Latitude.astype(float),
                             hovertext='Sample ID: ' + sub_data.Sample.astype(str) + '<br> Longitude (ºE): ' + sub_data.Longitude.astype(str) + '<br> Latitude (ºN): ' + sub_data.Latitude.astype(str),
                             hoverinfo="text", marker=dict(color='black', line=dict(color=sub_data.colors, width=2)), geojson="natural earth", showlegend=True,geo='geo')
                 layout = dict()
                 layout['geo2'] = dict(
                 projection_rotation=dict(lon=np.nanmean(sub_data.Longitude.astype(float)), lat=np.nanmean(sub_data.Latitude.astype(float)), roll=0),
                 center=dict(lon=np.nanmean(sub_data.Longitude.astype(float)), lat=np.nanmean(sub_data.Latitude.astype(float))),
                 projection_type='orthographic', lataxis_range=[-90, 90], lonaxis_range=[-180, 180],
                 domain=dict(x=[0.0, 0.67], y=[0.4, 0.95]), bgcolor='rgba(0,0,0,0)')
                 layout['geo'] = dict(lonaxis_range=[-180, 180], lataxis_range=[-80, 80],
                                  domain=dict(x=[0.0, 0.2], y=[0.4, 0.52]))
                 # go.Figure(data=[data_geo, data_geo], layout=layout)
                 fig.add_trace(data_geo, row=1, col=1).update_layout( go.Figure(data=[data_geo, data_geo], layout=layout).layout, margin={"r": 0, "t": 0, "l": 0, "b": 0})
                 fig.layout['geo'] = layout['geo']
                 data_geo.update({'geo': 'geo2', 'showlegend': False})  # Update geo trace before adding the zoomed
                 fig.add_trace(data_geo, row=1, col=1)
                 fig.data[len(fig.data) - 1]['geo'] = 'geo2'  # Need to update geo trace manually
                 fig.layout['geo2'] = layout['geo2']
             if len(subset_summary_df[subset_summary_df.Flag_GPScoordinatesonland == 1]) > 0:
                 # Definition of the zoomed (geo2) and inset (geo) maps
                 sub_data=subset_summary_df[subset_summary_df.Flag_GPScoordinatesonland == 1]
                 data_geo = dict(type='scattergeo',
                          name='GPS coordinates on land<br>(Flag_GPScoordinatesonland)',
                          lon=sub_data.Longitude.astype(float),
                          lat=sub_data.Latitude.astype(float),
                          hovertext='Sample ID: ' + sub_data.Sample.astype(str)+ '<br> Longitude (ºE): '+sub_data.Longitude.astype(str)+'<br> Latitude (ºN): '+sub_data.Latitude.astype(str),
                          hoverinfo="text", marker=dict(size=3, line=dict(color=sub_data.colors, width=2), color='black'), geojson="natural earth", showlegend=True,
                          geo='geo')

                 layout = dict()
                 layout['geo2'] = dict(projection_rotation=dict(lon=np.nanmean(sub_data.Longitude.astype(float)), lat=np.nanmean(sub_data.Latitude.astype(float)), roll=0),center=dict(lon=np.nanmean(sub_data.Longitude.astype(float)), lat=np.nanmean(sub_data.Latitude.astype(float))),projection_type='orthographic',lataxis_range=[-90,90], lonaxis_range=[-180, 180],
                 domain=dict(x=[0.0, 0.67], y=[0.4, 0.95]),bgcolor='rgba(0,0,0,0)')
                 layout['geo'] = dict(lonaxis_range=[-180,180],lataxis_range=[-80,80],
                 domain=dict(x=[0.0, 0.2], y=[0.4, 0.52]))
                 #go.Figure(data=[data_geo, data_geo], layout=layout)
                 fig.add_trace(data_geo, row=1, col=1).update_layout(go.Figure(data=[data_geo, data_geo], layout=layout).layout, margin={"r": 0, "t": 0, "l": 0, "b": 0})
                 fig.layout['geo'] = layout['geo']
                 data_geo.update({'geo': 'geo2', 'showlegend': False})  # Update geo trace before adding the zoomed
                 fig.add_trace(data_geo, row=1, col=1)
                 fig.data[len(fig.data) - 1]['geo'] = 'geo2'  # Need to update geo trace manually
                 fig.layout['geo2'] = layout['geo2']


             subset_summary_df['colors'] = np.where(subset_summary_df.Flag_dubiousGPScoordinates == 1, 'orange', 'black')
             if len(subset_summary_df[subset_summary_df.Flag_dubiousGPScoordinates == 1]) > 0:
                 sub_data = subset_summary_df[subset_summary_df.Flag_dubiousGPScoordinates == 1]
                 data_geo = dict(type='scattergeo',
                             name='Dubious GPS coordinates<br>(Flag_dubiousGPScoordinates)',
                             lon=sub_data.Longitude.astype(float),
                             lat=sub_data.Latitude.astype(float),
                             hovertext='Sample ID: ' + sub_data.Sample.astype( str) + '<br> Longitude (ºE): ' + sub_data.Longitude.astype(str) + '<br> Latitude (ºN): ' + sub_data.Latitude.astype(str),
                             hoverinfo="text", marker=dict(color='black', line=dict(color=sub_data.colors, width=2)),
                             geojson="natural earth", showlegend=True, geo='geo')
                 layout['geo2'] = dict(projection_rotation=dict(lon=np.nanmean(sub_data.Longitude.astype(float)), lat=np.nanmean(sub_data.Latitude.astype(float)), roll=0),
                                center=dict(lon=np.nanmean(sub_data.Longitude.astype(float)), lat=np.nanmean(sub_data.Latitude.astype(float))),
                                projection_type='orthographic', lataxis_range=[-90, 90], lonaxis_range=[-180, 180],
                                domain=dict(x=[0.0, 0.67], y=[0.4, 0.95]), bgcolor='rgba(0,0,0,0)')
                 fig.add_trace(data_geo, row=1, col=1)
                 fig.data[len(fig.data) - 1]['geo'] = 'geo2'  # Need to update geo trace manually

             # Scatterplot 1, top-right panel: Count per sample
             subset_summary_df['colors'] = np.where(subset_summary_df.Flag_count==1, 'rgba(0,0,255,0.4)', 'black') # Low ROI count
             fig.add_trace(go.Scatter(x=subset_summary_df.Sample,
                                  y=subset_summary_df.ROI_count,
                                  error_y=dict(type="data",array=subset_summary_df.Count_error,width=0, thickness=0.1,color=subset_summary_df.colors.values[0]),
                                  hovertext="Sample ID: " + subset_summary_df.Sample.astype(str), hoverinfo="none",
                                  marker=dict(size=3.5, color='black'),
                                  mode='markers',showlegend=False, visible=True), row=1, col=2)
             if len(subset_summary_df[subset_summary_df.Flag_count == 1]) > 0:
                 fig.add_trace(go.Scatter(name='Low ROI counts<br>(Flag_count)',x=subset_summary_df[subset_summary_df.Flag_count==1].Sample,
                             y=subset_summary_df[subset_summary_df.Flag_count==1].ROI_count,
                             error_y=dict(type="data",array=subset_summary_df[subset_summary_df.Flag_count==1].Count_error,width=0,thickness=0.1,color=subset_summary_df[subset_summary_df.Flag_count==1].colors.values[0]),
                             hoverinfo="none",
                             marker=dict(size=3.5,color='black', line=dict(color=subset_summary_df[subset_summary_df.Flag_count==1].colors, width=2)), mode='markers',
                             showlegend=True, visible=True), row=1, col=2)
             # Scatterplot 2, middle-right panel: Percentage of artefacts
             subset_summary_df['colors'] = np.where(subset_summary_df.Flag_artefacts == 1, 'rgba(212,85,0,0.6)','black')  # High percentage of artefacts
             fig.add_trace(go.Scatter(x=subset_summary_df.Sample,
                                      y=subset_summary_df.Artefacts_percentage,
                                      hovertext="Sample ID: " + subset_summary_df.Sample.astype(str) + '<br>% of artefacts: ' + np.round(100 * subset_summary_df.Artefacts_percentage, 1).astype(str) + '%',
                                      hoverinfo="text",
                                      marker=dict(size=4.5, color='black'), mode='markers', showlegend=False, visible=True), row=2, col=2)

             if len(subset_summary_df[subset_summary_df.Flag_artefacts == 1]) > 0:
                 fig.add_trace(go.Scatter(name='High percentage of artefacts<br>(Flag_artefacts)',
                             x=subset_summary_df[subset_summary_df.Flag_artefacts == 1].Sample,
                             y=subset_summary_df[subset_summary_df.Flag_artefacts == 1].Artefacts_percentage,
                             hoverinfo="none",
                             marker=dict(size=3.5,color='black', line=dict(color=subset_summary_df[subset_summary_df.Flag_artefacts == 1].colors, width=2)),
                             mode='markers',
                             showlegend=True, visible=True), row=2, col=2)


             subset_summary_df['colors'] = np.where(subset_summary_df.Flag_validation == 1, 'rgba(95,211,188,0.6)', 'black')
             fig.add_trace(go.Scatter(x=subset_summary_df.Sample,
                                  name='Percentage of validation',
                                  y=subset_summary_df.Validation_percentage,
                                  hovertext="Sample ID: " + subset_summary_df.Sample.astype(str)+'<br>% of validation: '+np.round(100*subset_summary_df.Validation_percentage,1).astype(str)+'%',
                                  hoverinfo="text", marker=dict(color='black'), mode='lines',
                                  showlegend=True, visible=True), row=2, col=2)

             if len(subset_summary_df[subset_summary_df.Flag_validation == 1]) > 0:
                 fig.add_trace(go.Scatter(name='Low percentage of validation<br>(Flag_validation)',x=subset_summary_df[subset_summary_df.Flag_validation == 1].Sample,
                                          y=subset_summary_df[subset_summary_df.Flag_validation == 1].Validation_percentage,
                                          hoverinfo="none",marker=dict(size=4.5, color='black', line=dict(color=subset_summary_df[subset_summary_df.Flag_validation == 1].colors,width=2)), mode='markers', showlegend=True, visible=True), row=2, col=2)


             # Table
             summary_df['Flag_missing']=summary_df.apply( lambda x: str(x.Flag_missing) + ' (' + x.Missing_field + ')' if x.Flag_missing == 1 else x.Flag_missing,axis=1)
             fig.add_trace(go.Table(header=dict(values=['Sample/Profile ID<br>']+[column+'<br>' for column in summary_df.columns if 'Flag' in column],align=np.repeat('center',1+len([column for column in summary_df.columns if 'Flag' in column])),
                                       line_color='rgba(255,255,255,0)',fill_color='rgba(255,255,255,1)'),columnwidth =[4,5,4,4,2,2,2.5,1.5,1],
                           cells=dict(values=summary_df[summary_df.Flag==1][['Sample_URL']+[column for column in summary_df.columns if 'Flag' in column]].T)), row=3, col=1)
             # Update subplot domains for small table
             if len(summary_df[summary_df.Flag==0][['Sample_URL']])<8: # scrolling threshold
                 fig.layout['geo2']['domain']['y'] = [ 0.02 + 0.42 - (0.384 / 8) * max([1, 8 - len(summary_df[summary_df.Flag == 1][['Sample_URL']])]), 0.95]  # fig.update_geos(domain={'y':[0.02+0.42-(0.384/8)*max([1,8-len(summary_df[summary_df.Flag==0][['Sample_URL']])]), 1.0]}) # Update map
                 fig.layout['geo2']['domain']['x'] = [0, 0.56]  # fig.update_geos(domain={'y':[0.02+0.42-(0.384/8)*max([1,8-len(summary_df[summary_df.Flag==0][['Sample_URL']])]), 1.0]}) # Update map
                 fig.layout['geo']['domain']['y'] = [0.02 + 0.32 - (0.384 / 8) * max([1, 8 - len(summary_df[summary_df.Flag == 1][['Sample_URL']])]) + 0.05, 0.02 + 0.32 - (0.384 / 8) * max([1, 8 - len(summary_df[summary_df.Flag == 1][['Sample_URL']])]) + 0.25]
                 fig.layout['geo']['domain']['x'] = [0, 0.12]
                 fig.update_yaxes(domain=[0.05 + (0.768 - (0.42 / 8) * max([1, 8 - len(summary_df[summary_df.Flag == 1][['Sample_URL']])])), 0.95], row=1, col=2)  # Update scatter 1
                 fig.update_yaxes(domain=[ 0.05 + 0.42 - (0.384 / 8) * max([1, 8 - len(summary_df[summary_df.Flag == 1][['Sample_URL']])]), 0.768 - (0.42 / 8) * max([1, 8 - len(summary_df[summary_df.Flag == 1][['Sample_URL']])])], row=2, col=2)  # Update scatter 3
                 fig.update_traces(domain={'y': [0, 0.42 - (0.384 / 8) * max([1, 8 - len(summary_df[summary_df.Flag == 1][['Sample_URL']])])]}, row=3, col=1)  # Update table

             title =  r'<a href="{}">{}</a>'.format(df_standardizer['Project_source'][project_id], 'Project ID:' + str(project_id))#'Cruise:' + str(df["Cruise"][0]) + "<br> Hover on datapoint to see sample ID" if str(df["Cruise"][0]) != 'nan' else 'Project ID:' + str(project_id) + "<br> Hover on datapoint to see sample ID"
             fig.update_layout(legend=dict(y=0.95,xanchor="left",yanchor='top',x=-0.02,bgcolor='rgba(255,255,255,0)'),
                margin=dict(l=0, r=30, t=30, b=30),
                title={'text': title,'xanchor': 'center', 'yanchor': 'top', 'x': 0.5},
                xaxis={'title': '', 'showticklabels':False, 'tickfont':dict(size=10), 'titlefont':dict(size=12)},
                xaxis2={'title': 'Sample ID', 'nticks': 3, 'tickangle': 0, 'tickfont':dict(size=10), 'titlefont':dict(size=12)},
                yaxis={"tickmode": "array", 'title': '# of ROI', 'tickfont':dict(size=10), 'titlefont':dict(size=12),
                       "tickvals": pd.to_numeric([f"{n:.1g}" for n in np.power(10,np.arange(np.floor(np.log10(subset_summary_df['ROI_count'].max()))-1, np.ceil(np.log10(subset_summary_df['ROI_count'].max()))+1,1 ))]),
                       'ticktext':['10<sup>{}</sup>'.format(int(exponent)) for exponent in np.arange(np.floor(np.log10(subset_summary_df['ROI_count'].max()))-1, np.ceil(np.log10(subset_summary_df['ROI_count'].max()))+1,1 )]},
                yaxis2={'title.standoff':7000,"tickmode": "array", 'title': 'Percentage of artefacts/validation', 'tickfont':dict(size=10), 'titlefont':dict(size=12)},
                hovermode="x")
             fig.update_yaxes( type="log", row=1, col=2)
             fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
             fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
             fig.for_each_xaxis(lambda x: x.update(showgrid=False))
             fig.for_each_yaxis(lambda x: x.update(showgrid=False))
             fig.write_html(path_to_report)
             print('\nSaving cruise report to', path_to_report, sep=' ')
         else:
             print('\nNo samples left after dropping GPS coordinates and sample ID NAs. Skipping project report')
    else:
        print('\nNo project datafile found.')

# Function (4): Perform EcoTaxa export files standardization and harmonization based on standardizer spreadsheets
#standardization_func(standardizer_path='~/GIT/PSSdb/raw/project_IFCB_standardizer.xlsx',project_id=2248)
def standardization_func(standardizer_path,project_id,plot='nbss',df_taxonomy=df_taxonomy):
    """
       Objective: This function uses the instrument-specific standardizer spreadsheet to standardize variable names and units. This should facilitate the size spectrum analysis of projects collected with UVP, IFCB, and Zooscan
           :param standardizer_path: Full path of the standardizer spreadsheet containing project ID of interest
           :param project_id: unique integer for project ID to be standardized
           :param plot: Specification of the information plotted in the interactive standardization report. Options include a diversity pie chart of the taxonomic annotations ("diversity") or the sample-specific Normalized Biovolume Size Spectrum ("nbss")
           :param df_taxonomy: A taxonomic look-up table containing the native taxonomic annotations and standard taxonomy from the World Register of Marine Species. The table will be updated with new annotations if needed
    :return: standardized project dataframe and interactive standardization report
    """

    # Open standardizer spreadsheet
    sheets = [sheet for sheet in pd.ExcelFile(path_to_project_list).sheet_names if 'metadata' not in sheet]
    df_projects = pd.concat(map(lambda sheet: pd.read_excel(path_to_project_list, sheet_name=sheet), sheets))
    df_projects['Project_portal'] = df_projects.Project_localpath.apply(lambda path: Path(path).expanduser().name)
    df_projects = df_projects[df_projects['PSSdb_access'] == True]
    sheet = 0 if df_projects[df_projects.Project_ID.astype(str).isin([str(project_id)])].Project_portal.unique()[ 0] != 'ecopart' else  df_projects[df_projects.Project_ID.astype(str).isin([str(project_id)])].Project_portal.unique()[ 0]  # Default is first sheet unless project is from Ecopart
    sheetname = 'Data' if (sheet == 0) and len(re.findall(r'UVP', df_projects[df_projects.Project_ID.astype(str).isin([str(project_id)])].Instrument.unique()[0])) == 0 else 'ecotaxa' if ( len(re.findall( r'UVP', df_projects[ df_projects.Project_ID.astype( str).isin( [ str(project_id)])].Instrument.unique()[ 0])) > 0) and ( sheet == 0) else sheet


    df_standardizer = pd.read_excel(standardizer_path,sheet_name=sheet,index_col=0)

    # Control. Standardization only works for IFCB, scanners, and UVP projects at the moment
    if project_id not in df_standardizer.index:
        print('\nProject ID is not included in standardizer. Quitting')
        return
    if df_standardizer.loc[project_id]['Instrument'] not in ['Zooscan','ZooCam','IFCB','UVP','Other','Scanner']:
        print('\nStandardization only applies to Zooscan, ZooCam, IFCB, UVP and other scanners projects. Quitting')
        return

    # Retrieve fields to standardize in standardizer (indicated by variable_field)
    columns_for_description =dict(zip([item.split(':',1)[0].capitalize() for item in df_standardizer.loc[project_id]['Sampling_description'].split(';')],[eval(re.sub(r'(\w+)', r'"\1"',item.split(':',1)[1])) for item in df_standardizer.loc[project_id]['Sampling_description'].split(';')])) if str(df_standardizer.loc[project_id]['Sampling_description'])!='nan' else df_standardizer.loc[project_id]['Sampling_description']
    fields_for_description = pd.Series([columns_for_description[key][index]  for key, value in columns_for_description.items()  if (type(value)==dict) for index,fields in columns_for_description[key].items() if ('field' in index)],[key.capitalize()  for key, value in columns_for_description.items()  if (type(value)==dict) for index,fields in columns_for_description[key].items() if ('field' in index)]) if str(df_standardizer.loc[project_id]['Sampling_description'])!='nan' else pd.Series({})
    units_for_description = pd.Series([columns_for_description[key][index]  for key, value in columns_for_description.items()  if (type(value)==dict) for index,fields in columns_for_description[key].items() if ('unit' in index)],[key.capitalize()  for key, value in columns_for_description.items()  if (type(value)==dict) for index,fields in columns_for_description[key].items() if ('unit' in index)]) if str(df_standardizer.loc[project_id]['Sampling_description'])!='nan' else pd.Series({})

    columns_of_interest =  [string for string in list(df_standardizer.columns) if "_field" in string]
    fields=fields_of_interest = dict(zip(columns_of_interest,df_standardizer.loc[project_id][columns_of_interest].values.flatten().tolist()))
    fields=dict((key.replace("_field", ""), value.split(","))  if isinstance(value, str) else (key.replace("_field", ""), "".split(",")) for key, value in fields.items() )
    fields_of_interest =  dict((key.replace("_field", ""), value) for key, value in fields_of_interest.items() if isinstance(value, str))
    fields_of_interest_list = list(( key, field) for key, value in fields_of_interest.items() for field in value.split(','))
    fields_of_interest_series =pd.Series(list(value[1] for value in fields_of_interest_list),list(value[0] for value in fields_of_interest_list))

    # Retrieve flagged samples/profiles. Deprecated: Standardization is now performed regardless of the QC flags
    """
    if str(df_standardizer.loc[project_id]['Flag_path'])!='nan':
          df_flagged=pd.read_table(df_standardizer.loc[project_id]['Flag_path'],sep=',')
          flagged_samples=df_flagged.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Sample'].tolist() if len(df_flagged.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)'))>0 else ['']
    else:
          flagged_samples =['']
    """
    flagged_samples = ['']
    path_to_data = Path(df_standardizer.at[project_id, "Project_localpath"]).expanduser()
    path_files_list=list(path_to_data.rglob('**/*_{}_*'.format(str(project_id)))) if  path_to_data.stem!=cfg['UVP_consolidation_subdir'] else list( path_to_data.glob("ecotaxa_export_"+str(project_id)+'_*'))
    path_files_list=[path for path in path_files_list if ('flag' not in str(path)) & ("._" not in str(path))]
    path_to_standard_dir = Path(cfg['raw_dir']).expanduser() / cfg['standardized_raw_subdir'] / path_to_data.stem / path_files_list[0].parent.stem
    path_to_standard_plot = list(Path(path_to_git / cfg['figures_subdir'] / cfg['figures_standardizer'] / path_to_data.stem / path_files_list[0].parent.stem).glob('standardized_project_{}.html'.format(str(project_id)))) if len(list(Path(path_to_git / cfg['figures_subdir'] / cfg['figures_standardizer'] / path_to_data.stem / path_files_list[0].parent.stem).glob('standardized_project_{}.html'.format(str(project_id))))) else [Path(path_to_git / cfg['figures_subdir'] / cfg['figures_standardizer'] / path_to_data.stem / path_files_list[0].parent.stem/'standardized_project_{}.html'.format(str(project_id)))]
    summary_df_standardized_all=pd.DataFrame({})
    df_nbss=pd.DataFrame({})
    if len(path_to_standard_plot):
        with open(str(path_to_standard_plot[-1])) as f:
            html = f.read()
        call_arg_str = re.findall(r'Plotly\.newPlot\((.*)\)', html)[1]
        call_args = json.loads(f'[{call_arg_str}]')
        plotly_json = {'data': call_args[1], 'layout': call_args[2]}
        #plotly.io.from_json(json.dumps(plotly_json))
        summary_df_standardized_all=pd.DataFrame({'Latitude':plotly_json['data'][0]['lat'],'Longitude':plotly_json['data'][0]['lon'],'Sample':pd.Series(plotly_json['data'][0]['hovertext']).str.replace(r'Sample ID: ','').values})
        summary_df_standardized_all=pd.merge(summary_df_standardized_all,pd.DataFrame({'Sample':plotly_json['layout']['updatemenus'][0]['buttons'][0]['args'][0]['x'][0],'Abundance':plotly_json['layout']['updatemenus'][0]['buttons'][0]['args'][0]['y'][0]}),how='outer',on='Sample') #plotly_json['layout']['updatemenus'][0]['buttons'][0]['label']
        summary_df_standardized_all = pd.merge(summary_df_standardized_all, pd.DataFrame({'Sample': plotly_json['layout']['updatemenus'][0]['buttons'][1]['args'][0]['x'][0], 'Average_diameter': plotly_json['layout']['updatemenus'][0]['buttons'][1]['args'][0]['y'][0]}), how='outer',on='Sample') #plotly_json['layout']['updatemenus'][0]['buttons'][1]['label']
        summary_df_standardized_all = pd.merge(summary_df_standardized_all, pd.DataFrame({'Sample': plotly_json['layout']['updatemenus'][0]['buttons'][2]['args'][0]['x'][0], 'Std_diameter': plotly_json['layout']['updatemenus'][0]['buttons'][2]['args'][0]['y'][0]}), how='outer',on='Sample') #plotly_json['layout']['updatemenus'][0]['buttons'][2]['label']
        if plot=='nbss':
            df_nbss=pd.concat(map(lambda dict_args:pd.DataFrame({'Sample':dict_args['legendgroup'],'NBSS':dict_args['y'],'size_class_mid':dict_args['x']}) if 'legendgroup' in dict_args.keys() else pd.Dataframe({}),plotly_json['data'][3:]))
            df_nbss['Group_index']=pd.Categorical(df_nbss.Sample,df_nbss.Sample.unique()).codes



    if len(path_files_list)>0: # Check for native format datafile
        columns_dict=dict(zip(list(fields_of_interest_series.values),list(fields_of_interest_series.index)))
        columns_dict['object_number'] = 'object_number'  # Adding this column for UVP consolidated projects
        dtypes_dict_all = dict(zip(['Cruise', 'Station', 'Profile', 'Sample', 'Latitude', 'Longitude', 'Depth_min', 'Depth_max', 'Sampling_date', 'Sampling_time', 'Volume_analyzed', 'ROI', 'Area', 'Pixel', 'Minor', 'Biovolume', 'ESD', 'Category', 'Annotation', 'Sampling_type', 'Sampling_lower_size', 'Sampling_upper_size', 'object_number', 'Sampling_description'],[str, str, str, str, float, float, float, float, str, str, float, str, int, float, float, float, float, str, str, str, float, float, int, str]))
        dtypes_dict = {fields_of_interest_series.loc[key]: value for key, value in dtypes_dict_all.items() if key in fields_of_interest_series.index}

        # Check for existing standardized file(s):
        path_to_standard_dir.mkdir(parents=True, exist_ok=True)
        path_to_standard_file = list(path_to_standard_dir.rglob('standardized_project_{}_*.csv'.format(str(project_id))))
        path_to_standard_plot[-1].parent.mkdir(parents=True, exist_ok=True)
        if len(path_to_standard_file):
            #dtypes_dict_all['Area'] = float  # Replacing the data type of Area after converting pixel to micrometer-metric
            #df_standardized_existing = pd.concat(map(lambda path:(columns:=pd.read_table(path,nrows=0).columns,pd.read_table(path,sep=",",dtype=dtypes_dict_all,usecols=['Sample','Category']))[-1],natsorted(path_to_standard_file)))#pd.concat(map(lambda path:(columns:=pd.read_table(path,nrows=0).columns,pd.read_table(path,sep=",",dtype=dtypes_dict_all))[-1],natsorted(path_to_standard_file)))
            remaining_raw_files =[path for path in path_files_list if not Path(path_to_standard_dir / 'standardized_project_{}_{}.csv'.format(project_id,str(path.stem)[str(path.stem).rfind('_' + str(project_id) + '_') + 2 + len( str(project_id)):len(str( path.stem))].replace( '_features', ''))).exists()]
            if len(remaining_raw_files):
                print('\nExisting standardized file(s) found at {}. Generating standardized file(s) for new samples'.format(path_to_standard_dir))
            else:
                print( '\nExisting standardized file(s) found at {}, but no additional samples.\nSkipping project after saving report (if missing)'.format(path_to_standard_dir))
                """
                if path_to_standard_plot.is_file() == False:
                    #summary_df_standardized = pd.DataFrame({})
                    #for profile in df_standardized_existing.reset_index(drop=True).Sample.unique():
                        #subset_df_standardized = df_standardized_existing[df_standardized_existing.Sample == profile]
                    summary_subset_df = df_standardized_existing[(df_standardized_existing.Longitude <= 180) & (df_standardized_existing.Longitude >= -180) & ( df_standardized_existing.Latitude <= 90) & ( df_standardized_existing.Latitude >= -90)].astype(dict(zip(['Sample', 'Station', 'Latitude', 'Longitude', 'Profile', 'Volume_imaged', 'Depth_min', 'Depth_max'],[str]*8))).groupby(['Sample', 'Station', 'Latitude', 'Longitude', 'Profile', 'Volume_imaged', 'Depth_min', 'Depth_max'], dropna=True).apply(lambda x: pd.Series({'Count': x.ROI_number.astype(float).sum(), 'Abundance': x.ROI_number.astype(float).sum() / (x.Volume_imaged.astype(float).unique()[0]),# individuals per liter
                             'Average_diameter': np.nanmean(2 * np.power(np.repeat(x.Area.astype(float) ,x.ROI_number.astype(float))/ np.pi, 0.5)) if 'Area' in df_standardized_existing.columns else np.nanmean( np.repeat(x.ESD.astype(float),x.ROI_number.astype(float)())),  # micrometer
                             'Std_diameter': np.nanstd( 2 * np.power(np.repeat(x.Area.astype(float) ,x.ROI_number.astype(float)) / np.pi, 0.5)) if 'Area' in df_standardized_existing.columns else np.nanstd(np.repeat(x.ESD.astype(float),x.ROI_number.astype(float)()))})).reset_index()
                    summary_subset_df = summary_subset_df.groupby(['Sample', 'Station', 'Latitude', 'Longitude', 'Profile', 'Volume_imaged'],dropna=True).apply(lambda x: pd.Series({'Abundance': np.nanmean(x.Abundance),  # individuals per liter
                                                 'Average_diameter': np.nanmean(x.Average_diameter),  # micrometer
                                                 'Std_diameter': np.nanmean(x.Std_diameter)})).reset_index()  # micrometer

                        #summary_df_standardized = pd.concat([summary_df_standardized, summary_subset_df], axis=0, ignore_index=True).reset_index(drop=True)
                    path_to_standard_plot.parent.mkdir(parents=True, exist_ok=True)

                    if ( (plot == 'diversity') and len(summary_df_standardized) > 0):
                        # Using report function defined below
                        fig = standardization_report_func(df_summary=summary_subset_df , df_standardized=df_standardized_existing.assign(Project_source=df_standardizer['Project_source'][project_id]), df_nbss=None, plot=plot)
                    elif ((plot == 'nbss') and len(summary_df_standardized) > 0):
                        if (df_standardized_existing.Instrument.unique()[0] == 'IFCB') & len(path_to_standard_file) > 1:  # Use multi-processing to compute IFCB dashboard nbss
                            nbss = pd.DataFrame({})
                            chunk = 1000
                            with ThreadPool() as pool:
                                # for result in pool.map(lambda path: (columns := pd.read_table(path, nrows=0).columns, pd.read_table(path, usecols=[header for  header in ['object_number']+list(df_standardizer.loc[project_id, [  column  for column in df_standardizer.columns if 'field' in column]].dropna().values) if columns.isin([header]).any()]).assign(Sample_localpath=str(path).replace(str(Path.home()),'~')).rename(columns=columns_dict))[ -1], project_path, chunksize=chunk):  # pool.map(lambda path: (columns := pd.read_table(path,sep=",", nrows=0).columns, pd.read_table(path,sep=",",usecols=[column for column in ['Project_ID','Instrument','Longitude','Latitude','Sample','Sampling_date','Sampling_time','Depth_min','Depth_max','Volume_imaged','ROI','ROI_number','Area','Category'] if column in columns],dtype=dtypes_dict_all))[-1],standardized_files, chunksize=chunk):
                                for result in pool.map(lambda sample: process_nbss_standardized_files(path=None,df=df_standardized_existing.query('Sample=="{}"'.format(sample)).assign(Project_source=df_standardizer['Project_source'][project_id],type_particles=lambda x: np.where( x.ROI.isna(), 'all', 'vignettes')), category=['type_particles'],depth_selection=False)[2], df_standardized_existing.Sample.unique(), chunksize=chunk):  # pool.map(lambda path: (columns := pd.read_table(path,sep=",", nrows=0).columns, pd.read_table(path,sep=",",usecols=[column for column in ['Project_ID','Instrument','Longitude','Latitude','Sample','Sampling_date','Sampling_time','Depth_min','Depth_max','Volume_imaged','ROI','ROI_number','Area','Category'] if column in columns],dtype=dtypes_dict_all))[-1],standardized_files, chunksize=chunk):
                                    nbss = pd.concat([nbss, result], axis=0)
                            nbss = nbss.reset_index(drop=True)

                        else:
                            nbss = process_nbss_standardized_files(path=None, df=df_standardized_existing.assign(Project_source=df_standardizer['Project_source'][project_id]), category=[],depth_selection=False)[2]
                        # Using report function defined below
                        fig = standardization_report_func(df_summary=summary_subset_df , df_standardized=df_standardized_existing.assign(Project_source=df_standardizer['Project_source'][project_id]), df_nbss=nbss, plot=plot)
                    fig.write_html(path_to_standard_plot)
                    print('\nSaving standardized export plot to', path_to_standard_plot, sep=' ')
                """
                return
        else:
            df_standardized_existing = pd.DataFrame()
            remaining_raw_files=natsorted(path_files_list)
        # Load export tsv file
        df_standardized_new=pd.DataFrame({})
        format='{desc}{bar}' if len(path_files_list)>1 else '{desc}'
        print('\nSaving standardized datafile(s) to', path_to_standard_dir, sep=' ')
        with tqdm(desc='', total=len(remaining_raw_files), bar_format=format, position=0, leave=True) as bar:
            for file in natsorted(remaining_raw_files):
                df = pd.concat(map(lambda path: (columns:=pd.read_table(path,nrows=0).columns,pd.read_table(path,usecols=[header for header in ['object_number']+list(fields_of_interest_series.values) if columns.isin([header]).any()]).rename(columns=columns_dict).assign(File_path=path).astype({key:value for key,value in  dtypes_dict_all.items() if key in columns}))[-1],[file])).reset_index(drop=True)
                old_columns = df.columns
                # Standardize column names
                #df.columns = [(list(fields_of_interest_series.index)[list(fields_of_interest_series.values).index(value)]) for value in  old_columns.drop('File_path')] + ['File_path'] # pd.MultiIndex.from_arrays([[list(fields_of_interest.keys())[list(fields_of_interest.values()).index(value)] for value in df.columns],df.columns])
                index_duplicated_fields = pd.Series(fields_of_interest.keys()).isin(list(df.columns)) == False
                if any(index_duplicated_fields):
                     # Append duplicated fields
                     df = pd.concat([df.reset_index(drop=True), pd.DataFrame(dict( zip(list(pd.Series(fields_of_interest.keys())[index_duplicated_fields]), list(map(lambda item: df.loc[:, columns_dict[item]].values.tolist() if columns_dict[item] in df.columns else pd.NA, list(pd.Series( fields_of_interest.values())[index_duplicated_fields])))))).reset_index(drop=True)], axis=1)

                # Load variables used to describe the sampling collection, method, refs, etc. (Last column)
                df_method= pd.concat(map(lambda path:(columns:=pd.read_table(path,nrows=0).columns,pd.read_table(path,usecols=[header for header in [fields_of_interest_series['Sample']]+fields_for_description.values.tolist() if columns.isin([header]).any()]).rename(columns=dict(zip([header for header in [fields_of_interest_series['Sample']]+fields_for_description.values.tolist() if columns.isin([header]).any()],[(['Sample']+fields_for_description.index.tolist())[index] for index,header in enumerate(pd.Series([fields_of_interest_series['Sample']]+fields_for_description.values.tolist())) if columns.isin([header]).any()]))))[-1],[file])).drop_duplicates().reset_index(drop=True)
                df_method[fields_for_description.index] = df_method[fields_for_description.index].apply(lambda x: PA(x, dtype=ureg[dict(units_for_description)[x.name]].units) if x.name in dict(units_for_description).keys() else x)
                df_method[fields_for_description.index]=df_method[fields_for_description.index].apply(lambda x:x.name+":"+x.astype(str))

                # Remove flagged samples/profiles
                df=df[df['Sample'].astype(str).isin(flagged_samples)==False].reset_index(drop=True)
                df_method=df_method[df_method['Sample'].astype(str).isin(flagged_samples)==False].reset_index(drop=True)
                if len(df)==0:
                    continue

                # Append method description
                additional_description =' '.join([':'.join([additional_key,columns_for_description[additional_key]]) for additional_key in columns_for_description.keys() if additional_key not in fields_for_description.index]) if type(columns_for_description)==dict else []
                if (len(fields_for_description.index)>0) or (len(additional_description) > 0):
                    df_method=pd.merge(df_method,df[['Sample']],how='left',on=['Sample'])
                    df=df.assign(Sampling_description=additional_description+' ' + df_method[fields_for_description.index].apply(' '.join, axis=1) if len(additional_description) > 0 else df_method[fields_for_description.index].apply(' '.join, axis=1) )
                else:
                    df=df.assign(Sampling_description=columns_for_description)

                # Set missing values to NA
                na = str(df_standardizer.loc[project_id]['NA_value']).split(';')
                # Convert known data types to ensure masking of missing values is correct
                df=df.astype({key:value for key,value in dict(zip(['Sample','Latitude','Longitude','Volume_analyzed','Pixel'],[str,float,float,float,float])).items() if key in df.columns})
                # Convert NAs
                df = df.mask(df.apply(lambda x: x.astype(str).isin([convert(value, df.dtypes[x.name]) for value in pd.Series(na) if is_float(value) == False])))
                columns_to_convert = [column for column in df.columns if column not in ['Area', 'Depth_min', 'Depth_max', 'Minor_axis', 'Major_axis', 'ESD', 'Biovolume']]
                df[columns_to_convert] =  df.mask(df.apply(lambda x: x.astype(str).isin([convert(value, df.dtypes[x.name]) for value in pd.Series(na)])))[columns_to_convert]

                # Convert datetime and longitude format
                if 'Sampling_date' in df.columns:
                    if is_float(df.at[0, 'Sampling_date']):
                        df.loc[df['Sampling_date'].astype(float).isna() == False, 'Sampling_date'] = df.loc[df['Sampling_date'].astype(float).isna() == False, 'Sampling_date'].astype(float).astype(int).astype(str)
                        df = df.loc[df['Sampling_date'].astype(float).isna() == False]
                if 'Sampling_time' in df.columns:
                    if is_float(df.at[0, 'Sampling_time']):
                        df.loc[df['Sampling_time'].astype(float).isna() == False, 'Sampling_time'] = df.loc[df['Sampling_time'].astype(float).isna() == False, 'Sampling_time'].astype(float).astype(int).astype(str)

                df['datetime'] = pd.to_datetime(df['Sampling_date'].astype(str) + ' ' + df['Sampling_time'].astype(str).str.zfill(6),format=' '.join(df_standardizer.loc[project_id][['Sampling_date_format', 'Sampling_time_format']]),utc=True) if all(pd.Series(['Sampling_date', 'Sampling_time']).isin(df.columns)) else pd.to_datetime(df['Sampling_date'].astype(str), format=df_standardizer.at[project_id, 'Sampling_date_format'], utc=True) if 'Sampling_date' in df.columns else pd.to_datetime(df['Sampling_time'].astype(str).str.zfill(6), format=df_standardizer.at[project_id, 'Sampling_time_format'],utc=True) if 'Sampling_time' in df.columns else pd.NaT
                df['Longitude'] = (df.Longitude + 180) % 360 - 180  # Converting all longitude to [-180,180] decimal degrees

                # Transform variables based on field units
                units=[string for string in list(df_standardizer.columns) if "_unit" in string]
                units_of_interest =  dict(zip(units,df_standardizer.loc[project_id][units].values.flatten().tolist()))
                # Test whether all units are defined:                             [re.sub(r'^-?[0-9]_times_10power-?[0-9]_','',str(unit)) for unit in list(units_of_interest.values())]
                if any([unit not in ureg for unit in list(units_of_interest.values()) if str(unit)!='nan']):
                    print('\nQuitting. Please fill out standardizer with units from the list below or define your custom unit in {} using known units:\n'.format(path_to_data.parent.parent.parent /'Scripts'/'units_def.txt'),full_list_units, sep='')
                    continue

                # Convert pixel_to_size ratio in units pixel per millimeter
                if str(ureg(units_of_interest['Pixel_unit'].split('_per_')[0]).to_base_units().units)=='pixel':
                    pixel_size_ratio =  PQ(list(df.Pixel/ureg(units_of_interest['Pixel_unit'].split('_per_')[1]).to('millimeter')), 'pixel/millimeter')
                if str(ureg(units_of_interest['Pixel_unit'].split('_per_')[0]).to_base_units().units)=='meter':
                    pixel_size_ratio = PQ(list(1/(df.Pixel*ureg(units_of_interest['Pixel_unit'].split('_per_')[0]).to('millimeter'))), 'pixel/millimeter')

                image_resize_factor = np.where(df_method['image_resize_factor'].isna(),8,df_method['image_resize_factor'].isna()) if 'image_resize_factor' in df_method.columns else 8
                blob_x_grow = np.where(df_method['blob_x_grow'].isna(),20,df_method['blob_x_grow']) if 'blob_x_grow' in df_method.columns else 20
                blob_y_grow = np.where(df_method['blob_y_grow'].isna(),5,df_method['blob_y_grow']) if 'blob_y_grow' in df_method.columns else 5
                min_blob_area = np.where(df['Sampling_lower_size'].isna(),2000,df['Sampling_lower_size']) if 'Sampling_lower_size' in df.columns else 2000

                # Convert pixel to millimeter for upper/lower size
                if units_of_interest['Sampling_upper_size_unit'] == 'pixel':
                        df['Sampling_upper_size'] = df['Sampling_upper_size'] / (pixel_size_ratio.magnitude)
                        units_of_interest['Sampling_upper_size_unit'] = 'millimeter'
                if units_of_interest['Sampling_lower_size_unit'] == 'pixel':
                        df['Sampling_lower_size'] = df['Sampling_lower_size'] / (pixel_size_ratio.magnitude ** 2)
                        units_of_interest['Sampling_lower_size_unit'] = 'millimeter'
                # Set lower and upper size imaging threshold based on camera resolution and settings if missing (in pixels). # Upper threshold based on longest camera pixel dimension / Lower threshold is based on area
                if (df_standardizer.loc[project_id]['Instrument']=='IFCB') and (Path(df_standardizer.loc[project_id]['Project_localpath']).stem==cfg['IFCB_dir']):
                    df['Sampling_lower_size'] = np.nan
                    units_of_interest['Sampling_lower_size_unit'] = 'micrometer'

                camera_resolution = {'IFCB': {'Sampling_lower_size': image_resize_factor* min_blob_area, 'Sampling_upper_size':1360}, # Equivalent to minimumBlobArea/(blobXgrowAmount x blobYgrowAmount). Info stored in hdr file. See IFCB user group post: https://groups.google.com/g/ifcb-user-group/c/JYEoiyWNnLU/m/nt3FplR8BQAJ
                                     'UVP5HD': {'Sampling_lower_size': 30, 'Sampling_upper_size': 2048},
                                     'UVP5SD': {'Sampling_lower_size': 30, 'Sampling_upper_size': 1280},
                                     'UVP6': {'Sampling_lower_size': 30, 'Sampling_upper_size': 2464},
                                     'Zooscan': {'Sampling_lower_size': 631, 'Sampling_upper_size': 22640}}

                if (df_standardizer['Project_source'][project_id] == 'https://ecotaxa.obs-vlfr.fr/prj/' + str(project_id)):
                    with ecotaxa_py_client.ApiClient(configuration) as api_client:
                        api_instance = projects_api.ProjectsApi(api_client)
                        api_response = api_instance.project_query(project_id)
                    instrument = api_response['instrument']
                else:
                    instrument = df_standardizer['Instrument'][project_id]

                if 'Sampling_upper_size' not in df.columns:
                    df['Sampling_upper_size'] = np.nan
                    units_of_interest['Sampling_upper_size_unit'] = 'micrometer'
                if 'Sampling_lower_size' not in df.columns:
                    df['Sampling_lower_size'] = np.nan
                    units_of_interest['Sampling_lower_size_unit'] = 'micrometer'

                if (len(df['Sampling_upper_size'].loc[np.isnan(df['Sampling_upper_size'])])>0) and (instrument in camera_resolution.keys()):
                    df['Sampling_upper_size'].loc[np.isnan(df['Sampling_upper_size'])] = (camera_resolution[instrument]['Sampling_upper_size']/pixel_size_ratio.magnitude[0])*ureg('millimeter').to(units_of_interest['Sampling_upper_size_unit'])
                if (len(df['Sampling_lower_size'].loc[np.isnan(df['Sampling_lower_size'])])>0) and (instrument in camera_resolution.keys()):
                    df['Sampling_lower_size'].loc[np.isnan(df['Sampling_lower_size'])] = (camera_resolution[instrument]['Sampling_lower_size'] / (pixel_size_ratio.magnitude[0]**2)) * ureg('millimeter').to(units_of_interest['Sampling_lower_size_unit']) if instrument !='Zooscan' else (2*(((camera_resolution[instrument]['Sampling_lower_size'] / (pixel_size_ratio.magnitude[0]**2))/np.pi)**0.5)) * ureg('millimeter').to(units_of_interest['Sampling_lower_size_unit'])


                # Update taxonomic annotations lookup table (df_taxonomy) with new annotations
                if 'Category' in df.columns:
                    new_categories =df.dropna(subset='Category').Category[df.dropna(subset='Category').Category.isin(list(df_taxonomy.Category))==False].unique()
                    if len(new_categories):
                        try:
                            df_taxonomy_new = pd.concat( map(lambda hierarchy: annotation_in_WORMS(hierarchy.replace("_"," ")).assign(Category=hierarchy), new_categories))
                            df_taxonomy_new['URL'] = df_taxonomy_new.WORMS_ID.apply( lambda id: 'https://www.marinespecies.org/aphia.php?p=taxdetails&id={}'.format( id.replace('urn:lsid:marinespecies.org:taxname:', '')) if len(id) else '')
                            df_taxonomy_new['Life_stage'] = df_taxonomy_new.Functional_group.apply(lambda group: ';'.join([ast.literal_eval(dict)['Life stage'] for dict in group.split(';') if len(group) > 0 and 'Life stage' in ast.literal_eval(dict).keys()]))
                            df_taxonomy_new['functional_group'] = df_taxonomy_new.Functional_group.apply(lambda group: ';'.join([ dict.replace('{', '').replace(', ', ' (').replace( '}', ')').replace( "'", "") if len( group) > 0 and len( ast.literal_eval( dict)) > 1 else dict.replace( '{', '').replace( ', ', ' (').replace( '}', '').replace( "'", "") if len( group) > 0 and len( ast.literal_eval(dict)) == 1 else ''  for dict  in group.split( ';')]))
                            df_taxonomy = pd.concat([df_taxonomy.reset_index(drop=True), df_taxonomy_new[[column for column in df_taxonomy.columns if column in df_taxonomy_new.columns]].reset_index(drop=True)],axis=0, ignore_index=True)
                            df_taxonomy = df_taxonomy.sort_values(['Type', 'Category'], ascending=[False, True]).reset_index( drop=True)
                            print('New taxonomic annotations found in the project. Updating the taxonomic annotations lookup table')
                            with pd.ExcelWriter(str(path_to_taxonomy), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                                df_taxonomy.to_excel(writer, sheet_name='Data', index=False)
                        except:
                            pass


                # Use pint units system to convert units in standardizer spreadsheet to standard units
                # (degree for latitude/longitude, meter for depth, multiple of micrometer for plankton size)

                df_standardized = df.assign(Instrument=df_standardizer.loc[project_id,'Instrument'],Project_ID=project_id,
                               Cruise=df.Cruise if 'Cruise' in df.columns else pd.NA,
                               Station=df.Station if 'Station' in df.columns else pd.NA,
                               Profile=df.Profile if 'Profile' in df.columns else pd.NA,
                               Sample=df.Sample if 'Sample' in df.columns else pd.NA,
                               Sampling_date=df.datetime.dt.strftime('%Y%m%d') if all(np.isnan(df.datetime)==False) else pd.NA,
                               Sampling_time=df.datetime.dt.strftime('%H%M%S') if all(np.isnan(df.datetime)==False) else pd.NA,
                               Pixel=pixel_size_ratio.magnitude[0] if 'Pixel' in df.columns else pd.NA, # in pixels per millimeter
                               Volume_analyzed=1000*(list(df.Volume_analyzed)*ureg(units_of_interest['Volume_analyzed_unit']).to(ureg.cubic_meter))if 'Volume_analyzed' in df.columns else pd.NA, # cubic decimeter
                               Latitude=list(df.Latitude)*ureg(units_of_interest['Latitude_unit']).to(ureg.degree) if 'Latitude' in df.columns else pd.NA, # degree decimal
                               Longitude=list(df.Longitude) * ureg(units_of_interest['Longitude_unit']).to(ureg.degree) if 'Longitude' in df.columns else pd.NA,  # degree decimal
                               Depth_min=list(df.Depth_min) * ureg(units_of_interest['Depth_min_unit']).to(ureg.meter) if 'Depth_min' in df.columns else pd.NA,  # meter
                               Depth_max=list(df.Depth_max) * ureg(units_of_interest['Depth_max_unit']).to(ureg.meter) if 'Depth_max' in df.columns else pd.NA,  # meter
                               Annotation=np.where(df.Annotation.apply(lambda x:'' if str(x)=='nan' else x).astype(str).str.len()>0,df.Annotation,'unclassified') if 'Annotation' in df.columns else 'predicted' if ('Annotation' not in df.columns) and ('Category' in df.columns) else 'unclassified',
                               Category=df.Category.astype(str).apply(lambda x:'' if str(x)=='nan' else x) if 'Category' in df.columns else '',
                               Area=1e06*((df.Area*[1*ureg(units_of_interest['Area_unit']).to(ureg.square_pixel).magnitude for ntimes in range(df[['Area']].shape[1])])/((np.vstack([pixel_size_ratio**2 for ntimes in range(df[['Area']].shape[1])]).T)[0,:])) if all(pd.Series(['Area','Pixel']).isin(df.columns)) else pd.NA, # square micrometers
                               Biovolume=1e09*(df.Biovolume*[1*ureg(units_of_interest['Biovolume_unit']).to(ureg.cubic_pixel).magnitude for ntimes in range(df[['Biovolume']].shape[1])]/((np.vstack([pixel_size_ratio**3 for ntimes in range(df[['Biovolume']].shape[1])]).T)[0,:])) if all(pd.Series(['Biovolume','Pixel']).isin(df.columns)) else pd.NA, # cubic micrometers
                               Minor_axis=1e03*(df.Minor_axis*[1*ureg(units_of_interest['Minor_axis_unit']).to(ureg.pixel).magnitude for ntimes in range(df[['Minor_axis']].shape[1])]/((np.vstack([pixel_size_ratio for ntimes in range(df[['Minor_axis']].shape[1])]).T)[0,:])) if all(pd.Series(['Minor_axis','Pixel']).isin(df.columns)) else pd.NA, #  micrometers
                               Major_axis=1e03 * (df.Major_axis * [1 * ureg(units_of_interest['Major_axis_unit']).to(ureg.pixel).magnitude for ntimes in range(df[['Major_axis']].shape[1])] / ((np.vstack([pixel_size_ratio for ntimes in range(df[['Major_axis']].shape[1])]).T)[0,:])) if all( pd.Series(['Major_axis', 'Pixel']).isin(df.columns)) else pd.NA, # micrometers
                               ESD=1e03*(df.ESD*[1*ureg(units_of_interest['ESD_unit']).to(ureg.pixel).magnitude for ntimes in range(df[['ESD']].shape[1])]/((np.vstack([pixel_size_ratio for ntimes in range(df[['ESD']].shape[1])]).T)[0,:])) if all(pd.Series(['ESD','Pixel']).isin(df.columns)) else pd.NA, # micrometers
                               Sampling_type=df.Sampling_type if 'Sampling_type' in df.columns else pd.NA,
                               Sampling_lower_size=list(df.Sampling_lower_size)*ureg(units_of_interest['Sampling_lower_size_unit']).to(ureg.micrometer) if 'Sampling_lower_size' in df.columns else pd.NA, # micrometer
                               Sampling_upper_size=list(df.Sampling_upper_size)*ureg(units_of_interest['Sampling_upper_size_unit']).to(ureg.micrometer) if 'Sampling_upper_size' in df.columns else pd.NA # micrometer
                               )


                # Set volume analyzed to 5 mL for IFCB projects
                if (df_standardizer.loc[project_id]['Instrument'] in ['IFCB']) and all(df_standardized['Volume_analyzed'].isna()):
                    df_standardized['Volume_analyzed']=PQ(5*ureg('milliliter').to('liter')).magnitude

                # Convert volume analyzed to volume imaged to account for samples dilution or fractionation
                if 'Dilution' in df.columns:
                    dilution_factor = df.Dilution#np.where(df.Dilution > 1, df.Dilution, 1 / df.Dilution)
                    dilution_factor[pd.Series(dilution_factor).isna()] = 1  # Replace NA values with 1
                else:
                    dilution_factor=1

                df_standardized = df_standardized.assign(Volume_imaged=df_standardized.Volume_analyzed/dilution_factor) # cubic decimeters

                # Save data and metadata sheet
                df_standardized=pd.merge(df_standardized,df[['File_path', 'Sample']].drop_duplicates(),how='left',on='Sample')
                df_standardized =df_standardized.astype({key:value for key,value in dict(zip(['Sample','Profile','Cruise','Station','Sampling_date','Sampling_time', 'Category', 'Annotation'],8*[str])).items() if key in df_standardized.columns})

                df_standardized['object_number']=df_standardized['object_number'] if 'object_number' in df_standardized.columns else np.repeat(1,len(df_standardized))
                df_standardized=df_standardized.rename(columns={'object_number':'ROI_number'})
                df_standardized=df_standardized[['Project_ID','Cruise','Instrument','Sampling_type', 'Station', 'Profile','Sample', 'Latitude', 'Longitude', 'Sampling_date', 'Sampling_time','Depth_min', 'Depth_max', 'Volume_analyzed', 'Volume_imaged', 'ROI','ROI_number', 'Annotation','Category', 'Minor_axis', 'Major_axis', 'ESD', 'Area', 'Biovolume','Pixel','Sampling_lower_size','Sampling_upper_size','Sampling_description']]
                df_standardized.to_csv( path_to_standard_dir / 'standardized_project_{}_{}.csv'.format(project_id,str(file.stem)[str(file.stem).rfind('_' + str(project_id) + '_') + 2 + len( str(project_id)):len(str( file.stem))].replace( '_features', '')), sep=",", index=False)
                df_standardized_new=pd.concat([df_standardized_new,df_standardized],axis=0)
                # Append summary
                df_standardized = df_standardized[(df_standardized.Longitude <= 180) & (df_standardized.Longitude >= -180) & (df_standardized.Latitude <= 90) & (df_standardized.Latitude >= -90)]

                if len(df_standardized):
                    df_standardized['Profile'] = df_standardized['Profile'].astype(str)  # Required to pass groupby if missing
                    df_standardized['Station'] = df_standardized['Station'].astype(str)  # Required to pass groupby if missing
                    df_standardized['Sample'] = df_standardized['Sample'].astype( str)  # Required to pass groupby if missing
                    group = ['Sample', 'Station', 'Latitude', 'Longitude', 'Profile'] if 'UVP' in instrument else ['Sample', 'Station', 'Latitude', 'Longitude', 'Profile', 'Volume_imaged', 'Depth_min', 'Depth_max']

                    df_standardized_summary = df_standardized.astype(dict(zip(group,[str]*len(group)))).groupby(group).apply(lambda x: pd.Series({'Cumulative_volume_imaged':x[['Sample','Volume_imaged']].drop_duplicates().Volume_imaged.sum(),'Depth_range_min': x.Depth_min.astype(float).min(),'Depth_range_max': x.Depth_max.astype( float).max()})).reset_index()
                    df_standardized = pd.merge(df_standardized, df_standardized_summary, how='left', on=group)
                    summary_df_standardized = df_standardized.groupby(['Sample', 'Station', 'Latitude', 'Longitude', 'Profile', 'Cumulative_volume_imaged','Depth_range_min', 'Depth_range_max'], dropna=True).apply(lambda x: pd.Series({'Count': x.ROI_number.sum(),
                         'Abundance': x.ROI_number.sum() / (x.Cumulative_volume_imaged.unique()[0]),# individuals per liter
                         'Average_diameter': np.nanmean(2 * np.power(np.repeat(x.Area, x.ROI_number) / np.pi,0.5)) if 'Area' in df_standardized.columns else np.nanmean( np.repeat(x.ESD, x.ROI_number)),  # micrometer
                         'Std_diameter': np.nanstd(2 * np.power(np.repeat(x.Area, x.ROI_number) / np.pi,0.5)) if 'Area' in df_standardized.columns else np.nanstd(np.repeat(x.ESD, x.ROI_number))})).reset_index()
                    summary_df_standardized = summary_df_standardized.groupby( ['Sample', 'Station', 'Latitude', 'Longitude', 'Profile'], dropna=True).apply(lambda x: pd.Series({'Abundance': np.nanmean(x.Abundance),  # individuals per liter
                                             'Average_diameter': np.nanmean(x.Average_diameter),  # micrometer
                                             'Std_diameter': np.nanmean(x.Std_diameter)})).reset_index()  # micrometer

                    summary_df_standardized_all=pd.concat([summary_df_standardized_all,summary_df_standardized],axis=0).reset_index(drop=True)
                    if plot=='nbss':
                        nbss_sample=pd.concat(map(lambda sample: process_nbss_standardized_files(path=None, df=df_standardized.astype({'Category':str,'Sampling_type':str}).query('Sample=="{}"'.format(sample)).assign(type_particles=lambda x: np.where(x.ROI.isna(), 'all', 'vignettes')), category=['type_particles'],depth_selection=False)[2],df_standardized.Sample.unique()))
                        df_nbss=pd.concat([df_nbss,nbss_sample],axis=0).reset_index(drop=True)
                bar.set_description('Saving standardized datafile(s)', refresh=True)
                ok = bar.update(n=1)



        # with pd.ExcelWriter(str(path_to_standard_file),engine="xlsxwriter") as writer:
        # df_standardized.to_csv(path_to_standard_file, sep="\t",index=False) # to_excel(writer, sheet_name='Data', index=False)
        # df_standardized_metadata.to_csv(path_to_metadata_standard, sep="\t",index=False)  # to_excel(writer, sheet_name='Metadata',index=False)
        df_standardized =df_standardized_new.reset_index(drop=True) #pd.concat([df_standardized_existing.reset_index(drop=True), df_standardized_new.reset_index(drop=True)], axis=0,ignore_index=True)
        df_standardized_metadata=pd.DataFrame({'Variables':df_standardized.columns,'Variable_types':df_standardized.dtypes,
        'Units/Values/Timezone':['','','','','','','','degree','degree','yyyymmdd (UTC)','hhmmss (UTC)','meter','meter','cubic_decimeter','cubic_decimeter','','','','']+['micrometer' for ntimes in range(df_standardized[['Minor_axis']].shape[1])]+['micrometer' for ntimes in range(df_standardized[['Major_axis']].shape[1])]+['micrometer' for ntimes in range(df_standardized[['ESD']].shape[1])]+['square_micrometer' for ntimes in range(df_standardized[['Area']].shape[1])]+['cubic_micrometer' for ntimes in range(df_standardized[['Biovolume']].shape[1])]+['pixel_per_millimeter','micrometer','micrometer',''],
        'Description':['Project ID','Project cruise','Instrument','Sampling type (e.g. platform, gear, strategy)','Station ID (native format)','Profile ID (native format)','Sample ID (native format)','Latitude','Longitude','Sampling date','Sampling time','Minimum sampling depth','Maximum sampling depth','Volume analyzed (not accounting for sample dilution and/or fractionation)','Volume imaged (accounting for sample dilution and/or fractionation)','Region of interest ID (native format)','Number of ROI for a given area estimate (equals to 1 for all scanner, IFCB projects and UVP large particles)','Status of ROI annotation (e.g. unclassified, predicted, validated)','ROI assigned taxonomy']+ ['Object minor ellipsoidal axis derived from '+ field if field!='' or len(fields['Minor_axis'])>1 else 'Object minor ellipsoidal axis' for field in fields['Minor_axis'] ]+ ['Object major ellipsoidal axis derived from '+ field if field!='' or len(fields['Major_axis'])>1 else 'Object major ellipsoidal axis' for field in fields['Major_axis'] ]+ ['Object equivalent spherical diameter derived from '+ field if field!='' or len(fields['ESD'])>1 else 'Object equivalent spherical diameter' for field in fields['ESD'] ]+ ['Object surface area derived from '+ field if field!='' or len(fields['Area'])>1 else 'Object surface area' for field in fields['Area'] ]+['Object biovolume derived from '+ field if field!='' or len(fields['Biovolume'])>1 else 'Object biovolume' for field in fields['Biovolume']]+['Pixel size used for size conversion','Smallest sampled size','Largest sampled size','Additional description of the sampling method or protocol']})
       # Create readme file:
        if not Path(path_to_standard_dir.parent.parent / 'README.txt').is_file():
            with open(str(Path(path_to_standard_dir.parent.parent / 'README.txt')), 'w') as file:
                file.write( "README file for project standardized files (First created on February 10, 2023):\n\nThis directory contains the standardized table(s) of accessible projects.\nEach table include the following variables:\n\n{}\n\nContact us: nmfs.pssdb@noaa.gov".format('\n'.join(list(df_standardized_metadata[['Variables', 'Units/Values/Timezone', 'Description']].apply(lambda x: str(x.Variables) + " (" + str(x['Units/Values/Timezone']).strip() + "): " + x.Description if len(x['Units/Values/Timezone']) else str(x.Variables) + ": " + x.Description, axis=1).values))))

        # Interactive plots
        if ((plot=='diversity') and len(summary_df_standardized)>0):
            # Using report function defined below
            fig=standardization_report_func(df_summary=summary_df_standardized_all,df_standardized=df_standardized.dropna(subset=['ROI']).assign(Project_source=df_standardizer['Project_source'][project_id]),df_nbss=None,plot=plot)
            fig.write_html(path_to_standard_plot)
            print('Saving standardized export plot to', path_to_standard_plot, sep=' ')
        elif ((plot=='nbss') and len(summary_df_standardized_all)>0):
            df_nbss=df_nbss.assign(Instrument=df_standardizer.at[project_id,'Instrument'],Project_ID=project_id)
            group=['Instrument', 'Project_ID', 'Sample']
            nbss=pd.merge(df_nbss.drop(columns=['Group_index']), df_nbss.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename( {'index': 'Group_index'}, axis='columns'), how='left', on=group)
            # Using report function defined below
            if (instrument == 'IFCB') & (len(path_files_list)> 1): # group dataset by 6 month intervals and save
                summary_df_standardized_all['Sample_datetime']=pd.to_datetime(summary_df_standardized_all.Sample.str[0:7],format='D%Y%m')
                summary_df_standardized_all['Datetime_bin'] = summary_df_standardized_all.Sample_datetime.dt.strftime("%Y-")+pd.cut(summary_df_standardized_all.Sample_datetime.dt.strftime("%m").astype(float),[0,6,12],labels=['01','06']).astype(str)
                nbss['Sample_datetime']=pd.to_datetime(nbss.Sample.str[0:7],format='D%Y%m')
                nbss['Datetime_bin'] = nbss.Sample_datetime.dt.strftime("%Y-") + pd.cut( nbss.Sample_datetime.dt.strftime("%m").astype(float), [0, 6, 12], labels=['01', '06']).astype(str)

                for increment,bin in enumerate(summary_df_standardized_all.Datetime_bin.unique()):
                    fig = standardization_report_func(df_summary=summary_df_standardized_all.query('Datetime_bin=="{}"'.format(bin)),df_standardized=pd.DataFrame({'Project_source' : df_standardizer['Project_source'][project_id],'Project_ID':project_id,'Instrument':df_standardizer['Instrument'][project_id]},index=[0]),df_nbss=nbss.query('Datetime_bin=="{}"'.format(bin)), plot=plot)
                    fig.write_html(Path(path_to_standard_plot[0].parent)/'standardized_project_{}_{}.html'.format(project_id,bin))

            else:
                fig = standardization_report_func(df_summary=summary_df_standardized_all,df_standardized=df_standardized.assign( Project_source=df_standardizer['Project_source'][project_id]),df_nbss=nbss, plot=plot)
                fig.write_html(Path(path_to_standard_plot[0].parent)/'standardized_project_{}.html'.format(project_id))
            print('\nSaving standardized export plot to', path_to_standard_plot[0].parent, sep=' ')

        else:
            print('\nNo samples left after dropping samples NAs. Skipping standardization report')
    else: print('\nExport file for project {} not found. Please run 1.export_projects.py to export file automatically or export manually'.format(str(project_id)))

def standardization_report_func(df_summary,df_standardized,df_nbss,plot='diversity'):

    if plot == 'diversity':
        subplot_specs = [[{"type": "scattergeo", "rowspan": 2}, {"type": "scatter"}], [None, {"type": "pie"}]]
        subplot_titles = ['Map of the projet samples:', '', 'Annotations diversity']
    if plot == 'nbss':
        subplot_specs = [[{"type": "scattergeo", "rowspan": 2}, {"type": "scatter"}], [None, {"type": "scatter"}]]
        subplot_titles = ['', '', '']

    fig = make_subplots(rows=2, cols=2,specs=subplot_specs, subplot_titles=subplot_titles, column_widths=[0.5, 0.5], row_heights=[0.3, 0.7], vertical_spacing=0.1)
    # Map, left panel
    if len(df_summary):
        data_geo = dict(type='scattergeo',
                        name='',
                        lon=df_summary.dropna(subset=['Latitude', 'Longitude', 'Sample']).Longitude,
                        lat=df_summary.dropna(subset=['Latitude', 'Longitude', 'Sample']).Latitude,
                        # size=summary_df_standardized.Abundance,
                        hovertext='Sample ID: ' + df_summary.dropna( subset=['Latitude', 'Longitude', 'Sample']).Sample,
                        hoverinfo="text",
                        marker=dict(color='black', size=0.5), geojson="natural earth", showlegend=False, geo='geo')
        layout = dict()
        layout['geo2'] = dict(projection_rotation=dict(
                lon=np.nanmean(df_summary.dropna(subset=['Latitude', 'Longitude', 'Sample']).Longitude),
                lat=np.nanmean(df_summary.dropna(subset=['Latitude', 'Longitude', 'Sample']).Latitude),
                roll=0),center=dict(lon=np.nanmean(df_summary.dropna(subset=['Latitude', 'Longitude', 'Sample']).Longitude),
                lat=np.nanmean(df_summary.dropna(subset=['Latitude', 'Longitude', 'Sample']).Latitude)),
            projection_type='orthographic', lataxis_range=[-90, 90], lonaxis_range=[-180, 180],
            domain=dict(x=[0.02, 0.45], y=[0.0, 0.9]), bgcolor='rgba(0,0,0,0)')
        layout['geo'] = dict(lonaxis_range=[-180, 180], lataxis_range=[-80, 80],domain=dict(x=[0.0, 0.1], y=[0.0, 0.12]))
        # go.Figure(data=[data_geo, data_geo], layout=layout)
        fig.add_trace(data_geo, row=1, col=1).update_layout(go.Figure(data=[data_geo, data_geo], layout=layout).layout)
        data_geo.update({'geo': 'geo2', 'showlegend': False, 'marker': dict(color='black', size=2.5)})  # Update geo trace before adding the zoomed
        fig.add_trace(data_geo, row=1, col=1)
        fig.data[1]['geo'] = 'geo2'  # Need to update geo trace manually

        # fig.add_trace(go.Scattergeo(name='',lon=summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Longitude, lat=summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Latitude, hovertext='Sample ID: ' + summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Sample, hoverinfo="text", marker=dict(color='black'), geojson="natural earth", showlegend=False),row=1, col=1).update_geos(projection_rotation=dict(lon=np.nanmean(summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Longitude), lat=np.nanmean(summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Latitude), roll=0),center=dict(lon=np.nanmean(summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Longitude), lat=np.nanmean(summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Latitude)),projection_type='orthographic')
        # Scatterplot 1, top-right panel
        cols = df_summary.columns
        fig.add_trace(go.Scatter(name='', x=df_summary.dropna(subset=['Sample', 'Abundance']).Sample,
                                 y=df_summary.dropna(subset=['Sample', 'Abundance']).Abundance,
                                 #  error_y=dict(type="data",array=summary_df_standardized.Std_diameter),
                                 hovertext="Sample ID: " + df_summary.dropna( subset=['Sample', 'Abundance']).Sample,
                                 hoverinfo="text",
                                 marker=dict(color='black'), mode='markers',
                                 showlegend=True, visible=True), row=1, col=2)

    # Scatterplot 2, bottom-right panel
    if plot=='nbss':
        colors = sequential_hcl(h=[-220, 20], c=[0, 90, 0], l=[95, 30], power=[1.5, 1.5], rev=True)(len(df_nbss['Sample'].unique()))
        for i, sample in enumerate(np.sort(df_nbss.Sample.unique())):
            subset_nbss = df_nbss[df_nbss['Sample'] == sample]  # nbss#
            for v in subset_nbss.Group_index.unique():  # subset_nbss.Volume_imaged.unique():
                nbss = subset_nbss[subset_nbss['Group_index'] == v].sort_values(['size_class_mid']).dropna()  # subset_nbss[subset_nbss['Volume_imaged'] == v].sort_values(['bin_gmean']).dropna()
                show = True if v == subset_nbss.Group_index.unique()[0] else False

                fig.add_trace(go.Scatter(y=nbss.NBSS, x=nbss.size_class_mid,
                         line=dict(color=colors[i]),
                         name=str(sample),
                         hovertext='Sample ID: ' + nbss.Sample.astype(str) + '<br>Size bin: ' + np.round( nbss['size_class_mid'].astype(float), 1).astype(str) + ' \u03BCm<br>NBSS: ' + np.round( nbss['NBSS'].astype(float), 6).astype(str) + ' \u03BCm\u00b3 dm\u207B\u00b3 \u03BCm\u207B\u00b3',
                         hoverinfo="text",
                         legendrank=i,legendgroup=str(sample),
                         mode='lines', showlegend=show, visible=True), row=2, col=2)

    # Scatterplot 2, bottom-right panel
    if plot == 'diversity':
        subset_df = pd.merge(df_standardized, df_taxonomy[['Category', 'EcoTaxa_hierarchy', 'Full_hierarchy', 'Rank', 'Type', 'Domain'] + [rank for rank in rank_categories if rank in df_taxonomy.columns]], left_on='Category', right_on='Category', how='left')
        columns_ID = ['Project_ID']  # ['sample_profileid','sample_stationid','object_lat','object_lon','Depth_bin_med']
        column_group = ['Genus']
        group_categories = df_taxonomy[df_taxonomy.EcoTaxa_hierarchy.isin(subset_df.EcoTaxa_hierarchy.tolist())][ [column for column in df_taxonomy.columns if column in rank_categories and column in np.array(rank_categories)[np.where(pd.Categorical(rank_categories, categories=rank_categories, ordered=True) <= column_group[0])[0].tolist()]]]
        group_categories = group_categories.sort_values(by=column_group[0])
        group_categories[group_categories.isnull()] = ''
        if len(subset_df.dropna(subset=['Full_hierarchy'])) > 0:
            summary_df_taxonomy = subset_df.dropna(subset=['Full_hierarchy']).groupby(by=columns_ID + ['Phylum'],observed=True).apply( lambda x: pd.Series({'ROI_count': x.ROI.count()})).reset_index().pivot(index=columns_ID,columns='Phylum', values='ROI_count').reset_index()
            melt_df_taxonomy = summary_df_taxonomy.melt( value_vars=[column for column in summary_df_taxonomy.columns if column not in columns_ID],id_vars=columns_ID, var_name='Group', value_name='Count').dropna()
            melt_df_taxonomy = pd.merge(melt_df_taxonomy, group_categories, left_on='Group', right_on='Phylum', how='left')
            melt_df_taxonomy.Group = pd.Categorical(melt_df_taxonomy.Group, categories=group_categories['Phylum'].unique(), ordered=True)
            melt_df_taxonomy.Group.cat.categories
            melt_df_taxonomy = melt_df_taxonomy.sort_values(by=['Group'] + columns_ID)
            levels = pd.Categorical(melt_df_taxonomy.columns[melt_df_taxonomy.columns.isin(rank_categories)].tolist(),categories=pd.Series(rank_categories)[pd.Series(rank_categories).isin(melt_df_taxonomy.columns.tolist())].tolist(),ordered=True)
            colors_dict = taxon_color_palette(dataframe=group_categories, levels=levels, palette=px.colors.qualitative.T10)
            levels_dict = dict(zip(levels, px.colors.sequential.Reds[0:len(levels)][::-1]))
            # fig=px.pie(melt_df_taxonomy,values='Count',names='Group',color='Group',hole=0.3,color_discrete_map=colors_dict).update_traces(sort=False,textinfo='percent+label').update_layout(title_text="Annotation diversity",annotations=[dict(text='% ROI', x=0.5, y=0.5, font_size=20, showarrow=False)])
            level_new_data = melt_df_taxonomy[columns_ID + ['Group', 'Count']].drop_duplicates().groupby(by=columns_ID).apply(lambda x: pd.Series({'Total': x.Count.sum(), 'Rank': 'Phylum'})).reset_index()
            fig.add_trace(go.Pie(name='Phylum', labels=level_new_data['Rank'], values=level_new_data['Total'], hoverinfo='none',
                       textinfo='none', direction='clockwise', hole=0.2, legendgroup=2, showlegend=True, sort=False,
                       marker=dict(colors=['rgba(255,255,255,0)'], line=dict(color=levels_dict['Phylum'], width=6)),
                       textposition='inside', insidetextorientation='radial'), row=2, col=2).update_layout(
                       annotations=[dict(text='% ROI', x=0.78, y=0.31, font_size=12, showarrow=False), dict(text='Pie chart of the annotations diversity ' + str((np.round(100 * subset_df['Type'].value_counts(normalize=True), 1).astype(str) + '%').to_dict()), x=0.7,y=-0.05, font_size=12, showarrow=False)])
            fig.add_trace(go.Pie(name='Taxonomic phylum:', labels=melt_df_taxonomy[['Group', 'Count']].drop_duplicates()['Group'],
                       values=melt_df_taxonomy[['Group', 'Count']].drop_duplicates()['Count'], hole=0.2,
                       direction='clockwise', legendgroup=1, sort=False, hoverinfo='percent+label',
                       textinfo='percent+label', textfont_size=10, textposition='none',marker=dict(colors=[color for taxon, color in colors_dict.items() if taxon in melt_df_taxonomy[['Group', 'Count']].drop_duplicates()['Group'].values.tolist()])), row=2, col=2)
            levels = df_taxonomy.columns[df_taxonomy.columns.isin(rank_categories)].tolist()
            for index, level in enumerate(levels):
                if level == 'Phylum':
                    continue
                else:
                    summary_df_taxonomy = subset_df.dropna(subset=['Full_hierarchy']).groupby( by=columns_ID + levels[0:index + 1], observed=True, dropna=False).apply(lambda x: pd.Series({'ROI_count': x.ROI.count()})).reset_index()
                    summary_df_taxonomy = summary_df_taxonomy[summary_df_taxonomy['Phylum'].isnull() == False]
                    summary_df_taxonomy.loc[:, levels[0:index + 1]] = summary_df_taxonomy.loc[:,levels[0:index + 1]].fillna('unassigned')
                    summary_df_taxonomy = pd.merge(summary_df_taxonomy,summary_df_taxonomy.groupby(by=columns_ID + [levels[index - 1]], observed=True, group_keys=False).apply( lambda x: pd.Series({'Total_ROI_count': x.ROI_count.sum()})).reset_index()[[levels[index - 1], 'Total_ROI_count']], on=levels[index - 1],how='left')
                    summary_df_taxonomy['Phylum'] = pd.Categorical(summary_df_taxonomy['Phylum'], categories=group_categories['Phylum'].unique(),ordered=True)
                    new_data = summary_df_taxonomy  # melt_df_taxonomy.groupby(level).apply(lambda x :pd.Series({'Count': sum(x.Count)})).reset_index()
                    new_data = new_data.sort_values(by=levels[0:index])
                    new_colors = [colors_dict[taxon] if taxon != 'unassigned' else '#FFFFFF' for taxon in new_data[ new_data.columns.tolist()[index + 1]]]  # [color for taxon,color in colors_dict.items() if taxon in new_data[new_data.columns.tolist()[2]].values.tolist()]
                    new_data['Hierarchy'] = new_data.apply(lambda x: '>'.join(x[levels[0:index + 1]][levels[0:index + 1]] + ' (' + levels[0:index + 1] + ')'), axis=1)
                    level_new_data = new_data.groupby(by=columns_ID).apply(lambda x: pd.Series({'Total': x.ROI_count.sum(), 'Rank': level})).reset_index()
                    fig.add_trace(go.Pie(name=level, labels=level_new_data['Rank'], values=level_new_data['Total'],
                                         hoverinfo='none', textinfo='none', direction='clockwise',
                                         hole=0.2 + index * (1 / (len(levels) + 1)), legendgroup=2, showlegend=True,
                                         sort=False,
                                         marker=dict(colors=['#FFFFFF'], line=dict(color=levels_dict[level], width=6)),
                                         textposition='inside', insidetextorientation='radial'), row=2, col=2)
                    fig.add_trace(go.Pie(name=level, labels=new_data['Hierarchy'], values=new_data['ROI_count'],
                                         hoverinfo='label+percent', textinfo='none', direction='clockwise',
                                         hole=0.2 + index * (1 / (len(levels) + 1)), showlegend=False, sort=False,
                                         marker=dict(colors=new_colors, line=dict(color='#000000', width=0)),
                                         textposition='inside', insidetextorientation='radial'), row=2, col=2)

    cols_labels = np.where(cols.isin(['Abundance']), 'Abundance (particles dm<sup>-3</sup>)', '')
    cols_labels = np.where(cols.isin(['Average_diameter']), u'Average equivalent circular diameter (\u03BCm)', cols_labels)
    cols_labels = np.where(cols.isin(['Std_diameter']), u'Std equivalent circular diameter (\u03BCm)', cols_labels)
    button_scatter1 = [dict(method="restyle",  # argument to change data
                        args=[{'x': [df_summary['Sample'], 'undefined'],
                               'y': [df_summary[cols[k]], 'undefined'],
                               'visible': [True, True, True]}, [2]],# visible should have the same length as number of subplots, second argument specifies which subplot is updated
                        label=cols_labels[k]) for k in np.where(np.isin(list(cols), ['Abundance', 'Average_diameter', 'Std_diameter']))[0]]
    title =  r'<a href="{}">{}</a>'.format(df_standardized['Project_source'].unique()[0], 'Project ID:' + str(df_standardized['Project_ID'].unique()[0]))
    fig.update_layout(updatemenus=list([dict(active=0, buttons=button_scatter1, x=0.87, y=1.1, xanchor='left', yanchor='top')]),
    title={'text': title, 'xanchor': 'center', 'yanchor': 'top', 'x': 0.5},
    xaxis={'title': 'Sample ID', 'nticks': 2, 'tickangle': 0, 'tickfont': dict(size=10), 'titlefont': dict(size=12)},
    xaxis2={'title': u'Equivalent circular diameter (\u03BCm)', 'tickfont': dict(size=10), 'titlefont': dict(size=12),
            "tickvals": np.sort(np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(0, np.ceil(np.log10(df_nbss['size_class_mid'].max())), 1)))).tolist() if plot == 'nbss' else [''],
            'ticktext': [size if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(0, np.ceil(np.log10(df_nbss['size_class_mid'].max())), 1))))] if plot == 'nbss' else ['']},
    yaxis={"tickmode": "array", 'title': 'Variable selected in <br> dropdown menu ', 'tickfont': dict(size=10),'titlefont': dict(size=12)},
    yaxis2={"tickmode": "array",'title': u'Normalized Biomass Size Spectrum <br> (\u03BCm\u00b3 dm\u207B\u00b3 \u03BCm\u207B\u00b3)','tickfont': dict(size=10), 'titlefont': dict(size=12)}, boxmode='group')
    fig.update_yaxes(type="log", row=1, col=2)
    fig.update_yaxes(type="log", row=2, col=2)
    fig.update_xaxes(type="log", row=2, col=2, ticks="inside")
    fig.for_each_xaxis(lambda x: x.update(showgrid=False))
    fig.for_each_yaxis(lambda x: x.update(showgrid=False, ticks="inside"))
    return fig

# Function (5): Update flag summaries to account for overruled sample(s)
def update_summaries(path_to_flags=path_to_git / cfg['dataset_subdir'] / cfg['flag_subdir'],path_to_summaries=path_to_git / cfg['dataset_subdir'] / cfg['summary_subdir']):
    """
        Objective: This function update the summary files(=project-specific table including a summary of individual sample source, path, locations, date, study area, flag) according to individual flag files(=project-specific table including a summary of individual sample source, path, overall flag, individual flags, and overruling boolean factor)
           :param path_to_flags: Full path of the subdirectory containing all flag datafiles
           :param path_to_summaries: Full path of the subdirectory containing all summary datafiles
        :return: updated project summary datafiles at path_to_summaries/[data source/ instrument]/Summary_project_[Project_ID].csv
    """

    flag_datafiles=list(Path(path_to_flags).expanduser().rglob('project_*_flags.csv'))
    summary_datafiles=list(Path(path_to_summaries).expanduser().rglob('Summary_project_*.csv'))
    remaining_flags = flag_datafiles
    with tqdm(desc='Updating project-specific summary table', total=len(summary_datafiles), bar_format='{desc}{bar}', position=0, leave=True) as bar:
        for path in summary_datafiles:
            df_summary=pd.read_table(path,sep=",",parse_dates=['Sample_datetime'],date_parser=lambda date: datetime.datetime.strptime(date,'%Y-%m-%dT%H%M%SZ'))
            percent=np.round(100*bar.n/len(summary_datafiles))
            bar.set_description("Updating project-specific summary table {} ({}%)".format(path.stem.replace('Summary_',''),percent) , refresh=True)
            path_flag=Path(str(path).replace( str(path_to_summaries) + (str(path).replace(str(path_to_summaries), '').replace(path.name, '')),str(path_to_flags) + os.sep+str(str(Path(str(path).replace(str(path_to_summaries),'').replace(path.name,''))).split(os.sep)[1])).replace('Summary_project_', '/project_').replace('.csv', '_flags.csv'))
            if path_flag.exists():
                df_flag=pd.read_table(path_flag,sep=",")
                df_flag['Sample_flag']=np.where(((df_flag.Flag==0) & (df_flag.Overrule==False)) | ((df_flag.Flag==1) & (df_flag.Overrule==True)),0,1)
                df_summary=pd.merge(df_summary,df_flag[['Sample','Sample_flag']],how='left',on='Sample') if 'Sample_flag' not in df_summary.columns else pd.merge(df_summary.drop(columns='Sample_flag'),df_flag[['Sample','Sample_flag']],how='left',on='Sample')
                df_summary=df_summary.sort_values(['Sample_datetime'])
                df_summary['Sample_datetime']=df_summary['Sample_datetime'].dt.strftime('%Y-%m-%dT%H%M%SZ').astype(str)
                df_summary.to_csv(path,index=False)
                remaining_flags.remove(path_flag)
            else:
                print('Flag table not found for {} project {}. Skipping'.format(str(path).replace(str(path_to_summaries),'').replace(path.name,''),path.stem.split('_')[2]))
            ok = bar.update(n=1)
    if len(remaining_flags):
        sheets_project_list = [sheet for sheet in pd.ExcelFile(path_to_project_list).sheet_names if 'metadata' not in sheet]
        columns_project_list = dict(map(lambda sheet: (sheet, pd.read_excel(path_to_project_list, sheet_name=sheet).columns.tolist()),sheets_project_list))
        df_projects_metadata = pd.read_excel(path_to_project_list, sheet_name='metadata')
        df_projects = pd.concat(map(lambda sheet: pd.read_excel(path_to_project_list, sheet_name=sheet).assign(Portal=sheet),sheets_project_list)).reset_index(drop=True)
        with tqdm(desc='Summary table not found', total=len(remaining_flags),bar_format='{desc}{bar}', position=0, leave=True) as bar:
            for path in remaining_flags:
                project_id=path.stem.replace('project_','').replace('_flags','')
                percent = np.round(100 * bar.n / len(remaining_flags))
                bar.set_description( "Generating project-specific summary table {} ({}%)".format(project_id, percent), refresh=True)
                df_project = df_projects[(df_projects.Project_ID.astype(str) == str(project_id)) & (df_projects.Portal.astype(str).isin(['ecotaxa', 'ifcb']))]
                quality_control_func( standardizer_path=path_to_flags.parent / 'project_{}_standardizer.xlsx'.format(re.split(r"\s+|\d+", df_project.Instrument.unique()[0])[0]), project_id=project_id, report_path=path_to_git / 'reports' ,validation_threshold=0.95)
                bar.update(n=1)



if __name__ == '__main__':
    print('Assisting standardizer completion. \nPlease enter the following information (attention: method sensitive to empty space and lower/upper case. Do not use quotation)\n')
    funcs_args=input('Press:\n 1 for printing all fields for Ecotaxa project using the API (downloading export files is not required)\n 1b for printing all possible units\n 2 for generating/saving interactive plot of raw ecotaxa export\n 3 for flagging project\n 4 for generating/saving/plotting project export standardization\n 5 for updating all summary files according to latest flag datafiles')

    if funcs_args == '1b':
        unitfile_args=input('Custom units full path:')
        if len(unitfile_args)>0:
            filling_standardizer_unit_func(custom_units_path=unitfile_args)
        else:
            filling_standardizer_unit_func()
        quit()

    standardizerfile_args = input('Standardizer spreadsheet full path:')
    # Example: ~/GIT/PSSdb/raw/project_Zooscan_standardizer.xlsx

    project_args=input('Unique project ID of interest (integer):')
    # Example: 377
    if funcs_args=='1':
        configfile_args = input('Configuration file (containing EcoTaxa authentication) full path:')
        # Example: ~/GIT/PSSdb/scripts/configuration_masterfile_pw.yaml
        if len(configfile_args) > 0 and len(standardizerfile_args) > 0:
            if len(project_args) == 0:
                filling_standardizer_field_func(config_path=configfile_args, standardizer_path=standardizerfile_args)
            else:
                filling_standardizer_field_func(config_path=configfile_args, standardizer_path=standardizerfile_args,project_id=int(project_args))
        else: print('Missing information. Please run again')


    if funcs_args=='2':
        if len(project_args) > 0 and len(standardizerfile_args) > 0:
            diagnostic_plots_func(standardizer_path=standardizerfile_args,project_id=project_args)
        else:
            print('Missing information. Please run again')

    if funcs_args == '3':
        if len(project_args) > 0 and len(standardizerfile_args) > 0:
             quality_control_func(standardizer_path=standardizerfile_args, project_id=project_args, report_path=path_to_git / 'reports')
        else:
             print('Missing information. Please run again')

    if funcs_args=='4':
        if len(project_args) > 0 and len(standardizerfile_args) > 0:
            standardization_func(standardizer_path=standardizerfile_args, project_id=project_args)
        else:
            print('Missing information. Please run again')

    if funcs_args=='5':
        update_summaries(path_to_flags=path_to_git / cfg['dataset_subdir'] / cfg['flag_subdir'],path_to_summaries=path_to_git / cfg['dataset_subdir'] / cfg['summary_subdir'])
