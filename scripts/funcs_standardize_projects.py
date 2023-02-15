##  Objective: This file contains one function to:
# (1) print the fields included in Ecotaxa projects using the API (downloading project is thus not necessary). This should facilitate the completion of the standardizer spreadsheets
# (1bis) print the list of standard units from pint module
# (2) save interactive plot of raw ecotaxa export file to identify outliers, missing values, and units. This should facilitate the completion of the standardizer spreadsheets
# (3) flag samples to be left out of standardized projects based on 5 criteria (missing necessary data/metadata, location on land, low ROI count, high number of bubbles and other artefacts, multiple size calibration factors) and generate corresponding interactive report
# (4) perform EcoTaxa export files standardization (standard column names and units) and harmonization (across instruments) based on standardizer spreadsheets
#  Running function (4) is required after project export in order to reduce (i.e. select only the raw columns necessary for further processing), harmonize (all the PSSdb projects now have similar columns, regardless of specific instrumentation and image acquisition/processing routines), and standardize (All columns follow standardized notations and units according to international database, including taxonomic annotations) project export files

## Requirements: Change path on
# l.69: configuration file with username/password information (Ecotaxa_API_pw.yaml)
# l.93: custom units definition file (units_def.txt)
# l.106: plankton taxonomic information (plankton_annotated_taxonomy.xlsx). File created/updated below

## TO DO: Perform standardization on other imaging sensor projects (at the moment function 3 only works for Zooscan, IFCB, and UVP)


## Python modules
# Path modules
from pathlib import Path # Handling of path object

# Config modules
import wget
import yaml # requires installation of PyYAML package

# Prompt for export confirmation
import sys

# Panda and numpy modules (used for dataframes and arays)
import pandas as pd
import numpy as np
from natsort import natsorted

# Globe data modules
import geopandas as gpd
path_to_zip = Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent / 'ne_10m_admin_0_countries.zip'

# Download high resolution countries shapefile
import requests
import shutil
if path_to_zip.with_suffix('').is_dir()==False:
  with requests.Session() as sess:
    url = 'https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip'
    rsp = sess.get(url, headers={'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/102.0.0.0 Safari/537.36',}, stream=True)
    with open(path_to_zip, "wb") as fd:
        for a_chunk in rsp.iter_content():  # Loop over content, i.e. eventual HTTP chunks
            # rsp.raise_for_status()
            fd.write(a_chunk)
    shutil.unpack_archive(path_to_zip, path_to_zip.with_suffix(''))  # Unzip export file
    path_to_zip.unlink(missing_ok=True)


from shapely.geometry import Polygon, mapping, Point, MultiPolygon
from shapely import wkt
from shapely.ops import unary_union

# Plot module
import plotly.express as px # Use pip install plotly==5.8.0
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotnine import * # Python equivalent to ggplot2. Use pip install plotnine. Do not execute in Pycharm (bug- no fix yet): https://youtrack.jetbrains.com/issue/PY-42390
from colorspace import sequential_hcl # Use: pip install git+https://github.com/retostauffer/python-colorspace

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

path_to_config=Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_to_git=Path(cfg['git_dir']).expanduser()
path_to_project_list=path_to_git / cfg['proj_list']

path_to_config_usr=Path('~/GIT/PSSdb/scripts/Ecotaxa_API_pw.yaml').expanduser()
with open(path_to_config_usr ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

with ecotaxa_py_client.ApiClient() as client:
    api = authentification_api.AuthentificationApi(client)
    token = api.login(LoginReq(username=cfg_pw['ecotaxa_user'],password=cfg_pw['ecotaxa_pass']))
configuration = ecotaxa_py_client.Configuration(host = "https://ecotaxa.obs-vlfr.fr/api",access_token=token, discard_unknown_keys=True)
configuration.verify_ssl=False

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
ureg.load_definitions(Path.home()/'GIT'/'PSSdb'/'scripts' /'units_def.txt')  # This text file is used to define custom units from standard units  (e.g. square_pixel etc.)
full_list_units = list(dict.fromkeys(sorted(dir(ureg))))  # list(dict.fromkeys(sorted(list(np.concatenate([dir(getattr(ureg.sys,system)) for system in dir(ureg.sys)]).flat))))
import re

# Poisson distribution
from scipy.stats import poisson # Estimate counts uncertainties assuming detection follows a Poisson distribution
import math

# Standardization of taxonomic annotations
from funcs_standardize_annotations import taxon_color_palette,annotation_in_WORMS
import ast
from tqdm import tqdm

path_to_taxonomy=Path(Path.home()/'GIT'/'PSSdb'/'ancillary' /'plankton_annotated_taxonomy.xlsx').expanduser()
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

# Functions here:

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
def filling_standardizer_unit_func(custom_units_path='~/GIT/PSSdb/scripts/units_def.txt'):
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
    path_to_plotdir = path_to_git / 'figures' / 'Standardizer'
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
def flag_func(dataframe):
    """
        Objective: This function assign flags to image samples contained in a dataframe based on multiple criteria.
        Criteria include missing data/metadata, anomalous GPS location, low ROI counts, presence of artefacts, multiple size calibration factors
        :param dataframe: dataframe whose samples need to be quality controlled
        :return: dataframe with flags
    """
    # Flag #1: Missing required data/metadata
    dataframe['Flag_missing']=0
    index_flag=pd.isnull(dataframe).any(1).to_numpy().nonzero()[0]
    dataframe.loc[index_flag,'Flag_missing']=1
    dataframe['Missing_field'] =''
    dataframe['Missing_field'] = pd.isnull(dataframe).apply(lambda x: ';'.join(x.index[x == True].tolist()), axis=1)

    # Flag #3: low ROI count
    if 'Sample' not in dataframe.columns and 'Profile' in dataframe.columns:
        dataframe['Sample']=dataframe['Profile']
    if 'Sample' in dataframe.columns:
        # Flag #2: anomalous GPS location
        dataframe['Flag_GPScoordinatesonland'] = 0
        world = gpd.read_file(path_to_zip.with_suffix('') / 'ne_10m_admin_0_countries.shp') #gpd.datasets.get_path("naturalearth_lowres")
        world_polygon = unary_union(world.geometry.tolist())
        gdf = gpd.GeoDataFrame(dataframe[['Sample','Longitude', 'Latitude']].drop_duplicates().dropna(), geometry=gpd.points_from_xy(dataframe[['Sample','Longitude', 'Latitude']].drop_duplicates().dropna().Longitude, dataframe[['Sample','Longitude', 'Latitude']].drop_duplicates().dropna().Latitude))
        # gdf=unary_union(gdf.geometry.tolist())
        gdf['Flag_GPScoordinatesonland']=list(map(lambda x:x.within(world_polygon), gdf.geometry.tolist()))
        dataframe['Flag_GPScoordinatesonland'] = np.where(dataframe['Sample'].isin(gdf[gdf['Flag_GPScoordinatesonland'] == True]['Sample'].tolist()), 1, 0)
        dataframe['Flag_dubiousGPScoordinates'] = 0
        dataframe['Flag_dubiousGPScoordinates'] = np.where(dataframe['Sample'].isin(gdf.get('Sample')[list(map(lambda x: x.contains(Point(0, 0)), gdf.geometry.tolist()))]), 1,dataframe['Flag_dubiousGPScoordinates'])

        # Calculate count uncertainties assuming Poisson distribution. Upper-lower count limit are based on 5% uncertainty
        summary=dataframe.groupby(['Sample']).apply(lambda x :pd.Series({'ROI_count':len(x.ROI),'Count_uncertainty':poisson.pmf(k=len(x.ROI),mu=len(x.ROI)),'Count_lower':poisson.ppf((0.05/2), mu=len(x.ROI)),'Count_upper':poisson.ppf(1-(0.05/2), mu=len(x.ROI))})).reset_index()
        # kwargs={'mapping':'x:Sample,y:ROI_count,ymin:Count_lower,ymax:Count_upper','scale_y_log10()':'','theme(axis_text_x=element_blank())':''}
        # ggplot_funcs(dataframe=summary,output_path='~/GIT/PSSdb_LOV_orga/Plots/EcoTaxa_3315_count_uncertainty.png',xlab='Sample ID',ylab='Number of pictures per sample',geom='geom_pointrange',**kwargs)
        dataframe['Flag_count'] = 0
        dataframe['Flag_count'] = np.where(dataframe['Sample'].isin(summary[summary.Count_uncertainty > 0.05].Sample.tolist()), 1,dataframe['Flag_count'])
        # dataframe'0' in dataframe['Flag_count'].astype('str')].Sample.unique()
        # Flag #4: Presence of artefacts (>20% of the total sample ROIs)
        dataframe['Flag_artefacts'] = 0
        if 'Category' in dataframe.columns:
            # Extract number of artefacts per samples
            dataframe.Category = np.where(dataframe.Category == '', pd.NA, dataframe.Category)
            summary = dataframe.dropna(subset=['Category']).groupby(['Sample']).apply(lambda x: pd.Series({'Artefacts_count': len(x[ x.Category.str.lower().apply( lambda annotation: len( re.findall( r'bead|bubble|artefact|artifact', annotation)) > 0)].ROI), 'Artefacts_percentage': len(x[x.Category.str.lower().apply( lambda annotation: len( re.findall( r'bead|bubble|artefact|artifact', annotation)) > 0)].ROI) / len( x.ROI)})).reset_index()
            dataframe['Flag_artefacts'] = np.where(dataframe['Sample'].isin(summary[summary.Artefacts_percentage > 0.2].Sample.tolist()), 1, dataframe['Flag_artefacts'])
        # dataframe['0' in dataframe['Flag_artefacts'].astype('str')].Sample.unique()
        # Flag #5: Multiple size calibration factors
        if 'Pixel' in dataframe.columns:
            dataframe['Flag_size']=0
            # Extract pixel size per sample/profile
            summary=dataframe.groupby(['Pixel']).apply(lambda x: pd.Series({'Sample_percentage':len(x.Sample.unique())//len(dataframe.Sample.unique()),'Profile_percentage':len(x.Profile.unique())/len(dataframe.Profile.unique()) if 'Profile' in x.columns else 1})).reset_index()
            summary=summary.sort_values(by=['Profile_percentage'],ascending=False).reset_index()
            pixel_list=summary.loc[1:,'Pixel'].astype('str').to_list() if len(summary.Pixel.values)>1 else ['inf']
            dataframe['Flag_size'] = np.where(dataframe['Pixel'].astype('str').isin(pixel_list),1, dataframe['Flag_size'])
            # dataframe['0' in dataframe['Flag_size'].astype('str')].Sample.unique()

    dataframe['Flag']=np.where(dataframe[[column for column in dataframe.columns if 'Flag_' in column]].sum(axis=1)==0,0,1) # Default is 0, aka no anomaly
    return dataframe

#filling_standardizer_flag_func(standardizer_path='~/GIT/PSSdb/raw/project_Zooscan_standardizer.xlsx',project_id=4023,report_path='~/GIT/PSSdb_LOV/Reports')
def filling_standardizer_flag_func(standardizer_path,project_id,report_path):
    """
    Objective: This function adds flags to image samples based on multiple criteria and automatically generate a report for a given project
    :param standardizer_path: Full path of the standardizer spreadsheet containing project ID of interest
    :param project_id: ID of the project to be flagged in the standardizer spreadsheet
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
    project_path_list = list(Path(df_standardizer['Project_localpath'][project_id]).expanduser().rglob("*_"+str(project_id)+'_*')) if Path(df_standardizer['Project_localpath'][project_id]).expanduser().stem!=cfg['UVP_consolidation_subdir'] else list(Path(df_standardizer['Project_localpath'][project_id]).expanduser().glob("*_"+str(project_id)+'_*'))


    if len([path for path in project_path_list if 'flag' not in str(path)])>0:
         project_path = natsorted([path for path in project_path_list if 'flag' not in str(path)])
         df = pd.concat(map(lambda path: (columns := pd.read_table(path, nrows=0).columns, pd.read_table(path, usecols=[header for  header in df_standardizer.loc[project_id, [  column  for column in df_standardizer.columns if 'field' in column]].dropna().values if columns.isin([header]).any()]))[ -1], project_path))
         old_columns = df.columns
         df.columns = [df_standardizer.loc[ project_id, [column for column in df_standardizer.columns if 'field' in column]].dropna()[ df_standardizer.loc[project_id, [column for column in df_standardizer.columns if 'field' in column]].dropna() == column].index[0].replace( '_field', '') for column in df.columns]
         index_duplicated_fields = df_standardizer.loc[project_id, [column for column in df_standardizer.columns if 'field' in column]].dropna().index.str.replace( "_field", "").isin(list(df.columns)) == False
         if any(index_duplicated_fields):
             # Append duplicated fields
             df = pd.concat([df.reset_index(drop=True), pd.DataFrame(dict(zip(list(df_standardizer.loc[project_id,[column for column in df_standardizer.columns if 'field' in column]].dropna().index.str.replace("_field","")[index_duplicated_fields]), list( map(lambda item: df.loc[:, df.columns[old_columns.isin([item])][0]].values.tolist(),list(df_standardizer.loc[project_id,[column for column in df_standardizer.columns if 'field' in column]].dropna().values[index_duplicated_fields])))))).reset_index( drop=True)], axis=1)
         # Skip existing samples if flag file exists
         if Path(flag_overrule_path).expanduser().is_file():
             df_flags = pd.read_csv(flag_overrule_path, sep=',')
             df = df[df.Sample.isin(list(df_flags.Sample)) == False]
             df_flagged_existing = pd.merge(df[df.Sample.isin(list(df_flags.Sample))], df_flags, how='left',on='Sample')
             if len(df):
                 print('Existing flags found at {}. Updating the file with new samples'.format(flag_overrule_path))
             else:
                 print('Existing flags found at {}, but no additional samples. Skipping project'.format(flag_overrule_path))
                 return
         else:
             df_flags = pd.DataFrame()
             df_flagged_existing = pd.DataFrame()

         # Set missing values to NA
         na = str(df_standardizer.loc[project_id]['NA_value']).split(';')
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

         df.Category = np.where(df.Category.isna(), '', df.Category) # Avoid flagging of small particles without annotations for UVP projects
         df = df.astype(dict(zip(['Sample', 'Latitude', 'Longitude', 'Volume_analyzed', 'Pixel'], [str, float, float, float, float])))
         df = df.mask(df.apply(lambda x: x.astype(str).isin([convert(value, df.dtypes[x.name]) for value in pd.Series(na) if is_float(value) == False])))
         columns_to_convert = [column for column in df.columns if column not in ['Area', 'Depth_min', 'Depth_max', 'Minor_axis', 'ESD', 'Biovolume']]
         df[columns_to_convert] = df.mask(df.apply(lambda x: x.astype(str).isin([convert(value, df.dtypes[x.name]) for value in pd.Series(na)])))[columns_to_convert]
         df.Longitude = (df.Longitude + 180) % 360 - 180  # Converting all longitude to [-180,180] decimal degrees

         # Append flags
         flagged_df=flag_func(df)
         # Override flag_missing for optional variables: Sampling_size, Dilution, Annotation, Category
         flagged_df['Missing_field'] = flagged_df['Missing_field'].apply(lambda x: '' if pd.Series((field in ['Sampling_upper_size', 'Sampling_lower_size', 'Category', 'Annotation', 'Dilution'] for field in x.split(';'))).all() else x)
         flagged_df.loc[flagged_df['Missing_field'] == '', 'Flag_missing'] = 0
         flagged_df['Flag'] = np.where(flagged_df[[column for column in flagged_df.columns if 'Flag_' in column]].sum(axis=1)==0,0,1)

         # Adding URL to check flagged samples
         # Query EcoTaxa API to retrieve samples ID
         if (df_standardizer['Project_source'][project_id]=='https://ecotaxa.obs-vlfr.fr/prj/'+str(project_id)):
             with ecotaxa_py_client.ApiClient(configuration) as api_client:
                 api_instance = samples_api.SamplesApi(api_client)
                 samples = api_instance.samples_search(project_ids=int(project_id), id_pattern='')  #
             flagged_df = pd.merge(flagged_df, pd.DataFrame({'Sample': pd.Series(map(lambda x: x['orig_id'], samples)).astype(str(flagged_df.dtypes['Sample'])),'Sample_ID': list(map(lambda x: x['sampleid'], samples))}),how='left', on='Sample')
         flagged_df['Sample_URL'] = flagged_df[['Sample', 'Sample_ID']].apply(lambda x: pd.Series({'Sample_URL': r'{}?taxo=&taxochild=&ipp=100&zoom=100&sortby=&magenabled=0&popupenabled=0&statusfilter=&samples={}&sortorder=asc&dispfield=&projid={}&pageoffset=0"'.format(df_standardizer['Project_source'][project_id],x.Sample_ID,project_id)}),axis=1) if df_standardizer['Project_source'][project_id]=='https://ecotaxa.obs-vlfr.fr/prj/'+str(project_id) else flagged_df[['Sample']].apply(lambda x: pd.Series({'Sample_URL': r'{}&bin={}'.format(df_standardizer['Project_source'][project_id],x.Sample)}),axis=1) if 'ifcb' in df_standardizer['Project_source'][project_id] else ''

         # Generating flags overruling datafile that will be read to filter samples out during standardization
         path_to_datafile = Path(cfg['raw_dir']).expanduser() / 'flags' / Path(df_standardizer['Project_localpath'][project_id]).stem
         path_to_datafile.mkdir(parents=True, exist_ok=True)
         report_path = Path(report_path).expanduser() / Path(df_standardizer['Project_localpath'][project_id]).stem
         report_path.mkdir(parents=True, exist_ok=True)
         report_filename = 'Report_project_' + str(project_id) + '.html'
         path_to_report = Path(report_path).expanduser() / report_filename
         # Create readme file:
         if not Path(path_to_datafile.parent / 'README.txt').is_file():
             with open(str(Path(path_to_datafile.parent / 'README.txt')), 'w') as file:
                 file.write( "README file for project control quality check (First created on February 10, 2023):\n\nThis directory contains a summary table of project control quality flags.\nCurrent flagging is done based on 6 critera:\n\n-GPS coordinates on land (Falg_GPScoordinatesonland): takes 1 if GPS coordinates correspond to a location on land (according to a 10m spatial resolution)\n\n-Dubious GPS coordinates (Flag_dubiousGPScoordinates): takes 1 if location is at 0 degrees latitude and longitude\n\n-Missing variables (Flag_missing): takes 1 if a variable required for project standardization is missing (check required variables on the project standardizer spreadsheets)\n\n-Low particles count per sample (Flag_count): takes 1 if the total number of particles detected in a given sample yields high uncertainty (>5%) assuming counts follow a Poisson distribution\n\n-High percentage of artefacts (Flag_artefacts): takes 1 if the percentage of artefacts in a given sample is superior to 20%\n\n-Multiple pixel-to-size calibration factors (Flag_size): takes 1 if the project include multiple size calibration factors\n\nEach flag is assigned a boolean factor which takes 0 if the given sample has not been flagged or 1 if it has been flagged.\nThe overall flag (Flag) is predicted based on the result of 1 or more flagged critera.\nThe overrule boolean factor allows to overrule this flag (reset False to True) if the data owner deemed the sample should be kept for further processing.\nFlagged samples will be discarded during standardization.\nThe flag table will be updated if a project have been updated with additional samples.\nInteractive reports of project control quality check can be found at: {}\n\nContact us: nmfs.pssdb@noaa.gov".format(report_path))
         if not Path(report_path.parent / 'README.txt').is_file():
             with open(str(Path(report_path.parent / 'README.txt')), 'w') as file:
                 file.write("README file for project control quality check report (First created on February 10, 2023):\n\nThis directory contains interactive reports of project control quality flags.\nThe bottom table include flagged samples that can be checked by clicking on the URL of the first column.\nCurrent flagging is done based on 6 critera:\n\n-GPS coordinates on land (Falg_GPScoordinatesonland): takes 1 if GPS coordinates correspond to a location on land (according to a 10m spatial resolution)\n\n-Dubious GPS coordinates (Flag_dubiousGPScoordinates): takes 1 if location is at 0 degrees latitude and longitude\n\n-Missing variables (Flag_missing): takes 1 if a variable required for project standardization is missing (check required variables on the project standardizer spreadsheets)\n\n-Low particles count per sample (Flag_count): takes 1 if the total number of particles detected in a given sample yields high uncertainty (>5%) assuming counts follow a Poisson distribution\n\n-High percentage of artefacts (Flag_artefacts): takes 1 if the percentage of artefacts in a given sample is superior to 20%\n\n-Multiple pixel-to-size calibration factors (Flag_size): takes 1 if the project include multiple size calibration factors\n\nEach flag is assigned a boolean factor which takes 0 if the given sample has not been flagged or 1 if it has been flagged.\nThe overall flag (Flag) is predicted based on the result of 1 or more flagged critera.\nFlagged samples will be discarded during standardization.\nThe flag report will be updated if a project have been updated with additional samples.\nTables of project control quality check can be found at: {}\n\nContact us: nmfs.pssdb@noaa.gov".format( path_to_datafile.parent))

         #if len(flag_overrule_path) == 0 or Path(flag_overrule_path).expanduser().is_file() == False:
         overrule_name ='project_{}_flags.csv'.format(str(project_id))
         print('Saving flags to', str(path_to_datafile / overrule_name),'\nSet Overrule to True if you wish to keep samples for further processing', sep=" ")
         overruled_df = flagged_df[['Sample', 'Sample_URL', 'Flag'] + [column for column in flagged_df.columns if('Flag_' in column) or ('Missing_field' in column)]].drop_duplicates()  # flagged_df[flagged_df['Flag']==0][['Sample','Flag']].drop_duplicates()
         overruled_df['Overrule'] = False
         overruled_df = overruled_df.sort_values(by=['Flag'], ascending=False)
         pd.concat([df_flags.reset_index(drop=True), overruled_df.reset_index(drop=True)], axis=0).sort_values( by=['Flag', 'Sample'], ascending=[False, True]).to_csv(str(path_to_datafile / overrule_name), sep=',',index=False)

         # Update standardizer spreadsheet with flagged samples path and save
         df_standardizer['Flag_path'] = df_standardizer['Flag_path'].astype(str)
         df_standardizer['Flag_path'][project_id] = str(path_to_datafile / overrule_name).replace(str(Path.home()), '~')
         print('Updating standardizer spreadsheet with path of flagged samples/profiles ID datafile')
         df_standardizer['Project_ID'] = df_standardizer.index
         df_standardizer = df_standardizer[['Project_ID'] + df_standardizer.columns[:-1].tolist()]
         with pd.ExcelWriter(str(standardizer_path), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                df_standardizer.to_excel(writer, sheet_name=sheetname, index=False)
                #df_standardizer_metadata.to_excel(writer, sheet_name='Metadata', index=False)

         # Generating interactive project report
         fig = make_subplots(subplot_titles=['','','','Flag: 0 (bad flag), 1 (no flag)'],rows=3, cols=2,specs=[[{"type": "scattergeo","rowspan": 2}, {"type": "scatter",'t':0.05}],
                                              [None, {"type": "scatter",'b':0.05}],
                                              [{"type":'table',"colspan": 2},None]]
                        ,row_heights=[0.3,0.3, 0.4], vertical_spacing=0.02)
         fig.layout.annotations[0].update(xanchor='left',yanchor='top',x=0.0,y=1.0,font={'size':12})
         flagged_df = pd.concat([df_flagged_existing.reset_index(drop=True), flagged_df.reset_index(drop=True)], axis=0)

         # Map, left column
         if 'Category' not in df.columns:
             flagged_df['Category']=''
             flagged_df['Annotation']='unclassified'
         flagged_df['Annotation'] = flagged_df['Annotation'] if 'Annotation' in flagged_df.columns else 'predicted' if ( 'Annotation' not in flagged_df.columns) and ('Category' in flagged_df.columns) else 'unclassified'
         # Drop small particles to calculate the percentage of validation/artefacts for UVP projects
         flagged_df['Category'] = np.where(flagged_df.Category == '', pd.NA, df.Category)

         summary_df = flagged_df.groupby(['Sample','Sample_URL',  'Latitude', 'Longitude']+[column for column in flagged_df.columns if ('Flag' in column) or ('Missing_field' in column)], dropna=False).apply(lambda x: pd.Series({'ROI_count': x.ROI.count(),'Count_error':np.diff(poisson.ppf([0.05/2,1-(0.05/2)], mu=len(x.ROI)))[0],'Validation_percentage':len(x[x['Annotation'].isin(['validated'])].ROI) / len(x.ROI),'Artefacts_percentage':len(x[x.Category.str.lower().apply(lambda annotation:len(re.findall(r'bead|bubble|artefact|artifact',annotation))>0)].ROI) / len(x.ROI)})).reset_index()
         summary_df['Sample_URL'] = flagged_df.groupby(['Sample','Sample_URL',  'Latitude', 'Longitude']+[column for column in flagged_df.columns if ('Flag' in column) or ('Missing_field' in column)], dropna=False).apply(lambda x: pd.Series({'ROI_count': x.ROI.count(),'Count_error':np.diff(poisson.ppf([0.05/2,1-(0.05/2)], mu=len(x.ROI)))[0],'Validation_percentage':len(x.dropna(subset=['Category'])[x.dropna(subset=['Category'])['Annotation'].isin(['validated'])].ROI) / len(x.dropna(subset=['Category']).ROI),'Artefacts_percentage':len(x.dropna(subset=['Category'])[x.dropna(subset=['Category']).Category.str.lower().apply(lambda annotation:len(re.findall(r'bead|bubble|artefact|artifact',annotation))>0)].ROI) / len(x.dropna(subset=['Category']).ROI)})).reset_index()
         subset_summary_df=summary_df.dropna(subset=['Latitude', 'Longitude', 'Sample','Sample_URL'])
         if len(subset_summary_df):
             subset_summary_df['colors']=np.where(subset_summary_df.Flag_GPScoordinatesonland==1, 'red','black')
             subset_summary_df.Sample = pd.Categorical(subset_summary_df.Sample, categories=subset_summary_df.Sample.unique(),ordered=True)
             #gdf = gpd.GeoDataFrame(subset_summary_df,geometry=gpd.points_from_xy(subset_summary_df.Longitude,subset_summary_df.Latitude)).set_index('Sample').set_crs(epsg=4326, inplace=True)
             if len(subset_summary_df[subset_summary_df.Flag_GPScoordinatesonland == 1]) > 0:
                 # Definition of the zoomed (geo2) and inset (geo) maps
                 sub_data=subset_summary_df[subset_summary_df.Flag_GPScoordinatesonland == 1]
                 data_geo = dict(type='scattergeo',
                          name='GPS coordinates on land<br>(Flag_GPScoordinatesonland)',
                          lon=sub_data.Longitude,
                          lat=sub_data.Latitude,
                          hovertext='Sample ID: ' + sub_data.Sample.astype(str)+ '<br> Longitude (ºE): '+sub_data.Longitude.astype(str)+'<br> Latitude (ºN): '+sub_data.Latitude.astype(str),
                          hoverinfo="text", marker=dict(size=3, color='black'), geojson="natural earth", showlegend=True,
                          geo='geo')

                 layout = dict()
                 layout['geo2'] = dict(projection_rotation=dict(lon=np.nanmean(sub_data.Longitude), lat=np.nanmean(sub_data.Latitude), roll=0),center=dict(lon=np.nanmean(sub_data.Longitude), lat=np.nanmean(sub_data.Latitude)),projection_type='orthographic',lataxis_range=[-90,90], lonaxis_range=[-180, 180],
                 domain=dict(x=[0.0, 0.67], y=[0.4, 0.95]),bgcolor='rgba(0,0,0,0)')
                 layout['geo'] = dict(lonaxis_range=[-180,180],lataxis_range=[-80,80],
                 domain=dict(x=[0.0, 0.2], y=[0.4, 0.52]))
                 #go.Figure(data=[data_geo, data_geo], layout=layout)
                 fig.add_trace(data_geo, row=1, col=1).update_layout(go.Figure(data=[data_geo, data_geo], layout=layout).layout, margin={"r": 0, "t": 0, "l": 0, "b": 0})
                 fig.layout['geo'] = layout['geo']
                 data_geo.update({'geo': 'geo2', 'showlegend': True})  # Update geo trace before adding the zoomed
                 fig.add_trace(data_geo, row=1, col=1)
                 fig.data[len(fig.data) - 1]['geo'] = 'geo2'  # Need to update geo trace manually
                 fig.layout['geo2'] = layout['geo2']

             if len(subset_summary_df[subset_summary_df.Flag_GPScoordinatesonland==0])>0:
                 sub_data = subset_summary_df[subset_summary_df.Flag_GPScoordinatesonland == 0]
                 data_geo =dict(type='scattergeo',
                             name='No flag',
                             lon=sub_data.Longitude,
                             lat=sub_data.Latitude,
                             hovertext='Sample ID: ' + sub_data.Sample.astype(str) + '<br> Longitude (ºE): ' + sub_data.Longitude.astype(str) + '<br> Latitude (ºN): ' + sub_data.Latitude.astype(str),
                             hoverinfo="text", marker=dict(color='black', line=dict(color=sub_data.colors, width=2)), geojson="natural earth", showlegend=True,geo='geo')
                 layout = dict()
                 layout['geo2'] = dict(
                 projection_rotation=dict(lon=np.nanmean(sub_data.Longitude), lat=np.nanmean(sub_data.Latitude), roll=0),
                 center=dict(lon=np.nanmean(sub_data.Longitude), lat=np.nanmean(sub_data.Latitude)),
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
             subset_summary_df['colors'] = np.where(subset_summary_df.Flag_dubiousGPScoordinates == 1, 'orange', 'black')
             if len(subset_summary_df[subset_summary_df.Flag_dubiousGPScoordinates == 1]) > 0:
                 sub_data = subset_summary_df[subset_summary_df.Flag_dubiousGPScoordinates == 1]
                 data_geo = dict(type='scattergeo',
                             name='Dubious GPS coordinates<br>(Flag_dubiousGPScoordinates)',
                             lon=sub_data.Longitude,
                             lat=sub_data.Latitude,
                             hovertext='Sample ID: ' + sub_data.Sample.astype( str) + '<br> Longitude (ºE): ' + sub_data.Longitude.astype(str) + '<br> Latitude (ºN): ' + sub_data.Latitude.astype(str),
                             hoverinfo="text", marker=dict(color='black', line=dict(color=sub_data.colors, width=2)),
                             geojson="natural earth", showlegend=True, geo='geo')
                 layout['geo2'] = dict(projection_rotation=dict(lon=np.nanmean(sub_data.Longitude), lat=np.nanmean(sub_data.Latitude), roll=0),
                                center=dict(lon=np.nanmean(sub_data.Longitude), lat=np.nanmean(sub_data.Latitude)),
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
                             hovertext="Sample ID: " + subset_summary_df[subset_summary_df.Flag_count==1].Sample.astype(str)+'<br>ROI imaged: '+subset_summary_df[subset_summary_df.Flag_count==1].ROI_count.astype(int).astype(str), hoverinfo="text",
                             marker=dict(size=3.5,color='black', line=dict(color=subset_summary_df[subset_summary_df.Flag_count==1].colors, width=2)), mode='markers',
                             showlegend=True, visible=True), row=1, col=2)
             if len(subset_summary_df[subset_summary_df.Flag_count==0]) > 0:
                 fig.add_trace(go.Scatter(x=subset_summary_df[subset_summary_df.Flag_count==0].Sample,
                                 y=subset_summary_df[subset_summary_df.Flag_count==0].ROI_count,
                                 error_y=dict(type="data", array=subset_summary_df[subset_summary_df.Flag_count==0].Count_error,width=0,thickness=0.1,color=subset_summary_df[subset_summary_df.Flag_count==0].colors.values[0]),
                                 hovertext="Sample ID: " + subset_summary_df[subset_summary_df.Flag_count==0].Sample.astype(str)+'<br>ROI imaged: '+subset_summary_df[subset_summary_df.Flag_count==0].ROI_count.astype(int).astype(str), hoverinfo="text",
                                 marker=dict(size=4.5,color='black', line=dict(color=subset_summary_df[subset_summary_df.Flag_count==0].colors,width=2)), mode='markers',showlegend=False, visible=True), row=1, col=2)
             # Scatterplot 2, middle-right panel: Percentage of artefacts
             subset_summary_df['colors'] = np.where(subset_summary_df.Flag_artefacts == 1, 'rgba(212,85,0,0.6)','black')  # High percentage of artefacts
             fig.add_trace(go.Scatter(x=subset_summary_df.Sample,
                                  name='Percentage of validation',
                                  y=subset_summary_df.Validation_percentage,
                                  hovertext="Sample ID: " + subset_summary_df.Sample.astype(str)+'<br>% of validation: '+np.round(100*subset_summary_df.Validation_percentage,1).astype(str)+'%',
                                  hoverinfo="text", marker=dict(color='black'), mode='lines',
                                  showlegend=True, visible=True), row=2, col=2)

             if len(subset_summary_df[subset_summary_df.Flag_artefacts == 1]) > 0:
                 fig.add_trace(go.Scatter(name='High percentage of artefacts<br>(Flag_artefacts)',
                             x=subset_summary_df[subset_summary_df.Flag_artefacts == 1].Sample,
                             y=subset_summary_df[subset_summary_df.Flag_artefacts == 1].Artefacts_percentage,
                             hovertext="Sample ID: " + subset_summary_df[subset_summary_df.Flag_artefacts == 1].Sample.astype(str)+'<br>% of artefacts: '+np.round(100*subset_summary_df[subset_summary_df.Flag_artefacts == 1].Artefacts_percentage,1).astype(str)+'%',
                             hoverinfo="text",
                             marker=dict(size=3.5,color='black', line=dict(color=subset_summary_df[subset_summary_df.Flag_artefacts == 1].colors, width=2)),
                             mode='markers',
                             showlegend=True, visible=True), row=2, col=2)

             if len(subset_summary_df[subset_summary_df.Flag_artefacts == 0]) > 0:
                fig.add_trace(go.Scatter(x=subset_summary_df[subset_summary_df.Flag_artefacts == 0].Sample,
                                     y=subset_summary_df[subset_summary_df.Flag_artefacts == 0].Artefacts_percentage,
                                     hovertext="Sample ID: " + subset_summary_df[subset_summary_df.Flag_artefacts == 0].Sample.astype(str)+'<br>% of artefacts: '+np.round(100*subset_summary_df[subset_summary_df.Flag_artefacts == 0].Artefacts_percentage,1).astype(str)+'%',
                                     hoverinfo="text",
                                     marker=dict(size=4.5,color='black', line=dict(color=subset_summary_df[subset_summary_df.Flag_artefacts == 0].colors, width=2)),
                                     mode='markers', showlegend=False, visible=True), row=2, col=2)
             # Table
             summary_df['Flag_missing']=summary_df.apply( lambda x: str(x.Flag_missing) + ' (' + x.Missing_field + ')' if x.Flag_missing == 1 else x.Flag_missing,axis=1)
             fig.add_trace(go.Table(header=dict(values=['Sample/Profile ID<br>']+[column+'<br>' for column in summary_df.columns if 'Flag' in column],align=np.repeat('center',1+len([column for column in summary_df.columns if 'Flag' in column])),
                                       line_color='rgba(255,255,255,0)',fill_color='rgba(255,255,255,1)'),
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

             title = r'<a href="{}">{}</a>'.format(df_standardizer['Project_source'][project_id], 'Cruise:' + ', '.join(df["Cruise"].unique())) if all(df["Cruise"].isna()) == False else r'<a href="{}">{}</a>'.format(df_standardizer['Project_source'][project_id], 'Project ID:' + str(project_id))#'Cruise:' + str(df["Cruise"][0]) + "<br> Hover on datapoint to see sample ID" if str(df["Cruise"][0]) != 'nan' else 'Project ID:' + str(project_id) + "<br> Hover on datapoint to see sample ID"
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
             print('Saving cruise report to', path_to_report, sep=' ')
         else:
             print('No samples left after dropping GPS coordinates and sample ID NAs. Skipping project report')
    else:
        print('No project datafile found.')

# Function (4): Perform EcoTaxa export files standardization and harmonization based on standardizer spreadsheets
#standardization_func(standardizer_path='~/GIT/PSSdb/raw/project_IFCB_standardizer.xlsx',project_id=2248)
def standardization_func(standardizer_path,project_id,plot='diversity',df_taxonomy=df_taxonomy):
    """
       Objective: This function uses the instrument-specific standardizer spreadsheet to standardize variable names and units
       This should facilitate the size spectrum analysis of projects collected with UVP, IFCB, and Zooscan
       :param standardizer_path: Full path of the standardizer spreadsheet containing project ID of interest
       :param project_id: unique integer for project ID to be standardized
       :return: standardized project dataframe
       """

    # Open standardizer spreadsheet
    sheets = [sheet for sheet in pd.ExcelFile(path_to_project_list).sheet_names if 'metadata' not in sheet]
    df_projects = pd.concat(map(lambda sheet: pd.read_excel(path_to_project_list, sheet_name=sheet), sheets))
    df_projects['Project_portal'] = df_projects.Project_localpath.apply(lambda path: Path(path).expanduser().name)
    df_projects = df_projects[df_projects['PSSdb_access'] == True]
    sheet = 0 if df_projects[df_projects.Project_ID.astype(str).isin([str(project_id)])].Project_portal.unique()[ 0] != 'ecopart' else  df_projects[df_projects.Project_ID.astype(str).isin([str(project_id)])].Project_portal.unique()[ 0]  # Default is first sheet unless project is from Ecopart
    sheetname = 'Data' if (sheet == 0) and len(re.findall(r'UVP', df_projects[df_projects.Project_ID.astype(str).isin([str(project_id)])].Instrument.unique()[0])) == 0 else 'ecotaxa' if ( len(re.findall( r'UVP', df_projects[ df_projects.Project_ID.astype( str).isin( [ str(project_id)])].Instrument.unique()[ 0])) > 0) and ( sheet == 0) else sheet


    df_standardizer = pd.read_excel(standardizer_path,sheet_name=sheet,index_col=0)

    # Control. Standardization only works for IFCB, Zooscan, and UVP projects at the moment
    if project_id not in df_standardizer.index:
        print('Project ID is not included in standardizer. Quitting')
        return
    if df_standardizer.loc[project_id]['Instrument'] not in ['Zooscan','IFCB','UVP']:
        print('Standardization only applies to Zooscan, IFCB, and UVP projects. Quitting')
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

    # Retrieve flagged samples/profiles
    if str(df_standardizer.loc[project_id]['Flag_path'])!='nan':
          df_flagged=pd.read_table(df_standardizer.loc[project_id]['Flag_path'],sep=',')
          flagged_samples=df_flagged.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)')['Sample'].tolist() if len(df_flagged.query('(Flag==1 & Overrule==False) or (Flag==0 & Overrule==True)'))>0 else ['']
    else:
          flagged_samples =['']

    path_to_data = Path(df_standardizer.at[project_id, "Project_localpath"]).expanduser()
    path_files_list=list(path_to_data.rglob('**/*_{}_*'.format(str(project_id)))) if  path_to_data.stem!=cfg['UVP_consolidation_subdir'] else list( path_to_data.glob("*_"+str(project_id)+'_*'))
    path_files_list=[path for path in path_files_list if '_flag' not in str(path)]
    if len(path_files_list)>0: # Check for native format datafile
        # Load export tsv file
        df = pd.concat(map(lambda path: (columns:=pd.read_table(path,nrows=0).columns,pd.read_table(path,usecols=[header for header in fields_of_interest_series.values if columns.isin([header]).any()]).assign(File_path=path))[-1],natsorted(path_files_list)))
        old_columns = df.columns
        # Standardize column names
        df.columns = [(list(fields_of_interest_series.index)[list(fields_of_interest_series.values).index(value)]) for value in  old_columns.drop('File_path')] + ['File_path'] # pd.MultiIndex.from_arrays([[list(fields_of_interest.keys())[list(fields_of_interest.values()).index(value)] for value in df.columns],df.columns])
        index_duplicated_fields = pd.Series(fields_of_interest.keys()).isin(list(df.columns)) == False
        if any(index_duplicated_fields):
            # Append duplicated fields
            df = pd.concat([df.reset_index(drop=True), pd.DataFrame(dict( zip(list(pd.Series(fields_of_interest.keys())[index_duplicated_fields]), list(map(lambda item: df.loc[:, df.columns[ old_columns.isin( [ fields_of_interest.get( item)])][ 0]].values.tolist(), list(pd.Series( fields_of_interest.keys())[index_duplicated_fields])))))).reset_index(drop=True)], axis=1)

        # Load variables used to describe the sampling collection, method, refs, etc. (Last column)
        df_method= pd.concat(map(lambda path:(columns:=pd.read_table(path,nrows=0).columns,pd.read_table(path,usecols=[header for header in [fields_of_interest_series['Sample']]+fields_for_description.values.tolist() if columns.isin([header]).any()]).rename(columns=dict(zip([header for header in [fields_of_interest_series['Sample']]+fields_for_description.values.tolist() if columns.isin([header]).any()],[(['Sample']+fields_for_description.index.tolist())[index] for index,header in enumerate(pd.Series([fields_of_interest_series['Sample']]+fields_for_description.values.tolist())) if columns.isin([header]).any()]))))[-1],natsorted(path_files_list)))
        df_method[fields_for_description.index] = df_method[fields_for_description.index].apply(lambda x: PA(x, dtype=ureg[dict(units_for_description)[x.name]].units) if x.name in dict(units_for_description).keys() else x)
        df_method[fields_for_description.index]=df_method[fields_for_description.index].apply(lambda x:x.name+":"+x.astype(str))

        # Check for existing standardized file(s):
        path_to_standard_dir = path_to_git / cfg['standardized_raw_subdir'] / path_to_data.stem / path_files_list[0].parent.stem
        path_to_standard_dir.mkdir(parents=True, exist_ok=True)
        path_to_standard_file = list(path_to_standard_dir.rglob('standardized_project_{}*.csv'.format(str(project_id))))
        if len(path_to_standard_file):
            df_standardized_existing = pd.concat(map(lambda path: (columns := pd.read_table(path, nrows=0).columns, pd.read_table(path, sep=","))[-1],natsorted(path_to_standard_file)))
            remaining_standard_samples = df[df.Sample.isin(list(df_standardized_existing.Sample.unique())) == False].Sample.unique()
            if len(remaining_standard_samples):
                print('Existing standardized file(s) found at {}. Generating standardized file(s) for new samples'.format(path_to_standard_dir))
                df=df[df.Sample.isin(list(remaining_standard_samples))].reset_index(drop=True)
                df_method = df_method[df_method.Sample.isin(list(remaining_standard_samples))].reset_index(drop=True)
            else:
                print('Existing standardized file(s) found at {}, but no additional samples. Skipping project'.format(path_to_standard_dir))
                return
        else:
            df_standardized_existing = pd.DataFrame()

        # Remove flagged samples/profiles
        df=df[df['Sample'].isin(flagged_samples)==False].reset_index(drop=True)
        df_method=df_method[df_method['Sample'].isin(flagged_samples)==False].reset_index(drop=True)
        if len(df)==0:
            print('No samples left after flagging. Skipping project ', project_id, sep='')
            return

        # Append method description
        additional_description =' '.join([':'.join([additional_key,columns_for_description[additional_key]]) for additional_key in columns_for_description.keys() if additional_key not in fields_for_description.index]) if type(columns_for_description)==dict else []
        if (len(fields_for_description.index)>0) or (len(additional_description) > 0):
            df=df.assign(Sampling_description=additional_description+' ' + df_method[fields_for_description.index].apply(' '.join, axis=1) if len(additional_description) > 0 else df_method[fields_for_description.index].apply(' '.join, axis=1) )
        else:
            df=df.assign(Sampling_description=columns_for_description)

        # Set missing values to NA
        na = str(df_standardizer.loc[project_id]['NA_value']).split(';')
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

        # Convert known data types to ensure masking of missing values is correct
        df=df.astype(dict(zip(['Sample','Latitude','Longitude','Volume_analyzed','Pixel'],[str,float,float,float,float])))
        # Convert NAs
        df = df.mask(df.apply(lambda x: x.astype(str).isin([convert(value, df.dtypes[x.name]) for value in pd.Series(na) if is_float(value) == False])))
        columns_to_convert = [column for column in df.columns if column not in ['Area', 'Depth_min', 'Depth_max', 'Minor_axis', 'ESD', 'Biovolume']]
        df[columns_to_convert] =  df.mask(df.apply(lambda x: x.astype(str).isin([convert(value, df.dtypes[x.name]) for value in pd.Series(na)])))[columns_to_convert]

        # Transform variables based on field units
        units=[string for string in list(df_standardizer.columns) if "_unit" in string]
        units_of_interest =  dict(zip(units,df_standardizer.loc[project_id][units].values.flatten().tolist()))
        # Test whether all units are defined:                             [re.sub(r'^-?[0-9]_times_10power-?[0-9]_','',str(unit)) for unit in list(units_of_interest.values())]
        if any([unit not in ureg for unit in list(units_of_interest.values()) if str(unit)!='nan']):
            print('Quitting. Please fill out standardizer with units from the list below or define your custom unit in {} using known units:\n'.format(path_to_data.parent.parent.parent /'Scripts'/'units_def.txt'),full_list_units, sep='')
            return

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

        if len(df['Sampling_upper_size'].loc[np.isnan(df['Sampling_upper_size'])])>0:
            df['Sampling_upper_size'].loc[np.isnan(df['Sampling_upper_size'])] = (camera_resolution[instrument]['Sampling_upper_size']/pixel_size_ratio.magnitude[0])*ureg('millimeter').to(units_of_interest['Sampling_upper_size_unit'])
        if len(df['Sampling_lower_size'].loc[np.isnan(df['Sampling_lower_size'])])>0:
            df['Sampling_lower_size'].loc[np.isnan(df['Sampling_lower_size'])] = (camera_resolution[instrument]['Sampling_lower_size'] / (pixel_size_ratio.magnitude[0]**2)) * ureg('millimeter').to(units_of_interest['Sampling_lower_size_unit']) if instrument !='Zooscan' else (2*(((camera_resolution[instrument]['Sampling_lower_size'] / (pixel_size_ratio.magnitude[0]**2))/np.pi)**0.5)) * ureg('millimeter').to(units_of_interest['Sampling_lower_size_unit'])

        # Use pint units system to convert units in standardizer spreadsheet to standard units
        # (degree for latitude/longitude, meter for depth, multiple of micrometer for plankton size)

        df.datetime = pd.to_datetime(df['Sampling_date'].astype(str) + ' ' + df['Sampling_time'].astype(str).str.zfill(6), format=' '.join(df_standardizer.loc[project_id][['Sampling_date_format', 'Sampling_time_format']]),utc=True) if all(pd.Series(['Sampling_date', 'Sampling_time']).isin(df.columns)) else pd.to_datetime(df['Sampling_date'].astype(str), format=df_standardizer.at[project_id, 'Sampling_date_format'],utc=True) if 'Sampling_date' in df.columns else pd.to_datetime(df['Sampling_time'].astype(str).str.zfill(6),format=df_standardizer.at[ project_id, 'Sampling_time_format'],utc=True) if 'Sampling_time' in df.columns else pd.NaT
        df.Longitude = (df.Longitude + 180) % 360 - 180  # Converting all longitude to [-180,180] decimal degrees

        # Update taxonomic annotations standards (df_taxonomy) with new annotations
        if 'Category' in df.columns:
            new_categories = df.dropna(subset=['Category']).Category[df.dropna(subset=['Category']).Category.isin(list(df_taxonomy.Category)) == False].unique()
            if len(new_categories):
                df_taxonomy_new = pd.concat( map(lambda hierarchy: annotation_in_WORMS(hierarchy.replace("_", " ")).assign(Category=hierarchy), new_categories))
                df_taxonomy_new['URL'] = df_taxonomy_new.WORMS_ID.apply( lambda id: 'https://www.marinespecies.org/aphia.php?p=taxdetails&id={}'.format( id.replace('urn:lsid:marinespecies.org:taxname:', '')) if len(id) else '')
                df_taxonomy_new['Life_stage'] = df_taxonomy_new.Functional_group.apply(lambda group: ';'.join([ast.literal_eval(dict)['Life stage'] for dict in group.split(';') if len(group) > 0 and 'Life stage' in ast.literal_eval(dict).keys()]))
                df_taxonomy_new['functional_group'] = df_taxonomy_new.Functional_group.apply(lambda group: ';'.join([ dict.replace('{', '').replace(', ', ' (').replace( '}', ')').replace( "'", "") if len( group) > 0 and len( ast.literal_eval( dict)) > 1 else dict.replace( '{', '').replace( ', ', ' (').replace( '}', '').replace( "'", "") if len( group) > 0 and len( ast.literal_eval(dict)) == 1 else ''  for dict  in group.split( ';')]))
                df_taxonomy = pd.concat([df_taxonomy.reset_index(drop=True), df_taxonomy_new[df_taxonomy.columns].reset_index(drop=True)],axis=0, ignore_index=True)
                df_taxonomy = df_taxonomy.sort_values(['Type', 'Category'], ascending=[False, True]).reset_index( drop=True)
                with pd.ExcelWriter(str(path_to_taxonomy), engine="openpyxl", mode="a",if_sheet_exists="replace") as writer:
                    df_taxonomy.to_excel(writer, sheet_name='Data', index=False)

        df_standardized = df.assign(Instrument=df_standardizer.loc[project_id,'Instrument'],Project_ID=project_id,
                       Station=df.Station if 'Station' in df.columns else pd.NA,
                       Profile=df.Profile if 'Profile' in df.columns else pd.NA,
                       Sample=df.Sample if 'Sample' in df.columns else pd.NA,
                       Sampling_date=df.datetime.dt.strftime('%Y%m%d') if all(np.isnan(df.datetime)==False) else pd.NA,
                       Sampling_time=df.datetime.dt.strftime('%H%M%S') if all(np.isnan(df.datetime)==False) else pd.NA,
                       Pixel=pixel_size_ratio if 'Pixel' in df.columns else pd.NA, # in pixels per millimeter
                       Volume_analyzed=1000*(list(df.Volume_analyzed)*ureg(units_of_interest['Volume_analyzed_unit']).to(ureg.cubic_meter))if 'Volume_analyzed' in df.columns else pd.NA, # cubic decimeter
                       Latitude=list(df.Latitude)*ureg(units_of_interest['Latitude_unit']).to(ureg.degree) if 'Latitude' in df.columns else pd.NA, # degree decimal
                       Longitude=list(df.Longitude) * ureg(units_of_interest['Longitude_unit']).to(ureg.degree) if 'Longitude' in df.columns else pd.NA,  # degree decimal
                       Depth_min=list(df.Depth_min) * ureg(units_of_interest['Depth_min_unit']).to(ureg.meter) if 'Depth_min' in df.columns else pd.NA,  # meter
                       Depth_max=list(df.Depth_max) * ureg(units_of_interest['Depth_max_unit']).to(ureg.meter) if 'Depth_max' in df.columns else pd.NA,  # meter
                       Annotation=np.where(df.Annotation.apply(lambda x:'' if str(x)=='nan' else x).astype(str).str.len()>0,df.Annotation,'unclassified') if 'Annotation' in df.columns else 'predicted' if ('Annotation' not in df.columns) and ('Category' in df.columns) else 'unclassified',
                       Category=df.Category.astype(str).apply(lambda x:'' if str(x)=='nan' else x) if 'Category' in df.columns else '',
                       Area=1e06*((df.Area*[1*ureg(units_of_interest['Area_unit']).to(ureg.square_pixel).magnitude for ntimes in range(df[['Area']].shape[1])])/((np.vstack([pixel_size_ratio**2 for ntimes in range(df[['Area']].shape[1])]).T)[0,:])) if all(pd.Series(['Area','Pixel']).isin(df.columns)) else pd.NA, # square micrometers
                       Biovolume=1e09*(df.Biovolume*[1*ureg(units_of_interest['Biovolume_unit']).to(ureg.cubic_pixel).magnitude for ntimes in range(df[['Biovolume']].shape[1])]/((np.vstack([pixel_size_ratio**3 for ntimes in range(df[['Biovolume']].shape[1])]).T)[0,:])) if all(pd.Series(['Biovolume','Pixel']).isin(df.columns)) else pd.NA, # cubic micrometers
                       Minor_axis=1e03*(df.Minor_axis*[1*ureg(units_of_interest['Minor_axis_unit']).to(ureg.pixel).magnitude for ntimes in range(df[['Minor_axis']].shape[1])]/((np.vstack([pixel_size_ratio for ntimes in range(df[['Minor_axis']].shape[1])]).T)[0,:])) if all(pd.Series(['Minor_axis','Pixel']).isin(df.columns)) else pd.NA, # cubic micrometers
                       ESD=1e03*(df.ESD*[1*ureg(units_of_interest['ESD_unit']).to(ureg.pixel).magnitude for ntimes in range(df[['ESD']].shape[1])]/((np.vstack([pixel_size_ratio for ntimes in range(df[['ESD']].shape[1])]).T)[0,:])) if all(pd.Series(['ESD','Pixel']).isin(df.columns)) else pd.NA, # micrometers
                       Sampling_type=df.Sampling_type if 'Sampling_type' in df.columns else pd.NA,
                       Sampling_lower_size=list(df.Sampling_lower_size)*ureg(units_of_interest['Sampling_lower_size_unit']).to(ureg.micrometer) if 'Sampling_lower_size' in df.columns else pd.NA, # micrometer
                       Sampling_upper_size=list(df.Sampling_upper_size)*ureg(units_of_interest['Sampling_upper_size_unit']).to(ureg.micrometer) if 'Sampling_upper_size' in df.columns else pd.NA # micrometer
                       )


        # Set volume analyzed to 5 mL for IFCB projects
        if df_standardizer.loc[project_id]['Instrument'] in ['IFCB'] and all(df_standardized[['Volume_analyzed']].isna()):
            df_standardized['Volume_analyzed']=PQ(5*ureg('milliliter').to('liter')).magnitude

        # Convert volume analyzed to volume imaged to account for samples dilution or fractionation
        if 'Dilution' in df.columns:
            dilution_factor = df.Dilution#np.where(df.Dilution > 1, df.Dilution, 1 / df.Dilution)
            dilution_factor[pd.Series(dilution_factor).isna()] = 1  # Replace NA values with 1
        else:
            dilution_factor=1

        df_standardized = df_standardized.assign(Volume_imaged=df_standardized.Volume_analyzed/dilution_factor) # cubic decimeters
        # Generate metadata sheet
        df_standardized=df_standardized[['Project_ID','Cruise','Instrument','Sampling_type', 'Station', 'Profile','Sample', 'Latitude', 'Longitude', 'Sampling_date', 'Sampling_time','Depth_min', 'Depth_max', 'Volume_analyzed', 'Volume_imaged',
                                    'ROI', 'Annotation','Category', 'Minor_axis', 'ESD', 'Area', 'Biovolume','Pixel','Sampling_lower_size','Sampling_upper_size','Sampling_description']]
        df_standardized_metadata=pd.DataFrame({'Variables':df_standardized.columns,'Variable_types':df_standardized.dtypes,
        'Units/Values/Timezone':['','','','','','','','degree','degree','yyyymmdd (UTC)','hhmmss (UTC)','meter','meter','cubic_decimeter','cubic_decimeter','','','']+['micrometer' for ntimes in range(df_standardized[['Minor_axis']].shape[1])]+['micrometer' for ntimes in range(df_standardized[['ESD']].shape[1])]+['square_micrometer' for ntimes in range(df_standardized[['Area']].shape[1])]+['cubic_micrometer' for ntimes in range(df_standardized[['Biovolume']].shape[1])]+['pixel_per_millimeter','micrometer','micrometer',''],
        'Description':['Project ID','Project cruise','Instrument','Sampling type (e.g. platform, gear, strategy)','Station ID (native format)','Profile ID (native format)','Sample ID (native format)','Latitude','Longitude','Sampling date','Sampling time','Minimum sampling depth','Maximum sampling depth','Volume analyzed (not accounting for sample dilution and/or fractionation)','Volume imaged (accounting for sample dilution and/or fractionation)','Region of interest ID (native format)','Status of ROI annotation (e.g. unclassified, predicted, validated)','ROI assigned taxonomy']+ ['Object minor ellipsoidal axis derived from '+ field if field!='' or len(fields['Minor_axis'])>1 else 'Object minor ellipsoidal axis' for field in fields['Minor_axis'] ]+ ['Object equivalent spherical diameter derived from '+ field if field!='' or len(fields['ESD'])>1 else 'Object equivalent spherical diameter' for field in fields['ESD'] ]+ ['Object surface area derived from '+ field if field!='' or len(fields['Area'])>1 else 'Object surface area' for field in fields['Area'] ]+['Object biovolume derived from '+ field if field!='' or len(fields['Biovolume'])>1 else 'Object biovolume' for field in fields['Biovolume']]+['Pixel size used for size conversion','Smallest sampled size','Largest sampled size','Additional description of the sampling method or protocol']})
        # Create readme file:
        if not Path(path_to_standard_dir.parent.parent / 'README.txt').is_file():
            with open(str(Path(path_to_standard_dir.parent.parent / 'README.txt')), 'w') as file:
                file.write( "README file for project standardized files (First created on February 10, 2023):\n\nThis directory contains the standardized table(s) of accessible projects.\nEach table include the following variables:\n\n{}\n\nContact us: nmfs.pssdb@noaa.gov".format('\n'.join(list(df_standardized_metadata[['Variables', 'Units/Values/Timezone', 'Description']].apply(lambda x: str(x.Variables) + " (" + str(x['Units/Values/Timezone']).strip() + "): " + x.Description if len(x['Units/Values/Timezone']) else str(x.Variables) + ": " + x.Description, axis=1).values))))

        # Save standardized dataframe
        path_to_standard_plot = path_to_git / 'figures' / 'standardizer' /path_to_data.stem /  path_files_list[0].parent.stem/ 'standardized_project_{}.html'.format(str(project_id))
        path_to_standard_plot.parent.mkdir(parents=True, exist_ok=True)
        print('Saving standardized datafile to', path_to_standard_dir,sep=' ')
        path_dict = df[['File_path', 'Sample']].drop_duplicates().groupby(['File_path'])['Sample'].apply(list).to_dict()
        df.groupby(['File_path']).apply(lambda x :df_standardized[df_standardized.Sample.isin(path_dict.get(x.File_path.unique()[0]))].to_csv(path_to_standard_dir/'standardized_project_{}_{}.csv'.format(project_id,str(x.File_path.unique()[0].stem)[str(x.File_path.unique()[0].stem).rfind('_'+str(project_id)+'_')+2+len(str(project_id)):len(str(x.File_path.unique()[0].stem))].replace('_features','')), sep=",",index=False))
        #with pd.ExcelWriter(str(path_to_standard_file),engine="xlsxwriter") as writer:
        #df_standardized.to_csv(path_to_standard_file, sep="\t",index=False) # to_excel(writer, sheet_name='Data', index=False)
        #df_standardized_metadata.to_csv(path_to_metadata_standard, sep="\t",index=False)  # to_excel(writer, sheet_name='Metadata',index=False)

        # Interactive plots
        df_standardized=pd.concat([df_standardized_existing.reset_index(drop=True),df_standardized.reset_index(drop=True)],axis=0,ignore_index=True)
        if plot=='diversity':
            subplot_specs = [[{"type": "scattergeo", "rowspan": 2}, {"type": "scatter"}], [None, {"type": "pie"}]]
            subplot_titles=['Map of the projet samples:','','Annotations diversity']
        fig = make_subplots(rows=2, cols=2,specs=subplot_specs,subplot_titles=subplot_titles,
                            column_widths=[0.5, 0.5],row_heights=[0.3,0.7], vertical_spacing=0.1)
        # Map, left panel
        df_standardized['Profile']= df_standardized['Profile'].astype(str) # Required to pass groupby if missing
        df_standardized['Station'] = df_standardized['Station'].astype(str)# Required to pass groupby if missing
        df_standardized['Sample'] = df_standardized['Sample'].astype(str)  # Required to pass groupby if missing
        df_standardized = df_standardized[(df_standardized.Longitude <= 180) & (df_standardized.Longitude >= -180) & ( df_standardized.Latitude <= 90) & (df_standardized.Latitude >= -90)]

        summary_df_standardized=df_standardized.groupby( ['Sample', 'Station', 'Latitude', 'Longitude', 'Profile', 'Volume_imaged','Depth_min','Depth_max'],dropna=True).apply(lambda x: pd.Series({'Count':x.ROI.count(),'Abundance': x.ROI.count() / ( x.Volume_analyzed.unique()[0]), # individuals per liter
                                                                                                                                                                       'Average_diameter':np.nanmean(2*np.power(x.Area/np.pi,0.5)) if 'Area' in df_standardized.columns else np.nanmean(x.ESD), # micrometer
                                                                                                                                                                       'Std_diameter':np.nanstd(2*np.power(x.Area/np.pi,0.5)) if 'Area' in df_standardized.columns else np.nanstd(x.ESD)})).reset_index()
        summary_df_standardized=summary_df_standardized.groupby(['Sample', 'Station', 'Latitude', 'Longitude', 'Profile', 'Volume_imaged'],dropna=True).apply(lambda x: pd.Series({'Abundance': np.nanmean(x.Abundance), # individuals per liter
                                                                                                                                     'Average_diameter':np.nanmean(x.Average_diameter), # micrometer
                                                                                                                                     'Std_diameter':np.nanmean(x.Std_diameter)})).reset_index() # micrometer

        if len(summary_df_standardized):
            data_geo = dict(type='scattergeo',
                        name='',
                        lon=summary_df_standardized.dropna(subset=['Latitude', 'Longitude', 'Sample']).Longitude,
                        lat=summary_df_standardized.dropna(subset=['Latitude', 'Longitude', 'Sample']).Latitude,
                        # size=summary_df_standardized.Abundance,
                        hovertext='Sample ID: ' + summary_df_standardized.dropna(
                            subset=['Latitude', 'Longitude', 'Sample']).Sample, hoverinfo="text",
                        marker=dict(color='black',size=0.5), geojson="natural earth", showlegend=False, geo='geo')
            layout = dict()
            layout['geo2'] = dict(
            projection_rotation=dict(lon=np.nanmean(summary_df_standardized.dropna(subset=['Latitude', 'Longitude', 'Sample']).Longitude), lat=np.nanmean(summary_df_standardized.dropna(subset=['Latitude', 'Longitude', 'Sample']).Latitude), roll=0),
            center=dict(lon=np.nanmean(summary_df_standardized.dropna(subset=['Latitude', 'Longitude', 'Sample']).Longitude), lat=np.nanmean(summary_df_standardized.dropna(subset=['Latitude', 'Longitude', 'Sample']).Latitude)),
            projection_type='orthographic',lataxis_range=[-90,90], lonaxis_range=[-180, 180],
            domain=dict(x=[0.0, 0.47], y=[0.0, 0.9]), bgcolor='rgba(0,0,0,0)')
            layout['geo'] = dict(lonaxis_range=[-180, 180], lataxis_range=[-80, 80], domain=dict(x=[0.0, 0.1], y=[0.0, 0.12]))

            # go.Figure(data=[data_geo, data_geo], layout=layout)
            fig.add_trace(data_geo, row=1, col=1).update_layout(go.Figure(data=[data_geo, data_geo], layout=layout).layout)
            data_geo.update({'geo': 'geo2', 'showlegend': False,'marker':dict(color='black',size=2.5)})  # Update geo trace before adding the zoomed
            fig.add_trace(data_geo, row=1, col=1)
            fig.data[1]['geo'] = 'geo2'  # Need to update geo trace manually

            #fig.add_trace(go.Scattergeo(name='',lon=summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Longitude, lat=summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Latitude, hovertext='Sample ID: ' + summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Sample, hoverinfo="text", marker=dict(color='black'), geojson="natural earth", showlegend=False),row=1, col=1).update_geos(projection_rotation=dict(lon=np.nanmean(summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Longitude), lat=np.nanmean(summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Latitude), roll=0),center=dict(lon=np.nanmean(summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Longitude), lat=np.nanmean(summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample']).Latitude)),projection_type='orthographic')
            # Scatterplot 1, top-right panel
            cols=summary_df_standardized.columns
            fig.add_trace(go.Scatter(name='',x=summary_df_standardized.dropna(subset=['Sample','Abundance']).Sample, y=summary_df_standardized.dropna(subset=['Sample','Abundance']).Abundance,
                               #  error_y=dict(type="data",array=summary_df_standardized.Std_diameter),
                                 hovertext="Sample ID: " + summary_df_standardized.dropna(subset=['Sample','Abundance']).Sample, hoverinfo="text",
                                 marker=dict(color='black'),mode='markers',
                                 showlegend=True, visible=True), row=1, col=2)
            # Scatterplot 2, bottom-right panel
            df_standardized=pd.merge(df_standardized,pd.DataFrame([2 * ((df_standardized[['Area']] / np.pi) ** 0.5).to_numpy() if 'Area' in df_standardized.columns else pd.NA][0],columns=['ECD' for ntimes in range(df_standardized[['Area']].shape[1])]),left_index=True,right_index=True )# in micrometer
            if plot=='diversity':
                subset_df = pd.merge(df_standardized, df_taxonomy[['Category', 'EcoTaxa_hierarchy', 'Full_hierarchy', 'Rank', 'Type','Domain']+[rank for rank in rank_categories if rank in df_taxonomy.columns]], left_on='Category', right_on='Category', how='left')
                columns_ID = ['Project_ID']  # ['sample_profileid','sample_stationid','object_lat','object_lon','Depth_bin_med']
                column_group = ['Genus']

                df_taxonomy[[column for column in df_taxonomy.columns if column in rank_categories and column in np.array(rank_categories)[np.where( pd.Categorical(rank_categories, categories=rank_categories, ordered=True) <= column_group[0])[0].tolist()]] + ['Functional_group']].drop_duplicates().dropna().groupby(column_group[0]).apply(lambda x: pd.Series({'group': x.Functional_group.value_counts().index[0]})).reset_index().sort_values(by='group')
                group_categories = df_taxonomy[df_taxonomy.EcoTaxa_hierarchy.isin(subset_df.EcoTaxa_hierarchy.tolist())][[column for column in df_taxonomy.columns if column in rank_categories and column in np.array(rank_categories)[np.where(pd.Categorical(rank_categories, categories=rank_categories, ordered=True) <= column_group[0])[ 0].tolist()]]]
                group_categories =group_categories.sort_values(by=column_group[0])
                group_categories[group_categories.isnull()] = ''
                if len(subset_df.dropna(subset=['Full_hierarchy'])) > 0:
                    summary_df_taxonomy =  subset_df.dropna(subset=['Full_hierarchy']).groupby(by=columns_ID + ['Phylum'], observed=True).apply( lambda x: pd.Series({'ROI_count': x.ROI.count()})).reset_index().pivot(index=columns_ID,columns='Phylum', values='ROI_count').reset_index()
                    melt_df_taxonomy = summary_df_taxonomy.melt(
                    value_vars=[column for column in summary_df_taxonomy.columns if column not in columns_ID],
                    id_vars=columns_ID, var_name='Group', value_name='Count').dropna()
                    melt_df_taxonomy = pd.merge(melt_df_taxonomy, group_categories, left_on='Group', right_on='Phylum',how='left')
                    melt_df_taxonomy.Group = pd.Categorical(melt_df_taxonomy.Group,categories=group_categories['Phylum'].unique(), ordered=True)
                    melt_df_taxonomy.Group.cat.categories
                    melt_df_taxonomy = melt_df_taxonomy.sort_values(by=['Group'] + columns_ID)
                    levels = pd.Categorical(melt_df_taxonomy.columns[melt_df_taxonomy.columns.isin(rank_categories)].tolist(),categories=pd.Series(rank_categories)[ pd.Series(rank_categories).isin(melt_df_taxonomy.columns.tolist())].tolist(), ordered=True)
                    colors_dict = taxon_color_palette(dataframe=group_categories, levels=levels,palette=px.colors.qualitative.T10)
                    levels_dict = dict(zip(levels, px.colors.sequential.Reds[0:len(levels)][::-1]))
                    # fig=px.pie(melt_df_taxonomy,values='Count',names='Group',color='Group',hole=0.3,color_discrete_map=colors_dict).update_traces(sort=False,textinfo='percent+label').update_layout(title_text="Annotation diversity",annotations=[dict(text='% ROI', x=0.5, y=0.5, font_size=20, showarrow=False)])
                    level_new_data = melt_df_taxonomy[columns_ID + ['Group', 'Count']].drop_duplicates().groupby(by=columns_ID).apply(lambda x: pd.Series({'Total': x.Count.sum(), 'Rank': 'Phylum'})).reset_index()
                    fig.add_trace(go.Pie(name='Phylum', labels=level_new_data['Rank'], values=level_new_data['Total'], hoverinfo='none',
                       textinfo='none', direction='clockwise', hole=0.2, legendgroup=2, showlegend=True, sort=False,
                       marker=dict(colors=['rgba(255,255,255,0)'], line=dict(color=levels_dict['Phylum'], width=6)),
                       textposition='inside', insidetextorientation='radial'), row=2, col=2).update_layout(annotations=[dict(text='% ROI', x=0.78,y=0.31, font_size=12, showarrow=False),dict(text='Pie chart of the annotations diversity '+str((np.round(100*subset_df['Type'].value_counts(normalize=True),1).astype(str)+ '%').to_dict()), x=0.7, y=-0.05, font_size=12, showarrow=False)])

                    fig.add_trace( go.Pie(name='Taxonomic phylum:', labels=melt_df_taxonomy[['Group', 'Count']].drop_duplicates()['Group'],
                       values=melt_df_taxonomy[['Group', 'Count']].drop_duplicates()['Count'], hole=0.2,
                       direction='clockwise', legendgroup=1, sort=False, hoverinfo='percent+label',
                       textinfo='percent+label', textfont_size=10, textposition='none',
                       marker=dict(colors=[color for taxon, color in colors_dict.items() if taxon in melt_df_taxonomy[['Group', 'Count']].drop_duplicates()['Group'].values.tolist()])), row=2, col=2)
                    levels = df_taxonomy.columns[df_taxonomy.columns.isin(rank_categories)].tolist()
                    for index, level in enumerate(levels):

                        if level == 'Phylum':
                            continue
                        else:
                            summary_df_taxonomy = subset_df.dropna(subset=['Full_hierarchy']).groupby(by=columns_ID + levels[0:index + 1], observed=True, dropna=False).apply(lambda x: pd.Series({'ROI_count': x.ROI.count()})).reset_index()
                            summary_df_taxonomy = summary_df_taxonomy[summary_df_taxonomy['Phylum'].isnull() == False]
                            summary_df_taxonomy.loc[:, levels[0:index + 1]] = summary_df_taxonomy.loc[:,levels[0:index + 1]].fillna('unassigned')
                            summary_df_taxonomy = pd.merge(summary_df_taxonomy,summary_df_taxonomy.groupby(by=columns_ID + [levels[index - 1]], observed=True, group_keys=False).apply( lambda x: pd.Series({'Total_ROI_count': x.ROI_count.sum()})).reset_index()[ [levels[index - 1], 'Total_ROI_count']], on=levels[index - 1],how='left')
                            summary_df_taxonomy['Phylum'] = pd.Categorical(summary_df_taxonomy['Phylum'],categories=group_categories['Phylum'].unique(),ordered=True)
                            new_data = summary_df_taxonomy  # melt_df_taxonomy.groupby(level).apply(lambda x :pd.Series({'Count': sum(x.Count)})).reset_index()
                            new_data = new_data.sort_values(by=levels[0:index])
                            new_colors = [colors_dict[taxon] if taxon != 'unassigned' else '#FFFFFF' for taxon in new_data[new_data.columns.tolist()[index + 1]]]  # [color for taxon,color in colors_dict.items() if taxon in new_data[new_data.columns.tolist()[2]].values.tolist()]
                            new_data['Hierarchy'] = new_data.apply(lambda x: '>'.join(x[levels[0:index + 1]][levels[0:index + 1]] + ' (' + levels[0:index + 1] + ')'), axis=1)
                            level_new_data = new_data.groupby(by=columns_ID).apply( lambda x: pd.Series({'Total': x.ROI_count.sum(), 'Rank': level})).reset_index()
                            fig.add_trace(go.Pie(name=level, labels=level_new_data['Rank'], values=level_new_data['Total'],
                               hoverinfo='none', textinfo='none', direction='clockwise',
                               hole=0.2 + index * (1 / (len(levels) + 1)), legendgroup=2, showlegend=True, sort=False,
                               marker=dict(colors=['#FFFFFF'], line=dict(color=levels_dict[level], width=6)),
                               textposition='inside', insidetextorientation='radial'), row=2, col=2)
                            fig.add_trace(go.Pie(name=level, labels=new_data['Hierarchy'], values=new_data['ROI_count'],
                                               hoverinfo='label+percent', textinfo='none', direction='clockwise',
                                               hole=0.2 + index * (1 / (len(levels) + 1)), showlegend=False, sort=False,
                                               marker=dict(colors=new_colors, line=dict(color='#000000', width=0)),
                                               textposition='inside', insidetextorientation='radial'), row=2, col=2)


            cols_labels=np.where(cols.isin(['Abundance']),'Abundance (particles dm<sup>-3</sup>)','')
            cols_labels = np.where(cols.isin(['Average_diameter']), u'Average equivalent circular diameter (\u03BCm)', cols_labels)
            cols_labels = np.where(cols.isin(['Std_diameter']), u'Std equivalent circular diameter (\u03BCm)', cols_labels)

            button_scatter1 = [dict(method="restyle",  # argument to change data
                                args=[{'x': [summary_df_standardized['Sample'], 'undefined'],
                                       'y': [summary_df_standardized[cols[k]], 'undefined'],
                                       'visible': [True, True, True]}, [2]],
                                # visible should have the same length as number of subplots, second argument specifies which subplot is updated
                                label=cols_labels[k]) for k in np.where(np.isin(list(cols), ['Abundance', 'Average_diameter', 'Std_diameter']))[0]]

            title =  r'<a href="{}">{}</a>'.format(df_standardizer['Project_source'][project_id], 'Cruise:' + ', '.join(df_standardized["Cruise"].unique())) if all(df_standardized["Cruise"].isna())==False else r'<a href="{}">{}</a>'.format(df_standardizer['Project_source'][project_id], 'Project ID:' + str(project_id))
            fig.update_layout(updatemenus=list([dict(active=0, buttons=button_scatter1, x=0.87, y=1.1, xanchor='left', yanchor='top')]),
            title={'text': title, 'xanchor': 'center', 'yanchor': 'top', 'x': 0.5},
            xaxis={'title': 'Sample ID', 'nticks': 2, 'tickangle': 0, 'tickfont': dict(size=10),'titlefont': dict(size=12)},
            yaxis={"tickmode": "array",'title': 'Variable selected in <br> dropdown menu ', 'tickfont': dict(size=10),'titlefont': dict(size=12)},
            boxmode='group')
            fig.update_yaxes(type="log", row=1, col=2)
            fig.update_yaxes(type="log", row=2, col=2)
            fig.update_xaxes(type="log", row=2, col=2,ticks="inside")
            fig.for_each_xaxis(lambda x: x.update(showgrid=False))
            fig.for_each_yaxis(lambda x: x.update(showgrid=False,ticks="inside"))
            fig.write_html(path_to_standard_plot)
            print('Saving standardized export plot to', path_to_standard_plot, sep=' ')
        else:
            print('No samples left after dropping sample NAs. Skipping project standardization plot')
    else: print('Export file for project {} not found. Please run 1.export_projects.py to export file automatically or export manually'.format(str(project_id)))


if __name__ == '__main__':
    print('Assisting standardizer completion. \nPlease enter the following information (attention: method sensitive to empty space and lower/upper case. Do not use quotation)\n')
    funcs_args=input('Press:\n 1 for printing all fields for Ecotaxa project using the API (downloading export files is not required)\n 1b for printing all possible units\n 2 for generating/saving interactive plot of raw ecotaxa export\n 3 for flagging project\n 4 for generating/saving/plotting project export standardization\n')

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
        # Example: ~/GIT/PSSdb/raw/Ecotaxa_API_pw.yaml
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
             filling_standardizer_flag_func(standardizer_path=standardizerfile_args, project_id=project_args)
        else:
             print('Missing information. Please run again')

    if funcs_args=='4':
        if len(project_args) > 0 and len(standardizerfile_args) > 0:
            standardization_func(standardizer_path=standardizerfile_args, project_id=project_args)
        else:
            print('Missing information. Please run again')