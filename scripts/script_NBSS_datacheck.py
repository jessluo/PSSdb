## Objective: This script generates a series of interactive plots to check the Normalized Biovolume Size Spectrum for individual datasets (=projects)

## Python modules:

# Data handling
import numpy as np
import pandas as pd
import re
from natsort import natsorted
from itertools import chain

# Multiprocessing mapping
from multiprocessing.pool import ThreadPool # Use: pip install multiprocess

# NBSS calculation and statistics
try:
    from funcs_NBS import *
except:
    from scripts.funcs_NBS import *

# Functions for plotting:
import colorspace  # Use: pip install git+https://github.com/retostauffer/python-colorspace
from colormap import rgb2hex
palette_depth=colorspace.sequential_hcl(name="Blues 2").colors()
palette_inst=colorspace.sequential_hcl(name="Purple-Blue").colors()
palette_latitude=colorspace.diverging_hcl(name="Purple-Brown").colors()
import plotly
import plotly.express as px # Use: pip install plotly==5.8.0
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "browser"
from plotly.subplots import make_subplots
import patchworklib as pw
pw.param["margin"] = 0.1
pw.param["dpi"] = 300
font_family='Helvetica'
import plotnine
from plotnine import *  # Python equivalent of ggplot2. Use shell command: pip install plotnine. Do not execute in Pycharm (bug- no fix yet): https://youtrack.jetbrains.com/issue/PY-42390
import geopandas as gpd
from shapely.geometry import Polygon, mapping
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in  np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)

## Workflow starts here:
path_to_git=Path('~/GIT/PSSdb').expanduser()
path_to_data=path_to_git / 'raw'
path_to_project_list=path_to_data/ "project_list_all.xlsx"
df_list=pd.concat(map(lambda x:pd.read_excel(path_to_project_list,sheet_name=x).assign(Portal=x),['ecotaxa','ecopart','ifcb']))
standardizer_files=list(path_to_data.glob('project_*_standardizer.xlsx'))
df_standardizer=pd.concat(map(lambda x:pd.read_excel(x),standardizer_files)).reset_index(drop=True)
df_standardizer['Project_files']=df_standardizer.apply(lambda x:';'.join([str(path) for path in list(path_to_standard_files.rglob('*_{}_*.tsv'.format(x['Project_ID'])))]),axis=1)
df_standardizer['Portal']=df_standardizer.Project_localpath.apply(lambda x: 'ecopart' if Path(x).stem.lower()=='ecopart' else 'ecotaxa' if (Path(x).stem.lower()=='ecotaxa') or (Path(x).stem.lower()=='ecotaxa_ecopart_consolidation')  else 'dashboard')
df_standardizer['Project_files']=df_standardizer.apply(lambda x:';'.join([str(path) for path in list(Path(x['Project_localpath']).expanduser().rglob('*raw_{}_metadata.tsv'.format(x['Project_ID'])))]) if x.Portal=='ecopart' else x.Project_files  ,axis=1)
df_standardizer=pd.merge(df_standardizer,df_list[df_list.PSSdb_access==True][['Project_ID','Project_title','Contact_name','Contact_email','Portal']],how='left',on=['Project_ID','Portal'])

path_to_standard_files=Path(path_to_data/ 'raw_standardized').expanduser()

# Define columns datatypes for standardized files
dtypes_dict_all = dict(
        zip(['Cruise', 'Station', 'Profile', 'Sample', 'Latitude', 'Longitude', 'Depth_min', 'Depth_max',
             'Sampling_date', 'Sampling_time', 'Volume_analyzed', 'Volume_imaged', 'ROI', 'Area', 'Pixel', 'Minor', 'Biovolume', 'ESD',
             'Category', 'Annotation', 'Sampling_type', 'Sampling_lower_size', 'Sampling_upper_size', 'ROI_number',
             'Sampling_description'],
            [str, str, str, str, float, float, float, float, str, str, float,float, str,float, float, float, float, float, str,
             str, str, float, float, int, str]))

for instrument in ['UVP','Zooscan','IFCB']:
    print("Generating interactive plot for all {} standardized datasets. Please wait".format(instrument))
    standardized_files={project:natsorted(list(Path(path_to_standard_files / Path(df_standardizer.loc[(df_standardizer.Project_ID==project) & (df_standardizer.Instrument==instrument),'Project_localpath'].values[0]).stem).rglob('standardized_project_{}_*'.format(project)))) for project in df_standardizer.loc[df_standardizer.Instrument==instrument,'Project_ID'] if len(list(Path(path_to_standard_files / Path(df_standardizer.loc[(df_standardizer.Project_ID==project) & (df_standardizer.Instrument==instrument),'Project_localpath'].values[0]).stem).rglob('standardized_project_{}_*'.format(project))))}
    standardized_files=natsorted(set(chain(*standardized_files.values())))
    # df=pd.concat(map(lambda path: (columns:=pd.read_table(path,dtype=dtypes_dict_all,sep=",").columns,pd.read_table(path,sep=",",usecols=[column for column in ['Project_ID','Instrument','Longitude','Latitude','Sample','Sampling_date','Sampling_time','Depth_min','Depth_max','Volume_imaged','ROI','ROI_number','Area','Category'] if column in columns],dtype=dtypes_dict_all))[-1],standardized_files))
    chunk = 1000
    df = pd.DataFrame()
    with ThreadPool() as pool:
        for result in pool.map(lambda path: (columns := pd.read_table(path,sep=",", nrows=0).columns, pd.read_table(path,sep=",",usecols=[column for column in ['Project_ID','Instrument','Longitude','Latitude','Sample','Sampling_date','Sampling_time','Depth_min','Depth_max','Volume_imaged','ROI','ROI_number','Area','Category'] if column in columns],dtype=dtypes_dict_all))[-1],standardized_files, chunksize=chunk):
            df = pd.concat([df, result], axis=0)
    df = df.reset_index(drop=True)
