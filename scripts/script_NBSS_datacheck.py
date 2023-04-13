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

# Gridding and NBSS calculation and statistics
try:
    from funcs_NBS import *
    from funcs_gridding import *
except:
    from scripts.funcs_NBS import *
    from scripts.funcs_gridding import *

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
path_to_standard_files=Path(path_to_data/ 'raw_standardized').expanduser()

df_list=pd.concat(map(lambda x:pd.read_excel(path_to_project_list,sheet_name=x).assign(Portal=x),['ecotaxa','ecopart','ifcb']))
standardizer_files=list(path_to_data.glob('project_*_standardizer.xlsx'))
df_standardizer=pd.concat(map(lambda x:pd.read_excel(x),standardizer_files)).reset_index(drop=True)
df_standardizer['Project_files']=df_standardizer.apply(lambda x:';'.join([str(path) for path in list(path_to_standard_files.rglob('*_{}_*.tsv'.format(x['Project_ID'])))]),axis=1)
df_standardizer['Portal']=df_standardizer.Project_localpath.apply(lambda x: 'ecopart' if Path(x).stem.lower()=='ecopart' else 'ecotaxa' if (Path(x).stem.lower()=='ecotaxa') or (Path(x).stem.lower()=='ecotaxa_ecopart_consolidation')  else 'dashboard')
df_standardizer['Project_files']=df_standardizer.apply(lambda x:';'.join([str(path) for path in list(Path(x['Project_localpath']).expanduser().rglob('*raw_{}_metadata.tsv'.format(x['Project_ID'])))]) if x.Portal=='ecopart' else x.Project_files  ,axis=1)
df_standardizer=pd.merge(df_standardizer,df_list[df_list.PSSdb_access==True][['Project_ID','Project_title','Contact_name','Contact_email','Portal']],how='left',on=['Project_ID','Portal'])

path_to_bins = path_to_data / 'ecopart_size_bins.tsv'
df_bins=pd.read_table(path_to_bins,sep=",")
from scipy import stats

bins=df_bins['biovol_um3']
df_bins = pd.DataFrame({'sizeClasses': pd.cut(bins.to_numpy(), bins.to_numpy(), include_lowest=True).categories.values,  # Define bin categories (cubic micrometers)
                        'range_size_bin': np.concatenate(np.diff((np.resize(np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])),(len(bins) - 1, 2))), axis=1)),  # cubic micrometers
                        'size_class_mid': stats.gmean(np.resize( np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])),(len(bins) - 1, 2)), axis=1)})  # Define geometrical mean of bin categories (cubic micrometers)

# Define columns datatypes for standardized files
dtypes_dict_all = dict(
        zip(['Cruise', 'Station', 'Profile', 'Sample', 'Latitude', 'Longitude', 'Depth_min', 'Depth_max',
             'Sampling_date', 'Sampling_time', 'Volume_analyzed', 'Volume_imaged', 'ROI', 'Area', 'Pixel', 'Minor', 'Biovolume', 'ESD',
             'Category', 'Annotation', 'Sampling_type', 'Sampling_lower_size', 'Sampling_upper_size', 'ROI_number',
             'Sampling_description'],
            [str, str, str, str, float, float, float, float, str, str, float,float, str,float, float, float, float, float, str,
             str, str, float, float, int, str]))
nbss_all=pd.DataFrame()
for instrument in ['UVP','Zooscan','IFCB']:
    print("Computing NBSS for all {} standardized datasets. Please wait".format(instrument))
    standardized_files={project:natsorted(list(Path(path_to_standard_files / Path(df_standardizer.loc[(df_standardizer.Project_ID==project) & (df_standardizer.Instrument==instrument),'Project_localpath'].values[0]).stem).rglob('standardized_project_{}_*'.format(project)))) for project in df_standardizer.loc[df_standardizer.Instrument==instrument,'Project_ID'] if len(list(Path(path_to_standard_files / Path(df_standardizer.loc[(df_standardizer.Project_ID==project) & (df_standardizer.Instrument==instrument),'Project_localpath'].values[0]).stem).rglob('standardized_project_{}_*'.format(project))))}
    standardized_files=natsorted(set(chain(*standardized_files.values())))
    df=pd.concat(map(lambda path: (columns:=pd.read_table(path,dtype=dtypes_dict_all,sep=",").columns,pd.read_table(path,sep=",",usecols=[column for column in ['Project_ID','Instrument','Longitude','Latitude','Sample','Sampling_date','Sampling_time','Depth_min','Depth_max','Volume_imaged','ROI','ROI_number','Area','Category'] if column in columns],dtype=dtypes_dict_all))[-1],standardized_files))
    # Read all datafiles
    """
    chunk = 1000
    df = pd.DataFrame()
    with ThreadPool() as pool:
        for result in pool.map(lambda path: (columns := pd.read_table(path,sep=",", nrows=0).columns, pd.read_table(path,sep=",",usecols=[column for column in ['Project_ID','Instrument','Longitude','Latitude','Sample','Sampling_date','Sampling_time','Depth_min','Depth_max','Volume_imaged','ROI','ROI_number','Area','Category'] if column in columns],dtype=dtypes_dict_all))[-1],standardized_files, chunksize=chunk):
            df = pd.concat([df, result], axis=0)
    """
    df = df.reset_index(drop=True)
    if (instrument in ['IFCB','Zooscan']) and ('ROI_number' not in df.columns):
        df=df.assign(ROI_number=1)
    # Assign biovolume, latitudinal/longitudinal/time bins as sample-specific values (so that NBSS computation is sample-specific)
    df=df.assign(Biovolume=(1 / 6) * np.pi * ((2 * ((df.Area / np.pi) ** 0.5))**3),Station_location=df.Project_ID.astype(str)+'_'+df.Sample.astype(str),date_bin=pd.to_datetime(df.Sampling_date+df.Sampling_time,format="%Y%m%d%H%M%S").astype(str),midLatBin=df.Latitude.astype(str), midLonBin=df.Longitude.astype(str))
    # Apply depth and artefacts filter
    df_depthselection=(df.drop_duplicates(['Project_ID','Sample','Depth_min','Depth_max'])[['Project_ID','Sample','Depth_min','Depth_max']]).groupby(['Project_ID','Sample','Depth_min','Depth_max']).apply(lambda x: pd.Series({'Depth_selected':((x.Depth_min<200) & (x.Depth_max<300)).values[0]})).reset_index()
    df_subset=df[(df.Sample.astype(str).isin(list(df_depthselection[df_depthselection.Depth_selected==True].Sample.astype(str)))) & (df.Category.str.lower().apply(lambda annotation:len(re.findall(r'bead|bubble|artefact|artifact|not-living',annotation))==0 if str(annotation)!='nan' else True))]
    # Assign sampling depth range
    df_subset_summary=(df_subset.drop_duplicates(['Project_ID','Sample','Depth_min','Depth_max'])[['Project_ID','Sample','Depth_min','Depth_max']]).groupby(['Project_ID','Sample','Depth_min','Depth_max']).apply(lambda x: pd.Series({'Min_obs_depth': x.Depth_min.min(), 'Max_obs_depth': x.Depth_max.max()})).reset_index()
    df_subset =pd.merge(df_subset,df_subset_summary,how='left',on=['Project_ID','Sample','Depth_min','Depth_max'])
    # Assign size bins and grouping index for each sample
    df_subset=df_subset.assign(sizeClasses=pd.cut(df_subset['Biovolume'],bins,include_lowest=True))
    df_subset =pd.merge(df_subset,df_bins,how='left',on=['sizeClasses'])
    group=['Project_ID','Sample', 'Latitude', 'Longitude', 'Volume_imaged','Min_obs_depth','Max_obs_depth']
    df_subset =pd.merge(df_subset , df_subset.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename({'index': 'Group_index'}, axis='columns'),how='left', on=group)
    # Compute cumulative volume per sample/profile
    df_subset = pd.merge(df_subset ,df_subset.groupby(group).apply(lambda x: pd.Series({'cumulative_volume': x[['Sample', 'Depth_min', 'Depth_max','Volume_imaged']].drop_duplicates().Volume_imaged.sum()})).reset_index(),how='left', on=group)
    # Compute NBSS without applying the thresholding
    nbss=df_subset.groupby(group +list(df_bins.columns)+ ['cumulative_volume','Group_index'],dropna=False).apply(lambda x: pd.Series({
            'NBSS_count':x.ROI_number.sum(),
            'NBSS': (sum(x.Biovolume*x.ROI_number) / x.cumulative_volume.unique() / ( (x.range_size_bin.unique())))[0]})).reset_index()
    nbss=nbss.assign(Instrument=instrument)
    # Merge to other instrument
    nbss_all=pd.concat([nbss_all,nbss])

# Generate ggplot
nbss_all['size_mid']=((nbss_all['size_class_mid']*6)/np.pi)**(1/3)
stat_nbss_all=nbss_all.groupby('size_mid').apply(lambda x:pd.Series({"y":np.nanmedian(x.NBSS),"ymin":np.nanquantile(x.NBSS,0.05),"ymax":np.nanquantile(x.NBSS,0.95)})).reset_index()
plot = (ggplot(nbss_all) +
        facet_wrap('~Instrument', nrow=1) +
         geom_line(mapping=aes(x='size_mid', y='NBSS', group='Group_index'), size=0.15, alpha=0.01) +
         geom_pointrange(data=stat_nbss_all,mapping=aes(x='size_mid', group='size_mid', ymin='ymin', ymax='ymax', y='y')) +
         labs(x='Equivalent circular diameter (µm)',y='Normalized Biovolume Size Spectrum \n ($\mu$m$^{3}$ dm$^{-3}$ $\mu$m$^{-3}$)', title='',colour='') +
         scale_x_log10(limits=[1e+0, 1e+05], breaks=[10 ** np.arange(0, 5, step=1, dtype=np.float)][0],
                              labels=['10$^{%s}$' % int(n) for n in np.arange(0, 5, step=1)]) +
        scale_y_log10(limits=[1e-5, 1e+5], breaks=[10 ** np.arange(-5, 5, step=2, dtype=np.float)][0],
                              labels=['10$^{%s}$' % int(n) for n in np.arange(-5, 5, step=2)]) +
        theme(axis_ticks_direction="inout", legend_direction='horizontal', legend_position='top',
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family=font_family, size=10),
                      legend_text=element_text(family=font_family, size=10),
                      axis_title=element_text(family=font_family, size=10),
                      axis_text_x=element_text(family=font_family, size=10),
                      axis_text_y=element_text(family=font_family, size=10, rotation=90),
                      plot_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(5, 3)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/NBSS_Instrument_all.png'.format(str(Path.home())), limitsize=False, dpi=600)

# Generate interactive plot
subplot_specs = [[{"type": "scattergeo"}, {"type": "scatter"}]]
subplot_titles = ['', '']
fig = make_subplots(rows=1, cols=2, specs=subplot_specs, subplot_titles=subplot_titles, column_widths=[0.5, 0.5])
# Map, left panel
nbss_all['color'] = pd.Categorical(nbss_all.Instrument, categories=['Zooscan', 'UVP', 'IFCB']).rename_categories({'UVP': 'rgba(255,255,255,100)', 'Zooscan': 'rgba(145, 111, 138, 100)','IFCB': 'rgba(211, 95, 95, 100)'}).astype(str)

for i, instrument in enumerate(nbss_all.Instrument.unique()):
    data_geo = dict(type='scattergeo', name=instrument,
                    lon=nbss_all[nbss_all.Instrument == instrument].dropna( subset=['Latitude', 'Longitude']).Longitude,
                    lat=nbss_all[nbss_all.Instrument == instrument].dropna(subset=['Latitude', 'Longitude']).Latitude,
                    hovertext='Project ID: ' + nbss_all[nbss_all.Instrument == instrument].Project_ID.astype(str),
                    hoverinfo="text",
                    marker=dict(color=nbss_all[nbss_all.Instrument == instrument].color, size=1.75), legendrank=i,
                    legendgroup=instrument,
                    geojson="natural earth", showlegend=True, geo='geo')
    layout = dict()
    layout['geo'] = dict(lonaxis_range=[-180, 180], domain=dict(x=[0.0, 0.4], y=[0.0, 1.0]), lataxis_range=[-80, 80])
    fig.add_trace(data_geo, row=1, col=1).update_layout(go.Figure(data=[data_geo, data_geo], layout=layout).layout)
    # Scatterplot 1, NBSS
    for v in nbss_all[nbss_all.Instrument == instrument].Group_index.unique():
        subset_nbss = nbss_all[nbss_all.Instrument == instrument][nbss_all[nbss_all.Instrument == instrument]['Group_index'] == v].sort_values(['size_mid']).dropna()
        show = True if v == nbss_all[nbss_all.Instrument == instrument].Group_index.unique()[0] else False
        fig.add_trace(go.Scatter(y=subset_nbss.NBSS, x=subset_nbss.size_mid,
                                 name=subset_nbss.Project_ID.astype(str).unique()[0],
                                 line=dict(color=subset_nbss.color.unique()[0], width=0.15),
                                 hovertext='Project ID: ' + subset_nbss.Project_ID.astype(str).unique()[0] + '<br> Sample ID: ' + subset_nbss.Sample.astype(str).unique()[0],
                                 hoverinfo="text",
                                 legendrank=i,
                                 legendgroup=instrument,
                                 mode='lines', showlegend=False, visible=True), row=1, col=2)

fig.update_layout( title={'text': 'Pelagic Size Structure Database (PSSdb)', 'xanchor': 'center', 'yanchor': 'top', 'x': 0.5},
    xaxis={'title': u'Equivalent circular diameter (\u03BCm)', 'tickfont': dict(size=10), 'titlefont': dict(size=12)},
    yaxis={"tickmode": "array",'title': u'Normalized Biovolume Size Spectrum <br> (\u03BCm\u00b3 m\u207B\u00b3 \u03BCm\u207B\u00b3)','tickfont': dict(size=10), 'titlefont': dict(size=12)}, boxmode='group')
fig.update_yaxes(type="log", row=1, col=2)
fig.update_xaxes(type="log", row=1, col=2, ticks="inside")
fig.update_yaxes(type="log", row=2, col=2)
fig.update_xaxes(type="log", row=2, col=2, ticks="inside")
fig.for_each_xaxis(lambda x: x.update(showgrid=False))
fig.for_each_yaxis(lambda x: x.update(showgrid=False, ticks="inside"))
fig.write_html('{}/GIT/PSSdb/figures/NBSS_Instrument_all.html'.format(str(Path.home())))