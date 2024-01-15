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
import statsmodels.api as sm

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

path_to_bins = path_to_git / 'ancillary' / 'ecopart_size_bins.tsv'
df_bins=pd.read_table(path_to_bins,sep=",")
from scipy import stats

bins=df_bins['ESD_um']
df_bins = pd.DataFrame({'sizeClasses': pd.cut(bins.to_numpy(), bins.to_numpy(), include_lowest=True).categories.values,  # Define bin categories (micrometers)
                        'range_size_bin': np.concatenate(np.diff((1 / 6) * np.pi * (np.resize( np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])), (len(bins) - 1, 2)) ** 3), axis=1)),  # cubic micrometers
                        'size_class_mid': stats.gmean(np.resize( np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])),(len(bins) - 1, 2)), axis=1)})  # Define geometrical mean of bin categories (micrometers)

# Define columns datatypes for standardized files
dtypes_dict_all = dict(
        zip(['Cruise', 'Station', 'Profile', 'Sample', 'Latitude', 'Longitude', 'Depth_min', 'Depth_max',
             'Sampling_date', 'Sampling_time', 'Volume_analyzed', 'Volume_imaged', 'ROI', 'Area', 'Pixel', 'Minor', 'Biovolume', 'ESD',
             'Category', 'Annotation', 'Sampling_type', 'Sampling_lower_size', 'Sampling_upper_size', 'ROI_number',
             'Sampling_description'],
            [str, str, str, str, float, float, float, float, str, str, float,float, str,float, float, float, float, float, str,
             str, str, float, float, int, str]))

def regress_nbss(nbss,threshold_count=0.2,threshold_size=0.2,n_bins=3):
    """
    Objective: This function computes the parameters (slope, intercept, coefficient of determination) of thresholded Normalized Biomass Size Spectrum using a log-linear regression.\n
    The function also returns the indices selected after thresholding, based on a threshold for count and size uncertainty (default is 20%) and a fixed number of empty bins (n=3)
    :param nbss: Sample-specific NBSS (2nd output, aka nbss_sum) returned by the function process_nbss_standardized_files. Attention, the dataframe should be sample-specific,otherwise it will not work. If it is not the case, use df.groupby(['Sample']).apply(lambda x: regress_nbss(x)).reset_index()
    :param threshold_count: Threshold for count uncertainty based on the number of ROI assigned to each size class. Default is 20% uncertainty. All size classes above this threshold will display empty NB (reset to NA) to account for these classes for the empty bins number
    :param threshold_size: Threshold for size uncertainty based on the instrument pixel size. Default is 20% uncertainty. All size classes above this threshold will display empty NB (reset to NA) to account for these classes for the empty bins number
    :param n_bins: Number of consecutive empty bins needed to identify the last size class for reliable regression. Default is 3 empty bins. All size classes greater than the last one recording this much empty bins  will be discarded as they are too close to the instrument detection limit.
    :return: Returns a dictionary including the total and selected abundance (NBSS integral), selected indices, and NBSS parameters
    """

    nbss = nbss.astype(dict(zip(['size_class_mid', 'NBSS'], [float] * 2))).sort_values(['size_class_mid']).reset_index().rename(columns={'index':'native_index'})

    nbss['Count_uncertainty'] = poisson.pmf(k=nbss['NBSS_count'],mu=nbss['NBSS_count']) if 'NBSS_count' in nbss.columns else 0
    nbss.loc[nbss.Count_uncertainty > threshold_count, 'NBSS'] = np.nan
    nbss['Size_uncertainty'] = norm.pdf( (1 / 6) * np.pi * (2 * (nbss.size_class_pixel.astype(float) / np.pi) ** 0.5) ** 3, loc=(1 / 6) * np.pi * (2 * (nbss.size_class_pixel.astype(float).min() / np.pi) ** 0.5) ** 3, scale=(1 / 6) * np.pi * (2 * (1 / np.pi) ** 0.5) ** 3) if 'size_class_pixel' in nbss.columns else 0
    nbss.loc[nbss.Size_uncertainty > threshold_size, 'NBSS'] = np.nan
    nbss['selected'] = (nbss.NBSS.index >= nbss.NBSS.index.to_list()[np.argmax(nbss.NBSS)])

    if any(nbss.NBSS.isnull().astype(int).groupby(nbss.NBSS.notnull().astype(int).cumsum()).sum()>=n_bins+1):
        nbss.loc[(nbss.NBSS.isnull().astype(int).groupby(nbss.NBSS.notnull().astype(int).cumsum()).sum().index[np.argwhere(nbss.NBSS.isnull().astype(int).groupby(nbss.NBSS.notnull().astype(int).cumsum()).sum().to_numpy()>=n_bins+1)[0][0]]+1):,'selected']=False
    selection=nbss.sort_values(['native_index']).selected.to_numpy()
    total_abundance=np.nansum(nbss.NBSS)
    nbss=nbss.drop(index=nbss[nbss.selected==False].index).dropna(subset=['NBSS'])
    selected_abundance=np.nansum(nbss.NBSS)
    try:
        y,x=nbss.NBSS.to_numpy(),((1/6)*np.pi*nbss.size_class_mid**3).to_numpy()
        x=sm.add_constant(x)
        wls_model =sm.WLS(y, x )
        reg =  wls_model.fit()
        return {'total_abundance':total_abundance,'selected_abundance':selected_abundance,'selection':selection,'slope':reg.params[1],'slope_sd':reg.bse[1],'intercept':reg.params[0],'intercept_sd':reg.bse[0],'R2':reg.rsquared}
    except:
        return {'total_abundance':total_abundance,'selected_abundance':selected_abundance,'selection': selection,'slope': pd.NA, 'slope_sd': pd.NA, 'intercept': pd.NA, 'intercept_sd': pd.NA,'R2': pd.NA}

def process_nbss_standardized_files(path=None,df=None,category=[],depth_selection=True):
    """
    Objective: This function computes the Normalized Biomass Size Spectrum (in \u03BCm\u00b3 dm\u207B\u00b3 \u03BCm\u207B\u00b3) for individual standardized files
    :param path: Full path to a project standardized file.
    :param df: Dataframe of standardized file if path is None
    :param category: A taxonomic group for group-specific NBSS computation. Default None will select all categories.
    :return: Returns the NBSS spectrum for individual sample
    """
    if (path) and (not df):
        columns=pd.read_table(path,sep=",", nrows=0).columns
        df=pd.read_table(path,sep=",",usecols=[column for column in ['Project_ID','Instrument','Longitude','Latitude','Station','Sample','Sampling_type','Sampling_date','Sampling_time','Depth_min','Depth_max','Volume_analyzed','Volume_imaged','ROI','ROI_number','Area','Category','Pixel','Sampling_lower_size','Sampling_upper_size'] if column in columns],dtype=dtypes_dict_all)
    instrument=df.Instrument.unique()[0]
    if (instrument in ['IFCB', 'Zooscan']) and ('ROI_number' not in df.columns):
        df = df.assign(ROI_number=1)
    # Assign biovolume, latitudinal/longitudinal/time bins as sample-specific values (so that NBSS computation is sample-specific)
    df = df.assign(Biovolume=(1 / 6) * np.pi * ((2 * ((df.Area / np.pi) ** 0.5)) ** 3),
                   ECD=2 * ((df.Area / np.pi) ** 0.5),
                   Station_location=df.Project_ID.astype(str) + '_' + df.Sample.astype(str),
                   Date_bin=pd.to_datetime(df.Sampling_date + df.Sampling_time, format="%Y%m%d%H%M%S").astype(str),
                   midLatBin=df.Latitude.astype(str), midLonBin=df.Longitude.astype(str),Sampling_type=df.Sampling_type if 'Sampling_type' in df.columns else 'nan')
    # Apply depth and artefacts filter, as well as category filter if non-null
    df_depthselection = (df.drop_duplicates(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'])[
        ['Project_ID', 'Sample', 'Depth_min', 'Depth_max']]).astype(dict(zip(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'],[str]*4))).groupby(
        ['Project_ID', 'Sample', 'Depth_min', 'Depth_max']).apply(lambda x: pd.Series({'Depth_selected': ((x.Depth_min.astype(float) < 200) & (x.Depth_max.astype(float) < 300)).values[0]})).reset_index()
    df=pd.merge(df.astype(dict(zip(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'],[str]*4))),df_depthselection,how='left',on=['Project_ID', 'Sample', 'Depth_min', 'Depth_max'])
    # Assign sampling depth range
    if instrument in ['UVP']:
        df = pd.merge(df.astype(dict(zip(['Project_ID', 'Sample'], [str] * 2))), df.astype(dict(zip(['Project_ID', 'Sample'], [str] * 2))).groupby(['Project_ID', 'Sample']).apply(lambda x: pd.Series({'Sampling_depth_range_min': x.Depth_min.astype(float).min(),'Sampling_depth_range_max': x.Depth_max.astype(float).max()})).reset_index(),how='left', on=['Project_ID', 'Sample'])
        df = df.rename(columns={'Depth_min': 'Depth_min_bin', 'Depth_max': 'Depth_max_bin','Sampling_depth_range_min': 'Depth_min','Sampling_depth_range_max': 'Depth_max'})

    df_summary = (df.drop_duplicates(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'])[['Project_ID', 'Sample', 'Depth_min', 'Depth_max']]).astype(dict(zip(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'], [str] * 4))).groupby(['Project_ID', 'Sample', 'Depth_min', 'Depth_max']).apply(lambda x: pd.Series({'Min_obs_depth': x.Depth_min.astype(float).min(), 'Max_obs_depth': x.Depth_max.astype(float).max()})).reset_index()
    df = pd.merge(df.astype(dict(zip(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'], [str] * 4))), df_summary, how='left', on=['Project_ID', 'Sample', 'Depth_min', 'Depth_max'])
    group = ['Instrument', 'Project_ID', 'Station', 'Sample', 'Date_bin', 'Latitude', 'Longitude', 'Min_obs_depth', 'Max_obs_depth', 'Sampling_lower_size', 'Sampling_upper_size']
    # Compute cumulative volume per sample/profile before dropping some observations
    df_volume = df.astype(dict(zip(group, [str] * len(group)))).groupby(group).apply(lambda x: pd.Series({'cumulative_volume':x[['Sample', 'Depth_min', 'Depth_max','Volume_imaged']].drop_duplicates().Volume_imaged.astype( float).sum()})).reset_index()
    df = pd.merge(df.astype(dict(zip(group, [str] * len(group)))), df_volume, how='left', on=group)

    df_subset = df[(df.Depth_selected == True) & (df.Sampling_type.str.lower().isin(['test', 'exp', 'junk', 'culture'])==False) & (df.Category.str.lower().apply(lambda annotation: len( re.findall(r'bead|bubble|artefact|artifact|not-living', annotation)) == 0 if str(annotation) != 'nan' else True))] if depth_selection else df[(df.Sampling_type.str.lower().isin(['test', 'exp', 'junk', 'culture'])==False) & (df.Category.str.lower().apply(lambda annotation: len( re.findall(r'bead|bubble|artefact|artifact|not-living', annotation)) == 0 if str(annotation) != 'nan' else True))]

    if len(df_subset)>0 & any(df_subset.ECD!=0):
        # Assign size bins and grouping index for each sample
        df_subset = df_subset.assign(sizeClasses=pd.cut(df_subset['ECD'], bins, include_lowest=True))
        df_subset = pd.merge(df_subset, df_bins, how='left', on=['sizeClasses'])
        # Assign a group index before groupby computation
        group=list(np.delete(group,pd.Series(group).isin(['Sampling_lower_size', 'Sampling_upper_size'])))
        df_subset = pd.merge(df_subset, df_subset.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left', on=group)

        #  Compute NBSS without applying the thresholding
        nbss_all = df_subset.astype(dict(zip(group+['sizeClasses']+ [ 'Group_index','Sampling_lower_size', 'Sampling_upper_size']+category,[str]*len(group+['sizeClasses']+ [ 'Group_index','Sampling_lower_size', 'Sampling_upper_size']+category)))).dropna(subset=['sizeClasses']).groupby(group + list(df_bins.columns) + [ 'Group_index','Sampling_lower_size', 'Sampling_upper_size']+category).apply(lambda x: pd.Series({'size_class_pixel':(x.Pixel.unique()[0]) * 1e-03 * x.size_class_mid.unique()[0],'NBSS_count': x.ROI_number.sum(), 'sum_biovolume': sum(x.Biovolume * x.ROI_number),'volume': x.cumulative_volume.unique()[0], 'NBSS': sum(x.groupby(['cumulative_volume']).apply( lambda y: (sum(y.Biovolume * y.ROI_number) / y.cumulative_volume.unique())[0])) /(x.range_size_bin.unique())[0]})).reset_index()
        #x=nbss_all.loc[list(nbss_all.groupby(group + list(df_bins.columns) + [ 'Group_index'], dropna=True).groups.values())[i]]
        nbss_sum=nbss_all.astype(dict(zip(group+['sizeClasses']+ [ 'Group_index','Sampling_lower_size', 'Sampling_upper_size']+category,[str]*len(group+['sizeClasses']+ [ 'Group_index','Sampling_lower_size', 'Sampling_upper_size']+category)))).groupby(group + list(df_bins.columns) + [ 'Group_index']+category, dropna=False).apply(
            lambda x: pd.Series({'size_class_pixel':np.nanmean(x.size_class_pixel),'volume':1 if instrument in ['Zooscan','Other','Scanner'] else x.volume.unique()[0],'NBSS_count':x.NBSS_count.sum(),'sum_biovolume':x.sum_biovolume.sum(),'NBSS':x.NBSS.sum()})).reset_index().sort_values(['Group_index','size_class_mid'])

        #nbss = df_subset.groupby(group + list(df_bins.columns) + [ 'Group_index'], dropna=False).apply(
        #lambda x: pd.Series({'NBSS_count': x.ROI_number.sum(),'sum_biovolume':sum(x.Biovolume * x.ROI_number),'volume':sum(x.cumulative_volume.unique()),'NBSS': sum(x.groupby(['cumulative_volume']).apply(lambda y:(sum(y.Biovolume * y.ROI_number) / y.cumulative_volume.unique())[0]*(y.cumulative_volume.unique()/sum(x.cumulative_volume.unique()))[0])) / (x.range_size_bin.unique())[ 0]})).reset_index()

    else:
        nbss_all=pd.DataFrame()
        nbss_sum = pd.DataFrame()

    return nbss_all,  nbss_sum

if __name__ == '__main__':
    nbss_all=pd.DataFrame()
    for instrument in ['UVP','Zooscan','IFCB']:
        print("Computing NBSS for all {} standardized datasets. Please wait".format(instrument))
        standardized_files={project:natsorted(list(Path(path_to_standard_files / Path(df_standardizer.loc[(df_standardizer.Project_ID==project) & (df_standardizer.Instrument==instrument),'Project_localpath'].values[0]).stem).rglob('standardized_project_{}_*'.format(project)))) for project in df_standardizer.loc[df_standardizer.Instrument==instrument,'Project_ID'] if len(list(Path(path_to_standard_files / Path(df_standardizer.loc[(df_standardizer.Project_ID==project) & (df_standardizer.Instrument==instrument),'Project_localpath'].values[0]).stem).rglob('standardized_project_{}_*'.format(project))))}
        standardized_files=natsorted(set(chain(*standardized_files.values())))
        #df=pd.concat(map(lambda path: (columns:=pd.read_table(path,dtype=dtypes_dict_all,sep=",").columns,pd.read_table(path,sep=",",usecols=[column for column in ['Project_ID','Instrument','Longitude','Latitude','Sample','Sampling_date','Sampling_time','Depth_min','Depth_max','Volume_imaged','ROI','ROI_number','Area','Category'] if column in columns],dtype=dtypes_dict_all))[-1],standardized_files))
        # Read all datafiles
        chunk = 1000
        nbss= pd.DataFrame() #nbss_all, nbss_max, nbss_sum, nbss_average=pd.concat(map(lambda path:process_nbss_standardized_files(path)[0],standardized_files))
        with ThreadPool() as pool:
            for result in pool.map(lambda path:process_nbss_standardized_files(path)[2],standardized_files, chunksize=chunk):#pool.map(lambda path: (columns := pd.read_table(path,sep=",", nrows=0).columns, pd.read_table(path,sep=",",usecols=[column for column in ['Project_ID','Instrument','Longitude','Latitude','Sample','Sampling_date','Sampling_time','Depth_min','Depth_max','Volume_imaged','ROI','ROI_number','Area','Category'] if column in columns],dtype=dtypes_dict_all))[-1],standardized_files, chunksize=chunk):
                nbss = pd.concat([nbss, result], axis=0)
        # Merge to other instrument
        nbss_all=pd.concat([nbss_all,nbss]).reset_index(drop=True)
    #Reset nbss_all Group_index:
    group = ['Instrument', 'Project_ID', 'Sample', 'Latitude', 'Longitude', 'volume', 'Min_obs_depth', 'Max_obs_depth']
    nbss_all = pd.merge(nbss_all.drop(columns=['Group_index']),nbss_all.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename( {'index': 'Group_index'}, axis='columns'), how='left', on=group)

    # Generate ggplot
    plot = (ggplot(nbss_all) +
        facet_wrap('~Instrument', nrow=1,scales='free') +
         geom_line(mapping=aes(x='size_class_mid', y='NBSS', group='Group_index'), size=0.01, alpha=0.6,colour='#{:02x}{:02x}{:02x}'.format(183, 200 , 190)) +
         #geom_pointrange(data=stat_nbss_all,mapping=aes(x='size_mid', group='size_mid', ymin='ymin', ymax='ymax', y='y'),alpha=0.1,colour='#{:02x}{:02x}{:02x}'.format(183, 200 , 190)) +
         stat_summary(mapping=aes(x='size_class_mid', y='NBSS', group='size_mid'),fun_data='median_hilow',geom='pointrange',alpha=0.3,size=0.05)+
         labs(x='Equivalent circular diameter (µm)',y='Normalized Biovolume Size Spectrum \n ($\mu$m$^{3}$ dm$^{-3}$ $\mu$m$^{-3}$)', title='',colour='') +
         scale_x_log10(limits=[1e+0, 1e+05], breaks=[10 ** np.arange(0, 5, step=1, dtype=np.float)][0],
                              labels=['10$^{%s}$' % int(n) for n in np.arange(0, 5, step=1)]) +
        scale_y_log10(limits=[1e-5, 1e+7], breaks=[10 ** np.arange(-5, 7, step=2, dtype=np.float)][0],
                              labels=['10$^{%s}$' % int(n) for n in np.arange(-5, 7, step=2)]) +
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
    plot[0].set_size_inches(7, 4)
    plot[0].savefig(fname='{}/GIT/PSSdb_LOV_orga/Plots/NBSS_Instrument_all.png'.format(str(Path.home())), limitsize=False, dpi=600)

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
        fig.add_trace(go.Scatter(y=subset_nbss.NBSS, x=subset_nbss.size_class_mid,
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
    fig.for_each_xaxis(lambda x: x.update(type="log",showgrid=False))
    fig.for_each_yaxis(lambda x: x.update(type="log",showgrid=False, ticks="inside"))
    fig.write_html('{}/GIT/PSSdb/figures/NBSS_Instrument_all.html'.format(str(Path.home())))