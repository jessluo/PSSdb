## Objective: This script generates all the summary statistics and plots used for PSSdb first-release data paper

## Modules:

# Modules for data and path handling:
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from wcmatch.pathlib import Path # Handling of path object
import os
import yaml
import math as m
path_to_config=Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg= yaml.safe_load(config_file)
try:
    from funcs_standardize_projects import *
except:
    from scripts.funcs_standardize_projects import *

# Functions for plotting:
from plotnine import *
import seaborn as sns
import matplotlib.pyplot as plt
palette=list((sns.color_palette("Blues",4).as_hex()))[0:1]+list((sns.color_palette("PuRd",4).as_hex()))
fontname = 'serif'
from scipy.stats import spearmanr
import geopandas as gpd
oceans = gpd.read_file(list(Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent.rglob('goas_v01.shp'))[0])
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)
theme_paper=theme(axis_ticks_direction="inout",
              panel_grid=element_blank(),
              axis_line = element_line(colour = "black"),
              panel_background=element_rect(fill='white'),
              panel_border=element_rect(color='black'),
              legend_title=element_text(family="serif", size=8),
              legend_position='top',
              legend_text=element_text(family="serif", size=8),
              axis_title=element_text(family="serif", size=8),
              axis_text_x=element_text(family="serif", size=8),
              axis_text_y=element_text(family="serif", size=8, rotation=90),
              plot_background=element_rect(fill='white'), strip_background=element_rect(fill='white'))
colors = { 'Scanner': 'red', 'UVP': 'blue','IFCB': 'green'}

## Workflow starts here:

path_to_datafile=(Path(cfg['git_dir']).expanduser()/ cfg['dataset_subdir']) / 'NBSS_data' / 'NBSS_ver_10_2023'
path_files= list(set(path_to_datafile.rglob('*_1b_*.csv')) - set(path_to_datafile.rglob('*Sensitivity_analysis/*')))
path_files_1a=list(set(path_to_datafile.rglob('*_1a_*.csv')) - set(path_to_datafile.rglob('*Sensitivity_analysis/*')))
df=pd.concat(map(lambda path: pd.read_table(path,sep=',').assign(Instrument=path.name[0:path.name.find('_')]),path_files)).drop_duplicates().reset_index(drop=True)
#df['intercept_NB_mean']=np.log10(np.exp(df.intercept_NB_mean))
df_cor=df[['NBSS_slope_mean','PSD_slope_mean','NBSS_intercept_mean','PSD_intercept_mean','NBSS_r2_mean','PSD_r2_mean']].corr(method='spearman')
grouping = ['year', 'month', 'latitude', 'longitude', 'Instrument']
df = pd.merge(df, df.drop_duplicates(subset=grouping, ignore_index=True)[grouping].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left', on=grouping)
df_stat=df.groupby(['Instrument'])[['NBSS_slope_mean','PSD_slope_mean','NBSS_intercept_mean','PSD_intercept_mean','NBSS_r2_mean','PSD_r2_mean']].describe()
df_1a = pd.concat(map(lambda path: pd.read_table(path,sep=',').assign(Instrument=path.name[0:path.name.find('_')]),path_files_1a)).drop_duplicates().reset_index(drop=True)
grouping = ['year', 'month', 'latitude', 'longitude', 'Instrument']
df_1a = pd.merge(df_1a, df_1a.drop_duplicates(subset=grouping, ignore_index=True)[grouping].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left', on=grouping)

# Generate maps and heatmap of the spatio-temporal coverage for PSSdb products 1


# Fig 2. Spatial coverage
df_summary=df.groupby(['Instrument','latitude','longitude']).agg({'Group_index':'count'}).reset_index().rename(columns={'Group_index':'count'})
plot=(ggplot(data=df_summary)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="longitude",y="latitude",fill='count'),shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1) +
  coord_cartesian(expand = 0)+
 geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), fill='black', color='black') +
 labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)') +
 scale_fill_gradientn(trans='log',colors=palette)+
theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(10.5,3)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Spatial_coverage_PSSdb.svg'.format(str(Path.home())), limitsize=False, dpi=600)

# Fig 2. Temporal coverage

df_summary=df.groupby(['Instrument','month','year']).agg({'Group_index':'count'}).reset_index().rename(columns={'Group_index':'count'})
plot=(ggplot(data=df_summary)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="month",y="year",fill='count'),size=5,shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1) +
 labs(x=r'Month', y=r'Year)') +
 scale_fill_gradientn(trans='log10',colors=palette)+scale_y_continuous(breaks=[1,4,7,10],labels=['Jan','Apr','Jul','Oct'])+
theme_paper).draw(show=False, return_ggplot=True)

plot[0].set_size_inches(6.5,2.8)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_2_Temporal_coverage_PSSdb.svg'.format(str(Path.home())), limitsize=False, dpi=600)

#Fig. 3a  Average NBSS and supplemental figure 2 Average PSD
# Normalized biovolume vs size
plot = (ggplot(data=df_1a)+
        geom_line(df_1a,aes(x='equivalent_circular_diameter_mean', y='normalized_biovolume_mean', color='Instrument',group='Group_index'), alpha=0.02, size = 0.1) +
        geom_point(aes(x='equivalent_circular_diameter_mean', y='normalized_biovolume_mean', color='Instrument'),size = 0.05, alpha=0.01, shape = 'o')+
        stat_summary(data=df_1a[df_1a.groupby(['Instrument']).equivalent_circular_diameter_mean.transform(lambda x: x.astype(str).isin(pd.Series(x.value_counts(normalize=True)[x.value_counts(normalize=True)>=np.quantile(x.value_counts(normalize=True),0.5)].index).astype(str)))],mapping=aes(x='equivalent_circular_diameter_mean', y='normalized_biovolume_mean', color='Instrument'),geom='line', fun_y=np.nanmedian, size = 1)+
        labs(y=r'Normalized Biovolume ($\mu$m$^{3}$ L$^{-1}$ $\mu$m$^{-3}$)', x=r'Equivalent circular diameter ($\mu$m)')+
        scale_color_manual(values = colors)+
        scale_y_log10(breaks=[10**np.arange(-5,7,step=2, dtype=np.float)][0],labels=['10$^{%s}$'% int(n) for n in np.arange(-5,7,step=2)])+
        scale_x_log10(breaks=[size  for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))],labels= [size if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))])+
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(4.5,3)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/fig_3_NBSS.pdf'.format(str(Path.home())), dpi=300)


#Normalized abundance vs size
plot = (ggplot(data=df_1a)+
        geom_line(df_1a,aes(x='equivalent_circular_diameter_mean', y='normalized_abundance_mean', color='Instrument',group='Group_index'), alpha=0.02, size = 0.1) +
        geom_point(aes(x='equivalent_circular_diameter_mean', y='normalized_abundance_mean', color='Instrument'),size = 0.05, alpha=0.01, shape = 'o')+
        stat_summary(data=df_1a[df_1a.groupby(['Instrument']).equivalent_circular_diameter_mean.transform(lambda x: x.astype(str).isin(pd.Series(x.value_counts(normalize=True)[x.value_counts(normalize=True)>=np.quantile(x.value_counts(normalize=True),0.5)].index).astype(str)))],mapping=aes(x='equivalent_circular_diameter_mean', y='normalized_abundance_mean', color='Instrument'),geom='line', fun_y=np.nanmean, size = 1)+
        labs(y=r'Normalized Abundance (particles L$^{-1}$ $\mu$m$^{-1}$)', x=r'Equivalent circular diameter ($\mu$m)')+
        scale_color_manual(values = colors)+
        scale_y_log10(breaks=[10**np.arange(-5,7,step=2, dtype=np.float)][0],labels=['10$^{%s}$'% int(n) for n in np.arange(-5,7,step=2)])+
        scale_x_log10(breaks=[size  for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))],labels= [size if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))])+
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(4.5,3)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/fig_3_PSD.pdf'.format(str(Path.home())), dpi=300)
plot = (ggplot(data=df_1a)+
        geom_line(df_1a,aes(x='biovolume_size_class', y='normalized_abundance_mean', color='Instrument',group='Group_index'), alpha=0.02, size = 0.1) +
        geom_point(aes(x='biovolume_size_class', y='normalized_abundance_mean', color='Instrument'),size = 0.05, alpha=0.01, shape = 'o')+
        stat_summary(data=df_1a[df_1a.groupby(['Instrument']).equivalent_circular_diameter_mean.transform(lambda x: x.astype(str).isin(pd.Series(x.value_counts(normalize=True)[x.value_counts(normalize=True)>=np.quantile(x.value_counts(normalize=True),0.5)].index).astype(str)))],mapping=aes(x='biovolume_size_class', y='normalized_abundance_mean', color='Instrument'),geom='line', fun_y=np.nanmean, size = 1)+
        labs(y=r'Normalized Abundance (particles L$^{-1}$ $\mu$m$^{-1}$)', x=r'Biovolume ($\mu$m$^{3}$)')+
        scale_color_manual(values = colors)+
        scale_y_log10(breaks=[10**np.arange(-5,7,step=2, dtype=np.float)][0],labels=['10$^{%s}$'% int(n) for n in np.arange(-5,7,step=2)])+
        scale_x_log10(labels=lambda bk: [f'{int(size):,}'  for size in bk])+
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(4.5,3)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/fig_3_PSD_biovolume.pdf'.format(str(Path.home())), dpi=100)


#Comparison between NBSS and PSD slope
slope, intercept, r_value, p_value, std_err = stats.linregress(x=df.dropna(subset=['NBSS_slope_mean','PSD_slope_mean'])['NBSS_slope_mean'],y=df.dropna(subset=['NBSS_slope_mean','PSD_slope_mean'])['PSD_slope_mean'])
plot = (ggplot(data=df)+
        geom_abline(slope=slope, intercept=intercept, alpha=0.1) +
        geom_point(aes(x='NBSS_slope_mean', y='PSD_slope_mean', color='Instrument'),size = 1, alpha=0.1, shape = 'o')+
        #geom_smooth(aes(x='NBSS_slope_mean', y='PSD_slope_mean'), method='lm', linetype='dashed', size = 0.5)+
        labs(y='PSD slope ( L$^{-1}$ $\mu$m$^{-1}$)', x=r'NBSS slope ( L$^{-1}$ $\mu$m$^{-3}$)')+
        annotate('text', label='y = '+str(np.round(slope, 2))+'x '+str(np.round(intercept, 2))+ ', R$^{2}$ = '+str(np.round(r_value, 3)),  x = np.nanquantile(df.NBSS_slope_mean,[0.95]),y = np.nanquantile(df.PSD_slope_mean,[0.95]),ha='right')+
        scale_color_manual(values = colors)+scale_x_continuous(limits=np.nanquantile(df.NBSS_slope_mean,[0.05,0.95]))+
        scale_y_continuous(limits=np.nanquantile(df.PSD_slope_mean,[0.05,0.95]))+
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(3,3)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/fig_3_slopes.svg'.format(str(Path.home())), dpi=600)
# Add density
plot = (ggplot(data=df.dropna(subset=['NBSS_slope_mean']))+
        geom_histogram(aes(x='NBSS_slope_mean', fill='Instrument'),color='#ffffff00', alpha=0.3)+
        labs(y='Count', x=r'NBSS slope ( L$^{-1}$ $\mu$m$^{-3}$)')+
        scale_fill_manual(values = colors)+#scale_x_continuous(limits=np.nanquantile(df.NBSS_slope_mean,[0.05,0.95]))+
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(3,1)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/fig_3_NBSSslopes_density.svg'.format(str(Path.home())), dpi=600)
plot = (ggplot(data=df)+
        geom_histogram(aes(x='PSD_slope_mean', fill='Instrument'),color='#ffffff00',size = 1, alpha=0.1)+
        labs(y='Density', x=r'PSD slope ( L$^{-1}$ $\mu$m$^{-1}$)')+
        scale_fill_manual(values = colors)+#scale_x_continuous(limits=np.nanquantile(df.PSD_slope_mean,[0.95,0.05]))+
        #coord_flip()+ flipping messes up the rendering
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(3,1)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/fig_3_PSDslopes_density.svg'.format(str(Path.home())), dpi=600)

#Comparison between NBSS and PSD intercept

slope, intercept, r_value, p_value, std_err = stats.linregress(x=df.dropna(subset=['NBSS_intercept_mean','PSD_intercept_mean'])['NBSS_intercept_mean'],y=df.dropna(subset=['NBSS_intercept_mean','PSD_intercept_mean'])['PSD_intercept_mean'])
plot = (ggplot(data=df)+
        geom_abline(slope=slope, intercept=intercept, alpha=0.1) +
        geom_point(aes(x='NBSS_intercept_mean', y='PSD_intercept_mean', color='Instrument'),size = 1, alpha=0.1, shape = 'o')+
        annotate('text', label='y = '+str(np.round(slope, 2))+'x + '+str(np.round(intercept, 2))+ ', R$^{2}$ = '+str(np.round(r_value, 3)),  x = np.nanquantile(df.NBSS_intercept_mean,[0.95]),y =  np.nanquantile(df.PSD_intercept_mean,[0.95]),ha='right')+
        labs(y='PSD intercept (particles L$^{-1}$ $\mu$m$^{-1}$)', x=r'NBSS intercept ($\mu$m$^{3}$ L$^{-1}$ $\mu$m$^{-3}$)')+
        scale_color_manual(values = colors)+
        scale_x_continuous(limits=np.nanquantile(df.NBSS_intercept_mean,[0.05,0.95]),labels=lambda bk:['10$^{%s}$'% int(size) for size in bk])+
        scale_y_continuous(limits=np.nanquantile(df.PSD_intercept_mean,[0.05,0.95]),labels=lambda bk:['10$^{%s}$'% int(size) for size in bk])+
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(3,3)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/fig_3_intercepts.svg'.format(str(Path.home())), dpi=600)
# Add density
plot = (ggplot(data=df)+
        geom_density(aes(x='NBSS_intercept_mean', fill='Instrument'),color='#ffffff00',size = 1, alpha=0.1)+
        #geom_smooth(aes(x='NBSS_slope_mean', y='PSD_slope_mean'), method='lm', linetype='dashed', size = 0.5)+
        labs(y='Density', x=r'NBSS intercept ($\mu$m$^{3}$ L$^{-1}$ $\mu$m$^{-3}$)')+
        scale_fill_manual(values = colors)+scale_x_continuous(limits=np.nanquantile(df.NBSS_intercept_mean,[0.05,0.95]))+
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(3,1)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/fig_3_NBSSintercepts_density.svg'.format(str(Path.home())), dpi=600)
plot = (ggplot(data=df)+
        geom_density(aes(x='PSD_intercept_mean', fill='Instrument'),color='#ffffff00',size = 1, alpha=0.1)+
        labs(y='Density', x=r'PSD intercept (particles L$^{-1}$ $\mu$m$^{-1}$)')+
        scale_fill_manual(values = colors)+scale_x_continuous(limits=np.nanquantile(df.PSD_intercept_mean,[0.95,0.05]))+
        #coord_flip()+ flipping messes up the rendering
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(3,1)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/fig_3_PSDintercepts_density.svg'.format(str(Path.home())), dpi=600)


# Fig 5. Latitudinal trends
plot=(ggplot(data=df.melt(id_vars=['Instrument','year', 'month', 'longitude','latitude', 'min_depth', 'max_depth'],value_vars=['NBSS_slope_mean',  'NBSS_intercept_mean','NBSS_r2_mean']))+
       facet_grid('Instrument~pd.Categorical(variable,["NBSS_slope_mean","NBSS_intercept_mean","NBSS_r2_mean"])',scales='free')+
       stat_summary( aes(x="latitude",y="value"),alpha=1,size=0.05) +
       labs(y=r'', x=r'Latitude ($^{\circ}$N)') +scale_x_continuous(limits=[-90,90], breaks = np.arange(-90, 100, 10))+
       scale_color_manual(values={'IFCB':'#{:02x}{:02x}{:02x}'.format(111, 145 , 111),'UVP':'#{:02x}{:02x}{:02x}'.format(147,167,172),'Scanner':'#{:02x}{:02x}{:02x}'.format(95,141,211)})+
       coord_flip()+
       theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(8,12)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Figure_5.pdf'.format(str(Path.home())),  dpi=600, bbox_inches='tight')


# Fig 6. Climatology per basin
df=pd.merge(df.astype(dict(zip(['longitude', 'latitude'],[str]*2))),df.astype(dict(zip(['longitude', 'latitude'],[str]*2))).drop_duplicates(subset=['longitude', 'latitude'], ignore_index=True)[['longitude', 'latitude']].reset_index().rename({'index': 'Group_index'}, axis='columns'),how='left',on=['longitude', 'latitude']).astype(dict(zip(['longitude', 'latitude'],[float]*2)))
df['Study_area']=df.ocean
#gdf = gpd.GeoDataFrame(df[['Group_index','longitude', 'latitude']].drop_duplicates().dropna(), geometry=gpd.points_from_xy( df[[ 'Group_index','longitude', 'latitude']].drop_duplicates().dropna().longitude,df[['Group_index','longitude', 'latitude']].drop_duplicates().dropna().latitude))
#df['Study_area'] = pd.merge(df, gpd.tools.sjoin_nearest(gdf, oceans, how='left')[['Group_index', 'name']], how='left',on='Group_index')['name'].astype(str)
data=df[(df.groupby(['Instrument','Study_area'])['month'].transform(lambda x: len(x.unique())>9)) &(df.Study_area.isin(['Indian Ocean','Mediterranean Region','North Atlantic Ocean','North Pacific Ocean','South Pacific Ocean']))].melt(id_vars=['Instrument', 'Group_index', 'Study_area','year', 'month', 'latitude', 'longitude', 'min_depth', 'max_depth'],value_vars=['NBSS_slope_mean',  'NBSS_intercept_mean','NBSS_r2_mean']).astype({'month':float,'value':float}).dropna(subset=['value'])
plot=(ggplot(data)+
    facet_grid('pd.Categorical(variable,categories=["NBSS_slope_mean",  "NBSS_intercept_mean","NBSS_r2_mean"])~pd.Categorical(Study_area,categories=["Indian Ocean","Mediterranean Region","North Atlantic Ocean","North Pacific Ocean","South Pacific Ocean"])',scales='free')+
    stat_summary( aes(x="month",y="value",color='Instrument'),alpha=0.4) +
    stat_summary(aes(x="month", y="value", color='Instrument'), alpha=1,geom='line') +
    stat_summary(aes(x="month", y="value", fill='Instrument'),geom='ribbon', alpha=0.3) +
labs(x=r'', y=r'') +
scale_color_manual(values=colors)+
scale_fill_manual(values=colors)+
scale_x_continuous(breaks=np.arange(1,13,3),labels=[calendar.month_abbr[month] for month in np.arange(1,13,3)])+
theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches( 8,6)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Figures_6.svg'.format(str(Path.home())),  dpi=600, bbox_inches='tight')



# Suppl Fig A1 (Sampling efforts):
path_to_project_list=Path(cfg['git_dir']).expanduser()/ cfg['proj_list']
df_list=pd.concat(map(lambda x:pd.read_excel(path_to_project_list,sheet_name=x).assign(Portal=x),['ecotaxa','ecopart','ifcb']))

path_to_datafile=(Path(cfg['raw_dir']).expanduser() /cfg['flag_subdir']).expanduser()
path_flags=list(Path(path_to_datafile).rglob('project_*_flags.csv'))

df_standardizer=pd.concat(map(lambda path: pd.read_table(path,sep=","),path_flags)).reset_index(drop=True)
# Select only standard files for the instrument of interest
path_files=[path for project in df_standardizer.index  for path in path_files if path in list(Path(df_standardizer.loc[project,'Project_localpath'].replace('raw',cfg['standardized_subdir'])).expanduser().rglob('standardized_project_{}_*'.format(df_standardizer.loc[project,'Project_ID'])))]
df=pd.concat(map(lambda path: pd.read_table(path,sep=',',usecols=['Instrument','Project_ID','Cruise','Sample','Latitude','Longitude','Sampling_date','Sampling_time','Depth_min','Depth_max','Sampling_description'], parse_dates={'Sampling_datetime': ['Sampling_date', 'Sampling_time']}),path_files)).drop_duplicates().reset_index(drop=True)
## Convert sampling description column
df_method=pd.concat(map (lambda desc:pd.DataFrame(ast.literal_eval('{"'+(' '.join(list(filter(lambda x:':' in x,desc.split(' '))))).replace(' ','","').replace(':','":"')+'"}') ,index=[0]).assign(Sampling_description=desc)  if str(desc)!='nan' else pd.DataFrame({'Sampling_description':desc},index=[0]),df.Sampling_description)).reset_index(drop=True)
df=pd.concat([df,df_method.drop(columns='Sampling_description')],axis=1)
# Apply depth and artefacts filter, as well as category filter if non-null
df_depthselection = (df.drop_duplicates(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'])[['Project_ID', 'Sample', 'Depth_min', 'Depth_max']]).astype(dict(zip(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'],[str]*4))).groupby(['Project_ID', 'Sample', 'Depth_min', 'Depth_max']).apply(lambda x: pd.Series({'Depth_selected': ((x.Depth_min.astype(float) < 200) & (x.Depth_max.astype(float) < 250)).values[0]})).reset_index()
df_depthselection.Depth_selected.value_counts(normalize=True)
df=pd.merge(df.astype(dict(zip(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'],[str]*4))),df_depthselection,how='left',on=['Project_ID', 'Sample', 'Depth_min', 'Depth_max']).assign(Depth_range=df[['Depth_min','Depth_max']].astype(dict(zip(['Depth_min','Depth_max'],2*[int]))).astype(dict(zip(['Depth_min','Depth_max'],2*[str]))).agg('-'.join, axis=1))
df.groupby(['Instrument']).apply(lambda x:pd.Series({'nb_samples':len(x.Sample.unique()),'Date_min':x.Sampling_datetime.dt.strftime('%Y %m').min(),'Date_max':x.Sampling_datetime.dt.strftime('%Y %m').max()})).reset_index()

# Group net depth range in 50m bins
if all(df.Instrument.isin(['Zooscan','Scanner'])):
    df['Depth_range']=df.Depth_range.apply(lambda depth:'-'.join([str(int(50 * round(float(depth.split('-')[0]) / 50))),str(int(50 * round(float(depth.split('-')[1]) / 50)))]))
if all(df.Instrument.isin(['IFCB'])):
    df['Depth_range']=df.Depth_range.apply(lambda depth:'-'.join([str(int(10 * round(float(depth.split('-')[0]) / 10))),str(int(10 * round(float(depth.split('-')[1]) / 10)))]))
if all(df.Instrument.isin(['UVP'])):
    df['Depth_min']=df.astype({'Depth_min':float}).groupby(['Instrument','Project_ID','Cruise','Sample','Latitude','Longitude','Sampling_datetime','Sampling_description'])['Depth_min'].transform(min)
    df['Depth_max']=df.astype({'Depth_max':float}).groupby(['Instrument','Project_ID','Cruise','Sample','Latitude','Longitude','Sampling_datetime','Sampling_description'])['Depth_max'].transform(max)
    df=df.drop(columns=['Depth_selected','Depth_range']).drop_duplicates().reset_index(drop=True)
    df=df.assign(Depth_range=df[['Depth_min','Depth_max']].astype(dict(zip(['Depth_min','Depth_max'],2*[float]))).astype(dict(zip(['Depth_min','Depth_max'],2*[int]))).astype(dict(zip(['Depth_min','Depth_max'],2*[str]))).agg('-'.join, axis=1),Depth_selected=df.groupby(['Instrument','Project_ID','Cruise','Sample','Latitude','Longitude','Sampling_datetime','Sampling_description']).apply(lambda x:  ((x.Depth_min.astype(float) < 200))).reset_index(['Instrument','Project_ID','Cruise','Sample','Latitude','Longitude','Sampling_datetime','Sampling_description'], drop=True))
    df['Depth_range']=df.Depth_range.astype(str).apply(lambda depth:'-'.join([str(int(100 * round(float(depth.split('-')[0]) / 100))),str(int(100 * round(float(depth.split('-')[1]) / 100)))]))
    df_depthselection = (df.drop_duplicates(['Project_ID', 'Sample', 'Depth_range'])[['Project_ID', 'Sample', 'Depth_range']]).astype(dict(zip(['Project_ID', 'Sample', 'Depth_range'],[str]*3))).reset_index(drop=True).groupby(['Project_ID', 'Sample', 'Depth_range']).apply(lambda x: pd.Series({'Depth_selected': ((float(x.Depth_range.unique()[0].split('-')[0]) < 200))})).reset_index()
    df_depthselection.Depth_selected.value_counts(normalize=True)
    df=pd.merge(df.astype(dict(zip(['Project_ID', 'Sample', 'Depth_range'],[str]*3))).drop(columns='Depth_selected'),df_depthselection,how='left',on=['Project_ID', 'Sample', 'Depth_range'])
    df=pd.merge(df.astype({'Project_ID':str}),df_list.astype({'Project_ID':str}).loc[df_list.Portal=='ecotaxa',['Project_ID','Instrument']],how='left',on='Project_ID')
df['Depth_range']=pd.Categorical(df['Depth_range'],categories=natsorted(df['Depth_range'].unique()))
df_summary=df.groupby(['Depth_selected','Depth_range']).apply(lambda x: pd.Series({'Sample_count':len(x.Sample)})).reset_index()

df_summary['Sample_percentage']=df_summary.groupby(['Depth_selected'])['Sample_count'].transform(lambda x:np.round(100*(x/np.nansum(x)),1))
df_summary['Sample_percentage']=np.where((df_summary.Sample_percentage<1) | (df_summary.Sample_percentage.isna()),'',df_summary.Sample_percentage.astype(str)+'%')
(100*(df_summary.groupby(['Depth_selected'])['Sample_count'].sum()/df_summary['Sample_count'].sum())).round(1)
plot=(ggplot(data=df_summary)+
  geom_col( aes(x="Depth_range",y='Sample_count',fill='Depth_selected'),width=1.5, position=position_dodge(width=0.9),color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,100),alpha=1) +
geom_text(df_summary,aes(x='Depth_range', y='Sample_count',label='Sample_percentage', group=1), position=position_dodge(width=0.9), size=8,va='bottom') +
      labs(x=r'Depth range (m)', y=r'# samples') +scale_fill_manual(values={True:'#{:02x}{:02x}{:02x}'.format(0, 0 , 0),False:'#{:02x}{:02x}{:02x}'.format(236, 236 , 236)},drop=True)+
theme(axis_ticks_direction="inout", legend_direction='horizontal', legend_position='top',
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=10),
                      legend_text=element_text(family="serif", size=10),
                      axis_title=element_text(family="serif", size=10),
                      axis_text_x=element_text(family="serif", size=8, rotation=90,ha='center'),
                      axis_text_y=element_text(family="serif", size=10, rotation=90),
                      plot_background=element_rect(fill='white'),strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)

plot[0].set_size_inches(9.2,4)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Sampling_depth_PSSdb.svg'.format(str(Path.home())), limitsize=False, dpi=600)

df_summary=df[(df.Depth_selected==True) &(df.Sampling_description.isna()==False)][['Cruise','Pmta_voltage','Pmta_threshold','Pmtb_voltage','Pmtb_threshold']].drop_duplicates().melt(id_vars=['Cruise'])
df_summary=df_summary.astype({'value':float}).groupby(['Cruise','variable']).agg({'value':['mean','std'] }).reset_index()
df_summary.columns=['Cruise','variable','mean','sd']

plot=(ggplot(data=df_summary.astype(dict(zip(['mean','sd'],2*[float]))))+
geom_pointrange(aes(x="variable", y='mean',ymax='mean+sd',ymin='mean-sd', color='Cruise', group='Cruise'), position=position_dodge(width=0.9), alpha=1) +
labs(x=r'', y=r'Instrument settings (mV)') +scale_color_manual(values={'EXPORTS':'#{:02x}{:02x}{:02x}'.format(145,111,111), 'NESLTER_broadscale':'#{:02x}{:02x}{:02x}{:02x}'.format(145,73,111,255), 'NESLTER_transect':'#{:02x}{:02x}{:02x}{:02x}'.format(145,73,111,200), 'SPIROPA':'#{:02x}{:02x}{:02x}{:02x}'.format(145,73,111,155), 'san-francisco-pier-17':'#{:02x}{:02x}{:02x}{:02x}'.format(145,170,186,255), 'santa-cruz-municipal-wharf':'#{:02x}{:02x}{:02x}{:02x}'.format(145,170,186,155)},drop=True)+
theme(axis_ticks_direction="inout", legend_direction='horizontal', legend_position='top',
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=10),
                      legend_text=element_text(family="serif", size=10),
                      axis_title=element_text(family="serif", size=10),
                      axis_text_x=element_text(family="serif", size=10),
                      axis_text_y=element_text(family="serif", size=10, rotation=90),
                      plot_background=element_rect(fill='white'),strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)

plot[0].set_size_inches(3.2,4)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Sampling_acq_PSSdb.svg'.format(str(Path.home())), limitsize=False, dpi=600)

plot=(ggplot(data=df[(df.Depth_selected==True) & (df.Net_mesh.isna()==False)])+
  facet_wrap('~Instrument', nrow=1)+
  geom_histogram( aes(x="Net_mesh"),color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1) +
geom_text(aes(x='Net_mesh', y=after_stat('count'),label=after_stat('prop*100'), group=1),stat='count', position=position_dodge(width=0.9), size=8,va='bottom', format_string='{:.1f}%') +
      labs(x=r'Net mesh (${\mu}$m)', y=r'# samples') +
theme(axis_ticks_direction="inout", legend_direction='horizontal', legend_position='top',
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=10),
                      legend_text=element_text(family="serif", size=10),
                      axis_title=element_text(family="serif", size=10),
                      axis_text_x=element_text(family="serif", size=10),
                      axis_text_y=element_text(family="serif", size=10, rotation=90),
                      plot_background=element_rect(fill='white'),strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)

plot[0].set_size_inches(4.5,4)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Sampling_acq_PSSdb.svg'.format(str(Path.home())), limitsize=False, dpi=600)

plot=(ggplot(data=df[(df.Depth_selected==True)])+
  facet_wrap('~Instrument', nrow=1)+
  geom_histogram( aes(x="Profiling_mode",group='Instrument'),color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1) +
geom_text(aes(x='Profiling_mode', y=after_stat('count'),label=after_stat('prop*100'), group=1),stat='count', position=position_dodge(width=0.9), size=8,va='bottom', format_string='{:.1f}%') +
      labs(x=r'Profiling_mode', y=r'# samples') +
theme(axis_ticks_direction="inout", legend_direction='horizontal', legend_position='top',
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=10),
                      legend_text=element_text(family="serif", size=10),
                      axis_title=element_text(family="serif", size=10),
                      axis_text_x=element_text(family="serif", size=10),
                      axis_text_y=element_text(family="serif", size=10, rotation=90),
                      plot_background=element_rect(fill='white'),strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)

plot[0].set_size_inches(3,4)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Sampling_acq_PSSdb.svg'.format(str(Path.home())), limitsize=False, dpi=600)




# extra figure for intercept values:


pal = {"IFCB":"tab:green",
           "UVP":"tab:gray",
           "Scanner":"tab:blue"}
ax=sns.barplot(data=df, x='Instrument', y='intercept_NB_mean', palette=pal, errorbar = 'sd', order= ['IFCB', 'Scanner', 'UVP'], width = 0.5)
sns.despine(top = True, right = True)
ax.set(xlabel = '', ylabel = 'Intercept \n' + r'log($\mu$m$^{3}$ dm$^{-3}$ $\mu$m$^{-3}$)')
plt.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/barplot_intercetpts.png'.format(str(Path.home())),  dpi=600)
plt.close()


# Fig 4. Intercept, slope, R2 global map

world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)
#df_summary=df.groupby(['Instrument','latitude','longitude']).agg({'N':'sum'}).reset_index().rename(columns={'N':'count'})
plot=(ggplot(data=df)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="longitude",y="latitude",fill='slope_NB_mean'),shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1, size=1) +
  coord_cartesian(expand = 0)+
 geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), fill='black', color='black') +
 labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)', color=r'Mean slope') +
 scale_fill_cmap(cmap_name='PuRd',limits=[-0.5, -2])+
#scale_fill_gradientn(trans='log10',colors=palette)+
theme(axis_ticks_direction="inout",
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=20),
                      legend_text=element_text(family="serif", size=15),
                      axis_title=element_text(family="serif", size=15),
                      axis_text_x=element_text(family="serif", size=15),
                      axis_text_y=element_text(family="serif", size=15, rotation=90),
                      plot_background=element_rect(fill='white'),strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)

plot[0].set_size_inches(13,3.1)
#plt.tight_layout()
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_4.1_slopes_NBSS.pdf'.format(str(Path.home())), dpi=600)
plt.close()

# Intercept: untransformed
df['intercept_unt'] = df['intercept_NB_mean'].apply(lambda x: math.exp(x))
Intercept_unt_summary=df.groupby(['Instrument']).apply(lambda x: pd.Series({'intercept_mean_unt': x.intercept_unt.mean(),
                                                                            'intercept_mean': x.intercept_NB_mean.mean(),
                                                                            'intercept_sd_unt': x.intercept_unt.std(),
                                                                             'intercept_sd': x.intercept_NB_mean.std()}))
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)
plot=(ggplot(data=df)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="longitude",y="latitude",fill='intercept_NB_mean'),shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1, size=1) +
  coord_cartesian(expand = 0)+
 geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), fill='black', color='black') +
 labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)', color=r'Mean Intercept') +
 scale_fill_cmap(cmap_name='PuRd', limits=[5,25]) +
 #scale_fill_gradientn(trans='log10',colors=palette)+
theme(axis_ticks_direction="inout",
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=20),
                      legend_text=element_text(family="serif", size=15),
                      axis_title=element_text(family="serif", size=15),
                      axis_text_x=element_text(family="serif", size=15),
                      axis_text_y=element_text(family="serif", size=15, rotation=90),
                      plot_background=element_rect(fill='white'),strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)

plot[0].set_size_inches(13,3.1)
plt.tight_layout()
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_4.2_intercepts_NBSS.pdf'.format(str(Path.home())), dpi=600)
plt.close()

#R2

world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)
#df_summary=df.groupby(['Instrument','latitude','longitude']).agg({'N':'sum'}).reset_index().rename(columns={'N':'count'})
plot=(ggplot(data=df)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="longitude",y="latitude",fill='r2_mean'),shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1, size=1) +
  coord_cartesian(expand = 0)+
 geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), fill='black', color='black') +
 labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)', color=r'Mean slope') +
 scale_fill_cmap(cmap_name='PuRd',limits=[0.85,0.99]) +
 #scale_fill_gradientn(trans='log10',colors=palette)+
theme(axis_ticks_direction="inout",
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=10),
                      legend_text=element_text(family="serif", size=15),
                      axis_title=element_text(family="serif", size=15),
                      axis_text_x=element_text(family="serif", size=15),
                      axis_text_y=element_text(family="serif", size=15, rotation=90),
                      plot_background=element_rect(fill='white'),strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)

plot[0].set_size_inches(13,3.1)
plt.tight_layout()
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_4.3_regress_coeff_NBSS.pdf'.format(str(Path.home())), dpi=600)
plt.close()

#Figure 5. Parameters by latitude
df['latitude_cat'] = pd.cut(df['latitude'], bins=list(np.arange(-90,100, 1))).apply(lambda x: np.round(x.mid))
graph_data_lat = df.groupby(['Instrument', 'latitude_cat']).apply(lambda x: pd.Series({'slope_mean':x.slope_NB_mean.mean(),
                                                                                      'slope_std':x.slope_NB_mean.std(),
                                                                                      'intercept_mean': x.intercept_NB_mean.mean(),
                                                                                      'intercept_std': x.intercept_NB_mean.std(),
                                                                                      'r2_mean': x.r2_NB_mean.mean(),
                                                                                      'r2_std': x.r2_NB_mean.std()})).reset_index()

colors = {'IFCB':'green', 'UVP':'gray', 'Scanner':'blue'}
#slope
plot = (ggplot(data=graph_data_lat)+
        geom_point(aes(x='latitude_cat', y='slope_mean', color='Instrument'))+
        geom_line(aes(x='latitude_cat', y='slope_mean', color='Instrument'))+
        geom_errorbar(aes(x ='latitude_cat', ymin = 'slope_mean-slope_std', ymax = 'slope_mean+slope_std', color = 'Instrument'))+
        coord_flip()+
        scale_color_manual(values = colors)+
        scale_x_continuous(breaks = np.arange(-90, 100, 10))+
        labs(y=r'Mean slope ( dm$^{-3}$ $\mu$m$^{-3}$)', x=r'Latitude ($^{\circ}$N)')+
        theme(axis_ticks_direction="inout",
              panel_grid=element_blank(),
              panel_border = element_blank(),
              axis_line = element_line(colour = "black"),
              panel_background=element_rect(fill='white'),
              #panel_border=element_rect(color='#222222'),
              legend_title=element_text(family="serif", size=10),
              legend_position=[0.5, 0.95],
              legend_text=element_text(family="serif", size=10),
              axis_title=element_text(family="serif", size=15),
              axis_text_x=element_text(family="serif", size=12),
              axis_text_y=element_text(family="serif", size=12, rotation=90),
              plot_background=element_rect(fill='white'), strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(5,8)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_5.1_slope_lat.pdf'.format(str(Path.home())), dpi=600)

#intercept
plot = (ggplot(data=graph_data_lat)+
        geom_point(aes(x='latitude_cat', y='intercept_mean', color='Instrument'))+
        geom_line(aes(x='latitude_cat', y='intercept_mean', color='Instrument'))+
        geom_errorbar(aes(x ='latitude_cat', ymin = 'intercept_mean-intercept_std', ymax = 'intercept_mean+intercept_std', color = 'Instrument'))+
        coord_flip()+
        scale_color_manual(values = colors)+
        scale_x_continuous(breaks = np.arange(-90, 100, 10))+
        labs(y=r'Mean Intercept ( $\mu$m$^{3}$ dm$^{-3}$ $\mu$m$^{-3}$)', x=r'Latitude ($^{\circ}$N)')+
        theme(axis_ticks_direction="inout",
              panel_grid=element_blank(),
              panel_border = element_blank(),
              axis_line = element_line(colour = "black"),
              panel_background=element_rect(fill='white'),
              #panel_border=element_rect(color='#222222'),
              legend_title=element_text(family="serif", size=10),
              legend_position=[0.5, 0.95],
              legend_text=element_text(family="serif", size=10),
              axis_title=element_text(family="serif", size=15),
              axis_text_x=element_text(family="serif", size=12),
              axis_text_y=element_text(family="serif", size=12, rotation=90),
              plot_background=element_rect(fill='white'), strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(5,8)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_5.2_intercept_lat.pdf'.format(str(Path.home())), dpi=600)

#R2
plot = (ggplot(data=graph_data_lat)+
        geom_point(aes(x='latitude_cat', y='r2_mean', color='Instrument'))+
        geom_line(aes(x='latitude_cat', y='r2_mean', color='Instrument'))+
        geom_errorbar(aes(x ='latitude_cat', ymin = 'r2_mean-r2_std', ymax = 'r2_mean+r2_std', color = 'Instrument'))+
        coord_flip()+
        scale_color_manual(values = colors)+
        scale_x_continuous(breaks = np.arange(-90, 100, 10))+
        labs(y=r'Mean R$^{2}$', x=r'Latitude ($^{\circ}$N)')+
        theme(axis_ticks_direction="inout",
              panel_grid=element_blank(),
              panel_border = element_blank(),
              axis_line = element_line(colour = "black"),
              panel_background=element_rect(fill='white'),
              #panel_border=element_rect(color='#222222'),
              legend_title=element_text(family="serif", size=10),
              legend_position=[0.5, 0.95],
              legend_text=element_text(family="serif", size=10),
              axis_title=element_text(family="serif", size=15),
              axis_text_x=element_text(family="serif", size=12),
              axis_text_y=element_text(family="serif", size=12, rotation=90),
              plot_background=element_rect(fill='white'), strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(5,8)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_5.3_R2_lat.pdf'.format(str(Path.home())), dpi=600)


##Figure 6. Climatology of parameters by ocean basin
# reasign variables
month_names = {1:'Jan', 2: 'Feb', 3:'Mar', 4:'Apr', 5:'May', 6: 'Jun', 7:'Jul', 8: 'Aug', 9:'Sep', 10:'Oct', 11: 'Nov', 12:'Dec'}
ocean_names = {'o00':'coastal', 'o07':'North_Pacific', 'o01': 'Arctic_Ocean', 'o22':'Mediterranean_Sea', 'o24':'Red_Sea', 'o05':'Indian_Ocean', 'o03':'South_Atlantic',
               'o04': 'Southern_Ocean', 'o06':'South_Pacific', 'o02': 'North_Atlantic', 'o21':'Baltic_Sea'}
df = df.replace({'month':month_names, 'ocean':ocean_names})
insufficient_data = ['Arctic_Ocean', 'Red_Sea', 'South_Atlantic', 'Southern_Ocean', 'Baltic_Sea', 'coastal']
df_clim = df[~df['ocean'].isin(insufficient_data)].reset_index()
pal = {"IFCB":"tab:green",
           "UVP":"tab:gray",
           "Scanner":"tab:blue"}


# Slope
fig, axs= plt.subplots(nrows=1, ncols=5, sharey=False, figsize=(30, 4))
ocean= list(df_clim['ocean'].unique())
for ocean, ax in zip(ocean, axs.ravel()):
    df_subset=df_clim[df_clim['ocean']== ocean]
    sns.pointplot(ax =ax, data=df_subset, x='month', y='slope_NB_mean',hue='Instrument', palette=pal, dodge= 0.3, errorbar = 'ci', order= ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    sns.despine(top = True, right = True)
    ax.set_ylim(-2, -0.4)
    ax.set_xticklabels(['Jan', '', 'Mar','', 'May','', 'Jul','', 'Sep','', 'Nov',''], fontsize = 15)
    ax.set_yticklabels([-2,-1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4] ,fontsize=15)
    ax.set_title(ocean.replace('_', ' '),fontsize = 18).set_fontname(fontname)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[3:], labels=labels[3:])
    ax.set(xlabel= None, ylabel = None)
axs[0].set_ylabel(r'Mean slope ( dm$^{-3}$ $\mu$m$^{-3}$)', fontsize = 17).set_fontname(fontname)
fig.set_figheight(4)
fig.set_figwidth(30)
plt.tight_layout
for ax in axs.flatten():
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname(fontname) for label in labels]
plt.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_6.1_slope_clim_basins.pdf'.format(str(Path.home())), dpi=600)
plt.close()



# Intercept
fig, axs= plt.subplots(nrows=1, ncols=5, sharey=False, figsize=(30, 4))
ocean= list(df_clim['ocean'].unique())
for ocean, ax in zip(ocean, axs.ravel()):
    df_subset=df_clim[df_clim['ocean']== ocean]
    sns.pointplot(ax =ax, data=df_subset, x='month', y='intercept_NB_mean',hue='Instrument', palette=pal, dodge= 0.3, errorbar = 'ci', order= ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    sns.despine(top = True, right = True)
    ax.set_ylim(0, 35)
    ax.set_xticklabels(['Jan', '', 'Mar','', 'May','', 'Jul','', 'Sep','', 'Nov',''], fontsize = 15)
    ax.set_yticklabels([0, 5,10,15,20,25,30,35] ,fontsize=15)
    ax.set_title(ocean.replace('_', ' '),fontsize = 18).set_fontname(fontname)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[3:], labels=labels[3:])
    ax.set(xlabel= None, ylabel = None)
axs[0].set_ylabel(r'Mean Intercept ( $\mu$m$^{3}$ dm$^{-3}$ $\mu$m$^{-3}$)', fontsize = 17).set_fontname(fontname)
fig.set_figheight(4)
fig.set_figwidth(30)
plt.tight_layout
for ax in axs.flatten():
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname(fontname) for label in labels]
plt.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_6.2_intercept_clim_basins.pdf'.format(str(Path.home())), dpi=600)
plt.close()


# R2
fig, axs= plt.subplots(nrows=1, ncols=5, sharey=False, figsize=(30, 4))
ocean= list(df_clim['ocean'].unique())
for ocean, ax in zip(ocean, axs.ravel()):
    df_subset=df_clim[df_clim['ocean']== ocean]
    sns.pointplot(ax=ax, data=df_subset, x='month', y='r2_NB_mean', hue='Instrument', palette=pal, dodge=0.3,errorbar='ci',order=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    sns.despine(top = True, right = True)
    ax.set_ylim(0.65, 1)
    ax.set_xticklabels(['Jan', '', 'Mar', '', 'May', '', 'Jul', '', 'Sep', '', 'Nov', ''], fontsize=15)
    ax.set_yticklabels([0.65,0.7,0.75, 0.8, 0.85, 0.9, 0.95, 1], fontsize=15)
    ax.set_title(ocean.replace('_', ' '), fontsize=18).set_fontname(fontname)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[3:], labels=labels[3:])
    ax.set(xlabel=None, ylabel=None)
axs[0].set_ylabel(r'Mean R$^{2}$', fontsize = 17).set_fontname(fontname)
fig.set_figheight(4)
fig.set_figwidth(30)
plt.tight_layout
for ax in axs.flatten():
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname(fontname) for label in labels]
plt.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_6.3_r2_clim_basins.pdf'.format(str(Path.home())), dpi=600)
plt.close()


## Figures of Sensitivity analysis (Supp figure 3)
path_to_datafile=(Path(cfg['git_dir']).expanduser()/ cfg['dataset_subdir']) / 'NBSS_data' / 'Sensitivity_analysis'
path_files=list(path_to_datafile.rglob('*_1b_*.csv'))
path_files_1a=list(path_to_datafile.rglob('*_1a_*.csv'))
df=pd.concat(map(lambda path: pd.read_table(path,sep=',').assign(Instrument=path.name.split('_')[0], Biovol_metric = path.name.split('_')[2]),path_files)).drop_duplicates().reset_index(drop=True)
colors = {'Biovolume-ellipsoid': 'teal', 'Biovolume-area': 'black', 'Biovolume-distance-map': 'slateblue'}
plot = (ggplot(data=df)+
        geom_violin(aes(x='Instrument', y='NBSS_slope_mean', color='Biovol_metric'),  position = position_dodge(width=1))+
        geom_point(aes(x='Instrument', y='NBSS_slope_mean', color='Biovol_metric'), position = position_dodge(width=1),size = 1, alpha=0.4, shape = 'o')+
        stat_summary(aes(x='Instrument', y='NBSS_slope_mean', color='Biovol_metric'),geom='point', fun_y=np.mean, shape=0, size = 5,  position = position_dodge(width=1))+
        labs(y=r'Mean slope ( dm$^{-3}$ $\mu$m$^{-3}$)', x='')+
        scale_color_manual(values = colors)+
        scale_fill_manual(values=['#00000000', '#00000000', '#00000000']) +
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(3,3)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Supp_fig_3.1_Biovol_sensitivity.pdf'.format(str(Path.home())), dpi=600)

plot = (ggplot(data=df)+
        geom_violin(aes(x='Instrument', y='NBSS_intercept_mean', color='Biovol_metric'),  position = position_dodge(width=1))+
        geom_point(aes(x='Instrument', y='NBSS_intercept_mean', color='Biovol_metric'), position = position_dodge(width=1),size = 1, alpha=0.4, shape = 'o')+
        stat_summary(aes(x='Instrument', y='NBSS_intercept_mean', color='Biovol_metric'),geom='point', fun_y=np.mean, shape=0, size = 5,  position = position_dodge(width=1))+
        labs(y=r'Mean Intercept ( $\mu$m$^{3}$ dm$^{-3}$ $\mu$m$^{-3}$)', x='')+
        scale_color_manual(values = colors)+
        scale_fill_manual(values=['#00000000', '#00000000', '#00000000']) +
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(3,3)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Supp_fig_3.2_Biovol_sensitivity.pdf'.format(str(Path.home())), dpi=600)

plot = (ggplot(data=df)+
        geom_violin(aes(x='Instrument', y='NBSS_r2_mean', color='Biovol_metric'),  position = position_dodge(width=1))+
        geom_point(aes(x='Instrument', y='NBSS_r2_mean', color='Biovol_metric'), position = position_dodge(width=1),size = 1, alpha=0.4, shape = 'o')+
        stat_summary(aes(x='Instrument', y='NBSS_r2_mean', color='Biovol_metric'),geom='point', fun_y=np.mean, shape=0, size = 5,  position = position_dodge(width=1))+
        labs(y=r'Mean R$^{2}$', x='')+
        scale_color_manual(values = colors)+
        scale_fill_manual(values=['#00000000', '#00000000', '#00000000']) +
        theme_paper).draw(show=False, return_ggplot=True)
plot[0].set_size_inches(3,3)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Supp_fig_3.3_Biovol_sensitivity.pdf'.format(str(Path.home())), dpi=600)
