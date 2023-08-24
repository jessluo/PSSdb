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
palette=list((sns.color_palette("PuRd",4).as_hex()))
fontname = 'serif'

## Workflow starts here:
# Generate stats for table 1:
path_to_project_list=Path(cfg['git_dir']).expanduser()/ cfg['proj_list']
df_list=pd.concat(map(lambda x:pd.read_excel(path_to_project_list,sheet_name=x).assign(Portal=x),['ecotaxa','ecopart','ifcb']))

path_to_datafile=(Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir']).expanduser()
path_standardizer=list(Path(path_to_datafile.parent).glob(['project_Zooscan*', 'project_Other*']))
path_standardizer=list(Path(path_to_datafile.parent).glob('project_UVP_*'))
path_standardizer=list(Path(path_to_datafile.parent).glob('project_IFCB_*'))

path_files=list(path_to_datafile.rglob('standardized_*.csv'))
df_standardizer=pd.concat(map(lambda path: pd.read_excel(path),path_standardizer)).reset_index(drop=True)
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
    df['Depth_range']=df.Depth_range.apply(lambda depth:'-'.join([str(int(50 * round(float(depth.split('-')[0]) / 50))),str(int(50 * round(float(depth.split('-')[1]) / 50)))]))

df['Depth_range']=pd.Categorical(df['Depth_range'],categories=natsorted(df['Depth_range'].unique()))
df_summary=df.groupby(['Depth_selected','Depth_range']).apply(lambda x: pd.Series({'Sample_count':len(x.Sample)})).reset_index()

df_summary['Sample_percentage']=df_summary.groupby(['Depth_selected'])['Sample_count'].transform(lambda x:np.round(100*(x/np.nansum(x)),1))
df_summary['Sample_percentage']=np.where((df_summary.Sample_percentage<1) | (df_summary.Sample_percentage.isna()),'',df_summary.Sample_percentage.astype(str)+'%')
(100*(df_summary.groupby(['Depth_selected'])['Sample_count'].sum()/df_summary['Sample_count'].sum())).round(1)
plot=(ggplot(data=df_summary)+
  geom_col( aes(x="Depth_range",y='Sample_count',fill='Depth_selected'),width=0.1, position=position_dodge(width=0.9),color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,100),alpha=1) +
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

df.groupby(['Instrument']).apply(lambda x:pd.Series({'nb_samples':len(x.Sample.unique()),'Date_min':x.Sampling_datetime.dt.strftime('%Y %m').min(),'Date_max':x.Sampling_datetime.dt.strftime('%Y %m').max()})).reset_index()


# Generate maps and heatmap of the spatio-temporal coverage for PSSdb products 1
path_to_datafile=(Path(cfg['git_dir']).expanduser()/ cfg['dataset_subdir']) / 'NBSS_data'
path_files=list(path_to_datafile.rglob('*_1b_*.csv'))
path_files_1a=list(path_to_datafile.rglob('*_1a_*.csv'))
df=pd.concat(map(lambda path: pd.read_table(path,sep=',').assign(Instrument=path.name[0:path.name.find('_')]),path_files)).drop_duplicates().reset_index(drop=True)
df_1a = pd.concat(map(lambda path: pd.read_table(path,sep=',').assign(Instrument=path.name[0:path.name.find('_')]),path_files_1a)).drop_duplicates().reset_index(drop=True)


# Fig 2. Spatial coverage
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)
df_summary=df.groupby(['Instrument','latitude','longitude']).agg({'N':'sum'}).reset_index().rename(columns={'N':'count'})
plot=(ggplot(data=df_summary)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="longitude",y="latitude",fill='count'),shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1) +
  coord_cartesian(expand = 0)+
 geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), fill='black', color='black') +
 labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)') +
 scale_fill_gradientn(trans='log10',colors=palette)+
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

plot[0].set_size_inches(6.5,2.8)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Spatial_coverage_PSSdb.svg'.format(str(Path.home())), limitsize=False, dpi=600)

# Fig 2. Temporal coverage

df_summary=df.groupby(['Instrument','month','year']).agg({'N':'sum'}).reset_index().rename(columns={'N':'count'})
plot=(ggplot(data=df_summary)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="month",y="year",fill='count'),size=5,shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1) +
 labs(x=r'Month', y=r'Year)') +
 scale_fill_gradientn(trans='log10',colors=palette)+scale_y_continuous(breaks=[1,4,7,10],labels=['Jan','Apr','Jul','Oct'])+
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

plot[0].set_size_inches(6.5,2.8)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_2_Temporal_coverage_PSSdb.svg'.format(str(Path.home())), limitsize=False, dpi=600)

#Fig. 3. Average NBSS
df_1a['log_biovolume_size_class'] = np.log10(df_1a['biovolume_size_class'])
df_1a['log_normalized_biovolume_mean'] = np.log10(df_1a['normalized_biovolume_mean'])
df_1a['log_diameter'] =df_1a['biovolume_size_class'].apply(lambda x: np.log10((x*6)/m.pi)**(1./3.))
df_1a['biovolume_size_class']= df_1a['biovolume_size_class'].astype(str)

graph_data = df_1a.groupby(['Instrument', 'log_biovolume_size_class']).apply(lambda x: pd.Series({'NB_mean':x.log_normalized_biovolume_mean.mean(), 'NB_std':x.log_normalized_biovolume_mean.std()})).reset_index()
colors = ['green', 'gray', 'blue']
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

for i, instrument in enumerate (['IFCB','UVP','Scanner']):
    data_plot = graph_data[graph_data['Instrument'] == instrument].reset_index()
    ax1.plot(data_plot.log_biovolume_size_class, data_plot.NB_mean, label=instrument, color= colors[i])
    ax1.fill_between(data_plot.log_biovolume_size_class, data_plot.NB_mean - data_plot.NB_std*2, data_plot.NB_mean + data_plot.NB_std*2, color=colors[i], alpha=0.2)

new_tick_values = np.array([10e1,  10e5,  10e9,  10e13])
new_tick_locations = np.array([1,  5,  9,  13])

def tick_function(X):
    V = np.round(np.log10(((X*6)/m.pi)**(1./3.))).astype(int)
    return V

#sns.despine(top = True, right = True)
ax1.set_yticks([-3, -1, 1, 3, 5])
ax1.spines['right'].set_visible(False)

ax1.set_ylabel( r'$log_{10}$ Normalized Biovolume ($\mu$m$^{3}$ dm$^{-3}$ $\Delta\mu$m$^{-3}$)')
ax1.set_xlabel(r'$log_{10}$ Biovolume ($\mu$m$^{3}$)')

ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(tick_function(new_tick_values))
ax2.spines['right'].set_visible(False)
ax2.set_xlabel(r'$log_{10}$ ESD ($\mu$m)')


plt.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_3_NBSS_summary.pdf'.format(str(Path.home())),  dpi=600)
plt.close()

##finding weird values in IFCB
IFCB_df = df_1a[df_1a.Instrument == 'IFCB']
IFCB_df['log_biovolume_size_class'] = np.log10(IFCB_df['biovolume_size_class'])
IFCB_df['log_normalized_biovolume_mean'] = np.log10(IFCB_df['normalized_biovolume_mean'])

IFCB_df['Sample_ID'] = IFCB_df.year.astype(str) +'_'+ IFCB_df.month.astype(str) +'_'+ IFCB_df.latitude.astype(str) + '_'+ IFCB_df.longitude.astype(str)
min_IFCB = IFCB_df[IFCB_df.log_biovolume_size_class == IFCB_df.log_biovolume_size_class.min()].reset_index()
min_IFCB_df = IFCB_df[IFCB_df.Sample_ID == min_IFCB.Sample_ID[0]].reset_index() # this data is from NESLTER broadscale

min_IFCB_df = min_IFCB_df.iloc[0:15]
sns.lineplot(data=min_IFCB_df , x='log_biovolume_size_class', y='log_normalized_biovolume_mean')


# Fig 4. Intercept, slope, R2 global map

world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)
#df_summary=df.groupby(['Instrument','latitude','longitude']).agg({'N':'sum'}).reset_index().rename(columns={'N':'count'})
plot=(ggplot(data=df)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="longitude",y="latitude",fill='slope_mean'),shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1, size=1) +
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
df['intercept_unt'] = df['intercept_mean'].apply(lambda x: math.exp(x))
Intercept_unt_summary=df.groupby(['Instrument']).apply(lambda x: pd.Series({'intercept_mean_unt': x.intercept_unt.mean(),
                                                                            'intercept_mean': x.intercept_mean.mean(),
                                                                            'intercept_sd_unt': x.intercept_unt.std(),
                                                                             'intercept_sd': x.intercept_mean.std()}))
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)
plot=(ggplot(data=df)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="longitude",y="latitude",fill='intercept_mean'),shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1, size=1) +
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
df['latitude_cat'] = pd.cut(df['latitude'], bins=list(np.arange(-90,100, 1))).apply(lambda x: np.round(x.mid))#.astype(str)
graph_data_lat = df.groupby(['Instrument', 'latitude_cat']).apply(lambda x: pd.Series({'slope_mean':x.slope_mean.mean(),
                                                                                      'slope_std':x.slope_mean.std(),
                                                                                      'intercept_mean': x.intercept_mean.mean(),
                                                                                      'intercept_std': x.intercept_mean.std(),
                                                                                      'r2_mean': x.r2_mean.mean(),
                                                                                      'r2_std': x.r2_mean.std()})).reset_index()

colors = {'IFCB':'green', 'UVP':'gray', 'Zooscan':'blue'}
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
           "Zooscan":"tab:blue"}


# Slope
fig, axs= plt.subplots(nrows=1, ncols=5, sharey=False, figsize=(30, 4))
ocean= list(df_clim['ocean'].unique())
for ocean, ax in zip(ocean, axs.ravel()):
    df_subset=df_clim[df_clim['ocean']== ocean]
    sns.pointplot(ax =ax, data=df_subset, x='month', y='slope_mean',hue='Instrument', palette=pal, dodge= 0.3, errorbar = 'ci', order= ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
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
    sns.pointplot(ax =ax, data=df_subset, x='month', y='intercept_mean',hue='Instrument', palette=pal, dodge= 0.3, errorbar = 'ci', order= ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
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
    sns.pointplot(ax=ax, data=df_subset, x='month', y='r2_mean', hue='Instrument', palette=pal, dodge=0.3,errorbar='ci',order=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
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