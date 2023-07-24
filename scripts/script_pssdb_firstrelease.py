## Objective: This script generates all the summary statistics and plots used for PSSdb first-release data paper

## Modules:

# Modules for data and path handling:
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from pathlib import Path  # Handling of path object
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


## Workflow starts here:
# Generate maps and heatmap of the spatio-temporal coverage for PSSdb products 1
path_to_datafile=(Path(cfg['git_dir']).expanduser()/ cfg['dataset_subdir']) / 'NBSS_data'
path_files=list(path_to_datafile.rglob('*_1b_*.csv'))
path_files_1a=list(path_to_datafile.rglob('*_1a_*.csv'))
df=pd.concat(map(lambda path: pd.read_table(path,sep=',').assign(Instrument=path.name[0:path.name.find('_')]),path_files)).drop_duplicates().reset_index(drop=True)
df_1a = pd.concat(map(lambda path: pd.read_table(path,sep=',').assign(Instrument=path.name[0:path.name.find('_')]),path_files_1a)).drop_duplicates().reset_index(drop=True)
fontname = 'serif'

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

for i, instrument in enumerate (['IFCB','UVP','Zooscan']):
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
min_IFCB = IFCB_df[IFCB_df.log_biovolume_size_class == IFCB_df.log_biovolume_size_class.min()]
# Fig 4. Intercept, slope, R2

world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)
#df_summary=df.groupby(['Instrument','latitude','longitude']).agg({'N':'sum'}).reset_index().rename(columns={'N':'count'})
plot=(ggplot(data=df)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="longitude",y="latitude",fill='slope_mean'),shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1, size=1) +
  coord_cartesian(expand = 0)+
 geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), fill='black', color='black') +
 labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)', color=r'Mean slope') +
 scale_fill_cmap(limits=[-0.5, -2]) +
 #scale_fill_gradientn(trans='log10',colors=palette)+
theme(axis_ticks_direction="inout",
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=10),
                      legend_text=element_text(family="serif", size=10),
                      axis_title=element_text(family="serif", size=10),
                      axis_text_x=element_text(family="serif", size=10),
                      axis_text_y=element_text(family="serif", size=10, rotation=90),
                      plot_background=element_rect(fill='white'),strip_background=element_rect(fill='white'))).draw(show=False, return_ggplot=True)

plot[0].set_size_inches(13,3.1)
plt.tight_layout()
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_4.1_slopes_NBSS.pdf'.format(str(Path.home())), dpi=600)
plt.close()

# Intercept: untransformed
#df['intercept_unt'] = df['intercept_mean'].apply(lambda x: m.e**x)

world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)
#df_summary=df.groupby(['Instrument','latitude','longitude']).agg({'N':'sum'}).reset_index().rename(columns={'N':'count'})
plot=(ggplot(data=df)+
  facet_wrap('~Instrument', nrow=1)+
  geom_point( aes(x="longitude",y="latitude",fill='intercept_mean'),shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1, size=1) +
  coord_cartesian(expand = 0)+
 geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), fill='black', color='black') +
 labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)', color=r'Mean Intercept') +
 scale_fill_cmap(limits=[5,25]) +
 #scale_fill_gradientn(trans='log10',colors=palette)+
theme(axis_ticks_direction="inout",
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=10),
                      legend_text=element_text(family="serif", size=10),
                      axis_title=element_text(family="serif", size=10),
                      axis_text_x=element_text(family="serif", size=10),
                      axis_text_y=element_text(family="serif", size=10, rotation=90),
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
 scale_fill_cmap(limits=[0.85,0.99]) +
 #scale_fill_gradientn(trans='log10',colors=palette)+
theme(axis_ticks_direction="inout",
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=10),
                      legend_text=element_text(family="serif", size=10),
                      axis_title=element_text(family="serif", size=10),
                      axis_text_x=element_text(family="serif", size=10),
                      axis_text_y=element_text(family="serif", size=10, rotation=90),
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



# Slope
fig, axs= plt.subplots(nrows=1, ncols=5, sharey=False, figsize=(25, 4))
ocean= list(df_clim['ocean'].unique())
for ocean, ax in zip(ocean, axs.ravel()):
    df_subset=df_clim[df_clim['ocean']== ocean]
    sns.boxplot(ax =ax, data=df_subset, x='month', y='slope_mean',color='skyblue', showfliers = False, order= ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    sns.despine(top = True, right = True)
    ax.set_ylim(-2, -0.4)
    ax.set_xticklabels(['Jan', '', 'Mar','', 'May','', 'Jul','', 'Sep','', 'Nov',''])
    ax.set_title(ocean.replace('_', ' ')).set_fontname(fontname)
    ax.set(xlabel= None, ylabel = None)
axs[0].set_ylabel(r'Mean slope ( dm$^{-3}$ $\mu$m$^{-3}$)').set_fontname(fontname)
plt.tight_layout
for ax in axs.flatten():
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname(fontname) for label in labels]
plt.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_6.1_slope_clim_basins.pdf'.format(str(Path.home())), dpi=600)
plt.close()



# Intercept
fig, axs= plt.subplots(nrows=1, ncols=5, sharey=False, figsize=(25, 4))
ocean= list(df_clim['ocean'].unique())
for ocean, ax in zip(ocean, axs.ravel()):
    df_subset=df_clim[df_clim['ocean']== ocean]
    sns.boxplot(ax =ax, data=df_subset, x='month', y='intercept_mean',color='skyblue', showfliers = False, order= ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    sns.despine(top = True, right = True)
    ax.set_ylim(2, 50)
    ax.set_xticklabels(['Jan', '', 'Mar','', 'May','', 'Jul','', 'Sep','', 'Nov',''])
    ax.set_title(ocean.replace('_', ' ')).set_fontname(fontname)
    #ax.set_title('')
    ax.set(xlabel= None, ylabel = None)
axs[0].set_ylabel(r'Mean Intercept ( $\mu$m$^{3}$ dm$^{-3}$ $\mu$m$^{-3}$)').set_fontname(fontname)
plt.tight_layout
for ax in axs.flatten():
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname(fontname) for label in labels]

plt.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_6.2_intercept_clim_basins.pdf'.format(str(Path.home())), dpi=600)
plt.close()


# R2
fig, axs= plt.subplots(nrows=1, ncols=5, sharey=False, figsize=(25, 4))
ocean= list(df_clim['ocean'].unique())
for ocean, ax in zip(ocean, axs.ravel()):
    df_subset=df_clim[df_clim['ocean']== ocean]
    sns.boxplot(ax =ax, data=df_subset, x='month', y='r2_mean',color='skyblue', showfliers = False, order= ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    sns.despine(top = True, right = True)
    ax.set_ylim(0.65, 1)
    ax.set_xticklabels(['Jan', '', 'Mar','', 'May','', 'Jul','', 'Sep','', 'Nov',''])
    ax.set_title(ocean.replace('_', ' ')).set_fontname(fontname)
    #ax.set_title('')
    ax.set(xlabel= None, ylabel = None)
axs[0].set_ylabel(r'Mean R$^{2}$').set_fontname(fontname)
plt.tight_layout
for ax in axs.flatten():
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontname(fontname) for label in labels]

plt.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Fig_6.3_r2_clim_basins.pdf'.format(str(Path.home())), dpi=600)
plt.close()