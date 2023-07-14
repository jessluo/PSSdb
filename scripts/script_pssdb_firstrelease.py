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
path_to_config=Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config ,'r') as config_file:
    cfg= yaml.safe_load(config_file)
from scripts.funcs_standardize_projects import *

# Functions for plotting:
import seaborn as sns
palette=list((sns.color_palette("PuRd",4).as_hex()))


## Workflow starts here:
# Generate maps and heatmap of the spatio-temporal coverage for PSSdb products 1
path_to_datafile=(Path(cfg['git_dir']).expanduser()/ cfg['dataset_subdir']) / 'NBSS_data'
path_files=list(path_to_datafile.rglob('*_1b_*.csv'))
df=pd.concat(map(lambda path: pd.read_table(path,sep=',').assign(Instrument=path.name[0:path.name.find('_')]),path_files)).drop_duplicates().reset_index(drop=True)

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

# Temporal coverage

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
plot[0].savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/Temporal_coverage_PSSdb.svg'.format(str(Path.home())), limitsize=False, dpi=600)
