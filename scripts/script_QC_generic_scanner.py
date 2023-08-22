## Objective: This script aims at checking the presence of all size fractions of generic scanner Ecotaxa projects and ensure SS estimates are consistent with ZooScan datasets

## Python modules:
# Data handling
import numpy as np
import pandas as pd
import re
import datetime
import scipy.optimize
from constrained_linear_regression import ConstrainedLinearRegression
from natsort import natsorted
from itertools import chain
from pathlib import Path
import matplotlib
#matplotlib.use('WebAgg')
# function to calculate NBSS from standardized files:
try:
    from script_NBSS_datacheck import *
except:
    from scripts.script_NBSS_datacheck import *

path_to_git=Path('~/GIT/PSSdb').expanduser()
path_to_data=path_to_git / 'raw'
path_to_project_list=path_to_data/ "project_list_all.xlsx"
path_to_standard_files=Path(path_to_data/ 'raw_standardized').expanduser()

standard_files=list(Path(path_to_standard_files / 'ecotaxa' / 'Other').rglob('standardized*'))
df_standardized=pd.concat(map(lambda path:pd.read_table(path,sep=",",keep_date_col=True,parse_dates={'Sampling_datetime':['Sampling_date','Sampling_time']}),standard_files)).drop_duplicates().reset_index(drop=True)
df_standardized=df_standardized.assign(Biovolume=(1 / 6) * np.pi * ((2 * ((df_standardized.Area / np.pi) ** 0.5)) ** 3),ECD=2 * ((df_standardized.Area / np.pi) ** 0.5))
# Apply depth and artefacts filter, as well as category filter if non-null
df_depthselection = (df_standardized.drop_duplicates(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'])[
        ['Project_ID', 'Sample', 'Depth_min', 'Depth_max']]).astype(dict(zip(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'],[str]*4))).groupby(
        ['Project_ID', 'Sample', 'Depth_min', 'Depth_max']).apply(lambda x: pd.Series({'Depth_selected': ((x.Depth_min.astype(float) < 200) & (x.Depth_max.astype(float) < 250)).values[0]})).reset_index()
df_standardized=pd.merge(df_standardized.astype(dict(zip(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'],[str]*4))),df_depthselection,how='left',on=['Project_ID', 'Sample', 'Depth_min', 'Depth_max'])
df_subset = df_standardized[(df_standardized.Depth_selected == True) & (df_standardized.Sampling_type.str.lower().isin(['test', 'exp', 'junk', 'culture'])==False) & (df_standardized.Category.str.lower().apply(lambda annotation: len( re.findall(r'bead|bubble|artefact|artifact|not-living', annotation)) == 0 if str(annotation) != 'nan' else True))]

df_subset_summary = (df_subset.drop_duplicates(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'])[['Project_ID', 'Sample', 'Depth_min', 'Depth_max']]).astype(dict(zip(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'], [str] * 4))).groupby(['Project_ID', 'Sample', 'Depth_min', 'Depth_max']).apply(lambda x: pd.Series({'Min_obs_depth': x.Depth_min.astype(float).min(), 'Max_obs_depth': x.Depth_max.astype(float).max()})).reset_index()
df_subset = pd.merge(df_subset.astype(dict(zip(['Project_ID', 'Sample', 'Depth_min', 'Depth_max'], [str] * 4))),df_subset_summary, how='left', on=['Project_ID', 'Sample', 'Depth_min', 'Depth_max'])
# Assign size bins and grouping index for each sample
df_subset = df_subset.assign(sizeClasses=pd.cut(df_subset['ECD'], bins, include_lowest=True))
df_subset = pd.merge(df_subset, df_bins, how='left', on=['sizeClasses'])
group = ['Instrument', 'Station', 'Sample', 'Sampling_datetime', 'Latitude', 'Longitude', 'Min_obs_depth', 'Max_obs_depth', 'Sampling_lower_size', 'Sampling_upper_size']
# Compute cumulative volume per sample/profile
df_subset_volume = df_subset.astype(dict(zip(group, [str] * len(group)))).groupby(group).apply(lambda x: pd.Series({'cumulative_volume': x[[ 'Sample','Depth_min','Depth_max','Volume_imaged']].drop_duplicates().Volume_imaged.astype(float).sum()})).reset_index()
df_subset = pd.merge(df_subset.astype(dict(zip(group, [str] * len(group)))), df_subset_volume, how='left', on=group)
# Assign a group index before groupby computation
group = list(np.delete(group, pd.Series(group).isin(['Sampling_lower_size', 'Sampling_upper_size'])))
df_subset = pd.merge(df_subset, df_subset.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left', on=group)
#  Compute NBSS without applying the thresholding: 1/ calculate fraction-specific NBSS 2/ Sum each fraction to obtain the total NBSS per net 3/ Integrate between 0-200m 4/ Average NBSS in subbins/bins
niter=0 # Change for uncertainty estimates
nbss_all = df_subset.groupby(group + list(df_bins.columns) + ['Group_index', 'Sampling_lower_size', 'Sampling_upper_size'],dropna=False).apply(lambda x: pd.Series({'size_class_pixel': (x.Pixel.unique()[0]) * 1e-03 * x.size_class_mid.unique()[0], 'NBSS_count': x.ROI_number.sum(),'sum_biovolume': sum(x.Biovolume * x.ROI_number), 'volume': x.cumulative_volume.unique()[0],'NBSS': sum(x.Biovolume * x.ROI_number) / x.cumulative_volume.unique()[0] / (x.range_size_bin.unique())[0],'NBSS_std': np.std(list(map(lambda count: np.random.normal(loc=np.random.choice(np.repeat(x.Biovolume, x.ROI_number), size=count, replace=True), scale=(1 / 3) * np.pi * (np.sqrt((x.Pixel.unique()[0]) * 1e-03 * x.size_class_mid.unique()[0]) ** 3) / ( (1e-03 * x.Pixel.unique()[0]) ** 3), size=None).sum() / x.cumulative_volume.unique()[0] / (x.range_size_bin.unique())[0], np.random.poisson(lam=sum(x.ROI_number), size=niter))))})).reset_index()
# Sum size fraction-specific NBSS
nbss_sum = nbss_all.groupby(group + list(df_bins.columns) + ['Group_index'] , dropna=False).apply(lambda x: pd.DataFrame(x[['NBSS_count', 'sum_biovolume', 'NBSS']].sum(axis='rows').to_dict(), index=[0]).assign(size_class_pixel=np.nanmean(x.size_class_pixel), NBSS_std=np.sqrt(sum(x.NBSS_std ** 2)),volume=1 if x.Instrument.isin(['Zooscan','Other']).any() else x.volume.unique()[0])).reset_index().sort_values(['Group_index', 'size_class_mid']).drop(columns='level_{}'.format(len(group + list(df_bins.columns) + ['Group_index'])))
nbss_sum =nbss_sum.assign( Latitude_subbin=pd.cut(nbss_sum.Latitude.astype(float),np.arange(-90,90,step=0.5),pd.IntervalIndex(pd.cut(np.arange(-90,90,step=0.5),np.arange(-90,90,step=0.5))).mid[1:].to_list()).astype(str),
                    Latitude_bin=pd.cut(nbss_sum.Latitude.astype(float),np.arange(-90,90,step=1),pd.IntervalIndex(pd.cut(np.arange(-90,90,step=1),np.arange(-90,90,step=1))).mid[1:].to_list()).astype(str),
                    Longitude_subbin=pd.cut(nbss_sum.Longitude.astype(float), np.arange(-180, 180, step=0.5),pd.IntervalIndex(pd.cut(np.arange(-180, 180, step=0.5), np.arange(-180, 180, step=0.5))).mid[ 1:].to_list()).astype(str),
                    Longitude_bin=pd.cut(nbss_sum.Longitude.astype(float), np.arange(-180, 180, step=1), pd.IntervalIndex(pd.cut(np.arange(-180, 180, step=1), np.arange(-180, 180, step=1))).mid[1:].to_list()).astype(str),
                    Date_subbin=pd.to_datetime(nbss_sum.Sampling_datetime).apply(lambda date:pd.to_datetime(date).strftime("%Y-{} %W").format(str(datetime.strptime(date.strftime('%Y %W-1'),'%Y %W-%w').month).zfill(2))),
                    Date_bin=pd.to_datetime(nbss_sum.Sampling_datetime).apply(lambda date:pd.to_datetime(date).strftime("%Y-{}").format(str(datetime.strptime(date.strftime('%Y %W-1'),'%Y %W-%w').month).zfill(2))))

# Integrate depth-specific NBSS over 0-200m. 1/ Drop sample in the grouping index and mix/max obs depth to combine epipelagic nets (0-100 and 100-200m) at the same location
group=['Instrument',  'Sampling_datetime','Date_subbin','Date_bin','Longitude','Longitude_subbin','Longitude_bin','Latitude','Latitude_subbin','Latitude_bin']
nbss_sum = pd.merge(nbss_sum.drop(columns='Group_index'), nbss_sum.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left', on=group)

nbss_sum_summary = nbss_sum.groupby(['Group_index']).apply(lambda x: pd.Series(
    {'Depth_range_min': min(x.drop_duplicates(['Min_obs_depth', 'Max_obs_depth']).Min_obs_depth.astype(float)),
     'Depth_range_max': max( x.drop_duplicates(['Min_obs_depth', 'Max_obs_depth']).Max_obs_depth.astype(float))})).reset_index()
nbss_sum = pd.merge(nbss_sum, nbss_sum_summary, how='left', on=['Group_index'])

nbss_int = nbss_sum.groupby(['Instrument','Group_index' ,'Sampling_datetime','Date_subbin','Date_bin','Longitude','Longitude_subbin','Longitude_bin','Latitude','Latitude_subbin','Latitude_bin'] + list( df_bins.columns)).apply(lambda x: pd.Series({'volume': x.volume.unique()[0], 'size_class_pixel': np.nanmean(x.size_class_pixel),
     'NBSS_count': sum(x.NBSS_count), 'NBSS_std': np.sqrt(sum(((( x.NBSS_std) ** 2) * ((np.where((x.Max_obs_depth.astype(float) - x.Min_obs_depth.astype(float)).values > 1,(x.Max_obs_depth.astype(float) - x.Min_obs_depth.astype(float)).values, 1) / np.where((x.Depth_range_max.astype(float).unique()[0] - x.Depth_range_min.astype(float).unique()[0]) > 1,(x.Depth_range_max.astype(float).unique()[0] - x.Depth_range_min.astype(float).unique()[0]), 1)) ** 2)))),
     'NBSS_int': sum(x.NBSS * np.where((x.Max_obs_depth.astype(float) - x.Min_obs_depth.astype(float)) > 1, (x.Max_obs_depth.astype(float) - x.Min_obs_depth.astype(float)), 1)) / np.where( (x.Depth_range_max.astype(float).unique()[0] - x.Depth_range_min.astype(float).unique()[0]) > 1,(x.Depth_range_max.astype(float).unique()[0] - x.Depth_range_min.astype(float).unique()[0]), 1),
      'NBSS_sum':sum(x.NBSS), 'Depth_range_min': x.Depth_range_min.unique()[0], 'Depth_range_max': x.Depth_range_max.unique()[0]})).reset_index()

# Average in subbins (Note that observations are cumulated, not averaged, for UVP and IFCb) and bins. Attention: on l.76, select either _sum or _int to generate NBSS
nbss_avg=nbss_int.groupby(['Instrument','Date_subbin','Date_bin','Longitude_subbin','Longitude_bin','Latitude_subbin','Latitude_bin','Depth_range_min','Depth_range_max','size_class_pixel'] + list( df_bins.columns)).apply(lambda x:pd.Series({'NBSS_count': np.nanmean(x.NBSS_count),'NBSS': np.nanmean(x.NBSS_sum)})).reset_index()
nbss_avg=nbss_avg.groupby(['Instrument','Date_bin','Longitude_bin','Latitude_bin','Depth_range_min','Depth_range_max','size_class_pixel'] + list( df_bins.columns)).apply(lambda x:pd.Series({'NBSS_count': np.nanmean(x.NBSS_count),'NBSS': np.nanmean(x.NBSS)})).reset_index()

# Apply thresholding to select unbiased NBSS
nbss_avg['selected']=nbss_avg.groupby(['Instrument','Date_bin','Longitude_bin','Latitude_bin','Depth_range_min','Depth_range_max'] ).apply(lambda x: pd.Series({'selected':(x.NBSS.index >= x.NBSS.index.to_list()[np.argmax(x.NBSS)]) & (x.NBSS.index <np.where(all(x.NBSS.drop(x.NBSS.index.to_list()[0:np.argmax(x.NBSS)]).isna()==False),x.index.max(),x.drop(x.NBSS.index.to_list()[0:np.argmax(x.NBSS)]).index[x.NBSS.drop(x.NBSS.index.to_list()[0:np.argmax(x.NBSS)]).isna()==False].max()))})).reset_index().explode('selected').reset_index(drop=True)['selected']

# Plot with non-generic scanner products
nbss_ref=pd.concat(map(lambda path:pd.read_csv(path,sep=',').assign(Instrument=(path.name)[0:(str(path.name).find('_'))]),list(Path('~/GIT/PSSdb/raw/NBSS_data/NBSS_ver_08_2023/NO_generic_scanner').expanduser().rglob('*_1a_*')))).reset_index(drop=True)
group=['Instrument','year','month','latitude','longitude']
nbss_ref= pd.merge(nbss_ref, nbss_ref.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left', on=group)
nbss_ref['biovolume_size_class']=np.round(nbss_ref['biovolume_size_class'],3)
nbss_summary=pd.merge(nbss_ref.groupby(['Instrument','biovolume_size_class']).apply(lambda x: pd.Series({'NB_5':np.nanquantile(x['normalized_biovolume_mean'],q=0.05),'NB_50':np.nanquantile(x['normalized_biovolume_mean'],q=0.5),'NB_95':np.nanquantile(x['normalized_biovolume_mean'],q=0.95)})).reset_index(),nbss_ref[['equivalent_circular_diameter_mean','biovolume_size_class']].round(4).drop_duplicates(),how='left',on=['biovolume_size_class'])
plot=(ggplot(data=nbss_ref)+
  geom_ribbon(nbss_summary,aes(x="equivalent_circular_diameter_mean",ymin="NB_5",ymax="NB_95",y='NB_50',fill='Instrument',group='Instrument'),color='#{:02x}{:02x}{:02x}{:02x}'.format(0,0,0,0),alpha=0.6) +
    geom_line(nbss_summary,aes(x="equivalent_circular_diameter_mean", y='NB_50', color='Instrument', fill='Instrument',group='Instrument'), alpha=1) +
    scale_fill_manual(values={'IFCB':'#{:02x}{:02x}{:02x}'.format(111, 145 , 111),'UVP':'#{:02x}{:02x}{:02x}'.format(147,167,172),'Scanner':'#{:02x}{:02x}{:02x}'.format(95,141,211)})+
 scale_color_manual(values={'IFCB': '#{:02x}{:02x}{:02x}'.format(111, 145, 111),'UVP': '#{:02x}{:02x}{:02x}'.format(147, 167, 172), 'Scanner': '#{:02x}{:02x}{:02x}'.format(95, 141, 211)}) +
stat_summary(nbss_avg[nbss_avg.selected==True], aes(x="size_class_mid", y="NBSS"), alpha=1, geom='line', color='black') +
stat_summary(nbss_avg[nbss_avg.selected==True], aes(x="size_class_mid", y="NBSS"), alpha=0.8, geom='ribbon', color='black') +
 labs(x=r'Equivalent circular diameter ($\mu$m)', y=r'Normalized Biovolume Size Spectrum ($\mu$m$^{3}$ dm$^{-3}$ $\mu$m$^{-3}$)') +
scale_y_log10(breaks=[10**np.arange(-5,7,step=2, dtype=np.float)][0],labels=['10$^{%s}$'% int(n) for n in np.arange(-5,7,step=2)])+
scale_x_log10(breaks=[size  for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))],labels= [size if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))])+
theme(axis_ticks_direction="inout", legend_direction='horizontal', legend_position='top',
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=10),
                      legend_text=element_text(family="serif", size=10),
                      axis_title=element_text(family="serif", size=10),
                      axis_text_x=element_text(family="serif", size=10),
                      axis_text_y=element_text(family="serif", size=10, rotation=90),
                      plot_background=element_rect(fill='white'),strip_background=element_rect(fill='white'))).draw(show=True, return_ggplot=True)
plot[0].savefig(fname='{}/GIT/PSSdb/figures/firstdata_paper/Generic_PSSdb_NBSS_int.png'.format(str(Path.home())), limitsize=False, dpi=600)


# Check size fractions
df_standardized=pd.concat(map(lambda path:pd.read_table(path,sep=",",usecols=['Project_ID','Sample','Latitude','Longitude','Sampling_lower_size','Sampling_upper_size']),standard_files)).drop_duplicates().reset_index(drop=True)
#df_standardized['Sample']=df_standardized.Sample.str.lower()
df_standardized.Sampling_lower_size.unique()
df_standardized.Sampling_lower_size=pd.Categorical(df_standardized.Sampling_lower_size,ordered=True,categories=[200,400,800])

df_check=df_standardized.groupby(['Sample','Latitude','Longitude']).apply(lambda x:pd.DataFrame(dict(zip(x.Sampling_lower_size.value_counts(sort=False).index,np.where(x.Sampling_lower_size.value_counts(sort=False)>0,[x.query('Sampling_lower_size=={}'.format(fraction)).Project_ID.astype(str).unique()[0] if len(x.query('Sampling_lower_size=={}'.format(fraction)).Project_ID.astype(str).unique()) else False for fraction in x.Sampling_lower_size.value_counts(sort=False).index],x.Sampling_lower_size.value_counts(sort=False)>0))),index=[0])).reset_index().drop(columns='level_3')
df_check['Complete']=np.where((df_check[df_standardized.Sampling_lower_size.cat.categories].astype(str)=='False').any(axis=1),False,True)
df_check[df_check.Complete==False].to_excel(Path('~/GIT/PSSdb/raw/QC_scanner_complete.xlsx'),index=False)
# Map of generic-scanner samples
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
world_polygon = pd.concat([pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon[0]).shape[0])}),pd.DataFrame(polygon[0], columns=['Longitude', 'Latitude'])], axis=1) if pd.DataFrame(polygon[0]).shape[1] > 1 else pd.concat([pd.DataFrame({'Country': np.repeat(country_polygon, pd.DataFrame(polygon).shape[0])}),pd.DataFrame(polygon, columns=['Longitude', 'Latitude'])], axis=1) for country, region in zip(world.name, world.geometry) for country_polygon, polygon in zip([str(country) + "_" + str(poly) for poly in np.arange(len(mapping(region)['coordinates'])).tolist()], mapping(region)['coordinates'])], axis=0)

plot=(ggplot(data=df_check)+
  geom_point( aes(x="Longitude",y="Latitude",fill='Complete'),shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(255, 255 , 255,0),alpha=1) +
  coord_cartesian(expand = 0)+
 geom_polygon(data=world_polygon, mapping=aes(x='Longitude', y='Latitude', group='Country'), fill='black', color='black') +
 labs(x=r'Longitude ($^{\circ}$E)', y=r'Latitude ($^{\circ}$N)') +
theme(axis_ticks_direction="inout", legend_direction='horizontal', legend_position='top',
                      panel_grid=element_blank(),
                      panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'),
                      legend_title=element_text(family="serif", size=10),
                      legend_text=element_text(family="serif", size=10),
                      axis_title=element_text(family="serif", size=10),
                      axis_text_x=element_text(family="serif", size=10),
                      axis_text_y=element_text(family="serif", size=10, rotation=90),
                      plot_background=element_rect(fill='white'),strip_background=element_rect(fill='white'))).draw(show=True, return_ggplot=True)


# Open flag files and save after flagging incomplete samples
df_sub=df_check[df_check.Complete==False]
path_to_flags=[path for path in list(Path(path_to_data / 'flags' / 'ecotaxa').glob("*_flags.csv")) if path.name.split("_")[1] in list(set([project for project in df_sub[df_standardized.Sampling_lower_size.cat.categories].values.flatten() if project!='False']))]

df_flags=pd.concat(map(lambda path: pd.read_csv(path).assign(Project_path=path),natsorted(path_to_flags))).reset_index(drop=True)
df_flags.loc[(df_flags.Sample.isin(df_sub.Sample.unique())) & (df_flags.Flag==0) & (df_flags.Overrule.astype(str)=='False'),'Overrule']=True
df_flags.loc[(df_flags.Sample.isin(df_sub.Sample.unique())) & (df_flags.Flag==1) & (df_flags.Overrule.astype(str)=='True'),'Overrule']=False
df_flags.groupby(['Project_path']).apply(lambda x:x.drop(columns=['Project_path']).to_csv(x.Project_path.unique()[0],index=False,sep=','))