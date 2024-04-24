## Objective: Generate a merged product that uses the method described in Soviadan et al.(in review) to select the maximum NBSS for a given sample

## Requirements: Computation of group-specific NBSS (Step 5)

## Python modules and path specifications:
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from pathlib import Path
import statistics as st
import os
import datetime, time
from natsort import natsorted
import shutil
from tqdm import tqdm
from functools import reduce
try:
    from funcs_NBS import *
    from funcs_standardize_projects import *
    from funcs_merged_products import *

except:
    from scripts.funcs_NBS import *
    from scripts.funcs_standardize_projects import *
    from scripts.funcs_merged_products import *

import seaborn as sns
def plot_unity(xdata, ydata, **kwargs):
    mn = min(xdata.min(), ydata.min())
    mx = max(xdata.max(), ydata.max())
    points = np.linspace(mn, mx, 100)
    plt.gca().plot(points, points, color='k', marker=None,
            linestyle='--', linewidth=1.0)
from scipy.stats import linregress
def r2(x, y, ax=None, **kws):
    ax = ax or plt.gca()
    data=pd.DataFrame({'x':x,'y':y}).dropna(subset=['x','y'])
    if len(data):
        slope, intercept, r_value, p_value, std_err = linregress(x=data.x, y=data.y)
        ax.annotate(f'$r^2 = {r_value ** 2:.2f}$\nEq: ${slope:.2f}x{intercept:+.2f}$',
                    xy=(.05, .95), xycoords=ax.transAxes, fontsize=8,
                    color='black', backgroundcolor='#FFFFFF99', ha='left', va='top')



## Workflow starts here:

# Load most recent NBSS taxa-specific products
path_to_products=natsorted(list((sorted(list(Path(path_to_git / cfg['dataset_subdir'] /cfg['product_subdir']).glob('NBSS_ver_*')),key = lambda fname:time.strptime(os.path.basename(fname),"NBSS_ver_%m_%Y"))[-1]).rglob('PFT/*-distribution_all_var_*'))) #rglob('PFT/*_1a_*')
grouping_factor=['date_bin','Station_location','midLatBin','midLonBin']#['ocean','year','month','latitude','longitude']
path_to_directory=Path(path_to_products[0].parent.parent / 'Merged')
path_to_directory.mkdir(exist_ok=True)
with tqdm(desc='Merging taxa-specific {} products', total=3, bar_format='{desc}{bar}', position=0, leave=True) as bar:
    for product in ['Size','Biomass','Weight']:
        bar.set_description("Merging taxa-specific {} products".format(product), refresh=True)
        # Select the associated files
        path_to_files=[path for path in path_to_products if product in  path.name]
        df_nbss=merge_products(path_to_products=path_to_files, grouping_factors=grouping_factor)
        df = df_nbss[df_nbss.n_instruments > 1]
        group = grouping_factor
        df = pd.merge(df, df.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename({'index': 'Group_station_index'}, axis='columns'), how='left', on=group)
        df =df.assign(Sample=df.midLatBin.astype(str) + '_' + df.midLonBin.astype(str) + "_" + df.date_bin.astype( str)).rename(columns={'midLonBin': 'Longitude', 'midLatBin': 'Latitude','Min_obs_depth':'min_depth','Max_obs_depth':'max_depth'}) #df.assign(Sample=df.latitude.astype(str) + '_' + df.longitude.astype(str) + "_" + df.year.astype( str) + df.month.astype(str).str.zfill(2)).rename(columns={'longitude': 'Longitude', 'latitude': 'Latitude'})
        ## Apply Soviadan et al. maximum selection method: https://doi.org/10.1101/2023.06.29.547051
        # Use class mid-point (biovolume, biomass, or weight) to assign the rest of the size classes information
        #data_bins =df_bins.drop(columns=['sizeClasses']).assign(sizeClasses=pd.cut(np.append(0,df_bins[dict_products_x[product]].to_numpy()),np.append(0,df_bins[dict_products_x[product]].to_numpy())).categories.values.astype(str))# df_bins.drop(columns=['sizeClasses']).assign(sizeClasses=pd.cut(np.append(0,df_bins[dict_products_x[product]].to_numpy()),np.append(0,df_bins[dict_products_x[product]].to_numpy())).categories.values.astype(str))
        df=pd.merge(df.drop(columns=[column for column in df.columns if column in df_bins.columns]).assign(sizeClasses=pd.cut(df[dict_products_x[product]],data_bins[dict_products_x[product]]).astype(str)),df_bins.drop(columns=['sizeClasses']).assign(sizeClasses=pd.cut(data_bins[dict_products_x[product]].to_numpy(),data_bins[dict_products_x[product]].to_numpy()).categories.values.astype(str)),how='left',on='sizeClasses') # data_bins
        # Use the maximum normalized abundance / biovolume / biomass / weight in each size class
        #x=df.loc[list(df.groupby(group+['PFT','sizeClasses']).groups.values())[3]]
        group=['Latitude', 'Longitude', 'date_bin','Group_station_index','n_instruments','merged_instruments']#['Latitude', 'Longitude', 'year', 'month','Group_station_index','n_instruments','merged_instruments']
        df.loc[(df.count_uncertainty>cfg['artefacts_threshold']) | (df.size_uncertainty>cfg['artefacts_threshold']),dict_products_y[product]]=np.nan
        df_selection=df.groupby(group+['PFT','sizeClasses']).apply(lambda x: pd.Series({**{'Validation_percentage':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].Validation_percentage.mean(),'min_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.min(),'max_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.max(),'selected_instrument':x.loc[(x[dict_products_y[product]]).idxmax(axis=0),'Instrument'].unique()[0]},**dict(zip(dict_products_y[product], (x[dict_products_y[product]]).max(axis=0)))}) if any((x[dict_products_y[product]].isna()==False).to_numpy()[0]) else None).reset_index()#df.groupby(group+['PFT','sizeClasses']).apply(lambda x: pd.Series({**{'min_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.min(),'max_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.max(),'selected_instrument':x.loc[(x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).idxmax(axis=0),'Instrument'].unique()[0]},**dict(zip(dict_products_y[product], ( x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).max(axis=0))),**dict(zip([variable.replace('_mean','_std') for variable in dict_products_y[product]],np.diag(x.loc[(x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).idxmax(axis=0),[variable.replace('_mean','_std') for variable in dict_products_y[product]]]))),**{'n':x.loc[(x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).idxmax(axis=0),'n'].unique()[0]}})).reset_index()#df.groupby(group+['PFT','sizeClasses']).apply(lambda x: pd.Series({**{'Validation_percentage':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].Validation_percentage.mean(),'min_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.min(),'max_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.max(),'selected_instrument':x.loc[(x[dict_products_y[product]]).idxmax(axis=0),'Instrument'].unique()[0]},**dict(zip(dict_products_y[product], (x[dict_products_y[product]]).max(axis=0))),**dict(zip([variable.replace('_mean','_std') for variable in dict_products_y[product]],np.diag(x.loc[(x[dict_products_y[product]]).idxmax(axis=0),[variable.replace('_mean','_std') for variable in dict_products_y[product]]]))),**{'n':x.loc[(x[dict_products_y[product]]).idxmax(axis=0),'n'].unique()[0]}})).reset_index()#df.groupby(group+['PFT','sizeClasses']).apply(lambda x: pd.Series({**{'min_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.min(),'max_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.max(),'selected_instrument':x.loc[(x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).idxmax(axis=0),'Instrument'].unique()[0]},**dict(zip(dict_products_y[product], ( x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).max(axis=0))),**dict(zip([variable.replace('_mean','_std') for variable in dict_products_y[product]],np.diag(x.loc[(x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).idxmax(axis=0),[variable.replace('_mean','_std') for variable in dict_products_y[product]]]))),**{'n':x.loc[(x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).idxmax(axis=0),'n'].unique()[0]}})).reset_index()
        df_selection=df_selection.dropna(subset='selected_instrument').reset_index(drop=True)
        # Plot the selection for each PFT
        plot = (ggplot(data=pd.merge(df_selection,df_bins,how='left',on='sizeClasses').sort_values(['Latitude', 'Longitude', 'Group_station_index','PFT','ECD_mid']).reset_index(drop=True).assign(latitude_bin=lambda x:pd.cut(x.Latitude.astype(float),bins=np.arange(-90,90,1),labels=np.arange(-90,90,1)[1:]).astype(float),PFT=lambda x: pd.Categorical(x.PFT.str.replace("_"," "),['Nanophytoplankton','Microphytoplankton','Mesophytoplankton','Ciliophora','Rhizaria','Crustaceans','Filter feeder gelatinous','Carnivorous gelatinous','Mollusca','Actinopterygii','Other zooplankton','Detritus']))) +
                facet_wrap('~PFT', nrow=2) +
                plotnine.stat_bin_2d( mapping=aes(x="ECD_mid", y="latitude_bin",group='latitude_bin', fill='selected_instrument'), bins=25,geom='point', shape='H', color='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50)).to_layer() +
                stat_summary(aes(x="ECD_mid", y="latitude_bin",group='latitude_bin', fill='selected_instrument'),geom='point', size=1.5, shape='H', color='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50), alpha=0.1) +
                scale_fill_manual(pal_instrument,limits=df_selection.selected_instrument.unique())+
                labs(fill='Selected instruments:',x=r'Equivalent circular diameter ($\mu$m)', y=r'Latitude ($\degree$N)') +
                scale_x_log10(breaks=[size  for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))],labels=lambda l: [int(size) if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in l])+
                theme_paper).draw(show=True)
        plot.set_size_inches(8.5, 3.8)
        plot.savefig(fname='{}/PSSdb_Merged-products_selection.svg'.format(str(path_to_directory)), dpi=300)
        df_instruments=df.assign(latitude_bin=lambda x: pd.cut(x.Latitude.astype(float), bins=np.arange(-90, 90, 1), labels=np.arange(-90, 90, 1)[1:]).astype(float),PFT=lambda x: pd.Categorical(x.PFT.str.replace("_", " "), ['Nanophytoplankton', 'Microphytoplankton', 'Mesophytoplankton', 'Ciliophora','Rhizaria', 'Crustaceans', 'Filter feeder gelatinous','Carnivorous gelatinous', 'Mollusca', 'Actinopterygii', 'Other zooplankton','Detritus'])).groupby(['latitude_bin','ECD_mid']).apply(lambda x: pd.Series(dict(zip(pd.Categorical(x.Instrument,df.Instrument.unique()).value_counts().index,(pd.Categorical(x.Instrument,df.Instrument.unique()).value_counts().values>=1).astype(int))))).reset_index().reset_index()
        plot = (ggplot(data=df_instruments.melt(id_vars=['latitude_bin','ECD_mid'],value_vars=df.Instrument.unique(),var_name='Instrument',value_name='value')) +
                stat_summary(aes(x="ECD_mid", y="Instrument", fill='value>0'),geom='point', size=3.5, shape='H', color='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50), alpha=1) +
                scale_fill_manual(values={True:'black',False:'#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50)},guide=False)+
                labs(fill='Selected instruments:', x=r'Equivalent circular diameter ($\mu$m)', y=r'') +
                scale_x_log10(breaks=[size for size in np.sort(np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))],labels=lambda l: [int(size) if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in l]) +theme_paper).draw(show=True)
        plot.set_size_inches(8.5, 1.8)
        plot.savefig(fname='{}/PSSdb_Merged-products_selection_Instruments_ECD.svg'.format(str(path_to_directory)), dpi=300)
        plot = (ggplot(data=df_instruments.melt(id_vars=['latitude_bin','ECD_mid'],value_vars=df.Instrument.unique(),var_name='Instrument',value_name='value')) +
                stat_summary(aes(x="latitude_bin", y="Instrument", fill='value>0'),geom='point', size=2.5, shape='H', color='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50), alpha=1) +
                coord_flip()+
                scale_fill_manual(values={True:'black',False:'#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50)},guide=False)+
                labs(fill='Selected instruments:', y=r'', x=r'Latitude ($\degree$N)') +
                theme_paper+theme(axis_text_x=element_text(rotation=45))).draw(show=True)
        plot.set_size_inches(1.2, 2.9)
        plot.savefig(fname='{}/PSSdb_Merged-products_selection_Instruments_latitude.svg'.format(str(path_to_directory)), dpi=300)

        #x=df_selection.loc[list(df_selection.groupby(['Latitude', 'Longitude', 'year', 'month','min_depth','max_depth','Group_station_index','n_instruments','merged_instruments','sizeClasses']).groups.values())[72]]
        df_merged=df_selection.groupby(['Latitude', 'Longitude', 'date_bin','min_depth','max_depth','Group_station_index','n_instruments','merged_instruments','sizeClasses']).apply(lambda x:pd.Series({**{'Validation_percentage':np.nanmean(x.Validation_percentage)},**dict(zip((pd.Series(dict_products_y[product])+'_living').to_list(),x.query('PFT!="Detritus"')[dict_products_y[product]].sum(axis=0))),**dict(zip((pd.Series(dict_products_y[product])+'_nonliving').to_list(),x.query('PFT=="Detritus"')[dict_products_y[product]].sum(axis=0)))})).reset_index()#df_selection.groupby(['Latitude', 'Longitude', 'year', 'month','min_depth','max_depth','Group_station_index','n_instruments','merged_instruments','sizeClasses']).apply(lambda x:pd.Series({**{'n':x.n.min()},**dict(zip((pd.Series(dict_products_y[product])+'_living').to_list(),x.query('PFT!="Detritus"')[dict_products_y[product]].sum(axis=0))),**dict(zip((pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_living').to_list(),((x.query('PFT!="Detritus"')[(pd.Series(dict_products_y[product]).str.replace('_mean','_std')).to_list()]**2).sum(axis=0))**0.5)),**dict(zip((pd.Series(dict_products_y[product])+'_nonliving').to_list(),x.query('PFT=="Detritus"')[dict_products_y[product]].sum(axis=0))),**dict(zip((pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_nonliving').to_list(),((x.query('PFT=="Detritus"')[(pd.Series(dict_products_y[product]).str.replace('_mean','_std')).to_list()]**2).sum(axis=0))**0.5))})).reset_index()
        df_merged=pd.merge(df_merged,df_bins,how='left',on='sizeClasses').sort_values(['Latitude', 'Longitude', 'date_bin','min_depth','max_depth','ECD_mid']).reset_index(drop=True)

        # Check merged products along with original data after generating final product
        df_merged['Sample']=df_merged.Latitude.astype(str) + '_' + df_merged.Longitude.astype(str) + "_" + df_merged.date_bin.astype( str)
        ## first, apply thresholding function, then average at bin level
        df_merged_thres=df_merged.groupby(['Sample','Latitude','Longitude','date_bin','min_depth','max_depth','Group_station_index','n_instruments','merged_instruments']).apply(lambda x: (threshold_func(binned_data=x.rename(columns=dict(zip((pd.Series(dict_products_y[product])+'_living').to_list(),(pd.Series(dict_products_y[product])).to_list()))).assign(count_uncertainty=0,size_uncertainty=0), empty_bins=3, threshold_count=0.2, threshold_size=0.2))).reset_index(drop=True)
        df_products_1a=NBSS_stats_func(df=df_merged_thres.rename(columns={'Latitude':'midLatBin','Longitude':'midLonBin','min_depth':'Min_obs_depth','max_depth':'Max_obs_depth'}).astype({'midLatBin':float,'midLonBin':float}), light_parsing=False, bin_loc=1, group_by='yyyymm')
        df_products_1a=pd.merge(df_products_1a,df_merged_thres.assign(latitude=lambda x:pd.cut(x.Latitude.astype(float),bins=np.arange(-90,91,1),labels=np.arange(-90,91,0.5)[1:-1:2]).astype(str),longitude=lambda x:pd.cut(x.Longitude.astype(float),bins=np.arange(-180,181,1),labels=np.arange(-180,181,0.5)[1:-1:2]).astype(str),year=lambda x:(x.date_bin.str[0:4]).astype(str),month=lambda x:(x.date_bin.str.split("_").str[-2]).astype(str)).astype(dict(zip(['latitude','longitude','year','month','biovolume_size_class','min_depth','max_depth','n_instruments','merged_instruments'],[str]*9))).groupby(['latitude','longitude','year','month','biovolume_size_class','min_depth','max_depth']).apply(lambda x: pd.Series({'n_instruments':len(x.merged_instruments.unique()[0].split("_")),'merged_instruments':'_'.join(natsorted(x.merged_instruments.unique()[0].split("_"))),'Validation_percentage':x.Validation_percentage.mean()})).reset_index().astype((df_products_1a.dtypes.loc[['latitude','longitude','year','month','min_depth','max_depth','biovolume_size_class']]).to_dict()),how='left',on=['latitude','longitude','year','month','min_depth','max_depth','biovolume_size_class'])
        df_products_1b=df_merged_thres.rename(columns={'Latitude':'midLatBin','Longitude':'midLonBin','min_depth':'Min_obs_depth','max_depth':'Max_obs_depth'}).astype(dict(zip(['midLatBin','midLonBin','date_bin','Min_obs_depth','Max_obs_depth','n_instruments','merged_instruments'],[str]*7))).groupby(['midLatBin','midLonBin','date_bin','Min_obs_depth','Max_obs_depth','n_instruments','merged_instruments']).apply(lambda x:linear_fit_func(x.astype({'midLatBin':float,'midLonBin':float,'Min_obs_depth':float,'Max_obs_depth':float}).assign(logSize=np.log10(x[dict_products_x[product]]),logNB=np.log10(x[dict_products_y[product][0]]),logECD=np.log10(x.ECD_mid),logPSD=np.log10(x.PSD) if 'PSD' in x.columns else np.nan), light_parsing = False, depth_parsing = False)).reset_index()
        df_products_1b=pd.merge(stats_linfit_func(df=df_products_1b, light_parsing = False, bin_loc = 1, group_by = 'yyyymm'),df_products_1b.assign(latitude=lambda x:pd.cut(x.midLatBin.astype(float),bins=np.arange(-90,91,1),labels=np.arange(-90,91,0.5)[1:-1:2]).astype(str),longitude=lambda x:pd.cut(x.midLonBin.astype(float),bins=np.arange(-180,181,1),labels=np.arange(-180,181,0.5)[1:-1:2]).astype(str),year=lambda x:(x.date_bin.str[0:4]).astype(str),month=lambda x:(x.date_bin.str.split("_").str[-1]).astype(str)).astype(dict(zip(['latitude','longitude','year','month','Min_obs_depth','Max_obs_depth'],[str]*6))).groupby(['latitude','longitude','year','month','Min_obs_depth','Max_obs_depth']).apply(lambda x: pd.Series({'n_instruments':len(x.merged_instruments.unique()[0].split("_")),'merged_instruments':'_'.join(natsorted(x.merged_instruments.unique()[0].split("_")))})).reset_index().rename(columns={'Min_obs_depth':'min_depth','Max_obs_depth':'max_depth'}).astype({'latitude':float,'longitude':float,'year':str,'month':str,'min_depth':float,'max_depth':float}),how='left',on=['latitude','longitude','year','month','min_depth','max_depth'])
        df_products_merged=pd.merge(df_products_1a,df_products_1b,how='left',on=['merged_instruments','year','month','latitude','longitude','min_depth','max_depth','n'])


        #x=df.loc[list(df.groupby(group+['Sample','sizeClasses']).groups.values())[3]]
        df_summary=df.groupby(group+['Sample','Instrument','sizeClasses']).apply(lambda x:pd.Series({**dict(zip((pd.Series(dict_products_y[product])+'_living').to_list(),x.query('PFT!="Detritus"')[dict_products_y[product]].sum(axis=0))),**dict(zip((pd.Series(dict_products_y[product])+'_nonliving').to_list(),x.query('PFT=="Detritus"')[dict_products_y[product]].sum(axis=0)))})).reset_index()#.apply(lambda x:pd.Series({**dict(zip((pd.Series(dict_products_y[product])+'_living').to_list(),x.query('PFT!="Detritus"')[dict_products_y[product]].sum(axis=0))),**dict(zip((pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_living').to_list(),((x.query('PFT!="Detritus"')[(pd.Series(dict_products_y[product]).str.replace('_mean','_std')).to_list()]**2).sum(axis=0))**0.5)),**dict(zip((pd.Series(dict_products_y[product])+'_nonliving').to_list(),x.query('PFT=="Detritus"')[dict_products_y[product]].sum(axis=0))),**dict(zip((pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_nonliving').to_list(),((x.query('PFT=="Detritus"')[(pd.Series(dict_products_y[product]).str.replace('_mean','_std')).to_list()]**2).sum(axis=0))**0.5))})).reset_index()
        df_summary=pd.merge(df_summary.assign(merged_instruments=df_summary.Instrument),df_merged[['Sample','min_depth','max_depth','sizeClasses','ECD_mid','biomass_mid']].drop_duplicates(),how='left',on=['Sample','sizeClasses'])
        df_summary=df_summary.dropna(subset=['min_depth', 'max_depth']).sort_values(['Sample', 'Latitude', 'Longitude', 'date_bin', 'min_depth', 'max_depth', 'Group_station_index', 'n_instruments', 'merged_instruments','ECD_mid']).reset_index(drop=True)
        df_summary_thres=df_summary.groupby(['Sample','Latitude','Longitude','date_bin','min_depth','max_depth','Group_station_index','n_instruments','merged_instruments']).apply(lambda x: (threshold_func(binned_data=x.rename(columns=dict(zip((pd.Series(dict_products_y[product])+'_living').to_list(),(pd.Series(dict_products_y[product])).to_list()))).assign(count_uncertainty=0,size_uncertainty=0), empty_bins=3, threshold_count=0.2, threshold_size=0.2))).reset_index(drop=True).rename(columns={'Latitude':'midLatBin','Longitude':'midLonBin','min_depth':'Min_obs_depth','max_depth':'Max_obs_depth'}).astype({'midLatBin':float,'midLonBin':float}).assign(size_class_mid=lambda x:(1/6)*np.pi*(x.ECD_mid**3))
        df_1a_original=df_summary_thres.groupby(['merged_instruments']).apply(lambda x: NBSS_stats_func(df=x, light_parsing=False, bin_loc=1, group_by='yyyymm')).reset_index().drop(columns='level_1')
        df_1b_original=df_summary_thres.astype(dict(zip(['midLatBin','midLonBin','date_bin','Min_obs_depth','Max_obs_depth','n_instruments','merged_instruments'],[str]*7))).groupby(['midLatBin','midLonBin','date_bin','Min_obs_depth','Max_obs_depth','n_instruments','merged_instruments']).apply(lambda x:linear_fit_func(x.astype({'midLatBin':float,'midLonBin':float,'Min_obs_depth':float,'Max_obs_depth':float}).assign(logSize=np.log10(x[dict_products_x[product]]),logNB=np.log10(x[dict_products_y[product][0]]),logECD=np.log10(x.ECD_mid),logPSD=np.log10(x.PSD) if 'PSD' in x.columns else np.nan), light_parsing = False, depth_parsing = False)).reset_index()
        df_1b_original=df_1b_original.groupby(['merged_instruments']).apply(lambda x: stats_linfit_func(df=x, light_parsing = False, bin_loc = 1, group_by = 'yyyymm')).reset_index().drop(columns='level_1')
        df_products_original=pd.merge(df_1a_original,df_1b_original,how='left',on=['merged_instruments','year','month','latitude','longitude','min_depth','max_depth','n'])

        df_all=pd.concat([df_products_merged[list(set(df_products_merged.columns) & set(df_products_original.columns))],df_products_original[list(set(df_products_merged.columns) & set(df_products_original.columns))]],axis=0).sort_values(['latitude','longitude','year','month','min_depth','max_depth','merged_instruments','equivalent_circular_diameter_mean']).reset_index(drop=True)[['latitude','longitude','year','month','min_depth','max_depth','merged_instruments','equivalent_circular_diameter_mean','biovolume_size_class','n','normalized_biovolume_mean','normalized_biovolume_std','normalized_abundance_mean','normalized_abundance_std','NBSS_slope_mean','NBSS_slope_std','NBSS_intercept_mean','NBSS_intercept_std','NBSS_r2_mean','NBSS_r2_std','PSD_slope_mean','PSD_slope_std','PSD_intercept_mean','PSD_intercept_std','PSD_r2_mean','PSD_r2_std','N_bins_min']]#pd.concat([df_merged[list(set(df_merged.columns) & set(df_summary.columns))],df_summary[list(set(df_merged.columns) & set(df_summary.columns))]],axis=0).sort_values(['Sample','min_depth','max_depth','n','merged_instruments','ECD_mid']).reset_index(drop=True)[group+['Sample','min_depth','max_depth','sizeClasses','ECD_mid','biomass_mid']+(pd.Series(dict_products_y[product])+'_living').to_list()+(pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_living').to_list()+(pd.Series(dict_products_y[product])+'_nonliving').to_list()+(pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_nonliving').to_list()]#pd.concat([df_merged[list(set(df_merged.columns) & set(df_summary.columns))],df_summary[list(set(df_merged.columns) & set(df_summary.columns))]],axis=0).sort_values(['Sample','min_depth','max_depth','merged_instruments','ECD_mid']).reset_index(drop=True)[group+['Sample','min_depth','max_depth','sizeClasses','ECD_mid','biomass_mid']+(pd.Series(dict_products_y[product])+'_living').to_list()+(pd.Series(dict_products_y[product])+'_nonliving').to_list()]#pd.concat([df_merged[list(set(df_merged.columns) & set(df_summary.columns))],df_summary[list(set(df_merged.columns) & set(df_summary.columns))]],axis=0).sort_values(['Sample','min_depth','max_depth','n','merged_instruments','ECD_mid']).reset_index(drop=True)[group+['Sample','min_depth','max_depth','sizeClasses','ECD_mid','biomass_mid']+(pd.Series(dict_products_y[product])+'_living').to_list()+(pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_living').to_list()+(pd.Series(dict_products_y[product])+'_nonliving').to_list()+(pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_nonliving').to_list()]
        group = ['Group_station_index','Sample', 'merged_instruments']
        df_all=pd.merge(df_all, df_all.drop_duplicates(subset=['year','month','latitude','longitude','min_depth','max_depth'], ignore_index=True)[['year','month','latitude','longitude','min_depth','max_depth']].reset_index().rename({'index': 'Group_station_index'}, axis='columns'), how='left', on=['year','month','latitude','longitude','min_depth','max_depth'])
        df_all['Sample']=df_all.latitude.astype(str) + '_' + df_all.longitude.astype(str) + "_" + df_all.year.astype( str)+ "_" + df_all.month.astype( str)
        df_all = pd.merge(df_all, df_all.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename( {'index': 'Group_index'}, axis='columns'), how='left', on=group)
        df_all=df_all.dropna(subset=['min_depth','max_depth']).reset_index(drop=True) # Some initial observations were dropped since there were above the count/size uncertainty threshold
        df_all['Group_index']=df_all.Group_station_index.astype(str)+"_"+df_all.merged_instruments.astype(str)
        df_all['Sample']=df_all.Sample.astype(str)+"_"+df_all.groupby(['Sample']).merged_instruments.transform(lambda x: x[x.str.len().idxmax()]  )
        df_all['Group_index'] = df_all.Group_index.astype(str) + "_slope" + df_all.NBSS_slope_mean.astype(float).round(2).astype(str) + "_intercept" + df_all.NBSS_intercept_mean.astype(float).round(2).astype(str)
        fig = standardization_report_func(df_summary=df_all.query('merged_instruments.str.contains("_")').assign(Sample=lambda x:x.latitude.astype(str)).rename(columns={'latitude':'Latitude','longitude':'Longitude'}).groupby(['Sample', 'Latitude', 'Longitude', 'min_depth', 'max_depth'], dropna=True).apply(lambda x: pd.Series({'Abundance': x['normalized_biovolume_mean'].astype(float).sum(),  # individuals per liter
                                                                                                                                                                                     'Average_diameter': np.nanmean( x['normalized_biovolume_mean'].astype(float)* x.equivalent_circular_diameter_mean.astype(float)),# micrometer
                                                                                                                                                                                     'Std_diameter': np.nanstd( x['normalized_biovolume_mean'].astype(float) * x.equivalent_circular_diameter_mean.astype(float))})).reset_index(),
            df_standardized=pd.DataFrame({}), df_nbss=df_all.dropna(subset=['NBSS_slope_mean']).rename( columns={'normalized_biovolume_mean': 'NBSS', 'equivalent_circular_diameter_mean': 'size_class_mid'}),plot='nbss')
        fig.write_html(Path(path_to_directory/'Merged_product_{}_living_multiinstruments.html'.format(product)))

        """
        # Compute slope and intercept (reference at 100 micrometers to be comparable) to append to Group_index
        #x=df_all.rename(columns={dict_products_y[product][0]+'_living': 'NBSS', 'ECD_mid': 'size_class_mid'}).loc[list(df_all.rename(columns={dict_products_y[product][0]+'_living': 'NBSS', 'ECD_mid': 'size_class_mid'}).groupby(['Sample','Group_index']).groups.values())[0]]
        df_stats=df_all.rename(columns={dict_products_y[product][0]+'_living': 'NBSS', 'ECD_mid': 'size_class_mid'}).groupby(['Sample','Group_index']).apply(lambda x: pd.DataFrame(regress_nbss(nbss=x.assign(size_class_mid=x.size_class_mid/100),threshold_count=0.2,threshold_size=0.2,n_bins=3),index=x.index))
        df_all=pd.concat([df_all,df_stats],axis=1)
        df_all.loc[df_all.slope.isna(),['slope','intercept']]='nan'
        df_all['Group_index']=df_all.Group_index.astype(str)+"_slope"+df_all.slope.astype(float).round(2).astype(str)+"_intercept"+df_all.intercept.astype(float).round(2).astype(str)
        
        fig = standardization_report_func(df_summary=df_all.query('merged_instruments.str.contains("_")').groupby(['Sample', 'Latitude', 'Longitude', 'min_depth', 'max_depth'], dropna=True).apply(lambda x: pd.Series({'Abundance': x[dict_products_y[product][0]+'_living'].astype(float).sum(),  # individuals per liter
                                                                                                                                                                                     'Average_diameter': np.nanmean( x[dict_products_y[product][0]+'_living'].astype(float)* x.ECD_mid.astype(float)),# micrometer
                                                                                                                                                                                     'Std_diameter': np.nanstd( x[dict_products_y[product][0]+'_living'].astype(float) * x.ECD_mid.astype(float))})).reset_index(),
            df_standardized=pd.DataFrame({}), df_nbss=df_all.query('selection==True').rename( columns={dict_products_y[product][0]+'_living': 'NBSS', 'ECD_mid': 'size_class_mid'}),plot='nbss')
        fig.write_html(Path(path_to_directory/'Merged_product_{}_living.html'.format(product)))
        fig = standardization_report_func(df_summary=df_all.query('merged_instruments.str.contains("_")').groupby(['Sample', 'Latitude', 'Longitude', 'min_depth', 'max_depth'], dropna=True).apply(lambda x: pd.Series({'Abundance': x[dict_products_y[product][0]+'_nonliving'].astype(float).sum(),  # individuals per liter
                                                                                                                                                                                     'Average_diameter': np.nanmean( x[dict_products_y[product][0]+'_nonliving'].astype(float)* x.ECD_mid.astype(float)),# micrometer
                                                                                                                                                                                     'Std_diameter': np.nanstd( x[dict_products_y[product][0]+'_nonliving'].astype(float) * x.ECD_mid.astype(float))})).reset_index(),
            df_standardized=pd.DataFrame({}), df_nbss=df_all.rename( columns={dict_products_y[product][0]+'_nonliving': 'NBSS', 'ECD_mid': 'size_class_mid'}),plot='nbss')
        fig.write_html(Path(path_to_directory/'Merged_product_{}_nonliving.html'.format(product)))
        """

        # Select only regression statistics to re-create the ESSD paper figures, append Global Oceans and Seas regions and save
        df_stats_summary=df_all.dropna(subset=['NBSS_slope_mean']).rename(columns={'latitude':'Latitude','longitude':'Longitude','NBSS_slope_mean':'slope','NBSS_intercept_mean':'intercept','NBSS_r2_mean':'r2'})[['Latitude','Longitude','Sample','year','month','merged_instruments','Group_index','Group_station_index','slope','intercept','r2']].drop_duplicates()#pd.merge(df_all[['Latitude','Longitude','Sample','year','month','merged_instruments','Group_index','Group_station_index']].drop_duplicates(),df_stats[['slope','intercept','R2']].drop_duplicates(),how='left',left_index=True, right_index=True)#pd.merge(df_all[['Latitude','Longitude','Sample','date_bin','merged_instruments','Group_index','Group_station_index']].drop_duplicates(),df_stats[['slope','intercept','R2']].drop_duplicates(),how='left',left_index=True, right_index=True)#pd.merge(df_all[['Latitude','Longitude','Sample','year','month','merged_instruments','Group_index','Group_station_index']].drop_duplicates(),df_stats[['slope','intercept','R2']].drop_duplicates(),how='left',left_index=True, right_index=True)
        gdf = gpd.GeoDataFrame(df_stats_summary[['Group_station_index', 'Longitude', 'Latitude']].drop_duplicates().dropna(),geometry=gpd.points_from_xy(df_stats_summary[['Group_station_index', 'Longitude', 'Latitude']].drop_duplicates().dropna().Longitude, df_stats_summary[['Group_station_index', 'Longitude', 'Latitude']].drop_duplicates().dropna().Latitude))
        df_stats_summary['Study_area'] = pd.merge(df_stats_summary,gpd.tools.sjoin(gdf, oceans, predicate="within", how='left')[['Group_station_index', 'name']], how='left', on='Group_station_index')['name'].astype(str)
        #df_stats_summary.to_csv(path_to_directory/'Merged_product_{}_stats_living.csv'.format(product),index=False)
        # Plot the selection for each PFT
        df_pivot=df_stats_summary.pivot_table(values=['slope','intercept', 'r2'],columns=['merged_instruments'],index=['Latitude', 'Longitude', 'Sample', 'year','month', 'Group_station_index']).reset_index()
        df_pivot.columns=['_'.join(col) if (type(col) is tuple) & (col[1]!='') else col[0] for col in df_pivot.columns.values]

        g = sns.pairplot(df_pivot[[column for column in df_pivot.columns if 'slope' in column]],kind='scatter',corner=True,diag_kind='kde',plot_kws=dict(marker="o",size=.3, color='black'),diag_kws=dict(color='black'))
        g.map_offdiag(plot_unity)
        g.map_offdiag(r2)
        g.fig.set_size_inches(10,10)
        g.savefig(fname='{}/PSSdb_Merged-products_slopes_comparison.svg'.format(str(path_to_directory)), dpi=300)

        plot = (ggplot(data=df_pivot) +
                geom_point(aes(x="slope_IFCB", y="slope_UVP",group='Sample'), size=1.5, shape='H', color='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50), alpha=0.1) +
                geom_abline(slope=1,intercept=0)+scale_x_continuous(limits=[min(df_pivot[['slope_IFCB','slope_UVP']].describe().loc['min']),max(df_pivot[['slope_IFCB','slope_UVP']].describe().loc['max'])])+scale_y_continuous(limits=[min(df_pivot[['slope_IFCB','slope_UVP']].describe().loc['min']),max(df_pivot[['slope_IFCB','slope_UVP']].describe().loc['max'])])+
                labs(x=r'IFCB spectral slope (L$^{-1}$ $\mu$m$^{-3}$)', y=r'UVP spectral slope (L$^{-1}$ $\mu$m$^{3}$)') +
                #scale_x_log10(breaks=[size  for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))],labels=lambda l: [int(size) if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in l])+
                theme_paper).draw(show=True)
        plot.set_size_inches(3,3)
        plot.savefig(fname='{}/PSSdb_Merged-products_slopes_comparison.svg'.format(str(path_to_directory)), dpi=300)
        plot = (ggplot(data=df_pivot) +
                geom_point(aes(x="intercept_IFCB", y="intercept_UVP", group='Sample'), size=1.5, shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50), alpha=0.1) +
                geom_abline(slope=1, intercept=0) +
                labs(x=r'IFCB spectral intercept ($\mu$m$^{3}$ L$^{-1}$ $\mu$m$^{-3}$)',  y=r'UVP spectral intercept ($\mu$m$^{3}$ L$^{-1}$ $\mu$m$^{-3}$)') +
                scale_y_continuous(limits=[min(df_pivot[['intercept_IFCB','intercept_UVP']].describe().loc['min']),max(df_pivot[['intercept_IFCB','intercept_UVP']].describe().loc['max'])],breaks=[size for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))], labels=lambda l: [int(size) if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else ''  for size in l]) +
                scale_x_continuous(limits=[min(df_pivot[['intercept_IFCB','intercept_UVP']].describe().loc['min']),max(df_pivot[['intercept_IFCB','intercept_UVP']].describe().loc['max'])],breaks=[size  for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))],labels=lambda l: [int(size) if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in l])+
                theme_paper).draw(show=True)
        plot.set_size_inches(3, 3)
        plot.savefig(fname='{}/PSSdb_Merged-products_intercepts_comparison.svg'.format(str(path_to_directory)), dpi=300)

        # Update progress bar and move to next product
        ok = bar.update(n=1)



#1) calculate NBSS products no adjustments, and without considering plankton functional groups
merged_raw_list = merge_taxa_products(grouping_factors= ['date_bin', 'Station_location', 'midLatBin', 'midLonBin'])





#2) merge PFT products calculated with different size metrics and calculate
merged_adjusted_raw_list = merge_adjust_taxa_products(grouping_factors= ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'PFT'])


#2) Calculate NBSS for each of these products


