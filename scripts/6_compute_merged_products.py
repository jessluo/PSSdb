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
import datetime
from natsort import natsorted
import shutil
from tqdm import tqdm
from functools import reduce
try:
    from funcs_merged_products import *
    from funcs_NBS import *
    from funcs_standardize_projects import *
except:
    from scripts.funcs_merged_products import *
    from scripts.funcs_NBS import *
    from scripts.funcs_standardize_projects import *


## Workflow starts here:

# Load most recent NBSS taxa-specific products
path_to_products=natsorted(list(natsorted(list(Path(path_to_git / cfg['dataset_subdir'] /cfg['product_subdir']).glob('NBSS_ver_*')))[0].rglob('PFT/*_1a_*')))
grouping_factor=['ocean','year','month','latitude','longitude']
path_to_directory=Path(path_to_products[0].parent.parent / 'Merged')
path_to_directory.mkdir(exist_ok=True)
with tqdm(desc='Merging taxa-specific {} products', total=3, bar_format='{desc}{bar}', position=0, leave=True) as bar:
    for product in ['Size','Biomass','Weight']:
        bar.set_description("Merging taxa-specific {} products".format(product), refresh=True)
        # Select the associated files
        path_to_files=[path for path in path_to_products if product in  path.name]
        df_nbss=merge_products(path_to_products=path_to_files, grouping_factors=['ocean', 'year', 'month', 'latitude', 'longitude'])
        df = df_nbss[df_nbss.n_instruments > 1]
        group = ['latitude', 'longitude', 'year', 'month']
        df = pd.merge(df, df.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename({'index': 'Group_station_index'}, axis='columns'), how='left', on=group)
        df = df.assign(Sample=df.latitude.astype(str) + '_' + df.longitude.astype(str) + "_" + df.year.astype( str) + df.month.astype(str).str.zfill(2)).rename(columns={'longitude': 'Longitude', 'latitude': 'Latitude'})
        ## Apply Soviadan et al. maximum selection method
        # Use class mid-point (biovolume, biomass, or weight) to assign the rest of the size classes information
        data_bins = df_bins.drop(columns=['sizeClasses']).assign(sizeClasses=pd.cut(np.append(0,df_bins[dict_products_x[product]].to_numpy()),np.append(0,df_bins[dict_products_x[product]].to_numpy())).categories.values.astype(str))
        df=pd.merge(df.drop(columns=[column for column in df.columns if column in data_bins.columns]).assign(sizeClasses=pd.cut(df[dict_products_x[product]],data_bins[dict_products_x[product]]).astype(str)),data_bins,how='left',on='sizeClasses')
        # Use the maximum normalized abundance / biovolume / biomass / weight in each size class
        #x=df.loc[list(df.groupby(group+['PFT','sizeClasses']).groups.values())[8]]
        group=['Latitude', 'Longitude', 'year', 'month','Group_station_index','n_instruments','merged_instruments']
        df_selection=df.groupby(group+['PFT','sizeClasses']).apply(lambda x: pd.Series({**{'min_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.min(),'max_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.max(),'selected_instrument':x.loc[(x[dict_products_y[product]]).idxmax(axis=0),'Instrument'].unique()[0]},**dict(zip(dict_products_y[product], (x[dict_products_y[product]]).max(axis=0))),**dict(zip([variable.replace('_mean','_std') for variable in dict_products_y[product]],np.diag(x.loc[(x[dict_products_y[product]]).idxmax(axis=0),[variable.replace('_mean','_std') for variable in dict_products_y[product]]]))),**{'n':x.loc[(x[dict_products_y[product]]).idxmax(axis=0),'n'].unique()[0]}})).reset_index()#df.groupby(group+['PFT','sizeClasses']).apply(lambda x: pd.Series({**{'min_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.min(),'max_depth':df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.max(),'selected_instrument':x.loc[(x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).idxmax(axis=0),'Instrument'].unique()[0]},**dict(zip(dict_products_y[product], ( x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).max(axis=0))),**dict(zip([variable.replace('_mean','_std') for variable in dict_products_y[product]],np.diag(x.loc[(x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).idxmax(axis=0),[variable.replace('_mean','_std') for variable in dict_products_y[product]]]))),**{'n':x.loc[(x[dict_products_y[product]]*np.repeat(np.where((x.max_depth.astype(float) - x.min_depth.astype(float)) > 1,(x.max_depth.astype(float) - x.min_depth.astype(float)), 1) / np.where( (df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()) > 1,(df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].max_depth.astype(float).max() - df[df[group].apply(tuple, axis=1).isin(x[group].apply(tuple, axis=1))].min_depth.astype(float).min()), 1),len(dict_products_y[product])).reshape(len(x),len(dict_products_y[product]))).idxmax(axis=0),'n'].unique()[0]}})).reset_index()
        # Plot the selection for each PFT
        plot = (ggplot(data=pd.merge(df_selection,data_bins,how='left',on='sizeClasses').sort_values(['Latitude', 'Longitude', 'Group_station_index','PFT','ECD_mid']).reset_index(drop=True).assign(latitude_bin=lambda x:pd.cut(x.Latitude.astype(float),bins=np.arange(-90,90,1),labels=np.arange(-90,90,1)[1:]).astype(float),PFT=lambda x: pd.Categorical(x.PFT.str.replace("_"," "),['Nanophytoplankton','Microphytoplankton','Mesophytoplankton','Ciliophora','Rhizaria','Crustaceans','Filter feeder gelatinous','Carnivorous gelatinous','Mollusca','Actinopterygii','Other zooplankton','Detritus']))) +
                facet_wrap('~PFT', nrow=2) +
                stat_summary(aes(x="ECD_mid", y="latitude_bin",group='latitude_bin', fill='selected_instrument'),geom='point', size=1.5, shape='H', color='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50), alpha=0.1) +
                labs(x=r'Equivalent circular diameter ($\mu$m)', y=r'Latitude ($\degree$N)') +
                scale_x_log10(breaks=[size  for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))],labels=lambda l: [int(size) if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in l])+
                theme_paper).draw(show=True)
        plot.set_size_inches(8.5, 3.8)
        plot.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/PSSdb_Merged-products_selection.svg'.format(str(Path.home())),limitsize=False, dpi=600)

        #x=df_selection.loc[list(df_selection.groupby(['Latitude', 'Longitude', 'year', 'month','min_depth','max_depth','Group_station_index','n_instruments','merged_instruments','sizeClasses']).groups.values())[72]]
        df_merged=df_selection.groupby(['Latitude', 'Longitude', 'year', 'month','min_depth','max_depth','Group_station_index','n_instruments','merged_instruments','sizeClasses']).apply(lambda x:pd.Series({**{'n':x.n.min()},**dict(zip((pd.Series(dict_products_y[product])+'_living').to_list(),x.query('PFT!="Detritus"')[dict_products_y[product]].sum(axis=0))),**dict(zip((pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_living').to_list(),((x.query('PFT!="Detritus"')[(pd.Series(dict_products_y[product]).str.replace('_mean','_std')).to_list()]**2).sum(axis=0))**0.5)),**dict(zip((pd.Series(dict_products_y[product])+'_nonliving').to_list(),x.query('PFT=="Detritus"')[dict_products_y[product]].sum(axis=0))),**dict(zip((pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_nonliving').to_list(),((x.query('PFT=="Detritus"')[(pd.Series(dict_products_y[product]).str.replace('_mean','_std')).to_list()]**2).sum(axis=0))**0.5))})).reset_index()
        df_merged=pd.merge(df_merged,data_bins,how='left',on='sizeClasses').sort_values(['Latitude', 'Longitude', 'year', 'month','min_depth','max_depth','ECD_mid']).reset_index(drop=True)

        # Check merged products along with original data
        df_merged['Sample']=df_merged.Latitude.astype(str) + '_' + df_merged.Longitude.astype(str) + "_" + df_merged.year.astype( str) + df_merged.month.astype(str).str.zfill(2)
        #x=df.loc[list(df.groupby(group+['Sample','sizeClasses']).groups.values())[3]]
        df_summary=df.groupby(group+['Sample','Instrument','sizeClasses']).apply(lambda x:pd.Series({**dict(zip((pd.Series(dict_products_y[product])+'_living').to_list(),x.query('PFT!="Detritus"')[dict_products_y[product]].sum(axis=0))),**dict(zip((pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_living').to_list(),((x.query('PFT!="Detritus"')[(pd.Series(dict_products_y[product]).str.replace('_mean','_std')).to_list()]**2).sum(axis=0))**0.5)),**dict(zip((pd.Series(dict_products_y[product])+'_nonliving').to_list(),x.query('PFT=="Detritus"')[dict_products_y[product]].sum(axis=0))),**dict(zip((pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_nonliving').to_list(),((x.query('PFT=="Detritus"')[(pd.Series(dict_products_y[product]).str.replace('_mean','_std')).to_list()]**2).sum(axis=0))**0.5))})).reset_index()
        df_summary=pd.merge(df_summary.assign(merged_instruments=df_summary.Instrument),df_merged[['Sample','min_depth','max_depth','n','sizeClasses','ECD_mid','biomass_mid']].drop_duplicates(),how='left',on=['Sample','sizeClasses'])
        df_all=pd.concat([df_merged[list(set(df_merged.columns) & set(df_summary.columns))],df_summary[list(set(df_merged.columns) & set(df_summary.columns))]],axis=0).sort_values(['Sample','min_depth','max_depth','n','merged_instruments','ECD_mid']).reset_index(drop=True)[group+['Sample','min_depth','max_depth','sizeClasses','ECD_mid','biomass_mid']+(pd.Series(dict_products_y[product])+'_living').to_list()+(pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_living').to_list()+(pd.Series(dict_products_y[product])+'_nonliving').to_list()+(pd.Series(dict_products_y[product]).str.replace('_mean','_std')+'_nonliving').to_list()]
        group = ['Group_station_index','Sample', 'merged_instruments']
        df_all = pd.merge(df_all, df_all.drop_duplicates(subset=group, ignore_index=True)[group].reset_index().rename( {'index': 'Group_index'}, axis='columns'), how='left', on=group)
        df_all['Group_index']=df_all.Group_index.astype(str)+"_"+df_all.merged_instruments.astype(str)
        df_all['Sample']=df_all.Sample.astype(str)+"_"+df_all.groupby(['Sample']).merged_instruments.transform(lambda x: x[x.str.len().idxmax()]  )
        # Compute slope and intercept (reference at 100 micrometers to be comparable) to append to Group_index
        #x=df_all.rename(columns={dict_products_y[product][0]+'_living': 'NBSS', 'ECD_mid': 'size_class_mid'}).loc[list(df_all.rename(columns={dict_products_y[product][0]+'_living': 'NBSS', 'ECD_mid': 'size_class_mid'}).groupby(['Sample','Group_index']).groups.values())[0]]
        df_stats=df_all.rename(columns={dict_products_y[product][0]+'_living': 'NBSS', 'ECD_mid': 'size_class_mid'}).groupby(['Sample','Group_index']).apply(lambda x: pd.DataFrame(regress_nbss(nbss=x.assign(size_class_mid=x.size_class_mid/100),threshold_count=0.2,threshold_size=0.2,n_bins=3),index=x.index))
        df_all=pd.concat([df_all,df_stats],axis=1)
        df_all.loc[df_all.slope.isna(),['slope','intercept']]='nan'
        df_all['Group_index']=df_all.Group_index.astype(str)+"_slope"+df_all.slope.astype(float).round(2).astype(str)+"_intercept"+df_all.intercept.astype(float).round(2).astype(str)

        fig = standardization_report_func(df_summary=df_all.query('merged_instruments.str.contains("_")').groupby(['Sample', 'Latitude', 'Longitude', 'min_depth', 'max_depth'], dropna=True).apply(lambda x: pd.Series({'Abundance': x[dict_products_y[product][0]+'_living'].astype(float).sum(),  # individuals per liter
                                                                                                                                                                                     'Average_diameter': np.nanmean( x[dict_products_y[product][0]+'_living'].astype(float)* x.ECD_mid.astype(float)),# micrometer
                                                                                                                                                                                     'Std_diameter': np.nanstd( x[dict_products_y[product][0]+'_living'].astype(float) * x.ECD_mid.astype(float))})).reset_index(),
            df_standardized=pd.DataFrame({}), df_nbss=df_all.rename( columns={dict_products_y[product][0]+'_living': 'NBSS', 'ECD_mid': 'size_class_mid'}),plot='nbss')
        fig.write_html(Path(path_to_directory/'Merged_product_{}_living.html'.format(product)))
        fig = standardization_report_func(df_summary=df_all.query('merged_instruments.str.contains("_")').groupby(['Sample', 'Latitude', 'Longitude', 'min_depth', 'max_depth'], dropna=True).apply(lambda x: pd.Series({'Abundance': x[dict_products_y[product][0]+'_nonliving'].astype(float).sum(),  # individuals per liter
                                                                                                                                                                                     'Average_diameter': np.nanmean( x[dict_products_y[product][0]+'_nonliving'].astype(float)* x.ECD_mid.astype(float)),# micrometer
                                                                                                                                                                                     'Std_diameter': np.nanstd( x[dict_products_y[product][0]+'_nonliving'].astype(float) * x.ECD_mid.astype(float))})).reset_index(),
            df_standardized=pd.DataFrame({}), df_nbss=df_all.rename( columns={dict_products_y[product][0]+'_nonliving': 'NBSS', 'ECD_mid': 'size_class_mid'}),plot='nbss')
        fig.write_html(Path(path_to_directory/'Merged_product_{}_nonliving.html'.format(product)))
        # Select only regression statistics to re-create the ESSD paper figures, append Global Oceans and Seas regions and save
        df_stats_summary=pd.merge(df_all[['Latitude','Longitude','Sample','year','month','merged_instruments','Group_index','Group_station_index']].drop_duplicates(),df_stats[['slope','intercept','R2']].drop_duplicates(),how='left',left_index=True, right_index=True)
        gdf = gpd.GeoDataFrame(df_stats_summary[['Group_station_index', 'Longitude', 'Latitude']].drop_duplicates().dropna(),geometry=gpd.points_from_xy(df_stats_summary[['Group_station_index', 'Longitude', 'Latitude']].drop_duplicates().dropna().Longitude, df_stats_summary[['Group_station_index', 'Longitude', 'Latitude']].drop_duplicates().dropna().Latitude))
        df_stats_summary['Study_area'] = pd.merge(df_stats_summary,gpd.tools.sjoin(gdf, oceans, predicate="within", how='left')[['Group_station_index', 'name']], how='left', on='Group_station_index')['name'].astype(str)
        #df_stats_summary.to_csv(path_to_directory/'Merged_product_{}_stats_living.csv'.format(product),index=False)
        # Plot the selection for each PFT
        df_pivot=df_stats_summary.pivot_table(values=['slope','intercept', 'R2'],columns=['merged_instruments'],index=['Latitude', 'Longitude', 'Sample', 'year', 'month', 'Group_station_index']).reset_index()
        df_pivot.columns=['_'.join(col) if (type(col) is tuple) & (col[1]!='') else col[0] for col in df_pivot.columns.values]

        plot = (ggplot(data=df_pivot) +
                coord_fixed()+
                geom_point(aes(x="slope_IFCB", y="slope_UVP",group='Sample'), size=1.5, shape='H', color='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50), alpha=0.1) +
                geom_abline(slope=1,intercept=0)+
                labs(x=r'IFCB spectral slope (L$^{-1}$ $\mu$m$^{-3}$)', y=r'UVP spectral slope (L$^{-1}$ $\mu$m$^{3}$)') +
                #scale_x_log10(breaks=[size  for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))],labels=lambda l: [int(size) if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in l])+
                theme_paper).draw(show=True)
        plot.set_size_inches(3,3)
        plot.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/PSSdb_Merged-products_slopes_comparison.svg'.format(str(Path.home())),limitsize=False, dpi=600)
        plot = (ggplot(data=df_pivot) +
                coord_equal() +
                geom_point(aes(x="10**intercept_IFCB", y="10**intercept_UVP", group='Sample'), size=1.5, shape='H',color='#{:02x}{:02x}{:02x}{:02x}'.format(0, 0, 0, 50), alpha=0.1) +
                geom_abline(slope=1, intercept=0) +
                labs(x=r'IFCB spectral intercept ($\mu$m$^{3}$ L$^{-1}$ $\mu$m$^{-3}$)',  y=r'UVP spectral intercept ($\mu$m$^{3}$ L$^{-1}$ $\mu$m$^{-3}$)') +
                scale_y_log10(breaks=[size for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))], labels=lambda l: [int(size) if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else ''  for size in l]) +
                scale_x_log10(breaks=[size  for size in np.sort( np.concatenate(np.arange(1, 10).reshape((9, 1)) * np.power(10, np.arange(1, 5, 1))))],labels=lambda l: [int(size) if (size / np.power(10, np.ceil(np.log10(size)))) == 1 else '' for size in l])+
                theme_paper).draw(show=True)
        plot.set_size_inches(6, 3)
        plot.savefig(fname='{}/GIT/PSSdb/figures/first_datapaper/PSSdb_Merged-products_intercepts_comparison.svg'.format( str(Path.home())), limitsize=False, dpi=600)

        # Update progress bar and move to next product
        ok = bar.update(n=1)



#1) calculate NBSS products no adjustments, and without considering plankton functional groups
merged_raw_list = merge_taxa_products(grouping_factors= ['date_bin', 'Station_location', 'midLatBin', 'midLonBin'])





#2) merge PFT products calculated with different size metrics and calculate
merged_adjusted_raw_list = merge_adjust_taxa_products(grouping_factors= ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'PFT'])


#2) Calculate NBSS for each of these products


