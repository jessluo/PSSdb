
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from glob import glob
try:
    from funcs_NBS import *
    from script_NBSS_datacheck import regress_nbss
except:
    from scripts.funcs_NBS import *
    from scripts.script_NBSS_datacheck import regress_nbss

path_to_git=Path('~/GIT/PSSdb').expanduser()
path_to_config = path_to_git /'scripts'/'configuration_masterfile.yaml'
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
from plotnine import *
theme_paper=theme(axis_ticks_direction="inout",
              panel_grid=element_blank(),
              axis_line = element_line(colour = "black"),
              panel_background=element_rect(fill='white'),
              panel_border=element_rect(color='black'),
              legend_title=element_text(family="serif", size=8),
              legend_position='top',
              legend_text=element_text(family="serif", size=8),
              strip_text=element_text(family="serif", size=8),
              axis_title=element_text(family="serif", size=8),
              axis_text_x=element_text(family="serif", size=8),
              axis_text_y=element_text(family="serif", size=8, rotation=90),
              plot_background=element_rect(fill='white'), strip_background=element_rect(fill='white'))
pal_instrument = { 'Scanner': '#6f918aff', 'UVP': 'red','IFCB': '#6f0058ff'}
import seaborn as sns

sns.set_context("paper",rc={"axes.labelsize": 8, 'axes.labelfamily': 'serif', "tick.labelsize": 8, 'tick.labelfamily': 'serif'})

def seaborn_plot_unity(xdata, ydata, **kwargs):
    mn = min(xdata.min(), ydata.min())
    mx = max(xdata.max(), ydata.max())
    points = np.linspace(mn, mx, 100)
    plt.gca().plot(points, points, color='k', marker=None,
            linestyle='--', linewidth=1.0)
from scipy.stats import linregress
def seaborn_r2(x, y,xy,color, ax=None, **kws):
    ax = ax or plt.gca()
    data=pd.DataFrame({'x':x,'y':y}).dropna(subset=['x','y'])
    if len(data):
        slope, intercept, r_value, p_value, std_err = linregress(x=data.x, y=data.y)
        ax.annotate(f'$r^2 = {r_value ** 2:.2f}$\nEq: ${slope:.2f}x{intercept:+.2f}$',
                    xy=xy, xycoords=ax.transAxes, fontsize=8,
                    color=color, backgroundcolor='#FFFFFF99', ha='left', va='top')

def seaborn_diag_func(data, label, color):
    for val in data.quantile([ .5]):
        plt.axvline(val, ls=':', color=color)
        plt.title(str(np.round(val,2)),size=8,family='serif', color=color,rotation=-90,y=0)

def seaborn_hide_axis(*args, **kwds):
    plt.gca().set_visible(False)


dtypes_products={'PFT':str,'Validation_percentage':float,'ocean':str,'year':str,'month':str,'n':float,'min_depth':float,'max_depth':float,'biomass_mid':float,'range_biomass_bin':float,'normalized_biomass_mean':float,'normalized_biomass_std':float,'biovolume_size_class':float,'normalized_biovolume_mean':float,'normalized_biovolume_std':float,'equivalent_circular_diameter_mean':float,'normalized_abundance_mean':float,'normalized_abundance_std':float}
dict_products_y={'Size':['NB','PSD'],'Biomass':['NB'],'Weight':['NB']}#{'Size':['normalized_biovolume_mean','normalized_abundance_mean'],'Biomass':'normalized_biomass_mean','Weight':'normalized_weight_mean'}
dict_products_x={'Size':'size_class_mid','Biomass':'biomass_mid','Weight':'biomass_mid'}#{'Size':'biovolume_size_class','Biomass':'biomass_mid','Weight':'biomass_mid'}

import geopandas as gpd
oceans = gpd.read_file(list(Path(gpd.datasets.get_path("naturalearth_lowres")).expanduser().parent.parent.rglob('goas_v01.shp'))[0])

# Append biomass bins
from scipy import stats
from natsort import natsorted
bins=np.power(2, np.arange(0, np.log2(100000) + 1 / 3, 1 / 3))  # Fixed size(um) bins used for UVP/EcoPart data. See https://ecopart.obs-vlfr.fr/. 1/3 ESD increments allow to bin particle of doubled biovolume in consecutive bins. np.exp(np.diff(np.log((1/6)*np.pi*(EcoPart_extended_bins**3))))
Ecopart_bins=pd.read_csv(Path(cfg['size_bins_path']).expanduser())
bins=Ecopart_bins.ESD_um.to_numpy()
path_to_allometry=path_to_git/cfg['allometry_lookup']

df_allometry=pd.read_excel(path_to_allometry)
df_allometry=df_allometry[df_allometry.Size_proxy=='Biovolume']

biomass_bins=(1e-6*(df_allometry.loc[df_allometry.Taxon=='Living','C_Intercept'].values[0]*(1e-09*(1 / 6) * np.pi * bins ** 3)**(df_allometry.loc[df_allometry.Taxon=='Living','C_Slope'].values[0])))#(1e-12*(df_allometry.loc[df_allometry.Taxon=='Protozoa','C_Intercept'].values[0]*((1 / 6) * np.pi * bins ** 3)**(df_allometry.loc[df_allometry.Taxon=='Protozoa','C_Slope'].values[0])))
df_bins = pd.DataFrame({'sizeClasses': pd.cut(Ecopart_bins.biovol_um3.to_numpy(), Ecopart_bins.biovol_um3.to_numpy()).categories.values.astype(str),  # Define bin categories (cubic micrometers)
                        'sizeClasses_ECD': pd.cut(bins, bins).categories.values,  # Define bin categories (um)
                        'size_range_ECD': np.diff(bins),  # Define width of individual bin categories (um)
                        'range_size_bin': np.concatenate(np.diff((1 / 6) * np.pi * (np.resize(np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])), (len(bins) - 1, 2)) ** 3), axis=1)),  # cubic micrometers
                        'ECD_mid': stats.gmean(np.resize( np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])), (len(bins) - 1, 2)), axis=1), # Define geometrical mean of bin categories (um)
                        'biovolume_size_class':stats.gmean((1 / 6) * np.pi * (np.resize(np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])), (len(bins) - 1, 2)) ** 3), axis=1), # Define geometrical mean of bin categories (cubic micrometer)
                        'size_class_mid':stats.gmean((1 / 6) * np.pi * (np.resize(np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])), (len(bins) - 1, 2)) ** 3), axis=1), # Define geometrical mean of bin categories (cubic micrometer)
                        'range_biomass_bin':np.concatenate(np.diff(np.resize(np.append(biomass_bins[0], np.append(np.repeat(biomass_bins[1:-1], repeats=2), biomass_bins[len(biomass_bins) - 1])), (len(biomass_bins) - 1, 2)), axis=1)), # in g
                        'biomass_mid':stats.gmean(np.resize( np.append(biomass_bins[0], np.append(np.repeat(biomass_bins[1:-1], repeats=2), biomass_bins[len(biomass_bins) - 1])), (len(biomass_bins) - 1, 2)), axis=1) #in g
})
data_bins=pd.DataFrame({'biovolume_size_class':(1 / 6) * np.pi *(bins**3),'size_class_mid':(1 / 6) * np.pi *(bins**3),'biomass_mid':biomass_bins})
def merge_products(path_to_products,grouping_factors= ['ocean','year','month','latitude','longitude']):
    df_nbss = pd.concat(map(lambda path: pd.read_csv(path, dtype=dtypes_products).assign(Instrument=path.name.split("_")[0]), path_to_products)).reset_index( drop=True)  # reduce(lambda left,right:pd.merge(left,right,how='outer',on=['PFT','min_depth','max_depth','Instrument','n']+grouping_factor),list(map(lambda path: pd.read_csv(path,dtype=dtypes_products).drop(columns=['Validation_percentage','QC_3std_dev','QC_min_n_size_bins','QC_R2']).assign(Instrument=path.name.split("_")[0]) if path.name.split("_")[2]!='Weight-distribution' else pd.read_csv(path,dtype=dtypes_products).rename(columns={'normalized_biomass_mean':'normalized_weight_mean','normalized_biomass_std':'normalized_weight_std'}).drop(columns=['Validation_percentage','QC_3std_dev','QC_min_n_size_bins','QC_R2']).assign(Instrument=path.name.split("_")[0]),path_to_products)) ).reset_index(drop=True)

    #df_nbss = pd.concat(map(lambda path: pd.read_csv(path, dtype=dtypes_products).drop( columns=['Validation_percentage', 'QC_3std_dev', 'QC_min_n_size_bins', 'QC_R2']).assign(Instrument=path.name.split("_")[0]) if path.name.split("_")[2] != 'Weight-distribution' else pd.read_csv(path,dtype=dtypes_products).rename( columns={'normalized_biomass_mean': 'normalized_weight_mean', 'normalized_biomass_std': 'normalized_weight_std'}).drop( columns=['Validation_percentage', 'QC_3std_dev', 'QC_min_n_size_bins', 'QC_R2']).assign(Instrument=path.name.split("_")[0]), path_to_products)).reset_index(drop=True)  # reduce(lambda left,right:pd.merge(left,right,how='outer',on=['PFT','min_depth','max_depth','Instrument','n']+grouping_factor),list(map(lambda path: pd.read_csv(path,dtype=dtypes_products).drop(columns=['Validation_percentage','QC_3std_dev','QC_min_n_size_bins','QC_R2']).assign(Instrument=path.name.split("_")[0]) if path.name.split("_")[2]!='Weight-distribution' else pd.read_csv(path,dtype=dtypes_products).rename(columns={'normalized_biomass_mean':'normalized_weight_mean','normalized_biomass_std':'normalized_weight_std'}).drop(columns=['Validation_percentage','QC_3std_dev','QC_min_n_size_bins','QC_R2']).assign(Instrument=path.name.split("_")[0]),path_to_products)) ).reset_index(drop=True)
    # Determine the bins with multiple instruments
    df_nbss = pd.merge(df_nbss.astype(dict(zip(grouping_factors, [str] * len(grouping_factors)))), df_nbss.astype(dict(zip(grouping_factors, [str] * len(grouping_factors)))).groupby( grouping_factors).apply(lambda x: pd.Series({'n_instruments': x.Instrument.nunique(),'merged_instruments': "_".join(natsorted( x.Instrument.unique()))})).reset_index(), how='left', on=grouping_factors)

    return df_nbss

def merge_taxa_products (path_to_directory,grouping_factors= ['date_bin', 'Station_location', 'midLatBin', 'midLonBin']):
    """
    Objective: read Scanner and UVP 1a products with taxonomic information, sum all biomass/biovolume by size to obtain a merged NBSS without adjustments.
    """
    file_list = [x for x in glob(str(Path(path_to_directory).expanduser()) + '/*distribution_all_var*.csv', recursive=True) if 'IFCB' not in x]
    merged_prod_path = Path(cfg['raw_dir']).expanduser() / 'NBSS_data' / 'PFT' / 'merged_products'
    if not os.path.exists(merged_prod_path):
        os.mkdir(merged_prod_path)
    for s in ['Biomass', 'Size', 'Weight']:
        UVP_df = pd.read_csv([x for x in file_list if 'UVP' in x and s in x][0], sep=',')
        print('UVP ' + s + 'length is ' + str(len(UVP_df)))
        Scanner_df = pd.read_csv([x for x in file_list if 'Scanner' in x and s in x][0], sep=',')
        print('Scanner ' + s + 'length is ' + str(len(Scanner_df)))
        df_merged = pd.concat([UVP_df, Scanner_df]).reset_index(drop=True)
        if s == 'Biomass':
            group_factors = grouping_factors +['biomassClasses']
            df_merged = df_merged.astype(dict(zip(['midLatBin', 'midLonBin','biomassClasses'], [str]*3))).groupby(group_factors).apply(lambda x: pd.Series({'Min_obs_depth': np.nanmin(x.Min_obs_depth),
                                                                                                                                    'Max_obs_depth': np.nanmax(x.Max_obs_depth),
                                                                                                                                    'biomass_mid':x.biomass_mid.unique()[0],
                                                                                                                                    'range_biomass_bin': x.range_biomass_bin.unique()[0],
                                                                                                                                    #'size_range_ECD': x.size_range_ECD.unique()[0],
                                                                                                                                    'biomass_mid': x.biomass_mid.unique()[0],
                                                                                                                                    'range_biomass_bin': x.range_biomass_bin.unique()[0],
                                                                                                                                    'Total_biomass': np.nansum(x.Total_biomass),
                                                                                                                                    'Biomass_mean': np.nanmean(x.Biomass_mean),
                                                                                                                                    #'size_class_pixel': np.nanmean(x.size_class_pixel),
                                                                                                                                    'ROI_number_sum':np.nansum(x.ROI_number_sum),
                                                                                                                                    'ROI_abundance_mean': np.nanmean(x.ROI_abundance_mean),
                                                                                                                                    'NB':np.nansum(x.NB),
                                                                                                                                    #'PSD':np.nansum(x.PSD),
                                                                                                                                    'count_uncertainty': np.nanmax(x.count_uncertainty),
                                                                                                                                    'size_uncertainty': np.nanmax(x.size_uncertainty)})).reset_index()
            del group_factors
            df_merged['logNB'] = np.log10(df_merged['NB'])
            df_merged['logSize'] = np.log10(df_merged['biomass_mid'].astype(float))
        elif s == 'Size':
            group_factors = grouping_factors + ['sizeClasses']
            df_merged = pd.concat([UVP_df, Scanner_df]).reset_index().groupby(group_factors).apply(lambda x: pd.Series({'Min_obs_depth': np.nanmin(x.Min_obs_depth),
                                                                                                                                    'Max_obs_depth': np.nanmax(x.Max_obs_depth),
                                                                                                                                    'size_class_mid': x.size_class_mid.unique()[0],
                                                                                                                                    'range_size_bin': x.range_size_bin.unique()[0],
                                                                                                                                    'ECD_mid': x.ECD_mid.unique()[0],
                                                                                                                                    'size_range_ECD':x.size_range_ECD.unique()[0],
                                                                                                                                    'Biovolume_mean': np.nanmean(x.Biovolume_mean),
                                                                                                                                    'size_class_pixel': np.nanmean(x.size_class_pixel),
                                                                                                                                    'ROI_number_sum': np.nansum(x.ROI_number_sum),
                                                                                                                                    'ROI_abundance_mean': np.nanmean(x.ROI_abundance_mean),
                                                                                                                                    'NB': np.nansum(x.NB),
                                                                                                                                    'PSD':np.nansum(x.PSD),
                                                                                                                                    'count_uncertainty': np.nanmax(x.count_uncertainty),
                                                                                                                                    'size_uncertainty': np.nanmax(x.size_uncertainty)})).reset_index()

            del group_factors
            df_merged['logNB'] = np.log10(df_merged['NB'])
            df_merged['logSize'] = np.log10(df_merged['size_class_mid'].astype(float))
            df_merged['logPSD'] = np.log10(df_merged['PSD'])
            df_merged['logECD'] = np.log10(df_merged['ECD_mid'].astype(float))

        elif s =='Weight':
            group_factors = grouping_factors + ['sizeClasses']
            df_merged = pd.concat([UVP_df, Scanner_df]).reset_index().groupby(group_factors).apply(lambda x: pd.Series({'Min_obs_depth': np.nanmin(x.Min_obs_depth),
                                                                                                                                    'Max_obs_depth': np.nanmax(x.Max_obs_depth),
                                                                                                                                    'size_class_mid':x.size_class_mid.unique()[0],
                                                                                                                                    'range_size_bin': x.range_size_bin.unique()[0],
                                                                                                                                    'ECD_mid': x.ECD_mid.unique()[0],
                                                                                                                                    'size_range_ECD': x.size_range_ECD.unique()[0],
                                                                                                                                    'biomass_mid': x.biomass_mid.unique()[0],
                                                                                                                                    'range_biomass_bin': x.biomass_mid.unique()[0],
                                                                                                                                    'Total_biomass': np.nansum(x.Total_biomass),
                                                                                                                                    'size_class_pixel': np.nanmean(x.size_class_pixel),
                                                                                                                                    'ROI_number_sum':np.nansum(x.ROI_number_sum),
                                                                                                                                    'ROI_abundance_mean': np.nanmean(x.ROI_abundance_mean),
                                                                                                                                    'NB':np.nansum(x.NB),
                                                                                                                                    'PSD':np.nansum(x.PSD),
                                                                                                                                    'count_uncertainty': np.nanmax(x.count_uncertainty),
                                                                                                                                    'size_uncertainty':np.nanmax(x.size_uncertainty)})).reset_index()
            del group_factors
            df_merged['logNB'] = np.log10(df_merged['NB'])
            df_merged['logSize'] = np.log10(df_merged['size_class_mid'].astype(float))
            df_merged['logPSD'] = np.log10(df_merged['PSD'])
            df_merged['logECD'] = np.log10(df_merged['ECD_mid'].astype(float))

        df_merged.to_csv(str(merged_prod_path) + '/' + s + '_Merged_Size-distribution_all_var_v' + currentYear + '-' + currentMonth + '.csv', index=False)
    merged_raw_list = glob(str(merged_prod_path) + '/*_Merged_*')
    return merged_raw_list


def merge_adjust_taxa_products (path_to_directory,grouping_factors= ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'PFT', 'sizeClasses']):
    """
    Objective: read Scanner and UVP 1a products with taxonomic information, correct and merge products based on Soviadan et al.
    """
    file_list = [x for x  in glob(str(Path(path_to_directory).expanduser()) + '/*distribution_all_var*.csv',recursive=True) if 'IFCB' not in x]
    merged_prod_path = Path(cfg['raw_dir']).expanduser() / 'NBSS_data' / 'PFT'/ 'merged_products'
    if not os.path.exists(merged_prod_path):
        os.mkdir(merged_prod_path)
    for s in ['Biomass', 'Size', 'Weight']:
        UVP_df = pd.read_csv([x for x in file_list if 'UVP' in x and s in x][0], sep = ',')
        print('UVP ' + s+ 'length is ' + str(len(UVP_df)))
        Scanner_df = pd.read_csv([x for x in file_list if 'Scanner' in x and s in x][0], sep = ',')
        print('Scanner ' + s + 'length is ' + str(len(Scanner_df)))
        #remove rows from these datasets that have a particle count <5
        #max_nb_list = pd.concat([UVP_df, Scanner_df]).groupby(grouping_factors, sort = False).apply(lambda x : pd.Series({'max_NB': np.nanmax(x.NB)})) # this only selects
        merged_df_test = pd.merge(UVP_df, Scanner_df, how= 'inner',  on=grouping_factors) # merged products will contain data only for grid cells with both instruments, thus how = inner
        merged_df_final = pd.DataFrame()
        print ('merging UVP and Scanner NBSS based on ' + s)
        for i in tqdm(range(0,len(merged_df_test))):
            for c in merged_df_test.columns:
                if np.nanmax([merged_df_test.NB_x[i], merged_df_test.NB_y[i]]) == merged_df_test[c].iloc[i]: # replace x and y with UVP and Scanner
                    col_name = c
            str_subset_cols = '_'+col_name.split('_')[-1]
            col_names = grouping_factors + [c for c in merged_df_test.columns if str_subset_cols in c]
            d={}
            for k in col_names:
                d[k]=merged_df_test[k].iloc[i]
            row_val =pd.DataFrame(d, index=[0])
            row_val.columns = row_val.columns.str.removesuffix(str_subset_cols)
            merged_df_final = pd.concat([merged_df_final, row_val])
        merged_df_final.to_csv(str(merged_prod_path) + '/' + s + '_Merged-Adjusted_Size-distribution_all_var_v' + currentYear + '-' + currentMonth + '.csv',index=False)
    merged_raw_list = glob(str(merged_prod_path) +'/*Merged-Adjusted*')
    return merged_raw_list


"""


IFCB_dict={}
#date_IFCB = []
for i in df_list_IFCB:
    df = pd.read_csv(i, sep='\t', header=0, index_col=[0])
    df['instrument'] = 'IFCB'
    st_list = list(set(df['Station_location']))
    date_list = list(set(df['date_bin']))
    #date_list = list(set(df['date_bin'])) # next three lines are to query the dates
    #date_IFCB.append([n for n in set(df['date_bin'])])
    #dates_IFCB = [val for sublist in date_IFCB for val in sublist]
    for st in st_list:
        IFCB_dict[st] = {}
        for d in date_list:
            IFCB_dict[st][d] = df[(df['Station_location']==st) & (df['date_bin']==d)]


Zooscan_dict={}
#date_Zooscan = []
for i in df_list_Zooscan:
    df = pd.read_csv(i, sep='\t', header=0, index_col=[0])
    df['instrument'] = 'Zooscan'
    st_list = list(set(df['Station_location']))
    date_list = list(set(df['date_bin']))
    #date_list = list(set(df['date_bin'])) # next three lines are to query the dates
    #date_Zooscan.append([n for n in set(df['date_bin'])])
    #dates_Zooscan = [val for sublist in date_Zooscan for val in sublist]
    for st in st_list:
        Zooscan_dict[st] = {}
        for d in date_list:
            Zooscan_dict[st][d] = df[(df['Station_location']==st) & (df['date_bin']==d)]

UVP_dict={}
#date_UVP = []
for i in df_list_UVP:
    df = pd.read_csv(i, sep='\t', header=0, index_col=[0])
    df['instrument'] = 'UVP'
    st_list = list(set(df['Station_location']))
    date_list = list(set(df['date_bin']))
    #date_list = list(set(df['date_bin'])) # next three lines are to query the dates
    #date_UVP.append([n for n in set(df['date_bin'])])
    #dates_UVP = [val for sublist in date_UVP for val in sublist]
    for st in st_list:
        UVP_dict[st] = {}
        for d in date_list:
            UVP_dict[st][d] = df[(df['Station_location']==st) & (df['date_bin']==d)]



# merge data based on location. 1/24/2023: there are no sites that have all three instruments. Only Zooscan + UVP and Zooscan + IFCB
merged_dfs={}
for i in Zooscan_dict:
    try:
        merged_st_df =  pd.concat([UVP_dict[i], Zooscan_dict[i]],  axis=0) #UVP_dict[i]],
        merged_st_df = merged_st_df.reset_index(drop=True)
        merged_dfs[i] = merged_st_df
    except:
        pass

# this is where I get the common stations across instruments
IFCB_keys = set(IFCB_dict.keys())
Zooscan_keys = set(Zooscan_dict.keys())
UVP_keys = set(UVP_dict.keys())
set1 = Zooscan_keys.intersection(IFCB_keys)
set2 = Zooscan_keys.intersection(UVP_keys)
common = set1.intersection(set2)
common


fig, ax1 = plt.subplots(1,1)
for i in merged_dfs.keys(): #['72.5_44.5', '78.5_79.5', '70.5_-53.5']
    ax1.plot(merged_dfs[i]['logSize'], merged_dfs[i]['logNBSS'], marker='.', label=merged_dfs[i]['instrument'])
ax1.set_ylabel('log (NBS)')
ax1.set_xlabel('log (Biovolume)')
ax1.legend()

# with seaborn:
for i in merged_dfs.keys(): #['72.5_44.5', '78.5_79.5', '70.5_-53.5']
    sns.lineplot(data = merged_dfs[i], x = 'logSize', y = 'logNBSS', marker='.', hue='instrument')

# separate instruments:

for st in IFCB_dict.keys(): #['72.5_44.5', '78.5_79.5', '70.5_-53.5']
    for d in IFCB_dict[st].keys():
        sns_plot = sns.lineplot(data = IFCB_dict[st][d], x = 'logSize', y = 'logNBSS', marker='.', color='red', linewidth = 0.5)

for st in Zooscan_dict.keys(): #['72.5_44.5', '78.5_79.5', '70.5_-53.5']
    for d in Zooscan_dict[st].keys():
        sns_plot = sns.lineplot(data = Zooscan_dict[st][d], x = 'logSize', y = 'logNBSS', marker='.', color='cyan', linewidth = 0.5)

for st in UVP_dict.keys(): #['72.5_44.5', '78.5_79.5', '70.5_-53.5']
    for d in UVP_dict[st].keys():
        sns_plot = sns.lineplot(data = UVP_dict[st][d], x = 'logSize', y = 'logNBSS', marker='.', color='green', linewidth = 0.5)

intercept = []
slope = []
r2 = []
lat = []
lon = []
fig, ax1 = plt.subplots(1,1)
for i in merged_dfs: #['72.5_44.5', '78.5_79.5', '70.5_-53.5']
    ax1.plot(merged_dfs[i]['logSize'], merged_dfs[i]['logNBSS'], marker='.', label=i)
    lr = LinearRegression()
    imp = SimpleImputer(missing_values=np.nan, strategy='most_frequent')
    x = merged_dfs[i]['logSize'].values.reshape(-1, 1)
    y = merged_dfs[i]['logNBSS'].values.reshape(-1, 1)
    comb = np.dstack((x, y)) ##  THIS IS NOT OK, we are filling up size bins with no data
    comb = np.squeeze(comb, axis=1)
    imp.fit(comb)
    comb = imp.transform(comb)
    x = comb[:, 0].reshape(-1, 1)
    y = comb[:, 1].reshape(-1, 1)
    lr.fit(x, y)
    ly_pred = lr.predict(x).reshape(-1, 1)
    intercept.append(list(lr.intercept_)) item for sublist in l for item in sublist
    slope.append(list(lr.coef_))
    r2.append(r2_score(y, ly_pred))
    lat.append(merged_dfs[i]['midLatBin'][0])
    lon.append(merged_dfs[i]['midLonBin'][0])
    #ax1.plot(x, ly_pred, color="black", linestyle = 'dashed')
intercept = [item for sublist in intercept for item in sublist]
slope = [item for sublist in slope for item in sublist]
ax1.set_ylabel('log (NBS)', fontsize=18)
ax1.set_xlabel('log (Biovolume)', fontsize=18)
ax1.autoscale(tight=True)
#ax1.legend()
plt.savefig('plot_all_NBSS_merged.pdf', dpi= 300)


# global map goes here


intercept_plot = [x*3 for x in intercept]
ax = plt.axes(projection=ccrs.PlateCarree())

plt.gca().coastlines('50m')
g1=ax.gridlines(draw_labels=True)
g1.xlines = False
g1.ylines = False


plt.scatter(lon, lat, label=None, c=slope, cmap='viridis', s=intercept_plot, linewidth=0, alpha=0.5, transform=ccrs.PlateCarree())
plt.gca().coastlines('50m')
ax.set_extent([-180, 180, -90, 90])
#ax.xlabel('longitude')
#ax.ylabel('latitude')
plt.colorbar(label='slope', orientation='horizontal', anchor=(0.5, 1))
#ax.clim(min(slope), max(slope))


labels = ["6", "9.5", "24.5"]

for n, area in enumerate([18, 28.5, 73.5]):
    plt.scatter([], [], c='k', alpha=0.3, s=area, label=labels[n], transform=ccrs.PlateCarree())

plt.legend(bbox_to_anchor=(0.75, 0), ncol = 3, scatterpoints=1, frameon=False,
           labelspacing=1, title='intercept')

"""