
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from glob import glob
try:
    from funcs_NBS import *
except:
    from scripts.funcs_NBS import *

path_to_git=Path('~/GIT/PSSdb').expanduser()
path_to_config = path_to_git /'scripts'/'configuration_masterfile.yaml'
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

def merge_taxa_products (grouping_factors= ['date_bin', 'Station_location', 'midLatBin', 'midLonBin']):
    """
    objective: read Scanner and UVP 1a products with taxonomic information, sum all biomass/biovolume by size to obtain a merged NBSS without adjustments.
    """
    currentMonth = str(datetime.datetime.now().month).rjust(2, '0')
    currentYear = str(datetime.datetime.now().year)
    file_list = [x for x in glob(str(Path(cfg['raw_dir']).expanduser() / 'NBSS_data' / 'PFT') + '/*distribution_all_var*.csv', recursive=True) if 'IFCB' not in x]
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


def merge_adjust_taxa_products (grouping_factors= ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'PFT', 'sizeClasses']):
    """
    objective: read Scanner and UVP 1a products with taxonomic information, correct and merge products based on Soviadan et al.
    """
    currentMonth = str(datetime.datetime.now().month).rjust(2, '0')
    currentYear = str(datetime.datetime.now().year)
    file_list = [x for x  in glob(str(Path(cfg['raw_dir']).expanduser() / 'NBSS_data' / 'PFT') + '/*distribution_all_var*.csv',recursive=True) if 'IFCB' not in x]
    merged_prod_path = Path(cfg['raw_dir']).expanduser() / 'NBSS_data' / 'PFT'/ 'merged_products'
    if not os.path.exists(merged_prod_path):
        os.mkdir(merged_prod_path)
    for s in ['Biomass', 'Size', 'Weight']:
        UVP_df = pd.read_csv([x for x in file_list if 'UVP' in x and s in x][0], sep = ',')
        print('UVP ' + s+ 'length is ' + str(len(UVP_df)))
        Scanner_df = pd.read_csv([x for x in file_list if 'Scanner' in x and s in x][0], sep = ',')
        print('Scanner ' + s + 'length is ' + str(len(Scanner_df)))
        #max_nb_list = pd.concat([UVP_df, Scanner_df]).groupby(grouping_factors, sort = False).apply(lambda x : pd.Series({'max_NB': np.nanmax(x.NB)})) # this only selects
        merged_df_test = pd.merge(UVP_df, Scanner_df, how= 'outer',  on=grouping_factors)
        merged_df_final = pd.DataFrame()
        print ('merging UVP and Scanner NBSS based on ' + s)
        for i in tqdm(range(0,len(merged_df_test))):
            for c in merged_df_test.columns:
                if np.nanmax([merged_df_test.NB_x[i], merged_df_test.NB_y[i]]) == merged_df_test[c].iloc[i]:
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
