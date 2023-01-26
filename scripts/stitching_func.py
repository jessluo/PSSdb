import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
from sklearn.impute import SimpleImputer
import cartopy.crs as ccrs
from funcs_read import *
import pandas as pd
import seaborn as sns

path_to_IFCB, df_list_IFCB =  df_list_func('IFCB', data_status = 'NBSS')

path_to_Zooscan, df_list_Zooscan =  df_list_func('Zooscan', data_status = 'NBSS')

path_to_UVP, df_list_UVP =  df_list_func('UVP', data_status = 'NBSS')



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
            IFCB_dict[st][d] = df[df['Station_location']==st & df['date_bin']==d]


Zooscan_dict = {}
date_Zooscan = []
for i in df_list_Zooscan:
    df = pd.read_csv(i, sep='\t', header=0, index_col=[0])
    df['instrument'] = 'Zooscan'
    st_list = list(set(df['Station_location']))
    date_list = list(set(df['date_bin'])) # next three lines are to query the dates
    date_Zooscan.append([n for n in set(df['date_bin'])])
    dates_Zooscan = [val for sublist in date_Zooscan for val in sublist]
    for n in st_list:
        Zooscan_dict[n] = df[df['Station_location'] == n]

UVP_dict = {}
date_UVP = []
for i in df_list_UVP:
    df = pd.read_csv(i, sep='\t', header=0, index_col=[0])
    df['instrument'] = 'UVP'
    st_list = list(set(df['Station_location']))
    date_list = list(set(df['date_bin'])) # next three lines are to query the dates
    date_UVP.append([n for n in set(df['date_bin'])])
    dates_UVP= [val for sublist in date_UVP for val in sublist]
    for n in st_list:
        UVP_dict[n] = df[df['Station_location'] == n]


min(dates_IFCB)
max(dates_IFCB)
min(dates_Zooscan)
max(dates_Zooscan)
min(dates_UVP)
max(dates_UVP)

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

for i in IFCB_dict.keys(): #['72.5_44.5', '78.5_79.5', '70.5_-53.5']
    sns_plot = sns.lineplot(data = IFCB_dict[i], x = 'logSize', y = 'logNBSS', marker='.', color='blue', linewidth = 0.5)

for i in Zooscan_dict.keys(): #['72.5_44.5', '78.5_79.5', '70.5_-53.5']
    sns_plot = sns.lineplot(data = Zooscan_dict[i], x = 'logSize', y = 'logNBSS', marker='.', color='red', linewidth = 0.5)

for i in UVP_dict.keys(): #['72.5_44.5', '78.5_79.5', '70.5_-53.5']
    sns_plot = sns.lineplot(data = UVP_dict[i], x = 'logSize', y = 'logNBSS', marker='.', color='grey', linewidth = 0.5)

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
