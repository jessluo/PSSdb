import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
from sklearn.impute import SimpleImputer
import cartopy.crs as ccrs

path_to_IFCB, df_list_IFCB =  df_list_func('IFCB', data_status = 'NBSS')

path_to_Zooscan, df_list_Zooscan =  df_list_func('Zooscan', data_status = 'NBSS')

IFCB_dict={}
for i in df_list_IFCB:
    df = pd.read_csv(i, sep='\t', header=0, index_col=[0])
    st_list = list(set(df['Station_location']))
    for n in st_list:
        IFCB_dict[n] = df[df['Station_location']==n]

Zooscan_dict = {}
for i in df_list_Zooscan:
    df = pd.read_csv(i, sep='\t', header=0, index_col=[0])
    st_list = list(set(df['Station_location']))
    for n in st_list:
        Zooscan_dict[n] = df[df['Station_location'] == n]

merged_dfs={}
for i in Zooscan_dict:
    try:
        merged_st_df =  pd.concat([IFCB_dict[i], Zooscan_dict[i]], axis=0)
        merged_st_df = merged_st_df.reset_index(drop=True)
        merged_dfs[i] = merged_st_df
    except:
        pass

fig, ax1 = plt.subplots(1,1)
for i in ['72.5_44.5', '78.5_79.5', '70.5_-53.5']: #['72.5_44.5', '78.5_79.5', '70.5_-53.5']
    ax1.plot(merged_dfs[i]['logSize'], merged_dfs[i]['logNBSS'], marker='.', label=i)
ax1.set_ylabel('log (NBS)')
ax1.set_xlabel('log (Biovolume)')
ax1.legend()

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
