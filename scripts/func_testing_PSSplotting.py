

# follow these to test functionality (this is to figure out why we are getting such horrible slopes)
path_to_data, ID = proj_id_list('IFCB', testing=True)
data_clean_all=pd.DataFrame()
for i in ID:
    data_clean = read_clean(path_to_data, i)
    data_clean_all = pd.concat([data_clean_all, data_clean], axis=0)

# function that calls on the standardizer here

data_clean_all['biovol_um3']= biovol_standardizer(data_clean_all['object_biovolume'],
                                     data_clean_all['process_pixel_to_micron'],  'IFCB')

data_clean_all['sizeClasses'], data_clean_all['range_size_bin'] = size_binning(data_clean_all['biovol_um3'])
data_clean_all['midDepthBin'] = depth_binning(data_clean_all['object_depth_min'])

data_clean_all['Station_ID'], data_clean_all['midLatBin'], data_clean_all['midLonBin'] = \
    station_binning(data_clean_all['object_lat'], data_clean_all['object_lon'])
# this bit was for generating random numbers and see if the equations were working properly
#data_clean_all['biovol_um3']=np.random.poisson(data_clean_all['biovol_um3'].mode(),
                                               #data_clean_all['biovol_um3'].count())

stats_biovol_SC = NB_SS(data_clean_all,
          station='Station_ID',
          depths='midDepthBin',
          size_range='range_size_bin',
          sizeClasses='sizeClasses',
          biovolume='biovol_um3',
          lat='midLatBin',
          lon='midLatBin',
          project_ID='proj_ID')

results_SS_all=pd.DataFrame()
for i in ID:
    results_SS = NB_SS_regress(stats_biovol_SC, station='Station_ID',
                               depths='midDepthBin', lat='midLatBin', lon='midLonBin',  ID=i)
    results_SS_all= pd.concat([results_SS_all, results_SS], axis=0)



#plotting of biovolume and NBSS here. this is to see how the data looks like

#first, all data
plt.plot(stats_biovol_SC['logSize'], stats_biovol_SC['logNBSS']) #def not linear :O

# NOTE: this can be useful to detect stations with problematic data
#now, parse by project and by stations:
stats_biovol_SC['proj_ID_min'] = stats_biovol_SC['proj_ID_min'].astype('int')
stats_biovol_SC['proj_ID_min'] = stats_biovol_SC['proj_ID_min'].astype('category')
ncols = 3
nrows = len(stats_biovol_SC.proj_ID_min.unique()) // ncols + (len(stats_biovol_SC.proj_ID_min.unique()) % ncols > 0)

for n, t in enumerate(stats_biovol_SC.proj_ID_min.unique()):
    ax = plt.subplot(nrows, ncols, n+1)
    data_plotting = stats_biovol_SC[(stats_biovol_SC['proj_ID_min'] == t)]
    for i in data_plotting.station_ID.unique():
        plt.plot(
            data_plotting.logSize[data_plotting['station_ID'] == i],
            data_plotting.logNBSS[data_plotting['station_ID'] == i])
    # for each taxonomic category and make a histogram with the counts of each biovolume bin
        ax.set_title(str(t).upper(), fontsize=10)
        ax.set_xlabel('Logbiovolume')
        ax.set_ylabel('LogNBSS')

        #what i notice from these graphs is that some stations have very little data
stats_biovol_SC.biovol_um3_count[(stats_biovol_SC['proj_ID_min']==3318)].min()





results_SS_all=pd.DataFrame()
for i in ID:
    results_SS = NB_SS_regress(stats_biovol_SC, i)
    results_SS_all = pd.concat([results_SS_all , results_SS], axis=0)
