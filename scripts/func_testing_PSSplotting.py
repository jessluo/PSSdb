

# follow these to test functionality (this is to figure out why we are getting such horrible slopes)
path_to_data, ID = proj_id_list('IFCB', testing=True)
data_clean_all=pd.DataFrame()
for i in ID:
    data_clean = read_clean(path_to_data, i)
    data_clean_all = pd.concat([data_clean_all, data_clean], axis=0)

data_clean_all = biovol(data_clean_all, 'IFCB')

# this bit was for generating random numbers and see if the equations were working properly
#data_clean_all['biovol_um3']=np.random.poisson(data_clean_all['biovol_um3'].mode(),
                                               #data_clean_all['biovol_um3'].count())

data_clean_all = binning(data_clean_all)
stats_biovol_SC = NB_SS(data_clean_all)

#plotting of biovolume and NBSS here. this is to see how the data looks like

#first, all data
plt.plot(stats_biovol_SC['logSize'], stats_biovol_SC['logNBSS']) #def not linear :O

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
