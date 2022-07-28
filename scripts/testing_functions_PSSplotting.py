## Plotting scripts


#1 Plot in one graph the data coming from 2 instruments. for this to work, merge needs to happen BEFORE calculating NBSS


#1) Check the (i) number of beads (ii)  and number of bubbles per time in



#plotting of biovolume and NBSS here. this is to see how the data looks like

#first, all data
plt.plot(stats_biovol_SC['logSize'], stats_biovol_SC['logNBSS']) #def not linear :O

# NOTE: this can be useful to detect stations with problematic data
#now, parse by project and by stations:
stats_biovol_SC['proj_ID_min'] = stats_biovol_SC['proj_ID_min'].astype('int')
stats_biovol_SC['Project_ID'] = stats_biovol_SC['Project_ID'].astype('category')
ncols = 3
nrows = len(stats_biovol_SC.Project_ID.unique()) // ncols + (len(stats_biovol_SC.Project_ID.unique()) % ncols > 0)


for n, t in enumerate(stats_biovol_SC.Project_ID.unique()):
    ax = plt.subplot(nrows, ncols, n+1)
    data_plotting = stats_biovol_SC[(stats_biovol_SC['Project_ID'] == t)]
    for i in data_plotting.Station_location.unique():
        plt.plot(
            data_plotting.logSize[data_plotting['Station_location'] == i],
            data_plotting.logNBSS[data_plotting['Station_location'] == i])
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

#use these to define names of columns and test each step of the functions
station = 'Station_ID'
depths = 'midDepthBin'
size_range = 'range_size_bin'
sizeClasses = 'sizeClasses'
biovolume = 'biovol_um3'
lat = 'midLatBin'
lon = 'midLonBin'
project_ID = 'proj_ID'


# histrogram to show the number of samples per size bin in all of the IFCB projects
plt.hist(x=stats_biovol_SC['biovol_um3_count'], bins= 5000)
plt.xlim(0, 500)
plt.xlabel(' number of samples in a Size Bin (including artifacts)')
plt.ylabel('frequency')


## Testing of standardization and merge of data from multiple instruments
# for a unified workflow, this should happen after running the script ecotaxa_API_export_all_data
import funcs_size_spectra

#standardize all downloaded files by instrument:
for i in ['IFCB', 'Zooscan']:
    mass_st_func(i)
    #bin all loaded files and calculate NBSS
    binning_all_func(i)
    #calculate NBSS for each individual file



