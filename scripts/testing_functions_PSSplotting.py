## Plotting scripts


from funcs_size_spectra import proj_id_list_func
# get the paths where the binned files are:
path_to_binned_IFCB, ID_list_IFCB =  proj_id_list_func('IFCB', standardized='binned', testing=False)

path_to_binned_Zooscan, ID_list_Zooscan =  proj_id_list_func('Zooscan', standardized='binned', testing=False)

# edit the list of projects, this is hardcoded and was obtained comparing dictionaries with the 'draft_funcs_merge'
IFCB_list_lgsample = [3302, 3309, 3321, 3324, 3325, 3334]
IFCB_list_smsample = [3147, 2248, 3289, 3296]
Zooscan_list = [378]

#read the IFCB tsv's, clean them and plot them:
path_to_IFCB = glob(str(path_to_binned_IFCB) + '/' + str(3309) + '*NBSS*') # for 3309, there is an issue in a single point where comparing NBSS between size bins does not work for filtering the shorter side of the spectrum
binned_data_IFCB = pd.read_csv(path_to_IFCB[0], sep='\t', header=0, index_col=[0])
binned_data_IFCB  = clean_lin_fit(binned_data_IFCB, instrument='IFCB', method = 'instrument_specific')
plt.plot(binned_data_IFCB ['logSize'], binned_data_IFCB ['logNBSS'])
plt.xlabel('log (Biovolume)')
plt.ylabel('log (NBS)')
plt.title('IFCB 3309 with lower threshold defined by instrument specific threshold')


#try to make a for loop to read all dataframes at the same time
path_to_Zooscan = glob(str(path_to_binned_Zooscan) + '/' + str(378) + '*NBSS*') # for 3309, there is an issue in a single point where comparing NBSS between size bins does not work for filtering the shorter side of the spectrum
IFCB_list_lgsample = [str(x) for x in IFCB_list_lgsample]
IFCB_dict = {}
for i in IFCB_list_lgsample:
    path_to_IFCB = glob(str(path_to_binned_IFCB) + '/' + i + '*NBSS*')
    IFCB_dict[i] = pd.read_csv(path_to_IFCB[0], sep='\t', header=0, index_col=[0])
    IFCB_dict[i] = clean_lin_fit(IFCB_dict[i], instrument='IFCB', method= 'MAX')
    binned_data_Zooscan = pd.read_csv(path_to_Zooscan[0], sep='\t', header=0, index_col=[0])
    subset_binned_Zooscan = binned_data_Zooscan.loc[(binned_data_Zooscan['Station_location'] == IFCB_dict[i].loc[0, 'Station_location'])]
    subset_binned_Zooscan = subset_binned_Zooscan.sort_values(by='range_size_bin')
    subset_binned_Zooscan = subset_binned_Zooscan.reset_index(drop=True)
    data_filt_Zooscan = clean_lin_fit(subset_binned_Zooscan, instrument='Zooscan', method= 'instrument_specific')
    IFCB_dict[i] = pd.concat([IFCB_dict[i], data_filt_Zooscan], axis=0)
    IFCB_dict[i] = IFCB_dict[i].reset_index(drop=True)

Station_list = []
for i in IFCB_list_lgsample:
    Station_list.append(IFCB_dict[i].loc[0, 'Station_location'])
    station = (IFCB_dict[i].loc[0, 'Station_location'])
    IFCB_dict[station] = IFCB_dict.pop(i)


fig, ax1 = plt.subplots(1,1)
for st in Station_list:
    ax1.plot(IFCB_dict[st]['logSize'], IFCB_dict[st]['logNBSS'], marker='.', label=st)
ax1.set_ylabel('log (NBS)')
ax1.set_xlabel('log (Biovolume)')
ax1.legend()
ax1.title.set_text('NBSS for IFCB (Max) and Zooscan (isntrument) for six stations sampled on  7/2013')

# read the Zooscan data
path_to_Zooscan = glob(str(path_to_binned_Zooscan) + '/' + str(378) + '*NBSS*')
binned_data_Zooscan = pd.read_csv(path_to_Zooscan[0], sep='\t', header=0, index_col=[0])
#subset zooscan data based on IFCB's station location
subset_binned_Zooscan = binned_data_Zooscan.loc[(binned_data_Zooscan['Station_location'] == binned_data_IFCB.loc[0, 'Station_location'])]
subset_binned_Zooscan = subset_binned_Zooscan.sort_values(by='range_size_bin')
subset_binned_Zooscan= subset_binned_Zooscan.reset_index(drop=True)
data_filt_Zooscan = clean_lin_fit(subset_binned_Zooscan, instrument='Zooscan', method = 'instrument_specific')
plt.plot(data_filt_Zooscan['logSize'], data_filt_Zooscan['logNBSS'])
plt.xlabel('log (Biovolume)')
plt.ylabel('log (NBS)')
plt.title('Zooscan 378, Station 79.5_66.5 with instrument specific threshold')


## MERGE both datasets:
merged_data = pd.concat([data_filt_IFCB, data_filt_Zooscan], axis=0)
merged_data = merged_data.reset_index(drop=True)
plt.plot(merged_data['logSize'], merged_data['logNBSS'])
#1 Plot in one graph the data coming from 2 instruments. for this to work, merge needs to happen BEFORE calculating NBSS


#1) Check the (i) number of beads (ii)  and number of bubbles per time in



#plotting of biovolume and NBSS here. this is to see how the data looks like

#first, all data
plt.plot(data_filt['logSize'], data_filt['logNBSS']) #def not linear :O

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



