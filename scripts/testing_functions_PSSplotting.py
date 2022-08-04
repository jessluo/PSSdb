## Plotting scripts
import pandas as pd

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

#plot some example IFCB's NBSS with/without bubbles, with/without beads and with/without artefacts
binning_NBS_func('IFCB', 'artefacts')
path_to_binned_IFCB, ID_list_IFCB =  proj_id_list_func('IFCB', standardized='binned', testing=False)
ID_list_IFCB = ['2248', '3147', '3289', '3295', '3296', '3302', '3307', '3309']
IFCB_dict ={}
for i in ID_list_IFCB:
    path_to_IFCB = glob(str(path_to_binned_IFCB) + '/' + i + '*NBSS*')
    IFCB_dict[i] = pd.read_csv(path_to_IFCB[0], sep='\t', header=0, index_col=[0])
    #IFCB_dict[i] = clean_lin_fit(IFCB_dict[i], instrument='IFCB', method= 'MAX')

fig, ax1 = plt.subplots(1,1)
for i in ID_list_IFCB:
    ax1.plot(IFCB_dict[i]['logSize'], IFCB_dict[i]['logNBSS'], marker='.', label=i)
ax1.set_ylabel('log (NBS)')
ax1.set_xlabel('log (Biovolume)')
ax1.legend()
ax1.title.set_text('NBSS for IFCB without bubbles and without artefacts')

#given the similarity of data with and without artefacts, We'll dig deeper into it
path_to_binned_IFCB, ID_list_IFCB =  proj_id_list_func('IFCB', standardized='Yes', testing=False)
IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
for element in IFCB_empty_rows:
    if element in ID_list_IFCB:
        ID_list_IFCB.remove(element)
ID_list_IFCB = [str(x) for x in ID_list_IFCB]
IFCB_dict = {}
for i in ID_list_IFCB:
    path_to_IFCB = glob(str(path_to_binned_IFCB) + '/standardized_export_' + i + '.tsv')
    IFCB_dict[i] = pd.read_csv(path_to_IFCB[0], sep='\t', header=0, index_col=[0])


df = IFCB_dict['3343']
df = proj_bin_func(df, 'IFCB','artefact')
df_wArt = IFCB_dict['3343']
df_wArt = proj_bin_func(df_wArt, 'IFCB',  'none')
IFCB_proj ={}
IFCB_proj['NoArt'] = df
IFCB_proj['Art'] = df_wArt
fig, ax1 = plt.subplots(1,1)
for i in ['NoArt', 'Art']:
    ax1.plot(IFCB_proj[i]['logSize'], IFCB_proj[i]['logNBSS'], marker='.', label=i)
ax1.set_ylabel('log (NBS)')
ax1.set_xlabel('log (Biovolume)')
ax1.legend()
ax1.title.set_text('NBSS for IFCB 3343 with and without artefacts')


#try to make a for loop to read all dataframes at the same time
def stitching_func(IFCB_proj_type= 'large_sample', method1= 'MAX', method2 = 'MAX'):
    path_to_binned_IFCB, ID_list_IFCB = proj_id_list_func('IFCB', standardized='binned', testing=False)
    path_to_binned_Zooscan, ID_list_Zooscan = proj_id_list_func('Zooscan', standardized='binned', testing=False)

    if IFCB_proj_type == 'large_sample':
        IFCB_list= [3302, 3309, 3321, 3324, 3325, 3334]
    elif IFCB_proj_type == 'small_sample':
        IFCB_list= [3147, 2248, 3289, 3296]
    IFCB_list= [str(x) for x in IFCB_list]
    path_to_Zooscan = glob(str(path_to_binned_Zooscan) + '/' + str(378) + '*NBSS*') # for 3309, there is an issue in a single point where comparing NBSS between size bins does not work for filtering the shorter side of the spectrum
    IFCB_list= [str(x) for x in IFCB_list]
    IFCB_dict = {}
    for i in IFCB_list:
        path_to_IFCB = glob(str(path_to_binned_IFCB) + '/' + i + '*NBSS*')
        IFCB_dict[i] = pd.read_csv(path_to_IFCB[0], sep='\t', header=0, index_col=[0])
        IFCB_dict[i] = clean_lin_fit(IFCB_dict[i], instrument='IFCB', method= method1)
        binned_data_Zooscan = pd.read_csv(path_to_Zooscan[0], sep='\t', header=0, index_col=[0])
        subset_binned_Zooscan = binned_data_Zooscan.loc[(binned_data_Zooscan['Station_location'] == IFCB_dict[i].loc[0, 'Station_location'])]
        subset_binned_Zooscan = subset_binned_Zooscan.sort_values(by='range_size_bin')
        subset_binned_Zooscan = subset_binned_Zooscan.reset_index(drop=True)
        data_filt_Zooscan = clean_lin_fit(subset_binned_Zooscan, instrument='Zooscan', method= method2)
        IFCB_dict[i] = pd.concat([IFCB_dict[i], data_filt_Zooscan], axis=0)
        IFCB_dict[i] = IFCB_dict[i].reset_index(drop=True)
    Station_list = []
    for i in IFCB_list:
        Station_list.append(IFCB_dict[i].loc[0, 'Station_location'])
        station = (IFCB_dict[i].loc[0, 'Station_location'])
        IFCB_dict[station] = IFCB_dict.pop(i)

    fig, ax1 = plt.subplots(1,1)
    for st in Station_list:
        ax1.plot(IFCB_dict[st]['logSize'], IFCB_dict[st]['logNBSS'], marker='.', label=st)
    ax1.set_ylabel('log (NBS)')
    ax1.set_xlabel('log (Biovolume)')
    ax1.legend()
    ax1.title.set_text('NBSS for IFCB ' + method1 + ' and Zooscan ' + method2  + ' for ' + IFCB_proj_type)

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



### EXPLORING IFCB DATA ####
# plot number of bubbles per project

#set path, project list, remove projects with empty spaces, and convert list to string
path_to_binned_IFCB, ID_list_IFCB =  proj_id_list_func('IFCB', standardized='Yes', testing=False)
IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
for element in IFCB_empty_rows:
    if element in ID_list_IFCB:
        ID_list_IFCB.remove(element)
ID_list_IFCB = [str(x) for x in ID_list_IFCB]

IFCB_dict = {}
bubble_count = []
for i in ID_list_IFCB:
    path_to_IFCB = glob(str(path_to_binned_IFCB) + '/standardized_export_' + i + '.tsv')
    IFCB_dict[i] = pd.read_csv(path_to_IFCB[0], sep='\t', header=0, index_col=[0])
    bubbly = IFCB_dict[i][IFCB_dict[i].Category == 'bubble'].Category.value_counts()
    if len(bubbly) == 0:
        bubble_count.append(0)
    else:
        bubble_count.append(bubbly['bubble'])

plt.bar(ID_list_IFCB, bubble_count)
plt.ylabel('# of bubbles')
plt.xlabel('IFCB project ID')


# number of beads per time, couldn't do this for all projects, so need to subset every

df_test['freq'] = df_test.groupby(['Category'])['date_bin'].transform(len)
plt.hist(x=df_test.loc[df_test['Category'] == 'bead'].date_bin, bins = 78) #bins = 78

# reduce dataframe and get counts of beads without changing time axis
col_list = ['date_bin', 'Category']
bead_count = []
for proj in ID_list_IFCB:
    IFCB_dict[proj]['date_bin'] = date_binning_func(IFCB_dict[proj]['Sampling_date'], IFCB_dict[proj]['Sampling_time'], group_by = 'None')
    IFCB_dict[proj] = IFCB_dict[proj][col_list]
    IFCB_dict[proj]['freq'] = IFCB_dict[proj].groupby(['Category'])['date_bin'].transform(len)
    beady = IFCB_dict[proj][IFCB_dict[proj].Category == 'bead'].Category.value_counts()
    if len(beady) == 0:
        bead_count.append(0)
    else:
        bead_count.append(beady['bead'])

# remove from IFCB projects the ones that dont have beads:
for n, p in enumerate(ID_list_IFCB):
    if bead_count[n] == 0:
        ID_list_IFCB.remove(p)

# try to do a subplot for each project NONE OF THIS WORKED
fig, ax_arr = plt.subplots(5, 7)
for n, i in enumerate(ID_list_IFCB):
    to_graph = IFCB_dict[i]
    ax_arr[n].hist(x=to_graph.loc[to_graph['Category'] == 'bead'].date_bin, bins=5)
    ax_arr[n].set_title(n)
    ax_arr[n].set_ylim([0, max(to_graph.freq)])

fig, ax_arr = plt.subplots(5, 7)
for n, i in enumerate(ID_list_IFCB):
    to_graph = IFCB_dict[i]
    for i in range (5):
        for j in range (7):
            to_graph.loc[to_graph['Category'] == 'bead'].date_bin.hist(ax=ax_arr[i,j])


to_graph = IFCB_dict['3296']
to_graph.loc[to_graph['Category'] == 'bead'].date_bin.hist()
plt.ylabel('# of beads')
plt.xlabel('time')
plt.title(' IFCB #3296')

#plot all the time
beads_all_time = pd.DataFrame()
for i in ID_list_IFCB:
    to_graph = IFCB_dict[i]
    beads_all_time = pd.concat([beads_all_time, to_graph], axis=0)

beads_all_time.loc[ beads_all_time['Category'] == 'bead'].date_bin.hist( bins= 200)
plt.ylabel('# of beads')
plt.xlabel('time')
plt.title('all Tara Oceans projects')





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



