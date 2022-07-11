

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
          lon='midLonBin',
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

#open files as dataframes and bin them by location and depth:
def binning_all(instrument, scale):
    """
    Objective: read into the standardized tsv files and bin the data by size (biovolume), station and depth
    :param instrument: the device used for image adcquisition. important since the path will change
    :return: tsv files with the same data but binned (would we like to delete the standardized data?)
    """
    path_to_data, ID = proj_id_list(instrument, standardized=True testing=True)#generate path and project ID's

    for i in ID:
        df = read_func(path_to_data, i, cat='Category')# get a dataframe for each project
        # there has to be a step here where biovolume is calculated
        df['Biovolume_um3'] = biovol_um3_func()##INCOMPLETE, CONTINUE HERE
        df['sizeClasses'], df['range_size_bin'] = size_binning_func(df['biovolume'])
        df['midDepthBin'] = depth_binning_func(data_clean_all['Depth_max'])
        df['Station_ID'], df['midLatBin'], df['midLonBin'] = \
            station_binning_func(df['Latitude'], df['Longitude'])

