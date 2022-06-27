#function that calls on
def PSS_func()
    import funcs_size_spectra
    path_to_data, ID = proj_id_list('IFCB', testing=True)
    data_clean_all=pd.DataFrame()
        for i in ID:
            data_clean = read(path_to_data, i)
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
