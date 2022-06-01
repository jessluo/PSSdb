

# follow these to test functionality (this is to figure out why we are getting such horrible slopes)
path_to_data, ID = proj_id_list('IFCB', testing=True)
data_clean_all=pd.DataFrame()
for i in ID:
    data_clean = read_clean(path_to_data, i)
    data_clean_all = pd.concat([data_clean_all, data_clean], axis=0)

data_clean_all = biovol(data_clean_all, 'IFCB')

data_clean_all['biovol_um3']=np.random.poisson(data_clean_all['biovol_um3'].mode(),
                                               data_clean_all['biovol_um3'].count())

data_clean_all = binning(data_clean_all)
stats_biovol_SC = NB_SS(data_clean_all)

results_SS_all=pd.DataFrame()
for i in ID:
    results_SS = NB_SS_regress(stats_biovol_SC, i)
    results_SS_all = pd.concat([results_SS_all , results_SS], axis=0)
