


## 7)  MERGER FUNCTION here: one way to make this (especially if it will have to be done globally) is to use as imput
# the date only
#step A: generate a dictionary with dates for the projects
#
def corr_proj_func(instrument, variable):
    """
    Objective: create a dictionary where each project has an associated list of dates
    :param instrument:
    :return:
    """
    #from funcs_size_spectra import proj_id_list_func
    path_to_binned, ID_list = proj_id_list_func(instrument, standardized='binned', testing=False)
    if instrument == 'IFCB':# this removal is only for testing,
        #since the standardizer should be able to deal with projects with empty rows
        IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]
        for element in IFCB_empty_rows:
            if element in ID_list:
                ID_list.remove(element)
    date_dict ={}
    # make a list for each binned file that has the dates
    for n, ID in enumerate(ID_list):
        path_to_file = glob(str(path_to_binned)+ '/'+ str(ID)+ '*binned*')
        binned_data = pd.read_csv(path_to_file[0], sep='\t', header=0, index_col=[0])
        date_dict[ID] = binned_data[variable].unique().astype('str').tolist()

    return date_dict




def merge_data_func(binning, variable,  instrument1 = 'IFCB', instrument2 = 'Zooscan'):
    """
    Objective: find STANDARDIZED, BINNED tsv files from different instruments and stitch them together.
    :param date: date of interest of the dataset (how close do we want it to be?)
    :return: a dataframe that contains all the size classes across instruments
    """
    date_dict_IFCB = corr_proj_func(instrument1, variable)
    date_dict_Zooscan = corr_proj_func(instrument2, variable)
    list_proj_IFCB = []
    list_proj_Zooscan = []
    for key in date_dict_IFCB:
        if binning in date_dict_IFCB[key]:
            list_proj_IFCB.append(key)
    for key in date_dict_Zooscan:
        if binning in date_dict_Zooscan[key]:
            list_proj_Zooscan.append(key)

    return list_proj_IFCB, list_proj_Zooscan

#another thing that can happen is that i subset the dictionary based on a list of pre defined conditions

new_dict_Zooscan = dict((k, station_dict_Zooscan[k]) for k in proj_Zooscan_dates)# this works, I cut down the stations
# that were sampled in a particular date
new_dict_IFCB = dict((k, station_dict_IFCB[k]) for k in proj_IFCB_dates)

for key in new_dict_Zooscan:
    if new_dict_IFCB[3147][0] in new_dict_Zooscan[key]:
        print(key)
