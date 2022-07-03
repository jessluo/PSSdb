#script to test standardization of all IFCB projects and coupling together with the PSS functions

#import functions to 1. get project ID lists and
from funcs_size_spectra import proj_id_list_func as id_list
#get the standardizer function
from funcs_export_standardizer import standardization_func as standardize

config_path = '/Users/mc4214/GIT/PSSdb/scripts/Ecotaxa_API_pw.yaml'

standardizer_path = '/Users/mc4214/GIT/PSSdb/raw/project_Zooscan_standardizer.xlsx'

path_to_data, ID_list = id_list('Zooscan', testing = False)
# this removal is only for testing, since the standardizer should be able to deal with projects with empty rows
IFCB_empty_rows = [3290, 3294, 3297, 3298, 3299, 3313, 3314, 3315, 3318, 3326, 3335, 3337]

for element in IFCB_empty_rows:
    if element in ID_list:
        ID_list.remove(element)

# this is working well, but when the function runs into projects with empty lat & lon it gets stuck.
# that's why they need to be removed

for i in ID_list:
    standardize(config_path, standardizer_path, i)