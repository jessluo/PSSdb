# gridding standardized files generated in 2.
import warnings
warnings.filterwarnings('ignore')
## load packages
import numpy as np
import pandas as pd
import os
try:
    from funcs_read import * # functions to list standardized and binned files
    from funcs_gridding import *
except:
    from scripts.funcs_read import * # functions to list standardized and binned files
    from scripts.funcs_gridding import *
from pathlib import Path
from glob import glob
import yaml  # requires installation of PyYAML package
import os
import shutil
from tqdm import tqdm

path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

## set up parameters for script:
cats = cfg['keep_cats']
if cats == 'Y':
    keep_cat= cfg['cats_to_keep']
elif cats == 'N':
    keep_cat = 'none'


st_increment = float(cfg['st_increment'])
date_group = cfg['date_group']
depth_binning = cfg['depth_binning']
if depth_binning == True:
    depth_bins = cfg['depth_bins']

## processing starts here
for instrument in ['Scanner', 'UVP', 'IFCB']:
    files_subset_df_list = pd.DataFrame()
    files_data = proj_id_list_func(instrument, data_status='standardized')  # generate path and project ID's
    test_1 = cfg['N_test']
    if test_1 == 0:
        files_data = files_data
    else:
        files_data = files_data[0:test_1]
    gridpath = Path(cfg['raw_dir']).expanduser() / 'gridded_data'
    if not os.path.exists(gridpath):
        os.mkdir(gridpath)

    dirpath = Path(cfg['raw_dir']).expanduser() / 'gridded_data' / instrument
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
    elif os.path.isdir(dirpath) and len(os.listdir(dirpath)) != 0:  # and  os.path.exists(path_download)
        replace = input('There is already gridded data in ' + str(dirpath) + ' do you want to replace the files? \n Y/N')
    if (len(os.listdir(dirpath)) == 0) or (replace == 'Y'):
        shutil.rmtree(dirpath)
        os.mkdir(dirpath)
        n_del_files = 0
        #spatial gridding starts here:
        for i in tqdm(files_data):
            filename= i.split('/')[-1]
            print('gridding  ' + i)
            df = pd.read_csv(i, header=0)
            df = df.dropna(subset=['Latitude']).reset_index(drop=True)
            if len(df) == 0:
                print('no data left after removing ROIS with no Lat-lon information in' + filename)
                n_del_files += 1
                continue
            df = df[df['Area'] != 0].reset_index(drop=True)
            #if instrument == 'UVP':
                #df = remove_UVP_noVal(df)
                #if len(df) == 0:
                    #print ('no data left after assessing validation status for ' + filename)
                    #n_del_files += 1
                    #continue
            if (instrument == 'Scanner') or (instrument == 'UVP'):
                df = depth_parsing_func(df, instrument)
                if len(df) == 0:
                    print('no data left after restricting depths to less than 200 meters in ' + filename)
                    n_del_files += 1
                    continue
            elif instrument == 'IFCB':
                df = depth_parsing_func(df, instrument)
                df = df.replace(np.nan, '')
                sampling_type_to_remove = ['test', 'exp', 'junk', 'culture']
                df = df[~df.Sampling_type.str.contains('|'.join(sampling_type_to_remove))]
                df = df.replace('', np.nan)
                df = df[df['Sample'].str.contains('D20190603T001443_IFCB115') == False].reset_index()
                if len(df) == 0:
                    print('no data left after restricting depths to less than 200 meters in ' + filename)
                    n_del_files += 1
                    continue
            df = biovol_func(df, instrument, keep_cat='none')
            df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(st_increment, df['Latitude'], df['Longitude'])
            df['Instrument']= instrument # necessary to set the same instrument name for all scanner types
            filename = filename.replace("standardized", "gridded")
            df.to_csv(str(dirpath) + '/' + instrument + '_' + filename, index=False)
        grid_list = group_gridded_files_func(instrument) # saving files in 1x1  lat/lon cells to facilitate computation

        # Temporal binning function starts here: since we consider sample size to assign time bins, data needs to be aggregated by each cell
        file_list = proj_id_list_func(instrument, data_status='gridded', big_grid=True)
        for cell in tqdm(grid_list):
            print('assinging temporal bin to large cell ' + cell)
            file_subset = [file for file in file_list if cell in file]
            df_gridded = pd.concat(map((lambda path: (pd.read_csv(path))), file_subset)).reset_index(drop=True)
            df_gridded_temp_binned = pd.DataFrame()
            #df_gridded_temp_binned = df_gridded.groupby(['Station_location']).apply(lambda x: date_binning_func(x, group_by=date_group)).reset_index
            for st in list(set(df_gridded['Station_location'])):
                # print(st)
                df_st_subset = df_gridded[df_gridded['Station_location'] == st].reset_index(drop=True)
                df_st_subset = date_binning_func(df_st_subset, group_by=date_group)
                if len(df)==0:
                    print('no data left after removing data without proper time stamp in ' + filename) # necessary because some projects don't have time info: '/Users/mc4214/GIT/PSSdb/raw/raw_standardized/ecotaxa/Zooscan/standardized_project_5785_20221103_1928.csv'
                    n_del_files += 1
                    continue
                df_st_subset['date_bin'] = df_st_subset['date_bin'].astype(str)
                df_gridded_temp_binned = pd.concat([df_gridded_temp_binned, df_st_subset])
                # separate data into months of year and save to avoid the creation of large files
                for m in list(set(df_gridded_temp_binned['date_bin'])):
                    df_st_t_subset = df_gridded_temp_binned[df_gridded_temp_binned['date_bin'] == m].reset_index(drop=True)
                    df_project_list_bins = [{'gridded_file_ID':cell +'_temp_binned_'+m, 'Project_ID_list': df_st_t_subset.Project_ID.unique().tolist(), 'Sample_list': df_st_t_subset.Sample.unique().tolist()}]
                    df_project_list_bins = pd.DataFrame(df_project_list_bins)
                    files_subset_df_list = pd.concat([files_subset_df_list, df_project_list_bins])
                    df_st_t_subset.to_csv(str(dirpath) +'/'+ instrument + '_gridded_' + cell +'_temp_binned_'+m+'.csv', index=False)

        files_subset_df_list .to_csv(str(dirpath) + '/' + instrument + '_' + 'Project_Sample_list_gridded.csv', index=False)

        for file in glob(str(Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir']) + '*/**/grid_N*' + instrument + '*.csv', recursive=True):
            os.remove(file)


        print(str(n_del_files) + ' files were not gridded due to missing Lat/Lon, deep water sample, time information or validation status < 95% for ' +instrument)
        if depth_binning == False:
            metadata_bins = pd.DataFrame(
                    {'Variables': ['Biovolume_area', 'Biovolume_ellipsoid', 'Biovolume_orig', 'date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'lat_grid', 'lon_grid', 'grid_id',  'date_grouping', 'light_cond'],
                     # 'light_cond',
                     'Variable_types': ['float64','float64','float64', 'int64', 'object', 'float64', 'float64','object','object','object','object','object'],  # 'str',
                     'Units/Values/Timezone': ['cubic_micrometer','cubic_micrometer','cubic_micrometer', date_group, 'lat_lon', 'degree', 'degree','lat-lat','lon-lon','lat-lat_lon-lon', 'month-year', ''],  # '',
                     'Description': ['Biovolume calculated as a spherical projection of Area in cubic micrometers',
                                     'Biovolume calculated as an ellipsoid following Dubois (2022)  in cubic micrometers',
                                     'Biovolume from IFCB dashboard, from ifcb-analysis  in cubic micrometers',
                                     'binned date information',
                                     'string that serves as an identifier of a single cell of a  0.5x0.5 degree spatial grid',
                                     'latitude of the center point of the 0.5x0.5 degree cell',
                                     'longitude of the center point of the 0.5x0.5 degree cell',
                                     'latitude label used to group data when binning by date',
                                     'longitude label used to group data when binning by date',
                                     'grid label used to group data when binning by date',
                                     'label used to read data in chunks during step4',
                                     'label used to group data by day or night']})
        elif depth_binning == True:
            df['midDepthBin'] = depth_binning_func(df['Depth_max'], depth_bins=depth_bins)
            metadata_bins = pd.DataFrame(
                    {'Variables': ['Biovolume_area', 'Biovolume_ellipsoid', 'Biovolume_orig', 'date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'lat_grid', 'lon_grid', 'grid_id',  'date_grouping', 'light_cond', 'midDepthBin'],
                     # 'light_cond',
                     'Variable_types': ['float64','float64','float64', 'int64', 'object', 'float64', 'float64','object','object','object','object','object', 'float64'],  # 'str',
                     'Units/Values/Timezone': ['cubic_micrometer','cubic_micrometer','cubic_micrometer', date_group, 'lat_lon', 'degree', 'degree','lat-lat','lon-lon','lat-lat_lon-lon', 'month-year', '', 'meters'],  # '',
                     'Description': ['Biovolume calculated as a spherical projection of Area in cubic micrometers',
                                     'Biovolume calculated as an ellipsoid following Dubois (2022)  in cubic micrometers',
                                     'Biovolume from IFCB dashboard from ifcb-analysis  in cubic micrometers',
                                     'binned date information',
                                     'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                     'latitude of the center point of the 0.5x0.5 degree cell',
                                     'longitude of the center point of the 0.5x0.5 degree cell',
                                     'latitude label used to group data when binning by date',
                                     'longitude label used to group data when binning by date'
                                     'grid label used to group data when binning by date',
                                     'label used to read data in chunks during step4',
                                     'label used to group data by day or night',
                                     'mid depth in a depth bin']})
        metadata_bins.to_csv(str(dirpath) + '/' + instrument + '_' + 'metadata_gridded.csv', index=False)
    elif replace == 'N':
        print('previously gridded files will be kept')


