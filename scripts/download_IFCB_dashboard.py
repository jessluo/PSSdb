# Process to download raw data from IFCB dashboards



import pandas as pd
import os
import shutil
import numpy as np
import yaml
from tqdm import tqdm
import time
import datetime as dt
from pathlib import Path
from funcs_IFCB_dashboard import *

# set paths:
path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
# open the metadata of the standardized files
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
all_files_data = str(Path(cfg['raw_dir']).expanduser()) + '/IFCB_dashboard_projects_list.xlsx'
IFCB_dashboard_data_path = str(Path(cfg['IFCB_dir']).expanduser())
if not os.path.exists(IFCB_dashboard_data_path):
    os.mkdir(IFCB_dashboard_data_path)


# generate a dataframe with the timeseries list for each dashboard:
timeseries_data =  pd.read_excel(all_files_data)

# subset timeseries based on downloading the ones for testing or all

testing = input (' Do you want to download all projects or just the tests? \n Enter all or tests ')
if testing == 'tests':
    timeseries_data = timeseries_data[timeseries_data['Project_test'] == True].reset_index(drop=True)

subset = input (' Do you want to get data for all the time series or a date range or the test dates? \n Enter all, range or tests ')
if subset == 'range':
    start_date= input(' enter starting date of the data as YYYYMMDD')
    start_date= int(start_date)
    end_date = input(' enter final date of the data as YYYYMMDD')
    end_date = int(end_date)
elif subset == 'all':
    start_date = 20000101
    end_date = 21000101
elif subset == 'tests':
    start_date = 20180825
    end_date = 20180825

dashboard = input (' Do you want to get data from WHOI or CALOOS dashboards? \n Enter WHOI or CALOOS')
if dashboard == 'CALOOS':
    timeseries_data = timeseries_data.loc[timeseries_data['dashboard_url'] == 'https://ifcb.caloos.org/']
if dashboard == 'WHOI':
    timeseries_data = timeseries_data.loc[timeseries_data['dashboard_url'] == 'https://ifcb-data.whoi.edu/']




# download starts here
for n in range (0, len(timeseries_data)):
    start = time.time()
    Project_source = timeseries_data.loc[n, 'dashboard_url']
    Project_ID = timeseries_data.loc[n, 'Project_ID']
    # define path for download
    path_download = IFCB_dashboard_data_path + '/' + Project_ID  ## NOTE: the first part of the directory needs to be changed for the GIT PSSdb project
    pathname = Project_source + Project_ID + '/'

    if os.path.isdir(path_download) and len(os.listdir(path_download) ) != 0: # and  os.path.exists(path_download)
        replace = input('There is already downloaded data in ' +pathname+ ' do you want to replace the files? \n Y/N')
        if replace == 'Y':
            print('Overwriting '+ Project_ID  +' file(s), please wait')
            shutil.rmtree(path_download)
            os.mkdir(path_download)
        elif replace == 'N':
            print('Skipping ' + Project_ID)
            continue
    elif not os.path.exists(path_download):
        os.mkdir(path_download)

    # PROCESSING HERE

    METRIC = 'temperature'  # needs to be defined based on Karina's script
    # use the file list to download data and store it locally
    start_get_df_list = time.time()
    file_info = get_df_list_IFCB(base_url=Project_source, dataset=Project_ID, startdate=start_date, enddate=end_date) # remember this fuction has the option to define an interval of dates
    elapsed_time_df_list = (time.time() - start_get_df_list)
    print('getting df_list took ' + str(elapsed_time_df_list) + ' seconds')
    last_file = file_info.loc[len(file_info.index) - 1, 'pid']
    # the loop uses the file_info dataframe to get the 'pid' (file name) and download features and class scores files
    print ('extracting features files, metadata and top 5 class scores for timeseries ' + timeseries_data.loc[n, 'Project_title'] +' stored in ' + pathname)
    file_numbers = 0
    df_concatenated = pd.DataFrame()
    for i in tqdm(file_info.loc[:, 'pid']): # timeit and progressbar can be used here
        try:
            # generate features files
            pid_id = str(i)
            #print(pid_id)
            features_filename = pathname + pid_id + '_features.csv'
            try:
                start_processing = time.time()
                start_df_metadata = time.time()
                features_df = df_from_url(features_filename)
                # obtain metadata
                met_dict = metadata_dict(dashboard=Project_source, pid_id=pid_id)
                elapsed_time_df_list = (time.time() - start_df_metadata)
                print('getting features and metadata took ' + str(elapsed_time_df_list) + ' seconds')
                features_clean = pd.DataFrame()

                # extract data of interest from features file:
                features_clean['area'] = features_df['Area']
                features_clean['biovolume'] = features_df['Biovolume']
                features_clean['equiv_diameter'] = features_df['EquivDiameter']
                features_clean['major_axis_length'] = features_df['MajorAxisLength']
                features_clean['minor_axis_length'] = features_df['MinorAxisLength']
                features_clean['solidity'] = features_df['Solidity']
                features_clean['summed_area'] = features_df['summedArea']
                features_clean['summed_biovolume'] = features_df['summedBiovolume']
                features_clean['summed_major_axis_length'] = features_df['summedMajorAxisLength']
                features_clean['summed_minor_axis_length'] = features_df['summedMinorAxisLength']

                #add metadata
                features_clean['project_ID'] = Project_ID
                features_clean['datetime'] = met_dict['datetime']
                #features_df['date'] = str(pd.to_datetime(met_dict['datetime']).date().strftime('%Y%m%d'))
                #features_df['time'] = str(pd.to_datetime(met_dict['datetime']).time().strftime('%H%M%S'))
                features_clean['pixel_to_micron'] = met_dict['scale']
                features_clean['latitude'] = met_dict['latitude']
                features_clean['longitude'] = met_dict['longitude']
                features_clean['depth'] = met_dict['depth']
                features_clean['vol_analyzed'] = met_dict['ml_analyzed']
                features_clean['sample_type'] = met_dict['sample_type']
                features_clean['sample_cruise'] = met_dict['cruise']
                features_clean['number_of_rois'] = met_dict['number_of_rois']
                features_clean['concentration'] = met_dict['concentration']

                # now generate dataframe for the class scores
                class_filename = pathname + str(i) + '_class_scores.csv'

                start_df_scores = time.time()
                class_df = df_from_url(class_filename)

                elapsed_time_scores = (time.time() - start_df_scores)
                print('getting class scores took ' + str(elapsed_time_scores) + ' seconds')

                features_clean['roi_id'] = class_df['pid']
                class_df = class_df.set_index('pid')
                class_df = class_df.astype(float)
                # create lists for top five classes and their scores
                class_1 = []
                class_1_score = []
                class_2 = []
                class_2_score = []
                class_3 = []
                class_3_score = []
                class_4 = []
                class_4_score = []
                class_5 = []
                class_5_score = []
                try:
                    for roi in class_df.index:
                        #print(roi)
                        index_top5 = np.argsort(-class_df.loc[roi].values)[:5]
                        #print(index_top5)
                        class_1.append(class_df.columns[index_top5[0]])
                        class_1_score.append(class_df.loc[roi, class_df.columns[index_top5[0]]])

                        class_2.append(class_df.columns[index_top5[1]])
                        class_2_score.append(class_df.loc[roi, class_df.columns[index_top5[1]]])

                        class_3.append(class_df.columns[index_top5[2]])
                        class_3_score.append(class_df.loc[roi, class_df.columns[index_top5[2]]])

                        class_4.append(class_df.columns[index_top5[3]])
                        class_4_score.append(class_df.loc[roi, class_df.columns[index_top5[3]]])

                        class_5.append(class_df.columns[index_top5[4]])
                        class_5_score.append(class_df.loc[roi, class_df.columns[index_top5[4]]])
                except:
                    print(i + ' does not have a class score')
                features_clean['class_1'] = class_1
                features_clean['class_1_score'] = class_1_score

                features_clean['class_2'] = class_2
                features_clean['class_2_score'] = class_2_score

                features_clean['class_3'] = class_3
                features_clean['class_3_score'] = class_3_score

                features_clean['class_4'] = class_4
                features_clean['class_4_score'] = class_4_score

                features_clean['class_5'] = class_5
                features_clean['class_5_score'] = class_5_score

                #replace with nan the class and class scores below 0.0001
                cols = ["class_1", "class_2", "class_3", "class_4", "class_5"]
                cols_scores = ["class_1_score", "class_2_score", "class_3_score", "class_4_score", "class_5_score"]
                for no, c in enumerate(cols_scores):
                    features_clean[c] = features_clean[c].mask(features_clean[c] < 0.0001, np.nan)
                    features_clean[cols[no]] = features_clean[cols[no]].mask(features_clean[c] < 0.0001, np.nan)

                features_clean['datetime'] = pd.to_datetime(features_clean['datetime'])
                year = features_clean.loc[1, 'datetime'].strftime('%Y')
                dataset_id = features_clean.loc[1, 'roi_id'].split('_')
                dataset_id = dataset_id[1]
                if Project_source == 'https://ifcb.caloos.org/':
                    dashboard_id = 'CALOOS'
                elif Project_source == 'https://ifcb-data.whoi.edu/':
                    dashboard_id = 'WHOI'

                df_concatenated = pd.concat([df_concatenated, features_clean], ignore_index=True)

                #if len(df_concatenated.index) <= 100000:
                    #df_concatenated = pd.concat([df_concatenated, features_df], ignore_index=True)
                    #pass
                if pid_id == last_file :
                    file_numbers = file_numbers + 1
                    print ('saving file # ' + str(file_numbers) + ' of '+ Project_ID)
                    df_concatenated.to_csv(path_download + '/' + dashboard_id +'_'+ dataset_id +'_'+ df_concatenated.loc[1, 'project_ID'] +'_'+ year + '_features_' + str(file_numbers) +'.tsv', sep='\t')
                    df_concatenated = pd.DataFrame()
                    #pass
                elif len(df_concatenated.index) > 500000:
                    file_numbers = file_numbers + 1
                    print ('saving file # ' + str(file_numbers) + ' of '+ Project_ID)
                    df_concatenated.to_csv(path_download + '/' + dashboard_id +'_'+ dataset_id +'_'+ df_concatenated.loc[1, 'project_ID'] +'_'+ year + '_features_' + str(file_numbers) +'.tsv', sep='\t')
                    df_concatenated = pd.DataFrame()

                elapsed_time_processing = (time.time() - start_processing)
                print('all processing took ' + str(elapsed_time_processing) + ' seconds')

                # class_df.to_csv(path_download + '/' + str(i) + '_class_scores.csv', sep='\t') 11/16/2022 decided not to keep class scores
                # print(str(i) + ' download done ')
            except:
                print('there is no features or class_scores files for ' + str(file_info.loc[i, 'pid']))
        except:
            pass
        elapsed_time_processing = (time.time() - start_processing)
        print('all processing took ' + str(elapsed_time_processing) + ' seconds')

    elapsed_time_fl = (time.time() - start)
    print(timeseries_data.loc[n, 'Project_title'] + ' download took ' + str((elapsed_time_fl/3600)) + ' hours')





