# Functions necessary to download raw files from CALOOS and WHOI dashboards

## create list of dates available for the project of interest
def get_df_list_IFCB ( base_url, dataset, startdate=20000101, enddate=21000101):
    """
    Objective: generate list of files for a dataset in the IFCB dashboard
    :param startdate: set start date to  select data from the project, in format YYYYMMDD
    :param enddate: set end  date to  select data from the project, in format YYYYMMDD
    :param base_url: WHOI or CALOOS url
    :param dataset: data set contained within the WHOI or CALOOS dashboard
    :return: data frame that contains the dataset names under the 'pid' column
    """
    import requests
    import json
    import pandas as pd
    link = base_url + 'api/list_bins?dataset=' + dataset # this path is for the data LIST only
    #metadata_link = "https://ifcb.caloos.org/timeline?dataset="
    r = requests.get(link)
    r_content = r.content
    r_content = json.loads(r_content.decode('utf-8'))
    all_file_info = pd.DataFrame.from_dict(r_content['data'], orient = 'columns')
    all_file_info['date_info'] = [int(sample[1:9]) for sample in all_file_info['pid']]
    # constrain date ranges for file download.
    #start_stamp=int(startdate.strftime('%Y%m%d'))
    #end_stamp=int(y['enddate'].strftime('%Y%m%d'))
    file_info = all_file_info.loc[all_file_info['date_info']>= startdate].reset_index(drop =True) # | all_file_info['date_info']<= end_stamp]
    file_info = file_info.loc[file_info['date_info']<= enddate].reset_index(drop =True)
    return file_info

# We adopted some code contributed by Joe Futrelle (WHOI): https://github.com/joefutrelle/pyifcb

# make function to get datasets from urls
def df_from_url(url):
    """
    Objective: use requests to build dataframes with ROI information from contents obtained from IFCB dashboards
    :param url: url for the dataset
    :return: a dataframe with the desired information (see outputs of IFCB dashboards)
    """
    import requests
    import pandas as pd
    data = requests.get(url)
    data = data.content
    data = data.decode('utf-8')
    data = data.split('\n')
    data_dict = {}
    for n, r in enumerate(data):
        roi_info = r.split(',')
        data_dict[n] = roi_info
    data_dict.popitem()
    df = pd.DataFrame.from_dict(data_dict, orient='index')
    df.columns = df.iloc[0]
    df.drop(df.index[0], inplace = True)
    df.columns = df.columns.str.replace('\r', '')
    if df.columns.__contains__(None):
        print('features file has not been generated')
    else:
        return df


def metadata_dict(dashboard, pid_id):
    """
    :param dashboard: url to either the WHOI or CALOOS dashboard
    :param pid: the identifier for the data in each time bin displayed on the timeline in the dashboard
    :return: dictionary with the metadata
    """
    import requests
    import re
    url = f'{dashboard}api/bin/{pid_id}?include_coordinates=False'
    r = requests.get(url)
    assert r.ok
    record = r.json()
    return {
        'scale': record['scale'],
        'datetime': record['timestamp_iso'],
        'previous_bin': record['previous_bin_id'],
        'next_bin': record['next_bin_id'],
        'latitude': record['lat'],
        'longitude': record['lng'],
        'depth': record['depth'],
        'instrument': record['instrument'],
        'number_of_rois': record['num_images'],
        'ml_analyzed': float(re.sub(' .*', '', record['ml_analyzed'])),
        'concentration': record['concentration'],
        'sample_type': record['sample_type'],
        'cruise': record['cruise'],
        'cast': record['cast'],
        'niskin': record['niskin'],
    }
