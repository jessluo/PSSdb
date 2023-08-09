


import os
from pathlib import Path
from glob import glob
import yaml
import pandas as pd
import calendar
import numpy as np
import datetime as dt
from tqdm import tqdm
import pytz # to define the timezone of the time info
from datetime import datetime
from astral.sun import sun #run pip install astral
# to assign day/night to an ROI, first create a LocationInfo element. See instructions: https://astral.readthedocs.io/en/latest/
from astral import LocationInfo
# an alternative to astral
from suncalc import get_position, get_times # pip install suncalc
#from dtki import d_n_calc # R. Kiko's module, see PSSdb documentation to install
try:
    from funcs_read import *
except:
    from scripts.funcs_read import *



# NOTE:  biovol_func_old deprecated as of 1/11/2023. This function allowed the user to choose which type of area estimation to use to
#get biovolume
def biovol_func_old(df, instrument, area_type= 'object_area', remove_cat='none'):
    """
    Objective: calculate biovolume (in cubic micrometers) of each object, only for UVP and Zooscan projects, following
    the volume of  a sphere. Also, determine which area will be used to calculate the biovolume.
    This function also removes Biovolumes= 0 in the IFCB files
    :param df: a STANDARDIZED dataframe that contains object's Area (should be in micrometers^2)
     and object's Minor_axis (should be in micrometers)
    :param area_type: name of the column that contains the area to be used to calculate biovolume
    :param geom_shape: geometric shape to be used to calculate biovolume, either 'sphere' or 'ellipse'
    :return: Biovolume in cubic micrometers
    """
    import math as m
    from pathlib import Path
    from glob import glob
    import yaml
    import pandas as pd
    # find the Area column in the metadata file and use its position to select the area
    ID = df.loc[0, 'Project_ID']
    path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    path_to_metadata = glob(str(Path(cfg['git_dir']).expanduser() / cfg['standardized_subdir'])+ '/**/*'+ str(ID)+ '*metadata*')
    metadata = pd.read_csv(path_to_metadata[0], sep='\t')
    if instrument == 'IFCB':# for IFCB projects, it needs to remove the column that is repeated. Also, remove bubbles
        if 'summed' in area_type:
            df = df.drop(df.columns[[19]], axis=1)
            df = df[df.Biovolume != 0]
            df = df[df.Category != 'bubble']
            df = df.reset_index(drop=True)
        else:
            df = df.drop(df.columns[[20]], axis=1)
            df = df[df.Biovolume != 0]
            df = df[df.Category != 'bubble']
            df = df.reset_index(drop=True)
            if 'summed' in area_type:
                ind = metadata.loc[metadata['Description'].str.contains('summed')].index[0]
            else:
                ind = metadata.loc[(metadata['Description'].str.contains('object_area') == True) &
                                   (metadata['Description'].str.contains('summed') == False)].index[0]
    if instrument == 'IFCB':
        if remove_cat == 'beads':
            df = df[df['Category'].str.contains("bead") == False]
            df = df.reset_index(drop=True)
        if remove_cat == 'artefact':
            df = df[df['Category'].str.contains('artefact') == False]
            df = df.reset_index(drop=True)
        elif remove_cat == 'none':
            df = df
            if 'summed' in area_type:
                ind = metadata.loc[metadata['Description'].str.contains('summed')].index[0]
            else:
                ind = metadata.loc[(metadata['Description'].str.contains('object_area') == True) &
                                   (metadata['Description'].str.contains('summed') == False)].index[0]

    elif instrument == 'Zooscan':
        df = df[df['Category'].str.contains('artefact') == False]
        df = df[df['Category'].str.contains('detritus') == False]
        df = df.reset_index(drop=True)
        ind = metadata.loc[(metadata['Description'].str.contains('object_area')==True) &
                           (metadata['Description'].str.contains('exc')== False)].index[0]
    for i in range(0, len(df)):
        r = m.sqrt((df.iloc[i, [ind]]/m.pi))
        #df.loc[i, 'ESD'] = r*2
        df.loc[i, 'Biovolume'] = (4/3) * m.pi * (r**3)

    return df

def remove_UVP_noVal(df):
    '''
    Objective: first, remove all UVP ROI's with status = unclassified. For the remaining ROI's with annotation,
    calculate the percentage of validated ROI's. if the sample has <95% validation, discard
    '''
    df['Category'].replace('', np.nan, inplace=True)
    df = df.dropna(subset = ['Category']).reset_index(drop=True)

    if len(df.dropna(subset=['Category'])) == 0:
        df = pd.DataFrame()
        return df
    else:
        val = df.groupby(['Sample']).apply(lambda x: pd.Series({'Validation_percentage':len(x.dropna(subset=['Category'])[x.dropna(subset=['Category'])['Annotation'].isin(['validated'])].ROI) /
                                                                            len(x.dropna(subset=['Category']).ROI)}))
        for s in df.Sample.unique():
            if val.loc[s, 'Validation_percentage'] < 0.95:
                index = df[df['Sample'] == s].index
                df.drop(index, inplace=True)
            else:
                df.loc[df['Sample']==s, 'Validation_percentage'] = val.loc[s, 'Validation_percentage']

        return df

def biovol_func(df, instrument, keep_cat='none'):
    """
    Objective: calculate biovolume (in cubic micrometers) of each object  following
    the volume of  a sphere.
        :param df: a STANDARDIZED dataframe (output of step2_standardize_projects)
    This function by default removes objects based on 'Category':
       instrument = 'IFCB': removes 'beads' and 'artefact'
       instrument = 'Zooscan': removes 'detritus' and 'artefact'
       instrument = 'UVP': removes 'detritus' and 'artefact'
    :param keep_cat: list of strings that determines if categories that should be removed are kept
    :return: the original dataframe Biovolume in cubic micrometers
    """
    import math as m

    path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    path_to_taxonomy = str(Path(cfg['git_dir']).expanduser()) + '/ancillary/plankton_annotated_taxonomy.xlsx'
    taxonomy_df = pd.read_excel(path_to_taxonomy)
    df['Category'] = df['Category'].fillna('')
    df['Cat_remove'] = pd.merge(df, taxonomy_df, how='left', on=['Category'])['EcoTaxa_hierarchy']
    df['Cat_remove'] = df['Cat_remove'].fillna('')
    for n, i in enumerate(df['Cat_remove']):
        if i == '':
            df.loc[n, 'Cat_remove'] = df.loc[n, 'Category']
        else:
            pass

    if instrument == 'IFCB':
        cat_remove = ['beads','detritus', 'artefact', 'bubble']
        if keep_cat == 'none':
            cat_remove = cat_remove
        else:
            cat_remove = [x for x in cat_remove if x not in keep_cat]
        for i in cat_remove:
            try:
                df = df[df['Cat_remove'].str.contains(i) == False]
                df = df.reset_index(drop=True)
            except:
                pass
    elif (instrument == 'Zooscan'):
        cat_remove = ['artefact','detritus', 'plastic']
        if keep_cat == 'none':
            cat_remove = cat_remove
        else:
            cat_remove = [x for x in cat_remove if x not in keep_cat]
        for i in cat_remove:
            try:
                df = df[df['Cat_remove'].str.contains(i) == False]#df[~df.Category.str.contains('|'.join(i))]
                df = df.reset_index(drop=True)
            except:
                pass
    elif (instrument == 'UVP'):
        cat_remove = ['artefact','detritus', 'plastic']
        if keep_cat == 'none':
            cat_remove = cat_remove
        else:
            cat_remove = [x for x in cat_remove if x not in keep_cat]
        for i in cat_remove:
            try:
                df = df[df['Cat_remove'].str.contains(i) == False]#df[~df.Category.str.contains('|'.join(i))]
                df = df.reset_index(drop=True)
            except:
                pass
    df['Biovolume'] = df['Area'].apply(lambda x: (4/3)* m.pi * (m.sqrt((x/m.pi)) **3) ) # convert area to ESD, then calculate biovolume
    df = df.drop(columns=['Cat_remove']).reset_index(drop=True)

    return df




def depth_binning_func(depth='Depth_min', depth_bins=[0, 25, 50, 100, 200, 500, 1000, 3000, 8000]):
    # create depth bins based on Jessica's suggestion (Jessica's suggestions as default)
    import pandas as pd
    depth_bin = pd.cut(x=depth, bins=depth_bins, include_lowest=True)
    # create a categorical variable with the middle depth of each bin
    midDepth_bin = []
    for i in range(0, len(depth)):
        midDepth_bin.append(depth_bin.iloc[i].mid)
    return midDepth_bin



    #6.3 by station (geographical location)
def gridding_func(st_increment, lat= 'Latitude', lon= 'Longitude'):
    import numpy as np
    import pandas as pd
    # create a categorical variable that will serve as the station ID based on binned lat and lon
    # 1) create arrays that will serve to bin lat and long with 1 degree increments
    lat_increments = np.arange(np.floor(min(lat) - st_increment),
                               np.ceil(max(lat) + st_increment), st_increment)
    lon_increments = np.arange(np.floor(min(lon) - st_increment),
                               np.ceil(max(lon) + st_increment), st_increment)
    # 2) create bins
    lat_bin = pd.cut(x=lat, bins=lat_increments, include_lowest=True)

    lon_bin = pd.cut(x=lon, bins=lon_increments, include_lowest=True)
    # 3) get middle lat after assigning each sample to a lat/lon bin
    midLat_bin = []
    for i in range(0, len(lat_bin)):
        midLat_bin.append(lat_bin.iloc[i].mid)
    midLon_bin = []
    for i in range(0, len(lon_bin)):
        midLon_bin.append(lon_bin.iloc[i].mid)
    # create a station ID
    Station_ID = []
    for i in range(0, len(midLon_bin)):
        Station_ID.append(str(midLat_bin[i]) + '_' + str(midLon_bin[i]))

    return Station_ID, midLat_bin, midLon_bin

def group_gridded_files_func(instrument, already_gridded= 'N'):
    """
    Objective: assign a label to each gridded dataset, so when the compilation of files occurs, it is only done with a group of data and prevents
    the creation of a huge file
    """
    #first, create series of numbers that break the globe into 15x15 degree cells:
    lat_left = np.arange(-90, 90+15, 15, dtype=int)
    lat_right = np.arange(-75,105 + 15, 15, dtype=int)
    lat_int = pd.IntervalIndex.from_arrays(lat_left, lat_right, closed='both')
    lon_left = np.arange(-180, 180+15, 15, dtype=int)
    lon_right = np.arange(-165, 195+15, 15, dtype=int)
    lon_int = pd.IntervalIndex.from_arrays(lon_left, lon_right, closed='both')
    grid_list = []
    # create a dictionary that has unique keys for each grid interval
    #grid_dict = {}
    #for x in range(0, len(lon_int)):
        #for y in range(0, len(lat_int)):
            #grid_dict[(str(x) + '_' + str(y))] = {'lon': lon_int[x], 'lat': lat_int[y]} ## however, this might not be necessary
    # try now to assing lat & lon grid numbers to the dataframe directly
    file_list = proj_id_list_func(instrument, data_status='gridded')
    for i in tqdm(file_list):
        if already_gridded == 'N':
            print('generating subsets of data for 15 degree grids from ' + i)
        else:
            print('extracting big grid labels for subsetting from ' + i)
        filename = i.split('/')[-1]
        dirpath = i.replace(filename, '')
        df = pd.read_csv(i, header = 0)
        df['lat_grid'] = None
        df['lon_grid'] = None
        for la in df['midLatBin'].unique():
            df.loc[df['midLatBin'] == la, 'lat_grid']= str(lat_int[lat_int.contains(la) == True][0].left) + '-' + str(lat_int[lat_int.contains(la) == True][0].right)
        for lo in df['midLonBin'].unique():
            df.loc[df['midLonBin'] == lo, 'lon_grid'] = str(lon_int[lon_int.contains(lo) == True][0].left) + '-' + str(lon_int[lon_int.contains(lo) == True][0].right)
        # combine these numbers to form a grid identifier
        df['grid_id'] = df.lat_grid.astype(str).str.cat(df.lon_grid.astype(str), sep='_')
        df_subgrouped = [x for _, x in df.groupby('grid_id')]
        for s in df_subgrouped:
            grid_list.append(str(s['grid_id'].unique()[0]))
            if already_gridded == 'N':
                s.to_csv(str(dirpath) + 'grid_N_'+str(s['grid_id'].unique()[0]) + '_'+ filename, index=False)
            #os.remove(i)
    grid_list_unique = [*set(grid_list)]
    grid_list_unique = ['N_' + s for s in grid_list_unique]
    return grid_list_unique

#def date_binning_func_old(date, time, lat, lon, group_by= 'yyyymm'): # consider adding a day/night column, this will have to consider latitude and month

    #import pandas as pd
    #import calendar
    #import numpy as np
    #import datetime as dt
    #import pytz # to define the timezone of the time info
    #from datetime import datetime
    #from astral.sun import sun #run pip install astral
    # we also need to define timezone based on lat/lon. See https://www.geeksforgeeks.org/get-time-zone-of-a-given-location-using-python/
    #from timezonefinder import TimezoneFinder # run pip install timezonefinder
    #obj = TimezoneFinder()
    # to assign day/night to an ROI, first create a LocationInfo element. See instructions: https://astral.readthedocs.io/en/latest/
    #from astral import LocationInfo

    #calendar.setfirstweekday(6)
    #date = date.astype(str)
    #time = time.astype(str)
    #dateTime = date+time
    #date_bin = pd.to_datetime(date, format='%Y%m%d')
    #year = pd.DatetimeIndex(date_bin).year
    #month = pd.DatetimeIndex(date_bin).month #zfill(2).
    #day = pd.DatetimeIndex(date_bin).day
    #week_of_year = pd.DatetimeIndex(date_bin).isocalendar().week #zfill(2).
    #week_bin = pd.Series()
    #week_bin = week_bin.str.cat[]
    #week_bin = str(year) + '_' + str(month).zfill(2) + '_wk'+str(week_of_year).zfill(2)

    #light_cond = []
    #create a merged timestamp with date and time, useful for getting day and night info, which is done in lines 189-217
    #for n, i in enumerate(Time):
        #try:
            #Time_bin = pd.to_datetime(i, format='%Y%m%d%H%M%S')
            #utc = pytz.UTC #this and next two lines: set the time info to UTC format. Necessary to compare it with the sunrise/sunset information
            #Time_bin=Time_bin.replace(tzinfo=utc)
            #l = LocationInfo()
            #l.name = 'station'
            #l.region = 'region'
            #l.timezone = obj.timezone_at(lng=lat[n], lat=lon[n])
            #l.latitude = lat[n]
            #l.longitude = lon[n]
            # use the created l object to get times of sunrise and sunset for that location
            #try:
                #s = sun(l.observer, date=Time_bin)
                #sunrise = s['sunrise']
                #sunset = s['sunset']
                #if sunrise <= Time_bin <= sunset:
                    #light_cond.append('day')
                #else:
                    #light_cond.append('night')
            #except: # this is necessary because at high latitudes we might have issues defining day/night i.e Tara polar oceans: ValueError: Sun never reaches 6 degrees below the horizon, at this location.
                #light_cond.append('high_lat')
        #except:
            #light_cond.append(float('nan'))
    #if group_by == 'yyyy':
        #date_bin = date_bin.dt.strftime('%Y')
    #elif group_by == 'yyyymm':
        #date_bin = date_bin.dt.strftime('%Y%m')
    #elif group_by == 'mm':
        #date_bin = date_bin.dt.strftime('%m')
    #elif group_by == 'week':
        #def get_week_of_month(year, month, day):
            #x = np.array(calendar.monthcalendar(year, month))
            #week_of_month = np.where(x == day)[0][0]
            #if week_of_month == 0:
                #week_of_month = week_of_month +1
            #elif week_of_month == 5:
                #week_of_month = week_of_month - 1
            #return (week_of_month)
        #day = date_bin.dt.day
        #month = date_bin.dt.month
        #year = date_bin.dt.year
        #for i in range(0, len(date_bin)):
            #week=get_week_of_month(year[i], month[i], day[i])
            #date_bin[i] = str(year[i]) + str(month[i]).zfill(2) + '_wk' + str(week)
    #elif group_by == 'None':
        #date_bin == date_bin

    #return date_bin, light_cond

def date_binning_func(df, group_by= 'yyyymm',day_night=False): # , ignore_high_lat=True, consider adding a day/night column, this will have to consider latitude and month
    """
    Objective: reduce the date information, so that the data can be binned by month, year, or month and year. Also create a column that assigns a 'day' or 'night' category to each ROI
    :param date: column of a  standardized dataframe containing date information ('Sampling_date' in standardized ecotaxa projects)
    with format yyyymmdd (UTC)
    :param time: column of a  standardized dataframe containing time ('Sampling_time' in standardized ecotaxa projects)
    with format HHMMSS  (UTC)
    :param lat:
    :param lon:
    :param group_by: range of time to bin the data, options are: 'year', 'month' and 'year_month'
    :return:
    """

    # we also need to define timezone based on lat/lon. See https://www.geeksforgeeks.org/get-time-zone-of-a-given-location-using-python/
    from timezonefinder import TimezoneFinder # run pip install timezonefinder
    obj = TimezoneFinder()
    df = df.dropna(subset=['Latitude', 'Longitude']).reset_index(drop=True)
    df['Sampling_time']=df['Sampling_time'].astype(str)
    df['Sampling_time'] = df['Sampling_time'].apply(lambda x: x.ljust(6, '0'))
    date = df['Sampling_date'].astype(str)

    #time = df['Sampling_time'].astype(str).apply(lambda x: ('0' + x) if len(x) < 6 else x)
    time = df['Sampling_time'].astype(str)
    lat = df['Latitude']
    lon = df['Longitude']
    date_bin = pd.to_datetime(date, format='%Y%m%d')
    year = np.char.array(pd.DatetimeIndex(date_bin).year.values)
    month = np.char.array(pd.DatetimeIndex(date_bin).month.values).zfill(2)
    week_of_year = np.char.array(pd.DatetimeIndex(date_bin).isocalendar().week.values).zfill(2) #zfill(2).

    if group_by == 'yyyy':
        df['date_bin'] = (year).astype(str)
    elif group_by == 'yyyymm':
        df['date_bin'] = (year + b'_' + month).astype(str)
    elif group_by == 'mm':
        df['date_bin'] = (month).astype(str)
    elif group_by == 'week':
        df['year_week'] = (year + b'_' + week_of_year).astype(str)
        df['month'] = (month).astype(str)
        week_dict = {key: None for key in df['year_week'].unique()}
        for i in df['year_week'].unique():
            week_dict[i] = list(df['month'].where(df['year_week'] == i).dropna().unique())  # .astype(int)
            if len(week_dict[i]) == 1:
                week_dict[i] = week_dict[i][0]
            else:
                week_dict[i] = list(map(int, week_dict[i]))
                week_dict[i] = str(int(np.round(np.mean(week_dict[i])))).zfill(2)
        df['month'] = df['year_week'].map(week_dict)
        mwk_df = df['year_week'].str.split("_", expand=True)
        df['date_bin'] = mwk_df[0] + '_' + df['month'] + '_wk'+ mwk_df[1]
        df = df.drop(['year_week', 'month'], axis=1).reset_index(drop=True)

    elif group_by == 'None':
        df['date_bin'] == str(date_bin)
    if day_night == True:
        time_bin = pd.to_datetime(time, format='%H%M%S')
        day = np.char.array(pd.DatetimeIndex(date_bin).day.values).zfill(2)
        hour = np.char.array(pd.DatetimeIndex(time_bin).hour.values).zfill(2)  # zfill(2).
        minute = np.char.array(pd.DatetimeIndex(time_bin).minute.values).zfill(2)
        second = np.char.array(pd.DatetimeIndex(time_bin).second.values).zfill(2)

        df['dateTime'] = (year + month + day + hour + minute + second).astype(str)
        df['dateTime']= df['dateTime'].apply(lambda x: pd.to_datetime(x, format='%Y%m%d%H%M%S'))
        #df['date_new'] = df['dateTime'].dt.strftime('%Y/%m/%d')
        df['date_new'] = pd.to_datetime(df['dateTime']).dt.date
        df['time_new'] = pd.to_datetime(df['dateTime']).dt.time
        light_cond=[]
        for i in range(0, len(df)):
            light_cond.append(d_n_calc(df.loc[i, 'date_new'], df.loc[i, 'time_new'], df.loc[i, 'Latitude'], df.loc[i, 'Longitude']))
        df['light_cond'] = light_cond
        df = df.drop(columns=['date_new', 'time_new']).reset_index(drop=True)


    return df
        # section below deprecated (5/9/2023)
        # dateTime= pd.to_datetime(all_time_info, format='%Y%m%d%H%M%S')
        #light_cond = []

        #create a merged timestamp with date and time, useful for getting day and night info, which is done in lines 189-217
       # if ignore_high_lat==False:
            #for n, t in enumerate(dateTime):
                #try:
                    #utc = pytz.UTC #this and next two lines: set the time info to UTC format. Necessary to compare it with the sunrise/sunset information
                    #t=t.replace(tzinfo=utc)
                    #l = LocationInfo()
                    #l.name = 'station'
                    #l.region = 'region'
                    #l.timezone = obj.timezone_at(lng=lon[n], lat=lat[n])
                    #l.latitude = lat[n]
                    #l.longitude = lon[n]
                    #use the created l object to get times of sunrise and sunset for that location
                    #try:
                        #s = sun(l.observer, date=t)
                        #sunrise = s['sunrise']
                        #sunset = s['sunset']
                        #if sunrise <= t <= sunset:
                            #light_cond.append('day')
                        #else:
                            #light_cond.append('night')
                    #except ValueError as e: # this is necessary because at high latitudes we might have issues defining day/night i.e Tara polar oceans: ValueError: Sun never reaches 6 degrees below the horizon, at this location.
                        #if e == 'Sun never reaches 6 degrees below the horizon, at this location':
                            #light_cond.append('high_lat')
                #except:
                    #print('error assigning day/night, check row' + str(n))
                    #light_cond.append(float('nan'))
        # an alternative way to get this info without worrying about Astral's error for polar latitudes is to use suncalc
        #else:
            #for n, t in enumerate(dateTime):
                #try:
                    #daylight_dict = get_times(t, lon[n], lat[n])
                    #if daylight_dict['sunrise'] <= t <= daylight_dict['sunset']:
                        #light_cond.append('day')
                    #else:
                        #light_cond.append('night')
                #except:
                    #light_cond.append(float('nan'))

        #df['light_cond']= light_cond

    #return df


