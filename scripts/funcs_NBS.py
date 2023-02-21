# functions to calculate normalized biomass size spectra, based on ecopart size bins and a logarithmic scale
# data needs to be gridded and temporally and spatially binned. proj_nbs_func summarizes all the steps





  #6.1 by size bins. Inputs: a column of the dataframe with biovolume, and the number of log bins to use.
    #returns two lists : the sizeClasses bins (categorical) and range_size_bins (float)
def size_binning_func(df_subset):
    """
    Objective: bin the dataframes (parsed by stations and/or dates and/or depths)
    :param biovolume: column of a dataframe that contains the biovolume (in cubic micrometers)
    :return: a binned dataframe by sizes, containing the range of the size bin for each
    """
    import numpy as np
    from pathlib import Path
    import pandas as pd
    import yaml
    # open dataframe with the size bins
    path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
    # open the metadata of the standardized files
    with open(path_to_config, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    path_to_bins = str(Path(cfg['raw_dir']).expanduser() / 'ecopart_size_bins.tsv')
    bins_df = pd.read_csv(path_to_bins)
    # create a categorical variable, which are bins defined by the numbers on the sizeClasses list
    # Assign bin to each data point based on biovolume, append bin range and bin class to the dataframe
    df_subset.loc[:, 'sizeClasses']= pd.cut(x=df_subset['Biovolume'], bins=bins_df['biovol_um3'], include_lowest=True)# size classes defined by biovolume
    range_size_bin = []
    for i in range(0, len(df_subset['sizeClasses'])):
        range_size_bin.append(df_subset.loc[i, 'sizeClasses'].length)
    df_subset.loc[:, 'range_size_bin'] = range_size_bin
    bins = df_subset['sizeClasses'].unique()
    bin_width = []# list that will contain the width of each size_bin
    for i in range(0, len(bins.categories)):
        bin_width.append(bins.categories[i].length)

    d = {'sizeClasses':bins.categories, 'range_size_bin': bin_width}
    df_bins = pd.DataFrame(d)
    # obtain the size of each size class bin and append it to the stats dataframe
    # unfortunately the Interval type produced by pd.cut is hard to work with. So they need to be modified

    return df_subset, df_bins


# 7)calculate x and y for NBSS, this includes adding total biovolume per size bin for each station and depth bin,
# inputs: a dataframe and the volume sampled of each (in cubic meters). other inputs are the column names
# from the original dataset that want to be retained (i.e;  proj ID, station ID, depth, size classes,
# range of size classes, biovolume, lat & lon ) plus the variables that will be used to group the data by
# stations, depths and size classes
# using 5 ml for IFCB
def NB_SS_func(df_binned, df_bins, ignore_depth = 'no', sample_id='Sample', dates='date_bin', station= 'Station_location', depths = 'midDepthBin',\
               min_depth = 'Min_obs_depth', max_depth= 'Max_obs_depth', light= 'light_cond', size_range= 'range_size_bin', sizeClasses= 'sizeClasses', biovolume='Biovolume',\
               lat='midLatBin', lon='midLonBin', vol_filtered='Volume_imaged'):
    """
    
    :param df_binned:
    :param df_bins:
    :param ignore_depth:
    :param dates:
    :param station:
    :param depths:
    :param size_range:
    :param sizeClasses:
    :param biovolume:
    :param lat:
    :param lon:
    :param vol_filtered:
    :return:
    """
    import numpy as np
    import pandas as pd
    if ignore_depth == 'no':
        # create a dataframe with summary statistics for each station, depth bin and size class
        # these column names should all be the same, since the input is a dataframe from the 'binning' and 'biovol' functions
        # group data by bins
        #test = df_binned.groupby([sizeClasses])

        NBS_biovol_df = df_binned.groupby([dates, station, light, depths, sizeClasses]).agg(
            {sample_id:'first', depths:'first', size_range:'first', lat: 'first', lon: 'first', vol_filtered: 'first', min_depth: 'first', max_depth:'first',
             biovolume:['sum', 'count', 'mean'] }).reset_index()
        df_vol = df_binned.groupby([dates, station, light, depths]).apply(lambda x: pd.Series({'cumulative_volume': x[['Sample', 'Min_obs_depth', 'Max_obs_depth','Volume_imaged']].drop_duplicates().Volume_imaged.sum()})).reset_index()

        # remove bins that have zero values
        # standardize by volume sample
    elif ignore_depth == 'yes':
        NBS_biovol_df = df_binned.groupby([dates, station, light, sizeClasses]).agg(
            {sample_id:'first', size_range:'first', lat: 'first', lon: 'first', vol_filtered: 'first', min_depth: 'first', max_depth:'first',
             biovolume:['sum', 'count' , 'mean']}).reset_index()
    # cumulative volume for 1x1 degree bin. For UVP: volume imaged is the volume of a picture. Consider sample ID and depths, to be abble to get the cumulative volume imaged for a grid
        df_vol = df_binned.groupby([dates, station, light]).apply(lambda x: pd.Series({ 'cumulative_volume': x[['Sample', 'Min_obs_depth', 'Max_obs_depth', 'Volume_imaged']].drop_duplicates().Volume_imaged.sum()})).reset_index() # don't understand why sampl

    NBS_biovol_df['cumulative_vol'] = df_vol.loc[0, 'cumulative_volume']


    NBS_biovol_df.columns = NBS_biovol_df.columns.map('_'.join).str.removesuffix("first")
    NBS_biovol_df.columns = NBS_biovol_df.columns.str.removesuffix("_")
    NBS_biovol_df = NBS_biovol_df[NBS_biovol_df['Biovolume_count'] != 0]
    NBS_biovol_df['NBSS'] = (NBS_biovol_df['Biovolume_sum']) / (NBS_biovol_df['cumulative_vol']) / (NBS_biovol_df[size_range])

    #two lines below deprecated as of 1/31/2023, use new biovolume:
    #NBS_biovol_df['NBSS'] = (NBS_biovol_df['Biovolume_sum']) / (np.nansum(list(set(NBS_biovol_df[vol_filtered])))) / (NBS_biovol_df[size_range])
    #NBS_biovol_df[vol_filtered] = np.nansum(list(set(NBS_biovol_df[vol_filtered]))) # this is in the case there are different volumes analyzed based on binning, they have to be summed

    # create two more columns with the parameters of normalized size spectra,:
    NBS_biovol_df['logNBSS'] = np.log(NBS_biovol_df['NBSS'])
    NBS_biovol_df['logSize'] = np.log(NBS_biovol_df[size_range])

    # now calculate the log biovolume for the df_bins dataframe
    df_bins['logSize'] = np.log(df_bins[size_range])
    # merge the two dataframes
    binned_NBS= pd.merge(df_bins, NBS_biovol_df, how='left', on=[sizeClasses, size_range, 'logSize'])
    # now fill the columns of date, station, lat/lon, project ID and volume
    for i in [dates, station, lat, lon, vol_filtered]:
        binned_NBS[i] = binned_NBS[i].fillna(NBS_biovol_df[i].unique()[0])
    # let's do the thresholding here:
    NBS_binned_thres = threshold_func(binned_NBS)

    return NBS_binned_thres

# 7)  create a function to remove data based on JO's comments and Haentjens THIS IS WORK IN PROGRESS:
# low threshold will be filtered because they are smaller than the next standardized size spectrum
# high threshold: obtain the highest size class, then iterate from highest to lowest and remove the large size bins
# that have less than an (instrument specific) threshold
# imputs: the data frame, and the column that contains the NBSS
def threshold_func(binned_data):
    """
    Objective: remove size bins based on Fabien's approach: (lower limit: max NBSS, upper limit: 3 NBSS with Nans)
    :param df: a binned dataframe with NBSS already calculated
    :return: a dataset without the bins that contain inaccurate NBSS
    """
    import numpy as np
    import pandas as pd

    #upper threshold based on max NBSS
    list_max_NBSS = binned_data.NBSS.tolist()
    binned_data_filt = binned_data.loc[list_max_NBSS.index(np.nanmax(list_max_NBSS)):len(binned_data)] #a ask if its better to just remove all data
    binned_data_filt = binned_data_filt.reset_index(drop=True)

    #lower threshold based on three consecutive size bins with nans
    for n, i in enumerate(np.isnan(binned_data_filt['NBSS'])):
        if (i == True) and (np.isnan(binned_data_filt.loc[n+1, 'NBSS']) == True) and (np.isnan(binned_data_filt.loc[n+2, 'NBSS']) == True):
            binned_data_filt = binned_data_filt.loc[0:n-1]
            break

    return binned_data_filt

def linear_fit_func(df1, depth_bin = False):
    """
    Objective: perform linear fit of NBSS. method of linear fit might change, current one is Least squares regression
    :param df: a dataframe with size bins and NBSS belonging to only one sample (one unique station, date and depth)
    from https://www.edureka.co/blog/least-square-regression/
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    #remove Nans
    df1 = df1.drop(columns=['NBSS_std'])
    df1 = df1.dropna()

    # extract x and y
    X = np.log(df1['range_size_bin'].values)#.reshape(-1, 1)
    Y = np.log(df1['NBSS_mean'].values)#.reshape(-1, 1)
    # Mean X and Y
    mean_x = np.mean(X)
    mean_y = np.mean(Y)
    # Total number of values
    n = len(X)
    # Model 1 -linear model : Y = m*X + b
    numer = 0 # numerator
    denom = 0 # denominator
    for i in range(n):
        numer += (X[i] - mean_x) * (Y[i] - mean_y)
        denom += (X[i] - mean_x) ** 2
    m = numer / denom # slope
    b = mean_y - (m * mean_x) # intercept
    #print("Coefficients")
    #print(m, c)
    # Calculating Root Mean Squares Error and R2 score
    rmse = 0 # root mean square error
    ss_tot = 0 # total sum of squares
    ss_res = 0 # total sum of squares of residuals
    for i in range(n):
        y_pred = b + m * X[i]
        rmse += (Y[i] - y_pred) ** 2
        ss_tot += (Y[i] - mean_y) ** 2
        ss_res += (Y[i] - y_pred) ** 2
    rmse = np.sqrt(rmse / n)
    #print("RMSE")
    #print(rmse)
    R2 = 1 - (ss_res / ss_tot)
    #print("R2 Score")
    #print(R2)

    # generate dataframe and append the results in the corresponding variable
    lin_fit = pd.DataFrame()
    if depth_bin == True:
        lin_fit.loc[0, 'depth'] = df1.loc[0, 'midDepthBin']
    lin_fit.loc[0, 'year'] = str(df1.loc[0, 'year'])
    lin_fit.loc[0, 'month'] = str(df1.loc[0, 'month'])

    lin_fit.loc[0, 'latitude'] = df1.loc[0, 'midLatBin']
    lin_fit.loc[0, 'longitude'] = df1.loc[0, 'midLonBin']
    lin_fit.loc[0, 'light_condition'] = df1.loc[0, 'light_cond']

    lin_fit.loc[0, 'min_depth'] = df1.loc[0, 'Min_obs_depth']
    lin_fit.loc[0, 'max_depth'] = df1.loc[0, 'Max_obs_depth']

    lin_fit.loc[0, 'slope'] = m
    lin_fit.loc[0, 'intercept'] = b
    lin_fit.loc[0, 'rmse'] = rmse
    lin_fit.loc[0, 'r2'] = R2
    lin_fit.loc[0, 'NBSS_sample_size'] = max(df1['NBSS_count'])

    return (lin_fit)

def parse_NBS_linfit_func(df, parse_by=['Station_location', 'date_bin', 'light_cond'], depth_bin=False):
    """
    Objective: parse out any df that come from files that include multiple stations, dates and depths.
    dataframe needs to be already binned. Steps are: 1) parse the data 2) bin ROIs by size  3) calculate the NBS
    (threshold included) and 4) create a dataframe with the Least squares regression result
    :param df: a standardized, binned dataframe
    :param parse_by: list of columns that will be used to parse the dataframe
    :return: a binned, tresholded NBS dataframe
    """
    import pandas as pd
    if depth_bin == True:
        parse_by.append('midDepthBin')
        df['midDepthBin'] = df['midDepthBin'].astype(str)
    df['date_bin'] = df['date_bin'].astype(str)
    NBSS_binned_all = pd.DataFrame()  # NBSS dataset
    for st in list(set(df[parse_by[0]])):
        # print(st)
        df_st_subset = df[df[parse_by[0]] == st].reset_index(drop=True)
        for t in list(set(df_st_subset[parse_by[1]])):
            # print(t)
            df_st_t_subset = df_st_subset[df_st_subset[parse_by[1]] == t].reset_index(drop=True)
            for l in list(set(df_st_subset[parse_by[2]])):
                df_st_t_l_subset = df_st_subset[df_st_subset[parse_by[2]] == l].reset_index(drop=True)
                if depth_bin == True:
                    for d in list(set(df_st_t_l_subset[parse_by[3]])):
                        df_st_t_l_d_subset = df_st_t_l_subset[df_st_t_l_subset[parse_by[3]] == d].reset_index(drop=True)
                        df_binned, df_bins = size_binning_func(df_st_t_l_d_subset)  # create size bins
                        NBS_biovol_df = NB_SS_func(df_binned, df_bins, ignore_depth='no', dates='date_bin',
                                                    sample_id='Sample',
                                                    station='Station_location',
                                                    depths='midDepthBin',
                                                    min_depth='Min_obs_depth', max_depth='Max_obs_depth',
                                                    size_range='range_size_bin',
                                                    sizeClasses='sizeClasses',
                                                    biovolume='Biovolume', lat='midLatBin', lon='midLonBin',
                                                    vol_filtered='Volume_imaged')  # calculate NBSS
                        #lin_fit = linear_fit_func(NBS_biovol_df, depth_bin=True)  2/20/2023: linear fit is going to be performed on the monthly 1 degree average of the NBSS
                        NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])
                else:
                    # print ('getting NBS for station ' + st + 'and date' + t)
                    df_binned, df_bins = size_binning_func(df_st_t_l_subset)  # create size bins
                    NBS_biovol_df = NB_SS_func(df_binned, df_bins, ignore_depth='yes', dates='date_bin',
                                            sample_id='Sample',
                                            station='Station_location',
                                            size_range='range_size_bin', sizeClasses='sizeClasses',
                                            min_depth='Min_obs_depth', max_depth='Max_obs_depth',
                                            biovolume='Biovolume', lat='midLatBin', lon='midLonBin',
                                            vol_filtered='Volume_imaged')  # calculate NBSS
                    # Create new dataset with slope, intercept, and R2 score for each linear fit
                    # lin_fit = linear_fit_func(NBS_biovol_df, depth_bin=True)  2/20/2023: linear fit is going to be performed on the monthly 1 degree average of the NBSS
                    NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])

    # after FULL DATAFRAME is assembled, get monthly, 1 degree bin average:
    NBS_avg, NBSS_1a = NBSS_stats_func(NBSS_binned_all)
    # Next, obntain the linear regression with the new, averaged NBSS
    lin_fit_data = pd.DataFrame()  # linear fit (least-squares) dataset

    for st in list(set(NBS_avg[parse_by[0]])):
        # print(st)
        NBS_avg_st_subset = NBS_avg[NBS_avg[parse_by[0]] == st].reset_index(drop=True)
        for t in list(set(NBS_avg_st_subset[parse_by[1]])):
            # print(t)
            NBS_avg_st_t_subset = NBS_avg_st_subset[NBS_avg_st_subset[parse_by[1]] == t].reset_index(drop=True)
            for l in list(set(NBS_avg_st_t_subset[parse_by[2]])):
                NBS_avg_st_t_l_subset = NBS_avg_st_t_subset[NBS_avg_st_t_subset[parse_by[2]] == l].reset_index(drop=True)

            if depth_bin == True:
                for d in list(set(NBS_avg_st_t_l_subset[parse_by[3]])):
                    NBS_avg_st_t_l_d_subset = NBS_avg_st_t_l_subset[NBS_avg_st_t_l_subset[parse_by[2]] == d].reset_index(drop=True)
                    lin_fit = linear_fit_func(NBS_avg_st_t_l_d_subset , depth_bin=True)
                    lin_fit_data = pd.concat([lin_fit_data, lin_fit])
            else:
                # print ('getting NBS for station ' + st + 'and date' + t 'and light condition' +l)
                lin_fit = linear_fit_func(NBS_avg_st_t_l_subset, depth_bin=False)
                lin_fit_data = pd.concat([lin_fit_data, lin_fit])

    return NBSS_1a, lin_fit_data


def NBSS_stats_func(df):
    """

    """
    from funcs_gridding import gridding_func
    import numpy as np
    bin_loc = input('Average NBSS by location? \n Enter Y/N ')
    if bin_loc == 'Y':
        st_increment = float(
            input('select the size of the spatial grid by which the NBSS will be averaged i.e for a 1x1 degree bin type 1 \n'))
        df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(st_increment, df['midLatBin'], df['midLonBin'])
    bin_date = input('Average NBSS by date? \n Enter Y/N ')
    if bin_date == 'Y':
        group_by = input(
            'input how will the NBSS will be averaged by date \n  yyyymm for month and year \n yyyy for year \n ')
        date_df = df['date_bin'].str.split("_", expand=True)
        if group_by == 'yyyy':
            df['date_bin'] = date_df[0]
        elif group_by == 'yyyymm':
            df['year_week'] = date_df[0] + '_' + date_df[2]
            df['month'] = date_df[1]
            week_dict = {key: None for key in df['year_week'].unique()}
            for i in df['year_week'].unique():
                week_dict[i] = list(df['month'].where(df['year_week'] == i).dropna().unique())  # .astype(int)
                if len(week_dict[i]) == 1:
                    week_dict[i] = week_dict[i][0]
                else:
                    week_dict[i] = list(map(int, week_dict[i]))
                    week_dict[i] = str(int(np.round(np.mean(week_dict[i])))).zfill(2)
            df['year'] = date_df[0]
            df['month'] = df['year_week'].map(week_dict)
            df['date_bin'] = date_df[0] + '_' + df['month']
    NBSS = df.filter(['range_size_bin', 'date_bin', 'month', 'year', 'Station_location', 'midLatBin', 'midLonBin','light_cond', 'Min_obs_depth', 'Max_obs_depth', 'NBSS'], axis=1)

    NBSS_avg = NBSS.groupby(['date_bin', 'Station_location', 'light_cond' 'range_size_bin']).agg(
            {'midLatBin':'first', 'midLonBin': 'first', 'year': 'first', 'month': 'first', 'Min_obs_depth':'first', 'Max_obs_depth':'first',
             'NBSS':[ 'count' , 'mean', 'std']}).reset_index()

    # generate a cleaner dataset (without grouping variables) to finally present 1a
    NBSS_1a=NBSS_avg
    NBSS_1a.columns = NBSS_1a.columns.map('_'.join).str.removesuffix("first")
    NBSS_1a.columns = NBSS_1a.columns.str.removesuffix("_")
    NBSS_1a= NBSS_1a.filter(['year', 'month', 'midLatBin', 'midLonBin', 'light_cond', 'Min_obs_depth', 'Max_obs_depth', 'range_size_bin','NBSS_mean', 'NBSS_std', 'NBSS_count'], axis=1)
    NBSS_1a= NBSS_1a.rename(columns={'midLatBin':'latitude', 'midLonBin': 'longitude', 'light_cond':'light_condition', 'Min_obs_depth':'min_depth', 'Max_obs_depth':'max_depth', 'range_size_bin': 'size_biovolume', 'NBSS_count':'bin_sample_size' })

    return NBSS_avg, NBSS_1a

def stats_linfit_func(df):
    """
    Objective: create summary statistics for the linear fit dataframes
    param df: a dataframe containing the results of the linear fits for each station, time and depth for a given imaging instrument
    """
    from funcs_gridding import gridding_func
    import numpy as np
    bin_loc =input('Group data by location? \n Enter Y/N ')
    if bin_loc == 'Y':
        st_increment = float(input('select the size of the spatial grid to group the data i.e for a 1x1 degree bin type 1 \n' ))
        df['Station_location'], df['midLatBin'], df['midLonBin'] = gridding_func(st_increment, df['Latitude'], df['Longitude'])
    bin_date =input('Group data by date? \n Enter Y/N ')
    if bin_date == 'Y':
        group_by= input('input how will the data be grouped by date \n  yyyymm for month and year \n yyyy for year \n ')
        date_df= df['Date'].str.split("_", expand=True)
        if group_by == 'yyyy':
            df['date_bin'] = date_df[0]
        elif group_by == 'yyyymm':
            df['year_week'] = date_df[0] + '_' + date_df[2]
            df['month'] = date_df[1]
            week_dict = {key:None for key in df['year_week'].unique()}
            for i in df['year_week'].unique():
                week_dict[i] =list(df['month'].where(df['year_week'] == i).dropna().unique()) #.astype(int)
                if len(week_dict[i]) == 1:
                    week_dict[i] = week_dict[i][0]
                else:
                    week_dict[i] = list(map(int, week_dict[i]))
                    week_dict[i] = str(int(np.round(np.mean(week_dict[i])))).zfill(2)
            df['year'] = date_df[0]
            df['month'] = df['year_week'].map(week_dict)
            df['date_bin'] = date_df[0] + '_' + df['month']
            df = df.drop(columns=['Date', 'year_week'])


    lin_fit_stats = df.groupby(['Station_location', 'date_bin']).agg(
        {'Project_ID': 'first', 'midLatBin': 'first', 'midLonBin': 'first', 'year':'first', 'month': 'first',
         'Slope': ['count', 'mean', 'std'], 'Intercept': ['count', 'mean', 'std'], 'R2': ['count', 'mean', 'std']}).reset_index()

    lin_fit_stats.columns = lin_fit_stats.columns.map('_'.join).str.removesuffix("first")
    lin_fit_stats.columns = lin_fit_stats.columns.str.removesuffix("_")
    lin_fit_stats= lin_fit_stats.rename(columns={'Slope_count': 'Sample_size'})
    lin_fit_stats= lin_fit_stats[lin_fit_stats.columns.drop(list(lin_fit_stats.filter(regex='count')))]

    return lin_fit_stats



# Plotting Values and Regression Line
#max_x = np.max(X) + 100
#min_x = np.min(X) - 100
# Calculating line values x and y
#x = np.linspace(min_x, max_x, 1000)
#y = c + m * x
# Ploting Line
#plt.plot(x, y, color='#58b970', label='Regression Line')
# Ploting Scatter Points
#plt.scatter(X, Y, c='#ef5423', label='Scatter Plot')
#plt.xlabel('logSize')
#plt.ylabel('logNBSS')
#plt.legend()
#plt.show()