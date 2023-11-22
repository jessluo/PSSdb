#script to change units for the Size spectra regression
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import yaml
from pathlib import Path
from glob import glob

path_to_config = Path('~/GIT/PSSdb/scripts/configuration_masterfile.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
# generate paths to standardized, gridded files, and the standardizer spreadsheets

path_to_1a_list= glob(str(Path(cfg['raw_dir']).expanduser() / 'NBSS_data') + '/*1a*')
NBSSpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data'

conversion_factor = np.float(input('default size units of PSSdb are in micrometers, please input the conversion factor you want to apply to change the units'))
new_units = input('input a string with the new units')

def linear_fit_func(df1, light_parsing = False, depth_parsing = False):
    """
    Objective: perform linear fit of NBSS. method of linear fit might change, current one is Least squares regression
    :param df: a dataframe with size bins and NBSS belonging to only one sample (one unique station, date and depth)
    from https://www.edureka.co/blog/least-square-regression/
    """
    df1['logNB'] = np.log10(df1['normalized_biovolume_mean'])
    df1['logSize'] = np.log10((df1['biovolume_size_class'].astype(float)*(conversion_factor**3)))

    df1['logPSD'] = np.log10(df1['normalized_abundance_mean'])
    df1['logECD'] = np.log10((df1['equivalent_circular_diameter_mean'].astype(float)*conversion_factor))
    # extract x and y
    def regression(df1, X_var, Y_var):
        '''
        :param Y_var: string that identifies the column header of the df1 dataframe to use as Y variable
        '''
        X = df1[X_var].values#.reshape(-1, 1)
        Y = df1[Y_var].values#.reshape(-1, 1)
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
        return m,b,rmse,R2
    m_NB, b_NB, rmse_NB, R2_NB = regression(df1, 'logSize','logNB')
    m_PSD, b_PSD, rmse_PSD, R2_PSD = regression(df1, 'logECD', 'logPSD')
    # generate dataframe and append the results in the corresponding variable
    lin_fit = pd.DataFrame()
    if depth_parsing == True:
        lin_fit.loc[0, 'depth'] = df1.loc[0, 'midDepthBin']

    lin_fit.loc[0, 'year'] = str(df1.loc[0, 'year'])
    lin_fit.loc[0, 'month'] = str(df1.loc[0, 'month'])
    lin_fit.loc[0, 'latitude'] = df1.loc[0, 'latitude']
    lin_fit.loc[0, 'longitude'] = df1.loc[0, 'longitude']
    if light_parsing==True:
        lin_fit.loc[0, 'light_condition'] = df1.loc[0, 'light_cond']
    lin_fit.loc[0, 'ocean'] = df1.loc[0, 'ocean']
    lin_fit.loc[0, 'min_depth'] = df1.loc[0, 'min_depth']
    lin_fit.loc[0, 'max_depth'] = df1.loc[0, 'max_depth']

    lin_fit.loc[0, 'NBSS_slope'] = m_NB
    lin_fit.loc[0, 'NBSS_intercept'] = b_NB
    lin_fit.loc[0, 'NBSS_rmse'] = rmse_NB
    lin_fit.loc[0, 'NBSS_r2'] = R2_NB

    lin_fit.loc[0, 'PSD_slope'] = m_PSD
    lin_fit.loc[0, 'PSD_intercept'] = b_PSD
    lin_fit.loc[0, 'PSD_rmse'] = rmse_PSD
    lin_fit.loc[0, 'PSD_r2'] = R2_PSD


    return (lin_fit)

grouping = ['year', 'month', 'latitude', 'longitude', 'ocean', 'min_depth', 'max_depth']
for i in path_to_1a_list:
    filename = i.split('/')[-1]
    instrument = filename.split('_')[0]
    version_date = filename.split('_')[-1]
    df = pd.read_csv(i, sep = ',')
    df_grouped = pd.merge(df, df.drop_duplicates(subset=grouping, ignore_index=True)[grouping].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left', on=grouping)
    lin_fit_full = pd.DataFrame()
    for g in df_grouped.Group_index.unique():
        df_subset = df_grouped[df_grouped.Group_index==g].reset_index(drop=True)
        lin_fit = linear_fit_func(df_subset)
        lin_fit_full = pd.concat([lin_fit_full, lin_fit]).reset_index(drop = True)
        lin_fit_full.to_csv(str(NBSSpath) + '/' + instrument + '_1b_Size-spectra-fit_'+new_units+'_'+version_date,index=False)


