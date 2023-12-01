# Objective: This step computes taxa-specific Normalized Biovolume Size Spectrum based on steps 0,1,2,3. Computation follows that of step 4, except that taxonomic groups are assigned to each ROI to compute taxa-specific products

# Python modules and path specifications:

# Standardization of taxonomic annotations
try:
    from funcs_standardize_annotations import *
except:
    from scripts.funcs_standardize_annotations import *
import ast
import yaml# requires installation of PyYAML package
from pathlib import Path
#open config file and extract parameters:
path_to_git=Path('~/GIT/PSSdb').expanduser()
path_to_config = path_to_git /'scripts'/'configuration_masterfile.yaml'
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_to_allometry=path_to_git/cfg['allometry_lookup']
path_to_taxonomy=path_to_git/cfg['annotations_lookup']
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import numpy as np

import statistics as st
import os
from datetime import datetime

from glob import glob
import shutil
# Config modules
try:
    from funcs_read import *
    from funcs_NBS import *
except:
    from scripts.funcs_read import *
    from scripts.funcs_NBS import *

from tqdm import tqdm
import matplotlib.pyplot as plt
# Configurations:
light_parsing = cfg['day_night']
depth_parsing = cfg['depth_binning']

bin_loc = cfg['st_increment_avg']
group_by= cfg['date_group_avg']
sensitivity = cfg['sensitivity']

# Workflow starts here
df_allometry=pd.read_excel(path_to_allometry)
df_taxonomy=pd.read_excel(path_to_taxonomy)
#1: Merge taxonomic and allometric lookup tables generated on step 2 to assign a taxonomic group for each ROI
confirmation = input("Starting the computation of class- and functional type-specifc size spectra.\nDo you wish to update the taxonomic lookup table ({})? Enter Y or N\nY if new annotations should be merged to the WORMS database (https://www.marinespecies.org/aphia.php?p=search)\n".format(str(path_to_taxonomy).replace(str(Path.home()),'~')))
if confirmation=='Y':
    # Add column to allometric look-up table to merge with annotations
    if 'Full_hierarchy' not in df_allometry.columns:
        new_categories=df_allometry.Taxon.unique()
        df_allometry_new=pd.concat(map(lambda hierarchy:annotation_in_WORMS(hierarchy.replace("_",'>')).assign(Category=hierarchy),new_categories))
        df_allometry_new['URL'] = df_allometry_new.WORMS_ID.apply( lambda id: 'https://www.marinespecies.org/aphia.php?p=taxdetails&id={}'.format( id.replace('urn:lsid:marinespecies.org:taxname:', '')) if len(id) else '')
        df_allometry =pd.merge(df_allometry,df_allometry_new,left_on='Taxon',right_on='Category',suffixes=['_biomass','_WORMS'])
        df_allometry = df_allometry.sort_values(['Type', 'Category'], ascending=[False, True]).reset_index(drop=True)
        print('Updating allometric groups with the taxonomic annotations lookup table')
        with pd.ExcelWriter(str(path_to_allometry), engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
            df_allometry.to_excel(writer,sheet_name='Data',  index=False)


    #Re-assign WORMS taxonomy for misclassified categories
    new_categories=df_taxonomy[(df_taxonomy.Full_hierarchy.astype(str)=='nan') & (df_taxonomy.PFT.astype(str)=='nan')].Category.unique() if 'PFT' in df_taxonomy.columns else df_taxonomy[(df_taxonomy.Full_hierarchy.astype(str)=='nan')].Category.unique()
    categories_dict={'Amy_Gony_Protoc':'Gonyaulacaceae','Centric':'Bacillariophyceae','Ciliates':'Ciliophora','Coccolithophore':'Coccolithophyceae','Copepod_nauplii':'Copepoda','Cryptophyte':'Cryptophyceae','Cyl_Nitz':'Bacillariaceae','Det_Cer_Lau':'Hemiaulaceae','Diatoms':'Bacillariophyceae','Guin_Dact':'Rhizosoleniaceae','Pennate':'Bacillariophyceae','Pronoctaluca':'Pronoctiluca','Rhiz_Prob':'Rhizaria','Scrip_Het':'Peridiniaceae','ciliate':'Ciliophora','coccolithophorid':'Coccolithophyceae','metaneuplii':'Copepoda','multiple-diatoms':'Bacillariophyceae','nauplii':'Copepoda','pennate':'Bacillariophyceae','pennate_Pseudo-nitszchia':'Pseudonitzschia','pennate_Thalassionema':'Thalassionema','pennate_morphotype1':'Bacillariophyceae','protist with spike':'Protozoa','tempChaetoceros contortus':'Chaetoceros','tempCoccosphaerales':'Coccosphaerales','tempCryptophyceae':'Cryptophyceae','tempKatodinium glaucum':'Katodinium','tempPrasinophyceae':'Prasinophyceae','Dinophyceae_morphotype_2':'Dinophyceae' }
    if len(new_categories):
        df_taxonomy_new = pd.concat( map(lambda hierarchy: annotation_in_WORMS(re.sub('[^A-Za-z0-9]+','>',hierarchy)).assign(Category=hierarchy) if hierarchy not in categories_dict.keys() else annotation_in_WORMS(re.sub('[^A-Za-z0-9]+','>',categories_dict[hierarchy])).assign(Category=hierarchy) , new_categories))
        df_taxonomy_new['URL'] = df_taxonomy_new.WORMS_ID.apply(lambda id: 'https://www.marinespecies.org/aphia.php?p=taxdetails&id={}'.format(id.replace('urn:lsid:marinespecies.org:taxname:', '')) if len(id) else '')
        df_taxonomy = pd.concat([df_taxonomy[df_taxonomy.Category.isin(df_taxonomy_new.Category)==False].reset_index(drop=True), df_taxonomy_new[[column for column in df_taxonomy.columns if column in df_taxonomy_new.columns]].reset_index(drop=True)],axis=0, ignore_index=True)
        df_taxonomy = df_taxonomy.sort_values(['Type', 'Category'], ascending=[False, True]).reset_index(drop=True)
        print('New taxonomic annotations found in the project. Updating the taxonomic annotations lookup table')
        with pd.ExcelWriter(str(path_to_taxonomy), engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
            df_taxonomy.to_excel(writer, sheet_name='Data', index=False)

    # Merge allometry and taxonomy to assign functional group
    if len(df_taxonomy[(df_taxonomy[['Taxon','PFT']].astype(str)=='nan').all(axis=1)]):
        df_taxonomy.loc[(df_taxonomy[['Taxon','PFT']].astype(str)=='nan').all(axis=1),['Taxon','PFT']]=df_taxonomy[(df_taxonomy[['Taxon','PFT']].astype(str)=='nan').all(axis=1)].apply(lambda x:df_allometry.dropna(subset=['Full_hierarchy']).loc[df_allometry.dropna(subset=['Full_hierarchy']).Full_hierarchy.loc[ [i for i, sub in df_allometry.dropna(subset=['Full_hierarchy']).Full_hierarchy.astype(str).items() if (sub in str(x.Full_hierarchy)) and len(sub)]].str.len().idxmax(),['Taxon','Functional_group_allometry']] if len([i for i, sub in df_allometry.dropna(subset=['Full_hierarchy']).Full_hierarchy.astype(str).items() if (sub in str(x.Full_hierarchy)) and len(sub)]) else df_allometry.dropna(subset=['Full_hierarchy']).loc[df_allometry.dropna(subset=['Full_hierarchy']).Full_hierarchy.loc[ [i for i, sub in df_allometry.dropna(subset=['Full_hierarchy']).Full_hierarchy.astype(str).items() if (sub in str(x.Full_hierarchy))]].str.len().idxmax(),['Taxon','Functional_group_allometry']] if len(df_allometry.dropna(subset=['Full_hierarchy']).Full_hierarchy.loc[ [i for i, sub in df_allometry.dropna(subset=['Full_hierarchy']).Full_hierarchy.astype(str).items() if (sub in str(x.Full_hierarchy))]]) else ['',''],axis=1)
        print('Updating the taxonomic annotations lookup table with allometric functional group')
        with pd.ExcelWriter(str(path_to_taxonomy), engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
            df_taxonomy.to_excel(writer, sheet_name='Data', index=False)
elif confirmation=='N':
    pass
else:
    print('Input not conformed. Quitting, please re-run the script and type optional input as it appears.')
    quit()
#2: Loop through instrument-specific gridded files and compute class- and functional types-specifc Normalized Biovolume Size Spectra
for instrument in os.listdir(Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir']):
    #first, break apart datasets by big global grids, to avoid making one HUGE file of all gridded datasets per instrument
    #get paths to gridded files
    file_list = glob(str(Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir']) + '*/**/' + instrument +'*_temp_binned_*.csv', recursive=True) #generate path and project ID's but ONLY for parsed data
    grid_list = [re.search('N_(.+?).csv', file) for file in file_list]#group_gridded_files_func(instrument, already_gridded='Y')
    currentMonth = str(datetime.datetime.now().month).rjust(2, '0')
    currentYear = str(datetime.datetime.now().year)
    NBSSpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data' / str('NBSS_ver_'+currentMonth+'_'+currentYear)


    biovol = 'Biovolume_area'
    print ('generating PSSdb biovolume products for ' +instrument+' based on ' + biovol)
    NBSS_binned_all,NBSS_binned_all_PFT = pd.DataFrame(),pd.DataFrame()  # NBSS dataset
    lin_fit_data,lin_fit_data_PFT = pd.DataFrame(),pd.DataFrame()
    for i in tqdm(grid_list):
        i = i.group(1)
        print ('grouping data and calculating Normalized Biovolume Size Specta for cell number ' + i)
        file_subset = [file for file in file_list if i in file]
        #  Read gridded standardized files and merge to taxonomy lookup table
        NBS_biovol_df= pd.concat(map((lambda path: (pd.read_csv(path))), file_subset)).reset_index(drop=True)
        NBS_biovol_df=pd.merge(NBS_biovol_df, df_taxonomy[['Category','Taxon','PFT']],how='left',on=['Category'])

        if instrument == 'IFCB':
            NBS_biovol_df = NBS_biovol_df[(~NBS_biovol_df['Sample'].str.contains('D20200204T143120_IFCB109')) | (~NBS_biovol_df['Sample'].str.contains('D20131110T014214_IFCB101')) | (~NBS_biovol_df['Sample'].str.contains('D20131110T105100_IFCB101'))| (~NBS_biovol_df['Sample'].str.contains('D20190825T074148_IFCB109'))].reset_index(drop=True)
        NBS_biovol_df, df_bins = size_binning_func(NBS_biovol_df, biovol)  # create size bins
        # Append Sieburth size class nomination to phytoplankton functional group
        NBS_biovol_df['PFT'] = np.where(NBS_biovol_df['PFT'] == 'Phytoplankton', pd.cut(NBS_biovol_df.ECD_mid,[2,20,200,20000],labels=['Pico','Nano','Micro']).astype(str)+NBS_biovol_df.PFT.astype(str).str.lower(), NBS_biovol_df.PFT)
        # Compute taxa- and function group-specific NBSS
        NBS_biovol_df, lin_fit = NB_SS_func(NBS_biovol_df, df_bins, biovol_estimate=biovol, sensitivity=sensitivity,group=['Taxon','PFT'])
        NBSS_binned_all = pd.concat([NBSS_binned_all, NBS_biovol_df])
        lin_fit_data = pd.concat([lin_fit_data, lin_fit])
        if len(NBS_biovol_df):
            # Compute functional type specific NBSS
            group=['date_bin', 'Station_location', 'midLatBin', 'midLonBin','Min_obs_depth', 'Max_obs_depth',  'PFT']
            NBS_biovol_df_PFT=NBS_biovol_df.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin','Min_obs_depth', 'Max_obs_depth', 'sizeClasses', 'size_class_mid','range_size_bin', 'ECD_mid', 'size_range_ECD', 'PFT']).apply(lambda x:pd.Series({'Validation_percentage':np.round(np.nansum((x.Validation_percentage.astype(float)*x.ROI_number_sum.astype(float))/x.ROI_number_sum.astype(float).sum()),2),'Biovolume_mean':np.nansum((x.Biovolume_mean.astype(float)*x.ROI_number_sum)/x.ROI_number_sum.sum()),'size_class_pixel':np.nanmean(x.size_class_pixel.astype(float)),'ROI_number_sum':x.ROI_number_sum.astype(float).sum(),'ROI_abundance_mean':np.nansum((x.ROI_abundance_mean.astype(float)*x.ROI_number_sum.astype(float))/x.ROI_number_sum.astype(float).sum()),'NB':x.NB.astype(float).sum(),'PSD':x.PSD.astype(float).sum(),'count_uncertainty':poisson.pmf(k=np.nansum((x.ROI_abundance_mean.astype(float)*x.ROI_number_sum.astype(float))/x.ROI_number_sum.astype(float).sum()), mu=np.nansum((x.ROI_abundance_mean.astype(float)*x.ROI_number_sum.astype(float))/x.ROI_number_sum.astype(float).sum())),'size_uncertainty':np.nanmean(x.size_uncertainty.astype(float)), 'logNB':np.log10(x.NB.astype(float).sum()),'logSize':np.log10(x.size_class_mid.astype(float).unique()[0]),'logPSD':np.log10(x.PSD.astype(float).sum()),'logECD':np.log10(x.ECD_mid.astype(float).unique()[0])})).reset_index().sort_values(['date_bin', 'Station_location', 'midLatBin', 'midLonBin','Min_obs_depth', 'Max_obs_depth','PFT','size_class_mid']).reset_index(drop=True)
            ## Add empty size classes and threshold
            NBS_biovol_df_PFT.sizeClasses = pd.Categorical(NBS_biovol_df_PFT.sizeClasses, df_bins.sizeClasses.astype(str),ordered=True)  # this generates a categorical index for the column, this index can have a different lenght that the actual number of unique values, especially here since the categories come from df_bins
            NBS_biovol_df_PFT = pd.merge(NBS_biovol_df_PFT,NBS_biovol_df_PFT.drop_duplicates(subset=group, ignore_index=True)[ group].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left',on=group)
            multiindex = pd.MultiIndex.from_product([list( NBS_biovol_df_PFT.astype({column: 'category' for column in ['Group_index', 'sizeClasses']})[ column].cat.categories) for column in ['Group_index', 'sizeClasses']],names=['Group_index', 'sizeClasses'])
            NBS_biovol_df_PFT = pd.merge(NBS_biovol_df_PFT.drop_duplicates(['Group_index'])[group+['Group_index', 'Validation_percentage']], NBS_biovol_df_PFT.set_index(['Group_index', 'sizeClasses']).reindex(multiindex,fill_value=pd.NA).reset_index().drop( columns= group + [ 'Validation_percentage', 'size_class_mid','range_size_bin', 'ECD_mid', 'size_range_ECD']), how='right', on=['Group_index']).sort_values(['Group_index']).reset_index(drop=True).sort_values(['date_bin', 'Station_location', 'midLatBin', 'midLonBin','Min_obs_depth', 'Max_obs_depth','PFT','sizeClasses']).reset_index(drop=True)
            NBS_biovol_df_PFT = pd.merge(NBS_biovol_df_PFT.astype({'sizeClasses':str}), df_bins[['sizeClasses','size_class_mid', 'range_size_bin','ECD_mid', 'size_range_ECD']].astype({'sizeClasses':str}).drop_duplicates(), how='left', on='sizeClasses').astype({'size_class_mid':float})[['date_bin', 'Station_location', 'midLatBin', 'midLonBin','Min_obs_depth', 'Max_obs_depth', 'sizeClasses', 'size_class_mid', 'range_size_bin', 'ECD_mid', 'size_range_ECD', 'PFT','Validation_percentage', 'Biovolume_mean', 'size_class_pixel','ROI_number_sum', 'ROI_abundance_mean', 'NB', 'PSD','count_uncertainty', 'size_uncertainty', 'logNB', 'logSize', 'logPSD','logECD']].sort_values(['date_bin', 'Station_location', 'midLatBin', 'midLonBin','Min_obs_depth', 'Max_obs_depth','PFT','size_class_mid']).reset_index(drop=True)

            NBS_biovol_df_PFT = NBS_biovol_df_PFT.groupby(group).apply(lambda x: threshold_func(x)).reset_index(drop=True)
            ## Perform linear regression
            lin_fit_PFT = NBS_biovol_df_PFT.groupby( group+ ['Validation_percentage'] ).apply(lambda x: linear_fit_func(x)).reset_index().drop(columns='level_' + str(len(group + ['Validation_percentage']))).sort_values(['date_bin', 'Station_location'] + group).reset_index(drop=True).sort_values(group).reset_index(drop=True)
            NBSS_binned_all_PFT = pd.concat([NBSS_binned_all_PFT, NBS_biovol_df_PFT])
            lin_fit_data_PFT = pd.concat([lin_fit_data_PFT, lin_fit_PFT])


    NBSS_raw = NBSS_binned_all.filter(['date_bin', 'midLatBin', 'midLonBin','Taxon','Validation_percentage', 'size_class_mid', 'ECD_mid', 'NB', 'PSD','Min_obs_depth', 'Max_obs_depth'], axis=1)
    NBSS_raw_PFT = NBSS_binned_all_PFT.filter(['date_bin', 'midLatBin', 'midLonBin','PFT','Validation_percentage', 'size_class_mid', 'ECD_mid', 'NB', 'PSD','Min_obs_depth', 'Max_obs_depth'], axis=1)
    if len(NBSS_raw):
        # Average at bin levels
        NBSS_1a_class = NBSS_raw.groupby(['Taxon','Validation_percentage']).apply(lambda x: NBSS_stats_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_2'])
        NBSS_1a_class =  ocean_label_func(NBSS_1a_class.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')

        lin_fit_1b_class= lin_fit_data.groupby(['Taxon','Validation_percentage']).apply(lambda x: stats_linfit_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_2'])
        lin_fit_1b_class = ocean_label_func(lin_fit_1b_class.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')

        Path(NBSSpath / 'Class').mkdir(exist_ok=True, parents=True)
        NBSS_binned_all.to_csv(NBSSpath / 'Class' / str(instrument + '_Size-distribution_all_var_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
        NBSS_1a_class.to_csv(NBSSpath /'Class' / str(instrument + '_1a_Size-distribution_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
        lin_fit_1b_class.to_csv(NBSSpath / 'Class' / str(instrument + '_1b_Size-spectra-fit_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
    if len(NBSS_raw_PFT):
        NBSS_1a_PFT = NBSS_raw_PFT.groupby(['Taxon','Validation_percentage']).apply(lambda x: NBSS_stats_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_2'])
        NBSS_1a_PFT =  ocean_label_func(NBSS_1a_PFT.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')

        lin_fit_1b_PFT= lin_fit_data_PFT.groupby(['Taxon','Validation_percentage']).apply(lambda x: stats_linfit_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_2'])
        lin_fit_1b_PFT = ocean_label_func(lin_fit_1b_PFT.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')

        # Save data levels
        Path(NBSSpath / 'PFT').mkdir(exist_ok=True, parents=True)
        NBSS_binned_all_PFT.to_csv(NBSSpath / 'PFT' / str(instrument + '_Size-distribution_all_var_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
        NBSS_1a_PFT.to_csv(NBSSpath /'PFT' / str(instrument + '_1a_Size-distribution_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
        lin_fit_1b_PFT.to_csv(NBSSpath / 'PFT' / str(instrument + '_1b_Size-spectra-fit_v' + currentYear + '-' + currentMonth + '.csv'),index=False)