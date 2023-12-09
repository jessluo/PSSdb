# Objective: This step computes taxa-specific Normalized Biovolume Size Spectrum based on steps 0,1,2,3. Computation follows that of step 4, except that taxonomic groups are assigned to each ROI to compute taxa-specific products

# Python modules and path specifications:
import warnings
warnings.filterwarnings('ignore')


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


# Plot module
import plotly.express as px # Use pip install plotly==5.8.0
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotnine import * # Python equivalent to ggplot2. Use pip install plotnine. Do not execute in Pycharm (bug- no fix yet): https://youtrack.jetbrains.com/issue/PY-42390
from colorspace import sequential_hcl # Use: pip install git+https://github.com/retostauffer/python-colorspace

def report_product(biovolume_spectra,biomass_spectra,taxo_group='Taxon'):
    subplot_specs = [[{"type": "scattergeo"}, {"type": "scatter"}], [{"type": "scattergeo"}, {"type": "scatter"}]]
    subplot_titles = ['Biovolume integral', '','Biomass integral', '']
    fig_class = make_subplots(rows=2, cols=2,specs=subplot_specs, subplot_titles=subplot_titles, column_widths=[0.6, 0.4], row_heights=[0.5, 0.5], vertical_spacing=0.1)
    fig_class.update_yaxes(type="log", row=1, col=2)
    fig_class.update_yaxes(type="log", row=2, col=2)
    fig_class.update_xaxes(type="log", row=2, col=2, ticks="inside")
    fig_class.update_xaxes(type="log", row=1, col=2, ticks="inside")
    fig_class.for_each_xaxis(lambda x: x.update(showgrid=False))
    fig_class.for_each_yaxis(lambda x: x.update(showgrid=False, ticks="inside"))

    group=['year','month','latitude','longitude']
    biovolume_spectra['Group_NBSS'] = biovolume_spectra.groupby(['year', 'month', 'latitude', 'longitude',taxo_group], dropna=False).normalized_biovolume_mean.transform(sum)
    biovolume_spectra['Total_NBSS'] = biovolume_spectra.groupby(['year', 'month', 'latitude', 'longitude'],dropna=False).normalized_biovolume_mean.transform(sum)

    # Add spectrum
    biovolume_spectra = pd.merge(biovolume_spectra,biovolume_spectra.drop_duplicates(subset=group + [taxo_group, 'min_depth', 'max_depth'], ignore_index=True)[group + [taxo_group, 'min_depth', 'max_depth']].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left',on=group + [taxo_group, 'min_depth', 'max_depth']) if 'Group_index' not in biovolume_spectra.columns else biovolume_spectra
    nbss = biovolume_spectra.groupby('Group_index').apply(lambda x: x.append( x.tail(1).assign(equivalent_circular_diameter_mean=np.nan, normalized_biovolume_mean=np.nan))).reset_index( drop=True)

    fig_class.add_trace(go.Scatter(y=nbss.normalized_biovolume_mean, x=nbss.equivalent_circular_diameter_mean,line=dict(color='grey', width=0.15),name='All taxa (biovolume)',mode='lines', showlegend=True, visible=True), row=1, col=2)
    fig_class.add_trace(go.Scatter(y=nbss[nbss[taxo_group]==nbss[taxo_group].unique()[0]].normalized_biovolume_mean, x=nbss[nbss[taxo_group]==nbss[taxo_group].unique()[0]].equivalent_circular_diameter_mean, line=dict(color='black', width=0.5), name='taxa_biovolume', mode='lines', showlegend=True,visible=True), row=1, col=2)

    biomass_spectra['Group_NBSS'] = biomass_spectra.groupby(['year', 'month', 'latitude', 'longitude', taxo_group], dropna=False).normalized_biomass_mean.transform(sum)
    biomass_spectra['Total_NBSS'] = biomass_spectra.groupby(['year', 'month', 'latitude', 'longitude'], dropna=False).normalized_biomass_mean.transform( sum)

    # Add spectrum
    biomass_spectra= pd.merge(biomass_spectra, biomass_spectra.drop_duplicates( subset=group + [taxo_group, 'min_depth', 'max_depth'], ignore_index=True)[ group + [taxo_group, 'min_depth', 'max_depth']].reset_index().rename({'index': 'Group_index'}, axis='columns'),how='left', on=group + [taxo_group, 'min_depth','max_depth']) if 'Group_index' not in biomass_spectra.columns else biomass_spectra
    nbss_biomass =biomass_spectra.groupby('Group_index').apply( lambda x: x.append(x.tail(1).assign(biomass_mid=np.nan, normalized_biomass_mean=np.nan))).reset_index(drop=True)

    fig_class.add_trace(go.Scatter(y=nbss_biomass.normalized_biomass_mean, x=nbss_biomass.biomass_mid,
                                   line=dict(color='gray', width=0.15), name='All taxa (biomass)', mode='lines', showlegend=True,
                                   visible=True), row=2, col=2)
    fig_class.add_trace(go.Scatter(y=nbss_biomass[nbss_biomass[taxo_group]==nbss_biomass[taxo_group].unique()[0]].normalized_biomass_mean, x=nbss_biomass[nbss_biomass[taxo_group]==nbss_biomass[taxo_group].unique()[0]].biomass_mid,
                                   line=dict(color='black', width=0.5), name='taxa_biomass', mode='lines', showlegend=True,
                                   visible=True), row=2, col=2)

    # Button should be defined before geo traces to fix bug
    button_scatter1 = [dict(method="update",  # argument to change data
                            args=[{'lon': [biovolume_spectra[biovolume_spectra[taxo_group] == taxa].groupby( ['latitude', 'longitude']).Group_NBSS.mean().reset_index().longitude, biomass_spectra[biomass_spectra[taxo_group] == taxa].groupby( ['latitude', 'longitude']).Group_NBSS.mean().reset_index().longitude],
                                   'lat':  [biovolume_spectra[biovolume_spectra[taxo_group] == taxa].groupby( ['latitude', 'longitude']).Group_NBSS.mean().reset_index().latitude, biomass_spectra[biomass_spectra[taxo_group] == taxa].groupby( ['latitude', 'longitude']).Group_NBSS.mean().reset_index().latitude],
                                  'x': [nbss_biomass.biomass_mid,nbss_biomass[nbss_biomass[taxo_group]==taxa].biomass_mid,nbss.equivalent_circular_diameter_mean,nbss[nbss[taxo_group]==taxa].equivalent_circular_diameter_mean], #
                                  'y':[nbss_biomass.normalized_biomass_mean,nbss_biomass[nbss_biomass[taxo_group]==taxa].normalized_biomass_mean,nbss.normalized_biovolume_mean,nbss[nbss[taxo_group]==taxa].normalized_biovolume_mean], #
                                   'xaxis':['x2','x2','x','x'],'yaxis':['y2','y2','y','y'],
                                   'marker.color':[np.log10(biovolume_spectra[biovolume_spectra[taxo_group] == taxa].groupby( ['latitude', 'longitude']).Group_NBSS.mean().reset_index()['Group_NBSS']),np.log10(biomass_spectra[biomass_spectra[taxo_group] == taxa].groupby( ['latitude', 'longitude']).Group_NBSS.mean().reset_index()[ 'Group_NBSS'])]             ,
                                  'line.color':['gray','black', 'gray','black','gray','gray'],
                                   'name':['all taxa (biomass)','taxa_biomass','all taxa (biovolume)','taxa_biovolume','',''],
                                   'visible': [True]*len(fig_class.data)},{'traces': [0,1, 2,3,4,5]}],# visible should have the same length as number of subplots, second argument specifies which subplot is updated
                            label=taxa) for taxa in list(set(biovolume_spectra[taxo_group].unique()).intersection(set(biomass_spectra[taxo_group].unique())))]

    fig_class.update_layout(updatemenus=list([dict(active=0, buttons=button_scatter1, x=0.87, y=1.1, xanchor='left', yanchor='top')]),
        xaxis={'title': u'Equivalent circular diameter (\u03BCm)', 'tickfont': dict(size=10),'titlefont': dict(size=12)},
        xaxis2={'title': u'Wet weight (g)', 'tickfont': dict(size=10), 'titlefont': dict(size=12)},
        yaxis={"tickmode": "array",'title': u'Normalized Biovolume Size Spectrum <br> (\u03BCm\u00b3 dm\u207B\u00b3 \u03BCm\u207B\u00b3)', 'tickfont': dict(size=10), 'titlefont': dict(size=12)},
        yaxis2={"tickmode": "array",'title': u'Normalized Biomass Size Spectrum <br> (g dm\u207B\u00b3 g\u207B\u00b9)','tickfont': dict(size=10), 'titlefont': dict(size=12)})
    df_summary = biovolume_spectra.groupby(['latitude', 'longitude']).Total_NBSS.mean().reset_index()
    data_geo = dict(type='scattergeo',name='',
                    lon=df_summary.dropna(subset=['latitude', 'longitude']).longitude,
                    lat=df_summary.dropna(subset=['latitude', 'longitude']).latitude,
                    marker=dict(color=np.log10(df_summary.dropna(subset=['latitude', 'longitude'])['Total_NBSS']),colorbar=dict(titleside="top", outlinecolor="rgba(68, 68, 68, 0)", x=0.15,ticks="inside", orientation='h', title='$log_{10} particles per L$',len=0.3, xanchor='left'), colorscale=px.colors.sequential.Teal_r + px.colors.sequential.Reds, opacity=0.7,size=2 * (1 + np.log10(df_summary.dropna(subset=['latitude', 'longitude']).Total_NBSS) - min(np.log10(df_summary.dropna(subset=['latitude', 'longitude']).Total_NBSS)))),
                    geojson="natural earth", showlegend=False, geo='geo')
    layout = dict()
    layout['geo'] = dict(lonaxis_range=[-180, 180], lataxis_range=[-80, 80])
    # go.Figure(data=data_geo, layout=layout)
    fig_class.add_trace(data_geo, row=1, col=1).update_layout(go.Figure(data=data_geo, layout=layout).layout)
    df_summary = biomass_spectra.groupby(['latitude', 'longitude']).Total_NBSS.mean().reset_index()
    data_geo = dict(type='scattergeo',name='',lon=df_summary.dropna(subset=['latitude', 'longitude']).longitude, lat=df_summary.dropna(subset=['latitude', 'longitude']).latitude,
                    marker=dict(color=np.log10(df_summary.dropna(subset=['latitude', 'longitude'])['Total_NBSS']),colorbar=dict(titleside="top", outlinecolor="rgba(68, 68, 68, 0)", x=0.15, y=-0.15,ticks="inside", orientation='h', title='$log_{10} particles per L$',len=0.3, xanchor='left'),colorscale=px.colors.sequential.Teal_r + px.colors.sequential.Reds, opacity=0.7,size=2 * (1 + np.log10( df_summary.dropna(subset=['latitude', 'longitude']).Total_NBSS) - min(np.log10(df_summary.dropna(subset=['latitude', 'longitude']).Total_NBSS)))), geojson="natural earth", showlegend=False, geo='geo')
    layout = dict()
    layout['geo2'] = dict(lonaxis_range=[-180, 180], lataxis_range=[-80, 80])
    # go.Figure(data=data_geo, layout=layout)
    data_geo.update({'geo': 'geo2'})

    fig_class.add_trace(data_geo, row=2, col=1)
    fig_class.layout['geo2']={'domain': {'x': [0.1, 0.5], 'y': [0.0, 0.45]}}
    fig_class.layout['geo']['domain']['x']=[0.1, 0.5]


    return fig_class
fig=report_product(biovolume_spectra=NBSS_1a_class, biomass_spectra=NBSS_1a_class_biomass, taxo_group='PFT')
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import numpy as np
# Multiprocessing mapping
from multiprocessing.pool import ThreadPool # Use: pip install multiprocess
chunk=1000
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

from natsort import natsorted
import pint # Use pip install pint
warnings.filterwarnings('ignore', module='pint')
import pint_pandas # Use pip install pint-pandas
PA=pint_pandas.PintArray
from pint import UnitRegistry # Use pip install pint to install
ureg=UnitRegistry()
PQ=ureg.Quantity
#ureg.default_format = '~' # Add this if you want to use abbreviated unit names.
ureg.load_definitions(Path(cfg['git_dir']).expanduser()/cfg['units_file'])  # This text file is used to define custom units from standard units  (e.g. square_pixel etc.)
full_list_units = list(dict.fromkeys(sorted(dir(ureg))))  # list(dict.fromkeys(sorted(list(np.concatenate([dir(getattr(ureg.sys,system)) for system in dir(ureg.sys)]).flat))))

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
df_allometry=df_allometry[df_allometry.Size_proxy=='Biovolume']
# Append biomass bins
from scipy import stats
bins=np.power(2, np.arange(0, np.log2(100000) + 1 / 3, 1 / 3))  # Fixed size(um) bins used for UVP/EcoPart data. See https://ecopart.obs-vlfr.fr/. 1/3 ESD increments allow to bin particle of doubled biovolume in consecutive bins. np.exp(np.diff(np.log((1/6)*np.pi*(EcoPart_extended_bins**3))))
Ecopart_bins=pd.read_csv(Path(cfg['size_bins_path']).expanduser())
bins=Ecopart_bins.ESD_um.to_numpy()
biomass_bins=(1e-6*(df_allometry.loc[df_allometry.Taxon=='Living','C_Intercept'].values[0]*(1e-09*(1 / 6) * np.pi * bins ** 3)**(df_allometry.loc[df_allometry.Taxon=='Living','C_Slope'].values[0])))#(1e-12*(df_allometry.loc[df_allometry.Taxon=='Protozoa','C_Intercept'].values[0]*((1 / 6) * np.pi * bins ** 3)**(df_allometry.loc[df_allometry.Taxon=='Protozoa','C_Slope'].values[0])))
data_bins = pd.DataFrame({'sizeClasses': pd.cut(Ecopart_bins.biovol_um3.to_numpy(), Ecopart_bins.biovol_um3.to_numpy()).categories.values.astype(str),  # Define bin categories (cubic micrometers)
                        'sizeClasses_ECD': pd.cut(bins, bins).categories.values,  # Define bin categories (um)
                        'size_range_ECD': np.diff(bins),  # Define width of individual bin categories (um)
                        'range_size_bin': np.concatenate(np.diff((1 / 6) * np.pi * (np.resize(np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])), (len(bins) - 1, 2)) ** 3), axis=1)),  # cubic micrometers
                        'ECD_mid': stats.gmean(np.resize( np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])), (len(bins) - 1, 2)), axis=1), # Define geometrical mean of bin categories (um)
                        'size_class_mid':stats.gmean((1 / 6) * np.pi * (np.resize(np.append(bins[0], np.append(np.repeat(bins[1:-1], repeats=2), bins[len(bins) - 1])), (len(bins) - 1, 2)) ** 3), axis=1), # Define geometrical mean of bin categories (um)
                        'range_biomass_bin':np.concatenate(np.diff(np.resize(np.append(biomass_bins[0], np.append(np.repeat(biomass_bins[1:-1], repeats=2), biomass_bins[len(biomass_bins) - 1])), (len(biomass_bins) - 1, 2)), axis=1)), # in g
                        'biomass_mid':stats.gmean(np.resize( np.append(biomass_bins[0], np.append(np.repeat(biomass_bins[1:-1], repeats=2), biomass_bins[len(biomass_bins) - 1])), (len(biomass_bins) - 1, 2)), axis=1) #in g
})

# Convert all size units to cubic micrometers and mass units to gram
df_allometry['C_Intercept']=df_allometry['C_Intercept']*(df_allometry.Size_unit.apply(lambda x: PQ(1,'cubic_micrometer').to(x).magnitude)).values/(df_allometry.Elemental_mass_unit.apply(lambda x: PQ(1,'gram').to(x).magnitude)).values
df_allometry['Size_unit']='cubic_micrometer'
df_allometry['Elemental_mass_unit']='gram'

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
df_taxonomy=pd.read_excel(path_to_taxonomy)
df_taxonomy['Taxon']=df_taxonomy.Taxon.str.replace(" ","_")
#2: Loop through instrument-specific gridded files and compute class- and functional types-specifc Normalized Biovolume Size Spectra
for instrument in natsorted(os.listdir(Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir']))[1:2]:
    #first, break apart datasets by big global grids, to avoid making one HUGE file of all gridded datasets per instrument
    #get paths to gridded files
    file_list = glob(str(Path(cfg['raw_dir']).expanduser() / cfg['gridded_subdir']) + '*/**/' + instrument +'*_temp_binned_*.csv', recursive=True) #generate path and project ID's but ONLY for parsed data
    grid_list = [re.search('N_(.+?).csv', file) for file in file_list]#group_gridded_files_func(instrument, already_gridded='Y')
    currentMonth = str(datetime.datetime.now().month).rjust(2, '0')
    currentYear = str(datetime.datetime.now().year)
    NBSSpath = Path(cfg['raw_dir']).expanduser() / 'NBSS_data' / str('NBSS_ver_'+currentMonth+'_'+currentYear)


    biovol = 'Biovolume_area'
    print ('\ngenerating PSSdb biovolume/biomass products for ' +instrument+' based on ' + biovol)
    NBSS_binned_all,NBSS_binned_all_PFT,NBSS_binned_all_biomass,NBSS_binned_all_biomass_PFT = pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()  # NBSS dataset
    lin_fit_data,lin_fit_data_PFT,lin_fit_data_biomass,lin_fit_data_biomass_PFT = pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
    # Compute Normalized Biovolume Size Spectra
    with ThreadPool() as pool:
        for result in pool.map(lambda path: NB_SS_func( pd.merge(pd.merge(pd.read_csv(path),df_taxonomy[['Category', 'Taxon','PFT']], how='left', on='Category').assign(diameter=lambda x: ((x[biovol]*6)/m.pi)**(1./3.),sizeClasses=lambda x: pd.cut(x[biovol], bins=Ecopart_bins.biovol_um3.to_numpy()).astype(str)).rename(columns={'diameter':'diameter_from_{}'.format(biovol)}),data_bins,how='left',on='sizeClasses').assign(PFT=lambda x: np.where(x['PFT'] == 'Phytoplankton', pd.cut(x.ECD_mid,[0,2,20,200,20000],labels=['Pico','Nano','Micro','Meso']).astype(str)+x.PFT.astype(str).str.lower(), x.PFT)),data_bins, biovol_estimate=biovol, sensitivity=sensitivity,group=['Taxon','PFT'],thresholding=False)[0], file_list, chunksize=chunk):
            NBSS_binned_all = pd.concat([NBSS_binned_all, result], axis=0)
    if len(NBSS_binned_all):
        NBSS_binned_all =NBSS_binned_all.reset_index(drop=True)
        group = ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon']
        ## Add empty size classes, threshold, and compute linear fit based on Taxon only
        NBSS_binned_all.sizeClasses = pd.Categorical( NBSS_binned_all.sizeClasses, data_bins.sizeClasses.astype(str), ordered=True)  # this generates a categorical index for the column, this index can have a different lenght that the actual number of unique values, especially here since the categories come from df_bins
        NBSS_binned_all = pd.merge(NBSS_binned_all, NBSS_binned_all.drop_duplicates(subset=group, ignore_index=True)[ group].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left', on=group)
        multiindex = pd.MultiIndex.from_product([list( NBSS_binned_all.astype({column: 'category' for column in ['Group_index', 'sizeClasses']})[ column].cat.categories) for column in ['Group_index', 'sizeClasses']], names=['Group_index', 'sizeClasses'])
        df_test=NBSS_binned_all.set_index(['Group_index', 'sizeClasses'])[NBSS_binned_all.set_index(['Group_index', 'sizeClasses']).index.duplicated(keep=False)]

        NBSS_binned_all_thres = pd.merge( NBSS_binned_all.drop_duplicates(['Group_index'])[group + ['Group_index', 'Validation_percentage']],NBSS_binned_all.set_index(['Group_index', 'sizeClasses']).reindex(multiindex, fill_value=pd.NA).reset_index().drop(columns=group + ['Validation_percentage', 'size_class_mid', 'range_size_bin', 'ECD_mid', 'size_range_ECD']), how='right', on=['Group_index']).sort_values( ['Group_index']).reset_index(drop=True).sort_values( ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon', 'sizeClasses']).reset_index(drop=True)
        NBSS_binned_all_thres =pd.merge(NBSS_binned_all_thres.astype({'sizeClasses': str}), data_bins[['sizeClasses', 'size_class_mid', 'range_size_bin', 'ECD_mid', 'size_range_ECD']].astype({'sizeClasses': str}), how='left', on='sizeClasses').astype({'size_class_mid': float}).sort_values( ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon', 'size_class_mid']).reset_index(drop=True)
        NBSS_binned_all=NBSS_binned_all_thres.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon', 'Validation_percentage']).apply(lambda x: threshold_func(x)).reset_index(drop=True).sort_values( ['date_bin', 'Station_location','midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon','size_class_mid']).reset_index(drop=True)
        lin_fit = NBSS_binned_all.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon', 'Validation_percentage']).apply(lambda x: linear_fit_func(x)).reset_index().drop(columns='level_' + str( len(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon'] + [ 'Validation_percentage']))).sort_values( ['Taxon','date_bin', 'Station_location','midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth']).reset_index(drop=True)

        ## Sum Taxon-specific NBSS to obtain PFT-specifc NBSS
        NBSS_binned_all_thres_PFT = NBSS_binned_all_thres.dropna(subset=['NB']).groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth','sizeClasses', 'size_class_mid', 'range_size_bin', 'ECD_mid', 'size_range_ECD', 'PFT']).apply(lambda x: pd.Series({ 'Validation_percentage': np.round(np.nansum((x.Validation_percentage.astype(float) * x.ROI_number_sum.astype( float)) / x.ROI_number_sum.astype(float).sum()), 2),'Biovolume_mean': np.nansum((x.Biovolume_mean.astype(float) * x.ROI_number_sum) / np.nansum(x.ROI_number_sum)),'size_class_pixel': np.nanmean(x.size_class_pixel.astype(float)), 'ROI_number_sum': np.nansum(x.ROI_number_sum.astype( float)),'ROI_abundance_mean': np.nansum(( x.ROI_abundance_mean.astype(float) * x.ROI_number_sum.astype(float)) / np.nansum(x.ROI_number_sum.astype(float))),'NB': np.nansum(x.NB.astype( float)), 'PSD': np.nansum(x.PSD.astype(float)), 'count_uncertainty': poisson.pmf( k=np.nansum(x.ROI_number_sum.astype( float)), mu=np.nansum(x.ROI_number_sum.astype( float))),'size_uncertainty': np.nanmean( x.size_uncertainty.astype(float)), 'logNB': np.log10( np.nansum(x.NB.astype( float))), 'logSize': np.log10( x.size_class_mid.astype( float).unique()[0]),'logPSD': np.log10( np.nansum(x.PSD.astype( float))), 'logECD': np.log10(x.ECD_mid.astype( float).unique()[0])})).reset_index().sort_values( ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT', 'size_class_mid']).reset_index(drop=True)
        NBSS_binned_all_PFT=NBSS_binned_all_thres_PFT.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT', 'Validation_percentage']).apply(lambda x: threshold_func(x)).reset_index(drop=True).sort_values( ['date_bin', 'Station_location','midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT','size_class_mid']).reset_index(drop=True)
        lin_fit_PFT = NBSS_binned_all_PFT.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT', 'Validation_percentage']).apply(lambda x: linear_fit_func(x)).reset_index().drop(columns='level_' + str( len(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT'] + [ 'Validation_percentage']))).sort_values( ['PFT','date_bin', 'Station_location','midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth']).reset_index(drop=True)
        ## Save subbin datafiles
        NBSS_raw = NBSS_binned_all.filter(['date_bin', 'midLatBin', 'midLonBin','Taxon','PFT','Validation_percentage', 'size_class_mid', 'ECD_mid', 'NB', 'PSD','Min_obs_depth', 'Max_obs_depth'], axis=1)
        NBSS_raw_PFT = NBSS_binned_all_PFT.filter(['date_bin', 'midLatBin', 'midLonBin','PFT','Validation_percentage', 'size_class_mid', 'ECD_mid', 'NB', 'PSD','Min_obs_depth', 'Max_obs_depth'], axis=1)
        Path(NBSSpath / 'Class').mkdir(exist_ok=True, parents=True)
        NBSS_binned_all.to_csv(NBSSpath / 'Class' / str(instrument + '_Size-distribution_all_var_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
        Path(NBSSpath / 'PFT').mkdir(exist_ok=True, parents=True)
        NBSS_binned_all_PFT.to_csv(NBSSpath / 'PFT' / str(instrument + '_Size-distribution_all_var_v' + currentYear + '-' + currentMonth + '.csv'),index=False)

    # Compute Normalized Biomass Size Spectra
    with ThreadPool() as pool:
        for result in pool.map(lambda path: NB_SS_func( pd.merge(pd.merge(pd.merge(pd.read_csv(path),df_taxonomy[['Category', 'Taxon','PFT']], how='left', on='Category').assign(diameter=lambda x: ((x[biovol]*6)/m.pi)**(1./3.),sizeClasses=lambda x: pd.cut(x[biovol], bins=Ecopart_bins.biovol_um3.to_numpy()).astype(str)).rename(columns={'diameter':'diameter_from_{}'.format(biovol)}),data_bins,how='left',on='sizeClasses'), df_allometry.query('Size_proxy=="Biovolume"')[['Taxon','Size_unit','Elemental_mass_unit','C_Intercept','C_Slope']],how='left',on=['Taxon']).assign(Biomass_area=lambda x: x['C_Intercept']*(x[biovol]**x['C_Slope']),PFT=lambda x: np.where(x['PFT'] == 'Phytoplankton', pd.cut(x.ECD_mid,[0,2,20,200,20000],labels=['Pico','Nano','Micro','Meso']).astype(str)+x.PFT.astype(str).str.lower(), x.PFT)).drop(columns=['range_size_bin']).assign(range_size_bin=lambda x: x.range_biomass_bin),data_bins.drop(columns=['size_class_mid']).assign(size_class_mid=data_bins.biomass_mid/data_bins.biomass_mid[0]), biovol_estimate='Biomass_area', sensitivity=sensitivity,group=['Taxon','PFT'],thresholding=False)[0], file_list, chunksize=chunk):
            NBSS_binned_all_biomass = pd.concat([NBSS_binned_all_biomass, result], axis=0)
    if len(NBSS_binned_all_biomass):
        NBSS_binned_all_biomass =NBSS_binned_all_biomass.rename(columns={'Biovolume_mean':'Biomass_mean'}).reset_index(drop=True)
        NBSS_binned_all_biomass = pd.merge(NBSS_binned_all_biomass.astype({'sizeClasses':str}),data_bins.astype({'sizeClasses':str})[['sizeClasses', 'range_biomass_bin', 'biomass_mid']], how='left', on=['sizeClasses'])

        group = ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon']
        ## Add empty size classes, threshold, and compute linear fit based on Taxon only
        NBSS_binned_all_biomass.sizeClasses = pd.Categorical( NBSS_binned_all_biomass.sizeClasses, data_bins.sizeClasses.astype(str), ordered=True)  # this generates a categorical index for the column, this index can have a different lenght that the actual number of unique values, especially here since the categories come from df_bins
        NBSS_binned_all_biomass = pd.merge(NBSS_binned_all_biomass,   NBSS_binned_all_biomass.drop_duplicates(subset=group, ignore_index=True)[ group].reset_index().rename({'index': 'Group_index'}, axis='columns'), how='left', on=group)
        multiindex = pd.MultiIndex.from_product([list( NBSS_binned_all_biomass.astype({column: 'category' for column in ['Group_index', 'sizeClasses']})[ column].cat.categories) for column in ['Group_index', 'sizeClasses']], names=['Group_index', 'sizeClasses'])
        NBSS_binned_all_biomass_thres = pd.merge( NBSS_binned_all_biomass.drop_duplicates(['Group_index'])[group + ['Group_index', 'Validation_percentage']],NBSS_binned_all_biomass.set_index(['Group_index', 'sizeClasses']).reindex(multiindex, fill_value=pd.NA).reset_index().drop(columns=group + ['Validation_percentage', 'size_class_mid', 'range_size_bin', 'ECD_mid', 'size_range_ECD']), how='right', on=['Group_index']).sort_values( ['Group_index']).reset_index(drop=True).sort_values( ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon', 'sizeClasses']).reset_index(drop=True)
        NBSS_binned_all_biomass_thres =pd.merge(NBSS_binned_all_biomass_thres.astype({'sizeClasses': str}), data_bins[['sizeClasses', 'size_class_mid', 'range_size_bin', 'ECD_mid', 'size_range_ECD']].astype({'sizeClasses': str}), how='left', on='sizeClasses').astype({'size_class_mid': float}).sort_values( ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon', 'size_class_mid']).reset_index(drop=True)
        NBSS_binned_all_biomass=NBSS_binned_all_biomass_thres.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon', 'Validation_percentage']).apply(lambda x: threshold_func(x)).reset_index(drop=True).sort_values( ['date_bin', 'Station_location','midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon','size_class_mid']).reset_index(drop=True)
        NBSS_binned_all_biomass['Total_biomass'] = NBSS_binned_all_biomass.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon', 'Validation_percentage'], dropna=False).apply(lambda x: pd.DataFrame({'Total_biomass': 1e+06 * np.nansum(x.NB * x.range_biomass_bin)}, index=list( x.index))).reset_index().Total_biomass.values  # in microgram per liters or milligram per cubic meters
        lin_fit_biomass = NBSS_binned_all_biomass.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon', 'Validation_percentage']).apply(lambda x: linear_fit_func(x.drop(columns=['logSize']).assign(logSize=np.log10(x.biomass_mid / data_bins.biomass_mid[0])))).reset_index().drop(columns='level_' + str( len(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'Taxon'] + [ 'Validation_percentage']))).sort_values( ['Taxon','date_bin', 'Station_location','midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth']).reset_index(drop=True)

        ## Sum Taxon-specific NBSS to obtain PFT-specifc NBSS
        NBSS_binned_all_thres_PFT_biomass = NBSS_binned_all_biomass_thres.dropna(subset=['NB']).groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth','sizeClasses', 'size_class_mid', 'range_size_bin', 'ECD_mid', 'size_range_ECD', 'biomass_mid', 'range_biomass_bin', 'PFT']).apply(lambda x: pd.Series({ 'Validation_percentage': np.round(np.nansum((x.Validation_percentage.astype(float) * x.ROI_number_sum.astype( float)) / x.ROI_number_sum.astype(float).sum()), 2),'Biomass_mean': np.nansum((x.Biomass_mean.astype(float) * x.ROI_number_sum) / np.nansum(x.ROI_number_sum)),'size_class_pixel': np.nanmean(x.size_class_pixel.astype(float)), 'ROI_number_sum': np.nansum(x.ROI_number_sum.astype( float)),'ROI_abundance_mean': np.nansum(( x.ROI_abundance_mean.astype(float) * x.ROI_number_sum.astype(float)) / np.nansum(x.ROI_number_sum.astype(float))),'NB': np.nansum(x.NB.astype( float)), 'PSD': np.nansum(x.PSD.astype(float)), 'count_uncertainty': poisson.pmf( k=np.nansum(x.ROI_number_sum.astype( float)), mu=np.nansum(x.ROI_number_sum.astype( float))),'size_uncertainty': np.nanmean( x.size_uncertainty.astype(float)), 'logNB': np.log10( np.nansum(x.NB.astype( float))), 'logSize': np.log10( x.size_class_mid.astype( float).unique()[0]),'logPSD': np.log10( np.nansum(x.PSD.astype( float))), 'logECD': np.log10(x.ECD_mid.astype( float).unique()[0])})).reset_index().sort_values( ['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT', 'size_class_mid']).reset_index(drop=True)
        NBSS_binned_all_PFT_biomass=NBSS_binned_all_thres_PFT_biomass.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT', 'Validation_percentage']).apply(lambda x: threshold_func(x)).reset_index(drop=True).sort_values( ['date_bin', 'Station_location','midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT','size_class_mid']).reset_index(drop=True)
        NBSS_binned_all_PFT_biomass['Total_biomass'] = NBSS_binned_all_PFT_biomass.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT', 'Validation_percentage'], dropna=False).apply(lambda x: pd.DataFrame({'Total_biomass': 1e+06 * np.nansum(x.NB * x.range_biomass_bin)}, index=list( x.index))).reset_index().Total_biomass.values  # in microgram per liters or milligram per cubic meters
        lin_fit_PFT_biomass = NBSS_binned_all_PFT_biomass.groupby(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT', 'Validation_percentage']).apply(lambda x: linear_fit_func(x.drop(columns=['logSize']).assign(logSize=np.log10(x.biomass_mid / data_bins.biomass_mid[0])))).reset_index().drop(columns='level_' + str( len(['date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth', 'PFT'] + [ 'Validation_percentage']))).sort_values( ['PFT','date_bin', 'Station_location','midLatBin', 'midLonBin', 'Min_obs_depth', 'Max_obs_depth']).reset_index(drop=True)
        ## Save subbin datafiles
        NBSS_raw_biomass = NBSS_binned_all_biomass.filter(['date_bin', 'midLatBin', 'midLonBin','Taxon','PFT','Validation_percentage', 'size_class_mid', 'ECD_mid','biomass_mid', 'NB', 'PSD','Min_obs_depth', 'Max_obs_depth','Total_biomass'], axis=1)
        NBSS_raw_biomass_PFT = NBSS_binned_all_PFT_biomass.filter(['date_bin', 'midLatBin', 'midLonBin','PFT','Validation_percentage', 'size_class_mid', 'ECD_mid','biomass_mid', 'NB', 'PSD','Min_obs_depth', 'Max_obs_depth','Total_biomass'], axis=1)
        Path(NBSSpath / 'Class').mkdir(exist_ok=True, parents=True)
        NBSS_binned_all_biomass.to_csv(NBSSpath / 'Class' / str(instrument + '_Biomass-distribution_all_var_v' + currentYear + '-' + currentMonth + '.csv'), index=False)
        Path(NBSSpath / 'PFT').mkdir(exist_ok=True, parents=True)
        NBSS_binned_all_PFT_biomass.to_csv(NBSSpath / 'PFT' / str(instrument + '_Biomass-distribution_all_var_v' + currentYear + '-' + currentMonth + '.csv'),index=False)

    # Saving final dataframes and interactive report
    if len(NBSS_raw):
        # Average at bin levels
        NBSS_1a_class = NBSS_raw.groupby(['Taxon','Validation_percentage']).apply(lambda x: NBSS_stats_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_2'])
        NBSS_1a_class =  ocean_label_func(NBSS_1a_class.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')
        NBSS_1a_class =NBSS_1a_class.dropna(subset=['Taxon']).sort_values(['year', 'month', 'latitude', 'longitude', 'Taxon', 'equivalent_circular_diameter_mean']).reset_index(drop=True)

        lin_fit_1b_class= lin_fit.groupby(['Taxon','Validation_percentage']).apply(lambda x: stats_linfit_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_2'])
        lin_fit_1b_class = ocean_label_func(lin_fit_1b_class.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')
        lin_fit_1b_class =lin_fit_1b_class.sort_values(['Taxon','year','month','latitude','longitude']).reset_index(drop=True)

        NBSS_1a_class.to_csv(NBSSpath /'Class' / str(instrument + '_1a_Size-distribution_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
        lin_fit_1b_class.to_csv(NBSSpath / 'Class' / str(instrument + '_1b_Size-spectra-fit_v' + currentYear + '-' + currentMonth + '.csv'),index=False)

    if len(NBSS_raw_PFT):
        NBSS_1a_PFT = NBSS_raw_PFT.groupby(['PFT','Validation_percentage']).apply(lambda x: NBSS_stats_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_2'])
        NBSS_1a_PFT =  ocean_label_func(NBSS_1a_PFT.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')
        NBSS_1a_PFT = NBSS_1a_PFT.dropna(subset=['PFT']).sort_values(['year', 'month', 'latitude', 'longitude', 'PFT', 'equivalent_circular_diameter_mean']).reset_index( drop=True)

        lin_fit_1b_PFT= lin_fit_PFT.groupby(['PFT','Validation_percentage']).apply(lambda x: stats_linfit_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_2'])
        lin_fit_1b_PFT = ocean_label_func(lin_fit_1b_PFT.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')
        lin_fit_1b_PFT =lin_fit_1b_PFT.dropna(subset=['PFT']).sort_values(['PFT', 'year', 'month', 'latitude', 'longitude']).reset_index(drop=True)

        NBSS_1a_PFT.to_csv(NBSSpath /'PFT' / str(instrument + '_1a_Size-distribution_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
        lin_fit_1b_PFT.to_csv(NBSSpath / 'PFT' / str(instrument + '_1b_Size-spectra-fit_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
    if len(NBSS_raw_biomass):
        # Average at bin levels
        NBSS_1a_class_biomass = NBSS_raw_biomass.groupby(['Taxon','Validation_percentage','biomass_mid']).apply(lambda x: NBSS_stats_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_3']).rename(columns={'normalized_biovolume_mean':'normalized_biomass_mean','normalized_biovolume_sd':'normalized_biomass_sd'})
        NBSS_1a_class_biomass =  ocean_label_func(NBSS_1a_class_biomass.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')
        NBSS_1a_class_biomass = NBSS_1a_class_biomass.dropna(subset=['Taxon']).sort_values( ['year', 'month', 'latitude', 'longitude', 'Taxon', 'equivalent_circular_diameter_mean']).reset_index(drop=True)

        lin_fit_1b_class_biomass= lin_fit_biomass.groupby(['Taxon','Validation_percentage']).apply(lambda x: stats_linfit_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_2'])
        lin_fit_1b_class_biomass = ocean_label_func(lin_fit_1b_class_biomass.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')
        lin_fit_1b_class_biomass =lin_fit_1b_class_biomass.dropna(subset=['Taxon']).sort_values(['Taxon', 'year', 'month', 'latitude', 'longitude']).reset_index(drop=True)

        NBSS_1a_class_biomass.to_csv(NBSSpath /'Class' / str(instrument + '_1a_Biomass-distribution_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
        lin_fit_1b_class_biomass.to_csv(NBSSpath / 'Class' / str(instrument + '_1b_Biomass-spectra-fit_v' + currentYear + '-' + currentMonth + '.csv'),index=False)



    if len(NBSS_raw_biomass_PFT):
        NBSS_1a_biomass_PFT = NBSS_raw_biomass_PFT.groupby(['PFT','Validation_percentage','biomass_mid']).apply(lambda x: NBSS_stats_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_3']).rename(columns={'normalized_biovolume_mean':'normalized_biomass_mean','normalized_biovolume_sd':'normalized_biomass_sd'})
        NBSS_1a_biomass_PFT =  ocean_label_func(NBSS_1a_biomass_PFT.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')
        NBSS_1a_biomass_PFT=NBSS_1a_biomass_PFT.sort_values(['year', 'month', 'latitude', 'longitude', 'PFT', 'equivalent_circular_diameter_mean']).reset_index(drop=True)

        lin_fit_1b_biomass_PFT= lin_fit_PFT_biomass.groupby(['PFT','Validation_percentage']).apply(lambda x: stats_linfit_func(x, bin_loc=bin_loc, group_by=group_by)).reset_index().drop(columns=['level_2'])
        lin_fit_1b_biomass_PFT = ocean_label_func(lin_fit_1b_biomass_PFT.assign(month_int=lambda x: x.astype({'month':int}).month,year_int=lambda x: x.year.astype(int)).sort_values(by=['year_int', 'month_int']).drop(['year_int', 'month_int'], axis=1), 'longitude', 'latitude')
        lin_fit_1b_biomass_PFT =lin_fit_1b_biomass_PFT.sort_values(['PFT', 'year', 'month', 'latitude', 'longitude']).reset_index(drop=True)

        # Save data levels
        NBSS_1a_biomass_PFT.to_csv(NBSSpath /'PFT' / str(instrument + '_1a_Biomass-distribution_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
        lin_fit_1b_biomass_PFT.to_csv(NBSSpath / 'PFT' / str(instrument + '_1b_Biomass-spectra-fit_v' + currentYear + '-' + currentMonth + '.csv'),index=False)
