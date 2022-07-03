##  Objective: This file contains one function to:
# (1) print the fields included in Ecotaxa projects using the API (download is not necessary). This should facilitate the completion of the standardizer spreadsheets
# (1bis) print the list of standard units from pint module
# (2) save interactive plot of raw ecotaxa export file to identify outliers, missing values, and units. This should facilitate the completion of the standardizer spreadsheets
# (3) perform EcoTaxa export files standardization (standard column names and units) and harmonization (across instruments) based on standardizer spreadsheets
#  Running function (3) is required after project export in order to reduce (i.e. select only the raw columns necessary for further processing), harmonize (all the PSSdb projects now have similar columns, regardless of specific instrumentation and image acquisition/processing routines) standardize (All columns follow standardized notations and units according to international database) Ecotaxa project export files

## Requirements: Run ecotaxa_API_project_list to generate templates of standardizer spreadsheets

## Useful documentations:
# ZooScan / Zooprocess manual : https://sites.google.com/view/piqv/piqv-manuals/ecotaxa-manuals?authuser=0

## TO DO: Perform standardization on UVP and ISIIS projects (at the moment function 3 only works for Zooscan and IFCB)


## Python modules
# Path modules
from pathlib import Path # Handling of path object

# Config modules
import yaml # requires installation of PyYAML package

# Prompt for export confirmation
import sys

# Panda and numpy modules (used for dataframes and arays)
import pandas as pd
import numpy as np

# Plot module
import plotly.express as px # Use pip install plotly==5.8.0
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotnine import * # Python equivalent to ggplot2. Use pip install plotnine. Do not execute in Pycharm (bug- no fix yet): https://youtrack.jetbrains.com/issue/PY-42390
from colorspace import sequential_hcl # Use: pip install git+https://github.com/retostauffer/python-colorspace

# Ecotaxa API modules. Install module via git first: pip install git+https://github.com/ecotaxa/ecotaxa_py_client.git
import ecotaxa_py_client
from ecotaxa_py_client.api import projects_api
from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.model.login_req import LoginReq
from ecotaxa_py_client.api import objects_api
from ecotaxa_py_client.model.project_filters import ProjectFilters

# Unit conversion
import pint # Use pip install pint
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings('ignore', module='pint')
import pint_pandas # Use pip install pint-pandas
PA=pint_pandas.PintArray
PQ=pint.quantity.Quantity
from pint import UnitRegistry # Use pip install pint to install
ureg=UnitRegistry()
#ureg.default_format = '~' # Add this if you want to use abbreviated unit names.
ureg.load_definitions('/Users/mc4214/GIT/PSSdb/scripts/units_def.txt')  # This text file is used to define custom units from standard units  (e.g. square_pixel etc.)
full_list_units = list(dict.fromkeys(sorted(dir(ureg))))  # list(dict.fromkeys(sorted(list(np.concatenate([dir(getattr(ureg.sys,system)) for system in dir(ureg.sys)]).flat))))
import re


# Functions here:

# Function (1): Use EcoTaxa API to search for project fields automatically without downloading export files
def filling_standardizer_field_func(config_path,standardizer_path,**project_id):
    """
    Objective: This function uses EcoTaxa API to print the free fields in accessible projects.
    This should facilitate the completion of projects standardizer spreadsheets.
    :param config_path: Full path of the configuration file containing Ecotaxa authentication info
    :param standardizer_path: Full path of the standardizer spreadsheet containing project ID of interest
    :param (optional) project_id: integer for specific project ID
    :return: free fields for all accessible projects or specific project ID
    """
    kwargs = project_id
    project_id = kwargs.get('project_id')

    # Open configuration file to get: EcoTaxa authentication info
    path_to_config_usr = Path(config_path).expanduser()
    with open(path_to_config_usr, 'r') as config_file:
        cfg_pw = yaml.safe_load(config_file)
    # Search visible and accessible projects based on instrument filter (Authentication in EcoTaxa required)
    with ecotaxa_py_client.ApiClient() as client:
        api = authentification_api.AuthentificationApi(client)
        token = api.login(LoginReq(
            username=cfg_pw['ecotaxa_user'],
            password=cfg_pw['ecotaxa_pass']
        ))

    configuration = ecotaxa_py_client.Configuration(host="https://ecotaxa.obs-vlfr.fr/api", access_token=token,discard_unknown_keys=True)
    with ecotaxa_py_client.ApiClient(configuration) as api_client:
        api_instance = projects_api.ProjectsApi(api_client)
        try:
            # Search Projects and return fields
           api_response_project_search_accessible = api_instance.search_projects(also_others=False,# Return accessible projects given account rights
                                                                                  title_filter='',# Optional title project filter
                                                                                  instrument_filter='',# All instruments
                                                                                  order_field='projid'# Sorting variable. Use instrument or projid
                                                                                  )

        except ecotaxa_py_client.ApiException as e:
            print("Exception when calling ProjectsApi->search_projects: %s\n" % e)
    # Step 2: Generate a dictionary for individual project with variables of interest (i.e. free fields)
    subset_df_instrument=pd.read_excel(standardizer_path, usecols=["Project_ID"])
    project_fields=dict(map(lambda x: (x.projid,["acq_" + str for str in list(x.acquisition_free_cols.keys())] +
                                      ["object_" + str for str in list(x.obj_free_cols.keys())] +
                                      ["sample_" + str for str in list(x.sample_free_cols.keys())] +
                                      ["process_" + str for str in list(x.process_free_cols.keys())]),
                     list(filter(lambda w: any(subset_df_instrument['Project_ID'].isin([w.projid])),api_response_project_search_accessible))))

    if not project_id is None:
        if int(project_id) in list(project_fields.keys()):
            project_fields = {int(project_id): project_fields[int(project_id)]}


    print(yaml.dump(project_fields))
    return project_fields

# Function (1b): Print possible units
def filling_standardizer_unit_func(custom_units_path='~/GIT/PSSdb/scripts/units_def.txt'):
    """
    Objective: This function print the units available to fill out the standardizer spreadsheets.
    :param custom_units_path: Full path of the custom unit definition file
    :return: units
    """
    ureg = UnitRegistry()
    # ureg.default_format = '~' # Add this if you want to use abbreviated unit names.
    ureg.load_definitions(Path(custom_units_path).expanduser())  # This text file is used to define custom units from standard units  (e.g. square_pixel etc.)
    full_list_units = list(dict.fromkeys(sorted( dir(ureg))))  # list(dict.fromkeys(sorted(list(np.concatenate([dir(getattr(ureg.sys,system)) for system in dir(ureg.sys)]).flat))))
    print('Use unit from the list below or define custom units in {}:\n'.format(custom_units_path),full_list_units,sep='')


# Function (2): Use EcoTaxa API to plot variable whose units need to be checked automatically without downloading export files
def diagnostic_plots_func(config_path,standardizer_path,project_id):
    """
       Objective: This function uses standardizer spreadsheets to plot variable whose units need to be checked interactively for a unique project.
       This should facilitate (1) datapoints flagging and filtering (2) the completion of projects standardizer spreadsheets.
       :param config_path: Full path of the configuration file containing path of Ecotaxa export files
       :param standardizer_path: Full path of the standardizer spreadsheet containing project ID of interest
       :param project_id: unique integer for project ID to be plotted
       :return: print interactive plots
       """

    # Open configuration file to get: EcoTaxa authentication info
    path_to_config_usr = Path(config_path).expanduser()
    with open(path_to_config_usr, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    path_to_git = Path(cfg['git_dir']).expanduser()

    df_instrument = pd.read_excel(standardizer_path,index_col=0)
    path_to_plotdir = path_to_git / 'figures'
    path_to_plotdir.mkdir(parents=True, exist_ok=True)

    path_to_data = path_to_git / Path(cfg['dataset_subdir']) / df_instrument.at[project_id,"Instrument"]
    columns_of_interest =   ["Sample_ID","Cruise_field"]+[string.replace("_unit", "_field") for string in list(df_instrument.columns) if "_unit" in string]
    fields_of_interest =dict(zip(columns_of_interest,["sample_id"]+df_instrument.loc[project_id][columns_of_interest[1:]].values.flatten().tolist()))
    fields_of_interest = dict((key.replace("_field",""),value.split(",")[0]) for key,value in fields_of_interest.items() if isinstance(value, str))

    df=pd.read_table(list(path_to_data.glob('ecotaxa_export_{}*'.format(str(project_id))))[0],usecols=fields_of_interest.values())
    old_columns=df.columns
    df.columns = [(list(fields_of_interest.keys())[list(fields_of_interest.values()).index(value)]) for value in old_columns] #pd.MultiIndex.from_arrays([[list(fields_of_interest.keys())[list(fields_of_interest.values()).index(value)] for value in df.columns],df.columns])
    spatial_bins=1 # Fixed spatial bining 1x1 degrees
    vertical_bins=5 # Fixed vertical bin 50m
    df=df.assign(Latitude_group=pd.cut(df.Latitude,np.arange(-90,90,step=spatial_bins)),
                 Longitude_group=pd.cut(df.Longitude,np.arange(-180,180,step=spatial_bins)),
                 Depth_med=np.add(df.Depth_min,df.Depth_max)/2,
                 Depth_group=pd.cut(np.add(df.Depth_min,df.Depth_max)/2,np.arange(0,6000,step=vertical_bins),labels=np.arange(0,6000,step=vertical_bins)[:-1]))
    df=pd.merge(df,df.groupby(['Sample_ID'],dropna=False).agg(objects=(list(fields_of_interest.keys())[9],'count')).reset_index(),on=['Sample_ID'])
    sample_df=df.groupby(['Latitude_group','Longitude_group','Depth_med'],observed=True,dropna=False).apply(lambda x:pd.Series({'n_objects':np.nanmean(x.objects),
                                                                                                                   'n_objects_upper':np.nanmean(x.objects)+np.nanstd(x.objects),
                                                                                                                   'n_objects_lower': np.nanmean(x.objects)-np.nanstd(x.objects),
                                                                                                                   'std_objects':np.nanstd(x.objects)
                                                                                                       })).reset_index()
    profiles_df=df.groupby(['Latitude','Longitude','Sample_ID','Depth_group'],observed=True,dropna=False).apply(lambda x:pd.Series({'n_objects':np.nanmean(x.objects)})).reset_index()
    df=pd.merge(df,sample_df,on=['Latitude_group','Longitude_group','Depth_med'])

    plot=(ggplot(sample_df) +
                geom_point(mapping=aes(x='n_objects', y='std_objects'),color='#222222') +
                geom_abline(mapping=aes(intercept=0,slope=0.5))+
                labs(title=r'Count uncertainties in spatial bins ({}x{}$^\circ$)'.format(str(spatial_bins),str(spatial_bins)), x='Number of objects(images) per sample in spatial bins',y='Std of the number of objects per sample in spatial bins') +
                scale_y_log10(labels=lambda y: list(map(lambda x: (("{0:.{1:d}e}".format(x, 0)).split("e"))[0] + "x10$^{" +(("{0:.{1:d}e}".format(x, 2)).split("e"))[1] + "}$", y))) +
                scale_x_log10(labels=lambda y: list(map(lambda x: (("{0:.{1:d}e}".format(x, 0)).split("e"))[0] + "x10$^{" +(("{0:.{1:d}e}".format(x, 2)).split("e"))[1] + "}$", y))) +
                theme(panel_grid=element_blank(), panel_background=element_rect(fill='white'),
                      panel_border=element_rect(color='#222222'), axis_title=element_text(size=10),
                      axis_text_x=element_text(size=10), axis_text_y=element_text(size=10, rotation=90),
                      plot_background=element_blank()))
    #print(plot)
    subset_df=df.drop_duplicates(subset=list(fields_of_interest.keys())[0:5])
    subset_df=subset_df.assign(Area=subset_df.Area if 'Area' in subset_df.columns else None,
                     Biovolume=subset_df.Biovolume if 'Biovolume' in subset_df.columns else None,
                     Minor_axis=subset_df.Minor_axis if 'Minor_axis' in subset_df.columns else None)

    subset_df =subset_df.sort_values('Sample_ID')
    #subset_df = df.groupby(by=list(fields_of_interest.keys())[0:8],dropna=False).apply(lambda x:pd.Series({'n_objects':len(x[list(fields_of_interest.keys())[9]]),'Volume_analyzed':x.Volume_analyzed.unique()[0] if 'Volume_analyzed' in x.columns else None,
    #                                                                                                       'Area':x.Area.unique()[0] if 'Area' in x.columns else None,'Minor':x.Minor_axis.unique()[0] if 'Minor_axis' in x.columns else None,'ESD':x.ESD.unique()[0] if 'ESD' in x.columns else None,'Biovolume':x.Biovolume.unique()[0] if 'Biovolume' in x.columns else None,'Pixel':x.Pixel.unique()[0] if 'Pixel' in x.columns else None})).reset_index()
    cols = subset_df.columns
    fig = make_subplots(rows=2, cols=2,#shared_xaxes='columns', # Use shared x/yaxes to fix axis
                       # subplot_titles=("","Number of objects", "Volume analyzed", "Pixel:size ratio", "Area","Minor ellipsoidal axis","Equivalent spherical diameter","Biovolume"),
                      #  specs=[[{"type": "scattergeo", "rowspan": 7}, {"type": "scatter"}],[None, {"type": "scatter"}],[None, {"type": "scatter"}],[None, {"type": "scatter"}],[None, {"type": "scatter"}],[None, {"type": "scatter"}],[None, {"type": "scatter"}]],
                        specs=[[{"type": "scattergeo", "rowspan": 2}, {"type": "scatter"}],[None, {"type": "scatter"}]],
                       column_widths=[0.5,0.5],vertical_spacing=0.07)
    # Map, left panel
    fig.add_trace(go.Scattergeo(lon=subset_df.Longitude, lat=subset_df.Latitude,
                                  hovertext='Sample ID: '+subset_df.Sample_ID, hoverinfo="lon+lat+text",
                                marker=dict(color='black'),geojson="natural earth",showlegend=False),
                   row=1,col=1)
    # Scatterplot 1, top-right panel
    fig.add_trace(go.Scatter(x=subset_df.Sample_ID, y=subset_df.n_objects,
                             hovertext="Sample ID: "+subset_df.Sample_ID,hoverinfo="text",marker=dict(color='black'),
                             showlegend=False,visible=True),row=1,col=2) # Use mode='markers' for points only
    # Scatterplot 2, bottom-right panel
    profiles_df=profiles_df.sort_values(['Depth_group','Sample_ID'])
    colors=sequential_hcl(h = [-220,20], c = [0, 90, 0], l = [95, 30], power = [1.5,1.5], rev = True)(len(profiles_df['Sample_ID'].unique()))

    fig.add_trace(go.Box(y=df.Depth_group, x=df.n_objects,
                         fillcolor='black',line=dict(color='#222222'),
                         hoverinfo='skip',
                         # error_y=dict(type="data",array=subset_df.std_objects),
                         # hovertext="Sample ID: " + subset_df.Sample_ID, hoverinfo="text",
                         notched=True,
                         showlegend=False, visible=True), row=2, col=2)
    fig.update_traces(orientation='h',row=2,col=2)

    for i,station in enumerate(np.sort(profiles_df.Sample_ID.unique())):
        subset_profiles= profiles_df[profiles_df['Sample_ID']==station].sort_values(['Depth_group'])
        fig.add_trace(go.Scatter(y=subset_profiles.Depth_group, x=subset_profiles.n_objects,
                                 marker=dict(color=colors[i]), line=dict(color=colors[i]),
                                 name="Sample "+str(station), hovertext=subset_profiles.Sample_ID, hoverinfo="text",
                                 legendrank=i,
                                 # line_shape='spline',
                                 # error_y=dict(type="data",array=subset_df.std_objects),
                                 # hovertext="Sample ID: " + subset_df.Sample_ID, hoverinfo="text",
                                 mode='markers',showlegend=True, visible="legendonly"), row=2, col=2)





    fig.update_yaxes(type="linear", autorange="reversed", row=2, col=2)

    button_scatter1 = [dict(method="restyle",  # argument to change data
                            args=[{'x': [subset_df['Sample_ID'], 'undefined'],
                                   'y': [subset_df[cols[k]], 'undefined'],
                                   'visible': [True, True, True]}, [1]], # visible should have the same length as number of subplots, second argument specifies which subplot is updated
                            label=cols[k]) for k in np.where(
        np.isin(list(cols), ['n_objects', 'Volume_analyzed', 'Area', 'Minor', 'ESD', 'Biovolume', 'Pixel']))[0]]


    fig.update_layout(updatemenus=list([dict(active=0, buttons=button_scatter1, x=0.87, y=1.1, xanchor='left', yanchor='top')]),
        title={'text': 'Cruise:' + subset_df["Cruise"][0] + "<br> Hover on datapoint to see sample ID",'xanchor': 'center', 'yanchor': 'top', 'x': 0.5},
        xaxis={'title': 'Sample ID','nticks': 2, 'tickangle': 0},xaxis2={'title': 'Number of objects imaged'},
        yaxis={'title': 'Variable selected in dropdown menu <br>(fill unit in standardizer)'},yaxis2={'title':'Depth (m)'},
        hovermode="x",
       boxmode='group')

    fig.for_each_xaxis(lambda x: x.update(showgrid=False))
    fig.for_each_yaxis(lambda x: x.update(showgrid=False))
    fig.update_yaxes(type="log", row=1, col=2)


    fig.show()
    path_to_plot=path_to_plotdir / "Plot_standardizer_{}_{}.html".format(df_instrument.at[project_id,"Instrument"],str(project_id))
    print("Saving plots to: ",path_to_plot,sep=" ")
    fig.write_html(path_to_plot)
    """  
    # Currently, API function get_object_set only accepts object related fields, not acquisition/process/sample fields 
    # Search visible and accessible projects based on instrument filter (Authentication in EcoTaxa required)
    with ecotaxa_py_client.ApiClient() as client:
        api = authentification_api.AuthentificationApi(client)
        token = api.login(LoginReq(
            username=cfg_pw['ecotaxa_user'],
            password=cfg_pw['ecotaxa_pass']
        ))
    configuration = ecotaxa_py_client.Configuration(host="https://ecotaxa.obs-vlfr.fr/api", access_token=token,
                                                        discard_unknown_keys=True)
    fields_of_interest = "obj.latitude,obj.longitude,obj.depth_min,obj.depth_max," + ','.join(['fre.' + item[item.find("_") + 1:] for item in fields_of_interest if isinstance(item, str) and 'object_' in item])

   # with ecotaxa_py_client.ApiClient(configuration) as api_client:
   #     api_instance =  objects_api.ObjectsApi(api_client)
   #     objects= api_instance.get_object_set(int(project_id),ProjectFilters(statusfilter="PVD"),fields=fields_of_interest)
   df=pd.DataFrame(objects['details'], columns=fields_of_interest.split(','))

   """

# Function (3): Perform EcoTaxa export files standardization and harmonization based on standardizer spreadsheets
def standardization_func(config_path,standardizer_path,project_id):
    """
       Objective: This function uses the instrument-specific standardizer spreadsheet to standardize export variable name and units
       This should facilitate the size spectrum analysis of projects collected with UVP, IFCB, Zooscan, and ISIIS
       :param standardizer_path: Full path of the standardizer spreadsheet containing project ID of interest
       :param project_id: unique integer for project ID to be standardized
       :return: standardized project dataframe
       """
    # Open configuration file to get: EcoTaxa authentication info
    path_to_config_usr = Path(config_path).expanduser()
    with open(path_to_config_usr, 'r') as config_file:
        cfg = yaml.safe_load(config_file)
    path_to_git = Path(cfg['git_dir']).expanduser()

    # Open standardizer spreadsheet
    df_instrument = pd.read_excel(standardizer_path,index_col=0)

    # Control. Standardization only works for Zooscan projects at the moment
    if project_id not in df_instrument.index:
        print('Project ID is not included in standardizer. Quitting')
        return
        if df_instrument.loc[project_id]['Instrument'] not in [ 'Zooscan','IFCB']:
            print('Standardization only applies to Zooscan and IFCB projects. Quitting')
            return

    # Retrieve fields to standardize in standardizer (indicated by variable_field)
    columns_of_interest = ['Object_ID', "Sample_ID"] + [string for string in list(df_instrument.columns) if"_field" in string]
    fields=fields_of_interest = dict(zip(columns_of_interest,['object_id', "sample_id"] + df_instrument.loc[project_id][columns_of_interest[2:]].values.flatten().tolist()))
    fields=dict((key.replace("_field", ""), value.split(","))  if isinstance(value, str) else (key.replace("_field", ""), "".split(",")) for key, value in fields.items() )
    fields_of_interest = dict((key.replace("_field", ""), value) for key, value in fields_of_interest.items() if isinstance(value, str))
    fields_of_interest_list = list(( key, field) for key, value in fields_of_interest.items() for field in value.split(','))
    fields_of_interest_series =pd.Series(list(value[1] for value in fields_of_interest_list),list(value[0] for value in fields_of_interest_list))

    # Retrieve index to skip
    index_NA=df_instrument.loc[project_id]['Missing_index']
    if str(index_NA)!='nan':
        index_NA=[int(id) for id in str(df_instrument.loc[project_id]['Missing_index']).split(",")]
    else:
        index_NA=None


    path_to_data = path_to_git / Path(cfg['dataset_subdir']) / df_instrument.at[project_id, "Instrument"]
    if len(list(path_to_data.glob('ecotaxa_export_{}*'.format(str(project_id)))))>0:

        # Load export tsv file
       # df = pd.read_table(list(path_to_data.glob('ecotaxa_export_{}*'.format(str(project_id))))[0],usecols=fields_of_interest.values(),skiprows=index_NA)
        df = pd.read_table(list(path_to_data.glob('ecotaxa_export_{}*'.format(str(project_id))))[0],usecols=fields_of_interest_series.values, skiprows=index_NA)
        old_columns = df.columns
        # Standardize column names
        df.columns = [(list(fields_of_interest_series.index)[list(fields_of_interest_series.values).index(value)]) for value in  old_columns]  # pd.MultiIndex.from_arrays([[list(fields_of_interest.keys())[list(fields_of_interest.values()).index(value)] for value in df.columns],df.columns])
        # Transform variables based on field units
        units= [string for string in list(df_instrument.columns) if "_unit" in string]
        units_of_interest = dict(zip(units,df_instrument.loc[project_id][units].values.flatten().tolist()))
       # Test whether all units are defined:                             [re.sub(r'^-?[0-9]_times_10power-?[0-9]_','',str(unit)) for unit in list(units_of_interest.values())]
        if any([unit not in ureg for unit in list(units_of_interest.values()) if str(unit)!='nan']):
            print('Quitting. Please fill out standardizer units with units from the list below or define your custom unit in {} using known units:\n'.format(path_to_git/'Scripts'/'units_def.txt'),full_list_units, sep='')
            return
        # Convert pixel_to_size ratio (pixel per millimeter)
        if str(ureg(units_of_interest['Pixel_unit'].split('_per_')[0]).to_base_units().units)=='pixel':
            pixel_size_ratio = PQ(list(df.Pixel/ureg(units_of_interest['Pixel_unit'].split('_per_')[1]).to('millimeter')), 'pixel/millimeter')
        if str(ureg(units_of_interest['Pixel_unit'].split('_per_')[0]).to_base_units().units)=='meter':
            pixel_size_ratio = PQ(list(1/(df.Pixel*ureg(units_of_interest['Pixel_unit'].split('_per_')[0]).to('millimeter'))), 'pixel/millimeter')
        # Use pint units system to convert units in standardizer spreadsheet to standard units
        # (degree for latitude/longitude, meter for depth, multiple of micrometer for plankton size)
        df_standardized = df.assign(Project_ID=project_id,
                       Station_ID=df.Station_ID if 'Station_ID' in df.columns else pd.NA,
                       Profile=df.Profile_ID if 'Profile_ID' in df.columns else pd.NA,
                       Sampling_date=df.Sampling_date if 'Sampling_date' in df.columns else pd.NA,
                       Sampling_time=df.Sampling_time if 'Sampling_time' in df.columns else pd.NA,
                       Pixel=pixel_size_ratio if 'Pixel' in df.columns else pd.NA,
                       Volume_analyzed=1000*(list(df.Volume_analyzed)*ureg(units_of_interest['Volume_analyzed_unit']).to(ureg.cubic_meter))if 'Volume_analyzed' in df.columns else pd.NA, # cubic dm
                       Latitude=list(df.Latitude)*ureg(units_of_interest['Latitude_unit']).to(ureg.degree) if 'Latitude' in df.columns else pd.NA, # degree decimal
                       Longitude=list(df.Longitude) * ureg(units_of_interest['Longitude_unit']).to(ureg.degree) if 'Longitude' in df.columns else pd.NA,  # degree decimal
                       Depth_min=list(df.Depth_min) * ureg(units_of_interest['Depth_min_unit']).to(ureg.meter) if 'Depth_min' in df.columns else pd.NA,  # meter
                       Depth_max=list(df.Depth_max) * ureg(units_of_interest['Depth_max_unit']).to(ureg.meter) if 'Depth_max' in df.columns else pd.NA,  # meter
                       Area=1e06*((df.Area*[1*ureg(units_of_interest['Area_unit']).to(ureg.square_pixel).magnitude for ntimes in range(df[['Area']].shape[1])])/((np.vstack([pixel_size_ratio**2 for ntimes in range(df[['Area']].shape[1])]).T)[0,:])) if all(pd.Series(['Area','Pixel']).isin(df.columns)) else pd.NA, # square micrometers
                       Biovolume=1e09*(df.Biovolume*[1*ureg(units_of_interest['Biovolume_unit']).to(ureg.cubic_pixel).magnitude for ntimes in range(df[['Biovolume']].shape[1])]/((np.vstack([pixel_size_ratio**3 for ntimes in range(df[['Biovolume']].shape[1])]).T)[0,:])) if all(pd.Series(['Biovolume','Pixel']).isin(df.columns)) else pd.NA, # cubic micrometers
                       Minor_axis=1e03*(df.Minor_axis*[1*ureg(units_of_interest['Minor_axis_unit']).to(ureg.pixel).magnitude for ntimes in range(df[['Minor_axis']].shape[1])]/((np.vstack([pixel_size_ratio for ntimes in range(df[['Minor_axis']].shape[1])]).T)[0,:])) if all(pd.Series(['Minor_axis','Pixel']).isin(df.columns)) else pd.NA, # cubic micrometers
                       ESD=1e03*(df.ESD*[1*ureg(units_of_interest['ESD_unit']).to(ureg.pixel).magnitude for ntimes in range(df[['ESD']].shape[1])]/((np.vstack([pixel_size_ratio for ntimes in range(df[['ESD']].shape[1])]).T)[0,:])) if all(pd.Series(['ESD','Pixel']).isin(df.columns)) else pd.NA) # micrometers
        # Assume volume analyzed (5 mL) for IFCB projects
        if df_instrument.loc[project_id]['Instrument'] in ['IFCB'] and all(df_standardized[['Volume_analyzed']]==pd.NA):
            df_standardized['Volume_analyzed']=PQ(5*ureg('milliliter').to('liter')).magnitude

        # Convert volume analyzed to volume imaged to account for samples dilution or fractionation
        if 'Dilution' in df.columns:
            dilution_factor = np.where(df.Dilution > 1, df.Dilution, 1 / df.Dilution)
        else:
            dilution_factor=1

        df_standardized = df_standardized.assign(Volume_imaged=df_standardized.Volume_analyzed/dilution_factor) # cubic meters
        # Generate metadata sheet
        df_standardized=df_standardized[['Project_ID','Cruise', 'Station_ID', 'Profile','Sample_ID' , 'Latitude', 'Longitude', 'Sampling_date', 'Sampling_time','Depth_min', 'Depth_max', 'Volume_analyzed', 'Volume_imaged',
                                    'Object_ID', 'Category', 'Minor_axis', 'ESD', 'Area', 'Biovolume']]
        df_standardized_metadata=pd.DataFrame({'Variables':df_standardized.columns,'Variable_types':df_standardized.dtypes,
        'Units/Values/Timezone':['','','','','','degree','degree','yyyymmdd (GMT)','hhmmss (GMT)','meter','meter','cubic_decimeter','cubic_decimeter','','']+['micrometer' for ntimes in range(df_standardized[['Minor_axis']].shape[1])]+['micrometer' for ntimes in range(df_standardized[['ESD']].shape[1])]+['square_micrometer' for ntimes in range(df_standardized[['Area']].shape[1])]+['cubic_micrometer' for ntimes in range(df_standardized[['Biovolume']].shape[1])],
        'Description':['Project ID in EcoTaxa','Project cruise','Station ID (native format)','Profile ID (native format)','Sample ID (native format)','Latitude','Longitude','Sampling date','Sampling time','Minimum sampling depth','Maximum sampling depth','Volume analyzed (not accounting for sample dilution and/or fractionation)','Volume imaged (accounting for sample dilution and/or fractionation)','Object ID (native format)','Object assigned taxonomy']+ ['Object minor ellipsoidal axis derived from '+ field if field!='' or len(fields['Minor_axis'])>1 else 'Object minor ellipsoidal axis' for field in fields['Minor_axis'] ]+ ['Object equivalent spherical diameter derived from '+ field if field!='' or len(fields['ESD'])>1 else 'Object equivalent spherical diameter' for field in fields['ESD'] ]+ ['Object surface area derived from '+ field if field!='' or len(fields['Area'])>1 else 'Object surface area' for field in fields['Area'] ]+['Object biovolume derived from '+ field if field!='' or len(fields['Biovolume'])>1 else 'Object biovolume' for field in fields['Biovolume']]})

        # Save standardized dataframe
        path_to_standard_dir = path_to_data.parent.parent / 'raw_standardized' / str(df_instrument.loc[project_id]['Instrument'])
        path_to_standard_dir.mkdir(parents=True, exist_ok=True)
        path_to_standard_file=path_to_standard_dir / 'standardized_export_{}.tsv'.format(str(project_id))
        path_to_metadata_standard=path_to_standard_dir / 'standardized_export_{}_metadata.tsv'.format(str(project_id))
        path_to_standard_plot = path_to_git / 'figures' /  'Standardized_export_{}.html'.format(str(project_id))

        print('Saving standardized export (data and metadata) to',path_to_standard_file,sep=' ')
        #with pd.ExcelWriter(str(path_to_standard_file),engine="xlsxwriter") as writer:
        df_standardized.to_csv(path_to_standard_file, sep="\t") # to_excel(writer, sheet_name='Data', index=False)
        df_standardized_metadata.to_csv(path_to_metadata_standard, sep="\t")  # to_excel(writer, sheet_name='Metadata',index=False)

        # Interactive plots
        fig = make_subplots(rows=2, cols=2,
                            specs=[[{"type": "scattergeo", "rowspan": 2}, {"type": "scatter"}],[None, {"type": "scatter"}]],
                            column_widths=[0.5, 0.5], vertical_spacing=0.07)
        # Map, left panel
        summary_df_standardized = df_standardized.groupby( ['Sample_ID', 'Station_ID', 'Latitude', 'Longitude', 'Profile', 'Volume_imaged'],dropna=False).apply(lambda x: pd.Series({'Abundance': x.Object_ID.count() / ( x.Volume_analyzed.unique()[0]), # individuals per liter
                                                                                                                                                                       'Average_diameter':np.nanmean(2*np.power(x.Area/np.pi,0.5)) if 'Area' in df_standardized.columns else np.nanmean(x.ESD), # micrometer
                                                                                                                                                                       'Std_diameter':np.nanstd(2*np.power(x.Area/np.pi,0.5)) if 'Area' in df_standardized.columns else np.nanstd(x.ESD)})).reset_index() # micrometer

        units=['','',' (degree)',' (degree)','',' (L)',' (particles/L)',' (µm)',' (µm)']
        fig.add_trace(go.Scattergeo(lon=summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample_ID']).Longitude, lat=summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample_ID']).Latitude,
                                    #size=summary_df_standardized.Abundance,
                                    hovertext='Sample ID: ' + summary_df_standardized.dropna(subset=['Latitude','Longitude','Sample_ID']).Sample_ID, hoverinfo="lon+lat+text",
                                    marker=dict(color='black'), geojson="natural earth", showlegend=False),row=1, col=1)
        # Scatterplot 1, top-right panel
        cols=summary_df_standardized.columns
        fig.add_trace(go.Scatter(x=summary_df_standardized.dropna(subset=['Sample_ID','Abundance']).Sample_ID, y=summary_df_standardized.dropna(subset=['Sample_ID','Abundance']).Abundance,
                               #  error_y=dict(type="data",array=summary_df_standardized.Std_diameter),
                                 hovertext="Sample ID: " + summary_df_standardized.dropna(subset=['Sample_ID','Abundance']).Sample_ID, hoverinfo="text",
                                 marker=dict(color='black'),mode='markers',
                                 showlegend=False, visible=True), row=1, col=2)

        button_scatter1 = [dict(method="restyle",  # argument to change data
                                args=[{'x': [summary_df_standardized['Sample_ID'], 'undefined'],
                                       'y': [summary_df_standardized[cols[k]], 'undefined'],
                                       'visible': [True, True, True]}, [1]],
                                # visible should have the same length as number of subplots, second argument specifies which subplot is updated
                                label=cols[k] + units[k]) for k in np.where(
            np.isin(list(cols), ['Abundance', 'Average_diameter', 'Std_diameter']))[0]]

        fig.update_layout(
            updatemenus=list([dict(active=0, buttons=button_scatter1, x=0.87, y=1.1, xanchor='left', yanchor='top')]),
            title={'text': 'Cruise:' + df_standardized["Cruise"][0] + "<br> Hover on datapoint to see sample ID",
                   'xanchor': 'center', 'yanchor': 'top', 'x': 0.5},
            xaxis={'title': 'Sample ID','nticks': 2, 'tickangle': 0},
            yaxis={"tickmode": "array",'title': 'Variable selected in dropdown menu'},
            hovermode="x",
            boxmode='group')
        fig.update_yaxes(type="log", row=1, col=2)

        fig.for_each_xaxis(lambda x: x.update(showgrid=False))
        fig.for_each_yaxis(lambda x: x.update(showgrid=False))
        fig.write_html(path_to_standard_plot)
        print('Saving standardized export plot to', path_to_standard_plot, sep=' ')
    else: print('Export file for project {} not found. Please run ecotaxa_API_export_data.py to export file with the API or export manually'.format(str(project_id)))


if __name__ == '__main__':
    print('Assisting standardizer completion. \nPlease enter the following information (attention: method sensitive to empty space and lower/upper case. Do not use quotation)\n')
    funcs_args=input('Press:\n 1 for printing all fields for Ecotaxa project using the API (downloading export files is not required)\n 1b for printing all possible units\n 2 for generating/saving interactive plot of raw ecotaxa export\n 3 for generating/saving/plotting ecotaxa export standardization\n')

    if funcs_args == '1b':
        unitfile_args=input('Custom units full path:')
        if len(unitfile_args)>0:
            filling_standardizer_unit_func(custom_units_path=unitfile_args)
        else:
            filling_standardizer_unit_func()
        quit()


    # Example: 3
    configfile_args = input('Configuration file (containing Ecotaxa authentication info for 1, storage info for 2/3) full path:')
    # Example: ~/GIT/PSSdb/scripts/Ecotaxa_API.yaml
    standardizerfile_args = input('Standardizer spreadsheet full path:')
    # Example: ~/GIT/PSSdb/raw/project_Zooscan_standardizer.xlsx
    # Example: ~/GIT/PSSdb/raw/project_IFCB_standardizer.xlsx
    project_args=input('Unique project ID of interest (integer):')
    # Example: 377
    if funcs_args=='1':
        if len(configfile_args) > 0 and len(standardizerfile_args) > 0:
            if len(project_args) == 0:
                filling_standardizer_field_func(config_path=configfile_args, standardizer_path=standardizerfile_args)
            else:
                filling_standardizer_field_func(config_path=configfile_args, standardizer_path=standardizerfile_args,project_id=int(project_args))
        else: print('Missing information. Please run again')


    if funcs_args=='2':
        if len(configfile_args) > 0 and len(standardizerfile_args) > 0:
            if len(project_args) > 0:
                diagnostic_plots_func(config_path=configfile_args, standardizer_path=standardizerfile_args,
                                      project_id=int(project_args))

        else:
            print('Missing information. Please run again')

    if funcs_args=='3':
        if len(configfile_args) > 0 and len(standardizerfile_args) > 0 and len(project_args) > 0:
            standardization_func(config_path=configfile_args, standardizer_path=standardizerfile_args, project_id=int(project_args))

        else:
            print('Missing information. Please run again')