## Objective: This script uses python URL scraping to search for taxonomic information of specific taxon (contained in EcoTaxa export tables) on the World Register of Marine Species (WoRMNS)

##Documentation: https://www.marinespecies.org/ (WORMS) and https://www.thepythoncode.com/article/extracting-and-submitting-web-page-forms-in-python (Python URL scraping)

## Python modules:

# Modules for data and path handling:
import pandas as pd
import numpy as np
from pathlib import Path  # Handling of path object
import os

# Modules for webpage handling:
import urllib3
import re
import string
import requests
from bs4 import BeautifulSoup
import urllib
from urllib.parse import urljoin
from requests_html import HTMLSession


# Modules for EcoTaxa authentication and query:
import ecotaxa_py_client
from ecotaxa_py_client.api import objects_api
from ecotaxa_py_client.api import projects_api
from ecotaxa_py_client.model.project_filters import ProjectFilters
from ecotaxa_py_client.api import authentification_api
from ecotaxa_py_client.model.login_req import LoginReq
from ecotaxa_py_client.api import taxonomy_tree_api
from ecotaxa_py_client.model.taxon_model import TaxonModel

import yaml # requires installation of PyYAML package
path_to_config_usr=Path('scripts/Ecotaxa_API_pw.yaml').expanduser()
  with open(path_to_config_usr ,'r') as config_file:
    cfg_pw = yaml.safe_load(config_file)

  with ecotaxa_py_client.ApiClient() as client:
    api = authentification_api.AuthentificationApi(client)
    token = api.login(LoginReq(username=cfg_pw['ecotaxa_user'],password=cfg_pw['ecotaxa_pass']))

configuration = ecotaxa_py_client.Configuration(host = "https://ecotaxa.obs-vlfr.fr/api",access_token=token, discard_unknown_keys=True)


## Useful functions:
# Retrieve hidden search/query (form) from URL
def get_all_forms(url):
    """Returns all form tags found on a web page's `url` """
    with  HTMLSession() as session:
        # GET request
        res = session.get(url)
        # for javascript driven website
        # res.html.render()
        soup = BeautifulSoup(res.html.html, "html.parser")

    return soup.find_all("form")

# Retrieve inputs required for the form
def get_form_details(form):
    """Returns the HTML details of a form,
    including action, method and list of form controls (inputs, etc)"""
    details = {}
    # get the form action (requested URL)
    action = form.attrs.get("action").lower()
    # get the form method (POST, GET, DELETE, etc)
    # if not specified, GET is the default in HTML
    method = form.attrs.get("method", "get").lower()
    # get all form inputs
    inputs = []
    for input_tag in form.find_all("input"):
        # get type of input form control
        input_type = input_tag.attrs.get("type", "text")
        # get name attribute
        input_name = input_tag.attrs.get("name")
        # get the default value of that input tag
        input_value =input_tag.attrs.get("value", "")
        # add everything to that list
        inputs.append({"type": input_type, "name": input_name, "value": input_value})
    details["action"] = action
    details["method"] = method
    details["inputs"] = inputs
    return details

# Fill the inputs for a form
def get_form_inputs(form):
    """Returns the inputs for a HTML form"""
    form_details=get_form_details(form)
    data = {}
    for input_tag in form_details["inputs"]:
        if input_tag["type"] == "hidden" or input_tag["type"] == "checkbox":
            # if it's hidden, use the default value
            data[input_tag["name"]] = input_tag["value"]
        elif input_tag["type"] == "select":
            for i, option in enumerate(input_tag["values"], start=1):
                # iterate over available select options
                if option == input_tag["value"]:
                    print(f"{i} # {option} (default)")
                else:
                    print(f"{i} # {option}")
            choice = input(f"Enter the option for the select field '{input_tag['name']}' (1-{i}): ")
            try:
                choice = int(choice)
            except:
                # choice invalid, take the default
                value = input_tag["value"]
            else:
                value = input_tag["values"][choice - 1]
            data[input_tag["name"]] = value
        elif input_tag["type"] != "submit":
            # all others except submit, prompt the user to set it
             value = input(f"Enter the value of the field '{input_tag['name']}' (type: {input_tag['type']}): ")
            data[input_tag["name"]] = value
    return data

# Return results (Rank, Phylum, Class, Genus, URL, Reference, Citation) of the form
def get_form_hierarchy(hierarchy,form):
 with  HTMLSession() as session:
    hierarchy_levels = hierarchy.split('>')
    df_hierarchy = pd.DataFrame({'EcoTaxa_hierarchy': hierarchy,'Full_hierarchy': '','Rank': '','Category': '','Domain': '','Phylum': '', 'Class': '','Order':'','Family':'','Genus': '','Functional_group': '','WORMS_ID': '','Reference': '','Citation': ''}, index=[0])
    # Loop through hierarchy levels starting with the last level
    for level in hierarchy_levels[::-1]:
        data = {}
        for input_tag in form["inputs"]:
            if input_tag["type"] == "hidden" or input_tag["type"] == "checkbox":
                # if it's hidden, use the default value
                data[input_tag["name"]] = input_tag["value"]
            elif input_tag["type"] != "submit":
                # all others except submit, prompt the user to set it
                # value = input(f"Enter the value of the field '{input_tag['name']}' (type: {input_tag['type']}): ")
                data[input_tag["name"]] = level
        data['marine']='0' # Re-assign value for marine checkbox to 0 to allow for non-marine species
        data['vOnly'] = '0'  # Re-assign value for extant checkbox to 0 to allow for unaccepted
        # print(data)
        url_search = urljoin(url, form["action"])
        # Post/get form with corresponding inputs
        if form["method"] == "post":
            res = session.post(url_search, data=data)
        elif form["method"] == "get":
            res = session.get(url_search, params=data)
        # Transform form query into xml table and save results in dataframe
        soup = BeautifulSoup(res.html.html, "lxml")
        if len(soup.find_all(class_="list-group-item"))==0: # If post result in a single page
            taxo_soup = BeautifulSoup(res.html.html, 'lxml')
            if len(taxo_soup.find_all(class_="alert alert-warning"))>0:
                continue # Skip current level if hierarchy not found
        else:
            search_response = pd.DataFrame({'Taxon': [item.getText()  for item in soup.find_all(class_="list-group-item") if any(re.findall(r'unaccepted|uncertain|unassessed',item.getText()))==False],
                 'Link': [item.find('a').get('href')  for item in soup.find_all(class_="list-group-item") if any(re.findall(r'unaccepted|uncertain|unassessed',item.getText()))==False]})

            if len(search_response)==0:
                     continue # Skip current level if hierarchy not found
            if search_response.Taxon[0].split(' ')[0].lower() != level.lower():
                     continue # Skip current level if hierarchy not found
            else:
                     connection_webpage = requests.get(urljoin('https://www.marinespecies.org', search_response.Link[0]),data=data)
                     taxo_soup = BeautifulSoup(connection_webpage.content, 'lxml')
                     if len(taxo_soup.find_all(class_="alert alert-warning"))>0:
                         res = session.post(urljoin('https://www.marinespecies.org', search_response.Link[0]),data=data)
                         soup = BeautifulSoup(res.html.html, "lxml")
                         taxo_soup = BeautifulSoup(res.html.html, 'lxml')

        fields=[item.getText() if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else '' for item in taxo_soup.find_all(class_="col-xs-12 col-sm-4 col-lg-2 control-label" )]
        Status=re.sub(r'[' + '\n' + '\xa0' + ']', '', taxo_soup.find_all(class_="leave_image_space")[fields.index('Status')].getText()) if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else ''
        if 'unaccepted' in Status:  # Re-assign hierarchy level with accepted name if status in unaccepted
            if 'Accepted Name' not in fields:
                continue
            else:
            level = re.sub(r'[' + '\n' + '\xa0' + ']', '', taxo_soup.find_all(class_="leave_image_space")[fields.index('Accepted Name')].getText()) if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else level
            data['tName']=level
            res = session.post(url_search, data=data)
            # Transform form query into xml table and save results in dataframe
            soup = BeautifulSoup(res.html.html, "lxml")
            search_response = pd.DataFrame({'Taxon': [item.getText() for item in soup.find_all(class_="list-group-item") if any(re.findall(r'unaccepted|uncertain|unassessed', item.getText())) == False],'Link': [item.find('a').get('href') for item in soup.find_all(class_="list-group-item") if any(re.findall(r'unaccepted|uncertain|unassessed',item.getText())) == False]})
            if len(search_response)==0:
                continue
            else:
            connection_webpage = requests.get(urljoin('https://www.marinespecies.org', search_response.Link[0]),data=data)
            taxo_soup = BeautifulSoup(connection_webpage.content, 'lxml')
            fields = [item.getText() if len(taxo_soup.find_all(class_="alert alert-warning")) == 0 else '' for item in taxo_soup.find_all(class_="col-xs-12 col-sm-4 col-lg-2 control-label")]


        script =taxo_soup.find_all('script')
        Functional_group=[item.getText()[item.getText().find('Functional group')+22:[sub.start() for sub in re.finditer('&nbsp',item.getText()) if sub.start()>item.getText().find('Functional group')+22][0]]  for item in script if 'Functional group' in item.getText()]
        Category = np.where(taxo_soup.find_all(class_="leave_image_space")[1].getText().split('\n')[1] == 'Biota','Living', 'NA').tolist() if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else ''
        dict_hierarchy = {re.sub(r'[' + string.punctuation + ']', '', level.split('\xa0')[1]): level.split('\xa0')[0] for level in taxo_soup.find_all(class_="leave_image_space")[1].getText().split('\n') if '\xa0' in level} if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else dict({'':''})
        full_hierarchy = '>'.join([level.split('\xa0')[0] + level.split('\xa0')[1] for level in taxo_soup.find_all(class_="leave_image_space")[1].getText().split('\n') if '\xa0' in level]) if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else ''
        Domain = dict_hierarchy['Kingdom'] if 'Kingdom' in dict_hierarchy.keys() else ''
        Genus = dict_hierarchy['Genus'] if 'Genus' in dict_hierarchy.keys() else ''  # [item.getText().split('\xa0')[0] for item in taxo_soup.findAll('li') if 'Genus' in item.getText()][0]
        Family=dict_hierarchy['Family'] if 'Family' in dict_hierarchy.keys() else '' # [item.getText().split('\xa0')[0] for item in taxo_soup.findAll('li') if 'Family' in item.getText()][0]
        Order=dict_hierarchy['Order'] if 'Order' in dict_hierarchy.keys() else '' # [item.getText().split('\xa0')[0] for item in taxo_soup.findAll('li') if 'Order' in item.getText()][0]
        Class = dict_hierarchy[ 'Class'] if 'Class' in dict_hierarchy.keys() else '' # [item.getText().split('\xa0')[0] for item in taxo_soup.findAll('li') if 'Class' in item.getText()][0]
        Phylum = dict_hierarchy['Phylum'] if 'Phylum' in dict_hierarchy.keys() else '' # [item.getText().split('\xa0')[0] for item in taxo_soup.findAll('li') if 'Phylum' in item.getText()][0]
        Rank = re.sub(r'[' + '\n' + ' ' + ']', '', taxo_soup.find_all(class_="leave_image_space")[fields.index('Rank')].getText()) if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else ''
        URL = re.sub(r'[' + '(' + ')' + ']', '', taxo_soup.find_all(class_="aphia_core_cursor-help")[fields.index('AphiaID')].getText()) if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else ''
        Original_reference = taxo_soup.find_all(class_="correctHTML")[0].getText() if len(taxo_soup.find_all(class_="correctHTML"))>0 else ''
        Taxonomic_citation = [re.sub(r'[' + '\n' + ' ' + ']', '', item.getText()) for item in taxo_soup.find_all(class_="col-xs-12 col-sm-8 col-lg-10 pull-left") if item.attrs['id'] == 'Citation'][0] if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else ''
        df_hierarchy = pd.DataFrame({'EcoTaxa_hierarchy': hierarchy,
                                         'Full_hierarchy':full_hierarchy,
                                         'Rank': Rank,
                                         'Category': Category,
                                         'Domain': Domain,
                                         'Phylum': Phylum,
                                         'Class': Class,
                                         'Order': Order,
                                         'Family':Family,
                                         'Genus': Genus,
                                         'Functional_group': Functional_group if len(Functional_group) > 0 else '',
                                         'WORMS_ID': URL,
                                         'Reference': Original_reference,
                                         'Citation': Taxonomic_citation
                                     }, index=[0])
        break
    return df_hierarchy


## Workflow starts here:
#Retrieve fields of interest (including taxonomic hierarchies) in EcoTaxa projects
# see: https://github.com/ecotaxa/ecotaxa_py_client/blob/main/docs/ObjectsApi.md#get_object_set
fields_of_interest = "obj.orig_id,txo.id,txo.display_name,img.imgid,img.file_name,img.thumb_file_name,obj.depth_min"; colnames=['ROI_ID','ROI_Annotation_ID','ROI_Annotation','Image_ID','Image_filename','Thumb_filename','Depth']
fields_of_interest = "txo.id,txo.display_name";colnames=['ROI_Annotation_ID','ROI_Annotation']
with ecotaxa_py_client.ApiClient(configuration) as api_client:
     api_project_instance=projects_api.ProjectsApi(api_client)
     project=api_project_instance.search_projects(also_others=False, # Return visible projects with access rights
                                                                   title_filter='', # Optional title project filter
                                                                   instrument_filter='', # All instruments
                                                                   order_field='projid' # Sorting variable. Use instrument or projid
                                                                   )
     project_list = list(map(lambda x: str(x.projid),project))
     api_object_instance =  objects_api.ObjectsApi(api_client)
     objects= api_object_instance.get_object_set(project_list[0],ProjectFilters(statusfilter="PVD"),fields=fields_of_interest)
     df=pd.DataFrame(objects.details,columns=colnames).drop_duplicates()
     for project in project_list[1:]:
         objects = api_object_instance.get_object_set(project, ProjectFilters(statusfilter="PVD"),fields=fields_of_interest)
         df = pd.concat([df,pd.DataFrame(objects.details, columns=colnames)],axis=0).drop_duplicates()
     # df['ROI_ID']=df['ROI_ID'].apply(lambda x: ''.join(re.split('([0-9]+)', x)[0:-2]) + re.split('([0-9]+)', x)[-2].zfill(6))
     api_taxo_instance = taxonomy_tree_api.TaxonomyTreeApi(api_client)
     df['EcoTaxa_hierarchy'] = ['>'.join(api_taxo_instance.query_taxa(int(id)).lineage[::-1]) for id in df.ROI_Annotation_ID]
     df = df.sort_values(by=['EcoTaxa_hierarchy'])

# Retrieve the query form
url='https://www.marinespecies.org/aphia.php?p=search'
form=get_all_forms(url)[1]
# Retrieve form's inputs and add specific taxonomic query based on current taxon categories
form_details=get_form_details(form)
# Loop through EcoTaxa hierarchy and merge to annotated taxonomy
df_hierarchy=pd.DataFrame({})
for hierarchy in df.EcoTaxa_hierarchy: #hierarchy=df.EcoTaxa_hierarchy[0]
    print(hierarchy)
    data = get_form_hierarchy(hierarchy, form=form_details)
    # Update annotated taxonomy
    df_hierarchy=pd.concat([df_hierarchy, data],axis=0)
df_hierarchy=pd.merge(df[['ROI_Annotation','EcoTaxa_hierarchy']],df_hierarchy,on='EcoTaxa_hierarchy',how='right')
df_hierarchy.columns=['Category']+df_hierarchy.columns[1:].tolist()
df_hierarchy_metadata=pd.DataFrame({'Variables':df_hierarchy.columns,'Variable_types':df_hierarchy.dtypes,'Description':['Region of interest (object) annotation category','Full hierarchy of annotation category in EcoTaxa','Full taxonomic hierarchy of annotation category in World Register of Marine Species (WORMS)','Taxonomic rank of the annotation category in WORMS','Category of annotation category in EcoTaxa (e.g. Living)','Taxonomic domain/kingdom of the annotation category','Taxonomic phylum of the annotation category','Taxonomic class of the annotation category','Taxonomic order of the annotation category','Taxonomic family of the annotation category','Taxonomic genus of the annotation category', 'Functional group of the annotation category','Unique ID of the annotation category in the WORMS database','Reference for the annotation category description','Citation for the annotation category in WORMS']})
path_to_data=Path(os.getcwd()).expanduser() / 'raw' / 'plankton_annotated_taxonomy.xlsx'
with pd.ExcelWriter(str(path_to_data),engine="xlsxwriter") as writer:
    df_hierarchy.to_excel(writer, sheet_name='Data', index=False)
    df_hierarchy_metadata.to_excel(writer, sheet_name='Metadata', index=False)
df=pd.read_excel(path_to_data)