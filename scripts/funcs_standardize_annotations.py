## Objective: This script uses python URL scraping to search for taxonomic information of specific taxon (contained in EcoTaxa export tables) on the World Register of Marine Species (WoRMNS)

##Documentation: https://www.marinespecies.org/ (WORMS) and https://www.thepythoncode.com/article/extracting-and-submitting-web-page-forms-in-python (Python URL scraping)

## Requirements: Downloading chrome and chromedriver (path on l. 37) for your system (https://chromedriver.storage.googleapis.com/index.html?path=103.0.5060.134/)

## Python modules:
# Modules for data and path handling:
import pandas as pd
import numpy as np
from pathlib import Path  # Handling of path object
import os

# Modules for webpage handling/scraping:
import urllib3
import re
import string
import requests
from bs4 import BeautifulSoup
import urllib
from urllib.parse import urljoin
from requests_html import HTMLSession
# Selenium for web scraping: https://stackoverflow.com/questions/29385156/invoking-onclick-event-with-beautifulsoup-python
# Download chromedriver at:https://chromedriver.storage.googleapis.com/index.html?path=103.0.5060.134/
# Go in System Preferences > Security and Privacy > General > Allow
from selenium import webdriver # Use: pip3 install -U selenium
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager # Use: pip3 install webdriver-manager
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
options=Options()
options.add_argument("--headless") # Hidden windows
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import os
chromedriver = '~/Downloads/chromedriver' # Local path to the chromedriver. Note: Deprecated
#os.environ["webdriver.chrome.driver"] = chromedriver
import time

# Modules for colors:
import colorsys
import seaborn as sns
from colormap import rgb2hex
import matplotlib
import plotly
import plotly.express as px # Use: pip install plotly==5.8.0

## Functions start here:

def annotation_in_WORMS(hierarchy):
    """
    Objective: this function uses python requests to search taxonomic annotation on the World Register of Marine Sepcies (WORMS, https://www.marinespecies.org/).
    Note: Only accepted names will be used.\n
    :param hierarchy: Single or hierarchical taxonomic annotation to standardize with WORMS. If hierarchical taxonomic annotation, use > to separate taxonomic ranks
    :return: dataframe with corresponding rank, domain, phylum, class, order, family, genus, functional group, Taxon url, reference, citation, and URL for the annotation
    """
    hierarchy_levels = hierarchy.split('>')
    df_hierarchy = pd.DataFrame({'EcoTaxa_hierarchy': hierarchy, 'Full_hierarchy': '', 'Rank': '', 'Type': '', 'Domain': '', 'Phylum': '','Class': '', 'Order': '', 'Family': '', 'Genus': '', 'Functional_group': '', 'WORMS_ID': '', 'Reference': '','Citation': ''}, index=[0])
    url = 'https://www.marinespecies.org/aphia.php?p=taxlist'

    for annotation in hierarchy_levels[::-1]:
        with requests.Session() as session:
            data={'tName':annotation,'marine':'0','fossil':'0'}
            # Turn marine only and extant only search off
            session.post('https://www.marinespecies.org/aphia.php?p=search',data=data)
            taxon_list=session.post(url=url,data=data)
            soup = BeautifulSoup(taxon_list.content, "lxml")
            if len(soup.find_all(class_="list-group-item")) == 0:  # If post results in a single page:
                if len(soup.find_all(class_="alert alert-info")) > 0:
                    continue  # Skip current level if annotation not found
                else:
                    taxo_soup=soup
                    attributes_soup=''
                    script=''
                    annotation_webpage_with_attributes = session.get(urljoin(taxon_list.url, '#attributes'))
                    if annotation_webpage_with_attributes.ok:#len(list(filter(None,[id.get('id') if id.get('id') == "aphia_attributes_group_show" else id.get('href') if id.get('href') == '#attributes' else None for id in soup.findAll('a')]))) > 0:
                        attributes_soup = BeautifulSoup(annotation_webpage_with_attributes.content, 'lxml')
                        script = attributes_soup.find_all('script')
                        if len(attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})):
                            if 'No attributes found on lower taxonomic level' in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText() or len([item.getText() for item in script if 'Functional group' in item.getText()])==0:
                                attributes_soup = ''
                                script = ''
            else: # If post result in a list of taxa
                search_response = pd.DataFrame({'Taxon': [item.getText() for item in soup.find_all(class_="list-group-item") if any(re.findall(r'unaccepted|uncertain|unassessed',item.getText())) == False],'Link': [item.find('a').get('href') for item in soup.find_all(class_="list-group-item") if any(re.findall(r'unaccepted|uncertain|unassessed',item.getText())) == False]})
                if len(search_response) == 0:
                    continue  # Skip current level if annotation not found
                if search_response.Taxon[0].split(' ')[0].lower() != annotation.lower():
                    continue  # Skip current level if search response different from level
                else:
                    annotation_webpage = session.get(urljoin(url, search_response.Link[0]), data=data)
                    taxo_soup = BeautifulSoup(annotation_webpage.content, 'lxml')
                    if len(taxo_soup.find_all(class_="alert alert-info")) > 0:
                        res = session.post(urljoin(url, search_response.Link[0]), data=data)
                        soup = BeautifulSoup(res.content, "lxml")
                        taxo_soup = BeautifulSoup(res.content, 'lxml')
                        attributes_soup = ''
                        script = ''
                        annotation_webpage_with_attributes = session.get(urljoin(annotation_webpage.url,'#attributes'))  # get_attributes_output(url='https://www.marinespecies.org/aphia.php?p=search&adv=1',id=annotation_webpage.url.split('id=')[-1])
                        if annotation_webpage_with_attributes.ok:#len(list(filter(None,[id.get('id') if id.get('id') == "aphia_attributes_group_show" else id.get('href') if id.get('href') == '#attributes' else None for id in taxo_soup.findAll('a')]))) > 0:  # taxo_soup.findAll('a',onclick=True)[0].get('id')=="aphia_attributes_group_show":
                            attributes_soup = BeautifulSoup(annotation_webpage_with_attributes.content, 'lxml')
                            script = attributes_soup.find_all('script')
                            if len(attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})):
                                if 'No attributes found on lower taxonomic level' in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText() or len([item.getText() for item in script if 'Functional group' in item.getText()])==0:
                                    attributes_soup = ''
                                    script = ''
                    else:
                        attributes_soup=''
                        script=''
                        annotation_webpage_with_attributes = session.get(urljoin(annotation_webpage.url, '#attributes'))  # get_attributes_output(url='https://www.marinespecies.org/aphia.php?p=search&adv=1',id=annotation_webpage.url.split('id=')[-1])
                        if annotation_webpage_with_attributes.ok:#len(list(filter(None,[id.get('id') if id.get('id') == "aphia_attributes_group_show" else id.get('href') if id.get('href') == '#attributes' else None for id in taxo_soup.findAll('a')]))) > 0:
                            attributes_soup=BeautifulSoup(annotation_webpage_with_attributes.content, 'lxml')
                            script = attributes_soup.find_all('script')
                            if len(attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})):
                                if 'No attributes found on lower taxonomic level' in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText() or len([item.getText() for item in script if 'Functional group' in item.getText()])==0:
                                    attributes_soup = ''
            fields = [item.getText() if len(taxo_soup.find_all(class_="alert alert-info")) == 0 else '' for item in taxo_soup.find_all(class_="col-xs-12 col-sm-4 col-lg-2 control-label")]
            Status = re.sub(r'[' + '\n' + '\xa0' + ']', '',taxo_soup.find_all(class_="leave_image_space")[fields.index('Status')].getText()) if len(taxo_soup.find_all(class_="alert alert-info")) == 0 else ''
            if Status=='' or len(re.findall(r'unaccepted', Status)) > 0: #|uncertain|unassessed
                if 'Accepted Name' not in fields:
                    continue
                else:
                    #annotation = re.sub(r'[' + '\n' + '\xa0' + ']', '', taxo_soup.find_all(class_="leave_image_space")[ fields.index('Accepted Name')].getText()) if len( taxo_soup.find_all(class_="alert alert-info")) == 0 and 'Accepted Name' in fields else annotation
                    #data['tName'] = annotation.split(' (')[0]
                    taxon_list=session.get(urljoin(url, taxo_soup.find_all(class_="leave_image_space")[ fields.index('Accepted Name')].find_all('a')[0]['href']), data=data.fromkeys(['marine','fossil']))#session.post(url=url,data=data)
                    # Transform form query into xml table and save results in dataframe
                    soup = BeautifulSoup(taxon_list.content, "lxml")
                    search_response = pd.DataFrame({'Taxon': [item.getText() for item in soup.find_all(class_="list-group-item") if any(re.findall(r'unaccepted|uncertain|unassessed',item.getText())) == False],
                                                    'Link': [item.find('a').get('href') for item in soup.find_all(class_="list-group-item") if any(re.findall(r'unaccepted|uncertain|unassessed',item.getText())) == False]})
                    if len(search_response) == 0:
                        if len(soup.find_all(class_="alert alert-info")) > 0:
                            continue
                        else:
                            taxo_soup = soup
                            attributes_soup = ''
                            script = ''
                            annotation_webpage_with_attributes = session.get(urljoin(taxon_list.url, '#attributes'))  # get_attributes_output(url='https://www.marinespecies.org/aphia.php?p=search&adv=1',id=annotation_webpage.url.split('id=')[-1])
                            if annotation_webpage_with_attributes.ok:#len(list(filter(None,[id.get('id') if id.get('id') == "aphia_attributes_group_show" else id.get('href') if id.get('href') == '#attributes' else None for id in soup.findAll('a')]))) > 0:
                                attributes_soup = BeautifulSoup(annotation_webpage_with_attributes.content, 'lxml')
                                script = attributes_soup.find_all('script')
                                if len(attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})):
                                    if 'No attributes found on lower taxonomic level' in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[ 0].getText() or len([item.getText() for item in script if 'Functional group' in item.getText()]) == 0:
                                        attributes_soup = ''
                                        script = ''
                    else:
                        annotation_webpage = session.get(urljoin(url, search_response.Link[0]), data=data)
                        taxo_soup = BeautifulSoup(annotation_webpage.content, 'lxml')
                        attributes_soup = ''
                        annotation_webpage_with_attributes = session.get(urljoin(annotation_webpage.url, '#attributes'))  # get_attributes_output(url='https://www.marinespecies.org/aphia.php?p=search&adv=1',id=annotation_webpage.url.split('id=')[-1])
                        if annotation_webpage_with_attributes.ok: #taxo_soup.findAll('a', onclick=True)[0].get('id') == "aphia_attributes_group_show":
                            attributes_soup = BeautifulSoup(annotation_webpage_with_attributes.content, 'lxml')
                            script = attributes_soup.find_all('script')
                            if len(attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})):
                                if 'No attributes found on lower taxonomic level' in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText() or len([item.getText() for item in script if 'Functional group' in item.getText()])==0:
                                    attributes_soup = ''
                    fields = [item.getText() if len(taxo_soup.find_all(class_="alert alert-info")) == 0 else '' for item in taxo_soup.find_all(class_="col-xs-12 col-sm-4 col-lg-2 control-label")]

            Functional_group =[item.getText()[substring.start():]  for item in script for substring in re.finditer('Functional group', item.getText()) if 'Functional group' in item.getText()]
            Functional_group = ';'.join([str({'Functional group':group[group.find('Functional group') + 22:[sub.start() for sub in re.finditer('&nbsp', group) if sub.start() > group.find('Functional group') + 22 and sub.start()>group.find('nodes')][0]].replace(group[[sub.start() for sub in re.finditer('&nbsp', group) if sub.start() > group.find('Functional group') + 22 ][0]:[sub.start() for sub in re.finditer('&nbsp', group) if sub.start() > group.find('nodes') + 22 ][0]],''),group.split('&nbsp;" ,state: "" ,nodes: [{ text: "<b>')[1].split('<\\/b> ')[0] :group.split('&nbsp;" ,state: "" ,nodes: [{ text: "<b>')[1].split('<\\/b> ')[1].split('&nbsp')[0]}) if group.find('&nbsp;" ,state: "" ,nodes: [{ text: "<b>')!=-1 else str({'Functional group':group[group.find('Functional group') + 22:[sub.start() for sub in re.finditer('&nbsp', group) if sub.start() > group.find('Functional group') + 22][0]]}) for group in Functional_group])
            Type = np.where(taxo_soup.find_all(class_="leave_image_space")[1].getText().split('\n')[1] == 'Biota','Living', 'NA').tolist() if len( taxo_soup.find_all(class_="alert alert-info")) == 0 else ''
            dict_hierarchy = {re.sub(r'[' + string.punctuation + ']', '', level.split('\xa0')[1]): level.split('\xa0')[0] for level in taxo_soup.find_all(class_="leave_image_space")[1].getText().split('\n') if '\xa0' in level} if len( taxo_soup.find_all(class_="alert alert-info")) == 0 else dict({'': ''})
            full_hierarchy = '>'.join([level.split('\xa0')[0] + level.split('\xa0')[1] for level in taxo_soup.find_all(class_="leave_image_space")[1].getText().split('\n') if '\xa0' in level]) if len(taxo_soup.find_all(class_="alert alert-info")) == 0 else ''
            Domain = dict_hierarchy['Kingdom'] if 'Kingdom' in dict_hierarchy.keys() else ''
            Genus = dict_hierarchy['Genus'] if 'Genus' in dict_hierarchy.keys() else ''
            Family = dict_hierarchy['Family'] if 'Family' in dict_hierarchy.keys() else ''
            Order = dict_hierarchy['Order'] if 'Order' in dict_hierarchy.keys() else ''
            Class = dict_hierarchy['Class'] if 'Class' in dict_hierarchy.keys() else ''
            Phylum = dict_hierarchy['Phylum'] if 'Phylum' in dict_hierarchy.keys() else ''
            Rank = re.sub(r'[' + '\n' + ' ' + ']', '',taxo_soup.find_all(class_="leave_image_space")[fields.index('Rank')].getText()) if len(taxo_soup.find_all(class_="alert alert-info")) == 0 else ''
            URL = re.sub(r'[' + '(' + ')' + ']', '',taxo_soup.find_all(class_="aphia_core_cursor-help")[fields.index('AphiaID')].getText()) if len(taxo_soup.find_all(class_="alert alert-info")) == 0 else ''
            Original_reference = taxo_soup.find_all(class_="correctHTML")[0].getText() if len(taxo_soup.find_all(class_="correctHTML")) > 0 else ''
            Taxonomic_citation = [re.sub(r'[' + '\n' + ' ' + ']', '', item.getText()) for item in taxo_soup.find_all(class_="col-xs-12 col-sm-8 col-lg-10 pull-left") if item.attrs['id'] == 'Citation'][0] if len(taxo_soup.find_all(class_="alert alert-info")) == 0 else ''
            df_hierarchy = pd.DataFrame({'EcoTaxa_hierarchy': hierarchy,'Full_hierarchy': full_hierarchy, 'Rank': Rank,'Type': Type,'Domain': Domain,'Phylum': Phylum, 'Class': Class,'Order': Order, 'Family': Family, 'Genus': Genus,'Functional_group': Functional_group if len(Functional_group) > 0 else '','WORMS_ID': URL,'Reference': Original_reference,'Citation': Taxonomic_citation}, index=[0])
            break
    return df_hierarchy


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

# get attributes from the World Register of Marine Species
def get_attributes_output(url,id): # get_attributes_output(url='https://www.marinespecies.org/aphia.php?p=search&adv=1',id='Trichodesmium')
    """ This functions returns the xml table of the html response after click, given the click input id"""
    attributes_soup=''
    with webdriver.Chrome(chromedriver,options=options,keep_alive=False) as driver:#(service=Service(ChromeDriverManager().install()),keep_alive=False,options=options.add_argument("--window-size=1920,1200")) as driver:
        driver.get(url) # 'https://www.marinespecies.org/aphia.php?p=search&adv=1'
        driver.find_element(by='id',value='tName').send_keys(id) # ID
        driver.find_elements(by='name', value='marine')[1].send_keys('(any)') # Marine taxa only
        driver.find_elements(by='name', value='fossil')[1].send_keys('(any)')  # Fossil taxa only
        driver.find_elements(By.XPATH, '//button')[-1].click() # Submit search queries
        driver.get(driver.current_url+'#attributes')

        if len(driver.find_elements(by='id', value='aphia_attributes_group_show')) > 0:
            driver.find_element(by='id', value='aphia_attributes_group_show').click()  # Click to show attributes children
            time.sleep(10)  # Wait 10 sec before getting the webpage and killing the session
            attributes_soup = BeautifulSoup(driver.page_source, 'lxml')
            driver.close()

    return attributes_soup

# Return results (Rank, Domain, Phylum, Class, Order, Family, Genus, Functional group, URL, Reference, Citation) of the form
def get_form_hierarchy(hierarchy,url,form):
 """
 Objective: This functions returns the taxonomic informations of a given hierarchy stored in the World Register of Marine Species database (url)
 Note: Functions deprecated. Use annotations_in_WORMS function instead.
 :param hierarchy: Full hierarchy to be searched on url. Use > to separate different levels.
 :param url: Path of the taxonomic database query (https://www.marinespecies.org/aphia.php?p=search)
 :param form: Form containing the input(s) to be posted on the url
 :return: dataframe with corresponding rank, domain, phylum, class, order, family, genus, functional group, Taxon url, reference, and citation
 """
 with HTMLSession() as session:
    hierarchy_levels = hierarchy.split('>')
    df_hierarchy = pd.DataFrame({'EcoTaxa_hierarchy': hierarchy,'Full_hierarchy': '','Rank': '','Type': '','Domain': '','Phylum': '', 'Class': '','Order':'','Family':'','Genus': '','Functional_group': '','WORMS_ID': '','Reference': '','Citation': ''}, index=[0])
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
        if len(soup.find_all(class_="list-group-item"))==0: # If post results in a single page
            taxo_soup = BeautifulSoup(res.html.html, 'lxml')
            if len(taxo_soup.find_all(class_="alert alert-warning"))>0:
                continue # Skip current level if hierarchy not found
            else:
                attributes_soup = ''
                if len([id.get('id') for id in taxo_soup.findAll('a',onclick=True) if id.get('id')=="aphia_attributes_group_show"])>0:#taxo_soup.findAll('a',onclick=True)[0].get('id')=="aphia_attributes_group_show":
                    attributes_soup = get_attributes_output(url='https://www.marinespecies.org/aphia.php?p=search&adv=1',id=res.url.split('id=')[-1])#get_click_output(url=res.url + '#attributes',id="aphia_attributes_group_show")  # BeautifulSoup(driver.page_source, 'lxml')
                    if 'No attributes found on lower taxonomic level' in  attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText() or 'Functional group' not in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText():
                        attributes_soup = ''


        else: # If post results in a list of taxa
            search_response = pd.DataFrame({'Taxon': [item.getText()  for item in soup.find_all(class_="list-group-item") if any(re.findall(r'unaccepted|uncertain|unassessed',item.getText()))==False],
                 'Link': [item.find('a').get('href')  for item in soup.find_all(class_="list-group-item") if any(re.findall(r'unaccepted|uncertain|unassessed',item.getText()))==False]})

            if len(search_response)==0:
                     continue # Skip current level if hierarchy not found
            if search_response.Taxon[0].split(' ')[0].lower() != level.lower():
                     continue # Skip current level if search response different from level
            else:
                     connection_webpage = requests.get(urljoin(url, search_response.Link[0]),data=data)
                     taxo_soup = BeautifulSoup(connection_webpage.content, 'lxml')
                     if len(taxo_soup.find_all(class_="alert alert-warning"))>0:
                         res = session.post(urljoin(url, search_response.Link[0]),data=data)
                         soup = BeautifulSoup(res.html.html, "lxml")
                         taxo_soup = BeautifulSoup(res.html.html, 'lxml')
                         attributes_soup = ''
                         if len([id.get('id') for id in taxo_soup.findAll('a',onclick=True) if id.get('id')=="aphia_attributes_group_show"])>0:#taxo_soup.findAll('a',onclick=True)[0].get('id')=="aphia_attributes_group_show":
                             attributes_soup =get_attributes_output(url='https://www.marinespecies.org/aphia.php?p=search&adv=1',id=connection_webpage.url.split('id=')[-1]) #get_click_output(url=connection_webpage.url + '#attributes',id="aphia_attributes_group_show")  # BeautifulSoup(driver.page_source, 'lxml')
                             if 'No attributes found on lower taxonomic level' in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText() or 'Functional group' not in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText():
                                 attributes_soup = ''

                     else:
                         attributes_soup =''
                         if len([id.get('id') for id in taxo_soup.findAll('a',onclick=True) if id.get('id')=="aphia_attributes_group_show"])>0:
                             attributes_soup = get_attributes_output(url='https://www.marinespecies.org/aphia.php?p=search&adv=1',id=connection_webpage.url.split('id=')[-1])#get_click_output(url=connection_webpage.url + '#attributes', id="aphia_attributes_group_show")  # BeautifulSoup(driver.page_source, 'lxml')
                             if 'No attributes found on lower taxonomic level' in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText() or 'Functional group' not in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText():
                                 attributes_soup = ''


        fields=[item.getText() if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else '' for item in taxo_soup.find_all(class_="col-xs-12 col-sm-4 col-lg-2 control-label" )]
        Status=re.sub(r'[' + '\n' + '\xa0' + ']', '', taxo_soup.find_all(class_="leave_image_space")[fields.index('Status')].getText()) if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else ''
        if len(re.findall(r'unaccepted|uncertain|unassessed', Status))>0:  # Re-assign hierarchy level with accepted name if status in unaccepted
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
                   connection_webpage = requests.get(urljoin(url, search_response.Link[0]),data=data)
                   taxo_soup = BeautifulSoup(connection_webpage.content, 'lxml')
                   attributes_soup = ''
                   if taxo_soup.findAll('a',onclick=True)[0].get('id')=="aphia_attributes_group_show":
                          attributes_soup = get_attributes_output(url='https://www.marinespecies.org/aphia.php?p=search&adv=1',id=connection_webpage.url.split('id=')[-1])#get_click_output(url=connection_webpage.url + '#attributes',id="aphia_attributes_group_show")  # BeautifulSoup(driver.page_source, 'lxml')
                          if 'No attributes found on lower taxonomic level' in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText() or 'Functional group' not in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText():
                                attributes_soup = ''
                   fields = [item.getText() if len(taxo_soup.find_all(class_="alert alert-warning")) == 0 else '' for item in taxo_soup.find_all(class_="col-xs-12 col-sm-4 col-lg-2 control-label")]

        script =taxo_soup.find_all('script')
        Functional_group =[item.getText()[item.getText().find('Functional group') + 22:[sub.start() for sub in re.finditer('&nbsp', item.getText()) if sub.start() > item.getText().find('Functional group') + 22][0]] for item in script if 'Functional group' in item.getText()] if len(attributes_soup)==0 or [item.replace('Attributes assigned at lower taxonomic level','') for item in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText().split('\n') if 'Functional group' in item ][-1].replace('Functional group', '').strip()=='not applicable'  else [item.replace('Attributes assigned at lower taxonomic level','') for item in attributes_soup.findAll(attrs={'id': 'aphia_attributes_group'})[0].getText().split('\n') if 'Functional group' in item ][-1].replace('Functional group', '').strip()
        Type = np.where(taxo_soup.find_all(class_="leave_image_space")[1].getText().split('\n')[1] == 'Biota','Living', 'NA').tolist() if len(taxo_soup.find_all(class_="alert alert-warning"))==0 else ''
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
                                         'Type': Type,
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

def taxon_color_palette(dataframe, levels, palette=px.colors.qualitative.T10):
  if not all(pd.Series(levels).isin(dataframe.columns.values)):
    return print(','.join(levels) + ' not found in dataframe')
  else:
    # Create/Remove colors from the palette to match the number of higher level categories
    try:
      pal = matplotlib.colors.LinearSegmentedColormap.from_list(name='palette', colors=palette,N=len(dataframe[levels.categories[0]].unique()))
    except ValueError:
      palette = [rgb2hex(int(re.sub("[^0-9]", "", rgb.split(',')[0])), int(re.sub("[^0-9]", "", rgb.split(',')[1])),
                         int(re.sub("[^0-9]", "", rgb.split(',')[2]))) for rgb in palette]
      pal = matplotlib.colors.LinearSegmentedColormap.from_list(name='palette', colors=palette, N=len(dataframe[levels.categories[0]].unique()))

    hex_palette = sns.color_palette(pal(range(len(dataframe[levels.categories[0]].unique()))), len(dataframe[levels.categories[0]].unique())).as_hex()  # sns.color_palette(colors,len(melt_df_taxonomy.Group.unique())).as_hex()
    # sns.palplot(hex_palette)
    colors_dict = {level: hex_palette[index].upper() for index, level in enumerate(dataframe[levels.categories[0]].unique())}
    for index, level in enumerate(levels):
      if index == 0:
        continue
      else:
        for category in dataframe[levels.categories[index]].unique():
          sub_level = dataframe[dataframe[levels.categories[index]] == category][levels].drop_duplicates()[[levels.categories[index - 1]]].values[0][0]
          new_pal = matplotlib.colors.LinearSegmentedColormap.from_list(name='subpalette',colors=[colors_dict[sub_level], '#FFFFFF'],
                                                                        N=len(dataframe[dataframe[levels.categories[index - 1]] == sub_level][levels].drop_duplicates()[ levels.categories[index]].unique()) + 1)
          hex_new_palette = sns.color_palette(new_pal(range(len(dataframe[dataframe[levels.categories[index - 1]] == sub_level][levels].drop_duplicates()[levels.categories[index]].unique()))), len(dataframe[dataframe[levels.categories[index - 1]] == sub_level][levels].drop_duplicates()[levels.categories[index]].unique())).as_hex()  # sns.color_palette(colors,len(melt_df_taxonomy.Group.unique())).as_hex()
          # sns.palplot(hex_new_palette)
          colors_dict.update({taxon: hex_new_palette[i].upper() for i, taxon in enumerate(dataframe[dataframe[levels.categories[index - 1]] == sub_level][levels].drop_duplicates()[levels.categories[index]].unique())})
  return colors_dict