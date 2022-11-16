#process to generate an IFCB dashboard project list

import requests
import re
import os
import yaml
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from funcs_IFCB_dashboard import *

# set path
# set paths:
path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
# open the metadata of the standardized files
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)
path_download = str(Path(cfg['raw_dir']).expanduser())

path_save = path_download + '/IFCB_dashboard_projects_list.xlsx'

if os.path.exists(path_save):
    exists = input('Warning! the IFCB_dashboard_projects_list exists. Do you want to update? \n counting the number of ROIs in a time series is very time consuming \n y/n')
    if exists == 'y':
        #os.remove(path_save)
        # predefined information of dashboards:
        # dashboards = ['CALOOS', 'WHOI']
        contacts = ['Clarissa Anderson', 'Heidi Sozik']
        emails = ['cra002@ucsd.edu', 'hsosik@whoi.edu']
        urls = ['https://ifcb.caloos.org/', 'https://ifcb-data.whoi.edu/']

        ##test projects
        test = ['santa-cruz-municipal-wharf', 'NAAMES']

        # lists to be filled with the project
        Project_ID = []
        Project_title = []
        dashboard_url = []
        Instrument = []
        Contact_name = []
        Contact_email = []
        Pssdb_access = []
        Objects_number = []
        Percentage_classified = []
        Percentage_validated = []
        Latest_update = []
        Project_test = []
        n_roi = 0

        for ds, url in enumerate(urls):
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
            keys_proj_numbers = [key for key, val in data_dict.items() if any('timeline?dataset=' in s for s in val)]
            for n in keys_proj_numbers:
                for i in data_dict[n]:
                    proj_info = re.search('timeline\?dataset=(.*)', i)
                    try:
                        l = proj_info.group(1).split('>')
                        l[0] = l[0].replace('\'', '').replace('\"', '')
                        l[1] = l[1].replace('</a', '')
                        if l[0] not in Project_ID:
                            Project_ID.append(l[0])
                            Project_title.append(l[1])
                            dashboard_url.append(urls[ds])
                            Instrument.append('IFCB')
                            Contact_name.append(contacts[ds])
                            Contact_email.append(emails[ds])
                            Pssdb_access.append('True')
                            Percentage_classified.append(100)
                            Percentage_validated.append(float('nan'))
                            if (l[0] == 'stearns-wharf' or l[0] == 'Llopiz_preserved_microzooplankton'):
                                Project_test.append('True')
                            else:
                                Project_test.append('False')
                            # step here to get ROI count for each projects:
                            try:
                                pid_df = get_df_list_IFCB(urls[ds], l[0])
                                print('counting the number of ROIs for '+ l[0])
                                for i in tqdm(pid_df.loc[:, 'pid']):
                                    n_roi = n_roi + roi_number(urls[ds], i)
                                Objects_number.append(n_roi)
                                Latest_update.append(pid_df.loc[len(pid_df) - 1, 'sample_time'])
                            except:
                                Objects_number.append(float('nan'))
                                Latest_update.append(float('nan'))

                    except:
                        pass

        projects_df = pd.DataFrame()

        projects_df['Project_ID'] = Project_ID
        projects_df['Project_title'] = Project_title
        projects_df['dashboard_url'] = dashboard_url
        projects_df['Instrument'] = Instrument
        projects_df['Contact_name'] = Contact_name
        projects_df['Contact_email'] = Contact_email
        projects_df['Pssdb_access'] = Pssdb_access
        projects_df['Objects_number'] = Objects_number
        projects_df['Percentage_classified'] = Percentage_classified
        projects_df['Percentage_validated'] = Percentage_validated
        projects_df['Latest_update'] = Latest_update
        projects_df['Project_test'] = Project_test

        projects_df = projects_df[projects_df['Latest_update'].notna()]
        os.remove(path_save)
        #projects_df.to_csv(path_save, sep='\t')
        with pd.ExcelWriter(path_save, engine="xlsxwriter") as writer: #saving as excel file not working
            projects_df.to_excel(writer, sheet_name='Data', index=False)
        print('IFCB_dashboard_projects_list.xlsx file saved in ' + path_download)

else:
    print('process finished')



