#process to generate an IFCB dashboard project list

import requests
import re
import pandas as pd




url = 'https://ifcb-data.whoi.edu/'

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

#def matchingKeys(dictionary, searchString):
keys_proj_numbers = [key for key,val in data_dict.items() if any('timeline?dataset=' in s for s in val)]


proj_index = []
proj_name = []

for n in keys_proj_numbers:
    for i in data_dict[n]:
        proj_info = re.search('timeline\?dataset=(.*)', i)
        try:
            l = proj_info.group(1).split('>')
            l[0] = l[0].replace('\'', '').replace('\"', '')
            l[1] = l[1].replace('</a', '')
            print(l[0])
            print(l[1])
            if l[0] not in proj_index:
                proj_index.append(l[0])
                proj_name.append(l[1])
        except:
            pass







metadata_bins = pd.DataFrame({'Variables':['Biovolume',  'date_bin', 'Station_location', 'midLatBin', 'midLonBin', 'midDepthBin'],
                                 'Variable_types': ['float64',  'int64', 'object', 'float64', 'float64', 'float64'],
                                 'Units/Values/Timezone': ['cubic_micrometer',  date_group, 'lat_lon', 'degree', 'degree', 'meters'],
                                 'Description': ['Biovolume calculated as a spherical projection of ' + area_type ,
                                                 'binned date information',
                                                 'string that serves as an identifier of a single cell of a  1x1 degree spatial grid',
                                                 'latitude of the center point of the 1x1 degree cell',
                                                 'longitude of the center point of the 1x1 degree cell',
                                                 'middle point within a depth bin']})
