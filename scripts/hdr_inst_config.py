# GOAL: to extract important instrument configuration that determines the 'sampling' restrictions

import pyifcb as ifcb
url = 'http://ifcb-data.whoi.edu/IFCB101_BigelowJun2016/D20160709T210837_IFCB101.html'

#with ifcb.open_url(url) as sample_bin:
    #lid = sample_bin.lid
    print('{} has {} image(s)'.format(lid, len(sample_bin.images)))