# check out this example for lat and lon

# https://www.oreilly.com/library/view/python-data-science/9781491912126/ch04.html

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

from pathlib import Path
from glob import glob
# Config modules
import yaml  # requires installation of PyYAML package


#look for coastline layer and add it to that map, UNABLE to import cartopy, see issue
import cartopy.crs as ccrs


path_to_config = Path('~/GIT/PSSdb/scripts/Ecotaxa_API.yaml').expanduser()
with open(path_to_config, 'r') as config_file:
    cfg = yaml.safe_load(config_file)

# find path to get the data with slopes and intercepts

path_to_data = Path(cfg['git_dir']).expanduser() / cfg['dataset_subdir']

filename = path_to_data / 'NPSS_IFCB.tsv'

df = pd.read_csv(filename, sep='\t', encoding_errors='ignore', header=0)
df = df.dropna()
# Extract the data we're interested in
lat, lon = df['station_lat'], df['station_lon']
slope, intercept = df['slope'], df['intercept']

intercept_t = [a * 10 for a in intercept]

# Scatter the points, using size and color but no label

#plt.axis(aspect='equal')

#plt.figure(figsize=(5, 3))
#set lat lon grid
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([min(lon)-5, max(lon)+5, min(lat)-5, max(lat)+5] )
plt.gca().coastlines('50m')
g1=ax.gridlines(draw_labels=True)
g1.xlines = False
g1.ylines = False
ax.set_title(''.join("'PSS results for IFCB  projects 3315, 3318, 3326'"))

plt.scatter(lon, lat, label=None, c=slope, cmap='seismic', s=intercept_t, linewidth=0, alpha=0.5, transform=ccrs.PlateCarree())
#ax.xlabel('longitude')
#ax.ylabel('latitude')
plt.colorbar(label='slope', orientation='horizontal', anchor=(0.5, -1))
#ax.clim(min(slope), max(slope))


labels = ["3", "8", "12.5"]

for n, area in enumerate([30, 80, 125]):
    plt.scatter([], [], c='k', alpha=0.3, s=area, label=labels[n], transform=ccrs.PlateCarree())

plt.legend( ncol = 3, scatterpoints=1, frameon=False,
           labelspacing=1, title='intercept')


figname = 'slopes_intercept_IFCB.pdf'

savepath = Path(cfg['git_dir']).expanduser() / 'figures' / figname

plt.savefig(savepath)
