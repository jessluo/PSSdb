
import requests
from bs4 import BeautifulSoup

url = "https://ifcb.caloos.org/timeline?dataset=santa-cruz-municipal-wharf&bin=D20160513T221410_IFCB104"
page = requests.get(url)

soup = BeautifulSoup(page.content, "html.parser")

vol_analyzed = soup.find('span', {'id': 'stat-ml-analyzed'}).text
depth = soup.find('span', {'id': 'stat-depth'}).text
Latitude = soup.find('span', {'id': 'stat-lat'}).text
Longitude = soup.find('span', {'id': 'stat-lon'}).text