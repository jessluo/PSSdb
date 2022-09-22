# Objective: This script uses web crawling library to retrieve IFCB metadata from an online dashboard

from pathlib import Path
# Module for web crawling with javascript inputs
# Selenium for web scraping: https://stackoverflow.com/questions/29385156/invoking-onclick-event-with-beautifulsoup-python
# Download chromedriver at:https://chromedriver.chromium.org/downloads
# Go in System Preferences > Security and Privacy > General > Allow
from selenium import webdriver # Use: pip3 install -U selenium
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager # Use: pip3 install webdriver-manager
from selenium.webdriver.chrome.options import Options
chromedriver = '{}/Downloads/chromedriver'.format(str(Path.home())) # Local path to the chromedriver. Note: Deprecated
options=Options()
options.add_argument("--headless") # Use this to keep the browser session closed

import re
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

import time
url = "https://ifcb.caloos.org/timeline?dataset=santa-cruz-municipal-wharf&bin=D20160513T221410_IFCB104"

with webdriver.Chrome(chromedriver, options=options,keep_alive=False) as driver:
    driver.get(url)
    time.sleep(3) # Time needed for the webpage to be fully loaded, in seconds
    vol_analyzed=float(re.sub("[^\d\.]", "",driver.find_element(by='id', value='stat-ml-analyzed').text)) if len(driver.find_element(by='id', value='stat-ml-analyzed').text)>0 else pd.NA
    depth = float(re.sub("[^\d\.]", "", driver.find_element(by='id', value='stat-depth').text)) if len(driver.find_element(by='id', value='stat-depth').text)>0 else pd.NA
    Latitude = float(driver.find_element(by='id', value='stat-lat').text) if len(driver.find_element(by='id', value='stat-lat').text)>0 else pd.NA
    Longitude = float(driver.find_element(by='id', value='stat-lon').text) if len(driver.find_element(by='id', value='stat-lon').text)>0 else pd.NA
    driver.quit()
