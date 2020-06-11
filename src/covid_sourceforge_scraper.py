"""
Author: Rishov S. Chatterjee

This script is designed for scraping the various test datasets used by Randhawa et. al for classifying the taxonomy
of a viral genome sequence to be made available.

Modules used:
- re
- requests
- BeautifulSoup
- os

All the code in this script is the intellectual property of City of Hope National Medical Center.

"""
import re
import requests
from bs4 import BeautifulSoup
from os import getcwd, system

# Getting the COVID19 Data on Web Page from Source Forge
web_page = requests.get("https://sourceforge.net/projects/mldsp-gui/files/COVID19Dataset/Test1/Riboviria/")

# Parsing the HTML Content of the HTTP GET Request Using BeautifulSoup's HTML Parser
soup = BeautifulSoup(web_page.content, 'html.parser')

# Finding all the Table Headers with the special header: files_name_h
fasta_files = soup.find_all('th', headers="files_name_h")

# Extracting all the fasta download links using a list comprehension
links = [list(link.children)[0].attrs['href'] for link in fasta_files]

# Data Path
fasta_path = getcwd() + "/data"

# Getting the Starting Index of Each / to Pinpoint File Name
print([m.start() for m in re.finditer("/", links[0])])

# Downloading All the Files to the data folder
for link in links:
    file_name = link[80:100]
    get_req = requests.get(link, allow_redirects=True)
    open(f"{fasta_path}/{file_name}", 'wb').write(get_req.content)
    
