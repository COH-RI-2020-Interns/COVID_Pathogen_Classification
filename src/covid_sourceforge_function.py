"""
Author: Rishov S. Chatterjee

This script is designed for scraping the various test datasets used by Randhawa et. al for classifying the taxonomy
of a viral genome sequence to be made available.

Modules used:
- re
- time
- requests
- BeautifulSoup
- os

All the code in this script is the intellectual property of City of Hope National Medical Center.

"""
import re
import time
import requests
from bs4 import BeautifulSoup
from os import getcwd, listdir, system, chdir

def fasta_scraper(test_name: str, family_name: str):

    # Creating the Test Name as a Folder
    fasta_path = getcwd() + "/data"

    if test_name not in listdir(fasta_path):
        chdir("data")
        system(f"mkdir {test_name}")
        chdir("..")

    # Getting the COVID19 Data on Web Page from Source Forge
    web_page = requests.get(f"https://sourceforge.net/projects/mldsp-gui/files/COVID19Dataset/{test_name}/{family_name}/")

    # Parsing the HTML Content of the HTTP GET Request Using BeautifulSoup's HTML Parser
    soup = BeautifulSoup(web_page.content, 'html.parser')

    # Finding all the Table Headers with the special header: files_name_h
    fasta_files = soup.find_all('th', headers="files_name_h")

    # Extracting all the fasta download links using a list comprehension
    links = [list(link.children)[0].attrs['href'] for link in fasta_files]

    # Creating a Folder to Store All Files Pertaining to Virus Family
    chdir(f"data/{test_name}")
    system(f"mkdir {family_name}")
    chdir("../..")

    # Getting the Starting Index of Each / to Pinpoint File Name
    slash_routes = [m.start() for m in re.finditer("/", links[0])]

    # Downloading All the Files to the data folder
    for i, link in enumerate(links):
        file_name = link[slash_routes[-2]+1:slash_routes[-1]]
        get_req = requests.get(link, allow_redirects=True)
        open(f"{fasta_path}/{test_name}/{family_name}/{file_name}", 'wb').write(get_req.content)
        print(f"{test_name}/{family_name} Fasta Sequence {i+1}: {file_name} successfully downloaded.")
        time.sleep(5)
    return "Thanks for using the scraper!"
