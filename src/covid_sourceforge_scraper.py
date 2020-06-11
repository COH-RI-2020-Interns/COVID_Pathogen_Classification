"""
Author: Rishov S. Chatterjee

This Python script is designed for scraping the various test datasets used by Randhawa et. al for classifying the taxonomy
of a viral genome sequence to be made available.

Modules Used:
- requests
- BeautifulSoup

Modular Functions Used:
- src/covid_sourceforge_function.py -> fasta_scraper()

All the code in this script is the intellectual property of City of Hope National Medical Center.

"""
import requests
from bs4 import BeautifulSoup
from os import chdir
chdir("src")
from covid_sourceforge_function import fasta_scraper
chdir("..")

# Get Randhawa's Main Web Page on Source Forge
url = "https://sourceforge.net/projects/mldsp-gui/files/COVID19Dataset/"
web_page = requests.get(url)

# Get the HTML Elements of the Page
soup = BeautifulSoup(web_page.content, 'html.parser')

# Get List of Test Names
test_names = [soup.find_all('tr')[i]['title'] for i in range(2,9)]
print(test_names)

test_name = input("Please type in a test name from the above choices: ")

# Get Virus Families in Test Page
try:
    virus_page = requests.get(f"{url}{test_name}")
except TypeError:
    print("Unable to retrieve web page. Unknown test name.")

# Get HTML Elements of New Request
virus_soup = BeautifulSoup(virus_page.content, 'html.parser')

# Get List of Virus Families in the Requested Page
folders = virus_soup.find_all('tr', {'class':'folder'})
virus_names = [virus_soup.find_all('tr')[i]['title'] for i in range(2, len(folders)+1)]
print(virus_names)

virus_name = input("Please type in a viral family name from the above choices: ")

# Run the Scraper
try:
    fasta_scraper(test_name, virus_name)
except (TypeError, RuntimeError):
    print("Not able to find files to download. Try typing in the virus name again.")
